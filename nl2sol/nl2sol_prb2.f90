program main

!*****************************************************************************80
!
!! MAIN is the main program for NL2SOL_PRB2.
!
!  Discussion:
!
!    NL2SOL_PRB2 runs NL2SOL on various test problems.
!
!    This main program calls NLTEST to run NL2SOL, the nonlinear
!    least-squares solver, on various test problems.
!
!    The test problems used are from the references.  Some additional
!    test problems were suggested by Jorge More.  Calls passing
!    these problems to NLTEST have been commented out (since
!    there are enough other problems), but not removed, since
!    they may be of interest to other researchers.
!
!  Modified:
!
!    02 April 2006
!
!  Reference:
!
!    K M Brown,
!    A Quadratically Convergent Newton-like Method Based upon
!    Gaussian Elimination,
!    SIAM Journal on Numerical Analysis,
!    Volume 6, pages 560-569, 1969.
!
!    John Dennis, David Gay and Roy Welsch,
!    Algorithm 573: An Adaptive Nonlinear Least-squares Algorithm,
!    ACM Transactions on Mathematical Software,
!    Volume 7, Number 3, pages 348-368, 1981.
!
!    Philip Gill and Walter Murray,
!    Algorithms for the Solution of the Non-linear Least-squares Problem,
!    SIAM Journal on Numerical Analysis,
!    Volume 15, Number 5, pages 977-991, 1978.
!
!    R R Meyer,
!    Theoretical and Computational Aspects of Nonlinear Regression,
!    in Nonlinear Programming,
!    edited by J B Rosen, O L Mangasarian, and K Ritter,
!    pages 465-486,
!    Academic Press, New York, 1970.
!
  implicit none

  integer, parameter :: prob_max = 60

  integer i
  character irc(prob_max)
  integer is(6,prob_max)
  integer iv(100)
  integer j
  integer jac
  character jtyp(2)
  integer k
  integer mxfcsv
  integer mxitsv
  integer n
  character ( len = 8 ) name(prob_max)
  integer nout
  integer nprob
  integer pu
  real rs(5,prob_max)
  logical rstart
  real v(5000)
  character ( len = 20 ) version
  integer xscal1
  integer xscal2
!
!  common storage with nltest.
!
  common /testcm/ v, rs, nout, nprob, is, iv
  common /testch/ name, irc

  data rstart /.false./
  data jtyp(1),jtyp(2)/' ','*'/

  version = 'NL2SOL version 2.2'

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NL2SOL_PRB2:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the NL2SOL library.'
!
!  Assign default values to IV and V.
!
  call dfault ( iv, v )
  nout = 6
!
!  Non-default values may now be assigned.
!
!  Activate the next line to turn off detailed summary printing
!
  if ( .false. ) then
    iv(21) = 0
  end if

  pu = iv(21)

  nprob = 0
!
!  Rosenbrock.
!
  n = 2
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,2,1,'Rosenbrock', jac, rstart, xscal1, xscal2 )
!
!  Helix.
!
  n = 3
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,3,2,'Helix', jac, rstart, xscal1, xscal2 )
!
!  Singular.
!
  n = 4
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,4,3,'Singular', jac, rstart, xscal1, xscal2 )
!
!  Woods.
!
  n = 7
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,4,4,'Woods', jac, rstart, xscal1, xscal2 )
!
!  Zangwill.
!
  n = 3
  jac = 1
  xscal1 = 1
  xscal2 = 1
  call nltest(n,3,5,'Zangwill', jac, rstart, xscal1, xscal2 )
!
!  Engvall.
!
  n = 5
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,3,6,'Engvall', jac, rstart, xscal1, xscal2 )
!
!  Branin.
!
  n = 2
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,2,7,'Branin', jac, rstart, xscal1, xscal2 )
!
!  Beale.
!
  n = 3
  jac = 1
  xscal1 = 1
  xscal2 = 2
  call nltest(n,2,8,'Beale', jac, rstart, xscal1, xscal2 )
!
!  Cragg and Levy.
!
  jac = 1
  xscal1 = 1
  xscal2 = 2
  call nltest(5,4,9,'Cragg', jac, rstart, xscal1, xscal2 )
!
!  Box.
!
  n = 10
  jac = 1
  xscal1 = 1
  xscal2 = 2
  call nltest(n,3,10,'Box', jac, rstart, xscal1, xscal2 )
!
!  Davidon 1.
!
  n = 15
  jac = 1
  mxfcsv = iv(17)
  mxitsv = iv(18)
  iv(17) = 20
  iv(18) = 15
  xscal1 = 1
  xscal2 = 1
  call nltest(n,15,11,'Davidon1', jac, rstart, xscal1, xscal2 )
!
!  Freudenstein and Roth.
!
  n = 2
  jac = 1
  iv(17) = mxfcsv
  iv(18) = mxitsv
  xscal1 = 1
  xscal2 = 3
  call nltest(n,2,12,'Freudenstein', jac, rstart, xscal1, xscal2 )
!
!  Watson6.
!
  n = 31
  jac = 1
  xscal1 = 1
  xscal2 = 1
  call nltest(n,6,13,'Watson6', jac, rstart, xscal1, xscal2 )
!
!  Watson9.
!
  n = 31
  jac = 1
  xscal1 = 1
  xscal2 = 1
  call nltest(n,9,14,'Watson9', jac, rstart, xscal1, xscal2 )
!
!  Watson12.
!
  n = 31
  jac = 1
  xscal1 = 1
  xscal2 = 1
  call nltest(n,12,15,'Watson12', jac, rstart, xscal1, xscal2 )
!
!  Watson20.
!
  n = 31
  jac = 1
  mxfcsv = iv(17)
  iv(17) = 20
  mxitsv = iv(18)
  iv(18) = 15
  xscal1 = 1
  xscal2 = 3
  call nltest(n,20,16,'Watson20', jac, rstart, xscal1, xscal2 )
!
!  Chebyquad.
!
  n = 8
  jac = 1
  iv(17) = mxfcsv
  iv(18) = mxitsv
  xscal1 = 1
  xscal2 = 2
  call nltest(n,8,17,'Chebyquad', jac, rstart, xscal1, xscal2 )
!
!  Brown and Dennis.
!
  n = 20
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,4,18,'Brown', jac, rstart, xscal1, xscal2 )
!
!  Bard.
!
  n = 15
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,3,19,'Bard', jac, rstart, xscal1, xscal2 )
!
!  Jennrich and Sampson.
!
  n = 10
  jac = 1
  xscal1 = 1
  xscal2 = 1
  call nltest(n,2,20,'Jennrich', jac, rstart, xscal1, xscal2 )
!
!  Kowalik and Osborne.
!
  n = 11
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,4,21,'Kowalik', jac, rstart, xscal1, xscal2 )
!
!  Osborne 1.
!
  n = 33
  jac = 1
  xscal1 = 1
  xscal2 = 1
  call nltest(n,5,22,'Osborne1', jac, rstart, xscal1, xscal2 )
!
!  Osborne 2.
!
  n = 65
  jac = 1
  xscal1 = 1
  xscal2 = 2
  call nltest(n,11,23,'Osborne2', jac, rstart, xscal1, xscal2 )
!
!  Madsen.
!
  n = 3
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,2,24,'Madsen', jac, rstart, xscal1, xscal2 )
!
!  Meyer.
!
  n = 16
  jac = 1
  xscal2 = 1
  iv(17) = 400
  iv(18) = 300
  xscal1 = 1
  xscal2 = 3
  call nltest(n,3,25,'Meyer', jac, rstart, xscal1, xscal2 )
!
!  Brown5.
!
  n = 5
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,5,26,'Brown5', jac, rstart, xscal1, xscal2 )
!
!  Brown10.
!
  n = 10
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,10,27,'Brown10', jac, rstart, xscal1, xscal2 )
!
!  Brown30.
!  Need to increase V and IV dimensions for this problem?
!  Even so, initial sum of squares overflows!
!
  n = 30
  jac = 1
  xscal1 = 1
  xscal2 = 3
! call nltest(n,30,28,'Brown30', jac, rstart, xscal1, xscal2 )
!
!  Brown40.
!
  n = 40
  jac = 1
  xscal1 = 1
  xscal2 = 3
! call nltest(n,40,29,'Brown40', jac, rstart, xscal1, xscal2 )
!
!  Bard+10.
!
  n = 15
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,3,30,'Bard+10', jac, rstart, xscal1, xscal2 )
!
!  Kowalik and Osborne + 10.
!
  n = 11
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,4,31,'Kowal+10', jac, rstart, xscal1, xscal2 )
!
!  Meyer + 10.
!
  n = 16
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,3,32,'Meyer+10', jac, rstart, xscal1, xscal2 )
!
!  Watson6 + 10.
!
  n = 31
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,6,33,'Watson6+10', jac, rstart, xscal1, xscal2 )
!
!  Watson9 + 10.
!
  n = 31
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,9,34,'Watson9+10', jac, rstart, xscal1, xscal2 )
!
!  Watson12 + 10.
!
  n = 31
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,12,35,'Watson12+10', jac, rstart, xscal1, xscal2 )
!
!  Watson20 + 10.
!
  n = 31
  jac = 1
  xscal1 = 1
  xscal2 = 3
  call nltest(n,20,36,'Watson20+10', jac, rstart, xscal1, xscal2 )
!
!  Repeat Rosenbrock with finite-difference jacobian.
!
  if ( .false. ) then
  n = 2
  jac = 2
  xscal1 = 1
  xscal2 = 1
  iv(17) = 50
  iv(18) = 40
  call nltest(n,2,1,'Rosenbrock', jac, rstart, xscal1, xscal2 )
  end if
!
!  Repeat Brown with finite-difference jacobian.
!
  n = 20
  jac = 2
  xscal1 = 1
  xscal2 = 3
  v(29) = max ( 1.0E-07, v(29) )
  iv(17) = 30
  iv(18) = 20
  call nltest(n,4,18,'Brown', jac, rstart, xscal1, xscal2 )
!
!  Print summary.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( version )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Summary of test runs.'
  write ( *, '(a)' ) '  * means a finite difference jacobian was used.'
  write ( *, '(a)' ) ' '

  if (mod(k,56) == 1) then
    write(pu, 110)
  end if

 110     format( &
         48h0 problem    n   p  niter   nf   ng  iv1  x0scal,5x, &
         'final f     preldf     nreldf     reldx'/)

  do k = 1, nprob


         j = is(6,k)
         write(pu,120) jtyp(j), name(k), &
                       (is(i,k), i=1,5), irc(k), (rs(i,k), i=1,5)
 120     format(1x,a1,a,2i4,i7,2i5,3x,a1,f9.1,e13.3,3e11.3)

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NL2SOL_PRB2:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine nltest ( n, p, nex, title, jac, rstart, xscal1, xscal2 )

!*****************************************************************************80
!
!! NLTEST calls NL2SOL for a given problem.
!
!  Modified:
!
!    28 March 2006
!
!  Parameters:
!
!    Input, integer N, the number of equations or functions.
!
!    Input, integer P, the number of variables.
!
!    Input, integer NEX, the index of the problem.
!
!    Input, character ( len = * ) TITLE, the title of the problem.
!
!    Input, integer JAC.
!    1, the jacobian is available through a subroutine.
!    2, the jacobian must be approximated by finite differences.
!
!    Input, logical RSTART, ?
!
!    Input, integer XSCAL1, XSCAL2, control the number of runs,
!    and the scaling used.  The problem will be solved for
!    IRUN = XSCAL1 to XSCAL2, and the scaling used will be 10**(IRUN-1).
!
  implicit none

  integer, parameter :: p_max = 40
  integer, parameter :: prob_max = 60

  integer i
  character irc(prob_max)
  integer irun
  integer iv(100)
  integer jac
  integer n
  character name(prob_max)*8
  integer nex
  integer p
  integer, parameter :: prunit = 21
  integer pu
  character rc(10)
  integer, parameter :: reldx = 17
  real rs(5,prob_max)
  logical rstart
  logical rstrt
  real t
  external testj
  external testr
  character ( len = * ) title
  integer uiparm(1)
  real urparm(1)
  real v(5000)
  character ( len = 20 ) version
  real x(p_max)
  real x0scal
  integer xscal1
  integer xscal2

  common /testcm/ v, rs, nout, nprob, is, iv
  common /testch/ name, irc
  integer is(6,prob_max), nout, nprob

  character*2 alg(2)
  character*1 jtyp(2)

  integer f, f0, nfcall, nfcov, ngcall
  integer ngcov, niter, nreduc, preduc

  parameter (f=10, f0=13, nfcall=6, nfcov=40, ngcall=30)
  parameter ( ngcov=41, niter=31, nreduc=6, preduc=7 )

  data alg(1),alg(2)/'ol','no'/
  data jtyp(1),jtyp(2)/' ','*'/

  data rc / &
    '.', '+', 'x', 'r', 'b', 'a', 's', 'f', 'e', 'i' /

  uiparm(1) = nex
  rstrt = rstart
  version = 'NL2SOL version 2.2'

  if ( .not. rstrt ) then
    pu = iv(prunit)
    if ( pu .ne. 0 ) then
      write ( pu, '(a)' ) ' '
      write ( pu, '(a)' ) ' '
      write ( pu, 10 ) alg(jac), title, version
    end if
10  format (11h ***** nl2s,a2,12h on problem ,a,6h *****,6x,a)
  end if

  do irun = xscal1, xscal2
!
!  Initialize the solution vector X.
!
    if ( .not. rstrt ) then

      iv(1) = 12
      x0scal = 10.0E+00**( irun - 1 )

      call xinit ( p, x, nex )

      do i = 1, p
        x(i) = x0scal * x(i)
      end do

    end if

    if ( jac == 1 ) then
      call nl2sol ( n, p, x, testr, testj, iv, v, uiparm, urparm, testr )
    else if ( jac == 2 ) then
      call nl2sno ( n, p, x, testr, iv, v, uiparm, urparm, testr )
    end if

    if ( .not. rstrt .and. nprob .lt. 50 ) then
      nprob = nprob + 1
    end if

    name(nprob) = title

    is(1,nprob) = n
    is(2,nprob) = p
    is(3,nprob) = iv(niter)
    is(4,nprob) = iv(nfcall) - iv(nfcov)
    is(5,nprob) = iv(ngcall) - iv(ngcov)

    i = iv(1)
    irc(nprob) = rc(i)
    is(6,nprob) = jac
    rs(1,nprob) = x0scal
    rs(2,nprob) = v(f)
    if ( v(f0) <= 0.0E+00 ) then
      rs(3,nprob) = 1.0E+00
    else
      rs(3,nprob) = v(preduc) / v(f0)
    end if
    if ( v(f0) <= 0.0E+00 ) then
      rs(4,nprob) = 1.0E+00
    else
      rs(4,nprob) = v(nreduc) / v(f0)
    end if
    rs(5,nprob) = v(reldx)

    rstrt = .false.

    if ( nout /= 0 ) then

      if (nprob == 1) then
        write(nout,50) version
      end if

 50   format(1h1,11x,a,10x,24hnl2sol test summary.....,10x, &
      32h(* = finite-difference jacobian)/ &
      48h0 problem    n   p  niter   nf   ng  iv1  x0scal,5x, &
      39hfinal f     preldf     nreldf     reldx/)
      write(nout,60) jtyp(jac), title, &
      (is(i,nprob),i=1,5),irc(nprob),(rs(i,nprob),i=1,5)
 60   format(1x,a1,a,2i4,i7,2i5,3x,a1,f9.1,e13.3,3e11.3)

    end if

  end do

  return
end
subroutine testj ( n, p, x, nfcall, jac, uiparm, urparm, ufparm )

!*****************************************************************************80
!
!! TESTJ evaluates the jacobian matrix.
!
!  Discussion:
!
!    This routine evaluates the jacobian matrix J for the various
!    test problems listed in the references.
!
!  Reference:
!
!    Philip Gill and Walter Murray,
!    Algorithms for the Solution of the Non-linear Least-squares Problem,
!    SIAM Journal on Numerical Analysis,
!    Volume 15, Number 5, pages 977-991, 1978.
!
!    R R Meyer,
!    Theoretical and Computational Aspects of Nonlinear Regression,
!    in Nonlinear Programming,
!    edited by J B Rosen, O L Mangasarian, and K Ritter,
!    pages 465-486,
!    Academic Press, New York, 1970.
!
!    K M Brown,
!    A Quadratically Convergent Newton-like Method Based upon
!    Gaussian Elimination,
!    SIAM Journal on Numerical Analysis,
!    Volume 6, pages 560-569, 1969.
!
!  Parameters:
!
!    Input, integer N, is the number of functions, the length of R,
!    and the number of rows used in J.
!
!    Input, integer P, the number of variables.
!
!    Input, real X(P), the point at which the jacobian is to be evaluated.
!
!    Input, integer NFCALL, is the invocation count of TESTR.  It is
!    not needed by this routine.
!
!    Output, real JAC(N,P), the jacobian matrix at X.
!
!    Input, integer UIPARM(1), contains the value of NEX, the index
!    of the problem being solved.
!
!    Input, real URPARM(*), is a user parameter vector, which is
!    not needed by this routine.
!
!    Input, external UFPARM, is the name of a user-chosen function,
!    which is not needed here.
!
  implicit none

  integer n
  integer p

  real e
  real expmin
  integer i
  integer j
  real jac(n,p)
  integer k
  integer nex
  integer nfcall
  real r2
  real t
  real theta
  real ti
  real tim1
  real tip1
  real tpi
  real tpim1
  real tpip1
  real, parameter :: twopi = 2.0E+00 * 3.141592653589793E+00
  real u
  external ufparm
  real uftolg
  integer uiparm(1)
  real ukow(11)
  real urparm(*)
  real v
  real w
  real x(p)
  real z

  save expmin
  save uftolg
  save ukow

  data expmin / 0.0E+00 /
  data uftolg / 0.0E+00 /
  data ukow / &
    4.0E+00,   2.0E+00,   1.0E+00, 0.5E+00,    0.25E+00, &
    0.167E+00, 0.125E+00, 0.1E+00, 0.0833E+00, 0.0714E+00, &
    0.0625E+00 /

  nex = uiparm(1)

  do j = 1, p
    do i = 1, n
      jac(i,j) = 0.0E+00
    end do
  end do
!
!  Rosenbrock.
!
  if ( nex == 1 ) then

    jac(1,1) = -20.0E+00 * x(1)
    jac(1,2) = 10.0E+00
    jac(2,1) = -1.0E+00
    jac(2,2) = 0.0E+00
!
!  Helix.
!
  else if ( nex == 2 ) then

    t = x(1)**2 + x(2)**2
    ti = 100.0E+00 / ( twopi * t )
    jac(1,1) = ti * x(2)
    t = 10.0E+00 / sqrt ( t )
    jac(2,1) = x(1) * t
    jac(3,1) = 0.0E+00
    jac(1,2) = -ti * x(1)
    jac(2,2) = x(2) * t
    jac(3,2) = 0.0E+00
    jac(1,3) = 10.0E+00
    jac(2,3) = 0.0E+00
    jac(3,3) = 1.0E+00
!
!  Singular.
!
  else if ( nex == 3 ) then

    jac(1,1) = 1.0E+00
    jac(1,2) = 10.0E+00
    jac(2,3) = sqrt ( 5.0E+00 )
    jac(2,4) = -jac(2,3)
    jac(3,2) = 2.0E+00 * ( x(2) - 2.0 * x(3) )
    jac(3,3) = -2.0E+00 * jac(3,2)
    jac(4,1) = sqrt( 40.0E+00 ) * ( x(1) - x(4) )
    jac(4,4) = -jac(4,1)
!
!  Woods.
!
  else if ( nex == 4 ) then

    jac(1,1) = -20.0E+00 * x(1)
    jac(1,2) = 10.0E+00
    jac(2,1) = -1.0E+00
    jac(3,4) = sqrt ( 90.0E+00 )
    jac(3,3) = -2.0E+00 * x(3) * jac(3,4)
    jac(4,3) = -1.0E+00
    jac(5,2) = sqrt ( 9.9E+00 )
    jac(5,4) = jac(5,2)
    jac(6,2) = sqrt ( 0.2 )
    jac(7,4) = jac(6,2)
!
!  Zangwill.
!
  else if ( nex == 5 ) then

    do k = 1, 3
      do i = 1, 3
        jac(i,k) = 1.0E+00
      end do
    end do
    jac(1,2) = -1.0E+00
    jac(2,1) = -1.0E+00
    jac(3,3) = -1.0E+00
!
!  Engvall.
!
  else if ( nex == 6 ) then

    jac(1,1) = 2.0E+00 * x(1)
    jac(1,2) = 2.0E+00 * x(2)
    jac(1,3) = 2.0E+00 * x(3)
    jac(2,1) = jac(1,1)
    jac(2,2) = jac(1,2)
    jac(2,3) = 2.0E+00 * ( x(3) - 2.0E+00 )
    jac(3,1) = 1.0E+00
    jac(3,2) = 1.0E+00
    jac(3,3) = 1.0E+00
    jac(4,1) = 1.0E+00
    jac(4,2) = 1.0E+00
    jac(4,3) = -1.0E+00
    t = 2.0E+00 * ( 5.0E+00 * x(3) - x(1) + 1.0E+00 )
    jac(5,1) = 3.0E+00 * x(1)**2 - t
    jac(5,2) = 6.0E+00 * x(2)
    jac(5,3) = 5.0E+00 * t
!
!  Branin.
!
  else if ( nex == 7 ) then

    jac(1,1) = 4.0E+00
    jac(1,2) = 4.0E+00
    jac(2,1) = 3.0E+00 + (x(1) - 2.0E+00 ) * ( 3.0E+00 * x(1) - 2.0E+00 * x(2) - 2.0E+00 ) + &
       x(2) * x(2)
    jac(2,2) = 1.0E+00 + 2.0E+00 * ( 2.0E+00 * x(1) - x(2) * x(2) ) - ( x(1) - x(2) )**2
!
!  Beale.
!
  else if ( nex == 8 ) then

    jac(1,1) = x(2) - 1.0E+00
    jac(1,2) = x(1)
    jac(2,1) = x(2)**2 - 1.0E+00
    jac(2,2) = 2.0E+00 * x(1) * x(2)
    jac(3,1) = x(2)**3 - 1.0E+00
    jac(3,2) = 3.0E+00 * x(1) * x(2)**2
!
!  Cragg and Levy.
!
  else if ( nex == 9 ) then

    t = exp ( x(1) )
    jac(1,2) = -2.0E+00 * ( t - x(2) )
    jac(1,1) = -t * jac(1,2)
    jac(2,2) = 30.0E+00 * ( x(2) - x(3) )**2
    jac(2,3) = -jac(2,2)
    jac(3,3) = 2.0E+00 * sin ( x(3) - x(4) ) / ( cos ( x(3) - x(4) ) )**3
    jac(3,4) = -jac(3,3)
    jac(4,1) = 4.0E+00 * x(1)**3
    jac(5,4) = 1.0E+00
!
!  Box.
!
  else if ( nex == 10 ) then

    if ( expmin == 0.0E+00 ) then
      expmin = 1.999E+00 * log ( tiny ( expmin ) )
    end if

    do i = 1, 10

      ti = - 0.1E+00 * real ( i, kind = 4 )

      t = x(1) * ti

      if ( t < expmin ) then
        e = 0.0E+00
      else
        e = exp ( t )
      end if

      jac(i,1) = ti * e

      t = x(2) * ti
      if ( t < expmin ) then
        e = 0.0E+00
      else
        e = exp ( t )
      end if

      jac(i,2) = -ti * e
      jac(i,3) = exp ( 10.0E+00 * ti ) - exp ( ti )

    end do
!
!  Davidon 1.
!
  else if ( nex == 11 ) then

    do i = 1, n-1
      ti = real ( i, kind = 4 )
      t = 1.0E+00
      do k = 1, p
        jac(i,k) = t
        t = t * ti
      end do
    end do

    jac(n,1) = 1.0E+00
    do k = 2, p
      jac(n,k) = 0.0E+00
    end do
!
!  Freudenstein and Roth.
!
  else if ( nex == 12 ) then

    jac(1,1) = 1.0E+00
    jac(1,2) = -2.0E+00 + x(2) * ( 10.0E+00 - 3.0E+00 * x(2) )
    jac(2,1) = 1.0E+00
    jac(2,2) = -14.0E+00 + x(2) * ( 2.0E+00 + 3.0E+00 * x(2) )
!
!  Watson6.
!  Watson9.
!  Watson12.
!  Watson20.
!
  else if ( nex == 13 .or. &
            nex == 14 .or. &
            nex == 15 .or. &
            nex == 16 ) then

    do i = 1, 29
      ti = real ( i, kind = 4 ) / 29.0E+00
      r2 = x(1)
      t = 1.0E+00
      do k = 2, p
        t = t * ti
        r2 = r2 + t * x(k)
      end do
      r2 = -2.0E+00 * r2
      jac(i,1) = r2
      t = 1.0E+00
      r2 = ti * r2
      do k = 2, p
        jac(i,k) = t * ( real ( k - 1, kind = 4 ) + r2 )
        t = t * ti
      end do
    end do

    jac(30,1) = 1.0E+00
    jac(31,1) = -2.0E+00 * x(1)
    jac(31,2) = 1.0E+00
!
!  Chebyquad.
!
  else if ( nex == 17 ) then

    do k = 1, n
      tim1 = -1.0E+00 / real ( n, kind = 4 )
      z = 2.0E+00 * x(k) - 1.0E+00
      ti = z * tim1
      tpim1 = 0.0E+00
      tpi = 2.0E+00 * tim1
      z = z + z
      do i = 1, n
        jac(i,k) = tpi
        tpip1 = 4.0E+00 * ti + z * tpi - tpim1
        tpim1 = tpi
        tpi = tpip1
        tip1 = z * ti - tim1
        tim1 = ti
        ti = tip1
      end do
    end do
!
!  Brown and Dennis.
!
  else if ( nex == 18 ) then

    do i = 1, n
      ti = 0.2E+00 * real ( i, kind = 4 )
      jac(i,1) = 2.0E+00 * ( x(1) + x(2) * ti - exp ( ti ) )
      jac(i,2) = ti * jac(i,1)
      t = sin ( ti )
      jac(i,3) = 2.0E+00 * ( x(3) + x(4) * t - cos ( ti ) )
      jac(i,4) = t * jac(i,3)
    end do
!
!  Bard.
!
  else if ( nex == 19 ) then

    do i = 1, 15
      jac(i,1) = -1.0E+00
      u = real ( i, kind = 4 )
      v = 16.0 - u
      w = min ( u, v )
      t = u / ( x(2) * v + x(3) * w )**2
      jac(i,2) = v * t
      jac(i,3) = w * t
    end do
!
!  Jennrich and Sampson.
!
  else if ( nex == 20 ) then

    do i = 1, 10
      ti = real ( i, kind = 4 )
      jac(i,1) = -ti * exp ( ti * x(1) )
      jac(i,2) = -ti * exp ( ti * x(2) )
    end do
!
!  Kowalik and Osborne.
!
  else if ( nex == 21 ) then

    do i = 1, 11
      t = -1.0E+00 / ( ukow(i)**2 + x(3) * ukow(i) + x(4) )
      jac(i,1) = t * ( ukow(i)**2 + x(2) * ukow(i) )
      jac(i,2) = x(1) * ukow(i) * t
      t = t * jac(i,1) * x(1)
      jac(i,3) = ukow(i) * t
      jac(i,4) = t
    end do
!
!  Osborne 1.
!
  else if ( nex == 22 ) then

    do i = 1, 33
      ti = 10.0E+00 * real ( 1 - i, kind = 4 )
      jac(i,1) = -1.0E+00
      jac(i,2) = -exp ( x(4) * ti )
      jac(i,3) = -exp ( x(5) * ti )
      jac(i,4) = ti * x(2) * jac(i,2)
      jac(i,5) = ti * x(3) * jac(i,3)
    end do
!
!  Osborne 2.
!
!  UFTOLG is a machine-dependent constant.  It is just slightly
!  larger than the log of the smallest positive machine number.
!
  else if ( nex == 23 ) then

    if ( uftolg == 0.0E+00 ) then
      uftolg = 1.999E+00 * log ( tiny ( uftolg ) )
    end if

    do i = 1, 65
      ti = real ( 1 - i, kind = 4 ) * 0.1E+00
      jac(i,1) = -exp ( x(5) * ti )
      jac(i,5) = x(1) * ti * jac(i,1)
      do k = 2, 4
        t = x(k + 7) + ti
        theta = -x(k+4) * t * t
        if ( theta <= uftolg ) then
          r2 = 0.0E+00
        else
          r2 = -exp ( theta )
        end if
        jac(i,k) = r2
        r2 = -t * r2 * x(k)
        jac(i,k+4) = r2 * t
        jac(i,k+7) = 2.0E+00 * x(k+4) * r2
      end do
    end do
!
!  Madsen.
!
  else if ( nex == 24 ) then

    jac(1,1) = 2.0E+00 * x(1) + x(2)
    jac(1,2) = 2.0E+00 * x(2) + x(1)
    jac(2,1) = cos ( x(1) )
    jac(2,2) = 0.0E+00
    jac(3,1) = 0.0E+00
    jac(3,2) = -sin ( x(2) )
!
!  Meyer.
!
  else if ( nex == 25 ) then

    do i = 1, 16
      ti = real ( 5 * i + 45, kind = 4 )
      u = ti + x(3)
      t = exp ( x(2) / u )
      jac(i,1) = t
      jac(i,2) = x(1) * t / u
      jac(i,3) = -x(1) * x(2) * t / ( u * u )
    end do
!
!  Brown5.
!  Brown10.
!  Brown30.
!  Brown40.
!
  else if ( nex == 26 .or. &
            nex == 27 .or. &
            nex == 28 .or. &
            nex == 29 ) then

    do k = 1, n
      do i = 1, n-1
        if ( i == k ) then
          jac(i,k) = 2.0E+00
        else
          jac(i,k) = 1.0E+00
        end if
      end do
    end do

    do k = 1, n
      t = 1.0
      do i = 1, n
        if ( i /= k ) then
          t = t * x(i)
        end if
      end do
      jac(n,k) = t
    end do
!
!  Bard + 10.
!
  else if ( nex == 30 ) then

    do i = 1, 15
      jac(i,1) = -1.0E+00
      u = real ( i, kind = 4 )
      v = 16.0E+00 - u
      w = min ( u, v )
      t = u / ( x(2) * v + x(3) * w )**2
      jac(i,2) = v * t
      jac(i,3) = w * t
    end do
!
!  Kowalik and Osborne + 10.
!
  else if ( nex == 31 ) then

    do i = 1, 11
      t = -1.0E+00 / ( ukow(i)**2 + x(3) * ukow(i) + x(4) )
      jac(i,1) = t * ( ukow(i)**2 + x(2) * ukow(i) )
      jac(i,2) = x(1) * ukow(i) * t
      t = t * jac(i,1) * x(1)
      jac(i,3) = ukow(i) * t
      jac(i,4) = t
    end do
!
!  Meyer + 10.
!
  else if ( nex == 32 ) then

    do i = 1, 16
      ti = real ( 5 * i + 45, kind = 4 )
      u = ti + x(3)
      t = exp ( x(2) / u )
      jac(i,1) = t
      jac(i,2) = x(1) * t / u
      jac(i,3) = -x(1) * x(2) * t / ( u * u )
    end do
!
!  Watson6 + 10.
!  Watson9 + 10.
!  Watson12 + 10.
!  Watson20 + 10.
!
  else if ( nex == 33 .or. &
            nex == 34 .or. &
            nex == 35 .or. &
            nex == 36 ) then

    do i = 1, 29
      ti = real ( i, kind = 4 ) / 29.0E+00
      r2 = x(1)
      t = 1.0E+00
      do k = 2, p
        t = t * ti
        r2 = r2 + t * x(k)
      end do
      r2 = -2.0E+00 * r2
      jac(i,1) = r2
      t = 1.0
      r2 = ti * r2
      do k = 2, p
        jac(i,k) = t * ( real ( k - 1, kind = 4 ) + r2 )
        t = t * ti
      end do
    end do

    jac(30,1) = 1.0E+00
    jac(31,1) = -2.0E+00 * x(1)
    jac(31,2) = 1.0E+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TESTJ - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem index NEX = ', nex
    stop

  end if

  return
end
subroutine testr ( n, p, x, nfcall, r, uiparm, urparm, ufparm )

!*****************************************************************************80
!
!! TESTR evaluates the residual function.
!
!  Discussion:
!
!    This routine evaluates the residual vector R for the various
!    test functions in the references, as well as for some variations
!    suggested by Jorge More in a private communication
!    on some of these test problems, for 30 <= NEX.
!
!  Modified:
!
!    28 March 2006
!
!  Reference:
!
!    Philip Gill and Walter Murray,
!    Algorithms for the Solution of the Non-linear Least-squares Problem,
!    SIAM Journal on Numerical Analysis,
!    Volume 15, Number 5, pages 977-991, 1978.
!
!    R R Meyer,
!    Theoretical and Computational Aspects of Nonlinear Regression,
!    in Nonlinear Programming,
!    edited by J B Rosen, O L Mangasarian, and K Ritter,
!    pages 465-486,
!    Academic Press, New York, 1970.
!
!    K M Brown,
!    A Quadratically Convergent Newton-like Method Based upon
!    Gaussian Elimination,
!    SIAM Journal on Numerical Analysis,
!    Volume 6, pages 560-569, 1969.
!
!  Parameters:
!
!    Input, integer N, the number of equations or functions.
!
!    Input, integer P, the number of variables.
!
!    Input, real X(P), the point at which the residual vector
!    is to be evaluated.
!
!    Input/output, integer NFCALL; on input, the invocation count
!    of this routine.  In exceptional cases, NFCALL may be reset to
!    -1 on output to indicate an error occurred which prevented the
!    evaluation of R.
!
!    Input, integer UIPARM(1), contains the value of NEX, the index
!    of the problem being solved.
!
!    Input, real URPARM(*), is a user parameter vector, which is
!    not needed by this routine.
!
!    Input, external UFPARM, is the name of a user-chosen function,
!    which is not needed here.
!
!    Output, real R(N), the residual vector at X.
!
  implicit none

  integer n
  integer p

  real e1
  real e2
  real expmax
  real expmin
  integer i
  integer j
  integer nex
  integer nfcall
  real r(n)
  real r1
  real r2
  real ri
  real t
  real t1
  real t2
  real theta
  real ti
  real tim1
  real tip1
  real, parameter :: twopi = 2.0E+00 * 3.141592653589793E+00
  real u
  external ufparm
  real uftolg
  integer uiparm(1)
  real ukow(11)
  real urparm(*)
  real v
  real w
  real x(p)
  real ybard(15)
  real ykow(11)
  real ymeyer(16)
  real yosb1(33)
  real yosb2(65)
  real z

  save expmax
  save expmin
  save uftolg
  save ukow
  save ybard
  save ykow
  save ymeyer
  save yosb1
  save yosb2

  data expmax / 0.0E+00 /
  data expmin / 0.0E+00 /
  data uftolg / 0.0E+00 /
  data ukow / &
    4.0E+00,   2.0E+00,   1.0E+00, 0.5E+00,    0.25E+00, &
    0.167E+00, 0.125E+00, 0.1E+00, 0.0833E+00, 0.0714E+00, &
    0.0625E+00 /
  data ybard / &
    0.14E+00, 0.18E+00, 0.22E+00, 0.25E+00, 0.29E+00, &
    0.32E+00, 0.35E+00, 0.39E+00, 0.37E+00, 0.58E+00, &
    0.73E+00, 0.96E+00, 1.34E+00, 2.10E+00, 4.39E+00 /
  data ykow / &
    1.957E-01, 1.947E-01, 1.735E-01, 1.600E-01, 8.44E-02, &
    6.27E-02,  4.56E-02,  3.42E-02,  3.23E-02,  2.35E-02, &
    2.46E-02 /
  data ymeyer / &
    3.478E+04, 2.861E+04, 2.365E+04, 1.963E+04, 1.637E+04, &
    1.372E+04, 1.154E+04, 9.744E+03, 8.261E+03, 7.030E+03, &
    6.005E+03, 5.147E+03, 4.427E+03, 3.820E+03, 3.307E+03, &
    2.872E+03 /
  data yosb1 / &
    8.44E-01, 9.08E-01, 9.32E-01, 9.36E-01, 9.25E-01, &
    9.08E-01, 8.81E-01, 8.50E-01, 8.18E-01, 7.84E-01, &
    7.51E-01, 7.18E-01, 6.85E-01, 6.58E-01, 6.28E-01, &
    6.03E-01, 5.80E-01, 5.58E-01, 5.38E-01, 5.22E-01, &
    5.06E-01, 4.90E-01, 4.78E-01, 4.67E-01, 4.57E-01, &
    4.48E-01, 4.38E-01, 4.31E-01, 4.24E-01, 4.20E-01, &
    4.14E-01, 4.11E-01, 4.06E-01 /
  data yosb2 / &
    1.366E+00, 1.191E+00, 1.112E+00, 1.013E+00, 9.91E-01, &
    8.85E-01,  8.31E-01,  8.47E-01,  7.86E-01,  7.25E-01, &
    7.46E-01,  6.79E-01,  6.08E-01,  6.55E-01,  6.16E-01, &
    6.06E-01,  6.02E-01,  6.26E-01,  6.51E-01,  7.24E-01, &
    6.49E-01,  6.49E-01,  6.94E-01,  6.44E-01,  6.24E-01, &
    6.61E-01,  6.12E-01,  5.58E-01,  5.33E-01,  4.95E-01, &
    5.00E-01,  4.23E-01,  3.95E-01,  3.75E-01,  3.72E-01, &
    3.91E-01,  3.96E-01,  4.05E-01,  4.28E-01,  4.29E-01, &
    5.23E-01,  5.62E-01,  6.07E-01,  6.53E-01,  6.72E-01, &
    7.08E-01,  6.33E-01,  6.68E-01,  6.45E-01,  6.32E-01, &
    5.91E-01,  5.59E-01,  5.97E-01,  6.25E-01,  7.39E-01, &
    7.10E-01,  7.29E-01,  7.20E-01,  6.36E-01,  5.81E-01, &
    4.28E-01,  2.92E-01,  1.62E-01,  9.8E-02,   5.4E-02 /

  nex = uiparm(1)
!
!  Rosenbrock.
!
  if ( nex == 1 ) then

    r(1) = 10.0E+00 * ( x(2) - x(1)**2 )
    r(2) = 1.0E+00 - x(1)
!
!  Helix.
!
  else if ( nex == 2 ) then

    theta = atan2 ( x(2), x(1) ) / twopi

    if ( x(1) <= 0.0E+00 .and. x(2) <= 0.0E+00 ) then
      theta = theta + 1.0E+00
    end if

    r(1) = 10.0E+00 * ( x(3) - 10.0E+00 * theta )
    r(2) = 10.0E+00 * ( sqrt ( x(1)**2 + x(2)**2 ) - 1.0E+00 )
    r(3) = x(3)
!
!  Singular.
!
  else if ( nex == 3 ) then

    r(1) = x(1) + 10.0E+00 * x(2)
    r(2) = sqrt ( 5.0E+00 ) * ( x(3) - x(4) )
    r(3) = ( x(2) - 2.0E+00 * x(3) )**2
    r(4) = sqrt ( 10.0E+00 ) * ( x(1) - x(4) )**2
!
!  Woods.
!
  else if ( nex == 4 ) then

    r(1) = 10.0E+00 * ( x(2) - x(1)**2 )
    r(2) = 1.0E+00 - x(1)
    r(3) = sqrt ( 90.0E+00 ) * ( x(4) - x(3)**2 )
    r(4) = 1.0E+00 - x(3)
    r(5) = sqrt ( 9.9E+00 ) * ( x(2) + x(4) - 2.0E+00 )
    t = sqrt ( 0.2E+00 )
    r(6) = t * ( x(2) - 1.0E+00 )
    r(7) = t * ( x(4) - 1.0E+00 )
!
!  Zangwill.
!
  else if ( nex == 5 ) then

    r(1) =  x(1) - x(2) + x(3)
    r(2) = -x(1) + x(2) + x(3)
    r(3) =  x(1) + x(2) - x(3)
!
!  Engvall.
!
  else if ( nex == 6 ) then

    r(1) = x(1)**2 + x(2)**2 + x(3)**2 - 1.0E+00
    r(2) = x(1)**2 + x(2)**2 + ( x(3) - 2.0E+00 )**2 - 1.0E+00
    r(3) = x(1) + x(2) + x(3) - 1.0E+00
    r(4) = x(1) + x(2) - x(3) + 1.0E+00
    r(5) = x(1)**3 + 3.0E+00 * x(2)**2 &
      + ( 5.0E+00 * x(3) - x(1) + 1.0E+00 )**2 - 36.0E+00
!
!  Branin.
!
  else if ( nex == 7 ) then

    r(1) = 4.0E+00 * ( x(1) + x(2) )
    r(2) = r(1) + ( x(1) - x(2) ) * ( ( x(1) - 2.0E+00 )**2 + &
           x(2)**2 - 1.0E+00 )
!
!  Beale.
!
  else if ( nex == 8 ) then

    r(1) = 1.5E+00   - x(1) * ( 1.0E+00 - x(2)    )
    r(2) = 2.25E+00  - x(1) * ( 1.0E+00 - x(2)**2 )
    r(3) = 2.625E+00 - x(1) * ( 1.0E+00 - x(2)**3 )
!
!  Cragg and Levy.
!
  else if ( nex == 9 ) then

    r(1) = ( exp ( x(1) ) - x(2) )**2
    r(2) = 10.0E+00 * ( x(2) - x(3) )**3
    r(3) = ( sin ( x(3 ) - x(4) ) / cos ( x(3) - x(4) ) )**2
    r(4) = x(1)**4
    r(5) = x(4) - 1.0E+00
!
!  Box.
!
  else if ( nex == 10 ) then

    if ( expmax == 0.0E+00 ) then
      expmax = 1.999E+00 * log ( huge ( expmax ) )
    end if

    if ( expmin == 0.0E+00 ) then
      expmin = 1.999E+00 * log ( tiny ( expmin ) )
    end if

    if ( min ( x(1), x(2), x(3) ) <= -expmax ) then
      nfcall = -1
      return
    end if

    do i = 1, 10

      ti = -0.1E+00 * real ( i, kind = 4 )

      t1 = ti * x(1)

      if ( t1 <= expmin ) then
        e1 = 0.0E+00
      else
        e1 = exp ( t1 )
      end if

      t2 = ti * x(2)

      if ( t2 <= expmin ) then
        e2 = 0.0E+00
      else
        e2 = exp ( t2 )
      end if

      r(i) = ( e1 - e2 ) - x(3) * ( exp ( ti ) - exp ( 10.0E+00 * ti ) )

   end do
!
!  Davidon 1.
!
  else if ( nex == 11 ) then

    do i = 1, n-1
      r1 = 0.0E+00
      ti = real ( i, kind = 4 )
      t = 1.0E+00
      do j = 1, p
        r1 = r1 + t * x(j)
        t = t * ti
      end do
      r(i) = r1
    end do
    r(n) = x(1) - 1.0E+00
!
!  Freudenstein and Roth.
!
  else if ( nex == 12 ) then

    r(1) = -13.0E+00 + x(1) -  2.0E+00 * x(2) + 5.0E+00 * x(2)**2 - x(2)**3
    r(2) = -29.0E+00 + x(1) - 14.0E+00 * x(2) +           x(2)**2 + x(2)**3
!
!  Watson6.
!  Watson9.
!  Watson12.
!  Watson20.
!
  else if ( nex == 13 .or. &
            nex == 14 .or. &
            nex == 15 .or. &
            nex == 16 ) then

    do i = 1, 29
      ti = real ( i, kind = 4 ) / 29.0E+00
      r1 = 0.0E+00
      r2 = x(1)
      t = 1.0E+00
      do j = 2, p
        r1 = r1 + real ( j - 1, kind = 4 ) * t * x(j)
        t = t * ti
        r2 = r2 + t * x(j)
      end do
      r(i) = r1 - r2 * r2 - 1.0E+00
    end do
    r(30) = x(1)
    r(31) = x(2) - x(1)**2 - 1.0E+00
!
!  Chebyquad.
!
  else if ( nex == 17 ) then

    r(1:n) = 0.0E+00

    do j = 1, n
      tim1 = 1.0E+00
      ti = 2.0E+00 * x(j) - 1.0E+00
      z = ti + ti
      do i = 1, n
        r(i) = r(i) + ti
        tip1 = z * ti - tim1
        tim1 = ti
        ti = tip1
      end do
    end do

    do i = 1, n
      ti = 0.0E+00
      if ( mod ( i, 2 ) == 0 ) then
        ti = -1.0E+00 / real ( i * i - 1, kind = 4 )
      end if
      r(i) = ti - r(i) / real ( n, kind = 4 )
    end do
!
!  Brown and Dennis.
!
  else if ( nex == 18 ) then

    do i = 1, n
      ti = 0.2E+00 * real ( i, kind = 4 )
      r(i) = ( x(1) + x(2) * ti - exp ( ti ) )**2 + &
             ( x(3) + x(4) * sin ( ti ) - cos ( ti ) )**2
    end do
!
!  Bard.
!
  else if ( nex == 19 ) then

    do i = 1, 15
      u = real ( i, kind = 4 )
      v = 16.0E+00 - u
      w = min ( u, v )
      r(i) = ybard(i) - ( x(1) + u / ( x(2) * v + x(3) * w ) )
    end do
!
!  Jennrich and Sampson.
!
  else if ( nex == 20 ) then

    do i = 1, 10
      ti = real ( i, kind = 4 )
      r(i) = 2.0E+00 + 2.0E+00 * ti - ( exp ( ti * x(1) ) + exp ( ti * x(2) ) )
    end do
!
!  Kowalik and Osborne.
!
  else if ( nex == 21 ) then

     do i = 1, 11
       r(i) = ykow(i) - x(1) * ( ukow(i)**2 + x(2) * ukow(i) ) &
         / ( ukow(i)**2 + x(3) * ukow(i) + x(4) )
     end do
!
!  Osborne 1.
!
  else if ( nex == 22 ) then

    do i = 1, 33
      ti = 10.0E+00 * real ( 1 - i, kind = 4 )
      r(i) = yosb1(i) - ( x(1) + x(2) * exp ( x(4) * ti ) + &
        x(3) * exp ( x(5) * ti ) )
    end do
!
!  Osborne 2.
!
!  UFTOLG is a machine-dependent constant.  It is just slightly
!  larger than the log of the smallest positive machine number.
!
  else if ( nex == 23 ) then

    if ( uftolg == 0.0E+00 ) then
      uftolg = 1.999E+00 * log ( tiny ( uftolg ) )
    end if

    do i = 1, 65
      ti = 0.1E+00 * real ( 1 - i, kind = 4 )
      ri = x(1) * exp ( x(5) * ti )
      do j = 2, 4
        theta = -x(j+4) * ( ti + x(j+7) )**2
        if ( theta <= uftolg ) then
          t = 0.0E+00
        else
          t = exp ( theta )
        end if
        ri = ri + x(j) * t
      end do
      r(i) = yosb2(i) - ri
    end do
!
!  Madsen.
!
  else if ( nex == 24 ) then

    r(1) = x(1)**2 + x(2)**2 + x(1) * x(2)
    r(2) = sin ( x(1) )
    r(3) = cos ( x(2) )
!
!  Meyer.
!
  else if ( nex == 25 ) then

    do i = 1, 16
      ti = real ( 5 * i + 45, kind = 4 )
      r(i) = x(1) * exp ( x(2) / ( ti + x(3) ) ) - ymeyer(i)
    end do
!
!  Brown5.
!  Brown10.
!  Brown30.
!  Brown40.
!
  else if ( nex == 26 .or. &
            nex == 27 .or. &
            nex == 28 .or. &
            nex == 29 ) then

    r(1:n-1) = x(1:n-1) + sum ( x(1:n) ) - real ( n + 1, kind = 4 )
    r(n) = product ( x(1:n) ) - 1.0E+00
!
!  Bard + 10.
!
  else if ( nex == 30 ) then

    do i = 1, 15
      u = real ( i, kind = 4 )
      v = 16.0E+00 - u
      w = min ( u, v )
      r(i) = ybard(i) - ( x(1) + u / ( x(2) * v + x(3) * w ) ) + 10.0E+00
   end do
!
!  Kowalik and Osborne + 10.
!
  else if ( nex == 31 ) then

    do i = 1, 11
      r(i) = ykow(i) - x(1) * ( ukow(i)**2 + x(2) * ukow(i) ) &
        / ( ukow(i)**2 + x(3) * ukow(i) + x(4) ) + 10.0E+00
    end do
!
!  Meyer + 10.
!
  else if ( nex == 32 ) then

    do i = 1, 16
      ti = real ( 5 * i + 45, kind = 4 )
      r(i) = x(1) * exp ( x(2) / ( ti + x(3) ) ) - ymeyer(i) + 10.0E+00
    end do
!
!  Watson6 + 10.
!  Watson9 + 10.
!  Watson12 + 10.
!  Watson20 + 10.
!
  else if ( nex == 33 .or. &
            nex == 34 .or. &
            nex == 35 .or. &
            nex == 36 ) then

    do i = 1, 29
      ti = real ( i, kind = 4 ) / 29.0E+00
      r1 = 0.0E+00
      r2 = x(1)
      t = 1.0E+00
      do j = 2, p
        r1 = r1 + real ( j - 1, kind = 4 ) * t * x(j)
        t = t * ti
        r2 = r2 + t*x(j)
      end do
      r(i) = r1 - r2 * r2 - 1.0E+00 + 10.0E+00
    end do
    r(30) = x(1) + 10.0E+00
    r(31) = x(2) - x(1)**2 - 1.0E+00 + 10.0E+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TESTR - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem index NEX = ', nex
    stop

  end if

  return
end
subroutine xinit ( p, x, nex )

!*****************************************************************************80
!
!! XINIT initializes the solution vector X.
!
!  Discussion:
!
!    This routine initializes the solution vector X.
!
!  Modified:
!
!    27 March 2006
!
!  Parameters:
!
!    Input, integer P, the  number of parameters.
!
!    Output, real X(P), the initial value for the solution vector.
!
!    Input, integer NEX, the test problem number, between 1 and 36.
!
  implicit none

  integer p

  integer i
  integer nex
  real pp1inv
  real x(p)
!
!  Rosenbrock.
!
  if ( nex == 1 ) then

    x(1) = -1.2E+00
    x(2) = 1.0E+00
!
!  Helix.
!
  else if ( nex == 2 ) then

    x(1) = -1.0E+00
    x(2) = 0.0E+00
    x(3) = 0.0E+00
!
!  Singular.
!
  else if ( nex == 3 ) then

    x(1) = 3.0E+00
    x(2) = -1.0E+00
    x(3) = 0.0E+00
    x(4) = 1.0E+00
!
!  Woods.
!
  else if ( nex == 4 ) then

    x(1) = -3.0E+00
    x(2) = -1.0E+00
    x(3) = -3.0E+00
    x(4) = -1.0E+00
!
!  Zangwill.
!
  else if ( nex == 5 ) then

    x(1) = 100.0E+00
    x(2) = -1.0E+00
    x(3) = 2.5E+00
!
!  Engvall.
!
  else if ( nex == 6 ) then

    x(1) = 1.0E+00
    x(2) = 2.0E+00
    x(3) = 0.0E+00
!
!  Branin.
!
  else if ( nex == 7 ) then

    x(1) = 2.0E+00
    x(2) = 0.0E+00
!
!  Beale.
!
  else if ( nex == 8 ) then

    x(1) = 0.1E+00
    x(2) = 0.1E+00
!
!  Cragg and Levy.
!
  else if ( nex == 9 ) then

    x(1) = 1.0E+00
    x(2) = 2.0E+00
    x(3) = 2.0E+00
    x(4) = 2.0E+00
!
!  Box.
!
  else if ( nex == 10 ) then

    x(1) = 0.0E+00
    x(2) = 10.0E+00
    x(3) = 20.0E+00
!
!  Davidon 1.
!
  else if ( nex == 11 ) then

    do i = 1, p
      x(i) = 0.0E+00
    end do
!
!  Freudenstein and Roth.
!
  else if ( nex == 12 ) then

    x(1) = 15.0E+00
    x(2) = -2.0E+00
!
!  Watson6.
!  Watson9.
!  Watson12.
!  Watson20.
!
  else if ( nex == 13 .or. &
            nex == 14 .or. &
            nex == 15 .or. &
            nex == 16 ) then

    do i = 1, p
      x(i) = 0.0E+00
    end do
!
!  Chebyquad.
!
  else if ( nex == 17 ) then

    do i = 1, p
      x(i) = real ( i, kind = 4 ) / real ( p + 1, kind = 4 )
    end do
!
!  Brown and Dennis.
!
  else if ( nex == 18 ) then

    x(1) = 25.0E+00
    x(2) = 5.0E+00
    x(3) = -5.0E+00
    x(4) = -1.0E+00
!
!  Bard.
!
  else if ( nex == 19 ) then

    x(1) = 1.0E+00
    x(2) = 1.0E+00
    x(3) = 1.0E+00
!
!  Jennrich and Sampson.
!
  else if ( nex == 20 ) then

    x(1) = 0.3E+00
    x(2) = 0.4E+00
!
!  Kowalik and Osborne.
!
  else if ( nex == 21 ) then

    x(1) = 0.25E+00
    x(2) = 0.39E+00
    x(3) = 0.415E+00
    x(4) = 0.39E+00
!
!  Osborne 1.
!
  else if ( nex == 22 ) then

    x(1) = 0.5E+00
    x(2) = 1.5E+00
    x(3) = -1.0E+00
    x(4) = 0.01E+00
    x(5) = 0.02E+00
!
!  Osborne 2.
!
  else if ( nex == 23 ) then

    x(1) = 1.3E+00
    x(2) = 0.65E+00
    x(3) = 0.65E+00
    x(4) = 0.7E+00
    x(5) = 0.60E+00
    x(6) = 3.0E+00
    x(7) = 5.0E+00
    x(8) = 7.0E+00
    x(9) = 2.0E+00
    x(10) = 4.5E+00
    x(11) = 5.5E+00
!
!  Madsen.
!
  else if ( nex == 24 ) then

    x(1) = 3.0E+00
    x(2) = 1.0E+00
!
!  Meyer.
!
  else if ( nex == 25 ) then

    x(1) = 0.02E+00
    x(2) = 4000.0E+00
    x(3) = 250.0E+00
!
!  Brown5.
!  Brown10.
!  Brown30.
!  Brown40.
!
  else if ( nex == 26 .or. &
            nex == 27 .or. &
            nex == 28 .or. &
            nex == 29 ) then

    x(1:p) = 0.5E+00
!
!  Bard + 10.
!
  else if ( nex == 30 ) then

    x(1) = 1.0E+00
    x(2) = 1.0E+00
    x(3) = 1.0E+00
!
!  Kowalik and Osborne + 10.
!
  else if ( nex == 31 ) then

    x(1) = 0.25E+00
    x(2) = 0.39E+00
    x(3) = 0.415E+00
    x(4) = 0.39E+00
!
!  Meyer + 10.
!
  else if ( nex == 32 ) then

    x(1) = 0.02E+00
    x(2) = 4000.0E+00
    x(3) = 250.0E+00
!
!  Watson6 + 10.
!  Watson9 + 10.
!  Watson12 + 10.
!  Watson20 + 10.
!
  else if ( nex == 33 .or. &
            nex == 34 .or. &
            nex == 35 .or. &
            nex == 36 ) then

    do i = 1, p
      x(i) = 0.0E+00
    end do
!
!  Illegal value.
!
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XINIT - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem index NEX = ', nex
    stop

  end if

  return
end
