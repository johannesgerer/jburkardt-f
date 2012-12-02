function approxres ( i, h, s, givens, ldg )

!*****************************************************************************80
!
!! APPROXRES approximates the residual using a Givens updating scheme.
!
!  Discussion:
!
!    The rotation matrix is formed using
!
!      [H(I),H(I+1)]'
!
!    with the intent of zeroing H(I+1), but here is applied to the 2x1
!    vector
!
!      [S(I), S(I+1)]'.
!
!  Reference:
!
!    Barrett, Berry, Chan, Demmel, Donato, Dongarra, Eijkhout, Pozo,
!      Romine, van der Vorst,
!    Templates for the Solution of Linear Systems: Building Blocks
!      for Iterative Methods,
!    SIAM, 1994.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer I, ...
!
!    Input, real H(*), ...
!
!    put, real S(*), ...
!
!    put, GIVENS(LDG,2), ...
!
!    Input, integer LDG, the leading dimension of GIVENS.
!
!    Output, real APPROXRES, the approximate residual.
!
  implicit none

  integer ldg

  real approxres
  real givens(ldg,2)
  real h(*)
  integer i
  real s(*)

  call givens_set ( h(i), h(i+1), givens(i,1), givens(i,2) )

  call rotvec ( s(i), s(i+1), givens(i,1), givens(i,2) )

  approxres = s(i+1)

  return
end
subroutine bicg ( n, b, x, work, iter, resid, matvec, matvect, &
  psolve, psolve_t, info, curpform )

!*****************************************************************************80
!
!! BICG implements the BiConjugate Gradient method.
!
!  Discussion:
!
!    Preconditioning is used.  The convergence test is:
!
!      norm(b-A*x) / norm(b) < TOL.
!
!  Reference:
!
!    Barrett, Berry, Chan, Demmel, Donato, Dongarra, Eijkhout, Pozo,
!      Romine, van der Vorst,
!    Templates for the Solution of Linear Systems: Building Blocks
!      for Iterative Methods,
!    SIAM, 1994.
!
!  Modified:
!
!    27 March 2001
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real B(N), the right hand side.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Workspace, real WORK(N,6).
!    Workspace for residual, direction vector, etc.
!    Note that Z and Q, and ZTLD and QTLD share workspace.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, actual number of iterations performed.
!
!    Input/output, real RESID.
!    On input, the allowable convergence measure for norm(b-A*x) / norm(b).
!    On output, the final value of this measure.
!
!  MATVEC  (external subroutine)
!          The user must provide a subroutine to perform the
!          matrix-vector product
!
!               y : = alpha*A*x+beta*y,
!
!          where alpha and beta are scalars, x and y are vectors,
!          and A is a matrix. Vector x must remain unchanged.
!          The solution is over-written on vector y.
!
!          The call is:
!
!             CALL MATVEC(ALPHA, X, BETA, Y)
!
!          The matrix is passed into the routine in a common block.
!
!  MATVECT  (external subroutine)
!          The user must provide a subroutine to perform the
!          matrix-vector product
!
!               y : = alpha*A'*x+beta*y,
!
!          where alpha and beta are scalars, x and y are vectors,
!          and A' is the tranpose of a matrix A. Vector x must remain
!          unchanged.
!          The solution is over-written on vector y.
!
!          The call is:
!
!             CALL MATVECT(ALPHA, X, BETA, Y)
!
!          The matrix is passed into the routine in a common block.
!
!  PSOLVE  (external subroutine)
!          The user must provide a subroutine to perform the
!          preconditioner solve routine for the linear system
!
!               M*x = b,
!
!          where x and b are vectors, and M a matrix. Vector b must
!          remain unchanged.
!          The solution is over-written on vector x.
!
!          The call is:
!
!             CALL PSOLVE(X, B)
!
!          The preconditioner is passed into the routine in a common block.
!
!  PSOLVET  (external subroutine)
!          The user must provide a subroutine to perform the
!          preconditioner solve routine for the linear system
!
!               M'*x = b,
!
!          where x and y are vectors, and M' is the tranpose of a
!          matrix M. Vector b must remain unchanged.
!          The solution is over-written on vector x.
!
!          The call is:
!
!             CALL PSOLVET(X, B)
!
!          The preconditioner is passed into the routine in a common block.
!
!  INFO    (output) integer
!
!        = 0: Successful exit. Iterated approximate solution returned.
!
!          >  0: Convergence to tolerance not achieved. This will be
!                set to the number of iterations performed.
!
!          <  0: Illegal input parameter, or breakdown occurred
!                during iteration.
!
!                Illegal parameter:
!
!                   -1: matrix dimension N < 0
!                   -3: Maximum number of iterations ITER < = 0.
!
!                BREAKDOWN: If parameters RHO or OMEGA become smaller
!                   than some tolerance, the program will terminate.
!                   Here we check against tolerance BREAKTOL.
!
!                  -10: RHO < BREAKTOL: RHO and RTLD have become
!                                       orthogonal.
!
!                  BREAKTOL is set in function GETBREAK.
!
  implicit none

  integer n

  real b(n)
  real bnrm2
  character ( len = 8 ) curpform
  integer job
  integer info
  integer iter
  integer ndx1
  integer ndx2
  real resid
  real sclr1
  real sclr2
  real snrm2
  real tol
  real work(n*6)
  real x(n)

  external matvec
  external matvect
  external psolve
  external psolve_t

  info = 0
!
!  Test the input parameters.
!
  if ( n < 1 ) then
    info = -1
    write ( *, * ) ' '
    write ( *, * ) 'BICG - Fatal error!'
    write ( *, * ) '  N is less than 1.'
    write ( *, * ) '  N = ',n
    stop
  end if

  if ( iter <= 0 ) then
    info = -3
    return
  end if
!
!  Stop test may need some indexing info from REVCOM
!  use the init call to send the request across. REVCOM
!  will note these requests, and everytime it asks for
!  stop test to be done, it will provide the indexing info.
!
  ndx1 = 1
  ndx2 = -1
  tol = resid
  bnrm2 = snrm2 ( n, b, 1 )
!
!  First time call always init.
!
  job = 1

  do

    call bicg_revcom ( n, b, x, work, iter, resid, info, ndx1, ndx2, &
      sclr1, sclr2, job )
!
!  -1: Termination.
!
    if ( job == -1 ) then

      exit
!
!  1: Compute WORK(NDX2) = SCLR1 * A * WORK(NDX1) + SCLR2 * WORK(NDX2).
!
    else if ( job == 1 ) then

      call matvec ( sclr1, work(ndx1), sclr2, work(ndx2) )
!
!  2: Compute WORK(NDX2) = SCLR1 * A' * WORK(NDX1) + SCLR2 * WORK(NDX2).
!
    else if ( job == 2 ) then

      call matvect ( sclr1, work(ndx1), sclr2, work(ndx2) )
!
!  3: Solve M * WORK(NDX1) = WORK(NDX2).
!
    else if ( job == 3 ) then

      call psolve ( n, work(ndx1), work(ndx2), curpform )
!
!  4: Solve M' * WORK(NDX1) = WORK(NDX2).
!
    else if ( job == 4 ) then

      call psolve_t ( n, work(ndx1), work(ndx2), curpform )
!
!  5: Compute WORK(NDX2) = SCLR1 * A * X + SCLR2 * WORK(NDX2).
!
    else if ( job == 5 ) then

      call matvec ( sclr1, x, sclr2, work(ndx2) )
!
!  6: Do a stopping test on the relative residual reduction.
!
    else if ( job == 6 ) then

      call stopb ( n, work(ndx1), bnrm2, resid, tol, info )

    end if

    job = 2

  end do

  return
end
subroutine bicg_revcom ( n, b, x, work, iter, resid, info, ndx1, ndx2, &
  sclr1, sclr2, job )

!*****************************************************************************80
!
!! BICG_REVCOM is controlled by BICG using reverse communication.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real B(N), the right hand side.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Workspace, real WORK(N,6).
!    Workspace for residual, direction vector, etc.
!    Note that Z and Q, and ZTLD and QTLD share workspace.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, actual number of iterations performed.
!
!    Input/output, real RESID.
!    On input, the allowable convergence measure for norm(b-A*x) / norm(b).
!    On output, the final value of this measure.
!
!  INFO    (output) integer
!
!        = 0: Successful exit. Iterated approximate solution returned.
!
!          >  0: Convergence to tolerance not achieved. This will be
!                set to the number of iterations performed.
!
!          <  0: Illegal input parameter, or breakdown occurred
!                during iteration.
!
!                Illegal parameter:
!
!                   -1: matrix dimension N < 0
!                   -3: Maximum number of iterations ITER < = 0.
!                   -5: Erroneous NDX1/NDX2 in INIT call.
!                   -6: Erroneous RLBL.
!
!                BREAKDOWN: If parameters RHO or OMEGA become smaller
!                   than some tolerance, the program will terminate.
!                   Here we check against tolerance BREAKTOL.
!
!                  -10: RHO < BREAKTOL: RHO and RTLD have become
!                                       orthogonal.
!
!                  BREAKTOL is set in function GETBREAK.
!
!  NDX1    (input/output) integer.
!  NDX2    On entry in INIT call contain indices required by interface
!          level for stopping test.
!          All other times, used as output, to indicate indices into
!          WORK[] for the MATVEC, PSOLVE done by the interface level.
!
!  SCLR1   (output) real.
!  SCLR2   Used to pass the scalars used in MATVEC.
!
!    Input/output, integer JOB.
!    Used to communicate job code between the two levels.
!
  implicit none

  integer n

  real alpha
  real b(n)
  real beta
  real bnrm2
  real sdot
  real snrm2
  real getbreak
  integer i
  integer job
  integer info
  integer iter
  integer maxit
  integer ndx1
  integer ndx2
  integer need1
  integer need2
  integer p
  integer ptld
  integer q
  integer qtld
  integer r
  real resid
  real rho
  real :: rho1 = 0.0E+00
  real rhotol
  integer :: rlbl = 0
  integer rtld
  real sclr1
  real sclr2
  real tol
  real work(n,6)
  real x(n)
  integer z
  integer ztld
!
!  G95 won't let us SAVE everything and SAVE specific things as well!
!
  save

  if ( job == 1 ) then
    go to 1
  else if ( job == 2 ) then
    if ( rlbl == 2) go to 2
    if ( rlbl == 3) go to 3
    if ( rlbl == 4) go to 4
    if ( rlbl == 5) go to 5
    if ( rlbl == 6) go to 6
    if ( rlbl == 7) go to 7
    info = -6
    go to 20
  end if

 1    continue

  info = 0
  maxit = iter
  tol = resid
!
!  Alias workspace columns.
!
  r = 1
  rtld = 2
  z = 3
  ztld = 4
  p = 5
  ptld = 6
  q = 3
  qtld = 4
!
!  Check if caller will need indexing info.
!
  if ( ndx1 == -1 ) then
    need1 = ndx1
  else if ( ndx1 == 1 ) then
    need1 = ((r-1)*n)+1
  else if ( ndx1 == 2 ) then
    need1 = ((rtld-1)*n)+1
  else if ( ndx1 == 3 ) then
    need1 = ((z-1)*n)+1
  else if ( ndx1 == 4 ) then
    need1 = ((ztld-1)*n)+1
  else if ( ndx1 == 5 ) then
    need1 = ((p-1)*n)+1
  else if ( ndx1 == 6 ) then
    need1 = ((ptld-1)*n)+1
  else if ( ndx1 == 7 ) then
    need1 = ((q-1)*n)+1
  else if ( ndx1 == 8 ) then
    need1 = ((qtld-1)*n)+1
  else
    info = -5
    go to 20
  end if

  if ( ndx2 == -1 ) then
    need2 = ndx2
  else if ( ndx2 == 1 ) then
    need2 = ((r-1)*n)+1
  else if ( ndx2 == 2 ) then
    need2 = ((rtld-1)*n)+1
  else if ( ndx2 == 3 ) then
    need2 = ((z-1)*n)+1
  else if ( ndx2 == 4 ) then
    need2 = ((ztld-1)*n)+1
  else if ( ndx2 == 5 ) then
    need2 = ((p-1)*n)+1
  else if ( ndx2 == 6 ) then
    need2 = ((ptld-1)*n)+1
  else if ( ndx2 == 7 ) then
    need2 = ((q-1)*n)+1
  else if ( ndx2 == 8 ) then
    need2 = ((qtld-1)*n)+1
  else
    info = -5
    go to 20
  end if
!
!  Set breakdown parameters.
!
  rhotol = getbreak()
!
!  Set the initial residual.
!
  work(1:n,r) = b(1:n)

  if ( snrm2(n,x,1) /= 0.0E+00 ) then
    sclr1 = - 1.0E+00
    sclr2 = 0.0E+00
    ndx1 = ((rtld-1)*n)+1
    ndx2 = ((r   -1)*n)+1
    rlbl = 2
    job = 5
    return
  end if

 2    continue

  if ( snrm2(n,work(1,r),1) <= tol ) then
    go to 30
  end if

  work(1:n,rtld) = work(1:n,r)

  bnrm2 = snrm2(n,b,1)

  if ( bnrm2 == 0.0E+00 ) then
    bnrm2 = 1.0E+00
  end if

  iter = 0

   10 continue
!
!  Perform BiConjugate Gradient iteration.
!
  iter = iter + 1
!
!  Compute direction vectors PK and PTLD.
!
  ndx1 = ((z-1)*n)+1
  ndx2 = ((r-1)*n)+1
  rlbl = 3
  job = 3
  return

 3    continue

  ndx1 = ((ztld-1)*n)+1
  ndx2 = ((rtld-1)*n)+1
  rlbl = 4
  job = 4
  return

 4    continue

  rho = sdot(n,work(1,z),1,work(1,rtld),1)

  if ( abs(rho) < rhotol) then
    go to 25
  end if

  if ( 1 < iter ) then

    beta = rho / rho1
    call saxpy ( n, beta, work(1,p), 1, work(1,z), 1 )
    call saxpy ( n, beta, work(1,ptld), 1, work(1,ztld), 1 )

    work(1:n,p) = work(1:n,z)
    work(1:n,ptld) = work(1:n,ztld)

  else

    do i = 1, n
      work(i,p) = work(i,z)
    end do

    do i = 1, n
      work(i,ptld) = work(i,ztld)
    end do

  end if

  sclr1 = 1.0E+00
  sclr2 = 0.0E+00
  ndx1 = ((p-1)*n)+1
  ndx2 = ((q-1)*n)+1
  rlbl = 5
  job = 1
  return

 5    continue

  sclr1 = 1.0E+00
  sclr2 = 0.0E+00
  ndx1 = ((ptld-1)*n)+1
  ndx2 = ((qtld-1)*n)+1
  rlbl = 6
  job = 2
  return

 6    continue

  alpha = rho / sdot ( n, work(1,ptld), 1, work(1,q), 1 )
!
!  Compute current solution vector x.
!
  call saxpy ( n, alpha, work(1,p), 1, x, 1 )
!
!  Compute residual vector RK, find norm,
!  then check for tolerance.
!
  call saxpy ( n, - alpha, work(1,q), 1, work(1,r), 1 )

  ndx1 = need1
  ndx2 = need2
  rlbl = 7
  job = 6
  return

 7    continue

  if ( info == 1 ) then
    go to 30
  end if

  if ( iter == maxit ) then
    info = 1
    go to 20
  end if

  call saxpy ( n, -alpha, work(1,qtld), 1, work(1,rtld), 1 )
  rho1 = rho

  go to 10

   20 continue
!
!  Iteration fails.
!
  rlbl = -1
  job = -1
  return

   25 continue
!
!  Set breakdown flag.
!
  info = -10
  rlbl = -1
  job = -1
  return

   30 continue
!
!  Iteration successful; return.
!
  info = 0
  rlbl = -1
  job = -1

  return
end
subroutine bicgstab ( n, b, x, work, iter, resid, matvec, psolve, &
  info, curpform )

!*****************************************************************************80
!
!! BICGSTAB implements the BiConjugate Gradient Stabilized method.
!
!  Discussion:
!
!    Preconditioning is used.  The convergence test is:
!
!      norm(b-A*x) / norm(b) < TOL.
!
!  Modified:
!
!    02 May 2000
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real B(N), the right hand side vector.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Workspace, real WORK(N,7)
!    Workspace for residual, direction vector, etc.
!    Note that vectors R and S shared the same workspace.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, actual number of iterations performed.
!
!    Input/output, real RESID.
!    On input, the allowable convergence measure for norm(b-A*x) / norm(b).
!    On output, the final value of this measure.
!
!  MATVEC  (external subroutine)
!          The user must provide a subroutine to perform the
!          matrix-vector product
!
!               y : = alpha*A*x+beta*y,
!
!          where alpha and beta are scalars, x and y are vectors,
!          and A is a matrix. Vector x must remain unchanged.
!          The solution is over-written on vector y.
!
!          The call is:
!
!             CALL MATVEC(ALPHA, X, BETA, Y)
!
!          The matrix is passed into the routine in a common block.
!
!  PSOLVE  (external subroutine)
!          The user must provide a subroutine to perform the
!          preconditioner solve routine for the linear system
!
!               M*x = b,
!
!          where x and b are vectors, and M a matrix. Vector b must
!          remain unchanged.
!          The solution is over-written on vector b.
!
!          The call is:
!
!             CALL PSOLVE(X, B)
!
!          The preconditioner is passed into the routine in a common block.
!
!  INFO    (output) integer
!
!        = 0: Successful exit. Iterated approximate solution returned.
!
!          >  0: Convergence to tolerance not achieved. This will be
!                set to the number of iterations performed.
!
!          <  0: Illegal input parameter, or breakdown occurred
!                during iteration.
!
!                Illegal parameter:
!
!                   -1: matrix dimension N < 0
!                   -3: Maximum number of iterations ITER < = 0.
!
!                BREAKDOWN: If parameters RHO or OMEGA become smaller
!                   than some tolerance, the program will terminate.
!                   Here we check against tolerance BREAKTOL.
!
!                  -10: RHO < BREAKTOL: RHO and RTLD have become
!                                       orthogonal.
!                  -11: OMEGA < BREAKTOL: S and T have become
!                                         orthogonal relative to T'*T.
!
!                  BREAKTOL is set in function GETBREAK.
!
  implicit none

  integer n

  real b(n)
  real bnrm2
  character ( len = 8 ) curpform
  integer job
  integer info
  integer iter
  integer ndx1
  integer ndx2
  real resid
  real sclr1
  real sclr2
  real snrm2
  real tol
  real work(n*7)
  real x(n)

  external matvec
  external psolve

  info = 0
!
!  Test the input parameters.
!
  if ( n < 1 ) then
    info = -1
    write ( *, * ) ' '
    write ( *, * ) 'BICGSTAB - Fatal error!'
    write ( *, * ) '  N is less than 1.'
    write ( *, * ) '  N = ', n
    stop
  end if

  if ( iter <= 0 ) then
    info = -3
    return
  end if
!
!  Stop test may need some indexing info from REVCOM
!  use the init call to send the request across. REVCOM
!  will note these requests, and everytime it asks for
!  stop test to be done, it will provide the indexing info.
!
!  1 == R;
!  2 == RTLD;
!  3 == P;
!  4 == V;
!  5 == T;
!  6 == PHAT;
!  7 == SHAT;
!  8 == S;
! -1 == ignore;
!  any other == error
!
  ndx1 = 1
  ndx2 = -1
  tol = resid
  bnrm2 = snrm2 ( n, b, 1 )
!
!  First time call always init.
!
  job = 1

  do

    call bicgstab_revcom ( n, b, x, work, iter, resid, info, ndx1, ndx2, &
      sclr1, sclr2, job )
!
!  -1: Terminate.
!
    if ( job == -1 ) then

      exit
!
!  1: Compute WORK(NDX2) = SCLR1 * A * WORK(NDX1) + SCLR2 * WORK(NDX2).
!
    else if ( job == 1 ) then

      call matvec ( sclr1, work(ndx1), sclr2, work(ndx2) )
!
!  2: Solve M * WORK(NDX1) = WORK(NDX2).
!
    else if ( job == 2 ) then

      call psolve ( n, work(ndx1), work(ndx2), curpform )
!
!  3: Compute WORK(NDX2) = SCLR1 * A * X + SCLR2 * WORK(NDX2).
!
    else if ( job == 3 ) then

      call matvec ( sclr1, x, sclr2, work(ndx2) )
!
!  4: Do a stopping test on the relative residual reduction.
!
    else if ( job == 4 ) then

      call stopb ( n, work(ndx1), bnrm2, resid, tol, info )

    end if

    job = 2

  end do

  return
end
subroutine bicgstab_revcom ( n, b, x, work, iter, resid, info, ndx1, &
  ndx2, sclr1, sclr2, job )

!*****************************************************************************80
!
!! BICGSTAB_REVCOM is controlled by BICGSTAB using reverse communication.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real B(N), the right hand side vector.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Workspace, real WORK(N,7).
!    Workspace for residual, direction vector, etc.
!    Note that vectors R and S shared the same workspace.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, actual number of iterations performed.
!
!    Input/output, real RESID.
!    On input, the allowable convergence measure for norm(b-A*x) / norm(b).
!    On output, the final value of this measure.
!
!  INFO    (output) integer
!
!        = 0: Successful exit. Iterated approximate solution returned.
!
!          >  0: Convergence to tolerance not achieved. This will be
!                set to the number of iterations performed.
!
!          <  0: Illegal input parameter, or breakdown occurred
!                during iteration.
!
!                Illegal parameter:
!
!                   -1: matrix dimension N < 0
!                   -3: Maximum number of iterations ITER < = 0.
!                   -5: Erroneous NDX1/NDX2 in INIT call.
!                   -6: Erroneous RLBL.
!
!                BREAKDOWN: If parameters RHO or OMEGA become smaller
!                   than some tolerance, the program will terminate.
!                   Here we check against tolerance BREAKTOL.
!
!                  -10: RHO < BREAKTOL: RHO and RTLD have become
!                                       orthogonal.
!                  -11: OMEGA < BREAKTOL: S and T have become
!                                         orthogonal relative to T'*T.
!
!                  BREAKTOL is set in function GETBREAK.
!
!  NDX1    (input/output) integer.
!  NDX2    On entry in INIT call contain indices required by interface
!          level for stopping test.
!          All other times, used as output, to indicate indices into
!          WORK[] for the MATVEC, PSOLVE done by the interface level.
!
!  SCLR1   (output) real.
!  SCLR2   Used to pass the scalars used in MATVEC.
!
!  JOB    (input/output) integer.
!          Used to communicate job code between the two levels.
!
  implicit none

  real, parameter :: one = 1.0E+00

  integer n

  real :: alpha = 0.0E+00
  real b(n)
  real beta
  real bnrm2
  real sdot
  real snrm2
  real getbreak
  integer i
  integer job
  integer info
  integer iter
  integer maxit
  integer ndx1
  integer ndx2
  integer need1
  integer need2
  real :: omega = 0.0E+00
  real omegatol
  integer p
  integer phat
  integer r
  real resid
  real rho
  real :: rho1 = 0.0E+00
  real rhotol
  integer :: rlbl = 0
  integer rtld
  integer s
  real sclr1
  real sclr2
  integer shat
  integer t
  real tol
  integer v
  real work(n,7)
  real x(n)

  save

  if ( job == 1 ) then

     go to 1
!
!  Handle a resumption.
!
  else if ( job == 2 ) then
     if ( rlbl == 2) go to 2
     if ( rlbl == 3) go to 3
     if ( rlbl == 4) go to 4
     if ( rlbl == 5) go to 5
     if ( rlbl == 6) go to 6
     if ( rlbl == 7) go to 7
     info = -6
     go to 20
  end if

 1    continue

  info = 0
  maxit = iter
  tol = resid
!
!  Alias workspace columns.
!
  r  = 1
  rtld = 2
  p  = 3
  v  = 4
  t  = 5
  phat = 6
  shat = 7
  s  = 1
!
!  Check if caller will need indexing info.
!
  if ( ndx1 == -1 ) then
    need1 = ndx1
  else if ( ndx1 == 1 ) then
    need1 = ((r-1)*n)+1
  else if ( ndx1 == 2 ) then
    need1 = ((rtld-1)*n)+1
  else if ( ndx1 == 3 ) then
    need1 = ((p-1)*n)+1
  else if ( ndx1 == 4 ) then
    need1 = ((v-1)*n)+1
  else if ( ndx1 == 5 ) then
    need1 = ((t-1)*n)+1
  else if ( ndx1 == 6 ) then
    need1 = ((phat-1)*n)+1
  else if ( ndx1 == 7 ) then
    need1 = ((shat-1)*n)+1
  else if ( ndx1 == 8 ) then
    need1 = ((s-1)*n)+1
  else
    write ( *, * ) ' '
    write ( *, * ) 'BICGSTAB_REVCOM - Fatal error!'
    write ( *, * ) '  Illegal value of NDX1 = ',ndx1
    stop
  end if

  if ( ndx2==-1 ) then
    need2 = ndx2
  else if ( ndx2 == 1 ) then
    need2 = ((r-1)*n)+1
  else if ( ndx2 == 2 ) then
    need2 = ((rtld-1)*n)+1
  else if ( ndx2 == 3 ) then
    need2 = ((p-1)*n)+1
  else if ( ndx2 == 4 ) then
    need2 = ((v-1)*n)+1
  else if ( ndx2 == 5 ) then
    need2 = ((t-1)*n)+1
  else if ( ndx2 == 6 ) then
    need2 = ((phat-1)*n)+1
  else if ( ndx2 == 7 ) then
    need2 = ((shat-1)*n)+1
  else if ( ndx2 == 8 ) then
    need2 = ((s-1)*n)+1
  else
    write ( *, * ) ' '
    write ( *, * ) 'BICGSTAB_REVCOM - Fatal error!'
    write ( *, * ) '  Illegal value of NDX2 = ',ndx2
    stop
  end if
!
!  Set parameter tolerances.
!
  rhotol = getbreak()
  omegatol = getbreak()
!
!  Set the initial residual.
!
  work(1:n,r) = b(1:n)

  if ( snrm2(n,x,1) /= 0.0E+00 ) then
    sclr1 = -1.0E+00
    sclr2 = 1.0E+00
    ndx1 = -1
    ndx2 = ((r-1)*n)+1
    rlbl = 2
    job = 3
    return
  end if

 2    continue

  if ( snrm2(n,work(1,r),1) <= tol ) then
    go to 30
  end if

  do i = 1, n
    work(i,rtld) = work(i,r)
  end do

  bnrm2 = snrm2(n,b,1)
  if ( bnrm2 == 0.0E+00 ) then
    bnrm2 = 1.0E+00
  end if

  iter = 0

   10 continue
!
!  Perform BiConjugate Gradient Stabilized iteration.
!
  iter = iter + 1
  rho = sdot(n,work(1,rtld),1,work(1,r),1)

  if ( abs(rho) < rhotol ) then
    go to 25
  end if
!
!  Compute vector P.
!
  if ( iter > 1 ) then
    beta = ( rho / rho1 ) * ( alpha / omega )
    call saxpy ( n, -omega, work(1,v), 1, work(1,p), 1 )
    call sscal ( n, beta, work(1,p), 1 )
    call saxpy ( n, one, work(1,r), 1, work(1,p), 1 )
  else
    do i = 1, n
      work(i,p) = work(i,r)
    end do
  end if
!
!  Compute direction adjusting vector PHAT and scalar ALPHA.
!
  ndx1 = ((phat-1)*n)+1
  ndx2 = ((p   -1)*n)+1
  rlbl = 3
  job = 2
  return

 3    continue

  ndx1 = ((phat-1)*n)+1
  ndx2 = ((v   -1)*n)+1
  sclr1 = 1.0E+00
  sclr2 = 0.0E+00
  rlbl = 4
  job = 1
  return

 4    continue

  alpha = rho / sdot(n,work(1,rtld),1,work(1,v),1)
!
!  Early check for tolerance.
!
  call saxpy(n,-alpha, work(1,v),1,work(1,r),1)

  do i = 1, n
    work(i,s) = work(i,r)
  end do

  if ( snrm2(n,work(1,s),1) <= tol ) then
    call saxpy(n,alpha, work(1,phat),1,x,1)
    resid = snrm2(n,work(1,s),1) / bnrm2
    go to 30
  end if
!
!  Compute stabilizer vector SHAT and scalar OMEGA.
!
  ndx1 = ((shat-1)*n)+1
  ndx2 = ((s   -1)*n)+1
  rlbl = 5
  job = 2
  return

 5    continue

  ndx1 = ((shat-1)*n)+1
  ndx2 = ((t   -1)*n)+1
  sclr1 = 1.0E+00
  sclr2 = 0.0E+00
  rlbl = 6
  job = 1
  return

 6    continue

  omega = sdot(n,work(1,t),1,work(1,s),1) / sdot(n,work(1,t),1,work(1,t),1)
!
!  Compute new solution approximation vector X.
!
  call saxpy(n,alpha, work(1,phat),1,x,1)
  call saxpy(n,omega, work(1,shat),1,x,1)
!
!  Compute residual R, check for tolerance.
!
  call saxpy(n,-omega, work(1,t),1,work(1,r),1)
  ndx1 = need1
  ndx2 = need2
  rlbl = 7
  job = 4
  return

 7    continue

  if ( info == 1 ) then
    go to 30
  end if

  if ( iter == maxit ) then
    info = 1
    go to 20
  end if

  if ( abs(omega) < omegatol ) then
    go to 25
  else
    rho1 = rho
    go to 10
  end if

   20 continue
!
!  Iteration fails.
!
  rlbl = -1
  job = -1
  return

   25 continue
!
!  Set breakdown flag.
!
  if ( abs ( rho ) < rhotol ) then
    info = -10
  else if ( abs ( omega ) < omegatol ) then
    info = -11
  end if

  rlbl = -1
  job = -1
  return

   30 continue
!
!  Iteration successful; return.
!
  info = 0
  rlbl = -1
  job = -1

  return
end
subroutine cg ( n, b, x, work, iter, resid, matvec, psolve, info, curpform )

!*****************************************************************************80
!
!! CG implements the Conjugate Gradient method.
!
!  Discussion:
!
!    Preconditioning is used.  The convergence test is:
!
!      norm(b-A*x) / norm(b) < TOL.
!
!  Reference:
!
!    Barrett, Berry, Chan, Demmel, Donato, Dongarra, Eijkhout, Pozo,
!      Romine, van der Vorst,
!    Templates for the Solution of Linear Systems: Building Blocks
!      for Iterative Methods,
!    SIAM, 1994.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real B(N), the right hand side vector.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Workspace, real WORK(N,*).
!    Workspace for residual, direction vector, etc.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, actual number of iterations performed.
!
!    Input/output, real RESID.
!    On input, the allowable convergence measure for norm(b-A*x) / norm(b).
!    On output, the final value of this measure.
!
!  MATVEC  (external subroutine)
!          The user must provide a subroutine to perform the
!          matrix-vector product
!
!               y : = alpha*A*x+beta*y,
!
!          where alpha and beta are scalars, x and y are vectors,
!          and A is a matrix. Vector x must remain unchanged.
!          The solution is over-written on vector y.
!
!          The call is:
!
!             CALL MATVEC(ALPHA, X, BETA, Y)
!
!          The matrix is passed into the routine in a common block.
!
!  PSOLVE  (external subroutine)
!          The user must provide a subroutine to perform the
!          preconditioner solve routine for the linear system
!
!               M*x = b,
!
!          where x and b are vectors, and M a matrix. Vector b must
!          remain unchanged.
!          The solution is over-written on vector x.
!
!          The call is:
!
!             CALL PSOLVE(X, B)
!
!          The preconditioner is passed into the routine in a common block.
!
!  INFO    (output) integer
!
  implicit none

  integer n

  real b(n)
  real bnrm2
  character ( len = 8 ) curpform
  integer job
  integer info
  integer iter
  integer ndx1
  integer ndx2
  real resid
  real sclr1
  real sclr2
  real snrm2
  real tol
  real x(n)
  real work(*)

  external matvec
  external psolve

  info = 0
!
!  Test the input parameters.
!
  if ( n < 1 ) then
    info = -1
    write ( *, * ) ' '
    write ( *, * ) 'CG - Fatal error!'
    write ( *, * ) '  N is less than 1.'
    write ( *, * ) '  N = ',n
    stop
  end if

  if ( iter <= 0 ) then
    info = -3
    return
  end if
!
!  Stop test may need some indexing info from REVCOM
!  use the init call to send the request across. REVCOM
!  will note these requests, and everytime it asks for
!  stop test to be done, it will provide the indexing info.
!
  ndx1 = 1
  ndx2 = -1
  tol = resid
  bnrm2 = snrm2 ( n, b, 1 )
!
!  First time call always init.
!
  job = 1

  do

    call cg_revcom ( n, b, x, work, iter, resid, info, ndx1, ndx2, sclr1, &
      sclr2, job )
!
!  -1: Terminate.
!
    if ( job == -1 ) then

      exit
!
!  1: Compute WORK(NDX2) = SCLR1 * A * WORK(NDX1) + SCLR2 * WORK(NDX2).
!
    else if ( job == 1 ) then

      call matvec ( sclr1, work(ndx1), sclr2, work(ndx2) )
!
!  2: Solve M * WORK(NDX1) = WORK(NDX2).
!
    else if ( job == 2 ) then

      call psolve ( n, work(ndx1), work(ndx2), curpform )
!
!  3: Compute WORK(NDX2) = SCLR1 * A * X + SCLR2 * WORK(NDX2).
!
    else if ( job == 3 ) then

      call matvec ( sclr1, x, sclr2, work(ndx2) )
!
!  4: Do a stopping test on the relative residual reduction.
!
    else if ( job == 4 ) then

      call stopb ( n, work(ndx1), bnrm2, resid, tol, info )

    end if

    job = 2

  end do

  return
end
subroutine cg_revcom ( n, b, x, work, iter, resid, info, ndx1, ndx2, &
  sclr1, sclr2, job )

!*****************************************************************************80
!
!! CG_REVCOM is controlled by CG using reverse communication.
!
!  Reference:
!
!    Barrett, Berry, Chan, Demmel, Donato, Dongarra, Eijkhout, Pozo,
!      Romine, van der Vorst,
!    Templates for the Solution of Linear Systems: Building Blocks
!      for Iterative Methods,
!    SIAM, 1994.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real B(N), the right hand side vector.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Workspace, real WORK(N,4).
!    Workspace for residual, direction vector, etc.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, actual number of iterations performed.
!
!    Input/output, real RESID.
!    On input, the allowable convergence measure for norm(b-A*x) / norm(b).
!    On output, the final value of this measure.
!
!  INFO    (output) integer
!
!        = 0: Successful exit. Iterated approximate solution returned.
!
!          >  0: Convergence to tolerance not achieved. This will be
!                set to the number of iterations performed.
!
!          <  0: Illegal input parameter.
!
!                   -1: matrix dimension N < 0
!                   -3: Maximum number of iterations ITER < = 0.
!                   -5: Erroneous NDX1/NDX2 in INIT call.
!                   -6: Erroneous RLBL.
!
!  NDX1    (input/output) integer.
!  NDX2    On entry in INIT call contain indices required by interface
!          level for stopping test.
!          All other times, used as output, to indicate indices into
!          WORK[] for the MATVEC, PSOLVE done by the interface level.
!
!  SCLR1   (output) real.
!  SCLR2   Used to pass the scalars used in MATVEC.
!
!  JOB    (input/output) integer.
!          Used to communicate job code between the two levels.
!
  implicit none

  integer n

  real alpha
  real b(n)
  real beta
  real sdot
  real snrm2
  integer i
  integer job
  integer info
  integer iter
  integer maxit
  integer ndx1
  integer ndx2
  integer need1
  integer need2
  integer p
  integer q
  integer r
  real resid
  real rho
  real :: rho1 = 0.0E+00
  integer :: rlbl = 0
  real sclr1
  real sclr2
  real tol
  real work(n*4)
  real x(n)
  integer z

  save

  if ( job == 1 ) then
     go to 1
  else if ( job == 2 ) then
     if ( rlbl == 2) go to 2
     if ( rlbl == 3) go to 3
     if ( rlbl == 4) go to 4
     if ( rlbl == 5) go to 5
     info = -6
     go to 20
  end if

 1    continue

  info = 0
  maxit = iter
  tol = resid
!
!  Alias workspace columns.
!
  r = 1
  z = 2
  p = 3
  q = 4
!
!  Check if caller will need indexing info.
!
  if ( ndx1 == -1 ) then
    need1 = ndx1
  else if ( ndx1 == 1 ) then
    need1 = ((r-1)*n)+1
  else if ( ndx1 == 2 ) then
    need1 = ((z-1)*n)+1
  else if ( ndx1 == 3 ) then
    need1 = ((p-1)*n)+1
  else if ( ndx1 == 4 ) then
    need1 = ((q-1)*n)+1
  else
    info = -5
    go to 20
  end if

  if ( ndx2 == -1 ) then
    need2 = ndx2
  else if ( ndx2 == 1 ) then
    need2 = ((r-1)*n)+1
  else if ( ndx2 == 2 ) then
    need2 = ((z-1)*n)+1
  else if ( ndx2 == 3 ) then
    need2 = ((p-1)*n)+1
  else if ( ndx2 == 4 ) then
    need2 = ((q-1)*n)+1
  else
    info = -5
    go to 20
  end if
!
!  Set initial residual.
!
  do i = 1, n
    work(i+(r-1)*n) = b(i)
  end do

  if ( snrm2(n,x,1) /= 0.0E+00 ) then
    sclr1 = -1.0E+00
    sclr2 = 1.0E+00
    ndx1 = -1
    ndx2 = ((r-1)*n)+1
    rlbl = 2
    job = 3
    return
  end if

2     continue

  if ( snrm2(n,work(1:n+(r-1)*n),1) < tol ) then
    go to 30
  end if

  iter = 0

   10 continue
!
!  Perform preconditioned conjugate gradient iteration.
!
  iter = iter+1
!
!  Preconditioner Solve.
!
  ndx1 = ((z-1)*n)+1
  ndx2 = ((r-1)*n)+1
  rlbl = 3
  job = 2
  return

 3    continue

  rho = sdot(n,work(1:n+(r-1)*n),1,work(1:n+(z-1)*n),1)
!
!  Compute direction vector P.
!
  if ( iter > 1 ) then
    beta = rho / rho1
    call saxpy(n,beta, work(1:n+(p-1)*n),1,work(1:n+(z-1)*n),1)
  end if

  work(1:n+(p-1)*n) = work(1:n+(z-1)*n)
!
!  Compute scalar ALPHA (save A*P to Q).
!
  ndx1 = ((p-1)*n)+1
  ndx2 = ((q-1)*n)+1
  sclr1 = 1.0E+00
  sclr2 = 0.0E+00
  rlbl = 4
  job = 1
  return

 4    continue

  alpha = rho / sdot(n,work(1:n+(p-1)*n),1,work(1:n+(q-1)*n),1)
!
!  Compute current solution vector X.
!
  call saxpy(n,alpha, work(1:n+(p-1)*n),1,x,1)
!
!  Compute residual vector R, find norm,
!  then check for tolerance.
!
  call saxpy(n,-alpha,  work(1:n+(q-1)*n),1,work(1:n+(r-1)*n),1)

  ndx1 = need1
  ndx2 = need2
  rlbl = 5
  job = 4
  return

 5    continue

  if ( info == 1) then
    go to 30
  end if

  if ( iter == maxit ) then
    info = 1
    go to 20
  end if

  rho1 = rho

  go to 10

   20 continue
!
!  Iteration fails.
!
  rlbl = -1
  job = -1
  return

   30 continue
!
!  Iteration successful; return.
!
  info = 0
  rlbl = -1
  job = -1

  return
end
subroutine cgs ( n, b, x, work, iter, resid, matvec, psolve, info, &
  curpform )

!*****************************************************************************80
!
!! CGS implements the Conjugate Gradient Squared method.
!
!  Discussion:
!
!    The routine seeks a solution vector X satisfying the linear system
!
!      A * X = B.
!
!    Preconditioning is used.
!
!    The iteration is judged to have converged if:
!
!      norm(b-A*x) / norm(b) < TOL.
!
!    For this particular algorithm, an initial guess too close to
!    the actual solution can result in divergence.
!
!  Reference:
!
!    Barrett, Berry, Chan, Demmel, Donato, Dongarra, Eijkhout, Pozo, Romine,
!      van der Vorst,
!    Templates for the Solution of Linear Systems: Building Blocks
!      for Iterative Methods,
!    SIAM, 1994.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real B(N), the right hand side.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Workspace, real WORK(N,7).
!    Workspace for residual, direction vector, etc.
!    Note that vectors PHAT and QHAT, and UHAT and VHAT share
!    the same workspace.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, the actual number of iterations performed.
!
!    Input/output, real RESID.
!    On input, the allowable convergence measure for norm(b-A*x) / norm(b).
!    On output, the final value of this measure.
!
!    Input, external MATVEC.
!    The user must provide a subroutine to perform the matrix-vector product
!
!      y : = alpha*A*x+beta*y,
!
!    where alpha and beta are scalars, x and y are vectors,
!    and A is a matrix. Vector x must remain unchanged.
!    The solution is over-written on vector y.
!
!    The call is:
!
!      CALL MATVEC ( ALPHA, X, BETA, Y )
!
!    The matrix is passed into the routine in a common block.
!
!    Input, external PSOLVE.
!    The user must provide a subroutine to perform the
!    preconditioner solve routine for the linear system
!
!      M*x = b,
!
!    where x and b are vectors, and M a matrix. Vector b must remain unchanged.
!    The solution is over-written on vector x.  The call is:
!
!      CALL PSOLVE(X, B)
!
!    The preconditioner is passed into the routine in a common block.
!
!    Output, integer INFO.
!
!    = 0: Successful exit.
!    > 0: Convergence not achieved. This will be set
!         to the number of iterations performed.
!    < 0: Illegal input parameter.
!    -1:  matrix dimension N < 0
!    -3:  Maximum number of iterations ITER < = 0.
!    -10: RHO < BREAKTOL: RHO and RTLD have become orthogonal.
!
  implicit none

  integer n

  real b(n)
  real bnrm2
  character ( len = 8 ) curpform
  integer job
  integer info
  integer iter
  integer ndx1
  integer ndx2
  real resid
  real sclr1
  real sclr2
  real snrm2
  real tol
  real work(n*7)
  real x(n)

  external matvec
  external psolve

  info = 0
!
!  Test the input parameters.
!
  if ( n < 1 ) then
    info = -1
    write ( *, * ) ' '
    write ( *, * ) 'CGS - Fatal error!'
    write ( *, * ) '  N is less than 1.'
    write ( *, * ) '  N = ', n
    stop
  end if

  if ( iter <= 0 ) then
    info = -3
    return
  end if
!
!  Initialization.
!
  ndx1 = 1
  ndx2 = -1
  tol = resid
  bnrm2 = snrm2 ( n, b, 1 )
!
!  Set JOB so that the problem is initialized.
!
  job = 1

  do

    call cgs_revcom ( n, b, x, work, iter, resid, info, ndx1, ndx2, &
      sclr1, sclr2, job )
!
!  -1: Terminate.
!
    if ( job == -1 ) then

      exit
!
!  1: Compute WORK(NDX2) = SCLR1 * A * WORK(NDX1) + SCLR2 * WORK(NDX2).
!
    else if ( job == 1 ) then

      call matvec ( sclr1, work(ndx1), sclr2, work(ndx2) )
!
!  2: Solve M * WORK(NDX1) = WORK(NDX2).
!
    else if ( job == 2 ) then

      call psolve ( n, work(ndx1), work(ndx2), curpform )
!
!  3: Compute WORK(NDX2) = SCLR1 * A * X + SCLR2 * WORK(NDX2).
!
    else if ( job == 3 ) then

      call matvec ( sclr1, x, sclr2, work(ndx2) )
!
!  4: Do a stopping test on the relative residual reduction.
!
    else if ( job == 4 ) then

      call stopb ( n, work(ndx1), bnrm2, resid, tol, info )

    end if
!
!  Set JOB so that problem continues.
!
    job = 2

  end do

  return
end
subroutine cgs_ge ( n, b, x, work, iter, resid, info, a, p )

!*****************************************************************************80
!
!! CGS_GE implements the Conjugate Gradient Squared method for GE matrices.
!
!  Discussion:
!
!    The routine seeks a solution vector X satisfying the linear system
!
!      A * X = B.
!
!    Preconditioning is used.
!
!    The iteration is judged to have converged if:
!
!      norm(b-A*x) / norm(b) < TOL.
!
!    For this particular algorithm, an initial guess too close to
!    the actual solution can result in divergence.
!
!  Reference:
!
!    Barrett, Berry, Chan, Demmel, Donato, Dongarra, Eijkhout, Pozo, Romine,
!      van der Vorst,
!    Templates for the Solution of Linear Systems: Building Blocks
!      for Iterative Methods,
!    SIAM, 1994.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real B(N), the right hand side.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Workspace, real WORK(N,7).
!    Workspace for residual, direction vector, etc.
!    Note that vectors PHAT and QHAT, and UHAT and VHAT share
!    the same workspace.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, the actual number of iterations performed.
!
!    Input/output, real RESID.
!    On input, the allowable convergence measure for norm(b-A*x) / norm(b).
!    On output, the final value of this measure.
!
!    Output, integer INFO.
!
!    = 0: Successful exit.
!    > 0: Convergence not achieved. This will be set
!         to the number of iterations performed.
!    < 0: Illegal input parameter.
!    -1:  matrix dimension N < 0
!    -3:  Maximum number of iterations ITER < = 0.
!    -10: RHO < BREAKTOL: RHO and RTLD have become orthogonal.
!
  implicit none

  integer n

  real a(n,n)
  real b(n)
  real bnrm2
  integer job
  integer info
  integer iter
  integer ndx1
  integer ndx2
  real p(n,n)
  real resid
  real sclr1
  real sclr2
  real snrm2
  real tol
  real work(n*7)
  real x(n)

  info = 0
!
!  Test the input parameters.
!
  if ( n < 1 ) then
    info = -1
    write ( *, * ) ' '
    write ( *, * ) 'CGS - Fatal error!'
    write ( *, * ) '  N is less than 1.'
    write ( *, * ) '  N = ', n
    stop
  end if

  if ( iter <= 0 ) then
    info = -3
    return
  end if
!
!  Initialization.
!
  ndx1 = 1
  ndx2 = -1
  tol = resid
  bnrm2 = snrm2 ( n, b, 1 )
!
!  Set JOB so that the problem is initialized.
!
  job = 1

  do

    call cgs_revcom ( n, b, x, work, iter, resid, info, ndx1, ndx2, &
      sclr1, sclr2, job )
!
!  -1: Terminate.
!
    if ( job == -1 ) then

      exit
!
!  1: Compute WORK(NDX2) = SCLR1 * A * WORK(NDX1) + SCLR2 * WORK(NDX2).
!
    else if ( job == 1 ) then

      call matvec_ge ( n, a, sclr1, work(ndx1), sclr2, work(ndx2), &
        work(ndx2) )
!
!  2: Solve M * WORK(NDX1) = WORK(NDX2).
!  THIS IS THE PART I HAVE TO FIX UP NOW.
!
    else if ( job == 2 ) then

      call psolve_none ( n, work(ndx1), work(ndx2) )
!
!  3: Compute WORK(NDX2) = SCLR1 * A * X + SCLR2 * WORK(NDX2).
!
    else if ( job == 3 ) then

      call matvec_ge ( n, a, sclr1, x, sclr2, work(ndx2), work(ndx2) )
!
!  4: Do a stopping test on the relative residual reduction.
!
    else if ( job == 4 ) then

      call stopb ( n, work(ndx1), bnrm2, resid, tol, info )

    end if
!
!  Set JOB so that problem continues.
!
    job = 2

  end do

  return
end
subroutine cgs_revcom ( n, b, x, work, iter, resid, info, ndx1, ndx2, &
  sclr1, sclr2, job )

!*****************************************************************************80
!
!! CGS_REVCOM is controlled by CGS using reverse communication.
!
!  Reference:
!
!    Barrett, Berry, Chan, Demmel, Donato, Dongarra, Eijkhout, Pozo,
!      Romine, van der Vorst,
!    Templates for the Solution of Linear Systems: Building Blocks
!      for Iterative Methods,
!    SIAM, 1994.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real B(N), the right hand side vector.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Workspace, real WORK(N,7).
!    Workspace for residual, direction vector, etc.
!    Note that vectors PHAT and QHAT, and UHAT and VHAT share
!    the same workspace.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, actual number of iterations performed.
!
!    Input/output, real RESID.
!    On input, the allowable convergence measure for norm(b-A*x) / norm(b).
!    On output, the final value of this measure.
!
!  INFO    (output) integer
!
!        = 0: Successful exit.
!          >  0: Convergence not achieved. This will be set
!                to the number of iterations performed.
!
!          <  0: Illegal input parameter, or breakdown occured
!                during iteration.
!                   -1: matrix dimension N < 0
!                   -3: Maximum number of iterations ITER < = 0.
!                   -5: Erroneous NDX1/NDX2 in INIT call.
!                   -6: Erroneous RLBL.
!                   -10: RHO < BREAKTOL: RHO and RTLD have become orthogonal.
!
!  NDX1    (input/output) integer.
!  NDX2    On entry in INIT call contain indices required by interface
!          level for stopping test.
!          All other times, used as output, to indicate indices into
!          WORK[] for the MATVEC, PSOLVE done by the interface level.
!
!  SCLR1   (output) real.
!  SCLR2   Used to pass the scalars used in MATVEC.
!
!  JOB    (input/output) integer.
!          Used to communicate job code between the two levels.
!
  implicit none

  real, parameter :: one = 1.0E+00

  integer n

  real alpha
  real b(n)
  real beta
  real bnrm2
  real sdot
  real snrm2
  real getbreak
  integer i
  integer job
  integer info
  integer iter
  integer maxit
  integer ndx1
  integer ndx2
  integer need1
  integer need2
  integer p
  integer phat
  integer q
  integer qhat
  integer r
  real resid
  real rho
  real, save :: rho1 = 0.0E+00
  real rhotol
  integer, save :: rlbl = 0
  integer rtld
  real sclr1
  real sclr2
  real tol
  integer u
  integer uhat
  integer vhat
  real work(n,*)
  real x(n)

  save

  if ( job == 1 ) then
     go to 1
  else if ( job == 2 ) then
     if ( rlbl == 2) go to 2
     if ( rlbl == 3) go to 3
     if ( rlbl == 4) go to 4
     if ( rlbl == 5) go to 5
     if ( rlbl == 6) go to 6
     if ( rlbl == 7) go to 7
     info = -6
     go to 20
  end if

 1    continue

  info = 0
  maxit = iter
  tol = resid
!
!  Alias workspace columns.
!
  r  = 1
  rtld = 2
  p  = 3
  phat = 4
  q  = 5
  qhat = 6
  u  = 6
  uhat = 7
  vhat = 7
!
!  Check if caller will need indexing info.
!
  if ( ndx1 == -1 ) then
    need1 = ndx1
  else if ( ndx1 == 1 ) then
    need1 = ((r-1)*n)+1
  else if ( ndx1 == 2 ) then
    need1 = ((rtld-1)*n)+1
  else if ( ndx1 == 3 ) then
    need1 = ((p-1)*n)+1
  else if ( ndx1 == 4 ) then
    need1 = ((phat-1)*n)+1
  else if ( ndx1 == 5 ) then
    need1 = ((q-1)*n)+1
  else if ( ndx1 == 6 ) then
    need1 = ((qhat-1)*n)+1
  else if ( ndx1 == 7 ) then
    need1 = ((u-1)*n)+1
  else if ( ndx1 == 8 ) then
    need1 = ((uhat-1)*n)+1
  else if ( ndx1 == 9 ) then
    need1 = ((vhat-1)*n)+1
  else
    info = -5
    go to 20
  end if

  if ( ndx2 == -1 ) then
    need2 = ndx2
  else if ( ndx2 == 1 ) then
    need2 = ((r-1)*n)+1
  else if ( ndx2 == 2 ) then
    need2 = ((rtld-1)*n)+1
  else if ( ndx2 == 3 ) then
    need2 = ((p-1)*n)+1
  else if ( ndx2 == 4 ) then
    need2 = ((phat-1)*n)+1
  else if ( ndx2 == 5 ) then
    need2 = ((q-1)*n)+1
  else if ( ndx2 == 6 ) then
    need2 = ((qhat-1)*n)+1
  else if ( ndx2 == 7 ) then
    need2 = ((u-1)*n)+1
  else if ( ndx2 == 8 ) then
    need2 = ((uhat-1)*n)+1
  else if ( ndx2 == 9 ) then
    need2 = ((vhat-1)*n)+1
  else
    info = -5
    go to 20
  end if
!
!  Set the breakdown tolerance parameter.
!
  rhotol = getbreak()
!
!  Set initial residual.
!
  do i = 1, n
    work(i,r) = b(i)
  end do

  if ( snrm2(n,x,1) /= 0.0E+00 ) then
    sclr1 = -1.0E+00
    sclr2 = 1.0E+00
    ndx1 = -1
    ndx2 = ((r-1)*n)+1
    rlbl = 2
    job = 3
    return
  end if

 2    continue

  if ( snrm2(n,work(1,r),1) <= tol ) then
    go to 30
  end if

  bnrm2 = snrm2 ( n, b, 1 )

  if ( bnrm2 == 0.0E+00 ) then
    bnrm2 = 1.0E+00
  end if
!
!  Choose RTLD such that initially, (R,RTLD) = RHO is not equal to 0.
!  Here we choose RTLD = R.
!
  do i = 1, n
    work(i,rtld) = work(i,r)
  end do

  iter = 0

   10 continue
!
!  Perform Conjugate Gradient Squared iteration.
!
  iter = iter + 1
  rho = sdot(n,work(1,rtld),1,work(1,r),1)

  if ( abs(rho) < rhotol ) then
    go to 25
  end if
!
!  Compute direction vectors U and P.
!
  do i = 1, n
    work(i,u) = work(i,r)
  end do

  if ( iter > 1 ) then
!
!  Compute U.
!
    beta = rho / rho1

    call saxpy ( n, beta, work(1,q), 1, work(1,u), 1 )
!
!  Compute P.
!
    call sscal(n,beta**2, work(1,p),1)
    call saxpy(n,beta, work(1,q),1,work(1,p),1)
    call saxpy(n,one, work(1,u),1,work(1,p),1)

  else

    do i = 1, n
      work(i,p) = work(i,r)
    end do

  end if
!
!  Compute direction adjusting scalar ALPHA.
!
  ndx1 = ((phat-1)*n)+1
  ndx2 = ((p   -1)*n)+1
  rlbl = 3
  job = 2
  return

 3    continue

  ndx1 = ((phat-1)*n)+1
  ndx2 = ((vhat-1)*n)+1
  sclr1 = 1.0E+00
  sclr2 = 0.0E+00
  rlbl = 4
  job = 1
  return

 4    continue

  alpha = rho / sdot(n,work(1,rtld),1,work(1,vhat),1)

  do i = 1, n
    work(i,q) = work(i,u)
  end do

  call saxpy(n,-alpha, work(1,vhat),1,work(1,q),1)
!
!  Compute direction adjusting vectORT UHAT.
!  PHAT is being used as temporary storage here.
!
  do i = 1, n
    work(i,phat) = work(i,q)
  end do

  call saxpy(n,one, work(1,u),1,work(1,phat),1)
  ndx1 = ((uhat-1)*n)+1
  ndx2 = ((phat-1)*n)+1
  rlbl = 5
  job = 2
  return

 5    continue
!
!  Compute new solution approximation vector X.
!
  call saxpy(n,alpha, work(1,uhat),1,x,1)
!
!  Compute residual R and check for tolerance.
!
  ndx1 = ((uhat-1)*n)+1
  ndx2 = ((qhat-1)*n)+1
  sclr1 = 1.0E+00
  sclr2 = 0.0E+00
  rlbl = 6
  job = 1
  return

 6    continue

  call saxpy(n,-alpha, work(1,qhat),1,work(1,r),1)
  ndx1 = need1
  ndx2 = need2
  rlbl = 7
  job = 4
  return

 7    continue

  if ( info == 1) then
    go to 30
  end if

  if ( iter == maxit ) then
    info = 1
    go to 20
  end if

  rho1 = rho

  go to 10

   20 continue
!
!  Iteration fails.
!
  rlbl = -1
  job = -1
  return

   25 continue
!
!  Set the breakdown flag.
!
  if ( abs ( rho ) < rhotol ) then
    info = - 10
  end if

   30 continue
!
!  Iteration successful; return.
!
  info = 0
  rlbl = -1
  job = -1

  return
end
subroutine cheby ( n, a, b, x, work, iter, resid, matvec, psolve, &
  info, curpform )

!*****************************************************************************80
!
!! CHEBY implements the Chebyshev method.
!
!  Discussion:
!
!    Preconditioning is used.  This version requires explicit knowledge
!    of the maximum and minimum eigenvalues.  These eigenvalues must
!    be real and positive, which is the case for the symmetric positive
!    definite system.
!
!    The convergence test is:
!
!      norm(b-A*x) / norm(b) < TOL.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real B(N), the right hand side vector.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Workspace, real WORK(N,3).
!    Workspace for residual, direction vector, etc.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, actual number of iterations performed.
!
!    Input/output, real RESID.
!    On input, the allowable convergence measure for norm(b-A*x) / norm(b).
!    On output, the final value of this measure.
!
!  MATVEC  (external subroutine)
!          The user must provide a subroutine to perform the
!          matrix-vector product
!
!               y : = alpha*A*x+beta*y,
!
!          where alpha and beta are scalars, x and y are vectors,
!          and A is a matrix. Vector x must remain unchanged.
!          The solution is over-written on vector y.
!
!          The call is:
!
!             CALL MATVEC(ALPHA, X, BETA, Y)
!
!          The matrix is passed into the routine in a common block.
!
!  PSOLVE  (external subroutine)
!          The user must provide a subroutine to perform the
!          preconditioner solve routine for the linear system
!
!               M*x = b,
!
!          where x and b are vectors, and M a matrix. Vector b must
!          remain unchanged.
!          The solution is over-written on vector x.
!
!          The call is:
!
!             CALL PSOLVE(X, B)
!
!          The preconditioner is passed into the routine in a common block.
!
!  INFO    (output) integer
!
!        = 0: Successful exit. Iterated approximate solution returned.
!
!          >  0: Convergence to tolerance not achieved. This will be
!                set to the number of iterations performed.
!
!          <  0: Illegal input parameter, or breakdown occurred
!                during iteration.
!
!                Illegal parameter:
!
!                   -1: matrix dimension N < 0
!                   -3: Maximum number of iterations ITER < = 0.
!
  implicit none

  integer n

  real a(n,n)
  real b(n)
  real bnrm2
  character ( len = 8 ) curpform
  integer job
  integer info
  integer iter
  integer ndx1
  integer ndx2
  real resid
  real sclr1
  real sclr2
  real snrm2
  real tol
  real work(n*3)
  real x(n)

  external matvec
  external psolve

  info = 0
!
!  Test the input parameters.
!
  if ( n < 1 ) then
    info = -1
    write ( *, * ) ' '
    write ( *, * ) 'CHEBY - Fatal error!'
    write ( *, * ) '  N is less than 1.'
    write ( *, * ) '  N = ',n
    stop
  end if

  if ( iter <= 0 ) then
    info = -3
    return
  end if
!
!  Stop test may need some indexing info from REVCOM
!  use the init call to send the request across. REVCOM
!  will note these requests, and everytime it asks for
!  stop test to be done, it will provide the indexing info.
!
  ndx1 = 1
  ndx2 = -1
  tol = resid
  bnrm2 = snrm2 ( n, b, 1 )
!
!  First time call always init.
!
  job = 1

  do

    call cheby_revcom ( n, a, b, x, work, iter, resid, info, ndx1, &
      ndx2, sclr1, sclr2, job, curpform )
!
!  -1: Terminate.
!
    if ( job == -1 ) then

      exit
!
!  1: Compute WORK(NDX2) = SCLR1 * A * WORK(NDX1) + SCLR2 * WORK(NDX2).
!
    else if ( job == 1 ) then

      call matvec ( sclr1, work(ndx1), sclr2, work(ndx2) )
!
!  2: Solve M * WORK(NDX1) = WORK(NDX2).
!
    else if ( job == 2 ) then

      call psolve ( n, work(ndx1), work(ndx2), curpform )
!
!  3: Compute WORK(NDX2) = SCLR1 * A * X + SCLR2 * WORK(NDX2).
!
    else if ( job == 3 ) then

      call matvec ( sclr1, x, sclr2, work(ndx2) )
!
!  4: Do a stopping test on the relative residual reduction.
!
    else if ( job == 4 ) then

      call stopb ( n, work(ndx1), bnrm2, resid, tol, info )

    end if

    job = 2

  end do

  return
end
subroutine cheby_revcom ( n, a, b, x, work, iter, resid, info, &
  ndx1, ndx2, sclr1, sclr2, job, curpform )

!*****************************************************************************80
!
!! CHEBY_REVCOM is controlled by CHEBY using reverse communication.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real B(N), the right hand side vector.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Workspace, real WORK(N,3).
!    Workspace for residual, direction vector, etc.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, actual number of iterations performed.
!
!    Input/output, real RESID.
!    On input, the allowable convergence measure for norm(b-A*x) / norm(b).
!    On output, the final value of this measure.
!
!  INFO    (output) integer
!
!        = 0: Successful exit. Iterated approximate solution returned.
!
!          >  0: Convergence to tolerance not achieved. This will be
!                set to the number of iterations performed.
!
!          <  0: Illegal input parameter, or breakdown occurred
!                during iteration.
!
!                Illegal parameter:
!
!                   -1: matrix dimension N < 0
!                   -3: Maximum number of iterations ITER < = 0.
!                   -5: Erroneous NDX1/NDX2 in INIT call.
!                   -6: Erroneous RLBL.
!
!  NDX1    (input/output) integer.
!  NDX2    On entry in INIT call contain indices required by interface
!          level for stopping test.
!          All other times, used as output, to indicate indices into
!          WORK[] for the MATVEC, PSOLVE done by the interface level.
!
!  SCLR1   (output) real.
!  SCLR2   Used to pass the scalars used in MATVEC.
!
!  JOB    (input/output) integer.
!          Used to communicate job code between the two levels.
!
  implicit none

  integer n

  real a(n,n)
  real, save :: alpha = 0.0E+00
  real b(n)
  real beta
  real bnrm2
  real c
  character ( len = 8 ) curpform
  real d
  real snrm2
  real eigmax
  real eigmin
  integer i
  integer info
  integer iter
  integer job
  integer maxit
  integer ndx1
  integer ndx2
  integer need1
  integer need2
  integer p
  integer r
  real resid
  integer, save :: rlbl = 0
  real sclr1
  real sclr2
  real tol
  real work(n,*)
  real x(n)
  integer z

  save

  if ( job == 1 ) then
    go to 1
  else if ( job == 2 ) then
    if ( rlbl == 2) go to 2
    if ( rlbl == 3) go to 3
    if ( rlbl == 4) go to 4
    if ( rlbl == 5) go to 5
    info = -6
    go to 20
  end if

 1    continue

  info = 0

  maxit = iter
  tol = resid
!
!  Get the extremal eigenvalues of A.
!
  call geteig ( n, a, work, eigmax, eigmin, curpform )
!
!  Alias workspace columns.
!
  r = 1
  p = 2
  z = 3
!
!  Check if caller will need indexing info.
!
  if ( ndx1 == -1 ) then
    need1 = ndx1
  else if ( ndx1 == 1 ) then
    need1 = ((r-1)*n)+1
  else if ( ndx1 == 2 ) then
    need1 = ((p-1)*n)+1
  else if ( ndx1 == 3 ) then
    need1 = ((z-1)*n)+1
  else
    info = -5
    go to 20
  end if

  if ( ndx2 == -1 ) then
    need2 = ndx2
  else if ( ndx2 == 1 ) then
    need2 = ((r-1)*n)+1
  else if ( ndx2 == 2 ) then
    need2 = ((p-1)*n)+1
  else if ( ndx2 == 3 ) then
    need2 = ((z-1)*n)+1
  else
    info = -5
    go to 20
  end if
!
!  Set initial residual.
!
  do i = 1, n
    work(i,r) = b(i)
  end do

  if ( snrm2(n,x,1) /= 0.0E+00 ) then
    ndx1 = -1
    ndx2 = ((r-1)*n)+1
    sclr1 = -1.0E+00
    sclr2 = 1.0E+00
    rlbl = 2
    job = 3
    return
  end if

 2    continue

  if ( snrm2(n,work(1,r),1) < tol ) then
    go to 30
  end if

  bnrm2 = snrm2(n,b,1)
  if ( bnrm2 == 0.0E+00 ) then
    bnrm2 = 1.0E+00
  end if
!
!  Initialize ellipse parameters.
!
  c = ( eigmax - eigmin ) / 2.0E+00
  d = ( eigmax + eigmin ) / 2.0E+00

  iter = 0

   10 continue
!
!  Perform Chebyshev iteration.
!
  iter = iter+1

  ndx1 = ((z-1)*n)+1
  ndx2 = ((r-1)*n)+1
  rlbl = 3
  job = 2
  return

 3    continue

  if ( iter > 1 ) then
    beta = ( ( c * alpha ) / 2.0E+00 )**2
    alpha = 1.0E+00 / (d-beta)
    call saxpy(n,beta, work(1,p),1,work(1,z),1)
  else
    alpha = 2.0E+00 / d
  end if

  do i = 1, n
    work(i,p) = work(i,z)
  end do
!
!  Compute new approximation vector X; check accuracy.
!
  call saxpy(n,alpha, work(1,p),1,x,1)
  ndx1 = ((p-1)*n)+1
  ndx2 = ((r-1)*n)+1
  sclr1 = -alpha
  sclr2 = 1.0E+00
  rlbl = 4
  job = 1
  return

 4    continue

  ndx1 = need1
  ndx2 = need2
  rlbl = 5
  job = 4
  return

 5    continue

  if ( info == 1 ) then
    go to 30
  end if

  if ( iter == maxit ) then
    info = 1
    go to 20
  end if

  go to 10

   20 continue
!
!  Iteration fails.
!
  rlbl = -1
  job = -1
  return

   30 continue
!
!  Iteration successful; return.
!
  info = 0
  rlbl = -1
  job = -1

  return
end
subroutine dif ( m, n, a )

!*****************************************************************************80
!
!! DIF returns the second difference matrix.
!
!  Example:
!
!    N = 5
!
!    2 -1  .  .  .
!   -1  2 -1  .  .
!    . -1  2 -1  .
!    .  . -1  2 -1
!    .  .  . -1  2
!
!  Rectangular Properties:
!
!    A is tridiagonal.
!
!    Because A is tridiagonal, it has property A (bipartite).
!
!    A is integral: int ( A ) = A.
!
!    A is Toeplitz: constant along diagonals.
!
!    A has property A (bipartite).
!
!  Square Properties:
!
!    A is symmetric: A' = A.
!
!    Because A is symmetric, it is normal.
!
!    Because A is normal, it is diagonalizable.
!
!    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
!
!    A is positive definite.
!
!    A is an M matrix.
!
!    A is weakly diagonally dominant, but not strictly diagonally dominant.
!
!    A has an LU factorization A = L * U, without pivoting.
!
!      The matrix L is lower bidiagonal with subdiagonal elements:
!
!        L(I+1,I) = -I/(I+1)
!
!      The matrix U is upper bidiagonal, with diagonal elements
!
!        U(I,I) = (I+1)/I
!
!      and superdiagonal elements which are all -1.
!
!    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
!
!      L(I,I) =    sqrt ( (I+1) / I )
!      L(I,I-1) = -sqrt ( (I-1) / I )
!
!    The eigenvalues are
!
!      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
!                = 4 SIN**2(I*PI/(2*N+2))
!
!    The corresponding eigenvector X(I) has entries
!
!       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
!
!    Simple linear systems:
!
!      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
!
!      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
!
!    det ( A ) = N + 1.
!
!    The value of the determinant can be seen by induction,
!    and expanding the determinant across the first row:
!
!      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
!                = 2 * N - (N-1)
!                = N + 1
!
!  Reference:
!
!    Robert Gregory and David Karney,
!    Example 3.18,
!    A Collection of Matrices for Testing Computational Algorithms,
!    Wiley, New York, 1969, page 45, QA263 G862.
!
!    Morris Newman and John Todd,
!    Example A8,
!    The evaluation of matrix inversion programs,
!    Journal of the Society for Industrial and Applied Mathematics,
!    Volume 6, Number 4, pages 466-476, 1958.
!
!    John Todd,
!    Example A8,
!    Basic Numerical Mathematics,
!    Volume 2: Numerical Algebra,
!    Birkhauser, Basel and Academic Press, New York, 1977, page 1.
!
!    Joan Westlake,
!    Test Matrix A15,
!    A Handbook of Numerical Matrix Inversion and Solution of Linear Equations,
!    John Wiley, 1968.
!
!  Modified:
!
!    01 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns of A.
!
!    Output, real A(M,N), the matrix.
!
  implicit none

  integer m
  integer n

  real a(m,n)
  integer i
  integer j

  do i = 1, m
    do j = 1, n

      if ( j == i-1 ) then
        a(i,j) = - 1.0E+00
      else if ( j == i ) then
        a(i,j) = 2.0E+00
      else if ( j == i+1 ) then
        a(i,j) = - 1.0E+00
      else
        a(i,j) = 0.0E+00
      end if

    end do
  end do

  return
end
subroutine elem_vec ( i, n, alpha, e )

!*****************************************************************************80
!
!! ELEM_VEC constructs the I-th elementary vector E, scaled by ALPHA.
!
!  Definition:
!
!    In R^N, the I-th elementary vector has a 1 in entry I, and 0's
!    in all other entries.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer I, specifies the location of the nonzero entry.
!
!    Input, integer N, the order of the vector.
!
!    Input, real ALPHA, the value of the I-th entry of the output vector.
!
!    Output, real E(N), ALPHA times the I-th elementary vector.
!
  implicit none

  integer n

  real alpha
  real e(n)
  integer i

  e(1:n) = 0.0E+00
  e(i) = alpha

  return
end
function getbreak ( )

!*****************************************************************************80
!
!! GETBREAK is supposed to allow the user to set certain tolerances.
!
!  Discussion:
!
!    As sketched, the routine is quite crude.  Many routines call
!    GETBREAK for values, but no attempt is made to signal which
!    value is being asked for.
!
!  Modified:
!
!    28 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real GETBREAK, ...
!
  implicit none

  real getbreak

  getbreak = 0.00001

  return
end
subroutine geteig_ge ( n, a, work, eigmax, eigmin, curpform )

!*****************************************************************************80
!
!! GETEIG_GE computes the eigenvalues of the iteration matrix for the GE format.
!
!  Discussion:
!
!    The GE format is the LINPACK/LAPACK general matrix storage mode,
!    in which a full N by N matrix is stored in an N by N array.
!
!    GETEIG_GE uses an LAPACK routine for computing all the
!    eigenvalues of the matrix A.
!
!    This is for testing only, as this is more expensive than a direct
!    dense solver for the linear system.
!
!  Modified:
!
!    24 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
  implicit none

  integer n

  real a(n,n)
  character ( len = 8 ) curpform
  real slamch
  real eigmax
  real eigmin
  integer i
  integer info
  integer j
  logical lsamen
  real m(1)
  real matnorm
  real work(n,n+2)

  common /matpre/ m

  save /matpre/
!
!  As the matrix A is overwritten in the following routine, we
!  copy it to temporary workspace.
!
  work(1:n,1:n) = a(1:n,1:n)

  if ( lsamen ( 3, curpform, 'jacobi' ) ) then
    do i = 1, n
      work(i,i) = work(i,i) / m(i)
    end do
  end if
!
!  Call LAPACK eigenvalue routine.
!
  call ssyev ( 'no_vec', 'upper', n, work, n, work(1,n+1), work(1,n+2), &
    (3*n)-1, info )

  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GETEIG_GE - Warning:'
    write ( *, * ) '  SSYEV could not compute all eigenvalues.'
    write ( *, * ) '  Setting eigmin/max to default values.'
    eigmin = slamch ( 'eps' )
    eigmax = matnorm ( n, a )
    return
  end if

  eigmin = work(1,n+1)
  eigmax = work(n,n+1)
!
!  The eigenvalues should be positive.
!
  if ( eigmin < 0.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GETEIG_GE - Warning:'
    write ( *, * ) '  Computed min eigenvalue < = 0: set to epsilon'
    eigmin = slamch ( 'eps' )
  end if

  if ( eigmax < eigmin ) then
    write ( *, * ) ' '
    write ( *, * ) 'GETEIG_GE - Warning:'
    write ( *, * ) '  Max eigenvalue < min: set to |A|'
    eigmax = matnorm ( n, a )
  end if

  return
end
subroutine givens_apply ( i, h, givens, ldg )

!*****************************************************************************80
!
!! GIVENS_APPLY applies a sequence of Givens rotations to a column of H.
!
!  Discussion:
!
!    The Givens parameters are stored so that the first I-2 Givens
!    rotation matrices are known.  The I-1st Givens rotation is computed
!    using BLAS 1 routine SROTG.  Each rotation is applied to the 2x1 vector
!
!      [ H(J), H(J+1) ]',
!
!    which results in
!
!      H(J+1) = 0.
!
!  Reference:
!
!    Barrett, Berry, Chan, Demmel, Donato, Dongarra, Eijkhout, Pozo,
!      Romine, van der Vorst,
!    Templates for the Solution of Linear Systems: Building Blocks
!      for Iterative Methods,
!    SIAM, 1994.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer I, indicates that rotation matrix I-1 is to be constructed.
!
!    Input/output, real H(*), the I-th column of the matrix.
!
!    Output, real GIVENS(LDG,2), contains a sequence of pairs of parameters
!    definining Givens rotations.  GIVENS(I,1) and GIVENS(I,2) represent
!    the cosine and sine of the I-th Givens rotation.
!
!    Input, integer LDG, the leading dimension of GIVENS.
!
  implicit none

  integer ldg

  real givens(ldg,2)
  real h(*)
  integer i
  integer j
!
!  Construct the I-1st rotation matrix.
!
  call givens_set ( h(i), h(i+1), givens(i,1), givens(i,2) )
!
!  Apply 1,...,I-1st rotation matrices to the I-th column of H.
!
  do j = 1, i-1
    call rotvec ( h(j), h(j+1), givens(i,1), givens(i,2) )
  end do

  return
end
subroutine givens_set ( a, b, c, s )

!*****************************************************************************80
!
!! GIVENS_SET computes Givens rotation parameters.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, real A, B, ...
!
!    Output, real C, S, ...
!
  implicit none

  real a
  real b
  real c
  real s
  real temp

  if ( b == 0.0E+00  ) then
    c = 1.0E+00
    s = 0.0E+00
  else if ( abs ( b ) > abs ( a ) ) then
    temp = - a / b
    s = 1.0E+00 / sqrt ( 1.0 + temp**2 )
    c = temp * s
  else
    temp = - b / a
    c = 1.0E+00 / sqrt ( 1.0 + temp**2 )
    s = temp * c
  end if

  return
end
subroutine gmres ( n, b, x, restrt, iter, resid, matvec, psolve, info, &
  curpform )

!*****************************************************************************80
!
!! GMRES implements the Generalized Minimal Residual method.
!
!  Discussion:
!
!    The algorithm used is the Generalized Minimal Residual iterative
!    method with preconditioning.
!
!    The convergence test is:
!
!      norm(b-A*x) / norm(b) < TOL.
!
!  Modified:
!
!    10 June 2004
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real B(N), the right hand side vector.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Input, integer RESTRT.
!    Restart parameter, <= N. This parameter controls the amount
!    of memory required for matrix WORK2.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, the actual number of iterations performed.
!
!    Input/output, real RESID.
!    On input, the allowable error tolerance.
!    On output, the norm of the residual vector if solution approximated
!    to tolerance, otherwise reset to input tolerance.
!
!    Input, external MATVEC, the name of the matrix-vector multiply routine.
!    The user must provide a subroutine to perform the matrix-vector
!    product A*x = y.  Vector x must remain unchanged. The solution is
!    overwritten on vector y.  The call is:
!
!      CALL MATVEC ( X, Y )
!
!    Input, external PSOLVE, the name of the preconditioning routine.
!
!    Output, integer INFO, an error indicator.
!    0, successful exit;
!    1, the maximum number of iterations were performed, but
!       convergence was not achieved.
!
!    Input, character*8 CURPFORM, ?
!
!  Local parameters:
!
!    Workspace, real WORK(N,5).
!    Note that if the initial guess is the zero vector, then
!    storing the initial residual is not necessary.
!
!    Workspace, real WORK2(RESTRT,2*RESTRT+2).
!    This workspace is used for constructing and storing the
!    upper Hessenberg matrix. The two extra columns are used to
!    store the Givens rotation matrices.
!
  implicit none

  integer restrt
  integer n

  real b(n)
  real bnrm2
  character ( len = 8 ) curpform
  integer job
  integer info
  integer iter
  external matvec
  integer ndx1
  integer ndx2
  external psolve
  real resid
  real sclr1
  real sclr2
  real snrm2
  real tol
  real work(n*5)
  real work2(restrt*(2*restrt+2))
  real x(n)

  info = 0
!
!  Check the input parameters.
!
  if ( n < 1 ) then
    info = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GMRES - Fatal error!'
    write ( *, '(a)' ) '  N is less than 1.'
    write ( *, '(a,i6)' ) '  N = ', n
    stop
  end if

  if ( iter <= 0 ) then
    info = -3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GMRES - Fatal error!'
    write ( *, '(a)' ) '  ITER is less than 1.'
    write ( *, '(a,i6)' ) '  ITER = ', iter
    stop
  end if
!
!  Stop test may need some indexing info from REVCOM.
!
!  Use the "INIT" call to send the request across.
!
!  REVCOM will note these requests, and everytime it asks for
!  stop test to be done, it will provide the indexing info.
!
!  Note: V, and GIV contain # of history vectors each.
!  To access the i'th vector in V, add i to V*OFSET 1< = i<=RESTRT
!  To access the i'th vector in GIV, add i to GIV*OFSET 1< = i<=RESTRT
!
!  1 == R;
!  2 == S;
!  3 == W;
!  4 == Y;
!  5 == AV;
!  6 == H;
!  7*OFSET+i = = V;
!  8*OFSET+i = =GIV;
!  -1 = = ignore;
!  any other = = error
!
  ndx1 = 1
  ndx2 = -1
  tol = resid
  bnrm2 = snrm2 ( n, b, 1 )
!
!  JOB = 1 on first call, to request initialization.
!
  job = 1

  do

    call gmres_revcom ( n, b, x, restrt, work, work2, iter, resid, &
      info, ndx1, ndx2, sclr1, sclr2, job )
!
!  -1: Terminate.
!
    if ( job == -1 ) then

      exit
!
!  1: Compute WORK(NDX2) = SCLR1 * A * X + SCLR2 * WORK(NDX2).
!
    else if ( job == 1 ) then

      call matvec ( sclr1, x, sclr2, work(ndx2) )
!
!  2: Solve M * WORK(NDX1) = WORK(NDX2).
!
    else if ( job == 2 ) then

      call psolve ( n, work(ndx1), work(ndx2), curpform )
!
!  3: Compute WORK(NDX2) = SCLR1 * A * WORK(NDX1) + SCLR2 * WORK(NDX2).
!
    else if ( job == 3 ) then

      call matvec ( sclr1, work2(ndx1), sclr2, work(ndx2) )
!
!  4: Do a stopping test on the relative residual reduction.
!
    else if ( job == 4 ) then

      call stopb ( n, work(ndx1), bnrm2, resid, tol, info )

    end if

    job = 2

  end do

  return
end
subroutine gmres_revcom ( n, b, x, restrt, work, work2, iter, &
  resid, info, ndx1, ndx2, sclr1, sclr2, job)

!*****************************************************************************80
!
!! GMRES_REVCOM is controlled by GMRES using reverse communication.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real B(N), the right hand side vector.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Input, integer RESTRT.
!    Restart parameter, < = N. This parameter controls the amount
!    of memory required for matrix WORK2.
!
!    Workspace, real WORK(N,5).
!    Note that if the initial guess is the zero vector, then
!    storing the initial residual is not necessary.
!
!    Workspace, real WORK2(RESTRT,2*RESTRT+2).
!    This workspace is used for constructing and storing the
!    upper Hessenberg matrix. The two extra columns are used to
!    store the Givens rotation matrices.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, actual number of iterations performed.
!
!    Input/output, real RESID.
!    On input, the allowable error tolerance.
!    On output, the norm of the residual vector if solution approximated
!    to tolerance, otherwise reset to input tolerance.
!
!  INFO    (output) integer
!        = 0:  successful exit
!        = 1:  maximum number of iterations performed;
!                 convergence not achieved.
!            -5: Erroneous NDX1/NDX2 in INIT call.
!            -6: Erroneous RLBL.
!
!  NDX1    (input/output) integer.
!  NDX2    On entry in INIT call contain indices required by interface
!          level for stopping test.
!          All other times, used as output, to indicate indices into
!          WORK[] for the MATVEC, PSOLVE done by the interface level.
!
!  SCLR1   (output) real.
!  SCLR2   Used to pass the scalars used in MATVEC.
!
!  JOB    (input/output) integer.
!          Used to communicate job code between the two levels.
!
  implicit none

  integer, parameter :: ofset = 1000

  integer restrt
  integer n

  real approxres
  integer av
  real b(n)
  real  bnrm2
  real snrm2
  integer giv
  integer h
  integer i
  integer istep
  integer job
  integer info
  integer iter
  integer maxit
  integer ndx1
  integer ndx2
  integer need1
  integer need2
  integer r
  real resid
  integer, save :: rlbl = 0
  real rnorm
  integer s
  real sclr1
  real sclr2
  real tol
  integer v
  integer w
  real work(n,5)
  real work2(restrt,*)
  real x(n)
  integer y

  save

  if ( job == 1 ) then
     go to 1
  else if ( job == 2 ) then
     if ( rlbl == 2) go to 2
     if ( rlbl == 3) go to 3
     if ( rlbl == 4) go to 4
     if ( rlbl == 5) go to 5
     if ( rlbl == 6) go to 6
     if ( rlbl == 7) go to 7
     info = -6
     go to 200
  end if

 1    continue

  info = 0
  maxit = iter
  tol = resid
!
!  Alias workspace columns.
!
  r = 1
  s = 2
  w = 3
  y = 4
  av = 4
  h = 1
  v = h + restrt
  giv = v + restrt
!
!  Check if caller will need indexing info.
!
  if ( ndx1 == -1 ) then
    need1 = ndx1
  else if ( ndx1 == 1 ) then
    need1 = ((r-1)*n) + 1
  else if ( ndx1 == 2 ) then
    need1 = ((s-1)*n) + 1
  else if ( ndx1 == 3 ) then
    need1 = ((w-1)*n) + 1
  else if ( ndx1 == 4 ) then
    need1 = ((y-1)*n) + 1
  else if ( ndx1 == 5 ) then
    need1 = ((av-1)*n) + 1
  else if ( ndx1 == 6 ) then
    need1 = ((h-1)*n) + 1
  else if ( (ndx1 > v*ofset) .and. (ndx1 <= v*ofset+restrt) ) then
    need1 = ((ndx1-v*ofset-1)*n) + 1
  else if ( (ndx1 > giv*ofset) .and. (ndx1 <= giv*ofset+restrt) ) then
    need1 = ((ndx1-giv*ofset-1)*n) + 1
  else
    info = -5
    go to 100
  end if

  if ( ndx2==-1 ) then
    need2 = ndx2
  else if ( ndx2 == 1 ) then
        need2 = ((r-1)*n) + 1
     else if ( ndx2 == 2 ) then
        need2 = ((s-1)*n) + 1
     else if ( ndx2 == 3 ) then
        need2 = ((w-1)*n) + 1
     else if ( ndx2 == 4 ) then
        need2 = ((y-1)*n) + 1
     else if ( ndx2 == 5 ) then
        need2 = ((av-1)*n) + 1
     else if ( ndx2 == 6 ) then
        need2 = ((h-1)*n) + 1
     else if ( (ndx2 > v*ofset) .and. (ndx2 <= v*ofset+restrt) ) then
        need2 = ((ndx2-v*ofset-1)*n) + 1
     else if ( (ndx2 > giv*ofset) .and. (ndx2 <= giv*ofset+restrt) ) then
        need2 = ((ndx2-giv*ofset-1)*n) + 1
     else
        info = -5
        go to 100
     end if
!
!  Set initial residual.
!
  work(1:n,r) = b(1:n)

  if ( snrm2 ( n, x, 1 ) /= 0.0E+00 ) then
     sclr1 = -1.0E+00
     sclr2 = 1.0E+00
     ndx1 = -1
     ndx2 = ((r-1)*n) + 1
     rlbl = 2
     job = 1
     return
  end if

 2    continue

  if ( snrm2(n,work(1,r),1) < tol ) then
    go to 200
  end if

  bnrm2 = snrm2(n,b,1)

  if ( bnrm2 == 0.0E+00 ) then
    bnrm2 = 1.0E+00
  end if

   10 continue

  iter = iter + 1
!
!  Construct the first column of V, and initialize S to the
!  elementary vector E1 scaled by RNORM.
!
  ndx1 = ((v-1)*n) + 1
  ndx2 = ((r-1)*n) + 1
  rlbl = 3
  job = 2
  return

 3    continue

  rnorm = snrm2 ( n, work2(1,v), 1 )

  do i = 1, n
    work2(i,v) = work2(i,v) * rnorm
  end do

  call elem_vec ( 1, n, rnorm, work(1,s) )

  istep = 1

99    continue
    ndx1 = ((v+istep-1-1)*restrt) + 1
    ndx2 = ((av   -1)*restrt) + 1
    sclr1 = 1.0E+00
    sclr2 = 0.0E+00
    rlbl = 4
    job = 3
    return

 4      continue

    ndx1 = ((w -1)*n) + 1
    ndx2 = ((av-1)*n) + 1
    rlbl = 5
    job = 2
    return

 5      continue
!
!  Construct ISTEP-th column of H so that it is orthnormal to
!  the previous ISTEP-1 columns.
!
    call orthoh ( istep, n, work2(1,istep+h-1), work2(1,v), restrt, work(1,w) )
!
!  Apply Givens rotations to the ISTEP-th column of H. This
!  effectively reduces the Hessenberg matrix to upper
!  triangular form during the RESTRT iterations.
!
    if ( istep > 1 ) then
      call givens_apply ( istep, work2(1,istep+h-1), work2(1,giv), restrt )
    end if
!
!  Approximate residual norm. Check tolerance. If okay, compute
!  final approximation vector X and quit.
!
    resid = approxres ( istep, work2(1,istep+h-1), work(istep,s), &
      work2(1,giv), restrt) / bnrm2

    if ( resid <= tol ) then
      call update ( istep, n, x, work2(1,h), restrt, work(1,y), work(1,s), &
        work2(1,v), restrt )
      go to 200
    end if

  if ( istep < restrt ) then
    istep = istep + 1
    go to 99
  end if
!
!  Compute current solution vector X.
!
  call update ( restrt, n, x, work2(1,h), restrt, work(1,y), work(1,s), &
    work2(1,v), restrt )
!
!  Compute residual vector R, find norm,
!  then check for tolerance.
!
  work(1:n,r) = b(1:n)

  ndx1 = -1
  ndx2 = ((r-1)*n) + 1
  sclr1 = -1.0E+00
  sclr2 = 1.0E+00
  rlbl = 6
  job = 1
  return

 6    continue

  work(i+1,s) = snrm2(n,work(1,r),1)
  ndx1 = need1
  ndx2 = need2
  rlbl = 7
  job = 4
  return

 7    continue

  if ( info == 1) then
    go to 200
  end if

  if ( iter == maxit ) then
    info = 1
    go to 100
  end if

  go to 10

  100 continue
!
!  Iteration fails.
!
  rlbl = -1
  job = -1
  return

  200 continue
!
!  Iteration successful; return.
!
  info = 0
  rlbl = -1
  job = -1

  return
end
subroutine jacobi_ge ( n, a, b, x, iter, restol, resid, info )

!*****************************************************************************80
!
!! JACOBI_GE implements the Jacobi method for a matrix in GE format.
!
!  Discussion:
!
!    The GE format is the LINPACK/LAPACK general matrix storage mode,
!    in which a full N by N matrix is stored in an N by N array.
!
!    The relative error is measured:
!
!      norm ( X - X_1 ) / norm ( X ).
!
!  Modified:
!
!    12 December 2002
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(N,N), the system matrix, store in LINPACK GE form.
!
!    Input, real B(N), the right hand side vector.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, actual number of iterations performed.
!
!    Input, real RESTOL, the allowable convergence measure for
!    norm(x-x_1) / norm(x).
!
!    Output, real RESID, the final value of the convergence measure.
!
!    Output, integer INFO, error flag.
!    0: Successful exit. Iterated approximate solution returned.
!    >0: Convergence to tolerance not achieved. This will be
!        set to the number of iterations performed.
!
  implicit none

  integer n

  real a(n,n)
  real b(n)
  real diag(n)
  integer job
  integer info
  integer iter
  real resid
  real restol
  real sclr1
  real sclr2
  real temp(n)
  real x(n)
  real x1(n)
  real xnrm2

  info = 0
  resid = 0.0E+00
!
!  Test the input parameters.
!
  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JACOBI_GE - Fatal error!'
    write ( *, '(a)' ) '  Input quantity N < = 0!'
    stop
  end if

  if ( iter <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JACOBI_GE - Fatal error!'
    write ( *, '(a)' ) '  Maximum number of iterations < = 0.'
    write ( *, '(a,i6)' ) '  ITER = ', iter
    stop
  end if
!
!  Form the matrix splitting inv(M) and N.
!
  call jacobi_split_ge ( n, a, diag )
!
!  The first call is always "INIT", that is, JOB = 1.
!
  job = 1

  do

    call jacobi_revcom ( n, b, x, diag, temp, x1, iter, info, sclr1, &
      sclr2, job )
!
!  -1: Reconstruct the matrix, and terminate.
!
    if ( job == -1 ) then

      call jacobi_recon_ge ( n, a, diag )

      exit
!
!  1: Compute TEMP = SCLR1 * A * X + SCLR2 * TEMP.
!
    else if ( job == 1 ) then

      call matvec_ge ( n, a, sclr1, x, sclr2, temp, temp )
!
!  2: Do a stopping test on the current residual relative to the norm of X.
!
    else if ( job == 2 ) then

      call stopx ( n, x1, x, xnrm2, resid, restol, info )

    end if

    job = 2

  end do

  return
end
subroutine jacobi_revcom ( n, b, x, diag, temp, x1, iter, info, sclr1, sclr2, &
  job )

!*****************************************************************************80
!
!! JACOBI_REVCOM is controlled by JACOBI using reverse communication.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real B(N), the right hand side vector.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Input, real DIAG(N), the inverse of the diagonal entries of
!    the system matrix.
!
!    Input/output, real TEMP(N)...
!
!    Input/output, real X1(N), ...
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, actual number of iterations performed.
!
!    Input/output, integer INFO.
!
!    Output, real SCLR1, SCLR2, used to pass the scalars used in MATVEC.
!
!    Input/output, integer JOB, a commnication code.
!    On input:
!      1, this is the first call for a given problem.
!      2, this is a repeated call for a given problem.
!    On output:
!     -1, ITER > MAXIT, or stopping criterion was satisfied.
!      1, request for SCLR1 * A * X + SCLR2 * TEMP;
!      2, request for stopping criterion check.
!
  implicit none

  integer n

  real b(n)
  real diag(n)
  integer job
  integer info
  integer iter
  integer, save :: maxit = 0
  integer, save :: rlbl = 0
  real sclr1
  real sclr2
  real temp(n)
  real x(n)
  real x1(n)
!
!  JOB = 1, first call.
!
  if ( job == 1 ) then

    info = 0
    maxit = iter

    iter = 0
!
!  JOB = 2, repeated call.
!  Compute the error and check for acceptable convergence.
!
  else

    if ( rlbl == 2 ) then

      x(1:n) = diag(1:n) * temp(1:n)

      x1(1:n) = x1(1:n) - x(1:n)

      rlbl = 3
      job = 2
      return

    else if ( rlbl == 3 ) then

      if ( info == 1 ) then
        info = 0
        rlbl = -1
        job = -1
        return
      end if

      if ( iter >= maxit ) then
        info = 1
        rlbl = -1
        job = -1
        return
      end if

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'JACOBI_REVCOM - Fatal error!'
      write ( *, '(a,i6)' ) '  Illegal value of RLBL = ', rlbl
      stop

    end if

  end if
!
!  Perform the next Jacobi iteration
!
  iter = iter + 1
!
!  Save the current approximation to X in X1.
!
  x1(1:n) = x(1:n)
!
!  Apply iteration; result is updated approximation vector X.
!
  temp(1:n) = b(1:n)
  sclr1 = 1.0E+00
  sclr2 = 1.0E+00

  rlbl = 2
  job = 1

  return
end
subroutine jacobi_split_ge ( n, a, diag )

!*****************************************************************************80
!
!! JACOBI_SPLIT_GE splits a GE matrix for the Jacobi algorithm.
!
!  Discussion:
!
!    The GE format is the LINPACK/LAPACK general matrix storage mode,
!    in which a full N by N matrix is stored in an N by N array.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real A(N,N).
!    On input, the matrix to be split.
!    On output, the split matrix.
!
!    Output, real DIAG(N), the inverse of the diagonal entries of
!    the system matrix.
!
  implicit none

  integer n

  real a(n,n)
  integer i
  integer j
  real diag(n)

  do i = 1, n
    diag(i) = 1.0E+00 / a(i,i)
  end do

  a(1:n,1:n) = -a(1:n,1:n)

  do i = 1, n
    a(i,i) = 0.0E+00
  end do

  return
end
subroutine jacobi_recon_ge ( n, a, diag )

!*****************************************************************************80
!
!! JACOBI_RECON_GE reconstitutes a split GE matrix for the Jacobi algorithm.
!
!
!  Discussion:
!
!    The GE format is the LINPACK/LAPACK general matrix storage mode,
!    in which a full N by N matrix is stored in an N by N array.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real A(N,N).
!    On input, the matrix to be reconstituted.
!    On output, the reconstituted matrix.
!
!    Input, real DIAG(N), the inverse of the diagonal entries of
!    the system matrix.
!
  implicit none

  integer n

  real a(n,n)
  real diag(n)
  integer i

  a(1:n,1:n) = -a(1:n,1:n)

  do i = 1, n
    a(i,i) = 1.0E+00 / diag(i)
  end do

  return
end
subroutine matvec_ge ( n, a, alpha, x, beta, y, z )

!*****************************************************************************80
!
!! MATVEC_GE computes Z := ALPHA * A * X + BETA * Y for a GE matrix.
!
!  Discussion:
!
!    A GE matrix is one that is stored in LINPACK/LAPACK "general"
!    storage mode.  That is, the N by N matrix is stored in an N by
!    N array.
!
!  Modified:
!
!    20 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer n

  real a(n,n)
  real alpha
  real beta
  integer i
  integer j
  real x(n)
  real y(n)
  real z(n)

  do i = 1, n
    z(i) = beta * y(i)
    do j = 1, n
      z(i) = z(i) + alpha * a(i,j) * x(j)
    end do
  end do

  return
end
subroutine matvec_gb2 ( alpha, x, beta, y )

!*****************************************************************************80
!
!! MATVEC_GB computes Z := ALPHA * A * X + BETA * Y for a GB matrix.
!
!  Discussion:
!
!    A GB matrix is one that is stored in LINPACK/LAPACK "general band"
!    storage mode.  That is, the N by N matrix is assumed to have
!    a lower bandwidth of ML, an upper bandwidth of MU, and is stored
!    HOW EXACTLY?  DO WE INCLUDE EXTRA PIVOTING ENTRIES???
!
!  MVGB performs the matrix-vector product
!
!    y : = alpha*A*x+beta*y,
!
!  where alpha and beta are scalars, X and Y are vectors,
!  and A is a matrix.  Vector X must remain unchanged.
!  The solution is overwritten on vector Y.
!
!  The matrix A is passed into the routine in a common block.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
  implicit none

  integer, parameter :: maxdim = 200
  integer, parameter :: maxdim2 = maxdim * maxdim

  integer n

  real a(maxdim2)
  real alpha
  real beta
  integer ml
  integer mu
  real x(n)
  real y(n)

  common /matdim/ n
  common /matban/ ml,mu
  common /matrix/ a

  save /matban/
  save /matdim/
  save /matrix/

  if ( n > maxdim ) then
    write ( *, * ) ' '
    write ( *, * ) 'MVGB - Fatal error!'
    write ( *, * ) '  Matrix size N = ',n
    write ( *, * ) '  Exceeds maximum legal size MAXDIM = ',maxdim
    stop
  end if

  call sgbmv ( 'notranspose', n, n, ml, mu, alpha, a, n, x, 1, beta, &
    y, 1 )

  return
end
subroutine mvge ( alpha, x, beta, y )

!*****************************************************************************80
!
!! MVGE is a version of the MATVEC routine which assumes that the
!  matrix is in dense format, and uses the BLAS DGEMV.
!
!  MVGE performs the matrix-vector product
!
!    Y : = alpha * A * x + beta * y,
!
!  where alpha and beta are scalars, x and y are vectors,
!  and A is a matrix.  Vector x must remain unchanged.
!  The solution is over-written on vector y.
!
!  The matrix A is passed into the routine in a common block.
!
!    The GE format is the LINPACK/LAPACK general matrix storage mode,
!    in which a full N by N matrix is stored in an N by N array.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
  implicit none

  integer, parameter :: maxdim = 200
  integer, parameter :: maxdim2 = maxdim*maxdim

  integer n

  real a(maxdim2)
  real alpha
  real beta
  real x(n)
  real y(n)

  common /matdim/ n
  common /matrix/ a

  save /matdim/
  save /matrix/

  if ( n > maxdim ) then
    write ( *, * ) ' '
    write ( *, * ) 'MVGE - Fatal error!'
    write ( *, * ) '  Matrix size N = ',n
    write ( *, * ) '  Exceeds maximum legal size MAXDIM = ',maxdim
    stop
  end if

  call sgemv ( 'notranspose', n, n, alpha, a, n, x, 1, beta, y, 1 )

  return
end
subroutine orthoh ( i, n, h, v, ldv, w )

!*****************************************************************************80
!
!! ORTHOH constructs the I-th column of the upper Hessenberg matrix H
!  using the Gram-Schmidt process on V and W.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
  implicit none

  integer i
  integer ldv
  integer n

  real snrm2
  real h(i+1)
  integer j
  integer k
  real v(ldv,i+1)
  real w(n)
  real wnorm

  do k = 1, i

    h(k) = dot_product ( w(1:n), v(1:n,k) )

    call saxpy ( n, -h(k), v(1,k), 1, w, 1 )

  end do

  wnorm = snrm2 ( n, w, 1 )
  h(i+1) = wnorm

  v(1:n,i+1) = w(1:n) / wnorm

  return
end
subroutine psolve ( n, x, b, curpform )

!*****************************************************************************80
!
!! PSOLVE calls the appropriate preconditioner solver.
!
!  Modified:
!
!    22 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Output, real X(N), the solution.
!
!    Input, real B(N), the right hand side.
!
!    Input, character ( len = 8 ) CURPFORM, specifies the preconditioner.
!
  implicit none

  integer n

  real b(n)
  character ( len = 8 ) curpform
  logical lsame
  real x(n)

  if ( lsame ( curpform, 'identity' ) ) then

    call psolve_none ( n, x, b )

  else if ( lsame ( curpform, 'jacobi' ) ) then

    call psolve_jacobi ( n, x, b )

  else
    write ( *, * ) ' '
    write ( *, * ) 'PSOLVE - Fatal error!'
    write ( *, * ) '  Unknown preconditioner' // trim ( curpform )
    stop
  end if

  return
end
subroutine psolve_jacobi ( n, x, b )

!*****************************************************************************80
!
!! PSOLVE_JACOBI ??
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Output, real X(N), the solution.
!
!    Input, real B(N), the right hand side.
!
  implicit none

  integer, parameter :: maxdim = 200

  integer n

  real b(n)
  real m(maxdim)
  real x(n)

  common /matpre/ m

  save /matpre/

  x(1:n) = b(1:n) / m(1:n)

  return
end
subroutine psolve_jacobi_trans ( n, x, b )

!*****************************************************************************80
!
!! PSOLVE_JACOBI_TRANS solves the linear system Mx = b where M is diagonal.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Output, real X(N), the solution.
!
!    Input, real B(N), the right hand side.
!
  implicit none

  integer n

  real b(n)
  real x(n)

  call psolve_jacobi ( n, x, b )

  return
end
subroutine psolve_none ( n, x, b )

!*****************************************************************************80
!
!! PSOLVE_NONE is for the unpreconditioned version, i.e. just does
!  a vector copy (B to X ) then returns.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Output, real X(N), the solution.
!
!    Input, real B(N), the right hand side.
!
  implicit none

  integer n

  real b(n)
  real x(n)

  x(1:n) = b(1:n)

  return
end
subroutine psolve_none_trans ( n, x, b )

!*****************************************************************************80
!
!! PSOLVE_NONE_TRANS is for the unpreconditioned version, i.e. just does
!  a vector copy (B to X ) then returns.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Output, real X(N), the solution.
!
!    Input, real B(N), the right hand side.
!
  implicit none
!
  integer n

  real b(n)
  real x(n)

  x(1:n) = b(1:n)

  return
end
subroutine psolve_q ( n, x, b, which, curpform )

!*****************************************************************************80
!
!! PSOLVE_Q is a solver for QMR which requires left preconditioning
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Output, real X(N), the solution.
!
!    Input, real B(N), the right hand side.
!
  implicit none

  integer n

  real b(n)
  character ( len = 8 ) curpform
  logical lsame
  character ( len = 4 ) which
  real x(n)

  if ( lsame ( curpform, 'identity' ) ) then

     call psolve_none ( n, x, b )

  else if ( lsame ( curpform, 'jacobi' ) ) then

     if ( lsame(which, 'left' ) ) then
        call psolve_jacobi ( n, x, b )
     else
        call psolve_none ( n, x, b )
     end if

  else

    write ( *, * ) ' '
    write ( *, * ) 'PSOLVE_Q - Fatal error!'
    write ( *, * ) '  Unknown preconditioner', curpform
    stop

  end if

  return
end
subroutine psolve_t ( n, x, b, curpform )

!*****************************************************************************80
!
!! PSOLVE_T calls the appropriate solver.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Output, real X(N), the solution.
!
!    Input, real B(N), the right hand side.
!
  implicit none

  integer n

  real b(n)
  character ( len = 8 ) curpform
  logical lsame
  real x(n)

  if ( lsame ( curpform, 'identity') ) then
    call psolve_none_trans ( n, x, b )
  else if ( lsame ( curpform, 'jacobi') ) then
    call psolve_jacobi_trans ( n, x, b )
  else
    write ( *, * ) ' '
    write ( *, * ) 'PSOLVE_T - Fatal error!'
    write ( *, * ) '  Unknown preconditioner:',curpform
    stop
  end if

  return
end
subroutine psolve_t_q ( n, x, b, which, curpform )

!*****************************************************************************80
!
!! PSOLVE_T_Q is a solver for QMR which requires right preconditioning.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the system.
!
!    Output, real X(N), the solution.
!
!    Input, real B(N), the right hand side.
!
  implicit none

  integer n

  real b(n)
  character ( len = 8 ) curpform
  logical lsame
  character ( len = 4 ) which
  real x(n)

  if ( lsame ( curpform, 'identity' ) ) then
    call psolve_none ( n, x, b )
  else if ( lsame ( curpform, 'jacobi' ) ) then
    if ( lsame ( which, 'left') ) then
      call psolve_jacobi ( n, x, b )
    else
      call psolve_none ( n, x, b )
    end if
  else
    write ( *, * ) ' '
    write ( *, * ) 'PSOLVE_T_Q - Fatal error!'
    write ( *, * ) '  Unknown preconditioner' // curpform
    stop
  end if

  return
end
subroutine qmr ( n, b, x, work, iter, resid, matvec, matvect, &
  psolve_q, psolve_t_q, info, curpform )

!*****************************************************************************80
!
!! QMR implements the Quasi-Minimal Residual method.
!
!  Discussion:
!
!    Preconditioning is used.  The convergence test is:
!
!      (norm(b-A*x) / norm(b)) < TOL.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real B(N), the right hand side.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Workspace, real WORK(N,11).
!    Workspace for residual, direction vector, etc.
!    Note that W and WTLD, Y and YTLD, and Z and ZTLD share
!    workspace.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, actual number of iterations performed.
!
!    Input/output, real RESID.
!    On input, the allowable convergence measure for norm(b-A*x) / norm(b).
!    On output, the final value of this measure.
!
!  MATVEC  (external subroutine)
!          The user must provide a subroutine to perform the
!          matrix-vector product
!
!               y : = alpha*A*x+beta*y,
!
!          where alpha and beta are scalars, x and y are vectors,
!          and A is a matrix. Vector x must remain unchanged.
!          The solution is over-written on vector y.
!
!          The call is:
!
!             CALL MATVEC(ALPHA, X, BETA, Y)
!
!          The matrix is passed into the routine in a common block.
!
!  MATVECT  (external subroutine)
!          The user must provide a subroutine to perform the
!          matrix-vector product
!
!               y : = alpha*A'*x+beta*y,
!
!          where alpha and beta are scalars, x and y are vectors,
!          and A' is the tranpose of a matrix A. Vector x must remain
!          unchanged.
!          The solution is over-written on vector y.
!
!          The call is:
!
!             CALL MATVECT(ALPHA, X, BETA, Y)
!
!          The matrix is passed into the routine in a common block.
!
!  PSOLVEQ  (external subroutine)
!          The user must provide a subroutine to perform the
!          preconditioner solve routine for the linear system
!
!               M*x = b,
!
!          where x and b are vectors, and M a matrix. As QMR uses left
!          and right preconditioning and the preconditioners are in
!          common, we must specify in the call which to use. Vector b
!          must remain unchanged.
!          The solution is over-written on vector x.
!
!          The call is:
!
!             CALL PSOLVEQ(X, B, 'LEFT')
!
!          The preconditioner is passed into the routine in a common block.
!
!  PSOLVETQ  (external subroutine)
!          The user must provide a subroutine to perform the
!          preconditioner solve routine for the linear system
!
!               M'*x = b,
!
!          where x and y are vectors, and M' is the tranpose of a
!          matrix M. As QMR uses left and right preconditioning and
!          the preconditioners are in common, we must specify in the
!          call which to use. Vector b must remain unchanged.
!          The solution is over-written on vector x.
!
!          The call is:
!
!             CALL PSOLVETQ(X, B, 'LEFT')
!
!          The preconditioner is passed into the routine in a common block.
!
!  INFO    (output) integer
!
!        = 0: Successful exit. Iterated approximate solution returned.
!
!          >  0: Convergence to tolerance not achieved. This will be
!                set to the number of iterations performed.
!
!          <  0: Illegal input parameter, or breakdown occurred
!                during iteration.
!
!                Illegal parameter:
!
!                   -1: matrix dimension N < 0
!                   -3: Maximum number of iterations ITER < = 0.
!
!                BREAKDOWN: If parameters RHO or OMEGA become smaller
!                   than some tolerance, the program will terminate.
!                   Here we check against tolerance BREAKTOL.
!
!                  -10: RHO   < BREAKTOL: RHO and RTLD have become
!                                         orthogonal.
!                  -11: BETA  < BREAKTOL: EPS too small in relation to DELTA.
!                                         Convergence has stalled.
!                  -12: GAMMA < BREAKTOL: THETA too large.
!                                         Convergence has stalled.
!                  -13: DELTA < BREAKTOL: Y and Z have become
!                                         orthogonal.
!                  -14: EPS   < BREAKTOL: Q and PTLD have become
!                                         orthogonal.
!                  -15: XI    < BREAKTOL: Z too small. Convergence has stalled.
!
!                  BREAKTOL is set in function GETBREAK.
!
  implicit none

  integer n

  real b(n)
  real bnrm2
  character ( len = 8 ) curpform
  integer job
  integer info
  integer iter
  integer ndx1
  integer ndx2
  real resid
  real sclr1
  real sclr2
  real snrm2
  real tol
  real work(n*11)
  real x(n)

  external matvec
  external matvect
  external psolve_q
  external psolve_t_q

  info = 0
!
!  Test the input parameters.
!
  if ( n < 1 ) then
    info = -1
    write ( *, * ) ' '
    write ( *, * ) 'QMR - Fatal error!'
    write ( *, * ) '  N is less than 1.'
    write ( *, * ) '  N = ',n
    stop
  end if

  if ( iter <= 0 ) then
    info = -3
    return
  end if
!
!  Stop test may need some indexing info from REVCOM
!  use the init call to send the request across. REVCOM
!  will note these requests, and everytime it asks for
!  stop test to be done, it will provide the indexing info.
!
!  1 = =R; 2 == D; 3 == P; 4 == PTLD; 5 == Q; 6 == S; 7 == V;
!  8 = =VTLD; 9 == W; 10 == WTLD; 11 == Y; 12 == YTLD; 13 == Z;
!  14 = =ZTLD; -1 == ignore; any other == error
!
  ndx1 = 1
  ndx2 = -1
  tol = resid
  bnrm2 = snrm2 ( n, b, 1 )
!
!  First time call always init.
!
  job = 1

  do

    call qmr_revcom ( n, b, x, work, iter, resid, info, ndx1, ndx2, &
      sclr1, sclr2, job, curpform )
!
!  -1: Terminate.
!
    if ( job == -1 ) then

      exit
!
!  1: Compute WORK(NDX2) = SCLR1 * A * WORK(NDX1) + SCLR2 * WORK(NDX2).
!
    else if ( job == 1 ) then

      call matvec ( sclr1, work(ndx1), sclr2, work(ndx2) )
!
!  2: Compute WORK(NDX2) = SCLR1 * A' * WORK(NDX1) + SCLR2 * WORK(NDX2).
!
    else if ( job == 2 ) then

      call matvect ( sclr1, work(ndx1), sclr2, work(ndx2) )
!
!  Call left preconditioned solve.
!
    else if ( job == 3 ) then

      call psolve_q ( n, work(ndx1), work(ndx2), 'left', curpform )
!
!  Call right preconditioned solve.
!
    else if ( job == 4 ) then

      call psolve_q ( n, work(ndx1), work(ndx2), 'right', curpform )
!
!  Call left preconditioned transpose solve.
!
    else if ( job == 5 ) then

      call psolve_t_q ( n, work(ndx1), work(ndx2), 'left', curpform )
!
!  Call right preconditioned transpose solve.
!
    else if ( job == 6 ) then

      call psolve_t_q ( n, work(ndx1), work(ndx2), 'right', curpform )
!
!  7: Compute WORK(NDX2) = SCLR1 * A * X + SCLR2 * WORK(NDX2).
!
    else if ( job == 7 ) then

      call matvec ( sclr1, x, sclr2, work(ndx2) )
!
!  8: Do a stopping test on the relative residual reduction.
!
    else if ( job == 8 ) then

      call stopb ( n, work(ndx1), bnrm2, resid, tol, info )

    end if

    job = 2

  end do

  return
end
subroutine qmr_revcom ( n, b, x, work, iter, resid, info, ndx1, ndx2, &
  sclr1, sclr2, job, curpform )

!*****************************************************************************80
!
!! QMR_REVCOM is controlled by QMR using reverse communication.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real B(N), the right hand side.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Workspace, real WORK(N,11).
!    Workspace for residual, direction vector, etc.
!    Note that W and WTLD, Y and YTLD, and Z and ZTLD share workspace.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, actual number of iterations performed.
!
!    Input/output, real RESID.
!    On input, the allowable convergence measure for norm(b-A*x) / norm(b).
!    On output, the final value of this measure.
!
!  INFO    (output) integer
!
!        = 0: Successful exit. Iterated approximate solution returned.
!            -5: Erroneous NDX1/NDX2 in INIT call.
!            -6: Erroneous RLBL.
!
!          >  0: Convergence to tolerance not achieved. This will be
!                set to the number of iterations performed.
!
!          <  0: Illegal input parameter, or breakdown occurred
!                during iteration.
!
!                Illegal parameter:
!
!                   -1: matrix dimension N < 0
!                   -3: Maximum number of iterations ITER < = 0.
!
!                BREAKDOWN: If parameters RHO or OMEGA become smaller
!                   than some tolerance, the program will terminate.
!                   Here we check against tolerance BREAKTOL.
!
!                  -10: RHO   < BREAKTOL: RHO and RTLD have become
!                                         orthogonal.
!                  -11: BETA  < BREAKTOL: EPS too small in relation to DELTA.
!                                         Convergence has stalled.
!                  -12: GAMMA < BREAKTOL: THETA too large.
!                                         Convergence has stalled.
!                  -13: DELTA < BREAKTOL: Y and Z have become
!                                         orthogonal.
!                  -14: EPS   < BREAKTOL: Q and PTLD have become
!                                         orthogonal.
!                  -15: XI    < BREAKTOL: Z too small. Convergence has stalled.
!
!                  BREAKTOL is set in function GETBREAK.
!
!  NDX1    (input/output) integer.
!  NDX2    On entry in INIT call contain indices required by interface
!          level for stopping test.
!          All other times, used as output, to indicate indices into
!          WORK[] for the MATVEC, PSOLVE done by the interface level.
!
!  SCLR1   (output) real.
!  SCLR2   Used to pass the scalars used in MATVEC.
!
!  JOB    (input/output) integer.
!          Used to communicate job code between the two levels.
!
  implicit none

  real, parameter :: one = 1.0E+00

  integer n

  real b(n)
  real beta
  real betatol
  real bnrm2
  real c1
  character ( len = 8 ) curpform
  integer d
  real sdot
  real delta
  real deltatol
  real snrm2
  real, save :: eps = 0.0E+00
  real epstol
  real eta
  real gamma
  real gamma1
  real gammatol
  real getbreak
  integer i
  integer info
  integer iter
  integer job
  integer maxit
  integer ndx1
  integer ndx2
  integer need1
  integer need2
  integer p
  integer ptld
  integer q
  integer r
  real resid
  real rho
  real rho1
  real rhotol
  integer, save :: rlbl = 0
  integer s
  real sclr1
  real sclr2
  real theta
  real theta1
  real tol
  integer v
  integer vtld
  integer w
  real work(n,*)
  integer wtld
  real x(n)
  real xi
  real xitol
  integer y
  integer ytld
  integer z
  integer ztld

  save

  if ( job == 1 ) then
     go to 1
  else if ( job == 2 ) then
     if ( rlbl == 2) go to 2
     if ( rlbl == 3) go to 3
     if ( rlbl == 4) go to 4
     if ( rlbl == 5) go to 5
     if ( rlbl == 6) go to 6
     if ( rlbl == 7) go to 7
     if ( rlbl == 8) go to 8
     if ( rlbl == 9) go to 9
     if ( rlbl == 10) go to 10
     info = -6
     go to 20
  end if

 1    continue

  info = 0
  maxit = iter
  tol = resid
!
!  Alias workspace columns.
!
  r   = 1
  d   = 2
  p   = 3
  ptld = 4
  q   = 5
  s   = 6
  v   = 7
  vtld = 8
  w   = 9
  wtld = 9
  y   = 10
  ytld = 10
  z   = 11
  ztld = 11
!
!  Check if caller will need indexing info.
!
  if ( ndx1 == -1 ) then
    need1 = ndx1
  else if ( ndx1 == 1 ) then
        need1 = ((r-1)*n)+1
     else if ( ndx1 == 2 ) then
        need1 = ((d-1)*n)+1
     else if ( ndx1 == 3 ) then
        need1 = ((p-1)*n)+1
     else if ( ndx1 == 4 ) then
        need1 = ((ptld-1)*n)+1
     else if ( ndx1 == 5 ) then
        need1 = ((q-1)*n)+1
     else if ( ndx1 == 6 ) then
        need1 = ((s-1)*n)+1
     else if ( ndx1 == 7 ) then
        need1 = ((v-1)*n)+1
     else if ( ndx1 == 8 ) then
        need1 = ((vtld-1)*n)+1
     else if ( ndx1 == 9 ) then
        need1 = ((w-1)*n)+1
     else if ( ndx1 == 10 ) then
        need1 = ((wtld-1)*n)+1
     else if ( ndx1 == 11 ) then
        need1 = ((y-1)*n)+1
     else if ( ndx1 == 12 ) then
        need1 = ((ytld-1)*n)+1
     else if ( ndx1 == 13 ) then
        need1 = ((z-1)*n)+1
     else if ( ndx1 == 14 ) then
        need1 = ((ztld-1)*n)+1
     else
        info = -5
        go to 20
     end if

  if ( ndx2==-1 ) then
    need2 = ndx2
  else if ( ndx2 == 1 ) then
        need2 = ((r-1)*n)+1
     else if ( ndx2 == 2 ) then
        need2 = ((d-1)*n)+1
     else if ( ndx2 == 3 ) then
        need2 = ((p-1)*n)+1
     else if ( ndx2 == 4 ) then
        need2 = ((ptld-1)*n)+1
     else if ( ndx2 == 5 ) then
        need2 = ((q-1)*n)+1
     else if ( ndx2 == 6 ) then
        need2 = ((s-1)*n)+1
     else if ( ndx2 == 7 ) then
        need2 = ((v-1)*n)+1
     else if ( ndx2 == 8 ) then
        need2 = ((vtld-1)*n)+1
     else if ( ndx2 == 9 ) then
        need2 = ((w-1)*n)+1
     else if ( ndx2 == 10 ) then
        need2 = ((wtld-1)*n)+1
     else if ( ndx2 == 11 ) then
        need2 = ((y-1)*n)+1
     else if ( ndx2 == 12 ) then
        need2 = ((ytld-1)*n)+1
     else if ( ndx2 == 13 ) then
        need2 = ((z-1)*n)+1
     else if ( ndx2 == 14 ) then
        need2 = ((ztld-1)*n)+1
     else
        info = -5
        go to 20
     end if
!
!  Set breakdown tolerances.
!
  rhotol = getbreak()
  betatol = getbreak()
  gammatol = getbreak()
  deltatol = getbreak()
  epstol = getbreak()
  xitol  = getbreak()
!
!  Set initial residual.
!
  do i = 1, n
    work(i,r) = b(i)
  end do

  if ( snrm2(n,x,1) /= 0.0E+00 ) then
    sclr1 = -1.0E+00
    sclr2 = 0.0E+00
    ndx1 = ((d-1)*n)+1
    ndx2 = ((r-1)*n)+1
    rlbl = 2
    job = 7
    return
  end if

 2    continue

  if ( snrm2(n,work(1,r),1) < tol ) then
    go to 30
  end if

  bnrm2 = snrm2(n,b,1)

  if ( bnrm2 == 0.0E+00 ) then
    bnrm2 = 1.0E+00
  end if

  work(1:n,vtld) = work(1:n,r)
  ndx1 = ((y   -1)*n)+1
  ndx2 = ((vtld-1)*n)+1
  rlbl = 3
  job = 3
  return

 3    continue

  rho = snrm2(n,work(1,y),1)
  work(1:n,wtld) = work(1:n,r)
  ndx1 = ((z   -1)*n)+1
  ndx2 = ((wtld-1)*n)+1
  rlbl = 4
  job = 6
  return

 4    continue

  xi = snrm2(n,work(1,z),1)

  gamma = 1.0E+00
  eta = -one
  theta = 0.0E+00

  iter = 0

   40 continue
!
!  Perform Preconditioned QMR iteration.
!
  iter = iter + 1

  if ( abs ( rho ) < rhotol .or. abs ( xi ) < xitol ) then
    go to 25
  end if

  do i = 1, n
    work(i,v) = work(i,vtld) / rho
  end do

  do i = 1, n
    work(i,y) = work(i,y) / rho
  end do

  do i = 1, n
    work(i,w) = work(i,wtld) / xi
  end do

  do i = 1, n
    work(i,z) = work(i,z) / xi
  end do

  delta = sdot ( n, work(1,z), 1, work(1,y), 1 )

  if ( abs ( delta ) < deltatol ) then
    go to 25
  end if

  ndx1 = ((ytld-1)*n)+1
  ndx2 = ((y   -1)*n)+1
  rlbl = 5
  job = 4
  return

 5    continue

  ndx1 = ((ztld-1)*n)+1
  ndx2 = ((z   -1)*n)+1
  rlbl = 6
  job = 5
  return

 6    continue

  if ( iter > 1 ) then

    c1 = - xi * delta / eps
    call saxpy ( n, c1, work(1,p), 1, work(1,ytld), 1 )

    do i = 1, n
      work(i,p) = work(i,ytld)
    end do

    call saxpy ( n, -(rho*delta / eps), work(1,q), 1, work(1,ztld), 1 )

  else

    do i = 1, n
      work(i,p) = work(i,ytld)
    end do

  end if

     do i = 1, n
       work(i,q) = work(i,ztld)
     end do

     sclr1 = 1.0E+00
     sclr2 = 0.0E+00
     ndx1 = ((p   -1)*n)+1
     ndx2 = ((ptld-1)*n)+1
     rlbl = 7
     job = 1
     return

 7       continue

     eps = sdot(n,work(1,q),1,work(1,ptld),1)
     if ( abs(eps) < epstol ) then
       go to 25
     end if

     beta = eps / delta
     if ( abs(beta) < betatol) then
       go to 25
     end if

     do i = 1, n
       work(i,vtld) = work(i,ptld)
     end do

     call saxpy(n,-beta, work(1,v),1,work(1,vtld),1)
     call psolve_q ( n, work(1,y), work(1,vtld), 'left', curpform )

     rho1 = rho
     rho = snrm2(n,work(1,y),1)

     do i = 1, n
       work(i,wtld) = work(i,w)
     end do

     sclr1 = 1.0E+00
     sclr2 = -beta
     ndx1 = ((q   -1)*n)+1
     ndx2 = ((wtld-1)*n)+1
     rlbl = 8
     job = 2
     return

 8       continue

     ndx1 = ((z   -1)*n)+1
     ndx2 = ((wtld-1)*n)+1
     rlbl = 9
     job = 6
     return

 9       continue

     xi = snrm2(n,work(1,z),1)

     gamma1 = gamma
     theta1 = theta

     theta = rho / (gamma1*abs(beta))
     gamma = 1.0E+00 / sqrt(1.0+theta**2)
     if ( abs(gamma) < gammatol ) then
       go to 25
     end if

     eta = -eta*rho1 * gamma**2 / (beta * gamma1**2)

     if ( iter > 1 ) then
        call sscal(n,(theta1*gamma)**2, work(1,d),1)
        call saxpy(n,eta, work(1,p),1,work(1,d),1)
        call sscal(n,(theta1*gamma)**2, work(1,s),1)
        call saxpy(n,eta, work(1,ptld),1,work(1,s),1)
     else

        do i = 1, n
          work(i,d) = eta * work(i,p)
        end do

        do i = 1, n
          work(i,s) = eta * work(i,ptld)
        end do

     end if
!
!  Compute current solution vector x.
!
     call saxpy(n,one, work(1,d),1,x,1)
!
!  Compute residual vector rk, find norm,
!  then check for tolerance.
!
     call saxpy(n,-one, work(1,s),1,work(1,r),1)

     ndx1 = need1
     ndx2 = need2
     rlbl = 10
     job = 8
     return

 10   continue

  if ( info == 1 ) then
    go to 30
  end if

  if ( iter == maxit ) then
    info = 1
    go to 20
  end if

  go to 40

   20 continue
!
!  Iteration fails.
!
  rlbl = -1
  job = -1

  return

   25 continue
!
!  Method breakdown.
!
  if ( abs(rho) < rhotol ) then
    info = -10
  else if ( abs(beta) < betatol ) then
    info = -11
  else if ( abs(gamma) < gammatol ) then
    info = -12
  else if ( abs(delta) < deltatol ) then
    info = -13
  else if ( abs(eps) < epstol ) then
    info = -14
  else if ( abs(xi) < xitol ) then
    info = -15
  end if

  rlbl = -1
  job = -1

  return

   30 continue
!
!  Iteration successful; return.
!
  info = 0
  rlbl = -1
  job = -1

  return
end
subroutine resid_ge ( a, b, n, rmax, x )

!*****************************************************************************80
!
!! RESID_GE computes the residual A*X-B when A is stored in GE format.
!
!  Discussion:
!
!    The GE format is the LINPACK/LAPACK general matrix storage mode,
!    in which a full N by N matrix is stored in an N by N array.
!
!  Modified:
!
!    28 April 2000
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, real A(N,N), the unfactored N by N matrix.
!
!    Input, real B(N), the right hand side.
!
!    Input, integer N, the order of the matrix.
!
!    Output, real RMAX, the maximum absolute residual.
!
!    Input, real X(N), an estimate for the solution.
!
  implicit none

  integer n

  real a(n,n)
  real b(n)
  integer i
  real rmax
  real x(n)

  rmax = 0.0E+00

  do i = 1, n

    rmax = max ( rmax, abs ( b(i) - dot_product ( a(i,1:n), x(1:n) ) ) )

  end do

  return
end
subroutine rotvec ( x, y, c, s )

!*****************************************************************************80
!
!! ROTVEC applies a Givens rotation to a vector (X,Y).
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input/output, real X, Y, the two entries of the vector.
!
!    Input, real C, S, the cosine and sine of an angle, which define
!    the Givens rotation.
!
  implicit none

  real c
  real s
  real x
  real x2
  real y
  real y2

  x2 = c * x - s * y
  y2 = s * x + c * y

  x = x2
  y = y2

  return
end
subroutine r4vec_print_some ( n, a, max_print, title )

!*****************************************************************************80
!
!! R4VEC_PRINT_SOME prints "some" of an R4VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Modified:
!
!    13 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, real A(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer n

  real a(n)
  integer i
  integer max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(a)' ) ' '
  end if

  if ( n <= max_print ) then

    if ( all ( a(1:n) == aint ( a(1:n) ) ) ) then
      do i = 1, n
        write ( *, '(i6,2x,i6)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:n) ) < 1000000.0E+00 ) ) then
      do i = 1, n
        write ( *, '(i6,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, n
        write ( *, '(i6,2x,g14.6)' ) i, a(i)
      end do
    end if

  else if ( max_print >= 3 ) then

    if ( all ( a(1:max_print-2) == aint ( a(1:max_print-2) ) ) ) then
      do i = 1, max_print-2
        write ( *, '(i6,2x,i6)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:max_print-2) ) < 1000000.0E+00 ) ) then
      do i = 1, max_print-2
        write ( *, '(i6,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, max_print-2
        write ( *, '(i6,2x,g14.6)' ) i, a(i)
      end do
    end if

    write ( *, '(a)' ) '......  ..............'
    i = n

    if ( a(n) == aint ( a(n) ) ) then
      write ( *, '(i6,2x,i6)' ) i, int ( a(i) )
    else if (  abs ( a(n) ) < 1000000.0E+00 ) then
      write ( *, '(i6,2x,f14.6)' ) i, a(i)
    else
      write ( *, '(i6,2x,g14.6)' ) i, a(i)
    end if

  else

    if ( all ( a(1:max_print-1) == aint ( a(1:max_print-1) ) ) ) then
      do i = 1, max_print-1
        write ( *, '(i6,2x,i6)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:max_print-1) ) < 1000000.0E+00 ) ) then
      do i = 1, max_print-1
        write ( *, '(i6,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, max_print-1
        write ( *, '(i6,2x,g14.6)' ) i, a(i)
      end do
    end if

    i = max_print

    if ( a(n) == aint ( a(n) ) ) then
      write ( *, '(i6,2x,i6,a)' ) i, int ( a(i) ), '...more entries...'
    else if (  abs ( a(n) ) < 1000000.0E+00 ) then
      write ( *, '(i6,2x,f14.6,a)' ) i, a(i), '...more entries...'
    else
      write ( *, '(i6,2x,g14.6,a)' ) i, a(i), '...more entries...'
    end if

  end if

  return
end
function samax ( n, x, incx )

!*****************************************************************************80
!
!! SAMAX returns the maximum absolute value of the entries in a vector.
!
!  Modified:
!
!    08 April 1999
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real X(*), the vector to be examined.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Output, real SAMAX, the maximum absolute value of an element of X.
!
  implicit none

  integer i
  integer incx
  integer ix
  integer n
  real samax
  real x(*)

  if ( n <= 0 ) then

    samax = 0.0E+00

  else if ( n == 1 ) then

    samax = abs ( x(1) )

  else if ( incx == 1 ) then

    samax = abs ( x(1) )

    do i = 2, n
      if ( abs ( x(i) ) > samax ) then
        samax = abs ( x(i) )
      end if
    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    samax = abs ( x(ix) )
    ix = ix + incx

    do i = 2, n
      if ( abs ( x(ix) ) > samax ) then
        samax = abs ( x(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
subroutine sgb_to_sge ( m, ml, mu, n, a1, a2 )

!*****************************************************************************80
!
!! SGB_TO_SGE converts a general band matrix to general matrix format.
!
!  Discussion:
!
!    LINPACK and LAPACK band storage requires that an extra ML
!    superdiagonals be supplied to allow for fillin during Gauss
!    elimination.  Even though a band matrix is described as
!    having an upper bandwidth of MU, it effectively has an
!    upper bandwidth of MU+ML.  This routine will copy nonzero
!    values it finds in these extra bands, so that both unfactored
!    and factored matrices can be handled.
!
!  Modified:
!
!    18 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrices.
!    M must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths of A1.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1.
!
!    Input, integer N, the number of columns of the matrices.
!    N must be positive.
!
!    Input, real A1(2*ML+MU+1,N), the M by N general band matrix.
!
!    Output, real A2(M,N), the M by N general matrix, which
!    contains the information given in A1.
!
  implicit none

  integer ml
  integer mu
  integer n

  real a1(2*ml+mu+1,n)
  real a2(m,n)
  integer i
  integer ierror
  integer j
  integer m
!
!  Check the dimensions.
!
  call sgb_check ( m, n, ml, mu, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGB_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for A1.'
    return
  end if

  call sge_check ( m, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGB_TO_SGE - Fatal error!'
    write ( *, * ) '  Illegal dimensions for A2.'
    return
  end if

  do i = 1, m
    do j = 1, n
      if ( i - ml <= j .and. j <= i + mu + ml ) then
        a2(i,j) = a1(ml+mu+1+i-j,j)
      else
        a2(i,j) = 0.0E+00
      end if
    end do
  end do

  return
end
subroutine sge_check ( m, n, ierror )

!*****************************************************************************80
!
!! SGE_CHECK checks the dimensions of a general matrix.
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Output, integer IERROR, reports whether any errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 2 if M is illegal;
!    IERROR = IERROR + 4 if N is illegal.
!
  implicit none

  integer ierror
  integer m
  integer n

  ierror = 0

  if ( m < 1 ) then
    ierror = ierror + 2
    write ( *, * ) ' '
    write ( *, * ) 'SGE_CHECK - Illegal M = ', m
  end if

  if ( n < 1 ) then
    ierror = ierror + 4
    write ( *, * ) ' '
    write ( *, * ) 'SGE_CHECK - Illegal N = ', n
  end if

  return
end
subroutine slt_sl ( a, n, b )

!*****************************************************************************80
!
!! SLT_SL solves a lower triangular system.
!
!  Discussion:
!
!    No factorization of the lower triangular matrix is required.
!
!  Modified:
!
!    22 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(N,N), the lower triangular matrix.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real B(N).
!
!    On input, the right hand side.
!    On output, the solution vector.
!
  implicit none

  integer n

  real a(n,n)
  real b(n)
  integer i
  integer j

  do j = 1, n
    b(j) = b(j) / a(j,j)
    do i = j + 1, n
      b(i) = b(i) - a(i,j) * b(j)
    end do
  end do

  return
end
function snrm2 ( n, x, incx )

!*****************************************************************************80
!
!! SNRM2 computes the Euclidean norm of a vector.
!
!  Discussion:
!
!    The original SNRM2 algorithm is accurate but written in a bizarre,
!    unreadable and obsolete format.  This version goes for clarity.
!
!  Modified:
!
!    01 June 2000
!
!  Reference:
!
!    Lawson, Hanson, Kincaid, Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real X(*), the vector whose norm is to be computed.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Output, real SNRM2, the Euclidean norm of X.
!
  implicit none

  integer i
  integer incx
  integer ix
  integer n
  real samax
  real snrm2
  real stemp
  real x(*)
  real xmax

  if ( n <= 0 ) then

    snrm2 = 0.0E+00

  else

    xmax = samax ( n, x, incx )

    if ( xmax == 0.0E+00 ) then

      snrm2 = 0.0E+00

    else

      if ( incx >= 0 ) then
        ix = 1
      else
        ix = ( - n + 1 ) * incx + 1
      end if

      stemp = 0.0E+00
      do i = 1, n
        stemp = stemp + ( x(ix) / xmax )**2
        ix = ix + incx
      end do

      snrm2 = xmax * sqrt ( stemp )

    end if

  end if

  return
end
subroutine sor_ge ( n, a, b, x, omega, iter, restol, resid, info )

!*****************************************************************************80
!
!! SOR_GE implements the Successive Over-Relaxation method for a GE matrix.
!
!  Discussion:
!
!    The GE format is the LINPACK/LAPACK general matrix storage mode,
!    in which a full N by N matrix is stored in an N by N array.
!
!    The matrix splitting:
!
!      N = the strict upper triangular portion of A, stored in WORK.
!      M = the lower triangular portion of A.
!
!    Relative error measured:
!
!      norm ( X - X_1 ) / norm ( X ).
!
!  Modified:
!
!    22 November 2000
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(N,N), the matrix.
!
!    Input, real B(N), the right hand side vector.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Input, real OMEGA, the relaxation parameter.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, actual number of iterations performed.
!
!    Input, real RESTOL, the allowable convergence measure for
!    norm(x-x_1) / norm(x).
!
!    Output, real RESID, the final value of the convergence measure.
!
!    Output, integer INFO, error flag.
!     0: Successful exit. Iterated approximate solution returned.
!    >0: Convergence to tolerance not achieved. INFO will be
!        set to the number of iterations performed.
!
  implicit none

  integer n

  real a(n,n)
  real a_copy(n,n)
  real b(n)
  integer job
  integer info
  integer iter
  integer ndx1
  integer ndx2
  real omega
  real resid
  real restol
  real sclr1
  real sclr2
  real work(n*(n+3))
  real x(n)
  real xnrm2

  info = 0
!
!  Test the input parameters.
!
  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SOR_GE - Fatal error!'
    write ( *, * ) '  N is less than 1.'
    write ( *, * ) '  N = ', n
    stop
  end if

  if ( iter <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SOR_GE - Fatal error!'
    write ( *, * ) '  ITER is less than 1.'
    stop
  end if

  if ( omega <= 0.0E+00 .or. omega >= 2.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SOR_GE - Fatal error!'
    write ( *, * ) '  0 < OMEGA < 2 is required.'
    write ( *, * ) '  Input OMEGA = ', omega
    stop
  end if
!
!  Stop test may need some indexing info from REVCOM
!  use the init call to send the request across. REVCOM
!  will note these requests, and everytime it asks for
!  stop test to be done, it will provide the indexing info.
!
!  -1 == ignore;
!  1 == X1;
!  2 == TEMP;
!  3 == MM;
!
  ndx1 = 1
  ndx2 = -1
!
!  First time call always init.
!
  job = 1

  do

    call sor_revcom ( n, a, b, x, omega, work, iter, restol, info, &
      ndx1, ndx2, sclr1, sclr2, job, a_copy )
!
!  -1: Terminate.
!
    if ( job == -1 ) then

      exit
!
!  1: Compute WORK(NDX2) = SCLR1 * A * X + SCLR2 * WORK(NDX2).
!
    else if ( job == 1 ) then

      call matvec_ge ( n, a, sclr1, x, sclr2, work(ndx2), work(ndx2) )
!
!  2: Solve L*x = y.
!
    else if ( job == 2 ) then

      call slt_sl ( work(ndx1), n, x )
!
!  3: Do a stopping test on the current residual relative to the norm of X.
!
    else if ( job == 3 ) then

      call stopx ( n, work(ndx1), x, xnrm2, resid, restol, info )

    end if

    job = 2

  end do

  return
end
subroutine sor_revcom ( n, a, b, x, omega, work, iter, restol, &
  info, ndx1, ndx2, sclr1, sclr2, job, a_copy )

!*****************************************************************************80
!
!! SOR_REVCOM is controlled by SOR using reverse communication.
!
!  Modified:
!
!    22 November 2000
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real B(N), the right hand side vector.
!
!    Input/output, real X(N), on input, an estimate for the solution.
!    On output, an improved estimate.
!
!    Input, real OMEGA, the relaxation parameter.
!
!    Input/workspace, real WORK(N*(N+3)).
!    The relaxation parameter, OMEGA, should be input in WORK(1).
!    The amount of workspace can be significantly reduced (to 2*N)
!    by customizing the matrix-vector product and BSOLVE.
!
!    Input/output, integer ITER.
!    On input, the maximum iterations to be performed.
!    On output, actual number of iterations performed.
!
!    Input, real RESTOL, the allowable convergence measure for
!    norm(x-x_1) / norm(x).
!
!    Output, real RESID, the final value of the convergence measure.
!
!  INFO    (output) integer
!
!        = 0: Successful exit. Iterated approximate solution returned.
!            -5: Erroneous NDX1/NDX2 in INIT call.
!            -6: Erroneous RLBL.
!
!          >  0: Convergence to tolerance not achieved. This will be
!                set to the number of iterations performed.
!
!          <  0: Illegal input parameter, or breakdown occurred
!                during iteration.
!
!                Illegal parameter:
!
!                   -1: matrix dimension N < 0
!                   -3: Maximum number of iterations ITER < = 0.
!                   -4: Relaxation parameter OMEGA not in interval (0,2).
!
!  NDX1    (input/output) integer.
!  NDX2    On entry in INIT call contain indices required by interface
!          level for stopping test.
!          All other times, used as output, to indicate indices into
!          WORK[] for the MATVEC, PSOLVE done by the interface level.
!
!  SCLR1   (output) real.
!  SCLR2   Used to pass the scalars used in MATVEC.
!
!  JOB    (input/output) integer.
!          Used to communicate job code between the two levels.
!
  implicit none

  integer n

  real a(n,n)
  real a_copy(n,n)
  real b(n)
  real bnrm2
  integer i
  integer job
  integer info
  integer iter
  integer maxit
  integer mm
  integer ndx1
  integer ndx2
  integer need1
  integer need2
  real omega
  real restol
  integer, save :: rlbl = 0
  real sclr1
  real sclr2
  real snrm2
  integer temp
  real work(n,n+3)
  real x(n)
  integer x1

  save

  if ( job == 1 ) then
    go to 1
  else if ( job == 2 ) then
     if ( rlbl == 2) go to 2
     if ( rlbl == 3) go to 3
     if ( rlbl == 4) go to 4
     if ( rlbl == 5) go to 5
     info = -6
     go to 20
  end if

 1    continue

  info = 0
  maxit = iter
!
!  Alias workspace columns.
!
  x1 = 1
  temp = 2
  mm = 3
!
!  Check if caller will need indexing info.
!
  if ( ndx1 == -1 ) then
    need1 = ndx1
  else if ( ndx1 == 1 ) then
    need1 = ((x1-1)*n)+1
  else if ( ndx1 == 2 ) then
    need1 = ((temp-1)*n)+1
  else if ( ndx1 == 3 ) then
    need1 = ((mm-1)*n)+1
  else
    info = -5
    go to 20
  end if

  if ( ndx2 == -1 ) then
    need2 = ndx2
  else if ( ndx2 == 1 ) then
    need2 = ((x1-1)*n)+1
  else if ( ndx2 == 2 ) then
    need2 = ((temp-1)*n)+1
  else if ( ndx2 == 3 ) then
    need2 = ((mm-1)*n)+1
  else
    info = -5
    go to 20
  end if
!
!  Compute initial residual for convergence criteria.
!
  do i = 1, n
    work(i,x1) = b(i)
  end do

  if ( snrm2(n,x,1) /= 0.0E+00 ) then
    ndx1 = -1
    ndx2 = ((x1  -1)*n)+1
    sclr1 = -1.0E+00
    sclr2 = 1.0E+00
    rlbl = 2
    job = 1
    return
  end if

 2    continue
  if ( snrm2(n,work(1,x1),1) < restol ) then
    info = 0
    go to 20
  end if

  bnrm2 = snrm2(n,b,1)
  if ( bnrm2 == 0.0E+00 ) then
    bnrm2 = 1.0E+00
  end if
!
!  Matrix A is set to N. WORK(1:N,1:N) is set to MM.
!
  call sor_split_ge ( omega, n, a, b, work(1,mm), n, a_copy )

  iter = 0

   10 continue
!
!  Perform an SOR iteration.
!
  iter = iter + 1
!
!  Save the current approximation to X in X1,
!

  do i = 1, n
    work(i,x1) = x(i)
  end do
!
!  Apply iteration; result is updated approximation vector X.
!
  do i = 1, n
    work(i,temp) = b(i)
  end do

  ndx1 = -1
  ndx2 = ((temp-1)*n)+1
  sclr1 = 1.0E+00
  sclr2 = 1.0E+00
  rlbl = 3
  job = 1

  return

 3    continue

  x(1:n) = work(1:n,temp)
  ndx1 = ((mm  -1)*n)+1
  ndx2 = -1
!
!  Prepare for return.
!
  rlbl = 4
  job = 2
  return

 4    continue
!
!  Compute error and check for acceptable convergence.
!
  work(1:n,x1) = work(1:n,x1) - x(1:n)

  ndx1 = need1
  ndx2 = need2
  rlbl = 5
  job = 3
  return

 5    continue
  if ( info == 1 ) then
    info = 0
    go to 20
  end if

  if ( iter == maxit ) then
    info = 1
    go to 20
  end if

  go to 10

   20 continue

  call sor_recon_ge ( omega, n, a, b, a_copy )

  rlbl = -1
  job = -1

  return
end
subroutine sor_split_ge ( omega, n, a, b, work, a_copy )

!*****************************************************************************80
!
!! SOR_SPLIT_GE splits a GE matrix for the SOR algorithm.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, real OMEGA, the relaxation parameter.
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(N,N), ...
!
!    Input/output, real B(N), on input, the right hand side.
!    On output, the modified right hand side.
!
!    ?, real WORK(N,*), ?
!
!    Output, real A_COPY(N,N), a copy of the A matrix.
!
  implicit none

  integer n

  real a(n,n)
  real a_copy(n,n)
  real b(n)
  integer i
  integer j
  real omega
  real work(n,*)
!
!  Save a copy of A.
!
  a_copy(1:n,1:n) = a(1:n,1:n)
!
!  SPLIT
!  Set the splitting matrix M.
!
!  Set NN and B.
!
!  Temporarily store the matrix A in order to reconstruct
!  the original matrix. Because the lower triangular portion
!  of A must be zeroed, this is the easiest way to deal with it.
!  This causes the requirement that WORK be N x (2N+3).
!
    do i = 1, n
      work(i,i) = a(i,i)
      do j = 1, i-1
        work(i,j) = omega * a(i,j)
      end do
    end do

    do i = 1, n
      a(i,i) = ( 1.0E+00 - omega ) * a(i,i)
      do j = i+1, n
        a(i,j) = - omega * a(i,j)
      end do
    end do

    do i = 2, n
      do j = 1, i-1
        a(i,j) = 0.0E+00
      end do
    end do

    b(1:n) = omega * b(1:n)

  return
end
subroutine sor_recon_ge ( omega, n, a, b, a_copy )

!*****************************************************************************80
!
!! SOR_RECON_GE reconstructs the matrix A and right hand side B after splitting.
!
!  Modified:
!
!    26 November 2000
!
!  Parameters:
!
!    Input, real OMEGA, the relaxation parameter.
!
!    Input, integer N, the order of the matrix.
!
!    Output, real A(N,N), ...
!
!    Input/output, real B(N), on input, the modified right hand side.
!    On output, the original right hand side.
!
!    Input, real A_COPY(N,N), a copy of the original A matrix.
!
  implicit none

  integer n

  real a(n,n)
  real a_copy(n,n)
  real b(n)
  real omega

  a(1:n,1:n) = a_copy(1:n,1:n)

  b(1:n) = b(1:n) / omega

  return
end
subroutine stopb ( n, r, bnrm2, resid, restol, info )

!*****************************************************************************80
!
!! STOPB computes the stopping criterion on B.
!
!  Modified:
!
!    21 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real R(N), the residual.
!
!    Input, real BNRM2, the euclidean norm of the right hand side.
!
!    Output, real RESID, the computed stopping measure,
!    the 2 norm of the residual divided by the 2 norm of the
!    right hand side vector B.
!
!    Input, real RESTOL, the allowable convergence measure.
!
!    Input/output, integer INFO.
!
!    On input, if INFO is -1,then BNRM2 will be calculated.
!
!    On exit, 1/0 depending on whether stopping criterion
!    was met or not.
!
  implicit none

  integer n

  real bnrm2
  real snrm2
  integer info
  real r(n)
  real resid
  real restol

  resid = snrm2 ( n, r, 1 )

  if ( bnrm2 /= 0.0E+00 ) then
    resid = resid / bnrm2
  end if

  if ( resid <= restol ) then
    info = 1
  else
    info = 0
  end if

  return
end
subroutine stopx ( n, r, x, xnrm2, resid, restol, info )

!*****************************************************************************80
!
!! STOPX computes the stopping criterion on X.
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real R(N), the residual, A * X - B.
!
!    Input, real X(N), the solution estimate.
!
!    Output, real XNRM2, the euclidean norm of X.
!
!    Output, real RESID, the computed stopping measure,
!    the euclidean norm of the residual divided by the euclidean
!    norm of the current solution estimate.
!
!    Input, real RESTOL, the convergence tolerance.
!
!    Output, integer INFO.
!    0, the convergence tolerance was not met;
!    1, the convergence tolerance was met.
!
  implicit none

  integer n

  real snrm2
  integer info
  real r(n)
  real resid
  real restol
  real x(n)
  real xnrm2

  xnrm2 = snrm2 ( n, x, 1 )

  if ( xnrm2 == 0.0E+00 ) then
    xnrm2 = 1.0E+00
  end if

  resid = snrm2 ( n, r, 1 ) / xnrm2

  if ( resid <= restol ) then
    info = 1
  else
    info = 0
  end if

  return
end
subroutine sut_mxv ( a, x, b, m, n )

!*****************************************************************************80
!
!! SUT_MXV computes A * x, where A is an upper triangular matrix.
!
!  Modified:
!
!    05 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(N,N), the M by N upper triangular matrix, stored
!    in LINPACK general matrix storage.
!
!    Input, real X(N), the vector to be multiplied by A.
!
!    Output, real B(M), the product A * x.
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
  implicit none

  integer m
  integer n

  real a(n,n)
  real b(m)
  integer i
  integer ierror
  integer j
  double precision temp
  real x(n)
!
!  Check the dimensions.
!
  call sge_check ( m, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SUT_MXV - Fatal error!'
    write ( *, * ) '  Illegal dimensions!'
    return
  end if

  do i = 1, m
    temp = 0.0E+00
    do j = i, n
      temp = temp + dble ( a(i,j) ) * dble ( x(j) )
    end do
    b(i) = sngl ( temp )
  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine update ( i, n, x, h, ldh, y, s, v, ldv )

!*****************************************************************************80
!
!! UPDATE ??
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, integer I, ...
!
!    Input, integer N, the order of the matrix.
!
!    Output, real X(N), ...
!
!    Input, real H(LDH,*), ...
!
!    Input, integer LDH, the leading dimension of H.
!
!    Workspace, real Y(*), ...
!
!    Input, real S(*), ...
!
!    Input, real V(LDV,*), ...
!
!    Input, integer LDV, the leading dimension of V.
!
  implicit none

  integer ldh
  integer ldv
  integer n

  real h(ldh,*)
  integer i
  integer j
  real s(*)
  real v(ldv,*)
  real x(n)
  real y(*)
!
!  Solve H*y = s for upper triangular H.
!
  y(1:i) = s(1:i)

  call strsv ( 'upper', 'notrans', 'nonunit', i, h, ldh, y, 1 )
!
!  Compute the current solution vector X.
!
  do j = 1, i
    call saxpy ( n, y(j), v(1,j), 1, x, 1 )
  end do

  return
end
subroutine vecgen ( form, n, a, b )

!*****************************************************************************80
!
!! VECGEN generates a vector of all ones, zeros, or the row sum
!  of the matrix A.  In the last case, if the vector is used as
!  the right hand side of a linear system,
!
!    A*x = b,
!
!  then the solution is
!
!    x = (1,1,...,1)
!
!  Modified:
!
!    19 August 1999
!
!  Parameters:
!
!    Input, character ( len = 4 ) FORM.
!    Indicates which output vector is desired.
!    'ONES', set B to (1,1,...,1)
!    'ZERO', set B to (0,0,...,0)
!    'SUMR', set B to A*(1,1,...,1)^T.
!
!    Input, integer N, the order of the matrix and the dimension of B.
!
!    Input, real A(N,N), the matrix used when
!    FORM = 'SUMR'.
!
!    Output, real B(N), the output vector.
!
  implicit none

  integer n

  real a(n,n)
  real b(n)
  character ( len = 4 ) form
  integer i
  integer j
  logical lsamen
  real tmp

  if ( lsamen ( 3, form, 'ones' ) ) then

    b(1:n) = 1.0E+00

  else if ( lsamen ( 3, form, 'zeros' ) ) then

    b(1:n) = 0.0E+00

  else if ( lsamen ( 3, form, 'sumrow' ) ) then

    do i = 1, n
      b(i) = sum ( a(i,1:n) )
    end do

  else

     write ( *, * ) ' '
     write ( *, * ) 'VECGEN - Fatal error!'
     write ( *, * ) '  Unrecognized option = ' // form
     stop

  end if

  return
end
