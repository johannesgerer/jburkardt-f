subroutine bndacc ( g, mdg, nb, ip, ir, mt, jt )

!*****************************************************************************80
!
!! BNDACC accumulates information for a banded least squares problem.
!
!  Discussion:
!
!    BNDACC is used for the accumulation phase.  Then BNDSOL should be called
!    to solve linear systems.
!
!    The user must set IR = 1 and IP = 1 before the first call
!    to BNDACC for a new case.
!
!    BNDACC is to be called once for each block of data [ C(I), B(I) ]
!    to be introduced into the problem.  For each block of data, the
!    user must assign values to MT and JT, and copy the MT by (NB+1)
!    array of data [ C(I), B(I) ] into rows IR through IR+MT-1 of
!    the working array G.
!
!  Modified:
!
!    21 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) G(MDG,NB+1), the banded matrix being
!    accumulated.
!
!    Input, integer ( kind = 4 ) MDG, the leading dimension of G.
!    MDG must be at least equal to the number of rows.
!
!    Input, integer ( kind = 4 ) NB, the bandwidth of G, not counting the diagonal.
!
!    Input/output, integer ( kind = 4 ) IP, must be set to 1 by the user before
!    the first call.  Thereafter, its value is controlled by the routine.
!
!    Input/output, integer ( kind = 4 ) IR, indicates the index of the first row
!    of G into which new data is to be placed.  The user sets IR to
!    1 before the first call, but does not alter it thereafter.
!    Instead, the program keeps this quantity up to date as more
!    data is read in.
!
!    Input, integer ( kind = 4 ) MT, is set by the user to indicate the number of
!    new rows of data being introduced by the current call.
!    MT must be at least 0.
!
!    Input, integer ( kind = 4 ) JT, set by the user to indicate the column of the
!    submatrix A(I) that is identified with the first column of C(I).
!    This means that JT must be at least 1.
!
  implicit none

  integer ( kind = 4 ) mdg
  integer ( kind = 4 ) nb

  real ( kind = 8 ) g(mdg,nb+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ie
  integer ( kind = 4 ) ig
  integer ( kind = 4 ) ig1
  integer ( kind = 4 ) ig2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jg
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kh
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mh
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) mt
  integer ( kind = 4 ) mu
  real ( kind = 8 ) rho

  if ( mt <= 0 ) then
    return
  end if

  if ( jt /= ip ) then

    if ( ir < jt ) then

      do i = 1, mt
        ig1 = jt+mt-i
        ig2 = ir+mt-i
        do j = 1, nb+1
          g(ig1,j) = g(ig2,j)
        end do
      end do

      do i = 1, jt - ir
        ig = ir+i-1
        g(ig,1:nb+1) = 0.0D+00
      end do

      ir = jt

    end if

    mu = min ( nb-1, ir-ip-1 )

    do l = 1, mu

      k = min ( l, jt-ip )

      ig = ip+l
      do i = l+1, nb
        g(ig,i-k) = g(ig,i)
      end do

      g(ig,nb+1-k:nb) = 0.0D+00

    end do

    ip = jt

  end if

  mh = ir + mt - ip
  kh = min ( nb+1, mh )
  mode = 1

  do i = 1, kh

    if ( i+1 <= nb ) then

      call h12 ( mode, i, max ( i+1, ir-ip+1 ), mh, g(ip,i), 1, rho, &
        g(ip,i+1), 1, mdg, nb+1-i)

    else

      call h12 ( mode, i, max ( i+1, ir-ip+1 ), mh, g(ip,i), 1, rho, &
        g(ip,1), 1, mdg, 0 )

    end if

  end do
  
  ir = ip + kh

  if ( nb + 1 <= kh ) then
    g(ip+kh-1,1:nb) = 0.0D+00
  end if

  return
end
subroutine bndsol ( mode, g, mdg, nb, ip, ir, x, n, rnorm )

!*****************************************************************************80
!
!! BNDSOL solves a banded least squares problem accumulated by BNDACC.
!
!  Modified:
!
!    21 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MODE, specifies what system is to be solved.
!    1, solve R * X = Y, where R and Y are in stored in G.
!    2, solve R' * X = Y, where R is in G, Y in X.
!    3, solve R * X = Y, where R is in G, Y is in X.
!
!    Input, real ( kind = 8 ) G(MDG,NB+1), the banded matrix.
!
!    Input, integer ( kind = 4 ) MDG, the leading dimension of G.
!    MDG must be at least equal to the number of rows.
!
!    Input, integer ( kind = 4 ) NB, the bandwidth of G, not counting the diagonal.
!
!    Input, integer ( kind = 4 ) IP, IR, must have the values returned by BNDACC
!    from its last call.
!
!    Output, real ( kind = 8 ) X(N), the computed solution.
!
!    Input, integer ( kind = 4 ) N, the dimension of the solution.  BNDSOL will
!    use just the N by N logical submatrix of G.
!
!    Output, real ( kind = 8 ) RNORM, is set if MODE was 1 on input.
!    Its value is the L2 norm of G(N+1:IR-1,NB+1).
!
  implicit none

  integer ( kind = 4 ) mdg
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb

  real ( kind = 8 ) g(mdg,nb+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ie
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jg
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mode
  real ( kind = 8 ) rnorm
  real ( kind = 8 ) rsq
  real ( kind = 8 ) s
  real ( kind = 8 ) x(n)

  rnorm = 0.0D+00

  if ( mode == 1 ) then

    x(1:n) = g(1:n,nb+1)

    rsq = sqrt ( sum ( g(n+1:ir-1,nb+1)**2 ) )

  end if

  if ( mode == 1 .or. mode == 3 ) then

    do ii = 1, n
      i = n+1-ii

      s = 0.0D+00
      l = max ( 0, i-ip )

      ie = min ( n+1-i, nb )
      do j = 2, ie
        jg = j+l
        ix = i-1+j
        s = s + g(i,jg) * x(ix)
      end do

      if ( g(i,l+1) == 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BNDSOL - Fatal error!'
        write ( *, '(a)' ) '  Zero diagonal term.'
        write ( *, '(a,i6)' ) '  MODE = ', mode
        write ( *, '(a,i6)' ) '  I = ', i
        write ( *, '(a,i6)' ) '  J = ', j
        write ( *, '(a,i6)' ) '  L = ', l
        stop
      end if

      x(i) = ( x(i) - s ) / g(i,l+1)

    end do

  else if ( mode == 2 ) then

    do j = 1, n

      s = 0.0D+00

      if ( j /= 1 ) then
        i1 = max ( 1, j-nb+1 )
        do i = i1, j-1
          l = j - i + 1 + max ( 0, i-ip )
          s = s + x(i) * g(i,l)
        end do
      end if

      l = max ( 0, j-ip )

      if ( g(j,l+1) == 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BNDSOL - Fatal error!'
        write ( *, '(a)' ) '  Zero diagonal term.'
        write ( *, '(a,i6)' ) '  I = ', i
        write ( *, '(a,i6)' ) '  J = ', j
        write ( *, '(a,i6)' ) '  L = ', l
        stop
      end if

      x(j) = ( x(j) - s ) / g(j,l+1)

    end do

  end if

  return
end
function diff ( x, y )

!*****************************************************************************80
!
!! DIFF is used in tests that depend on machine precision.
!
!  Modified:
!
!    21 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the numbers whose difference is desired.
!
!    Output, real ( kind = 8 ) DIFF, the difference of X and Y.
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  diff = x - y

  return
end
subroutine g1 ( a, b, c, s, sig )

!*****************************************************************************80
!
!! G1 computes an orthogonal rotation matrix.
!
!  Discussion:
!
!    The routine is given a vector
!
!      ( A )
!      ( B )
!
!    and computes quantities C and S so that
!
!      (  C  S )   ( A )   ( SQRT ( A^2 + B^2 ) )
!      ( -S  C ) * ( B ) = ( 0                  ) 
!
!  Modified:
!
!    21 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the entries of the vector.
!
!    Output, real ( kind = 8 ) C, S, the entries of the rotation matrix.
!
!    Output, real ( kind = 8 ) SIG, the value of SQRT ( A**2 + B**2 ).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) s
  real ( kind = 8 ) sig
  real ( kind = 8 ) xr
  real ( kind = 8 ) yr

  if ( abs ( b ) < abs ( a ) ) then
    xr = b / a
    yr = sqrt ( 1.0D+00 + xr**2 )
    c = sign ( 1.0D+00 / yr, a )
    s = c * xr
    sig = abs ( a ) * yr
  else if ( b /= 0.0D+00 ) then
    xr = a / b
    yr = sqrt ( 1.0D+00 + xr**2 )
    s = sign ( 1.0D+00 / yr, b )
    c = s * xr
    sig = abs ( b ) * yr
  else
    c = 0.0D+00
    s = 1.0D+00
    sig = 0.0D+00
  end if

  return
end
subroutine g2 ( c, s, x, y )

!*****************************************************************************80
!
!! G2 applies a rotation matrix to a vector (X,Y).
!
!  Discussion:
!
!    This routine is no longer used by other routines in the package,
!    but is included for completeness.
!
!  Modified:
!
!    21 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, S, the cosine and sine of the rotation.
!
!    Input/output, real ( kind = 8 ) X, Y, the components of the vector
!    to be rotated.
!
  implicit none

  real ( kind = 8 ) c
  real ( kind = 8 ) s
  real ( kind = 8 ) x
  real ( kind = 8 ) x_new
  real ( kind = 8 ) y
  real ( kind = 8 ) y_new

  x_new =   c * x + s * y
  y_new = - s * x + c * y

  x = x_new
  y = y_new

  return
end
function gen ( anoise )

!*****************************************************************************80
!
!! GEN generates numbers for construction of test cases.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANOISE, determines the level of
!    "noise" to be added to the data.
!
!    Output, real ( kind = 8 ) GEN, a random value with noise added.
!
  implicit none

  real ( kind = 8 ), save :: aj = 0.0D+00
  real ( kind = 8 ) anoise
  real ( kind = 8 ) gen
  integer ( kind = 4 ), save :: i = 5
  integer ( kind = 4 ), save :: j = 7
  integer ( kind = 4 ), parameter :: mi = 891
  integer ( kind = 4 ), parameter :: mj = 457

  if ( anoise < 0.0D+00 ) then
    i = 5
    j = 7
    aj = 0.0D+00
    gen = 0.0D+00
    return
  end if
!
!  The sequence of values of J is bounded between 1 and 996.
!  If initial j = 1,2,3,4,5,6,7,8, or 9, the period is 332
!
  if ( 0.0D+00 < anoise ) then
    j = j * mj
    j = j - 997 * ( j / 997 )
    aj = dble ( j - 498 )
  end if
!
!  The sequence of values of I is bounded between 1 and 999.
!  If initial i = 1,2,3,6,7, or 9,  the period will be 50.
!  If initial i = 4 or 8   the period will be 25.
!  If initial i = 5        the period will be 10.
!
  i = i * mi
  i = i - 1000 * ( i / 1000 )

  gen = dble ( i - 500 ) + aj * anoise

  return
end
subroutine h12 ( mode, lpivot, l1, m, u, iue, up, c, ice, icv, ncv )

!*****************************************************************************80
!
!! H12 constructs or applies a Householder transformation.
!
!  Discussion:
!
!    The symmetric orthogonal Householder transformation matrix
!    has the form
!
!      Q = I + u * u' / b
!
!  Modified:
!
!    21 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MODE, selects the algorithm.
!    1, algorithm H1 to construct and apply a Householder transformation, 
!    2, algorithm H2 to apply a previously constructed transformation.
!
!    Input, integer ( kind = 4 ) LPIVOT, is the index of the pivot element, which
!    should satisfy 1 <= LPIVOT < L1.
!
!    Input, integer ( kind = 4 ) L1, M.   If L1 <= M, the transformation will be 
!    constructed to zero out elements indexed from L1 through M.   
!    But if M < L1, the subroutine does an identity transformation.
!
!    Input/output, real ( kind = 8 ) U(IUE,M), integer IUE, real ( kind = 8 ) UP.
!    If MODE = 1:
!      On input, 
!        U contains the pivot vector;  
!        IUE is the storage increment between elements.
!      On output, 
!        U and UP contain quantities defining the vector U of the
!        Householder transformation.
!    If MODE = 2:
!      On input, U and UP should contain the values previously computed
!      with MODE = 1, which will not be modified during this call.
!
!    Input/output, real ( kind = 8 ) C(*).  On input, C contains a matrix 
!    which will be regarded as a set of vectors to which the
!    Householder transformation is to be applied.
!    On output, C contains the set of transformed vectors.
!
!    Input, integer ( kind = 4 ) ICE, the storage increment between successive elements 
!    of a single vector in C.
!
!    Input, integer ( kind = 4 ) ICV, the storage increment between the first element
!    of succesive vectors in C.
!
!    Input, integer ( kind = 4 ) NCV, the number of vectors in C to be transformed.
!    If NCV <= 0, no operations are done on C.
!
  implicit none

  integer ( kind = 4 ) iue
  integer ( kind = 4 ) m

  real ( kind = 8 ) b
  real ( kind = 8 ) c(*)
  real ( kind = 8 ) cl
  real ( kind = 8 ) clinv
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) ice
  integer ( kind = 4 ) icv
  integer ( kind = 4 ) incr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) lpivot
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) ncv
  real ( kind = 8 ) sm
  real ( kind = 8 ) u(iue,m)
  real ( kind = 8 ) up

  if ( lpivot <= 0 ) then
    return
  else if ( l1 <= lpivot ) then
    return
  else if ( m < l1 ) then
    return
  end if

  cl = abs ( u(1,lpivot) )
!
!  Construct the transformation.
!
  if ( mode == 1 ) then

    do j = l1, m
      cl = max ( cl, abs ( u(1,j) ) )
    end do

    if ( cl <= 0.0D+00 ) then
      return
    end if

    clinv = 1.0D+00 / cl
    sm = ( u(1,lpivot) * clinv )**2
    do j = l1, m
      sm = sm + ( u(1,j) * clinv )**2
    end do
    cl = cl * sqrt ( sm )

    if ( 0.0D+00 < u(1,lpivot) ) then
      cl = - cl
    end if

    up = u(1,lpivot) - cl
    u(1,lpivot) = cl

  else if ( mode == 2 ) then

    if ( cl <= 0.0D+00 ) then
      return
    end if

  end if
!
!  Apply the transformation I+U*U'/B to the NCV vectors in C.
!
  if ( ncv <= 0 ) then
    return
  end if

  b = up * u(1,lpivot)
!
!  B must be nonpositive.
!
  if ( 0.0D+00 <= b ) then
    return
  end if

  b = 1.0D+00 / b
  i2 = 1 - icv + ice * ( lpivot - 1 )
  incr = ice * ( l1 - lpivot )

  do j = 1, ncv

    i2 = i2 + icv
    i3 = i2 + incr
    i4 = i2 + incr

    sm = c(i2) * up
    do i = l1, m
      sm = sm + c(i3) * u(1,i)
      i3 = i3 + ice
    end do

    if ( sm /= 0.0D+00 ) then
      sm = sm * b
      c(i2) = c(i2) + sm * up
      do i = l1, m
        c(i4) = c(i4) + sm * u(1,i)
        i4 = i4 + ice
      end do
    end if

  end do

  return
end
subroutine hfti ( a, mda, m, n, b, mdb, nb, tau, krank, rnorm, h, g, ip )

!*****************************************************************************80
!
!! HFTI: Householder forward triangulation with column interchanges.
!
!  Modified:
!
!    21 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(MDA,N), the M by N matrix.
!
!    Input, integer ( kind = 4 ) MDA, the leading dimension of A.  
!    MDA must be at least M.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in A.
!
!    Input/output, real ( kind = 8 ) B(MDB,NB), on input, the NB right hand 
!    side vectors of length M, and on output, the NB least squares solution
!    vectors of length N.
!
!    Input, integer ( kind = 4 ) MDB, the leading dimension of B, which 
!    must be at least M.
!
!    Input, integer ( kind = 4 ) NB, the number of columns in B.
!
!    Input, real ( kind = 8 ) TAU, a tolerance used in the pseudorank test.
!
!    Output, integer ( kind = 4 ) KRANK, the numerically determined rank of the
!    matrix A.  
!
!    Output, real ( kind = 8 ) RNORM(NB), the Euclidean norms of the residual
!    vectors for the least squares problems.
!
!    Workspace, real ( kind = 8 ) H(N).
!
!    Workspace, real ( kind = 8 ) G(N).
!
!    Workspace, integer IP(N).
!
  implicit none

  integer ( kind = 4 ) mda
  integer ( kind = 4 ) mdb
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb

  real ( kind = 8 ) a(mda,n)
  real ( kind = 8 ) b(mdb,nb)
  real ( kind = 8 ) diff
  real ( kind = 8 ), parameter :: factor = 0.001D+00
  real ( kind = 8 ) g(n)
  real ( kind = 8 ) h(n)
  real ( kind = 8 ) hmax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ip(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jb
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) krank
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ldiag
  integer ( kind = 4 ) lmax
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mode
  real ( kind = 8 ) rnorm(nb)
  real ( kind = 8 ) sm
  real ( kind = 8 ) tau
  real ( kind = 8 ) tmp

  k = 0
  ldiag = min ( m, n )

  if ( ldiag <= 0 ) then
    krank = 0
    return
  end if

  hmax = 0.0D+00
  h(1:n) = 0.0D+00

  do j = 1, ldiag
!
!  Update the squared column lengths and find LMAX.
!  In the original code, LMAX was not set if J <= 1, causing errors.
!
    if ( j <= 1 ) then

      lmax = 1

    else

      lmax = j
      do l = j, n
        h(l) = h(l) - a(j-1,l)**2
        if ( h(lmax) < h(l) ) then
          lmax = l
        end if
      end do

    end if
!
!  Compute the squared column lengths and find LMAX.
!
    if ( diff ( hmax + factor * h(lmax), hmax ) <= 0.0D+00 ) then

      lmax = j

      do l = j, n
        h(l) = 0.0D+00
        do i = j, m
          h(l) = h(l) + a(i,l)**2
        end do
        if ( h(lmax) < h(l) ) then
          lmax = l
        end if
      end do

      hmax = h(lmax)

    end if
!
!  Do column interchanges if needed.
!
    ip(j) = lmax

    if ( ip(j) /= j ) then

      do i = 1, m
        tmp       = a(i,j)
        a(i,j)    = a(i,lmax)
        a(i,lmax) = tmp
      end do
      h(lmax) = h(j)

    end if
!
!  Compute the J-th transformation and apply it to A and B.
!
    mode = 1

    if ( j < n ) then
      call h12 ( mode, j, j+1, m, a(1,j), 1, h(j), a(1,j+1), 1, mda, n-j )
    else
      call h12 ( mode, j, j+1, m, a(1,j), 1, h(j), a(1,1), 1, mda, n-j )
    end if

    mode = 2

    call h12 ( mode, j, j+1, m, a(1,j), 1, h(j), b, 1, mdb, nb )

  end do
!
!  Determine the pseudorank K using the tolerance TAU.
!
  k = ldiag
  do j = 1, ldiag
    if ( abs ( a(j,j) ) <= tau ) then
      k = j - 1
      exit
    end if
  end do
!
!  Compute the norms of the residual vectors.
!
  do jb = 1, nb
    rnorm(jb) = sqrt ( sum ( b(k+1:m,jb)**2 ) )
  end do
!
!  Special for pseudorank 0.
!
  if ( k <= 0 ) then
    b(1:n,1:nb) = 0.0D+00
    krank = k
    return
  end if
!
!  If the pseudorank is less than N, compute the Householder
!  decomposition of the first K rows.
!
  if ( k /= n ) then
    mode = 1
    do ii = 1, k
      i = k+1-ii
      call h12 ( mode, i, k+1, n, a(i,1), mda, g(i), a, mda, 1, i-1 )
    end do
  end if

  do jb = 1, nb
!
!  Solve the K by K triangular system.
!
    do l = 1, k
      sm = 0.0D+00
      i = k+1-l
      do j = i+1, k
        sm = sm + a(i,j) * b(j,jb)
      end do
      b(i,jb) = ( b(i,jb) - sm ) / a(i,i)
    end do
!
!  Complete computation of solution vector.
!
    if ( k < n ) then

      b(k+1:n,jb) = 0.0D+00

      mode = 2
      do i = 1, k
        call h12 ( mode, i, k+1, n, a(i,1), mda, g(i), b(1,jb), 1, mdb, 1 )
      end do

    end if
!
!  Reorder the solution vector to compensate for the column interchanges.
!
    do jj = 1, ldiag
      j = ldiag+1-jj
      if ( ip(j) /= j) then
        l = ip(j)
        tmp     = b(l,jb)
        b(l,jb) = b(j,jb)
        b(j,jb) = tmp
      end if
    end do

  end do
!
!  The solution vectors X are now in the first N rows of B.
!
  krank = k

  return
end
subroutine ldp ( g, mdg, m, n, h, x, xnorm, w, indx, mode )

!*****************************************************************************80
!
!! LDP implements least distance programming algorithm.
!
!  Discussion:
!
!    The problem is to determine an N vector X with minimum norm,
!    subject to 
!
!      H <= G * X
!
!    where G is an M by N matrix, and H is an M vector.
!
!  Modified:
!
!    21 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) G(MDG,N), the M by N matrix.
!    There is no restriction on the rank of G.
!
!    Input, integer ( kind = 4 ) MDG, the leading dimension of G.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in G.
!    M may be greater or less than N.
!
!    Input, real ( kind = 8 ) H(M), the right hand side of the constraint.
!
!    Output, real ( kind = 8 ) X(N), the computed solution.
!
!    Output, real ( kind = 8 ) XNORM, the value of the norm of X.
!
!    Workspace, real ( kind = 8 ) W((N+1)*(M+2)+2*M), used for NNLS.
!    First (N+1)*M entries = matrix E,
!    Next N+1 entries = vector F,
!    Next N+1 entries = vector Z,
!    Next M entries = vector Y,
!    Next M entries = vector WDUAL.
!
!    Output, integer ( kind = 4 ) INDX(M), defines the sets P and Z as follows:
!    INDX(1:NSETP) = set P, INDX(NSETP+1:N) = set Z.    
!
!    Output, integer ( kind = 4 ) MODE, error flag.
!    1, no error.
!    2, bad dimension.
!    3, maximum number of iterations taken without convergence.
!    4, the inequality constraints H <= G * X are incompatible.
!
  implicit none

  integer ( kind = 4 ) mdg
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) diff
  real ( kind = 8 ) fac
  real ( kind = 8 ) g(mdg,n)
  real ( kind = 8 ) h(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(m)
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) iwdual
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) iz
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jf
  integer ( kind = 4 ) mode
  real ( kind = 8 ) rnorm
  real ( kind = 8 ) w((n+1)*(m+2)+2*m)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xnorm

  if ( n <= 0 ) then
    mode = 2
    return
  end if

  x(1:n) = 0.0D+00
  xnorm = 0.0D+00

  if ( m <= 0 ) then
    mode = 1
    return
  end if

  iw = 0
  do j = 1, m
    do i = 1, n
      iw = iw + 1
      w(iw) = g(j,i)
    end do
    iw = iw + 1
    w(iw) = h(j)
  end do
  jf = iw + 1
!
!  Store N zeros followed by a 1 into F.
!
  do i = 1, n
    iw = iw + 1
    w(iw) = 0.0D+00
  end do
  w(iw+1) = 1.0D+00

  iz = iw + 2
  iy = iz + n + 1
  iwdual = iy + m

  call nnls ( w, n+1, n+1, m, w(jf), w(iy), rnorm, w(iwdual), w(iz), &
    indx, mode )

  if ( mode /= 1 ) then
    return
  end if

  if ( rnorm <= 0.0D+00 ) then
    mode = 4
    return
  end if

  fac = 1.0D+00

  iw = iy - 1
  do i = 1, m
    iw = iw + 1
    fac = fac - h(i) * w(iw)
  end do

  if ( diff ( 1.0D+00 + fac, 1.0D+00 ) == 0.0D+00 ) then
    mode = 4
    return
  end if

  fac = 1.0D+00 / fac

  do j = 1, n
    iw = iy - 1
    do i = 1, m
      iw = iw + 1
      x(j) = x(j) + g(i,j) * w(iw)
    end do
  end do

  x(1:n) = x(1:n) * fac

  xnorm = sqrt ( sum ( x(1:n)**2 ) )

  mode = 1

  return
end
subroutine mfeout ( a, mda, m, n, names, mode, unit, width )

!*****************************************************************************80
!
!! MFEOUT labeled matrix output for use with singular value analysis.
!
!  Discussion:
!
!    This 1995 version has additional arguments, UNIT and WIDTH,
!    to support user options regarding the output unit and the width of
!    print lines.   It also allows the user to choose the length of
!    names in the NAMES array.
!
!  Modified:
!
!    21 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(MDA,N), the M by N matrix to be output.
!
!    Input, integer ( kind = 4 ) MDA, the leading dimension of A.  
!    MDA must be at least M.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in A.
!
!    Input, character ( len = * ) NAMES(M), an array of names.
!    If NAMES(1) contains only blanks, the rest of the NAMES
!    array will be ignored.
!
!    Input, integer ( kind = 4 ) MODE.
!    = 1  write header for V matrix and use an F format.
!    = 2  write header for candidate solutions and use P format.
!
!    Input, integer ( kind = 4 ) UNIT, the output unit.  If 0 <= UNIT, then output is
!    to Fortran unit UNIT.  If UNIT == -1, output is to unit "*".
!
!    Input, integer ( kind = 4 ) WIDTH, the width of the output lines.  Each output line 
!    from this subroutine will have at most
!      max ( 26, min ( 124, WIDTH ) ) 
!    characters plus one additional leading character for fortran "carriage
!    control".  The carriage control character will always be a blank.
!
  implicit none

  integer ( kind = 4 ) mda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mda,n)
  logical blknam
  character ( len = 26 ), dimension ( 2 ) :: fmt1 = (/ &
    '(/7x,00x,8(5x,a4,i4,1x)/)', &
    '(/7x,00x,8(2x,a4,i4,4x)/)' /)
  character ( len = 26 ), dimension ( 2 ) :: fmt2 = (/ &
    '(1x,i4,1x,a00,1x,4p8f14.0)', &
    '(1x,i4,1x,a00,1x,8g14.6  )' /)
  character ( len = 4 ), dimension ( 2 ) :: head = (/ ' Col', 'Soln' /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) kblock
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lennam
  integer ( kind = 4 ) m
  integer ( kind = 4 ) maxcol
  integer ( kind = 4 ) mode
  character ( len = * ) names(m)
  integer ( kind = 4 ) namsiz
  integer ( kind = 4 ) nblock
  logical star
  integer ( kind = 4 ) unit
  integer ( kind = 4 ) width

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  star = unit < 0
  lennam = len ( names(1) )
  blknam = names(1) == ' '
  namsiz = 1

  if ( .not. blknam) then
    do i = 1, m
      namsiz = max ( namsiz, len_trim ( names(i) ) )
    end do
  end if

  write ( fmt1(mode)(6:7), '(i2.2)' ) namsiz
  write ( fmt2(mode)(12:13), '(i2.2)' ) namsiz

  if ( star ) then
    write ( *, '(a)' ) ' '
    if ( mode == 1 ) then
      write ( *, '(a)' ) '  V-matrix of the SVD of A*D.'
      write ( *, '(a)' ) '  (Elements of V scaled up by a factor of 10**4)'
    else
      write ( *, '(a)' ) '  Sequence of candidate solutions, X'
    end if
  else
    write ( unit, '(a)' ) ' '
    if ( mode == 1 ) then
      write ( unit, '(a)' ) '  V-matrix of the SVD of A*D.'
      write ( unit, '(a)' ) '  (Elements of V scaled up by a factor of 10**4)'
    else
      write ( unit, '(a)' ) '  Sequence of candidate solutions, X'
    end if
  end if
!
!  With NAMSIZ characters allowed for the "name" and MAXCOL
!  columns of numbers, the total line width, exclusive of a
!  carriage control character, will be 6 + lennam + 14 * maxcol.
!
  maxcol = max ( 1, min ( 8, ( width - 6 - namsiz ) / 14 ) )

  nblock = ( n + maxcol - 1 ) / maxcol
  j2 = 0

  do kblock = 1, nblock

     j1 = j2 + 1
     j2 = min ( n, j2 + maxcol )

     if ( star ) then
       write (*,fmt1(mode)) (head(mode),j,j = j1,j2)
     else
       write (unit,fmt1(mode)) (head(mode),j,j = j1,j2)
     end if

     do i = 1,m
        if ( star ) then
           if ( blknam ) then
             write (*,fmt2(mode)) i,' ',(a(i,j),j = j1,j2)
           else
             write (*,fmt2(mode)) i,names(i),(a(i,j),j = j1,j2)
           end if
        else
           if ( blknam ) then
             write (unit,fmt2(mode)) i,' ',(a(i,j),j = j1,j2)
           else
             write (unit,fmt2(mode)) i,names(i),(a(i,j),j = j1,j2)
           end if
        end if
      end do

  end do

  return
end
subroutine nnls ( a, mda, m, n, b, x, rnorm, w, zz, indx, mode )

!*****************************************************************************80
!
!! NNLS implements the nonnegative least squares algorithm.
!
!  Discussion:
!
!    Given an M by N matrix A and an M vector B, compute an
!    N vector X that solves the least squares problem (that is,
!    minimizes the L2 norm of the residual)
!
!      A * X = B  
!
!    subject to 
!
!      0 <= X.
!
!    Householder transformations are applied, which essentially
!    multiply the system by an orthogonal matrix Q:
!
!      Q * A * X = Q * B
!
!  Modified:
!
!    21 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(MDA,N).  On input, the M by N matrix.
!    On output, the product matrix Q * A, where Q is an M by M orthogonal
!    matrix generated implicitly by this routine.
!
!    Input, integer ( kind = 4 ) MDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in A.
!
!    Input/output, real ( kind = 8 ) B(M).  On input, contains the right
!    hand side vector.  On output, contains Q * B.
!
!    Output, real ( kind = 8 ) X(N), the computed solution.
!
!    Output, real ( kind = 8 ) RNORM, the L2 or Euclidean norm of the residual.
!
!    Workspace, real ( kind = 8 ) W(N).  On output, the dual solution
!    vector.  W will satisfy W(I) = 0 for I in set P, and W(I) <= 0
!    for I in set Z.  
!
!    Workspace, real ( kind = 8 ) Z(M).
!
!    Output, integer ( kind = 4 ) INDX(N), defines the sets P and Z as follows:
!    INDX(1:NSETP) = set P, INDX(NSETP+1:N) = set Z.      
!
!    Output, integer ( kind = 4 ) MODE, error flag.
!    1, the solution has been computed successfully.
!    2, the dimensions of the problem are bad, M <= 0 or N <= 0.
!    3, iteration count exceeded.  More than 3*N iterations.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mda,n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) asave
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) cc
  real ( kind = 8 ) diff
  real ( kind = 8 ) dummy(1)
  real ( kind = 8 ), parameter :: factor = 0.01D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) itmax
  integer ( kind = 4 ) iz
  integer ( kind = 4 ) iz1
  integer ( kind = 4 ) iz2
  integer ( kind = 4 ) izmax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jz
  integer ( kind = 4 ) l
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) npp1
  integer ( kind = 4 ) nsetp
  real ( kind = 8 ) rnorm
  integer ( kind = 4 ) rtnkey
  real ( kind = 8 ) sm
  real ( kind = 8 ) ss
  real ( kind = 8 ) t
  real ( kind = 8 ) temp
  real ( kind = 8 ) unorm
  real ( kind = 8 ) up
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) wmax
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) ztest
  real ( kind = 8 ) zz(m)

  mode = 1

  if ( m <= 0 .or. n <= 0 ) then
    mode = 2
    return
  end if

  iter = 0
  itmax = 3 * n
!
!  Initialize INDX and X.
!
  x(1:n) = 0.0D+00

  do i = 1, n
    indx(i) = i
  end do

  iz2 = n
  iz1 = 1
  nsetp = 0
  npp1 = 1
!
!  The main loop begins here.
!
30 continue
!
!  Quit if all coefficients are already in the solution.
!  or if M columns of A have been triangularized.
!
  if ( iz2 < iz1 .or. m <= nsetp ) then
    go to 350
  end if
!
!  Compute components of the dual (negative gradient) vector W.
!
  do iz = iz1, iz2
    j = indx(iz)
    w(j) = dot_product ( b(npp1:m), a(npp1:m,j) )
  end do
!
!  Find the largest positive W.
!
60 continue

  wmax = 0.0D+00
  do iz = iz1, iz2
    j = indx(iz)
    if ( wmax < w(j) ) then
      wmax = w(j)
      izmax = iz
    end if
  end do
!
!  If WMAX <= 0, go to termination.
!  This indicates satisfaction of the Kuhn-Tucker conditions.
!
  if ( wmax <= 0.0D+00 ) then
    go to 350
  end if

  iz = izmax
  j = indx(iz)
!
!  The sign of W(j) is OK for J to be moved to set P.
!
!  Begin the transformation and check the new diagonal element to avoid
!  near linear dependence.
!
  asave = a(npp1,j)

  call h12 ( 1, npp1, npp1+1, m, a(1,j), 1, up, dummy, 1, 1, 0 )

  unorm = sqrt ( sum ( a(1:nsetp,j)**2 ) )

  if ( 0.0D+00 <  diff ( unorm + abs ( a(npp1,j) ) * factor, unorm ) ) then
!
!  Column J is sufficiently independent.  
!
!  Copy B into ZZ, update ZZ and solve for ZTEST, the proposed new 
!  value for X(J).
!
    zz(1:m) = b(1:m)

    call h12 ( 2, npp1, npp1+1, m, a(1,j), 1, up, zz, 1, 1, 1 )

    ztest = zz(npp1) / a(npp1,j)
!
!  See if ZTEST is positive.
!
    if ( 0.0D+00 < ztest ) then
      go to 140
    end if

  end if
!
!  Reject J as a candidate to be moved from set Z to set P.
!  Restore A(NPP1,J).  Set W(J) = 0.  Loop back to test the dual
!  coefficients again.
!
  a(npp1,j) = asave
  w(j) = 0.0D+00
  go to 60
!
!  The index J = INDX(IZ)  has been selected to be moved from
!  set Z to set P.    
!
!  Update B, update indices,  Apply Householder
!  transformations to columns in new set Z,  zero subdiagonal elements in
!  column J, set W(J) = 0.
!
140 continue

  b(1:m) = zz(1:m)
  indx(iz) = indx(iz1)
  indx(iz1) = j
  iz1 = iz1+1
  nsetp = npp1
  npp1 = npp1+1

  if ( iz1 <= iz2) then
    do jz = iz1,iz2
      jj = indx(jz)
      call h12 ( 2, nsetp, npp1, m, a(1,j), 1, up, a(1,jj), 1, m, 1 )
    end do
  end if

  if ( nsetp /= m) then
    a(npp1:m,j) = 0.0D+00
  end if

  w(j) = 0.0D+00
!
!  Solve the triangular system.
!  Store the solution temporarily in ZZ.
!
  rtnkey = 1
  go to 400
!
!  Secondary loop begins here.
!
!  Iteration counter.
!
210 continue

  iter = iter + 1

  if ( itmax < iter ) then
    mode = 3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NNLS - Warning!'
    write ( *, '(a)' ) '  Quitting on iteration count.'
    go to 350
  end if
!
!  See if all new constrained coefficients are feasible.
!  If not, compute alpha.
!
  alpha = 2.0D+00

  do ip = 1, nsetp
    l = indx(ip)
    if ( zz(ip) <= 0.0D+00 ) then
      t = -x(l) / ( zz(ip) - x(l) )
      if ( t < alpha ) then
        alpha = t
        jj = ip
      end if
    end if
  end do
!
!  If all new constrained coefficients are feasible, then ALPHA will
!  still be 2.  If so exit from secondary loop to main loop.
!
  if ( alpha == 2.0D+00 ) go to 330
!
!  Otherwise use ALPHA which will be between 0.0D+00 and 1.0 to
!  interpolate between the old X and the new ZZ.
!
  do ip = 1, nsetp
    l = indx(ip)
    x(l) = x(l) + alpha * ( zz(ip) - x(l) )
  end do
!
!  Modify A and B and the index arrays to move coefficient I
!  from set P to set Z.
!
  i = indx(jj)

260 continue

  x(i) = 0.0D+00

  if ( jj /= nsetp ) then

    jj = jj+1

    do j = jj, nsetp

      ii = indx(j)
      indx(j-1) = ii
      call g1 ( a(j-1,ii), a(j,ii), cc, ss, a(j-1,ii) )
!
!  Apply procedure g2 (cc,ss,a(j-1,l),a(j,l))
!
      a(j,ii) = 0.0D+00

      do l = 1, n
        if ( l /= ii ) then
          temp = a(j-1,l)
          a(j-1,l) =  cc * temp + ss * a(j,l)
          a(j,l)   = -ss * temp + cc * a(j,l)
        end if
      end do
!
!  Apply procedure g2 (cc,ss,b(j-1),b(j))
!
      temp = b(j-1)
      b(j-1) =  cc * temp + ss * b(j)
      b(j)   = -ss * temp + cc * b(j)

    end do

  end if

  npp1 = nsetp
  nsetp = nsetp-1
  iz1 = iz1-1
  indx(iz1) = i
!
!  See if the remaining coefficients in set P are feasible.  They should
!  be because of the way ALPHA was determined.
!  If any are infeasible, it is due to round-off error.  Any
!  that are nonpositive will be set to zero
!  and moved from set P to set Z.
!
  do jj = 1, nsetp
    i = indx(jj)
    if ( x(i) <= 0.0D+00 ) then
      go to 260
    end if
  end do
!
!  Copy B into ZZ, then solve again and loop back.
!
  zz(1:m) = b(1:m)

  rtnkey = 2
  go to 400
!
!  End of secondary loop.
!
330 continue

  do ip = 1, nsetp
    i = indx(ip)
    x(i) = zz(ip)
  end do
!
!  All new coefficients are positive.  Loop back to beginning.
!
  go to 30
!
!  End of the main loop
!
!  Come here for termination.
!  Compute the norm of the final residual vector.
!
350 continue

  sm = sum ( b(npp1:m)**2 )

  if ( m < npp1 ) then
    w(1:n) = 0.0D+00
  end if

  rnorm = sqrt ( sm )
  return
!
!  The following block of code is used as an internal subroutine
!  to solve the triangular system, putting the solution in ZZ.
!
400 continue

  do ip = nsetp, 1, -1

    if ( ip /= nsetp ) then
      do ii = 1, ip
        zz(ii) = zz(ii) - a(ii,jj) * zz(ip+1)
      end do
    end if

    jj = indx(ip)
    zz(ip) = zz(ip) / a(ip,jj)

  end do

  go to 210

end
subroutine qrbd ( ipass, q, e, nn, v, mdv, nrv, c, mdc, ncc )

!*****************************************************************************80
!
!! QRBD uses the QR algorithm for the singular values of a bidiagonal matrix.
!
!  Discussion:
!
!    The bidiagonal matrix D, with diagonal entries Q(1:N) and
!    superdiagonal entries E(2:N) is pre- and post-multiplied by
!    elementary rotation matrices R and P so that
!
!      Rk*...*R2*R1 * D * P1'*P2'*...*Pk' = diag(S1,...,Sn)
!
!  Modified:
!
!    21 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson.
!
!  Reference:
!
!    Gene Golub, Christian Reinsch,
!    Singular Value Decomposition and Least Squares Solutions,
!    Numerische Mathematik,
!    Volume 14, Number 5, April 1970, pages 403-420.
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IPASS, error flag.
!    1, no error was detected.
!    2, an error occurred.
!
!    Input/output, real ( kind = 8 ) Q(NN).  On input, the diagonal 
!    entries of the matrix.  On output, the singular values, which 
!    are nonnegative, listed in nonincreasing order.
!
!    Input/output, real ( kind = 8 ) E(NN).  On input, E(2:N) contains
!    the superdiagonal entries of the matrix.
!
!    Input, integer ( kind = 4 ) NN, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) V(MDV,NN), an NRV by NN matrix
!    whose input value is to be postmultiplied by the transformations
!    P1' * P2' * ... * PK.
!
!    Input, integer ( kind = 4 ) MDV, the leading dimension of V.
!
!    Input, integer ( kind = 4 ) NRV, the number of rows in V.
!
!    Input/output, real ( kind = 8 ) C(MDC,NCC), an NN by NCC matrix
!    whose input value is to be premultiplied by the transformations
!    Rm * ... * R2 * R1.
!
!    Input, integer ( kind = 4 ) MDC, the leading dimension of C, which must be
!    at least NN.
!
!    Input, integer ( kind = 4 ) NCC, the number of columns in C.
!
  implicit none

  integer ( kind = 4 ) mdc
  integer ( kind = 4 ) mdv
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nn

  real ( kind = 8 ) c(mdc,ncc)
  real ( kind = 8 ) cs
  real ( kind = 8 ) diff
  real ( kind = 8 ) dnorm
  real ( kind = 8 ) e(nn)
  real ( kind = 8 ) f
  logical fail
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  logical havers
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ipass
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n10
  integer ( kind = 4 ) nqrs
  integer ( kind = 4 ) nrv
  real ( kind = 8 ) q(nn)
  real ( kind = 8 ) small
  real ( kind = 8 ) sn
  real ( kind = 8 ) t
  real ( kind = 8 ) temp
  real ( kind = 8 ) v(mdv,nn)
  logical wntv
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  n = nn
  ipass = 1

  if ( n <= 0 ) then
    return
  end if

  n10 = 10 * n
  wntv = 0 < nrv 
  havers = 0 < ncc
  fail = .false.
  nqrs = 0
  e(1) = 0.0D+00

  dnorm = 0.0D+00
  do j = 1, n
    dnorm = max ( dnorm, abs ( q(j) ) + abs ( e(j) ) )
  end do

  do k = n, 1, -1
!
!  Test for splitting or rank deficiencies.
!  First make test for last diagonal term, Q(K), being small.
!
20  continue

    if ( k == 1 ) then
      go to 50
    end if

    if ( diff ( dnorm+q(k), dnorm ) /= 0.0D+00 ) then
      go to 50
    end if
!
!  Since Q(K) is small we will make a special pass to
!  transform E(k) to zero.
!
!  Transformation constructed to zero out position (I,K).
!
    cs = 0.0D+00
    sn = -1.0D+00

    do i = k-1, 1, -1

      f =     -sn * e(i+1)
      e(i+1) = cs * e(i+1)
      call g1 ( q(i), f, cs, sn, q(i) )
!
!  Accumulate right hand transformations in V.
!
      if ( wntv ) then
        do j = 1, nrv
          temp = v(j,i)
          v(j,i) =   cs * temp + sn * v(j,k)
          v(j,k)  = -sn * temp + cs * v(j,k)
        end do
      end if

    end do
!
!  The matrix is now bidiagonal, and of lower order
!  since E(K) is zero.
!
50  continue

    do l = k, 1, -1

      if ( diff ( dnorm + e(l), dnorm ) == 0.0D+00 ) then
        go to 100
      end if

      if ( diff ( dnorm + q(l-1), dnorm ) == 0.0D+00 ) then
        go to 70
      end if

    end do
!
!  This loop can't complete since E(1) = zero.
!
    go to 100
!
!  Cancellation of E(L), 1 < L.
!
70  continue

    cs = 0.0D+00
    sn = -1.0D+00

    do i = l, k

      f =   -sn * e(i)
      e(i) = cs * e(i)
      if ( diff ( dnorm + f, dnorm ) == 0.0D+00 ) then
        go to 100
      end if

      call g1 ( q(i), f, cs, sn, q(i) )
!
!  Accumulate left hand transformations in C.
!
      if ( havers ) then
        do j = 1, ncc
          temp = c(i,j)
          c(i,j)   =   cs * temp + sn * c(l-1,j)
          c(l-1,j)  = -sn * temp + cs * c(l-1,j)
        end do
      end if

    end do
!
!  Test for convergence.
!
100 continue

    z = q(k)

    if ( l == k ) then
      go to 170
    end if
!
!  Shift from bottom 2 by 2 minor of B'*B.
!
    x = q(l)
    y = q(k-1)
    g = e(k-1)
    h = e(k)
    f = ( (y-z) * (y+z) + (g-h) * (g+h) ) / ( 2.0D+00 * h * y )
    g = sqrt ( 1.0D+00 + f * f )

    if ( 0.0D+00 <= f ) then
      t = f + g
    else
      t = f - g
    end if

    f = ( (x-z) * (x+z) + h * ( y / t - h ) ) / x
!
!  Next QR sweep.
!
    cs = 1.0D+00
    sn = 1.0D+00

    do i = l+1, k

      g = e(i)
      y = q(i)
      h = sn * g
      g = cs * g
      call g1 ( f, h, cs, sn, e(i-1) )
      f =  x * cs + g * sn
      g = -x * sn + g * cs
      h = y * sn
      y = y * cs
!
!  Accumulate rotations from the right in V.
!
      if ( wntv ) then

        do j = 1, nrv
          temp = v(j,i-1)
          v(j,i-1) =   cs * temp + sn * v(j,i)
          v(j,i)    = -sn * temp + cs * v(j,i)
        end do

      end if

      call g1 ( f, h, cs, sn, q(i-1) )

      f =  cs * g + sn * y
      x = -sn * g + cs * y
!
!  Accumulate rotations from the left in C.
!
      if ( havers ) then
        do j = 1, ncc
          temp = c(i-1,j)
          c(i-1,j) =   cs * temp + sn * c(i,j)
          c(i,j)    = -sn * temp + cs * c(i,j)
        end do
      end if

    end do

    e(l) = 0.0D+00
    e(k) = f
    q(k) = x
    nqrs = nqrs + 1

    if ( nqrs <= n10 ) then
      go to 20
    end if
!
!  Return to 'test for splitting'.
!
    small = abs(e(k))
    i = k
!
!  If failure to converge set smallest magnitude
!  term in off-diagonal to zero.  Continue on.
!
    do j = l, k
      temp = abs ( e(j) )
      if ( temp /= 0.0D+00 ) then
        if ( temp < small ) then
          small = temp
          i = j
        end if
      end if
    end do

    e(i) = 0.0D+00
    nqrs = 0
    fail = .true.
    go to 20
!
!  Cutoff for convergence failure. NQRS will be 2*N usually.
!
170 continue

    if ( z < 0.0D+00 ) then
      q(k) = -z
      if ( wntv ) then
        v(1:nrv,k) = - v(1:nrv,k)
      end if
    end if

190 continue

  end do
!
!  Convergence. Q(K) is made nonnegative.
!
  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    if ( q(i-1) < q(i) ) then
      go to 220
    end if
  end do

  if ( fail ) then
    ipass = 2
  end if

  return
!
!  Every singular value is in order.
!
220 continue

  do i = 2, n

    t = q(i-1)
    k = i - 1

    do j = i, n
      if ( t < q(j) ) then
        t = q(j)
        k = j
      end if
    end do

    if ( k == i - 1 ) then
      cycle
    end if

    q(k) = q(i-1)
    q(i-1) = t

    if ( havers ) then
      do j = 1, ncc
        temp     = c(i-1,j)
        c(i-1,j) = c(k,j)
        c(k,j)   = temp
      end do
    end if

    if ( wntv ) then
      do j = 1, nrv
        temp     = v(j,i-1)
        v(j,i-1) = v(j,k)
        v(j,k)   = temp
      end do
    end if

  end do
!
!  End of ordering algorithm.
!
  if ( fail ) then 
    ipass = 2
  end if

  return
end
subroutine sva ( a, mda, m, n, mdata, b, sing, kpvec, names, iscale, d, work )

!*****************************************************************************80
!
!! SVA carries out a singular value analysis.
!
!  Discussion:
!
!    This routine computes the singular value decomposition of the matrix 
!    of a least squares problem, and produces a printed report.
!
!    This 1995 version differs from the original 1973 version by the
!    addition of the arguments kpvec() and work(), and by allowing user to
!    choose the length of names in names().
!
!    kpvec() allows the user to exercise options regarding printing.
!    work() provides 2*n locations of work space.  originally sing() was
!    required to have 3*n elements, of which the last 2*n were used for
!    work space.  now sing() only needs n elements.
!
!  Modified:
!
!    21 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(MDA,N).  On input, the M by N 
!    matrix of the least squares problem to be analyzed.  This could be 
!    a matrix obtained by preliminary orthogonal transformations
!    applied to the actual problem matrix which may have had more rows.  
!    See MDATA below.
!
!    Input, integer ( kind = 4 ) MDA, the first dimensioning parameter for A.
!    It must be the case that max ( M, N ) <= MDA.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in A.
!    N < M or M <= N is permitted.  0 < M and 0 < N are required.
!
!    Input, integer ( kind = 4 ) MDATA, the number of rows in the actual least squares
!    problem.  Generally M <= MDATA.  MDATA is used only in computing
!    statistics for the report and is not used as a loop count or array
!    dimension.
!
!    Input/output, real ( kind = 8 ) B(M).  On input, the right-side vector 
!    of the least squares problem.  On output, the vector G = U'*B, where U
!    comes from the singular value decomposition of A.  
!
!    Output, real ( kind = 8 ) SING(N), the singular values of A, in
!    descending order, in locations 1 through min ( M, N ).
!    If M < N, SING(M+1:N) = 0.
!
!    Input, integer ( kind = 4 ) KPVEC(4), selects print options.
!    kpvec(1) determines whether the rest of the array is to be used or ignored.
!    if kpvec(1) = 1, the contents of (kpvec(i), i=2,4)
!    will be used to set internal variables as follows:
!    prblk = kpvec(2)
!    unit  = kpvec(3)
!    width = kpvec(4)
!    if kpvec(1) = 0 default settings will be used.  the user
!    need not dimension kpvec() greater than 1.  the subr will
!    set prblk = 111111, unit = -1, and width = 69.
!
!     names()  [in]  names(j), for j = 1, ..., n, may contain a
!              name for the jth component of the solution
!              vector.  the declared length of the elements of the
!              names() array is not specifically limited, but a
!              greater length reduces the space available for columns
!              of the matris to be printed.
!              if names(1) contains only blank characters,
!              it will be assumed that no names have been provided,
!              and this subr will not access the names() array beyond
!              the first element.
!
!    Input, integer ( kind = 4 ) ISCALE, selects the column scaling option.
!    1, use the identity scaling and ignore the D array.
!    2, scale nonzero Columns to have unit euclidean length; store 
!       reciprocal lengths of original nonzero columns in D.
!    3, user supplies column scale factors in D.  This routine
!       multiplies column J by D(J) and removes the scaling
!       from the solution at the end.
!
!    Input/output, real ( kind = 8 ) D(N), a vector whose use depends
!    on the value of ISCALE.
!
!    Workspace, real ( kind = 8 ) WORK(2*N).
!
!  Local parameters:
!
!    Local PRBLK, the decimal representation of prblk must be
!    representable as at most 6 digits, each being 0 or 1.
!    the decimal digits will be interpreted as independant
!    on/off flags for the 6 possible blocks of printed output.
!    examples:  111111 selects all blocks, 0 suppresses all
!    printing,  101010 selects the 1st, 3rd, and 5th blocks, etc.
!    The six blocks are:
!    1. header, with size and scaling option parameters.
!    2. v-matrix.  amount of output depends on m and n.
!    3. singular values and related quantities.  amount of output depends on n.
!    4. listing of ynorm and rnorm and their logarithms.
!       amount of output depends on n.
!    5. levenberg-marquart analysis.
!    6. candidate solutions.  amount of output depends on m and n.
!
!    Local, integer UNIT, selects the output unit.  If 0 <= UNIT,
!    UNIT will be used as the output unit number.  If UNIT = -1, 
!    output will be written to the "*" output unit, i.e., the standard 
!    system output unit.  The calling program is responsible for opening
!    and/or closing the selected output unit if the host system requires 
!    these actions.
!
!    Local, integer WIDTH, determines the width of blocks in the output.
!    The default value is 79.  WIDTH determines the width of blocks 2, 3, 
!    and 6 of the output report.  If 95 <= WIDTH, then block 3 will use 
!    95(+1) columns; otherwise 69(+1) columns.
!    Blocks 2 and 6 are printed by subroutine MFEOUT.  These blocks generally 
!    use at most WIDTH(+1) columns, but will use more if the names are so 
!    long that more space is needed to print one name and one numeric column.  
!    The (+1)'s above are reminders that in all cases there is one extra 
!    initial column for FORTRAN "carriage control".  The carriage control 
!    character will always be a blank.
!    Blocks 1, 4, and 5 have fixed widths of 63(+1), 66(+1) and
!    66(+1), respectively.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mda,n)
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) a3
  real ( kind = 8 ) a4
  real ( kind = 8 ) alamb
  real ( kind = 8 ) aln10
  real ( kind = 8 ) b(m)
  logical blk(6)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) del
  real ( kind = 8 ) el
  real ( kind = 8 ) el2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ie
  integer ( kind = 4 ) ipass
  integer ( kind = 4 ) iscale
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kpvec(4)
  integer ( kind = 4 ) mdata
  integer ( kind = 4 ) minmn
  integer ( kind = 4 ) minmn1
  integer ( kind = 4 ) mpass
  character ( len = * ) names(n)
  logical narrow
  integer ( kind = 4 ) nsol
  real ( kind = 8 ) pcoef
  integer ( kind = 4 ) prblk
  real ( kind = 8 ) rl
  real ( kind = 8 ) rnorm
  real ( kind = 8 ) rs
  real ( kind = 8 ) sb
  real ( kind = 8 ) sing(n)
  real ( kind = 8 ) sl
  logical star
  integer ( kind = 4 ) unit
  integer ( kind = 4 ) width
  real ( kind = 8 ) work(2*n)
  real ( kind = 8 ) yl
  real ( kind = 8 ) ynorm
  real ( kind = 8 ) ys
  real ( kind = 8 ) ysq

  220 format (1x/' index  sing. val.    p coef   ', &
         '   reciprocal    g coef         g**2   ', &
         '   cumulative  scaled sqrt'/ &
     31x,'   sing. val.',26x, &
         '  sum of sqrs  of cum.s.s.')
  221 format (1x/' index  sing. val.    p coef   ', &
             '   reciprocal    g coef     scaled sqrt'/ &
         31x,'   sing. val.',13x,'  of cum.s.s.')
  222 format (1x/' index  sing. val.    g coef         g**2   ', &
             '   cumulative  scaled sqrt'/ &
         44x,'  sum of sqrs  of cum.s.s.')

  230 format (' ',4x,'0',64x,2g13.4)
  231 format (' ',4x,'0',51x, g13.4)
  232 format (' ',4x,'0',38x,2g13.4)

  240 format (' ',i5,g12.4,6g13.4)

  260 format (1x,' m  = ',i6,',   n  =',i4,',   mdata  =',i8)
  270 format (1x/' singular value analysis of the least squares', &
         ' problem,  a*x = b,'/' scaled as (a*d)*y = b.')
  280 format (1x/' scaling option number',i2,'.  d is a diagonal', &
  ' matrix with the following diagonal elements..'/(5x,10e12.4))
  290 format (1x/' scaling option number 1.   d is the identity matrix.'/1x)
  300 format (1x/' index',12x,'ynorm      rnorm',11x, &
  '      log10      log10'/ &
  45x,'      ynorm      rnorm'/1x)
  310 format (' ',i5,6x,2e11.3,11x,2f11.3)
  320 format (1x/ &
     ' norms of solution and residual vectors for a range of values'/ &
     ' of the levenberg-marquardt parameter, lambda.'// &
         '      lambda      ynorm      rnorm', &
          '      log10      log10      log10'/ &
      34x,'     lambda      ynorm      rnorm')
  330 format (1x, 3e11.3, 3f11.3)

  if ( m <= 0 .or. n <= 0 ) then
    return
  end if

  minmn = min ( m, n )
  minmn1 = minmn + 1

  if ( kpvec(1) == 0 ) then
    prblk = 111111
    unit = -1
    width = 79
  else
    prblk = kpvec(2)
    unit = kpvec(3)
    width = kpvec(4)
  end if

  star = unit < 0
!
!  Build logical array blk() by testing decimal digits of PRBLK.
!
  do i = 6, 1, -1
    j = prblk/10
    blk(i) = 0 < ( prblk - 10 * j )
    prblk = j
  end do
!
!  Optionally print the header and m, n, mdata
!
  if ( blk(1) ) then
    if ( star ) then
      write (*,270)
      write (*,260) m,n,mdata
    else
      write (unit,270)
      write (unit,260) m,n,mdata
    end if
  end if
!
!  Handle scaling as selected by ISCALE.
!
  if ( iscale == 1 ) then

    if ( blk(1) ) then
      if ( star ) then
        write (*,290)
      else
        write (unit,290)
      end if
    end if

  else
!
!  Apply column scaling to A.
!
     do j = 1, n
        a1 = d(j)
        if ( iscale <= 2 ) then
           sb = 0.0D+00
           do i = 1, m
             sb = sb + a(i,j)**2
           end do
           a1 = sqrt ( sb )
           if ( a1 == 0.0D+00 ) then
             a1 = 1.0D+00
           end if
           a1 = 1.0D+00 / a1
           d(j) = a1
        end if
        do i = 1, m
           a(i,j) = a(i,j)*a1
        end do
     end do

     if ( blk(1) ) then
       if ( star ) then
         write (*,280) iscale,(d(j),j = 1, n)
       else
         write (unit,280) iscale,(d(j),j = 1, n)
       end if
     end if

  end if
!
!  Compute the singular value decomposition of the scaled matrix.
!
  call svdrs ( a, mda, m, n, b, m, 1, sing )
!
!  Determine nsol.
!
  nsol = minmn

  do j = 1, minmn
    if ( sing(j) == 0.0D+00 ) then
      nsol = j-1
      exit
    end if
  end do
!
!  The array b() contains the vector g.
!  Compute cumulative sums of squares of components of
!  g and store them in work(i), i = 1,...,minmn+1
!
  sb = sum ( b(minmn1:m)**2 )

  work(minmn+1) = sb
  do j = minmn, 1, -1
    sb = sb + b(j)**2
    work(j) = sb
  end do
!
!  Print the V matrix.
!
  if ( blk(2) ) then
    call mfeout ( a, mda, n, n, names, 1, unit, width )
  end if
!
!  Replace v by d*v in the array A.
!
  if ( 1 < iscale ) then
    do i = 1, n
      a(i,1:n) = d(i) * a(i,1:n)
    end do
  end if

  if ( blk(3) ) then
!
!  Print singular values and other summary results.
!
!  Output will be done using one of two layouts.  the narrow
!  layout uses 69 columns + 1 for carriage control, and makes two passes
!  through the computation.
!
!  The wide layout uses 95 columns + 1 for carriage control, and makes
!  only one pass through the computation.
!
!  G is now in  b() array.  v now in a(,) array.
!
  narrow = width < 95

  if ( narrow ) then
    mpass = 2
  else
    mpass = 1
  end if

  do 170 ipass = 1, mpass
     if ( star ) then
        if ( narrow ) then
           if ( ipass == 1 ) then
             write(*,221)
           else
             write(*,222)
           end if
        else
           write (*,220)
        end if
     else
        if ( narrow ) then
           if ( ipass == 1 ) then
              write(unit,221)
           else
              write(unit,222)
           end if
        else
           write (unit,220)
        end if
     end if
!
!  The following statement converts from integer to floating-point.
!
  a3 = work(1)
  a4 = sqrt ( a3 / max ( 1, mdata ) )

  if ( star ) then
    if ( narrow ) then
      if ( ipass == 1 ) then
        write(*,231) a4
      else
        write(*,232) a3, a4
      end if
    else
      write (*,230) a3,a4
    end if
  else
    if ( narrow ) then
        if ( ipass == 1 ) then
           write(unit,231) a4
        else
           write(unit,232) a3, a4
        end if
     else
        write (unit,230) a3,a4
     end if
  end if

  do k = 1, minmn
     if ( sing(k) == 0.0D+00 ) then
        pcoef = 0.0D+00
        if ( star ) then
          write (*,240) k,sing(k)
        else
          write (unit,240) k,sing(k)
        end if
     else
        pcoef = b(k) / sing(k)
        a1 = 1.0D+00 / sing(k)
        a2 = b(k)**2
        a3 = work(k+1)
        a4 = sqrt ( a3 / max ( 1, mdata-k ) ) 
        if ( star ) then
           if ( narrow ) then
              if ( ipass == 1 ) then
                 write(*,240) k,sing(k),pcoef,a1,b(k),      a4
              else
                 write(*,240) k,sing(k),         b(k),a2,a3,a4
              end if
           else
              write (*,240) k,sing(k),pcoef,a1,b(k),a2,a3,a4
           end if
        else
           if ( narrow ) then
              if ( ipass == 1 ) then
                 write(unit,240) k,sing(k),pcoef,a1,b(k),      a4
              else
                 write(unit,240) k,sing(k),         b(k),a2,a3,a4
              end if
           else
              write (unit,240) k,sing(k),pcoef,a1,b(k),a2,a3,a4
           end if

        end if
     end if

  end do

  170 continue

  end if

  if ( blk(4) ) then
!
!  Compute and print values of ynorm, rnorm and their logarithms.
!
  if ( star ) then
    write (*,300)
  else
    write (unit,300)
  end if

  ysq = 0.0D+00

  do j = 0, nsol

     if ( j /= 0 ) then
       ysq = ysq + (b(j) / sing(j))**2
     end if

     ynorm = sqrt(ysq)
     rnorm = sqrt(work(j+1))

     yl = -1000.0D+00
     if ( 0.0D+00 < ynorm ) then
       yl = log10 ( ynorm )
     end if

     rl = -1000.0D+00
     if ( 0.0D+00 < rnorm ) then
       rl = log10 ( rnorm )
     end if

     if ( star ) then
       write (*,310) j,ynorm,rnorm,yl,rl
     else
       write (unit,310) j,ynorm,rnorm,yl,rl
     end if

  end do

  end if

  if ( blk(5) .and. sing(1) /= 0.0D+00 ) then
!
!  Compute values of xnorm and rnorm for a sequence of values of
!  the levenberg-marquardt parameter.
!
     el = log10 ( sing(1) ) + 1.0D+00
     el2 = log10 ( sing(nsol) ) - 1.0D+00
     del = ( el2 - el ) / 20.0D+00
     aln10 = log ( 10.0D+00 )

     if ( star ) then
       write (*,320)
     else
       write (unit,320)
     end if

     do ie = 1, 21
!
!  Compute alamb = 10.0D+00 **el
!
        alamb = exp(aln10*el)
        ys = 0.0D+00
        rs = work(nsol+1)

        do i = 1, minmn
           sl = sing(i)**2 + alamb**2
           ys = ys + (b(i)*sing(i)/sl)**2
           rs = rs + (b(i)*(alamb**2)/sl)**2
        end do

        ynorm = sqrt(ys)
        rnorm = sqrt(rs)
        rl = -1000.0D+00
        if ( 0.0D+00 < rnorm ) then
          rl = log10 ( rnorm )
        end if

        yl = -1000.0D+00
        if ( 0.0D+00 < ynorm ) then
          yl = log10 ( ynorm )
        end if

        if ( star ) then
          write (*,330) alamb,ynorm,rnorm,el,yl,rl
        else
          write (unit,330) alamb,ynorm,rnorm,el,yl,rl
        end if

        el = el + del

    end do

  end if
!
!  Compute and optionally print candidate solutions.
!
  do k = 1, nsol
    pcoef = b(k) / sing(k)
    do i = 1, n
      a(i,k) = a(i,k) * pcoef
      if ( 1 < k ) then
        a(i,k) = a(i,k) + a(i,k-1)
      end if
    end do
  end do

  if ( blk(6) .and. 1 <= nsol ) then
    call mfeout ( a, mda, n, nsol, names, 2, unit, width )
  end if

  return
end
subroutine svdrs ( a, mda, m, n, b, mdb, nb, s )

!*****************************************************************************80
!
!! SVDRS: singular value decomposition with a right side vector.
!
!  Discussion:
!
!    This subroutine computes the singular value decomposition of the
!    M by N matrix A and optionally applies the transformations
!    from the left to the NB column vectors of the M by NB matrix B.
!
!    Either N <= M or M < N is permitted.
!
!    The singular value decomposition of A is of the form
!
!      A  =  U * S * V'
!
!    where U is M by M orthogonal, S is M by N diagonal with the
!    diagonal terms nonnegative and ordered from large to small, and
!    V is N by N orthogonal.  Note that these matrices also satisfy
!
!      S = U' * A * V
!
!    The matrix V is returned in the leading N rows and
!    columns of the array A.
!
!    The singular values, i.e. the diagonal terms of the matrix S,
!    are returned in the array S.  If M < N, positions M+1
!    through N of S will be set to zero.
!
!    The product matrix G = U' * B overwrites the input matrix B.
!
!    If the user wishes to obtain a minimum length least squares
!    solution of the linear system
!
!      A * x ~ = ~ b
!
!    the solution X can be constructed, following use of this subroutine,
!    by computing the sum for I = 1, ..., R of the outer products
!
!      (column I of V) * 1/S(I) * (row I of G)
!
!    Here R denotes the pseudorank of A which the user may choose
!    in the range 0 through min ( M, N ) based on the sizes of the
!    singular values.
!
!    This code gives special treatment to rows and columns that are
!    entirely zero.  This causes certain zero singular values to appear as
!    exact zeros rather than as about macheps times the largest singular value.
!    It similarly cleans up the associated columns of U and V.
!
!  Algorithm:
!
!    1. Exchange columns of A to pack nonzero columns to the left.
!       Set N_COPY = number of nonzero columns.
!       Use locations A(1,N), A(1,N-1),..., A(1,N_COPY+1) to record the
!       column permutations.
!
!    2. Exchange rows of A to pack nonzero rows to the top.
!       Quit packing if find N_COPY nonzero rows.  Make same row exchanges
!       in B.  Set M_COPY so that all nonzero rows of the permuted A
!       are in first M_COPY rows.  If M_COPY <= N_COPY then all M_COPY rows are
!       nonzero.  If N_COPY < M_COPY then the first N_COPY rows are known
!       to be nonzero, and rows NCOPY+1 through M_COPY may be zero or nonzero.
!
!    3. Apply the original algorithm to the M_COPY by N_COPY problem.
!
!    4. Move permutation record from A(,) to S(I), I = N_COPY+1,...,N.
!
!    5. Build V up from N_COPY by N_COPY to N by N by placing ones on
!       the diagonal and zeros elsewhere.  This is only partly done
!       explicitly.  It is completed during step 6.
!
!    6. Exchange rows of V to compensate for column exchanges of step 2.
!
!    7. Place zeros in S(N_COPY+1:N)  to represent zero singular values.
!
!  Modified:
!
!    21 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson.
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(MDA,N).
!    On input, the M by N matrix A.
!    On output contains the N by N matrix V.
!
!    Input, integer ( kind = 4 ) MDA, the leading dimension of A.
!    MDA must be at least the maximum of M and N.
!
!    Input, integer ( kind = 4 ) M, the number of rows of A and B.
!    M must be at least 1.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A,
!    the number of rows and columns of the matrix V.
!    N must be greater than 0.
!
!    Input/output, real ( kind = 8 ) B(MDB,NB).  If 0 < NB, this array 
!    must contain an M by NB matrix on input and will contain the
!    M by NB product matrix, G = U' * B on output.
!
!    Input, integer ( kind = 4 ) LDB, the leading dimension of B, which must be
!    at least M.
!
!    Input, integer ( kind = 4 ) NB, the number of right hand sides; that is, the
!    number of columns of data in B.  If no right hand sides are being
!    supplied, set NB = 0.
!
!    Output, real S(N), the singular values of A, with the ordering
!    0 <= S(N) <= S(N-1) <= ... <= S(2) <= S(1).
!    If M < N the singular values indexed from M+1 through N will be zero.
!
!  Local parameters:
!
!    real ( kind = 8 ) WORK(N,2).
!    locations 1 thru N will hold the off-diagonal terms of
!    the bidiagonal matrix for subroutine QRBD.  Locations N+1
!    thru 2*N will save info from one call to the next of H12.
!
  implicit none

  integer ( kind = 4 ) mda
  integer ( kind = 4 ) mdb
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb

  real ( kind = 8 ) a(mda,n)
  real ( kind = 8 ) b(mdb,nb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipass
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_copy
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) n_copy
  integer ( kind = 4 ) ns
  integer ( kind = 4 ) nsp1
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) t
  real ( kind = 8 ) work(n,2)

  n_copy = n

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVDRS - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of M = ', m
    return
  end if

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVDRS - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of N = ', n
    return
  end if

  if ( mda < m .or. mda < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVDRS - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of MDA = ', mda
    return
  end if

  if ( mdb < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVDRS - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of MDB = ', mdb
    return
  end if
!
!  If column J is entirely zero, exchange it with column N.
!
  j = n_copy

  do while ( 1 <= j )

    if ( all ( a(1:m,j) == 0.0D+00 ) ) then

      if ( j /= n_copy ) then
        a(1:m,j) = a(1:m,n_copy)
      end if

      a(1,n_copy) = j
      n_copy = n_copy - 1

    end if

    j = j - 1

  end do

  ns = 0
!
!  If N_COPY = 0 then A is entirely zero and the SVD computation can be skipped.
!
  if ( n_copy /= 0 ) then
!
!  Pack nonzero rows to the top.
!  Quit packing if find N_COPY nonzero rows.
!
    i = 1
    m_copy = m

    do while ( i <= n_copy .and. i < m_copy ) 

      if ( a(i,i) == 0.0D+00 ) then

        if ( all ( a(i,1:n_copy) == 0.0D+00 ) ) then

          do j = 1, nb
            t           = b(i,j)
            b(i,j)      = b(m_copy,j)
            b(m_copy,j) = t
          end do

          a(i,1:n_copy) = a(m_copy,1:n_copy)

          if ( m_copy <= n_copy ) then
            a(m_copy,1:n_copy) = 0.0D+00
          end if

           m_copy = m_copy - 1

        else
          i = i + 1
        end if

      else
        i = i + 1
      end if      

    end do
!
!  Begin the SVD algorithm:
!
!  1: Reduce the matrix to upper bidiagonal form with Householder
!     transformations
!       H(N_COPY)...H(1) * A * Q(1)...Q(N_COPY-2) = (D',0)'
!     where D is upper bidiagonal.
!
!  2: Apply H(N_COPY)...H(1) to B.  Here H(N_COPY)...H(1)*B replaces
!     B in storage.
!
!  3: The matrix product W = Q(1)...Q(N_COPY-2) overwrites the first
!     N_COPY rows of A in storage.
!
!  4: An SVD for D is computed.  Here K rotations Ri and Pi are
!     computed so that
!       Rk...R1*D*P1'...Pk' = diag(S1,...,Sm_copy)
!     to working accuracy.  The S(I) are nonnegative and nonincreasing.
!     Here Rk...R1*B overwrites B in storage while
!     A*P1'...Pk' overwrites A in storage.
!
!  5: It follows that, with the proper definitions, U'*B overwrites B, 
!     while V overwrites the first N_COPY row and columns of A.
!
    l = min ( m_copy, n_copy )
!
!  The following loop reduces A to upper bidiagonal and
!  also applies the premultiplying transformations to B.
!
    do j = 1, l
!
!  Note:  In the original code, there are instances where the calls
!  to H12 may involve formally incorrect memory accesses.
!
!  For instance, the case J = N is possible, and in that case,
!  the reference A(1,J+1) is illegal.  However, while it seems as though
!  an illegal memory access is possible, no illegal memory access is
!  committed.  In this case, for instance, the fact that N_COPY - J 
!  will be 0 means that there will be no operations carried out on 
!  the memory that is "illegally" being pointed to.  Some attempt has
!  been made to clean up the code, to keep overcautious compilers and
!  debuggers from flagging these occurrences.
!
      if ( j < m_copy ) then

        mode = 1

        if ( j < n ) then

          call h12 ( mode, j, j+1, m_copy, a(1,j), 1, t, a(1,j+1), 1, mda, &
            n_copy-j )
!
!  If N <= J, replace pointer to A(1,J+1) by pointer to A(1,J).  It
!  doesn't matter because the memory will not be used in this case.
!
        else

          call h12 ( mode, j, j+1, m_copy, a(1,j), 1, t, a(1,j  ), 1, mda, &
            n_copy-j )

        end if

        mode = 2
        call h12 ( mode, j, j+1, m_copy, a(1,j), 1, t, b, 1, mdb, nb )

      end if

      if ( j < n_copy-1 ) then

        mode = 1
        call h12 ( mode, j+1, j+2, n_copy, a(j,1), mda, work(j,2), a(j+1,1), &
          mda, 1, m_copy-j )

      end if

    end do
!
!  Copy the bidiagonal matrix into S and WORK for QRBD.
!
    do j = 2, l
      s(j) = a(j,j)
      work(j,1) = a(j-1,j)
    end do
    s(1) = a(1,1)

    ns = n_copy

    if ( m_copy < n_copy ) then
      ns = m_copy + 1
      s(ns) = 0.0D+00
      work(ns,1) = a(m_copy,m_copy+1)
    end if
!
!  Construct the explicit N_COPY by N_COPY product matrix, 
!  W = Q1*Q2*...*QL*I in A.
!
    do i = n_copy, 1, -1

      if ( i <= min ( m_copy, n_copy-2 ) ) then
        mode = 2
        call h12 ( mode, i+1, i+2, n_copy, a(i,1), mda, work(i,2), a(1,i+1), &
          1, mda, n_copy-i )
      end if

      a(i,1:n_copy) = 0.0D+00
      a(i,i) = 1.0D+00

    end do
!
!  Compute the singular value decomposition of the bidiagonal matrix.
!
    call qrbd ( ipass, s(1), work(1,1), ns, a, mda, n_copy, b, mdb, nb )

    if ( ipass == 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SVDRS - Warning:'
      write ( *, '(a)' ) '  Full accuracy was not attained in the singular'
      write ( *, '(a)' ) '  value decomposition of the bidiagonal matrix.'
    end if

  end if

  s(ns+1:n_copy) = 0.0D+00
!
!  Extract the record of permutations and store zeros.
!
  s(n_copy+1:n) = a(1,n_copy+1:n)
  a(1:n_copy,n_copy+1:n) = 0.0D+00
!
!  Permute rows and set zero singular values.
!
  do k = n_copy+1, n

    i = int ( s(k) )

    s(k) = 0.0D+00

    a(k,1:n) = a(i,1:n)

    a(i,1:n) = 0.0D+00
    a(i,k) = 1.0D+00

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
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
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

  character ( len = 8  ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
