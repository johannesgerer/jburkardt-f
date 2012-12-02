subroutine daxpy ( n, da, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DAXPY computes constant times a vector plus a vector.
!
!  Discussion:
!
!    Uses unrolled loops for increments equal to one.
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in DX and DY.
!
!    Input, real ( kind = 8 ) DA, the multiplier of DX.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of DX.
!
!    Input/output, real ( kind = 8 ) DY(*), the second vector.
!    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive 
!    entries of DY.
!
  implicit none

  real ( kind = 8 ) da
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  if ( n <= 0 ) then
    return
  end if

  if ( da  == 0.0D+00 ) then
    return
  end if
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
  if ( incx /= 1 .or. incy /= 1 ) then

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      dy(iy) = dy(iy) + da * dx(ix)
      ix = ix + incx
      iy = iy + incy
    end do
!
!  Code for both increments equal to 1.
!
  else

    m = mod ( n, 4 )

    do i = 1, m
      dy(i) = dy(i) + da * dx(i)
    end do

    do i = m+1, n, 4
      dy(i  ) = dy(i  ) + da * dx(i  )
      dy(i+1) = dy(i+1) + da * dx(i+1)
      dy(i+2) = dy(i+2) + da * dx(i+2)
      dy(i+3) = dy(i+3) + da * dx(i+3)
    end do

  end if

  return
end
function ddot ( n, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DDOT forms the dot product of two vectors.
!
!  Discussion:
!
!    This routine uses unrolled loops for increments equal to one.
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive entries 
!    in X.
!
!    Input, real ( kind = 8 ) DY(*), the second vector.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive entries 
!    in Y.
!
!    Output, real DDOT, the sum of the product of the corresponding
!    entries of X and Y.
!
  implicit none

  real ( kind = 8 ) ddot
  real ( kind = 8 ) dtemp
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  ddot = 0.0D+00
  dtemp = 0.0D+00

  if ( n <= 0 ) then
    return
  end if
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
  if ( incx /= 1 .or. incy /= 1 ) then

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      dtemp = dtemp + dx(ix) * dy(iy)
      ix = ix + incx
      iy = iy + incy
    end do
!
!  Code for both increments equal to 1.
!
  else

    m = mod ( n, 5 )

    do i = 1, m
      dtemp = dtemp + dx(i) * dy(i)
    end do

    do i = m+1, n, 5

      dtemp = dtemp + dx(i  ) * dy(i  ) &
                    + dx(i+1) * dy(i+1) &
                    + dx(i+2) * dy(i+2) &
                    + dx(i+3) * dy(i+3) &
                    + dx(i+4) * dy(i+4)
    end do

  end if

  ddot = dtemp

  return
end
function dnrm2 ( n, x, incx )

!*****************************************************************************80
!
!! DNRM2 returns the euclidean norm of a vector.
!
!  Discussion:
!
!     DNRM2 ( X ) = sqrt ( X' * X )
!
!  Author:
!
!    Sven Hammarling
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(*), the vector whose norm is to be computed.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive entries 
!    of X.
!
!    Output, real ( kind = 8 ) DNRM2, the Euclidean norm of X.
!
  implicit none

  real ( kind = 8 ) absxi
  real ( kind = 8 ) dnrm2
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n
  real ( kind = 8 ) norm
  real ( kind = 8 ) scale
  real ( kind = 8 ) ssq
  real ( kind = 8 ) x(*)

  if ( n < 1 .or. incx < 1 ) then

    norm  = 0.0D+00

  else if ( n == 1 ) then

    norm  = abs ( x(1) )

  else

    scale = 0.0D+00
    ssq = 1.0D+00

    do ix = 1, 1 + ( n - 1 )*incx, incx
      if ( x(ix) /= 0.0D+00 ) then
        absxi = abs ( x(ix) )
        if ( scale < absxi ) then
          ssq = 1.0D+00 + ssq * ( scale / absxi )**2
          scale = absxi
        else
          ssq = ssq + ( absxi / scale )**2
        end if
      end if
    end do
    norm  = scale * sqrt( ssq )
  end if

  dnrm2 = norm

  return
end
subroutine drot ( n, x, incx, y, incy, c, s )

!*****************************************************************************80
!
!! DROT applies a plane rotation.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input/output, real ( kind = 8 ) X(*), one of the vectors to be rotated.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of X.
!
!    Input/output, real ( kind = 8 ) Y(*), one of the vectors to be rotated.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive
!    elements of Y.
!
!    Input, real ( kind = 8 ) C, S, parameters (presumably the cosine and
!    sine of some angle) that define a plane rotation.
!
  implicit none

  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) n
  real ( kind = 8 ) s
  real ( kind = 8 ) stemp
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  if ( n <= 0 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    do i = 1, n
      stemp = c * x(i) + s * y(i)
      y(i) = c * y(i) - s * x(i)
      x(i) = stemp
    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      stemp = c * x(ix) + s * y(iy)
      y(iy) = c * y(iy) - s * x(ix)
      x(ix) = stemp
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine drotg ( sa, sb, c, s )

!*****************************************************************************80
!
!! DROTG constructs a Givens plane rotation.
!
!  Discussion:
!
!    This routine uses double precision real arithmetic.
!
!    Given values A and B, this routine computes
!
!    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
!          = sign ( B ) if abs ( A ) <= abs ( B );
!
!    R     = SIGMA * ( A * A + B * B );
!
!    C = A / R if R is not 0
!      = 1     if R is 0;
!
!    S = B / R if R is not 0,
!        0     if R is 0.
!
!    The computed numbers then satisfy the equation
!
!    (  C  S ) ( A ) = ( R )
!    ( -S  C ) ( B ) = ( 0 )
!
!    The routine also computes
!
!    Z = S     if abs ( A ) > abs ( B ),
!      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
!      = 1     if C is 0.
!
!    The single value Z encodes C and S, and hence the rotation:
!
!    If Z = 1, set C = 0 and S = 1;
!    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
!    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Algorithm 539, 
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) SA, SB.  On input, SA and SB are the values
!    A and B.  On output, SA is overwritten with R, and SB is
!    overwritten with Z.
!
!    Output, real ( kind = 8 ) C, S, the cosine and sine of the
!    Givens rotation.
!
  implicit none

  real ( kind = 8 ) c
  real ( kind = 8 ) r
  real ( kind = 8 ) roe
  real ( kind = 8 ) s
  real ( kind = 8 ) sa
  real ( kind = 8 ) sb
  real ( kind = 8 ) scale
  real ( kind = 8 ) z

  if ( abs ( sb ) < abs ( sa ) ) then
    roe = sa
  else
    roe = sb
  end if

  scale = abs ( sa ) + abs ( sb )

  if ( scale == 0.0D+00 ) then
    c = 1.0D+00
    s = 0.0D+00
    r = 0.0D+00
  else
    r = scale * sqrt ( ( sa / scale )**2 + ( sb / scale )**2 )
    r = sign ( 1.0D+00, roe ) * r
    c = sa / r
    s = sb / r
  end if

  if ( 0.0D+00 < abs ( c ) .and. abs ( c ) <= s ) then
    z = 1.0D+00 / c
  else
    z = s
  end if

  sa = r
  sb = z

  return
end
subroutine dscal ( n, sa, x, incx )

!*****************************************************************************80
!
!! DSCAL scales a vector by a constant.
!
!  Modified:
!
!    08 April 1999
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) SA, the multiplier.
!
!    Input/output, real ( kind = 8 ) X(*), the vector to be scaled.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) sa
  real ( kind = 8 ) x(*)

  if ( n <= 0 ) then

  else if ( incx == 1 ) then

    m = mod ( n, 5 )

    x(1:m) = sa * x(1:m)

    do i = m+1, n, 5
      x(i)   = sa * x(i)
      x(i+1) = sa * x(i+1)
      x(i+2) = sa * x(i+2)
      x(i+3) = sa * x(i+3)
      x(i+4) = sa * x(i+4)
    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    do i = 1, n
      x(ix) = sa * x(ix)
      ix = ix + incx
    end do

  end if

  return
end
subroutine dsvdc ( a, lda, m, n, s, e, u, ldu, v, ldv, work, job, info )

!*****************************************************************************80
!
!! DSVDC computes the singular value decomposition of a real rectangular matrix.
!
!  Discussion:
!
!    This routine reduces an M by N matrix A to diagonal form by orthogonal
!    transformations U and V.  The diagonal elements S(I) are the singular
!    values of A.  The columns of U are the corresponding left singular
!    vectors, and the columns of V the right singular vectors.
!
!    The form of the singular value decomposition is then
!
!      A(MxN) = U(MxM) * S(MxN) * V(NxN)'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(LDA,N).  On input, the M by N
!    matrix whose singular value decomposition is to be computed.
!    On output, the matrix has been destroyed.  Depending on the user's
!    requests, the matrix may contain other useful information.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!    LDA must be at least N.
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A.
!
!    Output, real ( kind = 8 ) S(MM), where MM = max(M+1,N).  The first
!    min(M,N) entries of S contain the singular values of A arranged in
!    descending order of magnitude.
!
!    Output, real ( kind = 8 ) E(MM), where MM = max(M+1,N).  Ordinarily
!    contains zeros.  However see the discussion of INFO for exceptions.
!
!    Output, real ( kind = 8 ) U(LDU,K).  If JOBA = 1 then K = M;
!    if 2 <= JOBA, then K = min(M,N).  U contains the M by M matrix of
!    left singular vectors.  U is not referenced if JOBA = 0.  If M <= N
!    or if JOBA = 2, then U may be identified with A in the subroutine call.
!
!    Input, integer ( kind = 4 ) LDU, the leading dimension of the array U.
!    LDU must be at least M.
!
!    Output, real ( kind = 8 ) V(LDV,N), the N by N matrix of right singular
!    vectors.  V is not referenced if JOB is 0.  If N <= M, then V may be
!    identified with A in the subroutine call.
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of the array V.
!    LDV must be at least N.
!
!    Workspace, real ( kind = 8 ) WORK(M).
!
!    Input, integer ( kind = 4 ) JOB, controls the computation of the singular
!    vectors.  It has the decimal expansion AB with the following meaning:
!      A =  0, do not compute the left singular vectors.
!      A =  1, return the M left singular vectors in U.
!      A >= 2, return the first min(M,N) singular vectors in U.
!      B =  0, do not compute the right singular vectors.
!      B =  1, return the right singular vectors in V.
!
!    Output, integer ( kind = 4 ) INFO, status indicator.
!    The singular values (and their corresponding singular vectors)
!    S(INFO+1), S(INFO+2),...,S(MN) are correct.  Here MN = min ( M, N ).
!    Thus if INFO is 0, all the singular values and their vectors are
!    correct.  In any event, the matrix B = U' * A * V is the bidiagonal
!    matrix with the elements of S on its diagonal and the elements of E on
!    its superdiagonal.  Thus the singular values of A and B are the same.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cs
  real ( kind = 8 ) e(*)
  real ( kind = 8 ) el
  real ( kind = 8 ) emm1
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) jobu
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kase
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) lls
  integer ( kind = 4 ) ls
  integer ( kind = 4 ) lu
  integer ( kind = 4 ), parameter :: maxit = 30
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) mm1
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) nct
  integer ( kind = 4 ) nctp1
  integer ( kind = 4 ) ncu
  integer ( kind = 4 ) nrt
  integer ( kind = 4 ) nrtp1
  real ( kind = 8 ) s(*)
  real ( kind = 8 ) scale
  real ( kind = 8 ) ddot
  real ( kind = 8 ) shift
  real ( kind = 8 ) sl
  real ( kind = 8 ) sm
  real ( kind = 8 ) smm1
  real ( kind = 8 ) sn
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) test
  real ( kind = 8 ) u(ldu,m)
  real ( kind = 8 ) v(ldv,n)
  logical wantu
  logical wantv
  real ( kind = 8 ) work(m)
  real ( kind = 8 ) ztest
!
!  Determine what is to be computed.
!
  wantu = .false.
  wantv = .false.
  jobu = mod ( job, 100 ) / 10

  if ( 1 < jobu ) then
    ncu = min ( m, n )
  else
    ncu = m
  end if

  if ( jobu /= 0 ) then
    wantu = .true.
  end if

  if ( mod ( job, 10 ) /= 0 ) then
    wantv = .true.
  end if
!
!  Reduce A to bidiagonal form, storing the diagonal elements
!  in S and the super-diagonal elements in E.
!
  info = 0
  nct = min ( m-1, n )
  nrt = max ( 0, min ( m, n-2 ) )
  lu = max ( nct, nrt )

  do l = 1, lu
!
!  Compute the transformation for the L-th column and
!  place the L-th diagonal in S(L).
!
    if ( l <= nct ) then

      s(l) = dnrm2 ( m-l+1, a(l,l), 1 )

      if ( s(l) /= 0.0D+00 ) then
        if ( a(l,l) /= 0.0D+00 ) then
          s(l) = sign ( s(l), a(l,l) )
        end if
        call dscal ( m-l+1, 1.0D+00 / s(l), a(l,l), 1 )
        a(l,l) = 1.0D+00 + a(l,l)
      end if

      s(l) = -s(l)

    end if

    do j = l+1, n
!
!  Apply the transformation.
!
      if ( l <= nct .and. s(l) /= 0.0D+00 ) then
        t = -ddot ( m-l+1, a(l,l), 1, a(l,j), 1 ) / a(l,l)
        call daxpy ( m-l+1, t, a(l,l), 1, a(l,j), 1 )
      end if
!
!  Place the L-th row of A into E for the
!  subsequent calculation of the row transformation.
!
      e(j) = a(l,j)

    end do
!
!  Place the transformation in U for subsequent back multiplication.
!
    if ( wantu .and. l <= nct ) then
      u(l:m,l) = a(l:m,l)
    end if
!
!  Compute the L-th row transformation and place the
!  L-th superdiagonal in E(L).
!
    if ( l <= nrt ) then

      e(l) = dnrm2 ( n-l, e(l+1), 1 )

      if ( e(l) /= 0.0D+00 ) then
        if ( e(l+1) /= 0.0D+00 ) then
          e(l) = sign ( e(l), e(l+1) )
        end if
        call dscal ( n-l, 1.0D+00 / e(l), e(l+1), 1 )
        e(l+1) = 1.0D+00 + e(l+1)
      end if

      e(l) = -e(l)
!
!  Apply the transformation.
!
      if ( l + 1 <= m .and. e(l) /= 0.0D+00 ) then

        work(l+1:m) = 0.0D+00

        do j = l+1, n
          call daxpy ( m-l, e(j), a(l+1,j), 1, work(l+1), 1 )
        end do

        do j = l+1, n
          call daxpy ( m-l, -e(j)/e(l+1), work(l+1), 1, a(l+1,j), 1 )
        end do

      end if
!
!  Place the transformation in V for subsequent back multiplication.
!
      if ( wantv ) then
        v(l+1:n,l) = e(l+1:n)
      end if

    end if

  end do
!
!  Set up the final bidiagonal matrix of order MN.
!
  mn = min ( m + 1, n )
  nctp1 = nct + 1
  nrtp1 = nrt + 1

  if ( nct < n ) then
    s(nctp1) = a(nctp1,nctp1)
  end if

  if ( m < mn ) then
    s(mn) = 0.0D+00
  end if

  if ( nrtp1 < mn ) then
    e(nrtp1) = a(nrtp1,mn)
  end if

  e(mn) = 0.0D+00
!
!  If required, generate U.
!
  if ( wantu ) then

    u(1:m,nctp1:ncu) = 0.0D+00

    do j = nctp1, ncu
      u(j,j) = 1.0D+00
    end do

    do ll = 1, nct

      l = nct - ll + 1

      if ( s(l) /= 0.0D+00 ) then

        do j = l+1, ncu
          t = -ddot ( m-l+1, u(l,l), 1, u(l,j), 1 ) / u(l,l)
          call daxpy ( m-l+1, t, u(l,l), 1, u(l,j), 1 )
        end do

        u(l:m,l) = -u(l:m,l)
        u(l,l) = 1.0D+00 + u(l,l)
        u(1:l-1,l) = 0.0D+00

      else

        u(1:m,l) = 0.0D+00
        u(l,l) = 1.0D+00

      end if

    end do

  end if
!
!  If it is required, generate V.
!
  if ( wantv ) then

    do ll = 1, n

      l = n - ll + 1

      if ( l <= nrt .and. e(l) /= 0.0D+00 ) then

        do j = l + 1, n
          t = -ddot ( n-l, v(l+1,l), 1, v(l+1,j), 1 ) / v(l+1,l)
          call daxpy ( n-l, t, v(l+1,l), 1, v(l+1,j), 1 )
        end do

      end if

      v(1:n,l) = 0.0D+00
      v(l,l) = 1.0D+00

    end do

  end if
!
!  Main iteration loop for the singular values.
!
  mm = mn
  iter = 0

  do while ( 0 < mn )
!
!  If too many iterations have been performed, set flag and return.
!
    if ( maxit <= iter ) then
      info = mn
      return
    end if
!
!  This section of the program inspects for
!  negligible elements in the S and E arrays.
!
!  On completion the variables KASE and L are set as follows:
!
!  KASE = 1     if S(MN) and E(L-1) are negligible and L < MN
!  KASE = 2     if S(L) is negligible and L < MN
!  KASE = 3     if E(L-1) is negligible, L < MN, and
!               S(L), ..., S(MN) are not negligible (QR step).
!  KASE = 4     if E(MN-1) is negligible (convergence).
!
    do ll = 1, mn

      l = mn - ll

      if ( l == 0 ) then
        exit
      end if

      test = abs ( s(l) ) + abs ( s(l+1) )
      ztest = test + abs ( e(l) )

      if ( ztest == test ) then
        e(l) = 0.0D+00
        exit
      end if

    end do

    if ( l == mn - 1 ) then

      kase = 4

    else

      do lls = l + 1, mn + 1

        ls = mn - lls + l + 1

        if ( ls == l ) then
          exit
        end if

        test = 0.0D+00
        if ( ls /= mn ) then
          test = test + abs ( e(ls) )
        end if

        if ( ls /= l + 1 ) then
          test = test + abs ( e(ls-1) )
        end if

        ztest = test + abs ( s(ls) )

        if ( ztest == test ) then
          s(ls) = 0.0D+00
          exit
        end if

      end do

      if ( ls == l ) then
        kase = 3
      else if ( ls == mn ) then
        kase = 1
      else
        kase = 2
        l = ls
      end if

    end if

    l = l + 1
!
!  Deflate negligible S(MN).
!
    if ( kase == 1 ) then

      mm1 = mn - 1
      f = e(mn-1)
      e(mn-1) = 0.0D+00

      do kk = l, mm1

        k = mm1 - kk + l
        t1 = s(k)
        call drotg ( t1, f, cs, sn )
        s(k) = t1

        if ( k /= l ) then
          f = -sn * e(k-1)
          e(k-1) = cs * e(k-1)
        end if

        if ( wantv ) then
          call drot ( n, v(1,k), 1, v(1,mn), 1, cs, sn )
        end if

      end do
!
!  Split at negligible S(L).
!
    else if ( kase == 2 ) then

      f = e(l-1)
      e(l-1) = 0.0D+00

      do k = l, mn

        t1 = s(k)
        call drotg ( t1, f, cs, sn )
        s(k) = t1
        f = -sn * e(k)
        e(k) = cs * e(k)
        if ( wantu ) then
          call drot ( m, u(1,k), 1, u(1,l-1), 1, cs, sn )
        end if

      end do
!
!  Perform one QR step.
!
    else if ( kase == 3 ) then
!
!  Calculate the shift.
!
      scale = max ( abs ( s(mn) ), abs ( s(mn-1) ), abs ( e(mn-1) ), &
                    abs ( s(l) ), abs ( e(l) ) )

      sm = s(mn) / scale
      smm1 = s(mn-1) / scale
      emm1 = e(mn-1) / scale
      sl = s(l) / scale
      el = e(l) / scale
      b = ( ( smm1 + sm ) * ( smm1 - sm ) + emm1 * emm1 ) / 2.0D+00
      c = sm  * sm * emm1 * emm1
      shift = 0.0D+00

      if ( b /= 0.0D+00 .or. c /= 0.0D+00 ) then
        shift = sqrt ( b * b + c )
        if ( b < 0.0D+00 ) then
          shift = -shift
        end if
        shift = c / ( b + shift )
      end if

      f = ( sl + sm ) * ( sl - sm ) + shift
      g = sl * el
!
!  Chase zeros.
!
      mm1 = mn - 1

      do k = l, mm1

        call drotg ( f, g, cs, sn )

        if ( k /= l ) then
          e(k-1) = f
        end if

        f = cs * s(k) + sn * e(k)
        e(k) = cs * e(k) - sn * s(k)
        g = sn * s(k+1)
        s(k+1) = cs * s(k+1)

        if ( wantv ) then
          call drot ( n, v(1,k), 1, v(1,k+1), 1, cs, sn )
        end if

        call drotg ( f, g, cs, sn )
        s(k) = f
        f = cs * e(k) + sn * s(k+1)
        s(k+1) = -sn * e(k) + cs * s(k+1)
        g = sn * e(k+1)
        e(k+1) = cs * e(k+1)

        if ( wantu .and. k < m ) then
          call drot ( m, u(1,k), 1, u(1,k+1), 1, cs, sn )
        end if

      end do

      e(mn-1) = f
      iter = iter + 1
!
!  Convergence.
!
    else if ( kase == 4 ) then
!
!  Make the singular value nonnegative.
!
      if ( s(l) < 0.0D+00 ) then
        s(l) = -s(l)
        if ( wantv ) then
          v(1:n,l) = -v(1:n,l)
        end if
      end if
!
!  Order the singular value.
!
      do

        if ( l == mm ) then
          exit
        end if

        if ( s(l+1) <= s(l) ) then
          exit
        end if

        t = s(l)
        s(l) = s(l+1)
        s(l+1) = t

        if ( wantv .and. l < n ) then
          call dswap ( n, v(1,l), 1, v(1,l+1), 1 )
        end if

        if ( wantu .and. l < m ) then
          call dswap ( m, u(1,l), 1, u(1,l+1), 1 )
        end if

        l = l + 1

      end do

      iter = 0
      mn = mn - 1

    end if

  end do

  return
end
subroutine dswap ( n, x, incx, y, incy )

!*****************************************************************************80
!
!! DSWAP interchanges two vectors.
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input/output, real ( kind = 8 ) X(*), one of the vectors to swap.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive 
!    entries of X.
!
!    Input/output, real ( kind = 8 ) Y(*), one of the vectors to swap.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive 
!    elements of Y.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) temp
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  if ( n <= 0 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    m = mod ( n, 3 )

    do i = 1, m
      temp = x(i)
      x(i) = y(i)
      y(i) = temp
    end do

    do i = m + 1, n, 3

      temp = x(i)
      x(i) = y(i)
      y(i) = temp

      temp = x(i+1)
      x(i+1) = y(i+1)
      y(i+1) = temp

      temp = x(i+2)
      x(i+2) = y(i+2)
      y(i+2) = temp

    end do

  else

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      temp = x(ix)
      x(ix) = y(iy)
      y(iy) = temp
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) &
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end
subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to 
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine. 
!
!    It has been modified to produce the product Q' * Z, where Z is an input 
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
!    The changes consist (essentially) of applying the orthogonal 
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, real ( kind = 8 ) E(N), the subdiagonal entries of the 
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, real ( kind = 8 ) Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ), parameter :: itn = 30
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  real ( kind = 8 ) prec
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) z(n)

  prec = epsilon ( prec )

  if ( n == 1 ) then
    return
  end if

  e(n) = 0.0D+00

  do l = 1, n

    j = 0

    do

      do m = l, n

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        exit
      end if

      if ( itn <= j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IMTQLX - Fatal error!'
        write ( *, '(a)' ) '  Iteration limit exceeded.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  L = ', l
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  N = ', n
        stop
      end if

      j = j + 1
      g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
      r =  sqrt ( g * g + 1.0D+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0D+00
      c = 1.0D+00
      p = 0.0D+00
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r =  sqrt ( c * c + 1.0D+00 )
          e(i+1) = f * r
          s = 1.0D+00 / r
          c = c * s
        else
          s = f / g
          r =  sqrt ( s * s + 1.0D+00 )
          e(i+1) = g * r
          c = 1.0D+00 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0D+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0D+00

    end do

  end do
!
!  Sorting.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end if

  end do

  return
end
function r8_choose ( n, k )

!*****************************************************************************80
!
!! R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in R8 arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, K, are the values of N and K.
!
!    Output, real ( kind = 8 ) R8_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0.0D+00

  else if ( mn == 0 ) then

    value = 1.0D+00

  else

    mx = max ( k, n - k )
    value = real ( mx + 1, kind = 8 )

    do i = 2, mn
      value = ( value * real ( mx + i, kind = 8 ) ) / real ( i, kind = 8 )
    end do

  end if

  r8_choose = value

  return
end
function r8_sign ( x )

!*****************************************************************************80
!
!! R8_SIGN returns the sign of an R8.
!
!  Discussion:
!
!    value = -1 if X < 0;
!    value = +1 if X => 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number whose sign is desired.
!
!    Output, real ( kind = 8 ) R8_SIGN, the sign of X:
!
  implicit none

  real ( kind = 8 ) r8_sign
  real ( kind = 8 ) x

  if ( x < 0.0D+00 ) then
    r8_sign = -1.0D+00
  else
    r8_sign = +1.0D+00
  end if

  return
end
function r8vec_in_ab ( n, x, a, b )

!*****************************************************************************80
!
!! R8VEC_IN_AB is TRUE if the entries of an R8VEC are in the range [A,B].
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in X.
!
!    Input, real ( kind = 8 ) X(N), the vector.
!
!    Input, real A, B, the limits of the range.
!
!    Output, logical R8VEC_IN_AB, is TRUE if every entry of A is
!    between A and B.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical r8vec_in_ab
  real ( kind = 8 ) x(n)

  if ( any ( x(1:n) < a .or. b < x(1:n) ) ) then
    r8vec_in_ab = .false.
  else
    r8vec_in_ab = .true.
  end if

  return
end
subroutine r8vec_linspace ( n, a_first, a_last, a )

!*****************************************************************************80
!
!! R8VEC_LINSPACE creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A_FIRST, A_LAST, the first and last entries.
!
!    Output, real ( kind = 8 ) A(N), a vector of linearly spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_first
  real ( kind = 8 ) a_last
  integer ( kind = 4 ) i

  if ( n == 1 ) then

    a(1) = ( a_first + a_last ) / 2.0D+00

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = 8 ) * a_first &
             + real (     i - 1, kind = 8 ) * a_last ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT prints an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, a1(i), a2(i)
  end do

  return
end
subroutine svd_solve ( m, n, a, b, x )

!*****************************************************************************80
!
!! SVD_SOLVE solves a linear system in the least squares sense.
!
!  Discussion:
!
!    The vector X returned by this routine should always minimize the 
!    Euclidean norm of the residual ||A*x-b||.
!
!    If the matrix A does not have full column rank, then there are multiple
!    vectors that attain the minimum residual.  In that case, the vector
!    X returned by this routine is the unique such minimizer that has the 
!    the minimum possible Euclidean norm, that is, ||A*x-b|| and ||x||
!    are both minimized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2012
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, real ( kind = 8 ) B(M), the right hand side.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) a_copy(m,n)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) e(max(m+1,n))
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) job
  real ( kind = 8 ) sdiag(max(m+1,n))
  real ( kind = 8 ) smax
  real ( kind = 8 ) stol
  real ( kind = 8 ) sub(n)
  real ( kind = 8 ) u(m,m)
  real ( kind = 8 ) ub(m)
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
  real ( kind = 8 ) x(n)
!
!  Get the SVD.
!
  a_copy(1:m,1:n) = a(1:m,1:n)
  lda = m
  ldu = m
  ldv = n
  allocate ( work(1:m) )
  job = 11

  call dsvdc ( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, job, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVD_SOLVE - Failure!'
    write ( *, '(a)' ) '  The SVD could not be calculated.'
    write ( *, '(a)' ) '  LINPACK routine DSVDC returned a nonzero'
    write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
    stop
  end if

  ub(1:m) = matmul ( transpose ( u(1:m,1:m) ), b(1:m) )

  sub(1:n) = 0.0D+00
!
!  For singular problems, there may be tiny but nonzero singular values
!  that should be ignored.  This is a reasonable attempt to avoid such 
!  problems, although in general, the user might wish to control the tolerance.
!
  smax = maxval ( sdiag(1:n) )
  if ( smax <= epsilon ( smax ) ) then
    smax = 1.0D+00
  end if

  stol = epsilon ( smax ) * smax

  do i = 1, n
    if ( i <= m ) then
      if ( stol <= sdiag(i) ) then
        sub(i) = ub(i) / sdiag(i)
      end if
    end if
  end do

  x(1:n) = matmul ( v(1:n,1:n), sub(1:n) )

  return
end
function t_double_product_integral ( i, j )

!*****************************************************************************80
!
!! T_DOUBLE_PRODUCT_INTEGRAL: integral (-1<=x<=1) T(i,x)*T(j,x)/sqrt(1-x^2) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, the polynomial indices.
!    0 <= I, J.
!
!    Output, real ( kind = 8 ) T_DOUBLE_PRODUCT_INTEGRAL, the integral.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t_double_product_integral
  real ( kind = 8 ) value

  if ( i < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'T_DOUBLE_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= I is required.'
    stop
  end if

  if ( j < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'T_DOUBLE_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= J is required.'
    stop
  end if

  if ( i /= j ) then
    value = 0.0D+00
  elseif ( i == 0 ) then
    value = pi
  elseif ( 0 < i ) then
    value = pi / 2.0D+00
  end if

  t_double_product_integral = value

  return
end
function t_integral ( e )

!*****************************************************************************80
!
!! T_INTEGRAL: integral ( -1 <= x <= +1 ) x^e dx / sqrt ( 1 - x^2 ).
!
!  Discussion:
!
!    Set 
!      x = cos ( theta ), 
!      dx = - sin ( theta ) d theta = - sqrt ( 1 - x^2 ) d theta
!    to transform the integral to
!      integral ( 0 <= theta <= pi ) - ( cos ( theta ) )^e d theta
!    which becomes
!      0 if E is odd,
!      (1/2^e) * choose ( e, e/2 ) * pi if E is even.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) E, the exponent of X.
!    0 <= E.
!
!    Output, real ( kind = 8 ) T_INTEGRAL, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) e
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) t_integral
  real ( kind = 8 ) value

  if ( mod ( e, 2 ) == 1 ) then

    value = 0.0D+00

  else

    value = r8_choose ( e, e / 2 ) * pi / 2.0D+00 ** e

  end if

  t_integral = value

  return
end
subroutine t_polynomial ( m, n, x, v )

!*****************************************************************************80
!
!! T_POLYNOMIAL evaluates Chebyshev polynomials T(n,x).
!
!  Discussion:
!
!    Chebyshev polynomials are useful as a basis for representing the
!    approximation of functions since they are well conditioned, in the sense
!    that in the interval [-1,1] they each have maximum absolute value 1.
!    Hence an error in the value of a coefficient of the approximation, of
!    size epsilon, is exactly reflected in an error of size epsilon between
!    the computed approximation and the theoretical approximation.
!
!    Typical usage is as follows, where we assume for the moment
!    that the interval of approximation is [-1,1].  The value
!    of N is chosen, the highest polynomial to be used in the
!    approximation.  Then the function to be approximated is
!    evaluated at the N+1 points XJ which are the zeroes of the N+1-th
!    Chebyshev polynomial.  Let these values be denoted by F(XJ).
!
!    The coefficients of the approximation are now defined by
!
!      C(I) = 2/(N+1) * sum ( 1 <= J <= N+1 ) F(XJ) T(I,XJ)
!
!    except that C(0) is given a value which is half that assigned
!    to it by the above formula,
!
!    and the representation is
!
!    F(X) approximated by sum ( 0 <= J <= N ) C(J) T(J,X)
!
!    Now note that, again because of the fact that the Chebyshev polynomials
!    have maximum absolute value 1, if the higher order terms of the
!    coefficients C are small, then we have the option of truncating
!    the approximation by dropping these terms, and we will have an
!    exact value for maximum perturbation to the approximation that
!    this will cause.
!
!    It should be noted that typically the error in approximation
!    is dominated by the first neglected basis function (some multiple of
!    T(N+1,X) in the example above).  If this term were the exact error,
!    then we would have found the minimax polynomial, the approximating
!    polynomial of smallest maximum deviation from the original function.
!    The minimax polynomial is hard to compute, and another important
!    feature of the Chebyshev approximation is that it tends to behave
!    like the minimax polynomial while being easy to compute.
!
!    To evaluate a sum like 
!
!      sum ( 0 <= J <= N ) C(J) T(J,X), 
!
!    Clenshaw's recurrence formula is recommended instead of computing the
!    polynomial values, forming the products and summing.
!
!    Assuming that the coefficients C(J) have been computed
!    for J = 0 to N, then the coefficients of the representation of the
!    indefinite integral of the function may be computed by
!
!      B(I) = ( C(I-1) - C(I+1))/2*(I-1) for I=1 to N+1, 
!
!    with
! 
!      C(N+1)=0
!      B(0) arbitrary.  
!
!    Also, the coefficients of the representation of the derivative of the 
!    function may be computed by:
!
!      D(I) = D(I+2)+2*I*C(I) for I=N-1, N-2, ..., 0, 
!
!    with
!
!      D(N+1) = D(N)=0.
!
!    Some of the above may have to adjusted because of the irregularity of C(0).
!
!    The formula is:
!
!      T(N,X) = COS(N*ARCCOS(X))
!
!  Differential equation:
!
!    (1-X*X) Y'' - X Y' + N N Y = 0
!
!  First terms:
!
!    T(0,X) =  1
!    T(1,X) =  1 X
!    T(2,X) =  2 X^2 -   1
!    T(3,X) =  4 X^3 -   3 X
!    T(4,X) =  8 X^4 -   8 X^2 +  1
!    T(5,X) = 16 X^5 -  20 X^3 +  5 X
!    T(6,X) = 32 X^6 -  48 X^4 + 18 X^2 - 1
!    T(7,X) = 64 X^7 - 112 X^5 + 56 X^3 - 7 X
!
!  Inequality:
!
!    abs ( T(N,X) ) <= 1 for -1 <= X <= 1
!
!  Orthogonality:
!
!    For integration over [-1,1] with weight
!
!      W(X) = 1 / sqrt(1-X*X), 
!
!    if we write the inner product of T(I,X) and T(J,X) as
!
!      < T(I,X), T(J,X) > = integral ( -1 <= X <= 1 ) W(X) T(I,X) T(J,X) dX
!
!    then the result is:
!
!      < T(I,X), T(J,X) > = 0    if I /= J
!      < T(I,X), T(J,X) > = PI/2 if I == J /= 0
!      < T(I,X), T(J,X) > = PI   if I == J == 0
!
!    A discrete orthogonality relation is also satisfied at each of
!    the N zeroes of T(N,X):  sum ( 1 <= K <= N ) T(I,X) * T(J,X)
!                              = 0 if I /= J
!                              = N/2 if I == J /= 0
!                              = N if I == J == 0
!
!  Recursion:
!
!    T(0,X) = 1,
!    T(1,X) = X,
!    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
!
!    T'(N,X) = N * ( -X * T(N,X) + T(N-1,X) ) / ( 1 - X^2 )
!
!  Special values:
!
!    T(N,1) = 1
!    T(N,-1) = (-1)^N
!    T(2N,0) = (-1)^N
!    T(2N+1,0) = 0
!    T(N,X) = (-1)^N * T(N,-X)
!
!  Zeroes:
!
!    M-th zero of T(N,X) is X = cos((2*M-1)*PI/(2*N)), M = 1 to N.
!
!  Extrema:
!
!    M-th extremum of T(N,X) is X = cos(PI*M/N), M = 0 to N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest polynomial to compute.
!
!    Input, real ( kind = 8 ) X(1:M), the evaluation points.
!
!    Output, real ( kind = 8 ) V(1:M,0:N), the values of the polynomials.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) v(1:m,0:n)
  real ( kind = 8 ) x(1:m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = x(1:m)
 
  do j = 2, n
    v(1:m,j) = 2.0D+00 * x(1:m) * v(1:m,j-1) - v(1:m,j-2)
  end do
 
  return
end
subroutine t_polynomial_ab ( a, b, m, n, x, v )

!*****************************************************************************80
!
!! T_POLYNOMIAL_AB: evaluates Chebyshev polynomials T(n,x) in [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the domain of definition.
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest polynomial to compute.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!    It must be the case that A <= X(*) <= B.
!
!    Output, real ( kind = 8 ) V(M,N+1), the values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) j
  real ( kind = 8 ) v(1:m,0:n)
  real ( kind = 8 ) x(1:m)
  real ( kind = 8 ) y(1:m)

  y(1:m) = ( ( b - x(1:m)     )   &
           - (     x(1:m) - a ) ) &
           / ( b          - a )

  call t_polynomial ( m, n, x, v )
 
  return
end
subroutine t_polynomial_coefficients ( n, c )

!*****************************************************************************80
!
!! T_POLYNOMIAL_COEFFICIENTS: coefficients of the Chebyshev polynomial T(n,x).
!
!  First terms:
!
!    N/K     0     1      2      3       4     5      6    7      8    9   10
!
!     0      1
!     1      0     1
!     2     -1     0      2
!     3      0    -3      0      4
!     4      1     0     -8      0       8
!     5      0     5      0    -20       0    16
!     6     -1     0     18      0     -48     0     32
!     7      0    -7      0     56       0  -112      0    64
!
!  Recursion:
!
!    T(0,X) = 1,
!    T(1,X) = X,
!    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the Chebyshev T
!    polynomials.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  c(1,1) = 1.0D+00
 
  do i = 2, n
    c(i,0)     =                        - c(i-2,0)
    c(i,1:i-2) = 2.0D+00 * c(i-1,0:i-3) - c(i-2,1:i-2)
    c(i,  i-1) = 2.0D+00 * c(i-1,  i-2)
    c(i,  i  ) = 2.0D+00 * c(i-1,  i-1)
  end do
 
  return
end
function t_polynomial_value ( n, x )

!*****************************************************************************80
!
!! T_POLYNOMIAL_VALUE: returns the single value T(n,x).
!
!  Discussion:
!
!    In cases where calling T_POLYNOMIAL is inconvenient, because it returns
!    a vector of values for multiple arguments X, this simpler interface
!    may be appropriate.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Input, real ( kind = 8 ) X, the argument of the polynomial.
!
!    Output, real ( kind = 8 ) T_POLYNOMIAL_VALUE, the value of T(n,x).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) t_polynomial_value
  real ( kind = 8 ), allocatable :: vec(:)
  real ( kind = 8 ) x

  m = 1
  allocate ( vec(0:n) )

  call t_polynomial ( m, n, x, vec )

  t_polynomial_value = vec(n)

  deallocate ( vec )

  return
end
subroutine t_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! T_POLYNOMIAL_VALUES returns values of Chebyshev polynomials T(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      ChebyshevT[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1000000000000000D+01, &
     0.8000000000000000D+00, &
     0.2800000000000000D+00, &
    -0.3520000000000000D+00, &
    -0.8432000000000000D+00, &
    -0.9971200000000000D+00, &
    -0.7521920000000000D+00, &
    -0.2063872000000000D+00, &
     0.4219724800000000D+00, &
     0.8815431680000000D+00, &
     0.9884965888000000D+00, &
     0.7000513740800000D+00, &
     0.1315856097280000D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine t_polynomial_zeros ( n, z )

!*****************************************************************************80
!
!! T_POLYNOMIAL_ZEROS returns zeroes of the Chebyshev polynomial T(n,x).
!
!  Discussion:
!
!    The I-th zero of T(N,X) is cos((2*I-1)*PI/(2*N)), I = 1 to N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) Z(N), the zeroes of T(N,X).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) z(n)

  do i = 1, n
    angle = real ( 2 * i - 1, kind = 8 ) * pi / real ( 2 * n, kind = 8 );
    z(i) = cos ( angle );
  end do

  return
end
subroutine t_project_coefficients ( n, f, c )

!*****************************************************************************80
!
!! T_PROJECT_COEFFICIENTS: function projected onto Chebyshev polynomials T(n,x).
!
!  Discussion:
!
!    It is assumed that the interval of definition is -1 <= x <= +1.
!
!    Over this interval, f(x) will be well approximated by
!
!      f(x) approx sum ( 0 <= i <= n ) c(i) * T(i,x)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!
!    Input, external real ( kind = 8 ) function F ( X ), evaluates the function.
!
!    Output, real ( kind = 8 ) C(0:N), the projection coefficients of f(x) onto
!    T(0,x) through T(n,x).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n)
  real ( kind = 8 ) d(0:n)
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fac
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) total
  real ( kind = 8 ) y

  do k = 0, n
    y = cos ( pi * ( real ( k, kind = 8 ) + 0.5 ) / real ( n + 1, kind = 8 ) )
    d(k) = f ( y )
  end do

  fac = 2.0D+00 / real ( n + 1, kind = 8 )

  do j = 0, n
    total = 0.0D+00
    do k = 0, n
      total = total + d(k) * cos ( ( pi * real ( j, kind = 8 ) ) &
        * ( ( real ( k, kind = 8 ) + 0.5 ) / real ( n + 1, kind = 8 ) ) )
    end do
    c(j) = fac * total
  end do

  c(0) = c(0) / 2.0D+00

  return
end
subroutine t_project_coefficients_ab ( n, f, a, b, c )

!*****************************************************************************80
!
!! T_PROJECT_COEFFICIENTS_AB: function projected onto T(n,x) over [a,b].
!
!  Discussion:
!
!    It is assumed that the interval of definition is a <= x <= b.
!
!    Over this interval, f(x) will be well approximated by
!
!      f(x) approx sum ( 0 <= i <= n ) c(i) * T(i,(2x-a-b)/(b-a))
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!
!    Input, external real ( kind = 8 ) function F ( X ), evaluates the function.
!
!    Input, real ( kind = 8 ) A, B, the interval of definition.
!
!    Output, real ( kind = 8 ) C(0:N), the projection coefficients of f(x) onto
!    T(0,x) through T(n,x).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c(0:n)
  real ( kind = 8 ) d(0:n)
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fac
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) total
  real ( kind = 8 ) y

  do k = 0, n

    t = cos ( pi * ( real ( k, kind = 8 ) - 0.5 ) / real ( n + 1, kind = 8 ) )

    y = ( ( 1.0D+00 + t ) * b   &
        + ( 1.0D+00 - t ) * a ) &
        /   2.0D+00

    d(k) = f ( y )

  end do

  fac = 2.0D+00 / real ( n + 1, kind = 8 )

  do j = 0, n
    total = 0.0D+00
    do k = 0, n
      total = total + d(k) * cos ( ( pi * real ( j, kind = 8 ) ) &
        * ( ( real ( k, kind = 8 ) + 0.5 ) / real ( n + 1, kind = 8 ) ) )
    end do
    c(j) = fac * total
  end do

  c(0) = c(0) / 2.0D+00

  return
end
subroutine t_project_coefficients_data ( a, b, m, n, x, d, c )

!*****************************************************************************80
!
!! T_PROJECT_COEFFICIENTS_DATA: project data onto Chebyshev polynomials T(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the domain of definition.
!
!    Input, integer ( kind = 4 ) M, the number of data values.
!
!    Input, integer ( kind = 4 ) N, the desired order of the Chebyshev 
!    expansion.
!
!    Input, real ( kind = 8 ) X(M), the data abscissas.  These need not 
!    be sorted.  It must be the case that A <= X() <= B.
!
!    Input, real ( kind = 8 ) D(M), the data values.
!
!    Output, real ( kind = 8 ) C(0:N), the approximate Chebshev coefficients.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c(0:n)
  real ( kind = 8 ) d(m)
  logical r8vec_in_ab
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ) x(m)

  if ( .not. r8vec_in_ab ( m, x, a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' T_PROJECT_COEFFICIENTS_DATA- Fatal error!'
    write ( *, '(a)' ) '  Some X not in [A,B].'
    stop
  end if
!
!  Compute the M by N+1 Chebyshev Vandermonde matrix V.
!
  call t_polynomial_ab ( a, b, m, n, x, v )
!
!  Compute the least-squares solution C.
!
  call svd_solve ( m, n + 1, v, d, c )

  return
end
subroutine t_project_value ( m, n, x, c, v )

!*****************************************************************************80
!
!! T_PROJECT_VALUE evaluates an expansion in Chebyshev polynomials T(n,x).
!
!  Discussion:
!
!    The projection is assumed to be based on the interval [-1,+1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Input, real ( kind = 8 ) C(0:N), the expansion coefficients.
!
!    Output, real ( kind = 8 ) V(M), the value of the Chebyshev function.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n 

  real ( kind = 8 ) b0(m)
  real ( kind = 8 ) b1(m)
  real ( kind = 8 ) b2(m)
  real ( kind = 8 ) c(0:n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) v(m)
  real ( kind = 8 ) x(m)

  b1(1:m) = 0.0D+00
  b0(1:m) = 0.0D+00

  do j = n, 0, -1
    b2(1:m) = b1(1:m)
    b1(1:m) = b0(1:m)
    b0(1:m) = c(j) + 2.0D+00 * x(1:m) * b1(1:m) - b2(1:m)
  end do

  v(1:m) = 0.5D+00 * ( c(0) + b0(1:m) - b2(1:m) )

  return
end
subroutine t_project_value_ab ( m, n, x, c, a, b, v )

!*****************************************************************************80
!
!! T_PROJECT_VALUE_AB evaluates an expansion in Chebyshev polynomials T(n,x).
!
!  Discussion:
!
!    The projection is assumed to be based on the interval [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Input, real ( kind = 8 ) C(0:N), the expansion coefficients.
!
!    Input, real ( kind = 8 ) A, B, the interval of definition.
!
!    Output, real ( kind = 8 ) V(M), the value of the Chebyshev function.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n 

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) b0(m)
  real ( kind = 8 ) b1(m)
  real ( kind = 8 ) b2(m)
  real ( kind = 8 ) c(0:n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) v(m)
  real ( kind = 8 ) x(m)

  b1(1:m) = 0.0D+00
  b0(1:m) = 0.0D+00

  do j = n, 0, -1
    b2(1:m) = b1(1:m)
    b1(1:m) = b0(1:m)
    b0(1:m) = c(j) + 2.0D+00 / ( b - a ) * ( 2.0D+00 * x(1:m) - a - b ) &
      * b1(1:m) - b2(1:m)
  end do

  v(1:m) = 0.5D+00 * ( c(0) + b0(1:m) - b2(1:m) )

  return
end
subroutine t_quadrature_rule ( n, t, w )

!*****************************************************************************80
!
!! T_QUADRATURE_RULE: quadrature rule for T(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Output, real ( kind = 8 ) T(N), W(N), the points and weights of the rule.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bj(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) w(n)

  t(1:n) = 0.0D+00

  bj(1) = sqrt ( 0.5D+00 )
  bj(2:n) = 0.5D+00

  w(1) = sqrt ( pi )
  w(2:n) = 0.0D+00

  call imtqlx ( n, t, bj, w )

  w(1:n) = w(1:n)**2

  return
end
function t_triple_product_integral ( i, j, k )

!*****************************************************************************80
!
!! T_TRIPLE_PRODUCT_INTEGRAL: integral (-1<=x<=1) T(i,x)*T(j,x)*T(k,x)/sqrt(1-x^2) dx
!
!  Discussion:
!
!    
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Mason, David Handscomb,
!    Chebyshev Polynomials,
!    CRC Press, 2002,
!    ISBN: 0-8493-035509,
!    LC: QA404.5.M37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, K, the polynomial indices.
!    0 <= I, J.
!
!    Output, real ( kind = 8 ) T_TRIPLE_PRODUCT_INTEGRAL, the integral.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t_double_product_integral
  real ( kind = 8 ) t_triple_product_integral
  real ( kind = 8 ) value

  if ( i < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'T_TRIPLE_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= I is required.'
    stop
  end if

  if ( j < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'T_TRIPLE_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= J is required.'
    stop
  end if

  if ( k < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'T_TRIPLE_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= K is required.'
    stop
  end if

  value = 0.5D+00 * ( &
      t_double_product_integral (       i + j,   k ) + &
    + t_double_product_integral ( abs ( i - j ), k ) )

  t_triple_product_integral = value

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

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
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
function u_double_product_integral ( i, j )

!*****************************************************************************80
!
!! U_DOUBLE_PRODUCT_INTEGRAL: integral (-1<=x<=1) U(i,x)*U(j,x)*sqrt(1-x^2) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, the polynomial indices.
!    0 <= I, J.
!
!    Output, real ( kind = 8 ) U_DOUBLE_PRODUCT_INTEGRAL, the integral.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) u_double_product_integral
  real ( kind = 8 ) value

  if ( i < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'U_DOUBLE_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= I is required.'
    stop
  end if

  if ( j < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'U_DOUBLE_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= J is required.'
    stop
  end if

  if ( i /= j ) then
    value = 0.0D+00
  else
    value = pi / 2.0D+00
  end if

  u_double_product_integral = value

  return
end
function u_integral ( e )

!*****************************************************************************80
!
!! U_INTEGRAL: integral ( -1 <= x <= +1 ) x^e sqrt ( 1 - x^2 ) dx.
!
!  Discussion:
!
!     E    U_INTEGRAL
!    --    --------------
!     0         pi /    2   
!     2         pi /    8
!     4         pi /   16
!     6     5 * pi /  128
!     8     7 * pi /  256
!    10    21 * pi / 1024
!    12    33 * pi / 2048
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) E, the exponent of X.
!    0 <= E.
!
!    Output, real ( kind = 8 ) U_INTEGRAL, the value of the integral.
!
  implicit none

  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  integer ( kind = 4 ) e
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) u_integral
  real ( kind = 8 ) value

  if ( mod ( e, 2 ) == 1 ) then

    value = 0.0D+00

  else

    arg1 = 0.5D+00 * real ( 1 + e, kind = 8 )
    arg2 = 2.0D+00 + 0.5D+00 * real ( e, kind = 8 )
    value = 0.5D+00 * sqrt ( pi ) * gamma ( arg1 ) / gamma ( arg2 )

  end if

  u_integral = value

  return
end
subroutine u_polynomial ( m, n, x, v )

!*****************************************************************************80
!
!! U_POLYNOMIAL evaluates Chebyshev polynomials U(n,x).
!
!  Discussion:
!
!    The formula is:
!
!      If |X| <= 1, then
!
!        U(N,X) = sin ( (N+1) * arccos(X) ) / sqrt ( 1 - X^2 )
!               = sin ( (N+1) * arccos(X) ) / sin ( arccos(X) )
!
!      else
!
!        U(N,X) = sinh ( (N+1) * arccosh(X) ) / sinh ( arccosh(X) )
!
!  Differential equation:
!
!    (1-X*X) Y'' - 3 X Y' + N (N+2) Y = 0
!
!  First terms:
!
!    U(0,X) =   1
!    U(1,X) =   2 X
!    U(2,X) =   4 X^2 -   1
!    U(3,X) =   8 X^3 -   4 X
!    U(4,X) =  16 X^4 -  12 X^2 +  1
!    U(5,X) =  32 X^5 -  32 X^3 +  6 X
!    U(6,X) =  64 X^6 -  80 X^4 + 24 X^2 - 1
!    U(7,X) = 128 X^7 - 192 X^5 + 80 X^3 - 8X
!
!  Orthogonality:
!
!    For integration over [-1,1] with weight
!
!      W(X) = sqrt(1-X*X), 
!
!    we have
!
!      < U(I,X), U(J,X) > = integral ( -1 <= X <= 1 ) W(X) U(I,X) U(J,X) dX 
!
!    then the result is:
!
!      < U(I,X), U(J,X) >  =  0    if I /= J
!      < U(I,X), U(J,X) >  =  PI/2 if I == J
!
!  Recursion:
!
!    U(0,X) = 1,
!    U(1,X) = 2 * X,
!    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
!
!  Special values:
!
!    U(N,1) = N + 1
!    U(2N,0) = (-1)^N
!    U(2N+1,0) = 0
!    U(N,X) = (-1)^N * U(N,-X)
!
!  Zeroes:
!
!    M-th zero of U(N,X) is X = cos( M*PI/(N+1)), M = 1 to N
!
!  Extrema:
!
!    M-th extremum of U(N,X) is X = cos( M*PI/N), M = 0 to N
!
!  Norm:
!
!    Integral ( -1 <= X <= 1 ) ( 1 - X^2 ) * U(N,X)^2 dX = PI/2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest polynomial to compute.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) V(M,0:N), the values of the N+1 Chebyshev
!    polynomials.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = 2.0D+00 * x(1:m)

  do i = 2, n
    v(1:m,i) = 2.0D+00 * x(1:m) * v(1:m,i-1) - v(1:m,i-2)
  end do
 
  return
end
subroutine u_polynomial_ab ( a, b, m, n, x, v )

!*****************************************************************************80
!
!! U_POLYNOMIAL_AB: evaluates Chebyshev polynomials U(n,x) in [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the domain of definition.
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest polynomial to compute.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!    It must be the case that A <= X(*) <= B.
!
!    Output, real ( kind = 8 ) V(M,N+1), the values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) j
  real ( kind = 8 ) v(1:m,0:n)
  real ( kind = 8 ) x(1:m)
  real ( kind = 8 ) y(1:m)

  y(1:m) = ( ( b - x(1:m)     )   &
           - (     x(1:m) - a ) ) &
           / ( b          - a )

  call u_polynomial ( m, n, x, v )
 
  return
end
subroutine u_polynomial_coefficients ( n, c )

!*****************************************************************************80
!
!! U_POLYNOMIAL_COEFFICIENTS: coefficients of Chebyshev polynomials U(n,x).
!
!  First terms:
!
!    N/K     0     1      2      3       4     5      6    7      8    9   10
!
!     0      1
!     1      0     2
!     2     -1     0      4
!     3      0    -4      0      8
!     4      1     0    -12      0      16
!     5      0     6      0    -32       0    32
!     6     -1     0     24      0     -80     0     64
!     7      0    -8      0     80       0  -192      0   128
!
!  Recursion:
!
!    U(0,X) = 1,
!    U(1,X) = 2*X,
!    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  c(1,1) = 2.0D+00
 
  do i = 2, n
    c(i,0)     =                        - c(i-2,0)
    c(i,1:i-2) = 2.0D+00 * c(i-1,0:i-3) - c(i-2,1:i-2)
    c(i,  i-1) = 2.0D+00 * c(i-1,  i-2)
    c(i,  i  ) = 2.0D+00 * c(i-1,  i-1)
  end do
 
  return
end
subroutine u_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! U_POLYNOMIAL_VALUES returns values of Chebyshev polynomials U(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      ChebyshevU[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1000000000000000D+01, &
     0.1600000000000000D+01, &
     0.1560000000000000D+01, &
     0.8960000000000000D+00, &
    -0.1264000000000000D+00, &
    -0.1098240000000000D+01, &
    -0.1630784000000000D+01, &
    -0.1511014400000000D+01, &
    -0.7868390400000000D+00, &
     0.2520719360000000D+00, &
     0.1190154137600000D+01, &
     0.1652174684160000D+01, &
     0.1453325357056000D+01 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine u_polynomial_zeros ( n, z )

!*****************************************************************************80
!
!! U_POLYNOMIAL_ZEROS returns zeroes of Chebyshev polynomials U(n,x).
!
!  Discussion:
!
!    The I-th zero of U(N,X) is cos((I-1)*PI/(N-1)), I = 1 to N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) Z(N), the zeroes of U(N,X).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) z(n)

  do i = 1, n
    angle = real ( i, kind = 8 ) * pi / real ( n + 1, kind = 8 )
    z(i) = cos ( angle )
  end do

  return
end
subroutine u_quadrature_rule ( n, t, w )

!*****************************************************************************80
!
!! U_QUADRATURE_RULE: quadrature rule for U(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Output, real ( kind = 8 ) T(N), W(N), the points and weights of the rule.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bj(n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) w(n)

  t(1:n) = 0.0D+00

  bj(1:n) = 0.5D+00

  w(1) = sqrt ( pi / 2.0D+00 )
  w(2:n) = 0.0D+00

  call imtqlx ( n, t, bj, w )

  w(1:n) = w(1:n)**2

  return
end
function v_double_product_integral ( i, j )

!*****************************************************************************80
!
!! V_DOUBLE_PRODUCT_INTEGRAL: integral (-1<=x<=1) V(i,x)*V(j,x)*sqrt(1+x)/sqrt(1-x) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, the polynomial indices.
!    0 <= I, J.
!
!    Output, real ( kind = 8 ) V_DOUBLE_PRODUCT_INTEGRAL, the integral.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) v_double_product_integral
  real ( kind = 8 ) value

  if ( i < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'V_DOUBLE_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= I is required.'
    stop
  end if

  if ( j < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'V_DOUBLE_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= J is required.'
    stop
  end if

  if ( i /= j ) then
    value = 0.0D+00
  else
    value = pi
  end if

  v_double_product_integral = value

  return
end
subroutine v_polynomial ( m, n, x, v )

!*****************************************************************************80
!
!! V_POLYNOMIAL evaluates Chebyshev polynomials V(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest polynomial to compute.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) V(M,0:N), the values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = 2.0D+00 * x(1:m) - 1.0D+00

  do i = 2, n
    v(1:m,i) = 2.0D+00 * x(1:m) * v(1:m,i-1) - v(1:m,i-2)
  end do
 
  return
end
subroutine v_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! V_POLYNOMIAL_VALUES returns values of Chebyshev polynomials V(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      u = Sqrt[(x+1)/2],
!      ChebyshevT[2*n+1,u] / u
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     1.0000000000000000D+00, &
     0.6000000000000000D+00, &
    -0.0400000000000000D+00, &
    -0.6640000000000000D+00, &
    -1.0224000000000000D+00, &
    -0.9718400000000000D+00, &
    -0.5325440000000000D+00, &
     0.1197696000000000D+00, &
     0.7241753600000000D+00, &
     1.0389109760000000D+00, &
     0.9380822016000000D+00, &
     0.4620205465600000D+00, &
    -0.1988493271040000D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine v_polynomial_zeros ( n, z )

!*****************************************************************************80
!
!! V_POLYNOMIAL_ZEROS returns zeroes of Chebyshev polynomials V(n,x).
!
!  Discussion:
!
!    The I-th zero of U(N,X) is cos((I-1/2)*PI/(N+1/2)), I = 1 to N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) Z(N), the zeroes.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) z(n)

  do i = 1, n
    angle = real ( 2 * n - 2 * i + 1, kind = 8 ) * pi / real ( 2 * n + 1, kind = 8 )
    z(i) = cos ( angle )
  end do

  return
end
function w_double_product_integral ( i, j )

!*****************************************************************************80
!
!! W_DOUBLE_PRODUCT_INTEGRAL: integral (-1<=x<=1) W(i,x)*W(j,x)*sqrt(1-x)/sqrt(1+x) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, the polynomial indices.
!    0 <= I, J.
!
!    Output, real ( kind = 8 ) W_DOUBLE_PRODUCT_INTEGRAL, the integral.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) value
  real ( kind = 8 ) w_double_product_integral

  if ( i < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'W_DOUBLE_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= I is required.'
    stop
  end if

  if ( j < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'W_DOUBLE_PRODUCT_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  0 <= J is required.'
    stop
  end if

  if ( i /= j ) then
    value = 0.0D+00
  else
    value = pi
  end if

  w_double_product_integral = value

  return
end
subroutine w_polynomial ( m, n, x, v )

!*****************************************************************************80
!
!! W_POLYNOMIAL evaluates Chebyshev polynomials W(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest polynomial to compute.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) V(M,0:N), the values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = 2.0D+00 * x(1:m) + 1.0D+00

  do i = 2, n
    v(1:m,i) = 2.0D+00 * x(1:m) * v(1:m,i-1) - v(1:m,i-2)
  end do
 
  return
end
subroutine w_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! W_POLYNOMIAL_VALUES returns values of Chebyshev polynomials W(n,x).
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      u = Sqrt[(x+1)/2],
!      ChebyshevU[2*n,u]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 13

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     1.000000000000000D+00, &
     2.600000000000000D+00, &
     3.160000000000000D+00, &
     2.456000000000000D+00, &
     0.769600000000000D+00, &
    -1.224640000000000D+00, &
    -2.729024000000000D+00, &
    -3.141798400000000D+00, &
    -2.297853440000000D+00, &
    -0.534767104000000D+00, &
     1.442226073600000D+00, &
     2.842328821760000D+00, &
     3.105500041216000D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00, &
    0.8D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine w_polynomial_zeros ( n, z )

!*****************************************************************************80
!
!! W_POLYNOMIAL_ZEROS returns zeroes of Chebyshev polynomials W(n,x).
!
!  Discussion:
!
!    The I-th zero of U(N,X) is cos(I*PI/(N+1/2)), I = 1 to N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 April 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) Z(N), the zeroes.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) z(n)

  do i = 1, n
    angle = real ( 2 * ( n - i + 1 ), kind = 8 ) * pi / real ( 2 * n + 1, kind = 8 )
    z(i) = cos ( angle )
  end do

  return
end
