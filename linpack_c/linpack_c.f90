subroutine cchdc ( a, lda, p, work, ipvt, job, info )

!*****************************************************************************80
!
!! CCHDC: Cholesky decomposition of a Hermitian positive definite matrix.
!
!  Discussion:
!
!    A pivoting option allows the user to estimate the condition of a
!    Hermitian positive definite matrix or determine the rank of a
!    Hermitian positive semidefinite matrix.
!
!    For Hermitian positive definite matrices, INFO = P is the normal return.
!
!    For pivoting with Hermitian positive semidefinite matrices, INFO will
!    in general be less than P.  However, INFO may be greater than
!    the rank of A, since rounding error can cause an otherwise zero
!    element to be positive.  Indefinite systems will always cause
!    INFO to be less than P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2007
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) A(LDA,P).  On input, A contains the
!    matrix whose decomposition is to be computed.  Only the upper half of A
!    need be stored.  The lower part of the array A is not referenced.
!    On output, A contains in its upper half the Cholesky factor
!    of the matrix A as it has been permuted by pivoting.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) P, the order of the matrix.
!
!    Workspace, complex WORK(P).
!
!    Input/output, integer ( kind = 4 ) IPVT(P).  IPVT is not referenced if
!    JOB == 0.  On input, IPVT contains integers that control the selection of
!    the pivot elements, if pivoting has been requested.  Each diagonal element
!    A(K,K) is placed in one of three classes according to the input
!    value of IPVT(K):
!      IPVT(K) >  0, X(K) is an initial element.
!      IPVT(K) == 0, X(K) is a free element.
!      IPVT(K) <  0, X(K) is a final element.
!    Before the decomposition is computed, initial elements are moved by
!    symmetric row and column interchanges to the beginning of the array A
!    and final elements to the end.  Both initial and final elements
!    are frozen in place during the computation and only free elements
!    are moved.  At the K-th stage of the reduction, if A(K,K) is occupied
!    by a free element, it is interchanged with the largest free element
!    A(L,L) with K <= L.
!    On output, IPVT(K) contains the index of the diagonal element
!    of A that was moved into the J-th position, if pivoting was requested.
!
!    Input, integer ( kind = 4 ) JOB, specifies whether column pivoting is
!    to be done.
!    0, no pivoting is done.
!    nonzero, pivoting is done.
!
!    Output, integer ( kind = 4 ) INFO, contains the index of the last positive
!    diagonal element of the Cholesky factor.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) p

  complex ( kind = 4 ) a(lda,p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) ipvt(p)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kb
  integer ( kind = 4 ) l
  real    ( kind = 4 ) maxdia
  integer ( kind = 4 ) maxl
  logical              negk
  integer ( kind = 4 ) pl
  integer ( kind = 4 ) plp1
  integer ( kind = 4 ) pu
  logical              swapk
  complex ( kind = 4 ) t
  complex ( kind = 4 ) temp
  complex ( kind = 4 ) work(p)

  pl = 1
  pu = 0
  info = p

  if ( job /= 0 ) then
!
!  Pivoting has been requested.  Rearrange the elements according to IPVT.
!
    do k = 1, p

      swapk = ( 0 < ipvt(k) )
      negk = ( ipvt(k) < 0 )

      if ( negk ) then
        ipvt(k) = -k
      else
        ipvt(k) = k
      end if

      if ( swapk ) then

        if ( k /= pl ) then

          call cswap ( pl-1, a(1:pl-1,k), 1, a(1:pl-1,pl), 1 )

          t        = a(k,k)
          a(k,k)   = a(pl,pl)
          a(pl,pl) = t

          a(pl,k) = conjg ( a(pl,k) )
          plp1 = pl + 1

          do j = plp1, p

            if ( j < k ) then

              t       = conjg ( a(pl,j) )
              a(pl,j) = conjg ( a(j,k) )
              a(j,k)  = t

            else if ( j /= k ) then

              t       = a(k,j)
              a(k,j)  = a(pl,j)
              a(pl,j) = t

            end if

          end do

          ipvt(k) = ipvt(pl)
          ipvt(pl) = k

        end if

        pl = pl + 1

      end if

    end do

    pu = p

    do kb = pl, p

      k = p - kb + pl

      if ( ipvt(k) < 0 ) then

        ipvt(k) = -ipvt(k)

        if ( pu /= k ) then

          call cswap ( k-1, a(1:k-1,k), 1, a(1:k-1,pu), 1 )

          t        = a(k,k)
          a(k,k)   = a(pu,pu)
          a(pu,pu) = t

          a(k,pu) = conjg ( a(k,pu) )

          do j = k + 1, p

            if ( j < pu ) then

              t       = conjg ( a(k,j) )
              a(k,j)  = conjg ( a(j,pu) )
              a(j,pu) = t

            else if ( j /= pu ) then

              t       = a(k,j)
              a(k,j)  = a(pu,j)
              a(pu,j) = t

            end if

          end do

          i        = ipvt(k)
          ipvt(k)  = ipvt(pu)
          ipvt(pu) = i

        end if

        pu = pu - 1

      end if

    end do

  end if

  do k = 1, p
!
!  Reduction loop.
!
    maxdia = real ( a(k,k) )
    maxl = k
!
!  Determine the pivot element.
!
    if ( pl <= k .and. k < pu ) then
      do l = k+1, pu
        if ( maxdia < real ( a(l,l) ) ) then
          maxdia = real ( a(l,l) )
          maxl = l
        end if
      end do
    end if
!
!  Quit if the pivot element is not positive.
!
    if ( maxdia <= 0.0E+00 ) then
      info = k - 1
      return
    end if
!
!  Start the pivoting and update IPVT.
!
    if ( k /= maxl ) then

      call cswap ( k-1, a(1:k-1,k), 1, a(1:k-1,maxl), 1 )
      a(maxl,maxl) = a(k,k)
      a(k,k) = cmplx ( maxdia, 0.0E+00 )

      i          = ipvt(maxl)
      ipvt(maxl) = ipvt(k)
      ipvt(k)    = i

      a(k,maxl) = conjg ( a(k,maxl) )

    end if
!
!  Reduction step.  Pivoting is contained across the rows.
!
    work(k) = cmplx ( sqrt ( real ( a(k,k) ) ), 0.0E+00 )
    a(k,k) = work(k)

    do j = k+1, p

      if ( k /= maxl ) then

        if ( j < maxl ) then

          t         = conjg ( a(k,j) )
          a(k,j)    = conjg ( a(j,maxl) )
          a(j,maxl) = t

        else if ( j /= maxl ) then

          t         = a(k,j)
          a(k,j)    = a(maxl,j)
          a(maxl,j) = t

        end if

      end if

      a(k,j) = a(k,j) / work(k)
      work(j) = conjg ( a(k,j) )
      temp = -a(k,j)
      a(k+1:j,j) = a(k+1:j,j) + temp * work(k+1:j)

    end do

  end do

  return
end
subroutine cchdd ( r, ldr, p, x, z, ldz, nz, y, rho, c, s, info )

!*****************************************************************************80
!
!! CCHDD downdates an augmented Cholesky decomposition.
!
!  Discussion:
!
!    CCHDD downdates an augmented Cholesky decomposition or the
!    triangular factor of an augmented QR decomposition.
!    Specifically, given an upper triangular matrix R of order P,  a
!    row vector X, a column vector Z, and a scalar Y, CCHDD
!    determines a unitary matrix U and a scalar ZETA such that
!
!          ( R   Z  )     ( RR  ZZ )
!      U * (        )  =  (        ),
!          ( 0 ZETA )     (  X   Y )
!
!    where RR is upper triangular.  If R and Z have been obtained
!    from the factorization of a least squares problem, then
!    RR and ZZ are the factors corresponding to the problem
!    with the observation (X,Y) removed.  In this case, if RHO
!    is the norm of the residual vector, then the norm of
!    the residual vector of the downdated problem is
!      sqrt ( RHO**2 - ZETA**2 ).
!    CCHDD will simultaneously downdate several triplets (Z,Y,RHO)
!    along with R.
!
!    For a less terse description of what CCHDD does and how
!    it may be applied, see the LINPACK guide.
!
!    The matrix U is determined as the product U(1)*...*U(P)
!    where U(I) is a rotation in the (P+1,I)-plane of the
!    form
!
!      ( C(I)  -conjg ( S(I) ) )
!      (                       ).
!      ( S(I)           C(I)   )
!
!    The rotations are chosen so that C(I) is real.
!
!    The user is warned that a given downdating problem may
!    be impossible to accomplish or may produce
!    inaccurate results.  For example, this can happen
!    if X is near a vector whose removal will reduce the
!    rank of R.  Beware.
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
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) R(LDR,P); on input, the upper
!    triangular matrix that is to be downdated.  On output, the downdated
!    matrix.  The part of R below the diagonal is not referenced.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of R.  P <= LDR.
!
!    Input, integer ( kind = 4 ) P, the order of the matrix.
!
!    Input, complex ( kind = 4 ) X(P), the row vector that is to
!    be removed from R.
!
!    Input/output, complex ( kind = 4 ) Z(LDZ,NZ); on input, an array of NZ
!    P-vectors which are to be downdated along with R.  On output,
!    the downdated vectors.
!
!    Input, integer ( kind = 4 ) LDZ, the leading dimension of Z.  P <= LDZ.
!
!    Input, integer ( kind = 4 ) NZ, the number of vectors to be downdated.
!    NZ may be zero, in which case Z, Y, and R are not referenced.
!
!    Input, complex ( kind = 4 ) Y(NZ), the scalars for the downdating
!    of the vectors Z.
!
!    Input/output, real ( kind = 4 ) RHO(NZ).  On input, the norms of the residual
!    vectors that are to be downdated.  On output, the downdated norms.
!
!    Output, real ( kind = 4 ) C(P), the cosines of the transforming rotations.
!
!    Output, complex ( kind = 4 ) S(P), the sines of the transforming rotations.
!
!    Output, integer ( kind = 4 ) INFO:
!     0, if the entire downdating was successful.
!    -1, if R could not be downdated.  In this case, all quantities
!        are left unaltered.
!     1, if some RHO could not be downdated.  The offending RHO's are
!        set to -1.
!
  implicit none

  integer ( kind = 4 ) ldr
  integer ( kind = 4 ) ldz
  integer ( kind = 4 ) nz
  integer ( kind = 4 ) p

  real    ( kind = 4 ) a
  real    ( kind = 4 ) alpha
  real    ( kind = 4 ) azeta
  complex ( kind = 4 ) b
  real    ( kind = 4 ) c(p)
  complex ( kind = 4 ) cdotc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  real    ( kind = 4 ) norm
  complex ( kind = 4 ) r(ldr,p)
  real    ( kind = 4 ) rho(nz)
  complex ( kind = 4 ) s(p)
  real    ( kind = 4 ) scale
  real    ( kind = 4 ) scnrm2
  complex ( kind = 4 ) t
  complex ( kind = 4 ) x(p)
  complex ( kind = 4 ) xx
  complex ( kind = 4 ) y(nz)
  complex ( kind = 4 ) z(ldz,nz)
  complex ( kind = 4 ) zeta
!
!  Solve the system hermitian(R) * A = X, placing the result in S.
!
  info = 0
  s(1) = conjg ( x(1) ) / conjg ( r(1,1) )
!
!  The FORTRAN90 DOT_PRODUCT automatically conjugates the first vector.
!
  do j = 2, p
    s(j) = conjg ( x(j) ) - dot_product ( r(1:j-1,j), s(1:j-1) )
    s(j) = s(j) / conjg ( r(j,j) )
  end do

  norm = scnrm2 ( p, s, 1 )

  if ( 1.0E+00 <= norm ) then
    info = -1
    return
  end if

  alpha = sqrt ( 1.0E+00 - norm**2 )
!
!  Determine the transformations.
!
  do ii = 1, p
    i = p - ii + 1
    scale = alpha + abs ( s(i) )
    a = alpha / scale
    b = s(i) / scale
    norm = sqrt ( a**2 + ( real ( b ) )**2 + ( aimag ( b ) )**2 )
    c(i) = a / norm
    s(i) = conjg ( b ) / norm
    alpha = scale * norm
  end do
!
!  Apply the transformations to R.
!
  do j = 1, p
    xx = cmplx ( 0.0E+00, 0.0E+00 )
    do ii = 1, j
      i = j - ii + 1
      t = c(i) * xx + s(i) * r(i,j)
      r(i,j) = c(i) * r(i,j) - conjg ( s(i) ) * xx
      xx = t
    end do
  end do
!
!  If required, downdate Z and RHO.
!
  do j = 1, nz

    zeta = y(j)

    do i = 1, p
      z(i,j) = ( z(i,j) - conjg ( s(i) ) * zeta ) / c(i)
      zeta = c(i) * zeta - s(i) * z(i,j)
    end do

    azeta = abs ( zeta )

    if ( rho(j) < azeta ) then
      info = 1
      rho(j) = -1.0E+00
    else
      rho(j) = rho(j) * sqrt ( 1.0E+00 - ( azeta / rho(j) )**2 )
    end if

  end do

  return
end
subroutine cchex ( r, ldr, p, k, l, z, ldz, nz, c, s, job )

!*****************************************************************************80
!
!! CCHEX updates a Cholesky factorization.
!
!  Discussion:
!
!    CCHEX updates a Cholesky factorization
!
!      A = hermitian(R) * R
!
!    of a positive definite matrix A of order P under diagonal
!    permutations of the form
!
!      E' * A * E
!
!    where E is a permutation matrix.  Specifically, given
!    an upper triangular matrix R and a permutation matrix
!    E (which is specified by K, L, and JOB), CCHEX determines
!    a unitary matrix U such that
!
!      U * R * E = RR,
!
!    where RR is upper triangular.  At the user's option, the
!    transformation U will be multiplied into the array Z.
!
!    If A = hermitian(X)*X, so that R is the triangular part of the
!    QR factorization of X, then RR is the triangular part of the
!    QR factorization of X * E, that is, X with its columns permuted.
!
!    For a less terse description of what CCHEX does and how
!    it may be applied, see the LINPACK guide.
!
!    The matrix Q is determined as the product U(L-K)*...*U(1)
!    of plane rotations of the form
!
!      (    C(I)       S(I) )
!      (                    ) ,
!      ( -conjg(S(i))  C(I) )
!
!    where C(I) is real, the rows these rotations operate on
!    are described below.
!
!    There are two types of permutations, which are determined
!    by the value of job.
!
!    JOB = 1, right circular shift:
!    The columns are rearranged in the following order.
!
!      1, ..., K-1, L, K, K+1, ..., L-1, L+1, ..., P.
!
!    U is the product of L-K rotations U(I), where U(I)
!    acts in the (L-I,L-I+1)-plane.
!
!    JOB = 2, left circular shift:
!    The columns are rearranged in the following order
!
!      1, ..., K-1, K+1, K+2, ..., L, L, L+1, ..., P.
!
!    U is the product of L-K rotations U(I), where U(I)
!    acts in the (K+I-1,K+I)-plane.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) R(LDR,P); On input, the upper triangular factor
!    that is to be updated.  On output, the updated factor.  Elements
!    below the diagonal are not referenced.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of R, which is at least P.
!
!    Input, integer ( kind = 4 ) P, the order of the matrix.
!
!    Input, integer ( kind = 4 ) K, the first column to be permuted.
!
!    Input, integer ( kind = 4 ) L, the last column to be permuted.
!    L must be strictly greater than K.
!
!    Input/output, complex ( kind = 4 ) Z(LDZ,NZ); on input, an array of NZ
!    P-vectors into which the transformation U is multiplied.  On output, the
!    updated matrix.  Z is not referenced if NZ = 0.
!
!    Input, integer ( kind = 4 ) LDZ, the leading dimension of Z, which must
!    be at least P.
!
!    Input, integer ( kind = 4 ) NZ, the number of columns of the matrix Z.
!
!    Output, real ( kind = 4 ) C(P), the cosines of the transforming rotations.
!
!    Output, complex ( kind = 4 ) S(P), the sines of the transforming rotations.
!
!    Input, integer ( kind = 4 ) JOB, determines the type of permutation.
!    1, right circular shift.
!    2, left circular shift.
!
  implicit none

  integer ( kind = 4 ) ldr
  integer ( kind = 4 ) ldz
  integer ( kind = 4 ) nz
  integer ( kind = 4 ) p

  real    ( kind = 4 ) c(p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) il
  integer ( kind = 4 ) iu
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  complex ( kind = 4 ) r(ldr,p)
  complex ( kind = 4 ) s(p)
  complex ( kind = 4 ) t
  complex ( kind = 4 ) z(ldz,nz)

  if ( job == 1 ) then
!
!  Right circular shift.
!
!  Reorder the columns.
!
    do i = 1, l
      ii = l - i + 1
      s(i) = r(ii,l)
    end do

    do jj = k, l - 1
      j = l - 1 - jj + k
      r(1:j,j+1) = r(1:j,j)
      r(j+1,j+1) = cmplx ( 0.0E+00, 0.0E+00 )
    end do

    do i = 1, k - 1
      ii = l - i + 1
      r(i,k) = s(ii)
    end do
!
!  Calculate the rotations.
!
    t = s(1)
    do i = 1, l - k
      call crotg ( s(i+1), t, c(i), s(i) )
      t = s(i+1)
    end do

    r(k,k) = t
    do j = k+1, p
      il = max ( 1, l - j + 1 )
      do ii = il, l - k
        i = l - ii
        t = c(ii) * r(i,j) + s(ii) * r(i+1,j)
        r(i+1,j) = c(ii) * r(i+1,j) - conjg ( s(ii) ) * r(i,j)
        r(i,j) = t
      end do
    end do
!
!  If required, apply the transformations to Z.
!
    do j = 1, nz
      do ii = 1, l - k
        i = l - ii
        t = c(ii) * z(i,j) + s(ii) * z(i+1,j)
        z(i+1,j) = c(ii) * z(i+1,j) - conjg ( s(ii) ) * z(i,j)
        z(i,j) = t
      end do
    end do

  else
!
!  Left circular shift.
!
!  Reorder the columns.
!
    do i = 1, k
      ii = l - k + i
      s(ii) = r(i,k)
    end do

    do j = k, l - 1
      r(1:j,j) = r(1:j,j+1)
      jj = j - k + 1
      s(jj) = r(j+1,j+1)
    end do

    do i = 1, k
      ii = l - k + i
      r(i,l) = s(ii)
    end do

    r(k+1:l,l) = cmplx ( 0.0E+00, 0.0E+00 )
!
!  Reduction loop.
!
    do j = k, p
!
!  Apply the rotations.
!
      if ( j /= k ) then
        iu = min ( j - 1, l - 1 )
        do i = k, iu
          ii = i - k + 1
          t = c(ii) * r(i,j) + s(ii) * r(i+1,j)
          r(i+1,j) = c(ii) * r(i+1,j) - conjg ( s(ii) ) * r(i,j)
          r(i,j) = t
        end do
      end if

      if ( j < l ) then
        jj = j - k + 1
        t = s(jj)
        call crotg ( r(j,j), t, c(jj), s(jj) )
      end if

    end do
!
!  Apply the rotations to Z.
!
    do j = 1, nz
      do i = k, l - 1
        ii = i - k + 1
        t = c(ii) * z(i,j) + s(ii) * z(i+1,j)
        z(i+1,j) = c(ii) * z(i+1,j) - conjg ( s(ii) ) * z(i,j)
        z(i,j) = t
      end do
    end do

  end if

  return
end
subroutine cchud ( r, ldr, p, x, z, ldz, nz, y, rho, c, s )

!*****************************************************************************80
!
!! CCHUD updates an augmented Cholesky decomposition.
!
!  Discussion:
!
!    CCHUD updates an augmented Cholesky decomposition of the
!    triangular part of an augmented QR decomposition.  Specifically,
!    given an upper triangular matrix R of order P, a row vector
!    X, a column vector Z, and a scalar Y, CCHUD determines a
!    unitary matrix U and a scalar ZETA such that
!
!           ( R  Z )     ( RR   ZZ  )
!      U  * (      )  =  (          ),
!           ( X  Y )     (  0  ZETA )
!
!    where RR is upper triangular.  If R and Z have been
!    obtained from the factorization of a least squares
!    problem, then RR and ZZ are the factors corresponding to
!    the problem with the observation (X,Y) appended.  In this
!    case, if RHO is the norm of the residual vector, then the
!    norm of the residual vector of the updated problem is
!    sqrt ( RHO**2 + ZETA**2 ).  CCHUD will simultaneously update
!    several triplets (Z,Y,RHO).
!
!    For a less terse description of what CCHUD does and how
!    it may be applied see the LINPACK guide.
!
!    The matrix U is determined as the product U(P)*...*U(1),
!    where U(I) is a rotation in the (I,P+1) plane of the
!    form
!
!      (          C(I)    S(I) )
!      (                       ).
!      ( -conjg ( S(I) )  C(I) )
!
!    The rotations are chosen so that C(I) is real.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) R(LDR,P), the upper triangular matrix
!    that is to be updated.  The part of R below the diagonal is
!    not referenced.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of R.
!    P <= LDR.
!
!    Input, integer ( kind = 4 ) P, the order of the matrix.
!
!    Input, complex ( kind = 4 ) X(P), the row to be added to R.
!
!    Input/output, complex ( kind = 4 ) Z(LDZ,NZ), NZ P-vectors to
!    be updated with R.
!
!    Input, integer ( kind = 4 ) LDZ, the leading dimension of Z.
!    P <= LDZ.
!
!    Integer, integer NZ, the number of vectors to be updated.
!    NZ may be zero, in which case Z, Y, and RHO are not referenced.
!
!    Input, complex ( kind = 4 ) Y(NZ), the scalars for updating the vectors Z.
!
!    Input/output, real ( kind = 4 ) RHO(NZ); on input, the norms of the residual
!    vectors that are to be updated.  If RHO(J) is negative, it is
!    left unaltered.  On output, the updated values.
!
!    Output, real ( kind = 4 ) C(P). the cosines of the transforming rotations.
!
!    Output, complex ( kind = 4 ) S(P), the sines of the transforming rotations.
!
  implicit none

  integer ( kind = 4 ) ldr
  integer ( kind = 4 ) ldz
  integer ( kind = 4 ) nz
  integer ( kind = 4 ) p

  real    ( kind = 4 ) azeta
  real    ( kind = 4 ) c(p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  complex ( kind = 4 ) r(ldr,p)
  real    ( kind = 4 ) rho(nz)
  complex ( kind = 4 ) s(p)
  real    ( kind = 4 ) scale
  complex ( kind = 4 ) t
  complex ( kind = 4 ) x(p)
  complex ( kind = 4 ) xj
  complex ( kind = 4 ) y(nz)
  complex ( kind = 4 ) z(ldz,nz)
  complex ( kind = 4 ) zeta
!
!  Update R.
!
  do j = 1, p

    xj = x(j)
!
!  Apply the previous rotations.
!
    do i = 1, j - 1
      t =  c(i) * r(i,j) +         s(i)   * xj
      xj = c(i) * xj     - conjg ( s(i) ) * r(i,j)
      r(i,j) = t
    end do
!
!  Compute the next rotation.
!
    call crotg ( r(j,j), xj, c(j), s(j) )

  end do
!
!  If required, update Z and RHO.
!
  do j = 1, nz

    zeta = y(j)

    do i = 1, p
      t =    c(i) * z(i,j) +         s(i)   * zeta
      zeta = c(i) * zeta   - conjg ( s(i) ) * z(i,j)
      z(i,j) = t
    end do

    azeta = abs ( zeta )

    if ( azeta /= 0.0E+00 .and. 0.0E+00 <= rho(j) ) then
      scale = azeta + rho(j)
      rho(j) = scale * sqrt ( ( azeta / scale )**2 + ( rho(j) / scale )**2 )
    end if

  end do

  return
end
subroutine cgbco ( abd, lda, n, ml, mu, ipvt, rcond, z )

!*****************************************************************************80
!
!! CGBCO factors a complex band matrix and estimates its condition.
!
!  Discussion:
!
!    If RCOND is not needed, CGBFA is slightly faster.
!
!    To solve A*X = B, follow CGBCO by CGBSL.
!
!    To compute inverse(A)*C, follow CGBCO by CGBSL.
!
!    To compute determinant(A), follow CGBCO by CGBDI.
!
!  Band storage:
!
!    If A is a band matrix, the following program segment
!    will set up the input.
!
!      ml = (band width below the diagonal)
!      mu = (band width above the diagonal)
!      m = ml + mu + 1
!      do j = 1, n
!        i1 = max ( 1, j - mu )
!        i2 = min ( n, j + ml )
!        do i = i1, i2
!          k = i - j + m
!          abd(k,j) = a(i,j)
!        end do
!      end do
!
!    This uses rows ML+1 through 2*ML+MU+1 of ABD.
!    In addition, the first ML rows in ABD are used for
!    elements generated during the triangularization.
!    The total number of rows needed in ABD is 2*ML+MU+1.
!    The ML+MU by ML+MU upper left triangle and the
!    ML by ML lower right triangle are not referenced.
!
!  Example:
!
!    If the original matrix A is
!
!      11 12 13  0  0  0
!      21 22 23 24  0  0
!       0 32 33 34 35  0
!       0  0 43 44 45 46
!       0  0  0 54 55 56
!       0  0  0  0 65 66
!
!     Then N = 6, ML = 1, MU = 2, 5 <= LDA and ABD should contain
!
!       *  *  *  +  +  +
!       *  * 13 24 35 46
!       * 12 23 34 45 56
!      11 22 33 44 55 66
!      21 32 43 54 65  *
!
!    * = not used,
!    + = used for pivoting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) ABD(LDA,N), on input, contains the matrix in
!    band storage.  The columns of the matrix are stored in the columns
!    of ABD and the diagonals of the matrix are stored in rows ML+1
!    through 2*ML+MU+1 of ABD.  On output, an upper triangular matrix
!    in band storage and the multipliers which were used to obtain it.
!    The factorization can be written A = L*U where L is a product of
!    permutation and unit lower triangular matrices and U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ABD.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, the number of diagonals below the main diagonal.
!    0 <= ML < N.
!
!    Input, integer ( kind = 4 ) MU, the number of diagonals above the main diagonal.
!    0 <= MU < N.
!    More efficient if ML <= MU.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, real ( kind = 4 ) RCOND, an estimate of the reciprocal condition of A.
!    For the system A*X = B, relative perturbations in A and B of size
!    epsilon may cause relative perturbations in X of size (EPSILON/RCOND).
!    If RCOND is so small that the logical expression
!      1.0 + RCOND == 1.0
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate
!    underflows.
!
!    Workspace, complex Z(N), a work vector whose contents are usually
!    unimportant.  If A is close to a singular matrix, then Z is an
!    approximate null vector in the sense that
!      norm ( A * Z ) = RCOND * norm ( A ) * norm ( Z ).
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) abd(lda,n)
  real    ( kind = 4 ) anorm
  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) cdotc
  complex ( kind = 4 ) csign1
  complex ( kind = 4 ) ek
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ju
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lm
  integer ( kind = 4 ) lz
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) mu
  real    ( kind = 4 ) rcond
  real    ( kind = 4 ) s
  real    ( kind = 4 ) scasum
  real    ( kind = 4 ) sm
  complex ( kind = 4 ) t
  complex ( kind = 4 ) wk
  complex ( kind = 4 ) wkm
  real    ( kind = 4 ) ynorm
  complex ( kind = 4 ) z(n)
!
!  Compute 1-norm of A.
!
  anorm = 0.0E+00
  l = ml + 1
  is = l + mu

  do j = 1, n

    anorm = max ( anorm, scasum ( l, abd(is:is+l-1,j), 1 ) )

    if ( ml + 1 < is ) then
      is = is - 1
    end if

    if ( j <= mu ) then
      l = l + 1
    end if

    if ( n - ml <= j ) then
      l = l - 1
    end if

  end do
!
!  Factor
!
  call cgbfa ( abd, lda, n, ml, mu, ipvt, info )
!
!  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
!
!  Estimate = norm(Z)/norm(Y) where A*Z = Y and hermitian(A)*Y = E.
!
!  Hermitian(A) is the conjugate transpose of A.
!
!  The components of E are chosen to cause maximum local
!  growth in the elements of W where hermitian(U)*W = E.
!
!  The vectors are frequently rescaled to avoid overflow.
!
!  Solve hermitian(U) * W = E.
!
  ek = cmplx ( 1.0E+00, 0.0E+00 )

  z(1:n) = cmplx ( 0.0E+00, 0.0E+00 )

  m = ml + mu + 1
  ju = 0

  do k = 1, n

    if ( cabs1 ( z(k) ) /= 0.0E+00 ) then
      ek = csign1 ( ek, -z(k) )
    end if

    if ( cabs1 ( abd(m,k) ) < cabs1 ( ek - z(k) ) ) then
      s = cabs1 ( abd(m,k) ) / cabs1 ( ek - z(k) )
      z(1:n) = z(1:n) * s
      ek = cmplx ( s, 0.0E+00 ) * ek
    end if

    wk = ek - z(k)
    wkm = -ek - z(k)
    s = cabs1 ( wk )
    sm = cabs1 ( wkm )

    if ( cabs1 ( abd(m,k) ) /= 0.0E+00 ) then
      wk = wk / conjg ( abd(m,k) )
      wkm = wkm / conjg ( abd(m,k) )
    else
      wk = cmplx ( 1.0E+00, 0.0E+00 )
      wkm = cmplx ( 1.0E+00, 0.0E+00 )
    end if

    ju = min ( max ( ju, mu + ipvt(k) ), n )
    mm = m

    if ( k+1 <= ju ) then

      do j = k+1, ju
        mm = mm - 1
        sm = sm + cabs1 ( z(j) + wkm * conjg ( abd(mm,j) ) )
        z(j) = z(j) + wk * conjg ( abd(mm,j) )
        s = s + cabs1 ( z(j) )
      end do

      if ( s < sm ) then
        t = wkm - wk
        wk = wkm
        mm = m
        do j = k+1, ju
          mm = mm - 1
          z(j) = z(j) + t * conjg ( abd(mm,j) )
        end do
      end if

    end if

    z(k) = wk

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
!
!  Solve hermitian(L) * Y = W.
!
  do k = n, 1, -1

    lm = min ( ml, n - k )

    if ( k < n ) then
!     z(k) = z(k) + cdotc ( lm, abd(m+1:m+lm,k), 1, z(k+1:k+lm), 1 )
      z(k) = z(k) + dot_product ( abd(m+1:m+lm,k), z(k+1:k+lm) )
    end if

    if ( 1.0E+00 < cabs1 ( z(k) ) ) then
      s = 1.0E+00 / cabs1 ( z(k) )
      z(1:n) = z(1:n) * s
    end if

    l = ipvt(k)

    t    = z(l)
    z(l) = z(k)
    z(k) = t

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = 1.0E+00
!
!  Solve L * V = Y.
!
  do k = 1, n

    l = ipvt(k)

    t = z(l)
    z(l) = z(k)
    z(k) = t

    lm = min ( ml, n - k )

    if ( k < n ) then
      z(k+1:k+lm) = z(k+1:k+lm) + t * abd(m+1:m+lm,k)
    end if

    if ( 1.0E+00 < cabs1 ( z(k) ) ) then
      s = 1.0E+00 / cabs1 ( z(k) )
      z(1:n) = z(1:n) * s
      ynorm = s * ynorm
    end if

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm
!
!  Solve U * Z = W.
!
  do k = n, 1, -1

    if ( cabs1 ( abd(m,k) ) < cabs1 ( z(k) ) ) then
      s = cabs1 ( abd(m,k) ) / cabs1 ( z(k) )
      z(1:n) = z(1:n) * s
      ynorm = s * ynorm
    end if

    if ( cabs1 ( abd(m,k) ) /= 0.0E+00 ) then
      z(k) = z(k) / abd(m,k)
    else
      z(k) = cmplx ( 1.0E+00, 0.0E+00 )
    end if

    lm = min ( k, m ) - 1
    la = m - lm
    lz = k - lm
    t = -z(k)
    z(lz:lz+lm-1) = z(lz:lz+lm-1) + t * abd(la:la+lm-1,k)

  end do
!
!  Make ZNORM = 1.
!
  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm

  if ( anorm /= 0.0E+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0E+00
  end if

  return
end
subroutine cgbdi ( abd, lda, n, ml, mu, ipvt, det )

!*****************************************************************************80
!
!! CGBDI computes the determinant of a band matrix factored by CGBCO or CGBFA.
!
!  Discussion:
!
!    If the inverse is needed, use CGBSL N times.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input, complex ( kind = 4 ) ABD(LDA,N), the output from CGBCO or CGBFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, the number of diagonals below the main diagonal.
!
!    Input, integer ( kind = 4 ) MU, the number of diagonals above the main diagonal.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from CGBCO or CGBFA.
!
!    Output, complex ( kind = 4 ) DET(2), determinant of original matrix.
!    Determinant = DET(1) * 10.0**DET(2) with 1.0 <= cabs1 ( DET(1) ) < 10.0
!    or DET(1) = 0.0.  Also, DET(2) is strictly real.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) abd(lda,n)
  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu

  m = ml + mu + 1
  det(1) = cmplx ( 1.0E+00, 0.0E+00 )
  det(2) = cmplx ( 0.0E+00, 0.0E+00 )

  do i = 1, n

    if ( ipvt(i) /= i ) then
      det(1) = -det(1)
    end if

    det(1) = det(1) * abd(m,i)

    if ( cabs1 ( det(1) ) == 0.0E+00 ) then
      exit
    end if

    do while ( cabs1 ( det(1) ) < 1.0E+00 )
      det(1) = det(1) * cmplx ( 10.0E+00, 0.0E+00 )
      det(2) = det(2) - cmplx ( 1.0E+00, 0.0E+00 )
    end do

    do while ( 10.0E+00 <= cabs1 ( det(1) ) )
      det(1) = det(1) / cmplx ( 10.0E+00, 0.0E+00 )
      det(2) = det(2) + cmplx ( 1.0E+00, 0.0E+00 )
    end do

  end do

  return
end
subroutine cgbfa ( abd, lda, n, ml, mu, ipvt, info )

!*****************************************************************************80
!
!! CGBFA factors a complex band matrix by elimination.
!
!  Discussion:
!
!    CGBFA is usually called by CGBCO, but it can be called
!    directly with a saving in time if RCOND is not needed.
!
!  Band storage:
!
!    If A is a band matrix, the following program segment
!    will set up the input.
!
!      ml = (band width below the diagonal)
!      mu = (band width above the diagonal)
!      m = ml + mu + 1
!      do j = 1, n
!        i1 = max ( 1, j - mu )
!        i2 = min ( n, j + ml )
!        do i = i1, i2
!          k = i - j + m
!          abd(k,j) = a(i,j)
!        end do
!      end do
!
!    This uses rows ML+1 through 2*ML+MU+1 of ABD.
!    In addition, the first ML rows in ABD are used for
!    elements generated during the triangularization.
!    The total number of rows needed in ABD is 2*ML+MU+1.
!    The ML+MU by ML+MU upper left triangle and the
!    ML by ML lower right triangle are not referenced.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) ABD(LDA,N), on input, contains the matrix in
!    band storage.  The columns of the matrix are stored in the columns
!    of ABD and the diagonals of the matrix are stored in rows ML+1
!    through 2*ML+MU+1 of ABD.  On output, an upper triangular matrix
!    in band storage and the multipliers which were used to obtain it.
!    The factorization can be written A = L*U where L is a product of
!    permutation and unit lower triangular matrices and U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ABD.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, the number of diagonals below the main diagonal.
!    0 <= ML < N.
!
!    Input, integer ( kind = 4 ) MU, the number of diagonals above the main diagonal.
!    0 <= MU < N.  More efficient if ML <= MU.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO.
!    0, normal value.
!    K, if U(K,K) == 0.0.  This is not an error condition for this
!    subroutine, but it does indicate that CGBSL will divide by zero if
!    called.  Use RCOND in CGBCO for a reliable indication of singularity.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) abd(lda,n)
  real    ( kind = 4 ) cabs1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) icamax
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j0
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) ju
  integer ( kind = 4 ) jz
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lm
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) mu
  complex ( kind = 4 ) t

  m = ml + mu + 1
  info = 0
!
!  Zero initial fill-in columns.
!
  j0 = mu + 2
  j1 = min ( n, m ) - 1

  do jz = j0, j1
    i0 = m + 1 - jz
    do i = i0, ml
      abd(i,jz) = cmplx ( 0.0E+00, 0.0E+00 )
    end do
  end do

  jz = j1
  ju = 0
!
!  Gaussian elimination with partial pivoting.
!
  do k = 1, n-1
!
!  Zero next fill-in column
!
    jz = jz + 1
    if ( jz <= n ) then
      abd(1:ml,jz) = cmplx ( 0.0E+00, 0.0E+00 )
    end if
!
!  Find L = pivot index.
!
    lm = min ( ml, n - k )
    l = icamax ( lm+1, abd(m:m+lm,k), 1 ) + m - 1
    ipvt(k) = l + k - m
!
!  Zero pivot implies this column already triangularized.
!
    if ( cabs1 ( abd(l,k) ) == 0.0E+00 ) then
      info = k
      cycle
    end if
!
!  Interchange if necessary.
!
    if ( l /= m ) then
      t        = abd(l,k)
      abd(l,k) = abd(m,k)
      abd(m,k) = t
    end if
!
!  Compute multipliers.
!
    t = - cmplx ( 1.0E+00, 0.0E+00  ) / abd(m,k)
    abd(m+1:m+lm,k) = abd(m+1:m+lm,k) * t
!
!  Row elimination with column indexing.
!
    ju = min ( max ( ju, mu + ipvt(k) ), n )
    mm = m

    do j = k+1, ju
      l = l - 1
      mm = mm - 1
      t = abd(l,j)
      if ( l /= mm ) then
        abd(l,j) = abd(mm,j)
        abd(mm,j) = t
      end if
      abd(mm+1:mm+lm,j) = abd(mm+1:mm+lm,j) + t * abd(m+1:m+lm,k)
    end do

  end do

  ipvt(n) = n

  if ( cabs1 ( abd(m,n) ) == 0.0E+00 ) then
    info = n
  end if

  return
end
subroutine cgbsl ( abd, lda, n, ml, mu, ipvt, b, job )

!*****************************************************************************80
!
!! CGBSL solves a complex band system factored by CGBCO or CGBFA.
!
!  Discussion:
!
!    CGBSL can solve A * X = B or hermitan ( A ) * X = B.
!
!    A division by zero will occur if the input factor contains a
!    zero on the diagonal.  Technically this indicates singularity
!    but it is often caused by improper arguments or improper
!    setting of LDA.  It will not occur if the subroutines are
!    called correctly and if CGBCO has set 0.0 < RCOND
!    or CGBFA has set INFO = 0.
!
!    To compute inverse ( A ) * C where C is a matrix with P columns:
!
!      call cgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
!
!      if ( rcond is not too small ) then
!        do j = 1, p
!          call cgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
!        end do
!      end if
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input, complex ( kind = 4 ) ABD(LDA,N), the output from CGBCO or CGBFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ABD.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, the number of diagonals below the main diagonal.
!
!    Input, integer ( kind = 4 ) MU, the number of diagonals above the main diagonal.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from CGBCO or CGBFA.
!
!    Input/output, complex ( kind = 4 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
!    Input, integer ( kind = 4 ) JOB.
!    0, to solve A*x = b,
!    nonzero, to solve hermitian(A)*x = b, where hermitian(A) is the
!    conjugate transpose.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) abd(lda,n)
  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) cdotc
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) lm
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  complex ( kind = 4 ) t

  m = mu + ml + 1

  if ( job == 0 ) then
!
!  JOB = 0, solve A * X = B.
!
!  First solve L * Y = B.
!
    if ( ml /= 0 ) then

      do k = 1, n - 1

        lm = min ( ml, n - k )
        l = ipvt(k)
        t = b(l)

        if ( l /= k ) then
          b(l) = b(k)
          b(k) = t
        end if

        b(k+1:k+lm) = b(k+1:k+lm) + t * abd(m+1:m+lm,k)

      end do

    end if
!
!  Now solve U * X = Y.
!
    do k = n, 1, -1
      b(k) = b(k) / abd(m,k)
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm
      t = -b(k)
      b(lb:lb+lm-1) = b(lb:lb+lm-1) + t * abd(la:la+lm-1,k)
    end do

  else
!
!  JOB = nonzero, solve hermitian(A) * X = B.
!
!  First solve hermitian ( U ) * Y = B.
!
    do k = 1, n
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm
!     t = cdotc ( lm, abd(la:la+lm-1,k), 1, b(lb:lb+lm-1), 1 )
      t = dot_product ( abd(la:la+lm-1,k), b(lb:lb+lm-1) )
      b(k) = ( b(k) - t ) / conjg ( abd(m,k) )
    end do
!
!  Now solve hermitian ( L ) * X = Y.
!
    if ( ml /= 0 ) then

      do k = n-1, 1, -1

        lm = min ( ml, n - k )
!       b(k) = b(k) + cdotc ( lm, abd(m+1:m+lm,k), 1, b(k+1:k+lm), 1 )
        b(k) = b(k) + dot_product ( abd(m+1:m+lm,k), b(k+1:k+lm) )
        l = ipvt(k)

        if ( l /= k ) then
          t    = b(l)
          b(l) = b(k)
          b(k) = t
        end if

      end do

    end if

  end if

  return
end
subroutine cgeco ( a, lda, n, ipvt, rcond, z )

!*****************************************************************************80
!
!! CGECO factors a complex matrix and estimates its condition.
!
!  Discussion:
!
!    If RCOND is not needed, CGEFA is slightly faster.
!
!    To solve A*X = B, follow CGECO by CGESL.
!
!    To compute inverse(A)*C, follow CGECO by CGESL.
!
!    To compute determinant(A), follow CGECO by CGEDI.
!
!    To compute inverse(A), follow CGECO by CGEDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 )A(LDA,N), on input, the matrix to be factored.
!    On output, an upper triangular matrix and the multipliers which were
!    used to obtain it.  The factorization can be written A = L*U where
!    L is a product of permutation and unit lower triangular matrices
!    and U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, real ( kind = 4 ) RCOND, an estimate of the reciprocal condition of A.
!    For the system A*X = B, relative perturbations in A and B of size
!    EPSILON may cause relative perturbations in X of size (EPSILON/RCOND).
!    If RCOND is so small that the logical expression
!      1.0 + RCOND == 1.0
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate
!    underflows.
!
!    Workspace, complex Z(N), a work vector whose contents are usually
!    unimportant.  If A is close to a singular matrix, then Z is
!    an approximate null vector in the sense that
!      norm ( A * Z ) = RCOND * norm ( A ) * norm ( Z ).
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  real    ( kind = 4 ) anorm
  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) cdotc
  complex ( kind = 4 ) csign1
  complex ( kind = 4 ) ek
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real    ( kind = 4 ) rcond
  real    ( kind = 4 ) s
  real    ( kind = 4 ) scasum
  real    ( kind = 4 ) sm
  complex ( kind = 4 ) t
  complex ( kind = 4 ) wk
  complex ( kind = 4 ) wkm
  real    ( kind = 4 ) ynorm
  complex ( kind = 4 ) z(n)
!
!  Compute the 1-norm of A.
!
  anorm = 0.0E+00
  do j = 1, n
    anorm = max ( anorm, scasum ( n, a(1:n,j), 1 ) )
  end do
!
!  Factor.
!
  call cgefa ( a, lda, n, ipvt, info )
!
!  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
!
!  Estimate = norm(Z)/norm(Y) where A*Z = Y and hermitian(A)*Y = E.
!
!  Hermitian(A) is the conjugate transpose of A.
!
!  The components of E are chosen to cause maximum local
!  growth in the elements of W where hermitian(U)*W = E.
!
!  The vectors are frequently rescaled to avoid overflow.
!
!  Solve hermitian(U)*W = E.
!
  ek = cmplx ( 1.0E+00, 0.0E+00 )
  z(1:n) = cmplx ( 0.0E+00, 0.0E+00 )

  do k = 1, n

    if ( cabs1 ( z(k) ) /= 0.0E+00 ) then
      ek = csign1 ( ek, -z(k) )
    end if

    if ( cabs1 ( a(k,k) ) < cabs1 ( ek - z(k) ) ) then
      s = cabs1 ( a(k,k) ) / cabs1 ( ek - z(k) )
      z(1:n) = z(1:n) * s
      ek = cmplx ( s, 0.0E+00 ) * ek
    end if

    wk = ek - z(k)
    wkm = -ek - z(k)
    s = cabs1 ( wk )
    sm = cabs1 ( wkm )

    if ( cabs1 ( a(k,k) ) /= 0.0E+00 ) then
      wk = wk / conjg ( a(k,k) )
      wkm = wkm / conjg ( a(k,k) )
    else
      wk = cmplx ( 1.0E+00, 0.0E+00 )
      wkm = cmplx ( 1.0E+00, 0.0E+00 )
    end if

    do j = k+1, n
      sm = sm + cabs1 ( z(j) + wkm * conjg ( a(k,j) ) )
      z(j) = z(j) + wk * conjg ( a(k,j) )
      s = s + cabs1 ( z(j) )
    end do

    if ( s < sm ) then
      t = wkm - wk
      wk = wkm
      z(k+1:n) = z(k+1:n) + t * conjg ( a(k,k+1:n) )
    end if

    z(k) = wk

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
!
!  Solve hermitian(L) * Y = W.
!
  do k = n, 1, -1

    if ( k < n ) then
!     z(k) = z(k) + cdotc ( n-k, a(k+1:n,k), 1, z(k+1:n), 1 )
      z(k) = z(k) + dot_product ( a(k+1:n,k), z(k+1:n) )
    end if

    if ( 1.0E+00 < cabs1 ( z(k) ) ) then
      s = 1.0E+00 / cabs1 ( z(k) )
      z(1:n) = z(1:n) * s
    end if

    l = ipvt(k)

    t    = z(l)
    z(l) = z(k)
    z(k) = t

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = 1.0E+00
!
!  Solve L * V = Y.
!
  do k = 1, n

    l = ipvt(k)

    t = z(l)
    z(l) = z(k)
    z(k) = t

    if ( k < n ) then
      z(k+1:n) = z(k+1:n) + t * a(k+1:n,k)
    end if

    if ( 1.0E+00 < cabs1 ( z(k) ) ) then
      s = 1.0E+00 / cabs1 ( z(k) )
      z(1:n) = z(1:n) * s
      ynorm = s * ynorm
    end if

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm
!
!  Solve U * Z = V.
!
  do k = n, 1, -1

    if ( cabs1 ( a(k,k) ) < cabs1 ( z(k) ) ) then
      s = cabs1 ( a(k,k) ) / cabs1 ( z(k) )
      z(1:n) = z(1:n) * s
      ynorm = s * ynorm
    end if

    if ( cabs1 ( a(k,k) ) /= 0.0E+00 ) then
      z(k) = z(k) / a(k,k)
    else
      z(k) = cmplx ( 1.0E+00, 0.0E+00 )
    end if

    t = -z(k)
    z(1:k-1) = z(1:k-1) + t * a(1:k-1,k)

  end do
!
!  Make ZNORM = 1.
!
  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm

  if ( anorm /= 0.0E+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0E+00
  end if

  return
end
subroutine cgedi ( a, lda, n, ipvt, det, work, job )

!*****************************************************************************80
!
!! CGEDI computes the determinant and inverse of a matrix.
!
!  Discussion:
!
!    The matrix must have been factored by CGECO or CGEFA.
!
!    A division by zero will occur if the input factor contains
!    a zero on the diagonal and the inverse is requested.
!    It will not occur if the subroutines are called correctly
!    and if CGECO has set 0.0 < RCOND or CGEFA has set
!    INFO == 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 April 2007
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) A(LDA,N); on input, the factor information
!    from CGECO or CGEFA.  On output, the inverse matrix, if it
!    was requested,
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from CGECO or CGEFA.
!
!    Output, complex ( kind = 4 ) DET(2), the determinant of the original matrix,
!    if requested.  Otherwise not referenced.
!    Determinant = DET(1) * 10.0**DET(2) with
!    1.0 <= cabs1 ( DET(1) ) < 10.0 or DET(1) == 0.0.
!    Also, DET(2) is strictly real.
!
!    Workspace, complex WORK(N).
!
!    Input, integer ( kind = 4 ) JOB.
!    11, both determinant and inverse.
!    01, inverse only.
!    10, determinant only.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  complex ( kind = 4 ) t
  complex ( kind = 4 ) work(n)
!
!  Compute the determinant.
!
  if ( job / 10 /= 0 ) then

    det(1) = cmplx ( 1.0E+00, 0.0E+00 )
    det(2) = cmplx ( 0.0E+00, 0.0E+00 )

    do i = 1, n

      if ( ipvt(i) /= i ) then
        det(1) = -det(1)
      end if

      det(1) = det(1) * a(i,i)

      if ( cabs1 ( det(1) ) == 0.0E+00 ) then
        exit
      end if

      do while ( cabs1 ( det(1) ) < 1.0E+00 )
        det(1) = det(1) * cmplx ( 10.0E+00, 0.0E+00 )
        det(2) = det(2) - cmplx ( 1.0E+00, 0.0E+00 )
      end do

      do while ( 10.0E+00 <= cabs1 ( det(1) ) )
        det(1) = det(1) / cmplx ( 10.0E+00, 0.0E+00 )
        det(2) = det(2) + cmplx ( 1.0E+00, 0.0E+00 )
      end do

    end do

  end if
!
!  Compute inverse(U).
!
  if ( mod ( job, 10 ) /= 0 ) then

    do k = 1, n

      a(k,k) = cmplx ( 1.0E+00, 0.0E+00 ) / a(k,k)
      t = -a(k,k)
      a(1:k-1,k) = a(1:k-1,k) * t

      do j = k+1, n
        t = a(k,j)
        a(k,j) = cmplx ( 0.0E+00, 0.0E+00 )
        a(1:k,j) = a(1:k,j) + t * a(1:k,k)
      end do

    end do
!
!  Form inverse(U) * inverse(L).
!
    do k = n-1, 1, -1

      work(k+1:n) = a(k+1:n,k)
      a(k+1:n,k) = cmplx ( 0.0E+00, 0.0E+00 )

      do j = k+1, n
        t = work(j)
        a(1:n,k) = a(1:n,k) + t * a(1:n,j)
      end do

      l = ipvt(k)

      if ( l /= k ) then
        work(1:n) = a(1:n,k)
        a(1:n,k)  = a(1:n,l)
        a(1:n,l)  = work(1:n)
      end if

    end do

  end if

  return
end
subroutine cgefa ( a, lda, n, ipvt, info )

!*****************************************************************************80
!
!! CGEFA factors a complex matrix by Gaussian elimination.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) A(LDA,N); on input, the matrix to be factored.
!    On output, an upper triangular matrix and the multipliers which were
!    used to obtain it.  The factorization can be written A = L*U where
!    L is a product of permutation and unit lower triangular matrices and
!    U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO,
!    0, normal value.
!    K, if U(K,K) == 0.0.  This is not an error condition for this
!    subroutine, but it does indicate that CGESL or CGEDI will divide by zero
!    if called.  Use RCOND in CGECO for a reliable indication of singularity.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  real    ( kind = 4 ) cabs1
  integer ( kind = 4 ) icamax
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  complex ( kind = 4 ) t
!
!  Gaussian elimination with partial pivoting.
!
  info = 0

  do k = 1, n - 1
!
!  Find L = pivot index.
!
    l = icamax ( n-k+1, a(k:n,k), 1 ) + k - 1
    ipvt(k) = l
!
!  Zero pivot implies this column already triangularized.
!
    if ( cabs1 ( a(l,k) ) == 0.0E+00 ) then
      info = k
      cycle
    end if
!
!  Interchange if necessary.
!
    if ( l /= k ) then
      t      = a(l,k)
      a(l,k) = a(k,k)
      a(k,k) = t
    end if
!
!  Compute multipliers
!
    t = - cmplx ( 1.0E+00, 0.0E+00 ) / a(k,k)
    a(k+1:n,k) = a(k+1:n,k) * t
!
!  Row elimination with column indexing
!
    do j = k+1, n
      t = a(l,j)
      if ( l /= k ) then
        a(l,j) = a(k,j)
        a(k,j) = t
      end if
      a(k+1:n,j) = a(k+1:n,j) + t * a(k+1:n,k)
    end do

  end do

  ipvt(n) = n

  if ( cabs1 ( a(n,n) ) == 0.0E+00 ) then
    info = n
  end if

  return
end
subroutine cgesl ( a, lda, n, ipvt, b, job )

!*****************************************************************************80
!
!! CGESL solves a complex system factored by CGECO or CGEFA.
!
!  Discussion:
!
!    A division by zero will occur if the input factor contains a
!    zero on the diagonal.  Technically this indicates singularity
!    but it is often caused by improper arguments or improper
!    setting of LDA.  It will not occur if the subroutines are
!    called correctly and if CGECO has set 0.0 < RCOND
!    or CGEFA has set INFO == 0.
!
!    To compute inverse(A) * C where C is a matrix with P columns:
!
!      call cgeco(a,lda,n,ipvt,rcond,z)
!
!      if (rcond is not too small) then
!        do j = 1, p
!          call cgesl ( a, lda, n, ipvt, c(1,j), 0 )
!        end do
!      end if
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input, complex ( kind = 4 ) A(LDA,N), the factored matrix information,
!    as output from CGECO or CGEFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from CGECO or CGEFA.
!
!    Input/output, complex ( kind = 4 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
!    Input, integer ( kind = 4 ) JOB.
!    0, to solve A*X = B,
!    nonzero, to solve hermitian(A)*X = B where hermitian(A) is the
!    conjugate transpose.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) cdotc
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  complex ( kind = 4 ) t

  if ( job == 0 ) then
!
!  JOB = 0, solve A * X = B.
!
!  First solve L * Y = B.
!
    do k = 1, n-1
      l = ipvt(k)
      t = b(l)
      if ( l /= k ) then
        b(l) = b(k)
        b(k) = t
      end if
      b(k+1:n) = b(k+1:n) + t * a(k+1:n,k)
    end do
!
!  Now solve U * X = Y.
!
    do k = n, 1, -1
      b(k) = b(k) / a(k,k)
      t = -b(k)
      b(1:k-1) = b(1:k-1) + t * a(1:k-1,k)
    end do

  else
!
!  JOB nonzero, solve hermitian(A) * X = B.
!
!  First solve hermitian(U) * Y = B.
!
    do k = 1, n
!     t = cdotc ( k-1, a(1:k-1,k), 1, b(1:k-1), 1 )
      t = dot_product ( a(1:k-1,k), b(1:k-1) )
      b(k) = ( b(k) - t ) / conjg ( a(k,k) )
    end do
!
!  Now solve hermitian(L) * X = Y.
!
    do k = n-1, 1, -1
!     b(k) = b(k) + cdotc ( n-k, a(k+1:n,k), 1, b(k+1:n), 1 )
      b(k) = b(k) + dot_product ( a(k+1:n,k), b(k+1:n) )
      l = ipvt(k)
      if ( l /= k ) then
        t    = b(l)
        b(l) = b(k)
        b(k) = t
      end if
    end do

  end if

  return
end
subroutine cgtsl ( n, c, d, e, b, info )

!*****************************************************************************80
!
!! CGTSL solves a complex general tridiagonal system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, complex ( kind = 4 ) C(N); on input, the subdiagonal of the
!    tridiagonal matrix in entries C(2:N).  On output, C has been overwritten.
!
!    Input/output, complex ( kind = 4 ) D(N); on input, the diagonal of the
!    tridiagonal matrix.  On output, D has been overwritten.
!
!    Input/output, complex ( kind = 4 ) E(N); on input, the superdiagonal of the
!    tridiagonal matrix in entries E(1:N-1).  On output, E has been overwritten.
!
!    Input/output, complex ( kind = 4 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
!    Output, integer ( kind = 4 ) INFO.
!    0, normal value.
!    K, if the K-th element of the diagonal becomes exactly zero.  The
!    subroutine returns when this is detected.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) c(n)
  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) d(n)
  complex ( kind = 4 ) e(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) k
  complex ( kind = 4 ) t

  info = 0
  c(1) = d(1)

  if ( 1 <= n-1 ) then

    d(1) = e(1)
    e(1) = cmplx ( 0.0E+00, 0.0E+00 )
    e(n) = cmplx ( 0.0E+00, 0.0E+00 )

    do k = 1, n-1

      if ( cabs1 ( c(k) ) <= cabs1 ( c(k+1) ) ) then

        t      = b(k)
        b(k)   = b(k+1)
        b(k+1) = t

        t      = c(k)
        c(k)   = c(k+1)
        c(k+1) = t

        t      = d(k)
        d(k)   = d(k+1)
        d(k+1) = t

        t      = e(k)
        e(k)   = e(k+1)
        e(k+1) = t

      end if

      if ( cabs1 ( c(k) ) == 0.0E+00 ) then
        info = k
        return
      end if

      t = -c(k+1) / c(k)
      c(k+1) = d(k+1) + t * d(k)
      d(k+1) = e(k+1) + t * e(k)
      e(k+1) = cmplx ( 0.0E+00, 0.0E+00 )
      b(k+1) = b(k+1) + t * b(k)

    end do

  end if

  if ( cabs1 ( c(n) ) == 0.0E+00 ) then
    info = n
    return
  end if
!
!  Back solve.
!
  b(n) = b(n) / c(n)

  if ( 1 < n ) then

    b(n-1) = ( b(n-1) - d(n-1) * b(n) ) / c(n-1)

    do k = n-2, 1, -1
      b(k) = ( b(k) - d(k) * b(k+1) - e(k) * b(k+2) ) / c(k)
    end do

  end if

  return
end
subroutine chico ( a, lda, n, ipvt, rcond, z )

!*****************************************************************************80
!
!! CHICO factors a complex hermitian matrix and estimates its condition.
!
!  Discussion:
!
!    If RCOND is not needed, CHIFA is slightly faster.
!
!    To solve A*X = B, follow CHICO by CHISL.
!
!    To compute inverse(A)*C, follow CHICO by CHISL.
!
!    To compute inverse(A), follow CHICO by CHIDI.
!
!    To compute determinant(A), follow CHICO by CHIDI.
!
!    To compute inertia(A), follow CHICO by CHIDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2007
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) A(LDA,N); on input, the hermitian matrix to be
!    factored.  On output, a block diagonal matrix and the multipliers
!    which were used to obtain it.  The factorization can be written
!    A = U*D*hermitian(U) where U is a product of permutation and unit
!    upper triangular matrices, hermitian(U) is the conjugate transpose
!    of U, and D is block diagonal with 1 by 1 and 2 by 2 blocks.
!    Only the diagonal and upper triangle are used.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, real ( kind = 4 ) RCOND, an estimate of the reciprocal condition of
!    the matrix.  For the system A*X = B, relative perturbations in A and B
!    of size EPSILON may cause relative perturbations in X of size
!    (EPSILON/RCOND).  If RCOND is so small that the logical expression
!      1.0 + RCOND == 1.0
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate underflows.
!
!    Workspace, complex Z(N), a work vector whose contents are usually
!    unimportant.  If A is close to a singular matrix, then Z is an
!    approximate null vector in the sense that
!      norm(A*Z) = RCOND * norm(A) * norm(Z).
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  complex ( kind = 4 ) ak
  complex ( kind = 4 ) akm1
  real    ( kind = 4 ) anorm
  complex ( kind = 4 ) bk
  complex ( kind = 4 ) bkm1
  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) cdotc
  complex ( kind = 4 ) csign1
  complex ( kind = 4 ) denom
  complex ( kind = 4 ) ek
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp
  integer ( kind = 4 ) kps
  integer ( kind = 4 ) ks
  real    ( kind = 4 ) rcond
  real    ( kind = 4 ) s
  real    ( kind = 4 ) scasum
  complex ( kind = 4 ) t
  real    ( kind = 4 ) ynorm
  complex ( kind = 4 ) z(n)
!
!  Find norm of A using only upper half.
!
  do j = 1, n

    z(j) = cmplx ( scasum ( j, a(1:j,j), 1 ), 0.0E+00 )

    do i = 1, j - 1
      z(i) = cmplx ( real ( z(i) ) + cabs1 ( a(i,j) ), 0.0E+00 )
    end do

  end do

  anorm = 0.0E+00
  do j = 1, n
    anorm = max ( anorm, real ( z(j) ) )
  end do
!
!  Factor.
!
  call chifa ( a, lda, n, ipvt, info )
!
!  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
!
!  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
!
!  The components of E are chosen to cause maximum local
!  growth in the elements of W where U*D*W = E.
!
!  The vectors are frequently rescaled to avoid overflow.
!
!  Solve U*D*W = E.
!
  ek = cmplx ( 1.0E+00, 0.0E+00 )
  z(1:n) = cmplx ( 0.0E+00, 0.0E+00 )

  k = n

  do while ( 0 < k )

    if ( ipvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    kp = abs ( ipvt(k) )
    kps = k + 1 - ks

    if ( kp /= kps ) then
      t      = z(kps)
      z(kps) = z(kp)
      z(kp)  = t
    end if

    if ( cabs1 ( z(k) ) /= 0.0E+00 ) then
      ek = csign1 ( ek, z(k) )
    end if

    z(k) = z(k) + ek
    z(1:k-ks) = z(1:k-ks) + z(k) * a(1:k-ks,k)

    if ( ks /= 1 ) then

      if ( cabs1 ( z(k-1) ) /= 0.0E+00 ) then
        ek = csign1 ( ek, z(k-1) )
      end if

      z(k-1) = z(k-1) + ek
      z(1:k-ks) = z(1:k-ks) + z(k-1) * a(1:k-ks,k-1)

    end if

    if ( ks /= 2 ) then

      if ( cabs1 ( a(k,k) ) < cabs1 ( z(k) ) ) then
        s = cabs1 ( a(k,k) ) / cabs1 ( z(k) )
        z(1:n) = z(1:n) * s
        ek = cmplx ( s, 0.0E+00 ) * ek
      end if

      if ( cabs1 ( a(k,k) ) /= 0.0E+00 ) then
        z(k) = z(k) / a(k,k)
      else
        z(k) = cmplx ( 1.0E+00, 0.0E+00 )
      end if

    else

      ak = a(k,k) / conjg ( a(k-1,k) )
      akm1 = a(k-1,k-1) / a(k-1,k)
      bk = z(k) / conjg ( a(k-1,k) )
      bkm1 = z(k-1) / a(k-1,k)
      denom = ak * akm1 - cmplx ( 1.0E+00, 0.0E+00 )
      z(k) = ( akm1 * bk - bkm1 ) / denom
      z(k-1) = ( ak * bkm1 - bk ) / denom

    end if

    k = k - ks

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
!
!  Solve hermitian(U) * Y = W.
!
  k = 1

  do while ( k <= n )

    if ( ipvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= 1 ) then

!     z(k) = z(k) + cdotc ( k-1, a(1:k-1,k), 1, z(1:k-1), 1 )
      z(k) = z(k) + dot_product ( a(1:k-1,k), z(1:k-1) )

      if ( ks == 2 ) then
!       z(k+1) = z(k+1) + cdotc ( k-1, a(1:k-1,k+1), 1, z(1:k-1), 1 )
        z(k+1) = z(k+1) + dot_product ( a(1:k-1,k+1), z(1:k-1) )
      end if

      kp = abs ( ipvt(k) )

      if ( kp /= k ) then
        t     = z(k)
        z(k)  = z(kp)
        z(kp) = t
      end if

    end if

    k = k + ks

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = 1.0E+00
!
!  Solve U*D*V = Y.
!
  k = n

  do while ( 0 < k )

    if ( ipvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= ks ) then

      kp = abs ( ipvt(k) )
      kps = k + 1 - ks

      if ( kp /= kps ) then
        t      = z(kps)
        z(kps) = z(kp)
        z(kp)  = t
      end if

      z(1:k-ks) = z(1:k-ks) + z(k) * a(1:k-ks,k)

      if ( ks == 2 ) then
        z(1:k-ks) = z(1:k-ks) + z(k-1) * a(1:k-ks,k-1)
      end if

    end if

    if ( ks /= 2 ) then

      if ( cabs1 ( a(k,k) ) < cabs1 ( z(k) ) ) then
        s = cabs1 ( a(k,k) ) / cabs1 ( z(k) )
        z(1:n) = z(1:n) * s
        ynorm = s * ynorm
      end if

      if ( cabs1 ( a(k,k) ) /= 0.0E+00 ) then
        z(k) = z(k) / a(k,k)
      else
        z(k) = cmplx ( 1.0E+00, 0.0E+00 )
      end if

    else

      ak = a(k,k) / conjg ( a(k-1,k) )
      akm1 = a(k-1,k-1) / a(k-1,k)
      bk = z(k) / conjg ( a(k-1,k) )
      bkm1 = z(k-1) / a(k-1,k)
      denom = ak * akm1 - 1.0E+00
      z(k) = ( akm1 * bk - bkm1 ) / denom
      z(k-1) = ( ak * bkm1 - bk ) / denom

    end if

    k = k - ks

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm
!
!  Solve hermitian(U) * Z = V.
!
  k = 1

  do while ( k <= n )

    if ( ipvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= 1 ) then

!     z(k) = z(k) + cdotc ( k-1, a(1:k-1,k), 1, z(1:k-1), 1 )
      z(k) = z(k) + dot_product ( a(1:k-1,k), z(1:k-1) )

      if ( ks == 2 ) then
!       z(k+1) = z(k+1) + cdotc ( k-1, a(1:k-1,k+1), 1, z(1:k-1), 1 )
        z(k+1) = z(k+1) + dot_product ( a(1:k-1,k+1), z(1:k-1) )
      end if

      kp = abs ( ipvt(k) )

      if ( kp /= k ) then
        t     = z(k)
        z(k)  = z(kp)
        z(kp) = t
      end if

    end if

    k = k + ks

  end do
!
!  Make ZNORM = 1.
!
  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm

  if ( anorm /= 0.0E+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0E+00
  end if

  return
end
subroutine chidi ( a, lda, n, ipvt, det, inert, work, job )

!*****************************************************************************80
!
!! CHIDI computes the determinant and inverse of a matrix factored by CHIFA.
!
!  Discussion:
!
!    CHIDI computes the determinant, inertia (number of positive, zero,
!    and negative eigenvalues) and inverse of a complex hermitian matrix
!    using the factors from CHIFA.
!
!    A division by zero may occur if the inverse is requested
!    and CHICO has set RCOND == 0.0 or CHIFA has set INFO /= 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) A(LDA,N); on input, the factored matrix
!    from CHIFA.  On output, if the inverse was requested, A contains
!    the inverse matrix.  The strict lower triangle of A is never
!    referenced.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from CHIFA.
!
!    Output, real ( kind = 4 ) DET(2), the determinant of the original matrix.
!    Determinant = DET(1) * 10.0**DET(2) with 1.0 <= abs ( DET(1) ) < 10.0
!    or DET(1) = 0.0.
!
!    Output, integer ( kind = 4 ) INERT(3), the inertia of the original matrix.
!    INERT(1) = number of positive eigenvalues.
!    INERT(2) = number of negative eigenvalues.
!    INERT(3) = number of zero eigenvalues.
!
!    Workspace, complex WORK(N),
!
!    Input, integer ( kind = 4 ) JOB, has the decimal expansion ABC where:
!    if C /= 0, the inverse is computed,
!    if B /= 0, the determinant is computed,
!    if A /= 0, the inertia is computed.
!    For example, JOB = 111 gives all three.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  real    ( kind = 4 ) ak
  complex ( kind = 4 ) akkp1
  real    ( kind = 4 ) akp1
  complex ( kind = 4 ) cdotc
  real    ( kind = 4 ) d
  real    ( kind = 4 ) det(2)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) inert(3)
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jb
  integer ( kind = 4 ) k
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) ks
  integer ( kind = 4 ) kstep
  logical              nodet
  logical              noert
  logical              noinv
  real    ( kind = 4 ) t
  complex ( kind = 4 ) temp
  complex ( kind = 4 ) work(n)

  noinv = mod ( job,   10 )       == 0
  nodet = mod ( job,  100 ) /  10 == 0
  noert = mod ( job, 1000 ) / 100 == 0

  if ( .not. nodet .or. .not. noert ) then

    if ( .not. noert ) then
      inert(1:3) = 0
    end if

    if ( .not. nodet ) then
      det(1) = 1.0E+00
      det(2) = 0.0E+00
    end if

    t = 0.0E+00

    do k = 1, n

      d = real ( a(k,k) )
!
!  Check if 1 by 1.
!
      if ( ipvt(k) <= 0 ) then
!
!  2 by 2 block
!  Use DET = ( D / T * C - T ) * T, T = abs ( S )
!  to avoid underflow/overflow troubles.
!  Take two passes through scaling.  Use T for flag.
!
        if ( t == 0.0E+00 ) then
          t = abs ( a(k,k+1) )
          d = ( d / t ) * real ( a(k+1,k+1) ) - t
        else
          d = t
          t = 0.0E+00
        end if

      end if

      if ( .not. noert ) then
        if ( 0.0E+00 < d ) then
          inert(1) = inert(1) + 1
        else if ( d < 0.0E+00 ) then
          inert(2) = inert(2) + 1
        else if ( d == 0.0E+00 ) then
          inert(3) = inert(3) + 1
        end if
      end if

      if ( .not. nodet ) then

        det(1) = det(1) * d

        if ( det(1) /= 0.0E+00 ) then

          do while ( abs ( det(1) ) < 1.0E+00 )
            det(1) = det(1) * 10.0E+00
            det(2) = det(2) - 1.0E+00
          end do

          do while ( 10.0E+00 <= abs ( det(1) ) )
            det(1) = det(1) / 10.0E+00
            det(2) = det(2) + 1.0E+00
          end do

        end if

      end if

    end do

  end if
!
!  Compute inverse(A).
!
  if ( .not. noinv ) then

    k = 1

    do while ( k <= n )

      km1 = k - 1

      if ( 0 <= ipvt(k) ) then
!
!  1 by 1
!
        a(k,k) = cmplx ( 1.0E+00 / real ( a(k,k) ), 0.0E+00 )

        if ( 1 <= km1 ) then

          work(1:km1) = a(1:km1,k)

          do j = 1, km1
            a(j,k) = cdotc ( j, a(1:j,j), 1, work, 1 )
            a(1:j-1,k) = a(1:j-1,k) + work(j) * a(1:j-1,j)
          end do

          a(k,k) = a(k,k) + cmplx ( &
            real ( cdotc ( km1, work, 1, a(1:km1,k), 1 ) ), 0.0E+00 )

        end if

        kstep = 1

      else
!
!  2 by 2
!
        t = abs ( a(k,k+1) )
        ak = real ( a(k,k) ) / t
        akp1 = real ( a(k+1,k+1) ) / t
        akkp1 = a(k,k+1) / t
        d = t * ( ak * akp1 - 1.0E+00 )
        a(k,k) = cmplx ( akp1 / d, 0.0E+00 )
        a(k+1,k+1) = cmplx ( ak / d, 0.0E+00 )
        a(k,k+1) = -akkp1 / d

        if ( 1 <= km1 ) then

          work(1:km1) = a(1:km1,k+1)

          do j = 1, km1
            a(j,k+1) = cdotc ( j, a(1:j,j), 1, work, 1 )
            a(1:j-1,k+1) = a(1:j-1,k+1) + work(j) * a(1:j-1,j)
          end do

          a(k+1,k+1) = a(k+1,k+1) + cmplx ( &
            real ( cdotc ( km1, work, 1, a(1:km1,k+1), 1 ) ), 0.0E+00 )

          a(k,k+1) = a(k,k+1) + cdotc ( km1, a(1:km1,k), 1, a(1:km1,k+1), 1 )

          work(1:km1) = a(1:km1,k)

          do j = 1, km1
            a(j,k) = cdotc ( j, a(1:j,j), 1, work, 1 )
            a(1:j-1,k) = a(1:j-1,k) + work(j) * a(1:j-1,j)
          end do

          a(k,k) = a(k,k) + cmplx ( &
            real ( cdotc ( km1, work, 1, a(1:km1,k), 1 ) ), 0.0E+00 )

        end if

        kstep = 2

      end if
!
!  Swap
!
      ks = abs ( ipvt(k) )

      if ( ks /= k ) then

        call cswap ( ks, a(1:ks,ks), 1, a(1:ks,k), 1 )

        do j = k, ks, -1
          temp    = conjg ( a(j,k) )
          a(j,k)  = conjg ( a(ks,j) )
          a(ks,j) = temp
        end do

        if ( kstep /= 1 ) then
          temp      = a(ks,k+1)
          a(ks,k+1) = a(k,k+1)
          a(k,k+1)  = temp
        end if

      end if

      k = k + kstep

    end do

  end if

  return
end
subroutine chifa ( a, lda, n, ipvt, info )

!*****************************************************************************80
!
!! CHIFA factors a complex hermitian matrix.
!
!  Discussion:
!
!    CHIFA performs the factoring by elimination with symmetric pivoting.
!
!    To solve A*X = B, follow CHIFA by CHISL.
!
!    To compute inverse(A)*C, follow CHIFA by CHISL.
!
!    To compute determinant(A), follow CHIFA by CHIDI.
!
!    To compute inertia(A), follow CHIFA by CHIDI.
!
!    To compute inverse(A), follow CHIFA by CHIDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2007
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) A(LDA,N); on input, the hermitian matrix to be
!    factored.  On output, a block diagonal matrix and the multipliers which
!    were used to obtain it.  The factorization can be written
!    A = U*D*hermitian(U) where U is a product of permutation and unit upper
!    triangular matrices, hermitian(U) is the conjugate transpose of U, and
!    D is block diagonal with 1 by 1 and 2 by 2 blocks.  Only the diagonal
!    and upper triangle are used.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO.
!    0, normal value.
!    K, if the K-th pivot block is singular.  This is not an error condition
!    for this subroutine, but it does indicate that CHISL or CHIDI may
!    divide by zero if called.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  real    ( kind = 4 ) absakk
  complex ( kind = 4 ) ak
  complex ( kind = 4 ) akm1
  real    ( kind = 4 ) alpha
  complex ( kind = 4 ) bk
  complex ( kind = 4 ) bkm1
  real    ( kind = 4 ) cabs1
  real    ( kind = 4 ) colmax
  complex ( kind = 4 ) denom
  integer ( kind = 4 ) icamax
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) k
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) km2
  integer ( kind = 4 ) kstep
  complex ( kind = 4 ) mulk
  complex ( kind = 4 ) mulkm1
  real    ( kind = 4 ) rowmax
  logical              swap
  complex ( kind = 4 ) t
!
!  Initialize.
!
!  ALPHA is used in choosing pivot block size.
!
  alpha = ( 1.0E+00 + sqrt ( 17.0E+00 ) ) / 8.0E+00

  info = 0
!
!  Main loop on K, which goes from N to 1.
!
  k = n

  do
!
!  Leave the loop if K = 0 or K = 1.
!
    if ( k == 0 ) then
      exit
    end if

    if ( k == 1 ) then
      ipvt(1) = 1
      if ( cabs1 ( a(1,1) ) == 0.0E+00 ) then
        info = 1
      end if
      exit
    end if
!
!  This section of code determines the kind of
!  elimination to be performed.  When it is completed,
!  KSTEP will be set to the size of the pivot block, and
!  SWAP will be set to TRUE if an interchange is
!  required.
!
    km1 = k - 1
    absakk = cabs1 ( a(k,k) )
!
!  Determine the largest off-diagonal element in column K.
!
    imax = icamax ( k-1, a(1:k-1,k), 1 )
    colmax = cabs1 ( a(imax,k) )

    if ( alpha * colmax <= absakk ) then

      kstep = 1
      swap = .false.

    else
!
!  Determine the largest off-diagonal element in row IMAX.
!
      rowmax = 0.0E+00

      do j = imax + 1, k
        rowmax = max ( rowmax, cabs1 ( a(imax,j) ) )
      end do

      if ( imax /= 1 ) then
        jmax = icamax ( imax-1, a(1:imax-1,imax), 1 )
        rowmax = max ( rowmax, cabs1 ( a(jmax,imax) ) )
      end if

      if ( alpha * rowmax <= cabs1 ( a(imax,imax) )  ) then
        kstep = 1
        swap = .true.
      else if ( alpha * colmax * ( colmax / rowmax ) <= absakk ) then
        kstep = 1
        swap = .false.
      else
        kstep = 2
        swap = ( imax /= km1 )
      end if

    end if
!
!  Column K is zero.  Set INFO and iterate the loop.
!
    if ( max ( absakk, colmax ) == 0.0E+00 ) then
      ipvt(k) = k
      info = k
      k = k - kstep
      cycle
    end if

    if ( kstep /= 2 ) then
!
!  1 x 1 pivot block.
!
      if ( swap ) then

        call cswap ( imax, a(1:imax,imax), 1, a(1:imax,k), 1 )

        do jj = imax, k
          j = k + imax - jj
          t         = conjg ( a(j,k) )
          a(j,k)    = conjg ( a(imax,j) )
          a(imax,j) = t
        end do

      end if
!
!  Perform the elimination.
!
      do jj = 1, km1
        j = k - jj
        mulk = -a(j,k) / a(k,k)
        t = conjg ( mulk )
        a(1:j,j) = a(1:j,j) + t * a(1:j,k)
        a(j,j) = cmplx ( real ( a(j,j) ), 0.0E+00 )
        a(j,k) = mulk
      end do
!
!  Set the pivot array.
!
      ipvt(k) = k

      if ( swap ) then
        ipvt(k) = imax
      end if

    else
!
!  2 x 2 pivot block.
!
      if ( swap ) then

        call cswap ( imax, a(1:imax,imax), 1, a(1:imax,k-1), 1 )

        do jj = imax, km1
          j = km1 + imax - jj
          t         = conjg ( a(j,k-1) )
          a(j,k-1)  = conjg ( a(imax,j) )
          a(imax,j) = t
        end do

        t         = a(k-1,k)
        a(k-1,k)  = a(imax,k)
        a(imax,k) = t

      end if
!
!  Perform the elimination.
!
      km2 = k - 2

      if ( 0 < k - 2 ) then

        ak = a(k,k) / a(k-1,k)
        akm1 = a(k-1,k-1) / conjg ( a(k-1,k) )
        denom = 1.0E+00 - ak * akm1

        do jj = 1, k - 2

          j = km1 - jj
          bk = a(j,k) / a(k-1,k)
          bkm1 = a(j,k-1) / conjg ( a(k-1,k) )
          mulk = ( akm1 * bk - bkm1 ) / denom
          mulkm1 = ( ak * bkm1 - bk ) / denom
          t = conjg ( mulk )
          a(1:j,j) = a(1:j,j) + t * a(1:j,k)
          t = conjg ( mulkm1 )
          a(1:j,j) = a(1:j,j) + t * a(1:j,k-1)
          a(j,k) = mulk
          a(j,k-1) = mulkm1
          a(j,j) = cmplx ( real ( a(j,j) ), 0.0E+00 )

        end do

      end if
!
!  Set the pivot array.
!
      if ( swap ) then
        ipvt(k) = -imax
      else
        ipvt(k) = 1 - k
      end if

      ipvt(k-1) = ipvt(k)

    end if

    k = k - kstep

  end do

  return
end
subroutine chisl ( a, lda, n, ipvt, b )

!*****************************************************************************80
!
!! CHISL solves a complex hermitian system factored by CHIFA.
!
!  Discussion:
!
!    A division by zero may occur if CHICO has set RCOND == 0.0
!    or CHIFA has set INFO /= 0.
!
!    To compute inverse(A) * C where C is a matrix with P columns
!
!      call chifa(a,lda,n,ipvt,info)
!      if ( info == 0 ) then
!        do j = 1, p
!          call chisl(a,lda,n,ipvt,c(1,j))
!        end do
!      end if
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input, complex ( kind = 4 ) A(LDA,N), the output from CHIFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from CHIFA.
!
!    Input/output, complex ( kind = 4 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  complex ( kind = 4 ) ak
  complex ( kind = 4 ) akm1
  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) bk
  complex ( kind = 4 ) bkm1
  complex ( kind = 4 ) cdotc
  complex ( kind = 4 ) denom
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp
  complex ( kind = 4 ) t
!
!  Loop backward applying the transformations and D inverse to B.
!
  k = n

  do while ( 0 < k )
!
!  1 x 1 pivot block.
!
    if ( 0 <= ipvt(k) ) then

      if ( k /= 1 ) then

        kp = ipvt(k)

        if ( kp /= k ) then
          t     = b(k)
          b(k)  = b(kp)
          b(kp) = t
        end if

        b(1:k-1) = b(1:k-1) + a(1:k-1,k) * b(k)

      end if
!
!  Apply D inverse.
!
      b(k) = b(k) / a(k,k)
      k = k - 1
!
!  2 x 2 pivot block.
!
    else

      if ( k /= 2 ) then

        kp = abs ( ipvt(k) )

        if ( kp /= k - 1 ) then
          t      = b(k-1)
          b(k-1) = b(kp)
          b(kp)  = t
        end if

        b(1:k-2) = b(1:k-2) + b(k)   * a(1:k-2,k)
        b(1:k-2) = b(1:k-2) + b(k-1) * a(1:k-2,k-1)

      end if
!
!  Apply D inverse.
!
      ak = a(k,k) / conjg ( a(k-1,k) )
      akm1 = a(k-1,k-1) / a(k-1,k)
      bk = b(k) / conjg ( a(k-1,k) )
      bkm1 = b(k-1) / a(k-1,k)
      denom = ak * akm1 - 1.0E+00
      b(k) = ( akm1 * bk - bkm1 ) / denom
      b(k-1) = ( ak * bkm1 - bk ) / denom
      k = k - 2

    end if

  end do
!
!  Loop forward applying the transformations.
!
  k = 1
  do while ( k <= n )
!
!  1 x 1 pivot block.
!
    if ( 0 <= ipvt(k) ) then

      if ( k /= 1 ) then

        b(k) = b(k) + cdotc ( k-1, a(1:k-1,k), 1, b(1:k-1), 1 )
        kp = ipvt(k)

        if ( kp /= k ) then
          t     = b(k)
          b(k)  = b(kp)
          b(kp) = t
        end if

      end if

      k = k + 1

    else
!
!  2 x 2 pivot block.
!
      if ( k /= 1 ) then

        b(k) = b(k) + cdotc ( k-1, a(1:k-1,k), 1, b(1:k-1), 1 )
        b(k+1) = b(k+1) + cdotc ( k-1, a(1:k-1,k+1), 1, b(1:k-1), 1 )
        kp = abs ( ipvt(k) )

        if ( kp /= k ) then
          t     = b(k)
          b(k)  = b(kp)
          b(kp) = t
        end if

      end if

      k = k + 2

    end if

  end do

  return
end
subroutine chpco ( ap, n, ipvt, rcond, z )

!*****************************************************************************80
!
!! CHPCO factors a complex hermitian packed matrix and estimates its condition.
!
!  Discussion:
!
!    If RCOND is not needed, CHPFA is slightly faster.
!
!    To solve A*X = B, follow CHPCO by CHPSL.
!
!    To compute inverse(A)*C, follow CHPCO by CHPSL.
!
!    To compute inverse(A), follow CHPCO by CHPDI.
!
!    To compute determinant(A), follow CHPCO by CHPDI.
!
!    To compute inertia(A), follow CHPCO by CHPDI.
!
!  Packed storage:
!
!    The following program segment will pack the upper
!    triangle of a hermitian matrix.
!
!      k = 0
!      do j = 1, n
!        do i = 1, j
!          k = k + 1
!          ap(k) = a(i,j)
!        end do
!      end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) AP(N*(N+1)/2); on input, the packed form of a
!    hermitian matrix A.  The columns of the upper triangle are stored
!    sequentially in a one-dimensional array of length N*(N+1)/2.  On
!    output, a block diagonal matrix and the multipliers which were used
!    to obtain it stored in packed form.  The factorization can be written
!    A = U*D*hermitian(U) where U is a product of permutation and unit
!    upper triangular matrices, hermitian(U) is the conjugate transpose
!    of U, and D is block diagonal with 1 by 1 and 2 by 2 blocks.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, real ( kind = 4 ) RCOND, an estimate of the reciprocal condition of
!    the matrix.  For the system A*X = B, relative perturbations in A and B
!    of size EPSILON may cause relative perturbations in X of size
!    (EPSILON/RCOND).  If RCOND is so small that the logical expression
!      1.0 + RCOND == 1.0
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate underflows.
!
!    Workspace, complex Z(N), a work vector whose contents are usually
!    unimportant.  If A is close to a singular matrix, then Z is an
!    approximate null vector in the sense that
!      norm(A*Z) = RCOND * norm(A) * norm(Z).
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) ak
  complex ( kind = 4 ) akm1
  real    ( kind = 4 ) anorm
  complex ( kind = 4 ) ap((n*(n+1))/2)
  complex ( kind = 4 ) bk
  complex ( kind = 4 ) bkm1
  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) cdotc
  complex ( kind = 4 ) csign1
  complex ( kind = 4 ) denom
  complex ( kind = 4 ) ek
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ikm1
  integer ( kind = 4 ) ikp1
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) km1k
  integer ( kind = 4 ) km1km1
  integer ( kind = 4 ) kp
  integer ( kind = 4 ) kps
  integer ( kind = 4 ) ks
  real    ( kind = 4 ) rcond
  real    ( kind = 4 ) s
  real    ( kind = 4 ) scasum
  complex ( kind = 4 ) t
  real    ( kind = 4 ) ynorm
  complex ( kind = 4 ) z(n)
!
!  Find norm of A using only upper half.
!
  j1 = 1

  do j = 1, n

    z(j) = cmplx ( scasum ( j, ap(j1:j1+j-1), 1 ), 0.0E+00 )
    ij = j1
    j1 = j1 + j

    do i = 1, j-1
      z(i) = cmplx ( real ( z(i) ) + cabs1 ( ap(ij) ), 0.0E+00 )
      ij = ij + 1
    end do

  end do

  anorm = 0.0E+00
  do j = 1, n
    anorm = max ( anorm, real ( z(j) ) )
  end do
!
!  Factor.
!
  call chpfa ( ap, n, ipvt, info )
!
!  RCOND = 1/(norm(A) * (estimate of norm(inverse(A)))).
!
!  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
!
!  The components of E are chosen to cause maximum local
!  growth in the elements of W where U*D*W = E.
!
!  The vectors are frequently rescaled to avoid overflow.
!
!  Solve U*D*W = E.
!
  ek = cmplx ( 1.0E+00, 0.0E+00 )
  z(1:n) = cmplx ( 0.0E+00, 0.0E+00 )

  k = n
  ik = ( n * ( n - 1 ) ) / 2

  do while ( 0 < k )

    kk = ik + k
    ikm1 = ik - ( k - 1 )

    if ( ipvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    kp = abs ( ipvt(k) )
    kps = k + 1 - ks

    if ( kp /= kps ) then
      t      = z(kps)
      z(kps) = z(kp)
      z(kp)  = t
    end if

    if ( cabs1 ( z(k) ) /= 0.0E+00 ) then
      ek = csign1 ( ek, z(k) )
    end if

    z(k) = z(k) + ek
    z(1:k-ks) = z(1:k-ks) + z(k) * ap(ik+1:ik+k-ks)

    if ( ks /= 1 ) then

      if ( cabs1 ( z(k-1) ) /= 0.0E+00 ) then
        ek = csign1 ( ek, z(k-1) )
      end if

      z(k-1) = z(k-1) + ek
      z(1:k-ks) = z(1:k-ks) + z(k-1) * ap(ikm1+1:ikm1+k-ks)

    end if

    if ( ks /= 2 ) then

      if ( cabs1 ( ap(kk) ) < cabs1 ( z(k) ) ) then
        s = cabs1 ( ap(kk) ) / cabs1 ( z(k) )
        z(1:n) = z(1:n) * s
        ek = cmplx ( s, 0.0E+00 ) * ek
      end if

      if ( cabs1 ( ap(kk) ) /= 0.0E+00 ) then
        z(k) = z(k) / ap(kk)
      else
        z(k) = cmplx ( 1.0E+00, 0.0E+00 )
      end if

    else

      km1k = ik + k - 1
      km1km1 = ikm1 + k - 1
      ak = ap(kk) / conjg ( ap(km1k) )
      akm1 = ap(km1km1) / ap(km1k)
      bk = z(k) / conjg ( ap(km1k) )
      bkm1 = z(k-1) / ap(km1k)
      denom = ak * akm1 - 1.0E+00
      z(k) = ( akm1 * bk - bkm1 ) / denom
      z(k-1) = ( ak * bkm1 - bk ) / denom

    end if

    k = k - ks
    ik = ik - k

    if ( ks == 2 ) then
      ik = ik - ( k + 1 )
    end if

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
!
!  Solve hermitian(U) * Y = W.
!
  k = 1
  ik = 0

  do while ( k <= n )

    if ( ipvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= 1 ) then

      z(k) = z(k) + cdotc ( k-1, ap(ik+1:ik+k-1), 1, z(1:k-1), 1 )
      ikp1 = ik + k

      if ( ks == 2 ) then
        z(k+1) = z(k+1) + cdotc ( k-1, ap(ikp1+1:ikp1+k-1), 1, z(1:k-1), 1 )
      end if

      kp = abs ( ipvt(k) )

      if ( kp /= k ) then
        t     = z(k)
        z(k)  = z(kp)
        z(kp) = t
      end if

    end if

    ik = ik + k
    if ( ks == 2 ) then
      ik = ik + ( k + 1 )
    end if

    k = k + ks

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = 1.0E+00
!
!  Solve U*D*V = Y.
!
  k = n
  ik = ( n * ( n - 1 ) ) / 2

  do while ( 0 < k )

    kk = ik + k
    ikm1 = ik - ( k - 1 )

    if ( ipvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= ks ) then

      kp = abs ( ipvt(k) )
      kps = k + 1 - ks

      if ( kp /= kps ) then
        t      = z(kps)
        z(kps) = z(kp)
        z(kp)  = t
      end if

      z(1:k-ks) = z(1:k-ks) + z(k) * ap(ik+1:ik+k-ks)

      if ( ks == 2 ) then
        z(1:k-ks) = z(1:k-ks) + z(k-1) * ap(ikm1+1:ikm1+k-ks)
      end if

    end if

    if ( ks /= 2 ) then

      if ( cabs1 ( ap(kk) ) < cabs1 ( z(k) ) ) then
        s = cabs1 ( ap(kk) ) / cabs1 ( z(k) )
        z(1:n) = z(1:n) * s
        ynorm = s * ynorm
      end if

      if ( cabs1 ( ap(kk) ) /= 0.0E+00 ) then
        z(k) = z(k) / ap(kk)
      else
        z(k) = cmplx ( 1.0E+00, 0.0E+00 )
      end if

    else

      km1k = ik + k - 1
      km1km1 = ikm1 + k - 1
      ak = ap(kk) / conjg ( ap(km1k) )
      akm1 = ap(km1km1) / ap(km1k)
      bk = z(k) / conjg ( ap(km1k) )
      bkm1 = z(k-1) / ap(km1k)
      denom = ak * akm1 - 1.0E+00
      z(k) = ( akm1 * bk - bkm1 ) / denom
      z(k-1) = ( ak * bkm1 - bk ) / denom

    end if

    k = k - ks
    ik = ik - k

    if ( ks == 2 ) then
      ik = ik - ( k + 1 )
    end if

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm
!
!  Solve hermitian(U) * Z = V.
!
  k = 1
  ik = 0

  do while ( k <= n )

    if ( ipvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= 1 ) then

      z(k) = z(k) + cdotc ( k-1, ap(ik+1:ik+k-1), 1, z(1:k-1), 1 )
      ikp1 = ik + k

      if ( ks == 2 ) then
        z(k+1) = z(k+1) + cdotc ( k-1, ap(ikp1+1:ikp1+k-1), 1, z(1:k-1), 1 )
      end if

      kp = abs ( ipvt(k) )

      if ( kp /= k ) then
        t     = z(k)
        z(k)  = z(kp)
        z(kp) = t
      end if

    end if

    ik = ik + k

    if ( ks == 2 ) then
      ik = ik + ( k + 1 )
    end if

    k = k + ks

  end do
!
!  Make ZNORM = 1.
!
  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm

  if ( anorm /= 0.0E+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0E+00
  end if

  return
end
subroutine chpdi ( ap, n, ipvt, det, inert, work, job )

!*****************************************************************************80
!
!! CHPDI: determinant, inertia and inverse of a complex hermitian matrix.
!
!  Discussion:
!
!    The routine uses the factors from CHPFA.
!
!    The matrix is stored in packed form.
!
!    A division by zero will occur if the inverse is requested and CHPCO has
!    set RCOND == 0.0 or CHPFA has set INFO /= 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) AP(N*(N+1)/2); on input, the factored
!    matrix from CHPFA.  If the inverse was requested, then on output, AP
!    contains the upper triangle of the inverse of the original matrix, stored
!    in packed form.  The columns of the upper triangle are stored sequentially
!    in a one-dimensional array.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from CHPFA.
!
!    Output, real ( kind = 4 ) DET(2), if requested, the determinant of the
!    original matrix.
!    Determinant = DET(1) * 10.0**DET(2) with 1.0 <= abs ( DET(1) ) < 10.0
!    or DET(1) = 0.0.
!
!    Output, integer ( kind = 4 ) INERT(3), if requested, the inertia of the
!     original matrix.
!    INERT(1) = number of positive eigenvalues.
!    INERT(2) = number of negative eigenvalues.
!    INERT(3) = number of zero eigenvalues.
!
!    Workspace, complex WORK(N).
!
!    Input, integer ( kind = 4 ) JOB, has the decimal expansion ABC where:
!    if C /= 0, the inverse is computed,
!    if B /= 0, the determinant is computed,
!    if A /= 0, the inertia is computed.
!    For example, JOB = 111 gives all three.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 4 ) ak
  complex ( kind = 4 ) akkp1
  real    ( kind = 4 ) akp1
  complex ( kind = 4 ) ap((n*(n+1))/2)
  complex ( kind = 4 ) cdotc
  real    ( kind = 4 ) d
  real    ( kind = 4 ) det(2)
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ikp1
  integer ( kind = 4 ) iks
  integer ( kind = 4 ) inert(3)
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jb
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) jkp1
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) kkp1
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) ks
  integer ( kind = 4 ) ksj
  integer ( kind = 4 ) kskp1
  integer ( kind = 4 ) kstep
  logical              nodet
  logical              noert
  logical              noinv
  real    ( kind = 4 ) t
  complex ( kind = 4 ) temp
  complex ( kind = 4 ) work(n)

  noinv = mod ( job, 10 ) == 0
  nodet = mod ( job, 100 ) / 10 == 0
  noert = mod ( job, 1000 ) / 100 == 0

  if ( .not. nodet .or. .not. noert ) then

    if ( .not. noert ) then
      inert(1:3) = 0
    end if

    if ( .not. nodet ) then
      det(1) = 1.0E+00
      det(2) = 0.0E+00
    end if

    t = 0.0E+00
    ik = 0

    do k = 1, n

      kk = ik + k
      d = real ( ap(kk) )
!
!  Check if 1 by 1
!
      if ( ipvt(k) <= 0 ) then
!
!  2 by 2 block
!  Use DET (D  S; S  C)  =  ( D / T * C - T ) * T, T = abs ( S )
!  to avoid underflow/overflow troubles.
!  Take two passes through scaling.  Use T for flag.
!
        if ( t == 0.0E+00 ) then
          ikp1 = ik + k
          kkp1 = ikp1 + k
          t = abs ( ap(kkp1) )
          d = ( d / t ) * real ( ap(kkp1+1) ) - t
        else
          d = t
          t = 0.0E+00
        end if

      end if

      if ( .not. noert ) then

        if ( 0.0E+00 < d ) then
          inert(1) = inert(1) + 1
        else if ( d < 0.0E+00 ) then
          inert(2) = inert(2) + 1
        else if ( d == 0.0E+00 ) then
          inert(3) = inert(3) + 1
        end if

      end if

      if ( .not. nodet ) then

        det(1) = det(1) * d

        if ( det(1) /= 0.0E+00 ) then

          do while ( abs ( det(1) ) < 1.0E+00 )
            det(1) = det(1) * 10.0E+00
            det(2) = det(2) - 1.0E+00
          end do

          do while ( 10.0E+00 <= abs ( det(1) ) )
            det(1) = det(1) / 10.0E+00
            det(2) = det(2) + 1.0E+00
          end do

        end if

      end if

      ik = ik + k

    end do

  end if
!
!  Compute inverse(A).
!
  if ( .not. noinv ) then

    k = 1
    ik = 0

    do while ( k <= n )

      km1 = k - 1
      kk = ik + k
      ikp1 = ik + k
      kkp1 = ikp1 + k
!
!  1 by 1
!
      if ( 0 <= ipvt(k) ) then

        ap(kk) = cmplx ( 1.0E+00 / real ( ap(kk) ), 0.0E+00 )

        if ( 1 <= km1 ) then

          work(1:km1) = ap(ik+1:ik+km1)

          ij = 0
          do j = 1, km1
            jk = ik + j
            ap(jk) = cdotc ( j, ap(ij+1:ij+j), 1, work, 1 )
            ap(ik+1:ik+j-1) = ap(ik+1:ik+j-1) + work(j) * ap(ij+1:ij+j-1)
            ij = ij + j
          end do

          ap(kk) = ap(kk) + cmplx ( real ( &
            cdotc ( km1, work, 1, ap(ik+1:ik+km1), 1 ) ), 0.0E+00 )

        end if

        kstep = 1
!
!  2 by 2
!
      else

        t = abs ( ap(kkp1) )
        ak = real ( ap(kk) ) / t
        akp1 = real ( ap(kkp1+1) ) / t
        akkp1 = ap(kkp1) / t
        d = t * ( ak * akp1 - 1.0E+00 )
        ap(kk) = cmplx ( akp1 / d, 0.0E+00 )
        ap(kkp1+1) = cmplx ( ak / d, 0.0E+00 )
        ap(kkp1) = -akkp1 / d

        if ( 1 <= km1 ) then

          work(1:km1) = ap(ikp1+1:ikp1+km1)

          ij = 0
          do j = 1, km1
            jkp1 = ikp1 + j
            ap(jkp1) = cdotc ( j, ap(ij+1:ij+j), 1, work, 1 )
            ap(ikp1+1:ikp1+j-1) = ap(ikp1+1:ikp1+j-1) &
              + work(j) * ap(ij+1:ij+j-1)
            ij = ij + j
          end do

          ap(kkp1+1) = ap(kkp1+1) + cmplx ( real ( cdotc ( km1, work, 1, &
            ap(ikp1+1), 1 ) ), 0.0E+00 )

          ap(kkp1) = ap(kkp1) &
            + cdotc ( km1, ap(ik+1:ik+km1), 1, ap(ikp1+1:ikp1+km1), 1 )

          work(1:km1) = ap(ik+1:ik+km1)

          ij = 0

          do j = 1, km1
            jk = ik + j
            ap(jk) = cdotc ( j, ap(ij+1:ij+j), 1, work, 1 )
            ap(ik+1:ik+j-1) = ap(ik+1:ik+j-1) + work(j) * ap(ij+1:ij+j-1)
            ij = ij + j
          end do

          ap(kk) = ap(kk) &
            + cmplx ( real ( cdotc ( km1, work, 1, ap(ik+1:ik+km1), 1 ) ), &
            0.0E+00 )

        end if

        kstep = 2

      end if
!
!  Swap
!
      ks = abs ( ipvt(k) )

      if ( ks /= k ) then

        iks = ( ks * ( ks - 1 ) ) / 2
        call cswap ( ks, ap(iks+1:iks+ks), 1, ap(ik+1:ik+ks), 1 )
        ksj = ik + ks

        do jb = ks, k

          j = k + ks - jb
          jk = ik + j

          temp    = conjg ( ap(jk) )
          ap(jk)  = conjg ( ap(ksj) )
          ap(ksj) = temp

          ksj = ksj - ( j - 1 )

        end do

        if ( kstep /= 1 ) then

          kskp1 = ikp1 + ks

          temp      = ap(kskp1)
          ap(kskp1) = ap(kkp1)
          ap(kkp1)  = temp

        end if

      end if

      ik = ik + k

      if ( kstep == 2 ) then
        ik = ik + k + 1
      end if

      k = k + kstep

    end do

  end if

  return
end
subroutine chpfa ( ap, n, ipvt, info )

!*****************************************************************************80
!
!! CHPFA factors a complex hermitian packed matrix.
!
!  Discussion:
!
!    To solve A*X = B, follow CHPFA by CHPSL.
!
!    To compute inverse(A)*C, follow CHPFA by CHPSL.
!
!    To compute determinant(A), follow CHPFA by CHPDI.
!
!    To compute inertia(A), follow CHPFA by CHPDI.
!
!    To compute inverse(A), follow CHPFA by CHPDI.
!
!  Packed storage:
!
!    The following program segment will pack the upper
!    triangle of a hermitian matrix.
!
!      k = 0
!      do j = 1, n
!        do i = 1, j
!          k = k + 1
!          ap(k) = a(i,j)
!        end do
!      end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) AP(N*(N+1)/2); on input, the packed form
!    of a hermitian matrix.  The columns of the upper triangle are
!    stored sequentially in a one-dimensional array.  On output, a
!    block diagonal matrix and the multipliers which were used to
!    obtain it stored in packed form.  The factorization can be
!    written A = U*D*hermitian(U) where U is a product of permutation
!    and unit upper triangular matrices , hermitian(U) is the
!    conjugate transpose of U, and D is block diagonal with 1 by 1
!    and 2 by 2 blocks.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO.
!    0, normal value.
!    K, if the K-th pivot block is singular.  This is not an error condition
!    for this subroutine, but it does indicate that CHPSL or CHPDI may divide
!    by zero if called.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 4 ) absakk
  complex ( kind = 4 ) ak
  complex ( kind = 4 ) akm1
  real    ( kind = 4 ) alpha
  complex ( kind = 4 ) ap((n*(n+1))/2)
  complex ( kind = 4 ) bk
  complex ( kind = 4 ) bkm1
  real    ( kind = 4 ) cabs1
  real    ( kind = 4 ) colmax
  complex ( kind = 4 ) denom
  integer ( kind = 4 ) icamax
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ijj
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ikm1
  integer ( kind = 4 ) im
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imim
  integer ( kind = 4 ) imj
  integer ( kind = 4 ) imk
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) jkm1
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) jmim
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) km1k
  integer ( kind = 4 ) km1km1
  integer ( kind = 4 ) km2
  integer ( kind = 4 ) kstep
  complex ( kind = 4 ) mulk
  complex ( kind = 4 ) mulkm1
  real    ( kind = 4 ) rowmax
  logical              swap
  complex ( kind = 4 ) t
!
!  Initialize.
!
!  ALPHA is used in choosing pivot block size.
!
  alpha = ( 1.0E+00 + sqrt ( 17.0E+00 ) ) / 8.0E+00

  info = 0
!
!  Main loop on K, which goes from N to 1.
!
  k = n
  ik = ( n * ( n - 1 ) ) / 2

  do
!
!  Leave the loop if K = 0 or K = 1.
!
    if ( k == 0 ) then
      exit
    end if

    if ( k == 1 ) then
      ipvt(1) = 1
      if ( cabs1 ( ap(1) ) == 0.0E+00 ) then
        info = 1
      end if
      exit
    end if
!
!  This section of code determines the kind of
!  elimination to be performed.  When it is completed,
!  KSTEP will be set to the size of the pivot block, and
!  SWAP will be set to TRUE if an interchange is
!  required.
!
    km1 = k - 1
    kk = ik + k
    absakk = cabs1 ( ap(kk) )
!
!  Determine the largest off-diagonal element in column K.
!
    imax = icamax ( k-1, ap(ik+1:ik+k-1), 1 )
    imk = ik + imax
    colmax = cabs1 ( ap(imk) )

    if ( alpha * colmax <= absakk ) then

      kstep = 1
      swap = .false.
!
!  Determine the largest off-diagonal element in row IMAX.
!
    else

      rowmax = 0.0E+00
      im = ( imax * ( imax - 1 ) ) / 2
      imj = im + 2 * imax

      do j = imax + 1, k
        rowmax = max ( rowmax, cabs1 ( ap(imj) ) )
        imj = imj + j
      end do

      if ( imax /= 1 ) then
        jmax = icamax ( imax-1, ap(im+1:im+imax-1), 1 )
        jmim = jmax + im
        rowmax = max ( rowmax, cabs1 ( ap(jmim) ) )
      end if

      imim = imax + im

      if ( alpha * rowmax <= cabs1 ( ap(imim) ) ) then
        kstep = 1
        swap = .true.
      else if ( alpha * colmax * ( colmax / rowmax ) <= absakk ) then
        kstep = 1
        swap = .false.
      else
        kstep = 2
        swap = ( imax /= km1 )
      end if

    end if
!
!  Column K is zero.  Set INFO and iterate the loop.
!
    if ( max ( absakk, colmax ) == 0.0E+00 ) then
      ipvt(k) = k
      info = k
      ik = ik - ( k - 1 )
      if ( kstep == 2 ) then
        ik = ik - ( k - 2 )
      end if
      k = k - kstep
      cycle
    end if

    if ( kstep /= 2 ) then
!
!  1 x 1 pivot block.
!
      if ( swap ) then

        call cswap ( imax, ap(im+1:im+imax), 1, ap(ik+1:ik+imax), 1 )
        imj = ik + imax

        do jj = imax, k

          j = k + imax - jj
          jk = ik + j

          t       = conjg ( ap(jk) )
          ap(jk)  = conjg ( ap(imj) )
          ap(imj) = t

          imj = imj - ( j - 1 )

        end do

      end if
!
!  Perform the elimination.
!
      ij = ik - ( k - 1 )
      do jj = 1, km1
        j = k - jj
        jk = ik + j
        mulk = -ap(jk) / ap(kk)
        t = conjg ( mulk )
        ap(ij+1:ij+j) = ap(ij+1:ij+j) + t * ap(ik+1:ik+j)
        ijj = ij + j
        ap(ijj) = cmplx ( real ( ap(ijj) ), 0.0E+00 )
        ap(jk) = mulk
        ij = ij - ( j - 1 )
      end do
!
!  Set the pivot array.
!
      if ( swap ) then
        ipvt(k) = imax
      else
        ipvt(k) = k
      end if
!
!  2 x 2 pivot block.
!
    else

      km1k = ik + k - 1
      ikm1 = ik - ( k - 1 )

      if ( swap ) then

        call cswap ( imax, ap(im+1:im+imax), 1, ap(ikm1+1:ikm1+imax), 1 )
        imj = ikm1 + imax

        do jj = imax, km1

          j = km1 + imax - jj
          jkm1 = ikm1 + j

          t        = conjg ( ap(jkm1) )
          ap(jkm1) = conjg ( ap(imj) )
          ap(imj)  = t

          imj = imj - ( j - 1 )
        end do

        t        = ap(km1k)
        ap(km1k) = ap(imk)
        ap(imk)  = t

      end if
!
!  Perform the elimination.
!
      km2 = k - 2

      if ( km2 /= 0 ) then

        ak = ap(kk) / ap(km1k)
        km1km1 = ikm1 + k - 1
        akm1 = ap(km1km1) / conjg ( ap(km1k) )
        denom = 1.0E+00 - ak * akm1
        ij = ik - ( k - 1 ) - ( k - 2 )

        do jj = 1, km2
          j = km1 - jj
          jk = ik + j
          bk = ap(jk) / ap(km1k)
          jkm1 = ikm1 + j
          bkm1 = ap(jkm1) / conjg ( ap(km1k) )
          mulk = ( akm1 * bk - bkm1 ) / denom
          mulkm1 = ( ak * bkm1 - bk ) / denom
          t = conjg ( mulk )
          ap(ij+1:ij+j) = ap(ij+1:ij+j) + t * ap(ik+1:ik+j)
          t = conjg ( mulkm1 )
          ap(ij+1:ij+j) = ap(ij+1:ij+j) + t * ap(ikm1+1:ikm1+j)
          ap(jk) = mulk
          ap(jkm1) = mulkm1
          ijj = ij + j
          ap(ijj) = cmplx ( real ( ap(ijj) ), 0.0E+00 )
          ij = ij - ( j - 1 )
        end do

      end if
!
!  Set the pivot array.
!
      if ( swap ) then
        ipvt(k) = -imax
      else
        ipvt(k) = 1 - k
      end if

      ipvt(k-1) = ipvt(k)

    end if

    ik = ik - ( k - 1 )
    if ( kstep == 2 ) then
      ik = ik - ( k - 2 )
    end if

    k = k - kstep

  end do

  return
end
subroutine chpsl ( ap, n, ipvt, b )

!*****************************************************************************80
!
!! CHPSL solves a complex hermitian system factored by CHPFA.
!
!  Discussion:
!
!    A division by zero may occur if CHPCO set RCOND to 0.0
!    or CHPFA set INFO nonzero.
!
!    To compute
!
!      inverse ( A ) * C
!
!    where C is a matrix with P columns
!
!      call chpfa(ap,n,ipvt,info)
!
!      if ( info == 0 ) then
!        do j = 1, p
!          call chpsl(ap,n,ipvt,c(1,j))
!        end do
!      end if
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input, complex ( kind = 4 ) AP(N*(N+1)/2), the output from CHPFA.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from CHPFA.
!
!    Input/output, complex ( kind = 4 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) ak
  complex ( kind = 4 ) akm1
  complex ( kind = 4 ) ap((n*(n+1))/2)
  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) bk
  complex ( kind = 4 ) bkm1
  complex ( kind = 4 ) cdotc
  complex ( kind = 4 ) denom
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ikm1
  integer ( kind = 4 ) ikp1
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) km1k
  integer ( kind = 4 ) km1km1
  integer ( kind = 4 ) kp
  complex ( kind = 4 ) t
!
!  Loop backward applying the transformations and inverse ( D ) to B.
!
  k = n
  ik = ( n * ( n - 1 ) ) / 2

  do while ( 0 < k )

    kk = ik + k
!
!  1 x 1 pivot block.
!
    if ( 0 <= ipvt(k) ) then

      if ( k /= 1 ) then

        kp = ipvt(k)

        if ( kp /= k ) then
          t     = b(k)
          b(k)  = b(kp)
          b(kp) = t
        end if

        b(1:k-1) = b(1:k-1) + b(k) * ap(ik+1:ik+k-1)

      end if
!
!  Apply D inverse.
!
      b(k) = b(k) / ap(kk)
      k = k - 1
      ik = ik - k

    else
!
!  2 x 2 pivot block.
!
      ikm1 = ik - ( k - 1 )

      if ( k /= 2 ) then

        kp = abs ( ipvt(k) )

        if ( kp /= k - 1 ) then

          t      = b(k-1)
          b(k-1) = b(kp)
          b(kp)  = t

        end if

        b(1:k-2) = b(1:k-2) + b(k)   * ap(ik+1:ik+k-2)
        b(1:k-2) = b(1:k-2) + b(k-1) * ap(ikm1+1:ikm1+k-2)

      end if
!
!  Apply D inverse.
!
      km1k = ik + k - 1
      kk = ik + k
      ak = ap(kk) / conjg ( ap(km1k) )
      km1km1 = ikm1 + k - 1
      akm1 = ap(km1km1) / ap(km1k)
      bk = b(k) / conjg ( ap(km1k) )
      bkm1 = b(k-1) / ap(km1k)
      denom = ak * akm1 - 1.0E+00
      b(k) = ( akm1 * bk - bkm1 ) / denom
      b(k-1) = ( ak * bkm1 - bk ) / denom
      k = k - 2
      ik = ik - ( k + 1 ) - k

    end if

  end do
!
!  Loop forward applying the transformations.
!
  k = 1
  ik = 0

  do while ( k <= n )
!
!  1 x 1 pivot block.
!
    if ( 0 <= ipvt(k) ) then

      if ( k /= 1 ) then

        b(k) = b(k) + cdotc ( k-1, ap(ik+1:ik+k-1), 1, b(1:k-1), 1 )
        kp = ipvt(k)

        if ( kp /= k ) then
          t     = b(k)
          b(k)  = b(kp)
          b(kp) = t
        end if

      end if

      ik = ik + k
      k = k + 1
!
!  2 x 2 pivot block.
!
    else

      if ( k /= 1 ) then

        b(k) = b(k) + cdotc ( k-1, ap(ik+1:ik+k-1), 1, b(1:k-1), 1 )
        ikp1 = ik + k
        b(k+1) = b(k+1) + cdotc ( k-1, ap(ikp1+1:ikp1+k-1), 1, b(1:k-1), 1 )
        kp = abs ( ipvt(k) )

        if ( kp /= k ) then
          t     = b(k)
          b(k)  = b(kp)
          b(kp) = t
        end if

      end if

      ik = ik + k + k + 1
      k = k + 2

    end if

  end do

  return
end
subroutine cpbco ( abd, lda, n, m, rcond, z, info )

!*****************************************************************************80
!
!! CPBCO factors a complex hermitian positive definite band matrix.
!
!  Discussion:
!
!    The routine also estimates the condition number of the matrix.
!
!    If RCOND is not needed, CPBFA is slightly faster.
!
!    To solve A*X = B, follow CPBCO by CPBSL.
!
!    To compute inverse(A)*C, follow CPBCO by CPBSL.
!
!    To compute determinant(A), follow CPBCO by CPBDI.
!
!  Band storage:
!
!    If A is a hermitian positive definite band matrix,
!    the following program segment will set up the input.
!
!      m = (band width above diagonal)
!      do j = 1, n
!        i1 = max ( 1, j-m )
!        do i = i1, j
!          k = i-j+m+1
!          abd(k,j) = a(i,j)
!        end do
!      end do
!
!    This uses M+1 rows of A, except for the M by M
!    upper left triangle, which is ignored.
!
!  Example:
!
!    If the original matrix is
!
!      11 12 13  0  0  0
!      12 22 23 24  0  0
!      13 23 33 34 35  0
!       0 24 34 44 45 46
!       0  0 35 45 55 56
!       0  0  0 46 56 66
!
!    then N = 6, M = 2 and ABD should contain
!
!       *  * 13 24 35 46
!       * 12 23 34 45 56
!      11 22 33 44 55 66
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) ABD(LDA,N); on input, the matrix to be factored.
!    The columns of the upper triangle are stored in the columns of ABD,
!    and the diagonals of the upper triangle are stored in the rows of ABD.
!    On output, an upper triangular matrix R, stored in band form, so that
!    A = hermitian(R) * R.  If INFO /= 0, the factorization is not complete.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ABD.
!    LDA must be at least M+1.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) M, the number of diagonals above the main diagonal.
!    0 <= M < N.
!
!    Output, real ( kind = 4 ) RCOND, an estimate of the reciprocal condition of
!    the matrix.  For the system A*X = B, relative perturbations in A and B
!    of size EPSILON may cause relative perturbations in X of size
!    (EPSILON/RCOND).  If RCOND is so small that the logical expression
!      1.0 + RCOND == 1.0
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate underflows.
!
!    Workspace, complex Z(N), a work vector whose contents are usually
!    unimportant.  If A is singular to working precision, then Z is
!    an approximate null vector in the sense that
!    norm ( A * Z ) = RCOND * norm ( A ) * norm ( Z ).
!    If INFO /= 0, Z is unchanged.
!
!    Output, integer ( kind = 4 ) INFO.
!    0, for normal return.
!    K, signals an error condition.  The leading minor of order K is not
!    positive definite.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) abd(lda,n)
  real    ( kind = 4 ) anorm
  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) cdotc
  complex ( kind = 4 ) csign1
  complex ( kind = 4 ) ek
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) lm
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mu
  real    ( kind = 4 ) rcond
  real    ( kind = 4 ) s
  real    ( kind = 4 ) scasum
  real    ( kind = 4 ) sm
  complex ( kind = 4 ) t
  complex ( kind = 4 ) wk
  complex ( kind = 4 ) wkm
  real    ( kind = 4 ) ynorm
  complex ( kind = 4 ) z(n)

  rcond = 0.0E+00
!
!  Find the norm of A.
!
  do j = 1, n

    l = min ( j, m + 1 )
    mu = max ( m + 2 - j, 1 )
    z(j) = cmplx ( scasum ( l, abd(mu:mu+l-1,j), 1 ), 0.0E+00 )
    k = j - l

    do i = mu, m
      k = k + 1
      z(k) = cmplx ( real ( z(k) ) + cabs1 ( abd(i,j) ), 0.0E+00 )
    end do

  end do

  anorm = 0.0E+00
  do j = 1, n
    anorm = max ( anorm, real ( z(j) ) )
  end do
!
!  Factor.
!
  call cpbfa ( abd, lda, n, m, info )

  if ( info /= 0 ) then
    return
  end if
!
!  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
!
!  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
!
!  The components of E are chosen to cause maximum local
!  growth in the elements of W where hermitian(R)*W = E.
!
!  The vectors are frequently rescaled to avoid overflow.
!
!  Solve hermitian(R)*W = E.
!
  ek = cmplx ( 1.0E+00, 0.0E+00 )
  z(1:n) = cmplx ( 0.0E+00, 0.0E+00 )

  do k = 1, n

    if ( cabs1 ( z(k) ) /= 0.0E+00 ) then
      ek = csign1 ( ek, -z(k) )
    end if

    if ( real ( abd(m+1,k) ) < cabs1 ( ek - z(k) ) ) then
      s = real ( abd(m+1,k) ) / cabs1 ( ek - z(k) )
      z(1:n) = z(1:n) * s
      ek = cmplx ( s, 0.0E+00 ) * ek
    end if

    wk = ek - z(k)
    wkm = - ek - z(k)
    s = cabs1 ( wk )
    sm = cabs1 ( wkm )
    wk = wk / abd(m+1,k)
    wkm = wkm / abd(m+1,k)
    j2 = min ( k + m, n )
    i = m + 1

    if ( k+1 <= j2 ) then

      do j = k+1, j2
        i = i - 1
        sm = sm + cabs1 ( z(j) + wkm * conjg ( abd(i,j) ) )
        z(j) = z(j) + wk * conjg ( abd(i,j) )
        s = s + cabs1 ( z(j) )
      end do

      if ( s < sm ) then
        t = wkm - wk
        wk = wkm
        i = m + 1
        do j = k+1, j2
          i = i - 1
          z(j) = z(j) + t * conjg ( abd(i,j) )
        end do
      end if

    end if

    z(k) = wk

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
!
!  Solve R * Y = W.
!
  do k = n, 1, -1

    if ( real ( abd(m+1,k) ) < cabs1 ( z(k) ) ) then
      s = real ( abd(m+1,k) ) / cabs1 ( z(k) )
      z(1:n) = z(1:n) * s
    end if

    z(k) = z(k) / abd(m+1,k)
    lm = min ( k - 1, m )
    la = m + 1 - lm
    lb = k - lm
    t = -z(k)
    z(lb:lb+lm-1) = z(lb:lb+lm-1) + t * abd(la:la+lm-1,k)

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = 1.0E+00
!
!  Solve hermitian(R)*V = Y.
!
  do k = 1, n

    lm = min ( k - 1, m )
    la = m + 1 - lm
    lb = k - lm
    z(k) = z(k) - cdotc ( lm, abd(la:la+lm-1,k), 1, z(lb:lb+lm-1), 1 )

    if ( real ( abd(m+1,k) ) < cabs1 ( z(k) ) ) then
      s = real ( abd(m+1,k) ) / cabs1 ( z(k) )
      z(1:n) = z(1:n) * s
      ynorm = s * ynorm
    end if

    z(k) = z(k) / abd(m+1,k)

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm
!
!  Solve R * Z = W.
!
  do k = n, 1, -1

    if ( real ( abd(m+1,k) ) < cabs1 ( z(k) ) ) then
      s = real ( abd(m+1,k) ) / cabs1 ( z(k) )
      z(1:n) = z(1:n) * s
      ynorm = s * ynorm
    end if

    z(k) = z(k) / abd(m+1,k)
    lm = min ( k - 1, m )
    la = m + 1 - lm
    lb = k - lm
    t = -z(k)
    z(lb:lb+lm-1) = z(lb:lb+lm-1) + t * abd(la:la+lm-1,k)

  end do
!
!  Make ZNORM = 1.
!
  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm

  if ( anorm /= 0.0E+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0E+00
  end if

  return
end
subroutine cpbdi ( abd, lda, n, m, det )

!*****************************************************************************80
!
!! CPBDI gets the determinant of a hermitian positive definite band matrix.
!
!  Discussion:
!
!    CPBDI uses the factors computed by CPBCO or CPBFA.
!
!    If the inverse is needed, use CPBSL N times.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input, complex ( kind = 4 ) ABD(LDA,N), the output from CPBCO or CPBFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array ABD.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) M, the number of diagonals above the main diagonal.
!
!    Output, real ( kind = 4 ) DET(2), the determinant of the original matrix in the
!    form determinant = DET(1) * 10.0**DET(2) with 1.0 <= DET(1) < 10.0
!    or DET(1) == 0.0.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) abd(lda,n)
  real    ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m

  det(1) = 1.0E+00
  det(2) = 0.0E+00

  do i = 1, n

    det(1) = det(1) * real ( abd(m+1,i) )**2

    if ( det(1) == 0.0E+00 ) then
      exit
    end if

    do while ( det(1) < 1.0E+00 )
      det(1) = det(1) * 10.0E+00
      det(2) = det(2) - 1.0E+00
    end do

    do while ( 10.0E+00 <= det(1) )
      det(1) = det(1) / 10.0E+00
      det(2) = det(2) + 1.0E+00
    end do

  end do

  return
end
subroutine cpbfa ( abd, lda, n, m, info )

!*****************************************************************************80
!
!! CPBFA factors a complex hermitian positive definite band matrix.
!
!  Discussion:
!
!    CPBFA is usually called by CPBCO, but it can be called
!    directly with a saving in time if RCOND is not needed.
!
!  Band storage:
!
!    If A is a hermitian positive definite band matrix,
!    the following program segment will set up the input.
!
!      m = (band width above diagonal)
!      do j = 1, n
!        i1 = max ( 1, j-m )
!        do i = i1, j
!          k = i-j+m+1
!          abd(k,j) = a(i,j)
!        end do
!      end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) ABD(LDA,N); on input, the matrix to be factored.
!    The columns of the upper triangle are stored in the columns of ABD
!    and the diagonals of the upper triangle are stored in the rows of ABD.
!    On output, an upper triangular matrix R, stored in band form, so that
!    A = hermitian(R)*R.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ABD.
!    LDA must be at least M+1.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) M, the number of diagonals above the main diagonal.
!    0 <= M < N.
!
!    Output, integer ( kind = 4 ) INFO.
!    0, for normal return.
!    K, if the leading minor of order K is not positive definite.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) abd(lda,n)
  complex ( kind = 4 ) cdotc
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mu
  real    ( kind = 4 ) s
  complex ( kind = 4 ) t

  info = 0

  do j = 1, n

    s = 0.0E+00
    ik = m + 1
    jk = max ( j - m, 1 )
    mu = max ( m + 2 - j, 1 )

    do k = mu, m
      t = abd(k,j) - cdotc ( k-mu, abd(ik:ik+k-mu-1,jk), 1, abd(mu:mu+k-mu-1,j), 1 )
      t = t / abd(m+1,jk)
      abd(k,j) = t
      s = s + real ( t * conjg ( t ) )
      ik = ik - 1
      jk = jk + 1
    end do

    s = real ( abd(m+1,j) ) - s

    if ( s <= 0.0E+00 .or. aimag ( abd(m+1,j) ) /= 0.0E+00 ) then
      info = j
      exit
    end if

    abd(m+1,j) = cmplx ( sqrt ( s ), 0.0E+00 )

  end do

  return
end
subroutine cpbsl ( abd, lda, n, m, b )

!*****************************************************************************80
!
!! CPBSL solves a complex hermitian positive definite band system.
!
!  Discussion:
!
!    The system matrix must have been factored by CPBCO or CPBFA.
!
!    A division by zero will occur if the input factor contains
!    a zero on the diagonal.  Technically this indicates
!    singularity but it is usually caused by improper subroutine
!    arguments.  It will not occur if the subroutines are called
!    correctly and INFO == 0.
!
!    To compute inverse(A) * C where C is a matrix with P columns:
!
!      call cpbco(abd,lda,n,rcond,z,info)
!
!      if (rcond is too small .or. info /= 0) then
!        error
!      end if
!
!      do j = 1, p
!        call cpbsl(abd,lda,n,c(1,j))
!      end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input, complex ( kind = 4 ) ABD(LDA,N), the output from CPBCO or CPBFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ABD.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) M, the number of diagonals above the main diagonal.
!
!    Input/output, complex ( kind = 4 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) abd(lda,n)
  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) cdotc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) lm
  integer ( kind = 4 ) m
  complex ( kind = 4 ) t
!
!  Solve hermitian(R) * Y = B.
!
  do k = 1, n
    lm = min ( k - 1, m )
    la = m + 1 - lm
    lb = k - lm
    t = cdotc ( lm, abd(la:la+lm-1,k), 1, b(lb:lb+lm-1), 1 )
    b(k) = ( b(k) - t ) / abd(m+1,k)
  end do
!
!  Solve R * X = Y.
!
  do k = n, 1, -1
    lm = min ( k - 1, m )
    la = m + 1 - lm
    lb = k - lm
    b(k) = b(k) / abd(m+1,k)
    t = -b(k)
    b(lb:lb+lm-1) = b(lb:lb+lm-1) + t * abd(la:la+lm-1,k)
  end do

  return
end
subroutine cpoco ( a, lda, n, rcond, z, info )

!*****************************************************************************80
!
!! CPOCO factors a complex hermitian positive definite matrix.
!
!  Discussion:
!
!    The routine also estimates the condition of the matrix.
!
!    If RCOND is not needed, CPOFA is slightly faster.
!
!    To solve A*X = B, follow CPOCO by CPOSL.
!
!    To compute inverse(A)*C, follow CPOCO by CPOSL.
!
!    To compute determinant(A), follow CPOCO by CPODI.
!
!    To compute inverse(A), follow CPOCO by CPODI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) A(LDA,N); on input, the hermitian matrix to be
!    factored.  On output, an upper triangular matrix R so that
!      A = hermitian(R)*R
!    where hermitian(R) is the conjugate transpose.  The strict lower
!    triangle is unaltered.  If INFO /= 0, the factorization is not complete.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 4 ) RCOND, an estimate of the reciprocal condition of
!    the matrix.  For the system A*X = B, relative perturbations in A and B
!    of size EPSILON may cause relative perturbations in X of size
!    (EPSILON/RCOND).  If RCOND is so small that the logical expression
!      1.0 + RCOND == 1.0
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate underflows.
!
!    Workspace, complex Z(N), a work vector whose contents are usually
!    unimportant.  If A is close to a singular matrix, then Z is an
!    approximate null vector in the sense that
!      norm(A*Z) = RCOND * norm(A) * norm(Z).
!
!    Output, integer ( kind = 4 ) INFO.
!    0, for normal return.
!    K, signals an error condition.  The leading minor of order K is not
!    positive definite.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  real    ( kind = 4 ) anorm
  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) cdotc
  complex ( kind = 4 ) csign1
  complex ( kind = 4 ) ek
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp1
  real    ( kind = 4 ) rcond
  real    ( kind = 4 ) s
  real    ( kind = 4 ) scasum
  real    ( kind = 4 ) sm
  complex ( kind = 4 ) t
  complex ( kind = 4 ) wk
  complex ( kind = 4 ) wkm
  real    ( kind = 4 ) ynorm
  complex ( kind = 4 ) z(n)

  rcond = 0.0E+00
!
!  Find norm of A using only upper half.
!
  do j = 1, n

    z(j) = cmplx ( scasum ( j, a(1:j,j), 1 ), 0.0E+00 )

    do i = 1, j - 1
      z(i) = cmplx ( real ( z(i) ) + cabs1 ( a(i,j) ), 0.0E+00 )
    end do

  end do

  anorm = 0.0E+00
  do j = 1, n
    anorm = max ( anorm, real ( z(j) ) )
  end do
!
!  Factor.
!
  call cpofa ( a, lda, n, info )

  if ( info /= 0 ) then
    return
  end if
!
!  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
!
!  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
!
!  The components of E are chosen to cause maximum local
!  growth in the elements of W where hermitian(R)*W = E.
!
!  The vectors are frequently rescaled to avoid overflow.
!
!  Solve hermitian(R)*W = E.
!
  ek = cmplx ( 1.0E+00, 0.0E+00 )
  z(1:n) = cmplx ( 0.0E+00, 0.0E+00 )

  do k = 1, n

    if ( cabs1 ( z(k) ) /= 0.0E+00 ) then
      ek = csign1 ( ek, -z(k) )
    end if

    if ( real ( a(k,k) ) < cabs1 ( ek - z(k) ) ) then
      s = real ( a(k,k) ) / cabs1 ( ek - z(k) )
      z(1:n) = z(1:n) * s
      ek = cmplx ( s, 0.0E+00 ) * ek
    end if

    wk = ek - z(k)
    wkm = -ek - z(k)
    s = cabs1 ( wk )
    sm = cabs1 ( wkm )
    wk = wk / a(k,k)
    wkm = wkm / a(k,k)
    kp1 = k + 1

    if ( kp1 <= n ) then

      do j = kp1, n
        sm = sm + cabs1 ( z(j) + wkm * conjg ( a(k,j) ) )
        z(j) = z(j) + wk * conjg ( a(k,j) )
        s = s + cabs1 ( z(j) )
      end do

      if ( s < sm ) then
        t = wkm - wk
        wk = wkm
        z(kp1:n) = z(kp1:n) + t * conjg ( a(k,kp1:n) )
      end if

    end if

    z(k) = wk

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
!
!  Solve R * Y = W.
!
  do k = n, 1, -1

    if ( real ( a(k,k) ) < cabs1 ( z(k) ) ) then
      s = real ( a(k,k) ) / cabs1 ( z(k) )
      z(1:n) = z(1:n) * s
    end if

    z(k) = z(k) / a(k,k)
    t = -z(k)
    z(1:k-1) = z(1:k-1) + t * a(1:k-1,k)

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = 1.0E+00
!
!  Solve hermitian(R) * V = Y.
!
  do k = 1, n

    z(k) = z(k) - cdotc ( k-1, a(1:k-1,k), 1, z(1:k-1), 1 )

    if ( real ( a(k,k) ) < cabs1 ( z(k) ) ) then
      s = real ( a(k,k) ) / cabs1 ( z(k) )
      z(1:n) = z(1:n) * s
      ynorm = s * ynorm
    end if

    z(k) = z(k) / a(k,k)

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm
!
!  Solve R * Z = V.
!
  do k = n, 1, -1

    if ( real ( a(k,k) ) < cabs1 ( z(k) ) ) then
      s = real ( a(k,k) ) / cabs1 ( z(k) )
      z(1:n) = z(1:n) * s
      ynorm = s * ynorm
    end if

    z(k) = z(k) / a(k,k)
    t = -z(k)
    z(1:k-1) = z(1:k-1) + t * a(1:k-1,k)

  end do
!
!  Make ZNORM = 1.
!
  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm

  if ( anorm /= 0.0E+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0E+00
  end if

  return
end
subroutine cpodi ( a, lda, n, det, job )

!*****************************************************************************80
!
!! CPODI: determinant, inverse of a complex hermitian positive definite matrix.
!
!  Discussion:
!
!    The matrix is assumed to have been factored by CPOCO, CPOFA or CQRDC.
!
!    A division by zero will occur if the input factor contains
!    a zero on the diagonal and the inverse is requested.
!    It will not occur if the subroutines are called correctly
!    and if CPOCO or CPOFA has set INFO == 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) A(LDA,N); on input, the output A from CPOCO or
!    CPOFA, or the output X from CQRDC.  On output, if CPOCO or CPOFA was
!    used to factor A, then CPODI produces the upper half of inverse(A).
!    If CQRDC was used to decompose X, then CPODI produces the upper half
!    of inverse(hermitian(X)*X) where hermitian(X) is the conjugate transpose.
!    Elements of A below the diagonal are unchanged.
!    If the units digit of JOB is zero, A is unchanged.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 4 ) DET(2), if requested, the determinant of A or of
!    hermitian(X)*X.  Determinant = DET(1) * 10.0**DET(2) with
!    1.0 <= abs ( DET(1) ) < 10.0 or DET(1) = 0.0.
!
!    Input, integer ( kind = 4 ) JOB.
!    11, both determinant and inverse.
!    01, inverse only.
!    10, determinant only.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  real    ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  complex ( kind = 4 ) t
!
!  Compute determinant
!
  if ( ( job / 10 ) /= 0 ) then

    det(1) = 1.0E+00
    det(2) = 0.0E+00

    do i = 1, n

      det(1) = det(1) * real ( a(i,i) )**2

      if ( det(1) == 0.0E+00 ) then
        exit
      end if

      do while ( det(1) < 1.0E+00 )
        det(1) = det(1) * 10.0E+00
        det(2) = det(2) - 1.0E+00
      end do

      do while ( 10.0E+00 <= det(1) )
        det(1) = det(1) / 10.0E+00
        det(2) = det(2) + 1.0E+00
      end do

    end do

  end if
!
!  Compute inverse(R).
!
  if ( mod ( job, 10 ) /= 0 ) then

    do k = 1, n

      a(k,k) = cmplx ( 1.0E+00, 0.0E+00 ) / a(k,k)
      t = -a(k,k)
      a(1:k-1,k) = a(1:k-1,k) * t

      do j = k+1, n
        t = a(k,j)
        a(k,j) = cmplx ( 0.0E+00, 0.0E+00 )
        a(1:k,j) = a(1:k,j) + t * a(1:k,k)
      end do

    end do
!
!  Form inverse(R) * hermitian(inverse(R)).
!
    do j = 1, n

      do k = 1, j-1
        t = conjg ( a(k,j) )
        a(1:k,k) = a(1:k,k) + t * a(1:k,j)
      end do

      t = conjg ( a(j,j) )
      a(1:j,j) = a(1:j,j) * t

    end do

  end if

  return
end
subroutine cpofa ( a, lda, n, info )

!*****************************************************************************80
!
!! CPOFA factors a complex hermitian positive definite matrix.
!
!  Discussion:
!
!    CPOFA is usually called by CPOCO, but it can be called
!    directly with a saving in time if RCOND is not needed.
!    (time for CPOCO) = (1 + 18/N) * (time for CPOFA).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) A(LDA,N); On input, the hermitian matrix to be
!    factored.  On output, an upper triangular matrix R so that
!      A = hermitian(R)*R
!    where hermitian(R) is the conjugate transpose.  The strict lower
!    triangle is unaltered.  If INFO /= 0, the factorization is not
!    complete.  Only the diagonal and upper triangle are used.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) INFO.
!    0, for normal return.
!    K, signals an error condition.  The leading minor of order K is
!    not positive definite.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  complex ( kind = 4 ) cdotc
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 4 ) s
  complex ( kind = 4 ) t

  info = 0

  do j = 1, n

    s = 0.0E+00

    do k = 1, j-1
      t = a(k,j) - cdotc ( k-1, a(1:k-1,k), 1, a(1:k-1,j), 1 )
      t = t / a(k,k)
      a(k,j) = t
      s = s + real ( t * conjg ( t ) )
    end do

    s = real ( a(j,j) ) - s

    if ( s <= 0.0E+00 .or. aimag ( a(j,j) ) /= 0.0E+00 ) then
      info = j
      exit
    end if

    a(j,j) = cmplx ( sqrt ( s ), 0.0E+00 )

  end do

  return
end
subroutine cposl ( a, lda, n, b )

!*****************************************************************************80
!
!! CPOSL solves a complex hermitian positive definite system.
!
!  Discussion:
!
!    CPOSL uses the factors computed by CPOCO or CPOFA.
!
!    A division by zero will occur if the input factor contains
!    a zero on the diagonal.  Technically this indicates
!    singularity but it is usually caused by improper subroutine
!    arguments.  It will not occur if the subroutines are called
!    correctly and INFO == 0.
!
!    To compute inverse(A) * C where C is a matrix with  p  columns
!
!      call cpoco(a,lda,n,rcond,z,info)
!
!      if (rcond is too small or info is not zero ) then
!        error
!      end if
!
!      do j = 1, p
!        call cposl(a,lda,n,c(1,j))
!      end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input, complex ( kind = 4 ) A(LDA,N), the output from CPOCO or CPOFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, complex ( kind = 4 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) cdotc
  integer ( kind = 4 ) k
  complex ( kind = 4 ) t
!
!  Solve hermitian(R) * Y = B.
!
  do k = 1, n
    t = cdotc ( k-1, a(1:k-1,k), 1, b(1:k-1), 1 )
    b(k) = ( b(k) - t ) / a(k,k)
  end do
!
!  Solve R * X = Y.
!
  do k = n, 1, -1
    b(k) = b(k) / a(k,k)
    b(1:k-1) = b(1:k-1) - a(1:k-1,k) * b(k)
  end do

  return
end
subroutine cppco ( ap, n, rcond, z, info )

!*****************************************************************************80
!
!! CPPCO factors a complex hermitian positive definite matrix.
!
!  Discussion:
!
!    The routine also estimates the condition of the matrix.
!
!    The matrix is stored in packed form.
!
!    If RCOND is not needed, CPPFA is slightly faster.
!
!    To solve A*X = B, follow CPPCO by CPPSL.
!
!    To compute inverse(A)*C, follow CPPCO by CPPSL.
!
!    To compute determinant(A), follow CPPCO by CPPDI.
!
!    To compute inverse(A), follow CPPCO by CPPDI.
!
!  Packed storage:
!
!    The following program segment will pack the upper
!    triangle of a hermitian matrix.
!
!      k = 0
!      do j = 1, n
!        do i = 1, j
!          k = k + 1
!          ap(k) = a(i,j)
!        end do
!      end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 May 2007
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) AP(N*(N+1)/2); on input, the packed form of a
!    hermitian matrix.  The columns of the upper triangle are stored
!    sequentially in a one-dimensional array.  On output, an upper
!    triangular matrix R, stored in packed form, so that A = hermitian(R) * R.
!    If INFO /= 0 , the factorization is not complete.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 4 ) RCOND, an estimate of the reciprocal condition of
!    the matrix.  For the system A*X = B, relative perturbations in A and B
!    of size EPSILON may cause relative perturbations in X of size
!    (EPSILON/RCOND).  If RCOND is so small that the logical expression
!      1.0 + RCOND == 1.0
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate underflows.
!
!    Workspace, complex Z(N), a work vector whose contents are usually
!    unimportant.  If A is close to a singular matrix, then Z is an
!    approximate null vector in the sense that
!      norm(A*Z) = RCOND * norm(A) * norm(Z).
!
!    Output, integer ( kind = 4 ) INFO.
!    0, for normal return.
!    K, signals an error condition.  The leading minor of order K is not
!    positive definite.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 4 ) anorm
  complex ( kind = 4 ) ap((n*(n+1))/2)
  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) cdotc
  complex ( kind = 4 ) csign1
  complex ( kind = 4 ) ek
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kj
  integer ( kind = 4 ) kk
  real    ( kind = 4 ) rcond
  real    ( kind = 4 ) s
  real    ( kind = 4 ) scasum
  real    ( kind = 4 ) sm
  complex ( kind = 4 ) t
  complex ( kind = 4 ) wk
  complex ( kind = 4 ) wkm
  real    ( kind = 4 ) ynorm
  complex ( kind = 4 ) z(n)

  rcond = 0.0E+00
!
!  Find norm of A.
!
  j1 = 1

  do j = 1, n

    z(j) = cmplx ( scasum ( j, ap(j1:j1+j-1), 1 ), 0.0E+00 )
    ij = j1
    j1 = j1 + j

    do i = 1, j-1
      z(i) = cmplx ( real ( z(i) ) + cabs1 ( ap(ij) ), 0.0E+00 )
      ij = ij + 1
    end do

  end do

  anorm = maxval ( real ( z(1:n) ) )
!
!  Factor.
!
  call cppfa ( ap, n, info )

  if ( info /= 0 ) then
    return
  end if
!
!  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
!
!  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
!
!  The components of E are chosen to cause maximum local
!  growth in the elements of W where hermitian(R)*W = E.
!
!  The vectors are frequently rescaled to avoid overflow.
!
!  Solve hermitian(R)*W = E.
!
  ek = cmplx ( 1.0E+00, 0.0E+00 )
  z(1:n) = cmplx ( 0.0E+00, 0.0E+00 )

  kk = 0

  do k = 1, n

    kk = kk + k

    if ( cabs1 ( z(k) ) /= 0.0E+00 ) then
      ek = csign1 ( ek, -z(k) )
    end if

    if ( real ( ap(kk) ) < cabs1 ( ek - z(k) ) ) then
      s = real ( ap(kk) ) / cabs1 ( ek - z(k) )
      z(1:n) = z(1:n) * s
      ek = cmplx ( s, 0.0E+00 ) * ek
    end if

    wk = ek - z(k)
    wkm = -ek - z(k)
    s = cabs1 ( wk )
    sm = cabs1 ( wkm )
    wk = wk / ap(kk)
    wkm = wkm / ap(kk)
    kj = kk + k

    if ( k+1 <= n ) then

      do j = k+1, n
        sm = sm + cabs1 ( z(j) + wkm * conjg ( ap(kj) ) )
        z(j) = z(j) + wk * conjg ( ap(kj) )
        s = s + cabs1 ( z(j) )
        kj = kj + j
      end do

      if ( s < sm ) then
        t = wkm - wk
        wk = wkm
        kj = kk + k
        do j = k+1, n
          z(j) = z(j) + t * conjg ( ap(kj) )
          kj = kj + j
        end do
      end if

    end if

    z(k) = wk

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
!
!  Solve R * Y = W.
!
  do k = n, 1, -1

    if ( real ( ap(kk) ) < cabs1 ( z(k) ) ) then
      s = real ( ap(kk) ) / cabs1 ( z(k) )
      z(1:n) = z(1:n) * s
    end if

    z(k) = z(k) / ap(kk)
    kk = kk - k
    t = -z(k)
    z(1:k-1) = z(1:k-1) + t * ap(kk+1:kk+k-1)

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = 1.0E+00
!
!  Solve hermitian(R) * V = Y.
!
  do k = 1, n

    z(k) = z(k) - cdotc ( k-1, ap(kk+1:kk+k-1), 1, z(1:k-1), 1 )
    kk = kk + k

    if ( real ( ap(kk) ) < cabs1 ( z(k) ) ) then
      s = real ( ap(kk) ) / cabs1 ( z(k) )
      z(1:n) = z(1:n) * s
      ynorm = s * ynorm
    end if

    z(k) = z(k) / ap(kk)

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm
!
!  Solve R * Z = V.
!
  do k = n, 1, -1

    if ( real ( ap(kk) ) < cabs1 ( z(k) ) ) then
      s = real ( ap(kk) ) / cabs1 ( z(k) )
      z(1:n) = z(1:n) * s
      ynorm = s * ynorm
    end if

    z(k) = z(k) / ap(kk)
    kk = kk - k
    t = -z(k)
    z(1:k-1) = z(1:k-1) + t * ap(kk+1:kk+k-1)

  end do
!
!  Make ZNORM = 1.
!
  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm

  if ( anorm /= 0.0E+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0E+00
  end if

  return
end
subroutine cppdi ( ap, n, det, job )

!*****************************************************************************80
!
!! CPPDI: determinant, inverse of a complex hermitian positive definite matrix.
!
!  Discussion:
!
!    The matrix is assumed to have been factored by CPPCO or CPPFA.
!
!    A division by zero will occur if the input factor contains
!    a zero on the diagonal and the inverse is requested.
!    It will not occur if the subroutines are called correctly
!    and if CPOCO or CPOFA has set INFO == 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) A((N*(N+1))/2); on input, the output from CPPCO
!    or CPPFA.  On output, the upper triangular half of the inverse.
!    The strict lower triangle is unaltered.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 4 ) DET(2), the determinant of original matrix if requested.
!    Otherwise not referenced.  Determinant = DET(1) * 10.0**DET(2)
!    with 1.0 <= DET(1) < 10.0 or DET(1) == 0.0.
!
!    Input, integer ( kind = 4 ) JOB.
!    11, both determinant and inverse.
!    01, inverse only.
!    10, determinant only.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) ap((n*(n+1))/2)
  real    ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) kj
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) kp1
  complex ( kind = 4 ) t
!
!  Compute determinant.
!
  if ( ( job / 10 ) /= 0 ) then

    det(1) = 1.0E+00
    det(2) = 0.0E+00
    ii = 0

    do i = 1, n

      ii = ii + i
      det(1) = det(1) * real ( ap(ii) )**2

      if ( det(1) == 0.0E+00 ) then
        exit
      end if

      do while ( det(1) < 1.0E+00 )
        det(1) = det(1) * 10.0E+00
        det(2) = det(2) - 1.0E+00
      end do

      do while ( 10.0E+00 <= det(1) )
        det(1) = det(1) / 10.0E+00
        det(2) = det(2) + 1.0E+00
      end do

    end do

  end if
!
!  Compute inverse ( R ).
!
  if ( mod ( job, 10 ) /= 0 ) then

    kk = 0

    do k = 1, n

      k1 = kk + 1
      kk = kk + k
      ap(kk) = cmplx ( 1.0E+00, 0.0E+00 ) / ap(kk)
      t = -ap(kk)
      ap(k1:k1+k-2) = ap(k1:k1+k-2) * t
      kp1 = k + 1
      j1 = kk + 1
      kj = kk + k

      do j = kp1, n
        t = ap(kj)
        ap(kj) = cmplx ( 0.0E+00, 0.0E+00 )
        ap(j1:j1+k-1) = ap(j1:j1+k-1) + t * ap(k1:k1+k-1)
        j1 = j1 + j
        kj = kj + j
      end do

    end do
!
!  Form inverse ( R ) * hermitian ( inverse ( R ) ).
!
    jj = 0
    do j = 1, n
      j1 = jj + 1
      jj = jj + j
      k1 = 1
      kj = j1

      do k = 1, j-1
        t = conjg ( ap(kj) )
        ap(k1:k1+k-1) = ap(k1:k1+k-1) + t * ap(j1:j1+k-1)
        k1 = k1 + k
        kj = kj + 1
      end do

      t = conjg ( ap(jj) )
      ap(j1:j1+j-1) = ap(j1:j1+j-1) * t

    end do

  end if

  return
end
subroutine cppfa ( ap, n, info )

!*****************************************************************************80
!
!! CPPFA factors a complex hermitian positive definite packed matrix.
!
!  Discussion:
!
!    The following program segment will pack the upper triangle of a
!    hermitian matrix.
!
!      k = 0
!      do j = 1, n
!        do i = 1, j
!          k = k + 1
!          ap(k) = a(i,j)
!        end do
!      end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) AP(N*(N+1)/2); on input, the packed form
!    of a hermitian matrix A.  The columns of the upper triangle are
!    stored sequentially in a one-dimensional array.  On output, an
!    upper triangular matrix R, stored in packed form, so that
!      A = hermitian(R) * R.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) INFO.
!    0, for normal return.
!    K, if the leading minor of order K is not positive definite.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) ap((n*(n+1))/2)
  complex ( kind = 4 ) cdotc
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kj
  integer ( kind = 4 ) kk
  real    ( kind = 4 ) s
  complex ( kind = 4 ) t

  info = 0
  jj = 0

  do j = 1, n

    s = 0.0E+00
    kj = jj
    kk = 0
!
!  (Note that DOT_PRODUCT, when invoked with complex arguments, implicitly
!  conjugates one argument when forming the result.
!
    do k = 1, j-1
      kj = kj + 1
      t = ap(kj) - dot_product ( ap(kk+1:kk+k-1), ap(jj+1:jj+k-1) )
      kk = kk + k
      t = t / ap(kk)
      ap(kj) = t
      s = s + real ( t * conjg ( t ) )
    end do

    jj = jj + j
    s = real ( ap(jj) ) - s

    if ( s <= 0.0E+00 .or. aimag ( ap(jj) ) /= 0.0E+00 ) then
      info = j
      exit
    end if

    ap(jj) = cmplx ( sqrt ( s ), 0.0E+00 )

  end do

  return
end
subroutine cppsl ( ap, n, b )

!*****************************************************************************80
!
!! CPPSL solves a complex hermitian positive definite linear system.
!
!  Discussion:
!
!    The matrix is assumed to have been factored by CPPCO or CPPFA.
!
!    A division by zero will occur if the input factor contains
!    a zero on the diagonal.  Technically this indicates
!    singularity but it is usually caused by improper subroutine
!    arguments.  It will not occur if the subroutines are called
!    correctly and INFO == 0.
!
!    To compute inverse(A) * C where C is a matrix with P columns:
!
!      call cppco(ap,n,rcond,z,info)
!
!      if (rcond is too small .or. info /= 0) then
!        error
!      end if
!
!      do j = 1, p
!        call cppsl(ap,n,c(1,j))
!      end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input, complex ( kind = 4 ) AP(N*(N+1)/2), the output from CPPCO or CPPFA.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, complex ( kind = 4 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) ap((n*(n+1))/2)
  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) cdotc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  complex ( kind = 4 ) t

  kk = 0
  do k = 1, n
    t = cdotc ( k-1, ap(kk+1:kk+k-1), 1, b(1:k-1), 1 )
    kk = kk + k
    b(k) = ( b(k) - t ) / ap(kk)
  end do

  do k = n, 1, -1
    b(k) = b(k) / ap(kk)
    kk = kk - k
    t = -b(k)
    b(1:k-1) = b(1:k-1) + t * ap(kk+1:kk+k-1)
  end do

  return
end
subroutine cptsl ( n, d, e, b )

!*****************************************************************************80
!
!! CPTSL solves a Hermitian positive definite tridiagonal linear system.
!
!  Discussion;
!
!    The system does not have to be factored first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, complex ( kind = 4 ) D(N).  On input, the diagonal of the
!    matrix.  On output, this has been overwritten by other information.
!
!    Input/output, complex ( kind = 4 ) E(N).  On input, the superdiagonal entries of the
!    matrix in locations E(1:N-1).  On output, this has been overwritten
!    by other information.
!
!    Input/output, complex ( kind = 4 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) d(n)
  complex ( kind = 4 ) e(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kbm1
  integer ( kind = 4 ) ke
  integer ( kind = 4 ) kf
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) nm1
  integer ( kind = 4 ) nm1d2
  complex ( kind = 4 ) t1
  complex ( kind = 4 ) t2
!
!  Check for 1 x 1 case.
!
  if ( n == 1 ) then
    b(1) = b(1) / d(1)
    return
  end if

  nm1 = n - 1
  nm1d2 = ( n - 1 ) / 2

  if ( n /= 2 ) then

    kbm1 = n - 1
!
!  Zero top half of subdiagonal and bottom half of superdiagonal.
!
    do k = 1, nm1d2
      t1 = conjg ( e(k) ) / d(k)
      d(k+1) = d(k+1) - t1 * e(k)
      b(k+1) = b(k+1) - t1 * b(k)
      t2 = e(kbm1) / d(kbm1+1)
      d(kbm1) = d(kbm1) - t2 * conjg ( e(kbm1) )
      b(kbm1) = b(kbm1) - t2 * b(kbm1+1)
      kbm1 = kbm1 - 1
    end do

  end if

  kp1 = nm1d2 + 1
!
!  Clean up for possible 2 x 2 block at center.
!
  if ( mod ( n, 2 ) == 0 ) then
    t1 = conjg ( e(kp1) ) / d(kp1)
    d(kp1+1) = d(kp1+1) - t1 * e(kp1)
    b(kp1+1) = b(kp1+1) - t1 * b(kp1)
    kp1 = kp1 + 1
  end if
!
!  Back solve starting at the center, going towards the top and bottom.
!
  b(kp1) = b(kp1) / d(kp1)

  if ( n /= 2 ) then

    k = kp1 - 1
    ke = kp1 + nm1d2 - 1

    do kf = kp1, ke
      b(k) = ( b(k) - e(k) * b(k+1) ) / d(k)
      b(kf+1) = ( b(kf+1) - conjg ( e(kf) ) * b(kf) ) / d(kf+1)
      k = k - 1
    end do

  end if

  if ( mod ( n, 2 ) == 0 ) then
    b(1) = ( b(1) - e(1) * b(2) ) / d(1)
  end if

  return
end
subroutine cqrdc ( x, ldx, n, p, qraux, ipvt, work, job )

!*****************************************************************************80
!
!! CQRDC computes the QR factorization of an N by P complex matrix.
!
!  Discussion:
!
!    CQRDC uses Householder transformations to compute the QR factorization
!    of an N by P matrix X.  Column pivoting based on the 2-norms of the
!    reduced columns may be performed at the user's option.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) X(LDX,P); on input, the matrix whose
!    decomposition is to be computed.  On output, the upper triangle contains
!    the upper triangular matrix R of the QR factorization.  Below its diagonal,
!    X contains information from which the unitary part of the decomposition
!    can be recovered.  If pivoting has been requested, the decomposition is
!    not that of the original matrix X, but that of X with its columns
!    permuted as described by IPVT.
!
!    Input, integer ( kind = 4 ) LDX, the leading dimension of X.  N <= LDX.
!
!    Input, integer ( kind = 4 ) N, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) P, the number of columns in the matrix X.
!
!    Output, complex ( kind = 4 ) QRAUX(P), further information required to
!    recover the unitary part of the decomposition.
!
!    Input/output, integer ( kind = 4 ) IPVT(P); on input, integers that control
!    the selection of the pivot columns.  The K-th column X(K) of X is placed
!    in one of three classes according to the value of IPVT(K):
!      IPVT(K) > 0, then X(K) is an initial column.
!      IPVT(K) == 0, then X(K) is a free column.
!      IPVT(K) < 0, then X(K) is a final column.
!    Before the decomposition is computed, initial columns are moved to the
!    beginning of the array X and final columns to the end.  Both initial
!    and final columns are frozen in place during the computation and only
!    free columns are moved.  At the K-th stage of the reduction, if X(K)
!    is occupied by a free column it is interchanged with the free column
!    of largest reduced norm.
!    On output, IPVT(K) contains the index of the column of the
!    original matrix that has been interchanged into
!    the K-th column, if pivoting was requested.
!    IPVT is not referenced if JOB == 0.
!
!    Workspace, complex WORK(P).  WORK is not referenced if JOB == 0.
!
!    Input, integer ( kind = 4 ) JOB, initiates column pivoting.
!    0, no pivoting is done.
!    nonzero, pivoting is done.
!
  implicit none

  integer ( kind = 4 ) ldx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p

  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) cdotc
  complex ( kind = 4 ) csign2
  integer ( kind = 4 ) ipvt(p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) job
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lup
  integer ( kind = 4 ) maxj
  real    ( kind = 4 ) maxnrm
  logical              negj
  complex ( kind = 4 ) nrmxl
  integer ( kind = 4 ) pl
  integer ( kind = 4 ) pu
  complex ( kind = 4 ) qraux(p)
  real    ( kind = 4 ) scnrm2
  logical              swapj
  complex ( kind = 4 ) t
  real    ( kind = 4 ) tt
  complex ( kind = 4 ) work(p)
  complex ( kind = 4 ) x(ldx,p)

  pl = 1
  pu = 0

  if ( job /= 0 ) then
!
!  Pivoting has been requested.  Rearrange the columns according to IPVT.
!
    do j = 1, p

      swapj = ( 0 < ipvt(j) )
      negj = ( ipvt(j) < 0 )

      if ( negj ) then
        ipvt(j) = -j
      else
        ipvt(j) = j
      end if

      if ( swapj ) then

        if ( j /= pl ) then
          call cswap ( n, x(1:n,pl), 1, x(1:n,j), 1 )
        end if

        ipvt(j) = ipvt(pl)
        ipvt(pl) = j
        pl = pl + 1

      end if

    end do

    pu = p

    do jj = 1, p

      j = p - jj + 1

      if ( ipvt(j) < 0 ) then

        ipvt(j) = -ipvt(j)

        if ( j /= pu ) then

          call cswap ( n, x(1:n,pu), 1, x(1:n,j), 1 )

          i        = ipvt(pu)
          ipvt(pu) = ipvt(j)
          ipvt(j)  = i

        end if

        pu = pu - 1

      end if

    end do

  end if
!
!  Compute the norms of the free columns.
!
  do j = pl, pu
    qraux(j) = cmplx ( scnrm2 ( n, x(1:n,j), 1 ), 0.0E+00 )
    work(j) = qraux(j)
  end do
!
!  Perform the Householder reduction of X.
!
  lup = min ( n, p )

  do l = 1, lup
!
!  Locate the column of largest norm and bring it
!  into the pivot position.
!
    if ( pl <= l .and. l < pu ) then

      maxnrm = 0.0E+00
      maxj = l

      do j = l, pu
        if ( maxnrm < real ( qraux(j) ) ) then
          maxnrm = real ( qraux(j) )
          maxj = j
        end if
      end do

      if ( maxj /= l ) then
        call cswap ( n, x(1:n,l), 1, x(1:n,maxj), 1 )
        qraux(maxj) = qraux(l)
        work(maxj) = work(l)

        i          = ipvt(maxj)
        ipvt(maxj) = ipvt(l)
        ipvt(l)    = i

      end if

    end if

    qraux(l) = cmplx ( 0.0E+00, 0.0E+00 )

    if ( l /= n ) then
!
!  Compute the Householder transformation for column L.
!
      nrmxl = cmplx ( scnrm2 ( n-l+1, x(l:n,l), 1 ), 0.0E+00 )

      if ( cabs1 ( nrmxl ) /= 0.0E+00 ) then

        if ( cabs1 ( x(l,l) ) /= 0.0E+00 ) then
          nrmxl = csign2 ( nrmxl, x(l,l) )
        end if

        x(l:n,l) = x(l:n,l) / nrmxl
        x(l,l) = cmplx ( 1.0E+00, 0.0E+00 ) + x(l,l)
!
!  Apply the transformation to the remaining columns,
!  updating the norms.
!
        do j = l+1, p

          t = -cdotc ( n-l+1, x(l:n,l), 1, x(l:n,j), 1 ) / x(l,l)
          x(l:n,j) = x(l:n,j) + t * x(l:n,l)

          if ( j < pl .or. pu < j ) then
            cycle
          end if

          if ( cabs1 ( qraux(j) ) == 0.0E+00 ) then
            cycle
          end if

          tt = 1.0E+00 - ( abs ( x(l,j) ) / real ( qraux(j) ) )**2
          tt = max ( tt, 0.0E+00 )
          t = cmplx ( tt, 0.0E+00 )
          tt = 1.0E+00 &
            + 0.05E+00 * tt * ( real ( qraux(j) ) / real ( work(j) ) )**2

          if ( tt /= 1.0E+00 ) then
            qraux(j) = qraux(j) * sqrt ( t )
          else
            qraux(j) = cmplx ( scnrm2 ( n-l, x(l+1:n,j), 1 ), 0.0E+00 )
            work(j) = qraux(j)
          end if

        end do
!
!  Save the transformation.
!
        qraux(l) = x(l,l)
        x(l,l) = -nrmxl

      end if

    end if

  end do

  return
end
subroutine cqrsl ( x, ldx, n, k, qraux, y, qy, qty, b, rsd, xb, job, info )

!*****************************************************************************80
!
!! CQRSL solves, transforms or projects systems factored by CQRDC.
!
!  Discussion:
!
!    The routine applies the output of CQRDC to compute coordinate
!    transformations, projections, and least squares solutions.
!
!    For K <= min ( N, P ), let XK be the matrix
!
!      XK = ( X(IPVT(1)), X(IPVT(2)), ... ,X(IPVT(k)) )
!
!    formed from columnns IPVT(1), ... ,IPVT(K) of the original
!    N by P matrix X that was input to CQRDC (if no pivoting was
!    done, XK consists of the first K columns of X in their
!    original order).  CQRDC produces a factored unitary matrix Q
!    and an upper triangular matrix R such that
!
!      XK = Q * ( R )
!               ( 0 )
!
!    This information is contained in coded form in the arrays
!    X and QRAUX.
!
!    The parameters QY, QTY, B, RSD, and XB are not referenced
!    if their computation is not requested and in this case
!    can be replaced by dummy variables in the calling program.
!
!    To save storage, the user may in some cases use the same
!    array for different parameters in the calling sequence.  A
!    frequently occuring example is when one wishes to compute
!    any of B, RSD, or XB and does not need Y or QTY.  In this
!    case one may identify Y, QTY, and one of B, RSD, or XB, while
!    providing separate arrays for anything else that is to be
!    computed.  Thus the calling sequence
!
!      call cqrsl ( x, ldx, n, k, qraux, y, dum, y, b, y, dum, 110, info )
!
!    will result in the computation of B and RSD, with RSD
!    overwriting Y.  More generally, each item in the following
!    list contains groups of permissible identifications for
!    a single callinng sequence.
!
!    1. ( Y, QTY, B )   ( RSD )      ( XB )  ( QY )
!    2. ( Y, QTY, RSD ) ( B )        ( XB )  ( QY )
!    3. ( Y, QTY, XB )  ( B )        ( RSD ) ( QY )
!    4. ( Y, QY )       ( QTY, B )   ( RSD ) ( XB )
!    5. ( Y, QY )       ( QTY, RSD ) ( B )   ( XB )
!    6. ( Y, QY )       ( QTY, XB )  ( B )   ( RSD )
!
!    In any group the value returned in the array allocated to
!    the group corresponds to the last member of the group.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input, complex ( kind = 4 ) X(LDX,P), the output of CQRDC.
!
!    Input, integer ( kind = 4 ) LDX, the leading dimension of X.
!
!    Input, integer ( kind = 4 ) N, the number of rows of the matrix XK, which
!    must have the same value as N in CQRDC.
!
!    Input, integer ( kind = 4 ) K, the number of columns of the matrix XK.
!    K must not be greater than min ( N, P), where P is the same as in the
!    calling sequence to CQRDC.
!
!    Input, complex ( kind = 4 ) QRAUX(P), the auxiliary output from CQRDC.
!
!    Input, complex ( kind = 4 ) Y(N), a vector that is to be manipulated by CQRSL.
!
!    Output, complex ( kind = 4 ) QY(N), contains Q*Y, if it has been requested.
!
!    Output, complex ( kind = 4 ) QTY(N), contains hermitian(Q)*Y, if it has
!    been requested.  Here hermitian(Q) is the conjugate transpose
!    of the matrix Q.
!
!    Output, complex ( kind = 4 ) B(K), the solution of the least squares problem
!      minimize norm2 ( Y - XK * B ),
!    if it has been requested.  If pivoting was requested in CQRDC,
!    the J-th component of B will be associated with column IPVT(J)
!    of the original matrix X that was input into CQRDC.
!
!    Output, complex ( kind = 4 ) RSD(N), the least squares residual Y - XK*B,
!    if it has been requested.  RSD is also the orthogonal projection
!    of Y onto the orthogonal complement of the column space of XK.
!
!    Output, complex ( kind = 4 ) XB(N), the least squares approximation XK*N,
!    if its computation has been requested.  XB is also the orthogonal
!    projection of Y onto the column space of X.
!
!    Input, integer ( kind = 4 ) JOB, specifies what is to be computed.  JOB has
!    the decimal expansion ABCDE, meaning:
!    if A /= 0, compute QY.
!    if B, D, D, or E /= 0, compute QTY.
!    if C /= 0, compute B.
!    if D /= 0, compute RSD.
!    if E /= 0, compute XB.
!    A request to compute B, RSD, or XB automatically triggers the
!    computation of QTY, for which an array must be provided in the
!    calling sequence.
!
!    Output, integer ( kind = 4 ) INFO, is zero unless the computation of B has
!    been requested and R is exactly singular.  In this case, INFO is the
!    index of the first zero diagonal element of R and B is left unaltered.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) ldx
  integer ( kind = 4 ) n

  complex ( kind = 4 ) b(k)
  real    ( kind = 4 ) cabs1
  logical              cb
  complex ( kind = 4 ) cdotc
  logical              cqty
  logical              cqy
  logical              cr
  logical              cxb
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) job
  integer ( kind = 4 ) ju
  integer ( kind = 4 ) kp1
  complex ( kind = 4 ) qraux(*)
  complex ( kind = 4 ) qty(n)
  complex ( kind = 4 ) qy(n)
  complex ( kind = 4 ) rsd(n)
  complex ( kind = 4 ) t
  complex ( kind = 4 ) temp
  complex ( kind = 4 ) x(ldx,*)
  complex ( kind = 4 ) xb(n)
  complex ( kind = 4 ) y(n)

  info = 0
!
!  Determine what is to be computed.
!
  cqy =  (       job / 10000        /= 0 )
  cqty = ( mod ( job,  10000 )      /= 0 )
  cb =   ( mod ( job,   1000 ) /100 /= 0 )
  cr =   ( mod ( job,    100 ) / 10 /= 0 )
  cxb =  ( mod ( job,     10 )      /= 0 )

  ju = min ( k, n - 1 )
!
!  Special action when N=1.
!
  if ( ju == 0 ) then

    if ( cqy ) then
      qy(1) = y(1)
    end if

    if ( cqty ) then
      qty(1) = y(1)
    end if

    if ( cxb ) then
      xb(1) = y(1)
    end if

    if ( cb ) then
      if ( cabs1 ( x(1,1) ) == 0.0E+00 ) then
        info = 1
      else
        b(1) = y(1) / x(1,1)
      end if
    end if

    if ( cr ) then
      rsd(1) = cmplx ( 0.0E+00, 0.0E+00 )
    end if

    return

  end if
!
!  Set up to compute QY or QTY.
!
  if ( cqy ) then
    qy(1:n) = y(1:n)
  end if

  if ( cqty ) then
    qty(1:n) = y(1:n)
  end if
!
!  Compute QY.
!
  if ( cqy ) then

    do jj = 1, ju

      j = ju - jj + 1

      if ( cabs1 ( qraux(j) ) /= 0.0E+00 ) then
        temp = x(j,j)
        x(j,j) = qraux(j)
        t = -cdotc ( n-j+1, x(j:n,j), 1, qy(j:n), 1 ) / x(j,j)
        qy(j:n) = qy(j:n) + t * x(j:n,j)
        x(j,j) = temp
      end if

    end do

  end if
!
!  Compute hermitian ( A ) * Y.
!
  if ( cqty ) then
    do j = 1, ju
      if ( cabs1 ( qraux(j) ) /= 0.0E+00 ) then
        temp = x(j,j)
        x(j,j) = qraux(j)
        t = -cdotc ( n-j+1, x(j:n,j), 1, qty(j:n), 1 ) / x(j,j)
        qty(j:n) = qty(j:n) + t * x(j:n,j)
        x(j,j) = temp
      end if
    end do
  end if
!
!  Set up to compute B, RSD, or XB.
!
  if ( cb ) then
    b(1:k) = qty(1:k)
  end if

  kp1 = k + 1

  if ( cxb ) then
    xb(1:k) = qty(1:k)
  end if

  if ( cr .and. k < n ) then
    rsd(k+1:n) = qty(k+1:n)
  end if

  if ( cxb ) then
    xb(k+1:n) = cmplx ( 0.0E+00, 0.0E+00 )
  end if

  if ( cr ) then
    rsd(1:k) = cmplx ( 0.0E+00, 0.0E+00 )
  end if
!
!  Compute B.
!
  if ( cb ) then

    do jj = 1, k

      j = k - jj + 1

      if ( cabs1 ( x(j,j) ) == 0.0E+00 ) then
        info = j
        exit
      end if

      b(j) = b(j) / x(j,j)

      if ( j /= 1 ) then
        t = -b(j)
        b(1:j-1) = b(1:j-1) + t * x(1:j-1,j)
      end if

    end do

  end if

  if ( cr .or.  cxb ) then
!
!  Compute RSD or XB as required.
!
    do jj = 1, ju

      j = ju - jj + 1

      if ( cabs1 ( qraux(j) ) /= 0.0E+00 ) then

        temp = x(j,j)
        x(j,j) = qraux(j)

        if ( cr ) then
          t = -cdotc ( n-j+1, x(j:n,j), 1, rsd(j:n), 1 ) / x(j,j)
          rsd(j:n) = rsd(j:n) + t * x(j:n,j)
        end if

        if ( cxb ) then
          t = -cdotc ( n-j+1, x(j:n,j), 1, xb(j:n), 1 ) / x(j,j)
          xb(j:n) = xb(j:n) + t * x(j:n,j)
        end if

        x(j,j) = temp

      end if

    end do

  end if

  return
end
subroutine csico ( a, lda, n, ipvt, rcond, z )

!*****************************************************************************80
!
!! CSICO factors a complex symmetric matrix.
!
!  Discussion:
!
!    The factorization is done by symmetric pivoting.
!
!    The routine also estimates the condition of the matrix.
!
!    If RCOND is not needed, CSIFA is slightly faster.
!
!    To solve A*X = B, follow CSICO by CSISL.
!
!    To compute inverse(A)*C, follow CSICO by CSISL.
!
!    To compute inverse(A), follow CSICO by CSIDI.
!
!    To compute determinant(A), follow CSICO by CSIDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) A(LDA,N); on input, the symmetric matrix
!    to be factored.  On output, a block diagonal matrix and the multipliers
!    which were used to obtain it.  The factorization can be written A = U*D*U'
!    where U is a product of permutation and unit upper triangular matrices,
!    U' is the transpose of U, and D is block diagonal with 1 by 1 and
!    2 by 2 blocks.  Only the diagonal and upper triangle are used.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, real ( kind = 4 ) RCOND, an estimate of the reciprocal condition of
!    the matrix.  For the system A*X = B, relative perturbations in A and B
!    of size EPSILON may cause relative perturbations in X of size
!    (EPSILON/RCOND).  If RCOND is so small that the logical expression
!      1.0 + RCOND == 1.0
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate underflows.
!
!    Workspace, complex Z(N), a work vector whose contents are usually
!    unimportant.  If A is close to a singular matrix, then Z is an
!    approximate null vector in the sense that
!      norm(A*Z) = RCOND * norm(A) * norm(Z).
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  complex ( kind = 4 ) ak
  complex ( kind = 4 ) akm1
  real    ( kind = 4 ) anorm
  complex ( kind = 4 ) bk
  complex ( kind = 4 ) bkm1
  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) cdotu
  complex ( kind = 4 ) csign1
  complex ( kind = 4 ) denom
  complex ( kind = 4 ) ek
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp
  integer ( kind = 4 ) kps
  integer ( kind = 4 ) ks
  real    ( kind = 4 ) rcond
  real    ( kind = 4 ) s
  real    ( kind = 4 ) scasum
  complex ( kind = 4 ) t
  real    ( kind = 4 ) ynorm
  complex ( kind = 4 ) z(n)
!
!  Find norm of A using only upper half.
!
  do j = 1, n
    z(j) = cmplx ( scasum ( j, a(1:j,j), 1 ), 0.0E+00 )
    do i = 1, j-1
      z(i) = cmplx ( real ( z(i) ) + cabs1 ( a(i,j) ), 0.0E+00 )
    end do
  end do

  anorm = 0.0E+00
  do j = 1, n
    anorm = max ( anorm, real ( z(j) ) )
  end do
!
!  Factor.
!
  call csifa ( a, lda, n, ipvt, info )
!
!  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
!
!  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
!
!  The components of E are chosen to cause maximum local
!  growth in the elements of W where U*D*W = E.
!
!  The vectors are frequently rescaled to avoid overflow.
!
!  Solve U*D*W = E.
!
  ek = cmplx ( 1.0E+00, 0.0E+00 )
  z(1:n) = cmplx ( 0.0E+00, 0.0E+00 )

  k = n

  do while ( 0 < k )

    if ( ipvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    kp = abs ( ipvt(k) )
    kps = k + 1 - ks

    if ( kp /= kps ) then
      t      = z(kps)
      z(kps) = z(kp)
      z(kp)  = t
    end if

    if ( cabs1 ( z(k) ) /= 0.0E+00 ) then
      ek = csign1 ( ek, z(k) )
    end if

    z(k) = z(k) + ek
    z(1:k-ks) = z(1:k-ks) + z(k) * a(1:k-ks,k)

    if ( ks /= 1 ) then

      if ( cabs1 ( z(k-1) ) /= 0.0E+00 ) then
        ek = csign1 ( ek, z(k-1) )
      end if

      z(k-1) = z(k-1) + ek
      z(1:k-ks) = z(1:k-ks) + z(k-1) * a(1:k-ks,k-1)

    end if

    if ( ks /= 2 ) then

      if ( cabs1 ( a(k,k) ) < cabs1 ( z(k) ) ) then
        s = cabs1 ( a(k,k) ) / cabs1 ( z(k) )
        z(1:n) = z(1:n) * s
        ek = cmplx ( s, 0.0E+00 ) * ek
      end if

      if ( cabs1 ( a(k,k) ) /= 0.0E+00 ) then
        z(k) = z(k) / a(k,k)
      else
        z(k) = cmplx ( 1.0E+00, 0.0E+00 )
      end if

    else

      ak = a(k,k) / a(k-1,k)
      akm1 = a(k-1,k-1) / a(k-1,k)
      bk = z(k) / a(k-1,k)
      bkm1 = z(k-1) / a(k-1,k)
      denom = ak * akm1 - 1.0E+00
      z(k) = ( akm1 * bk - bkm1 ) / denom
      z(k-1) = ( ak * bkm1 - bk ) / denom

    end if

    k = k - ks

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
!
!  Solve U' * Y = W.
!
  k = 1

  do while ( k <= n )

    if ( ipvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= 1 ) then

      z(k) = z(k) + cdotu ( k-1, a(1:k-1,k), 1, z(1:k-1), 1 )

      if ( ks == 2 ) then
        z(k+1) = z(k+1) + cdotu ( k-1, a(1:k-1,k+1), 1, z(1:k-1), 1 )
      end if

      kp = abs ( ipvt(k) )

      if ( kp /= k ) then
        t     = z(k)
        z(k)  = z(kp)
        z(kp) = t
      end if

    end if

    k = k + ks

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = 1.0E+00
!
!  Solve U*D*V = Y.
!
  k = n

  do while ( 0 < k )

    if ( ipvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= ks ) then

      kp = abs ( ipvt(k) )
      kps = k + 1 - ks

      if ( kp /= kps ) then
        t      = z(kps)
        z(kps) = z(kp)
        z(kp)  = t
      end if

      z(1:k-ks) = z(1:k-ks) + z(k) * a(1:k-ks,k)

      if ( ks == 2 ) then
        z(1:k-ks) = z(1:k-ks) + z(k-1) * a(1:k-ks,k-1)
      end if

    end if

    if ( ks /= 2 ) then

      if ( cabs1 ( a(k,k) ) < cabs1 ( z(k) ) ) then
        s = cabs1 ( a(k,k) ) / cabs1 ( z(k) )
        z(1:n) = z(1:n) * s
        ynorm = s * ynorm
      end if

      if ( cabs1 ( a(k,k) ) /= 0.0E+00 ) then
        z(k) = z(k) / a(k,k)
      else
        z(k) = cmplx ( 1.0E+00, 0.0E+00 )
      end if

    else

      ak = a(k,k) / a(k-1,k)
      akm1 = a(k-1,k-1) / a(k-1,k)
      bk = z(k) / a(k-1,k)
      bkm1 = z(k-1) / a(k-1,k)
      denom = ak * akm1 - 1.0E+00
      z(k) = ( akm1 * bk - bkm1 ) / denom
      z(k-1) = ( ak * bkm1 - bk ) / denom

    end if

    k = k - ks

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm
!
!  Solve U' * Z = V.
!
  k = 1

  do while ( k <= n )

    if ( ipvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= 1 ) then

      z(k) = z(k) + cdotu ( k-1, a(1:k-1,k), 1, z(1:k-1), 1 )

      if ( ks == 2 ) then
        z(k+1) = z(k+1) + cdotu ( k-1, a(1:k-1,k+1), 1, z(1:k-1), 1 )
      end if

      kp = abs ( ipvt(k) )

      if ( kp /= k ) then
        t     = z(k)
        z(k)  = z(kp)
        z(kp) = t
      end if

    end if

    k = k + ks

  end do
!
!  Make ZNORM = 1.
!
  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm

  if ( anorm /= 0.0E+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0E+00
  end if

  return
end
subroutine csidi ( a, lda, n, ipvt, det, work, job )

!*****************************************************************************80
!
!! CSIDI computes the determinant and inverse of a matrix factored by CSIFA.
!
!  Discussion:
!
!    It is assumed the complex symmetric matrix has already been factored
!    by CSIFA.
!
!    A division by zero may occur if the inverse is requested
!    and CSICO set RCOND == 0.0 or CSIFA set INFO nonzero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) A(LDA,N); on input, the output from CSIFA.
!    If the inverse was requested, then on output, A contains the upper triangle
!    of the inverse of the original matrix.  The strict lower triangle
!    is never referenced.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from CSIFA.
!
!    Output, complex ( kind = 4 ) DET(2), if requested, the determinant of the matrix.
!    Determinant = DET(1) * 10.0**DET(2) with 1.0 <= abs ( DET(1) ) < 10.0
!    or DET(1) = 0.0.  Also, DET(2) is strictly real.
!
!    Workspace, complex WORK(N).
!
!    Input, integer ( kind = 4 ) JOB, has the decimal expansion AB where
!    if B /= 0, the inverse is computed,
!    if A /= 0, the determinant is computed,
!    For example, JOB = 11 gives both.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  complex ( kind = 4 ) ak
  complex ( kind = 4 ) akkp1
  complex ( kind = 4 ) akp1
  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) cdotu
  complex ( kind = 4 ) d
  complex ( kind = 4 ) det(2)
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jb
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) ks
  integer ( kind = 4 ) kstep
  logical              nodet
  logical              noinv
  complex ( kind = 4 ) t
  complex ( kind = 4 ) work(n)

  noinv = mod ( job, 10 ) == 0
  nodet = mod ( job, 100 ) / 10 == 0

  if ( .not. nodet ) then

    det(1) = cmplx ( 1.0E+00, 0.0E+00 )
    det(2) = cmplx ( 0.0E+00, 0.0E+00 )
    t = cmplx ( 0.0E+00, 0.0E+00 )

    do k = 1, n

      d = a(k,k)
!
!   2 by 2 block.
!   Use det ( D  T ) = ( D / T * C - T ) * T
!           ( T  C )
!   to avoid underflow/overflow troubles.
!   Take two passes through scaling.  Use T for flag.
!
      if ( ipvt(k) <= 0 ) then

        if ( cabs1 ( t ) == 0.0E+00 ) then
          t = a(k,k+1)
          d = ( d / t ) * a(k+1,k+1) - t
        else
          d = t
          t = cmplx ( 0.0E+00, 0.0E+00 )
        end if

      end if

      det(1) = det(1) * d

      if ( cabs1 ( det(1) ) /= 0.0E+00 ) then

        do while ( cabs1 ( det(1) ) < 1.0E+00 )
          det(1) = det(1) * cmplx ( 10.0E+00, 0.0E+00 )
          det(2) = det(2) - cmplx ( 1.0E+00, 0.0E+00 )
        end do

        do while ( 10.0E+00 <= cabs1 ( det(1) ) )
          det(1) = det(1) / cmplx ( 10.0E+00, 0.0E+00 )
          det(2) = det(2) + cmplx ( 1.0E+00, 0.0E+00 )
        end do

      end if

    end do

  end if
!
!  Compute inverse ( A ).
!
  if ( .not. noinv ) then

    k = 1

    do while ( k <= n )

      km1 = k - 1

!
!  1 by 1
!
      if ( 0 <= ipvt(k) ) then

        a(k,k) = cmplx ( 1.0E+00, 0.0E+00 ) / a(k,k)

        if ( 1 <= km1 ) then

          work(1:km1) = a(1:km1,k)

          do j = 1, km1
            a(j,k) = cdotu ( j, a(1:j,j), 1, work, 1 )
            a(1:j-1,k) = a(1:j-1,k) + work(j) * a(1:j-1,j)
          end do

          a(k,k) = a(k,k) + cdotu ( km1, work, 1, a(1:km1,k), 1 )

        end if

        kstep = 1
!
!  2 by 2
!
      else

        t = a(k,k+1)
        ak = a(k,k) / t
        akp1 = a(k+1,k+1) / t
        akkp1 = a(k,k+1) / t
        d = t * ( ak * akp1 - cmplx ( 1.0E+00, 0.0E+00 ) )
        a(k,k) = akp1 / d
        a(k+1,k+1) = ak / d
        a(k,k+1) = -akkp1 / d

        if ( 1 <= km1 ) then

          work(1:km1) = a(1:km1,k+1)

          do j = 1, km1
            a(j,k+1) = cdotu ( j, a(1:j,j), 1, work, 1 )
            a(1:j-1,k+1) = a(1:j-1,k+1) + work(j) * a(1:j-1,j)
          end do

          a(k+1,k+1) = a(k+1,k+1) + cdotu ( km1, work, 1, a(1:km1,k+1), 1 )
          a(k,k+1) = a(k,k+1) + cdotu ( km1, a(1:km1,k), 1, a(1:km1,k+1), 1 )

          work(1:km1) = a(1:km1,k)

          do j = 1, km1
            a(j,k) = cdotu ( j, a(1:j,j), 1, work, 1 )
            a(1:j-1,k) = a(1:j-1,k) + work(j) * a(1:j-1,j)
          end do

          a(k,k) = a(k,k) + cdotu ( km1, work, 1, a(1:km1,k), 1 )

        end if

        kstep = 2

      end if
!
!  Swap.
!
      ks = abs ( ipvt(k) )

      if ( ks /= k ) then

        call cswap ( ks, a(1:ks,ks), 1, a(1:ks,k), 1 )

        do jb = ks, k
          j = k + ks - jb
          t       = a(j,k)
          a(j,k)  = a(ks,j)
          a(ks,j) = t
        end do

        if ( kstep /= 1 ) then
          t         = a(ks,k+1)
          a(ks,k+1) = a(k,k+1)
          a(k,k+1)  = t
        end if

      end if

      k = k + kstep

    end do

  end if

  return
end
subroutine csifa ( a, lda, n, ipvt, info )

!*****************************************************************************80
!
!! CSIFA factors a complex symmetric matrix.
!
!  Discussion:
!
!    The factorization is accomplished by elimination with symmetric pivoting.
!
!    To solve A*X = B, follow CSIFA by CSISL.
!
!    To compute inverse(A)*C, follow CSIFA by CSISL.
!
!    To compute determinant(A), follow CSIFA by CSIDI.
!
!    To compute inverse(A), follow CSIFA by CSIDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2007
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) A(LDA,N); on input, the symmetric matrix to be
!    factored.  On output, a block diagonal matrix and the multipliers which
!    were used to obtain it.  The factorization can be written A = U*D*U'
!    where U is a product of permutation and unit upper triangular matrices,
!    U' is the transpose of U, and D is block diagonal with 1 by 1 and 2 by 2
!    blocks.  Only the diagonal and upper triangle are used.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO.
!    0, normal value.
!    K, if the K-th pivot block is singular.  This is not an error condition
!    for this subroutine, but it does indicate that CSISL or CSIDI may
!    divide by zero if called.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  real    ( kind = 4 ) absakk
  complex ( kind = 4 ) ak
  complex ( kind = 4 ) akm1
  real    ( kind = 4 ) alpha
  complex ( kind = 4 ) bk
  complex ( kind = 4 ) bkm1
  real    ( kind = 4 ) cabs1
  real    ( kind = 4 ) colmax
  complex ( kind = 4 ) denom
  integer ( kind = 4 ) icamax
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) k
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) km2
  integer ( kind = 4 ) kstep
  complex ( kind = 4 ) mulk
  complex ( kind = 4 ) mulkm1
  real    ( kind = 4 ) rowmax
  logical              swap
  complex ( kind = 4 ) t
!
!  Initialize.
!
!  ALPHA is used in choosing pivot block size.
!
  alpha = ( 1.0E+00 + sqrt ( 17.0E+00 ) ) / 8.0E+00

  info = 0
!
!  Main loop on K, which goes from N to 1.
!
  k = n

  do
!
!  Leave the loop if K = 0 or K = 1.
!
    if ( k == 0 ) then
      exit
    end if

    if ( k == 1 ) then
      ipvt(1) = 1
      if ( cabs1 ( a(1,1) ) == 0.0E+00 ) then
        info = 1
      end if
      exit
    end if
!
!  This section of code determines the kind of
!  elimination to be performed.  When it is completed,
!  KSTEP will be set to the size of the pivot block, and
!  SWAP will be set to TRUE if an interchange is
!  required.
!
    km1 = k - 1
    absakk = cabs1 ( a(k,k) )
!
!  Determine the largest off-diagonal element in column K.
!
    imax = icamax ( k-1, a(1:k-1,k), 1 )
    colmax = cabs1 ( a(imax,k) )

    if ( alpha * colmax < absakk ) then

      kstep = 1
      swap = .false.
!
!  Determine the largest off-diagonal element in row IMAX.
!
    else

      rowmax = 0.0E+00

      do j = imax + 1, k
        rowmax = max ( rowmax, cabs1 ( a(imax,j) ) )
      end do

      if ( imax /= 1 ) then
        jmax = icamax ( imax-1, a(1:imax-1,imax), 1 )
        rowmax = max ( rowmax, cabs1 ( a(jmax,imax) ) )
      end if

      if ( alpha * rowmax <= cabs1 ( a(imax,imax) ) ) then
        kstep = 1
        swap = .true.
      else if ( alpha * colmax * ( colmax / rowmax ) <= absakk ) then
        kstep = 1
        swap = .false.
      else
        kstep = 2
        swap = imax /= km1
      end if

    end if
!
!  Column K is zero.  Set INFO and iterate the loop.
!
    if ( max ( absakk, colmax ) == 0.0E+00 ) then
      ipvt(k) = k
      info = k
      k = k - kstep
      cycle
    end if

    if ( kstep /= 2 ) then
!
!  1 x 1 pivot block.
!
!  Perform an interchange.
!
      if ( swap ) then

        call cswap ( imax, a(1:imax,imax), 1, a(1:imax,k), 1 )

        do jj = imax, k
          j = k + imax - jj
          t         = a(j,k)
          a(j,k)    = a(imax,j)
          a(imax,j) = t
        end do

      end if
!
!  Perform the elimination.
!
      do jj = 1, km1
        j = k - jj
        mulk = -a(j,k) / a(k,k)
        t = mulk
        a(1:j,j) = a(1:j,j) + t * a(1:j,k)
        a(j,k) = mulk
      end do
!
!  Set the pivot array.
!
      if ( swap ) then
        ipvt(k) = imax
      else
        ipvt(k) = k
      end if
!
!  2 x 2 pivot block.
!
    else

      if ( swap ) then

        call cswap ( imax, a(1:imax,imax), 1, a(1:imax,k-1), 1 )

        do jj = imax, km1
          j = km1 + imax - jj
          t         = a(j,k-1)
          a(j,k-1)  = a(imax,j)
          a(imax,j) = t
        end do

        t         = a(k-1,k)
        a(k-1,k)  = a(imax,k)
        a(imax,k) = t

      end if
!
!  Perform the elimination.
!
      km2 = k - 2

      if ( km2 /= 0 ) then

        ak = a(k,k) / a(k-1,k)
        akm1 = a(k-1,k-1) / a(k-1,k)
        denom = 1.0E+00 - ak * akm1

        do jj = 1, km2
          j = km1 - jj
          bk = a(j,k) / a(k-1,k)
          bkm1 = a(j,k-1) / a(k-1,k)
          mulk = ( akm1 * bk - bkm1 ) / denom
          mulkm1 = ( ak * bkm1 - bk ) / denom
          t = mulk
          a(1:j,j) = a(1:j,j) + t * a(1:j,k)
          t = mulkm1
          a(1:j,j) = a(1:j,j) + t * a(1:j,k-1)
          a(j,k) = mulk
          a(j,k-1) = mulkm1
        end do

      end if
!
!  Set the pivot array.
!
      if ( swap ) then
        ipvt(k) = -imax
      else
        ipvt(k) = 1 - k
      end if

      ipvt(k-1) = ipvt(k)

    end if

    k = k - kstep

  end do

  return
end
subroutine csisl ( a, lda, n, ipvt, b )

!*****************************************************************************80
!
!! CSISL solves a complex symmetric system that was factored by CSIFA.
!
!  Discussion:
!
!    A division by zero may occur if CSICO has set RCOND == 0.0
!    or CSIFA has set INFO /= 0.
!
!    To compute inverse(A) * C where C is a matrix with P columns:
!
!      call csifa(a,lda,n,ipvt,info)
!
!      if ( info /= 0 ) then
!        error
!      end if
!
!      do j = 1, p
!        call csisl(a,lda,n,ipvt,c(1,j))
!      end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input, complex ( kind = 4 ) A(LDA,N), the output from CSICO or CSIFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from CSICO or CSIFA.
!
!    Input/output, complex ( kind = 4 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(lda,n)
  complex ( kind = 4 ) ak
  complex ( kind = 4 ) akm1
  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) bk
  complex ( kind = 4 ) bkm1
  complex ( kind = 4 ) cdotu
  complex ( kind = 4 ) denom
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp
  complex ( kind = 4 ) t
!
!  Loop backward applying the transformations and D inverse to B.
!
  k = n

  do while ( 0 < k )
!
!  1 x 1 pivot block.
!
    if ( 0 <= ipvt(k) ) then

      if ( k /= 1 ) then

        kp = ipvt(k)

        if ( kp /= k ) then
          t     = b(k)
          b(k)  = b(kp)
          b(kp) = t
        end if

        b(1:k-1) = b(1:k-1) + b(k) * a(1:k-1,k)

      end if

      b(k) = b(k) / a(k,k)
      k = k - 1
!
!  2 x 2 pivot block.
!
    else

      if ( k /= 2 ) then

        kp = abs ( ipvt(k) )

        if ( kp /= k - 1 ) then
          t      = b(k-1)
          b(k-1) = b(kp)
          b(kp)  = t
        end if

        b(1:k-2) = b(1:k-2) + b(k)   * a(1:k-2,k)
        b(1:k-2) = b(1:k-2) + b(k-1) * a(1:k-2,k-1)

      end if

      ak = a(k,k) / a(k-1,k)
      akm1 = a(k-1,k-1) / a(k-1,k)
      bk = b(k) / a(k-1,k)
      bkm1 = b(k-1) / a(k-1,k)
      denom = ak * akm1 - 1.0E+00
      b(k) = ( akm1 * bk - bkm1 ) / denom
      b(k-1) = ( ak * bkm1 - bk ) / denom
      k = k - 2

    end if

  end do
!
!  Loop forward applying the transformations.
!
  k = 1

  do while ( k <= n )

    if ( 0 <= ipvt(k) ) then
!
!  1 x 1 pivot block.
!
      if ( k /= 1 ) then

        b(k) = b(k) + cdotu ( k-1, a(1:k-1,k), 1, b(1:k-1), 1 )
        kp = ipvt(k)

        if ( kp /= k ) then
          t     = b(k)
          b(k)  = b(kp)
          b(kp) = t
        end if

      end if

      k = k + 1
!
!  2 x 2 pivot block.
!
    else

      if ( k /= 1 ) then

        b(k) = b(k) + cdotu ( k-1, a(1:k-1,k), 1, b(1:k-1), 1 )
        b(k+1) = b(k+1) + cdotu ( k-1, a(1:k-1,k+1), 1, b(1:k-1), 1 )
        kp = abs ( ipvt(k) )

        if ( kp /= k ) then
          t     = b(k)
          b(k)  = b(kp)
          b(kp) = t
        end if

      end if

      k = k + 2

    end if

  end do

  return
end
subroutine cspco ( ap, n, ipvt, rcond, z )

!*****************************************************************************80
!
!! CSPCO factors a complex symmetric matrix stored in packed form.
!
!  Discussion:
!
!    The routine also estimates the condition of the matrix.
!
!    If RCOND is not needed, CSPFA is slightly faster.
!
!    To solve A*X = B, follow CSPCO by CSPSL.
!
!    To compute inverse(A)*C, follow CSPCO by CSPSL.
!
!    To compute inverse(A), follow CSPCO by CSPDI.
!
!    To compute determinant(A), follow CSPCO by CSPDI.
!
!  Packed storage:
!
!    The following program segment will pack the upper
!    triangle of a symmetric matrix.
!
!      k = 0
!      do j = 1, n
!        do i = 1, j
!          k = k + 1
!          ap(k) = a(i,j)
!        end do
!      end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) AP(N*(N+1)/2); on input, the packed form of a
!    symmetric matrix.  The columns of the upper triangle are stored
!    sequentially in a one-dimensional array.  On output, a block diagonal
!    matrix and the multipliers which were used to obtain it, stored in packed
!    form.  The factorization can be written A = U*D*U' where U is a product
!    of permutation and unit upper triangular matrices, U' is the transpose
!    of U, and D is block diagonal with 1 by 1 and 2 by 2 blocks.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, real ( kind = 4 ) RCOND, an estimate of the reciprocal condition of
!    the matrix.  For the system A*X = B, relative perturbations in A and B
!    of size EPSILON may cause relative perturbations in X of size
!    (EPSILON/RCOND).  If RCOND is so small that the logical expression
!      1.0 + RCOND == 1.0
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate underflows.
!
!    Workspace, complex Z(N), a work vector whose contents are usually
!    unimportant.  If A is close to a singular matrix, then Z is an
!    approximate null vector in the sense that
!      norm(A*Z) = RCOND * norm(A) * norm(Z).
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) ak
  complex ( kind = 4 ) akm1
  real    ( kind = 4 ) anorm
  complex ( kind = 4 ) ap((n*(n+1))/2)
  complex ( kind = 4 ) bk
  complex ( kind = 4 ) bkm1
  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) cdotu
  complex ( kind = 4 ) csign1
  complex ( kind = 4 ) denom
  complex ( kind = 4 ) ek
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ikm1
  integer ( kind = 4 ) ikp1
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) km1k
  integer ( kind = 4 ) km1km1
  integer ( kind = 4 ) kp
  integer ( kind = 4 ) kps
  integer ( kind = 4 ) ks
  real    ( kind = 4 ) rcond
  real    ( kind = 4 ) s
  real    ( kind = 4 ) scasum
  complex ( kind = 4 ) t
  real    ( kind = 4 ) ynorm
  complex ( kind = 4 ) z(n)
!
!  Find norm of A using only upper half.
!
  j1 = 1

  do j = 1, n

    z(j) = cmplx ( scasum ( j, ap(j1:j1+j-1), 1 ), 0.0E+00 )
    ij = j1
    j1 = j1 + j

    do i = 1, j-1
      z(i) = cmplx ( real ( z(i) ) + cabs1 ( ap(ij) ), 0.0E+00 )
      ij = ij + 1
    end do

  end do

  anorm = 0.0E+00
  do j = 1, n
    anorm = max ( anorm, real ( z(j) ) )
  end do
!
!  Factor.
!
  call cspfa ( ap, n, ipvt, info )
!
!  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
!
!  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
!
!  The components of E are chosen to cause maximum local
!  growth in the elements of W where U*D*W = E.
!
!  The vectors are frequently rescaled to avoid overflow.
!
!  Solve U*D*W = E.
!
  ek = cmplx ( 1.0E+00, 0.0E+00 )
  z(1:n) = cmplx ( 0.0E+00, 0.0E+00 )

  k = n
  ik = ( n * ( n - 1 ) ) / 2

  do while ( 0 < k )

    kk = ik + k
    ikm1 = ik - ( k - 1 )

    if ( ipvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    kp = abs ( ipvt(k) )
    kps = k + 1 - ks

    if ( kp /= kps ) then
      t      = z(kps)
      z(kps) = z(kp)
      z(kp)  = t
    end if

    if ( cabs1 ( z(k) ) /= 0.0E+00 ) then
      ek = csign1 ( ek, z(k) )
    end if

    z(k) = z(k) + ek
    z(1:k-ks) = z(1:k-ks) + z(k) * ap(ik+1:ik+k-ks)

    if ( ks /= 1 ) then

      if ( cabs1 ( z(k-1) ) /= 0.0E+00 ) then
        ek = csign1 ( ek, z(k-1) )
      end if

      z(k-1) = z(k-1) + ek
      z(1:k-ks) = z(1:k-ks) + z(k-1) * ap(ikm1+1:ikm1+k-ks)

    end if

    if ( ks /= 2 ) then

      if ( cabs1 ( ap(kk) ) < cabs1 ( z(k) ) ) then
        s = cabs1 ( ap(kk) ) / cabs1 ( z(k) )
        z(1:n) = z(1:n) * s
        ek = cmplx ( s, 0.0E+00 ) * ek
      end if

      if ( cabs1 ( ap(kk) ) /= 0.0E+00 ) then
        z(k) = z(k) / ap(kk)
      else
        z(k) = cmplx ( 1.0E+00, 0.0E+00 )
      end if

    else

      km1k = ik + k - 1
      km1km1 = ikm1 + k - 1
      ak = ap(kk) / ap(km1k)
      akm1 = ap(km1km1) / ap(km1k)
      bk = z(k) / ap(km1k)
      bkm1 = z(k-1) / ap(km1k)
      denom = ak * akm1 - 1.0E+00
      z(k) = ( akm1 * bk - bkm1 ) / denom
      z(k-1) = ( ak * bkm1 - bk ) / denom

    end if

    k = k - ks
    ik = ik - k
    if ( ks == 2 ) then
      ik = ik - ( k + 1 )
    end if

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
!
!  Solve trans(U) * Y = W.
!
  k = 1
  ik = 0

  do while ( k <= n )

    if ( ipvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= 1 ) then

      z(k) = z(k) + cdotu ( k-1, ap(ik+1:ik+k-1), 1, z(1:k-1), 1 )
      ikp1 = ik + k

      if ( ks == 2 ) then
        z(k+1) = z(k+1) + cdotu ( k-1, ap(ikp1+1:ikp1+k-1), 1, z(1:k-1), 1 )
      end if

      kp = abs ( ipvt(k) )

      if ( kp /= k ) then
        t     = z(k)
        z(k)  = z(kp)
        z(kp) = t
      end if

    end if

    ik = ik + k
    if ( ks == 2 ) then
      ik = ik + ( k + 1 )
    end if

    k = k + ks

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = 1.0E+00
!
!  Solve U*D*V = Y.
!
  k = n
  ik = ( n * ( n - 1 ) ) / 2

  do while ( 0 < k )

    kk = ik + k
    ikm1 = ik - ( k - 1 )

    if ( ipvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= ks ) then

      kp = abs ( ipvt(k) )
      kps = k + 1 - ks

      if ( kp /= kps ) then
        t      = z(kps)
        z(kps) = z(kp)
        z(kp)  = t
      end if

      z(1:k-ks) = z(1:k-ks) + z(k) * ap(ik+1:ik+k-ks)

      if ( ks == 2 ) then
        z(1:k-ks) = z(1:k-ks) + z(k-1) * ap(ikm1+1:ikm1+k-ks)
      end if

    end if

    if ( ks /= 2 ) then

      if ( cabs1 ( ap(kk) ) < cabs1 ( z(k) ) ) then
        s = cabs1 ( ap(kk) ) / cabs1 ( z(k) )
        z(1:n) = z(1:n) * s
        ynorm = s * ynorm
      end if

      if ( cabs1 ( ap(kk) ) /= 0.0E+00 ) then
        z(k) = z(k) / ap(kk)
      else
        z(k) = cmplx ( 1.0E+00, 0.0E+00 )
      end if

    else

      km1k = ik + k - 1
      km1km1 = ikm1 + k - 1
      ak = ap(kk) / ap(km1k)
      akm1 = ap(km1km1) / ap(km1k)
      bk = z(k) / ap(km1k)
      bkm1 = z(k-1) / ap(km1k)
      denom = ak * akm1 - 1.0E+00
      z(k) = ( akm1 * bk - bkm1 ) / denom
      z(k-1) = ( ak * bkm1 - bk ) / denom

    end if

    k = k - ks
    ik = ik - k

    if ( ks == 2 ) then
      ik = ik - ( k + 1 )
    end if

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm
!
!  Solve U' * Z = V.
!
  k = 1
  ik = 0

  do while ( k <= n )

    if ( ipvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= 1 ) then

      z(k) = z(k) + cdotu ( k-1, ap(ik+1:ik+k-1), 1, z(1:k-1), 1 )
      ikp1 = ik + k

      if ( ks == 2 ) then
        z(k+1) = z(k+1) + cdotu ( k-1, ap(ikp1+1:ikp1+k-1), 1, z(1:k-1), 1 )
      end if

      kp = abs ( ipvt(k) )

      if ( kp /= k ) then
        t     = z(k)
        z(k)  = z(kp)
        z(kp) = t
      end if

    end if

    ik = ik + k

    if ( ks == 2 ) then
      ik = ik + ( k + 1 )
    end if

    k = k + ks

  end do
!
!  Make ZNORM = 1.
!
  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm

  if ( anorm /= 0.0E+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0E+00
  end if

  return
end
subroutine cspdi ( ap, n, ipvt, det, work, job )

!*****************************************************************************80
!
!! CSPDI sets the determinant and inverse of a complex symmetric packed matrix.
!
!  Discussion:
!
!    CSPDI uses the factors from CSPFA.
!
!    The matrix is stored in packed form.
!
!    A division by zero will occur if the inverse is requested and CSPCO has
!    set RCOND to 0.0 or CSPFA has set INFO nonzero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) AP(N*(N+1)/2); on input, the matrix factors
!    from CSPFA.  On output, if the inverse was requested, the upper
!    triangle of the inverse of the original matrix, stored in packed
!    form.  The columns of the upper triangle are stored sequentially
!    in a one-dimensional array.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from CSPFA.
!
!    Output, complex ( kind = 4 ) DET(2), the determinant of the original matrix.
!    Determinant = DET(1) * 10.0**DET(2) with 1.0 <= abs ( DET(1) ) < 10.0
!    or DET(1) = 0.0.  Also, DET(2) is strictly real.
!
!    Workspace, complex WORK(N).
!
!    Input, integer ( kind = 4 ) JOB, has the decimal expansion AB where
!    if B /= 0, the inverse is computed,
!    if A /= 0, the determinant is computed,
!    For example, JOB = 11 gives both.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) ak
  complex ( kind = 4 ) akkp1
  complex ( kind = 4 ) akp1
  complex ( kind = 4 ) ap((n*(n+1))/2)
  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) cdotu
  complex ( kind = 4 ) d
  complex ( kind = 4 ) det(2)
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ikp1
  integer ( kind = 4 ) iks
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jb
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) jkp1
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) kkp1
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) ks
  integer ( kind = 4 ) ksj
  integer ( kind = 4 ) kskp1
  integer ( kind = 4 ) kstep
  logical              nodet
  logical              noinv
  complex ( kind = 4 ) t
  complex ( kind = 4 ) work(n)

  noinv = mod ( job, 10 ) == 0
  nodet = mod ( job, 100 ) / 10 == 0

  if ( .not. nodet ) then

    det(1) = cmplx ( 1.0E+00, 0.0E+00 )
    det(2) = cmplx ( 0.0E+00, 0.0E+00 )
    t = cmplx ( 0.0E+00, 0.0E+00 )
    ik = 0

    do k = 1, n

      kk = ik + k
      d = ap(kk)
!
!  2 by 2 block
!  Use det (D  T)  =  ( D / T * C - T ) * T
!          (T  C)
!  to avoid underflow/overflow troubles.
!  Take two passes through scaling.  Use T for flag.
!
      if ( ipvt(k) <= 0 ) then

        if ( cabs1 ( t ) == 0.0E+00 ) then
          ikp1 = ik + k
          kkp1 = ikp1 + k
          t = ap(kkp1)
          d = ( d / t ) * ap(kkp1+1) - t
        else
          d = t
          t = cmplx ( 0.0E+00, 0.0E+00 )
        end if

      end if

      if ( .not. nodet ) then

        det(1) = det(1) * d

        if ( cabs1 ( det(1) ) /= 0.0E+00 ) then

          do while ( cabs1 ( det(1) ) < 1.0E+00 )
            det(1) = det(1) * cmplx ( 10.0E+00, 0.0E+00 )
            det(2) = det(2) - cmplx ( 1.0E+00, 0.0E+00 )
          end do

          do while ( 10.0E+00 <= cabs1 ( det(1) ) )
            det(1) = det(1) / cmplx ( 10.0E+00, 0.0E+00 )
            det(2) = det(2) + cmplx ( 1.0E+00, 0.0E+00 )
          end do

        end if

      end if

      ik = ik + k

    end do

  end if
!
!  Compute inverse ( A ).
!
  if ( .not. noinv ) then

    k = 1
    ik = 0

    do while ( k <= n )

      km1 = k - 1
      kk = ik + k
      ikp1 = ik + k

      if ( 0 <= ipvt(k) ) then
!
!  1 by 1
!
        ap(kk) = cmplx ( 1.0E+00, 0.0E+00 ) / ap(kk)

        if ( 1 <= km1 ) then

          work(1:km1) = ap(ik+1:ik+km1)
          ij = 0

          do j = 1, km1
            jk = ik + j
            ap(jk) = cdotu ( j, ap(ij+1:ij+j), 1, work, 1 )
            ap(ik+1:ik+j-1) = ap(ik+1:ik+j-1) + work(j) * ap(ij+1:ij+j-1)
            ij = ij + j
          end do

          ap(kk) = ap(kk) + cdotu ( km1, work, 1, ap(ik+1:ik+km1), 1 )

        end if

        kstep = 1
!
!  2 by 2
!
      else

        kkp1 = ikp1 + k
        t = ap(kkp1)
        ak = ap(kk) / t
        akp1 = ap(kkp1+1) / t
        akkp1 = ap(kkp1) / t
        d = t * ( ak * akp1 - cmplx ( 1.0E+00, 0.0E+00 ) )
        ap(kk) = akp1 / d
        ap(kkp1+1) = ak / d
        ap(kkp1) = -akkp1 / d

        if ( 1 <= km1 ) then

          work(1:km1) = ap(ikp1+1:ikp1+km1)
          ij = 0

          do j = 1, km1
            jkp1 = ikp1 + j
            ap(jkp1) = cdotu ( j, ap(ij+1:ij+j), 1, work, 1 )
            ap(ikp1+1:ikp1+j-1) = &
              ap(ikp1+1:ikp1+j-1) + work(j) * ap(ij+1:ij+j-1)
            ij = ij + j
          end do

          ap(kkp1+1) = ap(kkp1+1) &
            + cdotu ( km1, work, 1, ap(ikp1+1:ikp1+km1), 1 )

          ap(kkp1) = ap(kkp1) &
            + cdotu ( km1, ap(ik+1:ik+km1), 1, ap(ikp1+1:ikp1+km1), 1 )

          work(1:km1) = ap(ik+1:ik+km1)
          ij = 0

          do j = 1, km1
            jk = ik + j
            ap(jk) = cdotu ( j, ap(ij+1:ij+j), 1, work, 1 )
            ap(ik+1:ik+j-1) = ap(ik+1:ik+j-1) + work(j) * ap(ij+1:ij+j-1)
            ij = ij + j
          end do

          ap(kk) = ap(kk) + cdotu ( km1, work, 1, ap(ik+1:ik+km1), 1 )

        end if

        kstep = 2

      end if
!
!  Swap.
!
      ks = abs ( ipvt(k) )

      if ( ks /= k ) then

        iks = ( ks * ( ks - 1 ) ) / 2
        call cswap ( ks, ap(iks+1:iks+ks), 1, ap(ik+1:ik+ks), 1 )
        ksj = ik + ks

        do jb = ks, k
          j = k + ks - jb
          jk = ik + j
          t       = ap(jk)
          ap(jk)  = ap(ksj)
          ap(ksj) = t
          ksj = ksj - ( j - 1 )
        end do

        if ( kstep /= 1 ) then
          kskp1 = ikp1 + ks
          t         = ap(kskp1)
          ap(kskp1) = ap(kkp1)
          ap(kkp1)  = t
        end if

      end if

      ik = ik + k

      if ( kstep == 2 ) then
        ik = ik + k + 1
      end if

      k = k + kstep

    end do

  end if

  return
end
subroutine cspfa ( ap, n, ipvt, info )

!*****************************************************************************80
!
!! CSPFA factors a complex symmetric matrix stored in packed form.
!
!  Discussion:
!
!    The factorization is done by elimination with symmetric pivoting.
!
!    To solve A*X = B, follow CSPFA by CSPSL.
!
!    To compute inverse(A)*C, follow CSPFA by CSPSL.
!
!    To compute determinant(A), follow CSPFA by CSPDI.
!
!    To compute inverse(A), follow CSPFA by CSPDI.
!
!  Packed storage:
!
!    The following program segment will pack the upper
!    triangle of a symmetric matrix.
!
!      k = 0
!      do j = 1, n
!        do i = 1, j
!          k = k + 1
!          ap(k)  = a(i,j)
!        end do
!      end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) AP(N*(N+1)/2); On input, the packed form of a
!    symmetric matrix A.  The columns of the upper triangle are stored
!    sequentially in a one-dimensional array.  On output, a block diagonal
!    matrix and the multipliers which were used to obtain it stored in
!    packed form.  The factorization can be written A = U*D*U' where U
!    is a product of permutation and unit upper triangular matrices,
!    U' is the transpose of U, and D is block diagonal with 1 by 1 and
!    2 by 2 blocks.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO.
!    0, normal value.
!    K, if the K-th pivot block is singular.  This is not an error condition
!    for this subroutine, but it does indicate that CSPSL or CSPDI may
!    divide by zero if called.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 4 ) absakk
  complex ( kind = 4 ) ak
  complex ( kind = 4 ) akm1
  real    ( kind = 4 ) alpha
  complex ( kind = 4 ) ap((n*(n+1))/2)
  complex ( kind = 4 ) bk
  complex ( kind = 4 ) bkm1
  real    ( kind = 4 ) cabs1
  real    ( kind = 4 ) colmax
  complex ( kind = 4 ) denom
  integer ( kind = 4 ) icamax
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ijj
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ikm1
  integer ( kind = 4 ) im
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imim
  integer ( kind = 4 ) imj
  integer ( kind = 4 ) imk
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) jkm1
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) jmim
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) km1k
  integer ( kind = 4 ) km1km1
  integer ( kind = 4 ) km2
  integer ( kind = 4 ) kstep
  complex ( kind = 4 ) mulk
  complex ( kind = 4 ) mulkm1
  real    ( kind = 4 ) rowmax
  logical              swap
  complex ( kind = 4 ) t
!
!  Initialize.
!
!  ALPHA is used in choosing pivot block size.
!
  alpha = ( 1.0E+00 + sqrt ( 17.0E+00 ) ) / 8.0E+00

  info = 0
!
!  Main loop on K, which goes from N to 1.
!
  k = n
  ik = ( n * ( n - 1 ) ) / 2

  do
!
!  Leave the loop if K = 0 or K = 1.
!
    if ( k == 0 ) then
      exit
    end if

    if ( k == 1 ) then
      ipvt(1) = 1
      if ( cabs1 ( ap(1) ) == 0.0E+00 ) then
        info = 1
      end if
      exit
    end if
!
!  This section of code determines the kind of
!  elimination to be performed.  When it is completed,
!  KSTEP will be set to the size of the pivot block, and
!  SWAP will be set to .true. if an interchange is
!  required.
!
    km1 = k - 1
    kk = ik + k
    absakk = cabs1 ( ap(kk) )
!
!  Determine the largest off-diagonal element in column K.
!
    imax = icamax ( k-1, ap(ik+1:ik+k-1), 1 )
    imk = ik + imax
    colmax = cabs1 ( ap(imk) )

    if ( alpha * colmax <= absakk ) then

      kstep = 1
      swap = .false.
!
!  Determine the largest off-diagonal element in row IMAX.
!
    else

      rowmax = 0.0E+00
      im = ( imax * ( imax - 1 ) ) / 2
      imj = im + 2 * imax

      do j = imax + 1, k
        rowmax = max ( rowmax, cabs1 ( ap(imj) ) )
        imj = imj + j
      end do

      if ( imax /= 1 ) then
        jmax = icamax ( imax-1, ap(im+1:im+imax-1), 1 )
        jmim = jmax + im
        rowmax = max ( rowmax, cabs1 ( ap(jmim) ) )
      end if

      imim = imax + im

      if ( alpha * rowmax <= cabs1 ( ap(imim) ) ) then
        kstep = 1
        swap = .true.
      else if ( alpha * colmax * ( colmax / rowmax ) <= absakk ) then
        kstep = 1
        swap = .false.
      else
        kstep = 2
        swap = imax /= km1
      end if

    end if
!
!  Column K is zero.  Set INFO and iterate the loop.
!
    if ( max ( absakk, colmax ) == 0.0E+00 ) then
      ipvt(k) = k
      info = k
      ik = ik - ( k - 1 )
      if ( kstep == 2 ) then
        ik = ik - ( k - 2 )
      end if
      k = k - kstep
      cycle
    end if

    if ( kstep /= 2 ) then
!
!  1 x 1 pivot block.
!
      if ( swap ) then

        call cswap ( imax, ap(im+1:im+imax), 1, ap(ik+1:ik+imax), 1 )
        imj = ik + imax

        do jj = imax, k
          j = k + imax - jj
          jk = ik + j
          t       = ap(jk)
          ap(jk)  = ap(imj)
          ap(imj) = t
          imj = imj - ( j - 1 )
        end do

      end if
!
!  Perform the elimination.
!
      ij = ik - ( k - 1 )

      do jj = 1, km1
        j = k - jj
        jk = ik + j
        mulk = -ap(jk) / ap(kk)
        t = mulk
        ap(ij+1:ij+j) = ap(ij+1:ij+j) + t * ap(ik+1:ik+j)
        ijj = ij + j
        ap(jk) = mulk
        ij = ij - ( j - 1 )
      end do
!
!  Set the pivot array.
!
      if ( swap ) then
        ipvt(k) = imax
      else
        ipvt(k) = k
      end if
!
!  2 x 2 pivot block.
!
    else

      km1k = ik + k - 1
      ikm1 = ik - ( k - 1 )

      if ( swap ) then

        call cswap ( imax, ap(im+1:im+imax), 1, ap(ikm1+1:ikm1+imax), 1 )
        imj = ikm1 + imax

        do jj = imax, km1
          j = km1 + imax - jj
          jkm1 = ikm1 + j
          t        = ap(jkm1)
          ap(jkm1) = ap(imj)
          ap(imj)  = t
          imj = imj - ( j - 1 )
        end do

        t        = ap(km1k)
        ap(km1k) = ap(imk)
        ap(imk)  = t

      end if
!
!  Perform the elimination.
!
      km2 = k - 2

      if ( km2 /= 0 ) then

        ak = ap(kk) / ap(km1k)
        km1km1 = ikm1 + k - 1
        akm1 = ap(km1km1) / ap(km1k)
        denom = 1.0E+00 - ak * akm1
        ij = ik - ( k - 1 ) - ( k - 2 )

        do jj = 1, km2
          j = km1 - jj
          jk = ik + j
          bk = ap(jk) / ap(km1k)
          jkm1 = ikm1 + j
          bkm1 = ap(jkm1) / ap(km1k)
          mulk = ( akm1 * bk - bkm1 ) / denom
          mulkm1 = ( ak * bkm1 - bk ) / denom
          t = mulk
          ap(ij+1:ij+j) = ap(ij+1:ij+j) + t * ap(ik+1:ik+j)
          t = mulkm1
          ap(ij+1:ij+j) = ap(ij+1:ij+j) + t * ap(ikm1+1:ikm1+j)
          ap(jk) = mulk
          ap(jkm1) = mulkm1
          ijj = ij + j
          ij = ij - ( j - 1 )
        end do

      end if
!
!  Set the pivot array.
!
      if ( swap ) then
        ipvt(k) = -imax
      else
        ipvt(k) = 1 - k
      end if

      ipvt(k-1) = ipvt(k)

    end if

    ik = ik - ( k - 1 )

    if ( kstep == 2 ) then
      ik = ik - ( k - 2 )
    end if

    k = k - kstep

  end do

  return
end
subroutine cspsl ( ap, n, ipvt, b )

!*****************************************************************************80
!
!! CSPSL solves a complex symmetric system factored by CSPFA.
!
!  Discussion:
!
!    A division by zero may occur if CSPCO has set RCOND == 0.0
!    or CSPFA has set INFO /= 0.
!
!    To compute inverse(A) * C where C is a matrix with P columns:
!
!      call cspfa(ap,n,ipvt,info)
!
!      if (info /= 0) then
!        error
!      end if
!
!      do j = 1, p
!        call cspsl(ap,n,ipvt,c(1,j))
!      end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input, complex ( kind = 4 ) AP(N*(N+1)/2), the output from CSPFA.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from CSPFA.
!
!    Input/output, complex ( kind = 4 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) ak
  complex ( kind = 4 ) akm1
  complex ( kind = 4 ) ap((n*(n+1))/2)
  complex ( kind = 4 ) b(n)
  complex ( kind = 4 ) bk
  complex ( kind = 4 ) bkm1
  complex ( kind = 4 ) cdotu
  complex ( kind = 4 ) denom
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ikm1
  integer ( kind = 4 ) ikp1
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) km1k
  integer ( kind = 4 ) km1km1
  integer ( kind = 4 ) kp
  complex ( kind = 4 ) t
!
!  Loop backward applying the transformations and d inverse to b.
!
  k = n
  ik = ( n * ( n - 1 ) ) / 2

  do while ( 0 < k )

    kk = ik + k
    if ( 0 <= ipvt(k) ) then
!
!  1 x 1 pivot block.
!
      if ( k /= 1 ) then

        kp = ipvt(k)
        if ( kp /= k ) then
          t     = b(k)
          b(k)  = b(kp)
          b(kp) = t
        end if

        b(1:k-1) = b(1:k-1) + b(k) * ap(ik+1:ik+k-1)

      end if
!
!  Apply D inverse.
!
      b(k) = b(k) / ap(kk)
      k = k - 1
      ik = ik - k
!
!  2 x 2 pivot block.
!
    else

      ikm1 = ik - ( k - 1 )

      if ( k /= 2 ) then

        kp = abs ( ipvt(k) )

        if ( kp /= k - 1 ) then
          t      = b(k-1)
          b(k-1) = b(kp)
          b(kp)  = t
        end if

        b(1:k-2) = b(1:k-2) + b(k)   * ap(ik+1:ik+k-2)
        b(1:k-2) = b(1:k-2) + b(k-1) * ap(ikm1+1:ikm1+k-2)

      end if
!
!  Apply D inverse.
!
      km1k = ik + k - 1
      kk = ik + k
      ak = ap(kk) / ap(km1k)
      km1km1 = ikm1 + k - 1
      akm1 = ap(km1km1) / ap(km1k)
      bk = b(k) / ap(km1k)
      bkm1 = b(k-1) / ap(km1k)
      denom = ak * akm1 - 1.0E+00
      b(k) = ( akm1 * bk - bkm1 ) / denom
      b(k-1) = ( ak * bkm1 - bk ) / denom
      k = k - 2
      ik = ik - ( k + 1 ) - k

    end if

  end do
!
!  Loop forward applying the transformations.
!
  k = 1
  ik = 0

  do while ( k <= n )
!
!  1 x 1 pivot block.
!
    if ( 0 <= ipvt(k) ) then

      if ( k /= 1 ) then
        b(k) = b(k) + cdotu ( k-1, ap(ik+1:ik+k-1), 1, b(1:k-1), 1 )
        kp = ipvt(k)
        if ( kp /= k ) then
          t     = b(k)
          b(k)  = b(kp)
          b(kp) = t
        end if
      end if

      ik = ik + k
      k = k + 1
!
!  2 x 2 pivot block.
!
    else

      if ( k /= 1 ) then

        b(k) = b(k) + cdotu ( k-1, ap(ik+1:ik+k-1), 1, b(1:k-1), 1 )
        ikp1 = ik + k
        b(k+1) = b(k+1) + cdotu ( k-1, ap(ikp1+1:ikp1+k-1), 1, b(1:k-1), 1 )
        kp = abs ( ipvt(k) )

        if ( kp /= k ) then
          t     = b(k)
          b(k)  = b(kp)
          b(kp) = t
        end if

      end if

      ik = ik + k + k + 1
      k = k + 2

    end if

  end do

  return
end
subroutine csvdc ( x, ldx, n, p, s, e, u, ldu, v, ldv, work, job, info )

!*****************************************************************************80
!
!! CSVDC applies the singular value decompostion to an N by P matrix.
!
!  Discussion:
!
!    The routine reduces a complex N by P matrix X, by unitary transformations
!    U and V, to diagonal form.
!
!    The diagonal elements, S(I), are the singular values of Z.  The
!    columns of U are the corresponding left singular vectors,
!    and the columns of V the right singular vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) X(LDX,P); on input, the matrix whose
!    singular value decomposition is to be computed.  X is destroyed on output.
!
!    Input, integer ( kind = 4 ) LDX, the leading dimension of X.  N <= LDX.
!
!    Input, integer ( kind = 4 ) N, the number of rows of the matrix.
!
!    Input, integer ( kind = 4 ) P, the number of columns of the matrix X.
!
!    Output, complex ( kind = 4 ) S(MM), where MM = max ( N + 1, P ), the
!    first min ( N, P ) entries of S contain the singular values of X arranged
!    in descending order of magnitude.
!
!    Output, complex ( kind = 4 ) E(MM), where MM = max ( N + 1, P ),
!    ordinarily contains zeros on output.  However, see the discussion of INFO
!    for exceptions.
!
!    Output, complex ( kind = 4 ) U(LDU,K).  If JOBA == 1 then K == N; if
!    JOBA >= 2, then K == min ( N, P ).  U contains the matrix of left singular
!    vectors.  U is not referenced if JOBA == 0.  If N <= P or if JOBA > 2,
!    then U may be identified with X in the subroutine call.
!
!    Input, integer ( kind = 4 ) LDU, the leading dimension of U.  N <= LDU.
!
!    Output, complex ( kind = 4 ) V(LDV,P), if requested, the matrix of right
!    singular vectors.  If P <= N, V may be identified with X in the call.
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of V.  P <= LDV.
!
!    Workspace, complex WORK(N).
!
!    Input, integer ( kind = 4 ) JOB, controls the computation of the singular
!    vectors.  It has the decimal expansion AB meaning:
!    A =  0, do not compute the left singular vectors.
!    A =  1, return the N left singular vectors in U.
!    A >= 2, returns the first min ( N, P ) left singular vectors in U.
!    B =  0, do not compute the right singular vectors.
!    B =  1, return the right singular vectors in V.
!
!    Output, integer ( kind = 4 ) INFO.  The singular values and their
!    corresponding singular vectors, S(INFO+1), S(INFO+2),..., S(M) are correct.
!    Here M = min ( N, P ).  Thus if INFO == 0, all the singular values and
!    their vectors are correct.  In any event, the matrix
!      B = hermitian(U)*X*V
!    is the bidiagonal matrix with the elements of S on its diagonal
!    and the elements of E on its super-diagonal.  Hermitian(U)
!    is the conjugate-transpose of U.  Thus the singular values of X
!    and B are the same.
!
  implicit none

  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) ldx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p

  real    ( kind = 4 ) b
  real    ( kind = 4 ) c
  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) cdotc
  real    ( kind = 4 ) cs
  complex ( kind = 4 ) csign2
  complex ( kind = 4 ) e(*)
  real    ( kind = 4 ) el
  real    ( kind = 4 ) emm1
  real    ( kind = 4 ) f
  real    ( kind = 4 ) g
  integer ( kind = 4 ) i
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
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: maxit = 30
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) mm1
  integer ( kind = 4 ) mp1
  integer ( kind = 4 ) nct
  integer ( kind = 4 ) nctp1
  integer ( kind = 4 ) ncu
  integer ( kind = 4 ) nrt
  integer ( kind = 4 ) nrtp1
  complex ( kind = 4 ) r
  complex ( kind = 4 ) s(*)
  real    ( kind = 4 ) scnrm2
  real    ( kind = 4 ) scale
  real    ( kind = 4 ) shift
  real    ( kind = 4 ) sl
  real    ( kind = 4 ) sm
  real    ( kind = 4 ) smm1
  real    ( kind = 4 ) sn
  complex ( kind = 4 ) t
  real    ( kind = 4 ) t1
  real    ( kind = 4 ) test
  complex ( kind = 4 ) u(ldu,*)
  complex ( kind = 4 ) v(ldv,p)
  logical              wantu
  logical              wantv
  complex ( kind = 4 ) work(n)
  complex ( kind = 4 ) x(ldx,p)
  real    ( kind = 4 ) ztest
!
!  Determine what is to be computed.
!
  wantu = .false.
  wantv = .false.
  jobu = mod ( job, 100 ) / 10

  if ( 1 < jobu ) then
    ncu = min ( n, p )
  else
    ncu = n
  end if

  if ( jobu /= 0 ) then
    wantu = .true.
  end if

  if ( mod ( job, 10 ) /= 0 ) then
    wantv = .true.
  end if
!
!  Reduce X to bidiagonal form, storing the diagonal elements
!  in S and the super-diagonal elements in E.
!
  info = 0
  nct = min ( n - 1, p )
  nrt = max ( 0, min ( p - 2, n ) )
  lu = max ( nct, nrt )

  do l = 1, lu
!
!  Compute the transformation for the L-th column and
!  place the L-th diagonal in S(L).
!
    if ( l <= nct ) then

      s(l) = cmplx ( scnrm2 ( n-l+1, x(l:n,l), 1 ), 0.0E+00 )

      if ( cabs1 ( s(l) ) /= 0.0E+00 ) then

        if ( cabs1 ( x(l,l) ) /= 0.0E+00 ) then
          s(l) = csign2 ( s(l), x(l,l) )
        end if

        x(l:n,l) = x(l:n,l) / s(l)
        x(l,l) = cmplx ( 1.0E+00, 0.0E+00 ) + x(l,l)

      end if

      s(l) = -s(l)

    end if

    do j = l + 1, p

      if ( l <= nct ) then
        if ( cabs1 ( s(l) ) /= 0.0E+00 ) then
          t = -cdotc ( n-l+1, x(l:n,l), 1, x(l:n,j), 1 ) / x(l,l)
          x(l:n,j) = x(l:n,j) + t * x(l:n,l)
        end if
      end if
!
!  Place the L-th row of X into E for the
!  subsequent calculation of the row transformation.
!
      e(j) = conjg ( x(l,j) )

    end do
!
!  Place the transformation in U for subsequent back multiplication.
!
    if ( wantu .and. l <= nct ) then
      u(l:n,l) = x(l:n,l)
    end if

    if ( l <= nrt ) then
!
!  Compute the L-th row transformation and place the
!  L-th super-diagonal in E(L).
!
      e(l) = cmplx ( scnrm2 ( p-l, e(l+1:p), 1 ), 0.0E+00 )

      if ( cabs1 ( e(l) ) /= 0.0E+00 ) then

        if ( cabs1 ( e(l+1) ) /= 0.0E+00 ) then
          e(l) = csign2 ( e(l), e(l+1) )
        end if

        e(l+1:p) = e(l+1:p) / e(l)
        e(l+1) = cmplx ( 1.0E+00, 0.0E+00 ) + e(l+1)

      end if

      e(l) = -conjg ( e(l) )
!
!  Apply the transformation.
!
      if ( l+1 <= n .and. cabs1 ( e(l) ) /= 0.0E+00 ) then

        work(l+1:n) = cmplx ( 0.0E+00, 0.0E+00 )

        do j = l+1, p
          work(l+1:n) = work(l+1:n) + e(j) * x(l+1:n,j)
        end do

        do j = l+1, p
          x(l+1:n,j) = x(l+1:n,j) + conjg ( -e(j) / e(l+1) ) * work(l+1:n)
        end do

      end if
!
!  Place the transformation in V for subsequent back multiplication.
!
      if ( wantv ) then
        v(l+1:p,l) = e(l+1:p)
      end if

    end if

  end do
!
!  Set up the final bidiagonal matrix of order M.
!
  m = min ( p, n + 1 )
  nctp1 = nct + 1
  nrtp1 = nrt + 1

  if ( nct < p ) then
    s(nctp1) = x(nctp1,nctp1)
  end if

  if ( n < m ) then
    s(m) = cmplx ( 0.0E+00, 0.0E+00 )
  end if

  if ( nrtp1 < m ) then
    e(nrtp1) = x(nrtp1,m)
  end if

  e(m) = cmplx ( 0.0E+00, 0.0E+00 )
!
!  If required, generate U.
!
  if ( wantu ) then

    do j = nctp1, ncu
      u(1:n,j) = cmplx ( 0.0E+00, 0.0E+00 )
      u(j,j) = cmplx ( 1.0E+00, 0.0E+00 )
    end do

    do ll = 1, nct

      l = nct - ll + 1

      if ( cabs1 ( s(l) ) /= 0.0E+00 ) then

        do j = l+1, ncu
          t = -cdotc ( n-l+1, u(l:n,l), 1, u(l:n,j), 1 ) / u(l,l)
          u(l:n,j) = u(l:n,j) + t * u(l:n,l)
        end do

        u(l:n,l) = - u(l:n,l)
        u(l,l) = cmplx ( 1.0E+00, 0.0E+00 ) + u(l,l)
        u(1:l-1,l) = cmplx ( 0.0E+00, 0.0E+00 )

      else

        u(1:n,l) = cmplx ( 0.0E+00, 0.0E+00 )
        u(l,l) = cmplx ( 1.0E+00, 0.0E+00 )

      end if

    end do

  end if
!
!  If it is required, generate V.
!
  if ( wantv ) then

    do ll = 1, p

      l = p - ll + 1

      if ( l <= nrt ) then

        if ( cabs1 ( e(l) ) /= 0.0E+00 ) then
          do j = l+1, p
            t = -cdotc ( p-l, v(l+1:p,l), 1, v(l+1:p,j), 1 ) / v(l+1,l)
            v(l+1:p,j) = v(l+1:p,j) + t * v(l+1:p,l)
          end do
        end if

      end if

      v(1:p,l) = cmplx ( 0.0E+00, 0.0E+00 )
      v(l,l) = cmplx ( 1.0E+00, 0.0E+00 )

    end do

  end if
!
!  Transform S and E so that they are real.
!
  do i = 1, m

    if ( cabs1 ( s(i) ) /= 0.0E+00 ) then

      t = cmplx ( abs ( s(i) ), 0.0E+00 )
      r = s(i) / t
      s(i) = t

      if ( i < m ) then
        e(i) = e(i) / r
      end if

      if ( wantu ) then
        u(1:n,i) = u(1:n,i) * r
      end if

    end if

    if ( i == m ) then
      exit
    end if

    if ( cabs1 ( e(i) ) /= 0.0E+00 ) then

      t = cmplx ( abs ( e(i) ), 0.0E+00 )
      r = t / e(i)
      e(i) = t
      s(i+1) = s(i+1) * r

      if ( wantv ) then
        v(1:p,i+1) = v(1:p,i+1) * r
      end if

    end if

  end do
!
!  Main iteration loop for the singular values.
!
  mm = m
  iter = 0

  do
!
!  Quit if all the singular values have been found.
!
    if ( m == 0 ) then
      exit
    end if
!
!  If too many iterations have been performed, set flag and return.
!
    if ( maxit <= iter ) then
      info = m
      exit
    end if
!
!  This section of the program inspects for negligible elements in S and E.
!
!  On completion, the variables KASE and L are set as follows.
!
!  KASE = 1     if S(M) and E(L-1) are negligible and L < M
!  KASE = 2     if S(L) is negligible and L < M
!  KASE = 3     if E(L-1) is negligible, L < M, and
!               S(L), ..., S(M) are not negligible (QR step).
!  KASE = 4     if E(M-1) is negligible (convergence).
!
    do ll = 1, m

      l = m - ll

      if ( l == 0 ) then
        exit
      end if

      test = abs ( s(l) ) + abs ( s(l+1) )
      ztest = test + abs ( e(l) )

      if ( ztest == test ) then
        e(l) = cmplx ( 0.0E+00, 0.0E+00 )
        exit
      end if

    end do

    if ( l == m - 1 ) then

      kase = 4

    else

      mp1 = m + 1

      do lls = l+1, mp1

        ls = m - lls + l + 1

        if ( ls == l ) then
          exit
        end if

        test = 0.0E+00

        if ( ls /= m ) then
          test = test + abs ( e(ls) )
        end if

        if ( ls /= l + 1 ) then
          test = test + abs ( e(ls-1) )
        end if

        ztest = test + abs ( s(ls) )

        if ( ztest == test ) then
          s(ls) = cmplx ( 0.0E+00, 0.0E+00 )
          exit
        end if

      end do

      if ( ls == l ) then
        kase = 3
      else if ( ls == m ) then
        kase = 1
      else
        kase = 2
        l = ls
      end if

    end if

    l = l + 1
!
!  Deflate negligible S(M).
!
    if ( kase == 1 ) then

      mm1 = m - 1
      f = real ( e(m-1) )
      e(m-1) = cmplx ( 0.0E+00, 0.0E+00 )

      do kk = l, mm1

        k = mm1 - kk + l
        t1 = real ( s(k) )
        call srotg ( t1, f, cs, sn )
        s(k) = cmplx ( t1, 0.0E+00 )

        if ( k /= l ) then
          f = -sn * real ( e(k-1) )
          e(k-1) = cs * e(k-1)
        end if

        if ( wantv ) then
          call csrot ( p, v(1:p,k), 1, v(1:p,m), 1, cs, sn )
        end if

      end do
!
!  Split at negligible S(L).
!
    else if ( kase == 2 ) then

      f = real ( e(l-1) )
      e(l-1) = cmplx ( 0.0E+00, 0.0E+00 )

      do k = l, m

        t1 = real ( s(k) )
        call srotg ( t1, f, cs, sn )
        s(k) = cmplx ( t1, 0.0E+00 )
        f = -sn * real ( e(k) )
        e(k) = cs * e(k)

        if ( wantu ) then
          call csrot ( n, u(1:n,k), 1, u(1:n,l-1), 1, cs, sn )
        end if

      end do
!
!  Perform one QR step.
!
    else if ( kase == 3 ) then
!
!  Calculate the shift.
!
      scale = max ( abs ( s(m) ), abs ( s(m-1) ), abs ( e(m-1) ), &
        abs ( s(l) ), abs ( e(l) ) )

      sm = real ( s(m) ) / scale
      smm1 = real ( s(m-1) ) / scale
      emm1 = real ( e(m-1) ) / scale
      sl = real ( s(l) ) / scale
      el = real ( e(l) ) / scale
      b = ( ( smm1 + sm ) * ( smm1 - sm ) + emm1**2 ) / 2.0E+00
      c = ( sm * emm1 )**2
      shift = 0.0E+00

      if ( b /= 0.0E+00 .or. c /= 0.0E+00 ) then
        shift = sqrt ( b**2 + c )
        if ( b < 0.0E+00 ) then
          shift = -shift
        end if
        shift = c / ( b + shift )
      end if

      f = ( sl + sm ) * ( sl - sm ) + shift
      g = sl * el
!
!  Chase zeros.
!
      mm1 = m - 1

      do k = l, mm1

        call srotg ( f, g, cs, sn )

        if ( k /= l ) then
          e(k-1) = cmplx ( f, 0.0E+00 )
        end if

        f = cs * real ( s(k) ) + sn * real ( e(k) )
        e(k) = cs * e(k) - sn * s(k)
        g = sn * real ( s(k+1) )
        s(k+1) = cs * s(k+1)

        if ( wantv ) then
          call csrot ( p, v(1:p,k), 1, v(1:p,k+1), 1, cs, sn )
        end if

        call srotg ( f, g, cs, sn )
        s(k) = cmplx ( f, 0.0E+00 )
        f = cs * real ( e(k) ) + sn * real ( s(k+1) )
        s(k+1) = -sn * e(k) + cs * s(k+1)
        g = sn * real ( e(k+1) )
        e(k+1) = cs * e(k+1)

        if ( wantu .and. k < n ) then
          call csrot ( n, u(1:n,k), 1, u(1:n,k+1), 1, cs, sn )
        end if

      end do

      e(m-1) = cmplx ( f, 0.0E+00 )
      iter = iter + 1
!
!  Convergence.
!
    else if ( kase == 4 ) then
!
!  Make the singular value positive.
!
      if ( real ( s(l) ) < 0.0E+00 ) then
        s(l) = -s(l)
        if ( wantv ) then
          v(1:p,l) = - v(1:p,l)
        end if
      end if
!
!  Order the singular values.
!
      do while ( l /= mm )

        if ( real ( s(l+1) ) <= real ( s(l) ) ) then
          exit
        end if

        t      = s(l)
        s(l)   = s(l+1)
        s(l+1) = t

        if ( wantv .and. l < p ) then
          call cswap ( p, v(1:p,l), 1, v(1:p,l+1), 1 )
        end if

        if ( wantu .and. l < n ) then
          call cswap ( n, u(1:n,l), 1, u(1:n,l+1), 1 )
        end if

        l = l + 1

      end do

      iter = 0
      m = m - 1

    end if

  end do

  return
end
subroutine ctrco ( t, ldt, n, rcond, z, job )

!*****************************************************************************80
!
!! CTRCO estimates the condition of a complex triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2007
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input, complex ( kind = 4 ) T(LDT,N), the triangular matrix.  The zero
!    elements of the matrix are not referenced, and the corresponding
!    elements of the array can be used to store other information.
!
!    Input, integer ( kind = 4 ) LDT, the leading dimension of T.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 4 ) RCOND, an estimate of the reciprocal condition
!    of T.  For the system T*X = B, relative perturbations in T and B of size
!    EPSILON may cause relative perturbations in X of size (EPSILON/RCOND).
!    If RCOND is so small that the logical expression
!      1.0 + RCOND == 1.0
!    is true, then T may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate
!    underflows.
!
!    Workspace, complex Z(N), a work vector whose contents are usually
!    unimportant.  If T is close to a singular matrix, then Z is
!    an approximate null vector in the sense that
!      norm(A*Z) = RCOND * norm(A) * norm(Z).
!
!    Input, integer ( kind = 4 ) JOB, indicates if matrix is upper or lower
!    triangular.
!    0, lower triangular.
!    nonzero, upper triangular.
!
  implicit none

  integer ( kind = 4 ) ldt
  integer ( kind = 4 ) n

  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) csign1
  complex ( kind = 4 ) ek
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  logical              lower
  real    ( kind = 4 ) rcond
  real    ( kind = 4 ) s
  real    ( kind = 4 ) scasum
  real    ( kind = 4 ) sm
  complex ( kind = 4 ) t(ldt,n)
  real    ( kind = 4 ) tnorm
  complex ( kind = 4 ) w
  complex ( kind = 4 ) wk
  complex ( kind = 4 ) wkm
  real    ( kind = 4 ) ynorm
  complex ( kind = 4 ) z(n)

  lower = ( job == 0 )
!
!  Compute 1-norm of T
!
  tnorm = 0.0E+00

  do j = 1, n

    if ( lower ) then
      l = n + 1 - j
      i1 = j
    else
      l = j
      i1 = 1
    end if

    tnorm = max ( tnorm, scasum ( l, t(i1:i1+l-1,j), 1 ) )

  end do
!
!  RCOND = 1/(norm(T)*(estimate of norm(inverse(T)))).
!
!  Estimate = norm(Z)/norm(Y) where T*Z = Y and hermitian(T)*Y = E.
!
!  Hermitian(T) is the conjugate transpose of T.
!
!  The components of E are chosen to cause maximum local
!  growth in the elements of Y.
!
!  The vectors are frequently rescaled to avoid overflow.
!
!  Solve hermitian(T)*Y = E.
!
  ek = cmplx ( 1.0E+00, 0.0E+00 )
  z(1:n) = cmplx ( 0.0E+00, 0.0E+00 )

  do kk = 1, n

    if ( lower ) then
      k = n + 1 - kk
    else
      k = kk
    end if

    if ( cabs1 ( z(k) ) /= 0.0E+00 ) then
      ek = csign1 ( ek, -z(k) )
    end if

    if ( cabs1 ( t(k,k) ) < cabs1 ( ek - z(k) ) ) then
      s = cabs1 ( t(k,k) ) / cabs1 ( ek - z(k) )
      z(1:n) = z(1:n) * s
      ek = cmplx ( s, 0.0E+00 ) * ek
    end if

    wk = ek - z(k)
    wkm = - ek - z(k)
    s = cabs1 ( wk )
    sm = cabs1 ( wkm )

    if ( cabs1 ( t(k,k) ) /= 0.0E+00 ) then
      wk = wk / conjg ( t(k,k) )
      wkm = wkm / conjg ( t(k,k) )
    else
      wk = cmplx ( 1.0E+00, 0.0E+00 )
      wkm = cmplx ( 1.0E+00, 0.0E+00 )
    end if

    if ( kk /= n ) then

      if ( lower ) then
        j1 = 1
        j2 = k - 1
      else
        j1 = k + 1
        j2 = n
      end if

      do j = j1, j2
        sm = sm + cabs1 ( z(j) + wkm * conjg ( t(k,j) ) )
        z(j) = z(j) + wk * conjg ( t(k,j) )
        s = s + cabs1 ( z(j) )
      end do

      if ( s < sm ) then
        w = wkm - wk
        wk = wkm
        z(j1:j2) = z(j1:j2) + w * conjg ( t(k,j1:j2) )
      end if

    end if

    z(k) = wk

  end do

  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = 1.0E+00
!
!  Solve T*Z = Y.
!
  do kk = 1, n

    if ( lower ) then
      k = kk
    else
      k = n + 1 - kk
    end if

    if ( cabs1 ( t(k,k) ) < cabs1 ( z(k) ) ) then
      s = cabs1 ( t(k,k) ) / cabs1 ( z(k) )
      z(1:n) = z(1:n) * s
      ynorm = s * ynorm
    end if

    if ( cabs1 ( t(k,k) ) /= 0.0E+00 ) then
      z(k) = z(k) / t(k,k)
    else
      z(k) = cmplx ( 1.0E+00, 0.0E+00 )
    end if

    if ( lower ) then
      i1 = k + 1
    else
      i1 = 1
    end if

    if ( kk < n ) then
      w = -z(k)
      z(i1:i1+n-kk-1) = z(i1:i1+n-kk-1) + w * t(i1:i1+n-kk-1,k)
    end if

  end do
!
!  Make ZNORM = 1.
!
  s = 1.0E+00 / scasum ( n, z, 1 )
  z(1:n) = z(1:n) * s
  ynorm = s * ynorm

  if ( tnorm /= 0.0E+00 ) then
    rcond = ynorm / tnorm
  else
    rcond = 0.0E+00
  end if

  return
end
subroutine ctrdi ( t, ldt, n, det, job, info )

!*****************************************************************************80
!
!! CTRDI computes the determinant and inverse of a complex triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input/output, complex ( kind = 4 ) T(LDT,N), the triangular matrix.  The
!    zero elements of the matrix are not referenced, and the corresponding
!    elements of the array can be used to store other information.
!    On output, if an inverse was requested, then T has been overwritten
!    by its inverse.
!
!    Input, integer ( kind = 4 ) LDT, the leading dimension of T.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, complex ( kind = 4 ) DET(2), the determinant of the original matrix,
!    if requested.  Otherwise not referenced.
!    Determinant = DET(1) * 10.0**DET(2) with 1.0 <= cabs1 ( DET(1) ) < 10.0
!    or DET(1) == 0.0.  Also, DET(2) is strictly real.
!
!    Input, integer ( kind = 4 ) JOB.
!    010, no determinant,    inverse, matrix is lower triangular.
!    011, no determinant,    inverse, matrix is upper triangular.
!    100,    determinant, no inverse.
!    110,    determinant,    inverse, matrix is lower triangular.
!    111,    determinant,    inverse, matrix is upper triangular.
!
!    Output, integer ( kind = 4 ) INFO.
!    0, an inverse was requested and the matrix is nonsingular.
!    K, an inverse was requested, but the K-th diagonal element
!    of T is zero.
!
  implicit none

  integer ( kind = 4 ) ldt
  integer ( kind = 4 ) n

  real    ( kind = 4 ) cabs1
  complex ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  complex ( kind = 4 ) t(ldt,n)
  complex ( kind = 4 ) temp

  if ( ( job / 100 ) /= 0 ) then

    det(1) = cmplx ( 1.0E+00, 0.0E+00 )
    det(2) = cmplx ( 0.0E+00, 0.0E+00 )

    do i = 1, n

      det(1) = det(1) * t(i,i)

      if ( cabs1 ( det(1) ) == 0.0E+00 ) then
        exit
      end if

      do while ( cabs1 ( det(1) ) < 1.0E+00 )
        det(1) = det(1) * cmplx ( 10.0E+00, 0.0E+00 )
        det(2) = det(2) - cmplx ( 1.0E+00, 0.0E+00 )
      end do

      do while ( 10.0E+00 <= cabs1 ( det(1) ) )
        det(1) = det(1) / cmplx ( 10.0E+00, 0.0E+00 )
        det(2) = det(2) + cmplx ( 1.0E+00, 0.0E+00 )
      end do

    end do

  end if
!
!  Compute inverse of upper triangular matrix.
!
  if ( mod ( job / 10, 10 ) /= 0 ) then

    if ( mod ( job, 10 ) /= 0 ) then

      info = 0

      do k = 1, n

        if ( cabs1 ( t(k,k) ) == 0.0E+00 ) then
          info = k
          exit
        end if

        t(k,k) = cmplx ( 1.0E+00, 0.0E+00 ) / t(k,k)
        temp = -t(k,k)
        t(1:k-1,k) = t(1:k-1,k) * temp

        do j = k+1, n
          temp = t(k,j)
          t(k,j) = cmplx ( 0.0E+00, 0.0E+00 )
          t(1:k,j) = t(1:k,j) + temp * t(1:k,k)
        end do

      end do
!
!  Compute inverse of lower triangular matrix.
!
    else

      info = 0

      do k = n, 1, -1

        if ( cabs1 ( t(k,k) ) == 0.0E+00 ) then
          info = k
          exit
        end if

        t(k,k) = cmplx ( 1.0E+00, 0.0E+00 ) / t(k,k)

        if ( k /= n ) then
          temp = -t(k,k)
          t(k+1:n,k) = t(k+1:n,k) * temp
        end if

        do j = 1, k-1
          temp = t(k,j)
          t(k,j) = cmplx ( 0.0E+00, 0.0E+00 )
          t(k:n,j) = t(k:n,j) + temp * t(k:n,k)
        end do

      end do

    end if

  end if

  return
end
subroutine ctrsl ( t, ldt, n, b, job, info )

!*****************************************************************************80
!
!! CTRSL solves triangular systems T*X=B or Hermitian(T)*X=B.
!
!  Discussion:
!
!    Hermitian ( T ) denotes the conjugate transpose of the matrix T.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
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
!    Input, complex ( kind = 4 ) T(LDT,N), the matrix of the system.  The zero
!    elements of the matrix are not referenced, and the corresponding
!    elements of the array can be used to store other information.
!
!    Input, integer ( kind = 4 ) LDT, the leading dimension of T.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, complex ( kind = 4 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
!    Input, integer ( kind = 4 ) JOB, specifies what kind of system is to
!    be solved.
!    00, solve T*X=B, T lower triangular,
!    01, solve T*X=B, T upper triangular,
!    10, solve hermitian(T)*X=B, T lower triangular,
!    11, solve hermitian(T)*X=B, T upper triangular.
!
!    Output, integer ( kind = 4 ) INFO.
!    0, the system is nonsingular.
!    K, the index of the first zero diagonal element of T.
!
  implicit none

  integer ( kind = 4 ) ldt
  integer ( kind = 4 ) n

  complex ( kind = 4 ) b(n)
  real    ( kind = 4 ) cabs1
  integer ( kind = 4 ) case
  complex ( kind = 4 ) cdotc
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) job
  complex ( kind = 4 ) t(ldt,n)
  complex ( kind = 4 ) temp
!
!  Check for zero diagonal elements.
!
  do info = 1, n
    if ( cabs1 ( t(info,info) ) == 0.0E+00 ) then
      return
    end if
  end do

  info = 0
!
!  Determine the task and go to it.
!
  case = 1

  if ( mod ( job, 10 ) /= 0 ) then
    case = 2
  end if

  if ( mod ( job, 100 ) / 10 /= 0 ) then
    case = case + 2
  end if
!
!  Solve T * X = B for T lower triangular.
!
  if ( case == 1 ) then

    b(1) = b(1) / t(1,1)

    do j = 2, n
      temp = -b(j-1)
      b(j:n) = b(j:n) + temp * t(j:n,j-1)
      b(j) = b(j) / t(j,j)
    end do
!
!  Solve T * X = B for T upper triangular.
!
  else if ( case == 2 ) then

    b(n) = b(n) / t(n,n)

    do jj = 2, n
      j = n - jj + 1
      temp = -b(j+1)
      b(1:j) = b(1:j) + temp * t(1:j,j+1)
      b(j) = b(j) / t(j,j)
    end do
!
!  Solve hermitian(T) * X = B for T lower triangular.
!
  else if ( case == 3 ) then

    b(n) = b(n) / conjg ( t(n,n) )

    do jj = 2, n
      j = n - jj + 1
      b(j) = b(j) - cdotc ( jj-1, t(j+1:j+jj-1,j), 1, b(j+1:j+jj-1), 1 )
      b(j) = b(j) / conjg ( t(j,j) )
    end do
!
!  Solve hermitian(T) * X = B for T upper triangular.
!
  else if ( case == 4 ) then

    b(1) = b(1) / conjg ( t(1,1) )

    do j = 2, n
      b(j) = b(j) - cdotc ( j-1, t(1:j-1,j), 1, b(1:j-1), 1 )
      b(j) = b(j) / conjg ( t(j,j) )
    end do

  end if

  return
end
