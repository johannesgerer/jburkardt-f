subroutine dchdc ( a, lda, p, work, ipvt, job, info )

!*****************************************************************************80
!
!! DCHDC computes the Cholesky decomposition of a positive definite matrix.
!
!  Discussion:
!
!    A pivoting option allows the user to estimate the condition of a
!    positive definite matrix or determine the rank of a positive
!    semidefinite matrix.
!
!    For positive definite matrices, INFO = P is the normal return.
!
!    For pivoting with positive semidefinite matrices, INFO will
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
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) A(LDA,P).
!    On input, A contains the matrix whose decomposition is to
!    be computed.  Only the upper half of A need be stored.
!    The lower part of the array a is not referenced.
!    On output, A contains in its upper half the Cholesky factor
!    of the input matrix, as it has been permuted by pivoting.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) P, the order of the matrix.
!
!    Input, real ( kind = 8 ) WORK(P) is a work array.
!
!    Input/output, integer ( kind = 4 ) IPVT(P).
!    On input, IPVT contains integers that control the selection
!    of the pivot elements, if pivoting has been requested.
!    Each diagonal element A(K,K) is placed in one of three classes
!    according to the value of IPVT(K).
!
!      > 0, then X(K) is an initial element.
!      = 0, then X(K) is a free element.
!      < 0, then X(K) is a final element.
!
!    Before the decomposition is computed, initial elements are moved by
!    symmetric row and column interchanges to the beginning of the array A
!    and final elements to the end.  Both initial and final elements are
!    frozen in place during the computation and only free elements are moved.
!    At the K-th stage of the reduction, if A(K,K) is occupied by a free
!    element, it is interchanged with the largest free element A(L,L) with
!    K <= L.  IPVT is not referenced if JOB is 0.
!
!    On output, IPVT(J) contains the index of the diagonal element
!    of A that was moved into the J-th position, if pivoting was requested.
!
!    Input, integer ( kind = 4 ) JOB, initiates column pivoting.
!    0, no pivoting is done.
!    nonzero, pivoting is done.
!
!    Output, integer ( kind = 4 ) INFO, contains the index of the last positive
!    diagonal element of the Cholesky factor.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) p

  real ( kind = 8 ) a(lda,p)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) ipvt(p)
  integer ( kind = 4 ) jt
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) maxdia
  integer ( kind = 4 ) maxl
  logical negk
  integer ( kind = 4 ) pl
  integer ( kind = 4 ) pu
  logical swapk
  real ( kind = 8 ) temp
  real ( kind = 8 ) work(p)

  pl = 1
  pu = 0
  info = p

  if ( job /= 0 ) then
!
!  Pivoting has been requested.
!  Rearrange the the elements according to IPVT.
!
    do k = 1, p

      swapk = 0 < ipvt(k)

      negk = ipvt(k) < 0

      if ( negk ) then
        ipvt(k) = -k
      else
        ipvt(k) = k
      end if

      if ( swapk ) then

        if ( k /= pl ) then

          call dswap ( pl-1, a(1,k), 1, a(1,pl), 1 )

          temp = a(k,k)
          a(k,k) = a(pl,pl)
          a(pl,pl) = temp

          do j = pl+1, p

            if ( j < k ) then
              temp = a(pl,j)
              a(pl,j) = a(j,k)
              a(j,k) = temp
            else if ( k < j ) then
              temp = a(k,j)
              a(k,j) = a(pl,j)
              a(pl,j) = temp
            end if

          end do

          ipvt(k) = ipvt(pl)
          ipvt(pl) = k

        end if

        pl = pl + 1

      end if

    end do

    pu = p

    do k = p, pl, -1

      if ( ipvt(k) < 0 ) then

        ipvt(k) = -ipvt(k)

        if ( pu /= k ) then

          call dswap ( k-1, a(1,k), 1, a(1,pu), 1 )

          temp = a(k,k)
          a(k,k) = a(pu,pu)
          a(pu,pu) = temp

          do j = k+1, p

            if ( j < pu ) then
              temp = a(k,j)
              a(k,j) = a(j,pu)
              a(j,pu) = temp
            else if ( pu < j ) then
              temp = a(k,j)
              a(k,j) = a(pu,j)
              a(pu,j) = temp
            end if

          end do

          jt = ipvt(k)
          ipvt(k) = ipvt(pu)
          ipvt(pu) = jt

        end if

        pu = pu - 1

      end if

    end do

  end if

  do k = 1, p
!
!  Reduction loop.
!
    maxdia = a(k,k)
    maxl = k
!
!  Determine the pivot element.
!
    if ( pl <= k .and. k < pu ) then

      do l = k+1, pu
        if ( maxdia < a(l,l) ) then
          maxdia = a(l,l)
          maxl = l
        end if
      end do

    end if
!
!  Quit if the pivot element is not positive.
!
    if ( maxdia <= 0.0D+00 ) then
      info = k - 1
      return
    end if
!
!  Start the pivoting and update IPVT.
!
    if ( k /= maxl ) then
      call dswap ( k-1, a(1,k), 1, a(1,maxl), 1 )
      a(maxl,maxl) = a(k,k)
      a(k,k) = maxdia
      jp = ipvt(maxl)
      ipvt(maxl) = ipvt(k)
      ipvt(k) = jp
    end if
!
!  Reduction step.
!  Pivoting is contained across the rows.
!
    work(k) = sqrt ( a(k,k) )
    a(k,k) = work(k)

    do j = k + 1, p

      if ( k /= maxl ) then

        if ( j < maxl ) then
          temp = a(k,j)
          a(k,j) = a(j,maxl)
          a(j,maxl) = temp
        else if ( maxl < j ) then
          temp = a(k,j)
          a(k,j) = a(maxl,j)
          a(maxl,j) = temp
        end if

      end if

      a(k,j) = a(k,j) / work(k)
      work(j) = a(k,j)
      temp = -a(k,j)
      call daxpy ( j-k, temp, work(k+1), 1, a(k+1,j), 1 )

    end do

  end do

  return
end
subroutine dchdd ( r, ldr, p, x, z, ldz, nz, y, rho, c, s, info )

!*****************************************************************************80
!
!! DCHDD downdates an augmented Cholesky decomposition.
!
!  Discussion:
!
!    DCHDD can also downdate the triangular factor of an augmented QR
!    decomposition.
!
!    Specifically, given an upper triangular matrix R of order P, a
!    row vector X, a column vector Z, and a scalar Y, DCHDD
!    determines an orthogonal matrix U and a scalar ZETA such that
!
!          (R   Z )     (RR  ZZ)
!      U * (      )  =  (      ),
!          (0 ZETA)     ( X   Y)
!
!    where RR is upper triangular.
!
!    If R and Z have been obtained from the factorization of a least squares
!    problem, then RR and ZZ are the factors corresponding to the problem
!    with the observation (X,Y) removed.  In this case, if RHO
!    is the norm of the residual vector, then the norm of
!    the residual vector of the downdated problem is
!    sqrt ( RHO * RHO - ZETA * ZETA ). DCHDD will simultaneously downdate
!    several triplets (Z, Y, RHO) along with R.
!
!    For a less terse description of what DCHDD does and how
!    it may be applied, see the LINPACK guide.
!
!    The matrix U is determined as the product U(1)*...*U(P)
!    where U(I) is a rotation in the (P+1,I)-plane of the form
!
!      ( C(I)      -S(I)    )
!      (                    ).
!      ( S(I)       C(I)    )
!
!    The rotations are chosen so that C(I) is real.
!
!    The user is warned that a given downdating problem may be impossible
!    to accomplish or may produce inaccurate results.  For example, this
!    can happen if X is near a vector whose removal will reduce the
!    rank of R.  Beware.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) R(LDR,P), the upper triangular matrix that
!    is to be  downdated.  The part of R below the diagonal is not referenced.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of the array R.
!    LDR must be at least P.
!
!    Input, integer ( kind = 4 ) P, the order of the matrix R.
!
!    Input, real ( kind = 8 ) X(P), the row vector that is to be removed from R.
!
!    Input/output, real ( kind = 8 ) Z(LDZ,NZ), an array of NZ P-vectors
!    which are to be downdated along with R.
!
!    Input, integer ( kind = 4 ) LDZ, the leading dimension of the array Z.
!    LDZ must be at least P.
!
!    Input, integer ( kind = 4 ) NZ, the number of vectors to be downdated.
!    NZ may be zero, in which case Z, Y, and RHO are not referenced.
!
!    Input, real ( kind = 8 ) Y(NZ), the scalars for the downdating of
!    the vectors Z.
!
!    Input/output, real ( kind = 8 ) RHO(NZ), the norms of the residual vectors.
!    On output these have been changed along with R and Z.
!
!    Output, real ( kind = 8 ) C(P), S(P), the cosines and sines of the
!    transforming rotations.
!
!    Output, integer ( kind = 4 ) INFO, return flag.
!     0, the entire downdating was successful.
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

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) azeta
  real ( kind = 8 ) b
  real ( kind = 8 ) c(p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  real ( kind = 8 ) norm
  real ( kind = 8 ) r(ldr,p)
  real ( kind = 8 ) rho(nz)
  real ( kind = 8 ) s(p)
  real ( kind = 8 ) scale
  real ( kind = 8 ) ddot
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) t
  real ( kind = 8 ) x(p)
  real ( kind = 8 ) xx
  real ( kind = 8 ) y(nz)
  real ( kind = 8 ) z(ldz,nz)
  real ( kind = 8 ) zeta
!
!  Solve R' * A = X, placing the result in the array S.
!
  info = 0
  s(1) = x(1) / r(1,1)

  do j = 2, p
    s(j) = x(j) - ddot ( j-1, r(1,j), 1, s, 1 )
    s(j) = s(j) / r(j,j)
  end do

  norm = dnrm2 ( p, s, 1 )

  if ( 1.0D+00 <= norm ) then
    info = -1
    return
  end if

  alpha = sqrt ( 1.0D+00 - norm * norm )
!
!  Determine the transformations.
!
  do ii = 1, p
    i = p - ii + 1
    scale = alpha + abs ( s(i) )
    a = alpha / scale
    b = s(i) / scale
    norm = sqrt ( a * a + b * b )
    c(i) = a / norm
    s(i) = b / norm
    alpha = scale * norm
  end do
!
!  Apply the transformations to R.
!
  do j = 1, p
    xx = 0.0D+00
    do ii = 1, j
      i = j - ii + 1
      t = c(i) * xx + s(i) * r(i,j)
      r(i,j) = c(i) * r(i,j) - s(i) * xx
      xx = t
    end do
  end do
!
!  If required, downdate Z and RHO.
!
  do j = 1, nz

    zeta = y(j)
    do i = 1, p
      z(i,j) = ( z(i,j) - s(i) * zeta ) / c(i)
      zeta = c(i) * zeta - s(i) * z(i,j)
    end do

    azeta = abs ( zeta )

    if ( rho(j) < azeta ) then
      info = 1
      rho(j) = -1.0D+00
    else
      rho(j) = rho(j) * sqrt ( 1.0D+00 - ( azeta / rho(j) )**2 )
    end if

  end do

  return
end
subroutine dchex ( r, ldr, p, k, l, z, ldz, nz, c, s, job )

!*****************************************************************************80
!
!! DCHEX updates the Cholesky factorization of a positive definite matrix.
!
!  Discussion:
!
!    The factorization has the form
!
!      A = R' * R
!
!    where A is a positive definite matrix of order P.
!
!    The updating involves diagonal permutations of the form
!
!      E' * A * E
!
!    where E is a permutation matrix.  Specifically, given
!    an upper triangular matrix R and a permutation matrix
!    E (which is specified by K, L, and JOB), DCHEX determines
!    an orthogonal matrix U such that
!
!      U * R * E = RR,
!
!    where RR is upper triangular.  At the user's option, the
!    transformation U will be multiplied into the array Z.
!    If A = X'*X, so that R is the triangular part of the
!    QR factorization of X, then RR is the triangular part of the
!    QR factorization of X*E, that is, X with its columns permuted.
!
!    For a less terse description of what DCHEX does and how
!    it may be applied, see the LINPACK guide.
!
!    The matrix Q is determined as the product U(L-K)*...*U(1)
!    of plane rotations of the form
!
!      (    C(I)       S(I) )
!      (                    ),
!      (   -S(I)       C(I) )
!
!    where C(I) is real, the rows these rotations operate on
!    are described below.
!
!    There are two types of permutations, which are determined
!    by the value of JOB.
!
!    1, right circular shift.  The columns are rearranged in the order:
!
!         1,...,K-1,L,K,K+1,...,L-1,L+1,...,P.
!
!       U is the product of L-K rotations U(I), where U(I)
!       acts in the (L-I,L-I+1)-plane.
!
!    2, left circular shift: the columns are rearranged in the order
!
!         1,...,K-1,K+1,K+2,...,L,K,L+1,...,P.
!
!       U is the product of L-K rotations U(I), where U(I)
!       acts in the (K+I-1,K+I)-plane.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) R(LDR,P).  On input, the upper
!    triangular factor that is to be updated.  Elements of R below the
!    diagonal are not referenced.  On output, R has been updated.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of the array R.
!    LDR must be at least P.
!
!    Input, integer ( kind = 4 ) P, the order of the matrix R.
!
!    Input, integer ( kind = 4 ) K, the first column to be permuted.
!
!    Input, integer ( kind = 4 ) L, the last column to be permuted.
!    L must be strictly greater than K.
!
!    Input/output real ( kind = 8 ) Z(LDZ,NZ), an array of NZ P-vectors into
!    which the transformation U is multiplied.  Z is not referenced if NZ = 0.
!    On output, Z has been updated.
!
!    Input, integer ( kind = 4 ) LDZ, the leading dimension of the array Z.
!    LDZ must be at least P.
!
!    Input, integer ( kind = 4 ) NZ, the number of columns of the matrix Z.
!
!    Output, real ( kind = 8 ) C(P), S(P), the cosines and sines of the
!    transforming rotations.
!
!    Input, integer ( kind = 4 ) JOB, determines the type of permutation.
!    1, right circular shift.
!    2, left circular shift.
!
  implicit none

  integer ( kind = 4 ) ldr
  integer ( kind = 4 ) ldz
  integer ( kind = 4 ) p
  integer ( kind = 4 ) nz

  real ( kind = 8 ) c(p)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) il
  integer ( kind = 4 ) iu
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) r(ldr,p)
  real ( kind = 8 ) s(p)
  real ( kind = 8 ) t
  real ( kind = 8 ) z(ldz,nz)
!
!  Right circular shift.
!
  if ( job == 1 ) then
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
      r(j+1,j+1) = 0.0D+00
    end do

    do i = 1, k-1
      ii = l - i + 1
      r(i,k) = s(ii)
    end do
!
!  Calculate the rotations.
!
    t = s(1)
    do i = 1, l - k
      call drotg ( s(i+1), t, c(i), s(i) )
      t = s(i+1)
    end do

    r(k,k) = t

    do j = k+1, p
      il = max ( 1, l-j+1 )
      do ii = il, l - k
        i = l - ii
        t = c(ii) * r(i,j) + s(ii) * r(i+1,j)
        r(i+1,j) = c(ii) * r(i+1,j) - s(ii) * r(i,j)
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
        z(i+1,j) = c(ii) * z(i+1,j) - s(ii) * z(i,j)
        z(i,j) = t
      end do
    end do
!
!  Left circular shift.
!
  else
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

    r(k+1:l,l) = 0.0D+00
!
!  Reduction loop.
!
    do j = k, p
!
!  Apply the rotations.
!
      if ( j /= k ) then

        iu = min ( j-1, l-1 )

        do i = k, iu
          ii = i - k + 1
          t = c(ii) * r(i,j) + s(ii) * r(i+1,j)
          r(i+1,j) = c(ii) * r(i+1,j) - s(ii) * r(i,j)
          r(i,j) = t
        end do

      end if

      if ( j < l ) then
        jj = j - k + 1
        t = s(jj)
        call drotg ( r(j,j), t, c(jj), s(jj) )
      end if

    end do
!
!  Apply the rotations to Z.
!
    do j = 1, nz
      do i = k, l - 1
        ii = i - k + 1
        t = c(ii) * z(i,j) + s(ii) * z(i+1,j)
        z(i+1,j) = c(ii) * z(i+1,j) - s(ii) * z(i,j)
        z(i,j) = t
      end do
    end do

  end if

  return
end
subroutine dchud ( r, ldr, p, x, z, ldz, nz, y, rho, c, s )

!*****************************************************************************80
!
!! DCHUD updates an augmented Cholesky decomposition.
!
!  Discussion:
!
!    DCHUD can also update the triangular part of an augmented QR
!    decomposition.
!
!    Specifically, given an upper triangular matrix R of order P, a row vector
!    X, a column vector Z, and a scalar Y, DCHUD determines a unitary matrix
!    U and a scalar ZETA such that
!
!           (R  Z)     (RR   ZZ )
!      U  * (    )  =  (        ),
!           (X  Y)     ( 0  ZETA)
!
!    where RR is upper triangular.
!
!    If R and Z have been obtained from the factorization of a least squares
!    problem, then RR and ZZ are the factors corresponding to the problem
!    with the observation (X,Y) appended.  In this case, if RHO is the
!    norm of the residual vector, then the norm of the residual vector of
!    the updated problem is sqrt ( RHO * RHO + ZETA * ZETA ).  DCHUD will
!    simultaneously update several triplets (Z, Y, RHO).
!
!    For a less terse description of what DCHUD does and how
!    it may be applied, see the LINPACK guide.
!
!    The matrix U is determined as the product U(P)*...*U(1),
!    where U(I) is a rotation in the (I,P+1) plane of the form
!
!      (     C(I)      S(I) )
!      (                    ).
!      (    -S(I)      C(I) )
!
!    The rotations are chosen so that C(I) is real.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) R(LDR,P), the upper triangular matrix
!    to be updated.  The part of R below the diagonal is not referenced.
!    On output, the matrix has been updated.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of the array R.
!    LDR must be at least equal to P.
!
!    Input, integer ( kind = 4 ) P, the order of the matrix R.
!
!    Input, real ( kind = 8 ) X(P), the row to be added to R.
!
!    Input/output, real ( kind = 8 ) Z(LDZ,NZ), contains NZ P-vectors
!    to be updated with R.
!
!    Input, integer ( kind = 4 ) LDZ, the leading dimension of the array Z.
!    LDZ must be at least P.
!
!    Input, integer ( kind = 4 ) NZ, the number of vectors to be updated.  NZ may be
!    zero, in which case Z, Y, and RHO are not referenced.
!
!    Input, real ( kind = 8 ) Y(NZ), the scalars for updating the vectors Z.
!
!    Input/output, real ( kind = 8 ) RHO(NZ).  On input, the norms of the
!    residual vectors to be updated.  If RHO(J) is negative, it is left
!    unaltered.
!
!    Output, real ( kind = 8 ) C(P), S(P), the cosines and sines of the
!    transforming rotations.
!
  implicit none

  integer ( kind = 4 ) ldr
  integer ( kind = 4 ) ldz
  integer ( kind = 4 ) nz
  integer ( kind = 4 ) p

  real ( kind = 8 ) azeta
  real ( kind = 8 ) c(p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r(ldr,p)
  real ( kind = 8 ) rho(nz)
  real ( kind = 8 ) s(p)
  real ( kind = 8 ) scale
  real ( kind = 8 ) t
  real ( kind = 8 ) x(p)
  real ( kind = 8 ) xj
  real ( kind = 8 ) y(nz)
  real ( kind = 8 ) z(ldz,nz)
  real ( kind = 8 ) zeta
!
!  Update R.
!
  do j = 1, p

    xj = x(j)
!
!  Apply the previous rotations.
!
    do i = 1, j-1
      t = c(i) * r(i,j) + s(i) * xj
      xj = c(i) * xj - s(i) * r(i,j)
      r(i,j) = t
    end do
!
!  Compute the next rotation.
!
    call drotg ( r(j,j), xj, c(j), s(j) )

  end do
!
!  If required, update Z and RHO.
!
  do j = 1, nz

    zeta = y(j)

    do i = 1, p
      t =    c(i) * z(i,j) + s(i) * zeta
      zeta = c(i) * zeta   - s(i) * z(i,j)
      z(i,j) = t
    end do

    azeta = abs ( zeta )

    if ( azeta /= 0.0D+00 .and. 0.0D+00 <= rho(j) ) then
      scale = azeta + rho(j)
      rho(j) = scale * sqrt ( ( azeta / scale )**2 + ( rho(j) / scale )**2 )
    end if

  end do

  return
end
subroutine dgbco ( abd, lda, n, ml, mu, ipvt, rcond, z )

!*****************************************************************************80
!
!! DGBCO factors a real band matrix and estimates its condition.
!
!  Discussion:
!
!    If RCOND is not needed, DGBFA is slightly faster.
!
!    To solve A*X = B, follow DGBCO by DGBSL.
!
!    To compute inverse(A)*C, follow DGBCO by DGBSL.
!
!    To compute determinant(A), follow DGBCO by DGBDI.
!
!  Example:
!
!    If the original matrix is
!
!      11 12 13  0  0  0
!      21 22 23 24  0  0
!       0 32 33 34 35  0
!       0  0 43 44 45 46
!       0  0  0 54 55 56
!       0  0  0  0 65 66
!
!    then for proper band storage,
!
!      N = 6, ML = 1, MU = 2, 5 <= LDA and ABD should contain
!
!       *  *  *  +  +  +      * = not used
!       *  * 13 24 35 46      + = used for pivoting
!       * 12 23 34 45 56
!      11 22 33 44 55 66
!      21 32 43 54 65  *
!
!  Band storage:
!
!    If A is a band matrix, the following program segment
!    will set up the input.
!
!      ml = (band width below the diagonal)
!      mu = (band width above the diagonal)
!      m = ml + mu + 1
!
!      do j = 1, n
!        i1 = max ( 1, j-mu )
!        i2 = min ( n, j+ml )
!        do i = i1, i2
!          k = i - j + m
!          abd(k,j) = a(i,j)
!        end do
!      end do
!
!    This uses rows ML+1 through 2*ML+MU+1 of ABD.  In addition, the first
!    ML rows in ABD are used for elements generated during the
!    triangularization.  The total number of rows needed in ABD is
!    2*ML+MU+1.  The ML+MU by ML+MU upper left triangle and the ML by ML
!    lower right triangle are not referenced.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) ABD(LDA,N).  On input, the matrix in band
!    storage.  The columns of the matrix are stored in the columns of ABD and
!    the diagonals of the matrix are stored in rows ML+1 through 2*ML+MU+1
!    of ABD.  On output, an upper triangular matrix in band storage and
!    the multipliers which were used to obtain it.  The factorization can
!    be written A = L*U where L is a product of permutation and unit lower
!    triangular matrices and U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array ABD.
!    2*ML + MU + 1 <= LDA is required.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, MU, the number of diagonals below and above the
!    main diagonal.  0 <= ML < N, 0 <= MU < N.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, real ( kind = 8 ) RCOND, an estimate of the reciprocal condition
!    of A.  For the system A*X = B, relative perturbations in A and B of size
!    EPSILON may cause relative perturbations in X of size EPSILON/RCOND.
!    If RCOND is so small that the logical expression
!      1.0D+00 + RCOND == 1.0D+00
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate underflows.
!
!    Workspace, real ( kind = 8 ) Z(N), a work vector whose contents are
!    usually unimportant.  If A is close to a singular matrix, then Z is an
!    approximate null vector in the sense that
!      norm(A*Z) = RCOND * norm(A) * norm(Z).
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) abd(lda,n)
  real ( kind = 8 ) anorm
  real ( kind = 8 ) ek
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
  real ( kind = 8 ) rcond
  real ( kind = 8 ) s
  real ( kind = 8 ) dasum
  real ( kind = 8 ) ddot
  real ( kind = 8 ) sm
  real ( kind = 8 ) t
  real ( kind = 8 ) wk
  real ( kind = 8 ) wkm
  real ( kind = 8 ) ynorm
  real ( kind = 8 ) z(n)
!
!  Compute the 1-norm of A.
!
  anorm = 0.0D+00
  l = ml + 1
  is = l + mu

  do j = 1, n

    anorm = max ( anorm, dasum ( l, abd(is,j), 1 ) )

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
!  Factor.
!
  call dgbfa ( abd, lda, n, ml, mu, ipvt, info )
!
!  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
!
!  Estimate = norm(Z)/norm(Y) where  a*z = y  and A'*Y = E.
!
!  A' is the transpose of A.  The components of E are
!  chosen to cause maximum local growth in the elements of W where
!  U'*W = E.  The vectors are frequently rescaled to avoid
!  overflow.
!
!  Solve U' * W = E.
!
  ek = 1.0D+00
  z(1:n) = 0.0D+00
  m = ml + mu + 1
  ju = 0

  do k = 1, n

    if ( z(k) /= 0.0D+00 ) then
      ek = sign ( ek, -z(k) )
    end if

    if ( abs ( abd(m,k) ) < abs ( ek - z(k) ) ) then
      s = abs ( abd(m,k) ) / abs ( ek - z(k) )
      z(1:n) = s * z(1:n)
      ek = s * ek
    end if

    wk = ek - z(k)
    wkm = -ek - z(k)
    s = abs ( wk )
    sm = abs ( wkm )

    if ( abd(m,k) /= 0.0D+00 ) then
      wk = wk / abd(m,k)
      wkm = wkm / abd(m,k)
    else
      wk = 1.0D+00
      wkm = 1.0D+00
    end if

    ju = min ( max ( ju, mu + ipvt(k) ), n )
    mm = m

    if ( k + 1 <= ju ) then

      do j = k + 1, ju
        mm = mm - 1
        sm = sm + abs ( z(j) + wkm * abd(mm,j) )
        z(j) = z(j) + wk * abd(mm,j)
        s = s + abs ( z(j) )
      end do

      if ( s < sm ) then
        t = wkm - wk
        wk = wkm
        mm = m
        do j = k + 1, ju
          mm = mm - 1
          z(j) = z(j) + t * abd(mm,j)
        end do
      end if

    end if

    z(k) = wk

  end do

  z(1:n) = z(1:n) / dasum ( n, z, 1 )
!
!  Solve L' * Y = W.
!
  do k = n, 1, -1

    lm = min ( ml, n - k )

    if ( k < m ) then
      z(k) = z(k) + ddot ( lm, abd(m+1,k), 1, z(k+1), 1 )
    end if

    if ( 1.0D+00 < abs ( z(k) ) ) then
      s = 1.0D+00 / abs ( z(k) )
      z(1:n) = s * z(1:n)
    end if

    l = ipvt(k)
    t = z(l)
    z(l) = z(k)
    z(k) = t

  end do

  z(1:n) = z(1:n) / dasum ( n, z, 1 )
  ynorm = 1.0D+00
!
!  Solve L * V = Y.
!
  do k = 1, n

    l = ipvt(k)
    t = z(l)
    z(l) = z(k)
    z(k) = t
    lm = min ( ml, n-k )

    if ( k < n ) then
      call daxpy ( lm, t, abd(m+1,k), 1, z(k+1), 1 )
    end if

    if ( 1.0D+00 < abs ( z(k) ) ) then
      s = 1.0D+00 / abs ( z(k) )
      z(1:n) = s * z(1:n)
      ynorm = s * ynorm
    end if

  end do

  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm
!
!  Solve U * Z = W.
!
  do k = n, 1, -1

    if ( abs ( abd(m,k) ) < abs ( z(k) ) ) then
      s = abs ( abd(m,k) ) / abs ( z(k) )
      z(1:n) = s * z(1:n)
      ynorm = s * ynorm
    end if

    if ( abd(m,k) /= 0.0D+00 ) then
      z(k) = z(k) / abd(m,k)
    else
      z(k) = 1.0D+00
    end if

    lm = min ( k, m ) - 1
    la = m - lm
    lz = k - lm
    t = -z(k)
    call daxpy ( lm, t, abd(la,k), 1, z(lz), 1 )

  end do
!
!  Make ZNORM = 1.0.
!
  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm

  if ( anorm /= 0.0D+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0D+00
  end if

  return
end
subroutine dgbdi ( abd, lda, n, ml, mu, ipvt, det )

!*****************************************************************************80
!
!! DGBDI computes the determinant of a band matrix factored by DGBCO or DGBFA.
!
!  Discussion:
!
!    If the inverse is needed, use DGBSL N times.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input, real ( kind = 8 ) ABD(LDA,N), the LU factor
!    information from DGBCO or DGBFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array ABD.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, MU, the number of diagonals below and above the
!    main diagonal.  0 <= ML < N, 0 <= MU < N.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from DGBCO or DGBFA.
!
!    Output, real ( kind = 8 ) DET(2), the determinant of the original matrix,
!    if requested.
!      determinant = DET(1) * 10.0**DET(2)
!    with  1.0D+00 <= abs ( DET(1) ) < 10.0D+00 or DET(1) = 0.0D+00.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) abd(lda,n)
  real ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu

  m = ml + mu + 1
  det(1) = 1.0D+00
  det(2) = 0.0D+00

  do i = 1, n

    if ( ipvt(i) /= i ) then
      det(1) = -det(1)
    end if

    det(1) = det(1) * abd(m,i)

    if ( det(1) == 0.0D+00 ) then
      return
    end if

    do while ( abs ( det(1) ) < 1.0D+00 )
      det(1) = det(1) * 10.0D+00
      det(2) = det(2) - 1.0D+00
    end do

    do while ( 10.0D+00 <= abs ( det(1) ) )
      det(1) = det(1) / 10.0D+00
      det(2) = det(2) + 1.0D+00
    end do

  end do

  return
end
subroutine dgbfa ( abd, lda, n, ml, mu, ipvt, info )

!*****************************************************************************80
!
!! DGBFA factors a real band matrix by elimination.
!
!  Discussion:
!
!    DGBFA is usually called by DGBCO, but it can be called
!    directly with a saving in time if RCOND is not needed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) ABD(LDA,N).  On input, the matrix in band
!    storage.  The columns of the matrix are stored in the columns of ABD
!    and the diagonals of the matrix are stored in rows ML+1 through
!    2*ML+MU+1 of ABD.  On output, an upper triangular matrix in band storage
!    and the multipliers which were used to obtain it.  The factorization
!    can be written A = L*U where L is a product of permutation and unit lower
!    triangular matrices and U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array ABD.
!    2*ML + MU + 1 <= LDA is required.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, MU, the number of diagonals below and above the
!    main diagonal.  0 <= ML < N, 0 <= MU < N.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, normal value.
!    K, if U(K,K) == 0.0D+00.  This is not an error condition for this
!      subroutine, but it does indicate that DGBSL will divide by zero if
!      called.  Use RCOND in DGBCO for a reliable indication of singularity.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) abd(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) idamax
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
  real ( kind = 8 ) t

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
      abd(i,jz) = 0.0D+00
    end do
  end do

  jz = j1
  ju = 0
!
!  Gaussian elimination with partial pivoting.
!
  do k = 1, n-1
!
!  Zero out the next fill-in column.
!
    jz = jz + 1
    if ( jz <= n ) then
      abd(1:ml,jz) = 0.0D+00
    end if
!
!  Find L = pivot index.
!
    lm = min ( ml, n-k )
    l = idamax ( lm+1, abd(m,k), 1 ) + m - 1
    ipvt(k) = l + k - m
!
!  Zero pivot implies this column already triangularized.
!
    if ( abd(l,k) == 0.0D+00 ) then

      info = k
!
!  Interchange if necessary.
!
    else

      if ( l /= m ) then
        t = abd(l,k)
        abd(l,k) = abd(m,k)
        abd(m,k) = t
      end if
!
!  Compute multipliers.
!
      t = -1.0D+00 / abd(m,k)
      call dscal ( lm, t, abd(m+1,k), 1 )
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
        call daxpy ( lm, t, abd(m+1,k), 1, abd(mm+1,j), 1 )
      end do

    end if

  end do

  ipvt(n) = n

  if ( abd(m,n) == 0.0D+00 ) then
    info = n
  end if

  return
end
subroutine dgbsl ( abd, lda, n, ml, mu, ipvt, b, job )

!*****************************************************************************80
!
!! DGBSL solves a real banded system factored by DGBCO or DGBFA.
!
!  Discussion:
!
!    DGBSL can solve either A * X = B  or  A' * X = B.
!
!    A division by zero will occur if the input factor contains a
!    zero on the diagonal.  Technically this indicates singularity
!    but it is often caused by improper arguments or improper
!    setting of LDA.  It will not occur if the subroutines are
!    called correctly and if DGBCO has set 0.0 < RCOND
!    or DGBFA has set INFO == 0.
!
!    To compute inverse(A) * C  where C is a matrix with P columns:
!
!      call dgbco ( abd, lda, n, ml, mu, ipvt, rcond, z )
!
!      if ( rcond is too small ) then
!        exit
!      end if
!
!      do j = 1, p
!        call dgbsl ( abd, lda, n, ml, mu, ipvt, c(1,j), 0 )
!      end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input, real ( kind = 8 ) ABD(LDA,N), the output from DGBCO or DGBFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array ABD.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, MU, the number of diagonals below and above the
!    main diagonal.  0 <= ML < N, 0 <= MU < N.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from DGBCO or DGBFA.
!
!    Input/output, real ( kind = 8 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
!    Input, integer ( kind = 4 ) JOB, job choice.
!    0, solve A*X=B.
!    nonzero, solve A'*X=B.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) abd(lda,n)
  real ( kind = 8 ) b(n)
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
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t

  m = mu + ml + 1
!
!  JOB = 0, Solve A * x = b.
!
!  First solve L * y = b.
!
  if ( job == 0 ) then

    if ( 0 < ml ) then

      do k = 1, n-1
        lm = min ( ml, n-k )
        l = ipvt(k)
        t = b(l)
        if ( l /= k ) then
          b(l) = b(k)
          b(k) = t
        end if
        call daxpy ( lm, t, abd(m+1,k), 1, b(k+1), 1 )
      end do

    end if
!
!  Now solve U * x = y.
!
    do k = n, 1, -1
      b(k) = b(k) / abd(m,k)
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm
      t = -b(k)
      call daxpy ( lm, t, abd(la,k), 1, b(lb), 1 )
    end do
!
!  JOB nonzero, solve A' * x = b.
!
!  First solve U' * y = b.
!
  else

    do k = 1, n
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm
      t = ddot ( lm, abd(la,k), 1, b(lb), 1 )
      b(k) = ( b(k) - t ) / abd(m,k)
    end do
!
!  Now solve L' * x = y.
!
    if ( 0 < ml ) then

      do k = n - 1, 1, -1
        lm = min ( ml, n - k )
        b(k) = b(k) + ddot ( lm, abd(m+1,k), 1, b(k+1), 1 )
        l = ipvt(k)
        if ( l /= k ) then
          t = b(l)
          b(l) = b(k)
          b(k) = t
        end if
      end do

    end if

  end if

  return
end
subroutine dgeco ( a, lda, n, ipvt, rcond, z )

!*****************************************************************************80
!
!! DGECO factors a real matrix and estimates its condition number.
!
!  Discussion:
!
!    If RCOND is not needed, DGEFA is slightly faster.
!
!    To solve A * X = B, follow DGECO by DGESL.
!
!    To compute inverse ( A ) * C, follow DGECO by DGESL.
!
!    To compute determinant ( A ), follow DGECO by DGEDI.
!
!    To compute inverse ( A ), follow DGECO by DGEDI.
!
!    For the system A * X = B, relative perturbations in A and B
!    of size EPSILON may cause relative perturbations in X of size
!    EPSILON/RCOND.
!
!    If RCOND is so small that the logical expression
!      1.0D+00 + RCOND == 1.0D+00
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate
!    underflows.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) A(LDA,N).  On input, a matrix to be
!    factored.  On output, the LU factorization of the matrix.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, real ( kind = 8 ) RCOND, an estimate of the reciprocal
!    condition number of A.
!
!    Output, real ( kind = 8 ) Z(N), a work vector whose contents are usually
!    unimportant.  If A is close to a singular matrix, then Z is an
!    approximate null vector in the sense that
!      norm ( A * Z ) = RCOND * norm ( A ) * norm ( Z ).
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) anorm
  real ( kind = 8 ) ek
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) rcond
  real ( kind = 8 ) s
  real ( kind = 8 ) sm
  real ( kind = 8 ) t
  real ( kind = 8 ) wk
  real ( kind = 8 ) wkm
  real ( kind = 8 ) ynorm
  real ( kind = 8 ) z(n)
!
!  Compute the L1 norm of A.
!
  anorm = 0.0D+00
  do j = 1, n
    anorm = max ( anorm, sum ( abs ( a(1:n,j) ) ) )
  end do
!
!  Compute the LU factorization.
!
  call dgefa ( a, lda, n, ipvt, info )
!
!  RCOND = 1 / ( norm(A) * (estimate of norm(inverse(A))) )
!
!  estimate of norm(inverse(A)) = norm(Z) / norm(Y)
!
!  where
!    A * Z = Y
!  and
!    A' * Y = E
!
!  The components of E are chosen to cause maximum local growth in the
!  elements of W, where U'*W = E.  The vectors are frequently rescaled
!  to avoid overflow.
!
!  Solve U' * W = E.
!
  ek = 1.0D+00
  z(1:n) = 0.0D+00

  do k = 1, n

    if ( z(k) /= 0.0D+00 ) then
      ek = sign ( ek, -z(k) )
    end if

    if ( abs ( a(k,k) ) < abs ( ek - z(k) ) ) then
      s = abs ( a(k,k) ) / abs ( ek - z(k) )
      z(1:n) = s * z(1:n)
      ek = s * ek
    end if

    wk = ek - z(k)
    wkm = -ek - z(k)
    s = abs ( wk )
    sm = abs ( wkm )

    if ( a(k,k) /= 0.0D+00 ) then
      wk = wk / a(k,k)
      wkm = wkm / a(k,k)
    else
      wk = 1.0D+00
      wkm = 1.0D+00
    end if

    if ( k+1 <= n ) then

      do j = k+1, n
        sm = sm + abs ( z(j) + wkm * a(k,j) )
        z(j) = z(j) + wk * a(k,j)
        s = s + abs ( z(j) )
      end do

      if ( s < sm ) then
        t = wkm - wk
        wk = wkm
        z(k+1:n) = z(k+1:n) + t * a(k,k+1:n)
      end if

    end if

    z(k) = wk

  end do

  z(1:n) = z(1:n) / sum ( abs ( z(1:n) ) )
!
!  Solve L' * Y = W
!
  do k = n, 1, -1

    z(k) = z(k) + dot_product ( a(k+1:n,k), z(k+1:n) )

    if ( 1.0D+00 < abs ( z(k) ) ) then
      z(1:n) = z(1:n) / abs ( z(k) )
    end if

    l = ipvt(k)

    t    = z(l)
    z(l) = z(k)
    z(k) = t

  end do

  z(1:n) = z(1:n) / sum ( abs ( z(1:n) ) )

  ynorm = 1.0D+00
!
!  Solve L * V = Y.
!
  do k = 1, n

    l = ipvt(k)

    t    = z(l)
    z(l) = z(k)
    z(k) = t

    z(k+1:n) = z(k+1:n) + t * a(k+1:n,k)

    if ( 1.0D+00 < abs ( z(k) ) ) then
      ynorm = ynorm / abs ( z(k) )
      z(1:n) = z(1:n) / abs ( z(k) )
    end if

  end do

  s = sum ( abs ( z(1:n) ) )
  z(1:n) = z(1:n) / s
  ynorm = ynorm / s
!
!  Solve U * Z = V.
!
  do k = n, 1, -1

    if ( abs ( a(k,k) ) < abs ( z(k) ) ) then
      s = abs ( a(k,k) ) / abs ( z(k) )
      z(1:n) = s * z(1:n)
      ynorm = s * ynorm
    end if

    if ( a(k,k) /= 0.0D+00 ) then
      z(k) = z(k) / a(k,k)
    else
      z(k) = 1.0D+00
    end if

    z(1:k-1) = z(1:k-1) - z(k) * a(1:k-1,k)

  end do
!
!  Normalize Z in the L1 norm.
!
  s = 1.0D+00 / sum ( abs ( z(1:n) ) )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm

  if ( anorm /= 0.0D+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0D+00
  end if

  return
end
subroutine dgedi ( a, lda, n, ipvt, det, work, job )

!*****************************************************************************80
!
!! DGEDI computes the determinant and inverse of a matrix factored by DGECO or DGEFA.
!
!  Discussion:
!
!    A division by zero will occur if the input factor contains
!    a zero on the diagonal and the inverse is requested.
!    It will not occur if the subroutines are called correctly
!    and if DGECO has set 0.0 < RCOND or DGEFA has set INFO == 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) A(LDA,N), on input, the  LU factor
!    information, as output by DGECO or DGEFA.  On output, the inverse
!    matrix, if requested.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from DGECO or DGEFA.
!
!    Workspace, real ( kind = 8 ) WORK(N).
!
!    Output, real ( kind = 8 ) DET(2), the determinant of original matrix if
!    requested.  The determinant = DET(1) * 10.0**DET(2)
!    with  1.0D+00 <= abs ( DET(1) ) < 10.0D+00
!    or DET(1) == 0.0D+00.
!
!    Input, integer ( kind = 4 ) JOB, specifies what is to be computed.
!    11, both determinant and inverse.
!    01, inverse only.
!    10, determinant only.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) t
  real ( kind = 8 ) work(n)
!
!  Compute the determinant.
!
  if ( job / 10 /= 0 ) then

    det(1) = 1.0D+00
    det(2) = 0.0D+00

    do i = 1, n

      if ( ipvt(i) /= i ) then
        det(1) = -det(1)
      end if

      det(1) = det(1) * a(i,i)

      if ( det(1) == 0.0D+00 ) then
        exit
      end if

      do while ( abs ( det(1) ) < 1.0D+00 )
        det(1) = det(1) * 10.0D+00
        det(2) = det(2) - 1.0D+00
      end do

      do while ( 10.0D+00 <= abs ( det(1) ) )
        det(1) = det(1) / 10.0D+00
        det(2) = det(2) + 1.0D+00
      end do

    end do

  end if
!
!  Compute inverse(U).
!
  if ( mod ( job, 10 ) /= 0 ) then

    do k = 1, n

      a(k,k) = 1.0D+00 / a(k,k)
      t = -a(k,k)
      call dscal ( k-1, t, a(1,k), 1 )

      do j = k+1, n
        t = a(k,j)
        a(k,j) = 0.0D+00
        call daxpy ( k, t, a(1,k), 1, a(1,j), 1 )
      end do

    end do
!
!  Form inverse(U) * inverse(L).
!
    do k = n-1, 1, -1

      work(k+1:n) = a(k+1:n,k)

      a(k+1:n,k) = 0.0D+00

      do j = k+1, n
        t = work(j)
        call daxpy ( n, t, a(1,j), 1, a(1,k), 1 )
      end do

      l = ipvt(k)
      if ( l /= k ) then
        call dswap ( n, a(1,k), 1, a(1,l), 1 )
      end if

    end do

  end if

  return
end
subroutine dgefa ( a, lda, n, ipvt, info )

!*****************************************************************************80
!
!! DGEFA factors a real general matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2001
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
!    Input/output, real ( kind = 8 ) A(LDA,N).
!    On intput, the matrix to be factored.
!    On output, an upper triangular matrix and the multipliers used to obtain
!    it.  The factorization can be written A=L*U, where L is a product of
!    permutation and unit lower triangular matrices, and U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, singularity indicator.
!    0, normal value.
!    K, if U(K,K) == 0.  This is not an error condition for this subroutine,
!    but it does indicate that DGESL or DGEDI will divide by zero if called.
!    Use RCOND in DGECO for a reliable indication of singularity.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) t
!
!  Gaussian elimination with partial pivoting.
!
  info = 0

  do k = 1, n - 1
!
!  Find L = pivot index.
!
    l = idamax ( n-k+1, a(k,k), 1 ) + k - 1
    ipvt(k) = l
!
!  Zero pivot implies this column already triangularized.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      cycle
    end if
!
!  Interchange if necessary.
!
    if ( l /= k ) then
      t = a(l,k)
      a(l,k) = a(k,k)
      a(k,k) = t
    end if
!
!  Compute multipliers.
!
    t = -1.0D+00 / a(k,k)
    call dscal ( n-k, t, a(k+1,k), 1 )
!
!  Row elimination with column indexing.
!
    do j = k+1, n
      t = a(l,j)
      if ( l /= k ) then
        a(l,j) = a(k,j)
        a(k,j) = t
      end if
      call daxpy ( n-k, t, a(k+1,k), 1, a(k+1,j), 1 )
    end do

  end do

  ipvt(n) = n

  if ( a(n,n) == 0.0D+00 ) then
    info = n
  end if

  return
end
subroutine dgesl ( a, lda, n, ipvt, b, job )

!*****************************************************************************80
!
!! DGESL solves a real general linear system A * X = B.
!
!  Discussion:
!
!    DGESL can solve either of the systems A * X = B or A' * X = B.
!
!    The system matrix must have been factored by DGECO or DGEFA.
!
!    A division by zero will occur if the input factor contains a
!    zero on the diagonal.  Technically this indicates singularity
!    but it is often caused by improper arguments or improper
!    setting of LDA.  It will not occur if the subroutines are
!    called correctly and if DGECO has set 0.0 < RCOND
!    or DGEFA has set INFO == 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2001
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
!    Input, real ( kind = 8 ) A(LDA,N), the output from DGECO or DGEFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from DGECO or DGEFA.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer ( kind = 4 ) JOB.
!    0, solve A * X = B;
!    nonzero, solve A' * X = B.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t
!
!  Solve A * X = B.
!
  if ( job == 0 ) then

    do k = 1, n-1

      l = ipvt(k)
      t = b(l)

      if ( l /= k ) then
        b(l) = b(k)
        b(k) = t
      end if

      call daxpy ( n-k, t, a(k+1,k), 1, b(k+1), 1 )

    end do

    do k = n, 1, -1
      b(k) = b(k) / a(k,k)
      t = -b(k)
      call daxpy ( k-1, t, a(1,k), 1, b(1), 1 )
    end do

  else
!
!  Solve A' * X = B.
!
    do k = 1, n
      t = ddot ( k-1, a(1,k), 1, b(1), 1 )
      b(k) = ( b(k) - t ) / a(k,k)
    end do

    do k = n-1, 1, -1

      b(k) = b(k) + ddot ( n-k, a(k+1,k), 1, b(k+1), 1 )
      l = ipvt(k)

      if ( l /= k ) then
        t = b(l)
        b(l) = b(k)
        b(k) = t
      end if

    end do

  end if

  return
end
subroutine dgtsl ( n, c, d, e, b, info )

!*****************************************************************************80
!
!! DGTSL solves a general tridiagonal linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input, integer ( kind = 4 ) N, the order of the tridiagonal matrix.
!
!    Input/output, real ( kind = 8 ) C(N), contains the subdiagonal of the
!    tridiagonal matrix in entries C(2:N).  On output, C is destroyed.
!
!    Input/output, real ( kind = 8 ) D(N).  On input, the diagonal of the
!    matrix.  On output, D is destroyed.
!
!    Input/output, real ( kind = 8 ) E(N), contains the superdiagonal of the
!    tridiagonal matrix in entries E(1:N-1).  On output E is destroyed.
!
!    Input/output, real ( kind = 8 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, normal value.
!    K, the K-th element of the diagonal becomes exactly zero.  The
!       subroutine returns if this error condition is detected.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) k
  real ( kind = 8 ) t

  info = 0
  c(1) = d(1)

  if ( 2 <= n ) then

    d(1) = e(1)
    e(1) = 0.0D+00
    e(n) = 0.0D+00

    do k = 1, n - 1
!
!  Find the larger of the two rows.
!
      if ( abs ( c(k) ) <= abs ( c(k+1) ) ) then
!
!  Interchange rows.
!
        t = c(k+1)
        c(k+1) = c(k)
        c(k) = t

        t = d(k+1)
        d(k+1) = d(k)
        d(k) = t

        t = e(k+1)
        e(k+1) = e(k)
        e(k) = t

        t = b(k+1)
        b(k+1) = b(k)
        b(k) = t

      end if
!
!  Zero elements.
!
      if ( c(k) == 0.0D+00 ) then
        info = k
        return
      end if

      t = -c(k+1) / c(k)
      c(k+1) = d(k+1) + t * d(k)
      d(k+1) = e(k+1) + t * e(k)
      e(k+1) = 0.0D+00
      b(k+1) = b(k+1) + t * b(k)

    end do

  end if

  if ( c(n) == 0.0D+00 ) then
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
subroutine dpbco ( abd, lda, n, m, rcond, z, info )

!*****************************************************************************80
!
!! DPBCO factors a real symmetric positive definite banded matrix.
!
!  Discussion:
!
!    DPBCO also estimates the condition of the matrix.
!
!    If RCOND is not needed, DPBFA is slightly faster.
!
!    To solve A*X = B, follow DPBCO by DPBSL.
!
!    To compute inverse(A)*C, follow DPBCO by DPBSL.
!
!    To compute determinant(A), follow DPBCO by DPBDI.
!
!  Band storage:
!
!    If A is a symmetric positive definite band matrix, the following
!    program segment will set up the input.
!
!      m = (band width above diagonal)
!      do j = 1, n
!        i1 = max (1, j-m)
!        do i = i1, j
!          k = i-j+m+1
!          abd(k,j) = a(i,j)
!        end do
!      end do
!
!    This uses M + 1 rows of A, except for the M by M upper left triangle,
!    which is ignored.
!
!    For example, if the original matrix is
!
!      11 12 13  0  0  0
!      12 22 23 24  0  0
!      13 23 33 34 35  0
!       0 24 34 44 45 46
!       0  0 35 45 55 56
!       0  0  0 46 56 66
!
!    then N = 6, M = 2  and ABD should contain
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
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) ABD(LDA,N).  On input, the matrix to be
!    factored.  The columns of the upper triangle are stored in the columns
!    of ABD and the diagonals of the upper triangle are stored in the rows
!    of ABD.  On output, an upper triangular matrix R, stored in band form,
!    so that A = R'*R.  If INFO /= 0, the factorization is not complete.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array ABD.
!    M+1 <= LDA is required.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) M, the number of diagonals above the main diagonal.
!
!    Output, real ( kind = 8 ) RCOND, an estimate of the reciprocal condition
!    of A.  For the system A*X = B, relative perturbations in A and B of size
!    EPSILON may cause relative perturbations in X of size EPSILON/RCOND.
!    If RCOND is so small that the logical expression
!      1.0D+00 + RCOND == 1.0D+00
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate underflows.
!
!    Output, real ( kind = 8 ) Z(N), a work vector whose contents are usually
!    unimportant.  If A is singular to working precision, then Z is an
!    approximate null vector in the sense that
!      norm(A*Z) = RCOND * norm(A) * norm(Z).
!    If INFO /= 0, Z is unchanged.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, for normal return.
!    K, signals an error condition.  The leading minor of order K is not
!    positive definite.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) abd(lda,n)
  real ( kind = 8 ) anorm
  real ( kind = 8 ) ek
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) lm
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mu
  real ( kind = 8 ) rcond
  real ( kind = 8 ) s
  real ( kind = 8 ) dasum
  real ( kind = 8 ) ddot
  real ( kind = 8 ) sm
  real ( kind = 8 ) t
  real ( kind = 8 ) wk
  real ( kind = 8 ) wkm
  real ( kind = 8 ) ynorm
  real ( kind = 8 ) z(n)
!
!  Find the norm of A.
!
  do j = 1, n

    l = min ( j, m+1 )
    mu = max ( m+2-j, 1 )
    z(j) = dasum ( l, abd(mu,j), 1 )
    k = j - l
    do i = mu, m
      k = k + 1
      z(k) = z(k) + abs ( abd(i,j) )
    end do

  end do

  anorm = maxval ( z(1:n) )
!
!  Factor.
!
  call dpbfa ( abd, lda, n, m, info )

  if ( info /= 0 ) then
    return
  end if
!
!  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
!
!  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
!
!  The components of E are chosen to cause maximum local
!  growth in the elements of W where R'*W = E.
!
!  The vectors are frequently rescaled to avoid overflow.
!
!  Solve R' * W = E.
!
  ek = 1.0D+00
  z(1:n) = 0.0D+00

  do k = 1, n

    if ( z(k) /= 0.0D+00 ) then
      ek = sign ( ek, -z(k) )
    end if

    if ( abd(m+1,k) < abs ( ek - z(k) ) ) then
      s = abd(m+1,k) / abs ( ek - z(k) )
      z(1:n) = s * z(1:n)
      ek = s * ek
    end if

    wk = ek - z(k)
    wkm = -ek - z(k)
    s = abs ( wk )
    sm = abs ( wkm )
    wk = wk / abd(m+1,k)
    wkm = wkm / abd(m+1,k)
    j2 = min ( k+m, n )
    i = m + 1

    if ( k+1 <= j2 ) then

      do j = k+1, j2
        i = i - 1
        sm = sm + abs ( z(j) + wkm * abd(i,j) )
        z(j) = z(j) + wk * abd(i,j)
        s = s + abs ( z(j) )
      end do

      if ( s < sm ) then

        t = wkm - wk
        wk = wkm
        i = m + 1

        do j = k+1, j2
          i = i - 1
          z(j) = z(j) + t * abd(i,j)
        end do

      end if

    end if

    z(k) = wk

  end do

  z(1:n) = z(1:n) / dasum ( n, z, 1 )
!
!  Solve R * Y = W.
!
  do k = n, 1, -1

    if ( abd(m+1,k) < abs ( z(k) ) ) then
      s = abd(m+1,k) / abs ( z(k) )
      z(1:n) = s * z(1:n)
    end if

    z(k) = z(k) / abd(m+1,k)
    lm = min ( k-1, m )
    la = m + 1 - lm
    lb = k - lm
    t = -z(k)
    call daxpy ( lm, t, abd(la,k), 1, z(lb), 1 )

  end do

  z(1:n) = z(1:n) / dasum ( n, z, 1 )

  ynorm = 1.0D+00
!
!  Solve R' * V = Y.
!
  do k = 1, n

    lm = min ( k-1, m )
    la = m + 1 - lm
    lb = k - lm

    z(k) = z(k) - ddot ( lm, abd(la,k), 1, z(lb), 1 )

    if ( abd(m+1,k) < abs ( z(k) ) ) then
      s = abd(m+1,k) / abs ( z(k) )
      z(1:n) = s * z(1:n)
      ynorm = s * ynorm
    end if

    z(k) = z(k) / abd(m+1,k)

  end do

  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm
!
!  Solve R * Z = W.
!
  do k = n, 1, -1

    if ( abd(m+1,k) < abs ( z(k) ) ) then
      s = abd(m+1,k) / abs ( z(k) )
      z(1:n) = s * z(1:n)
      ynorm = s * ynorm
    end if

    z(k) = z(k) / abd(m+1,k)
    lm = min ( k-1, m )
    la = m + 1 - lm
    lb = k - lm
    t = -z(k)
    call daxpy ( lm, t, abd(la,k), 1, z(lb), 1 )

  end do
!
!  Make ZNORM = 1.0.
!
  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm

  if ( anorm /= 0.0D+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0D+00
  end if

  return
end
subroutine dpbdi ( abd, lda, n, m, det )

!*****************************************************************************80
!
!! DPBDI computes the determinant of a matrix factored by DPBCO or DPBFA.
!
!  Discussion:
!
!    If the inverse is needed, use DPBSL N times.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input, real ( kind = 8 ) ABD(LDA,N), the output from DPBCO or DPBFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array ABD.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) M, the number of diagonals above the main diagonal.
!
!    Output, real ( kind = 8 ) DET(2), the determinant of the original
!    matrix in the form
!      determinant = DET(1) * 10.0**DET(2)
!    with 1.0D+00 <= DET(1) < 10.0D+00 or DET(1) == 0.0D+00.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) abd(lda,n)
  real ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
!
!  Compute the determinant.
!
  det(1) = 1.0D+00
  det(2) = 0.0D+00

  do i = 1, n

    det(1) = det(1) * abd(m+1,i) * abd(m+1,i)

    if ( det(1) == 0.0D+00 ) then
      return
    end if

    do while ( det(1) < 1.0D+00 )
      det(1) = det(1) * 10.0D+00
      det(2) = det(2) - 1.0D+00
    end do

    do while ( 10.0D+00 <= det(1) )
      det(1) = det(1) / 10.0D+00
      det(2) = det(2) + 1.0D+00
    end do

  end do

  return
end
subroutine dpbfa ( abd, lda, n, m, info )

!*****************************************************************************80
!
!! DPBFA factors a real symmetric positive definite matrix stored in band form.
!
!  Discussion:
!
!    DPBFA is usually called by DPBCO, but it can be called
!    directly with a saving in time if RCOND is not needed.
!
!    If A is a symmetric positive definite band matrix,
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
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) ABD(LDA,N).  On input, the matrix to be
!    factored.  The columns of the upper triangle are stored in the columns
!    of ABD and the diagonals of the upper triangle are stored in the
!    rows of ABD.  On output, an upper triangular matrix R, stored in band
!    form, so that A = R' * R.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array ABD.
!    M+1 <= LDA is required.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) M, the number of diagonals above the main diagonal.
!
!    Output, integer ( kind = 4 ) INFO, error indicator.
!    0, for normal return.
!    K, if the leading minor of order K is not positive definite.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) abd(lda,n)
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mu
  real ( kind = 8 ) ddot
  real ( kind = 8 ) s
  real ( kind = 8 ) t

  do j = 1, n

    s = 0.0D+00
    ik = m + 1
    jk = max ( j - m, 1 )
    mu = max ( m + 2 - j, 1 )

    do k = mu, m
      t = abd(k,j) - ddot ( k-mu, abd(ik,jk), 1, abd(mu,j), 1 )
      t = t / abd(m+1,jk)
      abd(k,j) = t
      s = s + t * t
      ik = ik - 1
      jk = jk + 1
    end do

    s = abd(m+1,j) - s

    if ( s <= 0.0D+00 ) then
      info = j
      return
    end if

    abd(m+1,j) = sqrt ( s )

  end do

  info = 0

  return
end
subroutine dpbsl ( abd, lda, n, m, b )

!*****************************************************************************80
!
!! DPBSL solves a real SPD band system factored by DPBCO or DPBFA.
!
!  Discussion:
!
!    The matrix is assumed to be a symmetric positive definite (SPD)
!    band matrix.
!
!    To compute inverse(A) * C  where C is a matrix with P columns:
!
!      call dpbco ( abd, lda, n, rcond, z, info )
!
!      if ( rcond is too small .or. info /= 0) go to ...
!
!      do j = 1, p
!        call dpbsl ( abd, lda, n, c(1,j) )
!      end do
!
!    A division by zero will occur if the input factor contains
!    a zero on the diagonal.  Technically this indicates
!    singularity but it is usually caused by improper subroutine
!    arguments.  It will not occur if the subroutines are called
!    correctly and INFO == 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input, real ( kind = 8 ) ABD(LDA,N), the output from DPBCO or DPBFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array ABD.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) M, the number of diagonals above the main diagonal.
!
!    Input/output, real ( kind = 8 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) abd(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) lm
  integer ( kind = 4 ) m
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t
!
!  Solve R' * Y = B.
!
  do k = 1, n
    lm = min ( k-1, m )
    la = m + 1 - lm
    lb = k - lm
    t = ddot ( lm, abd(la,k), 1, b(lb), 1 )
    b(k) = ( b(k) - t ) / abd(m+1,k)
  end do
!
!  Solve R * X = Y.
!
  do k = n, 1, -1
    lm = min ( k-1, m )
    la = m + 1 - lm
    lb = k - lm
    b(k) = b(k) / abd(m+1,k)
    t = -b(k)
    call daxpy ( lm, t, abd(la,k), 1, b(lb), 1 )
  end do

  return
end
subroutine dpoco ( a, lda, n, rcond, z, info )

!*****************************************************************************80
!
!! DPOCO factors a real SPD matrix and estimates its condition.
!
!  Discussion:
!
!    If RCOND is not needed, DPOFA is slightly faster.
!
!    To solve A*X = B, follow DPOCO by DPOSL.
!
!    To compute inverse(A)*C, follow DPOCO by DPOSL.
!
!    To compute determinant(A), follow DPOCO by DPODI.
!
!    To compute inverse(A), follow DPOCO by DPODI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) A(LDA,N).  On input, the symmetric
!    matrix to be factored.  Only the diagonal and upper triangle are used.
!    On output, an upper triangular matrix R so that A = R'*R where R'
!    is the transpose.  The strict lower triangle is unaltered.
!    If INFO /= 0, the factorization is not complete.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) RCOND, an estimate of the reciprocal
!    condition of A.  For the system A*X = B, relative perturbations in
!    A and B of size EPSILON may cause relative perturbations in X of
!    size EPSILON/RCOND.  If RCOND is so small that the logical expression
!      1.0D+00 + RCOND == 1.0D+00
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate underflows.
!
!    Output, real ( kind = 8 ) Z(N), a work vector whose contents are usually
!    unimportant.  If A is close to a singular matrix, then Z is an
!    approximate null vector in the sense that
!      norm(A*Z) = RCOND * norm(A) * norm(Z).
!    If INFO /= 0, Z is unchanged.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, for normal return.
!    K, signals an error condition.  The leading minor of order K is not
!    positive definite.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) anorm
  real ( kind = 8 ) ek
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) rcond
  real ( kind = 8 ) s
  real ( kind = 8 ) dasum
  real ( kind = 8 ) ddot
  real ( kind = 8 ) sm
  real ( kind = 8 ) t
  real ( kind = 8 ) wk
  real ( kind = 8 ) wkm
  real ( kind = 8 ) ynorm
  real ( kind = 8 ) z(n)
!
!  Find norm of A using only upper half.
!
  do j = 1, n
    z(j) = dasum ( j, a(1,j), 1 )
    do i = 1, j-1
      z(i) = z(i) + abs ( a(i,j) )
    end do
  end do

  anorm = maxval ( z(1:n) )
!
!  Factor.
!
  call dpofa ( a, lda, n, info )

  if ( info /= 0 ) then
    return
  end if
!
!  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
!
!  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
!
!  The components of E are chosen to cause maximum local
!  growth in the elements of W where R'*W = E.
!
!  The vectors are frequently rescaled to avoid overflow.
!
!  Solve R' * W = E.
!
  ek = 1.0D+00
  z(1:n) = 0.0D+00

  do k = 1, n

    if ( z(k) /= 0.0D+00 ) then
      ek = sign ( ek, -z(k) )
    end if

    if ( a(k,k) < abs ( ek - z(k) ) ) then
      s = a(k,k) / abs ( ek - z(k) )
      z(1:n) = s * z(1:n)
      ek = s * ek
    end if

    wk = ek - z(k)
    wkm = -ek - z(k)
    s = abs ( wk )
    sm = abs ( wkm )
    wk = wk / a(k,k)
    wkm = wkm / a(k,k)

    if ( k + 1 <= n ) then

      do j = k+1, n
        sm = sm + abs ( z(j) + wkm * a(k,j) )
        z(j) = z(j) + wk * a(k,j)
        s = s + abs ( z(j) )
      end do

      if ( s < sm ) then
        t = wkm - wk
        wk = wkm
        z(k+1:n) = z(k+1:n) + t * a(k,k+1:n)
      end if

    end if

    z(k) = wk

  end do

  z(1:n) = z(1:n) / dasum ( n, z, 1 )
!
!  Solve R * Y = W.
!
  do k = n, 1, -1

    if ( a(k,k) < abs ( z(k) ) ) then
      s = a(k,k) / abs ( z(k) )
      z(1:n) = s * z(1:n)
    end if

    z(k) = z(k) / a(k,k)
    t = -z(k)
    call daxpy ( k-1, t, a(1,k), 1, z(1), 1 )

  end do

  z(1:n) = z(1:n) / dasum ( n, z, 1 )
  ynorm = 1.0D+00
!
!  Solve R' * V = Y.
!
  do k = 1, n

    z(k) = z(k) - ddot ( k-1, a(1,k), 1, z(1), 1 )

    if ( a(k,k) < abs ( z(k) ) ) then
      s = a(k,k) / abs ( z(k) )
      z(1:n) = s * z(1:n)
      ynorm = s * ynorm
    end if

    z(k) = z(k) / a(k,k)

  end do

  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm
!
!  Solve R * Z = V.
!
  do k = n, 1, -1

    if ( a(k,k) < abs ( z(k) ) ) then
      s = a(k,k) / abs ( z(k) )
      z(1:n) = s * z(1:n)
      ynorm = s * ynorm
    end if

    z(k) = z(k) / a(k,k)
    t = -z(k)
    call daxpy ( k-1, t, a(1,k), 1, z(1), 1 )

  end do
!
!  Make ZNORM = 1.0.
!
  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm

  if ( anorm /= 0.0D+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0D+00
  end if

  return
end
subroutine dpodi ( a, lda, n, det, job )

!*****************************************************************************80
!
!! DPODI computes the determinant and inverse of a certain matrix.
!
!  Discussion:
!
!    The matrix is real symmetric positive definite.
!    DPODI uses the factors computed by DPOCO, DPOFA or DQRDC.
!
!    A division by zero will occur if the input factor contains
!    a zero on the diagonal and the inverse is requested.
!    It will not occur if the subroutines are called correctly
!    and if DPOCO or DPOFA has set INFO == 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) A(LDA,N).  On input, the output A from
!    DPOCO or DPOFA, or the output X from DQRDC.  On output, if DPOCO or
!    DPOFA was used to factor A then DPODI produces the upper half of
!    inverse(A).  If DQRDC was used to decompose X then DPODI produces
!    the upper half of inverse(X'*X) where X' is the transpose.
!    Elements of A below the diagonal are unchanged.  If the units digit
!    of JOB is zero, A is unchanged.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, integer ( kind = 4 ) JOB, specifies the task.
!    11, both determinant and inverse.
!    01, inverse only.
!    10, determinant only.
!
!    Output, real ( kind = 8 ) DET(2), the determinant of A or of X'*X
!    if requested.
!      determinant = DET(1) * 10.0**DET(2)
!    with 1.0D+00 <= DET(1) < 10.0D+00 or DET(1) == 0.0D+00.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  real ( kind = 8 ) t
!
!  Compute the determinant.
!
  if ( job / 10 /= 0 ) then

    det(1) = 1.0D+00
    det(2) = 0.0D+00

    do i = 1, n

      det(1) = det(1) * a(i,i) * a(i,i)

      if ( det(1) == 0.0D+00 ) then
        exit
      end if

      do while ( det(1) < 1.0D+00 )
        det(1) = det(1) * 10.0D+00
        det(2) = det(2) - 1.0D+00
      end do

      do while ( 10.0D+00 <= det(1) )
        det(1) = det(1) / 10.0D+00
        det(2) = det(2) + 1.0D+00
      end do

    end do

  end if
!
!  Compute inverse(R).
!
  if ( mod ( job, 10 ) /= 0 ) then

    do k = 1, n

      a(k,k) = 1.0D+00 / a(k,k)
      t = -a(k,k)
      call dscal ( k-1, t, a(1,k), 1 )

      do j = k+1, n
        t = a(k,j)
        a(k,j) = 0.0D+00
        call daxpy ( k, t, a(1,k), 1, a(1,j), 1 )
      end do

    end do
!
!  Form inverse(R) * (inverse(R))'.
!
    do j = 1, n
      do k = 1, j-1
        t = a(k,j)
        call daxpy ( k, t, a(1,j), 1, a(1,k), 1 )
      end do
      t = a(j,j)
      call dscal ( j, t, a(1,j), 1 )
    end do

  end if

  return
end
subroutine dpofa ( a, lda, n, info )

!*****************************************************************************80
!
!! DPOFA factors a real symmetric positive definite matrix.
!
!  Discussion:
!
!    DPOFA is usually called by DPOCO, but it can be called
!    directly with a saving in time if RCOND is not needed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) A(LDA,N).  On input, the symmetric matrix
!    to be  factored.  Only the diagonal and upper triangle are used.
!    On output, an upper triangular matrix R so that A = R'*R
!    where R' is the transpose.  The strict lower triangle is unaltered.
!    If INFO /= 0, the factorization is not complete.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, for normal return.
!    K, signals an error condition.  The leading minor of order K is not
!    positive definite.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) s
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t

  do j = 1, n

    s = 0.0D+00

    do k = 1, j-1
      t = a(k,j) - ddot ( k-1, a(1,k), 1, a(1,j), 1 )
      t = t / a(k,k)
      a(k,j) = t
      s = s + t * t
    end do

    s = a(j,j) - s

    if ( s <= 0.0D+00 ) then
      info = j
      return
    end if

    a(j,j) = sqrt ( s )

  end do

  info = 0

  return
end
subroutine dposl ( a, lda, n, b )

!*****************************************************************************80
!
!! DPOSL solves a linear system factored by DPOCO or DPOFA.
!
!  Discussion:
!
!    To compute inverse(A) * C where C is a matrix with P columns:
!
!      call dpoco ( a, lda, n, rcond, z, info )
!
!      if ( rcond is not too small .and. info == 0 ) then
!        do j = 1, p
!          call dposl ( a, lda, n, c(1,j) )
!        end do
!      end if
!
!    A division by zero will occur if the input factor contains
!    a zero on the diagonal.  Technically this indicates
!    singularity but it is usually caused by improper subroutine
!    arguments.  It will not occur if the subroutines are called
!    correctly and INFO == 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input, real ( kind = 8 ) A(LDA,N), the output from DPOCO or DPOFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) k
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t
!
!  Solve R' * Y = B.
!
  do k = 1, n
    t = ddot ( k-1, a(1,k), 1, b(1), 1 )
    b(k) = ( b(k) - t ) / a(k,k)
  end do
!
!  Solve R * X = Y.
!
  do k = n, 1, -1
    b(k) = b(k) / a(k,k)
    t = -b(k)
    call daxpy ( k-1, t, a(1,k), 1, b(1), 1 )
  end do

  return
end
subroutine dppco ( ap, n, rcond, z, info )

!*****************************************************************************80
!
!! DPPCO factors a real symmetric positive definite matrix in packed form.
!
!  Discussion:
!
!    DPPCO also estimates the condition of the matrix.
!
!    If RCOND is not needed, DPPFA is slightly faster.
!
!    To solve A*X = B, follow DPPCO by DPPSL.
!
!    To compute inverse(A)*C, follow DPPCO by DPPSL.
!
!    To compute determinant(A), follow DPPCO by DPPDI.
!
!    To compute inverse(A), follow DPPCO by DPPDI.
!
!  Packed storage:
!
!    The following program segment will pack the upper triangle of
!    a symmetric matrix.
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
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) AP((N*(N+1))/2).  On input, the packed
!    form of a symmetric matrix A.  The columns of the upper triangle are
!    stored sequentially in a one-dimensional array.  On output, an upper
!    triangular matrix R, stored in packed form, so that A = R'*R.
!    If INFO /= 0, the factorization is not complete.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) RCOND, an estimate of the reciprocal condition
!    of A.  For the system A*X = B, relative perturbations in A and B of size
!    EPSILON may cause relative perturbations in X of size EPSILON/RCOND.
!    If RCOND is so small that the logical expression
!      1.0D+00 + RCOND == 1.0D+00
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate underflows.
!
!    Output, real ( kind = 8 ) Z(N), a work vector whose contents are usually
!    unimportant.  If A is singular to working precision, then Z is an
!    approximate null vector in the sense that
!      norm(A*Z) = RCOND * norm(A) * norm(Z).
!    If INFO /= 0, Z is unchanged.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, for normal return.
!    K, signals an error condition.  The leading minor of order K is
!    not positive definite.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) anorm
  real ( kind = 8 ) ap((n*(n+1))/2)
  real ( kind = 8 ) ek
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kj
  integer ( kind = 4 ) kk
  real ( kind = 8 ) rcond
  real ( kind = 8 ) s
  real ( kind = 8 ) dasum
  real ( kind = 8 ) ddot
  real ( kind = 8 ) sm
  real ( kind = 8 ) t
  real ( kind = 8 ) wk
  real ( kind = 8 ) wkm
  real ( kind = 8 ) ynorm
  real ( kind = 8 ) z(n)
!
!  Find the norm of A.
!
  j1 = 1
  do j = 1, n
    z(j) = dasum ( j, ap(j1), 1 )
    ij = j1
    j1 = j1 + j
    do i = 1, j-1
      z(i) = z(i) + abs ( ap(ij) )
      ij = ij + 1
    end do
  end do

  anorm = maxval ( z(1:n) )
!
!  Factor.
!
  call dppfa ( ap, n, info )

  if ( info /= 0 ) then
    return
  end if
!
!  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
!
!  Estimate = norm(Z)/norm(Y) where A * Z = Y and A * Y = E.
!
!  The components of E are chosen to cause maximum local
!  growth in the elements of W where R'*W = E.
!
!  The vectors are frequently rescaled to avoid overflow.
!
!  Solve R' * W = E.
!
  ek = 1.0D+00
  z(1:n) = 0.0D+00

  kk = 0

  do k = 1, n

    kk = kk + k

    if ( z(k) /= 0.0D+00 ) then
      ek = sign ( ek, -z(k) )
    end if

    if ( ap(kk) < abs ( ek - z(k) ) ) then
      s = ap(kk) / abs ( ek - z(k) )
      z(1:n) = s * z(1:n)
      ek = s * ek
    end if

    wk = ek - z(k)
    wkm = -ek - z(k)
    s = abs ( wk )
    sm = abs ( wkm )
    wk = wk / ap(kk)
    wkm = wkm / ap(kk)
    kj = kk + k

    if ( k + 1 <= n ) then

      do j = k + 1, n
        sm = sm + abs ( z(j) + wkm * ap(kj) )
        z(j) = z(j) + wk * ap(kj)
        s = s + abs ( z(j) )
        kj = kj + j
      end do

      if ( s < sm ) then

        t = wkm - wk
        wk = wkm
        kj = kk + k

        do j = k+1, n
          z(j) = z(j) + t * ap(kj)
          kj = kj + j
        end do

      end if

    end if

    z(k) = wk

  end do

  z(1:n) = z(1:n) / dasum ( n, z, 1 )
!
!  Solve R * Y = W.
!
  do k = n, 1, -1

    if ( ap(kk) < abs ( z(k) ) ) then
      s = ap(kk) / abs ( z(k) )
      z(1:n) = s * z(1:n)
    end if

    z(k) = z(k) / ap(kk)
    kk = kk - k
    t = -z(k)
    call daxpy ( k-1, t, ap(kk+1), 1, z(1), 1 )

  end do

  z(1:n) = z(1:n) / dasum ( n, z, 1 )

  ynorm = 1.0D+00
!
!  Solve R' * V = Y.
!
  do k = 1, n

    z(k) = z(k) - ddot ( k-1, ap(kk+1), 1, z(1), 1 )
    kk = kk + k

    if ( ap(kk) < abs ( z(k) ) ) then
      s = ap(kk) / abs ( z(k) )
      z(1:n) = s * z(1:n)
      ynorm = s * ynorm
    end if

    z(k) = z(k) / ap(kk)

  end do

  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm
!
!  Solve R * Z = V.
!
  do k = n, 1, -1

    if ( ap(kk) < abs ( z(k) ) ) then
      s = ap(kk) / abs ( z(k) )
      z(1:n) = s * z(1:n)
      ynorm = s * ynorm
    end if

    z(k) = z(k) / ap(kk)
    kk = kk - k
    t = -z(k)
    call daxpy ( k-1, t, ap(kk+1), 1, z(1), 1 )

  end do
!
!  Make ZNORM = 1.0.
!
  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm

  if ( anorm /= 0.0D+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0D+00
  end if

  return
end
subroutine dppdi ( ap, n, det, job )

!*****************************************************************************80
!
!! DPPDI: determinant and inverse of a matrix factored by DPPCO or DPPFA.
!
!  Discussion:
!
!    A division by zero will occur if the input factor contains
!    a zero on the diagonal and the inverse is requested.
!    It will not occur if the subroutines are called correctly
!    and if DPOCO or DPOFA has set INFO == 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) AP((N*(N+1))/2).  On input, the output
!    from DPPCO or DPPFA.  On output, the upper triangular half of the
!    inverse, if requested.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) DET(2), the determinant of the original matrix
!    if requested.
!      determinant = DET(1) * 10.0**DET(2)
!    with  1.0D+00 <= DET(1) < 10.0D+00 or DET(1) == 0.0D+00.
!
!    Input, integer ( kind = 4 ) JOB, job request.
!    11, both determinant and inverse.
!    01, inverse only.
!    10, determinant only.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) ap((n*(n+1))/2)
  real ( kind = 8 ) det(2)
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
  real ( kind = 8 ) t
!
!  Compute the determinant.
!
  if ( job / 10 /= 0 ) then

    det(1) = 1.0D+00
    det(2) = 0.0D+00
    ii = 0

    do i = 1, n

      ii = ii + i

      det(1) = det(1) * ap(ii) * ap(ii)

      if ( det(1) == 0.0D+00 ) then
        exit
      end if

      do while ( det(1) < 1.0D+00 )
        det(1) = det(1) * 10.0D+00
        det(2) = det(2) - 1.0D+00
      end do

      do while ( 10.0D+00 <= det(1) )
        det(1) = det(1) / 10.0D+00
        det(2) = det(2) + 1.0D+00
      end do

    end do

  end if
!
!  Compute inverse(R).
!
  if ( mod ( job, 10 ) /= 0 ) then

    kk = 0

    do k = 1, n

      k1 = kk + 1
      kk = kk + k
      ap(kk) = 1.0D+00 / ap(kk)
      t = -ap(kk)
      call dscal ( k-1, t, ap(k1), 1 )
      j1 = kk + 1
      kj = kk + k

      do j = k + 1, n
        t = ap(kj)
        ap(kj) = 0.0D+00
        call daxpy ( k, t, ap(k1), 1, ap(j1), 1 )
        j1 = j1 + j
        kj = kj + j
      end do

    end do
!
!  Form inverse(R) * (inverse(R))'.
!
    jj = 0

    do j = 1, n

      j1 = jj + 1
      jj = jj + j
      k1 = 1
      kj = j1

      do k = 1, j-1
        t = ap(kj)
        call daxpy ( k, t, ap(j1), 1, ap(k1), 1 )
        k1 = k1 + k
        kj = kj + 1
      end do

      t = ap(jj)
      call dscal ( j, t, ap(j1), 1 )

    end do

  end if

  return
end
subroutine dppfa ( ap, n, info )

!*****************************************************************************80
!
!! DPPFA factors a real symmetric positive definite matrix in packed form.
!
!  Discussion:
!
!    DPPFA is usually called by DPPCO, but it can be called
!    directly with a saving in time if RCOND is not needed.
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
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) AP((N*(N+1))/2).  On input, the packed
!    form of a symmetric matrix A.  The columns of the upper triangle are
!    stored sequentially in a one-dimensional array.  On output, an upper
!    triangular matrix R, stored in packed form, so that A = R'*R.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, for normal return.
!    K, if the leading minor of order K is not positive definite.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) ap((n*(n+1))/2)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kj
  integer ( kind = 4 ) kk
  real ( kind = 8 ) s
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t

  info = 0
  jj = 0

  do j = 1, n

    s = 0.0D+00
    kj = jj
    kk = 0

    do k = 1, j-1

      kj = kj + 1
      t = ap(kj) - ddot ( k-1, ap(kk+1), 1, ap(jj+1), 1 )
      kk = kk + k
      t = t / ap(kk)
      ap(kj) = t
      s = s + t * t

    end do

    jj = jj + j
    s = ap(jj) - s

    if ( s <= 0.0D+00 ) then
      info = j
      return
    end if

    ap(jj) = sqrt ( s )

  end do

  return
end
subroutine dppsl ( ap, n, b )

!*****************************************************************************80
!
!! DPPSL solves a real SPD system factored by DPPCO or DPPFA.
!
!  Discussion:
!
!    To compute inverse(A) * C where C is a matrix with P columns
!
!      call dppco ( ap, n, rcond, z, info )
!
!      if ( rcond is too small .or. info /= 0 ) then
!        exit
!      end if
!
!      do j = 1, p
!        call dppsl ( ap, n, c(1,j) )
!      end do
!
!    A division by zero will occur if the input factor contains
!    a zero on the diagonal.  Technically this indicates
!    singularity but it is usually caused by improper subroutine
!    arguments.  It will not occur if the subroutines are called
!    correctly and INFO == 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input, real ( kind = 8 ) AP((N*(N+1))/2), the output from DPPCO or DPPFA.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) ap((n*(n+1))/2)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t

  kk = 0

  do k = 1, n
    t = ddot ( k-1, ap(kk+1), 1, b(1), 1 )
    kk = kk + k
    b(k) = ( b(k) - t ) / ap(kk)
  end do

  do k = n, 1, -1
    b(k) = b(k) / ap(kk)
    kk = kk - k
    t = -b(k)
    call daxpy ( k-1, t, ap(kk+1), 1, b(1), 1 )
  end do

  return
end
subroutine dptsl ( n, d, e, b )

!*****************************************************************************80
!
!! DPTSL solves a positive definite tridiagonal linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N), on input the diagonal of the
!    tridiagonal matrix.  On output, D is destroyed.
!
!    Input, real ( kind = 8 ) E(N), the offdiagonal of the tridiagonal matrix in
!    entries E(1:N-1).
!
!    Input/output, real ( kind = 8 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b(n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kbm1
  integer ( kind = 4 ) ke
  integer ( kind = 4 ) kf
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) nm1d2
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
!
!  Check for 1 x 1 case.
!
  if ( n == 1 ) then
    b(1) = b(1) / d(1)
    return
  end if

  nm1d2 = ( n - 1 ) / 2

  if ( 2 < n ) then

    kbm1 = n - 1
!
!  Zero top half of subdiagonal and bottom half of superdiagonal.
!
    do k = 1, nm1d2
      t1 = e(k) / d(k)
      d(k+1) = d(k+1) - t1 * e(k)
      b(k+1) = b(k+1) - t1 * b(k)
      t2 = e(kbm1) / d(kbm1+1)
      d(kbm1) = d(kbm1) - t2 * e(kbm1)
      b(kbm1) = b(kbm1) - t2 * b(kbm1+1)
      kbm1 = kbm1 - 1
    end do

  end if

  kp1 = nm1d2 + 1
!
!  Clean up for possible 2 x 2 block at center.
!
  if ( mod ( n, 2 ) == 0 ) then
    t1 = e(kp1) / d(kp1)
    d(kp1+1) = d(kp1+1) - t1 * e(kp1)
    b(kp1+1) = b(kp1+1) - t1 * b(kp1)
    kp1 = kp1 + 1
  end if
!
!  Back solve starting at the center, going towards the top and bottom.
!
  b(kp1) = b(kp1) / d(kp1)

  if ( 2 < n ) then

    k = kp1 - 1
    ke = kp1 + nm1d2 - 1

    do kf = kp1, ke
      b(k) = ( b(k) - e(k) * b(k+1) ) / d(k)
      b(kf+1) = ( b(kf+1) - e(kf) * b(kf) ) / d(kf+1)
      k = k - 1
    end do

  end if

  if ( mod ( n, 2 ) == 0 ) then
    b(1) = ( b(1) - e(1) * b(2) ) / d(1)
  end if

  return
end
subroutine dqrdc ( a, lda, n, p, qraux, jpvt, work, job )

!*****************************************************************************80
!
!! DQRDC computes the QR factorization of a real rectangular matrix.
!
!  Discussion:
!
!    DQRDC uses Householder transformations.
!
!    Column pivoting based on the 2-norms of the reduced columns may be
!    performed at the user's option.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) A(LDA,P).  On input, the N by P matrix
!    whose decomposition is to be computed.  On output, A contains in
!    its upper triangle the upper triangular matrix R of the QR
!    factorization.  Below its diagonal A contains information from
!    which the orthogonal part of the decomposition can be recovered.
!    Note that if pivoting has been requested, the decomposition is not that
!    of the original matrix A but that of A with its columns permuted
!    as described by JPVT.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!    LDA must be at least N.
!
!    Input, integer ( kind = 4 ) N, the number of rows of the matrix A.
!
!    Input, integer ( kind = 4 ) P, the number of columns of the matrix A.
!
!    Output, real ( kind = 8 ) QRAUX(P), contains further information required
!    to recover the orthogonal part of the decomposition.
!
!    Input/output, integer ( kind = 4 ) JPVT(P).  On input, JPVT contains
!    integers that control the selection of the pivot columns.  The K-th
!    column A(*,K) of A is placed in one of three classes according to the
!    value of JPVT(K).
!      > 0, then A(K) is an initial column.
!      = 0, then A(K) is a free column.
!      < 0, then A(K) is a final column.
!    Before the decomposition is computed, initial columns are moved to
!    the beginning of the array A and final columns to the end.  Both
!    initial and final columns are frozen in place during the computation
!    and only free columns are moved.  At the K-th stage of the
!    reduction, if A(*,K) is occupied by a free column it is interchanged
!    with the free column of largest reduced norm.  JPVT is not referenced
!    if JOB == 0.  On output, JPVT(K) contains the index of the column of the
!    original matrix that has been interchanged into the K-th column, if
!    pivoting was requested.
!
!    Workspace, real ( kind = 8 ) WORK(P).  WORK is not referenced if JOB == 0.
!
!    Input, integer ( kind = 4 ) JOB, initiates column pivoting.
!    0, no pivoting is done.
!    nonzero, pivoting is done.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p

  real ( kind = 8 ) a(lda,p)
  integer ( kind = 4 ) jpvt(p)
  real ( kind = 8 ) qraux(p)
  real ( kind = 8 ) work(p)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lup
  integer ( kind = 4 ) maxj
  real ( kind = 8 ) maxnrm
  real ( kind = 8 ) nrmxl
  integer ( kind = 4 ) pl
  integer ( kind = 4 ) pu
  real ( kind = 8 ) ddot
  real ( kind = 8 ) dnrm2
  logical swapj
  real ( kind = 8 ) t
  real ( kind = 8 ) tt

  pl = 1
  pu = 0
!
!  If pivoting is requested, rearrange the columns.
!
  if ( job /= 0 ) then

    do j = 1, p

      swapj = 0 < jpvt(j)

      if ( jpvt(j) < 0 ) then
        jpvt(j) = -j
      else
        jpvt(j) = j
      end if

      if ( swapj ) then

        if ( j /= pl ) then
          call dswap ( n, a(1,pl), 1, a(1,j), 1 )
        end if

        jpvt(j) = jpvt(pl)
        jpvt(pl) = j
        pl = pl + 1

      end if

    end do

    pu = p

    do j = p, 1, -1

      if ( jpvt(j) < 0 ) then

        jpvt(j) = -jpvt(j)

        if ( j /= pu ) then
          call dswap ( n, a(1,pu), 1, a(1,j), 1 )
          jp = jpvt(pu)
          jpvt(pu) = jpvt(j)
          jpvt(j) = jp
        end if

        pu = pu - 1

      end if

    end do

  end if
!
!  Compute the norms of the free columns.
!
  do j = pl, pu
    qraux(j) = dnrm2 ( n, a(1,j), 1 )
  end do

  work(pl:pu) = qraux(pl:pu)
!
!  Perform the Householder reduction of A.
!
  lup = min ( n, p )

  do l = 1, lup
!
!  Bring the column of largest norm into the pivot position.
!
    if ( pl <= l .and. l < pu ) then

      maxnrm = 0.0D+00
      maxj = l
      do j = l, pu
        if ( maxnrm < qraux(j) ) then
          maxnrm = qraux(j)
          maxj = j
        end if
      end do

      if ( maxj /= l ) then
        call dswap ( n, a(1,l), 1, a(1,maxj), 1 )
        qraux(maxj) = qraux(l)
        work(maxj) = work(l)
        jp = jpvt(maxj)
        jpvt(maxj) = jpvt(l)
        jpvt(l) = jp
      end if

    end if
!
!  Compute the Householder transformation for column L.
!
    qraux(l) = 0.0D+00

    if ( l /= n ) then

      nrmxl = dnrm2 ( n-l+1, a(l,l), 1 )

      if ( nrmxl /= 0.0D+00 ) then

        if ( a(l,l) /= 0.0D+00 ) then
          nrmxl = sign ( nrmxl, a(l,l) )
        end if

        call dscal ( n-l+1, 1.0D+00 / nrmxl, a(l,l), 1 )
        a(l,l) = 1.0D+00 + a(l,l)
!
!  Apply the transformation to the remaining columns, updating the norms.
!
        do j = l + 1, p

          t = -ddot ( n-l+1, a(l,l), 1, a(l,j), 1 ) / a(l,l)
          call daxpy ( n-l+1, t, a(l,l), 1, a(l,j), 1 )

          if ( pl <= j .and. j <= pu ) then

            if ( qraux(j) /= 0.0D+00 ) then

              tt = 1.0D+00 - ( abs ( a(l,j) ) / qraux(j) )**2
              tt = max ( tt, 0.0D+00 )
              t = tt
              tt = 1.0D+00 + 0.05D+00 * tt * ( qraux(j) / work(j) )**2

              if ( tt /= 1.0D+00 ) then
                qraux(j) = qraux(j) * sqrt ( t )
              else
                qraux(j) = dnrm2 ( n-l, a(l+1,j), 1 )
                work(j) = qraux(j)
              end if

            end if

          end if

        end do
!
!  Save the transformation.
!
        qraux(l) = a(l,l)
        a(l,l) = -nrmxl

      end if

    end if

  end do

  return
end
subroutine dqrsl ( a, lda, n, k, qraux, y, qy, qty, b, rsd, ab, job, info )

!*****************************************************************************80
!
!! DQRSL computes transformations, projections, and least squares solutions.
!
!  Discussion:
!
!    DQRSL requires the output of DQRDC.
!
!    For K <= min(N,P), let AK be the matrix
!
!      AK = ( A(JPVT(1)), A(JPVT(2)), ..., A(JPVT(K)) )
!
!    formed from columns JPVT(1), ..., JPVT(K) of the original
!    N by P matrix A that was input to DQRDC.  If no pivoting was
!    done, AK consists of the first K columns of A in their
!    original order.  DQRDC produces a factored orthogonal matrix Q
!    and an upper triangular matrix R such that
!
!      AK = Q * (R)
!               (0)
!
!    This information is contained in coded form in the arrays
!    A and QRAUX.
!
!    The parameters QY, QTY, B, RSD, and AB are not referenced
!    if their computation is not requested and in this case
!    can be replaced by dummy variables in the calling program.
!    To save storage, the user may in some cases use the same
!    array for different parameters in the calling sequence.  A
!    frequently occuring example is when one wishes to compute
!    any of B, RSD, or AB and does not need Y or QTY.  In this
!    case one may identify Y, QTY, and one of B, RSD, or AB, while
!    providing separate arrays for anything else that is to be
!    computed.
!
!    Thus the calling sequence
!
!      call dqrsl ( a, lda, n, k, qraux, y, dum, y, b, y, dum, 110, info )
!
!    will result in the computation of B and RSD, with RSD
!    overwriting Y.  More generally, each item in the following
!    list contains groups of permissible identifications for
!    a single calling sequence.
!
!      1. (Y,QTY,B) (RSD) (AB) (QY)
!
!      2. (Y,QTY,RSD) (B) (AB) (QY)
!
!      3. (Y,QTY,AB) (B) (RSD) (QY)
!
!      4. (Y,QY) (QTY,B) (RSD) (AB)
!
!      5. (Y,QY) (QTY,RSD) (B) (AB)
!
!      6. (Y,QY) (QTY,AB) (B) (RSD)
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
!    17 May 2005
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
!    Input, real ( kind = 8 ) A(LDA,P), contains the output of DQRDC.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) N, the number of rows of the matrix AK.  It
!    must have the same value as N in DQRDC.
!
!    Input, integer ( kind = 4 ) K, the number of columns of the matrix AK.  K
!    must not be greater than min(N,P), where P is the same as in the
!    calling sequence to DQRDC.
!
!    Input, real ( kind = 8 ) QRAUX(P), the auxiliary output from DQRDC.
!
!    Input, real ( kind = 8 ) Y(N), a vector to be manipulated by DQRSL.
!
!    Output, real ( kind = 8 ) QY(N), contains Q * Y, if requested.
!
!    Output, real ( kind = 8 ) QTY(N), contains Q' * Y, if requested.
!
!    Output, real ( kind = 8 ) B(K), the solution of the least squares problem
!      minimize norm2 ( Y - AK * B),
!    if its computation has been requested.  Note that if pivoting was
!    requested in DQRDC, the J-th component of B will be associated with
!    column JPVT(J) of the original matrix A that was input into DQRDC.
!
!    Output, real ( kind = 8 ) RSD(N), the least squares residual Y - AK * B,
!    if its computation has been requested.  RSD is also the orthogonal
!    projection of Y onto the orthogonal complement of the column space
!    of AK.
!
!    Output, real ( kind = 8 ) AB(N), the least squares approximation Ak * B,
!    if its computation has been requested.  AB is also the orthogonal
!    projection of Y onto the column space of A.
!
!    Input, integer ( kind = 4 ) JOB, specifies what is to be computed.  JOB has
!    the decimal expansion ABCDE, with the following meaning:
!      if A /= 0, compute QY.
!      if B /= 0, compute QTY.
!      if C /= 0, compute QTY and B.
!      if D /= 0, compute QTY and RSD.
!      if E /= 0, compute QTY and AB.
!    Note that a request to compute B, RSD, or AB automatically triggers
!    the computation of QTY, for which an array must be provided in the
!    calling sequence.
!
!    Output, integer ( kind = 4 ) INFO, is zero unless the computation of B has
!    been requested and R is exactly singular.  In this case, INFO is the
!    index of the first zero diagonal element of R, and B is left unaltered.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,*)
  real ( kind = 8 ) ab(n)
  real ( kind = 8 ) b(k)
  logical cab
  logical cb
  logical cqty
  logical cqy
  logical cr
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) job
  integer ( kind = 4 ) ju
  real ( kind = 8 ) qraux(*)
  real ( kind = 8 ) qty(n)
  real ( kind = 8 ) qy(n)
  real ( kind = 8 ) rsd(n)
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t
  real ( kind = 8 ) temp
  real ( kind = 8 ) y(n)
!
!  Set info flag.
!
  info = 0
!
!  Determine what is to be computed.
!
  cqy =        job / 10000         /= 0
  cqty = mod ( job,  10000 )       /= 0
  cb =   mod ( job,   1000 ) / 100 /= 0
  cr =   mod ( job,    100 ) /  10 /= 0
  cab =  mod ( job,     10 )       /= 0

  ju = min ( k, n - 1 )
!
!  Special action when N = 1.
!
  if ( ju == 0 ) then

    if ( cqy ) then
      qy(1) = y(1)
    end if

    if ( cqty ) then
      qty(1) = y(1)
    end if

    if ( cab ) then
      ab(1) = y(1)
    end if

    if ( cb ) then

      if ( a(1,1) == 0.0D+00 ) then
        info = 1
      else
        b(1) = y(1) / a(1,1)
      end if

    end if

    if ( cr ) then
      rsd(1) = 0.0D+00
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

      if ( qraux(j) /= 0.0D+00 ) then
        temp = a(j,j)
        a(j,j) = qraux(j)
        t = -ddot ( n-j+1, a(j,j), 1, qy(j), 1 ) / a(j,j)
        call daxpy ( n-j+1, t, a(j,j), 1, qy(j), 1 )
        a(j,j) = temp
      end if

    end do

  end if
!
!  Compute Q'*Y.
!
  if ( cqty ) then

    do j = 1, ju
      if ( qraux(j) /= 0.0D+00 ) then
        temp = a(j,j)
        a(j,j) = qraux(j)
        t = -ddot ( n-j+1, a(j,j), 1, qty(j), 1 ) / a(j,j)
        call daxpy ( n-j+1, t, a(j,j), 1, qty(j), 1 )
        a(j,j) = temp
      end if
    end do

  end if
!
!  Set up to compute B, RSD, or AB.
!
  if ( cb ) then
    b(1:k) = qty(1:k)
  end if

  if ( cab ) then
    ab(1:k) = qty(1:k)
  end if

  if ( cr .and. k < n ) then
    rsd(k+1:n) = qty(k+1:n)
  end if

  if ( cab .and. k+1 <= n ) then
    ab(k+1:n) = 0.0D+00
  end if

  if ( cr ) then
    rsd(1:k) = 0.0D+00
  end if
!
!  Compute B.
!
  if ( cb ) then

    do jj = 1, k

      j = k - jj + 1

      if ( a(j,j) == 0.0D+00 ) then
        info = j
        exit
      end if

      b(j) = b(j) / a(j,j)

      if ( j /= 1 ) then
        t = -b(j)
        call daxpy ( j-1, t, a(1,j), 1, b, 1 )
      end if

    end do

  end if
!
!  Compute RSD or AB as required.
!
  if ( cr .or. cab ) then

    do jj = 1, ju

      j = ju - jj + 1

      if ( qraux(j) /= 0.0D+00 ) then

        temp = a(j,j)
        a(j,j) = qraux(j)

        if ( cr ) then
          t = -ddot ( n-j+1, a(j,j), 1, rsd(j), 1 ) / a(j,j)
          call daxpy ( n-j+1, t, a(j,j), 1, rsd(j), 1 )
        end if

        if ( cab ) then
          t = -ddot ( n-j+1, a(j,j), 1, ab(j), 1 ) / a(j,j)
          call daxpy ( n-j+1, t, a(j,j), 1, ab(j), 1 )
        end if

        a(j,j) = temp

      end if

    end do

  end if

  return
end
subroutine dsico ( a, lda, n, kpvt, rcond, z )

!*****************************************************************************80
!
!! DSICO factors a real symmetric matrix and estimates its condition.
!
!  Discussion:
!
!    If RCOND is not needed, DSIFA is slightly faster.
!
!    To solve A * X = B, follow DSICO by DSISL.
!
!    To compute inverse(A)*C, follow DSICO by DSISL.
!
!    To compute inverse(A), follow DSICO by DSIDI.
!
!    To compute determinant(A), follow DSICO by DSIDI.
!
!    To compute inertia(A), follow DSICO by DSIDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) A(LDA,N).  On input, the symmetric
!    matrix to be factored.  Only the diagonal and upper triangle are used.
!    On output, a block diagonal matrix and the multipliers which
!    were used to obtain it.  The factorization can be written A = U*D*U'
!    where U is a product of permutation and unit upper triangular
!    matrices, U' is the transpose of U, and D is block diagonal
!    with 1 by 1 and 2 by 2 blocks.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) KPVT(N), pivot indices.
!
!    Output, real ( kind = 8 ) RCOND, an estimate of the reciprocal condition
!    of A.  For the system A*X = B, relative perturbations in A and B of size
!    EPSILON may cause relative perturbations in X of size EPSILON/RCOND.
!    If RCOND is so small that the logical expression
!      1.0D+00 + RCOND == 1.0D+00
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate underflows.
!
!    Output, real ( kind = 8 ) Z(N), a work vector whose contents are usually
!    unimportant.  If A is close to a singular matrix, then Z is an
!    approximate null vector in the sense that
!      norm(A*Z) = RCOND * norm(A) * norm(Z).
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) ak
  real ( kind = 8 ) akm1
  real ( kind = 8 ) anorm
  real ( kind = 8 ) bk
  real ( kind = 8 ) bkm1
  real ( kind = 8 ) denom
  real ( kind = 8 ) ek
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp
  integer ( kind = 4 ) kps
  integer ( kind = 4 ) kpvt(n)
  integer ( kind = 4 ) ks
  real ( kind = 8 ) rcond
  real ( kind = 8 ) s
  real ( kind = 8 ) dasum
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t
  real ( kind = 8 ) ynorm
  real ( kind = 8 ) z(n)
!
!  Find the norm of A, using only entries in the upper half of the matrix.
!
  do j = 1, n
    z(j) = dasum ( j, a(1,j), 1 )
    do i = 1, j-1
      z(i) = z(i) + abs ( a(i,j) )
    end do
  end do

  anorm = maxval ( z(1:n) )
!
!  Factor.
!
  call dsifa ( a, lda, n, kpvt, info )
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
!  Solve U * D * W = E.
!
  ek = 1.0D+00
  z(1:n) = 0.0D+00
  k = n

  do while ( k /= 0 )

    if ( kpvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    kp = abs ( kpvt(k) )
    kps = k + 1 - ks

    if ( kp /= kps ) then
      t = z(kps)
      z(kps) = z(kp)
      z(kp) = t
    end if

    if ( z(k) /= 0.0D+00 ) then
      ek = sign ( ek, z(k) )
    end if

    z(k) = z(k) + ek
    call daxpy ( k-ks, z(k), a(1,k), 1, z(1), 1 )

    if ( ks /= 1 ) then
      if ( z(k-1) /= 0.0D+00 ) then
        ek = sign ( ek, z(k-1) )
      end if
      z(k-1) = z(k-1) + ek
      call daxpy ( k-ks, z(k-1), a(1,k-1), 1, z(1), 1 )
    end if

    if ( ks /= 2 ) then

      if ( abs ( a(k,k) ) < abs ( z(k) ) ) then
        s = abs ( a(k,k) ) / abs ( z(k) )
        z(1:n) = s * z(1:n)
        ek = s * ek
      end if

      if ( a(k,k) /= 0.0D+00 ) then
        z(k) = z(k) / a(k,k)
      else
        z(k) = 1.0D+00
      end if

    else

      ak = a(k,k) / a(k-1,k)
      akm1 = a(k-1,k-1) / a(k-1,k)
      bk = z(k) / a(k-1,k)
      bkm1 = z(k-1) / a(k-1,k)
      denom = ak * akm1 - 1.0D+00
      z(k) = ( akm1 * bk - bkm1 ) / denom
      z(k-1) = ( ak * bkm1 - bk ) / denom

    end if

    k = k - ks

  end do

  z(1:n) = z(1:n) / dasum ( n, z, 1 )
!
!  Solve U' * Y = W.
!
  k = 1

  do while ( k <= n )

    if ( kpvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= 1 ) then

      z(k) = z(k) + ddot ( k-1, a(1,k), 1, z(1), 1 )

      if ( ks == 2 ) then
        z(k+1) = z(k+1) + ddot ( k-1, a(1,k+1), 1, z(1), 1 )
      end if

      kp = abs ( kpvt(k) )

      if ( kp /= k ) then
        t = z(k)
        z(k) = z(kp)
        z(kp) = t
      end if

    end if

    k = k + ks

  end do

  z(1:n) = z(1:n) / dasum ( n, z, 1 )

  ynorm = 1.0D+00
!
!  Solve U * D * V = Y.
!
  k = n

  do while ( k /= 0 )

    if ( kpvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= ks ) then

      kp = abs ( kpvt(k) )
      kps = k + 1 - ks

      if ( kp /= kps ) then
        t = z(kps)
        z(kps) = z(kp)
        z(kp) = t
      end if

      call daxpy ( k-ks, z(k), a(1,k), 1, z(1), 1 )

      if ( ks == 2 ) then
        call daxpy ( k-ks, z(k-1), a(1,k-1), 1, z(1), 1 )
      end if

    end if

    if ( ks /= 2 ) then

      if ( abs ( a(k,k) ) < abs ( z(k) ) ) then
        s = abs ( a(k,k) ) / abs ( z(k) )
        z(1:n) = s * z(1:n)
        ynorm = s * ynorm
      end if

      if ( a(k,k) /= 0.0D+00 ) then
        z(k) = z(k) / a(k,k)
      else
        z(k) = 1.0D+00
      end if

    else

      ak = a(k,k) / a(k-1,k)
      akm1 = a(k-1,k-1) / a(k-1,k)
      bk = z(k) / a(k-1,k)
      bkm1 = z(k-1) / a(k-1,k)
      denom = ak * akm1 - 1.0D+00
      z(k) = ( akm1 * bk - bkm1 ) / denom
      z(k-1) = ( ak * bkm1 - bk ) / denom

    end if

    k = k - ks

  end do

  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm
!
!  Solve U' * Z = V.
!
  k = 1

  do while ( k <= n )

    if ( kpvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= 1 ) then

      z(k) = z(k) + ddot ( k-1, a(1,k), 1, z(1), 1 )
      if ( ks == 2 ) then
        z(k+1) = z(k+1) + ddot ( k-1, a(1,k+1), 1, z(1), 1 )
      end if
      kp = abs ( kpvt(k) )

      if ( kp /= k ) then
        t = z(k)
        z(k) = z(kp)
        z(kp) = t
      end if

    end if

    k = k + ks

  end do
!
!  Make ZNORM = 1.0.
!
  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm

  if ( anorm /= 0.0D+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0D+00
  end if

  return
end
subroutine dsidi ( a, lda, n, kpvt, det, inert, work, job )

!*****************************************************************************80
!
!! DSIDI: determinant, inertia and inverse of a real symmetric matrix.
!
!  Discussion:
!
!    DSIDI uses the factors from DSIFA.
!
!    A division by zero may occur if the inverse is requested
!    and DSICO has set RCOND == 0.0D+00 or DSIFA has set INFO /= 0.
!
!    Variables not requested by JOB are not used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) A(LDA,N).  On input, the output from DSIFA.
!    On output, the upper triangle of the inverse of the original matrix.
!    The strict lower triangle is never referenced.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) KPVT(N), the pivot vector from DSIFA.
!
!    Output, real ( kind = 8 ) DET(2), the determinant of the original matrix,
!    if requested.
!      determinant = DET(1) * 10.0**DET(2)
!    with 1.0D+00 <= abs ( DET(1) ) < 10.0D+00 or DET(1) = 0.0.
!
!    Output, integer ( kind = 4 ) INERT(3), the inertia of the original matrix.
!    INERT(1) = number of positive eigenvalues.
!    INERT(2) = number of negative eigenvalues.
!    INERT(3) = number of zero eigenvalues.
!
!    Workspace, real ( kind = 8 ) WORK(N).
!
!    Input, integer ( kind = 4 ) JOB, specifies the tasks.
!    JOB has the decimal expansion ABC where
!    If C /= 0, the inverse is computed,
!    If B /= 0, the determinant is computed,
!    If A /= 0, the inertia is computed.
!    For example, JOB = 111 gives all three.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) ak
  real ( kind = 8 ) akkp1
  real ( kind = 8 ) akp1
  real ( kind = 8 ) d
  real ( kind = 8 ) det(2)
  logical dodet
  logical doert
  logical doinv
  integer ( kind = 4 ) inert(3)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jb
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kpvt(n)
  integer ( kind = 4 ) ks
  integer ( kind = 4 ) kstep
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t
  real ( kind = 8 ) temp
  real ( kind = 8 ) work(n)

  doinv = mod ( job,   10 )       /= 0
  dodet = mod ( job,  100 ) /  10 /= 0
  doert = mod ( job, 1000 ) / 100 /= 0

  if ( dodet .or. doert ) then

    if ( doert ) then
      inert(1:3) = 0
    end if

    if ( dodet ) then
      det(1) = 1.0D+00
      det(2) = 0.0D+00
    end if

    t = 0.0D+00

    do k = 1, n

      d = a(k,k)
!
!  2 by 2 block.
!
!  use det (d  s)  =  (d/t * c - t) * t,  t = abs ( s )
!          (s  c)
!  to avoid underflow/overflow troubles.
!
!  Take two passes through scaling.  Use T for flag.
!
      if ( kpvt(k) <= 0 ) then

        if ( t == 0.0D+00 ) then
          t = abs ( a(k,k+1) )
          d = ( d / t ) * a(k+1,k+1) - t
        else
          d = t
          t = 0.0D+00
        end if

      end if

      if ( doert ) then
        if ( 0.0D+00 < d ) then
          inert(1) = inert(1) + 1
        else if ( d < 0.0D+00 ) then
          inert(2) = inert(2) + 1
        else if ( d == 0.0D+00 ) then
          inert(3) = inert(3) + 1
        end if
      end if

      if ( dodet ) then

        det(1) = det(1) * d

        if ( det(1) /= 0.0D+00 ) then

          do while ( abs ( det(1) ) < 1.0D+00 )
            det(1) = det(1) * 10.0D+00
            det(2) = det(2) - 1.0D+00
          end do

          do while ( 10.0D+00 <= abs ( det(1) ) )
            det(1) = det(1) / 10.0D+00
            det(2) = det(2) + 1.0D+00
          end do

        end if

      end if

    end do

  end if
!
!  Compute inverse(A).
!
  if ( doinv ) then

    k = 1

    do while ( k <= n )

      if ( 0 <= kpvt(k) ) then
!
!  1 by 1.
!
        a(k,k) = 1.0D+00 / a(k,k)

        if ( 2 <= k ) then

          work(1:k-1) = a(1:k-1,k)

          do j = 1, k-1
            a(j,k) = ddot ( j, a(1,j), 1, work, 1 )
            call daxpy ( j-1, work(j), a(1,j), 1, a(1,k), 1 )
          end do

          a(k,k) = a(k,k) + ddot ( k-1, work, 1, a(1,k), 1 )

        end if

        kstep = 1
!
!  2 by 2.
!
      else

        t = abs ( a(k,k+1) )
        ak = a(k,k) / t
        akp1 = a(k+1,k+1) / t
        akkp1 = a(k,k+1) / t
        d = t * ( ak * akp1 - 1.0D+00 )
        a(k,k) = akp1 / d
        a(k+1,k+1) = ak / d
        a(k,k+1) = -akkp1 / d

        if ( 2 <= k ) then

          work(1:k-1) = a(1:k-1,k+1)

          do j = 1, k-1
            a(j,k+1) = ddot ( j, a(1,j), 1, work, 1 )
            call daxpy ( j-1, work(j), a(1,j), 1, a(1,k+1), 1 )
          end do

          a(k+1,k+1) = a(k+1,k+1) + ddot ( k-1, work, 1, a(1,k+1), 1 )
          a(k,k+1) = a(k,k+1) + ddot ( k-1, a(1,k), 1, a(1,k+1), 1 )

          work(1:k-1) = a(1:k-1,k)

          do j = 1, k-1
            a(j,k) = ddot ( j, a(1,j), 1, work, 1 )
            call daxpy ( j-1, work(j), a(1,j), 1, a(1,k), 1 )
          end do

          a(k,k) = a(k,k) + ddot ( k-1, work, 1, a(1,k), 1 )

        end if

        kstep = 2

      end if
!
!  Swap.
!
      ks = abs ( kpvt(k) )

      if ( ks /= k ) then

        call dswap ( ks, a(1,ks), 1, a(1,k), 1 )

        do jb = ks, k
          j = k + ks - jb
          temp = a(j,k)
          a(j,k) = a(ks,j)
          a(ks,j) = temp
        end do

        if ( kstep /= 1 ) then
          temp = a(ks,k+1)
          a(ks,k+1) = a(k,k+1)
          a(k,k+1) = temp
        end if

      end if

      k = k + kstep

    end do

  end if

  return
end
subroutine dsifa ( a, lda, n, kpvt, info )

!*****************************************************************************80
!
!! DSIFA factors a real symmetric matrix.
!
!  Discussion:
!
!    To solve A*X = B, follow DSIFA by DSISL.
!
!    To compute inverse(A)*C, follow DSIFA by DSISL.
!
!    To compute determinant(A), follow DSIFA by DSIDI.
!
!    To compute inertia(A), follow DSIFA by DSIDI.
!
!    To compute inverse(A), follow DSIFA by DSIDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) A(LDA,N).  On input, the symmetric matrix
!    to be factored.  Only the diagonal and upper triangle are used.
!    On output, a block diagonal matrix and the multipliers which
!    were used to obtain it.  The factorization can be written A = U*D*U'
!    where U is a product of permutation and unit upper triangular
!    matrices, U' is the transpose of U, and D is block diagonal
!    with 1 by 1 and 2 by 2 blocks.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) KPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, normal value.
!    K, if the K-th pivot block is singular.  This is not an error
!    condition for this subroutine, but it does indicate that DSISL
!    or DSIDI may divide by zero if called.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) absakk
  real ( kind = 8 ) ak
  real ( kind = 8 ) akm1
  real ( kind = 8 ) alpha
  real ( kind = 8 ) bk
  real ( kind = 8 ) bkm1
  real ( kind = 8 ) colmax
  real ( kind = 8 ) denom
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imaxp1
  integer ( kind = 4 ) info
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kpvt(n)
  integer ( kind = 4 ) kstep
  real ( kind = 8 ) mulk
  real ( kind = 8 ) mulkm1
  real ( kind = 8 ) rowmax
  logical swap
  real ( kind = 8 ) t
!
!  ALPHA is used in choosing pivot block size.
!
  alpha = ( 1.0D+00 + sqrt ( 17.0D+00 ) ) / 8.0D+00

  info = 0
!
!  Main loop on K, which goes from N to 1.
!
  k = n

  do while ( 0 < k )

    if ( k == 1 ) then
      kpvt(1) = 1
      if ( a(1,1) == 0.0D+00 ) then
        info = 1
      end if
      return
    end if
!
!  This section of code determines the kind of
!  elimination to be performed.  When it is completed,
!  KSTEP will be set to the size of the pivot block, and
!  SWAP will be set to TRUE if an interchange is required.
!
    absakk = abs ( a(k,k) )
!
!  Determine the largest off-diagonal element in column K.
!
    imax = idamax ( k-1, a(1,k), 1 )
    colmax = abs ( a(imax,k) )

    if ( alpha * colmax <= absakk ) then

      kstep = 1
      swap = .false.
!
!  Determine the largest off-diagonal element in row IMAX.
!
    else

      rowmax = 0.0D+00
      imaxp1 = imax + 1
      do j = imax+1, k
        rowmax = max ( rowmax, abs ( a(imax,j) ) )
      end do

      if ( imax /= 1 ) then
        jmax = idamax ( imax-1, a(1,imax), 1 )
        rowmax = max ( rowmax, abs ( a(jmax,imax) ) )
      end if

      if ( alpha * rowmax <= abs ( a(imax,imax) ) ) then
        kstep = 1
        swap = .true.
      else if ( alpha * colmax * ( colmax / rowmax ) <= absakk ) then
        kstep = 1
        swap = .false.
      else
        kstep = 2
        swap = ( imax /= k-1 )
      end if

    end if
!
!  Column K is zero.
!  Set INFO and iterate the loop.
!
    if ( max ( absakk, colmax ) == 0.0D+00 ) then

      kpvt(k) = k
      info = k
!
!  1 x 1 pivot block.
!
!  Perform an interchange.
!
    else if ( kstep /= 2 ) then

      if ( swap ) then

        call dswap ( imax, a(1,imax), 1, a(1,k), 1 )

        do jj = imax, k
          j = k + imax - jj
          t = a(j,k)
          a(j,k) = a(imax,j)
          a(imax,j) = t
        end do

      end if
!
!  Perform the elimination.
!
      do jj = 1, k-1
        j = k - jj
        mulk = -a(j,k) / a(k,k)
        t = mulk
        call daxpy ( j, t, a(1,k), 1, a(1,j), 1 )
        a(j,k) = mulk
      end do
!
!  Set the pivot array.
!
      if ( swap ) then
        kpvt(k) = imax
      else
        kpvt(k) = k
      end if
!
!  2 x 2 pivot block.
!
!  Perform an interchange.
!
    else

      if ( swap ) then

        call dswap ( imax, a(1,imax), 1, a(1,k-1), 1 )

        do jj = imax, k-1
          j = k-1 + imax - jj
          t = a(j,k-1)
          a(j,k-1) = a(imax,j)
          a(imax,j) = t
        end do

        t = a(k-1,k)
        a(k-1,k) = a(imax,k)
        a(imax,k) = t

      end if
!
!  Perform the elimination.
!
      if ( k-2 /= 0 ) then

        ak = a(k,k) / a(k-1,k)
        akm1 = a(k-1,k-1) / a(k-1,k)
        denom = 1.0D+00 - ak * akm1

        do jj = 1, k-2

          j = k-1 - jj
          bk = a(j,k) / a(k-1,k)
          bkm1 = a(j,k-1) / a(k-1,k)
          mulk = ( akm1 * bk - bkm1 ) / denom
          mulkm1 = ( ak * bkm1 - bk ) / denom
          t = mulk
          call daxpy ( j, t, a(1,k), 1, a(1,j), 1 )
          t = mulkm1
          call daxpy ( j, t, a(1,k-1), 1, a(1,j), 1 )
          a(j,k) = mulk
          a(j,k-1) = mulkm1

        end do

      end if
!
!  Set the pivot array.
!
      if ( swap ) then
        kpvt(k) = -imax
      else
        kpvt(k) = 1 - k
      end if

      kpvt(k-1) = kpvt(k)

    end if

    k = k - kstep

  end do

  return
end
subroutine dsisl ( a, lda, n, kpvt, b )

!*****************************************************************************80
!
!! DSISL solves a real symmetric system factored by DSIFA.
!
!  Discussion:
!
!    To compute inverse(A) * C where C is a matrix with P columns
!
!      call dsifa ( a, lda, n, kpvt, info )
!
!      if ( info == 0 ) then
!        do j = 1, p
!          call dsisl ( a, lda, n, kpvt, c(1,j) )
!        end do
!      end if
!
!    A division by zero may occur if the inverse is requested
!    and DSICO has set RCOND == 0.0D+00 or DSIFA has set INFO /= 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input, real ( kind = 8 ) A(LDA,N), the output from DSIFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) KPVT(N), the pivot vector from DSIFA.
!
!    Input/output, real ( kind = 8 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) ak
  real ( kind = 8 ) akm1
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) bk
  real ( kind = 8 ) bkm1
  real ( kind = 8 ) denom
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp
  integer ( kind = 4 ) kpvt(n)
  real ( kind = 8 ) ddot
  real ( kind = 8 ) temp
!
!  Loop backward applying the transformations and D inverse to B.
!
  k = n

  do while ( 0 < k )

    if ( 0 <= kpvt(k) ) then
!
!  1 x 1 pivot block.
!
      if ( k /= 1 ) then

        kp = kpvt(k)
!
!  Interchange.
!
        if ( kp /= k ) then
          temp = b(k)
          b(k) = b(kp)
          b(kp) = temp
        end if
!
!  Apply the transformation.
!
        call daxpy ( k-1, b(k), a(1,k), 1, b(1), 1 )

      end if
!
!  Apply D inverse.
!
      b(k) = b(k) / a(k,k)
      k = k - 1

    else
!
!  2 x 2 pivot block.
!
      if ( k /= 2 ) then

        kp = abs ( kpvt(k) )
!
!  Interchange.
!
        if ( kp /= k-1 ) then
          temp = b(k-1)
          b(k-1) = b(kp)
          b(kp) = temp
        end if
!
!  Apply the transformation.
!
        call daxpy ( k-2, b(k), a(1,k), 1, b(1), 1 )
        call daxpy ( k-2, b(k-1), a(1,k-1), 1, b(1), 1 )

      end if
!
!  Apply D inverse.
!
      ak = a(k,k) / a(k-1,k)
      akm1 = a(k-1,k-1) / a(k-1,k)
      bk = b(k) / a(k-1,k)
      bkm1 = b(k-1) / a(k-1,k)
      denom = ak * akm1 - 1.0D+00
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

    if ( 0 <= kpvt(k) ) then
!
!  1 x 1 pivot block.
!
      if ( k /= 1 ) then
!
!  Apply the transformation.
!
        b(k) = b(k) + ddot ( k-1, a(1,k), 1, b(1), 1 )
        kp = kpvt(k)
!
!  Interchange.
!
        if ( kp /= k ) then
          temp = b(k)
          b(k) = b(kp)
          b(kp) = temp
        end if

      end if

      k = k + 1

    else
!
!  2 x 2 pivot block.
!
      if ( k /= 1 ) then
!
!  Apply the transformation.
!
        b(k) = b(k) + ddot ( k-1, a(1,k), 1, b(1), 1 )
        b(k+1) = b(k+1) + ddot ( k-1, a(1,k+1), 1, b(1), 1 )
        kp = abs ( kpvt(k) )
!
!  Interchange.
!
        if ( kp /= k ) then
          temp = b(k)
          b(k) = b(kp)
          b(kp) = temp
        end if

      end if

      k = k + 2

    end if

  end do

  return
end
subroutine dspco ( ap, n, kpvt, rcond, z )

!*****************************************************************************80
!
!! DSPCO factors a real symmetric matrix stored in packed form.
!
!  Discussion:
!
!    DSPCO uses elimination with symmetric pivoting and estimates
!    the condition of the matrix.
!
!    If RCOND is not needed, DSPFA is slightly faster.
!
!    To solve A*X = B, follow DSPCO by DSPSL.
!
!    To compute inverse(A)*C, follow DSPCO by DSPSL.
!
!    To compute inverse(A), follow DSPCO by DSPDI.
!
!    To compute determinant(A), follow DSPCO by DSPDI.
!
!    To compute inertia(A), follow DSPCO by DSPDI.
!
!  Packed storage:
!
!    The following program segment will pack the upper triangle of a
!    symmetric matrix.
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
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) AP((N*(N+1))/2).  On input, the packed form
!    of a symmetric matrix A.  The columns of the upper triangle are stored
!    sequentially in a one-dimensional array.  On output, a block diagonal
!    matrix and the multipliers which were used to obtain it, stored in
!    packed form.  The factorization can be written A = U*D*U'
!    where U is a product of permutation and unit upper triangular
!    matrices, U' is the transpose of U, and D is block diagonal
!    with 1 by 1 and 2 by 2 blocks.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) KPVT(N), the pivot indices.
!
!    Output, real ( kind = 8 ) RCOND, an estimate of the reciprocal condition
!    of A.  For the system A*X = B, relative perturbations in A and B of size
!    EPSILON may cause relative perturbations in X of size EPSILON/RCOND.
!    If RCOND is so small that the logical expression
!      1.0D+00 + RCOND == 1.0D+00
!    is true, then A may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate underflows.
!
!    Output, real ( kind = 8 ) Z(N) a work vector whose contents are usually
!    unimportant.  If A is close to a singular matrix, then Z is an
!    approximate null vector in the sense that
!      norm(A*Z) = RCOND * norm(A) * norm(Z).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) ak
  real ( kind = 8 ) akm1
  real ( kind = 8 ) anorm
  real ( kind = 8 ) ap((n*(n+1))/2)
  real ( kind = 8 ) bk
  real ( kind = 8 ) bkm1
  real ( kind = 8 ) denom
  real ( kind = 8 ) ek
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ikm1
  integer ( kind = 4 ) ikp1
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) km1k
  integer ( kind = 4 ) km1km1
  integer ( kind = 4 ) kp
  integer ( kind = 4 ) kps
  integer ( kind = 4 ) kpvt(n)
  integer ( kind = 4 ) ks
  real ( kind = 8 ) rcond
  real ( kind = 8 ) s
  real ( kind = 8 ) dasum
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t
  real ( kind = 8 ) ynorm
  real ( kind = 8 ) z(n)
!
!  Find norm of A using only upper half.
!
  j1 = 1
  do j = 1, n
    z(j) = dasum ( j, ap(j1), 1 )
    ij = j1
    j1 = j1 + j
    do i = 1, j-1
      z(i) = z(i) + abs ( ap(ij) )
      ij = ij + 1
    end do
  end do

  anorm = maxval ( z(1:n) )
!
!  Factor.
!
  call dspfa ( ap, n, kpvt, info )
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
!  Solve U * D * W = E.
!
  ek = 1.0D+00
  z(1:n) = 0.0D+00

  k = n
  ik = ( n * ( n - 1 ) ) / 2

  do while ( k /= 0 )

    kk = ik + k
    ikm1 = ik - ( k - 1 )

    if ( kpvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    kp = abs ( kpvt(k) )
    kps = k + 1 - ks

    if ( kp /= kps ) then
      t = z(kps)
      z(kps) = z(kp)
      z(kp) = t
    end if

    if ( z(k) /= 0.0D+00 ) then
      ek = sign ( ek, z(k) )
    end if

    z(k) = z(k) + ek
    call daxpy ( k-ks, z(k), ap(ik+1), 1, z(1), 1 )

    if ( ks /= 1 ) then
      if ( z(k-1) /= 0.0D+00 ) then
        ek = sign ( ek, z(k-1) )
      end if
      z(k-1) = z(k-1) + ek
      call daxpy ( k-ks, z(k-1), ap(ikm1+1), 1, z(1), 1 )
    end if

    if ( ks /= 2 ) then

      if ( abs ( ap(kk) ) < abs ( z(k) ) ) then
        s = abs ( ap(kk) ) / abs ( z(k) )
        z(1:n) = s * z(1:n)
        ek = s * ek
      end if

      if ( ap(kk) /= 0.0D+00 ) then
        z(k) = z(k) / ap(kk)
      else
        z(k) = 1.0D+00
      end if

    else

      km1k = ik + k - 1
      km1km1 = ikm1 + k - 1
      ak = ap(kk) / ap(km1k)
      akm1 = ap(km1km1) / ap(km1k)
      bk = z(k) / ap(km1k)
      bkm1 = z(k-1) / ap(km1k)
      denom = ak * akm1 - 1.0D+00
      z(k) = ( akm1 * bk - bkm1 ) / denom
      z(k-1) = ( ak * bkm1 - bk ) / denom

    end if

    k = k - ks
    ik = ik - k
    if ( ks == 2 ) then
      ik = ik - ( k + 1 )
    end if

  end do

  z(1:n) = z(1:n) / dasum ( n, z, 1 )
!
!  Solve U' * Y = W.
!
  k = 1
  ik = 0

  do while ( k <= n )

    if ( kpvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= 1 ) then

      z(k) = z(k) + ddot ( k-1, ap(ik+1), 1, z(1), 1 )
      ikp1 = ik + k

      if ( ks == 2 ) then
        z(k+1) = z(k+1) + ddot ( k-1, ap(ikp1+1), 1, z(1), 1 )
      end if

      kp = abs ( kpvt(k) )

      if ( kp /= k ) then
        t = z(k)
        z(k) = z(kp)
        z(kp) = t
      end if

    end if

    ik = ik + k
    if ( ks == 2 ) then
      ik = ik + ( k + 1 )
    end if
    k = k + ks

  end do

  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
  ynorm = 1.0D+00
!
!  Solve U * D * V = Y.
!
  k = n

  ik = ( n * ( n - 1 ) ) / 2

  do while ( 0 < k )

    kk = ik + k
    ikm1 = ik - ( k - 1 )

    if ( kpvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= ks ) then

      kp = abs ( kpvt(k) )
      kps = k + 1 - ks

      if ( kp /= kps ) then
        t = z(kps)
        z(kps) = z(kp)
        z(kp) = t
      end if

      call daxpy ( k-ks, z(k), ap(ik+1), 1, z(1), 1 )

      if ( ks == 2 ) then
        call daxpy ( k-ks, z(k-1), ap(ikm1+1), 1, z(1), 1 )
      end if

    end if

    if ( ks /= 2 ) then

      if ( abs ( ap(kk) ) < abs ( z(k) ) ) then
        s = abs ( ap(kk) ) / abs ( z(k) )
        z(1:n) = s * z(1:n)
        ynorm = s * ynorm
      end if

      if ( ap(kk) /= 0.0D+00 ) then
        z(k) = z(k) / ap(kk)
      else
        z(k) = 1.0D+00
      end if

    else

      km1k = ik + k - 1
      km1km1 = ikm1 + k - 1
      ak = ap(kk) / ap(km1k)
      akm1 = ap(km1km1) / ap(km1k)
      bk = z(k) / ap(km1k)
      bkm1 = z(k-1) / ap(km1k)
      denom = ak * akm1 - 1.0D+00
      z(k) = ( akm1 * bk - bkm1 ) / denom
      z(k-1) = ( ak * bkm1 - bk ) / denom

    end if

    k = k - ks
    ik = ik - k
    if ( ks == 2 ) then
      ik = ik - ( k + 1 )
    end if

  end do

  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm
!
!  Solve U' * Z = V.
!
  k = 1
  ik = 0

  do while ( k <= n )

    if ( kpvt(k) < 0 ) then
      ks = 2
    else
      ks = 1
    end if

    if ( k /= 1 ) then

      z(k) = z(k) + ddot ( k-1, ap(ik+1), 1, z(1), 1 )
      ikp1 = ik + k

      if ( ks == 2 ) then
        z(k+1) = z(k+1) + ddot ( k-1, ap(ikp1+1), 1, z(1), 1 )
      end if

      kp = abs ( kpvt(k) )

      if ( kp /= k ) then
        t = z(k)
        z(k) = z(kp)
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
!  Make ZNORM = 1.0.
!
  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm

  if ( anorm /= 0.0D+00 ) then
    rcond = ynorm / anorm
  else
    rcond = 0.0D+00
  end if

  return
end
subroutine dspdi ( ap, n, kpvt, det, inert, work, job )

!*****************************************************************************80
!
!! DSPDI: determinant, inertia and inverse of a real symmetric matrix.
!
!  Discussion:
!
!    DSPDI uses the factors from DSPFA, where the matrix is stored in
!    packed form.
!
!    A division by zero will occur if the inverse is requested
!    and DSPCO has set RCOND == 0.0D+00 or DSPFA has set INFO /= 0.
!
!    Variables not requested by JOB are not used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input/output, real ( kind = 8 ) AP((N*(N+1))/2).  On input, the output from
!    DSPFA.  On output, the upper triangle of the inverse of the original
!    matrix, stored in packed form, if requested.  The columns of the upper
!    triangle are stored sequentially in a one-dimensional array.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) KPVT(N), the pivot vector from DSPFA.
!
!    Output, real ( kind = 8 ) DET(2), the determinant of the original matrix,
!    if requested.
!      determinant = DET(1) * 10.0**DET(2)
!    with 1.0D+00 <= abs ( DET(1) ) < 10.0D+00 or DET(1) = 0.0.
!
!    Output, integer ( kind = 4 ) INERT(3), the inertia of the original matrix,
!    if requested.
!    INERT(1) = number of positive eigenvalues.
!    INERT(2) = number of negative eigenvalues.
!    INERT(3) = number of zero eigenvalues.
!
!    Workspace, real ( kind = 8 ) WORK(N).
!
!    Input, integer ( kind = 4 ) JOB, has the decimal expansion ABC where:
!      if A /= 0, the inertia is computed,
!      if B /= 0, the determinant is computed,
!      if C /= 0, the inverse is computed.
!    For example, JOB = 111  gives all three.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) ak
  real ( kind = 8 ) akkp1
  real ( kind = 8 ) akp1
  real ( kind = 8 ) ap((n*(n+1))/2)
  real ( kind = 8 ) d
  real ( kind = 8 ) det(2)
  logical dodet
  logical doert
  logical doinv
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ikp1
  integer ( kind = 4 ) iks
  integer ( kind = 4 ) inert(3)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jb
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) jkp1
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) kkp1
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) kpvt(n)
  integer ( kind = 4 ) ks
  integer ( kind = 4 ) ksj
  integer ( kind = 4 ) kskp1
  integer ( kind = 4 ) kstep
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t
  real ( kind = 8 ) temp
  real ( kind = 8 ) work(n)

  doinv = mod ( job,   10 )       /= 0
  dodet = mod ( job,  100 ) /  10 /= 0
  doert = mod ( job, 1000 ) / 100 /= 0

  if ( dodet .or. doert ) then

    if ( doert ) then
      inert(1:3) = 0
    end if

    if ( dodet ) then
      det(1:2) = (/ 1.0D+00, 0.0D+00 /)
    end if

    t = 0.0D+00
    ik = 0

    do k = 1, n

      kk = ik + k
      d = ap(kk)
!
!  2 by 2 block
!  use det (d  s)  =  (d/t * c - t) * t,  t = abs ( s )
!          (s  c)
!  to avoid underflow/overflow troubles.
!
!  Take two passes through scaling.  Use T for flag.
!
      if ( kpvt(k) <= 0 ) then

        if ( t == 0.0D+00 ) then
          ikp1 = ik + k
          kkp1 = ikp1 + k
          t = abs ( ap(kkp1) )
          d = ( d / t ) * ap(kkp1+1) - t
        else
          d = t
          t = 0.0D+00
        end if

      end if

      if ( doert ) then
        if ( 0.0D+00 < d ) then
          inert(1) = inert(1) + 1
        else if ( d < 0.0D+00 ) then
          inert(2) = inert(2) + 1
        else if ( d == 0.0D+00 ) then
          inert(3) = inert(3) + 1
        end if
      end if

      if ( dodet ) then

        det(1) = det(1) * d

        if ( det(1) /= 0.0D+00 ) then

          do while ( abs ( det(1) ) < 1.0D+00 )
            det(1) = det(1) * 10.0D+00
            det(2) = det(2) - 1.0D+00
          end do

          do while ( 10.0D+00 <= abs ( det(1) ) )
            det(1) = det(1) / 10.0D+00
            det(2) = det(2) + 1.0D+00
          end do

        end if

      end if

      ik = ik + k

    end do

  end if
!
!  Compute inverse(A).
!
  if ( doinv ) then

    k = 1
    ik = 0

    do while ( k <= n )

      km1 = k - 1
      kk = ik + k
      ikp1 = ik + k
      kkp1 = ikp1 + k

      if ( 0 <= kpvt(k) ) then
!
!  1 by 1.
!
        ap(kk) = 1.0D+00 / ap(kk)

        if ( 2 <= k ) then

          work(1:k-1) = ap(ik+1:ik+k-1)

          ij = 0

          do j = 1, k-1
            jk = ik + j
            ap(jk) = ddot ( j, ap(ij+1), 1, work, 1 )
            call daxpy ( j-1, work(j), ap(ij+1), 1, ap(ik+1), 1 )
            ij = ij + j
          end do

          ap(kk) = ap(kk) + ddot ( k-1, work, 1, ap(ik+1), 1 )

        end if

        kstep = 1

      else
!
!  2 by 2.
!
        t = abs ( ap(kkp1) )
        ak = ap(kk) / t
        akp1 = ap(kkp1+1) / t
        akkp1 = ap(kkp1) / t
        d = t * ( ak * akp1 - 1.0D+00 )
        ap(kk) = akp1 / d
        ap(kkp1+1) = ak / d
        ap(kkp1) = -akkp1 / d

        if ( 1 <= km1 ) then

          work(1:km1) = ap(ikp1+1:ikp1+km1)

          ij = 0

          do j = 1, km1
            jkp1 = ikp1 + j
            ap(jkp1) = ddot ( j, ap(ij+1), 1, work, 1 )
            call daxpy ( j-1, work(j), ap(ij+1), 1, ap(ikp1+1), 1 )
            ij = ij + j
          end do

          ap(kkp1+1) = ap(kkp1+1) + ddot ( km1, work, 1, ap(ikp1+1), 1 )
          ap(kkp1) = ap(kkp1) + ddot ( km1, ap(ik+1), 1, ap(ikp1+1), 1 )

          work(1:km1) = ap(ik+1:ik+km1)

          ij = 0

          do j = 1, km1
            jk = ik + j
            ap(jk) = ddot ( j, ap(ij+1), 1, work, 1 )
            call daxpy ( j-1, work(j), ap(ij+1), 1, ap(ik+1), 1 )
            ij = ij + j
          end do

          ap(kk) = ap(kk) + ddot ( km1, work, 1, ap(ik+1), 1 )

        end if

        kstep = 2

      end if
!
!  Swap.
!
      ks = abs ( kpvt(k) )

      if ( ks /= k ) then

        iks = ( ks * ( ks - 1 ) ) / 2
        call dswap ( ks, ap(iks+1), 1, ap(ik+1), 1 )
        ksj = ik + ks

        do jb = ks, k
          j = k + ks - jb
          jk = ik + j
          temp = ap(jk)
          ap(jk) = ap(ksj)
          ap(ksj) = temp
          ksj = ksj - ( j - 1 )
        end do

        if ( kstep /= 1 ) then
          kskp1 = ikp1 + ks
          temp = ap(kskp1)
          ap(kskp1) = ap(kkp1)
          ap(kkp1) = temp
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
subroutine dspfa ( ap, n, kpvt, info )

!*****************************************************************************80
!
!! DSPFA factors a real symmetric matrix stored in packed form.
!
!  Discussion:
!
!    To solve A*X = B, follow DSPFA by DSPSL.
!
!    To compute inverse(A)*C, follow DSPFA by DSPSL.
!
!    To compute determinant(A), follow DSPFA by DSPDI.
!
!    To compute inertia(A), follow DSPFA by DSPDI.
!
!    To compute inverse(A), follow DSPFA by DSPDI.
!
!  Packed storage:
!
!    The following program segment will pack the upper triangle of a
!    symmetric matrix.
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
!    16 May 2005
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
!    Input/output, real ( kind = 8 ) AP((N*(N+1))/2).  On input, the packed
!    form of a symmetric matrix A.  The columns of the upper triangle are stored
!    sequentially in a one-dimensional array.  On output, a block diagonal
!    matrix and the multipliers which were used to obtain it stored in
!    packed form.  The factorization can be written A = U*D*U' where U
!    is a product of permutation and unit upper triangular matrices, U'
!    is the transpose of U, and D is block diagonal with 1 by 1 and 2
!    by 2 blocks.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) KPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, normal value.
!    K, if the K-th pivot block is singular.  This is not an error
!    condition for this subroutine, but it does indicate that DSPSL or
!    DSPDI may divide by zero if called.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) absakk
  real ( kind = 8 ) ak
  real ( kind = 8 ) akm1
  real ( kind = 8 ) alpha
  real ( kind = 8 ) ap((n*(n+1))/2)
  real ( kind = 8 ) bk
  real ( kind = 8 ) bkm1
  real ( kind = 8 ) colmax
  real ( kind = 8 ) denom
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ikm1
  integer ( kind = 4 ) im
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imim
  integer ( kind = 4 ) imj
  integer ( kind = 4 ) imk
  integer ( kind = 4 ) info
  integer ( kind = 4 ) idamax
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
  integer ( kind = 4 ) kpvt(n)
  integer ( kind = 4 ) kstep
  real ( kind = 8 ) mulk
  real ( kind = 8 ) mulkm1
  real ( kind = 8 ) rowmax
  logical swap
  real ( kind = 8 ) t
!
!  ALPHA is used in choosing pivot block size.
!
  alpha = ( 1.0D+00 + sqrt ( 17.0D+00 ) ) / 8.0D+00

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
      kpvt(1) = 1
      if ( ap(1) == 0.0D+00 ) then
        info = 1
      end if
      exit
    end if
!
!  This section of code determines the kind of elimination to be performed.
!  When it is completed, KSTEP will be set to the size of the pivot block,
!  and SWAP will be set to TRUE if an interchange is required.
!
    km1 = k - 1
    kk = ik + k
    absakk = abs ( ap(kk) )
!
!  Determine the largest off-diagonal element in column K.
!
    imax = idamax ( k-1, ap(ik+1), 1 )
    imk = ik + imax
    colmax = abs ( ap(imk) )

    if ( alpha * colmax <= absakk ) then

      kstep = 1
      swap = .false.
!
!  Determine the largest off-diagonal element in row IMAX.
!
    else

      rowmax = 0.0D+00
      im = ( imax * ( imax - 1 ) ) / 2
      imj = im + 2 * imax

      do j = imax+1, k
        rowmax = max ( rowmax, abs ( ap(imj) ) )
        imj = imj + j
      end do

      if ( imax /= 1 ) then
        jmax = idamax ( imax-1, ap(im+1), 1 )
        jmim = jmax + im
        rowmax = max ( rowmax, abs ( ap(jmim) ) )
      end if

      imim = imax + im

      if ( alpha * rowmax <= abs ( ap(imim) ) ) then
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
    if ( max ( absakk, colmax ) == 0.0D+00 ) then

      kpvt(k) = k
      info = k

    else

      if ( kstep /= 2 ) then
!
!  1 x 1 pivot block.
!
        if ( swap ) then
!
!  Perform an interchange.
!
          call dswap ( imax, ap(im+1), 1, ap(ik+1), 1 )
          imj = ik + imax

          do jj = imax, k
            j = k + imax - jj
            jk = ik + j
            t = ap(jk)
            ap(jk) = ap(imj)
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
          call daxpy ( j, t, ap(ik+1), 1, ap(ij+1), 1 )
          ap(jk) = mulk
          ij = ij - ( j - 1 )
        end do
!
!  Set the pivot array.
!
        if ( swap ) then
          kpvt(k) = imax
        else
          kpvt(k) = k
        end if

      else
!
!  2 x 2 pivot block.
!
        km1k = ik + k - 1
        ikm1 = ik - ( k - 1 )
!
!  Perform an interchange.
!
        if ( swap ) then

          call dswap ( imax, ap(im+1), 1, ap(ikm1+1), 1 )
          imj = ikm1 + imax

          do jj = imax, km1
            j = km1 + imax - jj
            jkm1 = ikm1 + j
            t = ap(jkm1)
            ap(jkm1) = ap(imj)
            ap(imj) = t
            imj = imj - ( j - 1 )
          end do

          t = ap(km1k)
          ap(km1k) = ap(imk)
          ap(imk) = t

        end if
!
!  Perform the elimination.
!
        if ( k /= 2 ) then

          ak = ap(kk) / ap(km1k)
          km1km1 = ikm1 + k - 1
          akm1 = ap(km1km1) / ap(km1k)
          denom = 1.0D+00 - ak * akm1
          ij = ik - ( k - 1 ) - ( k - 2 )

          do jj = 1, k-2

            j = km1 - jj
            jk = ik + j
            bk = ap(jk) / ap(km1k)
            jkm1 = ikm1 + j
            bkm1 = ap(jkm1) / ap(km1k)
            mulk = ( akm1 * bk - bkm1 ) / denom
            mulkm1 = ( ak * bkm1 - bk ) / denom
            t = mulk
            call daxpy ( j, t, ap(ik+1), 1, ap(ij+1), 1 )
            t = mulkm1
            call daxpy ( j, t, ap(ikm1+1), 1, ap(ij+1), 1 )
            ap(jk) = mulk
            ap(jkm1) = mulkm1
            ij = ij - ( j - 1 )
          end do

        end if
!
!  Set the pivot array.
!
        if ( swap ) then
          kpvt(k) = -imax
        else
          kpvt(k) = 1 - k
        end if

        kpvt(k-1) = kpvt(k)

      end if

    end if

    ik = ik - ( k - 1 )
    if ( kstep == 2 ) then
      ik = ik - ( k - 2 )
    end if

    k = k - kstep

  end do

  return
end
subroutine dspsl ( ap, n, kpvt, b )

!*****************************************************************************80
!
!! DSPSL solves the real symmetric system factored by DSPFA.
!
!  Discussion:
!
!    To compute inverse(A) * C where C is a matrix with P columns:
!
!      call dspfa ( ap, n, kpvt, info )
!
!      if ( info /= 0 ) go to ...
!
!      do j = 1, p
!        call dspsl ( ap, n, kpvt, c(1,j) )
!      end do
!
!    A division by zero may occur if DSPCO has set RCOND == 0.0D+00
!    or DSPFA has set INFO /= 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input, real ( kind = 8 ) AP((N*(N+1))/2), the output from DSPFA.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) KPVT(N), the pivot vector from DSPFA.
!
!    Input/output, real ( kind = 8 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) ak
  real ( kind = 8 ) akm1
  real ( kind = 8 ) ap((n*(n+1))/2)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) bk
  real ( kind = 8 ) bkm1
  real ( kind = 8 ) denom
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ikm1
  integer ( kind = 4 ) ikp1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) km1k
  integer ( kind = 4 ) km1km1
  integer ( kind = 4 ) kp
  integer ( kind = 4 ) kpvt(n)
  real ( kind = 8 ) ddot
  real ( kind = 8 ) temp
!
!  Loop backward applying the transformations and D inverse to B.
!
  k = n
  ik = ( n * ( n - 1 ) ) / 2

  do while ( 0 < k )

    kk = ik + k

    if ( 0 <= kpvt(k) ) then
!
!  1 x 1 pivot block.
!
      if ( k /= 1 ) then

        kp = kpvt(k)
!
!  Interchange.
!
        if ( kp /= k ) then
          temp = b(k)
          b(k) = b(kp)
          b(kp) = temp
        end if
!
!  Apply the transformation.
!
        call daxpy ( k-1, b(k), ap(ik+1), 1, b(1), 1 )

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

        kp = abs ( kpvt(k) )
!
!  Interchange.
!
        if ( kp /= k-1 ) then
          temp = b(k-1)
          b(k-1) = b(kp)
          b(kp) = temp
        end if
!
!  Apply the transformation.
!
        call daxpy ( k-2, b(k), ap(ik+1), 1, b(1), 1 )
        call daxpy ( k-2, b(k-1), ap(ikm1+1), 1, b(1), 1 )

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
      denom = ak * akm1 - 1.0D+00
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

    if ( 0 <= kpvt(k) ) then
!
!  1 x 1 pivot block.
!
      if ( k /= 1 ) then
!
!  Apply the transformation.
!
        b(k) = b(k) + ddot ( k-1, ap(ik+1), 1, b(1), 1 )
        kp = kpvt(k)
!
!  Interchange.
!
        if ( kp /= k ) then
          temp = b(k)
          b(k) = b(kp)
          b(kp) = temp
        end if

      end if

      ik = ik + k
      k = k + 1

    else
!
!  2 x 2 pivot block.
!
      if ( k /= 1 ) then
!
!  Apply the transformation.
!
        b(k) = b(k) + ddot ( k-1, ap(ik+1), 1, b(1), 1 )
        ikp1 = ik + k
        b(k+1) = b(k+1) + ddot ( k-1, ap(ikp1+1), 1, b(1), 1 )
        kp = abs ( kpvt(k) )
!
!  Interchange.
!
        if ( kp /= k ) then
          temp = b(k)
          b(k) = b(kp)
          b(kp) = temp
        end if

      end if

      ik = ik + k + k + 1
      k = k + 2

    end if

  end do

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

        do j = l+1, n
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

      do lls = l+1, mn+1

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
subroutine dtrco ( t, ldt, n, rcond, z, job )

!*****************************************************************************80
!
!! DTRCO estimates the condition of a real triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
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
!    Input, real ( kind = 8 ) T(LDT,N), the triangular matrix.  The zero
!    elements of the matrix are not referenced, and the corresponding
!    elements of the array can be used to store other information.
!
!    Input, integer ( kind = 4 ) LDT, the leading dimension of the array T.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) RCOND, an estimate of the reciprocal condition
!    of T.  For the system T*X = B, relative perturbations in T and B of size
!    EPSILON may cause relative perturbations in X of size EPSILON/RCOND.
!    If RCOND is so small that the logical expression
!      1.0D+00 + RCOND == 1.0D+00
!    is true, then T may be singular to working precision.  In particular,
!    RCOND is zero if exact singularity is detected or the estimate underflows.
!
!    Output, real ( kind = 8 ) Z(N) a work vector whose contents are usually
!    unimportant.  If T is close to a singular matrix, then Z is an
!    approximate null vector in the sense that
!      norm(A*Z) = RCOND * norm(A) * norm(Z).
!
!    Input, integer ( kind = 4 ) JOB, indicates the shape of T:
!    0, T is lower triangular.
!    nonzero, T is upper triangular.
!
  implicit none

  integer ( kind = 4 ) ldt
  integer ( kind = 4 ) n

  real ( kind = 8 ) ek
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  logical lower
  real ( kind = 8 ) rcond
  real ( kind = 8 ) s
  real ( kind = 8 ) dasum
  real ( kind = 8 ) sm
  real ( kind = 8 ) t(ldt,n)
  real ( kind = 8 ) tnorm
  real ( kind = 8 ) w
  real ( kind = 8 ) wk
  real ( kind = 8 ) wkm
  real ( kind = 8 ) ynorm
  real ( kind = 8 ) z(n)

  lower = job == 0
!
!  Compute the 1-norm of T.
!
  tnorm = 0.0D+00

  do j = 1, n

    if ( lower ) then
      l = n + 1 - j
      i1 = j
    else
      l = j
      i1 = 1
    end if

    tnorm = max ( tnorm, dasum ( l, t(i1,j), 1 ) )

  end do
!
!  RCOND = 1/(norm(T)*(estimate of norm(inverse(T)))).
!
!  Estimate = norm(Z)/norm(Y) where T * Z = Y and T' * Y = E.
!
!  T' is the transpose of T.
!
!  The components of E are chosen to cause maximum local
!  growth in the elements of Y.
!
!  The vectors are frequently rescaled to avoid overflow.
!
!  Solve T' * Y = E.
!
  ek = 1.0D+00
  z(1:n) = 0.0D+00

  do kk = 1, n

    if ( lower ) then
      k = n + 1 - kk
    else
      k = kk
    end if

    if ( z(k) /= 0.0D+00 ) then
      ek = sign ( ek, -z(k) )
    end if

    if ( abs ( t(k,k) ) < abs ( ek - z(k) ) ) then
      s = abs ( t(k,k) ) / abs ( ek - z(k) )
      z(1:n) = s * z(1:n)
      ek = s * ek
    end if

    wk = ek - z(k)
    wkm = -ek - z(k)
    s = abs ( wk )
    sm = abs ( wkm )

    if ( t(k,k) /= 0.0D+00 ) then
      wk = wk / t(k,k)
      wkm = wkm / t(k,k)
    else
      wk = 1.0D+00
      wkm = 1.0D+00
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
        sm = sm + abs ( z(j) + wkm * t(k,j) )
        z(j) = z(j) + wk * t(k,j)
        s = s + abs ( z(j) )
      end do

      if ( s < sm ) then
        w = wkm - wk
        wk = wkm
        do j = j1, j2
          z(j) = z(j) + w * t(k,j)
        end do
      end if

    end if

    z(k) = wk

  end do

  z(1:n) = z(1:n) / dasum ( n, z, 1 )

  ynorm = 1.0D+00
!
!  Solve T * Z = Y.
!
  do kk = 1, n

    if ( lower ) then
      k = kk
    else
      k = n + 1 - kk
    end if

    if ( abs ( t(k,k) ) < abs ( z(k) ) ) then
      s = abs ( t(k,k) ) / abs ( z(k) )
      z(1:n) = s * z(1:n)
      ynorm = s * ynorm
    end if

    if ( t(k,k) /= 0.0D+00 ) then
      z(k) = z(k) / t(k,k)
    else
      z(k) = 1.0D+00
    end if

    if ( lower ) then
      i1 = k + 1
    else
      i1 = 1
    end if

    if ( kk < n ) then
      w = -z(k)
      call daxpy ( n-kk, w, t(i1,k), 1, z(i1), 1 )
    end if

  end do
!
!  Make ZNORM = 1.0.
!
  s = 1.0D+00 / dasum ( n, z, 1 )
  z(1:n) = s * z(1:n)
  ynorm = s * ynorm

  if ( tnorm /= 0.0D+00 ) then
    rcond = ynorm / tnorm
  else
    rcond = 0.0D+00
  end if

  return
end
subroutine dtrdi ( t, ldt, n, det, job, info )

!*****************************************************************************80
!
!! DTRDI computes the determinant and inverse of a real triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2001
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
!    Input/output, real ( kind = 8 ) T(LDT,N).
!    On input, T contains the triangular matrix.  The zero elements of the
!    matrix are not referenced, and the corresponding elements of the array
!    can be used to store other information.
!    On output, T contains the inverse matrix, if requested.
!
!    Input, integer ( kind = 4 ) LDT, the leading dimension of T.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, real ( kind = 8 ) DET(2), the determinant of the matrix, if
!    requested.  The determinant = DET(1) * 10.0**DET(2), with
!    1.0 <= abs ( DET(1) ) < 10.0, or DET(1) == 0.
!
!    Input, integer ( kind = 4 ) JOB, specifies the shape of T, and the task.
!    010, inverse of lower triangular matrix.
!    011, inverse of upper triangular matrix.
!    100, determinant only.
!    110, determinant and inverse of lower triangular.
!    111, determinant and inverse of upper triangular.
!
!    Output, integer ( kind = 4 ) INFO.
!    If the inverse was requested, then
!    0, if the system was nonsingular;
!    nonzero, if the system was singular.
!
  implicit none

  integer ( kind = 4 ) ldt
  integer ( kind = 4 ) n

  real ( kind = 8 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  real ( kind = 8 ) t(ldt,n)
  real ( kind = 8 ) temp

  info = 0
!
!  Determinant.
!
  if ( job / 100 /= 0 ) then

    det(1) = 1.0D+00
    det(2) = 0.0D+00

    do i = 1, n

      det(1) = det(1) * t(i,i)

      if ( det(1) == 0.0D+00 ) then
        exit
      end if

      do while ( abs ( det(1) ) < 1.0D+00 )
        det(1) = det(1) * 10.0D+00
        det(2) = det(2) - 1.0D+00
      end do

      do while ( 10.0D+00 <= abs ( det(1) ) )
        det(1) = det(1) / 10.0D+00
        det(2) = det(2) + 1.0D+00
      end do

    end do

  end if

  if ( mod ( job / 10, 10 ) == 0 ) then
    return
  end if
!
!  Inverse of an upper triangular matrix.
!
  if ( mod ( job, 10 ) /= 0 ) then

    info = 0

    do k = 1, n

      if ( t(k,k) == 0.0D+00 ) then
        info = k
        exit
      end if

      t(k,k) = 1.0D+00 / t(k,k)
      t(1:k-1,k) = -t(1:k-1,k) * t(k,k)

      do j = k + 1, n
        temp = t(k,j)
        t(k,j) = 0.0D+00
        call daxpy ( k, temp, t(1,k), 1, t(1,j), 1 )
      end do

    end do
!
!  Inverse of a lower triangular matrix.
!
  else

    info = 0

    do k = n, 1, -1

      if ( t(k,k) == 0.0D+00 ) then
        info = k
        exit
      end if

      t(k,k) = 1.0D+00 / t(k,k)
      temp = -t(k,k)

      if ( k /= n ) then
        call dscal ( n-k, temp, t(k+1,k), 1 )
      end if

      do j = 1, k-1
        temp = t(k,j)
        t(k,j) = 0.0D+00
        call daxpy ( n-k+1, temp, t(k,k), 1, t(k,j), 1 )
      end do

    end do

  end if

  return
end
subroutine dtrsl ( t, ldt, n, b, job, info )

!*****************************************************************************80
!
!! DTRSL solves triangular linear systems.
!
!  Discussion:
!
!    DTRSL can solve T * X = B or T' * X = B where T is a triangular
!    matrix of order N.
!
!    Here T' denotes the transpose of the matrix T.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 2005
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
!    Input, real ( kind = 8 ) T(LDT,N), the matrix of the system.  The zero
!    elements of the matrix are not referenced, and the corresponding
!    elements of the array can be used to store other information.
!
!    Input, integer ( kind = 4 ) LDT, the leading dimension of the array T.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
!    Input, integer ( kind = 4 ) JOB, specifies what kind of system is
!    to be solved:
!    00, solve T * X = B, T lower triangular,
!    01, solve T * X = B, T upper triangular,
!    10, solve T'* X = B, T lower triangular,
!    11, solve T'* X = B, T upper triangular.
!
!    Output, integer ( kind = 4 ) INFO, singularity indicator.
!    0, the system is nonsingular.
!    nonzero, the index of the first zero diagonal element of T.
!
  implicit none

  integer ( kind = 4 ) ldt
  integer ( kind = 4 ) n

  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) job
  integer ( kind = 4 ) task
  real ( kind = 8 ) t(ldt,n)
  real ( kind = 8 ) temp
!
!  Check for zero diagonal elements.
!
  do j = 1, n
    if ( t(j,j) == 0.0D+00 ) then
      info = j
      return
    end if
  end do

  info = 0
!
!  Determine the task and go to it.
!
  if ( mod ( job, 10 ) == 0 ) then
    task = 1
  else
    task = 2
  end if

  if ( mod ( job, 100 ) / 10 /= 0 ) then
    task = task + 2
  end if
!
!  Solve T * X = B for T lower triangular.
!
  if ( task == 1 ) then

    b(1) = b(1) / t(1,1)

    do j = 2, n
      temp = -b(j-1)
      call daxpy ( n-j+1, temp, t(j,j-1), 1, b(j), 1 )
      b(j) = b(j) / t(j,j)
    end do
!
!  Solve T * X = B for T upper triangular.
!
  else if ( task == 2 ) then

    b(n) = b(n) / t(n,n)

    do jj = 2, n
      j = n - jj + 1
      temp = -b(j+1)
      call daxpy ( j, temp, t(1,j+1), 1, b(1), 1 )
      b(j) = b(j) / t(j,j)
    end do
!
!  Solve T' * X = B for T lower triangular.
!
  else if ( task == 3 ) then

    b(n) = b(n) / t(n,n)

    do j = n-1, 1, -1
      b(j) = ( b(j) - dot_product ( t(j+1:n,j), b(j+1:n) ) ) / t(j,j)
    end do
!
!  Solve T' * X = B for T upper triangular.
!
  else if ( task == 4 ) then

    b(1) = b(1) / t(1,1)

    do j = 2, n
      b(j) = ( b(j) - dot_product ( t(1:j-1,j), b(1:j-1) ) ) / t(j,j)
    end do

  end if

  return
end
