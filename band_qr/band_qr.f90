subroutine dgbbqr2 ( m, n, ml, mu, a, lda, tau, work, info )

!*****************************************************************************80
!
!! DGBBQR2 QR factors an M by N band matrix in GB format, without blocking.
!
!  Discussion:
!
!    DGBBQR2 computes a QR factorization of a real M by N band matrix A
!    with lower band ML and upper band MU: A = Q * R.
!
!    A is stored as a packed band matrix.
!    Input matrix A has ML subdiagonals and MU superdiagonals.
!    Output matrix R has ML + MU superdiagonals.
!
!    Example of A with M = 8, N = 6, ML = 2, and MU = 1.
!    Left, input matrix; Right, output matrix.
!
!      x  x  0  0  0  0       r   r   r   r   0   0
!      x  x  x  0  0  0       v1  r   r   r   r   0
!      x  x  x  x  0  0       v1  v2  r   r   r   r
!      0  x  x  x  x  0       0   v2  v3  r   r   r
!      0  0  x  x  x  x       0   0   v3  v4  r   r
!      0  0  0  x  x  x       0   0   0   v4  v5  r
!      0  0  0  0  x  x       0   0   0   0   v5  v6
!      0  0  0  0  0  x       0   0   0   0   0   v6
!
!  Licensing:
!
!    Copyright (c) 2007, Universidad Jaume I de Castellon
!    All rights reserved.
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions 
!    are met:
!      * Redistributions of source code must retain the above copyright
!        notice, this list of conditions and the following disclaimer.
!      * Redistributions in binary form must reproduce the above copyright
!        notice, this list of conditions and the following disclaimer in the
!        documentation and/or other materials provided with the distribution.
!      * Neither the name of the <organization> nor the
!        names of its contributors may be used to endorse or promote 
!        products derived from this software without specific prior written 
!        permission.
!
!    This software is provided by <copyright holder> ''as is'' and any
!    express or implied warranties, including, but not limited to, the 
!    implied warranties of merchantability and fitness for a particular 
!    purpose are disclaimed.  In no event shall <copyright holder> be liable 
!    for any direct, indirect, incidental, special, exemplary, or 
!    consequential damages (including, but not limited to, procurement of 
!    substitute goods or services; loss of use, data, or profits; or 
!    business interruption) however caused and on any theory of liability, 
!    whether in contract, strict liability, or tort (including negligence 
!    or otherwise) arising in any way out of the use of this software, even 
!    if advised of the possibility of such damage.
!
!  Modified:
!
!    04 April 2010
!
!  Author:
!
!    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
!    Dept. de Ingenieria y Ciencia de Computadores, Univ. Jaume I
!    12.080 Castellon, Spain
!    {gquintan,remon,quintana}@icc.uji.es
!
!  Reference:
!
!    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
!    LAPACK-Style Codes for the QR Factorization of Banded Matrices,
!    To Appear.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    0 < M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    0 < N.
!
!    Input, integer ( kind = 4 ) ML, the number of nonzero subdiagonals.
!    0 < ML.
!
!    Input, integer ( kind = 4 ) MU, the number of nonzero superdiagonals.
!    0 < MU.
!
!    Input/output, real ( kind = 8 ) A(LDA,N).
!    A is stored as a packed matrix: dense A is M by N, whereas
!    band A is (2*ML+MU+1) by N.
!    On entry, the M by N matrix AB in band storage in rows ML+1 to
!    2*ML+MU+1; rows 1 to ML of the array need not be set.
!    This matrix has ML subdiagonals and MU superdiagonals with data.
!    On exit, the elements on and above the diagonal of the array
!    contain the upper band matrix R.
!    The output matrix R has ML+MU superdiagonals.
!    The elements within the band below the diagonal, with the array
!    TAU, represent the orthogonal matrix Q as a product of
!    elementary reflectors.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!    2*ML+MU+1 < LDA.
!
!    Output, real ( kind = 8 ) TAU(min(M,N)), the scalar factors of the 
!    elementary reflectors.
!
!    Workspace, real ( kind = 8 ) WORK(min(N,MU+ML)).
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0: successful exit
!    nonzero, if INFO = -I, the I-th argument had an illegal value.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) diag
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mh
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) nh
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) tau(*)
  real ( kind = 8 ) work(*)
  real ( kind = 8 ), parameter :: zero = 0.0D+00
!
!  Test the input arguments.
!
  info = 0

  if ( m < 0 ) then
    info = - 1
  else if ( n < 0 ) then
    info = - 2
  else if ( ml < 0 ) then
    info = - 3
  else if ( mu < 0 ) then
    info = - 4
  else if ( lda < 2 * ml + mu + 1 ) then
    info = - 6
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dgbbqr2', -info )
    return
  end if
!
!  Quick return if possible.
!
  mn = min ( m, n )
  if ( mn == 0 ) then
    return
  end if
!
!  Set rows 1:ML to zero.
!
  call dlaset ( 'all', ml, n, zero, zero, a, lda )

  do j = 1, mn

    mh = min ( ml + 1, m - j + 1 )
!
!  Generate reflector H(j) to annihilate A(j+1:j+mh-1,j).
!
    call dlarfg ( mh, a(ml+mu+1,j), a(ml+mu+1+min(1,ml),j), 1, tau(j) )

    nh = min ( n - j, mu + ml )
!
!  Apply reflector H(j) to rest of matrix: A(j:j+mh-1,j+1:j+nh).
!
    if ( 0 < nh ) then

      diag = a(ml+mu+1,j)

      a(ml+mu+1,j) = one

      call dlarf ( 'left', mh, nh, a(ml+mu+1,j), 1, tau(j), a(ml+mu,j+1), &
        lda - 1, work )

      a(ml+mu+1,j) = diag

    end if

  end do

  return
end
subroutine dgbbqrf ( nb, m, n, ml, mu, a, lda, tau, work, info )

!*****************************************************************************80
!
!! DGBBQRF QR factors an M by N band matrix in GB format, using blocking.
!
!  Discussion:
!
!    DGBBQRF computes a QR factorization of a real M by N band matrix A
!    with lower band ML and upper band MU: A = Q * R.
!
!    A is stored as a packed matrix: dense A is M by N, whereas
!    band A is (2*ML+MU+NB) by N.
!    Input matrix AB has ML subdiagonals and MU superdiagonals.
!    Output matrix R has ML+MU+NB-1 superdiagonals.
!
!    Example of A with M = 8, N = 6, ML = 2, MU = 1, and NB = 2.
!    Left, input matrix; Right, output matrix.
!
!      x  x  0  0  0  0       r   r   r   r   r   0
!      x  x  x  0  0  0       v1  r   r   r   r   r
!      x  x  x  x  0  0       v1  v2  r   r   r   r
!      0  x  x  x  x  0       0   v2  v3  r   r   r
!      0  0  x  x  x  x       0   0   v3  v4  r   r
!      0  0  0  x  x  x       0   0   0   v4  v5  r
!      0  0  0  0  x  x       0   0   0   0   v5  v6
!      0  0  0  0  0  x       0   0   0   0   0   v6
!
!  Licensing:
!
!    Copyright (c) 2007, Universidad Jaume I de Castellon
!    All rights reserved.
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions 
!    are met:
!      * Redistributions of source code must retain the above copyright
!        notice, this list of conditions and the following disclaimer.
!      * Redistributions in binary form must reproduce the above copyright
!        notice, this list of conditions and the following disclaimer in the
!        documentation and/or other materials provided with the distribution.
!      * Neither the name of the <organization> nor the
!        names of its contributors may be used to endorse or promote 
!        products derived from this software without specific prior written 
!        permission.
!
!    This software is provided by <copyright holder> ''as is'' and any
!    express or implied warranties, including, but not limited to, the 
!    implied warranties of merchantability and fitness for a particular 
!    purpose are disclaimed.  In no event shall <copyright holder> be liable 
!    for any direct, indirect, incidental, special, exemplary, or 
!    consequential damages (including, but not limited to, procurement of 
!    substitute goods or services; loss of use, data, or profits; or 
!    business interruption) however caused and on any theory of liability, 
!    whether in contract, strict liability, or tort (including negligence 
!    or otherwise) arising in any way out of the use of this software, even 
!    if advised of the possibility of such damage.
!
!  Modified:
!
!    04 April 2010
!
!  Author:
!
!    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
!    Dept. de Ingenieria y Ciencia de Computadores, Univ. Jaume I
!    12.080 Castellon, Spain
!    {gquintan,remon,quintana}@icc.uji.es
!
!  Reference:
!
!    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
!    LAPACK-Style Codes for the QR Factorization of Banded Matrices,
!    To Appear.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NB, the block size to use.
!    1 < NB.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the matrix.
!    0 < M.
!
!    Input, integer ( kind = 4 ) N, the number of columns in the matrix.
!    0 < N.
!
!    Input, integer ( kind = 4 ) ML, the number of nonzero subdiagonals.
!    0 < ML.
!
!    Input, integer ( kind = 4 ) MU, the number of nonzero superdiagonals.
!    0 < MU.
!
!    Input/output, real ( kind = 8 ) A(LDA,N).
!    A is stored as a packed matrix: dense A is M by N, whereas
!    band A is (2*ML+MU+NB) by N.
!    On entry, the M by N matrix A in band storage in rows ML+NB to
!    2*ML+MU+NB; rows 1 to ML+NB-1 of the array need not be set.
!    This matrix has ML subdiagonals and MU superdiagonals with data.
!    On exit, the elements on and above the diagonal of the array
!    contain the upper band matrix R.
!    The output matrix R has ML+MU+NB-1 superdiagonals.
!    The elements within the band below the diagonal, with the array
!    TAU, represent the orthogonal matrix Q as a product of
!    elementary reflectors.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!    max ( 1, M ) < LDA.
!
!    Output, real ( kind = 8 ) TAU(min(M,N)), the scalar factors of the 
!    elementary reflectors.
!
!    Workspace, real ( kind = 8 ) WORK(DIMWORK).
!    NB*NB + MIN(N,ML+MU)*NB + MIN(M,ML+NB)*NB < DIMWORK.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, no errors.
!    nonzero, if INFO = -I, the I-th argument had an illegal value.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) am
  integer ( kind = 4 ) aml
  integer ( kind = 4 ) amu
  integer ( kind = 4 ) an
  integer ( kind = 4 ) anb
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) info
  integer ( kind = 4 ) irwk
  integer ( kind = 4 ) it
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jb
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jlc
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) ncol
  real ( kind = 8 ) tau(*)
  real ( kind = 8 ) work(*)
  real ( kind = 8 ), parameter :: zero = 0.0D+00
!
!  Test the input arguments.
!
  info = 0

  if ( nb < 1 ) then
    info = -1
  else if ( m < 0 ) then
    info = -2
  else if ( n < 0 ) then
    info = -3
  else if ( ml < 0 ) then
    info = -4
  else if ( mu < 0 ) then
    info = -5
  else if ( lda < max ( 1, 2 * ml + mu + nb ) ) then
    info = -7
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dgbbqrf', -info )
    return
  end if
!
!  Quick return if possible.
!
  if ( m == 0 .or. n == 0 .or. ml == 0 ) then
    return
  end if
!
!  Adjust matrix sizes: AM, AN, AML, AMU, and ANB.
!
  am = min ( m, n + ml )
  an = min ( n, m + mu )
  aml = min ( m - 1, ml )
  amu = min ( n - 1, mu )
  anb = nb
  if ( aml < anb ) then
    anb = aml
  end if

  mn = min ( am, an )

  it = 1
  iv = it + anb * anb
  irwk = iv + min ( aml + anb, am ) * anb

  ncol = min ( an, am - aml )
!
!  Factorization of full band matrix A(:,1:ncol).
!
  do j = 1, ncol, anb

    jb = min ( anb, ncol - j + 1 )
!
!  Factorize block A(j:j+aml+jb-1,j:j+jb-1).
!
    call dgebqr2 ( aml + jb, jb, aml, amu, a(ml+mu+nb,j), lda - 1, &
      tau(j), work(irwk), info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DGBBQRF - Fatal error!'
      write ( *, '(a)' ) '  Error return from DGEBQR2.'
      write ( *, '(a,i8)' ) '  INFO = ', info
      stop
    end if

    if ( j + jb <= an ) then
!
!  WORK(IV) := Y (the lower part of Y is padded with zeros).
!
      do jj = 1, jb

        do ii = 1, aml + jj
          work(iv-1+(aml+jb)*(jj-1)+ii) = a(ml+mu+nb+ii-jj,j+jj-1)
        end do

        do ii = aml + jj + 1, aml + jb
          work(iv-1+(aml+jb)*(jj-1)+ii) = zero
        end do

      end do
!
!  Form the triangular factor T of the block reflector.
!
      call dlarft ( 'forward', 'columnwise', aml + jb, jb, work(iv), &
        aml + jb, tau(j), work(it), jb )
!
!  Apply block reflector to A(j:j+aml+jb-1,j+jb:an) from the left.
!
      jlc = min ( an, j + jb - 1 + aml + amu )

      call dlarfb ( 'left', 'transpose', 'forward', 'columnwise', &
        aml + jb, jlc - j - jb + 1, jb, work(iv), aml + jb, work(it), jb, &
        a(ml+mu+nb-jb,j+jb), lda - 1, work(irwk), jlc - j - jb + 1 )

    end if

  end do
!
!  Factorization of rectangular matrix A(:,ncol+1:mn).
!
  do j = ncol + 1, mn, anb

    jb = min ( anb, mn - j + 1 )
!
!  Factorize block A(j:am,j:j+jb-1).
!
    call dgeqr2 ( am - j + 1, jb, a(ml+mu+nb,j), lda - 1, tau(j), work(irwk), &
      info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DGBBQRF - Fatal error!'
      write ( *, '(a)' ) '  Error return from DGEQR2.'
      write ( *, '(a,i8)' ) '  INFO = ', info
      stop
    end if

    if ( j + jb <= an ) then
!
!  Form the triangular factor T of the block reflector.
!
      call dlarft ( 'forward', 'columnwise', am - j + 1, jb, a(ml+mu+nb,j), &
        lda - 1, tau(j), work(it), jb )
!
!  Apply block reflector to A(j:am,j+jb:an) from the left.
!
      call dlarfb ( 'left', 'transpose', 'forward', 'columnwise', &
        am - j + 1, an - j - jb + 1, jb, a(ml+mu+nb,j), lda - 1, work(it), &
        jb, a(ml+mu+nb-jb,j+jb), lda - 1, work(irwk), an - j - jb + 1 )

    end if

  end do

  return
end
subroutine dgbbqrs ( m, n, ml, mu, nrhs, a, lda, tau, b, ldb, work, lwork, &
  info )

!*****************************************************************************80
!
!! DGBBQRS solves A*X = B when A has been factored by DGBBQR2.
!
!  Discussion:
!
!    Solve the least squares problem
!      min || A*X - B ||
!    using the QR factorization
!      A = Q*R
!    computed by DGBBQR2.
!
!    Here, the matrix A is an M by N band matrix, which was stored in the
!    standard LINPACK/LAPACK "GB" format.
!
!    This matrix was QR-factored by DGBBQR2, which is able to take advantage
!    of the GB format, and the results of the factorization overwrote
!    A, with some extra information stored in TAU.
!
!    Now one or more linear systems A * x = b are to solved using
!    the QR procedure.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix A.
!    0 < M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A.
!    0 < N < M.
!
!    Input, integer ( kind = 4 ) ML, the number of nonzero subdiagonals.
!    0 < ML.
!
!    Input, integer ( kind = 4 ) MU, the number of nonzero superdiagonals.
!    0 < MU.
!
!    Input, integer ( kind = 4 ) NRHS, the number of columns of B.  
!    0 < NRHS.
!
!    Input, real ( kind = 8 ) A(LDA,N), part of the QR factorization
!    computed by DGBBQRF.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!    2*ML+MU+1 < LDA.
!
!    Input, real ( kind = 8 ) TAU(N), the scalar factor of the elementary
!    reflectors H, as computed by DGBBQRF.
!
!    Input/output, real ( kind = 8 ) B(LDB,NRHS).
!    On entry, the M by NRHS right hand side matrix B.
!    On exit, the N by NRHS solution matrix X.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of B.
!    M < LDB.
!
!    Workspace, real ( kind = 8 ) WORK(LWORK).
!
!    Input, integer ( kind = 4 ) LWORK, the length of the array WORK.  
!    LWORK must be at least NRHS, and should be at least NRHS*NB, where 
!    NB is the block size for this environment.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, no errors.
!    nonzero, if INFO = -I, the I-th argument had an illegal value.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldb
  integer ( kind = 4 ) lwork
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nrhs

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) aii
  real ( kind = 8 ) b(ldb,nrhs)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) m
  real ( kind = 8 ) tau(n)
  real ( kind = 8 ) temp
  real ( kind = 8 ) work(lwork)
!
!  B := Q' * B.
!
!  Q is a real orthogonal matrix of order M, 
!  the product of N elementary reflectors.
!
!  B is an M by NRHS vector.
! 
  do i = 1, n
!
!  H(I) is applied to B(I:M,1:NRHS).
!
    if ( tau(i) /= 0.0D+00 ) then

      aii = a(i-i+ml+mu+1,i)
      a(i-i+ml+mu+1,i) = 1.0D+00
!
!  work := B' * v.
!
      do j = 1, nrhs
        temp = 0.0D+00
        do i2 = i, min ( m, i + ml )
          temp = temp + b(i2,j) * a(i2-i+ml+mu+1,i)
        end do
        work(j) = temp
      end do
!
!  B := B - tau * v * work'.
!
      do j = 1, nrhs
        if ( work(j) /= 0.0D+00 )then
          temp = work(j)
          do i2 = i, min ( m, i + ml )
            b(i2,j) = b(i2,j) - tau(i) * a(i2-i+ml+mu+1,i) * temp
          end do
        end if
      end do

      a(i-i+ml+mu+1,i) = aii

    end if

  end do
!
!  Solve R*X = B(1:n,:).
!
  do j = 1, nrhs
    do k = n, 1, -1
      if ( b(k,j) /= 0.0D+00 )then
        b(k,j) = b(k,j) / a(k-k+ml+mu+1,k)
        do i = max ( 1, k - ml - mu ), k - 1
          b(i,j) = b(i,j) - b(k,j) * a(i-k+ml+mu+1,k) 
       end do
      end if
    end do
  end do

  return
end
subroutine dgebqr2 ( m, n, ml, mu, a, lda, tau, work, info )

!*****************************************************************************80
!
!! DGEBQR2 QR factors an M by N band matrix in GE format, with no blocking.
!
!  Discussion:
!
!    DGEBQR2 computes a QR factorization of a real M by N band matrix A
!    with lower band ML and upper band MU: A = Q * R.
!
!    A is stored as a general matrix.
!    Input matrix A has ML subdiagonals and MU superdiagonals.
!    Output matrix R has ML+MU superdiagonals.
!
!    Example of A with M = 8, N = 6, ML = 2, and MU = 1.
!    Left, input matrix; Right, output matrix.
!
!      x  x  0  0  0  0       r   r   r   r   0   0
!      x  x  x  0  0  0       v1  r   r   r   r   0
!      x  x  x  x  0  0       v1  v2  r   r   r   r
!      0  x  x  x  x  0       0   v2  v3  r   r   r
!      0  0  x  x  x  x       0   0   v3  v4  r   r
!      0  0  0  x  x  x       0   0   0   v4  v5  r
!      0  0  0  0  x  x       0   0   0   0   v5  v6
!      0  0  0  0  0  x       0   0   0   0   0   v6
!
!    The matrix Q is represented as a product of elementary reflectors
!
!      Q = H(1) H(2) . . . H(k), where k = min(m,n).
!
!    Each H(i) has the form
!
!      H(i) = I - tau * v * v'
!
!    where tau is a real scalar, and v is a real vector with
!    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!    and tau in TAU(i).
!
!  Licensing:
!
!    Copyright (c) 2007, Universidad Jaume I de Castellon
!    All rights reserved.
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions 
!    are met:
!      * Redistributions of source code must retain the above copyright
!        notice, this list of conditions and the following disclaimer.
!      * Redistributions in binary form must reproduce the above copyright
!        notice, this list of conditions and the following disclaimer in the
!        documentation and/or other materials provided with the distribution.
!      * Neither the name of the <organization> nor the
!        names of its contributors may be used to endorse or promote 
!        products derived from this software without specific prior written 
!        permission.
!
!    This software is provided by <copyright holder> ''as is'' and any
!    express or implied warranties, including, but not limited to, the 
!    implied warranties of merchantability and fitness for a particular 
!    purpose are disclaimed.  In no event shall <copyright holder> be liable 
!    for any direct, indirect, incidental, special, exemplary, or 
!    consequential damages (including, but not limited to, procurement of 
!    substitute goods or services; loss of use, data, or profits; or 
!    business interruption) however caused and on any theory of liability, 
!    whether in contract, strict liability, or tort (including negligence 
!    or otherwise) arising in any way out of the use of this software, even 
!    if advised of the possibility of such damage.
!
!  Modified:
!
!    04 April 2010
!
!  Author:
!
!    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
!    Dept. de Ingenieria y Ciencia de Computadores, Univ. Jaume I
!    12.080 Castellon, Spain
!    {gquintan,remon,quintana}@icc.uji.es
!
!  Reference:
!
!    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
!    LAPACK-Style Codes for the QR Factorization of Banded Matrices,
!    To Appear.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    0 < M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    0 < N.
!
!    Input, integer ( kind = 4 ) ML, the number of nonzero subdiagonals.
!    0 < ML.
!
!    Input, integer ( kind = 4 ) MU, the number of nonzero superdiagonals.
!    0 < MU.
!
!    Input/output, real ( kind = 8 )(LDA,N).
!    On entry, the M by N matrix A. It has ML subdiagonals and MU
!    superdiagonals.
!    On exit, the elements on and above the diagonal of the array
!    contain the min(M,N) by N upper band matrix R.
!    The output matrix R has ML+MU superdiagonals.
!    The elements within the band below the diagonal, with the array
!    TAU, represent the orthogonal matrix Q as a product of
!    elementary reflectors.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.  
!    max ( 1, M ) < LDA.
!
!    Output, real ( kind = 8 ) TAU(min(M,N)),
!    the scalar factors of the elementary reflectors.
!
!    Workspace, real ( kind = 8 ) WORK(min(N,MU+ML)).
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, no errors.
!    nonzero, if INFO = -I, the I-th argument had an illegal value.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) diag
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mh
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) nh
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) tau(*)
  real ( kind = 8 ) work(*)
!
!  Test the input arguments.
!
  info = 0

  if ( m < 0 ) then
    info = -1
  else if ( n < 0 ) then
    info = -2
  else if ( ml < 0 ) then
    info = -3
  else if ( mu < 0 ) then
    info = -4
  else if ( lda < max ( 1, m ) ) then
    info = -6
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dgebqr2', -info )
    return
  end if
!
!  Quick return if possible.
!
  mn = min ( m, n )
  if ( mn == 0 ) then
    return
  end if

  do j = 1, mn

    mh = min ( ml + 1, m - j + 1 )
!
!  Generate reflector H(j) to annihilate A(j+1:j+mh-1,j).
!
    call dlarfg ( mh, a(j,j), a(min(j+1,m),j), 1, tau(j) )

    nh = min ( n - j, mu + ml )
!
!  Apply reflector H(j) to rest of matrix: A(j:j+mh-1,j+1:j+nh).
!
    if ( 0 < nh ) then

      diag = a(j,j)

      a(j,j) = one

      call dlarf ( 'left', mh, nh, a(j,j), 1, tau(j), a(j,j+1), lda, work )

      a(j,j) = diag

    end if

  end do

  return
end
subroutine dgebqrf ( nb, m, n, ml, mu, a, lda, tau, work, info )

!*****************************************************************************80
!
!! DGEBQRF QR factors an M by N band matrix stored in GE format, with blocking.
!
!  Discussion:
!
!    DGEBQRF computes a QR factorization of a real M by N band matrix A
!    with lower band ML and upper band MU: A = Q * R.
!
!    A is stored as a general matrix.
!    Input matrix A has ML subdiagonals and MU superdiagonals.
!    Output matrix R has ML+MU+NB-1 superdiagonals.
!
!    Example of A with M = 8, N = 6, ML = 2, MU = 1, and NB = 2.
!    Left, input matrix; Right, output matrix.
!
!      x  x  0  0  0  0       r   r   r   r   r   0
!      x  x  x  0  0  0       v1  r   r   r   r   r
!      x  x  x  x  0  0       v1  v2  r   r   r   r
!      0  x  x  x  x  0       0   v2  v3  r   r   r
!      0  0  x  x  x  x       0   0   v3  v4  r   r
!      0  0  0  x  x  x       0   0   0   v4  v5  r
!      0  0  0  0  x  x       0   0   0   0   v5  v6
!      0  0  0  0  0  x       0   0   0   0   0   v6
!
!  Licensing:
!
!    Copyright (c) 2007, Universidad Jaume I de Castellon
!    All rights reserved.
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions 
!    are met:
!      * Redistributions of source code must retain the above copyright
!        notice, this list of conditions and the following disclaimer.
!      * Redistributions in binary form must reproduce the above copyright
!        notice, this list of conditions and the following disclaimer in the
!        documentation and/or other materials provided with the distribution.
!      * Neither the name of the <organization> nor the
!        names of its contributors may be used to endorse or promote 
!        products derived from this software without specific prior written 
!        permission.
!
!    This software is provided by <copyright holder> ''as is'' and any
!    express or implied warranties, including, but not limited to, the 
!    implied warranties of merchantability and fitness for a particular 
!    purpose are disclaimed.  In no event shall <copyright holder> be liable 
!    for any direct, indirect, incidental, special, exemplary, or 
!    consequential damages (including, but not limited to, procurement of 
!    substitute goods or services; loss of use, data, or profits; or 
!    business interruption) however caused and on any theory of liability, 
!    whether in contract, strict liability, or tort (including negligence 
!    or otherwise) arising in any way out of the use of this software, even 
!    if advised of the possibility of such damage.
!
!  Modified:
!
!    04 April 2010
!
!  Author:
!
!    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
!    Dept. de Ingenieria y Ciencia de Computadores, Univ. Jaume I
!    12.080 Castellon, Spain
!    {gquintan,remon,quintana}@icc.uji.es
!
!  Reference:
!
!    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
!    LAPACK-Style Codes for the QR Factorization of Banded Matrices,
!    To Appear.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NB, the block size.
!    1 < NB.
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    0 < M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    0 < N.
!
!    Input, integer ( kind = 4 ) ML, the number of nonzero subdiagonals.
!    0 < ML.
!
!    Input, integer ( kind = 4 ) MU, the number of nonzero superdiagonals.
!    0 < MU.
!
!    Input/output, real ( kind = 8 ) A(LDA,N).
!    On entry, the M by N matrix A. It has ML subdiagonals and MU
!    superdiagonals.
!    On exit, the elements on and above the diagonal of the array
!    contain the min(M,N) by N upper band matrix R.
!    The output matrix R has ML+MU+NB-1 superdiagonals.
!    The elements within the band below the diagonal, with the array
!    TAU, represent the orthogonal matrix Q as a product of
!    elementary reflectors.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.  
!    max ( 1, M ) < LDA.
!
!    Output, real ( kind = 8 ) TAU(min(M,N)), the scalar factors of the 
!    elementary reflectors.
!
!    Workspace, real ( kind = 8 ) WORK(DIMWORK).
!    NB*NB + min(N,ML+MU)*NB + min(M,ML+NB)*NB < DIMWORK.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, no errors.
!    nonzero, if INFO = -I, the I-th argument had an illegal value.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) am
  integer ( kind = 4 ) aml
  integer ( kind = 4 ) amu
  integer ( kind = 4 ) an
  integer ( kind = 4 ) anb
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) info
  integer ( kind = 4 ) irwk
  integer ( kind = 4 ) it
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jb
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jlc
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) ncol
  real ( kind = 8 ) tau(*)
  real ( kind = 8 ) work(*)
  real ( kind = 8 ), parameter :: zero = 0.0D+00
!
!  Test the input arguments.
!
  info = 0

  if ( nb < 1 ) then
    info = -1
  else if ( m < 0 ) then
    info = -2
  else if ( n < 0 ) then
    info = -3
  else if ( ml < 0 ) then
    info = -4
  else if ( mu < 0 ) then
    info = -5
  else if ( lda < max ( 1, m ) ) then
    info = -7
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dgebqrf', -info )
    return
  end if
!
!  Quick return if possible.
!
  if ( m == 0 .or. n == 0 .or. ml == 0 ) then
    return
  end if
!
!  Adjust matrix sizes: AM, AN, AML, AMU, and ANB.
!
  am = min ( m, n + ml )
  an = min ( n, m + mu )
  aml = min ( m - 1, ml )
  amu = min ( n - 1, mu )
  anb = nb
  if ( aml < anb ) then
    anb = aml
  end if

  mn = min ( am, an )

  it = 1
  iv = it + anb * anb
  irwk = iv + min ( aml + anb, am ) * anb

  ncol = min ( an, am - aml )
!
!  Factorization of full band matrix A(:,1:ncol).
!
  do j = 1, ncol, anb

    jb = min ( anb, ncol - j + 1 )
!
!  Factorize block A(j:j+aml+jb-1,j:j+jb-1).
!
    call dgebqr2 ( aml + jb, jb, aml, amu, a(j,j), lda, tau(j), work(irwk), &
      info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DGEBQRF - Fatal error!'
      write ( *, '(a)' ) '  Error return from DGEBQR2.'
      write ( *, '(a,i8)' ) '  INFO = ', info
      stop
    end if

    if ( j + jb <= an ) then
!
!  WORK(IV) := Y (the lower part of Y is padded with zeros).
!
      do jj = 1, jb

        do ii = 1, aml + jj
          work(iv-1+(aml+jb)*(jj-1)+ii) = a(j+ii-1,j+jj-1)
        end do

        do ii = aml + jj + 1, aml + jb
          work(iv-1+(aml+jb)*(jj-1)+ii) = zero
        end do
      end do
!
!  Form the triangular factor T of the block reflector.
!
      call dlarft ( 'forward', 'columnwise', aml + jb, jb, work(iv), &
        aml + jb, tau(j), work(it), jb )
!
!  Apply block reflector to A(j:j+aml+jb-1,j+jb:an) from the left.
!
      jlc = min ( an, j + jb - 1 + aml + amu )

      call dlarfb ( 'left', 'transpose', 'forward', 'columnwise', aml + jb, &
        jlc - j - jb + 1, jb, work(iv), aml + jb, work(it), jb, a(j,j+jb), &
        lda, work(irwk), jlc - j - jb + 1 )

    end if

  end do
!
!  Factorization of rectangular matrix A(:,ncol+1:mn).
!
  do j = ncol + 1, mn, anb

    jb = min ( anb, mn - j + 1 )
!
!  Factorize block A(j:am,j:j+jb-1).
!
    call dgeqr2 ( am - j + 1, jb, a(j,j), lda, tau(j), work(irwk), info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DGEBQRF - Fatal error!'
      write ( *, '(a)' ) '  Error return from DGEQR2.'
      write ( *, '(a,i8)' ) '  INFO = ', info
      stop
    end if

    if ( j + jb <= an ) then
!
!  Form the triangular factor T of the block reflector.
!
      call dlarft ( 'forward', 'columnwise', am - j + 1, jb, a(j,j), lda, &
        tau(j), work(it), jb )
!
!  Apply block reflector to A(j:am,j+jb:an) from the left.
!
      call dlarfb ( 'left', 'transpose', 'forward', 'columnwise', &
        am - j + 1, an - j - jb + 1, jb, a(j,j), lda, work(it), jb, &
        a(j,j+jb), lda, work(irwk), an - j - jb + 1 )

    end if

  end do

  return
end
subroutine dgeqrs ( m, n, nrhs, a, lda, tau, b, ldb, work, lwork, info )

!*****************************************************************************80
!
!! DGEQRS solves a linear system factored by DGEQRF.
!
!  Discussion:
!
!    This routine is not part of the main LAPACK release.  It only has
!    the status of a "testing" routine.  We need direct access to this 
!    routine to make comparisons in source code and results.
!
!    This routine solves the least squares problem
!      min || A*X - B ||
!    using the QR factorization
!      A = Q*R
!    computed by DGEQRF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2010
!
!  Author:
!
!    This version by John Burkardt.
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Third Edition,
!    SIAM, 1999,
!    ISBN: 0898714478,
!    LC: QA76.73.F25L36
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix A.  
!    0 < M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A.  
!    0 < N < M.
!
!    Input, integer ( kind = 4 ) NRHS, the number of columns of B.  
!    0 < NRHS.
!
!    Input, real ( kind = 8 ) A(LDA,N)
!    Details of the QR factorization of the original matrix A as
!    returned by DGEQRF.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.  
!    M < LDA.
!
!    Input, real ( kind = 8 ) TAU(N)
!    Details of the orthogonal matrix Q.
!
!    Input/output, real ( kind = 8 ) B(LDB,NRHS),
!    On entry, the M by NRHS right hand side matrix B.
!    On exit, the N by NRHS solution matrix X.
!
!    Input, integer ( kind = 4 ) LDB, the leading dimension of the array B. 
!    M < LDB.
!
!    Workspace, real ( kind = 8 ) WORK(LWORK).
!
!    Input, integer ( kind = 4 ) LWORK,
!    The length of the array WORK.  LWORK must be at least NRHS,
!    and should be at least NRHS*NB, where NB is the block size
!    for this environment.
!
!    Output, integer ( kind = 4 ) INFO,
!    = 0: successful exit
!    < 0: if INFO = -I, the I-th argument had an illegal value
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldb
  integer ( kind = 4 ) lwork
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nrhs

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(ldb,nrhs)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) m
  real ( kind = 8 ), parameter :: one = 1.0D+00
  real ( kind = 8 ) tau(n)
  real ( kind = 8 ) work(lwork)

  info = 0

  if ( m < 0 ) then
    info = - 1
  else if ( n < 0 .or. m < n ) then
    info = - 2
  else if ( nrhs < 0 ) then
    info = - 3
  else if ( lda < max ( 1, m ) ) then
    info = - 5
  else if ( ldb < max ( 1, m ) ) then
    info = - 8
  else if ( ( lwork < 1 .or. lwork < nrhs ) .and. 0 < m .and. 0 < n ) then
     info = - 10
  end if

  if ( info /= 0 ) then
     call xerbla ( 'DGEQRS', - info )
     return
  end if
!
!  Quick return if possible.
!
  if ( n == 0 .or. nrhs == 0 .or. m == 0 ) then
    return
  end if
!
!  Compute B := Q' * B.
!
  call dormqr ( 'Left', 'Transpose', m, nrhs, n, a, lda, tau, b, ldb, &
    work, lwork, info )
!
!  Solve R * X = B(1:n,:).
!
  call dtrsm ( 'Left', 'Upper', 'No transpose', 'Non-unit', n, nrhs, &
    one, a, lda, b, ldb )

  return
end
subroutine dgeqrs_two ( m, n, nrhs, a, lda, tau, b, ldb, work, lwork, info )

!*****************************************************************************80
!
!! DGEQRS_TWO solves A*X = B when A has been factored by DGEQRF or DGEBQR2.
!
!  Discussion:
!
!    This routine solve the least squares problem
!      min || A*X - B ||
!    using the QR factorization
!      A = Q*R
!    computed by DGEQRF.
!
!    This routine is a revised version of the LAPACK testing code
!    DGEQRS.  It replaces all the calls to LAPACK subroutines by
!    explicit code.
!
!    It is intended to be a guide for devising a corresponding
!    routine to handle band matrices which have been QR factored
!    by DGBBQR2 or by DGBBQRF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Third Edition,
!    SIAM, 1999,
!    ISBN: 0898714478,
!    LC: QA76.73.F25L36
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix A.  
!    0 < M.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix A.
!    0 < N < M.
!
!    Input, integer ( kind = 4 ) NRHS, the number of columns of B.  
!    0 < NRHS.
!
!    Input, real ( kind = 8 ) A(LDA,N).
!    Details of the QR factorization of the original matrix A as
!    returned by DGEQRF.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.  
!    M < LDA.
!
!    Input, real ( kind = 8 ) TAU(N), the scalar factor of the elementary
!    reflectors.
!
!    Input/output, real ( kind = 8 ) B(LDB,NRHS).
!    On entry, the M by NRHS right hand side matrix B.
!    On exit, the N by NRHS solution matrix X.
!
!    Input, integer ( kind = 4 ) LDB, the leading dimension of the array B. 
!    M < LDB.
!
!    Workspace, real ( kind = 8 ) WORK(LWORK).
!
!    Input, integer ( kind = 4 ) LWORK.
!    The length of the array WORK.  LWORK must be at least NRHS,
!    and should be at least NRHS*NB, where NB is the block size
!    for this environment.
!
!    Output, integer ( kind = 4 ) INFO.
!    = 0: successful exit
!    < 0: if INFO = -i, the i-th argument had an illegal value
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldb
  integer ( kind = 4 ) lwork
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nrhs

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) aii
  real ( kind = 8 ) b(ldb,nrhs)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) tau(n)
  real ( kind = 8 ) temp
  real ( kind = 8 ) work(lwork)
!
!  B := Q' * B.
!
!  Q is an orthogonal matrix of order M, the product of N elementary reflectors.
!
!  B is an M by NRHS vector.
! 
  do i = 1, n
!
!  H(I) is applied to B(I:M,1:NRHS).
!
    if ( tau(i) /= 0.0D+00 ) then

      aii = a(i,i)
      a(i,i) = 1.0D+00
!
!  work := B' * v.
!
      do j = 1, nrhs
        temp = 0.0D+00
        do i2 = i, m
          temp = temp + b(i2,j) * a(i2,i)
        end do
        work(j) = temp
      end do
!
!  B := B - tau * v * work'.
!
      do j = 1, nrhs
        if ( work(j) /= 0.0D+00 )then
          temp = work(j)
          do i2 = i, m
            b(i2,j) = b(i2,j) - tau(i) * a(i2,i) * temp
          end do
        end if
      end do

      a(i,i) = aii

    end if

  end do
!
!  Solve R*X = B(1:n,:).
!
  do j = 1, nrhs
    do k = n, 1, -1
      if ( b(k,j) /= 0.0D+00 )then
        b(k,j) = b(k,j) / a(k,k)
        do i = 1, k - 1
          b(i,j) = b(i,j) - b(k,j) * a(i,k)
        end do
      end if
    end do
  end do

  return
end
subroutine r8gb_print ( m, n, ml, mu, a, title )

!*****************************************************************************80
!
!! R8GB_PRINT prints an R8GB matrix.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1..
!
!    Input, real ( kind = 8 ) A(2*ML+MU+1,N), the R8GB matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  integer ( kind = 4 ) m
  character ( len = * ) title

  call r8gb_print_some ( m, n, ml, mu, a, 1, 1, m, n, title )

  return
end
subroutine r8gb_print_some ( m, n, ml, mu, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8GB_PRINT_SOME prints some of an R8GB matrix.
!
!  Discussion:
!
!    The R8GB storage format is for an M by N banded matrix, with lower 
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
!    extra superdiagonals, which may be required to store nonzero entries 
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.  
!
!    R8GB storage is used by LINPACK and LAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1..
!
!    Input, real ( kind = 8 ) A(2*ML+MU+1,N), the R8GB matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ), parameter :: incx = 5
  integer   ( kind = 4 ) ml
  integer   ( kind = 4 ) mu
  integer   ( kind = 4 ) n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  character ( len = 14 ) ctemp(incx)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2hi
  integer   ( kind = 4 ) i2lo
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j2
  integer   ( kind = 4 ) j2hi
  integer   ( kind = 4 ) j2lo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  integer   ( kind = 4 ) m
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - mu - ml )
    i2hi = min ( ihi, m )
    i2hi = min ( i2hi, j2hi + ml )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i < j - ml - mu  .or. j + ml < i ) then
          ctemp(j2) = '              '
        else
          write ( ctemp(j2), '(g14.6)' ) a(i-j+ml+mu+1,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine r8ge_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8GE_PRINT prints an R8 GE matrix.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * )  title

  call r8ge_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8ge_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8GE_PRINT_SOME prints some of an R8 GE matrix.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ), parameter :: incx = 5
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2hi
  integer   ( kind = 4 ) i2lo
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j2
  integer   ( kind = 4 ) j2hi
  integer   ( kind = 4 ) j2lo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

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

  integer   ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer   ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
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

  character ( len = 8 ) ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
