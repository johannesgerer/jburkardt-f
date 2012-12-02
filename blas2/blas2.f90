subroutine cgbmv ( trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy )

!*****************************************************************************80
!
!! CGBMV computes y := alpha * A * x + beta * y, A a complex band matrix.
!
!  Discussion:
!
!    CGBMV performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
!
!     y := alpha*conjg( A' )*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n band matrix, with kl sub-diagonals and ku super-diagonals.
!
!  Parameters:
!
!  TRANS  - character.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - integer.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  KL     - integer.
!           On entry, KL specifies the number of sub-diagonals of the
!           matrix A. KL must satisfy  0  <=  KL.
!           Unchanged on exit.
!
!  KU     - integer.
!           On entry, KU specifies the number of super-diagonals of the
!           matrix A. KU must satisfy  0  <=  KU.
!           Unchanged on exit.
!
!  ALPHA  - complex         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - complex          array of DIMENSION ( LDA, n ).
!           Before entry, the leading ( kl + ku + 1 ) by n part of the
!           array A must contain the matrix of coefficients, supplied
!           column by column, with the leading diagonal of the matrix in
!           row ( ku + 1 ) of the array, the first super-diagonal
!           starting at position 2 in row ku, the first sub-diagonal
!           starting at position 1 in row ( ku + 2 ), and so on.
!           Elements in the array A that do not correspond to elements
!           in the band matrix (such as the top left ku by ku triangle)
!           are not referenced.
!           The following program segment will transfer a band matrix
!           from conventional full matrix storage to band storage:
!
!                 do j = 1, n
!                    k = kU + 1 - J
!                    do i = max ( 1, j - KU ), min ( M, J + KL )
!                       a( K + I,j) = matrix(i,j)
!                    end do
!                 end do
!
!           Unchanged on exit.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( kl + ku + 1 ).
!           Unchanged on exit.
!
!  X      - complex          array of DIMENSION at least
!           ( 1 + ( n - 1 ) * abs( incx ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 ) * abs( incx ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
!  BETA   - complex         .
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - complex          array of DIMENSION at least
!           ( 1 + ( m - 1 ) * abs( incy ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 ) * abs( incy ) ) otherwise.
!           Before entry, the incremented array Y must contain the
!           vector y. On exit, Y is overwritten by the updated vector y.
!
!
!  incy   - integer.
!           On entry, incy specifies the increment for the elements of
!           Y. incy must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer lda

  complex alpha
  complex beta
  integer incx
  integer incy, kl, ku, m, n
  character        trans
  complex            a( lda, * ), x( * ), y( * )
  complex, parameter :: one = ( 1.0E+00, 0.0E+00 )
  complex, parameter :: zero = ( 0.0E+00, 0.0E+00 )
  complex            temp
  integer            i, info, ix, iy, j, jx, jy, k, kup1, kx, ky, &
                         lenx, leny
  logical            noconj
  logical            lsame
!
!  Test the input.
!
      info = 0
      if   ( .not.lsame ( trans, 'N' ) .and.  &
               .not.lsame ( trans, 'T' ) .and.  &
               .not.lsame ( trans, 'C' )      ) then
         info = 1
      else if ( m<0 ) then
         info = 2
      else if ( n<0 ) then
         info = 3
      else if ( kl<0 ) then
         info = 4
      else if ( ku<0 ) then
         info = 5
      else if ( lda<( kl + ku + 1 ) ) then
         info = 8
      else if ( incx == 0 ) then
         info = 10
      else if ( incy == 0 ) then
         info = 13
      end if

  if ( info /= 0 ) then
    call xerbla ( 'cgbmv ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( ( m == 0 ).or.( n == 0 ) .or.  &
    ( ( alpha == zero ) .and. ( beta == one ) ) ) then
    return
  end if

      noconj = lsame ( trans, 'T' )
!
!  Set LENX and LENY, the lengths of the vectors X and Y, and set
!  up the start points in X and Y.
!
      if ( lsame ( trans, 'N' ) ) then
        lenx = n
        leny = m
      else
        lenx = m
        leny = n
      end if

      if ( 0 < incx ) then
        kx = 1
      else
        kx = 1 - ( lenx - 1 ) * incx
      end if

      if ( 0 < incy ) then
        ky = 1
      else
        ky = 1 - ( leny - 1 ) * incy
      end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through the band part of A.
!
!  First form  y := beta*y.
!
      if ( beta /= one ) then
         if ( incy == 1 ) then
            if ( beta == zero ) then
               y(1:leny) = zero
            else
               y(1:leny) = beta * y(1:leny)
            end if
         else
            iy = ky
            if ( beta == zero ) then
               do i = 1, leny
                  y(iy) = zero
                  iy = iy + incy
               end do
            else
               do i = 1, leny
                  y(iy) = beta * y(iy)
                  iy = iy + incy
               end do
            end if
         end if
      end if

      if ( alpha == zero ) then
         return
      end if

      kup1 = ku + 1
      if ( lsame ( trans, 'N' ) ) then
!
!  Form  y := alpha*A*x + y.
!
         jx = kx
         if ( incy == 1 ) then
            do j = 1, n
               if ( x(jx) /= zero ) then
                  temp = alpha * x(jx)
                  k    = kup1 - j
                  do i = max ( 1, j - ku ), min ( m, j + kl )
                     y(i) = y(i) + temp * a( k + i,j)
                  end do
               end if
               jx = jx + incx
            end do
         else
            do j = 1, n
               if ( x(jx) /= zero ) then
                  temp = alpha * x(jx)
                  iy = ky
                  k    = kup1 - j
                  do i = max ( 1, j - ku ), min ( m, j + kl )
                     y(iy) = y(iy) + temp * a( k + i,j)
                     iy = iy + incy
                  end do
               end if
               jx = jx + incx
               if ( ku < j ) then
                  ky = ky + incy
               end if
            end do
         end if
      else
!
!  Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
!
         jy = ky
         if ( incx == 1 ) then
            do j = 1, n
               temp = zero
               k    = kup1 - j
               if ( noconj ) then
                  do i = max ( 1, j - ku ), min ( m, j + kl )
                     temp = temp + a( k + i,j) * x(i)
                  end do
               else
                  do i = max ( 1, j - ku ), min ( m, j + kl )
                     temp = temp + conjg ( a( k + i,j) ) * x(i)
                  end do
               end if
               y(jy) = y(jy) + alpha * temp
               jy = jy + incy
            end do
         else
            do j = 1, n
               temp = zero
               ix = kx
               k    = kup1 - j
               if ( noconj ) then
                  do i = max ( 1, j - ku ), min ( m, j + kl )
                     temp = temp + a( k + i,j) * x(ix)
                     ix = ix + incx
                  end do
               else
                  do i = max ( 1, j - ku ), min ( m, j + kl )
                     temp = temp + conjg ( a( k + i,j) ) * x(ix)
                     ix = ix + incx
                  end do
               end if
               y(jy) = y(jy) + alpha * temp
               jy = jy + incy
               if ( ku < j ) then
                  kx = kx + incx
               end if
            end do
         end if
      end if

  return
end
subroutine cgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

!*****************************************************************************80
!
!! CGEMV computes y := alpha * A * x + beta * y, A a general complex matrix.
!
!  Discussion:
!
!    CGEMV performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
!
!     y := alpha*conjg( A' )*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Parameters:
!
!  TRANS  - character.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - integer.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - complex         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - complex          array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - complex          array of DIMENSION at least
!           ( 1 + ( n - 1 ) * abs( incx ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 ) * abs( incx ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
!  BETA   - complex         .
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - complex          array of DIMENSION at least
!           ( 1 + ( m - 1 ) * abs( incy ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 ) * abs( incy ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  incy   - integer.
!           On entry, incy specifies the increment for the elements of
!           Y. incy must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer lda

  complex alpha
  complex beta
      integer            incx, incy, m, n
      character        trans
      complex            a( lda, * ), x( * ), y( * )
      complex, parameter :: one = ( 1.0E+00, 0.0E+00 )
      complex, parameter :: zero = ( 0.0E+00, 0.0E+00 )
      complex            temp
      integer            i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
      logical            noconj
      logical            lsame
!
!  Test the input.
!
      info = 0
      if   ( .not.lsame ( trans, 'N' ) .and.  &
               .not.lsame ( trans, 'T' ) .and.  &
               .not.lsame ( trans, 'C' )      ) then
         info = 1
      else if ( m < 0 ) then
         info = 2
      else if ( n < 0 ) then
         info = 3
      else if ( lda < max ( 1, m ) ) then
         info = 6
      else if ( incx == 0 ) then
         info = 8
      else if ( incy == 0 ) then
         info = 11
      end if

      if ( info /= 0 ) then
         call xerbla ( 'cgemv ', info )
         return
      end if
!
!  Quick return if possible.
!
      if ( ( m == 0 ).or.( n == 0 ) .or.  &
          ( ( alpha == zero ) .and. ( beta == one ) ) ) then
         return
      end if

      noconj = lsame ( trans, 'T' )
!
!  Set  lenx  and  leny, the lengths of the vectors x and y, and set
!  up the start points in  X  and  Y.
!
      if ( lsame ( trans, 'N' ) ) then
        lenx = n
        leny = m
      else
        lenx = m
        leny = n
      end if

      if ( 0 < incx ) then
        kx = 1
      else
        kx = 1 - ( lenx - 1 ) * incx
      end if

      if ( 0 < incy ) then
        ky = 1
      else
        ky = 1 - ( leny - 1 ) * incy
      end if
!
!  Start the operations.  In this version the elements of A are
!  accessed sequentially with one pass through A.
!
!  First form  y := beta*y.
!
      if ( beta /= one ) then
         if ( incy == 1 ) then
            if ( beta == zero ) then
               do i = 1, leny
                  y(i) = zero
               end do
            else
               do i = 1, leny
                  y(i) = beta * y(i)
               end do
            end if
         else
            iy = ky
            if ( beta == zero ) then
               do i = 1, leny
                  y(iy) = zero
                  iy = iy + incy
               end do
            else
               do i = 1, leny
                  y(iy) = beta * y(iy)
                  iy = iy + incy
               end do
            end if
         end if
      end if

      if ( alpha == zero ) then
         return
      end if

      if ( lsame ( trans, 'N' ) ) then
!
!  Form  y := alpha*A*x + y.
!
         jx = kx
         if ( incy == 1 ) then
            do j = 1, n
               if ( x(jx) /= zero ) then
                  temp = alpha * x(jx)
                  do i = 1, m
                     y(i) = y(i) + temp * a(i,j)
                  end do
               end if
               jx = jx + incx
            end do
         else
            do j = 1, n
               if ( x(jx) /= zero ) then
                  temp = alpha * x(jx)
                  iy = ky
                  do i = 1, m
                     y(iy) = y(iy) + temp * a(i,j)
                     iy = iy + incy
                  end do
               end if
               jx = jx + incx
            end do
         end if
      else
!
!  Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
!
         jy = ky
         if ( incx == 1 ) then
            do j = 1, n
               temp = zero
               if ( noconj ) then
                  do i = 1, m
                     temp = temp + a(i,j) * x(i)
                  end do
               else
                  do i = 1, m
                     temp = temp + conjg ( a(i,j) ) * x(i)
                  end do
               end if
               y(jy) = y(jy) + alpha * temp
               jy = jy + incy
            end do
         else
            do j = 1, n
               temp = zero
               ix = kx
               if ( noconj ) then
                  do i = 1, m
                     temp = temp + a(i,j) * x(ix)
                     ix = ix + incx
                  end do
               else
                  do i = 1, m
                     temp = temp + conjg ( a(i,j) ) * x(ix)
                     ix = ix + incx
                  end do
               end if
               y(jy) = y(jy) + alpha * temp
               jy = jy + incy
            end do
         end if
      end if

  return
end
subroutine cgerc ( m, n, alpha, x, incx, y, incy, a, lda )

!*****************************************************************************80
!
!! CGERC performs the rank 1 operation A := A + alpha * x * conjg ( y' ).
!
!  Discussion:
!
!    ALPHA is a scalar, x is an m element vector, y is an n element
!    vector and A is an m by n matrix.
!
!  Parameters:
!
!  M      - integer.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - complex         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - complex          array of dimension at least
!           ( 1 + ( m - 1 ) * abs( incx ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
!  Y      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incy ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  incy   - integer.
!           On entry, incy specifies the increment for the elements of
!           Y. incy must not be zero.
!           Unchanged on exit.
!
!  A      - complex          array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
  implicit none

  integer lda

  complex alpha
  integer incx
  integer incy
  integer m
  integer n
  complex temp

  complex            a( lda, * ), x( * ), y( * )
  complex, parameter :: zero = ( 0.0E+00, 0.0E+00 )
  integer i, info, ix, j, jy, kx
!
!  Test the input.
!
      info = 0
      if   ( m<0 ) then
         info = 1
      else if ( n<0 ) then
         info = 2
      else if ( incx == 0 ) then
         info = 5
      else if ( incy == 0 ) then
         info = 7
      else if ( lda<max ( 1, m ) ) then
         info = 9
      end if

      if ( info /= 0 ) then
         call xerbla ( 'cgerc ', info )
         return
      end if
!
!  Quick return if possible.
!
  if ( ( m == 0 ).or.( n == 0 ) .or. ( alpha == zero ) ) then
    return
  end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through A.
!
      if ( 0 < incy ) then
        jy = 1
      else
        jy = 1 - ( n - 1 )*incy
      end if

      if ( incx == 1 ) then
         do j = 1, n
            if ( y(jy) /= zero ) then
               temp = alpha * conjg ( y(jy) )
               do i = 1, m
                  a(i,j) = a(i,j) + x(i) * temp
               end do
            end if
            jy = jy + incy
         end do
      else
         if ( 0 < incx ) then
            kx = 1
         else
            kx = 1 - ( m - 1 ) * incx
         end if
         do j = 1, n
            if ( y(jy) /= zero ) then
               temp = alpha * conjg ( y(jy) )
               ix = kx
               do i = 1, m
                  a(i,j) = a(i,j) + x(ix) * temp
                  ix = ix + incx
               end do
            end if
            jy = jy + incy
         end do
      end if

  return
end
subroutine cgeru ( m, n, alpha, x, incx, y, incy, a, lda )

!*****************************************************************************80
!
!! CGERU performs the rank 1 operation A := A + alpha * x * y'.
!
!  Discussion:
!
!    alpha is a scalar, x is an m element vector, y is an n element
!    vector and A is an m by n matrix.
!
!  Parameters:
!
!  M      - integer.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - complex         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - complex          array of dimension at least
!           ( 1 + ( m - 1 ) * abs( incx ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
!  Y      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incy ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  incy   - integer.
!           On entry, incy specifies the increment for the elements of
!           Y. incy must not be zero.
!           Unchanged on exit.
!
!  A      - complex          array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
  implicit none

  integer lda

  complex alpha
      integer            incx, incy, m, n
      complex            a( lda, * ), x( * ), y( * )
      complex, parameter :: zero = ( 0.0E+00, 0.0E+00 )
      complex            temp
      integer            i, info, ix, j, jy, kx
!
!  Test the input.
!
      info = 0
      if   ( m < 0 ) then
         info = 1
      else if ( n < 0 ) then
         info = 2
      else if ( incx == 0 ) then
         info = 5
      else if ( incy == 0 ) then
         info = 7
      else if ( lda<max ( 1, m ) ) then
         info = 9
      end if

      if ( info /= 0 ) then
         call xerbla ( 'cgeru ', info )
         return
      end if
!
!  Quick return if possible.
!
      if ( ( m == 0 ).or.( n == 0 ) .or. ( alpha == zero ) ) then
         return
      end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through A.
!
      if ( 0 < incy ) then
        jy = 1
      else
        jy = 1 - ( n - 1 ) * incy
      end if

      if ( incx == 1 ) then

         do j = 1, n
            if ( y(jy) /= zero ) then
               temp = alpha * y(jy)
               do i = 1, m
                  a(i,j) = a(i,j) + x(i) * temp
               end do
            end if
            jy = jy + incy
         end do

      else

         if ( 0 < incx ) then
            kx = 1
         else
            kx = 1 - ( m - 1 ) * incx
         end if

         do j = 1, n
            if ( y(jy) /= zero ) then
               temp = alpha * y(jy)
               ix = kx
               do i = 1, m
                  a(i,j) = a(i,j) + x(ix) * temp
                  ix = ix + incx
               end do
            end if
            jy = jy + incy
         end do

      end if

  return
end
subroutine chbmv ( uplo, n, k, alpha, a, lda, x, incx, beta, y, incy )

!*****************************************************************************80
!
!! CHBMV performs the matrix-vector operation y := alpha * A * x + beta * y.
!
!  Discussion:
!
!    ALPHA and BETA are scalars, x and y are n element vectors and
!    A is an n by n hermitian band matrix, with k super-diagonals.
!
!  Parameters:
!
!  UPLO   - character.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the band matrix A is being supplied as
!           follows:
!
!              UPLO = 'U' or 'u'   The upper triangular part of A is
!                                  being supplied.
!
!              UPLO = 'L' or 'l'   The lower triangular part of A is
!                                  being supplied.
!
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  K      - integer.
!           On entry, K specifies the number of super-diagonals of the
!           matrix A. K must satisfy  0  <=  K.
!           Unchanged on exit.
!
!  ALPHA  - complex         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - complex          array of DIMENSION ( LDA, n ).
!           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!           by n part of the array A must contain the upper triangular
!           band part of the hermitian matrix, supplied column by
!           column, with the leading diagonal of the matrix in row
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. The top left k by k triangle
!           of the array A is not referenced.
!           The following program segment will transfer the upper
!           triangular part of a hermitian band matrix from conventional
!           full matrix storage to band storage:
!
!                 do j = 1, n
!                    m = k + 1 - J
!                    do i = max ( 1, j - K ), J
!                       a( M + I,j) = matrix(i,j)
!                    end do
!                 end do
!
!           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
!           by n part of the array A must contain the lower triangular
!           band part of the hermitian matrix, supplied column by
!           column, with the leading diagonal of the matrix in row 1 of
!           the array, the first sub-diagonal starting at position 1 in
!           row 2, and so on. The bottom right k by k triangle of the
!           array A is not referenced.
!           The following program segment will transfer the lower
!           triangular part of a hermitian band matrix from conventional
!           full matrix storage to band storage:
!
!                 do j = 1, n
!                    m = 1 - J
!                    do i = J, min ( N, j + K )
!                       a( M + I,j) = matrix(i,j)
!                    end do
!                 end do
!
!           Note that the imaginary parts of the diagonal elements need
!           not be set and are assumed to be zero.
!           Unchanged on exit.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( k + 1 ).
!           Unchanged on exit.
!
!  X      - complex          array of DIMENSION at least
!           ( 1 + ( n - 1 ) * abs( incx ) ).
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
!  BETA   - complex         .
!           On entry, BETA specifies the scalar beta.
!           Unchanged on exit.
!
!  Y      - complex          array of DIMENSION at least
!           ( 1 + ( n - 1 ) * abs( incy ) ).
!           Before entry, the incremented array Y must contain the
!           vector y. On exit, Y is overwritten by the updated vector y.
!
!  incy   - integer.
!           On entry, incy specifies the increment for the elements of
!           Y. incy must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer lda

  complex a( lda, * )
  complex alpha
  complex beta
  integer i
  integer incx
  integer incy
  integer info
  integer ix
  integer iy
  integer j
  integer jx
  integer jy
  integer k
  integer kplus1
  integer kx
  integer ky
  integer l
  logical lsame
  integer n
  complex, parameter :: one = ( 1.0E+00, 0.0E+00 )
  complex temp1
  complex temp2
  character uplo
  complex x( * )
  complex y( * )
  complex, parameter :: zero = ( 0.0E+00, 0.0E+00 )
!
!  Test the input.
!
  info = 0

  if ( .not. lsame ( uplo, 'U' ) .and. .not. lsame ( uplo, 'L' ) ) then
    info = 1
  else if ( n < 0 ) then
    info = 2
  else if ( k < 0 ) then
    info = 3
  else if ( lda < ( k + 1 ) ) then
    info = 6
  else if ( incx == 0 ) then
    info = 8
  else if ( incy == 0 ) then
    info = 11
  end if

  if ( info /= 0 ) then
    call xerbla ( 'chbmv ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( ( n == 0 ) .or. ( ( alpha == zero ) .and. ( beta == one ) ) ) then
    return
  end if
!
!  Set up the start points in  X  and  Y.
!
      if ( 0 < incx ) then
        kx = 1
      else
        kx = 1 - ( n - 1 ) * incx
      end if

      if ( 0 < incy ) then
        ky = 1
      else
        ky = 1 - ( n - 1 ) * incy
      end if
!
!  Start the operations. In this version the elements of the array A
!  are accessed sequentially with one pass through A.
!
!  First form  y := beta*y.
!
      if ( beta /= one ) then
         if ( incy == 1 ) then
            if ( beta == zero ) then
              y(1:n) = zero
            else
              y(1:n) = beta * y(1:n)
            end if
         else
            iy = ky
            if ( beta == zero ) then
               do i = 1, n
                  y(iy) = zero
                  iy = iy + incy
               end do
            else
               do i = 1, n
                  y(iy) = beta * y(iy)
                  iy = iy + incy
               end do
            end if
         end if
      end if

      if ( alpha == zero ) then
         return
      end if

      if ( lsame ( uplo, 'U' ) ) then
!
!  Form  y  when upper triangle of A is stored.
!
         kplus1 = k + 1
         if ( ( incx == 1 ) .and. ( incy == 1 ) ) then
            do j = 1, n
               temp1 = alpha * x(j)
               temp2 = zero
               l     = kplus1 - j
               do i = max ( 1, j - k ), j - 1
                  y(i) = y(i) + temp1 * a( l + i,j)
                  temp2 = temp2  + conjg ( a( l + i,j) ) * x(i)
               end do
               y(j) = y(j) + temp1 * real( a( kplus1,j) ) &
                               + alpha * temp2
            end do
         else
            jx = kx
            jy = ky
            do j = 1, n
               temp1 = alpha * x(jx)
               temp2 = zero
               ix = kx
               iy = ky
               l     = kplus1 - j
               do i = max ( 1, j - k ), j - 1
                  y(iy) = y(iy) + temp1 * a( l + i,j)
                  temp2   = temp2   + conjg ( a( l + i,j) ) * x(ix)
                  ix = ix + incx
                  iy = iy + incy
               end do
               y(jy) = y(jy) + temp1 * real( a( kplus1,j) ) &
                                 + alpha * temp2
               jx = jx + incx
               jy = jy + incy
               if ( k < j ) then
                  kx = kx + incx
                  ky = ky + incy
               end if
            end do
         end if
      else
!
!  Form  y  when lower triangle of A is stored.
!
         if ( ( incx == 1 ) .and. ( incy == 1 ) ) then
            do j = 1, n
               temp1  = alpha * x(j)
               temp2 = zero
               y(j) = y(j) + temp1 * real( a( 1,j) )
               l      = 1      - j
               do i = j + 1, min( n, j + k )
                  y(i) = y(i) + temp1 * a( l + i,j)
                  temp2 = temp2  + conjg ( a( l + i,j) ) * x(i)
               end do
               y(j) = y(j) + alpha * temp2
            end do
         else
            jx = kx
            jy = ky
            do j = 1, n
               temp1   = alpha * x(jx)
               temp2   = zero
               y(jy) = y(jy) + temp1 * real( a( 1,j) )
               l       = 1       - j
               ix = jx
               iy = jy
               do i = j + 1, min( n, j + k )
                  ix = ix + incx
                  iy = iy + incy
                  y(iy) = y(iy) + temp1 * a( l + i,j)
                  temp2   = temp2   + conjg ( a( l + i,j) ) * x(ix)
               end do
               y(jy) = y(jy) + alpha * temp2
               jx = jx + incx
               jy = jy + incy
            end do
         end if
      end if

  return
end
subroutine chemv ( uplo, n, alpha, a, lda, x, incx, beta, y, incy )

!*****************************************************************************80
!
!! CHEMV performs the matrix-vector operation y := alpha * A * x + beta * y.
!
!  Discussion:
!
!    alpha and beta are scalars, x and y are n element vectors and
!    A is an n by n hermitian matrix.
!
!  Parameters:
!
!  UPLO   - character.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - complex         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - complex          array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the hermitian matrix and the strictly
!           lower triangular part of A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the hermitian matrix and the strictly
!           upper triangular part of A is not referenced.
!           Note that the imaginary parts of the diagonal elements need
!           not be set and are assumed to be zero.
!           Unchanged on exit.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incx ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
!  BETA   - complex         .
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incy ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y. On exit, Y is overwritten by the updated
!           vector y.
!
!  incy   - integer.
!           On entry, incy specifies the increment for the elements of
!           Y. incy must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer lda

  complex a( lda, * )
  complex alpha
  complex beta
  integer incx
  integer incy
  integer n
  character uplo
  complex x( * )
  complex y( * )
      complex, parameter :: one = ( 1.0E+00, 0.0E+00 )
      complex, parameter :: zero = ( 0.0E+00, 0.0E+00 )
      complex            temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, kx, ky
      logical            lsame
!
!  Test the input.
!
      info = 0
      if   ( .not.lsame ( uplo, 'U' ) .and.  &
               .not.lsame ( uplo, 'L' )      ) then
         info = 1
      else if ( n < 0 ) then
         info = 2
      else if ( lda < max ( 1, n ) ) then
         info = 5
      else if ( incx == 0 ) then
         info = 7
      else if ( incy == 0 ) then
         info = 10
      end if
      if ( info/=0 ) then
         call xerbla ( 'chemv ', info )
         return
      end if
!
!  Quick return if possible.
!
      if ( ( n == 0 ).or.( ( alpha == zero ) .and. ( beta == one ) ) ) then
         return
      end if
!
!  Set up the start points in  X  and  Y.
!
      if ( 0 < incx ) then
         kx = 1
      else
         kx = 1 - ( n - 1 ) * incx
      end if

      if ( 0 < incy ) then
         ky = 1
      else
         ky = 1 - ( n - 1 )*incy
      end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through the triangular part
!  of A.
!
!     First form  y := beta*y.
!
      if ( beta/= one ) then
         if ( incy == 1 ) then
            if ( beta == zero ) then
               y(1:n) = zero
            else
               do i = 1, n
                  y(i) = beta * y(i)
               end do
            end if
         else
            iy = ky
            if ( beta == zero ) then
               do i = 1, n
                  y(iy) = zero
                  iy = iy + incy
               end do
            else
               do i = 1, n
                  y(iy) = beta * y(iy)
                  iy = iy + incy
               end do
            end if
         end if
      end if
      if ( alpha == zero ) &
         return
      if ( lsame ( uplo, 'U' ) ) then
!
!  Form  y  when A is stored in upper triangle.
!
         if ( ( incx == 1 ) .and. ( incy == 1 ) ) then
            do  j = 1, n
               temp1 = alpha * x(j)
               temp2 = zero
               do i = 1, j - 1
                  y(i) = y(i) + temp1 * a(i,j)
                  temp2 = temp2  + conjg ( a(i,j) ) * x(i)
               end do
               y(j) = y(j) + temp1 * real( a(j,j) ) + alpha * temp2
            end do
         else
            jx = kx
            jy = ky
            do j = 1, n
               temp1 = alpha * x(jx)
               temp2 = zero
               ix = kx
               iy = ky
               do i = 1, j - 1
                  y(iy) = y(iy) + temp1 * a(i,j)
                  temp2 = temp2 + conjg ( a(i,j) ) * x(ix)
                  ix = ix + incx
                  iy = iy + incy
               end do
               y(jy) = y(jy) + temp1 * real( a(j,j) ) + alpha * temp2
               jx = jx + incx
               jy = jy + incy
            end do
         end if
      else
!
!  Form  y  when A is stored in lower triangle.
!
         if ( ( incx == 1 ) .and. ( incy == 1 ) ) then
            do j = 1, n
               temp1  = alpha * x(j)
               temp2 = zero
               y(j) = y(j) + temp1 * real( a(j,j) )
               do i = j + 1, n
                  y(i) = y(i) + temp1 * a(i,j)
                  temp2 = temp2  + conjg ( a(i,j) ) * x(i)
               end do
               y(j) = y(j) + alpha * temp2
            end do
         else
            jx = kx
            jy = ky
            do j = 1, n
               temp1   = alpha * x(jx)
               temp2   = zero
               y(jy) = y(jy) + temp1 * real( a(j,j) )
               ix = jx
               iy = jy
               do i = j + 1, n
                  ix = ix + incx
                  iy = iy + incy
                  y(iy) = y(iy) + temp1 * a(i,j)
                  temp2   = temp2   + conjg ( a(i,j) ) * x(ix)
               end do
               y(jy) = y(jy) + alpha * temp2
               jx = jx + incx
               jy = jy + incy
            end do
         end if
      end if

  return
end
subroutine cher ( uplo, n, alpha, x, incx, a, lda )

!*****************************************************************************80
!
!! CHER performs the hermitian rank 1 operation A := A + alpha*x*conjg( x' ).
!
!  Discussion:
!
!    ALPHA is a real scalar, x is an n element vector and A is an
!    n by n hermitian matrix.
!
!  Parameters:
!
!  UPLO   - character.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incx ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
!  A      - complex          array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the hermitian matrix and the strictly
!           lower triangular part of A is not referenced. On exit, the
!           upper triangular part of the array A is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the hermitian matrix and the strictly
!           upper triangular part of A is not referenced. On exit, the
!           lower triangular part of the array A is overwritten by the
!           lower triangular part of the updated matrix.
!           Note that the imaginary parts of the diagonal elements need
!           not be set, they are assumed to be zero, and on exit they
!           are set to zero.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
  implicit none

  integer lda

  real alpha
      integer            incx, n
      character        uplo
      complex            a( lda, * ), x( * )
      complex, parameter :: zero = ( 0.0E+00, 0.0E+00 )
      complex            temp
      integer            i, info, ix, j, jx, kx
      logical            lsame
!
!  Test the input.
!
      info = 0
      if   ( .not.lsame ( uplo, 'U' ) .and.  &
               .not.lsame ( uplo, 'L' )      ) then
         info = 1
      else if ( n < 0 ) then
         info = 2
      else if ( incx == 0 ) then
         info = 5
      else if ( lda < max ( 1, n ) ) then
         info = 7
      end if
      if ( info/=0 ) then
         call xerbla ( 'cher  ', info )
         return
      end if
!
!  Quick return if possible.
!
      if ( n == 0 .or. alpha == real ( zero ) ) then
         return
      end if
!
!  Set the start point in X if the increment is not unity.
!
      if ( incx <= 0 ) then
         kx = 1 - ( n - 1 ) * incx
      else if ( incx/=1 ) then
         kx = 1
      end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through the triangular part
!  of A.
!
      if ( lsame ( uplo, 'U' ) ) then
!
!  Form  A  when A is stored in upper triangle.
!
         if ( incx == 1 ) then
            do j = 1, n
               if ( x(j) /= zero ) then
                  temp = alpha * conjg ( x(j) )
                  do i = 1, j - 1
                     a(i,j) = a(i,j) + x(i) * temp
                  end do
                  a(j,j) = real( a( j,j) ) + real( x(j) * temp )
               else
                  a(j,j) = real( a( j,j) )
               end if
            end do
         else
            jx = kx
            do j = 1, n
               if ( x(jx) /= zero ) then
                  temp = alpha * conjg ( x(jx) )
                  ix = kx
                  do i = 1, j - 1
                     a(i,j) = a(i,j) + x(ix) * temp
                     ix = ix + incx
                  end do
                  a(j,j) = real( a( j,j) ) + real ( x(jx) * temp )
               else
                  a(j,j) = real( a( j,j) )
               end if
               jx = jx + incx
            end do
         end if
      else
!
!  Form  A  when A is stored in lower triangle.
!
         if ( incx == 1 ) then
            do j = 1, n
               if ( x(j) /= zero ) then
                  temp = alpha * conjg ( x(j) )
                  a(j,j) = real( a( j,j) ) + real( temp * x(j) )
                  do i = j + 1, n
                     a(i,j) = a(i,j) + x(i) * temp
                  end do
               else
                  a(j,j) = real( a( j,j) )
               end if
            end do
         else
            jx = kx
            do j = 1, n
               if ( x(jx) /= zero ) then
                  temp = alpha * conjg ( x(jx) )
                  a(j,j) = real( a( j,j) ) + real( temp * x(jx) )
                  ix = jx
                  do i = j + 1, n
                     ix = ix + incx
                     a(i,j) = a(i,j) + x(ix) * temp
                  end do
               else
                  a(j,j) = real( a( j,j) )
               end if
               jx = jx + incx
            end do
         end if
      end if

  return
end
subroutine cher2 ( uplo, n, alpha, x, incx, y, incy, a, lda )

!*****************************************************************************80
!
!! CHER2 performs the hermitian rank 2 operation
!
!     A := alpha*x*conjg( y' ) + conjg( alpha )*y*conjg( x' ) + A,
!
!  where alpha is a scalar, x and y are n element vectors and A is an n
!  by n hermitian matrix.
!
!  Parameters:
!
!  UPLO   - character.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - complex         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incx ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
!  Y      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incy ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  incy   - integer.
!           On entry, incy specifies the increment for the elements of
!           Y. incy must not be zero.
!           Unchanged on exit.
!
!  A      - complex          array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the hermitian matrix and the strictly
!           lower triangular part of A is not referenced. On exit, the
!           upper triangular part of the array A is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the hermitian matrix and the strictly
!           upper triangular part of A is not referenced. On exit, the
!           lower triangular part of the array A is overwritten by the
!           lower triangular part of the updated matrix.
!           Note that the imaginary parts of the diagonal elements need
!           not be set, they are assumed to be zero, and on exit they
!           are set to zero.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
  implicit none

  integer lda

  complex alpha
      integer            incx, incy, n
      character        uplo
      complex            a( lda, * ), x( * ), y( * )
      complex, parameter :: zero = ( 0.0E+00, 0.0E+00 )
      complex            temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, kx, ky
      logical            lsame
!
!  Test the input.
!
      info = 0
      if   ( .not.lsame ( uplo, 'U' ) .and.  &
               .not.lsame ( uplo, 'L' )      ) then
         info = 1
      else if ( n<0 ) then
         info = 2
      else if ( incx == 0 ) then
         info = 5
      else if ( incy == 0 ) then
         info = 7
      else if ( lda<max ( 1, n ) ) then
         info = 9
      end if
      if ( info/=0 ) then
         call xerbla ( 'cher2 ', info )
         return
      end if
!
!  Quick return if possible.
!
      if ( ( n == 0 ).or.( alpha == zero ) ) then
         return
      end if
!
!  Set up the start points in X and Y if the increments are not both
!  unity.
!
      if ( ( incx/=1 ).or.( incy/=1 ) ) then
         if ( 0 < incx ) then
            kx = 1
         else
            kx = 1 - ( n - 1 ) * incx
         end if
         if ( 0 < incy ) then
            ky = 1
         else
            ky = 1 - ( n - 1 ) * incy
         end if
         jx = kx
         jy = ky
      end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through the triangular part
!  of A.
!
      if ( lsame ( uplo, 'U' ) ) then
!
!  Form  A  when A is stored in the upper triangle.
!
         if ( ( incx == 1 ) .and. ( incy == 1 ) ) then
            do j = 1, n
               if ( ( x(j) /= zero ).or.( y(j) /= zero ) ) then
                  temp1 = alpha * conjg ( y(j) )
                  temp2 = conjg ( alpha * x(j) )
                  do i = 1, j - 1
                     a(i,j) = a(i,j) + x(i) * temp1 + y(i) * temp2
                  end do
                  a(j,j) = real( a( j,j) ) + &
                              real( x(j) * temp1 + y(j) * temp2 )
               else
                  a(j,j) = real( a( j,j) )
               end if
            end do
         else
            do j = 1, n
               if ( ( x(jx) /= zero ).or.( y(jy) /= zero ) ) then
                  temp1 = alpha * conjg ( y(jy) )
                  temp2 = conjg ( alpha * x(jx) )
                  ix = kx
                  iy = ky
                  do i = 1, j - 1
                     a(i,j) = a(i,j) + x(ix) * temp1 + y(iy) * temp2
                     ix = ix + incx
                     iy = iy + incy
                  end do
                  a(j,j) = real( a( j,j) ) + &
                              real( x(jx) * temp1 + y(jy) * temp2 )
               else
                  a(j,j) = real( a( j,j) )
               end if
               jx = jx + incx
               jy = jy + incy
            end do
         end if
      else
!
!  Form  A  when A is stored in the lower triangle.
!
         if ( ( incx == 1 ) .and. ( incy == 1 ) ) then
            do j = 1, n
               if ( ( x(j) /= zero ).or.( y(j) /= zero ) ) then
                  temp1     = alpha * conjg ( y(j) )
                  temp2     = conjg ( alpha * x(j) )
                  a(j,j) = real( a( j,j) ) + &
                              real( x(j) * temp1 + y(j) * temp2 )
                  do i = j + 1, n
                     a(i,j) = a(i,j) + x(i) * temp1 + y(i) * temp2
                  end do
               else
                  a(j,j) = real( a( j,j) )
               end if
            end do
         else
            do j = 1, n
               if ( ( x(jx) /= zero ).or.( y(jy) /= zero ) ) then
                  temp1     = alpha * conjg ( y(jy) )
                  temp2     = conjg ( alpha * x(jx) )
                  a(j,j) = real( a( j,j) ) + &
                              real( x(jx) * temp1 + y(jy) * temp2 )
                  ix = jx
                  iy = jy
                  do i = j + 1, n
                     ix = ix + incx
                     iy = iy + incy
                     a(i,j) = a(i,j) + x(ix) * temp1 &
                                           + y(iy) * temp2
                  end do
               else
                  a(j,j) = real( a( j,j) )
               end if
               jx = jx + incx
               jy = jy + incy
            end do
         end if
      end if

  return
end
subroutine chpmv ( uplo, n, alpha, ap, x, incx, beta, y, incy )

!*****************************************************************************80
!
!! CHPMV performs the matrix-vector operation y := alpha*A*x + beta*y.
!
!  Discussion:
!
!    ALPHA and BETA are scalars, x and y are n element vectors and
!    A is an n by n hermitian matrix, supplied in packed form.
!
!  Parameters:
!
!  UPLO   - character.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the matrix A is supplied in the packed
!           array AP as follows:
!
!              UPLO = 'U' or 'u'   The upper triangular part of A is
!                                  supplied in AP.
!
!              UPLO = 'L' or 'l'   The lower triangular part of A is
!                                  supplied in AP.
!
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - complex         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  AP     - complex          array of DIMENSION at least
!           ( ( n*( n + 1 ) )/2 ).
!           Before entry with UPLO = 'U' or 'u', the array AP must
!           contain the upper triangular part of the hermitian matrix
!           packed sequentially, column by column, so that ap( 1 )
!           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )
!           and a( 2, 2 ) respectively, and so on.
!           Before entry with UPLO = 'L' or 'l', the array AP must
!           contain the lower triangular part of the hermitian matrix
!           packed sequentially, column by column, so that ap( 1 )
!           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )
!           and a( 3, 1 ) respectively, and so on.
!           Note that the imaginary parts of the diagonal elements need
!           not be set and are assumed to be zero.
!           Unchanged on exit.
!
!  X      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incx ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
!  BETA   - complex         .
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incy ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y. On exit, Y is overwritten by the updated
!           vector y.
!
!  incy   - integer.
!           On entry, incy specifies the increment for the elements of
!           Y. incy must not be zero.
!           Unchanged on exit.
!
!
  implicit none

  complex alpha
  complex beta
      integer            incx, incy, n
      character        uplo
      complex            ap( * ), x( * ), y( * )
      complex, parameter :: one = ( 1.0E+00, 0.0E+00 )
      complex, parameter :: zero = ( 0.0E+00, 0.0E+00 )
      complex            temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, k, kk, kx, ky
      logical            lsame
!
!  Test the input.
!
      info = 0
      if   ( .not.lsame ( uplo, 'U' ) .and.  &
               .not.lsame ( uplo, 'L' )      ) then
         info = 1
      else if ( n<0 ) then
         info = 2
      else if ( incx == 0 ) then
         info = 6
      else if ( incy == 0 ) then
         info = 9
      end if
      if ( info/=0 ) then
         call xerbla ( 'chpmv ', info )
         return
      end if
!
!     Quick return if possible.
!
      if ( ( n == 0 ).or.( ( alpha == zero ) .and. ( beta == one ) ) ) &
         return
!
!  Set up the start points in  X  and  Y.
!
      if ( 0 < incx ) then
         kx = 1
      else
         kx = 1 - ( n - 1 ) * incx
      end if
      if ( 0 < incy ) then
         ky = 1
      else
         ky = 1 - ( n - 1 )*incy
      end if
!
!  Start the operations. In this version the elements of the array AP
!  are accessed sequentially with one pass through AP.
!
!  First form  y := beta*y.
!
      if ( beta/= one ) then
         if ( incy == 1 ) then
            if ( beta == zero ) then
               y(1:n) = zero
            else
               y(1:n) = beta * y(1:n)
            end if
         else
            iy = ky
            if ( beta == zero ) then
               do i = 1, n
                  y(iy) = zero
                  iy = iy + incy
               end do
            else
               do i = 1, n
                  y(iy) = beta * y(iy)
                  iy = iy + incy
               end do
            end if
         end if
      end if
      if ( alpha == zero ) then
         return
      end if
      kk = 1
      if ( lsame ( uplo, 'U' ) ) then
!
!  Form Y when AP contains the upper triangle.
!
         if ( ( incx == 1 ) .and. ( incy == 1 ) ) then
            do j = 1, n
               temp1 = alpha * x(j)
               temp2 = zero
               k     = kk
               do i = 1, j - 1
                  y(i) = y(i) + temp1 * ap(k)
                  temp2 = temp2  + conjg ( ap(k) ) * x(i)
                  k = k + 1
               end do
               y(j) = y(j) + temp1 * real( ap( kk + j - 1 ) ) &
                               + alpha * temp2
               kk = kk + j
            end do
         else
            jx = kx
            jy = ky
            do j = 1, n
               temp1 = alpha * x(jx)
               temp2 = zero
               ix = kx
               iy = ky
               do k = kk, kk + j - 2
                  y(iy) = y(iy) + temp1 * ap(k)
                  temp2   = temp2   + conjg ( ap(k) ) * x(ix)
                  ix = ix + incx
                  iy = iy + incy
               end do
               y(jy) = y(jy) + temp1 * real( ap( kk + j - 1 ) ) &
                                 + alpha * temp2
               jx = jx + incx
               jy = jy + incy
               kk = kk + j
            end do
         end if
      else
!
!  Form Y when AP contains the lower triangle.
!
         if ( ( incx == 1 ) .and. ( incy == 1 ) ) then
            do j = 1, n
               temp1  = alpha * x(j)
               temp2 = zero
               y(j) = y(j) + temp1 * real( ap( kk ) )
               k = kk + 1
               do i = j + 1, n
                  y(i) = y(i) + temp1 * ap(k)
                  temp2 = temp2  + conjg ( ap(k) ) * x(i)
                  k      = k      + 1
               end do
               y(j) = y(j) + alpha * temp2
               kk = kk + ( n - j + 1 )
            end do
         else
            jx = kx
            jy = ky
            do j = 1, n
               temp1   = alpha * x(jx)
               temp2   = zero
               y(jy) = y(jy) + temp1 * real( ap( kk ) )
               ix = jx
               iy = jy
               do k = kk + 1, kk + n - j
                  ix = ix + incx
                  iy = iy + incy
                  y(iy) = y(iy) + temp1 * ap(k)
                  temp2   = temp2   + conjg ( ap(k) ) * x(ix)
               end do
               y(jy) = y(jy) + alpha * temp2
               jx = jx + incx
               jy = jy + incy
               kk = kk + ( n - j + 1 )
            end do
         end if
      end if

  return
end
subroutine chpr ( uplo, n, alpha, x, incx, ap )

!*****************************************************************************80
!
!! CHPR performs the hermitian rank 1 operation A := A + alpha*x*conjg( x' ).
!
!  Discussion:
!
!    ALPHA is a real scalar, x is an n element vector and A is an
!    n by n hermitian matrix, supplied in packed form.
!
!  Parameters:
!
!  UPLO   - character.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the matrix A is supplied in the packed
!           array AP as follows:
!
!              UPLO = 'U' or 'u'   The upper triangular part of A is
!                                  supplied in AP.
!
!              UPLO = 'L' or 'l'   The lower triangular part of A is
!                                  supplied in AP.
!
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incx ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
!  AP     - complex          array of DIMENSION at least
!           ( ( n*( n + 1 ) )/2 ).
!           Before entry with  UPLO = 'U' or 'u', the array AP must
!           contain the upper triangular part of the hermitian matrix
!           packed sequentially, column by column, so that ap( 1 )
!           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )
!           and a( 2, 2 ) respectively, and so on. On exit, the array
!           AP is overwritten by the upper triangular part of the
!           updated matrix.
!           Before entry with UPLO = 'L' or 'l', the array AP must
!           contain the lower triangular part of the hermitian matrix
!           packed sequentially, column by column, so that ap( 1 )
!           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )
!           and a( 3, 1 ) respectively, and so on. On exit, the array
!           AP is overwritten by the lower triangular part of the
!           updated matrix.
!           Note that the imaginary parts of the diagonal elements need
!           not be set, they are assumed to be zero, and on exit they
!           are set to zero.
!
  implicit none

  real alpha
      integer            incx, n
      character        uplo
      complex            ap( * ), x( * )
      complex, parameter :: zero = ( 0.0E+0, 0.0E+0 )
      complex            temp
      integer            i, info, ix, j, jx, k, kk, kx
      logical            lsame
!
!  Test the input.
!
      info = 0
      if   ( .not.lsame ( uplo, 'U' ) .and.  &
               .not.lsame ( uplo, 'L' )      ) then
         info = 1
      else if ( n<0 ) then
         info = 2
      else if ( incx == 0 ) then
         info = 5
      end if
      if ( info/=0 ) then
         call xerbla ( 'chpr  ', info )
         return
      end if
!
!  Quick return if possible.
!
      if ( ( n == 0 ).or.( alpha == real( zero ) ) ) &
         return
!
!  Set the start point in X if the increment is not unity.
!
      if ( incx <= 0 ) then
         kx = 1 - ( n - 1 ) * incx
      else if ( incx/=1 ) then
         kx = 1
      end if
!
!  Start the operations. In this version the elements of the array AP
!  are accessed sequentially with one pass through AP.
!
      kk = 1
      if ( lsame ( uplo, 'U' ) ) then
!
!  Form  A  when upper triangle is stored in AP.
!
         if ( incx == 1 ) then
            do j = 1, n
               if ( x(j) /= zero ) then
                  temp = alpha * conjg ( x(j) )
                  k    = kk
                  do i = 1, j - 1
                     ap(k) = ap(k) + x(i) * temp
                     k       = k       + 1
                  end do
                  ap( kk + j - 1 ) = real( ap( kk + j - 1 ) ) &
                                     + real( x(j) * temp )
               else
                  ap( kk + j - 1 ) = real( ap( kk + j - 1 ) )
               end if
               kk = kk + j
            end do
         else
            jx = kx
            do j = 1, n
               if ( x(jx) /= zero ) then
                  temp = alpha * conjg ( x(jx) )
                  ix = kx
                  do k = kk, kk + j - 2
                     ap(k) = ap(k) + x(ix) * temp
                     ix = ix + incx
                  end do
                  ap( kk + j - 1 ) = real( ap( kk + j - 1 ) ) &
                                     + real( x(jx) * temp )
               else
                  ap( kk + j - 1 ) = real( ap( kk + j - 1 ) )
               end if
               jx = jx + incx
               kk = kk + j
            end do
         end if
      else
!
!  Form A when lower triangle is stored in AP.
!
         if ( incx == 1 ) then
            do j = 1, n
               if ( x(j) /= zero ) then
                  temp = alpha * conjg ( x(j) )
                  ap( kk ) = real( ap( kk ) ) + real( temp * x(j) )
                  k        = kk               + 1
                  do i = j + 1, n
                     ap(k) = ap(k) + x(i) * temp
                     k       = k       + 1
                  end do
               else
                  ap( kk ) = real( ap( kk ) )
               end if
               kk = kk + n - j + 1
            end do
         else
            jx = kx
            do j = 1, n
               if ( x(jx) /= zero ) then
                  temp = alpha * conjg ( x(jx) )
                  ap( kk ) = real( ap( kk ) ) + real( temp * x(jx) )
                  ix = jx
                  do k = kk + 1, kk + n - j
                     ix = ix + incx
                     ap(k) = ap(k) + x(ix) * temp
                  end do
               else
                  ap( kk ) = real( ap( kk ) )
               end if
               jx = jx + incx
               kk = kk + n - j + 1
            end do
         end if
      end if

  return
end
subroutine chpr2 ( uplo, n, alpha, x, incx, y, incy, ap )

!*****************************************************************************80
!
!! CHPR2 performs the hermitian rank 2 operation
!
!     A := alpha*x*conjg( y' ) + conjg( alpha )*y*conjg( x' ) + A,
!
!  where alpha is a scalar, x and y are n element vectors and A is an
!  n by n hermitian matrix, supplied in packed form.
!
!  Parameters:
!
!  UPLO   - character.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the matrix A is supplied in the packed
!           array AP as follows:
!
!              UPLO = 'U' or 'u'   The upper triangular part of A is
!                                  supplied in AP.
!
!              UPLO = 'L' or 'l'   The lower triangular part of A is
!                                  supplied in AP.
!
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - complex         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incx ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
!  Y      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incy ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  incy   - integer.
!           On entry, incy specifies the increment for the elements of
!           Y. incy must not be zero.
!           Unchanged on exit.
!
!  AP     - complex          array of DIMENSION at least
!           ( ( n*( n + 1 ) )/2 ).
!           Before entry with  UPLO = 'U' or 'u', the array AP must
!           contain the upper triangular part of the hermitian matrix
!           packed sequentially, column by column, so that ap( 1 )
!           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )
!           and a( 2, 2 ) respectively, and so on. On exit, the array
!           AP is overwritten by the upper triangular part of the
!           updated matrix.
!           Before entry with UPLO = 'L' or 'l', the array AP must
!           contain the lower triangular part of the hermitian matrix
!           packed sequentially, column by column, so that ap( 1 )
!           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )
!           and a( 3, 1 ) respectively, and so on. On exit, the array
!           AP is overwritten by the lower triangular part of the
!           updated matrix.
!           Note that the imaginary parts of the diagonal elements need
!           not be set, they are assumed to be zero, and on exit they
!           are set to zero.
!
  implicit none

  complex alpha
      integer            incx, incy, n
      character        uplo
      complex            ap( * ), x( * ), y( * )
      complex, parameter :: zero = ( 0.0E+00, 0.0E+00 )
      complex            temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, k, kk, kx, ky
      logical            lsame
!
!  Test the input.
!
      info = 0
      if   ( .not.lsame ( uplo, 'U' ) .and.  &
               .not.lsame ( uplo, 'L' )      ) then
         info = 1
      else if ( n<0 ) then
         info = 2
      else if ( incx == 0 ) then
         info = 5
      else if ( incy == 0 ) then
         info = 7
      end if
      if ( info/=0 ) then
         call xerbla ( 'chpr2 ', info )
         return
      end if
!
!  Quick return if possible.
!
      if ( ( n == 0 ).or.( alpha == zero ) ) &
         return
!
!  Set up the start points in X and Y if the increments are not both unity.
!
      if ( ( incx/=1 ).or.( incy/=1 ) ) then
         if ( 0 < incx ) then
            kx = 1
         else
            kx = 1 - ( n - 1 ) * incx
         end if
         if ( 0 < incy ) then
            ky = 1
         else
            ky = 1 - ( n - 1 )*incy
         end if
         jx = kx
         jy = ky
      end if
!
!  Start the operations. In this version the elements of the array AP
!  are accessed sequentially with one pass through AP.
!
      kk = 1
      if ( lsame ( uplo, 'U' ) ) then
!
!  Form A when upper triangle is stored in AP.
!
         if ( ( incx == 1 ) .and. ( incy == 1 ) ) then
            do j = 1, n
               if ( ( x(j) /= zero ).or.( y(j) /= zero ) ) then
                  temp1 = alpha * conjg ( y(j) )
                  temp2 = conjg ( alpha * x(j) )
                  k     = kk
                  do i = 1, j - 1
                     ap(k) = ap(k) + x(i) * temp1 + y(i) * temp2
                     k       = k       + 1
                  end do
                  ap( kk + j - 1 ) = real( ap( kk + j - 1 ) ) + &
                                     real( x(j) * temp1 + y(j) * temp2 )
               else
                  ap( kk + j - 1 ) = real( ap( kk + j - 1 ) )
               end if
               kk = kk + j
            end do
         else
            do j = 1, n
               if ( ( x(jx) /= zero ).or.( y(jy) /= zero ) ) then
                  temp1 = alpha * conjg ( y(jy) )
                  temp2 = conjg ( alpha * x(jx) )
                  ix = kx
                  iy = ky
                  do k = kk, kk + j - 2
                     ap(k) = ap(k) + x(ix) * temp1 + y(iy) * temp2
                     ix = ix + incx
                     iy = iy + incy
                  end do
                  ap( kk + j - 1 ) = real( ap( kk + j - 1 ) ) + &
                                     real( x(jx) * temp1 + &
                                           y(jy) * temp2 )
               else
                  ap( kk + j - 1 ) = real( ap( kk + j - 1 ) )
               end if
               jx = jx + incx
               jy = jy + incy
               kk = kk + j
            end do
         end if
      else
!
!  Form A when lower triangle is stored in AP.
!
         if ( ( incx == 1 ) .and. ( incy == 1 ) ) then
            do j = 1, n
               if ( ( x(j) /= zero ).or.( y(j) /= zero ) ) then
                  temp1   = alpha * conjg ( y(j) )
                  temp2   = conjg ( alpha * x(j) )
                  ap( kk ) = real( ap( kk ) ) + &
                             real( x(j) * temp1 + y(j) * temp2 )
                  k        = kk               + 1
                  do i = j + 1, n
                     ap(k) = ap(k) + x(i) * temp1 + y(i) * temp2
                     k       = k       + 1
                  end do
               else
                  ap( kk ) = real( ap( kk ) )
               end if
               kk = kk + n - j + 1
            end do
         else
            do j = 1, n
               if ( ( x(jx) /= zero ).or.( y(jy) /= zero ) ) then
                  temp1    = alpha * conjg ( y(jy) )
                  temp2    = conjg ( alpha * x(jx) )
                  ap( kk ) = real( ap( kk ) ) + &
                             real( x(jx) * temp1 + y(jy) * temp2 )
                  ix = jx
                  iy = jy
                  do k = kk + 1, kk + n - j
                     ix = ix + incx
                     iy = iy + incy
                     ap(k) = ap(k) + x(ix) * temp1 + y(iy) * temp2
                  end do
               else
                  ap( kk ) = real( ap( kk ) )
               end if
               jx = jx + incx
               jy = jy + incy
               kk = kk + n - j + 1
            end do
         end if
      end if

  return
end
subroutine ctbmv ( uplo, trans, diag, n, k, a, lda, x, incx )

!*****************************************************************************80
!
!! CTBMV performs one of the matrix-vector operations
!
!     x := A*x,   or   x := A'*x,   or   x := conjg( A' )*x,
!
!  where x is an n element vector and  A is an n by n unit, or non-unit,
!  upper or lower triangular band matrix, with ( k + 1 ) diagonals.
!
!  Parameters:
!
!  UPLO   - character.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - character.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   x := A*x.
!
!              TRANS = 'T' or 't'   x := A'*x.
!
!              TRANS = 'C' or 'c'   x := conjg( A' )*x.
!
!           Unchanged on exit.
!
!  DIAG   - character.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  K      - integer.
!           On entry with UPLO = 'U' or 'u', K specifies the number of
!           super-diagonals of the matrix A.
!           On entry with UPLO = 'L' or 'l', K specifies the number of
!           sub-diagonals of the matrix A.
!           K must satisfy  0  <=  K.
!           Unchanged on exit.
!
!  A      - complex          array of DIMENSION ( LDA, n ).
!           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!           by n part of the array A must contain the upper triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. The top left k by k triangle
!           of the array A is not referenced.
!           The following program segment will transfer an upper
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 do j = 1, n
!                    m = k + 1 - J
!                    do i = max ( 1, j - K ), J
!                       a( M + I,j) = matrix(i,j)
!                    end do
!                 end do
!
!           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
!           by n part of the array A must contain the lower triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row 1 of
!           the array, the first sub-diagonal starting at position 1 in
!           row 2, and so on. The bottom right k by k triangle of the
!           array A is not referenced.
!           The following program segment will transfer a lower
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 do j = 1, n
!                    m = 1 - J
!                    do i = J, min ( N, j + K )
!                       a( M + I,j) = matrix(i,j)
!                    end do
!                 end do
!
!           Note that when DIAG = 'U' or 'u' the elements of the array A
!           corresponding to the diagonal elements of the matrix are not
!           referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( k + 1 ).
!           Unchanged on exit.
!
!  X      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incx ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer lda

  integer            incx, k, n
  character        diag, trans, uplo
  complex            a( lda, * ), x( * )
  complex, parameter :: zero = ( 0.0E+00, 0.0E+00 )
  complex            temp
  integer            i, info, ix, j, jx, kplus1, kx, l
  logical            noconj, nounit
  logical            lsame
!
!  Test the input.
!
      info = 0
      if   ( .not.lsame ( uplo , 'U' ) .and.  &
               .not.lsame ( uplo , 'L' )      ) then
         info = 1
      else if ( .not.lsame ( trans, 'N' ) .and.  &
               .not.lsame ( trans, 'T' ) .and.  &
               .not.lsame ( trans, 'C' )      ) then
         info = 2
      else if ( .not.lsame ( diag , 'U' ) .and.  &
               .not.lsame ( diag , 'N' )      ) then
         info = 3
      else if ( n<0 ) then
         info = 4
      else if ( k<0 ) then
         info = 5
      else if ( lda<( k + 1 ) ) then
         info = 7
      else if ( incx == 0 ) then
         info = 9
      end if
      if ( info/=0 ) then
         call xerbla ( 'ctbmv ', info )
         return
      end if
!
!  Quick return if possible.
!
      if ( n == 0 ) then
         return
      end if

      noconj = lsame ( trans, 'T' )
      nounit = lsame ( diag , 'N' )
!
!  Set up the start point in X if the increment is not unity. This
!  will be  ( N - 1 ) * incx   too small for descending loops.
!
      if ( incx <= 0 ) then
         kx = 1 - ( n - 1 ) * incx
      else if ( incx/=1 ) then
         kx = 1
      end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through A.
!
      if ( lsame ( trans, 'N' ) ) then
!
!  Form x := A*x.
!
         if ( lsame ( uplo, 'U' ) ) then
            kplus1 = k + 1
            if ( incx == 1 ) then
               do j = 1, n
                  if ( x(j) /= zero ) then
                     temp = x(j)
                     l    = kplus1 - j
                     do i = max ( 1, j - k ), j - 1
                        x(i) = x(i) + temp * a( l + i,j)
                     end do
                     if ( nounit ) then
                        x(j) = x(j) * a( kplus1,j)
                     end if
                  end if
               end do
            else
               jx = kx
               do j = 1, n
                  if ( x(jx) /= zero ) then
                     temp = x(jx)
                     ix = kx
                     l    = kplus1  - j
                     do i = max ( 1, j - k ), j - 1
                        x(ix) = x(ix) + temp * a( l + i,j)
                        ix = ix + incx
                     end do
                     if ( nounit ) then
                        x(jx) = x(jx) * a( kplus1,j)
                     end if
                  end if
                  jx = jx + incx
                  if ( k < j ) then
                     kx = kx + incx
                  end if
               end do
            end if
         else
            if ( incx == 1 ) then
               do j = n, 1, -1
                  if ( x(j) /= zero ) then
                     temp = x(j)
                     l    = 1      - j
                     do i = min ( n, j + k ), j + 1, -1
                        x(i) = x(i) + temp * a( l + i,j)
                     end do
                     if ( nounit ) then
                        x(j) = x(j) * a( 1,j)
                     end if
                  end if
               end do
            else
               kx = kx + ( n - 1 ) * incx
               jx = kx
               do j = n, 1, -1
                  if ( x(jx) /= zero ) then
                     temp = x(jx)
                     ix = kx
                     l    = 1       - j
                     do i = min ( n, j + k ), j + 1, -1
                        x(ix) = x(ix) + temp * a( l + i,j)
                        ix = ix - incx
                     end do
                     if ( nounit ) &
                        x(jx) = x(jx) * a( 1,j)
                  end if
                  jx = jx - incx
                  if ( k <= ( n - j ) ) then
                     kx = kx - incx
                  end if
               end do
            end if
         end if
      else
!
!  Form  x := A'*x  or  x := conjg( A' )*x.
!
         if ( lsame ( uplo, 'U' ) ) then
            kplus1 = k + 1
            if ( incx == 1 ) then
               do j = n, 1, -1
                  temp = x(j)
                  l    = kplus1 - j
                  if ( noconj ) then
                     if ( nounit ) then
                        temp = temp * a( kplus1,j)
                     end if
                     do i = j - 1, max( 1, j - k ), -1
                        temp = temp + a( l + i,j) * x(i)
                     end do
                  else
                     if ( nounit ) then
                        temp = temp * conjg ( a( kplus1,j) )
                     end if
                     do i = j - 1, max( 1, j - k ), -1
                        temp = temp + conjg ( a( l + i,j) ) * x(i)
                     end do
                  end if
                  x(j) = temp
               end do
            else
               kx = kx + ( n - 1 ) * incx
               jx = kx
               do j = n, 1, -1
                  temp = x(jx)
                  kx   = kx      - incx
                  ix = kx
                  l    = kplus1  - j
                  if ( noconj ) then
                     if ( nounit ) then
                        temp = temp * a( kplus1,j)
                     end if
                     do i = j - 1, max( 1, j - k ), -1
                        temp = temp + a( l + i,j) * x(ix)
                        ix = ix - incx
                     end do
                  else
                     if ( nounit ) then
                        temp = temp*conjg ( a( kplus1,j) )
                     end if
                     do i = j - 1, max( 1, j - k ), -1
                        temp = temp + conjg ( a( l + i,j) ) * x(ix)
                        ix = ix - incx
                     end do
                  end if
                  x(jx) = temp
                  jx = jx - incx
               end do
            end if
         else
            if ( incx == 1 ) then
               do j = 1, n
                  temp = x(j)
                  l    = 1      - j
                  if ( noconj ) then
                     if ( nounit ) then
                        temp = temp * a( 1,j)
                     end if
                     do i = j + 1, min( n, j + k )
                        temp = temp + a( l + i,j) * x(i)
                     end do
                  else
                     if ( nounit ) then
                        temp = temp * conjg ( a( 1,j) )
                     end if
                     do i = j + 1, min( n, j + k )
                        temp = temp + conjg ( a( l + i,j) ) * x(i)
                     end do
                  end if
                  x(j) = temp
               end do
            else
               jx = kx
               do j = 1, n
                  temp = x(jx)
                  kx   = kx      + incx
                  ix = kx
                  l    = 1       - j
                  if ( noconj ) then
                     if ( nounit ) &
                        temp = temp * a( 1,j)
                     do i = j + 1, min( n, j + k )
                        temp = temp + a( l + i,j) * x(ix)
                        ix = ix + incx
                     end do
                  else
                     if ( nounit ) &
                        temp = temp*conjg ( a( 1,j) )
                     do i = j + 1, min( n, j + k )
                        temp = temp + conjg ( a( l + i,j) ) * x(ix)
                        ix = ix + incx
                     end do
                  end if
                  x(jx) = temp
                  jx = jx + incx
               end do
            end if
         end if
      end if

  return
end
subroutine ctbsv ( uplo, trans, diag, n, k, a, lda, x, incx )

!*****************************************************************************80
!
!! CTBSV solves one of the systems of equations
!
!     A*x = b,   or   A'*x = b,   or   conjg( A' )*x = b,
!
!  where b and x are n element vectors and A is an n by n unit, or
!  non-unit, upper or lower triangular band matrix, with ( k + 1 )
!  diagonals.
!
!  No test for singularity or near-singularity is included in this
!  routine. Such tests must be performed before calling this routine.
!
!  Parameters:
!
!  UPLO   - character.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - character.
!           On entry, TRANS specifies the equations to be solved as
!           follows:
!
!              TRANS = 'N' or 'n'   A*x = b.
!
!              TRANS = 'T' or 't'   A'*x = b.
!
!              TRANS = 'C' or 'c'   conjg( A' )*x = b.
!
!           Unchanged on exit.
!
!  DIAG   - character.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  K      - integer.
!           On entry with UPLO = 'U' or 'u', K specifies the number of
!           super-diagonals of the matrix A.
!           On entry with UPLO = 'L' or 'l', K specifies the number of
!           sub-diagonals of the matrix A.
!           K must satisfy  0  <=  K.
!           Unchanged on exit.
!
!  A      - complex          array of DIMENSION ( LDA, n ).
!           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!           by n part of the array A must contain the upper triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. The top left k by k triangle
!           of the array A is not referenced.
!           The following program segment will transfer an upper
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 do j = 1, n
!                    m = k + 1 - J
!                    do i = max ( 1, j - K ), J
!                       a( M + I,j) = matrix(i,j)
!                    end do
!                 end do
!
!           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
!           by n part of the array A must contain the lower triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row 1 of
!           the array, the first sub-diagonal starting at position 1 in
!           row 2, and so on. The bottom right k by k triangle of the
!           array A is not referenced.
!           The following program segment will transfer a lower
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 do j = 1, n
!                    m = 1 - J
!                    do i = J, min ( N, j + K )
!                       a( M + I,j) = matrix(i,j)
!                    end do
!                 end do
!
!           Note that when DIAG = 'U' or 'u' the elements of the array A
!           corresponding to the diagonal elements of the matrix are not
!           referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( k + 1 ).
!           Unchanged on exit.
!
!  X      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incx ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer lda

      integer            incx, k, n
      character        diag, trans, uplo
      complex            a( lda, * ), x( * )
      complex, parameter :: zero = ( 0.0E+00, 0.0E+00 )
      complex            temp
      integer            i, info, ix, j, jx, kplus1, kx, l
      logical            noconj, nounit
      logical            lsame
!
!  Test the input.
!
      info = 0
      if   ( .not.lsame ( uplo , 'U' ) .and.  &
               .not.lsame ( uplo , 'L' )      ) then
         info = 1
      else if ( .not.lsame ( trans, 'N' ) .and.  &
               .not.lsame ( trans, 'T' ) .and.  &
               .not.lsame ( trans, 'C' )      ) then
         info = 2
      else if ( .not.lsame ( diag , 'U' ) .and.  &
               .not.lsame ( diag , 'N' )      ) then
         info = 3
      else if ( n<0 ) then
         info = 4
      else if ( k<0 ) then
         info = 5
      else if ( lda<( k + 1 ) ) then
         info = 7
      else if ( incx == 0 ) then
         info = 9
      end if
      if ( info/=0 ) then
         call xerbla ( 'ctbsv ', info )
         return
      end if
!
!  Quick return if possible.
!
      if ( n == 0 ) then
         return
      end if

      noconj = lsame ( trans, 'T' )
      nounit = lsame ( diag , 'N' )
!
!  Set up the start point in X if the increment is not unity. This
!  will be  ( N - 1 ) * incx  too small for descending loops.
!
      if ( incx <= 0 ) then
         kx = 1 - ( n - 1 ) * incx
      else if ( incx/=1 ) then
         kx = 1
      end if
!
!  Start the operations. In this version the elements of A are
!  accessed by sequentially with one pass through A.
!
      if ( lsame ( trans, 'N' ) ) then
!
!  Form  x := inv( A )*x.
!
         if ( lsame ( uplo, 'U' ) ) then
            kplus1 = k + 1
            if ( incx == 1 ) then
               do j = n, 1, -1
                  if ( x(j) /= zero ) then
                     l = kplus1 - j
                     if ( nounit ) then
                        x(j) = x(j)/a( kplus1,j)
                     end if
                     temp = x(j)
                     do i = j - 1, max ( 1, j - k ), -1
                        x(i) = x(i) - temp * a( l + i,j)
                     end do
                  end if
               end do
            else
               kx = kx + ( n - 1 ) * incx
               jx = kx
               do j = n, 1, -1
                  kx = kx - incx
                  if ( x(jx) /= zero ) then
                     ix = kx
                     l  = kplus1 - j
                     if ( nounit ) &
                        x(jx) = x(jx) / a( kplus1,j)
                     temp = x(jx)
                     do i = j - 1, max ( 1, j - k ), -1
                        x(ix) = x(ix) - temp * a( l + i,j)
                        ix = ix - incx
                     end do
                  end if
                  jx = jx - incx
               end do
            end if
         else
            if ( incx == 1 ) then
               do j = 1, n
                  if ( x(j) /= zero ) then
                     l = 1 - j
                     if ( nounit ) then
                        x(j) = x(j) / a( 1,j)
                     end if
                     temp = x(j)
                     do i = j + 1, min ( n, j + k )
                        x(i) = x(i) - temp * a( l + i,j)
                     end do
                  end if
               end do
            else
               jx = kx
               do j = 1, n
                  kx = kx + incx
                  if ( x(jx) /= zero ) then
                     ix = kx
                     l  = 1  - j
                     if ( nounit ) then
                        x(jx) = x(jx)/a( 1,j)
                     end if
                     temp = x(jx)
                     do i = j + 1, min ( n, j + k )
                        x(ix) = x(ix) - temp * a( l + i,j)
                        ix = ix + incx
                     end do
                  end if
                  jx = jx + incx
               end do
            end if
         end if
      else
!
!  Form  x := inv( A' )*x  or  x := inv( conjg( A') )*x.
!
         if ( lsame ( uplo, 'U' ) ) then
            kplus1 = k + 1
            if ( incx == 1 ) then
               do j = 1, n
                  temp = x(j)
                  l    = kplus1 - j
                  if ( noconj ) then
                     do i = max ( 1, j - k ), j - 1
                        temp = temp - a( l + i,j) * x(i)
                     end do
                     if ( nounit ) then
                        temp = temp/a( kplus1,j)
                     end if
                  else
                     do i = max ( 1, j - k ), j - 1
                        temp = temp - conjg ( a( l + i,j) ) * x(i)
                     end do
                     if ( nounit ) &
                        temp = temp/conjg ( a( kplus1,j) )
                  end if
                  x(j) = temp
               end do
            else
               jx = kx
               do j = 1, n
                  temp = x(jx)
                  ix = kx
                  l    = kplus1  - j
                  if ( noconj ) then
                     do  i = max ( 1, j - k ), j - 1
                        temp = temp - a( l + i,j) * x(ix)
                        ix = ix + incx
                     end do
                     if ( nounit ) &
                        temp = temp/a( kplus1,j)
                  else
                     do i = max ( 1, j - k ), j - 1
                        temp = temp - conjg ( a( l + i,j) ) * x(ix)
                        ix = ix + incx
                     end do
                     if ( nounit ) then
                        temp = temp/conjg ( a( kplus1,j) )
                     end if
                  end if
                  x(jx) = temp
                  jx = jx + incx
                  if ( l < j ) then
                     kx = kx + incx
                  end if
               end do
            end if
         else
            if ( incx == 1 ) then
               do j = n, 1, -1
                  temp = x(j)
                  l    = 1      - j
                  if ( noconj ) then
                     do i = min ( n, j + k ), j + 1, -1
                        temp = temp - a( l + i,j) * x(i)
                     end do
                     if ( nounit ) &
                        temp = temp/a( 1,j)
                  else
                     do i = min ( n, j + k ), j + 1, -1
                        temp = temp - conjg ( a( l + i,j) ) * x(i)
                     end do
                     if ( nounit ) &
                        temp = temp/conjg ( a( 1,j) )
                  end if
                  x(j) = temp
               end do
            else
               kx = kx + ( n - 1 ) * incx
               jx = kx
               do j = n, 1, -1
                  temp = x(jx)
                  ix = kx
                  l = 1 - j
                  if ( noconj ) then
                     do i = min ( n, j + k ), j + 1, -1
                        temp = temp - a( l + i,j) * x(ix)
                        ix = ix - incx
                     end do
                     if ( nounit ) &
                        temp = temp/a( 1,j)
                  else
                     do i = min ( n, j + k ), j + 1, -1
                        temp = temp - conjg ( a( l + i,j) ) * x(ix)
                        ix = ix - incx
                     end do
                     if ( nounit ) then
                        temp = temp/conjg ( a( 1,j) )
                     end if
                  end if
                  x(jx) = temp
                  jx = jx - incx
                  if ( k <= ( n - j ) ) then
                     kx = kx - incx
                  end if
               end do
            end if
         end if
      end if

  return
end
subroutine ctpmv ( uplo, trans, diag, n, ap, x, incx )

!*****************************************************************************80
!
!! CTPMV performs one of the matrix-vector operations
!
!     x := A*x,   or   x := A'*x,   or   x := conjg( A' )*x,
!
!  where x is an n element vector and  A is an n by n unit, or non-unit,
!  upper or lower triangular matrix, supplied in packed form.
!
!  Parameters:
!
!  UPLO   - character.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - character.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   x := A*x.
!
!              TRANS = 'T' or 't'   x := A'*x.
!
!              TRANS = 'C' or 'c'   x := conjg( A' )*x.
!
!           Unchanged on exit.
!
!  DIAG   - character.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  AP     - complex          array of DIMENSION at least
!           ( ( n*( n + 1 ) )/2 ).
!           Before entry with  UPLO = 'U' or 'u', the array AP must
!           contain the upper triangular matrix packed sequentially,
!           column by column, so that ap( 1 ) contains a( 1, 1 ),
!           ap( 2 ) and ap( 3 ) contain a( 1, 2 ) and a( 2, 2 )
!           respectively, and so on.
!           Before entry with UPLO = 'L' or 'l', the array AP must
!           contain the lower triangular matrix packed sequentially,
!           column by column, so that ap( 1 ) contains a( 1, 1 ),
!           ap( 2 ) and ap( 3 ) contain a( 2, 1 ) and a( 3, 1 )
!           respectively, and so on.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  X      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incx ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer incx
  integer n
      character        diag, trans, uplo
      complex            ap( * ), x( * )
      complex, parameter :: zero = ( 0.0E+00, 0.0E+00 )
      complex            temp
      integer            i, info, ix, j, jx, k, kk, kx
      logical            noconj, nounit
      logical            lsame
!
!  Test the input.
!
      info = 0
      if   ( .not.lsame ( uplo , 'U' ) .and.  &
               .not.lsame ( uplo , 'L' )      ) then
         info = 1
      else if ( .not.lsame ( trans, 'N' ) .and.  &
               .not.lsame ( trans, 'T' ) .and.  &
               .not.lsame ( trans, 'C' )      ) then
         info = 2
      else if ( .not.lsame ( diag , 'U' ) .and.  &
               .not.lsame ( diag , 'N' )      ) then
         info = 3
      else if ( n<0 ) then
         info = 4
      else if ( incx == 0 ) then
         info = 7
      end if
      if ( info /= 0 ) then
         call xerbla ( 'ctpmv ', info )
         return
      end if
!
!  Quick return if possible.
!
      if ( n == 0 ) then
         return
      end if

      noconj = lsame ( trans, 'T' )
      nounit = lsame ( diag , 'N' )
!
!  Set up the start point in X if the increment is not unity. This
!  will be  ( N - 1 ) * incx  too small for descending loops.
!
      if ( incx <= 0 ) then
         kx = 1 - ( n - 1 ) * incx
      else if ( incx /= 1 ) then
         kx = 1
      end if
!
!  Start the operations. In this version the elements of AP are
!  accessed sequentially with one pass through AP.
!
      if ( lsame ( trans, 'N' ) ) then
!
!  Form  x:= A*x.
!
         if ( lsame ( uplo, 'U' ) ) then
            kk = 1
            if ( incx == 1 ) then
               do j = 1, n
                  if ( x(j) /= zero ) then
                     temp = x(j)
                     k    = kk
                     do i = 1, j - 1
                        x(i) = x(i) + temp * ap(k)
                        k      = k      + 1
                     end do
                     if ( nounit ) then
                        x(j) = x(j) * ap( kk + j - 1 )
                     end if
                  end if
                  kk = kk + j
               end do
            else
               jx = kx
               do j = 1, n
                  if ( x(jx) /= zero ) then
                     temp = x(jx)
                     ix = kx
                     do k = kk, kk + j - 2
                        x(ix) = x(ix) + temp * ap(k)
                        ix = ix + incx
                     end do
                     if ( nounit ) then
                        x(jx) = x(jx) * ap( kk + j - 1 )
                     end if
                  end if
                  jx = jx + incx
                  kk = kk + j
               end do
            end if
         else
            kk = ( n*( n + 1 ) )/2
            if ( incx == 1 ) then
               do j = n, 1, -1
                  if ( x(j) /= zero ) then
                     temp = x(j)
                     k    = kk
                     do i = n, j + 1, -1
                        x(i) = x(i) + temp * ap(k)
                        k      = k      - 1
                     end do
                     if ( nounit ) &
                        x(j) = x(j) * ap( kk - n + j )
                  end if
                  kk = kk - ( n - j + 1 )
               end do
            else
               kx = kx + ( n - 1 ) * incx
               jx = kx
               do j = n, 1, -1
                  if ( x(jx) /= zero ) then
                     temp = x(jx)
                     ix = kx
                     do k = kk, kk - ( n - ( j + 1 ) ), -1
                        x(ix) = x(ix) + temp * ap(k)
                        ix = ix - incx
                     end do
                     if ( nounit ) &
                        x(jx) = x(jx) * ap( kk - n + j )
                  end if
                  jx = jx - incx
                  kk = kk - ( n - j + 1 )
               end do
            end if
         end if
      else
!
!  Form  x := A'*x  or  x := conjg( A' )*x.
!
         if ( lsame ( uplo, 'U' ) ) then
            kk = ( n*( n + 1 ) )/2
            if ( incx == 1 ) then
               do j = n, 1, -1
                  temp = x(j)
                  k    = kk     - 1
                  if ( noconj ) then
                     if ( nounit ) then
                        temp = temp * ap( kk )
                     end if
                     do i = j - 1, 1, -1
                        temp = temp + ap(k) * x(i)
                        k    = k    - 1
                     end do
                  else
                     if ( nounit ) then
                        temp = temp*conjg ( ap( kk ) )
                     end if
                     do i = j - 1, 1, -1
                        temp = temp + conjg ( ap(k) ) * x(i)
                        k    = k    - 1
                     end do
                  end if
                  x(j) = temp
                  kk     = kk   - j
               end do
            else
               jx = kx + ( n - 1 ) * incx
               do j = n, 1, -1
                  temp = x(jx)
                  ix = jx
                  if ( noconj ) then
                     if ( nounit ) then
                        temp = temp * ap( kk )
                     end if
                     do k = kk - 1, kk - j + 1, -1
                        ix = ix - incx
                        temp = temp + ap(k) * x(ix)
                     end do
                  else
                     if ( nounit ) &
                        temp = temp*conjg ( ap( kk ) )
                     do k = kk - 1, kk - j + 1, -1
                        ix = ix - incx
                        temp = temp + conjg ( ap(k) ) * x(ix)
                     end do
                  end if
                  x(jx) = temp
                  jx = jx - incx
                  kk      = kk   - j
               end do
            end if
         else
            kk = 1
            if ( incx == 1 ) then
               do j = 1, n
                  temp = x(j)
                  k    = kk     + 1
                  if ( noconj ) then
                     if ( nounit ) &
                        temp = temp * ap( kk )
                     do i = j + 1, n
                        temp = temp + ap(k) * x(i)
                        k    = k    + 1
                     end do
                  else
                     if ( nounit ) &
                        temp = temp * conjg ( ap( kk ) )
                     do i = j + 1, n
                        temp = temp + conjg ( ap(k) ) * x(i)
                        k    = k    + 1
                     end do
                  end if
                  x(j) = temp
                  kk     = kk   + ( n - j + 1 )
               end do
            else
               jx = kx
               do j = 1, n
                  temp = x(jx)
                  ix = jx
                  if ( noconj ) then
                     if ( nounit ) &
                        temp = temp * ap( kk )
                     do k = kk + 1, kk + n - j
                        ix = ix + incx
                        temp = temp + ap(k) * x(ix)
                     end do
                  else
                     if ( nounit ) &
                        temp = temp*conjg ( ap( kk ) )
                     do k = kk + 1, kk + n - j
                        ix = ix + incx
                        temp = temp + conjg ( ap(k) ) * x(ix)
                     end do
                  end if
                  x(jx) = temp
                  jx = jx + incx
                  kk      = kk   + ( n - j + 1 )
               end do
            end if
         end if
      end if

  return
end
subroutine ctpsv ( uplo, trans, diag, n, ap, x, incx )

!*****************************************************************************80
!
!! CTPSV solves one of the systems of equations
!
!     A*x = b,   or   A'*x = b,   or   conjg( A' )*x = b,
!
!  where b and x are n element vectors and A is an n by n unit, or
!  non-unit, upper or lower triangular matrix, supplied in packed form.
!
!  No test for singularity or near-singularity is included in this
!  routine. Such tests must be performed before calling this routine.
!
!  Parameters:
!
!  UPLO   - character.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - character.
!           On entry, TRANS specifies the equations to be solved as
!           follows:
!
!              TRANS = 'N' or 'n'   A*x = b.
!
!              TRANS = 'T' or 't'   A'*x = b.
!
!              TRANS = 'C' or 'c'   conjg( A' )*x = b.
!
!           Unchanged on exit.
!
!  DIAG   - character.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  AP     - complex          array of DIMENSION at least
!           ( ( n*( n + 1 ) )/2 ).
!           Before entry with  UPLO = 'U' or 'u', the array AP must
!           contain the upper triangular matrix packed sequentially,
!           column by column, so that ap( 1 ) contains a( 1, 1 ),
!           ap( 2 ) and ap( 3 ) contain a( 1, 2 ) and a( 2, 2 )
!           respectively, and so on.
!           Before entry with UPLO = 'L' or 'l', the array AP must
!           contain the lower triangular matrix packed sequentially,
!           column by column, so that ap( 1 ) contains a( 1, 1 ),
!           ap( 2 ) and ap( 3 ) contain a( 2, 1 ) and a( 3, 1 )
!           respectively, and so on.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  X      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incx ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer incx
  integer n
      character        diag, trans, uplo
      complex            ap( * ), x( * )
      complex, parameter :: zero = ( 0.0E+00, 0.0E+00 )
      complex            temp
      integer            i, info, ix, j, jx, k, kk, kx
      logical            noconj, nounit
      logical            lsame
!
!  Test the input.
!
      info = 0
      if   ( .not.lsame ( uplo , 'U' ) .and.  &
               .not.lsame ( uplo , 'L' )      ) then
         info = 1
      else if ( .not.lsame ( trans, 'N' ) .and.  &
               .not.lsame ( trans, 'T' ) .and.  &
               .not.lsame ( trans, 'C' )      ) then
         info = 2
      else if ( .not.lsame ( diag , 'U' ) .and.  &
               .not.lsame ( diag , 'N' )      ) then
         info = 3
      else if ( n<0 ) then
         info = 4
      else if ( incx == 0 ) then
         info = 7
      end if
      if ( info/=0 ) then
         call xerbla ( 'ctpsv ', info )
         return
      end if
!
!  Quick return if possible.
!
      if ( n == 0 ) &
         return

      noconj = lsame ( trans, 'T' )
      nounit = lsame ( diag , 'N' )
!
!  Set up the start point in X if the increment is not unity. This
!  will be  ( N - 1 ) * incx  too small for descending loops.
!
      if ( incx <= 0 ) then
         kx = 1 - ( n - 1 ) * incx
      else if ( incx/=1 ) then
         kx = 1
      end if
!
!  Start the operations. In this version the elements of AP are
!  accessed sequentially with one pass through AP.
!
      if ( lsame ( trans, 'N' ) ) then
!
!  Form  x := inv( A )*x.
!
         if ( lsame ( uplo, 'U' ) ) then
            kk = ( n*( n + 1 ) )/2
            if ( incx == 1 ) then
               do j = n, 1, -1
                  if ( x(j) /= zero ) then
                     if ( nounit ) &
                        x(j) = x(j)/ap( kk )
                     temp = x(j)
                     k    = kk     - 1
                     do i = j - 1, 1, -1
                        x(i) = x(i) - temp * ap(k)
                        k      = k      - 1
                     end do
                  end if
                  kk = kk - j
               end do
            else
               jx = kx + ( n - 1 ) * incx
               do j = n, 1, -1
                  if ( x(jx) /= zero ) then
                     if ( nounit ) &
                        x(jx) = x(jx)/ap( kk )
                     temp = x(jx)
                     ix = jx
                     do k = kk - 1, kk - j + 1, -1
                        ix = ix - incx
                        x(ix) = x(ix) - temp * ap(k)
                     end do
                  end if
                  jx = jx - incx
                  kk = kk - j
               end do
            end if
         else
            kk = 1
            if ( incx == 1 ) then
               do j = 1, n
                  if ( x(j) /= zero ) then
                     if ( nounit ) &
                        x(j) = x(j)/ap( kk )
                     temp = x(j)
                     k    = kk     + 1
                     do i = j + 1, n
                        x(i) = x(i) - temp * ap(k)
                        k      = k      + 1
                     end do
                  end if
                  kk = kk + ( n - j + 1 )
               end do
            else
               jx = kx
               do j = 1, n
                  if ( x(jx) /= zero ) then
                     if ( nounit ) &
                        x(jx) = x(jx)/ap( kk )
                     temp = x(jx)
                     ix = jx
                     do k = kk + 1, kk + n - j
                        ix = ix + incx
                        x(ix) = x(ix) - temp * ap(k)
                     end do
                  end if
                  jx = jx + incx
                  kk = kk + ( n - j + 1 )
               end do
            end if
         end if
      else
!
!  Form  x := inv( A' )*x  or  x := inv( conjg( A' ) )*x.
!
         if ( lsame ( uplo, 'U' ) ) then
            kk = 1
            if ( incx == 1 ) then
               do j = 1, n
                  temp = x(j)
                  k    = kk
                  if ( noconj ) then
                     do i = 1, j - 1
                        temp = temp - ap(k) * x(i)
                        k    = k    + 1
                     end do
                     if ( nounit ) then
                        temp = temp/ap( kk + j - 1 )
                     end if
                  else
                     do i = 1, j - 1
                        temp = temp - conjg ( ap(k) ) * x(i)
                        k    = k    + 1
                     end do
                     if ( nounit ) &
                        temp = temp/conjg ( ap( kk + j - 1 ) )
                  end if
                  x(j) = temp
                  kk     = kk   + j
               end do
            else
               jx = kx
               do j = 1, n
                  temp = x(jx)
                  ix = kx
                  if ( noconj ) then
                     do k = kk, kk + j - 2
                        temp = temp - ap(k) * x(ix)
                        ix = ix + incx
                     end do
                     if ( nounit ) then
                        temp = temp/ap( kk + j - 1 )
                     end if
                  else
                     do k = kk, kk + j - 2
                        temp = temp - conjg ( ap(k) ) * x(ix)
                        ix = ix + incx
                     end do
                     if ( nounit ) then
                        temp = temp / conjg ( ap( kk + j - 1 ) )
                     end if
                  end if
                  x(jx) = temp
                  jx = jx + incx
                  kk      = kk   + j
               end do
            end if
         else
            kk = ( n * ( n + 1 ) )/2
            if ( incx == 1 ) then
               do j = n, 1, -1
                  temp = x(j)
                  k    = kk
                  if ( noconj ) then
                     do i = n, j + 1, -1
                        temp = temp - ap(k) * x(i)
                        k    = k    - 1
                     end do
                     if ( nounit ) &
                        temp = temp / ap( kk - n + j )
                  else
                     do i = n, j + 1, -1
                        temp = temp - conjg ( ap(k) ) * x(i)
                        k    = k    - 1
                     end do
                     if ( nounit ) &
                        temp = temp / conjg ( ap( kk - n + j ) )
                  end if
                  x(j) = temp
                  kk     = kk   - ( n - j + 1 )
               end do
            else
               kx = kx + ( n - 1 ) * incx
               jx = kx
               do j = n, 1, -1
                  temp = x(jx)
                  ix = kx
                  if ( noconj ) then
                     do k = kk, kk - ( n - ( j + 1 ) ), -1
                        temp = temp - ap(k) * x(ix)
                        ix = ix - incx
                     end do
                     if ( nounit ) then
                        temp = temp/ap( kk - n + j )
                     end if
                  else
                     do k = kk, kk - ( n - ( j + 1 ) ), -1
                        temp = temp - conjg ( ap(k) ) * x(ix)
                        ix = ix - incx
                     end do
                     if ( nounit ) then
                        temp = temp/conjg ( ap( kk - n + j ) )
                     end if
                  end if
                  x(jx) = temp
                  jx = jx - incx
                  kk      = kk   - ( n - j + 1 )
               end do
            end if
         end if
      end if

  return
end
subroutine ctrmv ( uplo, trans, diag, n, a, lda, x, incx )

!*****************************************************************************80
!
!! CTRMV performs one of the matrix-vector operations
!
!     x := A*x,   or   x := A'*x,   or   x := conjg( A' )*x,
!
!  where x is an n element vector and  A is an n by n unit, or non-unit,
!  upper or lower triangular matrix.
!
!  Parameters:
!
!  UPLO   - character.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - character.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   x := A*x.
!
!              TRANS = 'T' or 't'   x := A'*x.
!
!              TRANS = 'C' or 'c'   x := conjg( A' )*x.
!
!           Unchanged on exit.
!
!  DIAG   - character.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  A      - complex          array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced either, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incx ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer lda

  integer incx
  integer n
      character        diag, trans, uplo
      complex            a( lda, * ), x( * )
      complex, parameter :: zero = ( 0.0E+00, 0.0E+00 )
      complex            temp
      integer            i, info, ix, j, jx, kx
      logical            noconj, nounit
      logical            lsame
!
!  Test the input.
!
      info = 0
      if   ( .not.lsame ( uplo , 'U' ) .and.  &
               .not.lsame ( uplo , 'L' )      ) then
         info = 1
      else if ( .not.lsame ( trans, 'N' ) .and.  &
               .not.lsame ( trans, 'T' ) .and.  &
               .not.lsame ( trans, 'C' )      ) then
         info = 2
      else if ( .not.lsame ( diag , 'U' ) .and.  &
               .not.lsame ( diag , 'N' )      ) then
         info = 3
      else if ( n<0 ) then
         info = 4
      else if ( lda<max ( 1, n ) ) then
         info = 6
      else if ( incx == 0 ) then
         info = 8
      end if
      if ( info/=0 ) then
         call xerbla ( 'ctrmv ', info )
         return
      end if
!
!  Quick return if possible.
!
      if ( n == 0 ) then
         return
      end if

      noconj = lsame ( trans, 'T' )
      nounit = lsame ( diag , 'N' )
!
!  Set up the start point in X if the increment is not unity. This
!  will be  ( N - 1 ) * incx  too small for descending loops.
!
      if ( incx <= 0 ) then
         kx = 1 - ( n - 1 ) * incx
      else if ( incx /= 1 ) then
         kx = 1
      end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through A.
!
      if ( lsame ( trans, 'N' ) ) then
!
!  Form  x := A*x.
!
         if ( lsame ( uplo, 'U' ) ) then
            if ( incx == 1 ) then
               do j = 1, n
                  if ( x(j) /= zero ) then
                     temp = x(j)
                     do i = 1, j - 1
                        x(i) = x(i) + temp * a(i,j)
                     end do
                     if ( nounit ) then
                        x(j) = x(j) * a(j,j)
                     end if
                  end if
               end do
            else
               jx = kx
               do j = 1, n
                  if ( x(jx) /= zero ) then
                     temp = x(jx)
                     ix = kx
                     do i = 1, j - 1
                        x(ix) = x(ix) + temp * a(i,j)
                        ix = ix + incx
                     end do
                     if ( nounit ) then
                        x(jx) = x(jx) * a(j,j)
                     end if
                  end if
                  jx = jx + incx
               end do
            end if
         else
            if ( incx == 1 ) then
               do j = n, 1, -1
                  if ( x(j) /= zero ) then
                     temp = x(j)
                     do i = n, j + 1, -1
                        x(i) = x(i) + temp * a(i,j)
                     end do
                     if ( nounit ) then
                        x(j) = x(j) * a(j,j)
                     end if
                  end if
               end do
            else
               kx = kx + ( n - 1 ) * incx
               jx = kx
               do j = n, 1, -1
                  if ( x(jx) /= zero ) then
                     temp = x(jx)
                     ix = kx
                     do i = n, j + 1, -1
                        x(ix) = x(ix) + temp * a(i,j)
                        ix = ix - incx
                     end do
                     if ( nounit ) then
                        x(jx) = x(jx) * a(j,j)
                     end if
                  end if
                  jx = jx - incx
               end do
            end if
         end if
      else
!
!  Form  x := A'*x  or  x := conjg( A' )*x.
!
         if ( lsame ( uplo, 'U' ) ) then
            if ( incx == 1 ) then
               do j = n, 1, -1
                  temp = x(j)
                  if ( noconj ) then
                     if ( nounit ) then
                        temp = temp * a(j,j)
                     end if
                     do i = j - 1, 1, -1
                        temp = temp + a(i,j) * x(i)
                     end do
                  else
                     if ( nounit ) then
                        temp = temp*conjg ( a(j,j) )
                     end if
                     do i = j - 1, 1, -1
                        temp = temp + conjg ( a(i,j) ) * x(i)
                     end do
                  end if
                  x(j) = temp
               end do
            else
               jx = kx + ( n - 1 ) * incx
               do j = n, 1, -1
                  temp = x(jx)
                  ix = jx
                  if ( noconj ) then
                     if ( nounit ) then
                        temp = temp * a(j,j)
                     end if
                     do i = j - 1, 1, -1
                        ix = ix - incx
                        temp = temp + a(i,j) * x(ix)
                     end do
                  else
                     if ( nounit ) then
                        temp = temp * conjg ( a(j,j) )
                     end if
                     do i = j - 1, 1, -1
                        ix = ix - incx
                        temp = temp + conjg ( a(i,j) ) * x(ix)
                     end do
                  end if
                  x(jx) = temp
                  jx = jx - incx
               end do
            end if
         else
            if ( incx == 1 ) then
               do j = 1, n
                  temp = x(j)
                  if ( noconj ) then
                     if ( nounit ) then
                        temp = temp * a(j,j)
                     end if
                     do i = j + 1, n
                        temp = temp + a(i,j) * x(i)
                     end do
                  else
                     if ( nounit ) then
                        temp = temp * conjg ( a(j,j) )
                     end if
                     do i = j + 1, n
                        temp = temp + conjg ( a(i,j) ) * x(i)
                     end do
                  end if
                  x(j) = temp
               end do
            else
               jx = kx
               do j = 1, n
                  temp = x(jx)
                  ix = jx
                  if ( noconj ) then
                     if ( nounit ) then
                        temp = temp * a(j,j)
                     end if
                     do i = j + 1, n
                        ix = ix + incx
                        temp = temp + a(i,j) * x(ix)
                     end do
                  else
                     if ( nounit ) then
                        temp = temp * conjg ( a(j,j) )
                     end if
                     do i = j + 1, n
                        ix = ix + incx
                        temp = temp + conjg ( a(i,j) ) * x(ix)
                     end do
                  end if
                  x(jx) = temp
                  jx = jx + incx
               end do
            end if
         end if
      end if

  return
end
subroutine ctrsv ( uplo, trans, diag, n, a, lda, x, incx )

!*****************************************************************************80
!
!! CTRSV solves one of the systems of equations
!
!     A*x = b,   or   A'*x = b,   or   conjg( A' )*x = b,
!
!  where b and x are n element vectors and A is an n by n unit, or
!  non-unit, upper or lower triangular matrix.
!
!  No test for singularity or near-singularity is included in this
!  routine. Such tests must be performed before calling this routine.
!
!  Parameters:
!
!  UPLO   - character.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - character.
!           On entry, TRANS specifies the equations to be solved as
!           follows:
!
!              TRANS = 'N' or 'n'   A*x = b.
!
!              TRANS = 'T' or 't'   A'*x = b.
!
!              TRANS = 'C' or 'c'   conjg( A' )*x = b.
!
!           Unchanged on exit.
!
!  DIAG   - character.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  A      - complex          array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced either, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - complex          array of dimension at least
!           ( 1 + ( n - 1 ) * abs( incx ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.
!
!  incx   - integer.
!           On entry, incx specifies the increment for the elements of
!           X. incx must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer lda

  integer incx
  integer n
      character        diag, trans, uplo
      complex            a( lda, * ), x( * )
      complex, parameter :: zero = ( 0.0E+00, 0.0E+00 )
      complex            temp
      integer            i, info, ix, j, jx, kx
      logical            noconj, nounit
      logical            lsame
!
!  Test the input.
!
      info = 0
      if   ( .not.lsame ( uplo , 'U' ) .and.  &
               .not.lsame ( uplo , 'L' )      ) then
         info = 1
      else if ( .not.lsame ( trans, 'N' ) .and.  &
               .not.lsame ( trans, 'T' ) .and.  &
               .not.lsame ( trans, 'C' )      ) then
         info = 2
      else if ( .not.lsame ( diag , 'U' ) .and.  &
               .not.lsame ( diag , 'N' )      ) then
         info = 3
      else if ( n<0 ) then
         info = 4
      else if ( lda<max ( 1, n ) ) then
         info = 6
      else if ( incx == 0 ) then
         info = 8
      end if
      if ( info/=0 ) then
         call xerbla ( 'ctrsv ', info )
         return
      end if
!
!  Quick return if possible.
!
      if ( n == 0 ) then
         return
      end if

      noconj = lsame ( trans, 'T' )
      nounit = lsame ( diag , 'N' )
!
!  Set up the start point in X if the increment is not unity. This
!  will be  ( N - 1 ) * incx  too small for descending loops.
!
      if ( incx <= 0 ) then
         kx = 1 - ( n - 1 ) * incx
      else if ( incx/=1 ) then
         kx = 1
      end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through A.
!
      if ( lsame ( trans, 'N' ) ) then
!
!  Form  x := inv( A )*x.
!
         if ( lsame ( uplo, 'U' ) ) then
            if ( incx == 1 ) then
               do j = n, 1, -1
                  if ( x(j) /= zero ) then
                     if ( nounit ) then
                        x(j) = x(j)/a(j,j)
                     end if
                     temp = x(j)
                     do i = j - 1, 1, -1
                        x(i) = x(i) - temp * a(i,j)
                     end do
                  end if
               end do
            else
               jx = kx + ( n - 1 ) * incx
               do j = n, 1, -1
                  if ( x(jx) /= zero ) then
                     if ( nounit ) &
                        x(jx) = x(jx)/a(j,j)
                     temp = x(jx)
                     ix = jx
                     do i = j - 1, 1, -1
                        ix = ix - incx
                        x(ix) = x(ix) - temp * a(i,j)
                     end do
                  end if
                  jx = jx - incx
               end do
            end if
         else
            if ( incx == 1 ) then
               do j = 1, n
                  if ( x(j) /= zero ) then
                     if ( nounit ) then
                        x(j) = x(j) / a(j,j)
                     end if
                     temp = x(j)
                     do i = j + 1, n
                        x(i) = x(i) - temp * a(i,j)
                     end do
                  end if
               end do
            else
               jx = kx
               do j = 1, n
                  if ( x(jx) /= zero ) then
                     if ( nounit ) then
                        x(jx) = x(jx) / a(j,j)
                     end if
                     temp = x(jx)
                     ix = jx
                     do i = j + 1, n
                        ix = ix + incx
                        x(ix) = x(ix) - temp * a(i,j)
                     end do
                  end if
                  jx = jx + incx
               end do
            end if
         end if
      else
!
!  Form  x := inv( A' )*x  or  x := inv( conjg( A' ) )*x.
!
         if ( lsame ( uplo, 'U' ) ) then
            if ( incx == 1 ) then
               do j = 1, n
                  temp = x(j)
                  if ( noconj ) then
                     do i = 1, j - 1
                        temp = temp - a(i,j) * x(i)
                     end do
                     if ( nounit ) then
                        temp = temp/a(j,j)
                     end if
                  else
                     do i = 1, j - 1
                        temp = temp - conjg ( a(i,j) ) * x(i)
                     end do
                     if ( nounit ) then
                        temp = temp/conjg ( a(j,j) )
                     end if
                  end if
                  x(j) = temp
               end do
            else
               jx = kx
               do j = 1, n
                  ix = kx
                  temp = x(jx)
                  if ( noconj ) then
                     do i = 1, j - 1
                        temp = temp - a(i,j) * x(ix)
                        ix = ix + incx
                     end do
                     if ( nounit ) then
                        temp = temp/a(j,j)
                     end if
                  else
                     do i = 1, j - 1
                        temp = temp - conjg ( a(i,j) ) * x(ix)
                        ix = ix + incx
                     end do
                     if ( nounit ) then
                        temp = temp/conjg ( a(j,j) )
                     end if
                  end if
                  x(jx) = temp
                  jx = jx + incx
               end do
            end if
         else
            if ( incx == 1 ) then
               do j = n, 1, -1
                  temp = x(j)
                  if ( noconj ) then
                     do i = n, j + 1, -1
                        temp = temp - a(i,j) * x(i)
                     end do
                     if ( nounit ) then
                        temp = temp/a(j,j)
                     end if
                  else
                     do i = n, j + 1, -1
                        temp = temp - conjg ( a(i,j) ) * x(i)
                     end do
                     if ( nounit ) then
                        temp = temp/conjg ( a(j,j) )
                     end if
                  end if
                  x(j) = temp
               end do
            else
               kx = kx + ( n - 1 ) * incx
               jx = kx
               do j = n, 1, -1
                  ix = kx
                  temp = x(jx)
                  if ( noconj ) then
                     do i = n, j + 1, -1
                        temp = temp - a(i,j) * x(ix)
                        ix = ix - incx
                     end do
                     if ( nounit ) then
                        temp = temp/a(j,j)
                     end if
                  else
                     do i = n, j + 1, -1
                        temp = temp - conjg ( a(i,j) ) * x(ix)
                        ix = ix - incx
                     end do
                     if ( nounit ) then
                        temp = temp/conjg ( a(j,j) )
                     end if
                  end if
                  x(jx) = temp
                  jx = jx - incx
               end do
            end if
         end if
      end if

  return
end
subroutine dgbmv ( trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy )

!*****************************************************************************80
!
!! DGBMV computes y: = alpha*A*x + beta*y for banded matrix A.
!
!  Discussion:
!
!    DGBMV performs one of the matrix-vector operations
!
!      y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!    where alpha and beta are scalars, x and y are vectors and A is an
!    m by n band matrix, with kl sub-diagonals and ku super-diagonals.
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    Jack Dongarra, Argonne National Lab.
!    Jeremy Du Croz, Nag Central Office.
!    Sven Hammarling, Nag Central Office.
!    Richard Hanson, Sandia National Labs.
!
!  Parameters:
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'N'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 'T'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'C'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  KL     - INTEGER.
!           On entry, KL specifies the number of sub-diagonals of the
!           matrix A. KL must satisfy  0 <= KL.
!           Unchanged on exit.
!
!  KU     - INTEGER.
!           On entry, KU specifies the number of super-diagonals of the
!           matrix A. KU must satisfy  0 <= KU.
!           Unchanged on exit.
!
!  ALPHA  - real ( kind = 8 ).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, n ).
!           Before entry, the leading ( kl + ku + 1 ) by n part of the
!           array A must contain the matrix of coefficients, supplied
!           column by column, with the leading diagonal of the matrix in
!           row ( ku + 1 ) of the array, the first super-diagonal
!           starting at position 2 in row ku, the first sub-diagonal
!           starting at position 1 in row ( ku + 2 ), and so on.
!           Elements in the array A that do not correspond to elements
!           in the band matrix (such as the top left ku by ku triangle)
!           are not referenced.
!           The following program segment will transfer a band matrix
!           from conventional full matrix storage to band storage:
!
!                 do J = 1, N
!                    K = KU + 1 - J
!                    do I = max ( 1, J - KU ), min ( M, J + KL )
!                       A( K + I, J ) = matrix( I, J )
!                    end do
!                 end do
!
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( kl + ku + 1 ).
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'N'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - real ( kind = 8 ).
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - real ( kind = 8 ) array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'N'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry, the incremented array Y must contain the
!           vector y. On exit, Y is overwritten by the updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer lda

  real ( kind = 8 ) a(lda,*)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer i
  integer incx
  integer incy
  integer info
  integer ix
  integer iy
  integer j
  integer jx
  integer jy
  integer k
  integer kl
  integer ku
  integer kup1
  integer kx
  integer ky
  integer lenx
  integer leny
  logical, external :: lsame
  integer m
  intrinsic max
  intrinsic min
  integer n
  real ( kind = 8 ) temp
  character trans
  real ( kind = 8 ) x(*)
  external xerbla
  real ( kind = 8 ) y(*)
!
!  Test the input parameters.
!
  info = 0
  if ( .not. lsame ( trans, 'N' ) .and.  &
       .not. lsame ( trans, 'T' ) .and.  &
       .not. lsame ( trans, 'C' )      ) then
    info = 1
  else if ( m < 0 ) then
    info = 2
  else if ( n < 0 ) then
    info = 3
  else if ( kl < 0 ) then
    info = 4
  else if ( ku < 0 ) then
    info = 5
  else if ( lda < ( kl + ku + 1 ) ) then
    info = 8
  else if ( incx == 0 ) then
    info = 10
  else if ( incy == 0 ) then
    info = 13
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dgbmv ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( ( m == 0 ) .or. &
       ( n == 0 ) .or. &
       ( ( alpha == 0.0d+00 ) .and. ( beta == 1.0d+00 ) ) ) then
    return
  end if
!
!  Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!  up the start points in  X  and  Y.
!
  if ( lsame ( trans, 'N' ) ) then
    lenx = n
    leny = m
  else
    lenx = m
    leny = n
  end if

  if ( 0 < incx ) then
    kx = 1
  else
    kx = 1 - ( lenx - 1 ) * incx
  end if

  if ( 0 < incy ) then
    ky = 1
  else
    ky = 1 - ( leny - 1 ) * incy
  end if
!
!  Start the operations.  In this version the elements of A are
!  accessed sequentially with one pass through the band part of A.
!
!  First form y := beta*y.
!
  if ( beta /= 1.0d+00 ) then
    if ( incy == 1 ) then
      if ( beta == 0.0d+00 ) then
        y(1:leny) = 0.0d+00
      else
        y(1:leny) = beta * y(1:leny)
      end if
    else
      iy = ky
      if ( beta == 0.0d+00 ) then
        do i = 1, leny
          y(iy) = 0.0d+00
          iy = iy + incy
        end do
      else
        do i = 1, leny
          y(iy) = beta * y(iy)
          iy = iy + incy
        end do
      end if
    end if
  end if

  if ( alpha == 0.0d+00 ) then
    return
  end if

  kup1 = ku + 1
  if ( lsame ( trans, 'N' ) ) then
!
!  Form  y := alpha*A*x + y.
!
    jx = kx
    if ( incy == 1 ) then
      do j = 1, n
        if ( x(jx) /= 0.0d+00 ) then
          temp = alpha * x(jx)
          k = kup1 - j
          do i = max ( 1, j - ku ), min ( m, j + kl )
            y(i) = y(i) + temp * a(k+i,j)
          end do
        end if
        jx = jx + incx
      end do
    else
      do j = 1, n
        if ( x(jx) /= 0.0d+00 ) then
          temp = alpha * x(jx)
          iy = ky
          k = kup1 - j
          do i = max ( 1, j - ku ), min ( m, j + kl )
            y(iy) = y(iy) + temp * a(k+i,j)
            iy = iy + incy
          end do
        end if
        jx = jx + incx
        if ( ku < j ) then
         ky = ky + incy
        end if

      end do
    end if
  else
!
!  Form y := alpha*A'*x + y.
!
    jy = ky
    if ( incx == 1 ) then
      do j = 1, n
        temp = 0.0d+00
        k = kup1 - j
        do i = max ( 1, j - ku ), min ( m, j + kl )
          temp = temp + a(k+i,j) * x(i)
        end do
        y(jy) = y(jy) + alpha * temp
        jy = jy + incy
      end do
    else
      do j = 1, n
        temp = 0.0d+00
        ix = kx
        k = kup1 - j
        do i = max ( 1, j - ku ), min ( m, j + kl )
          temp = temp + a(k+i,j) * x(ix)
          ix = ix + incx
        end do
        y(jy) = y(jy) + alpha * temp
        jy = jy + incy

        if ( ku < j ) then
          kx = kx + incx
        end if

      end do
    end if
  end if

  return
end
subroutine dgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

!*****************************************************************************80
!
!! DGEMV computes y := alpha * A * x + beta * y for general matrix A.
!
!  Discussion:
!
!    DGEMV performs one of the matrix-vector operations
!
!      y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!    where alpha and beta are scalars, x and y are vectors and A is an
!    m by n matrix.
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    Jack Dongarra, Argonne National Lab.
!    Jeremy Du Croz, Nag Central Office.
!    Sven Hammarling, Nag Central Office.
!    Richard Hanson, Sandia National Labs.
!
!  Parameters:
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'N'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 'T'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'C'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real ( kind = 8 ).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'N'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - real ( kind = 8 ).
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - real ( kind = 8 ) array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'N'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer lda

  real ( kind = 8 ) a(lda,*)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer i
  integer incx
  integer incy
  integer info
  integer ix
  integer iy
  integer j
  integer jx
  integer jy
  integer kx
  integer ky
  integer lenx
  integer leny
  logical, external :: lsame
  integer m
  intrinsic max
  integer n
  real ( kind = 8 ) temp
  character trans
  real ( kind = 8 ) x(*)
  external xerbla
  real ( kind = 8 ) y(*)
!
!  Test the input parameters.
!
  info = 0
  if ( .not. lsame ( trans, 'N' ) .and.  &
       .not. lsame ( trans, 'T' ) .and.  &
       .not. lsame ( trans, 'C' )      ) then
    info = 1
  else if ( m < 0 ) then
    info = 2
  else if ( n < 0 ) then
    info = 3
  else if ( lda < max ( 1, m ) ) then
    info = 6
  else if ( incx == 0 ) then
    info = 8
  else if ( incy == 0 ) then
    info = 11
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dgemv ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( ( m == 0 ) .or. &
       ( n == 0 ) .or. &
       ( ( alpha == 0.0d+00 ) .and. ( beta == 1.0d+00 ) ) ) then
   return
  end if
!
!  Set LENX and LENY, the lengths of the vectors x and y, and set
!  up the start points in X and Y.
!
  if ( lsame ( trans, 'N' ) ) then
    lenx = n
    leny = m
  else
    lenx = m
    leny = n
  end if

  if ( 0 < incx ) then
    kx = 1
  else
    kx = 1 - ( lenx - 1 ) * incx
  end if

  if ( 0 < incy ) then
    ky = 1
  else
    ky = 1 - ( leny - 1 ) * incy
  end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through A.
!
!  First form  y := beta*y.
!
  if ( beta /= 1.0d+00 ) then
    if ( incy == 1 ) then
      if ( beta == 0.0d+00 ) then
        y(1:leny) = 0.0d+00
      else
        y(1:leny) = beta * y(1:leny)
      end if
    else
      iy = ky
      if ( beta == 0.0d+00 ) then
        do i = 1, leny
          y(iy) = 0.0d+00
          iy = iy + incy
        end do
      else
        do i = 1, leny
          y(iy) = beta * y(iy)
          iy = iy + incy
        end do
      end if
    end if
  end if

  if ( alpha == 0.0d+00 ) then
   return
  end if

  if ( lsame ( trans, 'N' ) ) then
!
!  Form  y := alpha*A*x + y.
!
    jx = kx
    if ( incy == 1 ) then
      do j = 1, n
        if ( x(jx) /= 0.0d+00 ) then
          temp = alpha * x(jx)
          do i = 1, m
            y(i) = y(i) + temp * a(i,j)
          end do
        end if
        jx = jx + incx
      end do
    else
      do j = 1, n
        if ( x(jx) /= 0.0d+00 ) then
          temp = alpha * x(jx)
          iy = ky
          do i = 1, m
            y(iy) = y(iy) + temp * a(i,j)
            iy = iy + incy
          end do
        end if
        jx = jx + incx
      end do
    end if
  else
!
!  Form  y := alpha*A'*x + y.
!
    jy = ky
    if ( incx == 1 ) then
      do j = 1, n
        temp = 0.0d+00
        do i = 1, m
          temp = temp + a(i,j) * x(i)
        end do
        y(jy) = y(jy) + alpha * temp
        jy = jy + incy
      end do
    else
      do j = 1, n
        temp = 0.0d+00
        ix = kx
        do i = 1, m
          temp = temp + a(i,j) * x(ix)
          ix = ix + incx
        end do
        y(jy) = y(jy) + alpha * temp
        jy = jy + incy
      end do
    end if
  end if

  return
end
subroutine dger ( m, n, alpha, x, incx, y, incy, a, lda )

!*****************************************************************************80
!
!! DGER computes A := alpha*x*y' + A.
!
!  Discussion:
!
!    DGER performs the rank 1 operation
!
!      A := alpha*x*y' + A,
!
!    where alpha is a scalar, x is an m element vector, y is an n element
!    vector and A is an m by n matrix.
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    Jack Dongarra, Argonne National Lab.
!    Jeremy Du Croz, Nag Central Office.
!    Sven Hammarling, Nag Central Office.
!    Richard Hanson, Sandia National Labs.
!
!  Parameters:
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real ( kind = 8 ).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
  implicit none

  real ( kind = 8 )   alpha
  integer            incx, incy, lda, m, n
  real ( kind = 8 )   a( lda, * ), x( * ), y( * )
  real ( kind = 8 )   temp
  integer            i, info, ix, j, jy, kx
  external           xerbla
  intrinsic          max
!
!  Test the input parameters.
!
  info = 0
  if ( m < 0 ) then
    info = 1
  else if ( n < 0 ) then
    info = 2
  else if ( incx == 0 ) then
    info = 5
  else if ( incy == 0 ) then
    info = 7
  else if ( lda < max ( 1, m ) ) then
    info = 9
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dger  ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( ( m == 0 ) .or. ( n == 0 ) .or. ( alpha == 0.0d+00 ) ) then
    return
  end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through A.
!
  if ( 0 < incy ) then
    jy = 1
  else
    jy = 1 - ( n - 1 ) * incy
  end if

  if ( incx == 1 ) then
    do j = 1, n
      if ( y(jy) /= 0.0d+00 ) then
        temp = alpha * y(jy)
        do i = 1, m
          a(i,j) = a(i,j) + x(i) * temp
        end do
      end if
      jy = jy + incy
    end do
  else
    if ( 0 < incx ) then
      kx = 1
    else
      kx = 1 - ( m - 1 ) * incx
    end if
    do j = 1, n
      if ( y(jy) /= 0.0d+00 ) then
        temp = alpha * y(jy)
        ix = kx
        do i = 1, m
          a(i,j) = a(i,j) + x(ix) * temp
          ix = ix + incx
        end do
      end if
      jy = jy + incy
    end do
  end if

  return
end
subroutine dsbmv ( uplo, n, k, alpha, a, lda, x, incx, beta, y, incy )

!*****************************************************************************80
!
!! DSBMV computes y := alpha*A*x + beta*y for symmetric banded matrix A.
!
!  Discussion:
!
!    DSBMV performs the matrix-vector  operation
!
!      y := alpha*A*x + beta*y,
!
!    where alpha and beta are scalars, x and y are n element vectors and
!    A is an n by n symmetric band matrix, with k super-diagonals.
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    Jack Dongarra, Argonne National Lab.
!    Jeremy Du Croz, Nag Central Office.
!    Sven Hammarling, Nag Central Office.
!    Richard Hanson, Sandia National Labs.
!
!  Parameters:
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the band matrix A is being supplied as
!           follows:
!
!              UPLO = 'U' or 'U'   The upper triangular part of A is
!                                  being supplied.
!
!              UPLO = 'L' or 'L'   The lower triangular part of A is
!                                  being supplied.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry, K specifies the number of super-diagonals of the
!           matrix A. K must satisfy  0 <= K.
!           Unchanged on exit.
!
!  ALPHA  - real ( kind = 8 ).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, n ).
!           Before entry with UPLO = 'U' or 'U', the leading ( k + 1 )
!           by n part of the array A must contain the upper triangular
!           band part of the symmetric matrix, supplied column by
!           column, with the leading diagonal of the matrix in row
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. The top left k by k triangle
!           of the array A is not referenced.
!           The following program segment will transfer the upper
!           triangular part of a symmetric band matrix from conventional
!           full matrix storage to band storage:
!
!                 do J = 1, N
!                    M = K + 1 - J
!                    do I = max ( 1, J - K ), J
!                       A( M + I, J ) = matrix( I, J )
!                    end do
!                 end do
!
!           Before entry with UPLO = 'L' or 'L', the leading ( k + 1 )
!           by n part of the array A must contain the lower triangular
!           band part of the symmetric matrix, supplied column by
!           column, with the leading diagonal of the matrix in row 1 of
!           the array, the first sub-diagonal starting at position 1 in
!           row 2, and so on. The bottom right k by k triangle of the
!           array A is not referenced.
!           The following program segment will transfer the lower
!           triangular part of a symmetric band matrix from conventional
!           full matrix storage to band storage:
!
!                 do J = 1, N
!                    M = 1 - J
!                    do I = J, min ( N, J + K )
!                       A( M + I, J ) = matrix( I, J )
!                    end do
!                 end do
!
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( k + 1 ).
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - real ( kind = 8 ).
!           On entry, BETA specifies the scalar beta.
!           Unchanged on exit.
!
!  Y      - real ( kind = 8 ) array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the
!           vector y. On exit, Y is overwritten by the updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
  implicit none

  real ( kind = 8 )   alpha, beta
  integer            incx, incy, k, lda, n
  character       uplo
  real ( kind = 8 )   a( lda, * ), x( * ), y( * )
  real ( kind = 8 )   temp1, temp2
  integer            i, info, ix, iy, j, jx, jy, kplus1, kx, ky, l
  logical, external :: lsame
  external           xerbla
  intrinsic          max, min
!
!  Test the input parameters.
!
  info = 0
  if ( .not.  lsame ( uplo, 'U' ) .and.  &
       .not.  lsame ( uplo, 'L' )      ) then
    info = 1
  else if ( n < 0 ) then
    info = 2
  else if ( k < 0 ) then
    info = 3
  else if ( lda < ( k + 1 ) ) then
    info = 6
  else if ( incx == 0 ) then
    info = 8
  else if ( incy == 0 ) then
    info = 11
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dsbmv ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( ( n == 0 ) .or. ( ( alpha == 0.0d+00 ) .and. ( beta == 1.0d+00 ) ) ) then
    return
  end if
!
!  Set up the start points in  X  and  Y.
!
  if ( 0 < incx ) then
    kx = 1
  else
    kx = 1 - ( n - 1 ) * incx
  end if

  if ( 0 < incy ) then
    ky = 1
  else
    ky = 1 - ( n - 1 ) * incy
  end if
!
!  Start the operations. In this version the elements of the array A
!  are accessed sequentially with one pass through A.
!
!  First form y := beta*y.
!
  if ( beta /= 1.0d+00 ) then
    if ( incy == 1 ) then
      if ( beta == 0.0d+00 ) then
        do i = 1, n
          y(i) = 0.0d+00
        end do
      else
        do i = 1, n
          y(i) = beta * y(i)
        end do
      end if
    else
      iy = ky
      if ( beta == 0.0d+00 ) then
        do i = 1, n
          y(iy) = 0.0d+00
          iy = iy + incy
        end do
      else
        do i = 1, n
          y(iy) = beta * y(iy)
          iy = iy + incy
        end do
      end if
    end if
  end if

  if ( alpha == 0.0d+00 ) then
    return
  end if

  if ( lsame ( uplo, 'U' ) ) then
!
!  Form y when upper triangle of A is stored.
!
    kplus1 = k + 1
    if ( ( incx == 1 ).and.( incy == 1 ) ) then
      do j = 1, n
        temp1 = alpha * x(j)
        temp2 = 0.0d+00
        l = kplus1 - j
        do i = max ( 1, j - k ), j - 1
          y(i) = y(i) + temp1 * a(l+i,j)
          temp2 = temp2 + a(l+i,j) * x(i)
        end do
        y(j) = y(j) + temp1 * a(kplus1,j) + alpha * temp2
      end do
    else
      jx = kx
      jy = ky
      do j = 1, n
        temp1 = alpha * x(jx)
        temp2 = 0.0d+00
        ix = kx
        iy = ky
        l = kplus1 - j
        do i = max ( 1, j - k ), j - 1
          y(iy) = y(iy) + temp1 * a(l+i,j)
          temp2 = temp2 + a(l+i,j) * x(ix)
          ix = ix + incx
          iy = iy + incy
        end do
        y(jy) = y(jy) + temp1 * a(kplus1,j) + alpha * temp2
        jx = jx + incx
        jy = jy + incy
        if ( k < j ) then
          kx = kx + incx
          ky = ky + incy
        end if
      end do
    end if
  else
!
!  Form y when lower triangle of A is stored.
!
    if ( ( incx == 1 ).and.( incy == 1 ) ) then
      do j = 1, n
        temp1 = alpha * x(j)
        temp2 = 0.0d+00
        y(j) = y(j) + temp1 * a(1,j)
        l = 1 - j
        do i = j + 1, min ( n, j + k )
          y(i) = y(i) + temp1 * a(l+i,j)
          temp2 = temp2 + a(l+i,j) * x(i)
        end do
        y(j) = y(j) + alpha * temp2
      end do
    else
      jx = kx
      jy = ky
      do j = 1, n
        temp1 = alpha * x(jx)
        temp2 = 0.0d+00
        y(jy) = y(jy) + temp1 * a(1,j)
        l = 1 - j
        ix = jx
        iy = jy
        do i = j + 1, min ( n, j + k )
          ix = ix + incx
          iy = iy + incy
          y(iy) = y(iy) + temp1 * a(l+i,j)
          temp2 = temp2 + a(l+i,j) * x(ix)
        end do
        y(jy) = y(jy) + alpha * temp2
        jx = jx + incx
        jy = jy + incy
     end do
    end if
  end if

  return
end
subroutine dspmv ( uplo, n, alpha, ap, x, incx, beta, y, incy )

!*****************************************************************************80
!
!! DSPMV computes y := alpha*A*x + beta*y for symmetric packed matrix A.
!
!  Discussion:
!
!    DSPMV performs the matrix-vector operation
!
!      y := alpha*A*x + beta*y,
!
!    where alpha and beta are scalars, x and y are n element vectors and
!    A is an n by n symmetric matrix, supplied in packed form.
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    Jack Dongarra, Argonne National Lab.
!    Jeremy Du Croz, Nag Central Office.
!    Sven Hammarling, Nag Central Office.
!    Richard Hanson, Sandia National Labs.
!
!  Parameters:
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the matrix A is supplied in the packed
!           array AP as follows:
!
!              UPLO = 'U' or 'U'   The upper triangular part of A is
!                                  supplied in AP.
!
!              UPLO = 'L' or 'L'   The lower triangular part of A is
!                                  supplied in AP.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real ( kind = 8 ).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  AP     - real ( kind = 8 ) array of DIMENSION at least
!           ( ( n * ( n + 1 ) ) / 2 ).
!           Before entry with UPLO = 'U' or 'U', the array AP must
!           contain the upper triangular part of the symmetric matrix
!           packed sequentially, column by column, so that AP(1)
!           contains a(1,1), AP(2) and AP(3) contain a(1,2)
!           and a(2,2) respectively, and so on.
!           Before entry with UPLO = 'L' or 'L', the array AP must
!           contain the lower triangular part of the symmetric matrix
!           packed sequentially, column by column, so that AP(1)
!           contains a(1,1), AP(2) and AP(3) contain a(2,1)
!           and a(3,1) respectively, and so on.
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - real ( kind = 8 ).
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y. On exit, Y is overwritten by the updated
!           vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
  implicit none

  real ( kind = 8 )   alpha, beta
  integer            incx, incy, n
  character       uplo
  real ( kind = 8 )   ap( * ), x( * ), y( * )
  real ( kind = 8 )   temp1, temp2
  integer            i, info, ix, iy, j, jx, jy, k, kk, kx, ky
  logical, external :: lsame
  external           xerbla
!
!  Test the input parameters.
!
  info = 0
  if ( .not. lsame ( uplo, 'U' ) .and.  &
       .not. lsame ( uplo, 'L' )      ) then
    info = 1
  else if ( n < 0 ) then
    info = 2
  else if ( incx == 0 ) then
    info = 6
  else if ( incy == 0 ) then
    info = 9
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dspmv ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( ( n == 0 ) .or. ( ( alpha == 0.0d+00 ) .and. ( beta == 1.0d+00 ) ) ) then
    return
  end if
!
!  Set up the start points in X and Y.
!
  if ( 0 < incx ) then
    kx = 1
  else
    kx = 1 - ( n - 1 ) * incx
  end if

  if ( 0 < incy ) then
    ky = 1
  else
    ky = 1 - ( n - 1 ) * incy
  end if
!
!  Start the operations. In this version the elements of the array AP
!  are accessed sequentially with one pass through AP.
!
!  First form  y := beta*y.
!
  if ( beta /= 1.0d+00 ) then
    if ( incy == 1 ) then
      if ( beta == 0.0d+00 ) then
        do i = 1, n
          y(i) = 0.0d+00
        end do
      else
        do i = 1, n
          y(i) = beta * y(i)
        end do
      end if
    else
      iy = ky
      if ( beta == 0.0d+00 ) then
        do i = 1, n
          y(iy) = 0.0d+00
          iy = iy + incy
        end do
      else
        do i = 1, n
          y(iy) = beta * y(iy)
          iy = iy + incy
        end do
      end if
    end if
  end if

  if ( alpha == 0.0d+00 ) then
   return
  end if

  kk = 1
  if ( lsame ( uplo, 'U' ) ) then
!
!  Form y when AP contains the upper triangle.
!
    if ( ( incx == 1 ).and.( incy == 1 ) ) then
      do j = 1, n
        temp1 = alpha * x(j)
        temp2 = 0.0d+00
        k = kk
        do i = 1, j - 1
          y(i) = y(i) + temp1 * ap(k)
          temp2 = temp2 + ap(k) * x(i)
          k = k + 1
        end do
        y(j) = y(j) + temp1 * ap(kk+j-1) + alpha * temp2
        kk = kk + j
      end do
    else
      jx = kx
      jy = ky
      do j = 1, n
        temp1 = alpha * x(jx)
        temp2 = 0.0d+00
        ix = kx
        iy = ky
        do k = kk, kk + j - 2
          y(iy) = y(iy) + temp1 * ap(k)
          temp2 = temp2 + ap(k) * x(ix)
          ix = ix + incx
          iy = iy + incy
        end do
        y(jy) = y(jy) + temp1 * ap(kk+j-1) + alpha * temp2
        jx = jx + incx
        jy = jy + incy
        kk = kk + j
      end do
    end if
  else
!
!  Form y when AP contains the lower triangle.
!
    if ( ( incx == 1 ).and.( incy == 1 ) ) then
      do j = 1, n
        temp1 = alpha * x(j)
        temp2 = 0.0d+00
        y(j) = y(j) + temp1 * ap(kk)
        k = kk + 1
        do i = j + 1, n
          y(i) = y(i) + temp1 * ap(k)
          temp2 = temp2 + ap(k) * x(i)
          k = k + 1
        end do
        y(j) = y(j) + alpha * temp2
        kk = kk + ( n - j + 1 )
      end do
    else
      jx = kx
      jy = ky
      do j = 1, n
        temp1 = alpha * x(jx)
        temp2 = 0.0d+00
        y(jy) = y(jy) + temp1 * ap(kk)
        ix = jx
        iy = jy
        do k = kk + 1, kk + n - j
          ix = ix + incx
          iy = iy + incy
          y(iy) = y(iy) + temp1 * ap(k)
          temp2 = temp2 + ap(k) * x(ix)
        end do
        y(jy) = y(jy) + alpha * temp2
        jx = jx + incx
        jy = jy + incy
        kk = kk + ( n - j + 1 )
      end do
    end if
  end if

  return
end
subroutine dspr2 ( uplo, n, alpha, x, incx, y, incy, ap )

!*****************************************************************************80
!
!! DSPR2 sets A := alpha*x*y' + alpha*y*x' + A, A is a symmetric packed matrix.
!
!  Discussion:
!
!    DSPR2 performs the symmetric rank 2 operation
!
!      A := alpha*x*y' + alpha*y*x' + A,
!
!    where alpha is a scalar, x and y are n element vectors and A is an
!    n by n symmetric matrix, supplied in packed form.
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    Jack Dongarra, Argonne National Lab.
!    Jeremy Du Croz, Nag Central Office.
!    Sven Hammarling, Nag Central Office.
!    Richard Hanson, Sandia National Labs.
!
!  Parameters:
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the matrix A is supplied in the packed
!           array AP as follows:
!
!              UPLO = 'U' or 'U'   The upper triangular part of A is
!                                  supplied in AP.
!
!              UPLO = 'L' or 'L'   The lower triangular part of A is
!                                  supplied in AP.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real ( kind = 8 ).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  AP     - real ( kind = 8 ) array of DIMENSION at least
!           ( ( n * ( n + 1 ) ) / 2 ).
!           Before entry with  UPLO = 'U' or 'U', the array AP must
!           contain the upper triangular part of the symmetric matrix
!           packed sequentially, column by column, so that AP(1)
!           contains a(1,1), AP(2) and AP(3) contain a(1,2)
!           and a(2,2) respectively, and so on. On exit, the array
!           AP is overwritten by the upper triangular part of the
!           updated matrix.
!           Before entry with UPLO = 'L' or 'L', the array AP must
!           contain the lower triangular part of the symmetric matrix
!           packed sequentially, column by column, so that AP(1)
!           contains a(1,1), AP(2) and AP(3) contain a(2,1)
!           and a(3,1) respectively, and so on. On exit, the array
!           AP is overwritten by the lower triangular part of the
!           updated matrix.
!
  implicit none

  real ( kind = 8 )   alpha
  integer            incx, incy, n
  character       uplo
  real ( kind = 8 )   ap( * ), x( * ), y( * )
  real ( kind = 8 )   temp1, temp2
  integer            i, info, ix, iy, j, jx, jy, k, kk, kx, ky
  logical, external :: lsame
  external           xerbla
!
!  Test the input parameters.
!
  info = 0
  if     ( .not. lsame ( uplo, 'U' ) .and.  &
        .not. lsame ( uplo, 'L' )      ) then
    info = 1
  else if ( n < 0 ) then
    info = 2
  else if ( incx == 0 ) then
    info = 5
  else if ( incy == 0 ) then
    info = 7
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dspr2 ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( ( n == 0 ) .or. ( alpha == 0.0d+00 ) ) then
    return
  end if
!
!  Set up the start points in X and Y if the increments are not both
!  unity.
!
  if ( ( incx /= 1 ) .or. ( incy /= 1 ) ) then

    if ( 0 < incx ) then
      kx = 1
    else
      kx = 1 - ( n - 1 ) * incx
    end if

    if ( 0 < incy ) then
      ky = 1
    else
      ky = 1 - ( n - 1 ) * incy
    end if

    jx = kx
    jy = ky

  end if
!
!  Start the operations. In this version the elements of the array AP
!  are accessed sequentially with one pass through AP.
!
  kk = 1
  if ( lsame ( uplo, 'U' ) ) then
!
!  Form  A  when upper triangle is stored in AP.
!
    if ( ( incx == 1 ) .and. ( incy == 1 ) ) then
      do j = 1, n
        if ( ( x(j) /= 0.0d+00 ) .or. ( y(j) /= 0.0d+00 ) ) then
          temp1 = alpha * y(j)
          temp2 = alpha * x(j)
          k = kk
          do i = 1, j
            ap(k) = ap(k) + x(i) * temp1 + y(i) * temp2
            k = k + 1
          end do
        end if
        kk = kk + j
      end do
    else
      do j = 1, n
        if ( ( x(jx) /= 0.0d+00 ) .or. ( y(jy) /= 0.0d+00 ) ) then
          temp1 = alpha * y(jy)
          temp2 = alpha * x(jx)
          ix = kx
          iy = ky
          do k = kk, kk + j - 1
            ap(k) = ap(k) + x(ix) * temp1 + y(iy) * temp2
            ix = ix + incx
            iy = iy + incy
          end do
        end if
        jx = jx + incx
        jy = jy + incy
        kk = kk + j
      end do
    end if
  else
!
!  Form  A  when lower triangle is stored in AP.
!
    if ( ( incx == 1 ) .and. ( incy == 1 ) ) then
      do j = 1, n
        if ( ( x(j) /= 0.0d+00 ) .or. ( y(j) /= 0.0d+00 ) ) then
          temp1 = alpha * y(j)
          temp2 = alpha * x(j)
          k = kk
          do i = j, n
            ap(k) = ap(k) + x(i) * temp1 + y(i) * temp2
            k = k + 1
          end do
        end if
        kk = kk + n - j + 1
      end do
    else
      do j = 1, n
        if ( ( x(jx) /= 0.0d+00 ) .or. ( y(jy) /= 0.0d+00 ) ) then
          temp1 = alpha * y(jy)
          temp2 = alpha * x(jx)
          ix = jx
          iy = jy
          do k = kk, kk + n - j
            ap(k) = ap(k) + x(ix) * temp1 + y(iy) * temp2
            ix = ix + incx
            iy = iy + incy
          end do
        end if
        jx = jx + incx
        jy = jy + incy
        kk = kk + n - j + 1
      end do
    end if
  end if

  return
end
subroutine dspr ( uplo, n, alpha, x, incx, ap )

!*****************************************************************************80
!
!! DSPR computes A := alpha*x*x' + A, where A is a symmetric packed matrix.
!
!  Discussion:
!
!    DSPR performs the symmetric rank 1 operation
!
!      A := alpha*x*x' + A,
!
!    where alpha is a real scalar, x is an n element vector and A is an
!    n by n symmetric matrix, supplied in packed form.
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    Jack Dongarra, Argonne National Lab.
!    Jeremy Du Croz, Nag Central Office.
!    Sven Hammarling, Nag Central Office.
!    Richard Hanson, Sandia National Labs.
!
!  Parameters:
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the matrix A is supplied in the packed
!           array AP as follows:
!
!              UPLO = 'U' or 'U'   The upper triangular part of A is
!                                  supplied in AP.
!
!              UPLO = 'L' or 'L'   The lower triangular part of A is
!                                  supplied in AP.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real ( kind = 8 ).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  AP     - real ( kind = 8 ) array of DIMENSION at least
!           ( ( n * ( n + 1 ) ) / 2 ).
!           Before entry with  UPLO = 'U' or 'U', the array AP must
!           contain the upper triangular part of the symmetric matrix
!           packed sequentially, column by column, so that AP(1)
!           contains a(1,1), AP(2) and AP(3) contain a(1,2)
!           and a(2,2) respectively, and so on. On exit, the array
!           AP is overwritten by the upper triangular part of the
!           updated matrix.
!           Before entry with UPLO = 'L' or 'L', the array AP must
!           contain the lower triangular part of the symmetric matrix
!           packed sequentially, column by column, so that AP(1)
!           contains a(1,1), AP(2) and AP(3) contain a(2,1)
!           and a(3,1) respectively, and so on. On exit, the array
!           AP is overwritten by the lower triangular part of the
!           updated matrix.
!
  implicit none

  real ( kind = 8 )   alpha
  integer            incx, n
  character       uplo
  real ( kind = 8 )   ap( * ), x( * )
  real ( kind = 8 )   temp
  integer            i, info, ix, j, jx, k, kk, kx
  logical, external :: lsame
  external           xerbla
!
!  Test the input parameters.
!
  info = 0
  if     ( .not. lsame ( uplo, 'U' ) .and.  &
         .not. lsame ( uplo, 'L' )      ) then
    info = 1
  else if ( n < 0 ) then
    info = 2
  else if ( incx == 0 ) then
    info = 5
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dspr  ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( ( n == 0 ) .or. ( alpha == 0.0d+00 ) ) then
    return
  end if
!
!  Set the start point in X if the increment is not unity.
!
  if ( incx <= 0 ) then
    kx = 1 - ( n - 1 ) * incx
  else if ( incx /= 1 ) then
    kx = 1
  end if
!
!  Start the operations. In this version the elements of the array AP
!  are accessed sequentially with one pass through AP.
!
  kk = 1
  if ( lsame ( uplo, 'U' ) ) then
!
!  Form A when upper triangle is stored in AP.
!
    if ( incx == 1 ) then
      do j = 1, n
        if ( x(j) /= 0.0d+00 ) then
          temp = alpha * x(j)
          k = kk
          do i = 1, j
            ap(k) = ap(k) + x(i) * temp
            k = k + 1
          end do
        end if
        kk = kk + j
      end do
    else
      jx = kx
      do j = 1, n
        if ( x(jx) /= 0.0d+00 ) then
          temp = alpha * x(jx)
          ix = kx
          do k = kk, kk + j - 1
            ap(k) = ap(k) + x(ix) * temp
            ix = ix + incx
          end do
        end if
        jx = jx + incx
        kk = kk + j
      end do
    end if
  else
!
!  Form A when lower triangle is stored in AP.
!
    if ( incx == 1 ) then
      do j = 1, n
        if ( x(j) /= 0.0d+00 ) then
          temp = alpha * x(j)
          k = kk
          do i = j, n
            ap(k) = ap(k) + x(i) * temp
            k = k + 1
          end do
        end if
        kk = kk + n - j + 1
      end do
    else
      jx = kx
      do j = 1, n
        if ( x(jx) /= 0.0d+00 ) then
          temp = alpha * x(jx)
          ix = jx
          do k = kk, kk + n - j
            ap(k) = ap(k) + x(ix) * temp
            ix = ix + incx
          end do
        end if
        jx = jx + incx
        kk = kk + n - j + 1
      end do
    end if
  end if

  return
end
subroutine dsymv ( uplo, n, alpha, a, lda, x, incx, beta, y, incy )

!*****************************************************************************80
!
!! DSYMV computes y := alpha*A*x + beta*y for symmetric matrix A.
!
!  Discussion:
!
!    DSYMV performs the matrix-vector operation
!
!      y := alpha*A*x + beta*y,
!
!    where alpha and beta are scalars, x and y are n element vectors and
!    A is an n by n symmetric matrix.
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    Jack Dongarra, Argonne National Lab.
!    Jeremy Du Croz, Nag Central Office.
!    Sven Hammarling, Nag Central Office.
!    Richard Hanson, Sandia National Labs.
!
!  Parameters:
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'U'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'L'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real ( kind = 8 ).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'U', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced.
!           Before entry with UPLO = 'L' or 'L', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - real ( kind = 8 ).
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y. On exit, Y is overwritten by the updated
!           vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer lda

  real ( kind = 8 ) a(lda,*)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer i
  integer incx
  integer incy
  integer info
  integer ix
  integer iy
  integer j
  integer jx
  integer jy
  integer kx
  integer ky
  logical, external :: lsame
  intrinsic max
  integer n
  real ( kind = 8 ) temp1
  real ( kind = 8 ) temp2
  character uplo
  real ( kind = 8 ) x(*)
  external xerbla
  real ( kind = 8 ) y(*)
!
!  Test the input parameters.
!
  info = 0
  if ( .not. lsame ( uplo, 'U' ) .and.  &
       .not. lsame ( uplo, 'L' )      ) then
    info = 1
  else if ( n < 0 ) then
    info = 2
  else if ( lda < max ( 1, n ) ) then
    info = 5
  else if ( incx == 0 ) then
    info = 7
  else if ( incy == 0 ) then
    info = 10
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dsymv ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( ( n == 0 ) .or. ( ( alpha == 0.0d+00 ) .and. ( beta == 1.0d+00 ) ) ) then
    return
  end if
!
!  Set up the start points in X and Y.
!
  if ( 0 < incx ) then
    kx = 1
  else
    kx = 1 - ( n - 1 ) * incx
  end if

  if ( 0 < incy ) then
    ky = 1
  else
    ky = 1 - ( n - 1 ) * incy
  end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through the triangular part
!  of A.
!
!  First form  y := beta*y.
!
  if ( beta /= 1.0d+00 ) then
    if ( incy == 1 ) then
      if ( beta == 0.0d+00 ) then
        do i = 1, n
          y(i) = 0.0d+00
        end do
      else
        do i = 1, n
          y(i) = beta * y(i)
        end do
      end if
    else
      iy = ky
      if ( beta == 0.0d+00 ) then
        do i = 1, n
          y(iy) = 0.0d+00
          iy = iy + incy
        end do
      else
        do i = 1, n
          y(iy) = beta * y(iy)
          iy = iy + incy
        end do
      end if
    end if
  end if

  if ( alpha == 0.0d+00 ) then
   return
  end if

  if ( lsame ( uplo, 'U' ) ) then
!
!  Form  y  when A is stored in upper triangle.
!
    if ( ( incx == 1 ) .and. ( incy == 1 ) ) then
      do j = 1, n
        temp1 = alpha * x(j)
        temp2 = 0.0d+00
        do i = 1, j - 1
          y(i) = y(i) + temp1 * a(i,j)
          temp2 = temp2 + a(i,j) * x(i)
        end do
        y(j) = y(j) + temp1 * a(j,j) + alpha * temp2
      end do
    else
      jx = kx
      jy = ky
      do j = 1, n
        temp1 = alpha * x(jx)
        temp2 = 0.0d+00
        ix = kx
        iy = ky
        do i = 1, j - 1
          y(iy) = y(iy) + temp1 * a(i,j)
          temp2 = temp2 + a(i,j) * x(ix)
          ix = ix + incx
          iy = iy + incy
        end do
        y(jy) = y(jy) + temp1 * a(j,j) + alpha * temp2
        jx = jx + incx
        jy = jy + incy
      end do
    end if
  else
!
!  Form y when A is stored in lower triangle.
!
    if ( ( incx == 1 ).and.( incy == 1 ) ) then
      do j = 1, n
        temp1 = alpha * x(j)
        temp2 = 0.0d+00
        y(j) = y(j) + temp1 * a(j,j)
        do i = j + 1, n
          y(i) = y(i) + temp1 * a(i,j)
          temp2 = temp2 + a(i,j) * x(i)
        end do
        y(j) = y(j) + alpha * temp2
      end do
    else
      jx = kx
      jy = ky
      do j = 1, n
        temp1 = alpha * x(jx)
        temp2 = 0.0d+00
        y(jy) = y(jy) + temp1 * a(j,j)
        ix = jx
        iy = jy
        do i = j + 1, n
          ix = ix + incx
          iy = iy + incy
          y(iy) = y(iy) + temp1 * a(i,j)
          temp2 = temp2 + a(i,j) * x(ix)
        end do
        y(jy) = y(jy) + alpha * temp2
        jx = jx + incx
        jy = jy + incy
      end do
    end if
  end if

  return
end
subroutine dsyr2 ( uplo, n, alpha, x, incx, y, incy, a, lda )

!*****************************************************************************80
!
!! DSYR2 computes A := alpha*x*y' + alpha*y*x' + A, where A is symmetric.
!
!  Discussion:
!
!    DSYR2 performs the symmetric rank 2 operation
!
!      A := alpha*x*y' + alpha*y*x' + A,
!
!    where alpha is a scalar, x and y are n element vectors and A is an n
!    by n symmetric matrix.
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    Jack Dongarra, Argonne National Lab.
!    Jeremy Du Croz, Nag Central Office.
!    Sven Hammarling, Nag Central Office.
!    Richard Hanson, Sandia National Labs.
!
!  Parameters:
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'U'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'L'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real ( kind = 8 ).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'U', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced. On exit, the
!           upper triangular part of the array A is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry with UPLO = 'L' or 'L', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced. On exit, the
!           lower triangular part of the array A is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
  implicit none

  real ( kind = 8 )   alpha
  integer            incx, incy, lda, n
  character       uplo
  real ( kind = 8 )   a( lda, * ), x( * ), y( * )
  real ( kind = 8 )   temp1, temp2
  integer            i, info, ix, iy, j, jx, jy, kx, ky
  logical, external :: lsame
  external           xerbla
  intrinsic          max
!
!  Test the input parameters.
!
  info = 0
  if     ( .not. lsame ( uplo, 'U' ) .and.  &
         .not. lsame ( uplo, 'L' )      ) then
    info = 1
  else if ( n < 0 ) then
    info = 2
  else if ( incx == 0 ) then
    info = 5
  else if ( incy == 0 ) then
    info = 7
  else if ( lda < max ( 1, n ) ) then
    info = 9
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dsyr2 ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( ( n == 0 ) .or. ( alpha == 0.0d+00 ) ) then
    return
  end if
!
!  Set up the start points in X and Y if the increments are not both
!  unity.
!
  if ( ( incx /= 1 ) .or. ( incy /= 1 ) ) then
    if ( 0 < incx ) then
      kx = 1
    else
      kx = 1 - ( n - 1 ) * incx
    end if
    if ( 0 < incy ) then
      ky = 1
    else
      ky = 1 - ( n - 1 ) * incy
    end if
    jx = kx
    jy = ky
  end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through the triangular part
!  of A.
!
  if ( lsame ( uplo, 'U' ) ) then
!
!  Form A when A is stored in the upper triangle.
!
    if ( ( incx == 1 ) .and. ( incy == 1 ) ) then
      do j = 1, n
        if ( ( x(j) /= 0.0d+00 ) .or. ( y(j) /= 0.0d+00 ) ) then
          temp1 = alpha * y(j)
          temp2 = alpha * x(j)
          do i = 1, j
            a(i,j) = a(i,j) + x(i) * temp1 + y(i) * temp2
          end do
        end if
      end do
    else
      do j = 1, n
        if ( ( x(jx) /= 0.0d+00 ) .or. ( y(jy) /= 0.0d+00 ) ) then
          temp1 = alpha * y(jy)
          temp2 = alpha * x(jx)
          ix = kx
          iy = ky
          do i = 1, j
            a(i,j) = a(i,j) + x(ix) * temp1 + y(iy) * temp2
            ix = ix + incx
            iy = iy + incy
          end do
        end if
        jx = jx + incx
        jy = jy + incy
      end do
    end if
  else
!
!  Form A when A is stored in the lower triangle.
!
    if ( ( incx == 1 ) .and. ( incy == 1 ) ) then
      do j = 1, n
        if ( ( x(j) /= 0.0d+00 ) .or. ( y(j) /= 0.0d+00 ) ) then
          temp1 = alpha * y(j)
          temp2 = alpha * x(j)
          do i = j, n
            a(i,j) = a(i,j) + x(i) * temp1 + y(i) * temp2
          end do
        end if
      end do
    else
      do j = 1, n
        if ( ( x(jx) /= 0.0d+00 ) .or. ( y(jy) /= 0.0d+00 ) ) then
          temp1 = alpha * y(jy)
          temp2 = alpha * x(jx)
          ix = jx
          iy = jy
          do i = j, n
            a(i,j) = a(i,j) + x(ix) * temp1 + y(iy) * temp2
            ix = ix + incx
            iy = iy + incy
          end do
        end if
        jx = jx + incx
        jy = jy + incy
      end do
    end if
  end if

  return
end
subroutine dsyr ( uplo, n, alpha, x, incx, a, lda )

!*****************************************************************************80
!
!! DSYR computes A := alpha*x*x' + A where A is a symmetric matrix.
!
!  Discussion:
!
!    DSYR performs the symmetric rank 1 operation
!
!      A := alpha*x*x' + A,
!
!    where alpha is a real scalar, x is an n element vector and A is an
!    n by n symmetric matrix.
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    Jack Dongarra, Argonne National Lab.
!    Jeremy Du Croz, Nag Central Office.
!    Sven Hammarling, Nag Central Office.
!    Richard Hanson, Sandia National Labs.
!
!  Parameters:
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'U'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'L'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real ( kind = 8 ).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'U', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced. On exit, the
!           upper triangular part of the array A is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry with UPLO = 'L' or 'L', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced. On exit, the
!           lower triangular part of the array A is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
  implicit none

  real ( kind = 8 )   alpha
  integer            incx, lda, n
  character       uplo
  real ( kind = 8 )   a( lda, * ), x( * )
  real ( kind = 8 )   temp
  integer            i, info, ix, j, jx, kx
  logical, external :: lsame
  external           xerbla
  intrinsic          max
!
!  Test the input parameters.
!
  info = 0
  if     ( .not. lsame ( uplo, 'U' ) .and.  &
        .not. lsame ( uplo, 'L' )      ) then
    info = 1
  else if ( n < 0 ) then
    info = 2
  else if ( incx == 0 ) then
    info = 5
  else if ( lda < max ( 1, n ) ) then
    info = 7
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dsyr  ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( ( n == 0 ) .or. ( alpha == 0.0d+00 ) ) then
    return
  end if
!
!  Set the start point in X if the increment is not unity.
!
  if ( incx <= 0 ) then
    kx = 1 - ( n - 1 ) * incx
  else if ( incx /= 1 ) then
    kx = 1
  end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through the triangular part
!  of A.
!
  if ( lsame ( uplo, 'U' ) ) then
!
!  Form  A  when A is stored in upper triangle.
!
    if ( incx == 1 ) then
      do j = 1, n
        if ( x(j) /= 0.0d+00 ) then
          temp = alpha * x(j)
          do i = 1, j
            a(i,j) = a(i,j) + x(i) * temp
          end do
        end if
      end do
    else
      jx = kx
      do j = 1, n
        if ( x(jx) /= 0.0d+00 ) then
          temp = alpha * x(jx)
          ix = kx
          do i = 1, j
            a(i,j) = a(i,j) + x(ix) * temp
            ix = ix + incx
          end do
        end if
        jx = jx + incx
      end do
    end if
  else
!
!  Form A when A is stored in lower triangle.
!
    if ( incx == 1 ) then
      do j = 1, n
        if ( x(j) /= 0.0d+00 ) then
          temp = alpha * x(j)
          do i = j, n
            a(i,j) = a(i,j) + x(i) * temp
          end do
        end if
      end do
    else
      jx = kx
      do j = 1, n
        if ( x(jx) /= 0.0d+00 ) then
          temp = alpha * x(jx)
          ix = jx
          do i = j, n
            a(i,j) = a(i,j) + x(ix) * temp
            ix = ix + incx
          end do
        end if
        jx = jx + incx
      end do
    end if
  end if

  return
end
subroutine dtbmv ( uplo, trans, diag, n, k, a, lda, x, incx )

!*****************************************************************************80
!
!! DTBMV computes x = A*x or x = A'*x for a triangular band matrix A.
!
!  Discussion:
!
!    DTBMV performs one of the matrix-vector operations
!
!      x := A*x,   or   x := A'*x,
!
!    where x is an n element vector and  A is an n by n unit, or non-unit,
!    upper or lower triangular band matrix, with ( k + 1 ) diagonals.
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    Jack Dongarra, Argonne National Lab.
!    Jeremy Du Croz, Nag Central Office.
!    Sven Hammarling, Nag Central Office.
!    Richard Hanson, Sandia National Labs.
!
!  Parameters:
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'U'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'L'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'N'   x := A*x.
!
!              TRANS = 'T' or 'T'   x := A'*x.
!
!              TRANS = 'C' or 'C'   x := A'*x.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'U'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'N'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry with UPLO = 'U' or 'U', K specifies the number of
!           super-diagonals of the matrix A.
!           On entry with UPLO = 'L' or 'L', K specifies the number of
!           sub-diagonals of the matrix A.
!           K must satisfy  0 <= K.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, n ).
!           Before entry with UPLO = 'U' or 'U', the leading ( k + 1 )
!           by n part of the array A must contain the upper triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. The top left k by k triangle
!           of the array A is not referenced.
!           The following program segment will transfer an upper
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 do J = 1, N
!                    M = K + 1 - J
!                    do I = max ( 1, J - K ), J
!                       A( M + I, J ) = matrix( I, J )
!                    end do
!                 end do
!
!           Before entry with UPLO = 'L' or 'L', the leading ( k + 1 )
!           by n part of the array A must contain the lower triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row 1 of
!           the array, the first sub-diagonal starting at position 1 in
!           row 2, and so on. The bottom right k by k triangle of the
!           array A is not referenced.
!           The following program segment will transfer a lower
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 do J = 1, N
!                    M = 1 - J
!                    do I = J, min ( N, J + K )
!                       A( M + I, J ) = matrix( I, J )
!                    end do
!                 end do
!
!           Note that when DIAG = 'U' or 'U' the elements of the array A
!           corresponding to the diagonal elements of the matrix are not
!           referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( k + 1 ).
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer            incx, k, lda, n
  character       diag, trans, uplo
  real ( kind = 8 )   a( lda, * ), x( * )
  real ( kind = 8 )   temp
  integer            i, info, ix, j, jx, kplus1, kx, l
  logical            nounit
  logical, external :: lsame
  external           xerbla
  intrinsic          max, min
!
!  Test the input parameters.
!
  info = 0
  if ( .not. lsame ( uplo , 'U' ) .and.  &
       .not. lsame ( uplo , 'L' )      ) then
    info = 1
  else if ( .not. lsame ( trans, 'N' ) .and.  &
            .not. lsame ( trans, 'T' ) .and.  &
            .not. lsame ( trans, 'C' )      ) then
    info = 2
  else if ( .not. lsame ( diag , 'U' ) .and.  &
            .not. lsame ( diag , 'N' )      ) then
    info = 3
  else if ( n < 0 ) then
    info = 4
  else if ( k < 0 ) then
    info = 5
  else if ( lda < ( k + 1 ) ) then
    info = 7
  else if ( incx == 0 ) then
    info = 9
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dtbmv ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( n == 0 ) then
    return
  end if

  nounit = lsame ( diag, 'N' )
!
!  Set up the start point in X if the increment is not unity. This
!  will be ( N - 1 ) * INCX too small for descending loops.
!
  if ( incx <= 0 ) then
    kx = 1 - ( n - 1 ) * incx
  else if ( incx /= 1 ) then
    kx = 1
  end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through A.
!
  if ( lsame ( trans, 'N' ) ) then
!
!  Form x := A*x.
!
    if ( lsame ( uplo, 'U' ) ) then
      kplus1 = k + 1
      if ( incx == 1 ) then
        do j = 1, n
          if ( x(j) /= 0.0d+00 ) then
            temp = x(j)
            l = kplus1 - j
            do i = max ( 1, j - k ), j - 1
              x(i) = x(i) + temp * a(l+i,j)
            end do
            if ( nounit ) then
              x(j) = x(j) * a(kplus1,j)
            end if
          end if
        end do
      else
        jx = kx
        do j = 1, n
          if ( x(jx) /= 0.0d+00 ) then
            temp = x(jx)
            ix = kx
            l = kplus1 - j
            do i = max ( 1, j - k ), j - 1
              x(ix) = x(ix) + temp * a(l+i,j)
              ix = ix + incx
            end do
            if ( nounit ) then
              x(jx) = x(jx) * a(kplus1,j)
            end if
          end if
          jx = jx + incx
          if ( k < j ) then
           kx = kx + incx
          end if
        end do
      end if
    else
      if ( incx == 1 ) then
        do j = n, 1, -1
          if ( x(j) /= 0.0d+00 ) then
            temp = x(j)
            l = 1 - j
            do i = min ( n, j + k ), j + 1, -1
              x(i) = x(i) + temp * a(l+i,j)
            end do
            if ( nounit ) then
              x(j) = x(j) * a(1,j)
            end if
          end if
        end do
      else
        kx = kx + ( n - 1 ) * incx
        jx = kx
        do j = n, 1, -1
          if ( x(jx) /= 0.0d+00 ) then
            temp = x(jx)
            ix = kx
            l = 1 - j
            do i = min ( n, j + k ), j + 1, -1
              x(ix) = x(ix) + temp * a(l+i,j)
              ix = ix - incx
            end do
            if ( nounit ) then
              x(jx) = x(jx) * a(1,j)
            end if
          end if
          jx = jx - incx
          if ( k <= ( n - j ) ) then
            kx = kx - incx
          end if
        end do
      end if
    end if
  else
!
!  Form x := A'*x.
!
    if ( lsame ( uplo, 'U' ) ) then
      kplus1 = k + 1
      if ( incx == 1 ) then
        do j = n, 1, -1
          temp = x(j)
          l = kplus1 - j
          if ( nounit ) then
            temp = temp * a(kplus1,j)
          end if
          do i = j - 1, max ( 1, j - k ), -1
            temp = temp + a(l+i,j) * x(i)
          end do
          x(j) = temp
        end do
      else
        kx = kx + ( n - 1 ) * incx
        jx = kx
        do j = n, 1, -1
          temp = x(jx)
          kx = kx - incx
          ix = kx
          l = kplus1 - j
          if ( nounit ) then
            temp = temp * a(kplus1,j)
          end if
          do i = j - 1, max ( 1, j - k ), -1
            temp = temp + a(l+i,j) * x(ix)
            ix = ix - incx
          end do
          x(jx) = temp
          jx = jx - incx
        end do
      end if
    else
      if ( incx == 1 ) then
        do j = 1, n
          temp = x(j)
          l = 1 - j
          if ( nounit ) then
            temp = temp * a(1,j)
          end if
          do i = j + 1, min ( n, j + k )
            temp = temp + a(l+i,j) * x(i)
          end do
          x(j) = temp
        end do
      else
        jx = kx
        do j = 1, n
          temp = x(jx)
          kx = kx + incx
          ix = kx
          l = 1 - j
          if ( nounit ) then
            temp = temp * a(1,j)
          end if
          do i = j + 1, min ( n, j + k )
            temp = temp + a(l+i,j) * x(ix)
            ix = ix + incx
          end do
          x(jx) = temp
          jx = jx + incx
        end do
      end if
    end if
  end if

  return
end
subroutine dtbsv ( uplo, trans, diag, n, k, a, lda, x, incx )

!*****************************************************************************80
!
!! DTBSV solves A*x = b or A'*x = b for triangular band matrix A.
!
!  Discussion:
!
!    DTBSV solves one of the systems of equations
!
!      A*x = b,   or   A'*x = b,
!
!    where b and x are n element vectors and A is an n by n unit, or
!    non-unit, upper or lower triangular band matrix, with ( k + 1 )
!    diagonals.
!
!    No test for singularity or near-singularity is included in this
!    routine. Such tests must be performed before calling this routine.
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    Jack Dongarra, Argonne National Lab.
!    Jeremy Du Croz, Nag Central Office.
!    Sven Hammarling, Nag Central Office.
!    Richard Hanson, Sandia National Labs.
!
!  Parameters:
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'U'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'L'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the equations to be solved as
!           follows:
!
!              TRANS = 'N' or 'N'   A*x = b.
!
!              TRANS = 'T' or 'T'   A'*x = b.
!
!              TRANS = 'C' or 'C'   A'*x = b.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'U'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'N'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry with UPLO = 'U' or 'U', K specifies the number of
!           super-diagonals of the matrix A.
!           On entry with UPLO = 'L' or 'L', K specifies the number of
!           sub-diagonals of the matrix A.
!           K must satisfy  0 <= K.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, n ).
!           Before entry with UPLO = 'U' or 'U', the leading ( k + 1 )
!           by n part of the array A must contain the upper triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. The top left k by k triangle
!           of the array A is not referenced.
!           The following program segment will transfer an upper
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 do J = 1, N
!                    M = K + 1 - J
!                    do I = max ( 1, J - K ), J
!                       A( M + I, J ) = matrix( I, J )
!                    end do
!                 end do
!
!           Before entry with UPLO = 'L' or 'L', the leading ( k + 1 )
!           by n part of the array A must contain the lower triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row 1 of
!           the array, the first sub-diagonal starting at position 1 in
!           row 2, and so on. The bottom right k by k triangle of the
!           array A is not referenced.
!           The following program segment will transfer a lower
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 do J = 1, N
!                    M = 1 - J
!                    do I = J, min ( N, J + K )
!                       A( M + I, J ) = matrix( I, J )
!                    end do
!                 end do
!
!           Note that when DIAG = 'U' or 'U' the elements of the array A
!           corresponding to the diagonal elements of the matrix are not
!           referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( k + 1 ).
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer            incx, k, lda, n
  character       diag, trans, uplo
  real ( kind = 8 )   a( lda, * ), x( * )
  real ( kind = 8 )   temp
  integer            i, info, ix, j, jx, kplus1, kx, l
  logical            nounit
  logical, external :: lsame
  external           xerbla
  intrinsic          max, min
!
!  Test the input parameters.
!
  info = 0
  if     ( .not. lsame ( uplo , 'U' ) .and.  &
          .not. lsame ( uplo , 'L' )      ) then
    info = 1
  else if ( .not. lsame ( trans, 'N' ) .and.  &
          .not. lsame ( trans, 'T' ) .and.  &
          .not. lsame ( trans, 'C' )      ) then
    info = 2
  else if ( .not. lsame ( diag , 'U' ) .and.  &
         .not. lsame ( diag , 'N' )      ) then
    info = 3
  else if ( n < 0 ) then
    info = 4
  else if ( k < 0 ) then
    info = 5
  else if ( lda < ( k + 1 ) ) then
    info = 7
  else if ( incx == 0 ) then
    info = 9
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dtbsv ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( n == 0 ) then
    return
  end if

  nounit = lsame ( diag, 'N' )
!
!  Set up the start point in X if the increment is not unity. This
!  will be  ( N - 1 ) * INCX  too small for descending loops.
!
  if ( incx <= 0 ) then
    kx = 1 - ( n - 1 ) * incx
  else if ( incx /= 1 ) then
    kx = 1
  end if
!
!  Start the operations. In this version the elements of A are
!  accessed by sequentially with one pass through A.
!
  if ( lsame ( trans, 'N' ) ) then
!
!  Form  x := inv( A )*x.
!
    if ( lsame ( uplo, 'U' ) ) then
      kplus1 = k + 1
      if ( incx == 1 ) then
        do j = n, 1, -1
          if ( x(j) /= 0.0d+00 ) then
            l = kplus1 - j
            if ( nounit ) then
              x(j) = x(j) / a(kplus1,j)
            end if
            temp = x(j)
            do i = j - 1, max ( 1, j - k ), -1
              x(i) = x(i) - temp * a(l+i,j)
            end do
          end if
        end do
      else
        kx = kx + ( n - 1 ) * incx
        jx = kx
        do j = n, 1, -1
          kx = kx - incx
          if ( x(jx) /= 0.0d+00 ) then
            ix = kx
            l = kplus1 - j
            if ( nounit ) then
              x(jx) = x(jx) / a(kplus1,j)
            end if
            temp = x(jx)
            do i = j - 1, max ( 1, j - k ), -1
              x(ix) = x(ix) - temp * a(l+i,j)
              ix = ix - incx
            end do
          end if
          jx = jx - incx
        end do
      end if
    else
      if ( incx == 1 ) then
        do j = 1, n
          if ( x(j) /= 0.0d+00 ) then
            l = 1 - j
            if ( nounit ) then
              x(j) = x(j) / a(1,j)
            end if
            temp = x(j)
            do i = j + 1, min ( n, j + k )
              x(i) = x(i) - temp * a(l+i,j)
            end do
          end if
        end do
      else
        jx = kx
        do j = 1, n
          kx = kx + incx
          if ( x(jx) /= 0.0d+00 ) then
            ix = kx
            l = 1  - j
            if ( nounit ) then
              x(jx) = x(jx) / a(1,j)
            end if
            temp = x(jx)
            do i = j + 1, min ( n, j + k )
              x(ix) = x(ix) - temp * a(l+i,j)
              ix = ix + incx
            end do
          end if
          jx = jx + incx
        end do
      end if
    end if
  else
!
!  Form  x := inv( A')*x.
!
    if ( lsame ( uplo, 'U' ) ) then
      kplus1 = k + 1
      if ( incx == 1 ) then
        do j = 1, n
          temp = x(j)
          l = kplus1 - j
          do i = max ( 1, j - k ), j - 1
            temp = temp - a(l+i,j) * x(i)
          end do
          if ( nounit ) then
            temp = temp / a(kplus1,j)
          end if
          x(j) = temp
        end do
      else
        jx = kx
        do j = 1, n
          temp = x(jx)
          ix = kx
          l = kplus1  - j
          do i = max ( 1, j - k ), j - 1
            temp = temp - a(l+i,j) * x(ix)
            ix = ix + incx
          end do
          if ( nounit ) then
            temp = temp / a(kplus1,j)
          end if
          x(jx) = temp
          jx = jx + incx
          if ( k < j ) then
            kx = kx + incx
          end if
        end do
      end if
    else
      if ( incx == 1 ) then
        do j = n, 1, -1
          temp = x(j)
          l = 1      - j
          do i = min ( n, j + k ), j + 1, -1
            temp = temp - a(l+i,j) * x(i)
          end do
          if ( nounit ) then
            temp = temp / a(1,j)
          end if
          x(j) = temp
        end do
      else
        kx = kx + ( n - 1 ) * incx
        jx = kx
        do j = n, 1, -1
          temp = x(jx)
          ix = kx
          l = 1 - j
          do i = min ( n, j + k ), j + 1, -1
            temp = temp - a(l+i,j) * x(ix)
            ix = ix   - incx
          end do
          if ( nounit ) then
            temp = temp / a(1,j)
          end if
          x(jx) = temp
          jx = jx   - incx
          if ( k <= ( n - j ) ) then
            kx = kx - incx
          end if
        end do
      end if
    end if
  end if

  return
end
subroutine dtpmv ( uplo, trans, diag, n, ap, x, incx )

!*****************************************************************************80
!
!! DTPMV computes x := A*x or x = A'*x for a packed triangular matrix A.
!
!  Discussion:
!
!    DTPMV performs one of the matrix-vector operations
!
!      x := A*x,   or   x := A'*x,
!
!    where x is an n element vector and  A is an n by n unit, or non-unit,
!    upper or lower triangular matrix, supplied in packed form.
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    Jack Dongarra, Argonne National Lab.
!    Jeremy Du Croz, Nag Central Office.
!    Sven Hammarling, Nag Central Office.
!    Richard Hanson, Sandia National Labs.
!
!  Parameters:
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'U'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'L'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'N'   x := A*x.
!
!              TRANS = 'T' or 'T'   x := A'*x.
!
!              TRANS = 'C' or 'C'   x := A'*x.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'U'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'N'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  AP     - real ( kind = 8 ) array of DIMENSION at least
!           ( ( n * ( n + 1 ) ) / 2 ).
!           Before entry with  UPLO = 'U' or 'U', the array AP must
!           contain the upper triangular matrix packed sequentially,
!           column by column, so that AP(1) contains a(1,1),
!           AP(2) and AP(3) contain a(1,2) and a(2,2)
!           respectively, and so on.
!           Before entry with UPLO = 'L' or 'L', the array AP must
!           contain the lower triangular matrix packed sequentially,
!           column by column, so that AP(1) contains a(1,1),
!           AP(2) and AP(3) contain a(2,1) and a(3,1)
!           respectively, and so on.
!           Note that when  DIAG = 'U' or 'U', the diagonal elements of
!           A are not referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
  implicit none

  real ( kind = 8 ) ap( * )
  character diag
  integer i
  integer incx
  integer info
  integer ix
  integer j
  integer jx
  integer k
  integer kk
  integer kx
  logical, external :: lsame
  integer n
  logical nounit
  real ( kind = 8 ) temp
  character trans
  character uplo
  real ( kind = 8 ) x( * )
  external xerbla
!
!  Test the input parameters.
!
  info = 0
  if ( .not. lsame ( uplo , 'U' ) .and.  &
       .not. lsame ( uplo , 'L' )      ) then
    info = 1
  else if ( .not. lsame ( trans, 'N' ) .and.  &
            .not. lsame ( trans, 'T' ) .and.  &
            .not. lsame ( trans, 'C' )      ) then
    info = 2
  else if ( .not. lsame ( diag , 'U' ) .and.  &
            .not. lsame ( diag , 'N' )      ) then
    info = 3
  else if ( n < 0 ) then
    info = 4
  else if ( incx == 0 ) then
    info = 7
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dtpmv ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( n == 0 ) then
    return
  end if

  nounit = lsame ( diag, 'N' )
!
!  Set up the start point in X if the increment is not unity. This
!  will be  ( N - 1 ) * INCX too small for descending loops.
!
  if ( incx <= 0 ) then
    kx = 1 - ( n - 1 ) * incx
  else if ( incx /= 1 ) then
    kx = 1
  end if
!
!  Start the operations. In this version the elements of AP are
!  accessed sequentially with one pass through AP.
!
  if ( lsame ( trans, 'N' ) ) then
!
!  Form x:= A*x.
!
    if ( lsame ( uplo, 'U' ) ) then
      kk = 1
      if ( incx == 1 ) then
        do j = 1, n
          if ( x(j) /= 0.0d+00 ) then
            temp = x(j)
            k = kk
            do i = 1, j - 1
              x(i) = x(i) + temp * ap(k)
              k = k + 1
            end do
            if ( nounit ) then
              x(j) = x(j) * ap(kk+j-1)
            end if
          end if
          kk = kk + j
        end do
      else
        jx = kx
        do j = 1, n
          if ( x(jx) /= 0.0d+00 ) then
            temp = x(jx)
            ix = kx
            do k = kk, kk + j - 2
              x(ix) = x(ix) + temp * ap(k)
              ix = ix + incx
            end do
            if ( nounit ) then
              x(jx) = x(jx) * ap(kk+j-1)
            end if
          end if
          jx = jx + incx
          kk = kk + j
        end do
      end if
    else
      kk = ( n * ( n + 1 ) ) / 2
      if ( incx == 1 ) then
        do j = n, 1, -1
          if ( x(j) /= 0.0d+00 ) then
            temp = x(j)
            k = kk
            do i = n, j + 1, -1
              x(i) = x(i) + temp * ap(k)
              k = k - 1
            end do
            if ( nounit ) then
              x(j) = x(j) * ap(kk-n+j)
            end if
          end if
          kk = kk - ( n - j + 1 )
        end do
      else
        kx = kx + ( n - 1 ) * incx
        jx = kx
        do j = n, 1, -1
          if ( x(jx) /= 0.0d+00 ) then
            temp = x(jx)
            ix = kx
            do k = kk, kk - ( n - ( j + 1 ) ), -1
              x(ix) = x(ix) + temp * ap(k)
              ix = ix - incx
            end do
            if ( nounit ) then
              x(jx) = x(jx) * ap(kk-n+j)
            end if
          end if
          jx = jx - incx
          kk = kk - ( n - j + 1 )
        end do
      end if
    end if
  else
!
!  Form  x := A'*x.
!
    if ( lsame ( uplo, 'U' ) ) then
      kk = ( n * ( n + 1 ) ) / 2
      if ( incx == 1 ) then
        do j = n, 1, -1
          temp = x(j)
          if ( nounit ) then
            temp = temp * ap(kk)
          end if
          k = kk - 1
          do i = j - 1, 1, -1
            temp = temp + ap(k) * x(i)
            k = k    - 1
          end do
          x(j) = temp
          kk = kk - j
        end do
      else
        jx = kx + ( n - 1 ) * incx
        do j = n, 1, -1
          temp = x(jx)
          ix = jx
          if ( nounit ) then
            temp = temp * ap(kk)
          end if
          do k = kk - 1, kk - j + 1, -1
            ix = ix   - incx
            temp = temp + ap(k) * x(ix)
          end do
          x(jx) = temp
          jx = jx - incx
          kk = kk - j
        end do
      end if
    else
      kk = 1
      if ( incx == 1 ) then
        do j = 1, n
          temp = x(j)
          if ( nounit ) then
            temp = temp * ap(kk)
          end if
          k = kk + 1
          do i = j + 1, n
            temp = temp + ap(k) * x(i)
            k = k + 1
          end do
          x(j) = temp
          kk = kk + ( n - j + 1 )
        end do
      else
        jx = kx
        do j = 1, n
          temp = x(jx)
          ix = jx
          if ( nounit ) then
            temp = temp * ap(kk)
          end if
          do k = kk + 1, kk + n - j
            ix = ix + incx
            temp = temp + ap(k) * x(ix)
          end do
          x(jx) = temp
          jx = jx + incx
          kk = kk + ( n - j + 1 )
        end do
      end if
    end if
  end if

  return
end
subroutine dtpsv ( uplo, trans, diag, n, ap, x, incx )

!*****************************************************************************80
!
!! DTPSV solves A*x = b or A'*x = b for a triangular packed matrix A.
!
!  Discussion:
!
!    DTPSV solves one of the systems of equations
!
!      A*x = b,   or   A'*x = b,
!
!    where b and x are n element vectors and A is an n by n unit, or
!    non-unit, upper or lower triangular matrix, supplied in packed form.
!
!    No test for singularity or near-singularity is included in this
!    routine. Such tests must be performed before calling this routine.
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    Jack Dongarra, Argonne National Lab.
!    Jeremy Du Croz, Nag Central Office.
!    Sven Hammarling, Nag Central Office.
!    Richard Hanson, Sandia National Labs.
!
!  Parameters:
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'U'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'L'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the equations to be solved as
!           follows:
!
!              TRANS = 'N' or 'N'   A*x = b.
!
!              TRANS = 'T' or 'T'   A'*x = b.
!
!              TRANS = 'C' or 'C'   A'*x = b.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'U'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'N'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  AP     - real ( kind = 8 ) array of DIMENSION at least
!           ( ( n * ( n + 1 ) ) / 2 ).
!           Before entry with  UPLO = 'U' or 'U', the array AP must
!           contain the upper triangular matrix packed sequentially,
!           column by column, so that AP(1) contains a(1,1),
!           AP(2) and AP(3) contain a(1,2) and a(2,2)
!           respectively, and so on.
!           Before entry with UPLO = 'L' or 'L', the array AP must
!           contain the lower triangular matrix packed sequentially,
!           column by column, so that AP(1) contains a(1,1),
!           AP(2) and AP(3) contain a(2,1) and a(3,1)
!           respectively, and so on.
!           Note that when  DIAG = 'U' or 'U', the diagonal elements of
!           A are not referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer            incx, n
  character       diag, trans, uplo
  real ( kind = 8 )   ap( * ), x( * )
  real ( kind = 8 )   temp
  integer            i, info, ix, j, jx, k, kk, kx
  logical            nounit
  logical, external :: lsame
  external           xerbla
!
!  Test the input parameters.
!
  info = 0
  if     ( .not. lsame ( uplo , 'U' ) .and.  &
          .not. lsame ( uplo , 'L' )      ) then
    info = 1
  else if ( .not. lsame ( trans, 'N' ) .and.  &
          .not. lsame ( trans, 'T' ) .and.  &
          .not. lsame ( trans, 'C' )      ) then
    info = 2
  else if ( .not. lsame ( diag , 'U' ) .and.  &
         .not. lsame ( diag , 'N' )      ) then
    info = 3
  else if ( n < 0 ) then
    info = 4
  else if ( incx == 0 ) then
    info = 7
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dtpsv ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( n == 0 ) then
    return
  end if

  nounit = lsame ( diag, 'N' )
!
!  Set up the start point in X if the increment is not unity. This
!  will be  ( N - 1 ) * INCX  too small for descending loops.
!
  if ( incx <= 0 ) then
    kx = 1 - ( n - 1 ) * incx
  else if ( incx /= 1 ) then
    kx = 1
  end if
!
!  Start the operations. In this version the elements of AP are
!  accessed sequentially with one pass through AP.
!
  if ( lsame ( trans, 'N' ) ) then
!
!  Form  x := inv( A )*x.
!
    if ( lsame ( uplo, 'U' ) ) then
      kk = ( n * ( n + 1 ) ) / 2
      if ( incx == 1 ) then
        do j = n, 1, -1
          if ( x(j) /= 0.0d+00 ) then
            if ( nounit ) then
              x(j) = x(j) / ap(kk)
            end if
            temp = x(j)
            k = kk - 1
            do i = j - 1, 1, -1
              x(i) = x(i) - temp * ap(k)
              k = k - 1
            end do
          end if
          kk = kk - j
        end do
      else
        jx = kx + ( n - 1 ) * incx
        do j = n, 1, -1
          if ( x(jx) /= 0.0d+00 ) then
            if ( nounit ) then
              x(jx) = x(jx) / ap(kk)
            end if
            temp = x(jx)
            ix = jx
            do k = kk - 1, kk - j + 1, -1
              ix = ix - incx
              x(ix) = x(ix) - temp * ap(k)
            end do
          end if
          jx = jx - incx
          kk = kk - j
        end do
      end if
    else
      kk = 1
      if ( incx == 1 ) then
        do j = 1, n
          if ( x(j) /= 0.0d+00 ) then
            if ( nounit ) then
              x(j) = x(j) / ap(kk)
            end if
            temp = x(j)
            k = kk + 1
            do i = j + 1, n
              x(i) = x(i) - temp * ap(k)
              k = k + 1
            end do
          end if
          kk = kk + ( n - j + 1 )
        end do
      else
        jx = kx
        do j = 1, n
          if ( x(jx) /= 0.0d+00 ) then
            if ( nounit ) then
              x(jx) = x(jx) / ap(kk)
            end if
            temp = x(jx)
            ix = jx
            do k = kk + 1, kk + n - j
              ix = ix + incx
              x(ix) = x(ix) - temp * ap(k)
            end do
          end if
          jx = jx + incx
          kk = kk + ( n - j + 1 )
        end do
      end if
    end if
  else
!
!  Form  x := inv( A' )*x.
!
    if ( lsame ( uplo, 'U' ) ) then
      kk = 1
      if ( incx == 1 ) then
        do j = 1, n
          temp = x(j)
          k = kk
          do i = 1, j - 1
            temp = temp - ap(k) * x(i)
            k = k + 1
          end do
          if ( nounit ) then
            temp = temp / ap(kk+j-1)
          end if
          x(j) = temp
          kk = kk + j
        end do
      else
        jx = kx
        do j = 1, n
          temp = x(jx)
          ix = kx
          do k = kk, kk + j - 2
            temp = temp - ap(k) * x(ix)
            ix = ix + incx
          end do
          if ( nounit ) then
            temp = temp / ap(kk+j-1)
          end if
          x(jx) = temp
          jx = jx + incx
          kk = kk + j
        end do
      end if
    else
      kk = ( n * ( n + 1 ) ) / 2
      if ( incx == 1 ) then
        do j = n, 1, -1
          temp = x(j)
          k = kk
          do i = n, j + 1, -1
            temp = temp - ap(k) * x(i)
            k = k - 1
          end do
          if ( nounit ) then
            temp = temp / ap(kk-n+j)
          end if
          x(j) = temp
          kk = kk - ( n - j + 1 )
        end do
      else
        kx = kx + ( n - 1 ) * incx
        jx = kx
        do j = n, 1, -1
          temp = x(jx)
          ix = kx
          do k = kk, kk - ( n - ( j + 1 ) ), -1
            temp = temp - ap(k) * x(ix)
            ix = ix - incx
          end do
          if ( nounit ) then
            temp = temp / ap(kk-n+j)
          end if
          x(jx) = temp
          jx = jx - incx
          kk = kk - ( n - j + 1 )
        end do
      end if
    end if
  end if

  return
end
subroutine dtrmv ( uplo, trans, diag, n, a, lda, x, incx )

!*****************************************************************************80
!
!! DTRMV computes x: = A*x or x = A'*x for a triangular matrix A.
!
!  Discussion:
!
!    DTRMV performs one of the matrix-vector operations
!
!      x := A*x,   or   x := A'*x,
!
!    where x is an n element vector and  A is an n by n unit, or non-unit,
!    upper or lower triangular matrix.
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    Jack Dongarra, Argonne National Lab.
!    Jeremy Du Croz, Nag Central Office.
!    Sven Hammarling, Nag Central Office.
!    Richard Hanson, Sandia National Labs.
!
!  Parameters:
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'U'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'L'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'N'   x := A*x.
!
!              TRANS = 'T' or 'T'   x := A'*x.
!
!              TRANS = 'C' or 'C'   x := A'*x.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'U'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'N'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'U', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           A is not referenced.
!           Before entry with UPLO = 'L' or 'L', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'U', the diagonal elements of
!           A are not referenced either, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer            incx, lda, n
  character       diag, trans, uplo
  real ( kind = 8 )   a( lda, * ), x( * )
  real ( kind = 8 )   temp
  integer            i, info, ix, j, jx, kx
  logical            nounit
  logical, external :: lsame
  external           xerbla
  intrinsic          max
!
!  Test the input parameters.
!
  info = 0
  if  ( .not. lsame ( uplo , 'U' ) .and.  &
        .not. lsame ( uplo , 'L' )      ) then
    info = 1
  else if ( .not. lsame ( trans, 'N' ) .and.  &
            .not. lsame ( trans, 'T' ) .and.  &
            .not. lsame ( trans, 'C' )      ) then
    info = 2
  else if ( .not. lsame ( diag , 'U' ) .and.  &
            .not. lsame ( diag , 'N' )      ) then
    info = 3
  else if ( n < 0 ) then
    info = 4
  else if ( lda < max ( 1, n ) ) then
    info = 6
  else if ( incx == 0 ) then
    info = 8
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dtrmv ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( n == 0 ) then
    return
  end if

  nounit = lsame ( diag, 'N' )
!
!  Set up the start point in X if the increment is not unity. This
!  will be  ( N - 1 ) * INCX  too small for descending loops.
!
  if ( incx <= 0 ) then
    kx = 1 - ( n - 1 ) * incx
  else if ( incx /= 1 ) then
    kx = 1
  end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through A.
!
  if ( lsame ( trans, 'N' ) ) then
!
!  Form x := A*x.
!
    if ( lsame ( uplo, 'U' ) ) then
      if ( incx == 1 ) then
        do j = 1, n
          if ( x(j) /= 0.0d+00 ) then
            temp = x(j)
            do i = 1, j - 1
              x(i) = x(i) + temp * a(i,j)
            end do
            if ( nounit ) then
              x(j) = x(j) * a(j,j)
            end if
          end if
        end do
      else
        jx = kx
        do j = 1, n
          if ( x(jx) /= 0.0d+00 ) then
            temp = x(jx)
            ix = kx
            do i = 1, j - 1
              x(ix) = x(ix) + temp * a(i,j)
              ix = ix + incx
            end do
            if ( nounit ) then
              x(jx) = x(jx) * a(j,j)
            end if
          end if
          jx = jx + incx
        end do
      end if
    else
      if ( incx == 1 ) then
        do j = n, 1, -1
          if ( x(j) /= 0.0d+00 ) then
            temp = x(j)
            do i = n, j + 1, -1
              x(i) = x(i) + temp * a(i,j)
            end do
            if ( nounit ) then
              x(j) = x(j) * a(j,j)
            end if
          end if
        end do
      else
        kx = kx + ( n - 1 ) * incx
        jx = kx
        do j = n, 1, -1
          if ( x(jx) /= 0.0d+00 ) then
            temp = x(jx)
            ix = kx
            do i = n, j + 1, -1
              x(ix) = x(ix) + temp * a(i,j)
              ix = ix - incx
            end do
            if ( nounit ) then
              x(jx) = x(jx) * a(j,j)
            end if
          end if
          jx = jx - incx
        end do
      end if
    end if
  else
!
!  Form x := A'*x.
!
    if ( lsame ( uplo, 'U' ) ) then
      if ( incx == 1 ) then
        do j = n, 1, -1
          temp = x(j)
          if ( nounit ) then
            temp = temp * a(j,j)
          end if
          do i = j - 1, 1, -1
            temp = temp + a(i,j) * x(i)
          end do
          x(j) = temp
        end do
      else
        jx = kx + ( n - 1 ) * incx
        do j = n, 1, -1
          temp = x(jx)
          ix = jx
          if ( nounit ) then
            temp = temp * a(j,j)
          end if
          do i = j - 1, 1, -1
            ix = ix   - incx
            temp = temp + a(i,j) * x(ix)
          end do
          x(jx) = temp
          jx = jx - incx
        end do
      end if
    else
      if ( incx == 1 ) then
        do j = 1, n
          temp = x(j)
          if ( nounit ) then
            temp = temp * a(j,j)
          end if
          do i = j + 1, n
            temp = temp + a(i,j) * x(i)
          end do
          x(j) = temp
        end do
      else
        jx = kx
        do j = 1, n
          temp = x(jx)
          ix = jx
          if ( nounit ) then
            temp = temp * a(j,j)
          end if
          do i = j + 1, n
            ix = ix + incx
            temp = temp + a(i,j) * x(ix)
          end do
          x(jx) = temp
          jx = jx + incx
        end do
      end if
    end if
  end if

  return
end
subroutine dtrsv ( uplo, trans, diag, n, a, lda, x, incx )

!*****************************************************************************80
!
!! DTRSV solves A*x = b or A'*x = b for triangular matrix A.
!
!  Discussion:
!
!    DTRSV solves one of the systems of equations
!
!      A*x = b,   or   A'*x = b,
!
!    where b and x are n element vectors and A is an n by n unit, or
!    non-unit, upper or lower triangular matrix.
!
!    No test for singularity or near-singularity is included in this
!    routine. Such tests must be performed before calling this routine.
!
!  Modified:
!
!    09 June 2005
!
!  Author:
!
!    Jack Dongarra, Argonne National Lab.
!    Jeremy Du Croz, Nag Central Office.
!    Sven Hammarling, Nag Central Office.
!    Richard Hanson, Sandia National Labs.
!
!  Parameters:
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'U'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'L'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the equations to be solved as
!           follows:
!
!              TRANS = 'N' or 'N'   A*x = b.
!
!              TRANS = 'T' or 'T'   A'*x = b.
!
!              TRANS = 'C' or 'C'   A'*x = b.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'U'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'N'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'U', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           A is not referenced.
!           Before entry with UPLO = 'L' or 'L', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'U', the diagonal elements of
!           A are not referenced either, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
  implicit none

  integer            incx, lda, n
  character       diag, trans, uplo
  real ( kind = 8 )   a( lda, * ), x( * )
  real ( kind = 8 )   temp
  integer            i, info, ix, j, jx, kx
  logical            nounit
  logical, external :: lsame
  external           xerbla
  intrinsic          max
!
!  Test the input parameters.
!
  info = 0
  if ( .not. lsame ( uplo , 'U' ) .and.  &
       .not. lsame ( uplo , 'L' )      ) then
    info = 1
  else if ( .not. lsame ( trans, 'N' ) .and.  &
            .not. lsame ( trans, 'T' ) .and.  &
            .not. lsame ( trans, 'C' )      ) then
    info = 2
  else if ( .not. lsame ( diag , 'U' ) .and.  &
            .not. lsame ( diag , 'N' )      ) then
    info = 3
  else if ( n < 0 ) then
    info = 4
  else if ( lda < max ( 1, n ) ) then
    info = 6
  else if ( incx == 0 ) then
    info = 8
  end if

  if ( info /= 0 ) then
    call xerbla ( 'dtrsv ', info )
    return
  end if
!
!  Quick return if possible.
!
  if ( n == 0 ) then
    return
  end if

  nounit = lsame ( diag, 'N' )
!
!  Set up the start point in X if the increment is not unity. This
!  will be  ( N - 1 ) * INCX too small for descending loops.
!
  if ( incx <= 0 ) then
    kx = 1 - ( n - 1 ) * incx
  else if ( incx /= 1 ) then
    kx = 1
  end if
!
!  Start the operations.  In this version the elements of A are
!  accessed sequentially with one pass through A.
!
  if ( lsame ( trans, 'N' ) ) then
!
!  Form x := inv( A )*x.
!
    if ( lsame ( uplo, 'U' ) ) then
      if ( incx == 1 ) then
        do j = n, 1, -1
          if ( x(j) /= 0.0d+00 ) then
            if ( nounit ) then
              x(j) = x(j) / a(j,j)
            end if
            temp = x(j)
            do i = j - 1, 1, -1
              x(i) = x(i) - temp * a(i,j)
            end do
          end if
        end do
      else
        jx = kx + ( n - 1 ) * incx
        do j = n, 1, -1
          if ( x(jx) /= 0.0d+00 ) then
            if ( nounit ) then
              x(jx) = x(jx) / a(j,j)
            end if
            temp = x(jx)
            ix = jx
            do i = j - 1, 1, -1
              ix = ix - incx
              x(ix) = x(ix) - temp * a(i,j)
            end do
          end if
          jx = jx - incx
        end do
      end if
    else
      if ( incx == 1 ) then
        do j = 1, n
          if ( x(j) /= 0.0d+00 ) then
            if ( nounit ) then
              x(j) = x(j) / a(j,j)
            end if
            temp = x(j)
            do i = j + 1, n
              x(i) = x(i) - temp * a(i,j)
            end do
          end if
        end do
      else
        jx = kx
        do j = 1, n
          if ( x(jx) /= 0.0d+00 ) then
            if ( nounit ) then
              x(jx) = x(jx) / a(j,j)
            end if
            temp = x(jx)
            ix = jx
            do i = j + 1, n
              ix = ix + incx
              x(ix) = x(ix) - temp * a(i,j)
            end do
          end if
          jx = jx + incx
        end do
      end if
    end if
  else
!
!  Form  x := inv( A' )*x.
!
    if ( lsame ( uplo, 'U' ) ) then
      if ( incx == 1 ) then
        do j = 1, n
          temp = x(j)
          do i = 1, j - 1
            temp = temp - a(i,j) * x(i)
          end do
          if ( nounit ) then
            temp = temp / a(j,j)
          end if
          x(j) = temp
        end do
      else
        jx = kx
        do j = 1, n
          temp = x(jx)
          ix = kx
          do i = 1, j - 1
            temp = temp - a(i,j) * x(ix)
            ix = ix + incx
          end do
          if ( nounit ) then
            temp = temp / a(j,j)
          end if
          x(jx) = temp
          jx = jx + incx
        end do
      end if
    else
      if ( incx == 1 ) then
        do j = n, 1, -1
          temp = x(j)
          do i = n, j + 1, -1
            temp = temp - a(i,j) * x(i)
          end do
          if ( nounit ) then
            temp = temp / a(j,j)
          end if
          x(j) = temp
        end do
      else
        kx = kx + ( n - 1 ) * incx
        jx = kx
        do j = n, 1, -1
          temp = x(jx)
          ix = kx
          do i = n, j + 1, -1
            temp = temp - a(i,j) * x(ix)
            ix = ix - incx
          end do
          if ( nounit ) then
            temp = temp / a(j,j)
          end if
          x(jx) = temp
          jx = jx - incx
        end do
      end if
    end if
  end if

  return
end
