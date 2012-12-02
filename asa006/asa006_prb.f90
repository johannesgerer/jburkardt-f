program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA006_PRB.
!
!  Discussion:
!
!    ASA006_PRB calls the ASA006 routines.
!
!  Modified:
!
!    09 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA006_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA006 library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA006_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates the use of CHOLESKY.
!
!  Modified:
!
!    01 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 15

  real ( kind = 8 ) a((n_max*(n_max+1))/2)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nullty
  real ( kind = 8 ) u((n_max*(n_max+1))/2)
  real ( kind = 8 ) ufull(n_max,n_max)
  real ( kind = 8 ) utu

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  CHOLESKY computes the Cholesky factorization'
  write ( *, '(a)' ) '  of a positive definite symmetric matrix.'
  write ( *, '(a)' ) '  A compressed storage format is used.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we look at the matrix A which is'
  write ( *, '(a)' ) '  N+1 on the diagonal and'
  write ( *, '(a)' ) '  N   on the off diagonals.'

  do n = 1, n_max

    nn = ( n * ( n + 1 ) ) / 2
!
!  Set A to the lower triangle of the matrix which is N+1 on the diagonal
!  and N on the off diagonals.
!
    k = 0
    do i = 1, n
      do j = 1, i - 1
        k = k + 1
        a(k) = real ( n, kind = 8 )
      end do
      k = k + 1
      a(k) = real ( n + 1, kind = 8 )
    end do

    call cholesky ( a, n, nn, u, nullty, ifault )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Matrix order N = ', n
    write ( *, '(a,i8)' ) '  Maxtrix nullity NULLTY = ', nullty

    k = 0
    do j = 1, n
      do i = 1, j
        k = k + 1
        ufull(i,j) = u(k)
      end do
      do i = j + 1, n
        ufull(i,j) = 0.0D+00
      end do
    end do
!
!  Compute U' * U and compare to A.
!
    k = 0
    diff = 0.0D+00
    do i = 1, n
      do j = 1, i
        k = k + 1
        utu = 0.0D+00
        do l = 1, n
          utu = utu + ufull(l,i) * ufull(l,j)
        end do
        diff = diff + ( a(k) - utu ) * ( a(k) - utu )
      end do
    end do

    diff = sqrt ( diff )

    write ( *, '(a,g14.6)' ) '  RMS ( A - U''*U ) = ', diff
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 demonstrates the use of CHOLESKY.
!
!  Modified:
!
!    01 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 15

  real ( kind = 8 ) a((n_max*(n_max+1))/2)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nullty
  real ( kind = 8 ) u((n_max*(n_max+1))/2)
  real ( kind = 8 ) ufull(n_max,n_max)
  real ( kind = 8 ) utu

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  CHOLESKY computes the Cholesky factorization'
  write ( *, '(a)' ) '  of a positive definite symmetric matrix.'
  write ( *, '(a)' ) '  A compressed storage format is used.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we look at the Hilbert matrix'
  write ( *, '(a)' ) '  A(I,J) = 1/(I+J-1).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this particular matrix, we expect the'
  write ( *, '(a)' ) '  errors to grow rapidly.'

  do n = 1, n_max

    nn = ( n * ( n + 1 ) ) / 2
!
!  Set A to the Hilbert matrix.
!
    k = 0
    do i = 1, n
      do j = 1, i
        k = k + 1
        a(k) = 1.0D+00 / real ( i + j - 1, kind = 8 )
      end do
    end do

    call cholesky ( a, n, nn, u, nullty, ifault )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Matrix order N = ', n
    write ( *, '(a,i8)' ) '  Maxtrix nullity NULLTY = ', nullty

    k = 0
    do j = 1, n
      do i = 1, j
        k = k + 1
        ufull(i,j) = u(k)
      end do
      do i = j + 1, n
        ufull(i,j) = 0.0D+00
      end do
    end do
!
!  Compute U' * U and compare to A.
!
    k = 0
    diff = 0.0D+00
    do i = 1, n
      do j = 1, i
        k = k + 1
        utu = 0.0D+00
        do l = 1, n
          utu = utu + ufull(l,i) * ufull(l,j)
        end do
        diff = diff + ( a(k) - utu ) * ( a(k) - utu )
      end do
    end do

    diff = sqrt ( diff )

    write ( *, '(a,g14.6)' ) '  RMS ( A - U''*U ) = ', diff
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 demonstrates the use of SUBCHL.
!
!  Modified:
!
!    09 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 15
  integer ( kind = 4 ), parameter :: nn_max = (n_max*(n_max+1))/2

  real ( kind = 8 ) a(nn_max)
  integer ( kind = 4 ) b(n_max)
  real ( kind = 8 ) det
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nullty
  real ( kind = 8 ) u(nn_max)
  real ( kind = 8 ) ufull(n_max,n_max)
  real ( kind = 8 ) utu

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  SUBCHL computes the Cholesky factor'
  write ( *, '(a)' ) '  of a submatrix '
  write ( *, '(a)' ) '  of a positive definite symmetric matrix.'
  write ( *, '(a)' ) '  A compressed storage format is used.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we look at the Hilbert matrix'
  write ( *, '(a)' ) '  A(I,J) = 1/(I+J-1).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this particular matrix, we expect the'
  write ( *, '(a)' ) '  errors to grow rapidly.'
!
!  Set A to the N_MAX order Hilbert matrix.
!
  k = 0
  do i = 1, n_max
    do j = 1, i
      k = k + 1
      a(k) = 1.0D+00 / real ( i + j - 1, kind = 8 )
    end do
  end do

  do n = 1, n_max

    do i = 1, n
      b(i) = i
    end do

    call subchl ( a, b, n, u, nullty, ifault, nn_max, det )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Matrix order N = ', n
    write ( *, '(a,i8)' ) '  Maxtrix nullity NULLTY = ', nullty
    write ( *, '(a,g14.6)' ) '  Matrix determinant DET = ', det

    k = 0
    do j = 1, n
      do i = 1, j
        k = k + 1
        ufull(i,j) = u(k)
      end do
      do i = j + 1, n
        ufull(i,j) = 0.0D+00
      end do
    end do
!
!  Compute U' * U and compare to A.
!
    k = 0
    diff = 0.0D+00
    do i = 1, n
      do j = 1, i
        k = k + 1
        utu = 0.0D+00
        do l = 1, n
          utu = utu + ufull(l,i) * ufull(l,j)
        end do
        diff = diff + ( a(k) - utu )**2
      end do
    end do

    diff = sqrt ( diff )

    write ( *, '(a,g14.6)' ) '  RMS ( A - U''*U ) = ', diff
  end do

  return
end
