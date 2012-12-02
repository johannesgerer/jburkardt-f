program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA007_PRB.
!
!  Discussion:
!
!    ASA007_PRB calls the ASA007 routines.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA007_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA007 library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA007_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates the use of SYMINV.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 15

  real ( kind = 8 ) a((n_max*(n_max+1))/2)
  real ( kind = 8 ) afull(n_max,n_max)
  real ( kind = 8 ) c((n_max*(n_max+1))/2)
  real ( kind = 8 ) cfull(n_max,n_max)
  real ( kind = 8 ) cta
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nullty
  real ( kind = 8 ) w(n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  SYMINV computes the inverse of a positive'
  write ( *, '(a)' ) '  definite symmetric matrix.'
  write ( *, '(a)' ) '  A compressed storage format is used.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we look at the matrix A which is'
  write ( *, '(a)' ) '  N+1 on the diagonal and'
  write ( *, '(a)' ) '  N   on the off diagonals.'

  do n = 1, n_max
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

    call syminv ( a, n, c, w, nullty, ifault )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Matrix order N = ', n
    write ( *, '(a,i8)' ) '  Maxtrix nullity NULLTY = ', nullty

    k = 0
    do j = 1, n
      do i = 1, j - 1
        k = k + 1
        cfull(i,j) = c(k)
        cfull(j,i) = c(k)
      end do
      k = k + 1
      cfull(j,j) = c(k)
    end do

    k = 0
    do j = 1, n
      do i = 1, j - 1
        k = k + 1
        afull(i,j) = a(k)
        afull(j,i) = a(k)
      end do
      k = k + 1
      afull(j,j) = a(k)
    end do
!
!  Compute C * A - I.
!
    diff = 0.0D+00
    do i = 1, n
      do j = 1, n
        cta = 0.0D+00
        do k = 1, n
          cta = cta + cfull(i,k) * afull(k,j)
        end do
        if ( i .eq. j ) then
          diff = diff + ( 1.0D+00 - cta )**2
        else
          diff = diff + cta**2
        end if
      end do
    end do

    diff = sqrt ( diff )

    write ( *, '(a,g14.6)' ) '  RMS ( C * A - I ) = ', diff

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 demonstrates the use of SYMINV.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 15

  real ( kind = 8 ) a((n_max*(n_max+1))/2)
  real ( kind = 8 ) afull(n_max,n_max)
  real ( kind = 8 ) c((n_max*(n_max+1))/2)
  real ( kind = 8 ) cfull(n_max,n_max)
  real ( kind = 8 ) cta
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nullty
  real ( kind = 8 ) w(n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  SYMINV computes the inverse of a positive'
  write ( *, '(a)' ) '  definite symmetric matrix.'
  write ( *, '(a)' ) '  A compressed storage format is used.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we look at the Hilbert matrix'
  write ( *, '(a)' ) '  A(I,J) = 1/(I+J-1).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this particular matrix, we expect the'
  write ( *, '(a)' ) '  errors to grow rapidly.'

  do n = 1, n_max
!
!  Set A to the lower triangle of the matrix which is N+1 on the diagonal
!  and N on the off diagonals.
!
    k = 0
    do i = 1, n
      do j = 1, i
        k = k + 1
        a(k) = 1.0D+00 / real ( i + j - 1, kind = 8 )
      end do
    end do

    call syminv ( a, n, c, w, nullty, ifault )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Matrix order N = ', n
    write ( *, '(a,i8)' ) '  Maxtrix nullity NULLTY = ', nullty

    k = 0
    do j = 1, n
      do i = 1, j - 1
        k = k + 1
        cfull(i,j) = c(k)
        cfull(j,i) = c(k)
      end do
      k = k + 1
      cfull(j,j) = c(k)
    end do

    k = 0
    do j = 1, n
      do i = 1, j - 1
        k = k + 1
        afull(i,j) = a(k)
        afull(j,i) = a(k)
      end do
      k = k + 1
      afull(j,j) = a(k)
    end do
!
!  Compute C * A - I.
!
    diff = 0.0D+00
    do i = 1, n
      do j = 1, n
        cta = 0.0D+00
        do k = 1, n
          cta = cta + cfull(i,k) * afull(k,j)
        end do
        if ( i .eq. j ) then
          diff = diff + ( 1.0D+00 - cta )**2
        else
          diff = diff + cta**2
        end if
      end do
    end do

    diff = sqrt ( diff )

    write ( *, '(a,g14.6)' ) '  RMS ( C * A - I ) = ', diff

  end do

  return
end
