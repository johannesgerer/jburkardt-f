program main

!*****************************************************************************80
!
!! MAIN is a test program for exponential propagator using Arnoldi approach.
!
!  Discussion:
!
!    This main program is a very simple test using diagonal matrices
!    (Krylov subspace methods are blind to the structure of the matrix
!    except for symmetry). This provides a good way of testing the
!    accuracy of the method as well as the error estimates.
!
!  Modified:
!
!    02 July 2005
!
  implicit none

  integer, parameter :: nmax = 150
  integer, parameter :: ih0 = 60
  integer, parameter :: ndmx = 20

  integer, parameter :: nzmax = 7 * nmax

  real ( kind = 8 ) a(nzmax)
  real ( kind = 8 ) :: a0 = 0.0
  real ( kind = 8 ) :: b0 = 1.0
  real ( kind = 8 ) ddot
  real ( kind = 8 ) eps
  real ( kind = 8 ) h
  integer ioff(10)
  integer j
  integer k
  integer m
  integer n
  integer ndiag
  real ( kind = 8 ) t
  real ( kind = 8 ) tn
  real ( kind = 8 ) u(ih0*nmax)
  real ( kind = 8 ) w(nmax)
  real ( kind = 8 ) w1(nmax)
  real ( kind = 8 ) x(nmax)
  real ( kind = 8 ) y(nmax)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB12'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test EXPPRO, which computes the matrix exponential.'
!
!  Set dimension of matrix.
!
  n = 100
!
!  Define matrix.
!  A is a single diagonal matrix (ndiag = 1 and ioff(1) = 0 )
!
  ndiag = 1
  ioff(1) = 0
!
!  entries in the diagonal are uniformly distributed.
!
  h = 1.0 / real ( n + 1, kind = 8 )
  do j = 1, n
    a(j) = real ( j + 1, kind = 8 ) / real ( n + 1, kind = 8 )
  end do
!
!  Set a value for TN
!
  tn = 2.0

  eps = 0.0001

  m = 5
  write ( *, '(a,i6)' ) '  Dimension of Krylov subspace M = ', m
!
!  Define initial conditions: chosen so that solution = (1,1,1,1..1)^T
!
  do j = 1, n
    w(j) = exp ( a(j) * tn )
    w1(j) = w(j)
  end do

  call expprod ( n, m, eps, tn, u, w, x, y, a, ioff, ndiag )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First 10 components of final answer:'
  WRITE(*,*)(w(k),k=1,10)
  write ( *, '(a)' ) ' '

  do k = 1, n
    w1(k) = exp ( -a(k) * tn ) * w1(k)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First 10 components of exact solution '
  WRITE(*,*)(w1(k),k=1,10)
!
!  Compute actual 2-norm of error.
!
  t = 0.0
  do k = 1, n
    t = t + ( w1(k) - w(k) )**2
  end do
  t = sqrt ( t / ddot ( n, w, 1, w, 1 ) )

  write ( *, '(a)' ) ' '
  write ( *, * ) '  RMS error (approx-exact)=', t

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB12'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine oped ( n, x, y, diag, ioff, ndiag )

!*****************************************************************************80
!
!! OPED performs a matrix by vector multiplication.
!
!  Discussion:
!
!    The matrix is diagonally structured and stored in diagonal
!    format.
!
  implicit none

  integer n
  integer ndiag

  real ( kind = 8 ) diag(n,ndiag)
  integer i1
  integer i2
  integer io
  integer ioff(ndiag)
  integer j
  integer k
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  y(1:n) = 0.0

  do j = 1, ndiag
    io = ioff(j)
    i1 = max ( 1, 1 - io )
    i2 = min ( n, n - io )
    do k = i1, i2
      y(k) = y(k) + diag(k,j) * x(k+io)
    end do
  end do

  return
end
