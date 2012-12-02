program main

!*****************************************************************************80
!
!! MAIN is the main program for BVLS_PRB.
!
!  Discussion:
!
!    This program demonstrates the use of BVLS for solving least squares
!    problems with include bounds on the variables.
!
!  Modified:
!
!    19 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
  implicit none

  interface

    subroutine bvls ( a, b, bnd, x, rnorm, nsetp, w, index, ierr )
      real ( kind ( 1e0 ) ) a(:,:)
      real ( kind ( 1e0 ) ) b(:)
      real ( kind ( 1e0 ) ) bnd(:,:)
      real ( kind ( 1e0 ) ) x(:)
      real ( kind ( 1e0 ) ) rnorm
      integer nsetp
      real ( kind ( 1e0 ) ) w(:)
      integer index(:)
      integer ierr
    end subroutine

  end interface

  integer, parameter :: mm = 10
  integer, parameter :: nn = 10
  integer, parameter :: mxcase = 6
  integer, parameter :: jstep = 5

  real ( kind(1e0) ) a(mm,nn)
  real ( kind(1e0) ) a2(mm,nn)
  real ( kind(1e0) ) b(mm)
  real ( kind(1e0) ) b2(mm)
  real ( kind(1e0) ) bnd(2,nn)
  real ( kind(1e0) ) bndtab(2,nn,mxcase)
  real ( kind(1e0) ) d(nn)
  integer i
  integer icase
  integer ierr
  integer index(nn)
  integer j
  integer j1
  integer j2
  integer m
  integer, dimension(mxcase) :: mtab = (/ &
    2, 2, 4,  5, 10, 6 /)
  integer n
  integer nsetp
  integer, dimension(mxcase) :: ntab = (/ &
    2, 4, 2, 10,  5, 4 /)
  real ( kind(1e0) ) r(mm)
  real ( kind(1e0) ) rnorm
  real ( kind(1e0) ) rnorm2
  real ( kind(1e0) ) unbnd
  real ( kind(1e0) ) unbtab(mxcase)
  real ( kind(1e0) ) w(nn)
  real ( kind(1e0) ) x(nn)

  data unbtab / 5 * 1.0e6,  999.0e0 /
  data ((bndtab(i,j,1),i=1,2),j=1,2)/ 1.,2.,    3.,4.  /
  data ((bndtab(i,j,2),i=1,2),j=1,4)/ 0,10,  0,10,  0,10,  0,10/
  data ((bndtab(i,j,3),i=1,2),j=1,2)/ 0,100,   -100,100/
  data ((bndtab(i,j,4),i=1,2),j=1,10)/&
         0,0,   -.3994e0,-.3994e0,  -1,1,     -.3e0,-.2e0,    21,22,&
                -4,-3,  45,46,          100,101,  1.e6,1.e6,  -1,1/
  data ((bndtab(i,j,5),i=1,2),j=1,5)/&
                  0,1,  -1,0,  0,1,  .3e0,.4e0,  .048e0,.049e0/
  data ((bndtab(i,j,6),i=1,2),j=1,4)/&
           -100.,100.,  999.,999.,   999.,999.,   999.,999. /

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BVLS_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BVLS library.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Bounded Variables Least Squares.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  If the algorithm succeeds the solution vector, X(),'
  write ( *, '(a)' ) '  and the dual vector, W(), should be related as follows:'
  write ( *, '(a)' ) '    X(i) not at a bound      =>       W(i)  = 0'
  write ( *, '(a)' ) '    X(i) at its lower bound  =>       W(i) <= 0'
  write ( *, '(a)' ) '    X(i) at its upper bound  =>  0 <= W(i)'
  write ( *, '(a)' ) '  except that if an upper bound and lower bound are equal, then'
  write ( *, '(a)' ) '  the corresponding X(i) must take that value and W(i) may have'
  write ( *, '(a)' ) '  any value.'

  do icase = 1, mxcase

    m = mtab(icase)
    n = ntab(icase)
    unbnd = unbtab(icase)

    do j = 1, n
      bnd(1,j) = bndtab(1,j,icase)
      bnd(2,j) = bndtab(2,j,icase)
    end do

	where ( bnd(1,1:n) == unbnd ) bnd(1,1:n) = -huge(1e0)
	where ( bnd(2,1:n) == unbnd ) bnd(2,1:n) =  huge(1e0)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '----------'
    write ( *, '(a,i3)' ) 'Case ', icase
    write ( *, '(a)' ) '----------'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i5,a,i5,a,g17.5)') &
      '  M =', m,',   N =', n,',   UNBND =', unbnd

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Bounds:'
    write ( *, '(a)' ) ' '

    do j1 = 1, n, jstep
      j2 = min ( j1 - 1 + jstep, n )
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g14.6)' ) bnd(1,j1:j2)
      write ( *, '(2x,5g14.6)' ) bnd(2,j1:j2)
    end do

    call random_number ( harvest = b(1:m) )
    call random_number ( harvest = a(1:m,1:n) )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Matrix A:'
    write ( *, '(a)' ) ' '

    do j1 = 1, n, jstep
      j2 = min ( j1 - 1 + jstep, n )
      write ( *, '(a)' ) ' '
      do i = 1,m
        write ( *, '(2x,5g14.6)' ) a(i,j1:j2)
      end do
    end do

    b2(1:m) = b(1:m)
    a2(1:m,1:n) = a(1:m,1:n)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  RHS B:'
    write ( *, '(a)' ) ' '
    write ( *, '(2x,5g14.6)' ) b(1:m)

    call bvls ( a2(1:m,1:n), b2, bnd, x, rnorm, nsetp, w, index, ierr )

    if ( 0 < ierr ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Abnormal error flag, IERR = ', ierr
      stop
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) &
      '  After BVLS:  No. of components not at constraints =', nsetp
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Solution vector, X:'
    write ( *, '(a)' ) ' '
    write ( *, '(2x,5g14.6)' ) x(1:n)

    r(1:m) = b(1:m) - matmul ( a(1:m,1:n), x(1:n) )

    rnorm2 = sqrt ( dot_product ( r(1:m), r(1:m) ) )

    d(1:n) = matmul ( r(1:m), a(1:m,1:n) )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  R = B - A*X Computed by the driver:'
    write ( *, '(a)' ) ' '
    write ( *, '(2x,5g14.6)' ) r(1:m)
    write ( *, '(a)' ) ' '
    write ( *, '(a,g17.5)') '  RNORM2 computed by the driver =', rnorm2
    write ( *, '(a,g17.5)') '  RNORM computed by BVLS       = ', rnorm

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  W = (A**T)*R Computed by the driver:'
    write ( *, '(a)' ) ' '
    write ( *, '(2x,5g14.6)' ) d(1:n)
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Dual vector from BVLS, W() ='
    write ( *, '(a)' ) ' '
    write ( *, '(2x,5g14.6)' ) w(1:n)

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BVLS_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
