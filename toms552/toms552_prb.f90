program main

!*****************************************************************************80
!
!! MAIN is the main program for TOMS552_PRB.
!
!  Modified:
!
!    27 December 2007
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS552_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS552 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS552_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 is a sample driver for CL1.
!
!  Discussion:
!
!    this program solves a k by n overdetermined system
!
!      ax=b
!
!    in the l1 sense subject to l equality constraints
!
!      cx=d
!
!    and m inequality constraints
!
!      ex.le.f.
!
!    complete details of the parameters may be
!    found in the documentation of the subroutine.
!
!    the arrays are currently dimensioned to allow problems
!    for which k+l+m .le. 100, n .le. 10.
!
!    the program may be tested on the following data.
!
!     k = 8
!     l = 3
!     m = 2
!     n = 5
!
!      q = 2  0  1  3  1  7
!          7  4  4 15  7  4
!          9  4  7 20  6  7
!          2  2  1  5  3  4
!          9  3  2 14 10  0
!          4  5  0  9  9  4
!          4  4  9 17 -1  9
!          1  6  2  9  5  6
!          0  4  5  9 -1  5
!          3  2  7 12 -2  1
!          3  6 12 21 -3  6
!          0  3  6  9 -3  5
!          6  2  4 12  4  6
!
!     kode = 0
!     toler = 1.e-5
!     iter = 130
!
  dimension cu(2,110)
  integer iu(2,110)
  dimension q(102,12)
  dimension res(100)
  integer s(100)
  dimension x(12)

  data klmd, klm2d, nklmd, n2d /100,102,110,12/
!
! input data.
!
  read (5,99999) k, l, m, n, kode, toler, iter

  klm = k + l + m
  n1 = n + 1
  do i=1,klm
    read (5,99998) (q(i,j),j=1,n1)
     write (6,99994) (q(i,j),j=1,n1)
  end do

  call cl1(k, l, m, n, klmd, klm2d, nklmd, n2d, q, &
  kode, toler, iter, x, res, error, cu, iu, s)
!
! output kode, iteration count and error norm.
!
  write ( *, '(a)' ) ' '
  write (6,99997) kode, iter, error
!
! output solution vector.
!
  write ( *, '(a)' ) ' '
  write (6,99996) (i,x(i),i=1,n)
!
! output residual error at each point.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Residuals:'
  write ( *, '(a)' ) ' '
  write (6,99995) (i,res(i),i=1,klm)
  return
99999 format (5i3, e10.0, i3)
99998 format (8f3.0)
99997 format (16h kode,iter,error, 2i10, e18.7)
99996 format (4h sol, i5, e18.7)
99995 format (6h error, i5, e18.7)
99994 format (2h  , 8f5.0)
end
