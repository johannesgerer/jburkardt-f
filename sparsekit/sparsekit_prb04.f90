program main

!*****************************************************************************80
!
!! MAIN is the main program for SPARSEKIT_PRB04.
!
!  test suite for the unary routines.
!  tests some of the routines in the module unary. Still needs to tests
!  many other routines.
!
  integer, parameter :: nxmax = 10
  integer, parameter :: nmx = nxmax * nxmax
  integer, parameter :: nnzmax = 10 * nmx
  integer, parameter :: ndns = 20

  real ( kind = 8 ) a(nnzmax)
  real ( kind = 8 ) a1(nnzmax)
  real ( kind = 8 ) dns(ndns,ndns)
  integer i
  integer ia(nmx+1)
  integer ia1(nnzmax)
  integer iwk(nmx*2+1)
  integer iwork(nnzmax*2)
  integer ja(nnzmax)
  integer ja1(nnzmax)
  integer perm(16)
  integer, dimension ( 16 ) :: qperm = (/ &
    1, 3, 6, 8, 9, 11, 14, 16, 2, 4, 5, 7, 10, 12, 13, 15 /)
  real ( kind = 8 ) stencil(100)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB04'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A set of tests for SPARSEKIT.'
!
!  define correct permutation
!
  do i=1, 16
    perm(qperm(i)) = i
  end do
!
!  dimension of grid
!
      nx = 4
      ny = 4
      nz = 1
!
!  generate grid problem.
!
      call gen57pt(nx,ny,nz,a,ja,ia,iwk,stencil)
      n = nx*ny*nz
!
!  write out the matrix
!
      nnz = ia(n+1)-1
      WRITE(*,*) '-----------------------------------------'
      WRITE(*,*) '  +++  initial matrix in CSR format +++ '
      WRITE(*,*) '-----------------------------------------'
      call dump ( n, a, ja, ia, 6 )
!
!  call csrdns
!
      call csrdns(n,n,a,ja,ia,dns,ndns,ierr)
!
!  write it out as a dense matrix.
!
      WRITE(*,*) '-----------------------------------------'
      WRITE(*,*) '  +++ initial matrix in DENSE format+++ '
      WRITE(*,*) '-----------------------------------------'
      call dmpdns ( n, n, ndns, dns )
!
!  red black ordering
!
      job = 1
      call dperm (n,a,ja,ia,a1,ja1,ia1,perm,perm,job)
      nnz = ia(n+1)-1
      WRITE(*,*) '-----------------------------------------'
      WRITE(*,*) '  +++ red-black matrix in CSR format +++ '
      WRITE(*,*) '-----------------------------------------'
      call dump (n,a1,ja1,ia1,6)
!
!  sort matrix
!
      call csort (n,a1,ja1,ia1,iwork,.true.)
      nnz = ia(n+1)-1
      WRITE(*,*) '-----------------------------------------'
      WRITE(*,*) '  +++     matrix after sorting    +++ '
      WRITE(*,*) '-----------------------------------------'
      call dump (n,a1,ja1,ia1,6)
!
!  convert into dense format
!
       call csrdns(n, n, a1,ja1,ia1,dns,ndns,ierr)
       WRITE(*,*) '-----------------------------------------'
       WRITE(*,*) '  +++ red-black matrix in DENSE format+++ '
       WRITE(*,*) '-----------------------------------------'
       call dmpdns(n,n, ndns, dns)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB04'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine dmpdns ( nrow, ncol, ndns, dns )

!*****************************************************************************80
!
! this subroutine prints out a dense matrix in a simple format.
! the zero elements of the matrix are omitted. The format for the
! nonzero elements is f4.1, i.e., very little precision is provided.
!
! on entry
!
! nrow = row dimension of matrix
! ncol = column dimension of matrix
! ndns = first dimension of array dns.
! dns  = double dimensional array of size n x n containing the matrix
!
! on return
! ---------
! matrix will be printed out.
!
  integer ndns

  real ( kind = 8 ) dns(ndns,*)
  character ( len = 80 ) fmt
  integer i
  integer j
  integer j1
  integer j2
  integer last
  integer ncol
  integer nrow
!
! prints out a dense matrix -- without the zeros.
!
  write ( *, '(4x,16i4)' ) ( j, j = 1, ncol )

      fmt(1:5) = '    |'
      j1 = 6
      do j=1, ncol
        j2 = j1+4
        fmt(j1:j2) = '----'
        j1 = j2
      end do
      last = j1
      fmt(last:last) = '|'
      write ( *, '(a)' ) fmt(1:last)
!
! undo loop 1 ---
!
      j1 = 6
      do j=1,ncol
        j2 = j1+4
        fmt(j1:j2) = '   '
        j1 = j2
      end do

  do i=1, nrow
        j1 = 6
        write (fmt,101) i
 101    format(1x,i2,' |')
        do 3 j=1, ncol
          j2= j1+4
          if ( dns(i,j) /= 0.0 ) then
            write ( fmt(j1:j2), '(f4.1)' ) dns(i,j)
 102        format(f4.1)
            endif
          j1 = j2
 3        continue
        fmt(last:last) = '|'
        write ( *, '(a)' ) fmt(1:last)
  end do

      fmt(1:5) = '    |'
      j1 = 6
      do j=1, ncol
        j2 = j1+4
        fmt(j1:j2) = '----'
        j1 = j2
      end do
      fmt(last:last) = '|'
      write ( *, '(a)' ) fmt(1:last)

  return
end
function afun ( x, y, z )

!*****************************************************************************80
!
  real ( kind = 8 ) afun, x,y, z

  afun = -1.0D+00

  return
end
function bfun (x,y,z)

!*****************************************************************************80
!
  real ( kind = 8 ) bfun, x,y, z
  bfun = -1.0D+00
  return
end
function cfun (x,y,z)

!*****************************************************************************80
!
  real ( kind = 8 ) cfun, x,y, z
  cfun = -1.0D+00
  return
end
function dfun (x,y,z)

!*****************************************************************************80
!
  real ( kind = 8 ) dfun, x,y, z
  data gamma /100.0D+00/
  dfun = 10.0D+00
  return
end
function efun (x,y,z)

!*****************************************************************************80
!
  real ( kind = 8 ) efun, x,y, z
  data gamma /100.0D+00/
  efun = 0.0D+00
  return
end
function ffun (x,y,z)

!*****************************************************************************80
!
  real ( kind = 8 ) ffun
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  ffun = 0.0D+00

  return
end
function gfun (x,y,z)

!*****************************************************************************80
!
  real ( kind = 8 ) gfun, x,y, z
  gfun = 0.0D+00
  return
end
