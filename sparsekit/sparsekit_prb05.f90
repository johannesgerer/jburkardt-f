program main

!*****************************************************************************80
!
!! MAIN is a test suite for Part I : I/O routines.
!
! tests the following : gen5pt.f, prtmt, readmt, amd pltmt.
! 1) generates a 100 x 100 5pt matrix,
! 2) prints it with a given format in file 'first.mat'
! 3) reads the matrix from 'first.mat' using readmat
! 4) prints it again in file 'second.mat' in  a different format
! 5) makes 4 pic files to show the different options of pltmt.
!    these are in job0.pic, job01.pic, job10.pic, job11.pic 
!                          coded by Y. Saad, RIACS, 08/31/1989.
!
  implicit none

  integer, parameter :: nxmax = 20
  integer, parameter :: nmx = nxmax * nxmax

  real ( kind = 8 ) a(7*nmx)
  character ( len = 2 ) guesol
  integer i
  integer ia(nmx)
  integer iau(nmx)
  integer ierr
  integer ifmt
  integer iout
  integer j
  integer ja(7*nmx)
  integer job
  integer k
  character ( len = 8 ) key
  integer mode
  integer n
  integer ncol
  integer nmax
  integer nnz
  integer nrhs
  integer nrow
  integer nx
  integer ny
  integer nz
  integer nzmax
  real ( kind = 8 ) rhs(3*nmx)
  real ( kind = 8 ) stencil(7)
  character ( len = 72 ) title
  character ( len = 3 ) type

  call timestamp ( )

  WRITE(*,*)' '
  WRITE(*,*)'SPARSEKIT_PRB05'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  WRITE(*,*)'A set of tests for SPARSKIT'
  WRITE(*,*)' '

  open (unit=7,file='first.mat',STATUS='replace')
  open (unit=8,file='second.mat',STATUS='replace')
  open (unit=20,file='job00.pic',STATUS='replace')
  open (unit=21,file='job01.pic',STATUS='replace')
  open (unit=22,file='job10.pic',STATUS='replace')
  open (unit=23,file='job11.pic',STATUS='replace')
!
!  dimension of grid
!
  nx = 10
  ny = 10
  nz = 1
!
!  generate grid problem.
!
  call gen57pt (nx,ny,nz,a,ja,ia,iau,stencil)
!
!  create the Harwell-Boeing matrix. Start by defining title,
!  and type. them define format and print it.
!
  write (title,9) nx, ny
 9    format('Five-point matrix on a square region', &
    ' using a ',I2,' by ',I2,' grid *SPARSKIT*')
  key = 'Fivept10'
  type= 'RSA'
  n = nx*ny*nz
  ifmt = 5
  job = 3
  guesol = 'GX'
!
!  define a right hand side of ones, an initial guess of two's
!  and an exact solution of three's.
!
  do k=1, 3*n
    rhs(k) = real( 1 +  (k-1)/n )
  end do

  call prtmt (n,n,a,ja,ia,rhs,guesol,title,key,type,ifmt,job,7)
!
!  read it again in same matrix a, ja, ia
!
  nmax = nmx
  nzmax = 7*nmx
  do k=1, 3*n
    rhs(k) = 0.0
  end do
  job = 3

  rewind 7
  nrhs = 3*n

  call readmt (nmax,nzmax,job,7,a,ja,ia,rhs,nrhs,guesol, &
               nrow,ncol,nnz,title,key,type,ierr)
  WRITE(*,*)'ierr = ',ierr,' nrhs ', nrhs
!
!  matrix read.  print it again in a different format
!
  ifmt = 102
  ncol = nrow
  job = 3

  call prtmt (nrow,ncol,a,ja,ia,rhs,guesol,title,key,type, &
                          ifmt,job,8)
!
!  print four pic files
!
  mode = 0
  do i=1, 2
    do j=1, 2
      job = (i-1)*10 +j-1
      iout = 20+(i-1)*2+j-1
      call pltmt(nrow,ncol,mode,ja,ia,title,key,type,job,iout)
    end do
  end do

  CLOSE(UNIT=7)
  CLOSE(UNIT=8)
  CLOSE(UNIT=20)
  CLOSE(UNIT=21)
  CLOSE(UNIT=22)
  CLOSE(UNIT=23)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB05'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function afun (x,y,z)

!*****************************************************************************80
!
!! AFUN ???
!
  implicit none

  real ( kind = 8 ) afun, x,y, z

  afun = -1.0

  return
end
function bfun (x,y,z)

!*****************************************************************************80
!
!! BFUN ???
!
  implicit none

  real ( kind = 8 ) bfun, x,y, z

  bfun = -1.0

  return
end
function cfun (x,y,z)

!*****************************************************************************80
!
!! CFUN ???
!
  implicit none

  real ( kind = 8 ) cfun, x,y, z

  cfun = -1.0

  return
end
function dfun (x,y,z)

!*****************************************************************************80
!
!! DFUN ???
!
  implicit none

  real ( kind = 8 ) dfun, gamma, x,y, z

  data gamma /100.0/

  dfun = 10.0

  return
end
function efun (x,y,z)

!*****************************************************************************80
!
!! EFUN ???
!
  implicit none

  real ( kind = 8 ) efun, gamma, x,y, z

  data gamma /100.0/

  efun = 0.0

  return
end
function ffun (x,y,z)

!*****************************************************************************80
!
!! FFUN ???
!
  implicit none

  real ( kind = 8 ) ffun, x,y, z

  ffun = 0.0

  return
end
function gfun (x,y,z)

!*****************************************************************************80
!
!! GFUN ???
!
  implicit none

  real ( kind = 8 ) gfun, x,y, z

  gfun = 0.0

  return
end
