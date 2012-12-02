program main

!*****************************************************************************80
!
! main program for generating BLOCK 5 point and 7-point matrices in the
! Harwell-Boeing format.  Creates a file containing a
! harwell-boeing matrix.
!
  implicit none

  integer, parameter :: nxmax = 30
  integer, parameter :: nmx = nxmax * nxmax

  real ( kind = 8 ) a(9,5*nmx)
  real ( kind = 8 ) ao(45*nmx)
  character ( len = 2 ) guesol
  integer ia(nmx)
  integer iao(5*nmx)
  integer iau(nmx)
  integer ifmt
  integer iout
  integer ja(7*nmx)
  integer jao(15*nmx)
  integer job
  character ( len = 8 ) key
  character ( len = 50 ) matfil
  integer n
  integer na
  integer nfree
  integer nfree2
  integer nr
  integer nx
  integer ny
  integer nz
  real ( kind = 8 ) rhs(1)
  real ( kind = 8 ) stencil(7,100)
  character ( len = 72 ) title
  character ( len = 3 ) type

  call timestamp ( )

  WRITE(*,*)' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB10'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  WRITE(*,*)'This program demonstrates the use of GEN57BL'
  WRITE(*,*)'to generate a block sparse matrix derived from a 5 or'
  WRITE(*,*)'7 point stencil on an NX by NY by NZ grid.'
  WRITE(*,*)' '
  WRITE(*,*)'In this example, the block size is 3.'
  WRITE(*,*)' '
  WRITE(*,*)'The matrix is then stored in Harwell-Boeing format'
  WRITE(*,*)'in a file, using routine PRTMT.'
  WRITE(*,*)' '
!
!  NFREE is the block size.
!
  nx = 10
  ny = 10
  nz = 1
  nfree = 3
  na = 9
  call gen57bl ( nx, ny, nz, nfree, na, n, a, ja, ia, iau, stencil )
! 
!  Convert from BSR (block sparse row) to CSR (column sparse row) format.
!
  nfree2 = nfree * nfree
  nr = n / nfree
  call bsrcsr ( n, nfree, na, a, ja, ia, ao, jao, iao )
!
!  Set other data needed for Harwell Boeing format.
!
  guesol = 'NN'
  title = ' BLOCK 5-POINT TEST MATRIX FROM SPARSKIT               '
  type = 'RUA'
  key = 'BLOCK5PT'
  ifmt = 104
  job = 2
!
!  Store matrix in file, using Harwell Boeing format.
!
  matfil = 'test.mat'
  iout = 7

  open ( unit = iout, file = matfil, status = 'replace' )

  call prtmt ( n, n, ao, jao, iao, rhs, guesol, title, key, type, &
    ifmt, job, iout )

  close ( unit = iout )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB10'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine afunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x, y, z, coeff(100)

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0
    end do
    coeff((j-1)*nfree+j) = -1.0
  end do

  return
end
subroutine bfunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x, y, z, coeff(100)

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0
    end do
    coeff((j-1)*nfree+j) = -1.0
  end do

  return
end
subroutine cfunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  integer i
  integer j
  integer nfree
  real ( kind = 8 )x, y, z, coeff(100)

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0
    end do
    coeff((j-1)*nfree+j) = -1.0
  end do

  return
end
subroutine dfunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  integer i
  integer j
  integer nfree
  real ( kind = 8 )x, y, z, coeff(100)

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0
    end do
  end do

  return
end
subroutine efunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  integer i
  integer j
  integer nfree
  real ( kind = 8 )x, y, z, coeff(100)

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0
    end do
  end do

  return
end
subroutine ffunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x, y, z, coeff(100)

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0
    end do
  end do

  return
end
subroutine gfunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x, y, z, coeff(100)

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0
    end do
  end do

  return
end
