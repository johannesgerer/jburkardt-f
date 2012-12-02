program main

!*****************************************************************************80
!
!! MAIN tests all the subroutines in matvec.
!
!  Discussion:
!
!    This test generates matrices and transforms them in appropriate formats
!    and then calls the appropriate routines.
!
!  NMAX is the maximum number of rows in the matrix, plus 1.
!
  implicit none

  integer, parameter :: nmax=10*10*10+1
!
!  NZMAX is the maximum number of nonzero entries in the matrix.
!  There are at most 7 nonzeroes in the 3-D Laplacian operator.
!
  integer, parameter :: nzmax = 7 * nmax

  real ( kind = 8 ) a1(nzmax)
  real ( kind = 8 ) a2(nzmax)
  integer ia1(nmax)
  integer ia2(nmax)
  integer idiag
  integer, dimension ( 2 ) :: idim = (/ 4, 10 /)
  integer ierr
  integer ii
  integer ioff(10)
  integer iwk1(nmax)
  integer iwk2(nmax)
  integer j
  integer ja1(nzmax)
  integer ja2(nzmax)
  integer jad(nzmax)
  integer jdiag
  integer jj
  integer k
  integer n
  integer na
  integer ndiag
  integer nfree
  integer nlev
  integer nx
  integer ny
  integer nz
  real ( kind = 8 ) scal
  real ( kind = 8 ) stencil(100)
  real ( kind = 8 ) x(nmax)
  real ( kind = 8 ) y(nmax)
  real ( kind = 8 ) y0(nmax)
  real ( kind = 8 ) y1(nmax)

  call timestamp ( )

  WRITE(*,*)' '
  WRITE(*,*)'SPARSEKIT_PRB02'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  WRITE(*,*)'A set of tests for SPARSKIT.'
  WRITE(*,*)'These test the conversion routines, and the'
  WRITE(*,*)'lower and upper triangular solution routines.'
  WRITE(*,*)' '
!
!  ii loop controls the size of the grid.
!
  do ii = 1, 2

    WRITE(*,*) '---------------- ii ',ii,'--------------------'
    nfree = 1
    NA=NFREE*NFREE
    nx = idim(ii)
    ny = nx
!
!  jj loop corresponds to 2-D and 3-D problems.
!
    do jj=1, 2

           WRITE(*,*) '     ----------- jj ',jj,' -------------'
           nz = 1
           if (jj .eq. 2) nz = 10
!
!  call matrix generation routine to create matrix in CSR format.
!
            WRITE(*,*)'nx=',nx
            WRITE(*,*)'ny=',ny
            WRITE(*,*)'nz=',nz
            call gen57bl(nx,ny,nz,nfree,na,n,a1,ja1,ia1,ia2,stencil)
            WRITE(*,*)'The solution vector will have length N=',N
!
!  initialize the solution vector x
!
            do j=1, n
               x(j) = real ( j, kind = 8 )
            end do
!
!  Get exact answer in y0 for later checks.
!
            call amux(n,x,y0,a1,ja1,ia1)
!
!  Convert CSR format to ELL format, and test AMUXE.
!
            call csrell(n,a1,ja1,ia1,7,a2,jad,n,ndiag,ierr)
            call amuxe(n, x, y, n, ndiag, a2,jad)
            call errpr(n, y, y0,'amuxe ')
!
!  Convert CSR format to DIA format, and test AMUXD.
!
            idiag = 7
            call csrdia (n, idiag,10,a1, ja1, ia1, nmax, a2, &
                 ioff, a2, ja2, ia2, jad)
            call amuxd (n,x,y,a2,nmax,idiag,ioff)
            call errpr (n, y, y0,'amuxd ')
!
!  Convert CSR format to CSC format, and test ATMUX.
!
            call csrcsc(n,1,1,a1,ja1,ia1,a2,ja2,ia2)
            call atmux(n,x,y,a2,ja2,ia2)
            call errpr(n, y, y0,'atmux ')
!
!  Convert CSR format to JAD format, and test AMUXJ.
!
            call csrjad (n,a1,ja1,ia1, jdiag, jad, a2, ja2, ia2)
            call amuxj (n, x, y, jdiag, a2, ja2, ia2)
            call dvperm (n, y, jad)
            call errpr (n, y, y0,'amuxj ')
!
!  convert JAD format back to CSR, and test JADCSR and AMUX.
!
            call jadcsr (n, jdiag, a2, ja2, ia2, jad, a1, ja1, ia1)
            call amux (n, x, y, a1, ja1, ia1)
            call errpr (n, y, y0,'jadcsr')
!
!  triangular systems solutions
!
! TESTING LDSOL
!
            call getl (n, a1, ja1, ia1, a2, ja2, ia2)
            call amux (n,x,y0, a2, ja2, ia2)
            call atmux(n,x,y1, a2, ja2, ia2)
            call csrmsr (n, a2, ja2, ia2, a2, ja2, y, iwk2)
            do k=1,n
               a2(k) = 1.0/ a2(k)
            end do
            call ldsol (n, y, y0, a2, ja2)
            call errpr (n, x, y, 'ldsol ')
!
! TESTING LDSOLL
!
            call levels (n, ja2, ja2, nlev, jad, iwk1, iwk2)
            call ldsoll (n, y, y0, a2, ja2, nlev, jad, iwk1)
            call errpr (n, x, y,'ldsoll')
!
! TESTING UDSOLC
!
! here we take advantage of the fact that the MSR format for U
! is the MSC format for L
!
            call udsolc (n, y, y1, a2, ja2)
            call errpr (n, x, y,'udsolc')
!
! TESTING LSOL
!
! here we exploit the fact that with MSR format a, ja, ja is actually
! the correct data structure for the strict lower triangular part of
! the CSR format. First rescale matrix.
!
            scal = 0.1
            do k=ja2(1), ja2(n+1)-1
               a2(k)=a2(k)*scal
            end do
            call amux(n, x, y0, a2, ja2, ja2)
            do j=1,n
               y0(j) = x(j) + y0(j)
            end do
            call lsol (n, y, y0, a2, ja2, ja2)
            call errpr (n, x, y,'lsol  ')
!
! TESTING UDSOL
!
            call getu (n, a1, ja1, ia1, a2, ja2, ia2)
            call amux (n,x,y0, a2, ja2, ia2)
            call atmux(n,x,y1, a2, ja2, ia2)
            call csrmsr (n, a2, ja2, ia2, a2, ja2, y, jad)
            do k=1,n
               a2(k) = 1.0/ a2(k)
            end do
            call udsol (n, y, y0, a2, ja2)
            call errpr (n, x, y,'udsol ')
!
! TESTING LDSOLC
!
! here we take advantage of the fact that the MSR format for L
! is the MSC format for U
!
            call ldsolc (n, y, y1, a2, ja2)
            call errpr (n, x, y,'ldsolc')
!
! TESTING USOL
!
! here we exploit the fact that with MSR format a, ja, ja is actually
! the correct data structure for the strict lower triangular part of
! the CSR format. First rescale matrix.
!
            scal = 0.1
            do k=ja2(1), ja2(n+1)-1
               a2(k)=a2(k)*scal
            end do
            call amux(n, x, y1, a2, ja2, ja2)
            do j=1,n
               y1(j) = x(j) + y1(j)
            end do
            call usol (n, y, y1, a2, ja2, ja2)
            call errpr (n, x, y,'usol  ')
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB02'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine afunbl ( nfree, x, y, z, coeff )

!*****************************************************************************80
!
!! AFUNBL ???
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
  real ( kind = 8 ) x, y, z, coeff(100)

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
  real ( kind = 8 ) x, y, z, coeff(100)

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
  real ( kind = 8 ) x, y, z, coeff(100)

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
function afun ()

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) afun

  afun = 0.0

  return
end
function bfun ()

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) bfun

  bfun = 0.0

  return
end
function cfun ()

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) cfun

  cfun = 0.0

  return
end
function dfun ( )

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) dfun

  dfun = 0.0

  return
end
function efun ( )

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) efun

  efun = 0.0

  return
end
function ffun ( )

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) ffun

  ffun = 0.0

  return
end
function gfun ( )

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) gfun

  gfun = 0.0

  return
end
subroutine errpr ( n, y, y1, msg )

!*****************************************************************************80
!
  implicit none

  integer n

  integer k
  real ( kind = 8 ) y(n), y1(n), t
  character*6 msg

  t = 0.0
  do k=1,n
    t = t+(y(k)-y1(k))**2
  end do

  t = sqrt(t)
  WRITE(*,*)'RMS error in ',msg,' =', t

  return
end
