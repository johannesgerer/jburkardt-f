program main

!*****************************************************************************80
!
!! MAIN is the main program for SPARSEKIT_PRB01.
!
!  Modified:
!
!    17 August 2008
!
  implicit none

  integer, parameter :: nxmax = 30
  integer, parameter :: nmx = nxmax * nxmax
  integer, parameter :: nzmax = 7 * nmx

  real ( kind = 8 ) a(nzmax)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(nzmax)
  real ( kind = 8 ) beta
  real ( kind = 8 ) c(nzmax)
  integer ia(nmx+1)
  integer ib(nmx+1)
  integer ic(nzmax)
  integer ierr
  integer ifmt
  integer iw(nmx)
  integer j
  integer ja(nzmax)
  integer jb(nzmax)
  integer jc(nzmax)
  integer jj
  integer job
  integer k
  character ( len = 8 ) key
  integer n
  integer na
  integer nx
  integer ny
  integer nz
  real ( kind = 8 ) s
  real ( kind = 8 ) stencil(100)
  character ( len = 71 ) title
  character ( len = 3 ) type
  real ( kind = 8 ) x(nmx)
  real ( kind = 8 ) y(nmx)
  real ( kind = 8 ) y1(nmx)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB01'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A test problem for SPARSEKIT.'

  nx=10
  ny=10
  nz = 1
  na = 5*nmx
  call gen57pt(nx,ny,nz,a,ja,ia,iw,stencil)
  beta  = alpha
  alpha = 0.0D+00

  call gen57pt(ny,nx,nz,b,jb,ib,iw,stencil)
  n = nx*ny*nz

  s = 3.812D+00

  call aplsb1(n,n,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,ierr)
  if (ierr .ne. 0)WRITE(*,*)' ierr = ',ierr

  call dump(n,c,jc,ic,6)

  do k=1,n
    x(k) = real ( k, kind = 8 ) / real ( n, kind = 8 )
  end do

  call ope (n,x,y1,a,ja,ia)
  call ope (n,x,y,b,jb,ib)
  do j=1, n
    y1(j) = s*y(j) + y1(j)
  end do

  call ope (n,x,y,c,jc,ic)

  write (*,*) ' ------------ checking APLSB --------------'
  call ydfnorm(n,y1,y)

  type = '--------'
  title=' test matrix for blassm c = a+b '
  key = 'rua'

  ifmt = 103

  job = -1

  do jj=1,2
    call apmbt(n,n,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
    if (ierr .ne. 0) WRITE(*,*)' ierr = ',ierr
    call ope (n,x,y1,a,ja,ia)
    call opet (n,x,y,b,jb,ib)
    s = real ( job, kind = 8 )

    y1(1:n) = y1(1:n) + s * y(1:n)

    call ope (n,x,y,c,jc,ic)

    WRITE(*,*) '  '
    WRITE(*,*) ' ------------ checking APMBT---------------'
    WRITE(*,*) ' ------------ with JOB = ',job,' -------------'
    call ydfnorm(n,y1,y)

    job = job + 2
  end do

  type = '--------'
  title=' test matrix for blassm c = a+b^T '

  s = 0.1232445D+00
  call aplsbt(n,n,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,iw,ierr)

  if (ierr .ne. 0) WRITE(*,*)' ierr = ',ierr
  call ope (n,x,y1,a,ja,ia)
  call opet (n,x,y,b,jb,ib)

  y1(1:n) = y1(1:n) + s * y(1:n)

  call ope (n,x,y,c,jc,ic)

  WRITE(*,*) '  '
  WRITE(*,*) ' ------------ checking APLSBT---------------'
  call ydfnorm(n,y1,y)
!
!  Testing matrix products
!
  job = 1
  call amub (n,n,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)

  if (ierr .ne. 0) WRITE(*,*)' ierr = ',ierr
  call ope(n,x,y,b,jb,ib)
  call ope(n,y,y1,a,ja,ia)

  call ope(n,x,y,c,jc,ic)

  WRITE(*,*) '  '
  WRITE(*,*) ' ------------ checking AMUB  ---------------'
  call ydfnorm(n,y1,y)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB01'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function afun ( x, y, z )

!*****************************************************************************80
!
!! AFUN
!
  implicit none

  real ( kind = 8 ) afun
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  afun = -1.0D+00

  return
end
function bfun ( x, y, z )

!*****************************************************************************80
!
!! BFUN
!
  implicit none

  real ( kind = 8 ) bfun
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  bfun = -1.0D+00

  return
end
function cfun ( x, y, z )

!*****************************************************************************80
!
!! CFUN
!
  implicit none

  real ( kind = 8 ) cfun
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  cfun = -1.0D+00

  return
end
function dfun ( x, y, z )

!*****************************************************************************80
!
!! DFUN
!
  implicit none

  real ( kind = 8 ) dfun
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  dfun = 0.0D+00

  return
end
function efun ( x, y, z )

!*****************************************************************************80
!
!! EFUN
!
  implicit none

  real ( kind = 8 ) efun
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  efun = 0.0D+00

  return
end
function ffun ( x, y, z )

!*****************************************************************************80
!
!! FFUN
!
  implicit none

  real ( kind = 8 ) ffun
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  ffun = 0.0D+00

  return
end
function gfun ( x, y, z )

!*****************************************************************************80
!
!! GFUN
!
  implicit none

  real ( kind = 8 ) gfun
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  gfun = 0.0D+00

  return
end
subroutine ope ( n, x, y, a, ja, ia )
 
!*****************************************************************************80
!
!! OPE computes A * x for a sparse matrix A.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  integer i
  integer ia(n+1)
  integer ja(*)
  integer k1
  integer k2
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)
!
! sparse matrix * vector multiplication
!
  do i=1,n
    k1 = ia(i)
    k2 = ia(i+1) -1
    y(i) = dot_product ( a(k1:k2), x(ja(k1:k2)) )
  end do

  return
end
subroutine opet ( n, x, y, a, ja, ia )

!*****************************************************************************80
!
!! OPET computes A' * x for a sparse matrix A.
!
  implicit none

  real ( kind = 8 ) a(*)
  integer i
  integer ia(*)
  integer ja(*)
  integer k
  integer n
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)
!
! sparse matrix * vector multiplication
!
  y(1:n) = 0.0D+00

  do  i=1,n
    do k=ia(i), ia(i+1)-1
      y(ja(k)) = y(ja(k)) + x(i)*a(k)
    end do
  end do

  return
end
subroutine ydfnorm ( n, y1, y )

!*****************************************************************************80
!
!! YDFNORM prints the L2 norm of the difference of two vectors.
!
  implicit none

  integer n

  real ( kind = 8 ) t
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y1(n)

  t = sqrt ( sum ( ( y(1:n) - y1(1:n) )**2 ) )
  write(*,*) '2-norm of error (exact answer-tested answer)=',t

  return
end
subroutine dump0 ( n, a, ja, ia )

!*****************************************************************************80
!
!! DUMP0
!
  implicit none

  real ( kind = 8 ) a(*)
  integer i
  integer ia(*)
  integer ja(*)
  integer n
  integer k
  integer k1
  integer k2

  do i = 1, n
    write(*,100) i
    k1=ia(i)
    k2 = ia(i+1)-1
    write(*,101) (ja(k),k=k1,k2)
    write(*,102) (a(k),k=k1,k2)
  end do

 100  format ('row :',i2,20(2h -))
 101  format('     column indices:',10i5)
 102  format('             values:',10f5.1)

  return
end
subroutine afunbl ( nfree, x, y, z, coeff )

!*****************************************************************************80
!
!! AFUNBL
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0D+00
    end do
    coeff((j-1)*nfree+j) = -1.0D+00
  end do

  return
end
subroutine bfunbl ( nfree, x, y, z, coeff )

!*****************************************************************************80
!
!! BFUNBL
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0D+00
    end do
    coeff((j-1)*nfree+j) = -1.0D+00
  end do

  return
end
subroutine cfunbl ( nfree, x, y, z, coeff )

!*****************************************************************************80
!
!! CFUNBL
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0D+00
    end do
    coeff((j-1)*nfree+j) = -1.0D+00
  end do

  return
end
subroutine dfunbl ( nfree, x, y, z, coeff )

!*****************************************************************************80
!
!! DFUNBL 
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0D+00
    end do
  end do

  return
end
subroutine efunbl ( nfree, x, y, z, coeff )

!*****************************************************************************80
!
!! EFUNBL
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0D+00
    end do
  end do

  return
end
subroutine ffunbl ( nfree, x, y, z, coeff )

!*****************************************************************************80
!
!! FFUNBL
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0D+00
    end do
  end do

  return
end
subroutine gfunbl ( nfree, x, y, z, coeff )

!*****************************************************************************80
!
!! GFUNBL
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0D+00
    end do
  end do

  return
end
subroutine xyk ( nel, xyke, x, y, ijk, node )

!*****************************************************************************80
!
!!  XYK evaluates the material property function xyk
!
!  Discussion:
!
!    In this version of the routine, the matrix returned is the identity matrix.
!
  implicit none

  integer node

  integer ijk(node,*)
  integer nel
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xyke(2,2)
  real ( kind = 8 ) y(*)

  xyke(1,1) = 1.0D+00
  xyke(2,2) = 1.0D+00
  xyke(1,2) = 0.0D+00
  xyke(2,1) = 0.0D+00

  return
end
