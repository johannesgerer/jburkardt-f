program main

!*****************************************************************************80
!
!! MAIN is the test program for ilut preconditioned gmres.
!
!  Discussion:
!
!    This program generates a sparse matrix using
!    matgen and then solves a linear system with an
!    artificial right hand side.
!
  implicit none

  integer, parameter :: nmax = 10 * 10 + 1
!
!  NZMAX would normally be 5*NMAX, reflecting the fact that 5 entries
!  are nonzero in each row of a matrix derived from the Laplacian.
!  However, ILUT requires that we allow a certain amount of fillin
!  during a partial factorization.
!
  integer, parameter :: nzmax = 19 * nmax

  real ( kind = 8 ) a(nzmax)
  real ( kind = 8 ) au(nzmax)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) amax
  real ( kind = 8 ) eps
  real ( kind = 8 ) gammax
  real ( kind = 8 ) gammay
  integer i
  integer ia(nmax)
  integer ierr
  integer :: im = 10
  integer :: iout = 6
  integer iw(nmax,3)
  integer j
  integer ja(nzmax)
  integer jau(nzmax)
  integer ju(nmax)
  integer k
  integer lfil
  integer :: maxits = 100
  integer meth
  integer n
  integer nwk
  integer :: nx = 10
  integer :: ny = 10
  integer :: nz = 1
  real ( kind = 8 ) random
  real ( kind = 8 ) tol
  real ( kind = 8 ) vv(nmax,20)
  real ( kind = 8 ) x(nmax)
  real ( kind = 8 ) xran(nmax)
  real ( kind = 8 ) y(nmax)

  common /func/ gammax, gammay, alpha

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB13'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Test the preconditioners and iterative solvers'
  write ( *, '(a)' ) 'in SPARSKIT.'
!
!  The PDE to be discretized is :
!
! -Lap u + gammax exp (xy)delx u + gammay exp (-xy) dely u +alpha u
!
! where Lap = 2-D laplacean, delx = part. der. wrt x,
! dely = part. der. wrt y.
! gammax, gammay, and alpha are passed via the commun func.
!
! data for PDE:
!
  alpha = -60.0
  gammax = 10.0
  gammay = 10.0
!
!  data for preconditioner
!
  nwk = nzmax
!
!  data for pgmres
!
  eps = 1.0E-07
!
!  same initial guess for gmres
!
  do j = 1, nmax
    xran(j) = random()
  end do
!
!  call gen57 to generate matrix in compressed sparse row format
!
  call gen57pt ( nx, ny, nz, a, ja, ia, ju, x )
!
!  define N.
!
  n = nx * ny * nz
!
! test all different methods:
! ILU0, MILU0, ILUT and with different values of tol and lfil
! ( from cheaper to more expensive preconditioners)
! The more accurate the preconditioner the fewer iterations
! are required in pgmres, in general.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Methods 1-5 '
  write ( *, '(a)' ) ' '

  do meth = 1, 5

    if ( meth == 1 ) then

      write ( *, '(a)' ) '1. Try ILU(0) Preconditioner'
      call ilu0 (n, a, ja, ia, au, jau, ju, iw, ierr)

    else if ( meth == 2 ) then

      write ( *, '(a)' ) '2. Try MILU(0) Preconditioner'
      call milu0 (n, a, ja, ia, au, jau, ju, iw, ierr)

    else if ( meth == 3 ) then

      write ( *, '(a)' ) '3. Try ILUT Preconditioner'
      write ( *, '(a)' ) 'with tol = 0.001, lfil=1.'
      tol  = 0.001
      lfil = 1

      call ilut (n,a,ja,ia,lfil,tol,au,jau,ju,nwk, &
        vv,vv(1,2),iw,iw(1,2),iw(1,3),ierr)

    else if ( meth == 4 ) then

      write ( *, '(a)' ) '4. Try ILUT Preconditioner'
      write ( *, '(a)' ) 'with tol = 0.001, lfil=5.'
      tol = 0.001
      lfil = 5

      call ilut (n,a,ja,ia,lfil,tol,au,jau,ju,nwk, &
        vv,vv(1,2),iw,iw(1,2),iw(1,3),ierr)

    else if ( meth == 5 ) then

      write ( *, '(a)' ) '5. Try ILUT Preconditioner'
      write ( *, '(a)' ) 'with tol = .0001, lfil=7.'
      tol = 0.0001
      lfil = 7

      call ilut (n,a,ja,ia,lfil,tol,au,jau,ju,nwk, &
        vv,vv(1,2),iw,iw(1,2),iw(1,3),ierr)

    end if
!
!  Check that return was succesful
!
    write ( *, * ) ' Precon set-up returned with ierr ', ierr

    if ( ierr /= 0 ) then
      continue
    end if
!
!  Generate right hand side = A * (1,2,3,...n)**T
!
    do k = 1, n
      x(k) = real ( k, kind = 8 )
    end do

    call ope ( n, x, y, a, ja, ia )
!
!  Generate initial guess.
!
    do j = 1, n
      x(j) = xran(j)
    end do

    call pgmres (n, im, y, x, vv, eps, maxits, iout, &
      a, ja, ia, au, jau, ju, ierr)

    write ( *, '(a,i6)' ) ' pgmres returned with ierr = ',ierr

    amax = 0.0
    do i = 1, n
      amax = max ( amax, abs ( x(i) - real ( i, kind = 8 ) ) )
    end do

    write ( *, '(a)' ) ' '
    write (*,*) 'Maximum error in solution = ', amax

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB13'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function afun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) afun, x,y, z

  afun = -1.0

  return
end
function bfun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) bfun, x,y, z

  bfun = -1.0

  return
end
function cfun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) cfun, x,y, z

  cfun = -1.0

  return
end
function dfun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) dfun, x,y, z, gammax, gammay, alpha

  common /func/ gammax, gammay, alpha

  dfun = gammax*exp(x*y)

  return
end
function efun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) efun, x,y, z, gammax, gammay, alpha
  common /func/ gammax, gammay, alpha

  efun = gammay*exp(-x*y)

  return
end
function ffun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) ffun, x,y, z

  ffun = 0.0

  return
end
function gfun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) gfun, x,y, z, gammax, gammay, alpha

  common /func/ gammax, gammay, alpha
  gfun = alpha

  return
end
subroutine afunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) coeff(*)
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  return
end
subroutine bfunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) coeff(*)
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  return
end
subroutine cfunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) coeff(*)
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  return
end
subroutine dfunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) coeff(*)
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  return
end
subroutine efunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) coeff(*)
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  return
end
subroutine ffunbl (nfree,x,y,z,coeff)
 
!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) coeff(*)
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  return
end
subroutine gfunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) coeff(*)
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  return
end
function random ( )
      
!*****************************************************************************80
!
!  This routine was extracted from ELEFUNT.
!
  integer, save :: iy = 100001
  real ( kind = 8 ) random

  iy = iy * 125
  iy = iy - (iy/2796203) * 2796203
  random = real ( iy, kind = 8 ) / 2796203.0

  return
end
