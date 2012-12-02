program main

!*****************************************************************************80
!
!! MAIN is the main program for SPARSEKIT_PRB06.
!
!  Discussion:
!
!    SPARSEKIT_PRB06 demonstrates how to read a Harwell-Boeing 
!    sparse matrix file.
!
  implicit none

  integer, parameter :: nmax = 500
  integer, parameter :: nzmax = 7000

  real ( kind = 8 ) a(nzmax)
  real ( kind = 8 ) a1(nzmax)
  character ( len = 80 ) filnam
  character ( len = 2 ) guesol
  integer ia(nmax+1)
  integer ia1(nmax+1)
  integer ierr
  integer iin
  integer ios
  integer iout
  integer ja(nzmax)
  integer ja1(nzmax)
  integer job
  character ( len = 8 ) key
  integer ncol
  integer nnz
  integer nrhs
  integer nrow
  real ( kind = 8 ) rhs(1)
  character ( len = 72 ) title
  character ( len = 3 ) type
  logical valued

  iout = 6
  job = 2
  nrhs = 0

  call timestamp ( )

  write ( *, * ) ' '
  write ( *, * ) 'SPARSEKIT_PRB06'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, * ) '  This program demonstrates the use of the SPARSKIT'
  write ( *, * ) '  routines READMT and DINFO1 to read and report on'
  write ( *, * ) '  a sparse matrix stored in a file in the format'
  write ( *, * ) '  used by the Harwell-Boeing Sparse Matrix Collection'
  write ( *, * ) '  or "HBSMC".'
  write ( *, * ) ' '

  filnam = 'saylor_hb.txt'
  iin = 20

  open ( unit = iin, file = filnam, status = 'old', iostat = ios )

  if ( ios < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Error!'
    write ( *, '(a)' ) '  Unable to open file.'
    stop
  end if

  call readmt ( nmax, nzmax, job, iin, a, ja, ia, rhs, nrhs, &
    guesol, nrow, ncol, nnz, title, key, type, ierr )

  close ( unit = iin )
!
!  If not readable, return.
!
  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Error!'
    write ( *, '(a)' ) '  Unable to read matrix.'
    write ( *, '(a,i6)' ) '  READMT returned IERR = ', ierr
    stop
  end if

  valued = ( 2 <= job )

  call dinfo1 ( ncol, iout, a, ja, ia, valued, title, key, type, a1, ja1, ia1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB06'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
