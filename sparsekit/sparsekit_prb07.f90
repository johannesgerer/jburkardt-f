program main

!*****************************************************************************80
!
!! MAIN generates a Markov chain matrix to test eigenvalue routines.
!
!  The matrix A produced by this routine can also be used to generate
!  a related singular matrix, I-A.
!  The matrix models  simple random walk on a triangular grid.
!  see additional comments in subroutine.
!
!   will create a matrix in the HARWELL/BOEING format and put it in
!   the file markov.mat
!
  implicit none

  integer, parameter :: nmax = 5000
  integer, parameter :: nzmax= 4 * nmax

  real ( kind = 8 ) a(nzmax)
  integer ia(nmax+1)
  integer ifmt
  integer ios
  integer iout
  integer ja(nzmax)
  integer job
  character ( len = 8 ) key
  integer m
  integer n
  character ( len = 72 ) title
  character ( len = 3 ) type
  real ( kind = 8 ) rhs(1)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB07'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Generate a Markov chain matrix to test the'
  write ( *, '(a)' ) '  eigenvalue routine.'

  open ( unit = 11, file = 'markov.mat', status = 'replace', iostat = ios )

  if ( ios < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPARSEKIT_PRB07 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    stop
  end if
!
!  Set grid size - will not accept too large grids.
!
  m = 5

  if ( 2 * nmax < m * ( m + 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPARSEKIT_PRB07 - Fatal error!'
    write ( *, '(a)' ) '  M is too large - unable to produce matrix.'
    stop
  end if
!
!  Call the matrix generator.
!
  call markgen ( m, n, a, ja, ia )
!
!  Store result in file.
!
  title = ' Test matrix from SPARSKIT - Markov chain model           '
  key = 'randwk01'
  type = 'rua'
  iout = 11
  job = 2
  ifmt = 10

  call prtmt ( n, n, a, ja, ia, rhs, 'NN', title, key, type, ifmt, job, iout )

  close ( unit = iout )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix has been stored in Harwell/Boeing format'
  write ( *, '(a)' ) '  in the file "markov.mat".'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB07'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end

