program main

!*****************************************************************************80
!
!! SCALAPACK_PRB1 is a simple example of the use of ScaLAPACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 September 2002
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real a(5,5)
  real b(5)
  integer context
  integer desca(9)
  integer descb(9)
  integer i
  integer iam
  integer ihi
  integer ilo
  integer info
!
!  How exactly is IPIV dimensioned?
!
  integer ipiv(10)
  integer istore
  integer j
  integer jhi
  integer jlo
  integer jstore
  integer m
  integer mycol
  integer myrow
  integer n
  integer nbrhs
  integer npcol
  integer nprocs
  integer nprow
!
!  1. Set up the processes and define the grid.
!
!  Define the number of rows and columns in the processor grid.
!
  nprow = 2
  npcol = 2
!
!  BLACS_PINFO tells this process its number, and the total number of
!  processes.
!  
  call blacs_pinfo ( iam, nprocs )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) '  Process ', iam, ' is active.'
!
!  If no processes have been set up, then we must set up all NPROW * NPCOL
!  processes via a call to BLACS_SETUP.  
!
!  We only want one process to do the setup.  There is always a process
!  number 0, so that's the one we use to take care of this chore.
!
  if ( nprocs < 1 ) then

    if ( iam == 0 ) then
      nprocs = nprow * npcol
    end if

    call blacs_setup ( iam, nprocs )

  end if
!  
!  Ask BLACS_GET to return the default system context in ICTXT.
!
  call blacs_get ( -1, 0, context )
!
!  BLACS_GRIDINIT sets up the process grid of NPROW by NPCOL
!  processes, in row-major order.  
!
!  Note that we pass the default system context in CONTEXT, and
!  BLACS_GRIDINIT replaces it with the BLACS context for this grid.
!
  call blacs_gridinit ( context, 'Row-major', nprow, npcol )
!
!  BLACS_GRIDINFO tells each process the "position" it has been
!  assigned in the process grid.
!
!  CONTEXT is input, the context handle.
!  NPROW and NPCOL are output, the number of processor rows and columns.
!  MYROW and MYCOL are output, this process's row and column.
!
  call blacs_gridinfo ( context, nprow, npcol, myrow, mycol )

  if ( myrow == -1 ) then
    call blacs_exit ( 0 )
    stop
  end if
!
!  2. Define the matrix and right hand side.
!  
!  Set the array descriptors for A and B.
!
  m = 10
  n = 10

  call descinit ( desca, m, n, 5, 5, 0, 0, context, 5, info )

  call descinit ( descb, n, 1, 5, 1, 0, 0, context, 5, info )
!
!  Store the portion of A and B that belongs to this processor.
!
  do i = 1, 5
    do j = 1, 5
      a(i,j) = min ( myrow*5 + i, mycol*5 + j )
    end do
  end do

  if ( mycol == 0 ) then
    do i = 1, 5
      if ( myrow == 1 .and. i == 1 ) then
        b(i) = 1.0
      else
        b(i) = 0.0
      end if
    end do
  end if
!
!  3. Solve the system.
!  
  if ( myrow == 0 .and. mycol == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SCALAPACK_PRB1:'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) '  ScaLAPACK demonstration.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Solve A*X=B using PSGESV, for dense matrices.'
    write ( *, '(a)' ) '  The matrix A is defined as A(I,J) = min(I,J).'
  end if
!
!  Call the ScaLAPACK routine PSGESV to solve the linear system A*X=B.
!
  call psgesv ( n, 1, a, 1, 1, desca, ipiv, b, 1, 1, descb, info )

  if ( myrow == 0 .and. mycol == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  PSGESV returned INFO = ', info
  end if
!
!  4. Shut down the processes.
!  
!  Free the BLACS context.
!
  call blacs_gridexit ( context )
!
!  Break the connection of this process to the BLACS.
!
  call blacs_exit ( 0 )

  if ( myrow == 0 .and. mycol == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SCALAPACK_PRB1:'
    write ( *, '(a)' ) '  Normal end of execution.'
  end if

  stop
end
