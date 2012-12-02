program main

!*****************************************************************************80
!
!! SCALAPACK_PRB2 solves A*x=b using the PSGESV routine.
!
!  Discussion:
!
!    This file contains an example of the use of ScaLAPACK for
!    solving a linear system.
!
!    There are six processes, arranged in a two dimensional grid of
!    two rows and three columns.  Numbering of rows and columns begins
!    at zero.
!
!    The arrangement of the processes is suggested in this diagram:
!
!    +------------+------------+------------+
!    !  Row 0     |  Row 0     |  Row 0     |
!    !  Column 0  |  Column 1  |  Column 2  |
!    +------------+------------+------------+
!    !  Row 1     |  Row 1     |  Row 1     |
!    !  Column 0  |  Column 1  |  Column 2  |
!    +------------+------------+------------+
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Integer ICTXT, specifies the BLACS context handle identifying the
!    created process grid.  The context itself is global.
!
!    Integer NPROW, the number of process rows in the grid
!    to be created.
!
!    Integer NPCOL, the number of process columns in the grid
!    to be created.
!
  integer, parameter :: dlen_ = 9
  integer, parameter :: mb = 2
  integer, parameter :: nb = 2 
  integer, parameter :: mxllda = 5
  integer, parameter :: mxlldb = 5 
  integer, parameter :: mxlocr = 5
  integer, parameter :: mxlocc = 4

  real a(mxllda,mxlocc)
  real b(mxlldb)
  integer csrc
  integer desca(dlen_)
  integer descb(dlen_)
  integer i
  integer ia
  integer iam
  integer ib
  integer ictxt
  integer ihi
  integer ilo
  integer info
  integer ipiv(mxlocr+nb)
  integer istore
  integer j
  integer ja
  integer jb
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
  integer nrhs
  integer rsrc
  integer what
!
!  SET UP THE PROCESSES.
!
!  Define the number of rows and columns in the processor grid.
!
  nprow = 2
  npcol = 3
!
!  BLACS_PINFO tells this process its number, and the total number of
!  processes.
!  
  call blacs_pinfo ( iam, nprocs )
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
!  SET UP THE BLACS GRID.
!  
!  Now that we are sure the processes have been set up, we ask BLACS_GET 
!  to return the default system context in ICTXT.
!
  what = 0
  call blacs_get ( -1, what, ictxt )
!
!  We use BLACS_GRIDINIT to set up the process grid of NPROW by NPCOL
!  processes, in row-major order.  
!
!  Note that we pass the default system context in ICTXT, and
!  BLACS_GRIDINIT replaces it with the BLACS context for this grid.
!
  call blacs_gridinit ( ictxt, 'Row-major', nprow, npcol )
!
!  Now BLACS_GRIDINFO tells each process the "position" it has been
!  assigned in the process grid.
!
!  ICTXT is input, the context handle.
!  NPROW and NPCOL are output, the number of processor rows and columns.
!  MYROW and MYCOL are output, this process's row and column.
!
  call blacs_gridinfo ( ictxt, nprow, npcol, myrow, mycol )
!
!  If BLACS_GRIDINFO returns MYROW = -1 to this process, then this
!  process is not needed in the processor grid.
!
  if ( myrow == -1 ) then
    call blacs_exit ( 0 )
    stop
  end if
!
!  DEFINE THE MATRIX AND RIGHT HAND SIDE.
!  
!  Initialize the array descriptors for the matrix A and right
!  hand side B.
!
  m = 9
  n = 9

  rsrc = 0
  csrc = 0

  call descinit ( desca, m, n, mb, nb, rsrc, csrc, ictxt, mxllda, info )

  nrhs = 1
  nbrhs = 1

  call descinit ( descb, n, nrhs, nb, nbrhs, rsrc, csrc, ictxt, mxlldb, info )
!
!  Generate the the matrix A and right hand side B.
!  Each particular process only handles a portion of this information.
!
!  In the computations below, I and J are the row and column indices
!  of the full, logical matrix.  
!
!  ISTORE and JSTORE are the row and column indices of the entry into
!  which a particular matrix element will be stored.
!
!  JSTORE will record how many matrix columns we have stored in
!  this processor's local portion of the matrix.
!
  jstore = 0
!
!  Find the first and last columns in the next block of the matrix...
!
  do jlo = nb*mycol+1, n, nb*npcol

    jhi = min ( jlo+nb-1, n )
! 
!  ...and prepare to set up column J of that block.
!
    do j = jlo, jhi
!
!  Find the first and last rows in the next block of the matrix...
!
      do ilo = mb*myrow+1, m, mb*nprow

        ihi = min ( ilo+mb-1, m )
!
!  ...and now set up (some of) the elements A(I,J).  
!
!  This portion of the matrix goes into a new column of the storage vector.
!
        jstore = jstore + 1
        istore = 0

        do i = ilo, ihi
          istore = istore + 1
          a(istore,jstore) = min ( i, j )
        end do

      end do

    end do

  end do
!
!  The right hand side vector ( or matrix ) is only given to the
!  processors with MYCOL = 0.
!
  if ( mycol == 0 ) then

    jstore = 1
!
!  Find the first and last rows in the next block of the right hand side...
!
    do ilo = mb*myrow+1, m, mb*nprow

      ihi = min ( ilo+mb-1, m )
!
!  ...and now set up (some of) the elements B(I).  
!
      istore = 0

      do i = ilo, ihi

        istore = istore + 1

        if ( i == 1 ) then
          b(istore) = 1.0
        else
          b(istore) = 0.0
        end if

      end do

    end do

  end if
!
!  SOLVE THE SYSTEM.
!  
  if ( myrow == 0 .and. mycol == 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SCALAPACK_PRB2'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, * ) '  ScaLAPACK demonstration.'
    write ( *, * ) ' '
    write ( *, * ) '  Solve A*X=B using PSGESV, for dense matrices.'
    write ( *, * ) '  The matrix A is defined as A(I,J) = min(I,J).'
    write ( *, * ) '  A has M = ', m, ' rows and N = ', n, ' columns.'
    write ( *, * ) '  The row blocking is MB = ', mb
    write ( *, * ) '  The column blocking is NB = ', nb
    write ( *, * ) '  The number of processes is NPROCS = ', nprocs
    write ( *, * ) '  The process grid contains NPROW = ', nprow, ' rows'
    write ( *, * ) '  and NPCOL = ', npcol, ' columns.'
  end if
!
!  Call the ScaLAPACK routine PSGESV to solve the linear system A*X=B.
!
  ia = 1
  ja = 1

  ib = 1
  jb = 1

  call psgesv ( n, nrhs, a, ia, ja, desca, ipiv, b, ib, jb, descb, info )

  if ( myrow == 0 .and. mycol == 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  PSGESV returned INFO = ', info
  end if
!
!  SHUT DOWN THE PROCESSES.
!  
!  Free the BLACS context.
!
  call blacs_gridexit ( ictxt )
!
!  Break the connection of this process to the BLACS.
!
  call blacs_exit ( 0 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SCALAPACK_PRB2:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
