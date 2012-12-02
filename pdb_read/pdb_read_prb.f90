program main

!*****************************************************************************80
!
!! MAIN is the main program for PDB_READ_PRB.
!
!  Discussion:
!
!    PDB_READ_PRB tests the PDB_READ reading routine.
!
!    PDB_READ_PRB can be invoked from the command line, with the file
!    to be read included as a command line argument:
!
!      pdb_read_prb file.pdb > report.txt
!
!  Modified:
!
!    07 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: atom_max = 3000
  integer ( kind = 4 ), parameter :: maxres = 1000
  integer ( kind = 4 ), parameter :: mxpatm = 18

  integer ( kind = 4 ) atom_num
  real ( kind = 8 ) coord(atom_max,3)
  character ( len = 80 ) filepdb
  integer ( kind = 4 ) hiatom
  integer ( kind = 4 ) hires
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilen
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ipxfargc
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) numarg
  integer ( kind = 4 ) numchain
  integer ( kind = 4 ) numline
  integer ( kind = 4 ) numlost
  integer ( kind = 4 ) numres
  integer ( kind = 4 ) prtatm(maxres,mxpatm)
  character ( len = 3 ) resnam(atom_max)
  integer ( kind = 4 ) resnum(atom_max)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) zmax
  real ( kind = 8 ) zmin

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PDB_READ_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PDB_READ library.'
!
!  Get the number of command line arguments.
!
!  Old style:
!
  numarg = iargc ( )
!
!  New style:
!
! numarg = ipxfargc ( )

  if ( 1 <= numarg ) then

    iarg = 1
!
!  Old style:
!
    call getarg ( iarg, filepdb )
!
!  New style:
!
!   call pxfgetarg ( iarg, filepdb, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'PDB_READ_PRB - Fatal error!'
!     write ( *, '(a)' ) '  Could not read a command line argument.'
!     stop
!   end if

  else

    write ( *, '(a)' ) '  Enter the PDB file name:'
    read ( *, '(a)', iostat = ios ) filepdb

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PDB_READ_PRB - Fatal error!'
      write ( *, '(a)' ) '  Could not read the PDB file name.'
      stop
    end if

  end if

  write ( *, '(a)') 'Reading PDB file : ' // trim ( filepdb )
!
!  Open the PDB file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = filepdb, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PDB_READ_PRB - Fatal error!'
    write ( *, '(a)' ) '  Could not open the PDB file.'
    stop
  end if
!
!  Read the ATOM information from the file.
!
  call pdb_read_atom ( coord, hiatom, hires, input_unit, atom_max, maxres, &
    mxpatm, atom_num, numchain, numline, numlost, numres, prtatm, &
    resnam, resnum, xmax, xmin, ymax, ymin, zmax, zmin )
!
!  Close the PDB file.
!
  close ( unit = input_unit )
!
!  Print a summary of the information.
!
  call pdb_summary ( hiatom, hires, atom_max, maxres, mxpatm, &
    atom_num, numchain, numline, numlost, numres, xmax, xmin, &
    ymax, ymin, zmax, zmin )
!
!  Reset any data that is out of bounds.
!
  call pdb_check ( hires, hiatom, atom_max, maxres, atom_num, numres )
!
!  Print out the PRTATM array.
!
  call pdb_print_prtatm ( maxres, mxpatm, numres, prtatm )
!
!  Print out the COORD array.
!
  call pdb_print_coord ( coord, atom_max, atom_num, resnam, resnum )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PDB_READ_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
