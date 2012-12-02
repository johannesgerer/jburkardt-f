program main

!*****************************************************************************80
!
!! MAIN is the main program for PDB_TO_XYZ.
!
!  Discussion:
!
!    PDB_TO_XYZ writes the ATOM coordinates of a PDB file to an XYZ file.
!
!    The program may be invoked with both files specified on the command
!    line.
!
!      pdb_to_xyz file.pdb file.xyz
!
!    If either or both files are missing, the program will ask for them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) atom_num
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) numarg
  character ( len = 80 ) pdb_file_name
  integer ( kind = 4 ) pdb_file_line_num
  integer ( kind = 4 ) pdb_file_unit
  character ( len = 80 ) xyz_file_name
  integer ( kind = 4 ) xyz_file_unit

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PDB_TO_XYZ'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Write PDB Atom XYZ coordinates to an XYZ file.'
  write ( *, '(a)' ) ' '
!
!  Get the number of command line arguments.
!
  numarg = iargc ( )

  if ( 1 <= numarg ) then

    iarg = 1
    call getarg ( iarg, pdb_file_name )

  else

    write ( *, '(a)' ) '  Enter the PDB file name:'
    read ( *, '(a)', iostat = ios ) pdb_file_name

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PDB_TO_XYZ - Fatal error!'
      write ( *, '(a)' ) '  Could not read the PDB file name.'
      stop
    end if

  end if

  if ( 2 <= numarg ) then

    iarg = 2
    call getarg ( iarg, xyz_file_name )

  else

    write ( *, '(a)' ) '  Enter the XYZ file name:'
    read ( *, '(a)', iostat = ios ) xyz_file_name

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PDB_TO_XYZ - Fatal error!'
      write ( *, '(a)' ) '  Could not read the XYZ file name.'
      stop
    end if

  end if

  write ( *, '(a)') 'Reading PDB file "' // trim ( pdb_file_name ) // '".'
!
!  Open the PDB file.
!
  call get_unit ( pdb_file_unit )

  open ( unit = pdb_file_unit, file = pdb_file_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PDB_TO_XYZ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the PDB file.'
    stop
  end if
!
!  Open the XYZ file.
!
  call get_unit ( xyz_file_unit )

  open ( unit = xyz_file_unit, file = xyz_file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PDB_TO_XYZ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the XYZ file.'
    stop
  end if

  write ( xyz_file_unit, '(a)' ) '#  ' // trim ( xyz_file_name )
  write ( xyz_file_unit, '(a)' ) '#  created by PDB_TO_XYZ'
  write ( xyz_file_unit, '(a)' ) '#'
!
!  Read the ATOM coordinate information from the file.
!
  call pdb_atom_to_xyz ( pdb_file_unit, xyz_file_unit, pdb_file_line_num, &
    atom_num )

  write ( *, '(a)' ) ' '
  write ( *, '(2x,i6,a)' ) pdb_file_line_num, ' lines read from PDB file.'
  write ( *, '(2x,i6,a)' ) atom_num, ' atom records found in PDB file.'
!
!  Close the PDB file.
!
  close ( unit = pdb_file_unit )
!
!  Close the XYZ file.
!
  close ( unit = xyz_file_unit )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PDB_TO_XYZ'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine pdb_atom_to_xyz ( pdb_file_unit, xyz_file_unit, pdb_file_line_num, &
  atom_num )

!*****************************************************************************80
!
!! PDB_ATOM_TO_XYZ reads the ATOM records in a PDB file.
!
!  Discussion:
!
!    The PDB and XYZ files are presumed to have been opened by the user.
!
!  Format:
!
!    COLUMNS  DATA TYPE     FIELD       DEFINITION
!    --------------------------------------------------------------------------
!     1 -  6  Record name   "ATOM  "
!     7 - 11  Integer       serial      Atom serial number.
!    13 - 16  Atom          name        Atom name.
!    17       Character     altLoc      Alternate location indicator.
!    18 - 20  Residue name  resName     Residue name.
!    22       Character     chainID     Chain identifier.
!    23 - 26  Integer       resSeq      Residue sequence number.
!    27       AChar         iCode       Code for insertion of residues.
!    31 - 38  Real(8.3)     x           Orthogonal coordinates for X, Angstroms.
!    39 - 46  Real(8.3)     y           Orthogonal coordinates for Y, Angstroms.
!    47 - 54  Real(8.3)     z           Orthogonal coordinates for Z, Angstroms.
!    55 - 60  Real(6.2)     occupancy   Occupancy.
!    61 - 66  Real(6.2)     tempFactor  Temperature factor.
!    73 - 76  LString(4)    segID       Segment identifier, left-justified.
!    77 - 78  LString(2)    element     Element symbol, right-justified.
!    79 - 80  LString(2)    charge      Charge on the atom.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PDB_FILE_UNIT, the FORTRAN unit number associated with
!    the PDB file.
!
!    Input, integer ( kind = 4 ) XYZ_FILE_UNIT, the FORTRAN unit number associated with
!    the XYZ file.
!
!    Output, integer ( kind = 4 ) PDB_FILE_LINE_NUM, the number of lines read from the
!    PDB file.
!
!    Output, integer ( kind = 4 ) ATOM_NUM, the number of ATOM records read.
!
  implicit none

  integer ( kind = 4 ) atom_num
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) pdb_file_line_num
  integer ( kind = 4 ) pdb_file_unit
  logical s_eqi
  character ( len = 80 ) string
  real ( kind = 8 ) x
  integer ( kind = 4 ) xyz_file_unit
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  pdb_file_line_num = 0
  atom_num = 0

  do

    read ( pdb_file_unit, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      exit
    end if

    pdb_file_line_num = pdb_file_line_num + 1

    if ( s_eqi ( string(1:6), 'ENDMDL' )  ) then

      exit

    else if ( s_eqi ( string(1:4), 'ATOM' ) ) then

      read ( string, &
        '(6x,5x,1x,4x,1x,3x,1x,1x,4x,1x,3x,3f8.3,12x,6x,4x,2x,2x)', &
        iostat = ios ) x, y, z

      if ( ios /= 0 ) then
        exit
      end if

      write ( xyz_file_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, y, z

      atom_num = atom_num + 1

    end if

  end do

  return
end
function s_eqi ( strng1, strng2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( s1 /= s2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
