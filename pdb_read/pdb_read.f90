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
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
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
subroutine pdb_check ( hires, hiatom, atom_max, maxres, atom_num, numres )

!*****************************************************************************80
!
!! PDB_CHECK truncates the PDB size data to internal maxima, if neccessary.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) HIATOM, the maximum atom index encountered.
!
!    Input/output, integer ( kind = 4 ) HIRES, the maximum residue index encountered.
!
!    Input, integer ( kind = 4 ) ATOM_MAX, the maximum number of atoms.
!
!    Input, integer ( kind = 4 ) MAXRES, the maximum number of residues.
!
!    Input/output, integer ( kind = 4 ) ATOM_NUM, the number of atoms read in.
!
!    Input/output, integer ( kind = 4 ) NUMRES, the number of residuals read in. 
!
  implicit none

  integer ( kind = 4 ) atom_max
  integer ( kind = 4 ) atom_num
  integer ( kind = 4 ) hiatom
  integer ( kind = 4 ) hires
  integer ( kind = 4 ) maxres
  integer ( kind = 4 ) numres

  if ( maxres < hires ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PDB_CHECK: Warning!'
    write ( *, '(a)' ) '  The highest residual index exceeds the maximum.'
    hires = maxres
  end if

  if ( maxres < numres ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PDB_CHECK: Warning!'
    write ( *, '(a)' ) '  The residual total exceeds the maximum.'
    numres = maxres
  end if

  if ( atom_max < hiatom ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PDB_CHECK: Warning!'
    write ( *, '(a)' ) '  The highest atom index exceeds the maximum.'
    hiatom = atom_max
  end if

  if ( atom_max < atom_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PDB_CHECK: Warning!'
    write ( *, '(a)' ) '  The atom total index exceeds the maximum.'
    atom_num = atom_max
  end if

  return
end
subroutine pdb_init ( coord, hiatom, hires, atom_max, maxres, mxpatm, &
  atom_num, numchain, numline, numlost, numres, prtatm, resnam, &
  resnum, xmax, xmin, ymax, ymin, zmax, zmin )

!*****************************************************************************80
!
!! PDB_INIT initializes data associated with a given PDB file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) COORD(ATOM_MAX,3) the coordinates of atoms.
!
!    Output, integer ( kind = 4 ) HIATOM, the maximum atom index encountered.
!
!    Output, integer ( kind = 4 ) HIRES, the maximum residue index encountered.
!
!    Input, integer ( kind = 4 ) ATOM_MAX, the maximum number of atoms.
!
!    Input, integer ( kind = 4 ) MAXRES, the maximum number of residues.
!
!    Input, integer ( kind = 4 ) MXPATM, the maximum number of atoms that can be
!    associated with a residue.
!
!    Output, integer ( kind = 4 ) ATOM_NUM, the number of atoms read in.
!
!    Output, integer ( kind = 4 ) NUMCHAIN, the number of chains.
!
!    Output, integer ( kind = 4 ) NUMLINE, the number of lines of text in the file.
!
!    Output, integer ( kind = 4 ) NUMLOST, the number of atom items lost.
!
!    Output, integer ( kind = 4 ) NUMRES, the number of residuals read in. 
!
!    Output, integer ( kind = 4 ) PRTATM(MAXRES,MXPATM), contains the indices of the atoms
!    that make up each residue.  In particular, for the I-th residue, 
!    the following entries are significant:
!      PRTATM(I,1) - N, the nitrogen;
!      PRTATM(I,2) - CA, the C-Alpha;
!      PRTATM(I,3) - C, the carbon;
!      PRTATM(I,4) - O, the oxygen.
!      PRTATM(I,5) - CB, another carbon?
!      PRTATM(I,5+) - the other atoms.
!
!    Output, character ( len = 3 ) RESNAM(ATOM_MAX), the residue name for each
!    entry.
!
!    Output, integer ( kind = 4 ) RESNUM(ATOM_MAX), the index of the residue to which
!    the atom belongs.
!
!    Output, real ( kind = 8 ) XMAX, XMIN, YMAX, YMIN, ZMAX, ZMIN, the 
!    maximum and minimum X, Y and Z atomic coordinates.
!
  implicit none

  integer ( kind = 4 ) atom_max
  integer ( kind = 4 ) maxres
  integer ( kind = 4 ) mxpatm

  integer ( kind = 4 ) atom_num
  real ( kind = 8 ) coord(atom_max,3)
  integer ( kind = 4 ) hiatom
  integer ( kind = 4 ) hires
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
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

  atom_num = 0
  coord(1:atom_max,1:3) = 0.0D+00
  hiatom = 0
  hires = 0
  numchain = 0
  numline = 0
  numlost = 0
  numres = 0
  prtatm(1:maxres,1:mxpatm) = 0
  resnam(1:atom_max) = ' '
  resnum(1:atom_max) = 0
  xmax = 0.0D+00
  xmin = 0.0D+00
  ymax = 0.0D+00
  ymin = 0.0D+00
  zmax = 0.0D+00
  zmin = 0.0D+00

  return
end
subroutine pdb_print_coord ( coord, atom_max, atom_num, resnam, resnum )

!*****************************************************************************80
!
!! PDB_PRINT_COORD prints out the COORD array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) COORD(ATOM_MAX,3) the coordinates of atoms.
!
!    Input, integer ( kind = 4 ) ATOM_MAX, the maximum number of atoms.
!
!    Input, integer ( kind = 4 ) ATOM_NUM, the maximum atom serial number read in.
!
!    Input, character ( len = 3 ) RESNAM(ATOM_MAX), the residue name for
!    each entry.
!
!    Input, integer ( kind = 4 ) RESNUM(ATOM_MAX), the index of the residue to which
!    the atom belongs.
!
  implicit none

  integer ( kind = 4 ) atom_max

  integer ( kind = 4 ) atom
  integer ( kind = 4 ) atom_num
  real ( kind = 8 ) coord(atom_max,3)
  integer ( kind = 4 ) j
  character ( len = 3 ) resnam(atom_max)
  integer ( kind = 4 ) resnum(atom_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Atom  ResName   Residue   Atomic XYZ coordinates:'
  write ( *, '(a)' ) ' '
 
  do atom = 1, atom_num
    write ( *, '(1x,i4,6x,a3,4x,i4,2x,3g14.6)' ) &
      atom, resnam(atom), resnum(atom), coord(atom,1:3)
  end do

  return
end
subroutine pdb_print_prtatm ( maxres, mxpatm, numres, prtatm )

!*****************************************************************************80
!
!! PDB_PRINT_PRTATM prints out the PRTATM array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAXRES, the maximum number of residues.
!
!    Input, integer ( kind = 4 ) MXPATM, the maximum number of atoms that can be
!    associated with a residue.
!
!    Input, integer ( kind = 4 ) NUMRES, the number of residuals read in.  
!
!    Input, integer ( kind = 4 ) PRTATM(MAXRES,MXPATM), contains the indices of the atoms
!    that make up each residue.  In particular, for the I-th residue, 
!    the following entries are significant:
!      PRTATM(I,1) - N, the nitrogen;
!      PRTATM(I,2) - CA, the C-Alpha;
!      PRTATM(I,3) - C, the carbon;
!      PRTATM(I,4) - O, the oxygen.
!      PRTATM(I,5) - CB, another carbon?
!      PRTATM(I,5+) - the other atoms.
!
  implicit none

  integer ( kind = 4 ) maxres
  integer ( kind = 4 ) mxpatm

  integer ( kind = 4 ) ires
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) numres
  integer ( kind = 4 ) prtatm(maxres,mxpatm)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRTATM entries:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Residue, PRTATM(I,*)'
  write ( *, '(a)' ) ' '

  do ires = 1, numres
    jhi = 0
    do j = 1, mxpatm
      if ( prtatm(ires,j) /= 0 ) then
        jhi = j
      end if
    end do
    write ( *, '(1x,19i5)' ) ires, prtatm(ires,1:jhi)
  end do

  return
end
subroutine pdb_read_atom ( coord, hiatom, hires, input_unit, atom_max, &
  maxres, mxpatm, atom_num, numchain, numline, numlost, numres, &
  prtatm, resnam, resnum, xmax, xmin, ymax, ymin, zmax, zmin )

!*****************************************************************************80
!
!! PDB_READ_ATOM reads the ATOM records in a PDB file.
!
!  Discussion:
!
!    The PDB file is presumed to have been opened by the user.
!
!    PDB_READ builds the COORD and PRTATM arrays from the ATOM data.
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
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) COORD(ATOM_MAX,3) the coordinates of atoms.
!
!    Output, integer ( kind = 4 ) HIATOM, the maximum atom index encountered.
!
!    Output, integer ( kind = 4 ) HIRES, the maximum residue index encountered.
!
!    Input, integer ( kind = 4 ) INPUT_UNIT, the FORTRAN unit number associated with
!    the file.
!
!    Input, integer ( kind = 4 ) ATOM_MAX, the maximum number of atoms.
!
!    Input, integer ( kind = 4 ) MAXRES, the maximum number of residues.
!
!    Input, integer ( kind = 4 ) MXPATM, the maximum number of atoms that can be
!    associated with a residue.
!
!    Output, integer ( kind = 4 ) ATOM_NUM, the maximum atom serial number read in.
!
!    Output, integer ( kind = 4 ) NUMCHAIN, the number of chains.
!
!    Output, integer ( kind = 4 ) NUMLINE, the number of lines of text in the file.
!
!    Output, integer ( kind = 4 ) NUMLOST, the number of atom items lost.
!
!    Output, integer ( kind = 4 ) NUMRES, the number of residuals read in.  
!
!    Output, integer ( kind = 4 ) PRTATM(MAXRES,MXPATM), contains the indices of the atoms
!    that make up each residue.  In particular, for the I-th residue, 
!    the following entries are significant:
!      PRTATM(I,1) - N, the nitrogen;
!      PRTATM(I,2) - CA, the C-Alpha;
!      PRTATM(I,3) - C, the carbon;
!      PRTATM(I,4) - O, the oxygen.
!      PRTATM(I,5) - CB, another carbon?
!      PRTATM(I,5+) - the other atoms.
!
!    Output, character ( len = 3 ) RESNAM(ATOM_MAX), the residue name for each
!    entry.
!
!    Output, integer ( kind = 4 ) RESNUM(ATOM_MAX), the index of the residue to which
!    the atom belongs.
!
!    Output, real ( kind = 8 ) XMAX, XMIN, YMAX, YMIN, ZMAX, ZMIN, the 
!    maximum and minimum X, Y and Z atomic coordinates.
!
  implicit none

  integer ( kind = 4 ) atom_max
  integer ( kind = 4 ) maxres
  integer ( kind = 4 ) mxpatm

  character altloc
  integer ( kind = 4 ) atom
  character ( len = 4 ) atom_name
  integer ( kind = 4 ) atom_num
  character chains
  character ( len = 2 ) charge
  real ( kind = 8 ) coord(atom_max,3)
  character ( len = 2 ) element
  integer ( kind = 4 ) hiatom
  integer ( kind = 4 ) hires
  integer ( kind = 4 ) ibase
  character icode
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iset
  integer ( kind = 4 ) j
  logical s_eqi
  integer ( kind = 4 ) numchain
  integer ( kind = 4 ) numline
  integer ( kind = 4 ) numlost
  integer ( kind = 4 ) numres
  real ( kind = 8 ) occ
  integer ( kind = 4 ) prtatm(maxres,mxpatm)
  character prvchn
  integer ( kind = 4 ) prvnum
  character ( len = 3 ) prvres
  real ( kind = 8 ) r1
  character ( len = 3 ) resnam(atom_max)
  character ( len = 3 ) resname
  integer ( kind = 4 ) resno
  integer ( kind = 4 ) resnum(atom_max)
  real ( kind = 8 ) rfactr
  character ( len = 4 ) segid
  character ( len = 128 ) string
  real ( kind = 8 ) temp
  character ( len = 4 ) w1
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) z
  real ( kind = 8 ) zmax
  real ( kind = 8 ) zmin
!
!  Initialize PDB information.
!
  call pdb_init ( coord, hiatom, hires, atom_max, maxres, mxpatm, &
    atom_num, numchain, numline, numlost, numres, prtatm, resnam, &
    resnum, xmax, xmin, ymax, ymin, zmax, zmin )

  ibase = 0
  prvchn = ' '
  prvnum = 0
  prvres = ' '

  iset = 0

  do

    read ( input_unit, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      exit
    end if

    numline = numline + 1
 
    if ( s_eqi ( string(1:6), 'ENDMDL' )  ) then

      exit

    else if ( s_eqi ( string(1:4), 'ATOM' ) ) then

      read ( string, &
        '(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)', &
        iostat = ios ) &
        w1, atom, atom_name, altloc, resname, chains, resno, &
        icode, x, y, z, occ, temp, segid, element, charge

      if ( ios /= 0 ) then
        exit
      end if

      atom_num = atom_num + 1
 
      hiatom = max ( hiatom, atom )
!
!  Remove a possible initial blank in ATOM_NAME or RESNAME.
!
      if ( atom_name(1:1) == ' ' ) then
        atom_name = atom_name(2:)
      end if
 
      if ( resname(1:1) == ' ' ) then
        resname = resname(2:)
      end if
!
!  If necessary, increment the number of residues read.
!
      if ( resno /= prvnum .or. resname /= prvres .or. chains /= prvchn ) then

        prvnum = resno

        prvres = resname

        if ( chains /= prvchn .and. 1 < atom ) then
          ibase = numres
          numchain = numchain + 1
        end if

        prvchn = chains

        numres = numres + 1

      end if
!
!  Correct the residue index by accounting for the chain it is on.
!
      resno = resno + ibase
!
!  Keep track of the highest residue index.
!
      hires = max ( hires, resno )
!
!  For each atom, store the atomic coordinates, and the name and number 
!  of the residue to which the atom belongs.
!
      if ( 1 <= atom .and. atom <= atom_max ) then

        coord(atom,1) = x
        coord(atom,2) = y
        coord(atom,3) = z

        resnam(atom) = resname
        resnum(atom) = numres

      end if

      if ( iset == 0 ) then
        xmax = x
        xmin = x
        ymax = y
        ymin = y
        zmax = z
        zmin = z
        iset = 1
      else
        xmax = max ( xmax, x )
        xmin = min ( xmin, x )
        ymax = max ( ymax, y )
        ymin = min ( ymin, y )
        zmax = max ( zmax, z )
        zmin = min ( zmin, z )
      end if
!
!  For each residue, add the atomic index to the PRTATM database.
!
      if ( 1 <= numres .and. numres <= maxres ) then
 
        if ( atom_name == 'C' ) then

          prtatm(numres,3) = atom

        else if ( atom_name == 'CA ' ) then

          prtatm(numres,2) = atom

        else if ( atom_name == 'CB' ) then

          prtatm(numres,5) = atom

        else if ( atom_name == 'N' ) then

          prtatm(numres,1) = atom

        else if ( atom_name == 'O' ) then

          prtatm(numres,4) = atom

        else

          numlost = numlost + 1

          do j = 6, mxpatm
            if ( prtatm(numres,j) == 0 ) then
              prtatm(numres,j) = atom
              numlost = numlost - 1
              exit
            end if
          end do

        end if

      end if
 
    end if
  
  end do
 
  return
end
subroutine pdb_summary ( hiatom, hires, atom_max, maxres, mxpatm, atom_num, &
  numchain, numline, numlost, numres, xmax, xmin, ymax, ymin, zmax, zmin )

!*****************************************************************************80
!
!! PDB_SUMMARY prints a summary of the data read from a PDB file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HIATOM, the maximum atom index encountered.
!
!    Input, integer ( kind = 4 ) HIRES, the maximum residue index encountered.
!
!    Input, integer ( kind = 4 ) ATOM_MAX, the maximum number of atoms.
!
!    Input, integer ( kind = 4 ) MAXRES, the maximum number of residues.
!
!    Input, integer ( kind = 4 ) MXPATM, the maximum number of atoms that can be
!    associated with a residue.
!
!    Input, integer ( kind = 4 ) ATOM_NUM, the number of atoms read in.
!
!    Input, integer ( kind = 4 ) NUMCHAIN, the number of chains.
!
!    Input, integer ( kind = 4 ) NUMLINE, the number of lines of text in the file.
!
!    Input, integer ( kind = 4 ) NUMLOST, the number of atom items lost.
!
!    Input, integer ( kind = 4 ) NUMRES, the number of residuals read in.  
!
!    Input, real ( kind = 8 ) XMAX, XMIN, YMAX, YMIN, ZMAX, ZMIN, the 
!    maximum and minimum X, Y and Z atomic coordinates.
!
  implicit none

  integer ( kind = 4 ) atom_max
  integer ( kind = 4 ) atom_num
  integer ( kind = 4 ) hiatom
  integer ( kind = 4 ) hires
  integer ( kind = 4 ) maxres
  integer ( kind = 4 ) mxpatm
  integer ( kind = 4 ) numchain
  integer ( kind = 4 ) numline
  integer ( kind = 4 ) numlost
  integer ( kind = 4 ) numres
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) zmax
  real ( kind = 8 ) zmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PDB_SUMMARY:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Storage limits:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '    Maximum number of atoms         ', atom_max
  write ( *, '(a,i6)' ) '    Maximum number of residues      ', maxres
  write ( *, '(a,i6)' ) '    Maximum atoms per residue       ', mxpatm
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data from the file:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '    Number of text lines in file    ', numline
  write ( *, '(a,i6)' ) '    Highest residue index           ', hires
  write ( *, '(a,i6)' ) '    Number of residues              ', numres
  write ( *, '(a,i6)' ) '    Highest atom index              ', hiatom
  write ( *, '(a,i6)' ) '    Number of atoms                 ', atom_num
  write ( *, '(a,i6)' ) '    Number of chains                ', numchain
  write ( *, '(a,i6)' ) '    Number of atom items lost       ', numlost
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Atomic coordinate ranges, in Angstroms:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '    X: ', xmin, xmax
  write ( *, '(a,2g14.6)' ) '    Y: ', ymin, ymax
  write ( *, '(a,2g14.6)' ) '    Z: ', zmin, zmax

  return
end
subroutine pdb_write_atom ( atom_num, atom_xyz, pdb_file_unit )

!*****************************************************************************80
!
!! PDB_WRITE_ATOM writes the ATOM records in a PDB file.
!
!  Discussion:
!
!    This routine accepts a set of 3D coordinates, and writes them
!    out as plausible ATOM records of a PDB file.  
!
!    The reason for doing this is to allow the data to be viewed
!    by programs which can display configurations of atoms described by
!    a PDB file.
!
!    In order to write plausible ATOM records, all the fields must
!    be filled in.  This routine makes up plausible values for these
!    other fields.
!
!  PDB ATOM Record Format:
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
!    08 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ATOM_NUM, the maximum atom serial number read in.
!
!    Input, real ( kind = 8 ) ATOM_XYZ(3,ATOM_NUM) the coordinates of atoms.
!
!    Input, integer ( kind = 4 ) PDB_FILE_UNIT, the FORTRAN unit number associated with
!    the file.
!
  implicit none

  integer ( kind = 4 ) atom_num

  character altloc
  integer ( kind = 4 ) atom
  character ( len = 4 ) atom_name
  real ( kind = 8 ) atom_xyz(3,atom_num)
  character chains
  character ( len = 2 ) charge
  character ( len = 2 ) element
  character icode
  real ( kind = 8 ) occ
  integer ( kind = 4 ) pdb_file_unit
  character ( len = 3 ) resname
  integer ( kind = 4 ) resno
  character ( len = 4 ) segid
  real ( kind = 8 ) temp

  atom_name = 'H   '
  altloc = ' '
  resname = '   '
  chains = ' '
  resno = 0
  icode = ' '
  occ = 1.0D+00
  temp = 1.0D+00
  segid = '    '
  element = '  '
  charge = '  '

  do atom = 1, atom_num

    write ( pdb_file_unit, &
      '(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ), &
      'ATOM  ', atom, atom_name, altloc, resname, chains, resno, &
      icode, atom_xyz(1:3,atom), occ, temp, segid, element, charge

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
