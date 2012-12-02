program main

!*****************************************************************************80
!
!! MAIN is the main program for CRYSTAL_COORDINATES.
!
!  Discussion:
!
!    CRYSTAL_COORDINATES generates a set of FCC coordinates 
!    and writes them to a file.
!
!    The program may be used to specify the parameters for
!    a set of points in an FCC lattice, to compute the coordinates of
!    those points, and write them to a file.
!
!    The coordinate file created by this program can be read by the 
!    companion MD programs, and also by the BallRoom plotting program.
!    The points may be used as the initial locations of a set of
!    molecules in a molecular dynamics simulation.
!
!    This is meant to be an interactive program.  Prompts for the input
!    parameters are provided, as well as some suggestions for the two 
!    materials simulated by the example MD programs (Lennard-Jonesium and
!    Aluminum).
!
!  Modified:
!
!    27 October 2005
!
!  Author:
!
!    Furio Ercolessi, 
!    SISSA, Trieste, 
!    May 1995, revised May 1997
!
!  Reference:
!
!    Furio Ercolessi,
!    A Molecular Dynamics Primer.
!
  implicit none

  real ( kind = 8 ) :: alat
  real ( kind = 8 ) :: displac
  character ( len = 255 ) :: file_name
  integer ( kind = 4 ) :: nx
  integer ( kind = 4 ) :: ny
  integer ( kind = 4 ) :: nz

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CRYSTAL_COORDINATES:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Generate the coordinates of a set of points'
  write ( *, '(a)' ) '  in a 3D crystal with FCC symmetry.'
!
!  Get the parameters from the user.
!
  call read_parameters ( file_name, alat, nx, ny, nz, displac )
!
!  Print the parameters as a record.
!
  call print_parameters ( file_name, alat, nx, ny, nz, displac )
!
!  Generate the points.
!
  call generate_crystal ( file_name, alat, nx, ny, nz, displac )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The coordinate data was written to "' // &
    trim ( file_name ) // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CRYSTAL_COORDINATES:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine generate_crystal ( file_name, alat, nx, ny, nz, displac )

!*****************************************************************************80
!
!! GENERATE_CRYSTAL generates the atomic coordinates.
!
!  Discussion:
!
!    Note that there are 4 atoms in each FCC cell.
!
!  Modified:
!
!    27 October 2005
!
!  Author:
!
!    Furio Ercolessi, 
!    SISSA, Trieste, 
!    May 1995, revised May 1997
!
!  Reference:
!
!    Furio Ercolessi,
!    A Molecular Dynamics Primer.
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the 
!    output file to be created.
!
!    Input, real ( kind = 8 ) ALAT, the lattice spacing.
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the number of cells along the
!    X, Y and Z directions.
!
!    Input, real ( kind = 8 ) DISPLAC, the maximum random displacement.
!
  implicit none

  real ( kind = 8 ) :: alat
  real ( kind = 8 ), dimension(3) :: box_size
  real ( kind = 8 ) :: displac
  character ( len = * ) :: file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: nbase = 4
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz
  real ( kind = 8 ), dimension(3) :: rands
  real ( kind = 8 ), dimension(3,4), parameter :: rcell = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    0.5D+00, 0.5D+00, 0.0D+00, &
    0.0D+00, 0.5D+00, 0.5D+00, &
    0.5D+00, 0.0D+00, 0.5D+00 /), (/ 3, 4 /) )
  integer ( kind = 4 ) unit
  real ( kind = 8 ) :: x
  real ( kind = 8 ) :: y
  real ( kind = 8 ) :: z

  call get_unit ( unit )
!
!  Open file for coordinates
!
  open ( unit = unit, file = file_name, status = 'replace', &
    form = 'formatted', action = 'write', position = 'rewind', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GENERATE_CRYSTAL - Fatal error!'
    write ( *, '(a)' ) '  Cannot open the file.'
    stop
  end if
!
!  Define the number of particles N, and the box size along x,y,z.
!
  n = 4 * nx * ny * nz
  box_size(1) = nx * alat
  box_size(2) = ny * alat
  box_size(3) = nz * alat
!
!  The '%' signals to ignore this line to the BallRoom program producing
!  an image from the simulation
!
  write ( unit, '(a1,l2,i7,3e23.15)', iostat = ios ) &
    '%', .false., n, box_size(1:3)

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GENERATE_CRYSTAL - Fatal error!'
    write ( *, '(a)' ) '  Unable to write initial line of output file.'
    stop
  end if
!
!  Generate the atomic coordinates and write them to the file.
!
  do k = 0, nz - 1
    do j = 0, ny - 1
      do i = 0, nx - 1
        do l = 1, nbase

          call random_number ( rands )

          x = alat * ( real ( i, kind = 8 ) + rcell(1,l) ) &
            + 2.0D+00 * displac * ( rands(1) - 0.5D+00 )
          y = alat * ( real ( j, kind = 8 ) + rcell(2,l) ) &
            + 2.0D+00 * displac * ( rands(2) - 0.5D+00 )
          z = alat * ( real ( k, kind = 8 ) + rcell(3,l) ) &
            + 2.0D+00 * displac * ( rands(3) - 0.5D+00 )

          write ( unit, '(2x,e21.15,2x,e21.15,2x,e21.15)', iostat = ios ) &
            x, y, z

          if ( ios /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'GENERATE_CRYSTAL - Fatal error!'
            write ( *, '(a)' ) '  Unable to write initial line of output file.'
            stop
          end if

        end do
      end do
    end do
  end do

  endfile ( unit = unit )
  close ( unit = unit )

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
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

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
subroutine print_parameters ( file_name, alat, nx, ny, nz, displac )

!*****************************************************************************80
!
!! PRINT_PARAMETERS prints the user parameters.
!
!  Modified:
!
!    27 October 2005
!
!  Author:
!
!    Furio Ercolessi, 
!    SISSA, Trieste, 
!    May 1995, revised May 1997
!
!  Reference:
!
!    Furio Ercolessi,
!    A Molecular Dynamics Primer.
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the 
!    output file to be created.
!
!    Input, real ( kind = 8 ) ALAT, the lattice spacing.
!
!    Input, integer ( kind = 4 ) NX, NY, NZ, the number of cells along the
!    X, Y and Z directions.
!
!    Input, real ( kind = 8 ) DISPLAC, the maximum random displacement.
!
  implicit none

  real ( kind = 8 ) :: alat
  real ( kind = 8 ) :: displac
  character ( len = * ) :: file_name
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CRYSTAL read the following user parameters:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ALAT      = the cell spacing.'
  write ( *, '(a)' ) '  NX        = number of cells in the X direction.'
  write ( *, '(a)' ) '  NY        = number of cells in the Y direction.'
  write ( *, '(a)' ) '  NZ        = number of cells in the Z direction.'
  write ( *, '(a)' ) '  DISPLAC   = maximum random displacement.'
  write ( *, '(a)' ) '  FILE_NAME = the output file name'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  ALAT      = ', alat
  write ( *, '(a,i6)' ) '  NX        = ', nx
  write ( *, '(a,i6)' ) '  NY        = ', ny
  write ( *, '(a,i6)' ) '  NZ        = ', nz
  write ( *, '(a,g14.6)' ) '  DISPLAC   = ', displac
  write ( *, '(a)' ) '  FILE_NAME = "' // trim ( file_name ) // '".'

  return
end
subroutine read_parameters ( file_name, alat, nx, ny, nz, displac )

!*****************************************************************************80
!
!! READ_PARAMETERS obtains the parameters from the user.
!
!  Modified:
!
!    27 October 2005
!
!  Author:
!
!    Furio Ercolessi, 
!    SISSA, Trieste, 
!    May 1995, revised May 1997
!
!  Reference:
!
!    Furio Ercolessi,
!    A Molecular Dynamics Primer.
!
!  Parameters:
!
!    Output, character ( len = * ) FILE_NAME, the name of the 
!    output file to be created.
!
!    Output, real ( kind = 8 ) ALAT, the lattice spacing.
!
!    Output, integer ( kind = 4 ) NX, NY, NZ, the number of cells along the
!    X, Y and Z directions.
!
!    Output, real ( kind = 8 ) DISPLAC, the maximum random displacement.
!
  implicit none

  real ( kind = 8 ) :: alat
  real ( kind = 8 ), parameter :: cutoff_al = 5.55805441821810D+00
  real ( kind = 8 ), parameter :: cutoff_lj = 2.5D+00
  real ( kind = 8 ) :: displac
  character ( len = * ) :: file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) nz

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Reminder:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For the Lennard-Jones potential, truncated at 2.5, '
  write ( *, '(a)' ) '    the equilibrium lattice spacing is about'
  write ( *, '(a)' ) '    a_eq = 1.5496 Angstroms at P=0;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For the Al glue system,'
  write ( *, '(a)' ) '    the equilibrium lattice spacing is about'
  write ( *, '(a)' ) '    a_eq = 4.032 Angstroms at P=0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Do not deviate too much from the value that'
  write ( *, '(a)' ) '  corresponds to the system you will be using.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter the lattice spacing A: '
  read ( *, *, iostat = ios ) alat

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_PARAMETERS - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading ALAT.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reminder:'
  write ( *, '(a)' ) '  For the minimum image with the Lennard-Jones '
  write ( *, '(a)' ) '  potential, the minimum number of cells along'
  write ( *, '(a)' ) '  the X, Y, and Z directions is:'
  write ( *, '(i12)' ) int ( 2.0D+00 * cutoff_lj / alat ) + 1
  write ( *, '(a)' ) '  or, for the Al glue system:'
  write ( *, '(i12)' ) int ( 2.0D+00 * cutoff_al / alat ) + 1                    
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Enter the number of cells along the X direction: '
  read ( *, *, iostat = ios ) nx

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_PARAMETERS - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading NX.'
    stop
  end if

  write ( *, '(a)' ) &
    '  Enter the number of cells along the Y direction: '
  read ( *, *, iostat = ios ) ny

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_PARAMETERS - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading NY.'
    stop
  end if

  write ( *, '(a)' ) &
    '  Enter the number of cells along the Z direction: '
  read ( *, *, iostat = ios ) nz

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_PARAMETERS - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading NZ.'
    stop
  end if

  write ( *, '(a)' ) &
    '  Enter the maximum random displacement: '
  write ( *, '(a)' ) '  (A value of 0.05 is reasonable.)'
  read ( *, *, iostat = ios ) displac

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_PARAMETERS - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading DISPLAC.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Enter the name of the coordinate file: '

  read ( *, '(a)', iostat = ios ) file_name

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'READ_PARAMETERS - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading FILE_NAME.'
    stop
  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Modified:
!
!    06 August 2005
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
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
