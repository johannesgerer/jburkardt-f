program main

!*****************************************************************************80
!
!! MAIN is the main program for SHALLOW_WATER_1D.
!
!  Discussion:
!
!    SHALLOW_WATER_1D approximates the 1D shallow water equations.
!
!    This code can be considered a 1D version of Cleve Moler's shallow
!    water equation solver.
!
!    The version of the shallow water equations being solved here is in
!    conservative form, and omits the Coriolis force.  The state variables
!    are H (the height) and UH (the mass velocity).
!
!    The equations have the form
!
!      dH/dt + d UH/dx = 0
!
!      d UH/dt + d ( U^2 H + 1/2 g H^2 )/dx = 0
!
!    Here U is the ordinary velocity, U = UH/H, and g is the gravitational
!    acceleration.
!
!    The initial conditions are used to specify ( H, UH ) at an equally
!    spaced set of points, and then the Lax-Wendroff method is used to advance
!    the solution through a number of equally spaced points in time, with 
!    boundary conditions supplying the first and last spatial values.
!
!
!    Some input values will result in an unstable calculation that
!    quickly blows up.  This is related to the Courant-Friedrichs-Lewy
!    condition, which requires that DT be small enough, relative to DX and
!    the velocity, that information cannot cross an entire cell.
!
!    A "reasonable" set of input quantities is
!
!      shallow_water_1d 41 100 1.0 0.2 9.8
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Cleve Moler,
!    "The Shallow Water Equations",
!    Experiments with MATLAB.
!
!  Parameters:
!
!    Input, integer NX, the number of spatial nodes.
!
!    Input, integer NT, the number of times steps.
!
!    Input, real X_LENGTH, the length of the region.
!
!    Input, real T_LENGTH, the time extent.
!
!    Input, real G, the gravity constant.  G = 9.8 meters per second**2.
!
!    Output, real H_ARRAY(NX,NT+1), the height for all space and time points.
!
!    Output, real UH_ARRAY(NX,NT+1), the mass velocity for all space and time points.
!
!    Output, real X(NX), the X coordinates.
!
!    Output, real T(NT+1), the T coordinates.
!
  implicit none

  integer ( kind = 4 ) arg_num
  real ( kind = 8 ) dx
  real ( kind = 8 ) dt
  character ( len = 80 ) filename_h
  character ( len = 80 ) filename_t
  character ( len = 80 ) filename_uh
  character ( len = 80 ) filename_x
  real ( kind = 8 ) g
  real ( kind = 8 ), allocatable :: h(:)
  real ( kind = 8 ), allocatable :: h_array(:,:)
  real ( kind = 8 ), allocatable :: hm(:)
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) it
  integer ( kind = 4 ) last
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nx
  character ( len = 255 ) string
  real ( kind = 8 ), allocatable :: t(:)
  real ( kind = 8 ) t_length
  real ( kind = 8 ), allocatable :: uh(:)
  real ( kind = 8 ), allocatable :: uh_array(:,:)
  real ( kind = 8 ), allocatable :: uhm(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) x_length

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHALLOW_WATER_1D'
  write ( *, '(a)' ) '  FORTRAN90 version'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )

  if ( arg_num < 5 ) then
    write ( *, '(a)' ) ' '
  end if
!
!  Get the quadrature file root name:
!
  if ( arg_num < 1 ) then
    nx = 41
  else
    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, nx, ierror, last )
  end if
  write ( *, '(a,i6)' ) '  NX = ', nx

  if ( arg_num < 2 ) then
    nt = 100
  else
    iarg = 2
    call getarg ( iarg, string )
    call s_to_i4 ( string, nt, ierror, last )
  end if
  write ( *, '(a,i6)' ) '  NT = ', nt

  if ( arg_num < 3 ) then
    x_length = 1.0D+00
  else
    iarg = 3
    call getarg ( iarg, string )
    call s_to_r8 ( string, x_length, ierror, last )
  end if
  write ( *, '(a,g14.6)' ) '  X_LENGTH = ', x_length

  if ( arg_num < 4 ) then
    t_length = 0.2D+00
  else
    iarg = 3
    call getarg ( iarg, string )
    call s_to_r8 ( string, t_length, ierror, last )
  end if
  write ( *, '(a,g14.6)' ) '  T_LENGTH = ', t_length

  if ( arg_num < 5 ) then
    g = 9.8D+00
  else
    iarg = 5
    call getarg ( iarg, string )
    call s_to_r8 ( string, g, ierror, last )
  end if
  write ( *, '(a,g14.6)' ) '  G = ', g
!
!  Allocate space.
!
  allocate ( h(1:nx) )
  allocate ( h_array(1:nx,1:nt+1) )
  allocate ( hm(1:nx-1) )
  allocate ( t(1:nt+1) )
  allocate ( uh(1:nx) )
  allocate ( uh_array(1:nx,1:nt+1) )
  allocate ( uhm(1:nx-1) )
  allocate ( x(1:nx) )
!
!  Define the locations of the nodes and time steps and the spacing.
!
  call r8vec_linspace ( nx, 0.0D+00, x_length, x )
  call r8vec_linspace ( nt + 1, 0.0D+00, t_length, t )

  dx = x_length / real ( nx - 1, kind = 8 )
  dt = t_length / real ( nt, kind = 8 )

  write ( *, * ) 'DX = ', dx
  write ( *, * ) 'DT = ', dt
!
!  Apply the initial conditions.
!
  call initial_conditions ( nx, nt, x, t(1), h, uh )
!
!  Apply the boundary conditions.
!
  call boundary_conditions ( nx, nt, x, t(1), h, uh )
!
!  Store the first time step into H_ARRAY and UH_ARRAY.
!
  h_array(1:nx,1) = h(1:nx)
  uh_array(1:nx,1) = uh(1:nx)
!
!  Take NT more time steps.
!
  do it = 1, nt
!
!  Take a half time step, estimating H and UH at the NX-1 spatial midpoints.
!
    hm(1:nx-1) = ( h(1:nx-1) + h(2:nx) ) / 2.0D+00 &
      - ( dt / 2.0D+00 ) * ( uh(2:nx) - uh(1:nx-1) ) / dx

    uhm(1:nx-1) = ( uh(1:nx-1) + uh(2:nx) ) / 2.0D+00 &
      - ( dt / 2.0D+00 ) * ( &
        uh(2:nx)**2    / h(2:nx)   + 0.5D+00 * g * h(2:nx)**2 &
      - uh(1:nx-1)**2  / h(1:nx-1) - 0.5D+00 * g * h(1:nx-1)**2 ) / dx
!
!  Take a full time step, evaluating the derivative at the half time step,
!  to estimate the solution at the NX-2 nodes.
!
    h(2:nx-1) = h(2:nx-1) &
      - dt * ( uhm(2:nx-1) - uhm(1:nx-2) ) / dx

    uh(2:nx-1) = uh(2:nx-1) &
      - dt * ( &
        uhm(2:nx-1)**2  / hm(2:nx-1) + 0.5D+00 * g * hm(2:nx-1)**2 &
      - uhm(1:nx-2)**2  / hm(1:nx-2) - 0.5D+00 * g * hm(1:nx-2)**2 ) / dx
!
!  Update the boundary conditions.
!
    call boundary_conditions ( nx, nt, x, t(it+1), h, uh )
!
!  Copy data into the big arrays.
!
    h_array(1:nx,it+1) = h(1:nx)
    uh_array(1:nx,it+1) = uh(1:nx)
    
  end do
!
!  Write data to files.
!
  filename_x = 'sw1d_x.txt'
  filename_t = 'sw1d_t.txt'
  filename_h = 'sw1d_h.txt'
  filename_uh = 'sw1d_uh.txt'

  call r8vec_write ( filename_x, nx, x )
  call r8vec_write ( filename_t, nt + 1, t )
  call r8mat_write ( filename_h, nx, nt + 1, h_array )
  call r8mat_write ( filename_uh, nx, nt + 1, uh_array )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X  values saved in file "' // trim ( filename_x ) // '".'
  write ( *, '(a)' ) '  T  values saved in file "' // trim ( filename_t ) // '".'
  write ( *, '(a)' ) '  H  values saved in file "' // trim ( filename_h ) // '".'
  write ( *, '(a)' ) '  UH values saved in file "' // trim ( filename_uh ) // '".'
!
!  Free memory.
!
  deallocate ( h )
  deallocate ( h_array )
  deallocate ( hm )
  deallocate ( t )
  deallocate ( uh )
  deallocate ( uh_array )
  deallocate ( uhm )
  deallocate ( x )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHALLOW_WATER_1D:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine boundary_conditions ( nx, nt, x, t, h, uh )

!*****************************************************************************80
!
!! INITIAL_CONDITIONS sets the initial conditions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, the number of spatial nodes.
!
!    Input, integer ( kind = 4 ) NT, the number of times steps.
!
!    Input, real ( kind = 8 ) X(NX), the coordinates of the nodes.
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input/output, real ( kind = 8 ) H(NX), the height, with H(1) and H(NX) 
!    adjusted for boundary conditions.
!
!    Input/output, real ( kind = 8 ) UH(NX), the mass velocity, with UH(1) 
!    and UH(NX) adjusted for boundary conditions.
!
  implicit none

  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nx

  integer ( kind = 4 ) bc
  real ( kind = 8 ) h(nx)
  real ( kind = 8 ) t
  real ( kind = 8 ) uh(nx)
  real ( kind = 8 ) x(nx)

  bc = 1
!
!  Periodic boundary conditions on H and UH.
!
  if ( bc == 1 ) then
    h(1) = h(nx-1)
    h(nx) = h(2)
    uh(1) = uh(nx-1)
    uh(nx) = uh(2)
!
!  Free boundary conditions on H and UH.
!
  else if ( bc == 2 ) then
    h(1) = h(2)
    h(nx) = h(nx-1)
    uh(1) = uh(2)
    uh(nx) = uh(nx-1)
!
!  Reflective boundary conditions on UH, free boundary conditions on H.
!
  else if ( bc == 3 ) then
    h(1) = h(2)
    h(nx) = h(nx-1)
    uh(1) = - uh(2)
    uh(nx) = - uh(nx-1)
  end if

  return
end
subroutine ch_cap ( ch )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Discussion:
!
!    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions,
!    which guarantee the ASCII collating sequence.
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
!    Input/output, character CH, the character to capitalize.
!
  implicit none

  character ch
  integer ( kind = 4 ) itemp

  itemp = iachar ( ch )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    ch = achar ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Discussion:
!
!    CH_EQI ( 'A', 'a' ) is TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical ch_eqi

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( ch, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!     CH  DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character CH, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.
!    If CH was 'illegal', then DIGIT is -1.
!
  implicit none

  character ch
  integer ( kind = 4 ) digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then

    digit = iachar ( ch ) - 48

  else if ( ch == ' ' ) then

    digit = 0

  else

    digit = - 1

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
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2008
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
subroutine initial_conditions ( nx, nt, x, t, h, uh )

!*****************************************************************************80
!
!! INITIAL_CONDITIONS sets the initial conditions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, the number of spatial nodes.
!
!    Input, integer ( kind = 4 ) NT, the number of times steps.
!
!    Input, real ( kind = 8 ) X(NX), the coordinates of the nodes.
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Output, real ( kind = 8 ) H(NX), the initial height for all space.
!
!    Output, real ( kind = 8 ) UH(NX), the initial mass velocity for all space.
!
  implicit none

  integer ( kind = 4 ) nt
  integer ( kind = 4 ) nx

  real ( kind = 8 ) h(nx)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) uh(nx)
  real ( kind = 8 ) x(nx)

  h(1:nx) = 2.0D+00 + sin ( 2.0D+00 * pi * x(1:nx) )

  uh(1:nx) = 0.0D+00

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real ( kind = 8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine r8vec_linspace ( n, a_first, a_last, a )

!*****************************************************************************80
!
!! R8VEC_LINSPACE creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A_FIRST, A_LAST, the first and last entries.
!
!    Output, real ( kind = 8 ) A(N), a vector of linearly spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_first
  real ( kind = 8 ) a_last
  integer ( kind = 4 ) i

  if ( n == 1 ) then

    a(1) = ( a_first + a_last ) / 2.0D+00

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = 8 ) * a_first &
             + real (     i - 1, kind = 8 ) * a_last ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine r8vec_write ( output_filename, n, x )

!*****************************************************************************80
!
!! R8VEC_WRITE writes an R8VEC file.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the data.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
  real ( kind = 8 ) x(n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if

  if ( 0 < n ) then
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, '(2x,g24.16)' ) x(j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine s_to_i4 ( s, ival, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters of S 
!    used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) length
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        length = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    length = len_trim ( s )
  else
    ierror = 1
    length = 0
  end if

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) DVAL, the value read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical ch_eqi
  character c
  real ( kind = 8 ) dval
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) length
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s

  nchar = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( nchar < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Scientific notation exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to NCHAR.
!
  if ( iterm /= 1 .and. length+1 == nchar ) then
    length = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
        / real ( jbot, kind = 8 ) )
    end if
  end if

  dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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

