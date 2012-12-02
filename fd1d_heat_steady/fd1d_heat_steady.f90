subroutine fd1d_heat_steady ( n, a, b, ua, ub, k, f, x, u )

!*****************************************************************************80
!
!! FD1D_HEAT_STEADY solves the steady 1D heat equation.
!
!  Discussion:
!
!    This program seeks a solution of the steady heat equation:
!
!      - d/dx ( K(X) dUdx ) = F(X)
!
!    over the interval [A,B] with boundary conditions
!
!      U(A) = UA,
!      U(B) = UB.
!
!    The code uses the finite difference method to approximate the
!    second derivative in space.  This results in a sparse linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of grid points.
!
!    Input, real ( kind = 8 ) A, B, the interval endpoints.
!
!    Input, real ( kind = 8 ) UA, UB, the values prescribed for U at 
!    the endpoints.
!
!    Input, function K(X), evaluates the thermal conductance at the N
!    points X.  Set K(X) = 1 if you don't care about this coefficient.
!
!    Input, function F(X), evaluates the heat source term at the N 
!    points X.  Set F(X) = 0 if you don't want any heat sources.
!
!    Output, real ( kind = 8 ) X(N), the grid points.
!
!    Output, real ( kind = 8 ) U(N), the approximation to the solution at 
!    the grid points.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) dx
  real ( kind = 8 ), external :: f
  integer ( kind = 4 ) i
  real ( kind = 8 ), external :: k
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) tri(3,n)
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) ua
  real ( kind = 8 ) ub
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xm
  real ( kind = 8 ) xp

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_HEAT_STEADY'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Finite difference solution of'
  write ( *, '(a)' ) '  the steady 1D heat equation'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    - d/dx ( k(x) dUdx ) = F(x)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  for space interval A <= X <= B with boundary conditions'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    U(A) = UA'
  write ( *, '(a)' ) '    U(B) = UB'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A second order difference approximation is used.'
!
!  Set the X values.
!
  dx = ( b - a ) / real ( n - 1, kind = 8 )

  do i = 1, n
    x(i) = ( real ( n - i,     kind = 8 ) * a   &
           + real (     i - 1, kind = 8 ) * b ) &
           / real ( n -     1, kind = 8 )
  end do
!
!  Set up the tridiagonal matrix.
!
  tri(1,1) = 0.0D+00
  tri(2,1) = 1.0D+00
  tri(3,1) = 0.0D+00
  rhs(1) = ua

  do i = 2, n - 1

    xm = ( x(i-1) + x(i) ) / 2.0D+00
    xp = ( x(i) + x(i+1) ) / 2.0D+00

    tri(1,i) = - k ( xm )              / dx / dx
    tri(2,i) = ( k ( xm ) + k ( xp ) ) / dx / dx
    tri(3,i) =            - k ( xp )   / dx / dx

    rhs(i) = f ( x(i) )

  end do

  tri(1,n) = 0.0D+00
  tri(2,n) = 1.0D+00
  tri(3,n) = 0.0D+00
  rhs(n) = ub
!
!  Solve the linear system.
!
  call r83np_fs ( n, tri, rhs, u )

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
subroutine r83np_fs ( n, a, b, x )

!*****************************************************************************80
!
!! R83NP_FS factors and solves an R83NP system.
!
!  Discussion:
!
!    The R83NP storage format is used for a tridiagonal matrix.
!    The subdiagonal   is in entries (1,2:N), 
!    the diagonal      is in entries (2,1:N), 
!    the superdiagonal is in entries (3,1:N-1). 
!
!    This algorithm requires that each diagonal entry be nonzero.
!    It does not use pivoting, and so can fail on systems that
!    are actually nonsingular.
!
!    The "R83NP" format used for this routine is different from the R83 format.
!    Here, we insist that the nonzero entries
!    for a given row now appear in the corresponding column of the
!    packed array.
!
!  Example:
!
!    Here is how an R83NP matrix of order 5 would be stored:
!
!       *  A21 A32 A43 A54
!      A11 A22 A33 A44 A55
!      A12 A23 A34 A45  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input/output, real ( kind = 8 ) A(3,N).
!    On input, the tridiagonal matrix.
!    On output, the data in these vectors has been overwritten
!    by factorization information.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmult
!
!  The diagonal entries can't be zero.
!
  do i = 1, n
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83NP_FS - Fatal error!'
      write ( *, '(a,i8,a)' ) '  A(2,', i, ') = 0.'
      return
    end if
  end do

  x(1:n) = b(1:n)

  do i = 2, n
    a(2,i) = a(2,i) - a(3,i-1) * a(1,i) / a(2,i-1)
    x(i)   = x(i)   - x(i-1)   * a(1,i) / a(2,i-1)
  end do

  x(n) = x(n) / a(2,n)
  do i = n-1, 1, -1
    x(i) = ( x(i) - a(3,i) * x(i+1) ) / a(2,i)
  end do

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
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) j
  character ( len = * )  output_filename
  integer   ( kind = 4 ) output_status
  integer   ( kind = 4 ) output_unit
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
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
