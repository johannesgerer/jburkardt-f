subroutine fd1d_bvp ( n, a, aprime, c, f, x, u )

!*****************************************************************************80
!
!! FD1D_BVP solves a two point boundary value problem.
!
!  Discussion:
!
!    The program uses the finite difference method to solve a BVP
!    (boundary value problem) in one dimension.
!
!    The problem is defined on the region x(1) <= x <= x(n).
!
!    The following differential equation is imposed inside the region:
!
!      - d/dx a(x) du/dx + c(x) * u(x) = f(x)
!
!    where a(x), c(x), and f(x) are given functions.  We write out
!    the equation in full as
!
!      - a(x) * u''(x) - a'(x) * u'(x) + c(x) * u(x) = f(x)
!
!    At the boundaries, the following conditions are applied:
!
!      u(x(1)) = 0.0
!      u(x(n)) = 0.0
!
!    We replace the function U(X) by a vector of N values U associated
!    with the nodes.
!
!    The first and last values of U are determined by the boundary conditions.
!
!    At each interior node I, we write an equation to help us determine
!    U(I).  We do this by approximating the derivatives of U(X) by
!    finite differences.  Let us write XL, XM, and XR for X(I-1), X(I)
!    and X(I+1).  Similarly we have UL, UM, and UR.  Other quantities to
!    be evaluated at X(I) = XM will also be labeled with an M:
!
!      - AM * ( UL - 2 UM + UR ) / DX^2 - A'M * ( UL - UR ) / ( 2 * DX ) = FM
!
!    These N-2 linear equations for the unknown coefficients complete the
!    linear system and allow us to compute the finite difference approximation
!    to the solution of the BVP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes.
!
!    Input, function A ( X ), evaluates a(x);
!
!    Input, function APRIME ( X ), evaluates a'(x);
!
!    Input, function C ( X ), evaluates c(x);
!
!    Input, function F ( X ), evaluates f(x);
!
!    Input, real ( kind = 8 ) X(N), the mesh points, which may be
!    nonuniformly spaced.
!
!    Output, real ( kind = 8 ) U(N), the value of the finite difference
!    approximation to the solution.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), external :: a
  real ( kind = 8 ) am
  real ( kind = 8 ) apm
  real ( kind = 8 ), external :: aprime
  real ( kind = 8 ), external :: c
  real ( kind = 8 ) cm
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) fm
  integer ( kind = 4 ) i
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) tri(3,n)
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xm
!
!  Equation 1 is the left boundary condition, U(X(1)) = 0.0;
!
  tri(1,1) = 0.0D+00
  tri(2,1) = 1.0D+00
  tri(3,1) = 0.0D+00
  rhs(1) = 0.0D+00
!
!  Now gather the multipliers of U(I-1) to get the matrix entry A(I,I-1),
!  and so on.
!
!  The code has been modified to allow variable X spacing.
!
  do i = 2, n - 1

    xm  = x(i)
    am  = a ( xm )
    apm = aprime ( xm )
    cm  = c ( xm )
    fm  = f ( xm )

    tri(1,i) = - 2.0D+00 * am / ( x(i) - x(i-1) ) / ( x(i+1) - x(i-1) ) &
      + apm / ( x(i+1) - x(i-1) )

    tri(2,i) = + 2.0D+00 * am / ( x(i) - x(i-1) ) / ( x(i+1) - x(i) ) &
      + cm

    tri(3,i) = - 2.0D+00 * am / ( x(i+1) - x(i) ) / ( x(i+1) - x(i-1) ) &
      - apm / ( x(i+1) - x(i-1) )

    rhs(i)   = fm

  end do
!
!  Equation N is the right boundary condition, U(X(N)) = 0.0;
!
  tri(1,n) = 0.0D+00
  tri(2,n) = 1.0D+00
  tri(3,n) = 0.0D+00
  rhs(n) = 0.0D+00
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
  logical              lopen

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
!    Input, real ( kind = 8 ) TABLE(M,N), the data.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) j
  character ( len = * )  output_filename
  integer   ( kind = 4 ) output_status
  integer   ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real      ( kind = 8 ) table(m,n)
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
subroutine r8vec_even ( n, alo, ahi, a )

!*****************************************************************************80
!
!! R8VEC_EVEN returns an R8VEC of evenly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    If N is 1, then the midpoint is returned.
!
!    Otherwise, the two endpoints are returned, and N-2 evenly
!    spaced points between them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values.
!
!    Input, real ( kind = 8 ) ALO, AHI, the low and high values.
!
!    Output, real ( kind = 8 ) A(N), N evenly spaced values.
!    Normally, A(1) = ALO and A(N) = AHI.
!    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) ahi
  real ( kind = 8 ) alo
  integer ( kind = 4 ) i

  if ( n == 1 ) then

    a(1) = 0.5D+00 * ( alo + ahi )

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = 8 ) * alo   &
             + real (     i - 1, kind = 8 ) * ahi ) &
             / real ( n     - 1, kind = 8 )
    end do

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

  character ( len = 8  ) ampm
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
