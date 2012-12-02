program main

!*****************************************************************************80
!
!! MAIN is the main program for CCN_RULE.
!
!  Discussion:
!
!    This program computes a nested Clenshaw Curtis quadrature rule
!    and writes it to a file.
!
!    The user specifies:
!    * N, the number of points in the rule;
!    * A, the left endpoint;
!    * B, the right endpoint;
!    * FILENAME, the root name of the output files.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  integer ( kind = 4 ) arg_num
  real ( kind = 8 ) b
  character ( len = 255 ) filename
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) last
  integer ( kind = 4 ) n
  real ( kind = 8 ) r(2)
  character ( len = 255 ) string
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CCN_RULE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compute a nested Clenshaw Curtis rule for approximating'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Integral ( A <= x <= B ) f(x) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  of order N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The user specifies N, A, B and FILENAME.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  N is the number of points;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A is the left endpoint;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  B is the right endpoint;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FILENAME is used to generate 3 files:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    filename_w.txt - the weight file'
  write ( *, '(a)' ) '    filename_x.txt - the abscissa file.'
  write ( *, '(a)' ) '    filename_r.txt - the region file.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get N.
!
  if ( 1 <= arg_num ) then
    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, n, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the rule order N:'
    read ( *, * ) n
  end if
!
!  Get A.
!
  if ( 2 <= arg_num ) then
    iarg = 2
    call getarg ( iarg, string )
    call s_to_r8 ( string, a, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the left endpoint A:'
    read ( *, * ) a
  end if
!
!  Get B.
!
  if ( 3 <= arg_num ) then
    iarg = 3
    call getarg ( iarg, string )
    call s_to_r8 ( string, b, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the right endpoint B:'
    read ( *, * ) b
  end if
!
!  Get FILENAME.
!
  if ( 4 <= arg_num ) then
    iarg = 4
    call getarg ( iarg, filename )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter FILENAME, the "root name" of the quadrature files).'
    read ( *, '(a)' ) filename
  end if
!
!  Input summary.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N = ', n
  write ( *, '(a,g14.6)' ) '  A = ', a
  write ( *, '(a,g14.6)' ) '  B = ', b
  write ( *, '(a)'    ) '  FILENAME = "' // trim ( filename ) // '".'
!
!  Construct the rule.
!
  r(1) = a
  r(2) = b

  allocate ( w(n) )
  allocate ( x(n) )

  call ccn_compute_points ( n, x )

  x_min = -1.0D+00
  x_max = +1.0D+00
  call nc_compute ( n, x_min, x_max, x, w )
!
!  Rescale the rule.
!
  call rescale ( a, b, n, x, w )
!
!  Write the rule.
!
  call rule_write ( n, x, w, r, filename )
!
!  Free memory.
!
  deallocate ( w )
  deallocate ( x )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CCN_RULE:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine ccn_compute_points ( n, x )

!*****************************************************************************80
!
!! CCN_COMPUTE_POINTS: compute Clenshaw Curtis Nested points.
!
!  Discussion:
!
!    We want to compute the following sequence:
!
!    1/2,
!    0, 1
!    1/4, 3/4
!    1/8, 3/8, 5/8, 7/8,
!    1/16, 3/16, 5/16, 7/16, 9/16, 11/16, 13/16, 15/16, and so on.
!
!    But we'd prefer that the numbers in each row be regrouped in pairs
!    that are symmetric about 1/2, with the number above 1/2 coming first.
!    Thus, the last row might become:
!    (9/16, 7/16), (11/16, 5/16), ..., (15/16, 1/16).
!
!    Once we have our sequence, we apply the Chebyshev transformation
!    which maps [0,1] to [-1,+1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements to compute.
!
!    Output, real ( kind = 8 ) X(N), the elements of the sequence.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) td
  integer ( kind = 4 ) tu
  real ( kind = 8 ) x(n)
!
!  Handle first three entries specially.
!
  if ( 1 <= n ) then
    x(1) = 0.5D+00
  end if

  if ( 2 <= n ) then
    x(2) = 1.0D+00
  end if

  if ( 3 <= n ) then
    x(3) = 0.0D+00
  end if

  m = 3
  d = 2

  do while ( m < n )

    tu = d + 1
    td = d - 1

    k = min ( d, n - m )

    do i = 1, k
      if ( mod ( i, 2 ) == 1 ) then
        x(m+i) = real ( tu, kind = 8 ) / 2.0D+00 / real ( k, kind = 8 )
        tu = tu + 2
      else
        x(m+i) = real ( td, kind = 8 ) / 2.0D+00 / real ( k, kind = 8 )
        td = td - 2
      end if
    end do

    m = m + k
    d = d * 2

  end do
!
!  Apply the Chebyshev transformation.
!
  x(1:n) = cos ( x(1:n) * pi )

  x(1) = 0.0D+00

  if ( 2 <= n ) then
    x(2) = -1.0D+00
  end if

  if ( 3 <= n ) then
    x(3) = +1.0D+00
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

  character              ch
  integer   ( kind = 4 ) itemp

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
  logical   ch_eqi

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

  character              ch
  integer   ( kind = 4 ) digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then

    digit = iachar ( ch ) - 48

  else if ( ch == ' ' ) then

    digit = 0

  else

    digit = -1

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
!    15 January 2008
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
subroutine nc_compute ( n, x_min, x_max, x, w )

!*****************************************************************************80
!
!! NC_COMPUTE computes a Newton-Cotes quadrature rule.
!
!  Discussion:
!
!    The integral:
!
!      Integral ( X_MIN <= X <= X_MAX ) F(X) dx
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
!
!    For the CLOSED rule, the abscissas include the end points.
!    For an OPEN rule, the abscissas do not include the end points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input, real ( kind = 8 ) X_MIN, X_MAX, the endpoints of the interval.
!
!    Input, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) W(N), the weights.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) d(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) w(n)
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) x_max
  real    ( kind = 8 ) x_min
  real    ( kind = 8 ) yvala
  real    ( kind = 8 ) yvalb

  do i = 1, n
!
!  Compute the Lagrange basis polynomial which is 1 at X(I),
!  and zero at the other nodes.
!
    d(1:n) = 0.0D+00
    d(i) = 1.0D+00

    do j = 2, n
      do k = j, n
        d(n+j-k) = ( d(n+j-k-1) - d(n+j-k) ) / ( x(n+1-k) - x(n+j-k) )
      end do
    end do

    do j = 1, n - 1
      do k = 1, n - j
        d(n-k) = d(n-k) - x(n-k-j+1) * d(n-k+1)
      end do
    end do
!
!  Evaluate the antiderivative of the polynomial at the endpoints.
!
    yvala = d(n) / real ( n, kind = 8 )
    do j = n - 1, 1, -1
      yvala = yvala * x_min + d(j) / real ( j, kind = 8 )
    end do
    yvala = yvala * x_min

    yvalb = d(n) / real ( n, kind = 8 )
    do j = n - 1, 1, -1
      yvalb = yvalb * x_max + d(j) / real ( j, kind = 8 )
    end do
    yvalb = yvalb * x_max

    w(i) = yvalb - yvala

  end do

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
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
subroutine rescale ( a, b, n, x, w )

!*****************************************************************************80
!
!! RESCALE rescales a Legendre quadrature rule from [-1,+1] to [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin,
!    A fast algorithm for the calculation of the roots of special functions,
!    SIAM Journal on Scientific Computing,
!    Volume 29, Number 4, pages 1420-1438, 2007.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the new interval.
!
!    Input, integer ( kind = 4 ) N, the order.
!
!    Input/output, real ( kind = 8 ) X(N), on input, the abscissas for [-1,+1].
!    On output, the abscissas for [A,B].
!
!    Input/output, real ( kind = 8 ) W(N), on input, the weights for [-1,+1].
!    On output, the weights for [A,B].
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) w(n)
  real    ( kind = 8 ) x(n)

  x(1:n) = ( ( a + b ) + ( b - a ) * x(1:n) ) / 2.0D+00
  w(1:n) = ( b - a ) * w(1:n) / 2.0D+00

  return
end
subroutine rule_write ( order, x, w, r, filename )

!*****************************************************************************80
!
!! RULE_WRITE writes a quadrature rule to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Input, real ( kind = 8 ) X(ORDER), the abscissas.
!
!    Input, real ( kind = 8 ) W(ORDER), the weights.
!
!    Input, real ( kind = 8 ) R(2), defines the region.
!
!    Input, character ( len = * ) FILENAME, specifies the output.
!    'filename_w.txt', 'filename_x.txt', 'filename_r.txt' defining weights,
!    abscissas, and region.
!
  implicit none

  integer   ( kind = 4 )  order

  character ( len = * )   filename
  character ( len = 255 ) filename_r
  character ( len = 255 ) filename_w
  character ( len = 255 ) filename_x
  integer   ( kind = 4 )  i
  integer   ( kind = 4 )  kind
  real      ( kind = 8 )  r(2)
  real      ( kind = 8 )  w(order)
  real      ( kind = 8 )  x(order)

  filename_w = trim ( filename ) // '_w.txt'
  filename_x = trim ( filename ) // '_x.txt'
  filename_r = trim ( filename ) // '_r.txt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Creating quadrature files.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "Root" file name is   "' // trim ( filename ) // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Weight file will be   "' // trim ( filename_w ) // '".'
  write ( *, '(a)' ) '  Abscissa file will be "' // trim ( filename_x ) // '".'
  write ( *, '(a)' ) '  Region file will be   "' // trim ( filename_r ) // '".'

  call r8mat_write ( filename_w, 1, order, w )
  call r8mat_write ( filename_x, 1, order, x )
  call r8mat_write ( filename_r, 1, 2,     r )

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
!    15 January 2008
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

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) istate
  integer   ( kind = 4 ) ival
  integer   ( kind = 4 ) length
  character ( len = * )  s

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
!! S_TO_R8 reads an R8 value from a string.
!
!  Discussion:
!
!    An "R8" value is simply a real number to be stored as a
!    variable of type "real ( kind = 8 )".
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
!    12 January 2009
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

  character              c
  logical                ch_eqi
  real      ( kind = 8 ) dval
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ihave
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) iterm
  integer   ( kind = 4 ) jbot
  integer   ( kind = 4 ) jsgn
  integer   ( kind = 4 ) jtop
  integer   ( kind = 4 ) length
  integer   ( kind = 4 ) ndig
  real      ( kind = 8 ) rbot
  real      ( kind = 8 ) rexp
  real      ( kind = 8 ) rtop
  character ( len = * )  s
  integer   ( kind = 4 ) s_length
  character           :: TAB = achar ( 9 )

  s_length = len_trim ( s )

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

    if ( s_length < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' .or. c == TAB ) then

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
!  entire string, then we're done, and LENGTH is equal to S_LENGTH.
!
  if ( iterm /= 1 .and. length+1 == s_length ) then
    length = s_length
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
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 2008
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

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  character ( len = 8 )  date
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  character ( len = 10 ) time
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y
  character ( len = 5 )  zone

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
