program main

!*****************************************************************************80
!
!! MAIN is the main program for BASIS_COMPARE.
!
!  Discussion:
!
!    BASIS_COMPARE compares two basis sets.
!
!    In an M-dimensional space, a set of N orthonormal vectors U is provided,
!    which span some subspace.  We let Q be the M by N matrix whose columns
!    are the vectors U(1:N).
!
!    One or more M-vectors V are then considered.  For each V, it is
!    desired to know the proportion of the vector that lies
!    within and outside the subspace spanned by U.  V admits a decomposition
!    of the form
!
!      V = V1 + V2
!
!    where
!
!      V1 is a vector in the subspace spanned by U,
!      V2 is a vector orthogonal to the subspace spanned by U,
!      sqrt ( ||V1||^2 + ||V2||^2 ) = ||V|| = Euclidean length of V.
!
!    The values
!
!      ALPHA = ||V1|| / ||V||
!      BETA  = ||V2|| / ||V||
!
!    must lie between 0 and 1; their squares sum to 1, and indicate the portion
!    of V which lies within the U subspace, and the portion which
!    lies in the orthogonal complement.  Equivalently, the angle THETA
!    between V and the subspace is determined by
!
!      TAN ( THETA ) = BETA / ALPHA
!
!    To compute ALPHA, we realize that Q' * V produces the N by 1 vector
!    of projection coefficients of V onto each U(I).  Therefore, V1, the
!    orthogonal projection of V into U is Q * ( Q' * V ):
!
!      V1 =       Q * Q'   * V
!
!    and the "remainder" is:
!
!      V2 = ( I - Q * Q' ) * V
!
!    from which we can determine ALPHA and BETA.
!
!  Modified:
!
!    30 April 2003
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( : ) :: a
  real ( kind = 8 ) alpha
  real ( kind = 8 ), allocatable, dimension ( : ) :: b
  real ( kind = 8 ) beta
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: c
  real ( kind = 8 ) c_max
  character ( len = 100 ) file_name
  logical file_name_exist
  character ( len = 100 ) file_name_first
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_half
  integer ( kind = 4 ) n
  real ( kind = 8 ) theta
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: u
  real ( kind = 8 ), allocatable, dimension ( : ) :: v
  real ( kind = 8 ) v_norm
  real ( kind = 8 ), allocatable, dimension ( : ) :: v1
  real ( kind = 8 ) v1_norm
  real ( kind = 8 ), allocatable, dimension ( : ) :: v2
  real ( kind = 8 ) v2_norm

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_COMPARE:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Let U be a set of N vectors in M-dimensional space.'
  write ( *, '(a)' ) '  Regard U as a basis for some subspace.'
  write ( *, '(a)' ) '  Now, for one or more vectors V, determine'
  write ( *, '(a)' ) '  the orthogonal decomposition:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    V = V1 + V2'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    V1 lies in subspace(U), and'
  write ( *, '(a)' ) '    V2 is orthogonal to V1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Defining:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    ALPHA = ||V1|| / ||V||,'
  write ( *, '(a)' ) '    BETA  = ||V2|| / ||V||,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  we have ALPHA^2 + BETA^2 = 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The values ALPHA and BETA indicate'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    ALPHA: how much of V lies IN subspace(U), and'
  write ( *, '(a)' ) '    BETA:  how much lies in the orthogonal complement.'
  write ( *, '(a)' ) ' '
!
!  Get the name of the first file in the basis set U.
!
  call s_input ( 'Enter the name of the first file in the basis set U:', &
    file_name_first, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASIS_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Error reading user input.'
    stop
  end if
!
!  Count the number of lines in the file to get the spatial dimension M.
!
  call file_row_count ( file_name_first, m_half )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of data lines in file is ', m_half
  m = 2 * m_half
  write ( *, '(a,i6)' ) '  Spatial dimension M = ', m
!
!  Count the number of basis files to get N, the number of basis vectors.
!
  file_name = file_name_first

  n = 0

  do

    inquire ( file = file_name, exist = file_name_exist )

    if ( .not. file_name_exist ) then
      exit
    end if

    n = n + 1
    call file_name_inc ( file_name )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of basis vectors is N = ', n
!
!  Allocate U for the basis vectors, and V1 and V2 for the two components.
!
  allocate ( a(1:m_half) )
  allocate ( b(1:m_half) )
  allocate ( c(1:n,1:n) )
  allocate ( u(1:m,1:n) )
  allocate ( v(1:m) )
  allocate ( v1(1:m) )
  allocate ( v2(1:m) )
!
!  Read the basis vectors into the basis vector matrix U.
!
  file_name = file_name_first

  do j = 1, n
    call data_d2_read ( file_name, m_half, a, b, ierror )
    u(1:m-1:2,j) = a(1:m_half)
    u(2:m  :2,j) = b(1:m_half)
    call file_name_inc ( file_name )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum entry of |U(1:M,1:N)| = ', &
    maxval ( abs ( u(1:m,1:n) ) )
!
!  Verify that the columns of U are orthonormal.
!
  c(1:n,1:n) = matmul ( transpose ( u(1:m,1:n) ), u(1:m,1:n) )
  do i = 1, n
    c(i,i) = c(i,i) - 1.0D+00
  end do

  c_max = maxval ( abs ( c(1:n,1:n) ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum value of U'' * U - I = ', c_max
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (If the columns of U were perfectly orthonormal,'
  write ( *, '(a)' ) '  this would be exactly zero.)'
!
!  Get the name of the first file in the set of test vectors V.
!
  call s_input ( 'Enter the name of the first file in the test vector set:', &
    file_name_first, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASIS_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Error reading user input.'
    stop
  end if
!
!  For each test vector:
!
  file_name = file_name_first

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     ||V||         Alpha        Beta        Theta'
  write ( *, '(a)' ) '                                           (Degrees)'
  write ( *, '(a)' ) ' '

  do
!
!  Read in the data comprising vector V.
!
    call data_d2_read ( file_name, m_half, a, b, ierror )

    if ( ierror /= 0 ) then
      exit
    end if

    v(1:m-1:2) = a(1:m_half)
    v(2:m  :2) = b(1:m_half)
!
!  Compute V1 = U * U' * V
!
    v1(1:m) = matmul ( u(1:m,1:n), &
                matmul ( transpose ( u(1:m,1:n) ), v(1:m) )  &
              )
    v2(1:m) = v(1:m) - v1(1:m)

    v_norm  = sqrt ( dot_product ( v(1:m),  v(1:m)  ) )
    v1_norm = sqrt ( dot_product ( v1(1:m), v1(1:m) ) )
    v2_norm = sqrt ( dot_product ( v2(1:m), v2(1:m) ) )

    alpha = v1_norm / v_norm
    beta  = v2_norm / v_norm
    theta = atan2 ( beta, alpha ) * 180.0D+00 / 3.14159265D+00

    write ( *, '(4g14.6)' ) v_norm, alpha, beta, theta
!
!  Increment the file name.
!
    call file_name_inc ( file_name )

  end do
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( b )
  deallocate ( c )
  deallocate ( u )
  deallocate ( v )
  deallocate ( v1 )
  deallocate ( v2 )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_COMPARE:'
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
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Examples:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
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
function ch_is_digit ( c )

!*****************************************************************************80
!
!! CH_IS_DIGIT returns .TRUE. if a character is a decimal digit.
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, .TRUE. if C is a digit, .FALSE. otherwise.
!
  implicit none

  character c
  logical ch_is_digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
    ch_is_digit = .true.
  else
    ch_is_digit = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
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
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  
!    If C was 'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine data_d2_read ( file_name, n, x, y, ierror )

!*****************************************************************************80
!
!! DATA_D2_READ reads a data set of pairs of real numbers stored in a file.
!
!  Discussion:
!
!    The data set can be thought of as a real M by 2 array.
!
!    Each row of the array corresponds to one data "item".
!
!    The data is stored in a file, one row (pair of values) at a time.
!
!    Each row (pair of values) begins on a new line.
!
!    Blank lines and comment lines (beginning with '#') are ignored.
!
!    Individual data values should be separated by spaces or commas.
!
!  Modified:
!
!    22 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to read.
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Output, real ( kind = 8 ) X(N), Y(N), the data values.
!
!    Output, integer ( kind = 4 ) IERROR, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) n

  character ( len = * ) file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) last
  integer ( kind = 4 ) length
  character ( len = 80 ) line
  integer ( kind = 4 ) n2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  ierror = 0

  call get_unit ( input )

  open ( unit = input, file = file_name, iostat = ios, status = 'old' )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_D2_READ - Fatal error!'
    write ( *, '(a)' ) '  Error opening the file.'
    return
  end if

  x(1:n) = huge ( x(1) )
  y(1:n) = huge ( y(1) )

  n2 = 0

  do

    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DATA_D2_READ - Fatal error!'
      write ( *, '(a)' ) '  Error reading the file.'
      return
    end if

    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      n2 = n2 + 1

      last = 0
      call s_to_r8 ( line(last+1:), x(n2), ierror, length )

      if ( ierror /= 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DATA_D2_READ - Fatal error!'
        write ( *, '(a)' ) '  Error decoding information in the file.'
        return
      end if

      last = last + length

      call s_to_r8 ( line(last+1:), y(n2), ierror, length )

      if ( ierror /= 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DATA_D2_READ - Fatal error!'
        write ( *, '(a)' ) '  Error decoding information in the file.'
        return
      end if

      if ( n2 == n ) then
        exit
      end if

    end if

  end do

  close ( unit = input )

  return
end
subroutine digit_inc ( c )

!*****************************************************************************80
!
!! DIGIT_INC increments a decimal digit.
!
!  Example:
!
!    Input  Output
!    -----  ------
!    '0'    '1'
!    '1'    '2'
!    ...
!    '8'    '9'
!    '9'    '0'
!    'A'    'A'
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
!    Input/output, character C, a digit to be incremented.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  call ch_to_digit ( c, digit )

  if ( digit == -1 ) then
    return
  end if

  digit = digit + 1

  if ( digit == 10 ) then
    digit = 0
  end if

  call digit_to_ch ( digit, c )

  return
end
subroutine digit_to_ch ( digit, c )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
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
!    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC generates the next filename in a series.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected, and
!    if the name contains no digits, then nothing is done.
!
!  Examples:
!
!      Input          Output
!      -----          ------
!      a7to11.txt     a7to12.txt
!      a7to99.txt     a8to00.txt
!      a9to99.txt     a0to00.txt
!      cat.txt        cat.txt
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  logical ch_is_digit
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( file_name )

  do i = lens, 1, -1

    c = file_name(i:i)

    if ( ch_is_digit ( c ) ) then

      call digit_inc ( c )

      file_name(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do

  return
end
subroutine file_row_count ( file_in_name, row_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of row records in a file.
!
!  Discussion:
!
!    It does not count lines that are blank, or that begin with a
!    comment symbol '#'.
!
!  Modified:
!
!    06 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_IN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
!
  implicit none

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) comment_num
  character ( len = * ) file_in_name
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 100 ) line
  integer ( kind = 4 ) record_num
  integer ( kind = 4 ) row_num

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( file_in_name )
    stop
  end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( file_in_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = record_num
      exit
    end if

    record_num = record_num + 1

    if ( line(1:1) == '#' ) then
      comment_num = comment_num + 1
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      comment_num = comment_num + 1
      cycle
    end if

    row_num = row_num + 1

  end do

  close ( unit = file_in_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FILE_ROW_COUNT:'
  write ( *, '(a,i6)' ) '  Number of records:         ', record_num
  write ( *, '(a,i6)' ) '  Number of data records:    ', row_num
  write ( *, '(a,i6)' ) '  Number of comment records: ', comment_num

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
subroutine s_input ( string, value, ierror )

!*****************************************************************************80
!
!! S_INPUT prints a prompt string and reads a string from the user.
!
!  Discussion:
!
!    If the input line starts with a comment character ('#'), or is blank,
!    the routine ignores that line, and tries to read the next one.
!
!  Modified:
!
!    27 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the prompt string.
!
!    Output, character ( len = * ) VALUE, the value input by the user.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero if no error occurred.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = 80 ) line
  character ( len = * ) string
  character ( len = * ) value

  ierror = 0
  value = ' '
!
!  Write the prompt.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  do

    read ( *, '(a)', iostat = ierror ) line

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'S_INPUT: Fatal error!'
      write ( *, '(a)' ) '  Input error!'
      stop
    end if
!
!  If the line begins with a comment character, go back and read the next line.
!
    if ( line(1:1) == '#' ) then
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    value = line

    exit

  end do

  return
end
subroutine s_to_r8 ( s, r, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Examples:
!
!    S                 R
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
!  Modified:
!
!    12 February 2001
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
!    Output, real ( kind = 8 ) R, the real value that was read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character c
  logical ch_eqi
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
  real ( kind = 8 ) r
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0D+00
  length = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1
    c = s(length+1:length+1)
!
!  Blank or TAB character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( ihave > 1 ) then
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
        isgn = - 1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = - 1
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
      else if ( ihave >= 6 .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Exponent marker.
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
    else if ( ihave < 11 .and. lge ( c, '0' ) .and. lle ( c, '9' ) ) then

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
    if ( iterm == 1 .or. length+1 >= nchar ) then
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

    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0D+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0D+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

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
