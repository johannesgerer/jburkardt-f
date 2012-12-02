program main

!*****************************************************************************80
!
!! MAIN is the main program for VAN_DER_CORPUT_DATASET.
!
!  Discussion:
!
!    VAN_DER_CORPUT_DATASET generates a van der Corput dataset.
!
!  Usage:
!
!    van_der_corput_dataset base seed n
!
!    where
!
!    * BASE, the base of the sequence;
!    * SEED, the index of the first element to be computed;
!    * N, the number of points to generate.
!
!    The program generates the data and writes it to the file
!
!      van_der_corput_BASE_SEED_N.txt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Johannes van der Corput,
!    Verteilungsfunktionen I & II,
!    Nederl. Akad. Wetensch. Proc.,
!    Volume 38, 1935, pages 813-820, pages 1058-1066.
!
  implicit none

  integer ( kind = 4 ) arg_num
  integer ( kind = 4 ) base
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) last
  integer ( kind = 4 ) n
  character ( len = 255 ) output_filename
  real ( kind = 8 ), allocatable, dimension ( : ) :: r
  integer ( kind = 4 ) seed
  character ( len = 255 ) string

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VAN_DER_CORPUT_DATASET'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Generate a van der Corput dataset.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program requests input values from the user:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  * BASE, the base;'
  write ( *, '(a)' ) '  * SEED, the index of the first element;'
  write ( *, '(a)' ) '  * N, the number of elements to generate,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program generates the data and writes it to'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    van_der_corput_BASE_SEED_N.txt'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get the base BASE.
!
  if ( 1 <= arg_num ) then
    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, base, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter BASE, the van der Corput base,'
    write ( *, '(a)' ) '  often a prime number.'
    write ( *, '(a)' ) '  (Try "2" if you do not have a preference).'
    read ( *, * ) base
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)') '  User input BASE = ', base
!
!  Get the seed SEED.
!
  if ( 2 <= arg_num ) then
    iarg = 2
    call getarg ( iarg, string )
    call s_to_i4 ( string, seed, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter SEED, the index of the first element.'
    write ( *, '(a)' ) '  (Try "0" or "1" if you do not have a preference).'
    read ( *, * ) seed
  end if

  write ( *, '(a,i12)') '  User input SEED = ', seed
!
!  Get the number of points N.
!
  if ( 3 <= arg_num ) then
    iarg = 3
    call getarg ( iarg, string )
    call s_to_i4 ( string, n, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter N, the number of points to generate.'
    write ( *, '(a)' ) '  (Try "25" if you do not have a preference).'
    read ( *, * ) n
  end if

  write ( *, '(a,i12)') '  User input N = ', n
!
!  Generate the data.
!
  allocate ( r(1:n) )

  call i4_to_van_der_corput_sequence ( seed, base, n, r )
!
!  Write the data to a file.
!
  write ( output_filename, '(a,i12,a,i12,a,i12,a)' ) &
    'van_der_corput_', base, '_', seed, '_', n, '.txt'

  call s_blank_delete ( output_filename )

  call r8mat_write ( output_filename, 1, n, r )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data was written to the file "' &
     // trim ( output_filename ) // '".'
!
!  Deallocate memory.
!
  deallocate ( r )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VAN_DER_CORPUT_DATASET:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
subroutine i4_to_van_der_corput_sequence ( seed, base, n, r )

!*****************************************************************************80
!
!! I4_TO_VAN_DER_CORPUT_SEQUENCE: next N elements of a van der Corput sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Halton,
!    On the efficiency of certain quasi-random sequences of points
!    in evaluating multi-dimensional integrals,
!    Numerische Mathematik,
!    Volume 2, pages 84-90, 1960.
!
!    Johnannes van der Corput,
!    Verteilungsfunktionen I & II,
!    Nederl. Akad. Wetensch. Proc.,
!    Volume 38, 1935, pages 813-820, pages 1058-1066.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED, the seed or index of the desired element.
!    SEED should be nonnegative.
!    SEED = 0 is allowed, and returns R = 0.
!
!    Input, integer ( kind = 4 ) BASE, the van der Corput base, which is typically
!    a prime number.  BASE must be greater than 1.
!
!    Input, integer ( kind = 4 ) N, the number of elements desired.
!
!    Output, real ( kind = 8 ) R(N), the SEED-th through (SEED+N-1)-th
!    elements of the van der Corput sequence for base BASE.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) base
  real ( kind = 8 ) base_inv
  integer ( kind = 4 ) digit(n)
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed2(n)

  if ( base <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_VAN_DER_CORPUT_SEQUENCE - Fatal error!'
    write ( *, '(a)' ) '  The input base BASE is <= 1!'
    write ( *, '(a,i6)' ) '  BASE = ', base
    stop
  end if

  if ( seed < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_VAN_DER_CORPUT_SEQUENCE - Fatal error!'
    write ( *, '(a)' ) '  The input base SEED is < 0!'
    write ( *, '(a,i6)' ) '  SEED = ', seed
    stop
  end if
!
!  Set SEED2 = (/ SEED, SEED+1, SEED+2, ..., SEED+N-1 /)
!
  call i4vec_indicator ( n, seed2 )

  seed2(1:n) = seed2(1:n) + seed - 1

  base_inv = 1.0D+00 / real ( base, kind = 8 )

  r(1:n) = 0.0D+00

  do while ( any ( seed2(1:n) /= 0 ) )
    digit(1:n) = mod ( seed2(1:n), base )
    r(1:n) = r(1:n) + real ( digit(1:n), kind = 8 ) * base_inv
    base_inv = base_inv / real ( base, kind = 8 )
    seed2(1:n) = seed2(1:n) / base
  end do

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector A(I)=I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = i
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
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) get
  integer ( kind = 4 ) put
  integer ( kind = 4 ) nchar
  character ( len = * )  s
  character, parameter :: TAB = char ( 9 )

  put = 0
  nchar = len_trim ( s )

  do get = 1, nchar

    c = s(get:get)

    if ( c /= ' ' .and. c /= TAB ) then
      put = put + 1
      s(put:put) = c
    end if

  end do

  s(put+1:nchar) = ' '

  return
end
subroutine s_to_i4 ( s, value, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an integer value from a string.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Licensing:
!
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) VALUE, the integer value read from the string.
!    If the string is blank, then VALUE will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters
!    of S used to make the integer.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) length
  character ( len = * ) s
  integer ( kind = 4 ) state
  integer ( kind = 4 ) value

  value = 0
  ierror = 0
  length = 0

  state = 0
  isgn = 1

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  STATE = 0, haven't read anything.
!
    if ( state == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        state = 1
        isgn = -1
      else if ( c == '+' ) then
        state = 1
        isgn = +1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 1, have read the sign, expecting digits or spaces.
!
    else if ( state == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 2, have read at least one digit, expecting more.
!
    else if ( state == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

        value = 10 * value + iachar ( c ) - iachar ( '0' )

      else

        value = isgn * value
        ierror = 0
        length = i - 1
        return

      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( state == 2 ) then

    value = isgn * value
    ierror = 0
    length = len_trim ( s )

  else

    value = 0
    ierror = 1
    length = 0

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
