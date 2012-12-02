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
subroutine i4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_PRINT prints an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, integer ( kind = 4 ) A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) a(m,n)
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character ( len = * )  title

  ilo = 1
  ihi = m
  jlo = 1
  jhi = n

  call i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

  return
end
subroutine i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_PRINT_SOME prints some of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ), parameter :: incx = 10
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) a(m,n)
  character ( len = 8 )  ctemp(incx)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2hi
  integer   ( kind = 4 ) i2lo
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j2
  integer   ( kind = 4 ) j2hi
  integer   ( kind = 4 ) j2lo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8)' ) j
    end do

    write ( *, '(''  Col '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,a,10a8)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine i4mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! I4MAT_WRITE writes an I4MAT file.
!
!  Discussion:
!
!    An I4MAT is an array of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2009
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
!    Input, integer ( kind = 4 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) j
  character ( len = * )  output_filename
  integer   ( kind = 4 ) output_status
  integer   ( kind = 4 ) output_unit
  character ( len = 30 ) string
  integer   ( kind = 4 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a4)' ) '(', m, 'i10)'
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
subroutine pseudo_euclidean_distance ( n, xy, d )

!*****************************************************************************80
!
!! PSEUDO_EUCLIDEAN_DISTANCE computes the distance used for the ATT example.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 November 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of cities.
!
!    Input, integer ( kind = 4 ) XY(2,N), the coordinates of each city.
!
!    Output, integer ( kind = 4 ) D(N,N), the distance table.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) d(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) rij
  integer ( kind = 4 ) tij
  integer ( kind = 4 ) xy(2,n)

  do i = 1, n
    do j = 1, n
      rij = sqrt ( real ( ( xy(1,i) - xy(1,j) )**2 &
                        + ( xy(2,i) - xy(2,j) )**2 ) )
      tij = nint ( rij )
      if ( tij < rij ) then
        d(i,j) = tij + 1
      else
        d(i,j) = tij
      end if
    end do
  end do

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.  
!
!  Discussion:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
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
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  character              c1
  character              c2
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) lenc
  logical                s_eqi
  character ( len = *  ) s1
  integer   ( kind = 4 ) s1_length
  character ( len = *  ) s2
  integer   ( kind = 4 ) s2_length

  s1_length = len ( s1 )
  s2_length = len ( s2 )
  lenc = min ( s1_length, s2_length )
 
  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do
 
  do i = lenc + 1, s1_length
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do
 
  do i = lenc + 1, s2_length
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do
 
  s_eqi = .true.
 
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

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) length
  character ( len = * )  s
  integer   ( kind = 4 ) state
  character              :: TAB = achar ( 9 )
  integer   ( kind = 4 ) value

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

      if ( c == ' ' .or. c == TAB ) then

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

      if ( c == ' ' .or. c == TAB ) then

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
subroutine s_word_extract_first ( s, w )

!*****************************************************************************80
!
!! S_WORD_EXTRACT_FIRST extracts the first word from a string.
!
!  Discussion:
!
!    A "word" is a string of characters terminated by a blank or
!    the end of the string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string.  On output, the first
!    word has been removed, and the remaining string has been shifted left.
!
!    Output, character ( len = * ) W, the leading word of the string.
!
  implicit none

  integer   ( kind = 4 ) get1
  integer   ( kind = 4 ) get2
  character ( len = * )  s
  integer   ( kind = 4 ) s_length
  character ( len = * )  w

  w = ' '

  s_length = len_trim ( s )

  if ( s_length < 1 ) then
    return
  end if
!
!  Find the first nonblank.
!
  get1 = 0

  do

    get1 = get1 + 1

    if ( s_length < get1 ) then
      return
    end if

    if ( s(get1:get1) /= ' ' ) then
      exit
    end if

  end do
!
!  Look for the last contiguous nonblank.
!
  get2 = get1

  do

    if ( s_length <= get2 ) then
      exit
    end if

    if ( s(get2+1:get2+1) == ' ' ) then
      exit
    end if

    get2 = get2 + 1

  end do
!
!  Copy the word.
!
  w = s(get1:get2)
!
!  Shift the string.
!
  s(1:get2) = ' '
  s = adjustl ( s )

  return
end
subroutine text_split ( text, keyword, colon, field )

!*****************************************************************************80
!
!! TEXT_SPLIT uses a colon to locate the keyword and field in a line of text.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) TEXT, a line of text.
!
!    Output, character ( len = * ) KEYWORD, the portion of TEXT that
!    precedes the first colon.
!
!    Output, character ( len = 1 ) COLON, the colon, if it occurred.
!    Otherwise, a blank.
!
!    Output, character ( len = * ) FIELD< the portion of TEXT that follows
!    the first colon.
!
  implicit none

  character ( len = 1 )  colon
  integer   ( kind = 4 ) colon_loc
  character ( len = * )  field
  character ( len = * )  keyword
  character ( len = * )  text
  integer   ( kind = 4 ) text_end

  colon_loc = index ( text, ':' )

  if ( colon_loc == 0 ) then
    keyword = text
    colon = ' '
    field = ' '
    return
  end if

  keyword = text(1:colon_loc-1)
  keyword = adjustl ( keyword )

  colon = ':'

  text_end = len_trim ( text )
  field = text(colon_loc+1:text_end)
  field = adjustl ( field )

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

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
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
subroutine tsp_edge_weight_read ( tsp_filename, tsp_dimension, &
  tsp_edge_weight_type, tsp_edge_weight_format, tsp_edge_weight )

!*****************************************************************************80
!
!! TSP_EDGE_WEIGHT_READ reads the edge weights of a TSP file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) TSP_FILENAME, the TSP filename.
!
!    Input, integer ( kind = 4 ) TSP_DIMENSION, the number of cities.
!
!    Input, character ( len = * ) TSP_EDGE_WEIGHT_TYPE:
!    'EXPLICIT': the edge weights are given as an explicit list.
!
!    Input, character ( len = * ) TSP_EDGE_WEIGHT_FORMAT:
!    'FULL_MATRIX': the full distance matrix is listed.
!    'LOWER_DIAG_ROW': the diagonal and lower triangle are listed.
!
!    Output, integer ( kind = 4 ) TSP_EDGE_WEIGHT(TSP_DIMENSION,TSP_DIMENSION),
!    the edge weights.
!
  implicit none

  integer ( kind = 4 )    tsp_dimension

  logical                 done
  logical                 found_keyword
  integer ( kind = 4 )    i
  integer ( kind = 4 )    ierror
  integer ( kind = 4 )    j
  integer ( kind = 4 )    length
  integer ( kind = 4 )    line_num
  logical                 matrix_full
  logical                 s_eqi
  character ( len = 255 ) text
  integer ( kind = 4 )    tsp_edge_weight(tsp_dimension,tsp_dimension)
  character ( len = * )   tsp_edge_weight_format
  character ( len = * )   tsp_edge_weight_type
  character ( len = * )   tsp_filename
  integer ( kind = 4 ), allocatable :: tsp_node_coord(:,:)
  integer ( kind = 4 )    tsp_unit
  integer ( kind = 4 )    value
  character ( len = 255 ) word

  call get_unit ( tsp_unit )

  open ( file = tsp_filename, unit = tsp_unit, status = 'old' )

  found_keyword = .false.
  line_num = 0
  matrix_full = .false.

  if ( s_eqi ( tsp_edge_weight_type, 'EXPLICIT' ) ) then

    do
!
!  Read a line of text.
!
      read ( tsp_unit, '(a)' ) text

      line_num = line_num + 1

      if ( .not. found_keyword ) then

        if ( s_eqi ( text, 'EDGE_WEIGHT_SECTION' ) ) then
          found_keyword = .true.
          i = 1
          j = 1
        end if

      else
!
!  Reading a line of weights.  Read as many numbers as you can from this line.
!
        done = .true.

        do

          call s_word_extract_first ( text, word )

          call s_to_i4 ( word, value, ierror, length )

          if ( length <= 0 ) then
            exit
          end if

          if ( s_eqi ( tsp_edge_weight_format, 'FULL_MATRIX' ) ) then

            tsp_edge_weight(i,j) = value

            i = i + 1
            if ( tsp_dimension < i ) then
              i = 1
              j = j + 1
              if ( tsp_dimension < j ) then
                matrix_full = .true.
                exit
              end if
            end if

          else if ( s_eqi ( tsp_edge_weight_format, 'LOWER_DIAG_ROW' ) ) then

            tsp_edge_weight(i,j) = value
            tsp_edge_weight(j,i) = value

            j = j + 1
            if ( i < j ) then
              j = 1
              i = i + 1
              if ( tsp_dimension < i ) then
                matrix_full = .true.
                exit
              end if
            end if

          else

            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'TSP_EDGE_WEIGHT_READ - Fatal error!'
            write ( *, '(a)' ) '  Do not know how to handle format:'
            write ( *, '(a)' ) '  "' // trim ( tsp_edge_weight_format ) // '".'
            stop

          end if

        end do

        if ( matrix_full ) then
          exit 
        end if
        
      end if       

    end do
!
!  ATT Edge Weight Specification.
!
  else if ( s_eqi ( tsp_edge_weight_type, 'ATT' ) ) then

    allocate ( tsp_node_coord(2,tsp_dimension) )

    do
!
!  Read a line of text.
!
      read ( tsp_unit, '(a)' ) text

      line_num = line_num + 1

      if ( .not. found_keyword ) then

        if ( s_eqi ( text, 'NODE_COORD_SECTION' ) ) then
          found_keyword = .true.
          i = 1
          j = 1
        end if

      else
!
!  Expect to read a line of (index, X, Y).
!
        done = .true.
        call s_word_extract_first ( text, word )
        call s_to_i4 ( word, value, ierror, length )

        call s_word_extract_first ( text, word )
        call s_to_i4 ( word, value, ierror, length )
        tsp_node_coord(1,j) = value

        call s_word_extract_first ( text, word )
        call s_to_i4 ( word, value, ierror, length )
        tsp_node_coord(2,j) = value

        if ( tsp_dimension <= j ) then
          exit 
        end if

        j = j + 1
        
      end if       

    end do

    call pseudo_euclidean_distance ( tsp_dimension, tsp_node_coord, &
      tsp_edge_weight )

    deallocate ( tsp_node_coord )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TSP_EDGE_WEIGHT_READ - Fatal error!'
    write ( *, '(a)' ) '  Do not know how to handle given type:'
    write ( *, '(a)' ) '  "' // trim ( tsp_edge_weight_type ) // '".'
    stop

  end if

  close ( unit = tsp_unit )

  return
end
subroutine tsp_header_read ( tsp_filename, tsp_name, tsp_type, tsp_comment, &
  tsp_dimension, tsp_capacity, tsp_edge_weight_type, tsp_edge_weight_format, &
  tsp_edge_data_format, tsp_node_coord_type, tsp_display_data_type )

!*****************************************************************************80
!
!! TSP_HEADER_READ reads the header of a TSP file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) TSP_FILENAME, the TSP filename.
!
  implicit none

  character ( len = 1 )   colon
  integer ( kind = 4 )    colon_loc
  character ( len = 255 ) field
  integer ( kind = 4 )    ierror
  character ( len = 255 ) keyword
  integer ( kind = 4 )    leng
  logical                 s_eqi
  character ( len = 255 ) text
  integer ( kind = 4 )    tsp_capacity
  character ( len = * )   tsp_comment
  integer ( kind = 4 )    tsp_dimension
  character ( len = * )   tsp_display_data_type
  character ( len = * )   tsp_edge_data_format
  character ( len = * )   tsp_edge_weight_format
  character ( len = * )   tsp_edge_weight_type
  character ( len = * )   tsp_filename
  integer ( kind = 4 )    tsp_header_lines
  character ( len = * )   tsp_name
  character ( len = * )   tsp_node_coord_type
  character ( len = * )   tsp_type
  integer ( kind = 4 )    tsp_unit
  logical, parameter ::   verbose = .false.
!
!  Initialize.
!
  tsp_capacity = 0
  tsp_comment = ' '
  tsp_dimension = 0
  tsp_display_data_type = 'DEFAULT'
  tsp_edge_data_format = ' '
  tsp_edge_weight_format = ' '
  tsp_edge_weight_type = ' '
  tsp_name = ' '
  tsp_node_coord_type = 'DEFAULT'
  tsp_type = ' '

  call get_unit ( tsp_unit )

  open ( file = tsp_filename, unit = tsp_unit, status = 'old' )

  tsp_header_lines = 0

  do
!
!  Read a line of text.
!
    read ( tsp_unit, '(a)' ) text
!
!  Split text into "KEYWORD : FIELD"
!
    call text_split ( text, keyword, colon, field )

    if ( colon /= ':' ) then
      exit
    end if

    tsp_header_lines = tsp_header_lines + 1
!
!  Store the field in the variable indicated by the keyword.
!
    if ( s_eqi ( keyword, 'NAME' ) ) then
      tsp_name = field
    else if ( s_eqi ( keyword, 'TYPE' ) ) then
      tsp_type = field
    else if ( s_eqi ( keyword, 'COMMENT' ) ) then
      tsp_comment = field
    else if ( s_eqi ( keyword, 'DIMENSION' ) ) then
      call s_to_i4 ( field, tsp_dimension, ierror, leng )
    else if ( s_eqi ( keyword, 'CAPACITY' ) ) then
      call s_to_i4 ( field, tsp_capacity, ierror, leng )
    else if ( s_eqi ( keyword, 'EDGE_WEIGHT_TYPE' ) ) then
      tsp_edge_weight_type = field
    else if ( s_eqi ( keyword, 'EDGE_WEIGHT_FORMAT' ) ) then
      tsp_edge_weight_format = field
    else if ( s_eqi ( keyword, 'EDGE_DATA_FORMAT' ) ) then
      tsp_edge_data_format = field
    else if ( s_eqi ( keyword, 'NODE_COORD_TYPE' ) ) then
      tsp_node_coord_type = field
    else if ( s_eqi ( keyword, 'DISPLAY_DATA_TYPE' ) ) then
      tsp_display_data_type = field
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TSP_HEADER_READ - Fatal error!'
      write ( *, '(a)' ) '  Unexpected keyword "' // trim ( keyword ) // '".'
      stop
    end if

  end do

  close ( unit = tsp_unit )
!
!  In some cases, if a variable was not defined, it has a default value.
!
  if ( tsp_node_coord_type == 'DEFAULT' ) then
    tsp_node_coord_type = 'NO_COORDS'
  end if

  if ( tsp_display_data_type == 'DEFAULT' ) then
    if ( s_eqi ( tsp_node_coord_type, 'TWOD_COORDS' ) .or. &
         s_eqi ( tsp_node_coord_type, 'THREED_COORDS' ) ) then
      tsp_display_data_type = 'COORD_DISPLAY'
    else
      tsp_display_data_type = 'NO_DISPLAY'
    end if
  end if

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TSP_HEADER_READ:'
    write ( *, '(a)' ) '  Closed "' // trim ( tsp_filename ) // '".'
    write ( *, '(a,i8)' ) '  Number of header lines is ', tsp_header_lines
  end if

  return
end
