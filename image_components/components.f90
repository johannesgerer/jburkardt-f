subroutine file_column_count ( input_filename, column_num )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of COLUMN_NUM words,
!    separated by spaces.  There may also be some blank lines, and some
!    comment lines,
!    which have a "#" in column 1.
!
!    The routine tries to find the first non-comment non-blank line and
!    counts the number of words in that line.
!
!    If all lines are blanks or comments, it goes back and tries to analyze
!    a comment line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the file.
!
!    Output, integer ( kind = 4 ) COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer ( kind = 4 ) column_num
  logical got_one
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  character ( len = 255 ) line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = input_status )

  if ( input_status /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' &
      // trim ( input_filename ) // '" on unit ', input_unit
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    got_one = .true.
    exit

  end do

  if ( .not. got_one ) then

    rewind ( input_unit )

    do

      read ( input_unit, '(a)', iostat = input_status ) line

      if ( input_status /= 0 ) then
        exit
      end if

      if ( len_trim ( line ) == 0 ) then
        cycle
      end if

      got_one = .true.
      exit

    end do

  end if

  close ( unit = input_unit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    column_num = -1
    return
  end if

  call s_word_count ( line, column_num )

  return
end
subroutine file_row_count ( input_filename, row_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of row records in a file.
!
!  Discussion:
!
!    It does not count lines that are blank, or that begin with a
!    comment symbol '#'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
!
  implicit none

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) comment_num
  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  character ( len = 255 ) line
  integer ( kind = 4 ) record_num
  integer ( kind = 4 ) row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
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

  close ( unit = input_unit )

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
subroutine i4block_components ( l, m, n, a, component_num, c )

!*****************************************************************************80
!
!! I4BLOCK_COMPONENTS assigns contiguous nonzero pixels to a common component.
!
!  Discussion:
!
!    On input, the A array contains values of 0 or 1.
!
!    The 0 pixels are to be ignored.  The 1 pixels are to be grouped
!    into connected components.
!
!    The pixel A(I,J,K) is "connected" to the pixels:
!
!      A(I-1,J,  K  ),  A(I+1,J,  K  ),
!      A(I,  J-1,K  ),  A(I,  J+1,K  ),
!      A(I,  J,  K-1),  A(I,  J,  K+1),
!
!    so most pixels have 6 neighbors.
!
!    On output, COMPONENT_NUM reports the number of components of nonzero
!    data, and the array C contains the component assignment for
!    each nonzero pixel, and is 0 for zero pixels.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, M, N, the order of the array.
!
!    Input, integer ( kind = 4 ) A(L,M,N), the pixel array.
!
!    Output, integer ( kind = 4 ) COMPONENT_NUM, the number of components
!    of nonzero data.
!
!    Output, integer ( kind = 4 ) C(L,M,N), the component array.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(l,m,n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c(l,m,n)
  integer ( kind = 4 ) c1
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) component
  integer ( kind = 4 ) component_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) north
  integer ( kind = 4 ) p(0:l*m*n)
  integer ( kind = 4 ) q(0:l*m*n)
  integer ( kind = 4 ) up
  integer ( kind = 4 ) west
!
!  Initialization.
!
  c(1:l,1:m,1:n) = 0
  component_num = 0
!
!  P is simply used to store the component labels.  The dimension used
!  here is, of course, usually an absurd overestimate.
!
  do i = 0, l * m * n
    p(i) = i
  end do
!
!  "Read" the array one pixel at a time.  If a (nonzero) pixel's north or
!  west neighbor already has a label, the current pixel inherits it.
!  In case the labels disagree, we need to adjust the P array so we can
!  later deal with the fact that the two labels need to be merged.
!
  do i = 1, l

    do j = 1, m

      do k = 1, n

        if ( i == 1 ) then
          north = 0
        else
          north = c(i-1,j,k)
        end if

        if ( j == 1 ) then
          west = 0
        else
          west = c(i,j-1,k)
        end if

        if ( k == 1 ) then
          up = 0
        else
          up = c(i,j,k-1)
        end if

        if ( a(i,j,k) /= 0 ) then
!
!  New component?
!
          if ( north == 0 .and. west == 0 .and. up == 0 ) then
            component_num = component_num + 1
            c(i,j,k) = component_num
!
!  One predecessor is labeled.
!
          else if ( north /= 0 .and. west == 0 .and. up == 0 ) then
            c(i,j,k) = north
          else if ( north == 0 .and. west /= 0 .and. up == 0 ) then
            c(i,j,k) = west
          else if ( north == 0 .and. west == 0 .and. up /= 0 ) then
            c(i,j,k) = up
!
!  Two predecessors are labeled.
!
          else if ( north == 0 .and. west /= 0 .and. up /= 0 ) then
            c(i,j,k) = min ( west, up )
            c1 = min ( p(west), p(up) )
            p(west) = c1
            p(up) = c1
          else if ( north /= 0 .and. west == 0 .and. up /= 0 ) then
            c(i,j,k) = min ( north, up )
            c1 = min ( p(north), p(up) )
            p(north) = c1
            p(up) = c1
          else if ( north /= 0 .and. west /= 0 .and. up == 0 ) then
            c(i,j,k) = min ( north, west )
            c1 = min ( p(north), p(west) )
            p(north) = c1
            p(west) = c1
!
!  Three predecessors are labeled.
!
          else if ( north /= 0 .and. west /= 0 .and. up /= 0 ) then
            c(i,j,k) = min ( north, west, up )
            c1 = min ( p(north), p(west), p(up) )
            p(north) = c1
            p(west) = c1
            p(up) = c1
          end if

        end if

      end do

    end do

  end do
!
!  When a component has multiple labels, have the higher labels
!  point to the lowest one.
!
  do component = component_num, 1, -1
    b = component
    do while ( p(b) /= b )
      b = p(b)
    end do
    p(component) = b
  end do
!
!  Locate the minimum label for each component.
!  Assign these mininum labels new consecutive indices.
!
  q(0:component_num) = 0
  i = 0
  do component = 1, component_num
    if ( p(component) == component ) then
      i = i + 1
      q(component) = i
    end if
  end do

  component_num = i
!
!  Replace the labels by consecutive labels.
!
  do i = 1, l
    do j = 1, m
      do k = 1, n
        c(i,j,k) = q ( p ( c(i,j,k) ) )
      end do
    end do
  end do

  return
end
subroutine i4mat_components ( m, n, a, component_num, c )

!*****************************************************************************80
!
!! I4MAT_COMPONENTS assigns contiguous nonzero pixels to a common component.
!
!  Discussion:
!
!    On input, the A array contains values of 0 or 1.
!
!    The 0 pixels are to be ignored.  The 1 pixels are to be grouped
!    into connected components.
!
!    The pixel A(I,J) is "connected" to the pixels A(I-1,J), A(I+1,J),
!    A(I,J-1) and A(I,J+1), so most pixels have 4 neighbors.
!
!    (Another choice would be to assume that a pixel was connected
!    to the other 8 pixels in the 3x3 block containing it.)
!
!    On output, COMPONENT_NUM reports the number of components of nonzero
!    data, and the array C contains the component assignment for
!    each nonzero pixel, and is 0 for zero pixels.
!
!  Picture:
!
!    Input A:
!
!      0  2  0  0 17  0  3
!      0  0  3  0  1  0  4
!      1  0  4  8  8  0  7
!      3  0  6 45  0  0  0
!      3 17  0  5  9  2  5
!
!    Output:
!
!      COMPONENT_NUM = 4
!
!      C:
!
!      0  1  0  0  2  0  3
!      0  0  2  0  2  0  3
!      4  0  2  2  2  0  3
!      4  0  2  2  0  0  0
!      4  4  0  2  2  2  2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of the array.
!
!    Input, integer ( kind = 4 ) A(M,N), the pixel array.
!
!    Output, integer ( kind = 4 ) COMPONENT_NUM, the number of components
!    of nonzero data.
!
!    Output, integer ( kind = 4 ) C(M,N), the component array.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c(m,n)
  integer ( kind = 4 ) component
  integer ( kind = 4 ) component_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) north
  integer ( kind = 4 ) p(0:m*n)
  integer ( kind = 4 ) q(0:m*n)
  integer ( kind = 4 ) west
!
!  Initialization.
!
  c(1:m,1:n) = 0
  component_num = 0
!
!  P is simply used to store the component labels.  The dimension used
!  here is, of course, usually an absurd overestimate.
!
  do i = 0, m * n
    p(i) = i
  end do
!
!  "Read" the array one pixel at a time.  If a (nonzero) pixel's north or
!  west neighbor already has a label, the current pixel inherits it.
!  In case the labels disagree, we need to adjust the P array so we can
!  later deal with the fact that the two labels need to be merged.
!
  do i = 1, m

    do j = 1, n

      if ( i == 1 ) then
        north = 0
      else
        north = c(i-1,j)
      end if

      if ( j == 1 ) then
        west = 0
      else
        west = c(i,j-1)
      end if

      if ( a(i,j) /= 0 ) then

        if ( north == 0 ) then

          if ( west == 0 ) then
            component_num = component_num + 1
            c(i,j) = component_num
          else
            c(i,j) = west
          end if

        else if ( north /= 0 ) then

          if ( west == 0 .or. west == north ) then
            c(i,j) = north
          else
            c(i,j) = min ( north, west )
            if ( north < west ) then
              p(west) = north
            else
              p(north) = west
            end if
          end if

        end if

      end if

    end do

  end do
!
!  When a component has multiple labels, have the higher labels
!  point to the lowest one.
!
  do component = component_num, 1, -1
    b = component
    do while ( p(b) /= b )
      b = p(b)
    end do
    p(component) = b
  end do
!
!  Locate the minimum label for each component.
!  Assign these mininum labels new consecutive indices.
!
  q(0:component_num) = 0
  i = 0
  do component = 1, component_num
    if ( p(component) == component ) then
      i = i + 1
      q(component) = i
    end if
  end do

  component_num = i
!
!  Replace the labels by consecutive labels.
!
  do i = 1, m
    do j = 1, n
      c(i,j) = q ( p ( c(i,j) ) )
    end do
  end do

  return
end
subroutine i4mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! I4MAT_DATA_READ reads data from an I4MAT file.
!
!  Discussion:
!
!    An I4MAT is an array of I4's.
!
!    The file may contain more than N points, but this routine
!    will return after reading N points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, integer ( kind = 4 ) TABLE(M,N), the data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) j
  character ( len = 255 ) line
  integer ( kind = 4 ) table(m,n)
  integer ( kind = 4 ) x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      ierror = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I4MAT_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i8)' ) '  Number of values expected per line M = ', m
      write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_i4vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

  end do

  close ( unit = input_unit )

  return
end
subroutine i4mat_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! I4MAT_HEADER_READ reads the header from an I4MAT.
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
!    04 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) M, spatial dimension.
!
!    Output, integer ( kind = 4 ) N, the number of points.
!
  implicit none

  character ( len = * ) input_filename
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  call file_column_count ( input_filename, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  return
end
subroutine i4vec_components ( n, a, component_num, c )

!*****************************************************************************80
!
!! I4VEC_COMPONENTS assigns contiguous nonzero pixels to a common component.
!
!  Discussion:
!
!    This calculation is trivial compared to the 2D problem, and is included
!    primarily for comparison.
!
!    On input, the A array contains values of 0 or 1.
!
!    The 0 pixels are to be ignored.  The 1 pixels are to be grouped
!    into connected components.
!
!    The pixel A(I) is "connected" to the pixels A(I-1) and A(I+1).
!
!    On output, COMPONENT_NUM reports the number of components of nonzero
!    data, and the array C contains the component assignment for
!    each nonzero pixel, and is 0 for zero pixels.
!
!  Picture:
!
!    Input A:
!
!      0 0 1 2 4 0 0 4 0 0 0 8 9 9 1 2 3 0 0 5 0 1 6 0 0 0 4 0
!
!    Output:
!
!      COMPONENT_NUM = 6
!
!      C:
!
!      0 0 1 1 1 0 0 2 0 0 0 3 3 3 3 3 3 0 0 4 0 5 5 0 0 0 6 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the pixel array.
!
!    Output, integer ( kind = 4 ) COMPONENT_NUM, the number of components
!    of nonzero data.
!
!    Output, integer ( kind = 4 ) C(N), the component array.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) component_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) west
!
!  Initialization.
!
  c(1:n) = 0
  component_num = 0
!
!  "Read" the array one pixel at a time.  If a (nonzero) pixel's west neighbor
!  already has a label, the current pixel inherits it.  Otherwise, we have
!  begun a new component.
!
  west = 0

  do j = 1, n

    if ( a(j) /= 0 ) then
      if ( west == 0 ) then
        component_num = component_num + 1
      end if
      c(j) = component_num
    end if

    west = c(j)

  end do

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
!    28 June 2000
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
subroutine s_to_i4vec ( s, n, ivec, ierror )

!*****************************************************************************80
!
!! S_TO_I4VEC reads an I4VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer ( kind = 4 ) N, the number of values expected.
!
!    Output, integer ( kind = 4 ) IVEC(N), the values read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ivec(n)
  integer ( kind = 4 ) length
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_i4 ( s(ilo:), ivec(i), ierror, length )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + length

  end do

  return
end
subroutine s_word_count ( s, nword )

!*****************************************************************************80
!
!! S_WORD_COUNT counts the number of "words" in a string.
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
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, integer ( kind = 4 ) NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens
  integer ( kind = 4 ) nword
  character ( len = * ) s

  nword = 0
  lens = len ( s )

  if ( lens <= 0 ) then
    return
  end if

  blank = .true.

  do i = 1, lens

    if ( s(i:i) == ' ' ) then
      blank = .true.
    else if ( blank ) then
      nword = nword + 1
      blank = .false.
    end if

  end do

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
