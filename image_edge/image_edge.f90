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
subroutine getint ( done, ierror, inunit, ival, string )

!*****************************************************************************80
!
!! GETINT reads an integer from a file.
!
!  Discussion:
!
!    The file, or at least the part read by GETINT, is assumed to
!    contain nothing but integers.  These integers may be separated
!    by spaces, or appear on separate lines.  Comments, which begin
!    with "#" and extend to the end of the line, may appear anywhere.
!
!    Each time GETINT is called, it tries to read the next integer
!    it can find.  It remembers where it was in the current line
!    of text.
!
!    The user should open a text file on FORTRAN unit INUNIT,
!    set STRING = ' ' and DONE = TRUE.  The GETINT routine will take
!    care of reading in a new STRING as necessary, and extracting
!    as many integers as possible from the line of text before 
!    reading in the next line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, logical DONE.
!    On input, if this is the first call, or the user has changed
!    STRING, then set DONE = TRUE.
!    On output, if there is no more data to be read from STRING,
!    then DONE is TRUE.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred while trying to read the integer.
!
!    Input, integer ( kind = 4 ) INUNIT, the FORTRAN unit from which to read.
!
!    Output, integer ( kind = 4 ) IVAL, the integer that was read.
!
!    Input/output, character ( len = * ) STRING, the text of the most recently 
!    read line of the file.
!
  implicit none

  logical done
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) inunit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  character ( len = * ) string
  character ( len = 80 ) word

  do

    call word_next_rd ( string, word, done )

    if ( .not. done ) then
      exit
    end if

    read ( inunit, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      ierror = 1
      return
    end if

    i = index ( string, '#' )
    if ( i /= 0 ) then
      string(i:) = ' '
    end if

  end do

  call s_to_i4 ( word, ival, ierror, last )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GETINT - Fatal error!'
    write ( *, '(a)' ) '  Error trying to convert string to integer.'
    stop
  end if

  return
end
subroutine i4mat_histogram ( m, n, a, histo_num, histo_gram )

!*****************************************************************************80
!
!! I4MAT_HISTOGRAM computes a histogram of the elements of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is an array of I4's.
!
!    It is assumed that the entries in the vector A are nonnegative.
!    Only values between 0 and HISTO_NUM will be histogrammed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the order of A.
!
!    Input, integer ( kind = 4 ) A(M,N), the array to examine.  
!
!    Input, integer ( kind = 4 ) HISTO_NUM, the maximum value for which a
!    histogram entry will be computed.
!
!    Output, integer ( kind = 4 ) HISTO_GRAM(0:HISTO_NUM), contains the 
!    number of entries of A with the values of 0 through HISTO_NUM.
!
  implicit none

  integer ( kind = 4 ) histo_num
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) histo_gram(0:histo_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  histo_gram(0:histo_num) = 0

  do j = 1, n
    do i = 1, m

      if ( 0 <= a(i,j) .and. a(i,j) <= histo_num ) then
        histo_gram(a(i,j)) = histo_gram(a(i,j)) + 1
      end if

    end do
  end do

  return
end
subroutine news ( m, n, a, e )

!*****************************************************************************80
!
!! NEWS demonstrates the NEWS stencil for edge detection.
!
!  Discussion:
!
!    Given a black and white image A, which we regard as an M by N array
!    of pixels, we want to produce an array E of the same shape, which
!    contains information describing the location of edges.
!
!    A simple algorithm for trying to detect edges in an array that
!    represents an image is the NEWS scheme.  For each pixel A(C),
!    we consider its North, East, West, and South pixel neighbors.  The
!    indexing of arrays and images do not correspond, so we will use
!    these directions instead:
!
!             A(N)
!              |
!              |
!      A(W)---A(C)---A(E)
!              |
!              |
!             A(S)
!
!    Entry E(C) of the edge array will be computed by
!
!      E(C) = abs ( A(N) - A(S) ) + abs ( A(E) - A(W) )
!
!    Pixels of A that represent edges will tend to have high values
!    of E, while pixels that are interior to a region of roughly the
!    same shade will tend to have low values.
!
!    Thus, an edge detection scheme would use the NEWS stencil to
!    compute the E array, determine E_MAX, the maximum entry in E, 
!    choose some threshold value E_THRESH, and declare pixel A(I,J) 
!    to be associated with an edge whenever E(I,J) is greater than E_THRESH.
!
!    In this program, we demonstrate the NEWS stencil using a PGM
!    grayscale image of coins.  At the end, we use the edge information
!    to produce a color image in which the edges of the coins have been
!    outlined in red.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in 
!    the image.
!
!    Input, integer ( kind = 4 ) A(M,N), the gray scale image data, presumably
!    integers between 0 and 255.
!
!    Output, integer ( kind = 4 ) E(M,N), is 1 for each pixel that is part of 
!    an edge, and 0 otherwise.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) b(m+2,n+2)
  integer ( kind = 4 ) e(m,n)
  integer ( kind = 4 ) e_max
  integer ( kind = 4 ) thresh
!
!  For neatness, we add a border of zeros to the image,
!  then fill in the border by copying the nearby original values.
!  This will be our M+2 by N+2 data array B.
!
  b(2:m+1,2:n+1) = a(1:m,1:n)

  b(1,  2:n+1) = b(  2,2:n+1)
  b(m+2,2:n+1) = b(m+1,2:n+1)

  b(2:m+1,  1) = b(2:m+1,  2)
  b(2:m+1,n+2) = b(2:m+1,n+1)

  b(  1,  1) = ( b(  1,  2) + b(  2,  1) ) / 2
  b(m+2,  1) = ( b(m+2,  2) + b(m+1,  1) ) / 2
  b(  1,n+2) = ( b(  1,n+1) + b(  2,n+2) ) / 2
  b(m+2,n+2) = ( b(m+2,n+1) + b(m+1,n+2) ) / 2
!
!  Apply the NEWS Operator.  We don't process the boundary pixels.
!
!  The picture is:
!
!   |  0 +1  0 |     |  0  0   0 |
!   |  0  0  0 |  +  | -1  0  +1 |
!   |  0 -1  0 |     |  0  0   0 |
!
  e(1:m,1:n) = abs ( - b(1:m,  2:n+1) + b(3:m+2,2:n+1) ) &
             + abs ( - b(2:m+1,1:n  ) + b(2:m+1,3:n+2) )
!
!  Remap E so the largest value is 255.
!
  e_max = maxval ( e(1:m,1:n) )
!
!  Threshold the data.  Set the threshold to give enough detail 
!  to guess the coin denominations.
!
  thresh = e_max / 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NEWS:'
  write ( *, '(a,i4)' ) '  E_MAX = ', e_max
  write ( *, '(a,i4)' ) '  Using threshold value THRESH = ', thresh

  where ( e < thresh )
    e = 0
  elsewhere
    e = 1
  end where

  return
end
subroutine pbma_write ( file_out_name, row_num, col_num, b )

!*****************************************************************************80
!
!! PBMA_WRITE writes an ASCII portable bit map file.
!
!  Example:
!
!    P1
!    # feep.pbma created by PBMA_IO(PBMA_WRITE).
!    24 7
!    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
!    0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0
!    0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 1 0
!    0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 1 0
!    0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0
!    0 1 0 0 0 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 0 0 0 0
!    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
! 
!  Modified:
!
!    04 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the file 
!    to which the data should be written.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows 
!    and columns of data.
!
!    Input, integer ( kind = 4 ) B(ROW_NUM,COL_NUM), the bit value of each 
!    pixel.  These should be 0 or 1.
!
  implicit none

  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

  integer ( kind = 4 ) b(row_num,col_num)
  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios

  ierror = 0
!
!  Open the file.
!
  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PBMA_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    stop
  end if
!
!  Write the header.
!
  call pbma_write_header ( file_out_name, file_out_unit, row_num, col_num )
!
!  Write the data.
!
  call pbma_write_data ( file_out_unit, row_num, col_num, b )
!
!  Close the file.
!
  close ( unit = file_out_unit )

  return
end
subroutine pbma_write_data ( file_out_unit, row_num, col_num, b )

!*****************************************************************************80
!
!! PBMA_WRITE_DATA writes the data of a PBMA file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FILE_OUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows and 
!    columns of data.
!
!    Input, integer ( kind = 4 ) B(ROW_NUM,COL_NUM), the bit value of each 
!    pixel.  These should be 0 or 1.
!
  implicit none

  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

  integer ( kind = 4 ) b(row_num,col_num)
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
!
!  Write the header.
!
  do i = 1, row_num
    do jlo = 1, col_num, 60
      jhi = min ( jlo + 59, col_num )
      write ( file_out_unit, '(60i1)' ) b(i,jlo:jhi)
    end do
  end do

  return
end
subroutine pbma_write_header ( file_out_name, file_out_unit, row_num, col_num )

!*****************************************************************************80
!
!! PBMA_WRITE_HEADER writes the header of a PBMA file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) FILE_OUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows and 
!    columns of data.
!
  implicit none

  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  character ( len = 2 ) :: magic = 'P1'
  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num
!
!  Write the header.
!
  write ( file_out_unit, '(a2)' ) magic
  write ( file_out_unit, '(a)' ) '# ' // trim ( file_out_name ) &
    // ' created by PBMA_IO::PBMA_WRITE.F90.'
  write ( file_out_unit, '(i8,2x,i8)' ) col_num, row_num

  return
end
subroutine pgma_read_data ( file_in_unit, row_num, col_num, g )

!*****************************************************************************80
!
!! PGMA_READ_DATA reads the data in an ASCII PGM file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FILE_IN_UNIT, the unit number of the file.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows and 
!    columns of data.
!
!    Output, integer ( kind = 4 ) G(ROW_NUM,COL_NUM), the gray data.
!
  implicit none

  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

  logical done
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) g(row_num,col_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  character ( len = 80 ) string

  ierror = 0
  done = .true.
  string = ' '

  do i = 1, row_num
    do j = 1, col_num

      call getint ( done, ierror, file_in_unit, g(i,j), string )
      if ( ierror /= 0 ) then
        close ( unit = file_in_unit )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PGMA_READ_DATA - Fatal error!'
        write ( *, '(a)' ) '  Problem reading G data.'
        stop
      end if

    end do
  end do

  return
end
subroutine pgma_read_header ( file_in_unit, row_num, col_num, g_max )

!*****************************************************************************80
!
!! PGMA_READ_HEADER reads the header of an ASCII PGM file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FILE_IN_UNIT, the unit number of the file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows and 
!    columns of data.
!
!    Output, integer ( kind = 4 ) G_MAX, the maximum gray value.
!
  implicit none

  logical done
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 2 ) magic
  integer ( kind = 4 ) g_max
  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num
  logical s_eqi
  character ( len = 80 ) string
!
!  Read the first line of data, which must begin with the magic number.
!
  read ( file_in_unit, '(a)', iostat = ios ) magic

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
    write ( *, '(a)' ) '  End or error while reading file.'
    ierror = 2
    stop
  end if

  if ( .not. s_eqi ( magic, 'P2' ) ) then
    ierror = 3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error.'
    write ( *, '(a)' ) '  First two bytes are not magic number "P2".'
    write ( *, '(a)' ) '  First two bytes are: "' // magic // '".'
    stop
  end if
!
!  Now search for COL_NUM, ROW_NUM, and G_MAX.
!
  done = .true.
  string = ' '

  call getint ( done, ierror, file_in_unit, col_num, string )

  if ( ierror /= 0 ) then
    close ( unit = file_in_unit )
    ierror = 4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
    write ( *, '(a)' ) '  Problem reading COL_NUM.'
    stop
  end if

  call getint ( done, ierror, file_in_unit, row_num, string )

  if ( ierror /= 0 ) then
    ierror = 4
    close ( unit = file_in_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
    write ( *, '(a)' ) '  Problem reading ROW_NUM.'
    stop
  end if

  call getint ( done, ierror, file_in_unit, g_max, string )

  if ( ierror /= 0 ) then
    ierror = 4
    close ( unit = file_in_unit )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
    write ( *, '(a)' ) '  Problem reading G_MAX.'
    stop
  end if

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
subroutine s_to_i4 ( s, ival, ierror, last )

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
!    28 August 1999
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
!    If blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character that was
!    part of the representation of IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lens
  character ( len = * ) s

  ierror = 0
  istate = 0

  isgn = 1
  ival = 0

  lens = len ( s )

  i = 0

  do

    i = i + 1

    c = s(i:i)

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
        exit
      end if

    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        exit
      end if

    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        istate = 3
      end if

    end if
!
!  Continue or exit?
!
    if ( istate == 3 ) then
      ival = isgn * ival
      last = i - 1
      exit
    else if ( lens <= i ) then
      if ( istate == 2 ) then
        ival = isgn * ival
        last = lens
      else
        ierror = 1
        last = 0
      end if
      exit
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
subroutine word_next_rd ( line, word, done )

!*****************************************************************************80
!
!! WORD_NEXT_RD "reads" words from a string, one at a time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 June 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string, presumably containing
!    words separated by spaces.
!
!    Output, character ( len = * ) WORD.
!    If DONE is FALSE,
!      WORD contains the "next" word read from LINE.
!    Else
!      WORD is blank.
!
!    Input/output, logical DONE.
!    On input, on the first call, or with a fresh value of LINE,
!      set DONE to TRUE.
!    Else
!      leave it at the output value of the previous call.
!    On output, if a new nonblank word was extracted from LINE
!      DONE is FALSE
!    ELSE
!      DONE is TRUE.
!    If DONE is TRUE, then you need to provide a new LINE of data.
!
!  Local Parameters:
!
!    NEXT is the next location in LINE that should be searched.
!
  implicit none

  logical done
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) lenl
  character ( len = * ) line
  integer ( kind = 4 ), save :: next = 1
  character ( len = 1 ), parameter :: TAB = char(9)
  character ( len = * ) word

  lenl = len_trim ( line )

  if ( done ) then
    next = 1
    done = .false.
  end if
!
!  Beginning at index NEXT, search LINE for the next nonblank.
!
  ilo = next

  do
!
!  ...LINE(NEXT:LENL) is blank.  Return with WORD=' ', and DONE=TRUE.
!
    if ( lenl < ilo ) then
      word = ' '
      done = .true.
      next = lenl + 1
      return
    end if
!
!  ...If the current character is blank, skip to the next one.
!
    if ( line(ilo:ilo) /= ' ' .and. line(ilo:ilo) /= TAB ) then
      exit
    end if

    ilo = ilo + 1

  end do
!
!  To get here, ILO must be the index of the nonblank starting
!  character of the next word.
!
!  Now search for the LAST nonblank character.
!
  next = ilo + 1

  do

    if ( lenl < next ) then
      word = line(ilo:next-1)
      return
    end if

    if ( line(next:next) == ' ' .or. line(next:next) == TAB ) then
      exit
    end if

    next = next + 1

  end do

  word = line(ilo:next-1)

  return
end
