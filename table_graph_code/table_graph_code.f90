program main

!*****************************************************************************80
!
!! MAIN is the main program for TABLE_GRAPH_CODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Usage:
!
!    table_graph_code input_file_name
!
  implicit none

  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: a
  integer   ( kind = 4 ) arg_num
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: code
  integer   ( kind = 4 ) iarg
  integer   ( kind = 4 ) iargc
  character ( len = 255 ) :: input_file_name
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ), allocatable, dimension ( : ) :: order

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TABLE_GRAPH_CODE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  For integer data stored in a TABLE file,'
  write ( *, '(a)' ) '  ITABLE_HEADER_READ reads the header information'
  write ( *, '(a)' ) '  (about the dimensions of the data);'
  write ( *, '(a)' ) '  ITABLE_DATA_READ reads the data.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get the input file name.
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, input_file_name )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TABLE_GRAPH_CODE:'
    write ( *, '(a)' ) '  Please enter the name of the input file.'

    read ( *, '(a)' ) input_file_name

  end if
!
!  Read the input file headere.
!
  call itable_header_read (  input_file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the header of "' // trim ( input_file_name ) // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Row dimension M =     ', m
  write ( *, '(a,i8)' ) '  Column dimension N  = ', n

  if ( m /= n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Fatal error!'
    write ( *, '(a)' ) '  The input matrix is not square.'
    return
  end if
!
!  Allocate data.
!
  allocate ( a(1:n,1:n) )
  allocate ( code(1:n,1:n) )
  allocate ( order(1:n) )
!
!  Read the input file data.
!
  call itable_data_read (  input_file_name, n, n, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the data in "' // trim ( input_file_name ) // '".'

  call i4mat_print_some ( n, n, a, 1, 1, 5, 5, &
    '  5 x 5 portion of data read from file:' )
!
!  Compute the matrix code.
!
  call mg_code_back ( a, n, code, order )
!
!  Print the results.
!
  call node_order_print ( n, order, '  Node ordering:' )

  call mg_code_print ( n, code, '  The order-dependent code:' )
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( code )
  deallocate ( order )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TABLE_GRAPH_CODE'
  write ( *, '(a)' ) '  Normal end of execution.'

  call timestamp ( )

  stop
end
subroutine itable_data_read ( input_file_name, m, n, table )

!*****************************************************************************80
!
!! ITABLE_DATA_READ reads data from an ITABLE file.
!
!  Discussion:
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
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) ( kind = 4 ) N, the number of points.  
!
!    Output, integer ( kind = 4 ) ( kind = 4 ) TABLE(M,N), the table data.
!
  implicit none

  integer ( kind = 4 ) ( kind = 4 ) m
  integer ( kind = 4 ) ( kind = 4 ) n

  integer ( kind = 4 ) ( kind = 4 ) ierror
  character ( len = * ) input_file_name
  integer ( kind = 4 ) ( kind = 4 ) input_status
  integer ( kind = 4 ) ( kind = 4 ) input_unit
  integer ( kind = 4 ) ( kind = 4 ) j
  character ( len = 255 ) line
  integer ( kind = 4 ) ( kind = 4 ) table(m,n)
  integer ( kind = 4 ) ( kind = 4 ) x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ITABLE_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_file_name ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      ierror = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ITABLE_DATA_READ - Fatal error!'
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
subroutine file_column_count ( input_file_name, column_num )

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
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) ( kind = 4 ) COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer ( kind = 4 ) ( kind = 4 ) column_num
  logical got_one
  character ( len = * ) input_file_name
  integer ( kind = 4 ) ( kind = 4 ) input_status
  integer ( kind = 4 ) ( kind = 4 ) input_unit
  character ( len = 256 ) line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = input_status )

  if ( input_status /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' &
      // trim ( input_file_name ) // '" on unit ', input_unit
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
subroutine file_row_count ( input_file_name, row_num )

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
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) ( kind = 4 ) ROW_NUM, the number of rows found.
!
  implicit none

  integer ( kind = 4 ) ( kind = 4 ) bad_num
  integer ( kind = 4 ) ( kind = 4 ) comment_num
  integer ( kind = 4 ) ( kind = 4 ) ierror
  character ( len = * ) input_file_name
  integer ( kind = 4 ) ( kind = 4 ) input_status
  integer ( kind = 4 ) ( kind = 4 ) input_unit
  character ( len = 100 ) line
  integer ( kind = 4 ) ( kind = 4 ) record_num
  integer ( kind = 4 ) ( kind = 4 ) row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_file_name ) // '" on unit ', input_unit
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
!    Output, integer ( kind = 4 ) ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) ( kind = 4 ) i
  integer ( kind = 4 ) ( kind = 4 ) ios
  integer ( kind = 4 ) ( kind = 4 ) iunit
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
subroutine i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_PRINT_SOME prints some of an integer matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 November 2003
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
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 10
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = 7 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7)') j
    end do

    write ( *, '(''  Col '',10a7)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(i7)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a7)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine i4vec_backtrack ( n, x, indx, k, nstack, stack, maxstack, ncan )

!*****************************************************************************80
!
!! I4VEC_BACKTRACK supervises a backtrack search for an integer vector.
!
!  Discussion:
!
!    The routine tries to construct an integer vector one index at a time,
!    using possible candidates as supplied by the user.
!
!    At any time, the partially constructed vector may be discovered to be
!    unsatisfactory, but the routine records information about where the
!    last arbitrary choice was made, so that the search can be
!    carried out efficiently, rather than starting out all over again.
!
!    First, call the routine with INDX = 0 so it can initialize itself.
!
!    Now, on each return from the routine, if INDX is:
!      1, you've just been handed a complete candidate vector;
!         Admire it, analyze it, do what you like.
!      2, please determine suitable candidates for position X(K).
!         Return the number of candidates in NCAN(K), adding each
!         candidate to the end of STACK, and increasing NSTACK.
!      3, you're done.  Stop calling the routine;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of positions to be filled in the vector.
!
!    Input/output, integer ( kind = 4 ) X(N), the partial or complete candidate vector.
!
!    Input/output, integer ( kind = 4 ) INDX, a communication flag.
!    On input,
!      0 to start a search.
!    On output:
!      1, a complete output vector has been determined and returned in X(1:N);
!      2, candidates are needed for position X(K);
!      3, no more possible vectors exist.
!
!    Output, integer ( kind = 4 ) K, if INDX=2, the current vector index being considered.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!
!    Input, integer ( kind = 4 ) STACK(MAXSTACK), a list of all current candidates for
!    all positions 1 through K.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the maximum length of the stack.
!
!    Input/output, integer ( kind = 4 ) NCAN(N), lists the current number of candidates for
!    positions 1 through K.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) indx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ncan(n)
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) stack(maxstack)
  integer ( kind = 4 ) x(n)
!
!  If this is the first call, request a candidate for position 1.
!
  if ( indx == 0 ) then
    k = 1
    nstack = 0
    indx = 2
    return
  end if
!
!  Examine the stack.
!
  do
!
!  If there are candidates for position K, take the first available
!  one off the stack, and increment K.
!
!  This may cause K to reach the desired value of N, in which case
!  we need to signal the user that a complete set of candidates
!  is being returned.
!
    if ( 0 < ncan(k) ) then

      x(k) = stack(nstack)
      nstack = nstack - 1

      ncan(k) = ncan(k) - 1

      if ( k /= n ) then
        k = k + 1
        indx = 2
      else
        indx = 1
      end if

      exit
!
!  If there are no candidates for position K, then decrement K.
!  If K is still positive, repeat the examination of the stack.
!
    else

      k = k - 1

      if ( k <= 0 ) then
        indx = 3
        exit
      end if

    end if

  end do

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an integer vector to the indicator vector.
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
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an integer vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) big
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  big = maxval ( abs ( a(1:n) ) )

  write ( *, '(a)' ) ' '
  if ( big < 1000 ) then
    do i = 1, n
      write ( *, '(2x,i6,2x,i4)' ) i, a(i)
    end do
  else if ( big < 1000000 ) then
    do i = 1, n
      write ( *, '(2x,i6,2x,i7)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(2x,i6,2x,i12)' ) i, a(i)
    end do
  end if

  return
end
subroutine itable_header_read ( input_file_name, m, n )

!*****************************************************************************80
!
!! ITABLE_HEADER_READ reads the header from an integer table file.
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
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) ( kind = 4 ) M, spatial dimension.
!
!    Output, integer ( kind = 4 ) ( kind = 4 ) N, the number of points. 
!
  implicit none

  character ( len = * ) input_file_name
  integer ( kind = 4 ) ( kind = 4 ) m
  integer ( kind = 4 ) ( kind = 4 ) n

  call file_column_count ( input_file_name, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ITABLE_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_file_name ) // '".'
    stop
  end if

  call file_row_count ( input_file_name, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ITABLE_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_file_name ) // '".'
    stop
  end if

  return
end
subroutine mg_code_back ( adj, nnode, code, order )

!*****************************************************************************80
!
!! MG_CODE_BACK computes a multigraph code via backtracking.
!
!  Discussion:
!
!    The code is the "largest" order code (in the lexicographic
!    sense) over all possible node orderings.  The backtracking method
!    organizes the search of all possible node orderings so that if
!    a partial node ordering is sure to generate an order code
!    less than the best so far, then all the orderings that begin with
!    this partial ordering are immediately dropped from consideration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) bestorder(nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) index
  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nopen
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) nswap
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ), allocatable, dimension ( : ) :: stack

  if ( nnode <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)') 'MG_CODE_BACK - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    stop
  end if

  maxstack = 10 * nnode
  allocate ( stack(1:maxstack) )

  nstack = 0
  stack(1) = 0

  ncomp = 0
  nswap = 0
!
!  Start with the identity ordering.
!
  call i4vec_indicator ( nnode, order )
!
!  Compute the corresponding code.
!
  call mg_order_code ( adj, nnode, nnode, code, order )
!
!  Take this ordering and code as the best so far.
!
  bestorder(1:nnode) = order(1:nnode)
  bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)
!
!  Now consider all possible orderings, and their codes.
!
  index = 0

  do

    call i4vec_backtrack ( nnode, order, index, nopen, nstack, stack, &
      maxstack, ncan )
!
!  If the backtrack routine has returned a complete candidate ordering, then
!  compute the resulting code, and see it it is better
!  then our current best.  Then go back for the next backtrack search.
!
    if ( index == 1 ) then

      call mg_order_code ( adj, nnode, nnode, code, order )

      call mg_code_compare ( bestcode, code, nnode, nnode, result )

      ncomp = ncomp + 1

      if ( result == -1 ) then

        nswap = nswap + 1

        bestorder(1:nnode) = order(1:nnode)
        bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)

      end if
!
!  If the backtrack routine has returned a partial reordering,
!  supply candidates for the next item in the ordering.
!
    else if ( index == 2 ) then

      call mg_code_cand ( adj, bestcode, code, nnode, ncomp, nopen, &
        order, stack, maxstack, nstack, ncan )

    else

      exit

    end if

  end do
!
!  Once we have examined all possibilites, we are done.
!
!  Set the output ordering to the best ordering, and the output
!  code to the corresponding best code.
!
  order(1:nnode) = bestorder(1:nnode)
  code(1:nnode,1:nnode) = bestcode(1:nnode,1:nnode)

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MG_CODE_BACK:'
    write ( *, '(a,i8)' ) '  Comparisons: ', ncomp
    write ( *, '(a,i8)' ) '  Swaps:       ', nswap
  end if

  deallocate ( stack )

  return
end
subroutine mg_code_cand ( adj, bestcode, code, nnode, ncomp, nopen, &
  order, stack, maxstack, nstack, ncan )

!*****************************************************************************80
!
!! MG_CODE_CAND finds candidates for a maximal multigraph code ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) BESTCODE(NNODE,NNODE), the best code so far.
!
!    Workspace, integer CODE(NNODE,NNODE).
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input/output, integer ( kind = 4 ) NCOMP, the number of code comparisons.
!    This routine updates NCOMP by 1 each time the routine is called.
!
!    Input, integer ( kind = 4 ) NOPEN, identifies the first open position in ORDER.
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), contains in entries 1 through NOPEN-1
!    the elements of the current partial list.
!
!    Input/output, integer ( kind = 4 ) STACK(MAXSTACK), used to store the new candidates.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the maximum size of the STACK array.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!    On output, NSTACK has been increased by the number of
!    candidates found.
!
!    Input/output, integer ( kind = 4 ) NCAN(NNODE), the number of candidates
!    for each position.
!
  implicit none

  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) max_adj(nnode)
  integer ( kind = 4 ) max_adj_ni
  integer ( kind = 4 ) max_max_adj
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nfree
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) nopen
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) stack(maxstack)

  if ( nopen < 1 .or. nnode < nopen ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MG_CODE_CAND - Fatal error 1!'
    write ( *, '(a)' ) '  1 <= NOPEN <= NNODE should be true, but'
    write ( *, '(a,i8)' ) '  NOPEN = ', nopen
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    stop
  end if

  ncan(nopen) = 0
!
!  If we have fixed at least one entry of the list, then compare
!  the code of the current incomplete ordering to the
!  code induced by the corresponding part of the current
!  best ordering.
!
!  If the current incomplete ordering is actually LESS than the
!  current best, then bail out with zero candidates.
!
  if ( 1 < nopen ) then

    call mg_order_code ( adj, nnode, nopen-1, code, order )

    call mg_code_compare ( bestcode, code, nnode, nopen-1, result )

    ncomp = ncomp + 1

    if ( result == +1 ) then
      ncan(nopen) = 0
      return
    end if

  end if
!
!  Get a list of those nodes which have not been used yet.
!
  nfree = nnode + 1 - nopen
  call perm_free ( order, nopen-1, ifree, nfree )
!
!  To find candidates, we consider all the ordered nodes.
!  We find the lowest ordered node that has unordered neighbors.
!  We take the unordered neighbors with maximal adjacency.
!
  ncan(nopen) = 0
!
!  For each ordered node NI...
!
  do i = 1, nopen-1

    ni = order(i)
!
!  ...note the maximum adjacency from NI to any unordered node NJ...
!
    max_adj_ni = 0
    do j = 1, nfree
      nj = ifree(j)
      max_adj_ni = max ( max_adj_ni, adj(ni,nj) )
    end do
!
!   ...and take as candidates all unordered nodes NJ with maximal
!   adjacency from NI.
!
    if ( 0 < max_adj_ni ) then

      do j = 1, nfree

        nj = ifree(j)

        if ( adj(ni,nj) == max_adj_ni ) then

          ncan(nopen) = ncan(nopen) + 1
          nstack = nstack + 1

          if ( maxstack < nstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'MG_CODE_CAND - Fatal error 2!'
            write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
            write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
            write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
            stop
          end if

          stack(nstack) = nj

        end if

      end do

      return

    end if

  end do
!
!  If we get here, no unordered nodes are connected in any way to
!  ordered nodes.  This can happen in two ways:
!
!  * NOPEN = 1;
!  * The graph is disconnected;
!
!  For each free node, compute the maximum adjacency to any other free node.
!  Take the maximum of this value over all free nodes.
!  Candidates are free nodes with this maximum value.
!
  max_max_adj = 0

  do i = 1, nfree

    ni = ifree(i)

    max_adj(i) = 0
    do j = 1, nfree
      nj = ifree(j)
      if ( ni /= nj ) then
        max_adj(i) = max ( max_adj(i), adj(ni,nj) )
      end if
    end do

    max_max_adj = max ( max_max_adj, max_adj(i) )

  end do
!
!  Now go back and compute the maximum adjacency for each free node.
!  Nodes which achieve the maximum are added to the candidate list.
!
  ncan(nopen) = 0

  do i = 1, nfree

    if ( max_adj(i) == max_max_adj ) then

      ncan(nopen) = ncan(nopen) + 1
      nstack = nstack + 1

      if ( maxstack < nstack ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MG_CODE_CAND - Fatal error 2!'
        write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
        write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
        write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
        stop
      end if

      stack(nstack) = ifree(i)

    end if

  end do

  return
end
subroutine mg_code_compare ( code1, code2, nnode, npart, result )

!*****************************************************************************80
!
!! MG_CODE_COMPARE compares two partial multigraph codes.
!
!  Discussion:
!
!    CODE1 = CODE2 if every digit of both codes is equal.  
!
!    Otherwise, traverse the entries in a funny diagonal way, suggested
!    by this diagram for the first 10 entries:
!
!       *  1  2  4  7
!       *  *  3  5  8
!       *  *  *  6  9
!       *  *  *  * 10
!       *  *  *  *  *
!
!    As we do that, we examine the corresponding digits of the two codes.  
!    For the first entry, (I,J), where the two codes differ, we say:
!
!      if ( CODE1(I,J) < CODE2(I,J) ) then we say
!        CODE1 < CODE2
!      else
!        CODE2 < CODE1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) CODE1(NNODE,NNODE), CODE2(NNODE,NNODE), codes to compare.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPART, specifies the portion of the codes
!    to compare.  NPART should be between 1 and NNODE.  
!
!    If NPART = NNODE, then the full codes are compared.
!
!    If NPART < NNODE, then only entries corresponding to I and J
!    both less than or equal to NPART will be compared.
!
!    Output, integer ( kind = 4 ) RESULT:
!    -1, CODE1 < CODE2;
!     0, CODE1 = CODE2;
!    +1, CODE2 < CODE1.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code1(nnode,nnode)
  integer ( kind = 4 ) code2(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) result

  do j = 2, npart

    do i = 1, j-1

      if ( code1(i,j) < code2(i,j) ) then

        result = - 1
        return

      else if ( code2(i,j) < code1(i,j) ) then

        result = + 1
        return

      end if

    end do

  end do

  result = 0

  return
end
subroutine mg_code_print ( nnode, code, title )

!*****************************************************************************80
!
!! MG_CODE_PRINT prints out a multigraph code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  character ( len = 80 ) string
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode
    string(i:i) = '.'
  end do

  do i = 1, nnode

    write ( *, '(2x,a,80i1)' ) string(1:i), code(i,i+1:nnode)

  end do

  return
end
subroutine mg_order_code ( adj, nnode, npart, code, order )

!*****************************************************************************80
!
!! MG_ORDER_CODE returns the multigraph code for a specific node ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!    
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPART, the number of nodes to consider.  
!    NPART should be between 1 and NNODE.  
!
!    If NPART is NNODE, then the full code is returned.  
!
!    If NPART is less than NNODE, then the code is computed as 
!    though only NPART nodes existed, namely, those specified in the 
!    first NPART entries of order.  This option is provided so that
!    the routine can compute the portion of a code determined
!    by an incomplete ordering of the nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code for this ordering.  
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.  ORDER(1)
!    is the first node, and so on.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) order(nnode)

  do i = 1, nnode

    if ( i <= npart ) then

      ni = order(i)

      if ( order(i) < 1 .or. nnode < order(i) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MG_ORDER_CODE - Fatal error!'
        write ( *, '(a)' ) '  ORDER is not a proper permutation.'
        stop
      end if

    else
      ni = 0
    end if

    do j = i + 1, nnode

      if ( j <= npart ) then

        nj = order(j)

        if ( order(j) < 1 .or. nnode < order(j) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'MG_ORDER_CODE - Fatal error!'
          write ( *, '(a)' ) '  ORDER is not a proper permutation.'
          stop
        end if

      else
        nj = 0
      end if

      if ( ni == 0 .or. nj == 0 ) then

        code(i,j) = 0
        code(j,i) = 0

      else if ( ( ni < nj .and. adj(ni,nj) /= 0 ) .or. &
                ( nj < ni .and. adj(nj,ni) /= 0 ) ) then

        code(i,j) = adj(ni,nj)
        code(j,i) = adj(nj,ni)

      else

        code(i,j) = 0
        code(j,i) = 0

      end if

    end do
  end do

  return
end
subroutine node_order_print ( nnode, order, title )

!*****************************************************************************80
!
!! NODE_ORDER_PRINT prints out a node ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), the node ordering.  ORDER(1) is the label
!    of the node which is to be taken as the first node, and so on.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) order(nnode)
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if
 
  inc = 15

  do ilo = 1, nnode, inc

    ihi = min ( ilo + inc - 1, nnode )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a,2x,15i4)' ) 'Order:', ( i, i = ilo, ihi )
    write ( *, '(2x,a,2x,15i4)' ) 'Label:', order(ilo:ihi)

  end do

  return
end
subroutine perm_free ( ipart, npart, ifree, nfree )

!*****************************************************************************80
!
!! PERM_FREE reports the number of unused items in a partial permutation.
!
!  Discussion:
!
!    It is assumed that the N objects being permuted are the integers
!    from 1 to N, and that IPART contains a "partial" permutation, that
!    is, the NPART entries of IPART represent the beginning of a
!    permutation of all N items.
!
!    The routine returns in IFREE the items that have not been used yet.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IPART(NPART), the partial permutation, which should
!    contain, at most once, some of the integers between 1 and
!    NPART+NFREE.
!
!    Input, integer ( kind = 4 ) NPART, the number of entries in IPART.  NPART may be 0.
!
!    Output, integer ( kind = 4 ) IFREE(NFREE), the integers between 1 and NPART+NFREE
!    that were not used in IPART.
!
!    Input, integer ( kind = 4 ) NFREE, the number of integers that have not been
!    used in IPART.  This is simply N - NPART.  NFREE may be zero.
!
  implicit none

  integer ( kind = 4 ) nfree
  integer ( kind = 4 ) npart

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree(nfree)
  integer ( kind = 4 ) ipart(npart)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) match
  integer ( kind = 4 ) n

  n = npart + nfree

  if ( npart < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
    write ( *, '(a)' ) '  NPART < 0.'
    write ( *, '(a,i8)' ) '  NPART = ', npart
    stop

  else if ( npart == 0 ) then

    call i4vec_indicator ( n, ifree )

  else if ( nfree < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
    write ( *, '(a)' ) '  NFREE < 0.'
    write ( *, '(a,i8)' ) '  NFREE = ', nfree
    stop

  else if ( nfree == 0 ) then

    return

  else

    k = 0

    do i = 1, n

      match = 0

      do j = 1, npart
        if ( ipart(j) == i ) then
          match = j
          exit
        end if
      end do

      if ( match == 0 ) then

        k = k + 1

        if ( nfree < k ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)'    ) 'PERM_FREE - Fatal error!'
          write ( *, '(a)'    ) '  The partial permutation is illegal.'
          write ( *, '(a)'    ) '  It should contain, at most once, some of'
          write ( *, '(a,i8)' ) '  the integers between 1 and N = ', n
          write ( *, '(a)'    ) '  The number of integers that have not'
          write ( *, '(a,i8)' ) '  been used is at least K = ', k
          write ( *, '(a,i8)' ) '  This value should be exactly NFREE = ', &
            nfree
          call i4vec_print ( npart, ipart, '  The partial permutation:' )
          stop
        end if

        ifree(k) = i

      end if

    end do

  end if

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
!    Output, integer ( kind = 4 ) ( kind = 4 ) IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) ( kind = 4 ) LENGTH, the number of characters of S 
!    used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) ( kind = 4 ) i
  integer ( kind = 4 ) ( kind = 4 ) ierror
  integer ( kind = 4 ) ( kind = 4 ) isgn
  integer ( kind = 4 ) ( kind = 4 ) istate
  integer ( kind = 4 ) ( kind = 4 ) ival
  integer ( kind = 4 ) ( kind = 4 ) length
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
!    Input, integer ( kind = 4 ) ( kind = 4 ) N, the number of values expected.
!
!    Output, integer ( kind = 4 ) ( kind = 4 ) IVEC(N), the values read from the string.
!
!    Output, integer ( kind = 4 ) ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer ( kind = 4 ) ( kind = 4 ) n

  integer ( kind = 4 ) ( kind = 4 ) i
  integer ( kind = 4 ) ( kind = 4 ) ierror
  integer ( kind = 4 ) ( kind = 4 ) ilo
  integer ( kind = 4 ) ( kind = 4 ) ivec(n)
  integer ( kind = 4 ) ( kind = 4 ) length
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
!    Output, integer ( kind = 4 ) ( kind = 4 ) NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer ( kind = 4 ) ( kind = 4 ) i
  integer ( kind = 4 ) ( kind = 4 ) lens
  integer ( kind = 4 ) ( kind = 4 ) nword
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
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2001
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
  integer ( kind = 4 ) ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) ( kind = 4 ) h
  integer ( kind = 4 ) ( kind = 4 ) m
  integer ( kind = 4 ) ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) ( kind = 4 ) n
  integer ( kind = 4 ) ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) ( kind = 4 ) values(8)
  integer ( kind = 4 ) ( kind = 4 ) y
  character ( len = 5 ) zone

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
