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
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf
!    FORTRAN90 version by John Burkardt
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
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements
!    I and J. (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

  return
end
subroutine st_data_print ( nrow, ncol, nnzero, row, col, a, title )

!*****************************************************************************80
!
!! ST_DATA_PRINT prints the data in an ST file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the assumed number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the assumed number of columns in the matrix.
!
!    Input, integer ( kind = 4 ) NNZERO, the assumed number of nonzeros in the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NNZERO), COL(NNZERO), the 0-based row and column
!    index of a nonzero matrix entry.
!
!    Input, real ( kind = 8 ) A(NNZERO), the value of a nonzero matrix entry.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) nnzero

  real ( kind = 8 ) a(nnzero)
  integer ( kind = 4 ) col(nnzero)
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) row(nnzero)
  character ( len = * ) title

  call st_data_print_some ( 0, nrow-1, 0, ncol-1, nnzero, row, col, a, title )

  return
end
subroutine st_data_print_some ( row1, row2, col1, col2, nnzero, row, col, a, &
  title )

!*****************************************************************************80
!
!! ST_DATA_PRINT_SOME prints some of the data in an ST file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ROW1, ROW2, the first and last rows to print.
!
!    Input, integer ( kind = 4 ) COL1, COL2, the first and last columns to print.
!
!    Input, integer ( kind = 4 ) NNZERO, the assumed number of nonzeros in the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NNZERO), COL(NNZERO), the 0-based row and column
!    index of a nonzero matrix entry.
!
!    Input, real ( kind = 8 ) A(NNZERO), the value of a nonzero matrix entry.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) nnzero

  real ( kind = 8 ) a(nnzero)
  integer ( kind = 4 ) col(nnzero)
  integer ( kind = 4 ) col1
  integer ( kind = 4 ) col2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) row(nnzero)
  integer ( kind = 4 ) row1
  integer ( kind = 4 ) row2
  character ( len = *  ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim(title)
  write ( *, '(a)' ) ' '

  do k = 1, nnzero
    if ( row1 <= row(k) .and. row(k) <= row2 .and. &
         col1 <= col(k) .and. col(k) <= col2 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8,2x,g16.8)' ) k, row(k), col(k), a(k)
    end if
  end do

  return
end
subroutine st_data_read ( input_filename, nrow, ncol, nnzero, row, col, a )

!*****************************************************************************80
!
!! ST_DATA_READ reads the data of an ST file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the ST file.
!
!    Input, integer ( kind = 4 ) NROW, the assumed number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the assumed number of columns in the matrix.
!
!    Input, integer ( kind = 4 ) NNZERO, the assumed number of nonzeros in the matrix.
!
!    Output, integer ( kind = 4 ) ROW(NNZERO), COL(NNZERO), the 0-based row and column
!    index of a nonzero matrix entry.
!
!    Output, real ( kind = 8 ) A(NNZERO), the value of a nonzero matrix entry.
!
  implicit none

  integer ( kind = 4 ) nnzero

  real ( kind = 8 ) a(nnzero)
  real ( kind = 8 ) aij
  integer ( kind = 4 ) col(nnzero)
  integer ( kind = 4 ) i
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) row(nnzero)

  call get_unit ( input_unit )

  open ( file = input_filename, unit = input_unit, status = 'old', iostat = ios )

  do k = 1, nnzero

    read ( input_unit, *, iostat = ios ) i, j, aij

    if ( ios /= 0 ) then
      exit
    end if

    row(k) = i
    col(k) = j
    a(k) = aij

  end do

  close ( unit = input_unit )

  return
end
subroutine st_data_write ( output_filename, nrow, ncol, nnzero, row, col, a )

!*****************************************************************************80
!
!! ST_DATA_WRITE writes the data of an ST file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the ST file.
!
!    Input, integer ( kind = 4 ) NROW, the assumed number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the assumed number of columns in the matrix.
!
!    Input, integer ( kind = 4 ) NNZERO, the assumed number of nonzeros in the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NNZERO), COL(NNZERO), the 0-based row and column
!    index of a nonzero matrix entry.
!
!    Input, real ( kind = 8 ) A(NNZERO), the value of a nonzero matrix entry.
!
  implicit none

  integer ( kind = 4 ) nnzero

  real ( kind = 8 ) a(nnzero)
  real ( kind = 8 ) aij
  integer ( kind = 4 ) col(nnzero)
  integer ( kind = 4 ) i
  character ( len = * )  output_filename
  integer ( kind = 4 ) output_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) row(nnzero)

  call get_unit ( output_unit )

  open ( file = output_filename, unit = output_unit, status = 'replace', &
    iostat = ios )

  do k = 1, nnzero

    write ( output_unit, '(2x,i8,2x,i8,2x,g16.8)', iostat = ios ) &
      row(k), col(k), a(k)

  end do

  close ( unit = output_unit )

  return
end
subroutine st_header_print ( nrow, ncol, nnzero )

!*****************************************************************************80
!
!! ST_HEADER_PRINT prints the header of an ST file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the assumed number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the assumed number of columns in the matrix.
!
!    Input, integer ( kind = 4 ) NNZERO, the assumed number of nonzeros in the matrix.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrow

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ST header information:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of rows =     ', nrow
  write ( *, '(a,i8)' ) '  Number of columns =  ', ncol
  write ( *, '(a,i8)' ) '  Number of nonzeros = ', nnzero

  return
end
subroutine st_header_read ( input_filename, nrow, ncol, nnzero )

!*****************************************************************************80
!
!! ST_HEADER_READ reads the header of an ST file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the ST file.
!
!    Output, integer ( kind = 4 ) NROW, the assumed number of rows in the matrix.
!
!    Output, integer ( kind = 4 ) NCOL, the assumed number of columns in the matrix.
!
!    Output, integer ( kind = 4 ) NNZERO, the assumed number of nonzeros
!    in the matrix.
!
  implicit none

  real ( kind = 8 ) aij
  integer ( kind = 4 ) i
  character ( len = * )  input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrow

  call get_unit ( input_unit )

  open ( file = input_filename, unit = input_unit, status = 'old', iostat = ios )

  nnzero = 0
  nrow = 0
  ncol = 0

  do

    read ( input_unit, *, iostat = ios ) i, j, aij

    if ( ios /= 0 ) then
      exit
    end if

    nnzero = nnzero + 1
    nrow = max ( nrow, i + 1 )
    ncol = max ( ncol, j + 1 )

  end do

  close ( unit = input_unit )

  return
end
subroutine st_header_write ( output_filename, nrow, ncol, nnzero )

!*****************************************************************************80
!
!! ST_HEADER_WRITE writes the header of an ST file.
!
!  Discussion:
!
!    There is no header for an ST file.  But for completeness, we supply
!    this dummy routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the ST file.
!
!    Input, integer ( kind = 4 ) NROW, the assumed number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the assumed number of columns in the matrix.
!
!    Input, integer ( kind = 4 ) NNZERO, the assumed number of nonzeros in the matrix.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnzero
  integer ( kind = 4 ) nrow
  character ( len = * ) output_filename

  return
end
subroutine st_rebase ( base1, base2, nnzero, indx )

!*****************************************************************************80
!
!! ST_REBASE changes the base of an index array.
!
!  Discussion:
!
!    Both the ROW and COL arrays are ordinarily 0-based in the ST format.
!    FORTRAN and MATLAB expect 1-based indices.
!
!    To convert ROW and COL from 0-based to 1-based form, call this routine
!    with BASE1=0, BASE2=1.
!
!    If ROW and COL from FORTRAN or MATLAB are to be converted to ST format,
!    call this routine with BASE1=1 and BASE2=0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) BASE1, BASE2, the old and new bases.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeros in the matrix.
!
!    Input/output, integer ( kind = 4 ) INDX(NNZERO), the index vector
!    to be rebased.
!
  implicit none

  integer ( kind = 4 ) nnzero

  integer ( kind = 4 ) base1
  integer ( kind = 4 ) base2
  integer ( kind = 4 ) indx(nnzero)

  indx(1:nnzero) = indx(1:nnzero) - base1 + base2

  return
end
subroutine st_sort_a ( nrow, ncol, nnzero, row, col, a )

!*****************************************************************************80
!
!! ST_SORT_A sorts the entries of an ST matrix by column.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns in the matrix.
!
!    Input, integer ( kind = 4 ) NNZERO, the number of nonzeros in the matrix.
!
!    Input/output, integer ( kind = 4 ) ROW(NNZERO), COL(NNZERO), the
!    0-based row and column index of a nonzero matrix entry.
!
!    Input/output, real ( kind = 8 ) A(NNZERO), the nonzero matrix entries.
!
  implicit none

  integer ( kind = 4 ) nnzero

  real ( kind = 8 ) a(nnzero)
  real ( kind = 8 ) aij
  integer ( kind = 4 ) cij
  integer ( kind = 4 ) col(nnzero)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) rij
  integer ( kind = 4 ) row(nnzero)
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( nnzero, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      rij    = row(i)
      row(i) = row(j)
      row(j) = rij

      cij    = col(i)
      col(i) = col(j)
      col(j) = cij

      aij  = a(i)
      a(i) = a(j)
      a(j) = aij
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      if ( col(i) == col(j) ) then

        if ( row(i) < row(j) ) then
          isgn = - 1
        else if ( row(i) == row(j) ) then
          isgn = 0
        else if ( row(j) < row(i) ) then
          isgn = + 1
        end if

      else if ( col(i) < col(j) ) then

        isgn = - 1

      else if ( col(j) < col(i) ) then

        isgn = + 1

      end if

    else if ( indx == 0 ) then

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
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2005
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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
