subroutine ch_cap ( ch )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
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

  itemp = ichar ( ch )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    ch = char ( itemp - 32 )
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
subroutine mm_comment_print ( comment )

!*****************************************************************************80
!
!! MM_COMMENT_PRINT prints a comment from a Matrix Market file.
!
!  Discussion:
!
!    Comment lines begin with a '%' character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMENT, the comment.
!
  implicit none

  character ( len = * ) comment

  write ( *, '(a)' ) trim ( comment )

  return
end
subroutine mm_comment_read ( input_unit, comment )

!*****************************************************************************80
!
!! MM_COMMENT_READ reads a comment from a Matrix Market file.
!
!  Discussion:
!
!    This routine simply reads one line from the file.  Comment
!    lines begin with a '%' character, and must occur between the
!    header line and the size line.  It is up to the user to
!    examine the line returned by this routine, and, if it does not
!    begin with a comment character, then the user may do any of:
!    *) terminate the processing, because you simply wanted to see comments;
!    *) backspace the file, so the line can be processed by MM_SIZE_READ_FILE;
!    *) simply pass the line returned by this routine directly
!       to MM_SIZE_READ_STRING.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT_UNIT, the input unit identifier.
!
!    Output, character ( len = * ) COMMENT, the next line of the file.
!
  implicit none

  character ( len = * ) comment
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios

  read ( input_unit, '(a)', iostat = ios ) comment

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MM_COMMENT_READ'
    write ( *, '(a)' ) '  Unexpected end of file while reading comments.'
    stop
  end if

  return
end
subroutine mm_comment_write ( output_unit, comment )

!*****************************************************************************80
!
!! MM_COMMENT_WRITE writes a comment to a Matrix Market file.
!
!  Discussion:
!
!    Comments may be written AFTER the header line and BEFORE the size line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output unit identifier.
!
!    Input, character ( len = * ) COMMENT, the comment to be written.
!    This routine will prepend a "%  " comment marker to the string.
!
  implicit none

  character ( len = * ) comment
  integer ( kind = 4 ) output_unit

  if ( len_trim ( comment ) == 0 ) then
    write ( output_unit, '(a)' ) '%'
  else
    write ( output_unit, '(a)' ) '%  ' // trim ( comment )
  end if

  return
end
subroutine mm_file_read ( input_unit, id, type, rep, field, symm, nrow, &
  ncol, nnz, nnzmax, indx, jndx, ival, rval, dval, cval )

!*****************************************************************************80
!
!! MM_FILE_READ reads data from a Matrix Market file.
!
!  Discussion:
!
!    The data may be either sparse coordinate format, or dense array format.
!
!    The unit input_unit must be open, and the file will be rewound on return.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2003
!
!  Author:
!
!    Original FORTRAN77 version by Karin A. Remington, NIST ACMD
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT_UNIT, the input unit identifier.
!
!    Output, character ( len = 14 ) ID, the Matrix Market identifier.
!    This value must be '%%MatrixMarket'.
!
!    Output, character ( len = 6 ) TYPE, the Matrix Market type.
!    This value must be 'matrix'.
!
!    Output, character ( len = 10 ) REP, the Matrix Market 'representation'
!    indicator.  Possible values include:
!    'coordinate'   (for sparse data)
!    'array'        (for dense data)
!    'elemental'    (to be added)
!
!    Output, character ( len = 7 ) FIELD, the Matrix Market 'field'.
!    Possible values include:
!    'real'
!    'double'
!    'complex'
!    'integer'
!    'pattern'
!
!    Output, character ( len = 19 ) SYMM, the Matrix Market symmetry.
!    Possible values include:
!    'symmetric'
!    'hermitian'
!    'skew-symmetric'
!    'general'
!
!    Output, integer ( kind = 4 ) NROW, the number of rows in the matrix.
!
!    Output, integer ( kind = 4 ) NCOL, the number of columns in the matrix.
!
!    Output, integer ( kind = 4 ) NNZ, the number of nonzero entries required
!    to store the matrix, if REP = 'coordinate'.
!
!    Input, integer ( kind = 4 ) NNZMAX, the maximum dimension of the appropriate
!    data array.
!
!    Output, integer ( kind = 4 ) INDX(NNZ), the row indices for coordinate format.
!    Not used if REP is 'array'.
!
!    Output, integer ( kind = 4 ) JNDX(NNZ), the column indices for coordinate
!    format.  Not used if REP is 'array'.
!
!    Output, integer ( kind = 4 ) IVAL(NNZ), matrix values, if FIELD is 'integer'.
!
!    Output, real ( kind = 4 ) RVAL(NNZ), matrix values, if FIELD is 'real'.
!
!    Output, real ( kind = 8 ) DVAL(NNZ), matrix values, if FIELD is 'double'.
!
!    Output, complex ( kind = 4 ) CVAL(NNZ), matrix values, if FIELD is 'complex'.
!
  implicit none

  complex ( kind = 4 ) cval(*)
  real ( kind = 8 ) dval(*)
  character ( len = 7 ) field
  character ( len = 14 ) id
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ival(*)
  integer ( kind = 4 ) jndx(*)
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnz
  integer ( kind = 4 ) nnzmax
  integer ( kind = 4 ) nrow
  character ( len = 10 ) rep
  real ( kind = 4 ) rval(*)
  character ( len = 19 ) symm
  character ( len = 1024 ) tmp1
  character ( len = 6 ) type
!
!  Read and check the header line.
!
  call mm_header_read ( input_unit, id, type, rep, field, symm )

  call mm_header_check ( id, type, rep, field, symm )
!
!  Read through the comment lines:
!
  do

    call mm_comment_read ( input_unit, tmp1 )

    if ( tmp1(1:1) /= '%' ) then
      exit
    end if

  end do
!
!  Current line is not a comment.
!
!  We can either process the string directly by MM_SIZE_READ_STRING,
!  or backspace the file and have MM_SIZE_READ_FILE handle it.
!
  if ( .true. ) then

    call mm_size_read_string ( tmp1, rep, symm, nrow, ncol, nnz )

  else

    backspace ( input_unit )

    call mm_size_read ( input_unit, rep, symm, nrow, ncol, nnz )

  end if
!
!  Make sure there is enough storage space.
!
  if ( nnzmax < nnz ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MM_FILE_READ - Fatal error!'
    write ( *, '(a)' ) '  Insufficient memory to read matrix.'
    write ( *, '(a)' ) '  Number of entries to be read = ', nnz
    write ( *, '(a)' ) '  Number of entries supplied = ',nnzmax
    stop
  end if
!
!  Read the data values.
!
  call mm_values_read ( input_unit, rep, field, nnz, indx, jndx, &
    ival, rval, dval, cval )

  return
end
subroutine mm_file_write ( output_unit, id, type, rep, field, symm, nrow, &
  ncol, nnz, indx, jndx, ival, rval, dval, cval )

!*****************************************************************************80
!
!! MM_FILE_WRITE writes data to a Matrix Market file.
!
!  Discussion:
!
!    The data may be either sparse coordinate format, or dense array format.
!
!    The unit OUTPUT_UNIT must be open.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2003
!
!  Author:
!
!    Original FORTRAN77 version by Karin A. Remington, NIST ACMD
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output unit identifier.
!
!    Input, character ( len = 14 ) ID, the Matrix Market identifier.
!    This value must be '%%MatrixMarket'.
!
!    Input, character ( len = 6 ) TYPE, the Matrix Market type.
!    This value must be 'matrix'.
!
!    Input, character ( len = 10 ) REP, the Matrix Market 'representation'
!    indicator.  Possible values include:
!    'coordinate'   (for sparse data)
!    'array'        (for dense data)
!    'elemental'    (to be added)
!
!    Input, character ( len = 7 ) FIELD, the Matrix Market 'field'.
!    Possible values include:
!    'real'
!    'double'
!    'complex'
!    'integer'
!    'pattern' (for REP = 'coordinate' only)
!
!    Input, character ( len = 19 ) SYMM, the Matrix Market symmetry.
!    Possible values include:
!    'symmetric'
!    'hermitian'
!    'skew-symmetric'
!    'general'
!
!    Input, integer ( kind = 4 ) NROW, the number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns in the matrix.
!
!    Input, integer ( kind = 4 ) NNZ, the number of nonzero entries required
!    to store the matrix, if REP = 'coordinate'.
!
!    Input, integer ( kind = 4 ) INDX(NNZ), the row indices for coordinate format.
!    Not used if REP is 'array'.
!
!    Input, integer ( kind = 4 ) JNDX(NNZ), the column indices for coordinate
!    format.  Not used if REP is 'array'.
!
!    Input, integer ( kind = 4 ) IVAL(NNZ), matrix values, if FIELD is 'integer'.
!
!    Input, real ( kind = 4 ) RVAL(NNZ), matrix values, if FIELD is 'real'.
!
!    Input, real ( kind = 8 ) DVAL(NNZ), matrix values, if FIELD is 'double'.
!
!    Input, complex ( kind = 4 ) CVAL(NNZ), matrix values, if FIELD is 'complex'.
!
  implicit none

  complex ( kind = 4 ) cval(*)
  real ( kind = 8 ) dval(*)
  character ( len = * ) field
  character ( len = 14 ) id
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) ival(*)
  integer ( kind = 4 ) jndx(*)
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnz
  integer ( kind = 4 ) nnz2
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) output_unit
  character ( len = * ) rep
  real ( kind = 4 ) rval(*)
  logical s_eqi
  character ( len = * ) symm
  character ( len = 6 ) type
!
!  Test the header values.
!
  call mm_header_check ( id, type, rep, field, symm )
!
!  Write the header line.
!
  call mm_header_write ( output_unit, id, type, rep, field, symm )
!
!  Write some comment lines.
!
  call mm_comment_write ( output_unit, ' ' )
  call mm_comment_write ( output_unit, &
    'This file created by MM_FILE_WRITE of MM_IO.F90.' )
  call mm_comment_write ( output_unit, ' ' )
!
!  Write the size line.
!
  call mm_size_write ( output_unit, rep, nrow, ncol, nnz )
!
!  Determine NNZ where necessary.
!
  call mm_nnz_set ( rep, symm, nrow, ncol, nnz2 )
!
!  Write the data.
!
  if ( s_eqi ( rep, 'coordinate' ) ) then
    call mm_values_write ( output_unit, rep, field, nnz, indx, jndx, &
      ival, rval, dval, cval )
  else
    call mm_values_write ( output_unit, rep, field, nnz2, indx, jndx, &
      ival, rval, dval, cval )
  end if

  return
end
function mm_header_check ( id, type, rep, field, symm )

!*****************************************************************************80
!
!! MM_HEADER_CHECK checks the header strings for a Matrix Market file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 14 ) ID, the Matrix Market identifier.
!    This value must be '%%MatrixMarket'.
!
!    Input, character ( len = 6 ) TYPE, the Matrix Market type.
!    This value must be 'matrix'.
!
!    Input, character ( len = 10 ) REP, the Matrix Market 'representation'
!    indicator.  Possible values include:
!    'coordinate'   (for sparse data)
!    'array'        (for dense data)
!    'elemental'    (to be added)
!
!    Input, character ( len = 7 ) FIELD, the Matrix Market 'field'.
!    Possible values include:
!    'real'
!    'double'
!    'complex'
!    'integer'
!    'pattern' (for REP = 'coordinate' only)
!
!    Input, character ( len = 19 ) SYMM, the Matrix Market symmetry.
!    Possible values include:
!    'symmetric'
!    'hermitian'
!    'skew-symmetric'
!    'general'
!
!    Output, logical MM_HEADER_CHECK, is TRUE if the header is OK.
!
  implicit none

  character ( len = * ) field
  character ( len = 14 ) id
  logical mm_header_check
  character ( len = * ) rep
  logical s_eqi
  logical s_neqi
  character ( len = * ) symm
  character ( len = 6 ) type
!
!  Test the input qualifiers.
!
  if ( s_neqi ( id, '%%MatrixMarket' ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MM_HEADER_CHECK - Fatal error!'
    write ( *, '(a)' ) '  The value of ID was illegal:'
    write ( *, '(a)' ) '    "' // trim ( id ) // '".'
    write ( *, '(a)' ) '  Legal values are:'
    write ( *, '(a)' ) '    "%%MatrixMarket"'
    mm_header_check = .false.
    return
  end if

  if ( s_neqi ( type, 'matrix' ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MM_HEADER_CHECK - Fatal error!'
    write ( *, '(a)' ) '  The value of TYPE was illegal:'
    write ( *, '(a)' ) '    "' // trim ( type ) // '".'
    write ( *, '(a)' ) '  Legal values are:'
    write ( *, '(a)' ) '    "matrix"'
    mm_header_check = .false.
    return
  end if

  if ( &
    s_neqi ( rep, 'coordinate' ) .and. &
    s_neqi ( rep, 'array'      )       &
  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MM_HEADER_CHECK - Fatal error!'
    write ( *, '(a)' ) '  The value of REP was illegal:'
    write ( *, '(a)' ) '    "' // trim ( rep ) // '".'
    write ( *, '(a)' ) '  Legal values are:'
    write ( *, '(a)' ) '    "array"'
    write ( *, '(a)' ) '    "coordinate"'
    mm_header_check = .false.
    return
  end if

  if ( s_eqi ( rep, 'coordinate' ) ) then

    if ( &
      s_neqi ( field, 'integer' ) .and. &
      s_neqi ( field, 'real'    ) .and. &
      s_neqi ( field, 'double'  ) .and. &
      s_neqi ( field, 'complex' ) .and. &
      s_neqi ( field, 'pattern' )       &
    ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MM_HEADER_CHECK - Fatal error!'
      write ( *, '(a)' ) '  The value of FIELD was illegal:'
      write ( *, '(a)' ) '    "' // trim ( field ) // '".'
      mm_header_check = .false.
      return
    end if

  else if ( s_eqi ( rep, 'array' ) ) then

    if ( &
      s_neqi ( field, 'integer' ) .and. &
      s_neqi ( field, 'real'    ) .and. &
      s_neqi ( field, 'double'  ) .and. &
      s_neqi ( field, 'complex' )       &
    ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MM_HEADER_CHECK - Fatal error!'
      write ( *, '(a)' ) '  The value of FIELD was illegal:'
      write ( *, '(a)' ) '    "' // trim ( field ) // '".'
      mm_header_check = .false.
      return
    end if

  end if

  if ( &
    s_neqi ( symm, 'general'        ) .and. &
    s_neqi ( symm, 'symmetric'      ) .and. &
    s_neqi ( symm, 'hermitian'      ) .and. &
    s_neqi ( symm, 'skew-symmetric' )       &
  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MM_HEADER_CHECK - Fatal error!'
    write ( *, '(a)' ) '  The value of SYMM was illegal:'
    write ( *, '(a)' ) '    "' // trim ( symm ) // '".'
    mm_header_check = .false.
    return
  end if

  mm_header_check = .true.

  return
end
subroutine mm_header_print ( input_file, id, type, rep, field, symm )

!*****************************************************************************80
!
!! MM_HEADER_PRINT prints header information from a Matrix Market file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE, the input file name.
!
!    Input, character ( len = 14 ) ID, the Matrix Market identifier.
!    This value must be '%%MatrixMarket'.
!
!    Input, character ( len = 6 ) TYPE, the Matrix Market type.
!    This value must be 'matrix'.
!
!    Input, character ( len = 10 ) REP, the Matrix Market 'representation'
!    indicator.  Possible values include:
!    'coordinate'   (for sparse data)
!    'array'        (for dense data)
!    'elemental'    (to be added)
!
!    Input, character ( len = 7 ) FIELD, the Matrix Market 'field'.
!    Possible values include:
!    'real'
!    'double'
!    'complex'
!    'integer'
!    'pattern'
!
!    Input, character ( len = 19 ) SYMM, the Matrix Market symmetry.
!    Possible values include:
!    'symmetric'
!    'hermitian'
!    'skew-symmetric'
!    'general'
!
  implicit none

  character ( len = 14 ) id
  character ( len = 7 ) field
  character ( len = * ) input_file
  character ( len = 10 ) rep
  character ( len = 19 ) symm
  character ( len = 6 ) type

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MM_HEADER_PRINT:'
  write ( *, '(a)' ) '  Header information from Matrix Market file "' &
    // trim ( input_file ) // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix Market ID = "' &
    // trim ( id ) // '".'
  write ( *, '(a)' ) '    "%%MatrixMarket" is only allowed value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix Market TYPE = "' &
    // trim ( type ) // '".'
  write ( *, '(a)' ) '    "matrix" is only allowed value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Representation type REP = "' &
    // trim ( rep ) // '".'
  write ( *, '(a)' ) '    "coordinate" for sparse matrices,'
  write ( *, '(a)' ) '    "array"      for dense matrices,'
  write ( *, '(a)' ) '    "elemental"  for unassembled finite element matrices.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Numeric FIELD = "' &
    // trim ( field ) // '".'
  write ( *, '(a)' ) '    "integer" for integer values,'
  write ( *, '(a)' ) '    "real"    for real values,'
  write ( *, '(a)' ) '    "double"  for double precision real values,'
  write ( *, '(a)' ) '    "complex" for complex values,'
  write ( *, '(a)' ) '    "pattern" for nonzero pattern only.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Symmetry type SYMM = "' &
    // trim ( symm ) // '".'
  write ( *, '(a)' ) '    "general"         no symmetry,'
  write ( *, '(a)' ) '    "symmetric"       A(I,J) = A(J,I),'
  write ( *, '(a)' ) '                      input only lower triangle.'
  write ( *, '(a)' ) '    "skew-symmetric"  A(I,J) = - A(J,I),'
  write ( *, '(a)' ) '                      input only strict lower triangle.'
  write ( *, '(a)' ) '    "Hermitian"       A(I,J) = A*(J,I),'
  write ( *, '(a)' ) '                      input only lower triangle.'

  return
end
subroutine mm_header_read ( input_unit, id, type, rep, field, symm )

!*****************************************************************************80
!
!! MM_HEADER_READ reads the header line from a Matrix Market file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 2004
!
!  Author:
!
!    Original FORTRAN77 version by Karin A. Remington, NIST ACMD
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT_UNIT, the input unit identifier.
!
!    Output, character ( len = 14 ) ID, the Matrix Market identifier.
!    This value must be '%%MatrixMarket'.
!
!    Output, character ( len = 6 ) TYPE, the Matrix Market type.
!    This value must be 'matrix'.
!
!    Output, character ( len = 10 ) REP, the Matrix Market 'representation'
!    indicator.  Possible values include:
!    'coordinate'   (for sparse data)
!    'array'        (for dense data)
!    'elemental'    (to be added)
!
!    Output, character ( len = 7 ) FIELD, the Matrix Market 'field'.
!    Possible values include:
!    'real'
!    'double'
!    'complex'
!    'integer'
!    'pattern'
!
!    Output, character ( len = 19 ) SYMM, the Matrix Market symmetry.
!    Possible values include:
!    'symmetric'
!    'hermitian'
!    'skew-symmetric'
!    'general'
!
  implicit none

  logical done
  character ( len = 7 ) field
  character ( len = 14 ) id
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 10 ) rep
  character ( len = 1024 ) s
  character ( len = 19 ) symm
  character ( len = 6 ) type
!
!  Read the header line.
!
  read ( input_unit, '(a)', iostat = ios ) s

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MM_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  I/O error reading the header line.'
    stop
  end if

  if ( len_trim ( s ) <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MM_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  The header line does not contain any data.'
    stop
  end if
!
!  Parse the blank-delimited words from the header line:
!
  done = .true.

  call s_w_next ( s, id, done )

  if ( done ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MM_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  The header line does not contain enough data.'
    write ( *, '(a)' ) '  Could not read the ID field.'
    stop
  end if

  call s_w_next ( s, type, done )

  if ( done ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MM_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  The header line does not contain enough data.'
    write ( *, '(a)' ) '  Could not read the TYPE field.'
    stop
  end if

  call s_w_next ( s, rep, done )

  if ( done ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MM_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  The header line does not contain enough data.'
    write ( *, '(a)' ) '  Could not read the REP field.'
    stop
  end if

  call s_w_next ( s, field, done )

  if ( done ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MM_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  The header line does not contain enough data.'
    write ( *, '(a)' ) '  Could not read the FIELD field.'
    stop
  end if

  call s_w_next ( s, symm, done )

  if ( done ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MM_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  The header line does not contain enough data.'
    write ( *, '(a)' ) '  Could not read the SYMM field.'
    stop
  end if

  return
end
subroutine mm_header_write ( output_unit, id, type, rep, field, symm )

!*****************************************************************************80
!
!! MM_HEADER_WRITE prints header information to a Matrix Market file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output unit identifier.
!
!    Input, character ( len = 14 ) ID, the Matrix Market identifier.
!    This value must be '%%MatrixMarket'.
!
!    Input, character ( len = 6 ) TYPE, the Matrix Market type.
!    This value must be 'matrix'.
!
!    Input, character ( len = 10 ) REP, the Matrix Market 'representation'
!    indicator.  Possible values include:
!    'coordinate'   (for sparse data)
!    'array'        (for dense data)
!    'elemental'    (to be added)
!
!    Input, character ( len = 7 ) FIELD, the Matrix Market 'field'.
!    Possible values include:
!    'real'
!    'double'
!    'complex'
!    'integer'
!    'pattern'
!
!    Input, character ( len = 19 ) SYMM, the Matrix Market symmetry.
!    Possible values include:
!    'symmetric'
!    'hermitian'
!    'skew-symmetric'
!    'general'
!
  implicit none

  character ( len = 7 ) field
  character ( len = 14 ) id
  integer ( kind = 4 ) output_unit
  character ( len = 10 ) rep
  character ( len = 19 ) symm
  character ( len = 6 ) type

  write ( output_unit, '(a,1x,a,1x,a,1x,a,1x,a)' ) &
    trim ( id ), trim ( type ), trim ( rep ), trim ( field ), trim ( symm )

  return
end
subroutine mm_nnz_set ( rep, symm, nrow, ncol, nnz )

!*****************************************************************************80
!
!! MM_NNZ_SET sets the value of NNZ for the ARRAY representation.
!
!  Discussion:
!
!    If the representation is not "ARRAY", then NNZ is returned as 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 10 ) REP, the Matrix Market 'representation'
!    indicator.  Possible values include:
!    'coordinate'   (for sparse data)
!    'array'        (for dense data)
!    'elemental'    (to be added)
!
!    Input, character ( len = 19 ) SYMM, the Matrix Market symmetry.
!    Possible values include:
!    'symmetric'
!    'hermitian'
!    'skew-symmetric'
!    'general'
!
!    Input, integer ( kind = 4 ) NROW, the number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns in the matrix.
!
!    Output, integer ( kind = 4 ) NNZ, the number of nonzero entries required
!    to store the matrix, if REP = 'coordinate'.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nnz
  character ( len = 10 ) rep
  logical s_eqi
  character ( len = 19 ) symm

  nnz = 0

  if ( s_eqi ( rep, 'coordinate' ) ) then

  else if ( s_eqi ( rep, 'array' ) ) then

    if ( s_eqi ( symm, 'general' ) ) then
      nnz = nrow * ncol
    else if ( s_eqi ( symm, 'symmetric' ) .or. &
              s_eqi ( symm, 'hermitian' ) &
    ) then
      nnz = ( nrow * ncol - nrow ) / 2 + nrow
    else if ( s_eqi ( symm, 'skew-symmetric' ) ) then
      nnz = ( nrow * ncol - nrow ) / 2
    end if

  end if

  return
end
subroutine mm_size_print ( input_file, rep, symm, nrow, ncol, nnz )

!*****************************************************************************80
!
!! MM_SIZE_PRINT prints size information from a Matrix Market file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE, the input file name.
!
!    Input, character ( len = 10 ) REP, the Matrix Market 'representation'
!    indicator.  Possible values include:
!    'coordinate'   (for sparse data)
!    'array'        (for dense data)
!    'elemental'    (to be added)
!
!    Input, character ( len = 19 ) SYMM, the Matrix Market symmetry.
!    Possible values include:
!    'symmetric'
!    'hermitian'
!    'skew-symmetric'
!    'general'
!
!    Input, integer ( kind = 4 ) NROW, the number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns in the matrix.
!
!    Input, integer ( kind = 4 ) NNZ, the number of nonzero entries required
!    to store the matrix, if REP = 'coordinate'.
!
  implicit none

  character ( len = * ) input_file
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnz
  integer ( kind = 4 ) nnz2
  integer ( kind = 4 ) nrow
  character ( len = 10 ) rep
  logical s_eqi
  character ( len = 17 ) symm

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MM_SIZE_PRINT:'
  write ( *, '(a)' ) '  Size information from Matrix Market file "' &
    // trim ( input_file ) // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of rows    NROW = ', nrow
  write ( *, '(a,i8)' ) '  Number of columns NCOL = ', ncol

  if ( s_eqi ( rep, 'coordinate' ) ) then
    write ( *, '(a,i8)' ) '  Declared number of nonzeros NNZ = ', nnz
  else if ( s_eqi ( rep, 'array' ) ) then
    write ( *, '(a,i8)' ) '  Dummy number of nonzeros NNZ = ', nnz
    call mm_nnz_set ( rep, symm, nrow, ncol, nnz2 )
    write ( *, '(a,i8)' ) '  Inferred number of nonzeros NNZ2 = ', nnz2
  end if

  return
end
subroutine mm_size_read_file ( input_unit, rep, symm, nrow, ncol, nnz )

!*****************************************************************************80
!
!! MM_SIZE_READ_FILE reads size information from a Matrix Market file.
!
!  Discussion:
!
!    This routine assumes that the next line of the file is the size
!    record.  However, a Matrix Market file may contain comment lines
!    between the first record and the size record.  In that case, you
!    might prefer using MM_COMMENT_READ to read each line until you
!    find the first noncomment (comment lines must begin with the
!    percent character) and then passing that line to MM_SIZE_READ_STRING.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT_FILE, the input file unit.
!
!    Input, character ( len = 10 ) REP, the Matrix Market 'representation'
!    indicator.  Possible values include:
!    'coordinate'   (for sparse data)
!    'array'        (for dense data)
!    'elemental'    (to be added)
!
!    Input, character ( len = 19 ) SYMM, the Matrix Market symmetry.
!    Possible values include:
!    'symmetric'
!    'hermitian'
!    'skew-symmetric'
!    'general'
!
!    Output, integer NROW, the number of rows in the matrix.
!
!    Output, integer NCOL, the number of columns in the matrix.
!
!    Output, integer NNZ, the number of nonzero entries required to store
!    the matrix, if REP = 'coordinate'.
!
  implicit none

  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnz
  integer ( kind = 4 ) nrow
  character ( len = 10 ) rep
  logical s_eqi
  character ( len = 19 ) symm

  if ( s_eqi ( rep, 'coordinate' ) ) then

    read ( input_unit, *, iostat = ios ) nrow, ncol, nnz

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MM_SIZE_READ_FILE - Fatal error!'
      write ( *, '(a)' ) '  I/O error while reading size information.'
      write ( *, '(a)' ) '  Expecting "NROW NCOL NNZ" values.'
      stop
    end if

  else if ( s_eqi ( rep, 'array' ) ) then

    read ( input_unit, *, iostat = ios ) nrow, ncol

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MM_SIZE_READ_FILE - Fatal error!'
      write ( *, '(a)' ) '  I/O error while reading size information.'
      write ( *, '(a)' ) '  Expecting "NROW NCOL" values.'
      stop
    end if

    call mm_nnz_set ( rep, symm, nrow, ncol, nnz )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MM_SIZE_READ_FILE - Fatal error!'
    write ( *, '(a)' ) '  Unexpected value of REP.'
    stop

  end if

  return
end
subroutine mm_size_read_string ( string, rep, symm, nrow, ncol, nnz )

!*****************************************************************************80
!
!! MM_SIZE_READ_STRING reads size information from a string.
!
!  Discussion:
!
!    We presume that the string contains the size record from a
!    Matrix Market file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the string containing the size record.
!
!    Input, character ( len = 10 ) REP, the Matrix Market 'representation'
!    indicator.  Possible values include:
!    'coordinate'   (for sparse data)
!    'array'        (for dense data)
!    'elemental'    (to be added)
!
!    Input, character ( len = 19 ) SYMM, the Matrix Market symmetry.
!    Possible values include:
!    'symmetric'
!    'hermitian'
!    'skew-symmetric'
!    'general'
!
!    Output, integer ( kind = 4 ) NROW, the number of rows in the matrix.
!
!    Output, integer ( kind = 4 ) NCOL, the number of columns in the matrix.
!
!    Output, integer ( kind = 4 ) NNZ, the number of nonzero entries required
!    to store the matrix, if REP = 'coordinate'.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnz
  integer ( kind = 4 ) nrow
  character ( len = 10 ) rep
  logical s_eqi
  character ( len = * ) string
  character ( len = 19 ) symm
  integer ( kind = 4 ) value(3)

  if ( s_eqi ( rep, 'coordinate' ) ) then

    call s_to_i4vec ( string, 3, value, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MM_SIZE_READ_STRING - Fatal error!'
      write ( *, '(a)' ) '  I/O error while reading size information.'
      write ( *, '(a)' ) '  Expecting "NROW NCOL NNZ" values.'
      stop
    end if

    nrow = value(1)
    ncol = value(2)
    nnz = value(3)

  else if ( s_eqi ( rep, 'array' ) ) then

    call s_to_i4vec ( string, 2, value, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MM_SIZE_READ_STRING - Fatal error!'
      write ( *, '(a)' ) '  I/O error while reading size information.'
      write ( *, '(a)' ) '  Expecting "NROW NCOL" values.'
      stop
    end if

    nrow = value(1)
    ncol = value(2)

    call mm_nnz_set ( rep, symm, nrow, ncol, nnz )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MM_SIZE_READ_STRING - Fatal error!'
    write ( *, '(a)' ) '  Unexpected value of REP.'
    stop

  end if

  return
end
subroutine mm_size_write ( output_unit, rep, nrow, ncol, nnz )

!*****************************************************************************80
!
!! MM_SIZE_WRITE writes size information to a Matrix Market file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the input file unit.
!
!    Input, character ( len = 10 ) REP, the Matrix Market 'representation'
!    indicator.  Possible values include:
!    'coordinate'   (for sparse data)
!    'array'        (for dense data)
!    'elemental'    (to be added)
!
!    Input, integer ( kind = 4 ) NROW, the number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns in the matrix.
!
!    Input, integer ( kind = 4 ) NNZ, the number of nonzero entries required
!    to store the matrix, if REP = 'coordinate'.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnz
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) output_unit
  character ( len = 10 ) rep
  logical s_eqi

  if ( s_eqi ( rep, 'coordinate' ) ) then
    write ( output_unit, * ) nrow, ncol, nnz
  else if ( s_eqi ( rep, 'array' ) ) then
    write ( output_unit, * ) nrow, ncol
  end if

  return
end
subroutine mm_values_print ( rep, field, nnz, indx, jndx, ival, rval, &
  dval, cval )

!*****************************************************************************80
!
!! MM_VALUES_PRINT prints the matrix values of a Matrix Market file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 10 ) REP, the Matrix Market 'representation'
!    indicator.  Possible values include:
!    'coordinate'   (for sparse data)
!    'array'        (for dense data)
!    'elemental'    (to be added)
!
!    Input, character ( len = 7 ) FIELD, the Matrix Market 'field'.
!    Possible values include:
!    'real'
!    'double'
!    'complex'
!    'integer'
!    'pattern'
!
!    Input, integer ( kind = 4 ) NNZ, the number of nonzero entries required to
!    store the matrix, if REP = 'coordinate'.
!
!    Input, integer ( kind = 4 ) INDX(NNZ), the row indices for coordinate format.
!    Not used if REP is 'array'.
!
!    Input, integer ( kind = 4 ) JNDX(NNZ), the column indices for coordinate
!    format.  Not used if REP is 'array'.
!
!    Input, integer ( kind = 4 ) IVAL(NNZ), matrix values, if FIELD is 'integer'.
!
!    Input, real ( kind = 4 ) RVAL(NNZ), matrix values, if FIELD is 'real'.
!
!    Input, real ( kind = 8 ) DVAL(NNZ), matrix values, if FIELD is 'double'.
!
!    Input, complex ( kind = 4 ) CVAL(NNZ), matrix values, if FIELD is 'complex'.
!
  implicit none

  integer nnz

  complex ( kind = 4 ) cval(nnz)
  real ( kind = 8 ) dval(nnz)
  character ( len = 7 ) field
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) indx(nnz)
  integer ( kind = 4 ) ival(nnz)
  integer ( kind = 4 ) jndx(nnz)
  character ( len = 10 ) rep
  real ( kind = 4 ) rval(nnz)

  ilo = 1
  ihi = nnz

  call mm_values_print_some ( rep, field, nnz, indx, jndx, ival, rval, &
    dval, cval, ilo, ihi )

  return
end
subroutine mm_values_print_some ( rep, field, nnz, indx, jndx, ival, rval, &
  dval, cval, ilo, ihi )

!*****************************************************************************80
!
!! MM_VALUES_PRINT_SOME prints some matrix values of a Matrix Market file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 10 ) REP, the Matrix Market 'representation'
!    indicator.  Possible values include:
!    'coordinate'   (for sparse data)
!    'array'        (for dense data)
!    'elemental'    (to be added)
!
!    Input, character ( len = 7 ) FIELD, the Matrix Market 'field'.
!    Possible values include:
!    'real'
!    'double'
!    'complex'
!    'integer'
!    'pattern'
!
!    Input, integer ( kind = 4 ) NNZ, the number of nonzero entries required to
!    store the matrix, if REP = 'coordinate'.
!
!    Input, integer ( kind = 4 ) INDX(NNZ), the row indices for coordinate format.
!    Not used if REP is 'array'.
!
!    Input, integer ( kind = 4 ) JNDX(NNZ), the column indices for coordinate
!    format.  Not used if REP is 'array'.
!
!    Input, integer ( kind = 4 ) IVAL(NNZ), matrix values, if FIELD is 'integer'.
!
!    Input, real ( kind = 4 ) RVAL(NNZ), matrix values, if FIELD is 'real'.
!
!    Input, real ( kind = 8 ) DVAL(NNZ), matrix values, if FIELD is 'double'.
!
!    Input, complex ( kind = 4 ) CVAL(NNZ), matrix values, if FIELD is 'complex'.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the minimum and maximum indices of
!    the data to print.
!
  implicit none

  integer nnz

  complex ( kind = 4 ) cval(nnz)
  real ( kind = 8 ) dval(nnz)
  character ( len = 7 ) field
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) indx(nnz)
  integer ( kind = 4 ) ival(nnz)
  integer ( kind = 4 ) jndx(nnz)
  character ( len = 10 ) rep
  real ( kind = 4 ) rval(nnz)
  logical s_eqi

  if ( s_eqi ( rep, 'coordinate' ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Sparse array storage by coordinate.'
    write ( *, '(a,i8,a,i8)' ) '  Listing entries ', ilo, ' through ', ihi
    write ( *, '(a)' ) ' '
  else if ( s_eqi ( rep, 'array' ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Dense array storage of all elements.'
    write ( *, '(a,i8,a,i8)' ) '  Listing entries ', ilo, ' through ', ihi
    write ( *, '(a)' ) ' '
  end if

  do i = max ( ilo, 1 ), min ( ihi, nnz )

    if ( s_eqi ( rep, 'coordinate' ) ) then

      if ( s_eqi ( field, 'integer' ) ) then
        write ( *, '(i4,2x,3i8)' ) i, indx(i), jndx(i), ival(i)
      else if ( s_eqi ( field, 'real' ) ) then
        write ( *, '(i4,2x,2i8,g14.6)' ) i, indx(i), jndx(i), rval(i)
      else if ( s_eqi ( field, 'double' ) ) then
        write ( *, '(i4,2x,2i8,g14.6)' ) i, indx(i), jndx(i), dval(i)
      else if ( s_eqi ( field, 'complex' ) ) then
        write ( *, '(i4,2x,2i8,2g14.6)' ) i, indx(i), jndx(i), cval(i)
      else if ( s_eqi ( field, 'pattern' ) ) then
        write ( *, '(i4,2x,2i8)' ) i, indx(i), jndx(i)
      end if

    else if ( s_eqi ( rep, 'array' ) ) then

      if ( s_eqi ( field, 'integer' ) ) then
        write ( *, '(i4,2x,i8)' ) i, ival(i)
      else if ( s_eqi ( field, 'real' ) ) then
        write ( *, '(i4,2x,g14.6)' ) i, rval(i)
      else if ( s_eqi ( field, 'double' ) ) then
        write ( *, '(i4,2x,g14.6)' ) i, dval(i)
      else if ( s_eqi ( field, 'complex' ) ) then
        write ( *, '(i4,2x,2g14.6)' ) i, cval(i)
      end if

    end if

  end do

  return
end
subroutine mm_values_read ( input_unit, rep, field, nnz, indx, jndx, &
  ival, rval, dval, cval )

!*****************************************************************************80
!
!! MM_VALUES_READ reads matrix values from a Matrix Market file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT_UNIT, the input unit identifier.
!
!    Input, character ( len = 10 ) REP, the Matrix Market 'representation'
!    indicator.  Possible values include:
!    'coordinate'   (for sparse data)
!    'array'        (for dense data)
!    'elemental'    (to be added)
!
!    Input, character ( len = 7 ) FIELD, the Matrix Market 'field'.
!    Possible values include:
!    'real'
!    'double'
!    'complex'
!    'integer'
!    'pattern'
!
!    Input, integer ( kind = 4 ) NNZ, the number of nonzero entries required
!    to store the matrix, if REP = 'coordinate'.
!
!    Output, integer ( kind = 4 ) INDX(NNZ), the row indices for coordinate
!    format.  Not used if REP is 'array'.
!
!    Output, integer ( kind = 4 ) JNDX(NNZ), the column indices for coordinate
!    format.  Not used if REP is 'array'.
!
!    Output, integer ( kind = 4 ) IVAL(NNZ), matrix values, if FIELD is 'integer'.
!
!    Output, real ( kind = 4 ) RVAL(NNZ), matrix values, if FIELD is 'real'.
!
!    Output, real ( kind = 8 ) DVAL(NNZ), matrix values, if FIELD is 'double'.
!
!    Output, complex ( kind = 4 ) CVAL(NNZ), matrix values, if FIELD is 'complex'.
!
  implicit none

  integer ( kind = 4 ) nnz

  complex ( kind = 4 ) cval(nnz)
  real ( kind = 8 ) dval(nnz)
  character ( len = 7 ) field
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(nnz)
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ival(nnz)
  integer ( kind = 4 ) jndx(nnz)
  character ( len = 10 ) rep
  real ( kind = 4 ) rval(nnz)
  logical s_eqi

  if ( s_eqi ( rep, 'coordinate' ) ) then

    if ( s_eqi ( field, 'integer' ) ) then
      do i = 1, nnz
        read ( input_unit, *, iostat = ios ) indx(i), jndx(i), ival(i)
        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'MM_VALUES_READ - Fatal error!'
          write ( *, '(a)' ) '  Error or end of file on value field ', i
          stop
        end if
      end do
    else if ( s_eqi ( field, 'real' ) ) then
      do i = 1, nnz
        read ( input_unit, *, iostat = ios ) indx(i), jndx(i), rval(i)
        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'MM_VALUES_READ - Fatal error!'
          write ( *, '(a)' ) '  Error or end of file on value field ', i
          stop
        end if
      end do
    else if ( s_eqi ( field, 'double' ) ) then
      do i = 1, nnz
        read ( input_unit, *, iostat = ios ) indx(i), jndx(i), dval(i)
        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'MM_VALUES_READ - Fatal error!'
          write ( *, '(a)' ) '  Error or end of file on value field ', i
          stop
        end if
      end do
    else if ( s_eqi ( field, 'complex' ) ) then
      do i = 1, nnz
        read ( input_unit, *, iostat = ios ) indx(i), jndx(i), cval(i)
        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'MM_VALUES_READ - Fatal error!'
          write ( *, '(a)' ) '  Error or end of file on value field ', i
          stop
        end if
      end do
    else if ( s_eqi ( field, 'pattern' ) ) then
      do i = 1, nnz
        read ( input_unit, *, iostat = ios ) indx(i), jndx(i)
        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'MM_VALUES_READ - Fatal error!'
          write ( *, '(a)' ) '  Error or end of file on value field ', i
          stop
        end if
      end do
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MM_VALUES_READ - Fatal error!'
      write ( *, '(a)' ) '  Illegal value of FIELD.'
      stop
    end if

  else if ( s_eqi ( rep == 'array' ) ) then

    if ( s_eqi ( field, 'integer' ) ) then
      do i = 1, nnz
        read ( input_unit, *, iostat = ios ) ival(i)
        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'MM_VALUES_READ - Fatal error!'
          write ( *, '(a)' ) '  Error or end of file on value field ', i
          stop
        end if
      end do
    else if ( s_eqi ( field, 'real' ) ) then
      do i = 1, nnz
        read ( input_unit, *, iostat = ios ) rval(i)
        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'MM_VALUES_READ - Fatal error!'
          write ( *, '(a)' ) '  Error or end of file on value field ', i
          stop
        end if
      end do
    else if ( s_eqi ( field, 'double' ) ) then
      do i = 1, nnz
        read ( input_unit, *, iostat = ios ) dval(i)
        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'MM_VALUES_READ - Fatal error!'
          write ( *, '(a)' ) '  Error or end of file on value field ', i
          stop
        end if
      end do
    else if ( s_eqi ( field, 'complex' ) ) then
      do i = 1, nnz
        read ( input_unit, *, iostat = ios ) cval(i)
        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'MM_VALUES_READ - Fatal error!'
          write ( *, '(a)' ) '  Error or end of file on value field ', i
          stop
        end if
      end do
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MM_VALUES_READ - Fatal error!'
      write ( *, '(a)' ) '  Illegal value of FIELD.'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MM_VALUES_READ - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of REP.'
    stop

  end if

  return
end
subroutine mm_values_write ( output_unit, rep, field, nnz, indx, jndx, &
  ival, rval, dval, cval )

!*****************************************************************************80
!
!! MM_VALUES_WRITE writes matrix values to a Matrix Market file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output unit identifier.
!
!    Input, character ( len = 10 ) REP, the Matrix Market 'representation'
!    indicator.  Possible values include:
!    'coordinate'   (for sparse data)
!    'array'        (for dense data)
!    'elemental'    (to be added)
!
!    Input, character ( len = 7 ) FIELD, the Matrix Market 'field'.
!    Possible values include:
!    'real'
!    'double'
!    'complex'
!    'integer'
!    'pattern'
!
!    Input, integer ( kind = 4 ) NNZ, the number of nonzero entries required
!    to store the matrix, if REP = 'coordinate'.
!
!    Input, integer ( kind = 4 ) INDX(NNZ), the row indices for coordinate
!    format.  Not used if REP is 'array'.
!
!    Input, integer ( kind = 4 ) JNDX(NNZ), the column indices for coordinate
!    format.  Not used if REP is 'array'.
!
!    Input, integer ( kind = 4 ) IVAL(NNZ), matrix values, if FIELD is 'integer'.
!
!    Input, real ( kind = 4 ) RVAL(NNZ), matrix values, if FIELD is 'real'.
!
!    Input, real ( kind = 8 ) DVAL(NNZ), matrix values, if FIELD is 'double'.
!
!    Input, complex ( kind = 4 ) CVAL(NNZ), matrix values, if FIELD is 'complex'.
!
  implicit none

  integer ( kind = 4 ) nnz

  complex ( kind = 4 ) cval(nnz)
  real ( kind = 8 ) dval(nnz)
  character ( len = 7 ) field
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(nnz)
  integer ( kind = 4 ) ival(nnz)
  integer ( kind = 4 ) jndx(nnz)
  integer ( kind = 4 ) output_unit
  character ( len = 10 ) rep
  real ( kind = 4 ) rval(nnz)
  logical s_eqi

  if ( s_eqi ( rep, 'coordinate' ) ) then

    if ( s_eqi ( field, 'integer' ) ) then
      do i = 1, nnz
        write ( output_unit, * ) indx(i), jndx(i), ival(i)
      end do
    else if ( s_eqi ( field, 'real' ) ) then
      do i = 1, nnz
        write ( output_unit, * ) indx(i), jndx(i), rval(i)
      end do
    else if ( s_eqi ( field, 'double' ) ) then
      do i = 1, nnz
        write ( output_unit, * ) indx(i), jndx(i), dval(i)
      end do
    else if ( s_eqi ( field, 'complex' ) ) then
      do i = 1, nnz
        write ( output_unit, * ) indx(i), jndx(i), cval(i)
      end do
    else if ( s_eqi ( field, 'pattern' ) ) then
      do i = 1, nnz
        write ( output_unit, * ) indx(i), jndx(i)
      end do
    end if

  else if ( s_eqi ( rep, 'array' ) ) then

    if ( s_eqi ( field, 'integer' ) ) then
      do i = 1, nnz
        write ( output_unit, * ) ival(i)
      end do
    else if ( s_eqi ( field, 'real' ) ) then
      do i = 1, nnz
        write ( output_unit, * ) rval(i)
      end do
    else if ( s_eqi ( field, 'double' ) ) then
      do i = 1, nnz
        write ( output_unit, * ) dval(i)
      end do
    else if ( s_eqi ( field, 'complex' ) ) then
      do i = 1, nnz
        write ( output_unit, * ) cval(i)
      end do
    end if

  end if

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
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

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character ( len = * ) s1
  character ( len = * ) s2

  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )

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

  do i = lenc + 1, len1
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
function s_neqi ( s1, s2 )

!*****************************************************************************80
!
!! S_NEQI compares two strings for non-equality, ignoring case.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_NEQI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_neqi
  character ( len = * ) s1
  character ( len = * ) s2

  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )

  s_neqi = .true.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc+1, len1
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc+1, len2
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_neqi = .false.

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
subroutine s_w_next ( s, word, done )

!*****************************************************************************80
!
!! S_W_NEXT "reads" words from a string, one at a time.
!
!  Special cases:
!
!    The following characters are considered to be a single word,
!    whether surrounded by spaces or not:
!
!      " ( ) { } [ ]
!
!    Also, if there is a trailing comma on the word, it is stripped off.
!    This is to facilitate the reading of lists.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string, presumably containing words
!    separated by spaces.
!
!    Output, character ( len = * ) WORD.
!
!    If DONE is FALSE, then WORD contains the "next" word read.
!    If DONE is TRUE, then WORD is blank, because there was no more to read.
!
!    Input/output, logical DONE.
!
!    On input with a fresh string, set DONE to TRUE.
!
!    On output, the routine sets DONE:
!      FALSE if another word was read,
!      TRUE if no more words could be read.
!
  implicit none

  logical done
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ), save :: lenc = 0
  integer ( kind = 4 ), save :: next = 1
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
  character ( len = * ) word
!
!  We "remember" LENC and NEXT from the previous call.
!
!  An input value of DONE = TRUE signals a new line of text to examine.
!
  if ( done ) then

    next = 1
    done = .false.
    lenc = len_trim ( s )

    if ( lenc <= 0 ) then
      done = .true.
      word = ' '
      return
    end if

  end if
!
!  Beginning at index NEXT, search the string for the next nonblank,
!  which signals the beginning of a word.
!
  ilo = next
!
!  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
!
  do

    if ( lenc < ilo ) then
      word = ' '
      done = .true.
      next = lenc + 1
      return
    end if
!
!  If the current character is blank, skip to the next one.
!
    if ( s(ilo:ilo) /= ' ' .and. s(ilo:ilo) /= TAB ) then
      exit
    end if

    ilo = ilo + 1

  end do
!
!  ILO is the index of the next nonblank character in the string.
!
!  If this initial nonblank is a special character,
!  then that's the whole word as far as we're concerned,
!  so return immediately.
!
  if ( s(ilo:ilo) == '"' .or. &
       s(ilo:ilo) == '(' .or. &
       s(ilo:ilo) == ')' .or. &
       s(ilo:ilo) == '{' .or. &
       s(ilo:ilo) == '}' .or. &
       s(ilo:ilo) == '[' .or. &
       s(ilo:ilo) == ']' ) then

    word = s(ilo:ilo)
    next = ilo + 1
    return

  end if
!
!  Now search for the last contiguous character that is not a
!  blank, TAB, or special character.
!
  next = ilo + 1

  do while ( next <= lenc )

    if ( s(next:next) == ' ' ) then
      exit
    else if ( s(next:next) == TAB ) then
      exit
    else if ( s(next:next) == '"' ) then
      exit
    else if ( s(next:next) == '(' ) then
      exit
    else if ( s(next:next) == ')' ) then
      exit
    else if ( s(next:next) == '{' ) then
      exit
    else if ( s(next:next) == '}' ) then
      exit
    else if ( s(next:next) == '[' ) then
      exit
    else if ( s(next:next) == ']' ) then
      exit
    end if

    next = next + 1

  end do
!
!  Ignore a trailing comma.
!
  if ( s(next-1:next-1) == ',' ) then
    word = s(ilo:next-2)
  else
    word = s(ilo:next-1)
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
