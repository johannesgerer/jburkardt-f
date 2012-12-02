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
function ch_is_space ( ch )

!*****************************************************************************80
!
!! CH_IS_SPACE is TRUE if a character is a whitespace character.
!
!  Discussion:
!
!    Instead of CHAR, we now use the ACHAR function, which
!    guarantees the ASCII collating sequence.
!
!    A whitespace character is a space, a form feed, a newline,
!    a carriage return, a tab, or a vertical tab.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character CH, a character to check.
!
!    Output, logical CH_IS_SPACE is TRUE if the character is a whitespace
!    character.
!
  implicit none

  character ch
  logical ch_is_space

  if ( ch == ' ' ) then
    ch_is_space = .true.
  else if ( ch == achar ( 12 ) ) then
    ch_is_space = .true.
  else if ( ch == achar ( 10 ) ) then
    ch_is_space = .true.
  else if ( ch == achar ( 13 ) ) then
    ch_is_space = .true.
  else if ( ch == achar ( 9 ) ) then
    ch_is_space = .true.
  else if ( ch == achar ( 11 ) ) then
    ch_is_space = .true.
  else
    ch_is_space = .false.
  end if

  return
end
subroutine cnf_data_read ( cnf_file_name, v_num, c_num, l_num, l_c_num, l_val )

!*****************************************************************************80
!
!! CNF_DATA_READ reads the data of a CNF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CNF_FILE_NAME, the name of the CNF file.
!
!    Input, integer ( kind = 4 ) V_NUM, the number of variables.
!
!    Input, integer ( kind = 4 ) C_NUM, the number of clauses.
!
!    Input, integer ( kind = 4 ) L_NUM, the number of signed literals.
!
!    Output, integer ( kind = 4 ) L_C_NUM(C_NUM), the number of signed
!    literals occuring in each clause.
!
!    Output, integer ( kind = 4 ) L_VAL(L_NUM), a list of all the signed
!    literals in all the clauses, ordered by clause.
!
  implicit none

  integer ( kind = 4  ) c_num
  integer ( kind = 4  ) l_num
  integer ( kind = 4  ) v_num

  integer ( kind = 4  ) c_num2
  logical ch_eqi
  logical ch_is_space
  character ( len = *   ) cnf_file_name
  integer ( kind = 4  ) cnf_file_status
  integer ( kind = 4  ) cnf_file_unit
  integer ( kind = 4  ) ierror
  integer ( kind = 4  ) l_c_num(c_num)
  integer ( kind = 4  ) l_c_num2
  integer ( kind = 4  ) l_num2
  integer ( kind = 4  ) l_val(l_num)
  integer ( kind = 4  ) l_val2
  integer ( kind = 4  ) length
  character ( len = 255 ) line
  logical                 s_eqi
  integer ( kind = 4  ) v_num2

  character ( len = 20  ) word

  call get_unit ( cnf_file_unit )

  open ( unit = cnf_file_unit, file = cnf_file_name, status = 'old', &
    iostat = cnf_file_status )

  if ( cnf_file_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open file!'
    stop
  end if
!
!  Read lines until you find one that is not blank and does not begin
!  with a "c".  This should be the header line.
!
  line = ' '

  do

    read ( cnf_file_unit, '(a)', iostat = cnf_file_status ) line

    if ( cnf_file_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading the file.'
      stop
    end if

    if ( line(1:1) == 'c' .or. line(1:1) == 'C' ) then
      cycle
    end if

    if ( 0 < len_trim ( line ) ) then
      exit
    end if

  end do
!
!  We expect to be reading the line "p cnf V_NUM C_NUM"
!
  if ( .not. ch_eqi ( line(1:1), 'p' ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  First non-comment non-blank line does not start'
    write ( *, '(a)' ) '  with "p " marker.'
    stop
  end if

  if ( .not. ch_is_space ( line(2:2) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Character after "p" must be whitespace.'
    stop
  end if
!
!  Remove the first two characters and shift left to first nonblank.
!
  line(1:1) = ' '
  line(2:2) = ' '
  line = adjustl ( line )
!
!  Expect the string 'CNF'
!
  if ( .not. s_eqi ( line(1:3), 'cnf' ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  First non-comment non-blank line does not start'
    write ( *, '(a)' ) '  with "p cnf" marker.'
    stop
  end if

  if ( .not. ch_is_space ( line(4:4) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Character after "p cnf" must be whitespace.'
    stop
  end if
!
!  Remove the first four characters and shift left.
!
  line(1:4) = ' '
  line = adjustl ( line )
!
!  Extract the next word, which is the number of variables.
!  You can compare this to V_NUM for an extra check.
!
  call s_word_extract_first ( line, word )

  if ( len_trim ( word ) <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Unexpected End of input.'
    stop
  end if

  call s_to_i4 ( word, v_num2, ierror, length )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Unexpected End of input.'
    stop
  end if
!
!  Extract the next word, which is the number of clauses.
!  You can compare this to C_NUM for an extra check.
!
  call s_word_extract_first ( line, word )

  if ( len_trim ( word ) == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Unexpected End of input.'
    stop
  end if

  call s_to_i4 ( word, c_num2, ierror, length )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Unexpected End of input.'
    stop
  end if
!
!  Read remaining lines, counting the literals, ignoring occurrences of '0'.
!
  l_num2 = 0
  c_num2 = 0
  l_c_num2 = 0
  line = ' '

  do

    read ( cnf_file_unit, '(a)', iostat = cnf_file_status ) line

    if ( cnf_file_status /= 0 ) then
      exit
    end if

    if ( line(1:1) == 'c' ) then
      cycle
    end if

    if ( len_trim ( line ) < 0 ) then
      exit
    end if

    do

      call s_word_extract_first ( line, word )

      if ( len_trim ( word ) <= 0 ) then
        exit
      end if

      call s_to_i4 ( word, l_val2, ierror, length )

      if ( ierror /= 0 ) then
        exit
      end if

      if ( l_val2 /= 0 ) then
        l_num2 = l_num2 + 1
        l_val(l_num2) = l_val2
        l_c_num2 = l_c_num2 + 1
      else
        c_num2 = c_num2 + 1
        l_c_num(c_num2) = l_c_num2
        l_c_num2 = 0
      end if

    end do

  end do
!
!  At the end:
!
!    C_NUM2 should equal C_NUM,
!    L_NUM2 should equal L_NUM.
!
!  Close file and return.
!
  close ( unit = cnf_file_unit )

  return
end
subroutine cnf_data_write ( c_num, l_num, l_c_num, l_val, output_unit )

!*****************************************************************************80
!
!! CNF_DATA_WRITE writes data to a CNF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) C_NUM, the number of clauses.
!
!    Input, integer ( kind = 4 ) L_NUM, the total number of signed literals.
!
!    Input, integer ( kind = 4 ) L_C_NUM(C_NUM), the number of signed
!    literals occuring in each clause.
!
!    Input, integer ( kind = 4 ) L_VAL(L_NUM), a list of all the signed
!    literals in all the clauses, ordered by clause.
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output unit.
!
  implicit none

  integer ( kind = 4 ) c_num
  integer ( kind = 4 ) l_num

  integer ( kind = 4 ) c
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_c
  integer ( kind = 4 ) l_c_num(c_num)
  integer ( kind = 4 ) l_val(l_num)
  integer ( kind = 4 ) output_unit
  character ( len = 80 ) string

  l = 0
  string = ' '

  do c = 1, c_num

    i1 = 1
    i2 = 10

    do l_c = 1, l_c_num(c)

      l = l + 1

      write ( string(i1:i2), '(1x,i7)' ) l_val(l)

      i1 = i1 + 10
      i2 = i2 + 10

      if ( mod ( l_c, 10 ) == 0 ) then
        call s_blanks_delete ( string )
        write ( output_unit, '(a)' ) string(1:len_trim(string))
        string = ' '
      end if

    end do

    string(i2+1:i2+2) = ' 0'
    call s_blanks_delete ( string )
    write ( output_unit, '(a)' ) string(1:len_trim(string))
    string = ' '

  end do

  return
end
function cnf_evaluate ( v_num, c_num, l_num, l_c_num, l_val, v_val )

!*****************************************************************************80
!
!! CNF_EVALUATE evaluates a formula in CNF form.
!
!  Discussion:
!
!    The formula is in conjunctive normal form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) V_NUM, the number of variables.
!
!    Input, integer ( kind = 4 ) C_NUM, the number of clauses.
!
!    Input, integer ( kind = 4 ) L_NUM, the total number of signed literals.
!
!    Input, integer ( kind = 4 ) L_C_NUM(C_NUM), the number of signed
!    literals occuring in each clause.
!
!    Input, integer ( kind = 4 ) L_VAL(L_NUM), a list of all the signed
!    literals in all the clauses, ordered by clause.
!
!    Input, logical V_VAL(V_NUM), the values assigned to the variables.
!
!    Output, logical CNF_EVALUATE, the value of the CNF formula for the
!    given variable values.
!
  implicit none

  integer ( kind = 4 ) c_num
  logical              cnf_evaluate
  integer ( kind = 4 ) l_num
  integer ( kind = 4 ) v_num

  integer ( kind = 4 ) c
  logical              c_val
  logical              f_val
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_c
  integer ( kind = 4 ) l_c_num(c_num)
  integer ( kind = 4 ) l_val(l_num)
  logical              s_val
  integer ( kind = 4 ) v_index
  logical              v_val(v_num)
  logical              value

  f_val = .true.

  l = 0

  do c = 1, c_num
!
!  The clause is false unless some signed literal is true.
!
    c_val = .false.
    do l_c = 1, l_c_num(c)
      l = l + 1
      s_val = ( 0 < l_val(l) )
      v_index = abs ( l_val(l) )
!
!  The signed literal is true if the sign "equals" the value.
!  Note that we CAN'T exit the loop because we need to run out the
!  L index!
!
      if ( v_val(v_index) .eqv. s_val ) then
        c_val = .true.
      end if
    end do
!
!  The formula is false if any clause is false.
!
    if ( .not. c_val ) then
      f_val = .false.
      exit
    end if

  end do

  cnf_evaluate = f_val

  return
end
subroutine cnf_header_read ( cnf_file_name, v_num, c_num, l_num )

!*****************************************************************************80
!
!! CNF_HEADER_READ reads the header of a CNF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 June 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CNF_FILE_NAME, the name of the CNF file.
!
!    Output, integer ( kind = 4 ) V_NUM, the number of variables.
!
!    Output, integer ( kind = 4 ) C_NUM, the number of clauses.
!
!    Output, integer ( kind = 4 ) L_NUM, the number of signed literals.
!
  implicit none

  integer ( kind = 4  ) c_num
  logical                 ch_eqi
  logical                 ch_is_space
  character ( len = *   ) cnf_file_name
  integer ( kind = 4  ) cnf_file_status
  integer ( kind = 4  ) cnf_file_unit
  integer ( kind = 4  ) ierror
  integer ( kind = 4  ) l_num
  integer ( kind = 4  ) l_val
  integer ( kind = 4  ) length
  character ( len = 255 ) line
  logical                 s_eqi
  integer ( kind = 4  ) v_num
  character ( len = 20  ) word

  call get_unit ( cnf_file_unit )

  open ( unit = cnf_file_unit, file = cnf_file_name, status = 'old', &
    iostat = cnf_file_status )

  if ( cnf_file_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open file!'
    stop
  end if
!
!  Read lines until you find one that is not blank and does not begin
!  with a "c".  This should be the header line.
!
  line = ' '

  do

    read ( cnf_file_unit, '(a)', iostat = cnf_file_status ) line

    if ( cnf_file_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading the file.'
      stop
    end if

    if ( line(1:1) == 'c' .or. line(1:1) == 'C' ) then
      cycle
    end if

    if ( 0 < len_trim ( line ) ) then
      exit
    end if

  end do
!
!  We expect to be reading the line "p cnf V_NUM C_NUM"
!
  if ( .not. ch_eqi ( line(1:1), 'p' ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  First non-comment non-blank line does not start'
    write ( *, '(a)' ) '  with "p " marker.'
    stop
  end if

  if ( .not. ch_is_space ( line(2:2) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  Character after "p" must be whitespace.'
    stop
  end if
!
!  Remove the first two characters and shift left.
!
  line(1:2) = ' '
  line = adjustl ( line )
!
!  Expect the string 'CNF'
!
  if ( .not. s_eqi ( line(1:3), 'cnf' ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  First non-comment non-blank line does not start'
    write ( *, '(a)' ) '  with "p cnf" marker.'
    stop
  end if

  if ( .not. ch_is_space ( line(4:4) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  Character after "p cnf" must be whitespace.'
    stop
  end if
!
!  Remove the first four characters and shift left.
!
  line(1:4) = ' '
  line = adjustl ( line )
!
!  Extract the next word, which is the number of variables.
!
  call s_word_extract_first ( line, word )

  if ( len_trim ( word ) <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  Unexpected End of input.'
    stop
  end if

  call s_to_i4 ( word, v_num, ierror, length )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  Unexpected End of input.'
    stop
  end if
!
!  Extract the next word, which is the number of clauses.
!
  call s_word_extract_first ( line, word )

  if ( len_trim ( word ) == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  Unexpected End of input.'
    stop
  end if

  call s_to_i4 ( word, c_num, ierror, length )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  Unexpected End of input.'
    stop
  end if
!
!  Read remaining lines, counting the literals, ignoring occurrences of '0'.
!
  l_num = 0
  line = ' '

  do

    read ( cnf_file_unit, '(a)', iostat = cnf_file_status ) line

    if ( cnf_file_status /= 0 ) then
      exit
    end if

    if ( line(1:1) == 'c' ) then
      cycle
    end if

    if ( len_trim ( line ) < 0 ) then
      exit
    end if

    do

      call s_word_extract_first ( line, word )

      if ( len_trim ( word ) <= 0 ) then
        exit
      end if

      call s_to_i4 ( word, l_val, ierror, length )

      if ( ierror /= 0 ) then
        exit
      end if

      if ( l_val /= 0 ) then
        l_num = l_num + 1
      end if

    end do

  end do
!
!  Close file and return.
!
  close ( unit = cnf_file_unit )

  return
end
subroutine cnf_header_write ( v_num, c_num, output_name, output_unit )

!*****************************************************************************80
!
!! CNF_HEADER_WRITE writes the header for a CNF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) V_NUM, the number of variables.
!
!    Input, integer ( kind = 4 ) C_NUM, the number of clauses.
!
!    Input, character ( len = * ) OUTPUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) OUTPUT_UNIT, the output unit.
!
  implicit none

  integer ( kind = 4 ) c_num
  character ( len = *  ) output_name
  integer ( kind = 4 ) output_unit
  character ( len = 80 ) string
  integer ( kind = 4 ) v_num

  write ( output_unit, '(a)' ) 'c ' // trim ( output_name )
  write ( output_unit, '(a)' ) 'c'
  write ( string, '(a,1x,i7,1x,i7)' ) 'p cnf', v_num, c_num
  call s_blanks_delete ( string )
  write ( output_unit, '(a)' ) string(1:len_trim(string))

  return
end
subroutine cnf_print ( v_num, c_num, l_num, l_c_num, l_val )

!*****************************************************************************80
!
!! CNF_PRINT prints CNF information.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) V_NUM, the number of variables.
!
!    Input, integer ( kind = 4 ) C_NUM, the number of clauses.
!
!    Input, integer ( kind = 4 ) L_NUM, the total number of signed literals.
!
!    Input, integer ( kind = 4 ) L_C_NUM(C_NUM), the number of signed
!    literals occuring in each clause.
!
!    Input, integer ( kind = 4 ) L_VAL(L_NUM), a list of all the signed
!    literals in all the clauses, ordered by clause.
!
  implicit none

  integer ( kind = 4 ) c_num
  integer ( kind = 4 ) l_num

  integer ( kind = 4 ) c
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l_c
  integer ( kind = 4 ) l_c_num(c_num)
  integer ( kind = 4 ) l_val(l_num)
  integer ( kind = 4 ) v_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CNF data printout:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of variables       V_NUM  = ', v_num
  write ( *, '(a,i8)' ) '  The number of clauses         C_NUM  = ', c_num
  write ( *, '(a,i8)' ) '  The number of signed literals L_NUM  = ', l_num
  l = 0
  do c = 1, c_num
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a,i8,a)' ) &
      '  Clause ', c, ' includes ', l_c_num(c), ' signed literals'
    do l_c = 1, l_c_num(c)
      l = l + 1
      write ( *, '(i4)' ) l_val(l)
    end do
  end do

  return
end
subroutine cnf_write ( v_num, c_num, l_num, l_c_num, l_val, output_name )

!*****************************************************************************80
!
!! CNF_WRITE writes the header and data of a CNF file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) V_NUM, the number of variables.
!
!    Input, integer ( kind = 4 ) C_NUM, the number of clauses.
!
!    Input, integer ( kind = 4 ) L_NUM, the total number of signed literals.
!
!    Input, integer ( kind = 4 ) L_C_NUM(C_NUM), the number of signed
!    literals occuring in each clause.
!
!    Input, integer ( kind = 4 ) L_VAL(L_NUM), a list of all the signed
!    literals in all the clauses, ordered by clause.
!
!    Input, character ( len = * ) OUTPUT_NAME, the name of the output file.
!
  implicit none

  integer ( kind = 4 ) c_num
  integer ( kind = 4 ) l_num

  integer ( kind = 4 ) l_c_num(c_num)
  integer ( kind = 4 ) l_val(l_num)
  character ( len = *  ) output_name
  integer ( kind = 4 ) output_unit
  integer ( kind = 4 ) v_num

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_name, status = 'replace' )

  call cnf_header_write ( v_num, c_num, output_name, output_unit )

  call cnf_data_write ( c_num, l_num, l_c_num, l_val, output_unit )

  close ( unit = output_unit )

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
  logical ( kind = 4 ) lopen

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
subroutine lvec_next ( n, lvec )

!*****************************************************************************80
!
!! LVEC_NEXT generates the next logical vector.
!
!  Discussion:
!
!    Let "0" represent FALSE and "1" represent TRUE.
!    Then the vectors have the order
!    (0,0,...,0),
!    (0,0,...,1),
!    ...
!    (1,1,...,1)
!
!    and the "next" vector after (1,1,...,1) is (0,0,...,0).  That is,
!    we allow wrap around.
!
!  Example:
!
!    N = 3
!
!    Input      Output
!    -----      ------
!    0 0 0  =>  0 0 1
!    0 0 1  =>  0 1 0
!    0 1 0  =>  0 1 1
!    0 1 1  =>  1 0 0
!    1 0 0  =>  1 0 1
!    1 0 1  =>  1 1 0
!    1 1 0  =>  1 1 1
!    1 1 1  =>  0 0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input/output, logical LVEC(N), on output, the successor to the
!    input vector.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical              lvec(n)

  do i = n, 1, -1

    if ( .not. lvec(i) ) then
      lvec(i) = .true.
      return
    end if
    lvec(i) = .false.
  end do

  return
end
subroutine s_blanks_delete ( s )

!*****************************************************************************80
!
!! S_BLANKS_DELETE replaces consecutive blanks by one blank.
!
!  Discussion:
!
!    Thanks to Bill Richmond for pointing out a programming flaw which
!    meant that, as characters were slid to the left through multiple
!    blanks, their original images were not blanked out.  This problem
!    is easiest resolved by using a copy of the string.
!
!    The remaining characters are left justified and right padded with blanks.
!    TAB characters are converted to spaces.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 September 2004
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

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character              newchr
  character              oldchr
  character ( len = *  ) s
  character ( len = len ( s ) ) s_copy
  integer ( kind = 4 ) s_length
  character, parameter :: TAB = achar ( 9 )

  s_length = len ( s )

  j = 0
  s_copy(1:s_length) = s(1:s_length)
  s(1:s_length) = ' '

  newchr = ' '

  do i = 1, s_length

    oldchr = newchr
    newchr = s_copy(i:i)

    if ( newchr == TAB ) then
      newchr = ' '
    end if

    if ( oldchr /= ' ' .or. newchr /= ' ' ) then
      j = j + 1
      s(j:j) = newchr
    end if

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
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lenc
  logical                s_eqi
  character ( len = *  ) s1
  integer ( kind = 4 ) s1_length
  character ( len = *  ) s2
  integer ( kind = 4 ) s2_length

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
!    13 January 2006
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

  integer ( kind = 4 ) get1
  integer ( kind = 4 ) get2
  character ( len = * ) s
  integer ( kind = 4 ) s_length
  character ( len = * ) w

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

  character ( len = 8  ) ampm
  integer ( kind = 4 ) d
  character ( len = 8  ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 ) time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5  ) zone

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
