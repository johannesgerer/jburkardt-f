program main

!*****************************************************************************80
!
!! MAIN is the main program for MATMAN.
!
!  Discussion:
!
!    MATMAN is a program for interactive linear algebra demonstrations.
!
!    This is version 1.63.
!
!  Journal:
!
!    09 October 2000:   
!    CHANGE calls VALUE_READ so I can add a column easily too.
!    COL_APP allows me to append a column.
!
!    29 September 2000: 
!    Adding IDENTITY operation which appends I.
!    Added modular arithmetic option.
!    Added RANDOM.
!
!    26 September 2000: 
!    Might almost be working.  Now either implement
!    modular arithmetic, or saving of pre and post multipliers.
!
!    25 September 2000: 
!    Add integer arithmetic option, decide what to do about division later.
!
!    24 September 2000: 
!    Throw out superfluous items!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 2000
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: maxcol = 30
  integer, parameter :: maxrow = 16

  real ( kind = 4 ) a(maxrow,maxcol)
  integer ( kind = 4 ) a_int(maxrow,maxcol)
  integer ( kind = 4 ) base
  logical ch_is_digit
  character ( len = 20 ) command
  integer ( kind = 4 ) i_temp
  integer ( kind = 4 ) iabot(maxrow,maxcol)
  integer ( kind = 4 ) iatop(maxrow,maxcol)
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) icol1
  integer ( kind = 4 ) icol2
  integer ( kind = 4 ) ictop(maxrow,maxcol)
  integer ( kind = 4 ) idbot
  integer ( kind = 4 ) idtop
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) iprint
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) irow1
  integer ( kind = 4 ) irow2
  character isay
  integer ( kind = 4 ) isbot
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jform
  integer ( kind = 4 ) jlo
  character ( len = 80 ) line
  character ( len = 80 ) line2
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  character ( len = 80 ) prompt
  logical s_eqi
  integer ( kind = 4 ) seed
  real ( kind = 4 ) sval

  call timestamp ( )
!
!  Say hello
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATMAN'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Version number 1.63'
  write ( *, '(a)' ) '  Last modified on 13 January 2002.'
  write ( *, '(a)' ) '  An interactive program to perform elementary'
  write ( *, '(a)' ) '  row and column operations on a matrix.'
  write ( *, '(a)' ) '  Developed by Charles Cullen and John Burkardt.'
!
!  Initializations.
!
  call init ( a, a_int, base, iabot, iatop, ierror, iform, &
    line, maxrow, maxcol, nrow, ncol )
!
!  Get the next command from the user.
!
  do
 
    iprint = 0
 
    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MATMAN could not carry out your command:'
      write ( *, '(a)' ) '"' // trim ( command ) // '"'
      ierror = 0
    end if

    write ( *, '(a)' ) ' '
    line = ' '
    prompt = 'command ( H for help )'
!
!  No check for terminators.
!
    iterm = 0
    call s_read ( line2, line, prompt, ierror, iterm )

    line = line2
 
    if ( line2 == ' ' ) then
      cycle
    end if
 
    if ( ierror /= 0 ) then
      ierror = 0
      command = 'Q'
    end if
!
!  Check to see if the command is an ERO, in a special format.
!
    call s_blank_delete ( line2 )

    if ( s_eqi ( line2, 'R' ) ) then

    else if ( s_eqi ( line2(1:1), 'R' ) .and. &
      ( line2(2:2) == ' ' .or. ch_is_digit ( line2(2:2) ) ) ) then

      call row_op_check ( command, ierror, line2 )

      if ( ierror /= 0 ) then
        cycle
      end if
 
      line = line2

    else if ( s_eqi ( line2(1:3), 'ROW' ) .and. &
      ( line2(4:4) == ' ' .or. ch_is_digit ( line2(4:4) ) ) ) then

      call row_op_check ( command, ierror, line2 )
      if ( ierror /= 0 ) then
        cycle
      end if

      line = line2
!
!  Check to see if the command is an ECO, in a special format.
!
    else if ( s_eqi ( line2, 'C' ) ) then

    else if ( s_eqi ( line2(1:1), 'C' ) .and. &
      ( line2(2:2) == ' ' .or. ch_is_digit ( line2(2:2) ) ) ) then

      call col_op_check ( command, ierror, line2 )

      if ( ierror /= 0 ) then
        cycle
      end if

      line = line2

    else if ( s_eqi ( line2(1:3), 'COL' ) .and. &
      ( line2(4:4) == ' ' .or. ch_is_digit ( line2(4:4) ) ) ) then

      call col_op_check ( command, ierror, line2 )
      if ( ierror /= 0 ) then
        cycle
      end if
      line = line2

    else if ( s_eqi ( line2(1:6), 'COLUMN' ) .and. &
      ( line2(7:7) == ' ' .or. ch_is_digit ( line2(7:7) ) ) ) then

      call col_op_check ( command, ierror, line2 )
      if ( ierror /= 0 ) then
        cycle
      end if
      line = line2
!
!  Maybe this is all I need in order to get my change command.
!
    else if ( s_eqi(line2(1:2), 'A(' ) ) then

      command = 'A('
      line = line2(3:)

    else

      command = ' '

    end if
!
!  If the command was not an ERO that had to be translated, 
!  read it the regular way.
!
!  Blank, slash, comma, semicolon, equals terminate the command.
!  The rest of the user input stays in LINE.
!
    if ( command == ' ' ) then

      iterm = 1
      call s_read ( command, line, prompt, ierror, iterm )

      if ( ierror /= 0 ) then
        ierror = 0
        command = 'Q'
      end if

    end if
!
!  A(I,J)=X
!
    if ( s_eqi ( command(1:2), 'A(' ) ) then
 
      call change ( base, iform, nrow, ncol, a, a_int, iatop, iabot, ierror )
 
      iprint = 1
!
!  B=Set up sample problem.
!
    else if ( s_eqi ( command, 'B' ) ) then

      call size_set ( ierror, line, maxrow, nrow, 'rows' )
      call size_set ( ierror, line, maxcol, ncol, 'columns' )

      if ( s_eqi ( isay,'S' ) ) then

        if ( maxcol < ncol + 1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Your problem is too big!'
          cycle
        end if
        ncol = ncol + 1

        call problem_solve_random ( a, a_int, iatop, iabot, base, iform, &
          nrow, ncol, seed )

      end if
!
!  BASE = set base for modular arithmetic.
!
    else if ( s_eqi ( command, 'BASE' ) ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Enter a positive integer to use as the base'
      write ( *, '(a)' ) 'for modular arithmetic.'

      read ( *, * ) base

      jform = 4

      call form ( a, a_int, iatop, iabot, base, iform, jform, nrow, ncol )

      iprint = 1 
!
!  CHECK = Echelon Form.
!
    else if ( s_eqi ( command, 'CHECK' ) ) then
 
      call la_opt ( a, a_int, iabot, iatop, base, iform, nrow, ncol )
!
!  COL_ADD=Add a multiple of one column to another.
!
    else if ( s_eqi ( command, 'COL_ADD' ) ) then
 
      call col_add_param ( ierror, base, iform, icol1, icol2, istop, isbot, &
        line, ncol, sval )
 
      if ( ierror /= 0 ) then
        cycle
      end if

      call col_add ( a, a_int, iatop, iabot, ierror, base, iform, icol1, icol2, &
        nrow, ncol, sval, istop, isbot )
 
      iprint = 1
!
!  COL_APP = Append column
!
    else if ( s_eqi ( command(1:7), 'COL_APP' ) ) then

      ncol = ncol + 1

      call col_app ( a, a_int, iatop, iabot, base, ierror, iform, line, &
        nrow, ncol )

      if ( ierror /= 0 ) then
        ncol = ncol - 1
      else
        iprint = 1
      end if
!
!  COL_AUTO=Automatic column reduction.
!
    else if ( s_eqi ( command, 'COL_AUTO' ) ) then
 
      call col_auto ( a, a_int, iatop, iabot, ierror, base, iform, nrow, ncol )
 
      iprint = 1
!
!  COL_DIV = Divide column by scalar.
!
    else if ( s_eqi ( command, 'COL_DIV' ) ) then
 
      call col_div_param ( icol, ierror, base, iform, isbot, istop, line, sval )

      call col_div ( a, a_int, iatop, iabot, ierror, base, iform, icol, &
        nrow, ncol, sval, istop, isbot )

      if ( ierror == 0 ) then
        iprint = 1
      end if
!
!  COL_MUL=Multiply column by scalar.
!
    else if ( s_eqi ( command, 'COL_MUL' ) ) then
 
      call col_mul_param ( ierror, base, iform, icol, istop, isbot, line, sval )
 
      if ( ierror /= 0 ) then
        cycle
      end if

      call col_mul ( a, a_int, iatop, iabot, ierror, base, iform, icol, &
        nrow, ncol, sval, istop, isbot )
 
      iprint = 1
!
!  COL_SWAP = Swap columns I and J.
!
    else if ( s_eqi ( command, 'COL_SWAP' ) ) then
 
      prompt = 'column I, column J.'
      call i4_read ( icol1, line, prompt, ierror )
      if ( ierror /= 0 ) then
        cycle
      end if
 
      call i4_read ( icol2, line, prompt, ierror )
      if ( ierror /= 0 ) then
        cycle
      end if
 
      call col_swap ( a, a_int, iatop, iabot, ierror, iform, icol1, icol2, &
        nrow, ncol )
 
      iprint = 1
!
!  DEC_DIGIT=Set number of digits.
!
    else if ( s_eqi ( command, 'DEC_DIGIT' ) ) then
 
      call dec_digit_set ( ierror, line )
!
!  DECimal = use decimal arithmetic
!
    else if ( s_eqi ( command(1:3), 'DEC' ) ) then
 
      jform = 2
 
      call form ( a, a_int, iatop, iabot, base, iform, jform, nrow, ncol )

      iprint = 1
!
!  E=Enter problem definition
!
    else if ( s_eqi ( command, 'E' ) ) then

      call size_set ( ierror, line, maxrow, nrow, 'rows' )
      call size_set ( ierror, line, maxcol, ncol, 'columns' )
 
      call la_inp1 ( nrow, ncol, a, a_int, iabot, iatop, ierror, base, iform, &
        line, 1, nrow, 1, ncol )
   
      if ( ierror /= 0 ) then
        cycle
      end if
 
      iprint = 1
!
!  H=Help.
!
    else if ( s_eqi ( command(1:1), 'H' ) ) then
 
      call help
!
!  I4_BIG = Set maximum integer.
!
    else if ( s_eqi ( command, 'I4_BIG' ) ) then
 
      prompt = 'maximum integer for rational representations.'
 
      call i4_read ( i_temp, line, prompt, ierror )

      call i4_data ( 'SET', 'I4_BIG', i_temp )
!
!  IDENTITY
!
    else if ( s_eqi ( command(1:2), 'ID' ) ) then

      if ( nrow <= 0 ) then
        call size_set ( ierror, line, maxrow, nrow, 'rows' )
        ncol = 0
        if ( ierror /= 0 ) then
          cycle
        end if
      end if

      if ( maxcol < ncol + nrow ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Not enough room to append the identity matrix.'
        cycle
      end if

      call mat_append_identity ( nrow, ncol, a, a_int, iatop, iabot, iform )

      iprint = 1
!
!  INIT = Initializations.
!
    else if ( s_eqi ( command(1:3), 'INIT' ) ) then

      call init ( a, a_int, base, iabot, iatop, ierror, iform, &
        line, maxrow, maxcol, nrow, ncol )
!
!  INT = Integer arithmetic
!
    else if ( s_eqi ( command(1:3), 'INT' ) ) then

      jform = 3

      call form ( a, a_int, iatop, iabot, base, iform, jform, nrow, ncol )

      iprint = 1
!
!  Q = Quit (after confirmation).
!  QUit = QUIT NOW!
!
    else if ( s_eqi ( command(1:1), 'Q' ) ) then
 
      if ( s_eqi ( command(2:2), 'U' ) ) then
        exit
      end if

      line = ' '
      prompt = '"Y" to confirm you want to quit.'
      iterm = 0

      call s_read ( isay, line, prompt, ierror, iterm )

      if ( ierror /= 0 ) then
        exit
      end if
 
      if ( s_eqi ( isay, 'Y' ) ) then
        exit
      end if
!
!  RANDOM = Enter problem definition
!
    else if ( s_eqi ( command(1:3), 'RAN' ) ) then

      call size_set ( ierror, line, maxrow, nrow, 'rows' )
      call size_set ( ierror, line, maxcol, ncol, 'columns' )

      call mat_random ( nrow, ncol, base, iform, seed, a, a_int, iatop, iabot )
 
      iprint = 1
!
!  RATional = use rational arithmetic
!
    else if ( s_eqi ( command(1:3), 'RAT' ) ) then
 
      jform = 0
 
      call form ( a, a_int, iatop, iabot, base, iform, jform, nrow, ncol )
 
      iprint = 1
!
!  REAl = use real arithmetic
!
    else if ( s_eqi ( command(1:3), 'REA' ) ) then
 
      jform = 1
 
      call form ( a, a_int, iatop, iabot, base, iform, jform, nrow, ncol )
 
      iprint = 1
!
!  ROW_ADD=Add a multiple of one row to another.
!
    else if ( s_eqi ( command, 'ROW_ADD' ) ) then
 
      call row_add_param ( ierror, base, iform, irow1, irow2, istop, isbot, &
        line, nrow, sval )
 
      if ( ierror /= 0 ) then
        cycle
      end if

      call row_add ( a, a_int, iatop, iabot, ierror, base, iform, irow1, irow2, &
        nrow, ncol, sval, istop, isbot )
 
      iprint = 1
!
!  ROW_AUTO=Automatic row reduction.
!
    else if ( s_eqi ( command, 'ROW_AUTO' ) ) then
 
      call row_auto ( a, a_int, iatop, iabot, ierror, base, iform, nrow, ncol )
 
      iprint = 1
!
!  ROW_DIV = Divide row by scalar.
!
    else if ( s_eqi ( command, 'ROW_DIV' ) ) then
 
      call row_div_param ( ierror, base, iform, irow, isbot, istop, line, sval )

      call row_div ( a, a_int, iatop, iabot, ierror, base, iform, irow, &
        nrow, ncol, sval, istop, isbot )
 
      if ( ierror == 0 ) then
        iprint = 1
      end if
!
!  ROW_MUL=Multiply row by scalar.
!
    else if ( s_eqi ( command, 'ROW_MUL' ) ) then

      call row_mul_param ( ierror, base, iform, irow, istop, isbot, line, sval )
 
      if ( ierror /= 0 ) then
        cycle
      end if

      call row_mul ( a, a_int, iatop, iabot, ierror, base, iform, irow, &
        nrow, ncol, sval, istop, isbot )
 
      iprint = 1
!
!  ROW_SWAP = Interchange rows I and J.
!
    else if ( s_eqi ( command, 'ROW_SWAP' ) ) then
 
      prompt = 'row I, row J.'
      call i4_read ( irow1, line, prompt, ierror )
      if ( ierror /= 0 ) then
        cycle
      end if
 
      call i4_read ( irow2, line, prompt, ierror )
      if ( ierror /= 0 ) then
        cycle
      end if
 
      call row_swap ( a, a_int, iatop, iabot, ierror, iform, &
        irow1, irow2, nrow, ncol )
 
      iprint = 1
!
!  #: comment, ignore
!
    else if ( command(1:1) == '#' ) then
!
!  Blank: type out the current matrix
!
    else if ( len_trim ( command ) == 0 ) then

      iprint = 1
!
!  Unrecognized command.
!
    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Ignoring the command "' // trim ( command ) // '".'
 
    end if
!
!  After certain operations, print out the matrix.
!
    if ( ierror == 0 .and. iprint == 1 ) then

      call mat_print ( a, a_int, iabot, iatop, iform, nrow, ncol, &
        'The current matrix:' )
 
      iprint = 0
 
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATMAN'
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
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 August 1999
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

  logical ch_eqi
  character c1
  character c2
  character cc1
  character cc2

  cc1 = c1
  cc2 = c2

  call ch_cap ( cc1 )
  call ch_cap ( cc2 )

  if ( cc1 == cc2 ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
function ch_is_alpha ( c )

!*****************************************************************************80
!
!! CH_IS_ALPHA returns TRUE if C is an alphabetic character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, a character to check.
!
!    Output, logical CH_IS_ALPHA is TRUE if C is an alphabetic character.
!
  implicit none

  character c
  logical ch_is_alpha

  if ( ( lle ( 'a', c ) .and. lle ( c, 'z' ) ) .or. &
       ( lle ( 'A', c ) .and. lle ( c, 'Z' ) ) ) then
    ch_is_alpha = .true.
  else
    ch_is_alpha = .false.
  end if

  return
end
function ch_is_digit ( c )

!*****************************************************************************80
!
!! CH_IS_DIGIT returns .TRUE. if a character is a decimal digit.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
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
subroutine change ( base, iform, nrow, ncol, a, a_int, iatop, iabot, ierror )

!*****************************************************************************80
!
!! CHANGE allows the user to change an entry in the array.
!
!  Discussion:
!
!    The function expects to process a user command of the form:
!
!      A(5,12) = 17
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 November 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
!    Input/output, real ( kind = 4 ) A(NROW,NCOL), the current real matrix.
!
!    Input/output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input/output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL),
!    the current rational or decimal matrix.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 for an error.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  integer ( kind = 4 ) base
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 22 ) chrtmp3
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) isbot
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) ival
  character ( len = 80 ) line
  character ( len = 80 ) prompt
  real ( kind = 4 ) rval

  prompt = 'row I, column J, new value S.'
!
!  Get the row number.
!
  call i4_read ( irow, line, prompt, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHANGE - Error!'
    write ( *, '(a)' ) '  I4_READ returned error flag!'   
    return
  end if
 
  if ( irow < 1 .or. nrow < irow ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHANGE - Error!'
    write ( *, '(a)' ) '  Illegal row value!'
    ierror = 1
    return
  end if
!
!  Get the column number.
!
  if ( line(1:1) == ',' ) then
    line(1:1) = ' '
  end if

  call i4_read ( icol, line, prompt, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHANGE - Error!'
    write ( *, '(a)' ) '  I4_READ returned error flag!'
    return
  end if
 
  if ( icol < 1 .or. ncol < icol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHANGE - Error!'
    write ( *, '(a)' ) '  Illegal column value!'    
    ierror = 1
    return
  end if
!
!  Read the value.
!
  if ( line(1:2) == ')=' ) then
    line(1:2) = ' '
  end if

  call value_read ( irow, icol, base, iform, nrow, ncol, line, a, a_int, &
    iatop, iabot, ierror )

  call i4_to_s_left ( irow, chrtmp1 )
  call i4_to_s_left ( icol, chrtmp2 )

  if ( iform == 0 ) then
    call rat_to_s_left ( iatop(irow,icol), iabot(irow,icol), chrtmp3 )
  else if ( iform == 1 ) then
    call r4_to_s_left ( a(irow,icol), chrtmp3 )
  else if ( iform == 2 ) then
    call dec_to_s_left ( iatop(irow,icol), iabot(irow,icol), chrtmp3 )
  else if ( iform == 3 ) then
    call i4_to_s_left ( a_int(irow,icol), chrtmp3 )
  else if ( iform == 4 ) then
    call i4_to_s_left ( a_int(irow,icol), chrtmp3 )
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A(' // trim ( chrtmp1 ) // ',' // trim ( chrtmp2 ) // &
    ') = ' // trim ( chrtmp3 )
  
  return
end
subroutine chrctf ( s, itop, ibot, ierror, lchar )

!*****************************************************************************80
!
!! CHRCTF reads an integer or rational fraction from a string.
!
!  Discussion:
!
!    The integer may be in real format, for example '2.25'.  The routine
!    returns ITOP and IBOT.  If the input number is an integer, ITOP
!    equals that integer, and IBOT is 1.  But in the case of 2.25,
!    the program would return ITOP = 225, IBOT = 100.
!
!    Legal input is:
!
!      blanks,
!      initial sign,
!      blanks,
!      integer ( kind = 4 ) part,
!      decimal point,
!      fraction part,
!      'E' or 'e' or 'D' or 'd', exponent marker,
!      exponent sign,
!      exponent integer part,
!      blanks,
!      final comma or semicolon.
!
!    with most quantities optional.
!
!  Example:
!
!    S               ITOP      IBOT
!
!    '1'               1         1
!    '     1   '       1         1
!    '1A'              1         1
!    '12,34,56'        12        1
!    '  34 7'          34        1
!    '-1E2ABCD'        -100      1
!    '-1X2ABCD'        -1        1
!    ' 2E-1'           2         10
!    '23.45'           2345      100
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate when no more characters
!    can be read to form a legal integer.  Blanks, commas,
!    or other nonnumeric data will, in particular, cause
!    the conversion to halt.
!
!    Output, integer ( kind = 4 ) ITOP, the integer read from the string,
!    assuming that no negative exponents or fractional parts
!    were used.  Otherwise, the 'integer' is ITOP/IBOT.
!
!    Output, integer ( kind = 4 ) IBOT, the integer divisor required to
!    represent numbers which are in real format or have a
!    negative exponent.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0 if no errors,
!    Value of IHAVE when error occurred otherwise.
!
!    Output, integer ( kind = 4 ) LCHAR, number of characters read from
!    the string to form the number.
!
  implicit none

  logical ch_eqi
  character c
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  character ( len = * ) s

  nchar = len_trim ( s )

  ierror = 0
  lchar = - 1
  isgn = 1
  itop = 0
  ibot = 1
  jsgn = 1
  jtop = 0
  ihave = 1
  iterm = 0

  do while ( lchar < nchar )

    lchar = lchar + 1
    c = s(lchar+1:lchar+1)
!
!  Blank.
!
    if ( c == ' ' ) then

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
        lchar = lchar + 1
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
    else if ( lge ( c, '0' ) .and. lle ( c, '9' ) .and. ihave < 11 ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        itop = 10 * itop + ndig
      else if ( ihave == 5 ) then
        itop = 10 * itop + ndig
        ibot = 10 * ibot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if

    if ( iterm == 1 ) then
      exit
    end if

  end do

  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHRCTF - Serious error!'
    write ( *, '(a,a)' ) '  Illegal input:', trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jsgn == 1 ) then
    itop = itop * 10**jtop
  else
    ibot = ibot * 10**jtop
  end if

  itop = isgn * itop

  return
end
subroutine chrctg ( string, itop, ibot, ierror, lchar )

!*****************************************************************************80
!
!! CHRCTG reads an integer, decimal fraction or a ratio from a string.
!
!  Discussion:
!
!    CHRCTG returns an equivalent ratio (ITOP/IBOT).
!
!    If the input number is an integer, ITOP equals that integer, and
!    IBOT is 1.   But in the case of 2.25, the program would return
!    ITOP = 225, IBOT = 100.
!
!    A ratio is either
!      a number
!    or
!      a number, "/", a number.
!
!    A "number" is defined as:
!
!      blanks,
!      initial sign,
!      integer ( kind = 4 ) part,
!      decimal point,
!      fraction part,
!      E,
!      exponent sign,
!      exponent integer part,
!      blanks,
!      final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    STRING            ITOP      IBOT
!
!    '1'               1         1
!    '     1   '       1         1
!    '1A'              1         1
!    '12,34,56'        12        1
!    '  34 7'          34        1
!    '-1E2ABCD'        -100      1
!    '-1X2ABCD'        -1        1
!    ' 2E-1'           2         10
!    '23.45'           2345      100
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate when no more characters
!    can be read to form a legal integer.  Blanks, commas,
!    or other nonnumeric data will, in particular, cause
!    the conversion to halt.
!
!    Output, integer ( kind = 4 ) ITOP, the integer read from the string,
!    assuming that no negative exponents or fractional parts
!    were used.  Otherwise, the 'integer' is ITOP/IBOT.
!
!    Output, integer ( kind = 4 ) IBOT, the integer divisor required to
!    represent numbers which are in decimal format or have a
!    negative exponent.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 for an error.
!
!    Output, integer ( kind = 4 ) LCHAR, the number of characters read from
!    STRING to form the number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ibotb
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) itopb
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) lchar2
  integer ( kind = 4 ) nchar
  character ( len = * ) string

  itop = 0
  ibot = 1
  lchar = 0
 
  call chrctf ( string, itop, ibot, ierror, lchar )

  if ( ierror /= 0 ) then
    return
  end if
!
!  The number is represented as a fraction.
!  If the next nonblank character is "/", then read another number.
!
  nchar = len_trim ( string )
 
  do i = lchar+1, nchar-1
 
    if ( string(i:i) == '/' ) then
 
      call chrctf ( string(i+1:), itopb, ibotb, ierror, lchar2 )

      if ( ierror /= 0 ) then
        return
      end if
 
      itop = itop * ibotb
      ibot = ibot * itopb
 
      itemp = i4_gcd ( itop, ibot )
 
      itop = itop / itemp
      ibot = ibot / itemp
 
      lchar = i + lchar2
 
      return
 
    else if ( string(i:i) /= ' ' ) then
 
      return
 
    end if
 
  end do
 
  return
end
subroutine chrinp ( ierror, line, prompt )

!*****************************************************************************80
!
!! CHRINP requests new input if the LINE buffer is empty.
!
!  Discussion:
!
!    CHRINP checks to see whether there is any more information in
!    the buffer array LINE.  If so, it simply updates the prompt
!    and returns.  Otherwise, it prints the prompt string out,
!    reads the input from the user, and reprints the prompt and
!    the user input on those I/O units where it is appropriate.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 for an error.
!
!    Input/output, character ( len = 80 ) LINE.
!
!    On input, LINE may contain information that the calling
!    program can use, or LINE may be empty.
!
!    On output, LINE is unchanged if it contained information
!    on input.  But if the input LINE was empty, then the
!    output LINE contains whatever information the user typed.
!
!    Input/output, character ( len = 80 ) PROMPT.
!    On input, the prompt string to be printed.
!    On output, PROMPT has been blanked out, up to the first comma.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) icomma
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) lchar
  character ( len = 80 ) line
  character ( len = 80 ) prompt

  ierror = 0
 
  if ( line == ' ' ) then
 
    write ( *, '(a)' ) 'Enter ' // trim ( prompt )
 
    read ( *, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CHRINP - Warning!'
      write ( *, '(a)' ) '  End of input!'
      stop
    end if

    call s_blanks_delete ( line )

  end if
!
!  If item was read, remove item from PROMPT list.
!
  if ( line /= ' ' ) then

    icomma = index ( prompt, ',' )

    if ( 0 < icomma .and. icomma < 80 .and. &
         prompt(icomma+1:icomma+1) == ' ' ) then

      icomma = icomma + 1

    end if

    call s_chop ( prompt, 1, icomma )

  end if
 
  return
end
subroutine col_add ( a, a_int, iatop, iabot, ierror, base, iform, icol1, &
  icol2, nrow, ncol, sval, istop, isbot )

!*****************************************************************************80
!
!! COL_ADD adds a multiple of one column to another.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(NROW,NCOL), the current matrix.
!
!    Input/output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input/output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Output, integer ( kind = 4 ) IERROR, 0 = no error, 1 = an error.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) ICOL1, the column which is to be modified.
!
!    Input, integer ( kind = 4 ) ICOL2, the column which is to be multiplied by
!    a given value and added to row ICOL1.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
!    Input, real ( kind = 4 ) SVAL, the multiplier to use if real arithmetic is employed.
!
!    Input, integer ( kind = 4 ) ISTOP, ISBOT, the fractional or decimal
!    multiplier to use.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  integer ( kind = 4 ) base
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 22 ) chrtmp3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) icol1
  integer ( kind = 4 ) icol2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) isbot
  integer ( kind = 4 ) isbot2
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) istop2
  integer ( kind = 4 ) itop
  real ( kind = 4 ) sval
!
!  Return immediately if the multiplier is zero.
!
  if ( ( iform == 0 .and. istop == 0 ) .or. &
       ( iform == 1 .and. sval == 0.0E+00 ) .or. &
       ( iform == 2 .and. istop == 0 ) .or. &
       ( iform == 3 .and. istop == 0 ) .or. &
       ( iform == 4 .and. istop == 0 ) ) then
      return
    end if
!
!  Carry out the operation.
!
  if ( iform == 0 ) then
 
    do i = 1, nrow
 
      call rat_mul ( isbot2, isbot, iabot(i,icol2), istop2, istop, &
        iatop(i,icol2), ierror )
 
      if ( ierror /= 0 ) then
        return
      end if

      call rat_add ( ibot, iabot(i,icol1), isbot2, itop, iatop(i,icol1), &
        istop2 )
 
      iatop(i,icol1) = itop
      iabot(i,icol1) = ibot
 
    end do
 
  else if ( iform == 1 ) then
 
    a(1:nrow,icol1) = a(1:nrow,icol1) + sval * a(1:nrow,icol2)

  else if ( iform == 2 ) then
 
    do i = 1, nrow
 
      call dec_mul ( isbot2, isbot, iabot(i,icol2), istop2, istop, &
        iatop(i,icol2) )
 
      call dec_add ( ibot, iabot(i,icol1), isbot2, itop, iatop(i,icol1), &
        istop2 )
 
      iatop(i,icol1) = itop
      iabot(i,icol1) = ibot
 
    end do

  else if ( iform == 3 ) then
 
    a_int(1:nrow,icol1) = a_int(1:nrow,icol1) + istop * a_int(1:nrow,icol2)

  else if ( iform == 4 ) then

    a_int(1:nrow,icol1) = &
      mod ( a_int(1:nrow,icol1) + istop * a_int(1:nrow,icol2), base )

  end if
!
!  Print out a message.
!
  if ( iform == 0 ) then
 
    if ( istop == isbot ) then
      chrtmp3 = '+'
    else if ( istop == - isbot ) then
      chrtmp3 = '-'
    else
      call rat_to_s_left ( istop, isbot, chrtmp3 )
    end if
 
  else if ( iform == 1 ) then
 
    if ( sval == 1.0E+00 ) then
      chrtmp3 = '+'
    else if ( sval == - 1.0E+00 ) then
      chrtmp3 = '-'
    else
      call r4_to_s_left ( sval, chrtmp3 )
    end if
 
  else if ( iform == 2 ) then

    if ( istop == 1 .and. isbot == 0 ) then 
      chrtmp3 = '+'
    else if ( istop == - 1 .and. isbot == 0 ) then
      chrtmp3 = '-'
    else
      call dec_to_s_left ( istop, isbot, chrtmp3 )
    end if
 
  else if ( iform == 3 ) then
 
    if ( istop == 1 ) then
      chrtmp3 = '+'
    else if ( istop == - 1 ) then
      chrtmp3 = '-'
    else
      call i4_to_s_left ( istop, chrtmp3 )
    end if

  else if ( iform == 4 ) then
 
    if ( istop == 1 ) then
      chrtmp3 = '+'
    else if ( istop == - 1 ) then
      chrtmp3 = '-'
    else
      call i4_to_s_left ( istop, chrtmp3 )
    end if

  end if
 
  call i4_to_s_left ( icol1, chrtmp1 )

  call i4_to_s_left ( icol2, chrtmp2 )

  if ( chrtmp3 == '-' .or. chrtmp3 == '+' ) then

    write ( *, '(a)' ) '  ECO: Col ' // trim ( chrtmp1 ) // ' <=  ' // 'Col ' // &
      trim ( chrtmp1 ) // ' ' // chrtmp3(1:1) // ' Col ' // trim ( chrtmp2 )

  else if ( chrtmp3(1:1) == '-' .or. chrtmp3(1:1) == '+' ) then

    write ( *, '(a)' ) '  ECO: Col ' // trim ( chrtmp1 ) // ' <=  ' // 'Col ' // &
      trim ( chrtmp1 ) // ' ' // chrtmp3(1:1) // ' ' // trim ( chrtmp3(2: ) ) // &
      ' Col ' // trim ( chrtmp2 )

  else

    write ( *, '(a)' ) '  ECO: Col ' // trim ( chrtmp1 ) // ' <=  ' // 'Col ' // &
      trim ( chrtmp1 ) // ' + ' // trim ( chrtmp3 ) // ' Col ' // trim ( chrtmp2 )

  end if
 
  return
end
subroutine col_add_param ( ierror, base, iform, icol1, icol2, istop, isbot, &
  line, ncol, sval )

!*****************************************************************************80
!
!! COL_ADD_PARAM gets and checks the column add parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IERROR, 0 = no error, 1 = an error.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) ICOL1, the column to which the multiple is to be added.
!
!    Input, integer ( kind = 4 ) ICOL2, the column which is to be multiplied and
!    added to another row.
!
!    Output, integer ( kind = 4 ) ISTOP, ISBOT, the parts of the rational
!    or decimal fraction of the multiplier, if that is the
!    arithmetic being used.
!
!    Workspace, character ( len = 80 ) LINE, used to hold the user's input.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns in the matrix.
!
!    Output, real ( kind = 4 ) SVAL, the multiplier, if real arithmetic is used.
!
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) icol1
  integer ( kind = 4 ) icol2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) isbot
  integer ( kind = 4 ) istop
  character ( len = 80 ) line
  integer ( kind = 4 ) ncol
  character ( len = 80 ) prompt
  real ( kind = 4 ) sval

  prompt = 'multiplier S, column I to add, target column J.'
!
!  Get the multiplier, SVAL or ISTOP/ISBOT.
!
  if ( iform == 0 ) then
    call rat_read ( istop, isbot, line, prompt, ierror )
  else if ( iform == 1 ) then
    call r4_read ( sval, line, prompt, ierror )
  else if ( iform == 2 ) then
    call dec_read ( istop, isbot, line, prompt, ierror )
    if ( ierror == 0 ) then
      call dec_round ( istop, isbot )
    end if
  else if ( iform == 3 ) then
    call i4_read ( istop, line, prompt, ierror )
  else if ( iform == 4 ) then
    call i4_read ( istop, line, prompt, ierror )
    istop = mod ( istop, base )
  end if

  if ( ierror /= 0 ) then
    return
  end if
!
!  Get the column to add.
!
  call i4_read ( icol2, line, prompt, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Error reading column index from line:'
    write ( *, '(a)' ) '"' // trim ( line ) // '"'
    return
  end if

  if ( icol2 < 1 .or. ncol < icol2 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Error!  Column index was not acceptable!'
    return
  end if
!
!  Get the column to which we are adding.
!
  call i4_read ( icol1, line, prompt, ierror )
  if ( ierror /= 0 ) return
 
  if ( icol1 < 1 .or. ncol < icol1 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Error!  The column index was not acceptable!'
    return
  end if
!
!  Make sure the columns are different.
!
  if ( icol1 == icol2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Error!  The columns should not be the same!'
    ierror = 1
    return
  end if
 
  ierror = 0
 
  return
end
subroutine col_app ( a, a_int, iatop, iabot, base, ierror, iform, line, &
  nrow, ncol )

!*****************************************************************************80
!
!! COL_APP allows the user to append a column to the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(NROW,NCOL), the real matrix.
!
!    Input/output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input/output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL), 
!    the rat/dec matrix.
!
!    Input, integer ( kind = 4 ) BASE, the base of modular arithmetic.
!
!    Output, integer ( kind = 4 ) IERROR, 0 = no error, 1 = an error.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Workspace, character ( len = 80 ) LINE, used to hold the user's input.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) irow
  character ( len = 80 ) line

  icol = ncol

  do irow = 1, nrow

    call value_read ( irow, icol, base, iform, nrow, ncol, line, a, a_int, &
      iatop, iabot, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COL_APP - Error!'
      write ( *, '(a)' ) '  VALUE_READ returned an error.'
      return
    end if

  end do

  return
end
subroutine col_auto ( a, a_int, iatop, iabot, ierror, base, iform, nrow, ncol )

!*****************************************************************************80
!
!! COL_AUTO automatically column reduces the current matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(NROW,NCOL), the current matrix.
!
!    Input/output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input/output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL).
!    IATOP and IABOT represent the rational or decimal matrix
!    to which elementary row operations will be applied.
!
!    Output, integer ( kind = 4 ) IERROR, 0 an error occurred, 1 no error occurred.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  real ( kind = 4 ) amax
  real ( kind = 4 ) atemp
  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) isbot
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) kcol
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lcol
  real ( kind = 4 ) sval

  ierror = 0
  kcol = 0
 
  do j = 1, ncol
 
    jcol = j
 
    do i = 1, nrow
 
      irow = i
!
!  In row IROW, seek the column between ICOL and NCOL with
!  maximum nonzero entry AMAX.
!
      imax = 0
      amax = 0.0E+00
 
      do kcol = jcol, ncol
 
        if ( iform == 0 ) then
          call rat_to_r4 ( iatop(irow,kcol), iabot(irow,kcol), atemp )
        else if ( iform == 1 ) then
          atemp = a(irow,kcol)
        else if ( iform == 2 ) then
          call dec_to_r4 ( iatop(irow,kcol), iabot(irow,kcol), atemp )
        else if ( iform == 3 ) then
          atemp = real ( a_int(irow,kcol) )
        else if ( iform == 4 ) then
          atemp = real ( a_int(irow,kcol) )
        end if
 
        atemp = abs ( atemp )
 
        if ( amax < atemp ) then
          amax = atemp
          imax = kcol
        end if
 
      end do
 
      if ( imax /= 0 ) then
        kcol = imax
        exit
      end if
 
    end do
 
    if ( kcol == 0 ) then
      return
    end if
!
!  Interchange the JCOL-th and the pivot columns.
!
    if ( kcol /= jcol ) then
      call col_swap ( a, a_int, iatop, iabot, ierror, iform, kcol, jcol, &
        nrow, ncol )
    end if
!
!  Divide the pivot column by A(IROW,JCOL) so that A(IROW,JCOL) = 1.
!
    if ( iform == 0 ) then
 
      istop = iatop(irow,jcol)
      isbot = iabot(irow,jcol)
 
    else if ( iform == 1 ) then
 
      sval = a(irow,jcol)
 
    else if ( iform == 2 ) then
 
      istop = iatop(irow,jcol)
      isbot = iabot(irow,jcol)

    else if ( iform == 3 ) then

      istop = a_int(irow,jcol)

    else if ( iform == 4 ) then

      istop = a_int(irow,jcol) 

    end if
 
    call col_div ( a, a_int, iatop, iabot, ierror, base, iform, jcol, &
      nrow, ncol, sval, istop, isbot )
!
!  Annihilate A(IROW,L) for L not equal to JCOL.
!
    do l = 1, nrow
 
      lcol = l
 
      if ( lcol /= jcol ) then
 
        if ( iform == 0 ) then

          if ( iatop(irow,lcol) /= 0 ) then

            istop = - iatop(irow,lcol)
            isbot = iabot(irow,lcol)

            call col_add ( a, a_int, iatop, iabot, ierror, base, iform, &
              lcol, jcol, nrow, ncol, sval, istop, isbot )

            iatop(irow,lcol) = 0
            iabot(irow,lcol) = 1

          end if

        else if ( iform == 1 ) then

          if ( a(irow,lcol) /= 0.0E+00 ) then

            sval = - a(irow,lcol)

            call col_add ( a, a_int, iatop, iabot, ierror, base, iform, &
              lcol, jcol, nrow, ncol, sval, istop, isbot )

            a(irow,lcol) = 0.0E+00

          end if

        else if ( iform == 2 ) then

          if ( iatop(irow,lcol) /= 0 ) then

            istop = - iatop(irow,lcol)
            isbot = iabot(irow,lcol)

            call col_add ( a, a_int, iatop, iabot, ierror, base, iform, &
              lcol, jcol, nrow, ncol, sval, istop, isbot )

            iatop(irow,lcol) = 0
            iabot(irow,lcol) = 0

          end if

        else if ( iform == 3 ) then

          if ( a_int(irow,lcol) /= 0 ) then

            istop = - a_int(irow,lcol)

            call col_add ( a, a_int, iatop, iabot, ierror, base, iform, &
              lcol, jcol, nrow, ncol, sval, istop, isbot )

            a_int(irow,lcol) = 0

          end if

        else if ( iform == 4 ) then

          if ( a_int(irow,lcol) /= 0 ) then

            istop = - a_int(irow,lcol)

            call col_add ( a, a_int, iatop, iabot, ierror, base, iform, &
              lcol, jcol, nrow, ncol, sval, istop, isbot )

            a_int(irow,lcol) = 0

          end if
        end if
 
      end if
 
    end do
 
  end do
 
  return
end
subroutine col_div ( a, a_int, iatop, iabot, ierror, base, iform, icol, &
  nrow, ncol, sval, istop, isbot )

!*****************************************************************************80
!
!! COL_DIV divides a column of the matrix by a scale factor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(NROW,NCOL), the current matrix.
!
!    Input/output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input/output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Output, integer ( kind = 4 ) IERROR, 0 an error occurred, 1 no error occurred.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) ICOL, the column to be divided.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
!    Input, real ( kind = 4 ) SVAL, the real divisor.
!
!    Input, integer ( kind = 4 ) ISTOP, ISBOT, the fractional or decimal divisor.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  integer ( kind = 4 ) base
  character ( len = 10 ) chrtmp1
  character ( len = 24 ) chrtmp2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) isbot
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) itop
  character ( len = 3 ) op
  real ( kind = 4 ) sval
!
!  Make sure that the column number is legal.
!
  if ( icol < 1 .or. ncol < icol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Error!  The column number is out of range!'
    ierror = 1
    return
  end if
!
!  Check for an illegal divisor of 0.
!
  if ( ( iform == 0 .and. istop == 0 ) .or. &
       ( iform == 1 .and. sval == 0.0E+00 ) .or. &
       ( iform == 2 .and. istop == 0 ) .or. &
       ( iform == 3 .and. istop == 0 ) .or. &
       ( iform == 4 .and. istop == 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COL_DIV - Error!'
    write ( *, '(a)' ) '  It is illegal to divide by 0!'  
    ierror = 1
    return
  end if
!
!  Check for a pointless division by 1.
!
  if ( ( iform == 0 .and. istop == isbot ) .or. &
       ( iform == 1 .and. sval == 1.0E+00 ) .or. &
       ( iform == 2 .and. istop == 1 .and. isbot == 0 ) .or. &
       ( iform == 3 .and. istop == 1 ) .or. &
       ( iform == 4 .and. istop == 1 ) ) then
    return
  end if
!
!  Carry out the division.
!
  if ( iform == 0 ) then
 
    do i = 1, nrow
 
      call rat_mul ( ibot, iabot(i,icol), istop, itop, iatop(i,icol), isbot, &
        ierror )
 
      if ( ierror /= 0 ) then
        return
      end if
 
      iatop(i,icol) = itop
      iabot(i,icol) = ibot
 
    end do
 
  else if ( iform == 1 ) then
    a(1:nrow,icol) = a(1:nrow,icol) / sval
  else if ( iform == 2 ) then
 
    do i = 1, nrow
 
      call dec_div ( iatop(i,icol), iabot(i,icol), istop, isbot, &
        iatop(i,icol), iabot(i,icol), ierror )
 
    end do
 
  else if ( iform == 3 ) then
    a_int(1:nrow,icol) = a_int(1:nrow,icol) / istop
  else if ( iform == 4 ) then
    a_int(1:nrow,icol) = a_int(1:nrow,icol) / istop
  end if
!
!  Print out a statement about what has been done.
!
  if ( iform == 0 ) then
 
    if ( isbot == 1 ) then

      call i4_to_s_left ( istop, chrtmp2 )
      op = ' / '

    else

      call rat_to_s_left ( isbot, istop, chrtmp2 )
      op = ' * '

    end if
 
  else if ( iform == 1 ) then
 
    call r4_to_s_left ( sval, chrtmp2 )
    op = ' / '

  else if ( iform == 2 ) then
 
    call dec_to_s_left ( istop, isbot, chrtmp2 )
    op = ' / '

  else if ( iform == 3 ) then
 
    call i4_to_s_left ( istop, chrtmp2 )
    op = ' / '

  else if ( iform == 4 ) then
 
    call i4_to_s_left ( istop, chrtmp2 )
    op = ' / '

  end if

  call i4_to_s_left ( icol, chrtmp1 )

  write ( *, '(a)' ) '  ECO: Col ' // trim ( chrtmp1 ) // ' <=  Col ' // &
    trim ( chrtmp1 ) // op // trim ( chrtmp2 )
 
  return
end
subroutine col_div_param ( icol, ierror, base, iform, isbot, istop, line, sval )

!*****************************************************************************80
!
!! COL_DIV_PARAM gets and checks the column divide parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IERROR, 0 an error occurred, 1 no error occurred.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Workspace, character ( len = 80 ) LINE, used to hold the user's input.
!
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) isbot
  integer ( kind = 4 ) istop
  character ( len = 80 ) line
  character ( len = 80 ) prompt
  real ( kind = 4 ) sval

  prompt = 'column I, divisor S.'
!
!  Read the column number to be divided.
!
  call i4_read ( icol, line, prompt, ierror )
  if ( ierror /= 0 ) then
    return
  end if
!
!  Read the divisor, either RVAL or ISTOP/ISBOT.
!
  if ( iform == 0 ) then
    call rat_read ( istop, isbot, line, prompt, ierror )
  else if ( iform == 1 ) then
    call r4_read ( sval, line, prompt, ierror )
  else if ( iform == 2 ) then
    call dec_read ( istop, isbot, line, prompt, ierror )
    if ( ierror /= 0 ) then
      call dec_round ( istop, isbot )
    end if
  else if ( iform == 3 ) then
    call i4_read ( istop, line, prompt, ierror )
  else if ( iform == 4 ) then
    call i4_read ( istop, line, prompt, ierror )
    istop = mod ( istop, base )
  end if

  return
end
subroutine col_mul ( a, a_int, iatop, iabot, ierror, base, iform, icol, &
  nrow, ncol, sval, istop, isbot )

!*****************************************************************************80
!
!! COL_MUL multiplies a column of the matrix by a scale factor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(NROW,NCOL), the current matrix.
!
!    Input/output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input/output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Output, integer ( kind = 4 ) IERROR, 0 an error occurred, 1 no error occurred.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) ICOL, the column that is to be multiplied.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
!    Input, real ( kind = 4 ) SVAL, the real row multiplier.
!
!    Input, integer ( kind = 4 ) ISTOP, ISBOT, the decimal or fractional row multiplier.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  integer ( kind = 4 ) base
  character ( len = 10 ) chrtmp1
  character ( len = 22 ) chrtmp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) isbot
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) itop
  real ( kind = 4 ) sval
!
!  Make sure column number is OK.
!
  if ( icol < 1 .or. ncol < icol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Error!  The column number is out of range!'
    ierror = 1
    return
  end if
!
!  For rational arithmetic, make sure bottom of scale factor
!  is not 0.
!
  if ( iform == 0 .and. isbot == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Error!  Illegal 0 divisor in multiplier!'     
    ierror = 1
    return
  end if
!
!  Check for multiplication by 0.
!
  if ( ( iform == 0 .and. istop == 0 ) .or. &
       ( iform == 1 .and. sval == 0.0E+00 ) .or. &
       ( iform == 2 .and. istop == 0 ) .or. &
       ( iform == 3 .and. istop == 0 ) .or. &
       ( iform == 4 .and. istop == 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Warning - Multiplication by zero is not an ERO.'
    ierror = 1
    return
  end if
!
!  Check for multiplication by 1.
!
  if ( ( iform == 0 .and. istop == isbot ) .or. &
       ( iform == 1 .and. sval == 1.0E+00 ) .or. &
       ( iform == 2 .and. istop == 1 .and. isbot == 0 ) .or. &
       ( iform == 3 .and. istop == 1 ) .or. &
       ( iform == 4 .and. istop == 1 ) ) then
      return
  end if
!
!  Carry out the multiplication.
!
  if ( iform == 0 ) then
 
    do i = 1, nrow
 
      call rat_mul ( ibot, iabot(i,icol), isbot, itop, iatop(i,icol), istop, &
        ierror )

      if ( ierror /= 0 ) then
        return
      end if

      iatop(i,icol) = itop
      iabot(i,icol) = ibot
 
    end do
 
  else if ( iform == 1 ) then
 
    a(1:nrow,icol) = sval * a(1:nrow,icol)
 
  else if ( iform == 2 ) then
 
    do i = 1, nrow
 
      call dec_mul ( ibot, iabot(i,icol), isbot, itop, iatop(i,icol), istop )

      if ( ierror /= 0 ) return
 
      iatop(i,icol) = itop
      iabot(i,icol) = ibot
 
    end do

  else if ( iform == 3 ) then
 
    a_int(1:nrow,icol) = istop * a_int(1:nrow,icol)

  else if ( iform == 3 ) then
 
    a_int(1:nrow,icol) = mod ( istop * a_int(1:nrow,icol), base )
 
  end if
!
!  Confirm the operation.
!
  if ( iform == 0 ) then
 
    if ( istop == - isbot ) then
      chrtmp = '-'
    else
      call rat_to_s_left ( istop, isbot, chrtmp )
    end if
 
  else if ( iform == 1 ) then
 
    if ( sval == - 1.0E+00 ) then
      chrtmp = '-'
    else
      call r4_to_s_left ( sval, chrtmp )
    end if
 
  else if ( iform == 2 ) then

    if ( istop == - 1 .and. isbot == 0 ) then
      chrtmp = '-'
    else
      call dec_to_s_left ( istop, isbot, chrtmp )
    end if

  else if ( iform == 3 ) then
 
    if ( istop == - 1 ) then
      chrtmp = '-'
    else
      call i4_to_s_left ( istop, chrtmp )
    end if

  else if ( iform == 4 ) then
 
    if ( istop == - 1 ) then
      chrtmp = '-'
    else
      call i4_to_s_left ( istop, chrtmp )
    end if

  end if
 
  call i4_to_s_left ( icol, chrtmp1 )

  write ( *, '(a)' ) '  ECO: Col ' // trim ( chrtmp1 ) // ' <= ' // trim ( chrtmp ) // &
    ' Col ' // trim ( chrtmp1 )
  
  return
end
subroutine col_mul_param ( ierror, base, iform, icol, istop, isbot, line, rval )

!*****************************************************************************80
!
!! COL_MUL_PARAM gets and checks the column multiply parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IERROR, 0 an error occurred, 1 no error occurred.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Output, integer ( kind = 4 ) ICOL, the column to be multiplied.
!
!    Output, integer ( kind = 4 ) ISTOP, ISBOT, the multiplier to use for
!    fractional or decimal arithmetic.
!
!    Workspace, character ( len = 80 ) LINE, used to hold the user's input.
!
!    Output, real ( kind = 4 ) RVAL, the multiplier to use for real arithmetic.
!
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) isbot
  integer ( kind = 4 ) istop
  character ( len = 80 ) line
  character ( len = 80 ) prompt
  real ( kind = 4 ) rval

  prompt = 'column I, multiplier S.'
!
!  Read the column number to be multiplied.
!
  call i4_read ( icol, line, prompt, ierror )
  if ( ierror /= 0 ) return
!
!  Read the multiplier, either RVAL or ISTOP/ISBOT.
!
  if ( iform == 0 ) then
    call rat_read ( istop, isbot, line, prompt, ierror )
  else if ( iform == 1 ) then
    call r4_read ( rval, line, prompt, ierror )
  else if ( iform == 2 ) then
    call dec_read ( istop, isbot, line, prompt, ierror )
    if ( ierror == 0 ) then   
      call dec_round ( istop, isbot )
    end if
  else if ( iform == 3 ) then
    call i4_read ( istop, line, prompt, ierror )
  else if ( iform == 4 ) then
    call i4_read ( istop, line, prompt, ierror )
    istop = mod ( istop, base )
  end if
 
  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COL_MUL_PARAM - Error!'
    write ( *, '(a)' ) '  Your input could not be understood.'
    return
  end if

  return
end
subroutine col_op_check ( command, ierror, line2 )

!*****************************************************************************80
!
!! COL_OP_CHECK checks for commands given in the form of ECO's.
!
!  Discussion:
!
!    The form of the elementary column operation commands includes:
!
!    The column interchange command:
!      CI1 <=> CI2
!    Note that this will fail if user types "C I1 <=> C I2"
!
!    The scalar multiply command:
!      CI1 <= S * CI1
!    with or without the "*".
!
!    The scalar divide command:
!      CI1 <= CI1 / S
!
!    The add row command:
!      CI1 <= CI1 + S *CI2
!    or
!      CI1 <= S * CI2 + CI1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = 4 ) COMMAND.
!    If the routine decides that the user has input an ERO in the
!    natural format, then COMMAND contains the necessary
!    one letter MATMAN command to carry out the ERO.
!
!    Output, integer ( kind = 4 ) IERROR, 0 an error occurred, 1 no error occurred.
!
!    Workspace, character ( len = 80 ) LINE2, a copy of the user input in LINE.
!
  implicit none

  character ( len = 22 ) chrtmp
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 20 ) command
  integer ( kind = 4 ) icol1
  integer ( kind = 4 ) icol2
  integer ( kind = 4 ) icol3
  integer ( kind = 4 ) idbot2
  integer ( kind = 4 ) idbot3
  integer ( kind = 4 ) idtop2
  integer ( kind = 4 ) idtop3
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isbot2
  integer ( kind = 4 ) isbot3
  integer ( kind = 4 ) istop2
  integer ( kind = 4 ) istop3
  integer ( kind = 4 ) lchar
  logical ldiv
  character ( len = 80 ) line2
  character ( len = 80 ) string

  command = ' '
!
!  1. Remove all blanks from the line, and capitalize it.
!
  call s_blank_delete ( line2 )
  call s_cap ( line2 )
!
!  2. Is the first character a "C" or "COL" or "COLUMN"?
!
  if ( line2(1:1) /= 'C' ) then
    return
  end if
 
  if ( line2(1:6) == 'COLUMN' ) then
    call s_chop ( line2, 1, 6 )
  else if ( line2(1:3) == 'COL' ) then
    call s_chop ( line2, 1, 3 )
  else
    call s_chop ( line2, 1, 1 )
  end if
!
!  3. The next item should be a column number, ICOL1.
!
  call s_to_i4 ( line2, icol1, ierror, lchar )
 
  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Your ECO command could not be understood.'
    write ( *, '(a)' ) '  The first column number "R1" did not make sense.'
    return
  end if
 
  call s_chop ( line2, 1, lchar )
!
!  4. Check for the column interchange string "=", "<>", "<=>" or "<->".
!
  if ( line2(1:2) == '<>' ) then
    string = '<>'
  else if ( line2(1:3) == '<=>' ) then
    string = '<=>'
  else if ( line2(1:3) == '<->' ) then
    string = '<->'
  else if ( line2(1:2) == '<=' ) then
    string = '<='
  else if ( line2(1:2) == '<-' ) then
    string = '<-'
  else if ( line2(1:2) == '=>' ) then
    string = '=>'
  else if ( line2(1:2) == '->' ) then
    string = '->'
  else if ( line2(1:1) == '=' ) then
    string = '='
  else if ( line2(1:2) == ':=' ) then
    string = ':='
  else
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Your ECO command could not be understood.'
    write ( *, '(a)' ) '  The assignment symbol <=> was missing.'
    return
  end if
 
  lchar = len_trim ( string )
 
  call s_chop ( line2, 1, lchar )
!
!  5. The next quantity could be an explicit signed scalar, S2,
!     or an implicit +-1.
!
  if ( line2(1:1) == 'C' ) then

    istop2 = 1.0E+00
    isbot2 = 1.0E+00

  else

    if ( line2(1:2) == '+C' ) then
      istop2 = 1.0E+00
      isbot2 = 1.0E+00
      call s_chop ( line2, 1, 1 )
    else if ( line2(1:2) == '-C' ) then
      istop2 = - 1.0E+00
      isbot2 = 1.0E+00
      call s_chop ( line2, 1, 1 )
    else
      call chrctg ( line2, istop2, isbot2, ierror, lchar )
      call s_chop ( line2, 1, lchar )
 
      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Your ECO command could not be understood.'
        write ( *, '(a)' ) '  The multiplier S2 did not make sense.'
        ierror = 1
        return
      end if
 
    end if
  end if
!
!  6. Is the next character an optional "*"?
!
  if ( line2(1:1) == '*' ) then
    call s_chop ( line2, 1, 1 )
  end if
!
!  7. Is the next character a "C"?
!
  if ( line2(1:6) == 'COLUMN' ) then
    call s_chop ( line2, 1, 6 )
  else if ( line2(1:3) == 'COL' ) then
    call s_chop ( line2, 1, 3 )
  else if ( line2(1:1) == 'C' ) then
    call s_chop ( line2, 1, 1 )
  else
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Your ECO command could not be understood.'
    write ( *, '(a)' ) '  Could not find the second column index.'
    
    return
  end if
!
!  8. The next item should be a column number, ICOL2.
!
  call s_to_i4 ( line2, icol2, ierror, lchar )
 
  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Your ECO command could not be understood.'
    write ( *, '(a)' ) '  The index of the second column "C2" did not make sense.'
    ierror = 1
    return
  end if

  call s_chop ( line2, 1, lchar )
!
!  9. If there's nothing more, this must be an interchange
!     or a scaling.  Form the equivalent MATMAN command.
!
  if ( line2 == ' ' ) then
 
    if ( icol1 == icol2 ) then

      command = 'COL_MUL'

      if ( isbot2 < 0 ) then
        isbot2 = - isbot2
        istop2 = - istop2
      end if

      call rat_to_s_left ( istop2, isbot2, chrtmp )
      call i4_to_s_left ( icol1, chrtmp1 )
      line2 = chrtmp1 // ' ' // chrtmp
      call s_blanks_delete ( line2 )
      return
    end if
 
    if ( istop2 == 1 .and. isbot2 == 1 ) then
      command = 'COL_SWAP'
      call i4_to_s_left ( icol1, chrtmp1 )
      call i4_to_s_left ( icol2, chrtmp2 )
      line2 = chrtmp1 // ' ' // chrtmp2
      call s_blanks_delete ( line2 )
      return
    end if
 
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Your ECO command could not be understood.'
    write ( *, '(a)' ) '  A MULTIPLY command must have C1 and C2 the same.'
    write ( *, '(a)' ) '  An INTERCHANGE command cannot have a multiplier.'
    return
  end if
!
!  10. Is the next quantity a '/', or perhaps a '*'?
!
  ldiv = .false.
 
  if ( line2(1:1) == '/' ) then
 
    ldiv = .true.
    call s_chop ( line2, 1, 1 )
    call chrctg ( line2, idtop2, idbot2, ierror, lchar )
 
    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Your ECO command could not be understood.'
      write ( *, '(a)' ) '  The divisor of column 2 did not make sense.'
      return
    end if
 
    istop2 = istop2 * idbot2
    isbot2 = isbot2 * idtop2
 
    if ( icol1 == icol2 ) then

      if ( ldiv ) then
        command = 'COL_DIV'
        call i4_swap ( istop2, isbot2 )
      else
        command = 'COL_MUL'
      end if

      if ( isbot2 < 0 ) then
        isbot2 = - isbot2
        istop2 = - istop2
      end if

      call i4_to_s_left ( icol1, chrtmp1 )
      call rat_to_s_left ( istop2, isbot2, chrtmp )
      line2 = chrtmp1 // ' ' // chrtmp
      call s_blanks_delete ( line2 )
      return

    end if
 
  else if ( line2(1:1) == '*' ) then
 
    call s_chop ( line2, 1, 1 )
    call chrctg ( line2, idtop2, idbot2, ierror, lchar )
 
    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Your ECO command could not be understood.'
      write ( *, '(a)' ) '  The multiplier of column 2 did not make sense.'
      return
    end if
 
    istop2 = istop2 * idtop2
    isbot2 = isbot2 * idbot2
 
    if ( icol1 == icol2 ) then

      command = 'COL_MUL'

      if ( isbot2 < 0 ) then
        isbot2 = - isbot2
        istop2 = - istop2
      end if

      call i4_to_s_left ( icol1, chrtmp1 )
      call rat_to_s_left ( istop2, isbot2, chrtmp )
      line2 = chrtmp1 // ' ' // chrtmp
      call s_blanks_delete ( line2 )
      return

    end if
 
  end if
!
!  11. Is the next quantity a scalar, S3?
!
  if ( line2(1:2) == '+C' ) then
 
    istop3 = 1.0E+00
    isbot3 = 1.0E+00
    call s_chop ( line2, 1, 1) 
 
  else if ( line2(1:2) == '-C' ) then
 
    istop3 = - 1.0E+00
    isbot3 = 1.0E+00
    call s_chop ( line2, 1, 1 )
 
  else
 
    call chrctg ( line2, istop3, isbot3, ierror, lchar )
 
    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Your ECO command could not be understood.'
      write ( *, '(a)' ) '  The multiplier S2 did not make sense.'
      ierror = 1
      return
    end if
 
    call s_chop ( line2, 1, lchar )
 
  end if
!
!  12. Is the next quantity an optional "*"?
!
  if ( line2(1:1) == '*' ) then
    call s_chop ( line2, 1, 1 )
  end if
!
!  13. Is the next quantity a "C"?
!
  if ( line2(1:6) == 'COLUMN' ) then
    call s_chop ( line2, 1, 6 )
  else if ( line2(1:3) == 'COL' ) then
    call s_chop ( line2, 1, 3 )
  else if ( line2(1:1) == 'C' ) then
    call s_chop ( line2, 1, 1) 
  else
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Your ECO command could not be understood.'
    write ( *, '(a)' ) '  The "C" marking the third column was misplaced.'
    return
  end if
!
!  14. The next item should be a column number, ICOL3.
!
  call s_to_i4 ( line2, icol3, ierror, lchar )
 
  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Your ECO command could not be understood.'
    write ( *, '(a)' ) '  The third column number "C3" did not make sense.'
    ierror = 1
    return
  end if
 
  call s_chop ( line2, 1, lchar )
!
!  15. Is the next quantity a '/', or perhaps a '*'?
!
  if ( line2(1:1) == '/' ) then
 
    call s_chop ( line2, 1, 1 )
    call chrctg ( line2, idtop3, idbot3, ierror, lchar )
 
    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Your ECO command could not be understood.'
      write ( *, '(a)' ) '  The divisor of column 3 did not make sense.'
      return
    end if
 
    istop3 = istop3 * idbot3
    isbot3 = isbot3 * idtop3
 
  else if ( line2(1:1) == '*' ) then
 
    call s_chop ( line2, 1, 1 )
    call chrctg ( line2, idtop3, idbot3, ierror, lchar )
 
    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Your ECO command could not be understood.'
      write ( *, '(a)' ) '  The multiplier of column 3 did not make sense.'
      return
    end if
 
    istop3 = istop3 * idtop3
    isbot3 = isbot3 * idbot3
 
  end if
!
!  16. Form the equivalent MATMAN command.
!
  if ( icol1 == icol2 ) then

    command = 'COL_ADD'

    if ( isbot3 < 0 ) then
      isbot3 = - isbot3
      istop3 = - istop3
    end if

    call i4_to_s_left ( icol3, chrtmp1 )
    call i4_to_s_left ( icol1, chrtmp2 )
    call rat_to_s_left ( istop3, isbot3, chrtmp )
    line2 =  chrtmp // ' ' // chrtmp1 // ' ' // chrtmp2
    call s_blanks_delete ( line2 )

  else if ( icol1 == icol3 ) then

    command = 'COL_ADD'

    if ( isbot2 < 0 ) then
      isbot2 = - isbot2
      istop2 = - istop2
    end if

    call rat_to_s_left ( istop2, isbot2, chrtmp )
    call i4_to_s_left ( icol2, chrtmp1 )
    call i4_to_s_left ( icol1, chrtmp2 )
    line2 = chrtmp // ' ' // chrtmp1 // ' ' // chrtmp2
    call s_blanks_delete ( line2 )

  else
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Your ECO command could not be understood.'
    write ( *, '(a)' ) '  C2 or C3 must equal C1 in an ECO command.'
  end if
 
  return
end
subroutine col_swap ( a, a_int, iatop, iabot, ierror, iform, icol1, icol2, &
  nrow, ncol )

!*****************************************************************************80
!
!! COL_SWAP swaps two columns of a matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(NROW,NCOL), the current matrix.
!
!    Input/output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input/output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Output, integer ( kind = 4 ) IERROR, 0 an error occurred, 1 no error occurred.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) ICOL1, ICOL2, the rows to swap.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol1
  integer ( kind = 4 ) icol2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
!
!  Skip out if the two columns are the same.
!
  if ( icol1 == icol2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COL_SWAP - Error!'
    write ( *, '(a)' ) '  You have asked to swap a column with itself!'
    return
  end if
!
!  Refuse to continue if a row is out of bounds.
!
  if ( ( icol1 < 1 .or. ncol < icol1 ) .or. &
       ( icol2 < 1 .or. ncol < icol2 ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COL_SWAP - Error!'
    write ( *, '(a)' ) '  One of the columns is illegal!'
    return
  end if
!
!  Swap the columns.
!
  do i = 1, nrow
 
    if ( iform == 0 ) then
 
      call i4_swap ( iatop(i,icol1), iatop(i,icol2) )
      call i4_swap ( iabot(i,icol1), iabot(i,icol2) )
 
    else if ( iform == 1 ) then
 
      call r4_swap ( a(i,icol1), a(i,icol2) )
 
    else if ( iform == 2 ) then
 
      call i4_swap ( iatop(i,icol1), iatop(i,icol2) )
      call i4_swap ( iabot(i,icol1), iabot(i,icol2) )
 
    else if ( iform == 3 ) then

      call i4_swap ( a_int(i,icol1), a_int(i,icol2) )

    end if
 
  end do
 
  call i4_to_s_left ( icol1, chrtmp1 )
  call i4_to_s_left ( icol2, chrtmp2 )

  write ( *, '(a)' ) '  ECO: Col ' // trim ( chrtmp1 ) // ' <=> Col ' &
    // trim ( chrtmp2 )

  return
end
subroutine dec_add ( ibot, ibot1, ibot2, itop, itop1, itop2 )

!*****************************************************************************80
!
!! DEC_ADD adds two decimal quantities.
!
!  Discussion:
!
!    The routine computes
!  
!      ITOP * 10**IBOT = ITOP1 * 10**IBOT1 + ITOP2 * 10**IBOT2
!
!    while trying to avoid integer overflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IBOT, the exponent of the result.
!
!    Input, integer ( kind = 4 ) IBOT1, IBOT2, the exponents of the numbers to be added.
!
!    Output, integer ( kind = 4 ) ITOP, the coefficient of the result.
!
!    Input, integer ( kind = 4 ) ITOP1, ITOP2, the coefficients of the numbers to be added.
!
  implicit none

  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ibot1
  integer ( kind = 4 ) ibot2
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) itop1
  integer ( kind = 4 ) itop2
  integer ( kind = 4 ) jtop1
  integer ( kind = 4 ) jtop2

  if ( itop1 == 0 ) then
    itop = itop2
    ibot = ibot2
    return
  else if ( itop2 == 0 ) then
    itop = itop1
    ibot = ibot1
    return
  else if ( ibot1 == ibot2 ) then
    itop = itop1 + itop2
    ibot = ibot1
    call dec_round ( itop, ibot )
    return
  end if
!
!  Line up the exponents.
!
  jtop1 = itop1
  jtop2 = itop2
 
  if ( ibot1 < ibot2 ) then
    jtop2 = jtop2 * 10**(ibot2-ibot1)
  else
    jtop1 = jtop1 * 10**(ibot1-ibot2)
  end if
!
!  Add the coefficients.
!
  itop = jtop1 + jtop2
  ibot = min ( ibot1, ibot2 )
!
!  Clean up the result.
!
  call dec_round ( itop, ibot )
 
  return
end
subroutine dec_digit_set ( ierror, line )

!*****************************************************************************80
!
!! DEC_DIGIT_SET allows the user to specify the number of decimal digits.
!
!  Discussion:
!
!    DEC_DIGIT is
!
!    * the number of digits used when converting a real number
!      to a fraction using the "FI" or "FD" command;
!
!    * the maximum number of digits in a decimal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IERROR, 0 an error occurred, 1 no error occurred.
!
!    Workspace, character ( len = 80 ) LINE, used to hold the user's input.
!
  implicit none

  character ( len = 10 ) chrtmp1
  integer ( kind = 4 ) dec_digit_max
  integer ( kind = 4 ) dec_digit
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) itemp
  character ( len = 80 ) line
  character ( len = 80 ) prompt

  dec_digit_max = 0
  call i4_data ( 'GET', 'DEC_DIGIT_MAX', dec_digit_max )

  write ( *, '(a)' ) 'How many decimal places should be used in '
  write ( *, '(a)' ) 'converting real results to a decimal?'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' 1 means 123.45 becomes 1 * 10**2'
  write ( *, '(a)' ) ' 2 means 123.45 becomes 12 * 10**1'
  write ( *, '(a)' ) ' 3 means 123.45 becomes 123'
  write ( *, '(a)' ) 'and so on.'
  
  call i4_to_s_left ( dec_digit_max, chrtmp1 )
  prompt = 'number of decimals (1 to ' // chrtmp1 // ').'
  call s_blanks_delete ( prompt )
 
  call i4_read ( itemp, line, prompt, ierror )
 
  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Your choice was not acceptable!'
    return
  end if
!
!  Absolutely do not let DEC_DIGIT be less than 1.
!
  if ( itemp < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The number of decimals must be positive!'
    ierror = 1
    return
  end if
!
!  Allow user to exceed the maximum, with a warning.
!
  if ( dec_digit_max < itemp ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Warning!'
    write ( *, '(a)' ) '  Your choice is larger than the recommended maximum!'
    write ( *, '(a)' ) '  which is ' // trim ( chrtmp1 )
    write ( *, '(a)' ) '  It is possible that calculations will break down'
    write ( *, '(a)' ) '  at any time!  Be careful!'
  end if
  
  dec_digit = itemp
  call i4_data ( 'SET', 'DEC_DIGIT', dec_digit )

  call i4_to_s_left ( dec_digit, chrtmp1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The number of decimal digits will now be ' // chrtmp1
  
  return
end
subroutine dec_div ( itop1, ibot1, itop2, ibot2, itop, ibot, ierror )

!*****************************************************************************80
!
!! DEC_DIV divides two decimal values.
!
!  Discussion:
!
!    A decimal quantity is stored as
!
!      (ITOP,IBOT) 
!
!    representing the value
!
!      ITOP * 10 ** IBOT.
!
!    The routine computes 
!
!      ITOP * 10**IBOT = (ITOP1 * 10**IBOT1) / (ITOP2 * 10**IBOT2)
!
!                      = (ITOP1/ITOP2) * 10**(IBOT1-IBOT2)
!
!    while avoiding integer overflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ITOP1, IBOT1, the numerator.
!
!    Input, integer ( kind = 4 ) ITOP2, IBOT2, the denominator.
!
!    Output, integer ( kind = 4 ) ITOP, IBOT, the result.
!
!    Output, integer ( kind = 4 ) IERROR, 0 an error occurred, 1 no error occurred.
!
  implicit none

  real    ( kind = 8 ) dval
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ibot1
  integer ( kind = 4 ) ibot2
  integer ( kind = 4 ) ibot3
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) itop1
  integer ( kind = 4 ) itop2
  integer ( kind = 4 ) itop3
!
!  First special case, top fraction is 0.
!
  if ( itop1 == 0 ) then
    itop = 0
    ibot = 0
    return
  end if
!
!  First error, bottom of fraction is 0.
!
  if ( itop2 == 0 ) then
    ierror = 1
    itop = 0
    ibot = 0
    return
  end if
!
!  Second special case, result is 1.
!
  if ( itop1 == itop2 .and. ibot1 == ibot2 ) then
    itop = 1
    ibot = 0
    return
  end if
!
!  Third special case, result is power of 10.
!
  if ( itop1 == itop2 ) then
    itop = 1
    ibot = ibot1 - ibot2
    return
  end if
!
!  Fourth special case: ITOP1/ITOP2 is exact.
!
  if ( ( itop1 / itop2 ) * itop2 == itop1 ) then
    itop = itop1 / itop2
    ibot = ibot1 - ibot2
    return
  end if
!
!  General case.
!
  dval = real ( itop1, kind = 8 ) / real ( itop2, kind = 8 )
 
  call r8_to_dec ( dval, itop3, ibot3 )
 
  itop = itop3
  ibot = ibot3 + ibot1 - ibot2
 
  return
end
subroutine dec_mul ( ibot, ibot1, ibot2, itop, itop1, itop2 )

!*****************************************************************************80
!
!! DEC_MUL multiplies two decimals.
!
!  Discussion:
!
!    The routine computes
!
!      ITOP * 10**IBOT = (ITOP1 * 10**IBOT1) * (ITOP2 * 10**IBOT2)
!                      = (ITOP1*ITOP2) * 10**(IBOT1+IBOT2)
!
!    while avoiding integer overflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IBOT, the exponent of the result.
!
!    Input, integer ( kind = 4 ) IBOT1, the exponent of the first factor.
!
!    Input, integer ( kind = 4 ) IBOT2, the exponent of the second factor.
!
!    Output, integer ( kind = 4 ) ITOP, the coefficient of the result.
!
!    Input, integer ( kind = 4 ) ITOP1, the coefficient of the first factor.
!
!    Input, integer ( kind = 4 ) ITOP2, the coefficient of the second factor.
!
  implicit none

  integer ( kind = 4 ) i4_big
  real    ( kind = 8 ) dval
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ibot1
  integer ( kind = 4 ) ibot2
  integer ( kind = 4 ) ibot3
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) itop1
  integer ( kind = 4 ) itop2
  integer ( kind = 4 ) itop3
  real    ( kind = 4 ) rmax
  real    ( kind = 4 ) temp

  i4_big = 0
  call i4_data ( 'GET', 'I4_BIG', i4_big )
  rmax = real ( i4_big )
!
!  The result is zero if either ITOP1 or ITOP2 is zero.
!
  if ( itop1 == 0 .or. itop2 == 0 ) then
    itop = 0
    ibot = 0
    return
  end if
!
!  The result is simple if either ITOP1 or ITOP2 is one.
!
  if ( itop1 == 1 .or. itop2 == 1 ) then
    itop = itop1 * itop2
    ibot = ibot1 + ibot2
    return
  end if
 
  temp = log ( real ( abs ( itop1 ) ) ) + log ( real ( abs ( itop2 ) ) )
 
  if ( temp < log ( rmax ) ) then
 
    itop = itop1 * itop2
    ibot = ibot1 + ibot2
 
  else
 
    dval = real ( itop1, kind = 8 ) * real ( itop2, kind = 8 )
 
    call r8_to_dec ( dval, itop3, ibot3 )
 
    itop = itop3
    ibot = ibot3 + ( ibot1 + ibot2 )
 
  end if
!
!  Clean up the result.
!
  call dec_round ( itop, ibot )
 
  return
end
subroutine dec_print ( iatop, iabot, nrow, ncol, title )

!*****************************************************************************80
!
!! DEC_PRINT prints out decimal vectors and matrices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL).
!    IATOP and IABOT represent the decimal matrix.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
!    Input, character ( len = * ) TITLE, a label for the object being printed.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  character ( len = 22 ) chrtmp
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 10 ) chrtmp3
  character ( len = 40 ) format2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) ichi
  integer ( kind = 4 ) iclo
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) izhi
  integer ( kind = 4 ) izlo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) jmin
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) lenc
  integer, parameter :: ncolum = 80
  integer ( kind = 4 ) npline
  character ( len = 100 ) output
  character ( len = * ) title
!
!  Figure out how wide we must make each column.
!
  imax = 0
  jmax = 0
 
  do i = 1, nrow
    do j = 1, ncol
 
      call dec_to_s_left ( iatop(i,j), iabot(i,j), chrtmp )
      lenc = len_trim ( chrtmp )
      jmax = max ( jmax, lenc )
 
    end do
  end do
 
  kmax = 2 + imax + 1 + jmax
  npline = ncolum / kmax
!
!  Set up the format for the heading.
!
    call i4_to_s_left ( npline, chrtmp2 )
    call i4_to_s_left ( kmax, chrtmp3 )
    format2 = '(' // chrtmp2 // 'i' // chrtmp3 // ')'
 
  call s_blank_delete ( format2 )
 
  do jmin = 1, ncol, npline
 
    jmax = min ( jmin+npline-1, ncol )
 
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(a)' ) ' '
 
    if ( 1 < jmin .or. jmax < ncol ) then
      write ( output, * ) 'Columns ', jmin, ' to ', jmax
      call s_blanks_delete ( output )
      write ( *, '(a)' ) trim ( output )
      write ( *, '(a)' ) ' '
    end if
 
    do i = 1, nrow
 
      output = ' '
 
      do j = jmin, jmax
        klo = 4 + (j-jmin) * kmax + 1
        khi = 4 + (j-jmin) * kmax + kmax
        call dec_to_s_left ( iatop(i,j), iabot(i,j), chrtmp )
        output(klo:khi) = adjustr ( chrtmp(1:kmax) )
      end do
 
      write ( *, '(a)' ) trim ( output )

    end do
 
  end do
 
  return
end
subroutine dec_read ( itop, ibot, line, prompt, ierror )

!*****************************************************************************80
!
!! DEC_READ reads a decimal, rational or integer, and returns a decimal fraction.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ITOP, IBOT, represents the decimal fraction.
!
!    Workspace, character ( len = 80 ) LINE, used to hold the user's input.
!
!    Input/output, character ( len = 80 ) PROMPT, the prompt string.
!
!    Output, integer ( kind = 4 ) IERROR, 0 an error occurred, 1 no error occurred.
!
  implicit none

  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ibot1
  integer ( kind = 4 ) ibot2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) itop1
  integer ( kind = 4 ) itop2
  integer ( kind = 4 ) lchar
  character ( len = 80 ) line
  character ( len = 80 ) prompt

  ierror = 0
  itop = 0
  ibot = 0
 
  do
 
    call chrinp ( ierror, line, prompt )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DEC_READ - Fatal error!'
      write ( *, '(a)' ) '  CHRINP returned error flag!'
      return
    end if
 
    if ( line /= ' ' ) then
      exit
    end if
 
  end do

  call s_to_dec ( line, itop1, ibot1, lchar )

  if ( len ( line ) <= lchar ) then
    itop = itop1
    ibot = ibot1
  else if ( line(lchar+1:lchar+1) /= '/' ) then
    itop = itop1
    ibot = ibot1
  else
    call s_chop ( line, 1, lchar+1 )
    call s_to_dec ( line, itop2, ibot2, lchar )
    call dec_div ( itop1, ibot1, itop2, ibot2, itop, ibot, ierror )
  end if
 
  call s_chop ( line, 1, lchar )

  return
end
subroutine dec_round ( itop, ibot )

!*****************************************************************************80
!
!! DEC_ROUND rounds a decimal fraction to a given number of digits.
!
!  Discussion:
!
!    The routine takes an arbitrary decimal fraction represented by
!
!      ITOP * 10**IBOT
!
!    and makes sure that ITOP has no more than the allowed number of
!    decimal digits.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) ITOP, IBOT, the coefficient and exponent
!    of a decimal fraction.  On return, ITOP has no more than 
!    the allowed number of decimal digits.
!
  implicit none

  integer ( kind = 4 ) dec_digit
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) itop

  if ( itop == 0 ) then
    ibot = 0
    return
  end if

  dec_digit = 0
  call i4_data ( 'GET', 'DEC_DIGIT', dec_digit )
  
  do while ( 10**dec_digit <= abs ( itop ) )
    itop = itop / 10
    ibot = ibot + 1
  end do
  
  do while ( ( itop / 10 ) * 10 == itop )
    itop = itop / 10
    ibot = ibot + 1
  end do
 
  return
end
subroutine dec_to_r4 ( itop, ibot, a )

!*****************************************************************************
!
!! DEC_TO_R4 converts a decimal ITOP * 10**IBOT to an R4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ITOP, IBOT, the coefficient and exponent
!    of the decimal value.
!
!    Output, real ( kind = 4 ) A, the equivalent real value.
!
  implicit none

  real ( kind = 4 ) a
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) itop

  a = itop * 10.0E+00**ibot
 
  return
end
subroutine dec_to_rat ( iatop, iabot )

!*****************************************************************************80
!
!! DEC_TO_RAT converts a decimal to a rational representation.
!
!  Discussion:
!
!    On input, a value is represented as IATOP * 10**IABOT.
!
!    On output, approximately the same value is represented as IATOP / IABOT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) IATOP, IABOT.
!    On input, these quantities represent the value IATOP * 10 ** IABOT.
!    On output, these quantities represent the value IATOP / IABOT.
!
  implicit none

  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) iabot
  integer ( kind = 4 ) iatop
  integer ( kind = 4 ) itmp

  if ( 0 <= iabot ) then
    iatop = iatop * 10**iabot
    iabot = 1
  else
    iabot = 10**(-iabot)
    itmp = i4_gcd ( iatop, iabot )
    iatop = iatop / itmp
    iabot = iabot / itmp
  end if
 
  return
end
subroutine dec_to_s_left ( ival, jval, s )

!*****************************************************************************80
!
!! DEC_TO_S_LEFT returns a left-justified representation of IVAL * 10**JVAL.
!
!  Example:
!
!    IVAL     JVAL       S
!    ----     ----       ------
!       0        0       0
!      21        3       21000
!      -3        0       -3
!     147       -2       14.7
!      16       -5       0.00016
!      34       30       Inf
!     123      -21       0.0000000000000000012
!      34      -30       0.0E+00
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, JVAL, integers which represent the decimal.
!
!    Output, character ( len = * ) S, the representation of the value.
!    The string is 'Inf' or '0.0' if the value was too large
!    or small to represent with a fixed point format.
!
  implicit none

  character ( len = 22 ) chrrep
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iget1
  integer ( kind = 4 ) iget2
  integer ( kind = 4 ) iput1
  integer ( kind = 4 ) iput2
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jval
  integer ( kind = 4 ) maxdigit
  integer ( kind = 4 ) ndigit
  integer ( kind = 4 ) nleft
  character ( len = * ) s

  s = ' '

  if ( ival == 0 ) then
    s = '0'
    return
  end if

  maxdigit = len ( s )
!
!  Store a representation of IVAL in CHRREP.
!
  write ( chrrep, '(i22)' ) ival
  call s_blank_delete ( chrrep )
  ndigit = len_trim ( chrrep )
!
!  Overflow if JVAL is positive, and MAXDIGIT < NDIGIT + JVAL.
!
  if ( 0 < jval ) then
    if ( maxdigit < ndigit + jval ) then
      s = 'Inf'
      return
    end if
  end if
!
!  Underflow if JVAL is negative, and 3 + NDIGIT - JVAL > MAXDIGIT.
!
  if ( jval < 0 ) then
    if ( 0 < ival ) then
      if ( maxdigit < 3 - ndigit - jval ) then
        s = '0.0'
        return
      end if
    else
      if ( maxdigit < 5 - ndigit - jval ) then
        s = '0.0'
        return
      end if
    end if
  end if
!
!  If JVAL is nonnegative, insert trailing zeros.
!
  if ( 0 <= jval ) then

    s(1:ndigit) = chrrep(1:ndigit)

    do i = ndigit+1, ndigit+jval
      s(i:i) = '0'
    end do

  else if ( jval < 0 ) then

    iput2 = 0
    iget2 = 0
!
!  Sign.
!
    if ( ival < 0 ) then
      iput1 = 1
      iput2 = 1
      iget2 = 1
      s(iput1:iput2) = '-'
      ndigit = ndigit - 1
    end if
!
!  Digits of the integral part.
!
    if ( 0 < ndigit + jval ) then
      iput1 = iput2 + 1
      iput2 = iput1 + ndigit + jval -1
      iget1 = iget2 + 1
      iget2 = iget1 + ndigit+jval - 1
      s(iput1:iput2) = chrrep(iget1:iget2)
    else
      iput1 = iput2 + 1
      iput2 = iput1
      s(iput1:iput2) = '0'
    end if
!
!  Decimal point.
!
    iput1 = iput2 + 1
    iput2 = iput1
    s(iput1:iput2) = '.'
!
!  Leading zeroes.
!
    do i = 1, - jval - ndigit
      iput1 = iput2 + 1
      iput2 = iput1
      s(iput1:iput2) = '0'
    end do

    nleft = min ( -jval, ndigit )
    nleft = min ( nleft, maxdigit - iput2 )
    iput1 = iput2 + 1
    iput2 = iput1 + nleft - 1
    iget1 = iget2 + 1
    iget2 = iget1 + nleft - 1
    s(iput1:iput2) = chrrep(iget1:iget2)

  end if

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
subroutine form ( a, a_int, iatop, iabot, base, iform, jform, nrow, ncol )

!*****************************************************************************80
!
!! FORM converts from one arithmetic form to another.
!
!  Discussion:
!
!    On input, IFORM contains a code for the current arithmetic
!    form, and JFORM contains the code for the new form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(NROW,NCOL), the current real matrix.
!
!    Input/output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input/output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL),
!    the current fractional or decimal matrix.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input/output, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) JFORM, the arithmetic to be converted to.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jform
  real ( kind = 4 ) r

  if ( iform == jform ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  You are already using the arithmetic type that'
    write ( *, '(a)' ) '  you have requested.'
    return
  end if

!
!  Convert the matrix data.
!
  do i = 1, nrow
    do j = 1, ncol
      if ( jform == 0 ) then

        if ( iform == 1 ) then
          call r4_to_rat ( a(i,j), iatop(i,j), iabot(i,j) )
        else if ( iform == 2 ) then
          call dec_to_rat ( iatop(i,j), iabot(i,j) )
        else if ( iform == 3 ) then
          iatop(i,j) = a_int(i,j)
          iabot(i,j) = 1
        else if ( iform == 4 ) then
          iatop(i,j) = a_int(i,j)
          iabot(i,j) = 1
        end if

      else if ( jform == 1 ) then
 
        if ( iform == 0 ) then
          call rat_to_r4 ( iatop(i,j), iabot(i,j), a(i,j) )
        else if ( iform == 2 ) then
          call dec_to_r4 ( iatop(i,j), iabot(i,j), a(i,j) )
        else if ( iform == 3 ) then
          a(i,j) = real ( a_int(i,j) )
        else if ( iform == 4 ) then
          a(i,j) = real ( a_int(i,j) )
        end if
 
      else if ( jform == 2 ) then
 
        if ( iform == 0 ) then
          call rat_to_dec ( iatop(i,j), iabot(i,j), ierror )
        else if ( iform == 1 ) then
          call r4_to_dec ( a(i,j), iatop(i,j), iabot(i,j) ) 
        else if ( iform == 3 ) then
          iatop(i,j) = a_int(i,j)
          iabot(i,j) = 0
        else if ( iform == 4 ) then
          iatop(i,j) = a_int(i,j)
          iabot(i,j) = 0
        end if

      else if ( jform == 3 ) then

        if ( iform == 0 ) then
          call rat_to_r4 ( iatop(i,j), iabot(i,j), r )
          a_int(i,j) = nint ( r )
        else if ( iform == 1 ) then
          a_int(i,j) = nint ( a(i,j) )
        else if ( iform == 2 ) then
          call dec_to_r4 ( iatop(i,j), iabot(i,j), r )
          a_int(i,j) = nint ( r )
        else if ( iform == 4 ) then
          a_int(i,j) = a_int(i,j)
        end if

      else if ( jform == 4 ) then

        if ( iform == 0 ) then
          call rat_to_r4 ( iatop(i,j), iabot(i,j), r )
          a_int(i,j) = mod ( nint ( r ), base )
        else if ( iform == 1 ) then
          a_int(i,j) = mod ( nint ( a(i,j) ), base )
        else if ( iform == 2 ) then
          call dec_to_r4 ( iatop(i,j), iabot(i,j), r )
          a_int(i,j) = mod ( nint ( r ), base )
        else if ( iform == 3 ) then
          a_int(i,j) = mod ( a_int(i,j), base )
        end if

      end if

    end do
  end do
!
!  Update the arithmetic form.
!
  iform = jform
 
  return
end
subroutine help

!*****************************************************************************80
!
!! HELP prints out a brief list of the available commands.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 September 2000
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATMAN commands:'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) 'A(I,J)=S  Set matrix entry to S.'
  write ( *, '(a)' ) 'CHECK     Check matrix for reduced row echelon form.'
  write ( *, '(a)' ) 'COL_APP   Append a column to the matrix.'
  write ( *, '(a)' ) 'COL_AUTO  Automatic column reduction.'
  write ( *, '(a)' ) 'DEC       Use decimal arithmetic.'
  write ( *, '(a)' ) 'DEC_DIGIT Set the number of decimal digits.'
  write ( *, '(a)' ) 'E         Enter matrix with I rows and J columns.'
  write ( *, '(a)' ) 'H         for help.'
  write ( *, '(a)' ) 'I4_BIG    Set size of largest integer for fractions.'
  write ( *, '(a)' ) 'ID        Append the identity matrix.'
  write ( *, '(a)' ) 'INIT      Initialize data.'
  write ( *, '(a)' ) 'INT       Use integer arithmetic.'
  write ( *, '(a)' ) 'Q         Quit.'
  write ( *, '(a)' ) 'RAN       Set up a problem with random data.'
  write ( *, '(a)' ) 'RAT       Use rational arithmetic.'
  write ( *, '(a)' ) 'REAL      Use real arithmetic.'
  write ( *, '(a)' ) 'ROW_AUTO  Automatic row reduction.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ERO''s and ECO''s:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R1 <=>  R2          C1 <=>   C2'
  write ( *, '(a)' ) 'R1 <= S R1          C1 <=  S C1'
  write ( *, '(a)' ) 'R1 <=   R1 + S R2   C1 <=    C1 + S C2'

  return
end
subroutine i4_data ( op, var, ival )

!*****************************************************************************80
!
!! I4_DATA stores and retrieves common data items.
!
!  Discussion:
!
!    This routine works like a sort of COMMON block.  It stores or returns
!    the values of certain variables.  Thus, it allows routines
!    to "communicate" without having to have data passed up and
!    down the calling tree in argument lists.
!
!    The variables stored by this version of the routine are:
!
!    'DEC_DIGIT', the number of digits stored for decimals;
!    'DEC_DIGIT_MAX', the maximum number of digits stored for decimals;
!    'I4_BIG', the biggest integer to use in calculations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OP, describes the operation to be done.
!    'SET' means set a value.
!    'INC' means increment a value (and return its new value)
!    'GET' means get a value.
!
!    Input, character ( len = * ) VAR, the name of the variable.
!
!    Input/output, integer ( kind = 4 ) IVAL.
!    If OP is 'SET', then the variable named in VAR is set to the
!    value IVAL.
!    If OP is 'GET', then the value of IVAL is set to the value of
!    the variable named in VAR.
!    If OP is 'INC', then the value of IVAL is incremented by 1,
!    and its new value is returned in VAR.
!
  implicit none

  integer, save :: dec_digit = 7
  integer, save :: dec_digit_max = 7
  integer, save :: i4_big = 2147483647
  integer ( kind = 4 ) ival
  character ( len = * ) op
  logical s_eqi
  character ( len = * ) var

  if ( s_eqi ( op, 'SET' ) ) then

    if ( s_eqi ( var, 'DEC_DIGIT' ) ) then
      dec_digit = ival
    else if ( s_eqi ( var, 'DEC_DIGIT_MAX' ) ) then
      dec_digit_max = ival
    else if ( s_eqi ( var, 'I4_BIG' ) ) then
      i4_big = ival
    end if

  else if ( s_eqi ( op, 'GET' ) ) then

    if ( s_eqi ( var, 'DEC_DIGIT' ) ) then
      ival = dec_digit
    else if ( s_eqi ( var, 'DEC_DIGIT_MAX' ) ) then
      ival = dec_digit_max
    else if ( s_eqi ( var, 'I4_BIG' ) ) then
      ival = i4_big
    end if

  else if ( s_eqi ( op, 'INC' ) ) then

    if ( s_eqi ( var, 'DEC_DIGIT' ) ) then
      dec_digit = dec_digit + 1
      ival = dec_digit
    else if ( s_eqi ( var, 'DEC_DIGIT_MAX' ) ) then
      dec_digit_max = dec_digit_max + 1
      ival = dec_digit_max
    else if ( s_eqi ( var, 'I4_BIG' ) ) then
      i4_big = i4_big + 1
      ival = i4_big
    end if

  end if
 
  return
end
function i4_gcd ( i, j )

!*****************************************************************************80
!
!! I4_GCD finds the greatest common divisor of I and J.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, two numbers whose greatest common divisor
!    is desired.
!
!    Output, integer ( kind = 4 ) I4_GCD, the greatest common divisor of I and J.
!
!    Note that only the absolute values of I and J are
!    considered, so that the result is always nonnegative.
!
!    If I or J is 0, I4_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
!
!    If I and J have no common factor, I4_GCD is returned as 1.
!
!    Otherwise, using the Euclidean algorithm, I4_GCD is the
!    largest common factor of I and J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j

  i4_gcd = 1
!
!  Return immediately if either I or J is zero.
!
  if ( i == 0 ) then
    i4_gcd = max ( 1, abs ( j ) )
    return
  else if ( j == 0 ) then
    i4_gcd = max ( 1, abs ( i ) )
    return
  end if
!
!  Set IP to the larger of I and J, IQ to the smaller.
!  This way, we can alter IP and IQ as we go.
!
  ip = max ( abs ( i ), abs ( j ) )
  iq = min ( abs ( i ), abs ( j ) )
!
!  Carry out the Euclidean algorithm.
!
  do

    ir = mod ( ip, iq )

    if ( ir == 0 ) then
      exit
    end if

    ip = iq
    iq = ir

  end do

  i4_gcd = iq

  return
end
subroutine i4_print ( a_int, nrow, ncol, title )

!*****************************************************************************80
!
!! I4_PRINT prints out integer vectors and matrices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
!    Input, character ( len = * ) TITLE, a label for the object being printed.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ) a_int(nrow,ncol)
  logical allint
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 10 ) chrtmp3
  character ( len = 40 ) format1
  character ( len = 40 ) format2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ichi
  integer ( kind = 4 ) iclo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) izhi
  integer ( kind = 4 ) izlo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) jmin
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) kmin
  integer, parameter :: ncolum = 80
  integer ( kind = 4 ) npline
  character ( len = 100 ) output
  character ( len = * ) title
!
!  Figure out how many numbers we can fit in NCOLUM columns.
!
  kmin = 11
  kmax = kmin

  do i = 1, nrow
    do j = 1, ncol

      do while ( 10.0E+00**(kmax-kmin) <= abs ( a_int(i,j) ) )
        kmax = kmax + 1
      end do

    end do
  end do

  npline = ncolum / kmax
!
!  If all integers, cut down KMAX, the width of each number,
!  and update NPLINE, the number of numbers we can print on one line.
!
    kmax = kmax - 7
    npline = ncolum / kmax

    call i4_to_s_left ( npline, chrtmp2 )
    call i4_to_s_left ( kmax, chrtmp3 )

    format1 = '(' // chrtmp2 // 'i' // chrtmp3 // ')'

  call s_blank_delete ( format1 )
!
!  The second format is for ...
!
    format2 = '(' // chrtmp2 // 'i' // chrtmp3 // ')'

  call s_blank_delete ( format2 )
!
!  Now print the data.
!
  do jmin = 1, ncol, npline

    jmax = min ( jmin + npline - 1, ncol )

    write ( *, '(a)' ) ' '

    if ( jmin == 1 ) then
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
    end if

    if ( 1 < jmin .or. jmax < ncol ) then
      write ( output, * ) 'Columns ', jmin, ' to ', jmax
      call s_blanks_delete ( output )
      write ( *, '(a)' ) trim ( output )
      write ( *, '(a)' ) ' '
    end if

    do i = 1, nrow
      write ( *, format1 ) a_int(i,jmin:jmax)
    end do

  end do

  return
end
subroutine i4_random ( ilo, ihi, i )

!*****************************************************************************80
!
!! I4_RANDOM returns a random integer in a given range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ILO, IHI, the minimum and maximum acceptable values.
!
!    Output, integer ( kind = 4 ) I, the randomly chosen integer.
!
  implicit none

  logical, save :: seed = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  real ( kind = 4 ) r
  real ( kind = 4 ) rhi
  real ( kind = 4 ) rlo

  if ( .not. seed ) then
    call random_seed ( )
    seed = .true.
  end if
!
!  Pick a random number in (0,1).
!
  call random_number ( harvest = r )
!
!  Set a real interval [RLO,RHI] which contains the integers [ILO,IHI],
!  each with a "neighborhood" of width 1.
!
  rlo = real ( ilo ) - 0.5E+00
  rhi = real ( ihi ) + 0.5E+00
!
!  Set I to the integer that is nearest the scaled value of R.
!
  i = nint ( ( 1.0E+00 - r ) * rlo + r * rhi )
!
!  In case of oddball events at the boundary, enforce the limits.
!
  i = max ( i, ilo )
  i = min ( i, ihi )

  return
end
subroutine i4_read ( intval, line, prompt, ierror )

!*****************************************************************************80
!
!! I4_READ reads an integer from the input buffer.
!
!  Discussion:
!
!    The routine accepts LINE which contains input and a PROMPT line.  
!    If LINE is empy, the PROMPT will be printed and LINE read from the 
!    In either case, the integer INTVAL will be read from LINE,
!    beginning at character 1 and ending at the first comma, slash,
!    blank, or the end of LINE.
!
!    The PROMPT should consist of a string of names of data items,
!    separated by commas, with the current one first.
!
!    The program will print 'ENTER' PROMPT and after reading LINE
!    will strip the characters corresponding to INTVAL from LINE,
!    and the part of PROMPT up to the first comma, leaving LINE and
!    PROMPT ready for another call to I4_READ, S_READ, RAT_READ or
!    R4_READ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) INTVAL, the integer that was read.
!
!    Workspace, character ( len = 80 ) LINE, used to hold the user's input.
!
!    Input, character ( len = 80 ) PROMPT, the prompt string.
!
!    Output, integer ( kind = 4 ) IERROR, 0 an error occurred, 1 no error occurred.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) lchar
  character ( len = 80 ) line
  character ( len = 80 ) prompt

  ierror = 0
  intval = 0
!
!  Retrieve a likely character string from input.
!
  do

    call chrinp ( ierror, line, prompt )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I_READ - Fatal error!'
      write ( *, '(a)' ) '  CHRINP returned error flag!'
      return
    end if

    if ( line /= ' ' ) then
      exit
    end if

  end do
!
!  Convert the character string to an integer.
!
  call s_to_i4 ( line, intval, ierror, lchar )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I_READ - Fatal error!'
    write ( *, '(a)' ) '  S_TO_I4 returned error flag!'
    
    return
  end if
!
!  Remove the character string from the input line.
!
  if ( lchar < len ( line ) ) then
    if ( line(lchar+1:lchar+1) == ',' ) then
      lchar = lchar + 1
    end if
  end if

  call s_chop ( line, 1, lchar )
 
  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP switches two I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
subroutine i4_to_s_left ( intval, s )

!*****************************************************************************80
!
!! I4_TO_S_LEFT converts an I4 to a left-justified string.
!
!  Example:
!
!    Assume that S is 6 characters long:
!
!    INTVAL  S
!
!         1  1
!        -1  -1
!         0  0
!      1952  1952
!    123456  123456
!   1234567  ******  <-- Not enough room!
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
!    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be left-justified.  If there is not enough space,
!    the string will be filled with stars.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) ival
  character ( len = * ) s

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  ival = intval
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = - ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  The absolute value of the integer goes into S(ILO:IHI).
!
  ipos = ihi
!
!  Find the last digit of IVAL, strip it off, and stick it into the string.
!
  do

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do i = 1, ihi
        s(i:i) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )

    s(ipos:ipos) = c
    ipos = ipos - 1

    if ( ival == 0 ) then
      exit
    end if

  end do
!
!  Shift the string to the left.
!
  s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
  s(ilo+ihi-ipos:ihi) = ' '

  return
end
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real    ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end
subroutine init ( a, a_int, base, iabot, iatop, ierror, iform, &
  line, maxrow, maxcol, nrow, ncol )

!*****************************************************************************80
!
!! INIT initializes the data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 4 ) A(MAXROW,MAXCOL), the current matrix.
!
!    Input/output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Output, integer ( kind = 4 ) IABOT(MAXROW,MAXCOL), IATOP(MAXROW,MAXCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Output, integer ( kind = 4 ) IERROR.
!    The error flag, which is initialized to zero by this routine.
!
!    Output, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Output, character ( len = 80 ) LINE,
!    a buffer used to hold the user's input.
!
!    Input, integer ( kind = 4 ) MAXCOL, the maximum number of matrix columns.
!
!    Input, integer ( kind = 4 ) MAXROW, the maximum number of matrix rows.
!
!    Output, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
  implicit none

  integer ( kind = 4 ) maxcol
  integer ( kind = 4 ) maxrow

  real ( kind = 4 ) a(maxrow,maxcol)
  integer ( kind = 4 ) a_int(maxrow,maxcol)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) dec_digit
  integer ( kind = 4 ) dec_digit_max
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_big
  integer ( kind = 4 ) iabot(maxrow,maxcol)
  integer ( kind = 4 ) iatop(maxrow,maxcol)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  character ( len = * ) line
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  a(1,1) = 1.0E+00
  a_int(1,1) = 1
  iatop(1,1) = 1
  iabot(1,1) = 1
 
  base = 0
  ierror = 0
  iform = 0
  line = ' '

  dec_digit_max = 7
  call i4_data ( 'SET', 'DEC_DIGIT_MAX', dec_digit_max )

  dec_digit = 4
  call i4_data ( 'SET', 'DEC_DIGIT', dec_digit )

  i4_big = huge ( i4_big )
  call i4_data ( 'SET', 'I4_BIG', i4_big )

  ncol = 0
  nrow = 0
 
  return
end
subroutine la_inp1 ( nrow, ncol, a, a_int, iabot, iatop, ierror, base, iform, &
  line, row1, row2, col1, col2 )

!*****************************************************************************80
!
!! LA_INP1 accepts the values of the entries of a matrix from the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
!    Output, real ( kind = 4 ) A(NROW,NCOL), the current matrix.
!
!    Output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Output, integer ( kind = 4 ) IERROR, 0 an error occurred, 1 no error occurred.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Workspace, character ( len = 80 ) LINE, used to hold the user's input.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  integer ( kind = 4 ) base
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 10 ) chrtmp3
  integer ( kind = 4 ) col1
  integer ( kind = 4 ) col2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) j
  character ( len = 80 ) line
  character ( len = 80 ) prompt
  integer ( kind = 4 ) row1
  integer ( kind = 4 ) row2
!
!  Enter a single value
!
  if ( row1 == row2 .and. col1 == col2 ) then
 
    call i4_to_s_left ( row1, chrtmp1 )
    call i4_to_s_left ( col1, chrtmp2 )

    prompt = 'A(' // trim ( chrtmp1 ) // ',' // trim ( chrtmp2 ) //')'
 
    if ( iform == 0 ) then
      call rat_read ( iatop(row1,col1), iabot(row1,col1), line, prompt, ierror )
    else if ( iform == 1 ) then
      call r4_read ( a(row1,col1), line, prompt, ierror )
    else if ( iform == 2 ) then
      call dec_read ( iatop(row1,col1), iabot(row1,col1), line, prompt, ierror )
      call dec_round ( iatop(row1,col1), iabot(row1,col1) )
    else if ( iform == 3 ) then
      call i4_read ( a_int(row1,col1), line, prompt, ierror )
    else if ( iform == 4 ) then
      call i4_read ( a_int(row1,col1), line, prompt, ierror )
      a_int(row1,col1) = mod ( a_int(row1,col1), base )
    end if

    if ( ierror /= 0 ) then
      return
    end if
!
!  Enter a single column.
!
  else if ( col1 == col2 ) then
 
    do i = row1, row2
 
      call i4_to_s_left ( i, chrtmp1 )
      call i4_to_s_left ( row2, chrtmp2 )
      call i4_to_s_left ( col1, chrtmp3 )
      prompt = 'entries ' // chrtmp1 // ' to ' // chrtmp2 // ' of column ' // &
        chrtmp3
      call s_blanks_delete ( prompt )
 
      if ( iform == 0 ) then
 
        call rat_read ( iatop(i,col1), iabot(i,col1), line, prompt, ierror )
 
      else if ( iform == 1 ) then
 
        call r4_read ( a(i,col1), line, prompt, ierror )
 
      else if ( iform == 2 ) then
 
        call dec_read ( iatop(i,col1), iabot(i,col1), line, prompt, ierror )
        call dec_round ( iatop(i,col1), iabot(i,col1) )
 
      else if ( iform == 3 ) then

        call i4_read ( a_int(i,col1), line, prompt, ierror )
      else if ( iform == 4 ) then
        call i4_read ( a_int(i,col1), line, prompt, ierror )
        a_int(i,col1) = mod ( a_int(i,col1), base )
      end if

      if ( ierror /= 0 ) then
        return
      end if
 
    end do
!
!  Enter one or more partial rows of at least 2 entries.
!
  else
 
    do i = row1, row2
      do j = col1, col2
 
        call i4_to_s_left ( j, chrtmp1 )
        call i4_to_s_left ( col2, chrtmp2 )
        call i4_to_s_left ( i, chrtmp3 )
        prompt = 'entries ' // chrtmp1 // ' to ' // chrtmp2 // ' of row ' // &
          chrtmp3
        call s_blanks_delete ( prompt )
 
        if ( iform == 0 ) then
 
          call rat_read ( iatop(i,j), iabot(i,j), line, prompt, ierror )
 
        else if ( iform == 1 ) then
 
          call r4_read ( a(i,j), line, prompt, ierror )
 
        else if ( iform == 2 ) then
 
          call dec_read ( iatop(i,j), iabot(i,j), line, prompt, ierror )
          call dec_round ( iatop(i,j), iabot(i,j) )

        else if ( iform == 3 ) then

          call i4_read ( a_int(i,j), line, prompt, ierror )

        else if ( iform == 3 ) then

          call i4_read ( a_int(i,j), line, prompt, ierror )
          a_int(i,j) = mod ( a_int(i,j), base )
 
        end if

        if ( ierror /= 0 ) then
          return
        end if
 
      end do
    end do
 
  end if
 
  return
end
subroutine la_opt ( a, a_int, iabot, iatop, base, iform, nrow, ncol )

!*****************************************************************************80
!
!! LA_OPT checks for row echelon or reduced row echelon form.
!
!  Discussion:
!
!    A matrix is in row echelon form if:
!
!    * The first nonzero entry in each row is 1.
!
!    * The leading 1 in a given row occurs in a column to
!      the right of the leading 1 in the previous row.
!
!    * Rows which are entirely zero must occur last.
!
!    The matrix is in reduced row echelon form if, in addition to
!    the first three conditions, it also satisfies:
!
!    * Each column containing a leading 1 has no other nonzero entries.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(NROW,NCOL), the current matrix.
!
!    Input, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  integer ( kind = 4 ) base
  character ( len = 22 ) chrtmp
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) izer
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lead
  integer ( kind = 4 ) leadp

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Checking the matrix for row echelon form...'
!
!  Check rule 1.
!
  do i = 1, nrow
    do j = 1, ncol
 
      if ( iform == 0 ) then
 
        if ( iatop(i,j) == 0 ) then
          cycle
        else if ( iatop(i,j) == iabot(i,j) ) then
          exit
        end if
 
      else if ( iform == 1 ) then
 
        if ( a(i,j) == 0.0E+00 ) then
          cycle
        else if ( a(i,j) == 1.0E+00 ) then
          exit
        end if
 
      else if ( iform == 2 ) then
 
        if ( iatop(i,j) == 0 ) then
          cycle
        else if ( iatop(i,j) == 1 .and. iabot(i,j) == 0 ) then
          exit
        end if

      else if ( iform == 3 ) then
 
        if ( a_int(i,j) == 0 ) then
          cycle
        else if ( a_int(i,j) == 1 ) then
          exit
        end if

      else if ( iform == 4 ) then
 
        if ( a_int(i,j) == 0 ) then
          cycle
        else if ( a_int(i,j) == 1 ) then
          exit
        end if

      end if
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix is NOT in row echelon form.'
      
      call i4_to_s_left ( i, chrtmp1 )
      write ( *, '(a)' ) '  The first nonzero entry in row ' // chrtmp1
      
      call i4_to_s_left ( j, chrtmp1 )
      write ( *, '(a)' ) '  which occurs in column ' // chrtmp1
      
      if ( iform == 0 ) then
        call rat_to_s_left ( iatop(i,j), iabot(i,j), chrtmp )
      else if ( iform == 1 ) then
        call r4_to_s_left ( a(i,j), chrtmp )
      else if ( iform == 2 ) then
        call dec_to_s_left ( iatop(i,j), iabot(i,j), chrtmp )
      else if ( iform == 3 ) then
        call i4_to_s_left ( a_int(i,j), chrtmp )
      else if ( iform == 4 ) then
        call i4_to_s_left ( a_int(i,j), chrtmp )
      end if
 
      write ( *, '(a)' ) '  is ' // trim ( chrtmp ) // ' rather than 1.'
      
      return
 
    end do
 
  end do
!
!  Check rule 2.
!
  lead = 0
 
  do i = 1, nrow
    do j = 1, ncol
 
      if ( iform == 0 ) then
 
        if ( iatop(i,j) == 0 ) then
          cycle
        else if ( iatop(i,j) == iabot(i,j) ) then
          leadp = lead
          lead = j
          if ( leadp < lead ) then
            exit
          end if
        end if
 
      else if ( iform == 1 ) then
 
        if ( a(i,j) == 0.0E+00 ) then
          cycle
        else if ( a(i,j) == 1.0E+00 ) then
          leadp = lead
          lead = j
          if ( leadp < lead ) then
            exit
          end if
        end if
 
      else if ( iform == 2 ) then
 
        if ( iatop(i,j) == 0 ) then
          cycle
        else if ( iatop(i,j) == 1 .and. iabot(i,j) == 0 ) then
          leadp = lead
          lead = j
          if ( leadp < lead ) then
            exit
          end if
        end if
 
      else if ( iform == 3 ) then
 
        if ( a_int(i,j) == 0 ) then
          cycle
        else if ( a_int(i,j) == 1 ) then
          leadp = lead
          lead = j
          if ( leadp < lead ) then
            exit
          end if
        end if

      else if ( iform == 4 ) then
 
        if ( a_int(i,j) == 0 ) then
          cycle
        else if ( a_int(i,j) == 1 ) then
          leadp = lead
          lead = j
          if ( leadp < lead ) then
            exit
          end if
        end if

      end if
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix is NOT in row echelon form.'
      
      call i4_to_s_left ( i, chrtmp1 )
      write ( *, '(a)' ) '  The first 1 in row ' // trim ( chrtmp1 ) // ' does'
      
      call i4_to_s_left ( i-1, chrtmp1 )
      write ( *, '(a)' ) '  NOT occur to the right of the first 1 in row ' // &
        trim ( chrtmp1 ) // '.'
      
      return
 
    end do
 
  end do
!
!  Check rule 3.
!
  izer = 0
 
  do i = 1, nrow
 
    if ( izer == 0 ) then
 
      izer = i

      do j = 1, ncol
 
        if ( ( iform == 0 .and. iatop(i,j) /= 0 ) .or. &
             ( iform == 1 .and. a(i,j) /= 0.0E+00 ) .or. &
             ( iform == 2 .and. iatop(i,j) /= 0 ) .or. &
             ( iform == 3 .and. a_int(i,j) /= 0 ) .or. &
             ( iform == 4 .and. a_int(i,j) /= 0 ) ) then
          izer = 0
          exit
        end if
 
      end do
 
    else
 
      do j = 1, ncol
 
        if ( ( iform == 0 .and. iatop(i,j) /= 0 ) .or. &
             ( iform == 1 .and. a(i,j) /= 0.0E+00 ) .or. &
             ( iform == 2 .and. iatop(i,j) /= 0 ) .or. &
             ( iform == 3 .and. a_int(i,j) /= 0 ) .or. &
             ( iform == 4 .and. a_int(i,j) /= 0 ) ) then
 
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  The matrix is NOT in row echelon form.'
          
          call i4_to_s_left ( izer, chrtmp1 )
          write ( *, '(a)' ) '  Row ' // trim ( chrtmp1 ) // ' is entirely zero.'
          call i4_to_s_left ( i, chrtmp1 )
          write ( *, '(a)' ) '  Row ' // trim ( chrtmp1 ) // ' occurs later, and has'
          write ( *, '(a)' ) '  nonzero entries in it!'
          
          return

        end if
 
      end do
 
    end if
 
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The matrix is in row echelon form.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now checking to see if the matrix is in '
  write ( *, '(a)' ) '  REDUCED row echelon form...'
!
!  Check rule 4.
!
  do i = 1, nrow
 
    do j = 1, ncol
!
!  We know first nonzero in this row will be 1.
!
      if ( ( iform == 0 .and. iatop(i,j) == 0 ) .or. &
           ( iform == 1 .and. a(i,j) == 0.0E+00 ) .or. &
           ( iform == 2 .and. iatop(i,j) == 0 ) .or. &
           ( iform == 3 .and. a_int(i,j) == 0 ) .or. &
           ( iform == 4 .and. a_int(i,j) == 0 ) ) then
        cycle
      end if
!
!  The leading 1 of this row is entry (i,j).
!
      do ii = 1, nrow
 
        if ( ii /= i ) then
 
          if ( ( iform == 0 .and. iatop(ii,j) == 0 ) .or. &
               ( iform == 1 .and. a(ii,j) == 0.0E+00 ) .or. &
               ( iform == 2 .and. iatop(ii,j) == 0 ) .or. &
               ( iform == 3 .and. a_int(ii,j) == 0 ) .or. &
               ( iform == 4 .and. a_int(ii,j) == 0 ) ) then
            cycle
          end if
 
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  The matrix is NOT in reduced row echelon form.'
          
          call i4_to_s_left ( i, chrtmp1 )

          call i4_to_s_left ( j, chrtmp2 )

          write ( *, '(a)' ) '  Row ' // trim ( chrtmp1 ) //  &
            ' has its leading 1 in column ' // trim ( chrtmp2 ) // '.'
          write ( *, '(a)' ) '  This means that all other entries of that ' // &
            'column should be zero.'
 
          if ( iform == 0 ) then
            call rat_to_s_left ( iatop(ii,j), iabot(ii,j), chrtmp )
          else if ( iform == 1 ) then
            call r4_to_s_left ( a(ii,j), chrtmp )
          else if ( iform == 2 ) then
            call dec_to_s_left ( iatop(ii,j), iabot(ii,j), chrtmp )
          else if ( iform == 3 ) then
            call i4_to_s_left ( a_int(ii,j), chrtmp )
          else if ( iform == 4 ) then
            call i4_to_s_left ( a_int(ii,j), chrtmp )
          end if

          call i4_to_s_left ( ii, chrtmp1 )

          write ( *, '(a)' ) '  But the entry in row ' // trim ( chrtmp1 ) &
           // ' is ' // trim ( chrtmp )
          
          return

        end if
 
      end do
 
      exit
 
    end do

  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The matrix is in reduced row echelon form.'
 
  return
end
subroutine mat_append_identity ( nrow, ncol, a, a_int, iatop, iabot, iform )

!*****************************************************************************80
!
!! MAT_APPEND_IDENTITY appends the identity matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns of
!    the matrix.  On output, NCOL has been increased to NCOL+NROW.
!
!    Input/output, real ( kind = 4 ) A(NROW,NCOL), the current matrix.
!
!    Input/output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input/output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
  implicit none

  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,*)
  integer ( kind = 4 ) a_int(nrow,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iabot(nrow,*)
  integer ( kind = 4 ) iatop(nrow,*)
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jform
  integer ( kind = 4 ) ncol

  do i = 1, nrow
    do j = 1, nrow
      if ( i == j ) then
        if ( iform == 0 ) then
          iatop(i,ncol+j) = 1
          iabot(i,ncol+j) = 1
        else if ( iform == 1 ) then
          a(i,ncol+j) = 1.0E+00
        else if ( iform == 2 ) then
          iatop(i,ncol+j) = 1
          iabot(i,ncol+j) = 0
        else if ( iform == 3 ) then
          a_int(i,ncol+j) = 1
        end if
      else
        if ( iform == 0 ) then
          iatop(i,ncol+j) = 0
          iabot(i,ncol+j) = 1
        else if ( iform == 1 ) then
          a(i,ncol+j) = 0.0E+00
        else if ( iform == 2 ) then
          iatop(i,ncol+j) = 0
          iabot(i,ncol+j) = 0
        else if ( iform == 3 ) then
          a_int(i,ncol+j) = 0
        end if
      end if
    end do
  end do

  ncol = ncol + nrow

  return
end
subroutine mat_print ( a, a_int, iabot, iatop, iform, nrow, ncol, title )

!*****************************************************************************80
!
!! MAT_PRINT prints out the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(NROW,NCOL), the current matrix.
!
!    Input, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input, integer ( kind = 4 ) IABOT(NROW,NCOL), IATOP(NROW,NCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
!    Input, character ( len = * ) TITLE, the title of the object to be printed.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) iform
  character ( len = * ) title

  if ( iform == 0 ) then
    call rat_print ( iatop, iabot, nrow, ncol, title )
  else if ( iform == 1 ) then
    call r4_print ( a, nrow, ncol, title )
  else if ( iform == 2 ) then
    call dec_print ( iatop, iabot, nrow, ncol, title )
  else if ( iform == 3 ) then
    call i4_print ( a_int, nrow, ncol, title )
  else if ( iform == 4 ) then
    call i4_print ( a_int, nrow, ncol, title )
  end if
 
  return
end
subroutine mat_random ( nrow, ncol, base, iform, seed, a, a_int, iatop, iabot )

!*****************************************************************************80
!
!! MAT_RANDOM sets up a random problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 4 ) A(NROW,NCOL), the current matrix.
!
!    Output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  if ( iform /= 4 ) then
    ilo = -10
    ihi = 10
  else
    ilo = 0
    ihi = base - 1
  end if

  do i = 1, nrow
    do j = 1, ncol

      ival = i4_uniform ( ilo, ihi, seed )

      if ( iform == 0 ) then
        iatop(i,j) = ival
        iabot(i,j) = 1
      else if ( iform == 1 ) then
        a(i,j) = real ( ival )
      else if ( iform == 2 ) then
        iatop(i,j) = ival
        iabot(i,j) = 0
      else if ( iform == 3 ) then
        a_int(i,j) = ival
      else if ( iform == 4 ) then
        a_int(i,j) = ival
      end if

    end do
  end do

  return
end
subroutine problem_solve_random ( a, a_int, iatop, iabot, base, iform, &
  nrow, ncol, seed )

!*****************************************************************************80
!
!! PROBLEM_SOLVE_RANDOM sets up a random linear system problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 4 ) A(NROW,NCOL), the current matrix.
!
!    Output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ival2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jform
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) sol(ncol-1)

  ilo = -10
  ihi = 10

  do i = 1, nrow
    do j = 1, ncol-1
      ival = i4_uniform ( ilo, ihi, seed )
      a_int(i,j) = ival
    end do
  end do

  do j = 1, ncol-1
    sol(j) = i4_uniform ( ilo, ihi, seed )
  end do

  do i = 1, nrow

    ival = 0
    do j = 1, ncol-1
      ival = ival + iatop(i,j) * sol(j)
    end do

    iatop(i,ncol) = ival
 
  end do

  jform = iform
  iform = 3

  call form ( a, a_int, iatop, iabot, base, iform, jform, nrow, ncol )

  return
end
subroutine r4_print ( a, nrow, ncol, title )

!*****************************************************************************80
!
!! R4_PRINT prints out real vectors and matrices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(NROW,NCOL), the current matrix.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
!    Input, character ( len = * ) TITLE, a label for the object being printed.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  logical allint
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 10 ) chrtmp3
  character ( len = 40 ) format1
  character ( len = 40 ) format2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ichi
  integer ( kind = 4 ) iclo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) izhi
  integer ( kind = 4 ) izlo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) jmin
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) kmin
  integer, parameter :: ncolum = 80
  integer ( kind = 4 ) npline
  character ( len = 100 ) output
  character ( len = * ) title
!
!  Figure out how many numbers we can fit in NCOLUM columns.
!
  kmin = 11
  kmax = kmin
 
  do i = 1, nrow
    do j = 1, ncol 
  
      do while ( 10.0E+00**(kmax-kmin) < abs ( a(i,j) ) )
        kmax = kmax + 1
      end do
 
    end do
  end do
 
  npline = ncolum / kmax
!
!  Check to see if the matrix entries are all integers.
!
  allint = .true.
 
  do i = 1, nrow
    do j = 1, ncol
 
      if ( a(i,j) /= real ( int ( a(i,j) ) ) ) then
        allint = .false.
      end if

    end do
  end do
!
!  If all integers, cut down KMAX, the width of each number,
!  and update NPLINE, the number of numbers we can print on one line.
!
  if ( allint ) then

    kmax = kmax - 7
    npline = ncolum / kmax
 
    call i4_to_s_left ( npline, chrtmp2 )
    call i4_to_s_left ( kmax, chrtmp3 )

    format1 = '(' // chrtmp2 // 'f' // chrtmp3 // '.0)'
!
!  If nonintegral entries, print 7 decimals.
!
  else
 
    call i4_to_s_left ( npline, chrtmp2 )
    call i4_to_s_left ( kmax, chrtmp3 )

    format1 = '(' // chrtmp2 // 'f' // chrtmp3 // '.7)'
 
  end if
 
  call s_blank_delete ( format1 )
!
!  The second format is for ...
!
    format2 = '(' // chrtmp2 // 'i' // chrtmp3 // ')'

  call s_blank_delete ( format2 )
!
!  Now print the data.
!
  do jmin = 1, ncol, npline
 
    jmax = min ( jmin + npline - 1, ncol )
 
    write ( *, '(a)' ) ' '
 
    if ( jmin == 1 ) then
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
    end if

    if ( 1 < jmin .or. jmax < ncol ) then
      write ( output, * ) 'Columns ', jmin, ' to ', jmax
      call s_blanks_delete ( output )
      write ( *, '(a)' ) trim ( output )
      write ( *, '(a)' ) ' '    
    end if
 
    do i = 1, nrow
      write ( *, format1 ) a(i,jmin:jmax)
    end do
 
  end do
 
  return
end
subroutine r4_read ( rval, line, prompt, ierror )

!*****************************************************************************80
!
!! R4_READ "reads" a real value from a line of text.
!
!  Discussion:
!
!    The routine accepts a line of characters which may contain some
!    user input.  If not, it prints out the PROMPT and reads new
!    information into LINE, seeking to find a real number RVAL to
!    return.
!
!    The routine will accept integers, decimals, and ratios of the
!    form R1/R2.  Real numbers may be in scientific notation, as
!    +12.34E-56.78
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 4 ) RVAL, the real value found in LINE.
!
!    Workspace, character ( len = 80 ) LINE, used to hold the user's input.
!
!    Workspace, character ( len = 80 ) PROMPT, the prompt string.
!
!    Output, integer ( kind = 4 ) IERROR, 0 an error occurred, 1 no error occurred.
!
  implicit none

  real ( kind = 4 ) bot
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) lchar
  character ( len = 80 ) line
  character ( len = 80 ) prompt
  real ( kind = 4 ) rval
  real ( kind = 4 ) top

  ierror = 0
  rval = 0.0E+00
  top = 0.0E+00
  bot = 1.0E+00
!
!  Read a character string.
!
  do
 
    call chrinp ( ierror, line, prompt )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R_READ - Fatal error!'
      write ( *, '(a)' ) '  CHRINP returned error flag!'
      return
    end if

    if ( line /= ' ' ) then
      exit
    end if

  end do
!
!  Convert the character string to a decimal value, TOP.
!
  call s_to_r4 ( line, top, ierror, lchar )
!
!  If we haven't used up all our characters,
!  and if the next character is '/',
!  then the user means to input the value as a ratio,
!  so prepare to read BOT as well.
!
  if ( lchar+1 < len ( line ) ) then
    if ( line(lchar+1:lchar+1) == '/' ) then
      lchar = lchar + 1
      call s_chop ( line, 1, lchar )
      call s_to_r4 ( line, bot, ierror, lchar )
      if ( bot == 0.0E+00 ) then
        bot = 1.0E+00
      end if
    end if
  end if
!
!  Set the value of RVAL.
!
  rval = top / bot
!
!  Chop out the characters that were used.
!
  if ( lchar < len ( line ) ) then
    if ( line(lchar+1:lchar+1) == ',' ) then
      lchar = lchar + 1
    end if
  end if

  call s_chop ( line, 1, lchar )
 
  return
end
subroutine r4_swap ( x, y )

!*****************************************************************************80
!
!! R4_SWAP switches two real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 4 ) z
  real ( kind = 4 ) x
  real ( kind = 4 ) y

  z = x
  x = y
  y = z

  return
end
subroutine r4_to_dec ( rval, itop, ibot )

!*****************************************************************************80
!
!! R4_TO_DEC converts a real value to a decimal fraction form.
!
!  Discussion:
!
!    The routine is given RVAL, and computes ITOP and IBOT, so that 
!    approximately:
!
!      RVAL = ITOP * 10 ** IBOT
!
!    However, only DEC_DIGIT digits of RVAL are used in constructing the 
!    representation, where DEC_DIGIT is the maximum number of decimal
!    digits used in the decimal representation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) RVAL, the real number to be converted.
!
!    Output, integer ( kind = 4 ) ITOP, IBOT, the approximate decimal representation.
!
!    ITOP is an integer, strictly between -10**DEC_DIGIT and 10**DEC_DIGIT.
!    IBOT is an integer exponent of 10.
!
  implicit none

  integer ( kind = 4 ) dec_digit
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) itop
  real ( kind = 4 ) rtop
  real ( kind = 4 ) rval
  real ( kind = 4 ) ten1
  real ( kind = 4 ) ten2
!
!  Special cases.
!
  if ( rval == 0.0E+00 ) then
    itop = 0
    ibot = 0
    return
  end if

  dec_digit = 0
  call i4_data ( 'GET', 'DEC_DIGIT', dec_digit )
!
!  Factor RVAL = RTOP * 10**IBOT
!
  rtop = rval
  ibot = 0
!
!  Now normalize so that 10**(DEC_DIGIT-1) <= ABS(RTOP) < 10**(DEC_DIGIT)
!
  ten1 = 10.0E+00**( dec_digit - 1 )
  ten2 = 10.0E+00**dec_digit
  
  do while ( abs ( rtop ) < ten1 )
    rtop = rtop * 10.0E+00
    ibot = ibot - 1
  end do

  do while ( ten2 <= abs ( rtop ) )
    rtop = rtop / 10.0E+00
    ibot = ibot + 1
  end do
!
!  ITOP is the integer part of RTOP, rounded.
!
  itop = nint ( rtop )
!
!  Now divide out any factors of ten from ITOP.
!
  if ( itop /= 0 ) then

    do while ( mod ( itop, 10 ) == 0 )
      itop = itop / 10
      ibot = ibot + 1
    end do

  end if
 
  return
end
subroutine r4_to_rat ( a, iatop, iabot )

!*****************************************************************************80
!
!! R4_TO_RAT converts a real value to a rational value.  
!
!  Discussion:
!
!    The rational value (IATOP/IABOT) is essentially computed by truncating 
!    the decimal representation of the real value after a given number of
!    decimal digits.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, the real value to be converted.
!
!    Output, integer ( kind = 4 ) IATOP, IABOT, the numerator and denominator
!    of the rational value that approximates A.
!
  implicit none

  real ( kind = 4 ) a
  integer ( kind = 4 ) dec_digit
  real ( kind = 4 ) factor
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) iabot
  integer ( kind = 4 ) iatop
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) jfac

  dec_digit = 0
  call i4_data ( 'GET', 'DEC_DIGIT', dec_digit )

  factor = 10.0E+00**dec_digit
 
  if ( 0 < dec_digit ) then
    ifac = 10**dec_digit
    jfac = 1
  else
    ifac = 1
    jfac = 10**( - dec_digit )
  end if
 
  itop = nint ( a * factor ) * jfac
  ibot = ifac
!
!  Factor out the greatest common factor.
!
  itemp = i4_gcd ( itop, ibot )
 
  iatop = itop / itemp
  iabot = ibot / itemp
 
  return
end
subroutine r4_to_s_left ( rval, s )

!*****************************************************************************80
!
!! R4_TO_S_LEFT represents a real using 14 left_justified characters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) RVAL, a real number.
!
!    Output, character ( len = * ) S, a left-justified character variable 
!    containing the representation of RVAL.
!
  implicit none

  character ( len = 14 ) chrtmp
  integer ( kind = 4 ) i
  real ( kind = 4 ) rval
  character ( len = * ) s
!
!  We can't seem to write directly into the string because of compiler
!  quibbles.
!
  if ( real ( int ( rval ) ) == rval .and. abs ( rval ) < 1.0E+13 ) then
 
    write ( chrtmp, '(i14)' ) int ( rval )
 
  else
 
    write ( chrtmp, '(g14.6)' ) rval
 
  end if
 
  do i = 1, len ( chrtmp )
    if ( chrtmp(i:i) /= ' ' ) then
      s = chrtmp(i:)
      return
    end if
  end do

  s = ' '

  return
end
subroutine r8_to_dec ( dval, itop, ibot )

!*****************************************************************************80
!
!! R8_TO_DEC converts an R8 to a decimal representation.
!
!  Discussion:
!
!    Given the R8 value DVAL, DBLDEC computes integers ITOP and 
!    IBOT so that it is approximately true that:
!
!      DVAL = ITOP * 10 ** IBOT
!
!    Only DEC_DIGIT digits of DVAL are used in constructing the 
!    representation, where DEC_DIGIT is the number of digits in the
!    decimal representation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DVAL, the value whose decimal representation 
!    is desired.
!
!    Output, integer ( kind = 4 ) ITOP, IBOT, the approximate decimal representation of DVAL.
!    ITOP is an integer, strictly between -10**DEC_DIGIT and 10**DEC_DIGIT.
!    IBOT is an integer exponent of 10.
!
  implicit none

  integer ( kind = 4 ) dec_digit
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) itop
  real    ( kind = 8 ) dtop
  real    ( kind = 8 ) dval
  real    ( kind = 4 ) ten1
  real    ( kind = 4 ) ten2
!
!  Special cases.
!
  if ( dval == 0.0D+00 ) then
    itop = 0
    ibot = 0
    return
  end if

  dec_digit = 0
  call i4_data ( 'GET', 'DEC_DIGIT', dec_digit )
!
!  Factor DVAL = DTOP * 10**IBOT
!
  dtop = dval
  ibot = 0
!
!  Now normalize so that 10**(DEC_DIGIT-1) <= ABS(DTOP) < 10**(DEC_DIGIT)
!
  ten1 = 10.0E+00**( dec_digit - 1 )
  ten2 = 10.0E+00**dec_digit
  
  do while ( abs ( dtop ) < ten1 )
    dtop = dtop * 10.0D+00
    ibot = ibot - 1
  end do

  do while ( ten2 <= abs ( dtop ) )
    dtop = dtop / 10.0D+00
    ibot = ibot + 1
  end do
!
!  ITOP is the integer part of DTOP, rounded.
!
  itop = nint ( dtop )
!
!  Now divide out any factors of ten from ITOP.
!
  if ( itop /= 0 ) then
 
    do while ( 10 * ( itop / 10 ) == itop )
      itop = itop / 10
      ibot = ibot + 1
    end do

  end if
 
  return
end
subroutine rat_add ( ibot, ibot1, ibot2, itop, itop1, itop2 )

!*****************************************************************************80
!
!! RAT_ADD adds two rational values.
!
!  Discussion:
!
!    ITOP / IBOT = ( ITOP1 / IBOT1 ) + ( ITOP2 / IBOT2 )
!
!    If numeric overflow would occur, the computation is done in
!    real ( kind = 4 ) arithmetic and then converted back to a rational value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IBOT, the denominator of the result.
!
!    Input, integer ( kind = 4 ) IBOT1, IBOT2, the denominators of the
!    two rational values to be added.
!
!    Output, integer ( kind = 4 ) ITOP, the numerator of the result.
!
!    Input, integer ( kind = 4 ) ITOP1, ITOP2, the numerators of the
!    two rational values to be added.
!
  implicit none

  integer ( kind = 4 ) i4_big
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ibot1
  integer ( kind = 4 ) ibot2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) itop1
  integer ( kind = 4 ) itop2
  integer ( kind = 4 ) jbot1
  integer ( kind = 4 ) jbot2
  integer ( kind = 4 ) jbot3
  integer ( kind = 4 ) jtop1
  integer ( kind = 4 ) jtop2
  real ( kind = 4 ) rbot
  real ( kind = 4 ) rmax
  real ( kind = 4 ) rtop1
  real ( kind = 4 ) rtop2
  real ( kind = 4 ) rtop3
  real ( kind = 4 ) rval
  real ( kind = 4 ) rval1
  real ( kind = 4 ) rval2

  ierror = 0
 
  i4_big = 0
  call i4_data ( 'GET', 'I4_BIG', i4_big )

  rmax = real ( i4_big )

  if ( itop1 == 0 ) then
    itop = itop2
    ibot = ibot2
    return
  else if ( itop2 == 0 ) then
    itop = itop1
    ibot = ibot1
    return
  end if
!
!  Make copies of the input arguments, since we will change them.
!
  jbot1 = ibot1
  jbot2 = ibot2
  jtop1 = itop1
  jtop2 = itop2
!
!  Compute the greatest common factor of the two denominators,
!  and factor it out.
!
  jbot3 = i4_gcd ( jbot1, jbot2 )
  jbot1 = jbot1 / jbot3
  jbot2 = jbot2 / jbot3
!
!  The fraction may now be formally written as:
!
!    (jtop1*jbot2 + jtop2*jbot1) / (jbot1*jbot2*jbot3)
!
!  Check the tops for overflow.
!
  rtop1 = real ( jtop1 ) * real ( jbot2 )
 
  if ( rmax < abs ( rtop1 ) ) then
    ierror = 1
    itop = 0
  else
    jtop1 = jtop1 * jbot2
  end if
 
  rtop2 = real ( jtop2 ) * real ( jbot1 )
 
  if ( rmax < abs ( rtop2 ) ) then
    ierror = 2
    itop = 0
  else
    jtop2 = jtop2 * jbot1
  end if
 
  rtop3 = real ( jtop1 ) + real ( jtop2 )
 
  if ( rmax < abs ( rtop3 ) ) then
    ierror = 3
    itop = 0
  else
    itop = jtop1 + jtop2
  end if
!
!  Check the bottom for overflow.
!
  rbot = real ( jbot1 ) * real ( jbot2 ) * real ( jbot3 )
 
  if ( rmax < abs ( rbot ) ) then
    ierror = 4
    ibot = 1
  else
    ibot = jbot1 * jbot2 * jbot3
  end if
!
!  If there was potential overflow, then do the computation in
!  real ( kind = 4 ) arithmetic and convert back.
!
  if ( ierror /= 0 ) then
    ierror = 0
    rval1 = real ( itop1 ) / real ( ibot1 )
    rval2 = real ( itop2 ) / real ( ibot2 )
    rval = rval1 + rval2
    call r4_to_rat ( rval, itop, ibot )
  end if
!
!  Put the fraction in lowest terms.
!
  itemp = i4_gcd ( itop, ibot )
  itop = itop / itemp
  ibot = ibot / itemp
!
!  Sign of bottom should be positive.
!
  if ( ibot < 0 ) then
    ibot = - ibot
    itop = - itop
  end if
 
  return
end
subroutine rat_mul ( ibot, ibot1, ibot2, itop, itop1, itop2, ierror )

!*****************************************************************************80
!
!! RAT_MUL multiplies two fractions.
!
!  Discussion:
!
!    ITOP / IBOT = ( ITOP1 / IBOT1 ) * ( ITOP2 / IBOT2 ).
!
!    If numeric overflow would occur, the computation is done in
!    real ( kind = 4 ) arithmetic and then converted back to a rational value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IBOT, the denominator of the result.
!
!    Input, integer ( kind = 4 ) IBOT1, IBOT2, the denominators of the
!    two rational values to be multiplied.
!
!    Output, integer ( kind = 4 ) ITOP, the numerator of the result.
!
!    Input, integer ( kind = 4 ) ITOP1, ITOP2, the numerators of the
!    two rational values to be multiplied.
!
!    Output, integer ( kind = 4 ) IERROR, 0 an error occurred, 1 no error occurred.
!
  implicit none

  integer ( kind = 4 ) i4_big
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ibot1
  integer ( kind = 4 ) ibot2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) itop1
  integer ( kind = 4 ) itop2
  integer ( kind = 4 ) jbot1
  integer ( kind = 4 ) jbot2
  integer ( kind = 4 ) jtop1
  integer ( kind = 4 ) jtop2
  real ( kind = 4 ) rbot
  real ( kind = 4 ) rmax
  real ( kind = 4 ) rtop
  real ( kind = 4 ) temp

  ierror = 0
 
  i4_big = 0
  call i4_data ( 'GET', 'I4_BIG', i4_big )

  rmax = real ( i4_big )

  if ( itop1 == 0 .or. itop2 == 0 ) then
    itop = 0
    ibot = 1
    return
  end if

  if ( ibot1 == 0 .or. ibot2 == 0 ) then
    ierror = 1
    itop = 0
    ibot = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_MUL - Fatal error!'
    write ( *, '(a)' ) '  A rational fraction has a zero denominator!'
    return
  end if
!
!  Make copies of the input arguments, since we will change them.
!
  jbot1 = ibot1
  jbot2 = ibot2
  jtop1 = itop1
  jtop2 = itop2
!
!  Get rid of all common factors in top and bottom.
!
  itemp = i4_gcd ( jtop1, jbot1 )
  jtop1 = jtop1 / itemp
  jbot1 = jbot1 / itemp
  itemp = i4_gcd ( jtop1, jbot2 )
  jtop1 = jtop1 / itemp
  jbot2 = jbot2 / itemp
  itemp = i4_gcd ( jtop2, jbot1 )
  jtop2 = jtop2 / itemp
  jbot1 = jbot1 / itemp
  itemp = i4_gcd ( jtop2, jbot2 )
  jtop2 = jtop2 / itemp
  jbot2 = jbot2 / itemp
!
!  The fraction (ITOP1*ITOP2)/(IBOT1*IBOT2) is in lowest terms.
!
!  Check the top ITOP1*ITOP2 for overflow.
!
  rtop = real ( jtop1 ) * real ( jtop2 )
 
  if ( rmax < abs ( rtop ) ) then
    ierror = 1
  else
    itop = jtop1 * jtop2
  end if
!
!  Check the bottom IBOT1*IBOT2 for overflow.
!
  rbot = real ( jbot1 ) * real ( jbot2 )
 
  if ( rmax < abs ( rbot ) ) then
    ierror = 2
  else
    ibot = jbot1 * jbot2
  end if
!
!  if there was an overflow, then compute RTOP / RBOT and convert
!  back to rational.
!
  if ( ierror == 1 .or. ierror == 2 ) then
    ierror = 0
    temp = rtop / rbot
    call r4_to_rat ( temp, itop, ibot )
  end if
!
!  Sign of bottom should be positive.
!
  if ( ibot < 0 ) then
    ibot = - ibot
    itop = - itop
  end if
!
!  The fraction is ITOP/IBOT with no loss of accuracy (unless there
!  was overflow).
!
  return
end
subroutine rat_print ( iatop, iabot, nrow, ncol, title )

!*****************************************************************************80
!
!! RAT_PRINT prints out rational vectors or matrices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
!    Input, character ( len = * ) TITLE, a label for the object being printed.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 10 ) chrtmp3
  character ( len = 40 ) format1
  character ( len = 40 ) format2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) ichi
  integer ( kind = 4 ) iclo
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) ione
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) izhi
  integer ( kind = 4 ) izlo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) jmin
  integer ( kind = 4 ) kmax
  integer, parameter :: ncolum = 80
  integer ( kind = 4 ) none
  integer ( kind = 4 ) npline
  character ( len = 100 ) output
  character ( len = * ) title
!
!  Figure out how many rationals we can get in NCOLUM columns.
!
  kmax = 3
 
  do i = 1, nrow
    do j = 1, ncol
 
      itemp = abs ( iatop(i,j) )
 
      do while ( 10**(kmax-2) <= itemp )
        kmax = kmax + 1
      end do
 
      itemp = abs ( iabot(i,j) )
 
      do while ( 10**(kmax-2) < itemp )
        kmax = kmax + 1
      end do
 
    end do
  end do
 
  kmax = kmax + 1
  npline = ncolum / kmax
!
!  Create the formats.
!
  call i4_to_s_left ( npline, chrtmp2 )
  call i4_to_s_left ( kmax, chrtmp3 )

  format1 = '(' // chrtmp2 // 'i' // chrtmp3 // ')'
 
  call s_blank_delete ( format1 )
 
  format2 = '(' // chrtmp2 // 'i' // chrtmp3 // ')'
 
  call s_blank_delete ( format2 )
 
  do jmin = 1, ncol, npline
 
    jmax = min ( jmin+npline-1, ncol )

    write ( *, '(a)' ) ' '
    
    if ( jmin == 1 ) then
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
    end if
 
    if ( 1 < jmin .or. jmax < ncol ) then
      write ( output, * ) 'Columns ', jmin, ' to ', jmax
      call s_blanks_delete ( output )
      write ( *, '(a)' ) trim ( output )
      write ( *, '(a)' ) ' '  
    end if
 
    do i = 1, nrow
 
      write ( *, format1 ) iatop(i,jmin:jmax)
      write ( output, format1 ) iabot(i,jmin:jmax)
!
!  Delete each denominator that is 1.  If all are 1, don't
!  even print out the line.
!
      none = 0
 
      do j = jmin, jmax
 
        if ( iabot(i,j) == 1 ) then
          ione = (j-jmin+1) * kmax
          output(ione:ione) = ' '
        else
          none = 1
        end if
 
      end do
 
      write ( *, '(a)' ) trim ( output )

      if ( jmax == ncol .and. i == nrow ) then
      else
        write ( *, '(a)' ) ' '
      end if
 
    end do
 
  end do
 
  return
end
subroutine rat_read ( itop, ibot, line, prompt, ierror )

!*****************************************************************************80
!
!! RAT_READ reads a rational value, expressed as integer, decimal or fraction.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ITOP, IBOT, the top and bottom of the
!    fraction that was read.
!
!    Workspace, character ( len = 80 ) LINE, used to hold the user's input.
!
!    Workspace, character ( len = 80 ) PROMPT.
!
!    Output, integer ( kind = 4 ) IERROR, 0 an error occurred, 1 no error occurred.
!
  implicit none

  integer   ( kind = 4 ) i4_gcd
  integer   ( kind = 4 ) ibot
  integer   ( kind = 4 ) ibot1
  integer   ( kind = 4 ) ibot2
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) itemp
  integer   ( kind = 4 ) itop
  integer   ( kind = 4 ) itop1
  integer   ( kind = 4 ) itop2
  integer   ( kind = 4 ) lchar
  character ( len = 80 ) line
  character ( len = 80 ) prompt

  ierror = 0
  itop = 0
  ibot = 1
 
  do
 
    call chrinp ( ierror, line, prompt )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RAT_READ - Fatal error!'
      write ( *, '(a)' ) '  CHRINP returned error flag!'
      return
    end if
 
    if ( line /= ' ' ) then
      exit
    end if

  end do
 
  call chrctf ( line, itop1, ibot1, ierror, lchar )
 
  if ( len ( line ) <= lchar ) then
    itop = itop1
    ibot = ibot1
  else if ( line(lchar+1:lchar+1) /= '/' ) then
    itop = itop1
    ibot = ibot1
  else
    lchar = lchar + 1
    call s_chop ( line, 1, lchar )
    call chrctf ( line, itop2, ibot2, ierror, lchar )
    itop = itop1 * ibot2
    ibot = ibot1 * itop2
  end if
 
  call s_chop ( line, 1, lchar )
!
!  Make sure fraction is in lowest terms.
!
  itemp = i4_gcd ( itop, ibot )
  itop = itop / itemp
  ibot = ibot / itemp
 
  return
end
subroutine rat_to_dec ( iatop, iabot, ierror )

!*****************************************************************************80
!
!! RAT_TO_DEC converts a rational value to a decimal value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) IATOP, IABOT.
!
!    On input, the rational value (IATOP/IABOT) to be converted.
!
!    On output, the rational decimal value IATOP * 10**IABOT.
!
!    Output, integer ( kind = 4 ) IERROR, 0 an error occurred, 1 no error occurred.
!
  implicit none

  real    ( kind = 8 ) dval
  integer ( kind = 4 ) iabot
  integer ( kind = 4 ) iatop
  integer ( kind = 4 ) ierror

  if ( iabot /= 0 ) then

    ierror = 0

    dval = real ( iatop, kind = 8 ) / real ( iabot, kind = 8 )
 
    call r8_to_dec ( dval, iatop, iabot )

  else

    ierror = 1
    iatop = 0
    iabot = 0

  end if
 
  return
end
subroutine rat_to_r4 ( iatop, iabot, a )

!*****************************************************************************80
!
!! RAT_TO_R4 converts rational values to R4s.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IATOP, IABOT, the rational quantity
!    (IATOP/IABOT) that is to be converted.
!
!    Output, real ( kind = 4 ) A, the value of the rational quantity.
!
  implicit none

  real    ( kind = 4 ) a
  integer ( kind = 4 ) iabot
  integer ( kind = 4 ) iatop

  if ( iabot == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RATREL - Warning!'
    write ( *, '(a)' ) '  The input fraction had a zero denominator.'
    a = 0.0E+00
  else
    a = real ( iatop ) / real ( iabot )
  end if
 
  return
end
subroutine rat_to_s_left ( ival, jval, string )

!*****************************************************************************80
!
!! RAT_TO_S_LEFT returns a left-justified representation of IVAL/JVAL.
!
!  Discussion:
!
!    If the ratio is negative, a minus sign precedes IVAL.
!    A slash separates IVAL and JVAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, JVAL, the two integers whose
!    ratio IVAL/JVAL is to be represented.
!
!    Note that if IVAL is nonzero and JVAL is 0, STRING will
!    be returned as "Inf" or "-Inf" (Infinity), and if both
!    IVAL and JVAL are zero, STRING will be returned as "NaN"
!    (Not-a-Number).
!
!    Output, character ( len = 22 ) STRING, a left-justified string
!    containing the representation of IVAL/JVAL.
!
  implicit none

  integer   ( kind = 4 ) ival
  integer   ( kind = 4 ) ival2
  integer   ( kind = 4 ) jval
  integer   ( kind = 4 ) jval2
  character ( len = * ) string
!
!  Take care of simple cases right away.
!
  if ( ival == 0 ) then
 
    if ( jval /= 0 ) then
      string = '0'
    else
      string = 'NaN'
    end if
 
  else if ( jval == 0 ) then
 
    if ( 0 < ival ) then
      string = 'Inf'
    else
      string = '-Inf'
    end if
!
!  Make copies of IVAL and JVAL.
!
  else
 
    ival2 = ival
    jval2 = jval
 
    if ( jval2 == 1 ) then
      write ( string, '(i11)' ) ival2
    else
      write ( string, '(i11, ''/'', i10)' ) ival2, jval2
    end if
 
    call s_blank_delete ( string )
 
  end if
 
  return
end
subroutine row_add ( a, a_int, iatop, iabot, ierror, base, iform, irow1, &
  irow2, nrow, ncol, sval, istop, isbot )

!*****************************************************************************80
!
!! ROW_ADD adds a multiple of one row to another.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(NROW,NCOL), the current matrix.
!
!    Input/output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input/output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 for an error.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) IROW1, the row which is to be modified.
!
!    Input, integer ( kind = 4 ) IROW2, the row which is to be multiplied by
!    a given value and added to row IROW1.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
!    Input, real ( kind = 4 ) SVAL, the multiplier to use if real arithmetic is employed.
!
!    Input, integer ( kind = 4 ) ISTOP, ISBOT, the fractional or decimal
!    multiplier to use.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  integer ( kind = 4 ) base
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 22 ) chrtmp3
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) irow1
  integer ( kind = 4 ) irow2
  integer ( kind = 4 ) isbot
  integer ( kind = 4 ) isbot2
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) istop2
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) j
  real ( kind = 4 ) sval
!
!  Return immediately if the multiplier is zero.
!
  if ( ( iform == 0 .and. istop == 0 ) .or. &
       ( iform == 1 .and. sval == 0.0E+00 ) .or. &
       ( iform == 2 .and. istop == 0 ) .or. &
       ( iform == 3 .and. istop == 0 ) .or. &
       ( iform == 4 .and. istop == 0 ) ) then
      return
    end if
!
!  Carry out the operation.
!
  if ( iform == 0 ) then
 
    do j = 1, ncol
 
      call rat_mul ( isbot2, isbot, iabot(irow2,j), istop2, istop, &
        iatop(irow2,j), ierror )
 
      if ( ierror /= 0 ) then
        return
      end if

      call rat_add ( ibot, iabot(irow1,j), isbot2, itop, iatop(irow1,j), istop2 )
 
      iatop(irow1,j) = itop
      iabot(irow1,j) = ibot
 
    end do
 
  else if ( iform == 1 ) then
 
    a(irow1,1:ncol) = a(irow1,1:ncol) + sval * a(irow2,1:ncol)
 
  else if ( iform == 2 ) then
 
    do j = 1, ncol
 
      call dec_mul ( isbot2, isbot, iabot(irow2,j), istop2, istop, &
        iatop(irow2,j) )
 
      call dec_add ( ibot, iabot(irow1,j), isbot2, itop, &
        iatop(irow1,j), istop2 )
 
      iatop(irow1,j) = itop
      iabot(irow1,j) = ibot
 
    end do
 
  else if ( iform == 3 ) then
 
    a_int(irow1,1:ncol) = a_int(irow1,1:ncol) + istop * a_int(irow2,1:ncol)

  else if ( iform == 4 ) then
 
    a_int(irow1,1:ncol) = &
      mod ( a_int(irow1,1:ncol) + istop * a_int(irow2,1:ncol), base )
 
  end if
!
!  Print out a message.
!
  if ( iform == 0 ) then
 
    if ( istop == isbot ) then
      chrtmp3 = '+'
    else if ( istop == - isbot ) then
      chrtmp3 = '-'
    else
      call rat_to_s_left ( istop, isbot, chrtmp3 )
    end if
 
  else if ( iform == 1 ) then
 
    if ( sval == 1.0E+00 ) then
      chrtmp3 = '+'
    else if ( sval == - 1.0E+00 ) then
      chrtmp3 = '-'
    else
      call r4_to_s_left ( sval, chrtmp3 )
    end if
 
  else if ( iform == 2 ) then

    if ( istop == 1 .and. isbot == 0 ) then 
      chrtmp3 = '+'
    else if ( istop == - 1 .and. isbot == 0 ) then
      chrtmp3 = '-'
    else
      call dec_to_s_left ( istop, isbot, chrtmp3 )
    end if

  else if ( iform == 3 ) then
 
    if ( istop == 1 ) then
      chrtmp3 = '+'
    else if ( istop == - 1 ) then
      chrtmp3 = '-'
    else
      call i4_to_s_left ( istop, chrtmp3 )
    end if

  else if ( iform == 4 ) then
 
    if ( istop == 1 ) then
      chrtmp3 = '+'
    else if ( istop == - 1 ) then
      chrtmp3 = '-'
    else
      call i4_to_s_left ( istop, chrtmp3 )
    end if

  end if
 
  call i4_to_s_left ( irow1, chrtmp1 )
  call i4_to_s_left ( irow2, chrtmp2 )

  if ( chrtmp3 == '-' .or. chrtmp3 == '+' ) then

    write ( *, '(a)' ) '  ERO: Row ' // trim ( chrtmp1 ) // ' <=  ' // 'Row ' // &
      trim ( chrtmp1 ) // ' ' // chrtmp3(1:1) // ' Row ' // trim ( chrtmp2 )

  else if ( chrtmp3(1:1) == '-' .or. chrtmp3(1:1) == '+' ) then

    write ( *, '(a)' ) '  ERO: Row ' // trim ( chrtmp1 ) // ' <=  ' // 'Row ' // &
      trim ( chrtmp1 ) // ' ' // chrtmp3(1:1) // ' ' // trim ( chrtmp3 ) // &
      ' Row ' // trim ( chrtmp2 )

  else

    write ( *, '(a)' ) '  ERO: Row ' // trim ( chrtmp1 ) // ' <=  ' // 'Row ' // &
      trim ( chrtmp1 ) // ' + ' // trim ( chrtmp3 ) // ' Row ' // &
      trim ( chrtmp2 )

  end if

  return
end
subroutine row_add_param ( ierror, base, iform, irow1, irow2, istop, isbot, &
  line, nrow, sval )

!*****************************************************************************80
!
!! ROW_ADD_PARAM gets and checks the row add parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 for an error.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) IROW1, the row to which the multiple is to be added.
!
!    Input, integer ( kind = 4 ) IROW2, the row which is to be multiplied and
!    added to another row.
!
!    Output, integer ( kind = 4 ) ISTOP, ISBOT, the parts of the rational
!    or decimal fraction of the multiplier, if that is the
!    arithmetic being used.
!
!    Workspace, character ( len = 80 ) LINE, used to hold the user's input.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows in the matrix.
!
!    Output, real ( kind = 4 ) SVAL, the multiplier, if real arithmetic is used.
!
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) irow1
  integer ( kind = 4 ) irow2
  integer ( kind = 4 ) isbot
  integer ( kind = 4 ) istop
  character ( len = 80 ) line
  integer ( kind = 4 ) nrow
  character ( len = 80 ) prompt
  real ( kind = 4 ) sval

  prompt = 'multiplier S, row I to add, target row J.'
!
!  Get the multiplier, SVAL or ISTOP/ISBOT.
!
  if ( iform == 0 ) then
 
    call rat_read ( istop, isbot, line, prompt, ierror )
 
  else if ( iform == 1 ) then
 
    call r4_read ( sval, line, prompt, ierror )
 
  else if ( iform == 2 ) then
 
    call dec_read ( istop, isbot, line, prompt, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ROW_ADD_PARAM - Error!'
      write ( *, '(a)' ) '  DEC_READ returns error code.'
      return
    end if
 
    call dec_round ( istop, isbot )
 
  else if ( iform == 3 ) then
    call i4_read ( istop, line, prompt, ierror )
  else if ( iform == 4 ) then
    call i4_read ( istop, line, prompt, ierror )
    istop = mod ( istop, base )
  end if

  if ( ierror /= 0 ) then
    return
  end if
!
!  Get the row to add, IROW2.
!
  call i4_read ( irow2, line, prompt, ierror )

  if ( ierror /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Error reading row index from line:'
    write ( *, '(a)' ) '"' // trim ( line ) // '"'
    return
  end if

  if ( irow2 < 1 .or. nrow < irow2 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Error!  Row index was not acceptable!'
    return
  end if
!
!  Get the row to which we are adding, IROW1.
!
  call i4_read ( irow1, line, prompt, ierror )
  if ( ierror /= 0 ) return
 
  if ( irow1 < 1 .or. nrow < irow1 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Error!  The row index was not acceptable!'
    return
  end if
!
!  Make sure the rows are different.
!
  if ( irow1 == irow2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Error!  The rows should not be the same!'
    ierror = 1
    return
  end if
 
  ierror = 0
 
  return
end
subroutine row_auto ( a, a_int, iatop, iabot, ierror, base, iform, nrow, ncol )

!*****************************************************************************80
!
!! ROW_AUTO automatically row reduces the current matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(NROW,NCOL), the current matrix.
!
!    Input/output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input/output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL).
!    IATOP and IABOT represent the rational or decimal matrix
!    to which elementary row operations will be applied.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 for an error.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  real ( kind = 4 ) amax
  real ( kind = 4 ) atemp
  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) isbot
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) krow
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lrow
  real ( kind = 4 ) sval

  ierror = 0
  krow = 0
 
  do i = 1, nrow
 
    irow = i
 
    do j = 1, ncol
 
      jcol = j
!
!  In column JCOL, seek the row between IROW and NROW with
!  maximum nonzero entry AMAX.
!
      imax = 0
      amax = 0.0E+00
 
      do krow = irow, nrow
 
        if ( iform == 0 ) then
          call rat_to_r4 ( iatop(krow,jcol), iabot(krow,jcol), atemp )
        else if ( iform == 1 ) then
          atemp = a(krow,jcol)
        else if ( iform == 2 ) then
          call dec_to_r4 ( iatop(krow,jcol), iabot(krow,jcol), atemp )
        else if ( iform == 3 ) then
          atemp = real ( a_int(krow,jcol) )
        else if ( iform == 4 ) then
          atemp = real ( a_int(krow,jcol) )
        end if
 
        atemp = abs ( atemp )
 
        if ( amax < atemp ) then
          amax = atemp
          imax = krow
        end if
 
      end do
 
      if ( imax /= 0 ) then
        krow = imax
        exit
      end if
 
    end do
 
    if ( krow == 0 ) then
      return
    end if
!
!  Interchange the IROW-th and the pivot rows.
!
    if ( krow /= irow ) then
      call row_swap ( a, a_int, iatop, iabot, ierror, iform, &
        krow, irow, nrow, ncol )
    end if
!
!  Divide the pivot row by A(IROW,JCOL) so that A(IROW,JCOL) = 1.
!
    if ( iform == 0 ) then
      istop = iatop(irow,jcol)
      isbot = iabot(irow,jcol)
    else if ( iform == 1 ) then
      sval = a(irow,jcol)
    else if ( iform == 2 ) then
      istop = iatop(irow,jcol)
      isbot = iabot(irow,jcol)
    else if ( iform == 3 ) then
      istop = a_int(irow,jcol)
    else if ( iform == 4 ) then
      istop = a_int(irow,jcol)
    end if
 
    call row_div ( a, a_int, iatop, iabot, ierror, base, iform, irow, &
      nrow, ncol, sval, istop, isbot )
!
!  Annihilate A(L,JCOL) for L not equal to IROW.
!
    do l = 1, nrow
 
      lrow = l
 
      if ( lrow /= irow ) then
 
        if ( iform == 0 ) then

          if ( iatop(lrow,jcol) /= 0 ) then

            istop = - iatop(lrow,jcol)
            isbot = iabot(lrow,jcol)

            call row_add ( a, a_int, iatop, iabot, ierror, base, iform, &
              lrow, irow, nrow, ncol, sval, istop, isbot )

            iatop(lrow,jcol) = 0
            iabot(lrow,jcol) = 1

          end if

        else if ( iform == 1 ) then

          if ( a(lrow,jcol) /= 0.0E+00 ) then

            sval = - a(lrow,jcol)

            call row_add ( a, a_int, iatop, iabot, ierror, base, iform, &
              lrow, irow, nrow, ncol, sval, istop, isbot )

            a(lrow,jcol) = 0.0E+00

          end if

        else if ( iform == 2 ) then

          if ( iatop(lrow,jcol) /= 0 ) then

            istop = - iatop(lrow,jcol)
            isbot = iabot(lrow,jcol)

            call row_add ( a, a_int, iatop, iabot, ierror, base, iform, &
              lrow, irow, nrow, ncol, sval, istop, isbot )

            iatop(lrow,jcol) = 0
            iabot(lrow,jcol) = 0

          end if

        else if ( iform == 3 ) then

          if ( a_int(lrow,jcol) /= 0 ) then

            istop = - a_int(lrow,jcol)

            call row_add ( a, a_int, iatop, iabot, ierror, base, iform, &
              lrow, irow, nrow, ncol, sval, istop, isbot )

            a_int(lrow,jcol) = 0

          end if

        else if ( iform == 4 ) then

          if ( a_int(lrow,jcol) /= 0 ) then

            istop = - a_int(lrow,jcol)

            call row_add ( a, a_int, iatop, iabot, ierror, base, iform, &
              lrow, irow, nrow, ncol, sval, istop, isbot )

            a_int(lrow,jcol) = 0

          end if

        end if
 
      end if
 
    end do
 
  end do
 
  return
end
subroutine row_div ( a, a_int, iatop, iabot, ierror, base, iform, irow, &
  nrow, ncol, sval, istop, isbot )

!*****************************************************************************80
!
!! ROW_DIV divides row IROW of the A matrix by a scale factor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(NROW,NCOL), the current matrix.
!
!    Input/output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input/output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL).
!    IATOP and IABOT represent the current rational or decimal
!    matrix.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 for an error.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) IROW, the row to be divided.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
!    Input, real ( kind = 4 ) SVAL, the real divisor.
!
!    Input, integer ( kind = 4 ) ISTOP, ISBOT, the fractional or decimal divisor.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  integer ( kind = 4 ) base
  character ( len = 10 ) chrtmp1
  character ( len = 24 )  chrtmp2
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) isbot
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) j
  character ( len = 3 ) op
  real ( kind = 4 ) sval
!
!  Make sure that the row number is legal.
!
  if ( irow < 1 .or. nrow < irow ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Error!  The row number is out of range!'
    ierror = 1
    return
  end if
!
!  Check for an illegal divisor of 0.
!
  if ( ( iform == 0 .and. istop == 0 ) .or. &
       ( iform == 1 .and. sval == 0.0E+00 ) .or. &
       ( iform == 2 .and. istop == 0 ) .or. &
       ( iform == 3 .and. istop == 0 ) .or. &
       ( iform == 4 .and. istop == 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROW_DIV - Error!'
    write ( *, '(a)' ) '  Requested division by zero.'
    ierror = 1
    return
  end if
!
!  Check for a pointless divisor of 1.
!
  if ( ( iform == 0 .and. istop == isbot ) .or. &
       ( iform == 1 .and. sval == 1.0E+00 ) .or.&
       ( iform == 2 .and. istop == 1 .and. isbot == 0 ) .or. &
       ( iform == 3 .and. istop == 1 ) .or. &
       ( iform == 4 .and. istop == 1 ) ) then
    return
  end if
!
!  Carry out the division.
!
  if ( iform == 0 ) then
 
    do j = 1, ncol
 
      call rat_mul ( ibot, iabot(irow,j), istop, itop, iatop(irow,j), isbot, &
        ierror )
 
      if ( ierror /= 0 ) then
        return
      end if
 
      iatop(irow,j) = itop
      iabot(irow,j) = ibot
 
    end do
 
  else if ( iform == 1 ) then
 
    a(irow,1:ncol) = a(irow,1:ncol) / sval
 
  else if ( iform == 2 ) then
 
    do j = 1, ncol
 
      call dec_div ( iatop(irow,j), iabot(irow,j), istop, isbot, &
        iatop(irow,j), iabot(irow,j), ierror )
 
    end do

  else if ( iform == 3 ) then
 
    a_int(irow,1:ncol) = a_int(irow,1:ncol) / istop

   else if ( iform == 4 ) then
 
    a_int(irow,1:ncol) = a_int(irow,1:ncol) / istop
 
  end if
!
!  Print out a statement about what has been done.
!
  if ( iform == 0 ) then
 
    if ( isbot == 1 ) then

      call i4_to_s_left ( istop, chrtmp2 )
      op = ' / '

    else

      call rat_to_s_left ( isbot, istop, chrtmp2 )
      op = ' * '

    end if
 
  else if ( iform == 1 ) then
 
    call r4_to_s_left ( sval, chrtmp2 )
    op = ' / '

  else if ( iform == 2 ) then
 
    call dec_to_s_left ( istop, isbot, chrtmp2 )
    op = ' / '

  else if ( iform == 3 ) then
 
    call i4_to_s_left ( istop, chrtmp2 )
    op = ' / '

  else if ( iform == 4 ) then
 
    call i4_to_s_left ( istop, chrtmp2 )
    op = ' / '
  end if

  call i4_to_s_left ( irow, chrtmp1 )

  write ( *, '(a)' ) '  ERO: Row ' // trim ( chrtmp1 ) // ' <=  Row ' // &
    trim ( chrtmp1 ) // op // trim ( chrtmp2 )
 
  return
end
subroutine row_div_param ( ierror, base, iform, irow, isbot, istop, line, sval )

!*****************************************************************************80
!
!! ROW_DIV_PARAM gets and checks the row divide parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 for an error.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Output, integer ( kind = 4 ) IROW, the row to be divided.
!
!    Workspace, character ( len = 80 ) LINE, used to hold the user's input.
!
!    Output, real ( kind = 4 ) SVAL, the real divisor.
!
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) isbot
  integer ( kind = 4 ) istop
  character ( len = 80 ) line
  character ( len = 80 ) prompt
  real ( kind = 4 ) sval

  prompt = 'row I, divisor S.'
!
!  Read the row number to be divided.
!
  call i4_read ( irow, line, prompt, ierror )
  if ( ierror /= 0 ) then
    return
  end if
!
!  Read the divisor, either RVAL or ISTOP/ISBOT.
!
  if ( iform == 0 ) then
 
    call rat_read ( istop, isbot, line, prompt, ierror )
 
  else if ( iform == 1 ) then
 
    call r4_read ( sval, line, prompt, ierror )
 
  else if ( iform == 2 ) then
 
    call dec_read ( istop, isbot, line, prompt, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ROW_DIV_PARAM - Fatal error!'
      write ( *, '(a)' ) '  DEC_READ returned error flag.'
      return
    end if

    call dec_round ( istop, isbot )
 
  else if ( iform == 3 ) then
    call i4_read ( istop, line, prompt, ierror )
  else if ( iform == 4 ) then
    call i4_read ( istop, line, prompt, ierror ) 
  end if
 
  if ( ierror /= 0 ) then
    return
  end if

  return
end
subroutine row_mul ( a, a_int, iatop, iabot, ierror, base, iform, irow, &
  nrow, ncol, sval, istop, isbot )

!*****************************************************************************80
!
!! ROW_MUL multiplies a row of the matrix by a scale factor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(NROW,NCOL), the current matrix.
!
!    Input/output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input/output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 for an error.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) IROW, the row that is to be multiplied.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
!    Input, real ( kind = 4 ) SVAL, the real row multiplier.
!
!    Input, integer ( kind = 4 ) ISTOP, ISBOT, the decimal or fractional row multiplier.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  integer ( kind = 4 ) base
  character ( len = 10 ) chrtmp1
  character ( len = 22 ) chrtmp
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) isbot
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) j
  real ( kind = 4 ) sval
!
!  Make sure row number is OK.
!
  if ( irow < 1 .or. nrow < irow ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Error!  The row number is out of range!'
    ierror = 1
    return
  end if
!
!  For rational arithmetic, make sure bottom of scale factor
!  is not 0.
!
  if ( iform == 0 .and. isbot == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Error!  Illegal 0 divisor in multiplier!'
    ierror = 1
    return
  end if
!
!  Check for multiplication by 0.
!
  if ( ( iform == 0 .and. istop == 0 ) .or. &
       ( iform == 1 .and. sval == 0.0E+00 ) .or. &
       ( iform == 2 .and. istop == 0 ) .or. &
       ( iform == 3 .and. istop == 0 ) .or. &
       ( iform == 4 .and. istop == 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Warning - Multiplication by zero is not an ERO.'
    ierror = 1
    return
  end if
!
!  Check for multiplication by 1.
!
  if ( ( iform == 0 .and. istop == isbot ) .or. &
       ( iform == 1 .and. sval == 1.0E+00 ) .or. &
       ( iform == 2 .and. istop == 1 .and. isbot == 0 ) .or. &
       ( iform == 3 .and. istop == 1 ) .or. &
       ( iform == 4 .and. istop == 1 ) ) then
    return
  end if
!
!  Carry out the multiplication.
!
  if ( iform == 0 ) then
 
    do j = 1, ncol
 
      call rat_mul ( ibot, iabot(irow,j), isbot, itop, iatop(irow,j), istop, &
        ierror )

      if ( ierror /= 0 ) then
        return
      end if

      iatop(irow,j) = itop
      iabot(irow,j) = ibot
 
    end do
 
  else if ( iform == 1 ) then
 
    a(irow,1:ncol) = sval * a(irow,1:ncol)
 
  else if ( iform == 2 ) then
 
    do j = 1, ncol
 
      call dec_mul ( ibot, iabot(irow,j), isbot, itop, iatop(irow,j), istop )

      if ( ierror /= 0 ) return
 
      iatop(irow,j) = itop
      iabot(irow,j) = ibot
 
    end do
 
  else if ( iform == 3 ) then

    a_int(irow,1:ncol) = istop * a_int(irow,1:ncol)

  else if ( iform == 4 ) then

    a_int(irow,1:ncol) = mod ( istop * a_int(irow,1:ncol), base )

  end if
!
!  Confirm the operation.
!
  if ( iform == 0 ) then
 
    if ( istop == - isbot ) then
      chrtmp = '-'
    else
      call rat_to_s_left ( istop, isbot, chrtmp )
    end if
 
  else if ( iform == 1 ) then
 
    if ( sval == - 1.0E+00 ) then
      chrtmp = '-'
    else
      call r4_to_s_left ( sval, chrtmp )
    end if
 
  else if ( iform == 2 ) then

    if ( istop == - 1 .and. isbot == 0 ) then
      chrtmp = '-'
    else
      call dec_to_s_left ( istop, isbot, chrtmp )
    end if

  else if ( iform == 3 ) then
 
    if ( istop == - 1 ) then
      chrtmp = '-'
    else
      call i4_to_s_left ( istop, chrtmp )
    end if

  else if ( iform == 4 ) then
 
    if ( istop == - 1 ) then
      chrtmp = '-'
    else
      call i4_to_s_left ( istop, chrtmp )
    end if

  end if
 
  call i4_to_s_left ( irow, chrtmp1 )

  write ( *, '(a)' ) '  ERO: Row ' // trim ( chrtmp1 ) // ' <= ' // &
    trim ( chrtmp ) // ' Row ' // trim ( chrtmp1 )

  return
end
subroutine row_mul_param ( ierror, base, iform, irow, istop, isbot, line, rval )

!*****************************************************************************80
!
!! ROW_MUL_PARAM gets and checks the row multiply parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 for an error.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Output, integer ( kind = 4 ) IROW, the row to be multiplied.
!
!    Output, integer ( kind = 4 ) ISTOP, ISBOT, the multiplier to use for
!    fractional or decimal arithmetic.
!
!    Workspace, character ( len = 80 ) LINE, used to hold the user's input.
!
!    Output, real ( kind = 4 ) RVAL, the multiplier to use for real arithmetic.
!
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) isbot
  integer ( kind = 4 ) istop
  character ( len = 80 ) line
  character ( len = 80 ) prompt
  real ( kind = 4 ) rval

  prompt = 'row I, multiplier S.'
!
!  Read the row number to be multiplied.
!
  call i4_read ( irow, line, prompt, ierror )
  if ( ierror /= 0 ) return
!
!  Read the multiplier, either RVAL or ISTOP/ISBOT.
!
  if ( iform == 0 ) then
 
    call rat_read ( istop, isbot, line, prompt, ierror )
 
  else if ( iform == 1 ) then
 
    call r4_read ( rval, line, prompt, ierror )
 
  else if ( iform == 2 ) then
 
    call dec_read ( istop, isbot, line, prompt, ierror )
 
    call dec_round ( istop, isbot )

  else if ( iform == 3 ) then

    call i4_read ( istop, line, prompt, ierror )

  else if ( iform == 4 ) then

    call i4_read ( istop, line, prompt, ierror )
    istop = mod ( istop, base ) 
  end if
 
  return
end
subroutine row_op_check ( command, ierror, line2 )

!*****************************************************************************80
!
!! ROW_OP_CHECK checks for commands given in the form of ERO's.
!
!  Discussion:
!
!    The form of the elementary row operation commands includes:
!
!    The row interchange command:
!      RI1 <=> RI2
!    Note that this will fail if user types "R I1 <=> R I2"
!
!    The scalar multiply command:
!      RI1 <= S * RI1
!    with or without the "*".
!
!    The scalar divide command:
!      RI1 <= RI1 / S
!
!    The add row command:
!      RI1 <= RI1 + S * RI2
!    or
!      RI1 <= S * RI2 + RI1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = 4 ) COMMAND.
!    If the routine decides that the user has input an ERO in the
!    natural format, then COMMAND contains the necessary
!    one letter MATMAN command to carry out the ERO.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 for an error.
!
!    Workspace, character ( len = 80 ) LINE2, a copy of the user input in LINE.
!
  implicit none

  character ( len = 22 ) chrtmp
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 20 ) command
  integer ( kind = 4 ) idbot2
  integer ( kind = 4 ) idbot3
  integer ( kind = 4 ) idtop2
  integer ( kind = 4 ) idtop3
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) irow1
  integer ( kind = 4 ) irow2
  integer ( kind = 4 ) irow3
  integer ( kind = 4 ) isbot2
  integer ( kind = 4 ) isbot3
  integer ( kind = 4 ) istop2
  integer ( kind = 4 ) istop3
  integer ( kind = 4 ) lchar
  logical ldiv
  character ( len = 80 ) line2
  character ( len = 80 ) string

  command = ' '
!
!  1. Remove all blanks from the line, and capitalize it.
!
  call s_blank_delete ( line2 )
  call s_cap ( line2 )
!
!  2. Is the first character an "R" or "ROW"?
!
  if ( line2(1:1) /= 'R' ) then
    return
  end if
 
  if ( line2(1:3) == 'ROW' ) then
    call s_chop ( line2, 1, 3 )
  else
    call s_chop ( line2, 1, 1 )
  end if
!
!  3. The next item should be a row number, IROW1.
!
  call s_to_i4 ( line2, irow1, ierror, lchar )
 
  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Your ERO command could not be understood.'
    write ( *, '(a)' ) '  The first row number "R1" did not make sense.'
    return
  end if
 
  call s_chop ( line2, 1, lchar )
!
!  4. Check for the row interchange string "=", "<>", "<=>" or "<->".
!
  if ( line2(1:2) == '<>' ) then
    string = '<>'
  else if ( line2(1:3) == '<=>' ) then
    string = '<=>'
  else if ( line2(1:3) == '<->' ) then
    string = '<->'
  else if ( line2(1:2) == '<=' ) then
    string = '<='
  else if ( line2(1:2) == '<-' ) then
    string = '<-'
  else if ( line2(1:2) == '=>' ) then
    string = '=>'
  else if ( line2(1:2) == '->' ) then
    string = '->'
  else if ( line2(1:1) == '=' ) then
    string = '='
  else if ( line2(1:2) == ':=' ) then
    string = ':='
  else
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Your ERO command could not be understood.'
    write ( *, '(a)' ) '  The assignment symbol <=> was missing.'
    return
  end if
 
  lchar = len_trim ( string )
 
  call s_chop ( line2, 1, lchar )
!
!  5. The next quantity could be an explicit signed scalar, S2,
!     or an implicit +-1.
!
  if ( line2(1:1) == 'R' ) then

    istop2 = 1.0E+00
    isbot2 = 1.0E+00

  else

    if ( line2(1:2) == '+R' ) then
      istop2 = 1.0E+00
      isbot2 = 1.0E+00
      call s_chop ( line2, 1, 1 )
    else if ( line2(1:2) == '-R' ) then
      istop2 = - 1.0E+00
      isbot2 = 1.0E+00
      call s_chop ( line2, 1, 1 )
    else
      call chrctg ( line2, istop2, isbot2, ierror, lchar )
      call s_chop ( line2, 1, lchar )
 
      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Your ERO command could not be understood.'
        write ( *, '(a)' ) '  The multiplier S2 did not make sense.'
        ierror = 1
        return
      end if
 
    end if
  end if
!
!  6. Is the next character an optional "*"?
!
  if ( line2(1:1) == '*' ) then
    call s_chop ( line2, 1, 1 )
  end if
!
!  7. Is the next character an "R"?
!
  if ( line2(1:3) == 'ROW' ) then
    call s_chop ( line2, 1, 3 )
  else if ( line2(1:1) == 'R' ) then
    call s_chop ( line2, 1, 1 )
  else
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Your ERO command could not be understood.'
    write ( *, '(a)' ) '  Could not find the second row index.'
    return
  end if
!
!  8. The next item should be a row number, IROW2.
!
  call s_to_i4 ( line2, irow2, ierror, lchar )
 
  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Your ERO command could not be understood.'
    write ( *, '(a)' ) '  The second row number "R2" did not make sense.'
    ierror = 1
    return
  else
    call s_chop ( line2, 1, lchar )
  end if
!
!  9. If there's nothing more, this must be an interchange
!     or a scaling.  Form the equivalent MATMAN command.
!
  if ( line2 == ' ' ) then
 
    if ( irow1 == irow2 ) then

      command = 'ROW_MUL'

      if ( isbot2 < 0 ) then
        isbot2 = - isbot2
        istop2 = - istop2
      end if

      call rat_to_s_left ( istop2, isbot2, chrtmp )
      call i4_to_s_left ( irow1, chrtmp1 )
      line2 = chrtmp1 // ' ' // chrtmp
      call s_blanks_delete ( line2 )
      return
    end if
 
    if ( istop2 == 1 .and. isbot2 == 1 ) then
      command = 'ROW_SWAP'
      call i4_to_s_left ( irow1, chrtmp1 )
      call i4_to_s_left ( irow2, chrtmp2 )
      line2 = chrtmp1 // ' ' // chrtmp2
      call s_blanks_delete ( line2 )
      return
    end if
 
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Your ERO command could not be understood.'
    write ( *, '(a)' ) '  A MULTIPLY command must have R1 and R2 the same.'
    write ( *, '(a)' ) '  An INTERCHANGE command cannot have a multiplier.'
    return
  end if
!
!  10. Is the next quantity a '/', or perhaps a '*'?
!
  ldiv = .false.
 
  if ( line2(1:1) == '/' ) then
 
    ldiv = .true.
    call s_chop ( line2, 1, 1 )
    call chrctg ( line2, idtop2, idbot2, ierror, lchar )
 
    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Your ERO command could not be understood.'
      write ( *, '(a)' ) '  The divisor of row 2 did not make sense.'
      return
    end if
 
    istop2 = istop2 * idbot2
    isbot2 = isbot2 * idtop2
 
    if ( irow1 == irow2 ) then

      if ( ldiv ) then
        command = 'ROW_DIV'
        call i4_swap ( istop2, isbot2 )
      else
        command = 'ROW_MUL'
      end if

      if ( isbot2 < 0 ) then
        isbot2 = - isbot2
        istop2 = - istop2
      end if

      call i4_to_s_left ( irow1, chrtmp1 )
      call rat_to_s_left ( istop2, isbot2, chrtmp )
      line2 = chrtmp1 // ' ' // chrtmp
      call s_blanks_delete ( line2 )
      return

    end if
 
  else if ( line2(1:1) == '*' ) then
 
    call s_chop ( line2, 1, 1 )
    call chrctg ( line2, idtop2, idbot2, ierror, lchar )
 
    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Your ERO command could not be understood.'
      write ( *, '(a)' ) '  The multiplier of row 2 did not make sense.'
      return
    end if
 
    istop2 = istop2 * idtop2
    isbot2 = isbot2 * idbot2
 
    if ( irow1 == irow2 ) then

      command = 'ROW_MUL'

      if ( isbot2 < 0 ) then
        isbot2 = - isbot2
        istop2 = - istop2
      end if

      call i4_to_s_left ( irow1, chrtmp1 )
      call rat_to_s_left ( istop2, isbot2, chrtmp )
      line2 = chrtmp1 // ' ' // chrtmp
      call s_blanks_delete ( line2 )
      return

    end if
 
  end if
!
!  11. Is the next quantity a scalar, S3?
!
  if ( line2(1:2) == '+R' ) then
 
    istop3 = 1.0E+00
    isbot3 = 1.0E+00
    call s_chop ( line2, 1, 1) 
 
  else if ( line2(1:2) == '-R' ) then
 
    istop3 = - 1.0E+00
    isbot3 = 1.0E+00
    call s_chop ( line2, 1, 1 )
 
  else
 
    call chrctg ( line2, istop3, isbot3, ierror, lchar )
 
    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Your ERO command could not be understood.'
      write ( *, '(a)' ) '  The multiplier S2 did not make sense.'
      ierror = 1
      return
    end if
 
    call s_chop ( line2, 1, lchar )
 
  end if
!
!  12. Is the next quantity an optional "*"?
!
  if ( line2(1:1) == '*' ) then
    call s_chop ( line2, 1, 1 )
  end if
!
!  13. Is the next quantity an "R" or ROW?
!
  if ( line2(1:3) == 'ROW' ) then
    call s_chop ( line2, 1, 3 )
  else if ( line2(1:1) == 'R' ) then
    call s_chop ( line2, 1, 1) 
  else
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Your ERO command could not be understood.'
    write ( *, '(a)' ) '  The "R" marking the third row was misplaced.'
    return
  end if
!
!  14. The next item should be a row number, IROW3.
!
  call s_to_i4 ( line2, irow3, ierror, lchar )
 
  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Your ERO command could not be understood.'
    write ( *, '(a)' ) '  The third row number "R3" did not make sense.'
    ierror = 1
    return
  end if
 
  call s_chop ( line2, 1, lchar )
!
!  15. Is the next quantity a '/', or perhaps a '*'?
!
  if ( line2(1:1) == '/' ) then
 
    call s_chop ( line2, 1, 1 )
    call chrctg ( line2, idtop3, idbot3, ierror, lchar )
 
    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Your ERO command could not be understood.'
      write ( *, '(a)' ) '  The divisor of R3 did not make sense.'
      return
    end if
 
    istop3 = istop3 * idbot3
    isbot3 = isbot3 * idtop3
 
  else if ( line2(1:1) == '*' ) then
 
    call s_chop ( line2, 1, 1 )
    call chrctg ( line2, idtop3, idbot3, ierror, lchar )
 
    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Your ERO command could not be understood.'
      write ( *, '(a)' ) '  The multiplier of R3 did not make sense.'
      return
    end if
 
    istop3 = istop3 * idtop3
    isbot3 = isbot3 * idbot3
 
  end if
!
!  16. Form the equivalent MATMAN ADD command.
!
  if ( irow1 == irow2 ) then

    command = 'ROW_ADD'

    if ( isbot3 < 0 ) then
      isbot3 = - isbot3
      istop3 = - istop3
    end if

    call i4_to_s_left ( irow3, chrtmp1 )
    call i4_to_s_left ( irow1, chrtmp2 )
    call rat_to_s_left ( istop3, isbot3, chrtmp )
    line2 =  chrtmp // ' ' // chrtmp1 // ' ' // chrtmp2
    call s_blanks_delete ( line2 )

  else if ( irow1 == irow3 ) then

    command = 'ROW_ADD'

    if ( isbot2 < 0 ) then
      isbot2 = - isbot2
      istop2 = - istop2
    end if

    call rat_to_s_left ( istop2, isbot2, chrtmp )
    call i4_to_s_left ( irow2, chrtmp1 )
    call i4_to_s_left ( irow1, chrtmp2 )
    line2 = chrtmp // ' ' // chrtmp1 // ' ' // chrtmp2
    call s_blanks_delete ( line2 )

  else
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Your ERO command could not be understood.'
    write ( *, '(a)' ) '  R2 or R3 must equal R1 in an ERO command.'
  end if
 
  return
end
subroutine row_swap ( a, a_int, iatop, iabot, ierror, iform, irow1, &
  irow2, nrow, ncol )

!*****************************************************************************80
!
!! ROW_SWAP swaps two rows of a matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) A(NROW,NCOL), the current matrix.
!
!    Input/output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input/output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL).
!    IATOP and IABOT represent the current rational or decimal matrix.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 for an error.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) IROW1, IROW2, the rows to be swapped.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) irow1
  integer ( kind = 4 ) irow2
  integer ( kind = 4 ) j
  character ( len = 100 ) output
!
!  Skip out if the two rows are the same.
!
  if ( irow1 == irow2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROW_SWAP - Note:'
    write ( *, '(a)' ) '  You have asked to swap a row with itself!'
    return
  end if
!
!  Refuse to continue if a row is out of bounds.
!
  if ( ( irow1 < 1 .or. nrow < irow1 ) .or. &
       ( irow2 < 1 .or. nrow < irow2 ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ROW_SWAP - Error!'
    write ( *, '(a)' ) '  One of the rows is illegal!'
    return
  end if
!
!  Swap the rows.
!
  do j = 1, ncol
 
    if ( iform == 0 ) then
 
      call i4_swap ( iatop(irow1,j), iatop(irow2,j) )
      call i4_swap ( iabot(irow1,j), iabot(irow2,j) )
 
    else if ( iform == 1 ) then
 
      call r4_swap ( a(irow1,j), a(irow2,j) )
 
    else if ( iform == 2 ) then
 
      call i4_swap ( iatop(irow1,j), iatop(irow2,j) )
      call i4_swap ( iabot(irow1,j), iabot(irow2,j) )
 
    else if ( iform == 3 ) then
 
      call i4_swap ( a_int(irow1,j), a_int(irow2,j) )
 
    end if
 
  end do
 
  write ( output, * ) '  ERO: Row ', irow1, ' <=> Row ', irow2
  call s_blanks_delete ( output )
  write ( *, '(a)' ) trim ( output )
 
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
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  iput = 0

  do iget = 1, len ( s )

    c = s(iget:iget)

    if ( c /= ' ' .and. c /= TAB ) then
      iput = iput + 1
      s(iput:iput) = c
    end if

  end do

  s(iput+1:) = ' '

  return
end
subroutine s_blanks_delete ( string )

!*****************************************************************************80
!
!! S_BLANKS_DELETE replaces consecutive blanks by one blank.
!
!  Discussion:
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
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) STRING, the string to be transformed.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nchar
  character newchr
  character oldchr
  character ( len = * ) string
  character, parameter :: TAB = char ( 9 )

  nchar = len ( string )
  j = 0
  newchr = ' '

  do i = 1, nchar

    oldchr = newchr
    newchr = string(i:i)

    if ( newchr == TAB ) then
      newchr = ' '
    end if

    string(i:i) = ' '

    if ( oldchr /= ' ' .or. newchr /= ' ' ) then
      j = j + 1
      string(j:j) = newchr
    end if

  end do

  return
end
subroutine s_cap ( string )

!*****************************************************************************80
!
!! S_CAP replaces any lowercase letters by uppercase ones in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) STRING, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nchar
  character ( len = * ) string

  nchar = len ( string )

  do i = 1, nchar

    c = string(i:i)
    call ch_cap ( c )
    string(i:i) = c

  end do

  return
end
subroutine s_chop ( s, ilo, ihi )

!*****************************************************************************80
!
!! S_CHOP "chops out" a portion of a string, and closes up the hole.
!
!  Example:
!
!    S = 'Fred is not a jerk!'
!
!    call s_chop ( S, 9, 12 )
!
!    S = 'Fred is a jerk!    '
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the locations of the first and last
!    characters to be removed.
!
  implicit none

  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ihi2
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ilo2
  integer ( kind = 4 ) lens
  character ( len = * ) s

  lens = len ( s )

  ilo2 = max ( ilo, 1 )
  ihi2 = min ( ihi, lens )

  if ( ihi2 < ilo2 ) then
    return
  end if

  s(ilo2:lens+ilo2-ihi2-1) = s(ihi2+1:lens)
  s(lens+ilo2-ihi2:lens) = ' '

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
subroutine s_read ( string, line, prompt, ierror, iterm )

!*****************************************************************************80
!
!! S_READ extracts a character string from the input buffer.
!
!  Discussion:
!
!    The routine accepts an input LINE and a PROMPT line.
!
!    If the LINE is empty, the PROMPT is printed and user input
!
!    In either case, enough characters are read from LINE to fill
!    STRING and the positions read are removed.
!
!    PROMPT is also updated.  On satisfactory input of STRING,
!    everything in PROMPT up to and including the first comma is removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) STRING.
!    The user's response to the PROMPT, as read from LINE.
!
!    Input/output, character ( len = 80 ) LINE.
!    A buffer containing the user's input.
!
!    Input/output, character ( len = 80 ) PROMPT.
!    On input, a prompt string that will be printed if necessary.
!
!    On output, if STRING has been read, then PROMPT is cleared out
!    up to, and including, the first comma.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 for an error.
!
!    Input, integer ( kind = 4 ) ITERM,
!    0, No check for terminators.
!    1, Blank, slash, comma, semicolon, equals, greater or
!       lesser signs terminate input.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) lchar
  character ( len = 80 ) line
  integer ( kind = 4 ) nchar
  character null
  character ( len = 80 ) prompt
  character ( len = * ) string

  ierror = 0
  null = char(0)
  string = ' '

  call chrinp ( ierror, line, prompt )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Null input acceptable for character input only.
!
  if ( line == ' ' ) return
  lchar = 0
  nchar = len ( string )
 
  do i = 1, nchar
 
    if ( lchar /= 0 ) then
      exit
    end if
 
    c = line(i:i)
 
    if ( iterm == 1 ) then
 
      if ( c == ' ' .or. c == null .or. c == '/' .or. c == ',' .or. &
           c == ';' .or. c == '=' ) then
        lchar = i
      end if
 
    end if
 
    if ( lchar == 0 ) then
      string(i:i) = c
    end if
 
  end do
!
!  Chop out the character positions that have been used.
!
  if ( lchar == 0 ) then
    lchar = nchar
  end if

  call s_chop ( line, 1, lchar )
!
!  Force the string to be flush left by removing leading blanks.
!
  line = adjustl ( line )

  return
end
subroutine s_to_dec ( s, itop, ibot, length )

!*****************************************************************************80
!
!! S_TO_DEC reads a number from a string, returning a decimal result.
!
!  Discussion:
!
!    The integer may be in real format, for example '2.25'.  It
!    returns ITOP and IBOT.  If the input number is an integer, ITOP
!    equals that integer, and IBOT is 1.  But in the case of 2.25,
!    the program would return ITOP = 225, IBOT = 100.
!
!    Legal input is
!
!          blanks,
!       2  initial sign,
!          blanks,
!       3  whole number,
!       4  decimal point,
!       5  fraction,
!       6  'E' or 'e' or 'D' or 'd', exponent marker,
!       7  exponent sign,
!       8  exponent,
!          blanks
!       9  comma or semicolon
!      10end of information
!
!  Example:
!
!    S                 ITOP      IBOT     Length  Meaning
!
!    '1'                  1         0          1        1
!    '     1   '          1         0          6        1
!    '1A'                 1         0          1        1
!    '12,34,56'          12         0          3       12
!    '  34 7'            34         0          4       34
!    '-1E2ABCD'          -1         2          4     -100
!    '-1X2ABCD'          -1         0          2       -1
!    ' 2E-1'              2        -1          5      0.2
!    '23.45'           2345        -2          5    23.45
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading begins at position 1 and
!    terminate when no more characters
!    can be read to form a legal integer.  Blanks, commas,
!    or other nonnumeric data will, in particular, cause
!    the conversion to halt.
!
!    Output, integer ( kind = 4 ) ITOP, the integer read from the string,
!    assuming that no negative exponents or fractional parts
!    were used.  Otherwise, the 'integer' is ITOP/IBOT.
!
!    Output, integer ( kind = 4 ) IBOT, the integer divisor required to
!    represent numbers which are in real format or have a
!    negative exponent.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters used.
!
  implicit none

  logical ch_is_digit
  character c
  integer ( kind = 4 ) digit
  integer ( kind = 4 ) exponent
  integer ( kind = 4 ) exponent_sign
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) length
  integer ( kind = 4 ) mantissa_sign
  character ( len = * ) s
  logical s_eqi

  itop = 0
  ibot = 0

  if ( len ( s ) <= 0 ) then
    length = 0
    return
  end if

  length = - 1
  exponent_sign = 0
  mantissa_sign = 1
  exponent = 0
  ihave = 1
  iterm = 0
!
!  Consider the next character in the string.
!
  do

    length = length + 1
    c = s(length+1:length+1)
!
!  Blank.
!
    if ( c == ' ' ) then

      if ( ihave == 1 ) then

      else if ( ihave == 2 ) then

      else
        iterm = 1
      end if
!
!  Comma or semicolon.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 9
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        mantissa_sign = - 1
      else if ( ihave == 6 ) then
        ihave = 7
        exponent_sign = -1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
        mantissa_sign = +1
      else if ( ihave == 6 ) then
        ihave = 7
        exponent_sign = +1
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else
        iterm = 1
      end if
!
!  Exponent marker.
!
    else if ( s_eqi ( c, 'E' ) .or. s_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if ( ch_is_digit ( c ) ) then

      if ( ihave <= 3 ) then
        ihave = 3
        call ch_to_digit ( c, digit )
        itop = 10 * itop + digit
      else if ( ihave <= 5 ) then
        ihave = 5
        call ch_to_digit ( c, digit )
        itop = 10 * itop + digit
        ibot = ibot - 1
      else if ( ihave <= 8 ) then
        ihave = 8
        call ch_to_digit ( c, digit )
        exponent = 10 * exponent + digit
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'S_TO_DEC: Fatal error!'
        write ( *, '(a,i6)' ) '  IHAVE = ', ihave
        stop
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if

    if ( iterm == 1 ) then
      exit
    end if

    if ( len ( s ) <= length + 1 ) then
      length = len ( s )
      exit
    end if

  end do
!
!  Number seems to have terminated.
!  Have we got a legal number?
!
  if ( ihave == 1 ) then
    return
  else if ( ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_DEC - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) trim ( s )
    return
  end if
!
!  Normalize.
!
  if ( 0 < itop ) then

    do while ( mod ( itop, 10 ) == 0 )
      itop = itop / 10
      ibot = ibot + 1
    end do

  end if
!
!  Consolidate the number in the form ITOP * 10**IBOT.
!
  ibot = ibot + exponent_sign * exponent
  itop = mantissa_sign * itop

  if ( itop == 0 ) then
   ibot = 0
  end if

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
!    Output, integer ( kind = 4 ) LAST, the last character of S used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
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
        last = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    last = len_trim ( s )
  else
    ierror = 1
    last = 0
  end if

  return
end
subroutine s_to_r4 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R4 reads an R4 from a string.
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
!  Example:
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Output, real ( kind = 4 ) R, the real value that was read from the string.
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
!    Output, integer ( kind = 4 ) LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical ch_eqi
  character c
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 4 ) r
  real ( kind = 4 ) rbot
  real ( kind = 4 ) rexp
  real ( kind = 4 ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0E+00
  lchar = - 1
  isgn = 1
  rtop = 0.0E+00
  rbot = 1.0E+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    lchar = lchar + 1
    c = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
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
        lchar = lchar + 1
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
      else if ( 6 <= ihave .and. ihave <= 8 ) then
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
        rtop = 10.0E+00 * rtop + real ( ndig )
      else if ( ihave == 5 ) then
        rtop = 10.0E+00 * rtop + real ( ndig )
        rbot = 10.0E+00 * rbot
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
    if ( iterm == 1 .or. nchar <= lchar + 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
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
    rexp = 1.0E+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0E+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0E+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

  return
end
subroutine size_set ( ierror, line, maxrow, nrow, title )

!*****************************************************************************80
!
!! SIZE_SET allows the user to specify the size of a new problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 for an error.
!
!    Workspace, character ( len = 80 ) LINE, used to hold the user's input.
!
!    Input, integer ( kind = 4 ) MAXROW, the maximum number of matrix rows.
!
!    Output, integer ( kind = 4 ) NROW, the number of rows in the matrix.
!
!    Input, character ( len = * ) TITLE, contains 'rows' or 'columns'.
!
  implicit none

  integer ( kind = 4 ) maxrow

  integer ( kind = 4 ) ierror
  character ( len = 80 ) line
  integer ( kind = 4 ) nrow
  character ( len = 80 ) prompt
  character ( len = * ) title

  ierror = 0
  nrow = 0
 
  prompt = 'number of ' // trim ( title )

  call i4_read ( nrow, line, prompt, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SIZE_SET - Fatal error!'
    write ( *, '(a)' ) '  I4_READ returned error flag.'
    return
  end if
 
  if ( nrow < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SIZE_SET - Error!'
    write ( *, '(a)' ) '  Your value is less than 1!'
    ierror = 1
    return
  end if

  if ( maxrow < nrow ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SIZE_SET - Error!'
    write ( *, '(a,i6)' ) '  Your value is greater than ', maxrow
    ierror = 1
    return
  end if

  write ( *, '(a,i6)' ) 'Set number of ' // trim ( title ) // ' to ', nrow

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
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
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
subroutine value_read ( irow, icol, base, iform, nrow, ncol, line, a, a_int, &
  iatop, iabot, ierror )

!*****************************************************************************80
!
!! VALUE_READ allows the user to specify one entry in the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IROW, ICOL, the row and column of the entry.
!
!    Input, integer ( kind = 4 ) BASE, the base for modular arithmetic.
!
!    Input, integer ( kind = 4 ) IFORM, 0/1/2/3/4, rat/real/dec/int/mod arithmetic.
!
!    Input, integer ( kind = 4 ) BASE, the base of modular arithmetic.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in the matrix.
!
!    Workspace, character ( len = 80 ) LINE, used to hold the user's input.
!
!    Input/output, real ( kind = 4 ) A(NROW,NCOL), the real matrix.
!
!    Input/output, integer ( kind = 4 ) A_INT(NROW,NCOL), the current 
!    integer or integer mod matrix.
!
!    Input/output, integer ( kind = 4 ) IATOP(NROW,NCOL), IABOT(NROW,NCOL), the 
!    rat/dec matrix.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 for an error.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  real ( kind = 4 ) a(nrow,ncol)
  integer ( kind = 4 ) a_int(nrow,ncol)
  integer ( kind = 4 ) base
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) iabot(nrow,ncol)
  integer ( kind = 4 ) iatop(nrow,ncol)
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iform
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) isbot
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) ival
  character ( len = 80 ) line
  character ( len = 80 ) prompt
  real ( kind = 4 ) rval

  call i4_to_s_left ( irow, chrtmp1 )
  call i4_to_s_left ( icol, chrtmp2 )

  prompt = 'value of A(' // trim ( chrtmp1 ) // ',' // trim ( chrtmp2 ) // ')'

  if ( iform == 0 ) then

    call rat_read ( istop, isbot, line, prompt, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'VALUE_READ - Error!'
      write ( *, '(a)' ) '  RAT_READ returned error flag!'
      return
    end if

    itemp = i4_gcd ( istop, isbot )
    iatop(irow,icol) = istop / itemp
    iabot(irow,icol) = isbot / itemp

  else if ( iform == 1 ) then

    call r4_read ( rval, line, prompt, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'VALUE_READ - Error!'
      write ( *, '(a)' ) '  r4_read returned error flag!'
      return
    end if

    a(irow,icol) = rval

  else if ( iform == 2 ) then

    call dec_read ( istop, isbot, line, prompt, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'VALUE_READ - Error!'
      write ( *, '(a)' ) '  DEC_READ returned error flag!'
      return
    end if

    call dec_round ( istop, isbot )

    iatop(irow,icol) = istop
    iabot(irow,icol) = isbot

   else if ( iform == 3 ) then

    call i4_read ( ival, line, prompt, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'VALUE_READ - Error!'
      write ( *, '(a)' ) '  I4_READ returned error flag!'
      return
    end if

    a_int(irow,icol) = ival

   else if ( iform == 4 ) then

    call i4_read ( ival, line, prompt, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'VALUE_READ - Error!'
      write ( *, '(a)' ) '  I4_READ returned error flag!'
      return
    end if

    ival = mod ( ival, base )
    a_int(irow,icol) = ival

  end if

  return
end
