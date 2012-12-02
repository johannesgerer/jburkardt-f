program main

!*****************************************************************************80
!
!! MAIN is the main program for ISING_2D_SIMULATION.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: c1
  integer ( kind = 4 ) iterations
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8) , dimension ( 1:5 ) :: prob
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) step
  real ( kind = 8 ) thresh

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ISING_2D_SIMULATION'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Monte Carlo simulation of a 2D Ising model.'
!
!  Get input from commandline or user.
!
  call get_input ( m, n, iterations, thresh, seed )

  prob = (/ 0.98D+00, 0.85D+00, 0.50D+00, 0.15D+00, 0.02D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of rows is M = ', m
  write ( *, '(a,i8)' ) '  The number of columns is N = ', n
  write ( *, '(a,i8)' ) '  The number of iterations taken is ITERATIONS = ', iterations
  write ( *, '(a,g14.6)' ) '  The threshhold THRESH = ', thresh
  write ( *, '(a,i12)' ) '  The seed SEED = ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The transition probability table, based on the number of'
  write ( *, '(a)' ) '  neighbors with the same spin.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      1         2         3         4         5'
  write ( *, '(a)' ) ' '
  write ( *, '(7f10.4)' ) prob(1:5)
!
!  Initialize the system.
!
  allocate ( c1(m,n) )

  call ising_2d_initialize ( m, n, thresh, seed, c1 )
!
!  Do the simulation.
!
  call transition ( m, n, iterations, prob, seed, c1 )
!
!  Free memory.
!
  deallocate ( c1 )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ISING_2D_SIMULATION'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
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
subroutine ch_to_digit ( ch, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!     CH  DIGIT
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
!    Input, character CH, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.
!    If CH was 'illegal', then DIGIT is -1.
!
  implicit none

  character ch
  integer ( kind = 4 ) digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then

    digit = iachar ( ch ) - 48

  else if ( ch == ' ' ) then

    digit = 0

  else

    digit = - 1

  end if

  return
end
subroutine get_input ( m, n, iterations, thresh, seed )

!*****************************************************************************80
!
!! GET_INPUT gets input parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Output, integer ( kind = 4 ) ITERATIONS, the number of iterations.
!
!    Output, real ( kind = 8 ) THRESH, the threshhold.
!
!    Output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
  implicit none

  integer ( kind = 4 ) arg_num
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iterations
  integer ( kind = 4 ) last
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  character ( len = 80 ) string
  real ( kind = 8 ) thresh

  arg_num = iargc ( )

  if ( 1 <= arg_num ) then
    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, m, ierror, last )
  else if ( .true. ) then
    m = 10
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter M, the number of rows:'
    read ( *, * ) m
  end if

  if ( 2 <= arg_num ) then
    iarg = 2
    call getarg ( iarg, string )
    call s_to_i4 ( string, n, ierror, last )
  else if ( .true. ) then
    n = 10
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter N, the number of columns:'
    read ( *, * ) n
  end if

  if ( 3 <= arg_num ) then
    iarg = 3
    call getarg ( iarg, string )
    call s_to_i4 ( string, iterations, ierror, last )
  else if ( .true. ) then
    iterations = 15
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the number of iterations to take.'
    read ( *, * ) iterations
  end if

  if ( 4 <= arg_num ) then
    iarg = 4
    call getarg ( iarg, string )
    call s_to_r8 ( string, thresh, ierror, last )
  else if ( .true. ) then
    thresh = 0.50D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the threshhold.'
    read ( *, * ) thresh
  end if

  if ( 5 <= arg_num ) then
    iarg = 5
    call getarg ( iarg, string )
    call s_to_i4 ( string, seed, ierror, last )
  else if ( .true. ) then
    seed = 123456789
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the random number seed.'
    read ( *, * ) seed
  end if

  return
end
subroutine ising_2d_agree ( m, n, c1, c5 )

!*****************************************************************************80
!
!! ISING_2D_AGREE returns the number of neighbors agreeing with each cell.
!
!  Discussion:
!
!    The count includes the cell itself, so it is between 1 and 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of cells in each 
!    spatial dimension.
!
!    Input, integer ( kind = 4 ) C1(M,N), an array of 1's and -1's.
!
!    Output, integer ( kind = 4 ) C5(M,N), the number of neighbors 
!    that agree.  1, 2, 3, 4, or 5.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) c1(m,n)
  integer ( kind = 4 ) c5(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  c5 = c1 &
     + cshift ( c1, -1,  1 ) &
     + cshift ( c1, +1,  1 ) &
     + cshift ( c1, -1,  2 ) &
     + cshift ( c1, +1,  2 )

  do j = 1, n
    do i = 1, m
      if ( 0 < c1(i,j) ) then
        c5(i,j) = ( 5 + c5(i,j) ) / 2
      else
        c5(i,j) = ( 5 - c5(i,j) ) / 2
      end if
    end do
  end do

  return
end
subroutine ising_2d_initialize ( m, n, thresh, seed, c1 )

!*****************************************************************************80
!
!! ISING_2D_INITIALIZE initializes the Ising array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) THRESH, the threshhold.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, integer ( kind = 4 ) C1(M,N), the initial Ising array.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) c1(m,n)
  real ( kind = 8 ) r(m,n)
  real ( kind = 8 ) thresh
  integer ( kind = 4 ) seed

  call r8mat_uniform_01 ( m, n, seed, r )

  where ( r <= thresh )
    c1 = -1
  elsewhere
    c1 = +1
  end where

  return
end
subroutine ising_2d_stats ( step, m, n, c1 )

!*****************************************************************************80
!
!! ISING_2D_STATS prints information about the current step.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) STEP, the step number.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) C1(M,N), the current state of the system.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) c1(m,n)
  integer ( kind = 4 ) pos_count
  real ( kind = 8 ) pos_percent
  integer ( kind = 4 ) step
  integer ( kind = 4 ) neg_count
  real ( kind = 8 ) neg_percent

  if ( step == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Step     Positives       Negatives'
    write ( *, '(a)' ) '             #    %          #    %'
    write ( *, '(a)' ) ' '
  end if

  pos_count = sum ( c1**2 + c1 ) / 2
  neg_count = m * n - pos_count
  pos_percent = real ( 100 * pos_count, kind = 8 ) / real ( m * n, kind = 8 )
  neg_percent = real ( 100 * neg_count, kind = 8 ) / real ( m * n, kind = 8 )

  write ( *, '(2x,i4,2x,i6,2x,f6.2,2x,i6,2x,f6.2)' ) &
    step, pos_count, pos_percent, neg_count, neg_percent 

  return
end
subroutine neighbor_2d_stats ( step, m, n, c1, c5 )

!*****************************************************************************80
!
!! NEIGHBOR_2D_STATS prints neighbor statistics about the current step.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) STEP, the step number.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) C1(M,N), the current state of the system.
!
!    Input, integer ( kind = 4 ) C5(M,N), the number of agreeable neighbors.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) c1(m,n)
  integer ( kind = 4 ) c5(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) stats(-5:5)
  integer ( kind = 4 ) step

  if ( step == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Step     Neighborhood Charge:'
    write ( *, '(a)' ) &
      '           -5    -4    -3    -2    -1    +1    +2    +3    +4    +5'
    write ( *, '(a)' ) ' '
  end if
  
  stats(-5:5) = 0
  do j = 1, n
    do i = 1, m
      stats(c5(i,j)) = stats(c5(i,j)) + 1
    end do
  end do
  write ( *, '(2x,i4,10(2x,i4))' ) step, stats(-5:-1), stats(1:5)
  
  return
end
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
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
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

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

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) length
  character ( len = * ) s
  integer ( kind = 4 ) state
  character :: TAB = achar ( 9 )
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
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 value from a string.
!
!  Discussion:
!
!    An "R8" value is simply a real number to be stored as a
!    variable of type "real ( kind = 8 )".
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
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
!    12 January 2009
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
!    Output, real ( kind = 8 ) DVAL, the value read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character c
  logical ch_eqi
  real ( kind = 8 ) dval
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) length
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * )  s
  integer ( kind = 4 ) s_length
  character  :: TAB = achar ( 9 )

  s_length = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( s_length < length + 1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
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
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
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
!  Scientific notation exponent marker.
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
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

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
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to S_LENGTH.
!
  if ( iterm /= 1 .and. length + 1 == s_length ) then
    length = s_length
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
        / real ( jbot, kind = 8 ) )
    end if
  end if

  dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

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
subroutine transition ( m, n, iterations, prob, seed, c1 )

!*****************************************************************************80
!
!! TRANSITION carries out a Monte Carlo simulation of a 3D Ising model.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) ITERATIONS, the number of iterations.
!
!    Input, real ( kind = 8 ) PROB(1:5).  PROB(I) represents the probability 
!    that the spin of a given cell will be reversed, given that it has I 
!    immediate neighbors (including itself) with spin the same as its own.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number 
!    generator.
!
!    Input/output, integer ( kind = 4 ) C1(M,N).  On input, the current 
!    state of the system.  On output, the state of the system after 
!    the iterations.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) c1(m,n)
  integer ( kind = 4 ) c5(m,n)
  integer ( kind = 4 ) iterations
  integer ( kind = 4 ) j
  real ( kind = 8 ) prob(1:5)
  integer ( kind = 4 ) step
  real ( kind = 8 ) r(m,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) threshhold(m,n)

  step = 0
  call ising_2d_stats ( step, m, n, c1 )

  do step = 1, iterations
!
!  C5 contains 1 through 5, the number of cells that agree with the center cell.
!
    call ising_2d_agree ( m, n, c1, c5 )

    if ( .false. ) then
      call neighbor_2d_stats ( step, m, n, c1, c5 )
    end if
!
!  Determine the chances of flipping cell (I,J).
!
    do j = 1, 5
      where ( c5 == j ) 
        threshhold = prob(j)
      end where
    end do

    call r8mat_uniform_01 ( m, n, seed, r )

    where ( r < threshhold ) 
      c1 = - c1
    endwhere

    call ising_2d_stats ( step, m, n, c1 )

  end do

  return
end
