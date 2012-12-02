program main

!*****************************************************************************80
!
!! MAIN is the main program for TABLE_ORTHONORMALIZE.
!
!  Discussion:
!
!    TABLE_ORTHONORMALIZE orthonormalizes the vectors in a TABLE file.
!
!    TABLE_ORTHONORMALIZE is a FORTRAN90 program that can read in a
!    set of vectors and orthonormalize them.
!
!    For instance, the program can read three (or in general N) 
!    files, each containing one of the M-dimensional vectors
!    V1, V2 and V3.  The program then outputs three new vectors,
!    Q1, Q2 and Q3, with the properties that 
!
!    * The Euclidean norm of Q(I) is 1, for all I;
!    * The dot product < Q(I), 0 < Q(J), for I and J distinct;
!    *  The vectors Q span the same space as the original vectors V.
!
!    The algorithm collects the N input vectors V into an M by N 
!    matrix.  It then calls the appropriate LAPACK routines to compute 
!    the QR factorization:
!
!       V = Q * R
!
!    and outputs the individual columns of Q as the new vectors.
!
!    Some effort has been made to relieve the user of the tedious
!    necessity of specifying how many input files there are (the number
!    N), and what their dimensionality is (the number M).
!    It is assumed that the input data files have "consecutive" names,
!    such as file01.txt, file02.txt, ... and so on.  
!    This is enough to figure out N.  And the program counts 
!    the number of lines in the first file, and assumes this value
!    of M may safely be applied to all the files.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 February 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  character ( len = 255 ) input_name
  logical input_name_exist
  character ( len = 255 ) input_name_first
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m_half
  integer ( kind = 4 ) n
  integer ( kind = 4 ) numarg
  character ( len = 255 ) output_name
  character ( len = 255 ) output_name_first
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: v
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: vtv
  real ( kind = 8 ), allocatable, dimension ( : ) :: vx
  real ( kind = 8 ), allocatable, dimension ( : ) :: vy

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TABLE_ORTHONORMALIZE:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Orthonormalize a set of N vectors in M space.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Given a set of vectors V1, V2, ..., VN,'
  write ( *, '(a)' ) '  produce a second set Q1, Q2, ..., QN,'
  write ( *, '(a)' ) '  with the properties that:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  * each Q has unit Euclidean norm;'
  write ( *, '(a)' ) '  * distinct Q(I) and Q(J) are orthogonal;'
  write ( *, '(a)' ) '  * the Q''s span the same space as the V''s.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Each V is assumed to be stored in a separate file.'
  write ( *, '(a)' ) '  The file names are assumed to be "numerically" ordered'
  write ( *, '(a)' ) '  such as "v_01.txt", "v_02.txt" and so on.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The output of the program will be a similar family'
  write ( *, '(a)' ) '  of files containing the Q''s.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  To invoke the program, simply list the FIRST input'
  write ( *, '(a)' ) '  file for V and output file for Q, such as:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    table_orthonormalize  v_01.txt  q_01.txt'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  OK, here we go.'

  numarg = iargc ( )
!
!  Get the first input filename.
!
  if ( numarg < 1 ) then

    call s_input ( 'Enter the name of the first input file', &
      input_name_first, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TABLE_ORTHONORMALIZE - Fatal error!'
      write ( *, '(a)' ) '  Error reading user input.'
      stop
    end if

  else

    iarg = 1
    call getarg ( iarg, input_name_first )

  end if
!
!  Get the first output filename.
!
  if ( numarg < 2 ) then

    call s_input ( 'Enter the name of the first output file', &
      output_name_first, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TABLE_ORTHONORMALIZE - Fatal error!'
      write ( *, '(a)' ) '  Error reading user input.'
      stop
    end if

  else

    iarg = 2
    call getarg ( iarg, output_name_first )

  end if
!
!  Count the number of lines in the file to get the spatial dimension M.
!
  call xy_read_header ( input_name_first, m_half )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of data lines in file is ', m_half
  m = 2 * m_half
  write ( *, '(a,i8)' ) '  Spatial dimension M = ', m
!
!  Count the number of basis files to get N, the number of basis vectors.
!
  input_name = input_name_first

  n = 0

  do

    inquire ( file = input_name, exist = input_name_exist )

    if ( .not. input_name_exist ) then
      exit
    end if

    n = n + 1
    call file_name_inc ( input_name )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of basis vectors is N = ', n
!
!  Knowing M and N, we can allocate space.
!
  allocate ( v(m,n) )
  allocate ( vx(1:m_half) )
  allocate ( vy(1:m_half) )
!
!  Read the basis vectors into the basis vector matrix U.
!
  input_name = input_name_first

  do j = 1, n
    call xy_read_data ( input_name, m_half, vx, vy )
    v(1:m-1:2,j) = vx(1:m_half)
    v(2:m  :2,j) = vy(1:m_half)
    call file_name_inc ( input_name )
  end do
!
!  Compute the Q matrix.
!
  call q_factor ( m, n, v )
!
!  Perform orthonormality test.
!
  allocate ( vtv(n,n) )

  vtv = matmul ( transpose ( v ), v )

  call r8mat_print ( n, n, vtv, '  V''V (should be I)' )

  deallocate ( vtv )
!
!  Output the Q vectors.
!
  output_name = output_name_first

  do j = 1, n
    vx(1:m_half) = v(1:m-1:2,j)
    vy(1:m_half) = v(2:m  :2,j)
    call xy_write ( output_name, m_half, vx, vy, ierror )
    call file_name_inc ( output_name )
  end do
!
!  Free memory.
!
  deallocate ( v )
  deallocate ( vx )
  deallocate ( vy )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TABLE_ORTHONORMALIZE:'
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
function ch_is_digit ( c )

!*****************************************************************************80
!
!! CH_IS_DIGIT returns TRUE if a character is a decimal digit.
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
!    Output, logical CH_IS_DIGIT, TRUE if C is a digit, FALSE otherwise.
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
!! CH_TO_DIGIT returns the integer ( kind = 4 ) value of a base 10 digit.
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
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer ( kind = 4 ) value.  If C was
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
subroutine digit_inc ( c )

!*****************************************************************************80
!
!! DIGIT_INC increments a decimal digit.
!
!  Example:
!
!    Input  Output
!    -----  ------
!    '0'    '1'
!    '1'    '2'
!    ...
!    '8'    '9'
!    '9'    '0'
!    'A'    'A'
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
!    Input/output, character C, a digit to be incremented.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  call ch_to_digit ( c, digit )

  if ( digit == -1 ) then
    return
  end if

  digit = digit + 1

  if ( digit == 10 ) then
    digit = 0
  end if

  call digit_to_ch ( digit, c )

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
subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC generates the next filename in a series.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected, and
!    if the name contains no digits, then nothing is done.
!
!  Example:
!
!      Input          Output
!      -----          ------
!      a7to11.txt     a7to12.txt
!      a7to99.txt     a8to00.txt
!      a9to99.txt     a0to00.txt
!      cat.txt        cat.txt
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
!    Input/output, character ( len = * ) FILE_NAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  logical ch_is_digit
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( file_name )

  do i = lens, 1, -1

    c = file_name(i:i)

    if ( ch_is_digit ( c ) ) then

      call digit_inc ( c )

      file_name(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do

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
subroutine q_factor ( m, n, a )

!*****************************************************************************80
!
!! Q_FACTOR determines the Q factor of a matrix.
!
!  Discussion:
!
!    DGEQRF implicitly computes the QR factorization of an M by N 
!    matrix A:
!
!      A(MxN) = Q(MxK) * R(KxN)
!
!    where K = min ( M, N ).  For our purposes, it should always
!    be the case that N < M, so that K = N.
!
!    DORGQR explicitly forms the Q matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Anderson, Bai, Bischof, Demmel, Dongarra, Du Croz, Greenbaum,
!    Hammarling, McKenney, Ostrouchov, Sorensen,
!    LAPACK User's Guide,
!    Second Edition,
!    SIAM, 1995,
!    QA76.73.F25L36
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in 
!    the matrix.  It must be the case that N < M.
!
!    Input/output, real ( kind = 8 ) A(M,N).  On input, a matrix to analyze.
!    On output, the factor Q in the QR factorization of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) lwork
  real ( kind = 8 ) tau(n)
  real ( kind = 8 ) work(n)

  if ( m < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Q_FACTOR - Fatal error!'
    write ( *, '(a)' ) '  M < N.'
    write ( *, '(a,i8)' ) '  M = ', m
    write ( *, '(a,i8)' ) '  N = ', n
    stop
  end if
!
!  Compute the QR factorization.
!
  k = n
  lda = m
  lwork = n

  call dgeqrf ( m, n, a, lda, tau, work, lwork, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Q_FACTOR - Warning:'
    write ( *, '(a,i8)' ) '  DGEQRF returned nonzero INFO = ', info
    stop
  end if

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The diagonal of the R matrix gives some'
    write ( *, '(a)' ) '  information about the condition of the matrix.'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,i8,2x,g14.6)' ) i, a(i,i)
    end do
  end if
!
!  Construct Q explicitly.
!
  call dorgqr ( m, n, k, a, lda, tau, work, lwork, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Q_FACTOR - Warning:'
    write ( *, '(a,i8)' ) '  DORGQR returned nonzero INFO = ', info
    stop
  end if

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 2004
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
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
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
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine s_input ( string, value, ierror )

!*****************************************************************************80
!
!! S_INPUT prints a prompt string and reads a string from the user.
!
!  Discussion:
!
!    If the input line starts with a comment character ('#'), or is blank,
!    the routine ignores that line, and tries to read the next one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the prompt string.
!
!    Output, character ( len = * ) VALUE, the value input by the user.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, 
!    which is zero if no error occurred.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = 80 ) line
  character ( len = * ) string
  character ( len = * ) value

  ierror = 0
  value = ' '
!
!  Write the prompt.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  do

    read ( *, '(a)', iostat = ierror ) line

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'S_INPUT: Fatal error!'
      write ( *, '(a)' ) '  Input error!'
      stop
    end if
!
!  If the line begins with a comment character, go back and read the next line.
!
    if ( line(1:1) == '#' ) then
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    value = line

    exit

  end do

  return
end
subroutine s_to_r8 ( s, r, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer ( kind = 4 ) part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer ( kind = 4 ) part,
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
!    Output, double precision R, the value that was read from the string.
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
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character c
  logical ch_eqi
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) length
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) r
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0D+00
  length = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1
    c = s(length+1:length+1)
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
        length = length + 1
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
        rtop = 10.0D+00 * rtop + real ( ndig )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig )
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
    if ( iterm == 1 .or. nchar <= length+1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to NCHAR.
!
  if ( iterm /= 1 .and. length+1 == nchar ) then
    length = nchar
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
    rexp = 1.0D+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0D+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0D+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

  return
end
subroutine s_to_r8vec ( s, n, rvec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2001
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
!    Output, real ( kind = 8 ) RVEC(N), the values read from the string.
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
  integer ( kind = 4 ) length
  real ( kind = 8 ) rvec(n)
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, length )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + length

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
subroutine xy_read_data ( file_in_name, point_num, x, y )

!*****************************************************************************80
!
!! XY_READ_DATA reads data from an XY file.
!
!  Discussion:
!
!    The number of points in the file can be determined by calling
!    XY_READ_HEADER first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_IN_NAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.  The program
!    will stop reading data once POINT_NUM values have been read.
!
!    Output, real ( kind = 8 ) X(POINT_NUM), Y(POINT_NUM), the XY 
!    point coordinates.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) point_num

  character ( len = * ) file_in_name
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 100 ) line
  real ( kind = 8 ) temp(dim_num)
  real ( kind = 8 ) x(point_num)
  real ( kind = 8 ) y(point_num)

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XY_READ_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( file_in_name )
    stop
  end if

  i = 0

  do while ( i < point_num )

    read ( file_in_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = i
      exit
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( line, dim_num, temp, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    i = i + 1

    x(i) = temp(1)
    y(i) = temp(2)

  end do

  close ( unit = file_in_unit )

  return
end
subroutine xy_read_header ( file_in_name, point_num )

!*****************************************************************************80
!
!! XY_READ_HEADER determines the number of pairs of data in an XY file.
!
!  Discussion:
!
!    This routine assumes that the file contains exactly three kinds of
!    records:
!
!    COMMENTS which begin with a '#' character in column 1;
!    BLANKS which contain nothing but 'whitespace';
!    XY coordinates, which each contain one pair of real values.
!
!    The routine ignores comments and blank lines and returns
!    the number of lines containing XY coordinates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_IN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the number of points in the file.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) point_num

  character ( len = * ) file_in_name
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 100 ) line
  real ( kind = 8 ) temp(dim_num)

  point_num = 0

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XY_READ_HEADER - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file "' // &
      trim ( file_in_name ) // '".'
    stop
  end if

  do

    read ( file_in_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = 1
      exit
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( line, dim_num, temp, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    point_num = point_num + 1

  end do

  close ( unit = file_in_unit )

  return
end
subroutine xy_write ( file_out_name, point_num, x, y, ierror )

!*****************************************************************************80
!
!! XY_WRITE writes an XY file.
!
!  Example:
!
!    # my_file.xy
!    # created by ORTHONORMALIZE::XY_WRITE.
!    #
!    #  Number of points = 5
!    #
!    0.0  0.0
!    1.0  2.0
!    3.0  5.0
!    2.0  1.0
!    8.0  7.5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2003
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
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(POINT_NUM), Y(POINT_NUM), the XY coordinates.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the file could not be opened.
!
  implicit none

  integer ( kind = 4 ) point_num

  logical, parameter :: debug = .false.
  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  real ( kind = 8 ) x(point_num)
  real ( kind = 8 ) y(point_num)

  ierror = 0
!
!  Open the file.
!
  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XY_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' // &
      trim ( file_out_name ) // '".'
    ierror = 1
    return
  end if
!
!  Write the header.
!
  call xy_write_header ( file_out_name, file_out_unit, point_num )
!
!  Write the data.
!
  call xy_write_data ( file_out_unit, point_num, x, y )
!
!  Close the file.
!
  close ( unit = file_out_unit )
!
!  Report
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'XY_WRITE - Note:'
    write ( *, '(a)' ) '  The data was written.'
    write ( *, '(a,i8)' ) '  Number of points =    ', point_num
  end if

  return
end
subroutine xy_write_data ( file_out_unit, point_num, x, y )

!*****************************************************************************80
!
!! XY_WRITE_DATA writes the data of an XY file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FILE_OUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) X(POINT_NUM), Y(POINT_NUM), the point coordinates.
!
  implicit none

  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(point_num)
  real ( kind = 8 ) y(point_num)

  do j = 1, point_num
    write ( file_out_unit, '(2f24.16)' ) x(j), y(j)
  end do

  return
end
subroutine xy_write_header ( file_out_name, file_out_unit, point_num )

!*****************************************************************************80
!
!! XY_WRITE_HEADER writes the header of an XY file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2003
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
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
  implicit none

  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) point_num
!
!  Write the header.
!
  write ( file_out_unit, '(a)'    ) '#  ' // trim ( file_out_name )
  write ( file_out_unit, '(a)'    ) '#  created by ORTHONORMALIZE::XY_WRITE'
  write ( file_out_unit, '(a)' )    '#'
  write ( file_out_unit, '(a,i8)' ) '#  Number of points = ', point_num
  write ( file_out_unit, '(a)' )    '#'

  return
end

