program main

!*****************************************************************************80
!
!! MAIN is the main program for NINT_EXACTNESS_MIXED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real      ( kind = 8 ), allocatable, dimension ( : ) :: alpha
  integer   ( kind = 4 ) arg_num
  real      ( kind = 8 ), allocatable, dimension ( : ) :: beta
  integer   ( kind = 4 ) degree
  integer   ( kind = 4 ) degree_max
  integer   ( kind = 4 ) dim
  integer   ( kind = 4 ) dim_num
  integer   ( kind = 4 ) dim_num2
  integer   ( kind = 4 ), allocatable, dimension ( : ) :: expon
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) iarg
  integer   ( kind = 4 ) iargc
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) last
  logical                more
  integer   ( kind = 4 ) point_num
  integer   ( kind = 4 ) point_num2
  real      ( kind = 8 ) quad_error
  character ( len = 255 ) quad_filename
  character ( len = 255 ) quad_a_filename
  character ( len = 255 ) quad_b_filename
  character ( len = 255 ) quad_r_filename
  character ( len = 255 ) quad_w_filename
  character ( len = 255 ) quad_x_filename
  real      ( kind = 8 ) r8_huge
  integer   ( kind = 4 ), allocatable, dimension ( : ) :: rule
  character ( len = 255 ) string
  integer   ( kind = 4 ) t
  real      ( kind = 8 ), allocatable, dimension ( : ) :: weight
  real      ( kind = 8 ), allocatable, dimension ( :, : ) :: x
  real      ( kind = 8 ), allocatable, dimension ( :, : ) :: x_range

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NINT_EXACTNESS_MIXED'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Investigate the polynomial exactness of'
  write ( *, '(a)' ) '  a multidimensional quadrature rule'
  write ( *, '(a)' ) '  for a region R = R1 x R2 x ... x RM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Individual rules may be for:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Legendre:'
  write ( *, '(a)' ) '  region: [-1,+1]'
  write ( *, '(a)' ) '  weight: w(x)=1'
  write ( *, '(a)' ) '  rules: Gauss-Legendre, Clenshaw-Curtis, Fejer2, Gauss-Patterson'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Jacobi:'
  write ( *, '(a)' ) '  region: [-1,+1]'
  write ( *, '(a)' ) '  weight: w(x)=(1-x)^alpha (1+x)^beta'
  write ( *, '(a)' ) '  rules: Gauss-Jacobi'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Laguerre:'
  write ( *, '(a)' ) '  region: [0,+oo)'
  write ( *, '(a)' ) '  weight: w(x)=exp(-x)'
  write ( *, '(a)' ) '  rules: Gauss-Laguerre'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Generalized Laguerre:'
  write ( *, '(a)' ) '  region: [0,+oo)'
  write ( *, '(a)' ) '  weight: w(x)=x^alpha exp(-x)'
  write ( *, '(a)' ) '  rules: Generalized Gauss-Laguerre'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Hermite:'
  write ( *, '(a)' ) '  region: (-oo,+o)'
  write ( *, '(a)' ) '  weight: w(x)=exp(-x*x)'
  write ( *, '(a)' ) '  rules: Gauss-Hermite'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Generalized Hermite:'
  write ( *, '(a)' ) '  region: (-oo,+oo)'
  write ( *, '(a)' ) '  weight: w(x)=|x|^alpha exp(-x*x)'
  write ( *, '(a)' ) '  rules: generalized Gauss-Hermite'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get the quadrature file root name:
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, quad_filename )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NINT_EXACTNESS_MIXED:'
    write ( *, '(a)' ) '  Enter the "root" name of the quadrature files.'

    read ( *, '(a)' ) quad_filename

  end if
!
!  Create the names of:
!    the quadrature A file;
!    the quadrature B file;
!    the quadrature R file;
!    the quadrature W file;
!    the quadrature X file,
!
  quad_a_filename = trim ( quad_filename ) // '_a.txt'
  quad_b_filename = trim ( quad_filename ) // '_b.txt'
  quad_r_filename = trim ( quad_filename ) // '_r.txt'
  quad_w_filename = trim ( quad_filename ) // '_w.txt'
  quad_x_filename = trim ( quad_filename ) // '_x.txt'
!
!  The second command line argument is the maximum degree.
!
  if ( 2 <= arg_num ) then

    iarg = 2
    call getarg ( iarg, string )
    call s_to_i4 ( string, degree_max, ierror, last )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NINT_EXACTNESS_MIXED:'
    write ( *, '(a)' ) '  Please enter the maximum total degree to check.'

    read ( *, * ) degree_max

  end if
!
!  Summarize the input.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NINT_EXACTNESS: User input:'
  write ( *, '(a)' ) '  Quadrature rule A file = "' // trim ( quad_a_filename ) // '".'
  write ( *, '(a)' ) '  Quadrature rule B file = "' // trim ( quad_b_filename ) // '".'
  write ( *, '(a)' ) '  Quadrature rule R file = "' // trim ( quad_r_filename ) // '".'
  write ( *, '(a)' ) '  Quadrature rule W file = "' // trim ( quad_w_filename ) // '".'
  write ( *, '(a)' ) '  Quadrature rule X file = "' // trim ( quad_x_filename ) // '".'
  write ( *, '(a,i8)' ) '  Maximum total degree to check = ', degree_max
!
!  Read the X file.
!
  call r8mat_header_read ( quad_x_filename, dim_num, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points  = ', point_num

  allocate ( x(dim_num,point_num) )

  call r8mat_data_read ( quad_x_filename, dim_num, point_num, x )
!
!  Read the W file.
!
  call r8mat_header_read ( quad_w_filename, dim_num2, point_num2 )

  if ( dim_num2 /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NINT_EXACTNESS_MIXED - Fatal error!'
    write ( *, '(a)' ) '  The quadrature weight file should have exactly'
    write ( *, '(a)' ) '  one value on each line.'
    stop
  end if

  if ( point_num2 /= point_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NINT_EXACTNESS_MIXED - Fatal error!'
    write ( *, '(a)' ) '  The quadrature weight file should have exactly'
    write ( *, '(a)' ) '  the same number of lines as the abscissa file.'
    stop
  end if

  allocate ( weight(point_num) )

  call r8mat_data_read ( quad_w_filename, 1, point_num, weight )
!
!  Read the R file.
!
  call r8mat_header_read ( quad_r_filename, dim_num2, point_num2 )

  if ( dim_num2 /= dim_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NINT_EXACTNESS_MIXED - Fatal error!'
    write ( *, '(a)' ) '  The quadrature region file should have the'
    write ( *, '(a)' ) '  same number of values on each line as the'
    write ( *, '(a)' ) '  abscissa file does.'
    stop
  end if

  if ( point_num2 /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NINT_EXACTNESS_MIXED - Fatal error!'
    write ( *, '(a)' ) '  The quadrature region file should have 2 lines.'
    stop
  end if

  allocate ( x_range(dim_num,2) )

  call r8mat_data_read ( quad_r_filename, dim_num, 2, x_range )
!
!  Read the A file.
!
  call r8mat_header_read ( quad_a_filename, dim_num2, point_num2 )

  if ( dim_num2 /= dim_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NINT_EXACTNESS_MIXED - Fatal error!'
    write ( *, '(a)' ) '  The quadrature A file should have the'
    write ( *, '(a)' ) '  same number of values on each line as the'
    write ( *, '(a)' ) '  abscissa file does.'
    stop
  end if

  if ( point_num2 /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NINT_EXACTNESS_MIXED - Fatal error!'
    write ( *, '(a)' ) '  The quadrature A file should have 1 line.'
    stop
  end if

  allocate ( alpha(dim_num) )

  call r8mat_data_read ( quad_a_filename, dim_num, 1, alpha )
!
!  Read the B file.
!
  call r8mat_header_read ( quad_b_filename, dim_num2, point_num2 )

  if ( dim_num2 /= dim_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NINT_EXACTNESS_MIXED - Fatal error!'
    write ( *, '(a)' ) '  The quadrature B file should have the'
    write ( *, '(a)' ) '  same number of values on each line as the'
    write ( *, '(a)' ) '  abscissa file does.'
    stop
  end if

  if ( point_num2 /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NINT_EXACTNESS_MIXED - Fatal error!'
    write ( *, '(a)' ) '  The quadrature B file should have 1 line.'
    stop
  end if

  allocate ( beta(dim_num) )

  call r8mat_data_read ( quad_b_filename, dim_num, 1, beta )
!
!  Try to determine the rule types.
!
  allocate ( rule(1:dim_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Analysis of integration region:'
  write ( *, '(a)' ) ' '

  do dim = 1, dim_num

    if ( x_range(dim,1) == -1.0D+00 .and. &
         x_range(dim,2) == +1.0D+00 ) then

      if ( alpha(dim) == 0.0 .and. beta(dim) == 0.0 ) then
        rule(dim) = 1
        write ( *, '(2x,i4,2x,a)' ) dim, '  Gauss Legendre'
      else
        rule(dim) = 2
        write ( *, '(2x,i4,2x,a,g14.6,a,g14.6)' ) &
          dim, '  Gauss Jacobi, ALPHA = ', alpha(dim), ' BETA = ', beta(dim)
      end if

    else if ( x_range(dim,1) == 0.0D+00 .and. &
              x_range(dim,2) == r8_huge ( ) ) then

      if ( alpha(dim) == 0.0 ) then
        rule(dim) = 3
        write ( *, '(2x,i4,2x,a)' ) dim, '  Gauss Laguerre.'
      else
        rule(dim) = 4
        write ( *, '(2x,i4,2x,a,g14.6)' ) &
          dim, '  Generalized Gauss Laguerre, ALPHA = ', alpha(dim)
      end if

    else if ( x_range(dim,1) == - r8_huge ( ) .and. &
              x_range(dim,2) == + r8_huge ( ) ) then

      if ( alpha(dim) == 0.0 ) then
        rule(dim) = 5
        write ( *, '(2x,i4,2x,a)' ) dim, '  Gauss Hermite.'
      else
        rule(dim) = 6
        write ( *, '(2x,i4,2x,a,g14.6)' ) &
          dim, '  Generalized Gauss Hermite, ALPHA = ', alpha(dim)
      end if

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NINT_EXACTNESS_MIXED - Fatal error!'
      write ( *, '(a)' ) '  Did not recognize region component.'
      write ( *, '(a,i8)' ) '  Dimension DIM = ', dim
      write ( *, '(a,g14.6)' ) '  A = ', x_range(dim,1)
      write ( *, '(a,g14.6)' ) '  B = ', x_range(dim,2)
      stop

    end if

  end do
!
!  Explore the monomials.
!
  allocate ( expon(1:dim_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          Error          Degree  Exponents'
  write ( *, '(a)' ) ' '

  do degree = 0, degree_max

    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( degree, dim_num, expon, more, h, t )

      call monomial_quadrature ( dim_num, point_num, rule, alpha, beta, &
        expon, weight, x, quad_error )

      write ( *, '(2x,f24.16,3x,i2,4x,10i3)' ) &
        quad_error, degree, expon(1:dim_num)

      if ( .not. more ) then
        exit
      end if

    end do

    write ( *, '(a)' ) ' '

  end do

  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( expon )
  deallocate ( rule )
  deallocate ( weight )
  deallocate ( x )
  deallocate ( x_range )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NINT_EXACTNESS_MIXED:'
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

  logical ch_eqi
  character c1
  character c1_cap
  character c2
  character c2_cap

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
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  
!    If C was 'illegal', then DIGIT is -1.
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
subroutine comp_next ( n, k, a, more, h, t )

!*****************************************************************************80
!
!! COMP_NEXT computes the compositions of the integer N into K parts.
!
!  Discussion:
!
!    A composition of the integer N into K parts is an ordered sequence
!    of K nonnegative integers which sum to N.  The compositions (1,2,1)
!    and (1,1,2) are considered to be distinct.
!
!    The routine computes one composition on each call until there are no more.
!    For instance, one composition of 6 into 3 parts is
!    3+2+1, another would be 6+0+0.
!
!    On the first call to this routine, set MORE = FALSE.  The routine
!    will compute the first element in the sequence of compositions, and
!    return it, as well as setting MORE = TRUE.  If more compositions
!    are desired, call again, and again.  Each time, the routine will
!    return with a new composition.
!
!    However, when the LAST composition in the sequence is computed 
!    and returned, the routine will reset MORE to FALSE, signaling that
!    the end of the sequence has been reached.
!
!    This routine originally used a SAVE statement to maintain the
!    variables H and T.  I have decided that it is safer
!    to pass these variables as arguments, even though the user should
!    never alter them.  This allows this routine to safely shuffle
!    between several ongoing calculations.
!
!
!    There are 28 compositions of 6 into three parts.  This routine will
!    produce those compositions in the following order:
!
!     I         A
!     -     ---------
!     1     6   0   0
!     2     5   1   0
!     3     4   2   0
!     4     3   3   0
!     5     2   4   0
!     6     1   5   0
!     7     0   6   0
!     8     5   0   1
!     9     4   1   1
!    10     3   2   1
!    11     2   3   1
!    12     1   4   1
!    13     0   5   1
!    14     4   0   2
!    15     3   1   2
!    16     2   2   2
!    17     1   3   2
!    18     0   4   2
!    19     3   0   3
!    20     2   1   3
!    21     1   2   3
!    22     0   3   3
!    23     2   0   4
!    24     1   1   4
!    25     0   2   4
!    26     1   0   5
!    27     0   1   5
!    28     0   0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 July 2008
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer whose compositions are desired.
!
!    Input, integer ( kind = 4 ) K, the number of parts in the composition.
!
!    Input/output, integer ( kind = 4 ) A(K), the parts of the composition.
!
!    Input/output, logical MORE, set by the user to start the computation,
!    and by the routine to terminate it.
!
!    Input/output, integer ( kind = 4 )  H, T, two internal parameters needed for the
!    computation.  The user should allocate space for these in the calling
!    program, include them in the calling sequence, but never alter them!
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) h
  logical more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t
!
!  The first computation.
!
  if ( .not. more ) then

    t = n
    h = 0
    a(1) = n
    a(2:k) = 0
!
!  The next computation.
!
  else

    if ( 1 < t ) then
      h = 0
    end if

    h = h + 1
    t = a(h)
    a(h) = 0
    a(1) = t - 1
    a(h+1) = a(h+1) + 1

  end if
!
!  This is the last element of the sequence if all the
!  items are in the last slot.
!
  more = ( a(k) /= n )

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
!    Output, integer ( kind = 4 ) COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer ( kind = 4 ) column_num
  logical got_one
  character ( len = * ) input_file_name
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  character ( len = 255 ) line
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
!    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
!
  implicit none

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) comment_num
  integer ( kind = 4 ) ierror
  character ( len = * ) input_file_name
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  character ( len = 255 ) line
  integer ( kind = 4 ) record_num
  integer ( kind = 4 ) row_num

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
function monomial_integral_generalized_hermite ( expon, alpha )

!*****************************************************************************80
!
!! MONOMIAL_INTEGRAL_GENERALIZED_HERMITE evaluates a 1D monomial generalized Hermite integral.
!
!  Discussion:
!
!    The integral being computed is
!
!      integral ( -oo < x < +oo ) x^n |x|^alpha exp(-x*x) dx 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent of the monomial.
!    0 <= EXPON.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of |x| in the weight.
!    -1.0 < ALPHA.
!
!    Output, real ( kind = 8 ) MONOMIAL_INTEGRAL_GENERALIZED_HERMITE, 
!    the value of the integral.
!
  implicit none

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) arg
  integer ( kind = 4 ) expon
  real    ( kind = 8 ) monomial_integral_generalized_hermite
  real    ( kind = 8 ) r8_gamma
  real    ( kind = 8 ) r8_huge
  real    ( kind = 8 ) value

  if ( mod ( expon, 2 ) == 1 ) then

    value = 0.0D+00

  else

    arg = alpha + real ( expon, kind = 8 )

    if ( arg <= - 1.0D+00 ) then
      value = - r8_huge ( )
    else
      arg = ( arg + 1.0D+00 ) / 2.0D+00
      value = r8_gamma ( arg )
    end if

  end if

  monomial_integral_generalized_hermite = value

  return
end
function monomial_integral_generalized_laguerre ( expon, alpha )

!*****************************************************************************80
!
!! MONOMIAL_INTEGRAL_GENERALIZED_LAGUERRE evaluates a 1D monomial generalized Laguerre integral.
!
!  Discussion:
!
!    The integral being computed is
!
!      integral ( 0 <= x < +oo ) x^n x^alpha exp(-x) dx 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent of the monomial.
!    0 <= EXPON.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of x in the weight.
!    -1.0 < ALPHA.
!
!    Output, real ( kind = 8 ) MONOMIAL_INTEGRAL_GENERALIZED_LAGUERRE, 
!    the value of the integral.
!
  implicit none

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) arg
  integer ( kind = 4 ) expon
  real    ( kind = 8 ) monomial_integral_generalized_laguerre
  real    ( kind = 8 ) r8_gamma
  real    ( kind = 8 ) value

  arg = alpha + real ( expon + 1, kind = 8 )

  value = r8_gamma ( arg )
  monomial_integral_generalized_laguerre = value

  return
end
function monomial_integral_hermite ( expon )

!*****************************************************************************80
!
!! MONOMIAL_INTEGRAL_HERMITE evaluates a 1D monomial Hermite integral.
!
!  Discussion:
!
!    H(n) = Integral ( -oo < x < +oo ) x^n exp(-x*x) dx
!
!    H(n) is 0 for n odd.
!
!    H(n) = (n-1)!! * sqrt(pi) / 2^(n/2) for n even.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.  
!    0 <= EXPON.
!
!    Output, real ( kind = 8 ) MONOMIAL_INTEGRAL_HERMITE, 
!    the value of the integral.
!
  implicit none

  integer ( kind = 4 ) expon
  real    ( kind = 8 ) monomial_integral_hermite
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) r8_factorial2
  real    ( kind = 8 ) r8_huge
  real    ( kind = 8 ) value

  if ( expon < 0 ) then

    value = - r8_huge ( )

  else if ( mod ( expon, 2 ) == 1 ) then

    value = 0.0D+00

  else

    value = r8_factorial2 ( expon - 1 ) * sqrt ( pi ) / 2.0D+00**( expon / 2 )

  end if

  monomial_integral_hermite = value

  return
end
function monomial_integral_jacobi ( expon, alpha, beta )

!*****************************************************************************80
!
!! MONOMIAL_INTEGRAL_JACOBI evaluates the integral of a monomial with Jacobi weight.
!
!  Discussion:
!
!    VALUE = Integral ( -1 <= X <= +1 ) x^EXPON (1-x)^ALPHA (1+x)^BETA dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of (1-X) in the weight factor.
!
!    Input, real ( kind = 8 ) BETA, the exponent of (1+X) in the weight factor.
!
!    Output, real ( kind = 8 ) MONOMIAL_INTEGRAL_JACOBI, the value of the integral.
!
  implicit none

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) arg1
  real    ( kind = 8 ) arg2
  real    ( kind = 8 ) arg3
  real    ( kind = 8 ) arg4
  real    ( kind = 8 ) beta
  real    ( kind = 8 ) c
  integer ( kind = 4 ) expon
  real    ( kind = 8 ) monomial_integral_jacobi
  real    ( kind = 8 ) r8_gamma
  real    ( kind = 8 ) s
  real    ( kind = 8 ) value
  real    ( kind = 8 ) value1
  real    ( kind = 8 ) value2

  c = real ( expon, kind = 8 )

  if ( mod ( expon, 2 ) == 0 ) then
    s = +1.0D+00
  else
    s = -1.0D+00
  end if

  arg1 = - alpha
  arg2 =   1.0D+00 + c
  arg3 =   2.0D+00 + beta + c
  arg4 = - 1.0D+00

  call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value1 )

  arg1 = - beta
  arg2 =   1.0D+00 + c
  arg3 =   2.0D+00 + alpha + c
  arg4 = - 1.0D+00

  call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value2 )

  value = r8_gamma ( 1.0D+00 + c ) * ( &
      s * r8_gamma ( 1.0D+00 + beta  ) * value1 &
    / r8_gamma ( 2.0D+00 + beta  + c ) &
    +     r8_gamma ( 1.0D+00 + alpha ) * value2 &
    / r8_gamma ( 2.0D+00 + alpha + c ) )

  monomial_integral_jacobi = value

  return
end
function monomial_integral_laguerre ( expon )

!*****************************************************************************80
!
!! MONOMIAL_INTEGRAL_LAGUERRE evaluates a 1D monomial Laguerre integral.
!
!  Discussion:
!
!    The integral being computed is
!
!      integral ( 0 <= x < +oo ) x^n * exp ( -x ) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!    0 <= EXPON.
!
!    Output, real ( kind = 8 ) MONOMIAL_INTEGRAL_LAGUERRE, 
!    the value of the integral.
!
  implicit none

  integer ( kind = 4 ) expon
  real    ( kind = 8 ) monomial_integral_laguerre
  real    ( kind = 8 ) r8_factorial

  monomial_integral_laguerre = r8_factorial ( expon )

  return
end
function monomial_integral_legendre ( expon )

!*****************************************************************************80
!
!! MONOMIAL_INTEGRAL_LEGENDRE evaluates a 1D monomial Legendre integral.
!
!  Discussion:
!
!    The integral being computed is
!
!      integral ( -1 <= x < +1 ) x^n dx 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXPON, the exponent.
!    0 <= EXPON.
!
!    Output, real ( kind = 8 ) MONOMIAL_INTEGRAL_LEGENDRE, 
!    the value of the integral.
!
  implicit none

  integer ( kind = 4 ) expon
  real    ( kind = 8 ) monomial_integral_legendre
  real    ( kind = 8 ) value

  if ( mod ( expon, 2 ) == 1 ) then
    value = 0.0D+00
  else if ( mod ( expon, 2 ) == 0 ) then
    value = 2.0D+00 / real ( expon + 1, kind = 8 )
  end if

  monomial_integral_legendre = value

  return
end
function monomial_integral_mixed ( dim_num, rule, alpha, beta, expon )

!*****************************************************************************80
!
!! MONOMIAL_INTEGRAL_MIXED evaluates a multi-D monomial mixed integral.
!
!  Discussion:
!
!    This routine evaluates the integral, over a multidimensional region.
!    of a monomial of the form
!
!      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
!
!    where the exponents are nonnegative integers.  Note that
!    if the combination 0^0 is encountered, it should be treated
!    as 1.
!
!    The integration is carried out in a region that is a direct product
!    of 1D factors that may be of various types,
!    and the integration includes the weight functions associated with
!    the 1D factors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) RULE(DIM_NUM), the component rules.
!    1, Gauss-Legendre rule on [-1,+1];
!    2, Gauss-Jacobi rule on [-1,+1];
!    3, Gauss-Laguerre rule on [0,+oo);
!    4, Generalized Gauss-Laguerre rule on [0,+oo);
!    5, Gauss-Hermite rule on (-oo,+oo);
!    6, Generalized Gauss-Hermite rule on (-oo,+oo).
!
!    Input, real ( kind = 8 ) ALPHA(DIM_NUM), BETA(DIM_NUM), parameters that
!    may be needed for Jacobi, Generalized-Laguerre, or Generalized Hermite rules.
!
!    Input, integer ( kind = 4 ) EXPON(DIM_NUM), the exponents.
!
!    Output, real ( kind = 8 ) MONOMIAL_INTEGRAL_MIXED, 
!    the value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real    ( kind = 8 ) alpha(dim_num)
  real    ( kind = 8 ) beta(dim_num)
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) expon(dim_num)
  real    ( kind = 8 ) monomial_integral_generalized_hermite
  real    ( kind = 8 ) monomial_integral_generalized_laguerre
  real    ( kind = 8 ) monomial_integral_hermite
  real    ( kind = 8 ) monomial_integral_jacobi
  real    ( kind = 8 ) monomial_integral_legendre
  real    ( kind = 8 ) monomial_integral_laguerre
  real    ( kind = 8 ) monomial_integral_mixed
  integer ( kind = 4 ) rule(dim_num)
  real    ( kind = 8 ) value

  value = 1.0D+00

  do dim = 1, dim_num

    if ( rule(dim) == 1 ) then
      value = value * monomial_integral_legendre ( expon(dim) )
    else if ( rule(dim) == 2 ) then
      value = value &
        * monomial_integral_jacobi ( expon(dim), alpha(dim), beta(dim) )
    else if ( rule(dim) == 3 ) then
      value = value * monomial_integral_laguerre ( expon(dim) )
    else if ( rule(dim) == 4 ) then
      value = value &
        * monomial_integral_generalized_laguerre ( expon(dim), alpha(dim) )
    else if ( rule(dim) == 5 ) then
      value = value * monomial_integral_hermite ( expon(dim) )
    else if ( rule(dim) == 6 ) then
      value = value &
        * monomial_integral_generalized_hermite ( expon(dim), alpha(dim) )
    end if

  end do

  monomial_integral_mixed = value

  return
end
subroutine monomial_quadrature ( dim_num, point_num, rule, alpha, beta, &
  expon, weight, x, quad_error )

!*****************************************************************************80
!
!! MONOMIAL_QUADRATURE applies a quadrature rule to a monomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points in the rule.
!
!    Input, integer ( kind = 4 ) RULE(DIM_NUM), the component rules.
!    1, Gauss-Legendre rule on [-1,+1];
!    2, Gauss-Jacobi rule on [-1,+1];
!    3, Gauss-Laguerre rule on [0,+oo);
!    4, Generalized Gauss-Laguerre rule on [0,+oo);
!    5, Gauss-Hermite rule on (-oo,+oo);
!    6, Generalized Gauss-Hermite rule on (-oo,+oo).
!
!    Input, real ( kind = 8 ) ALPHA(DIM_NUM), BETA(DIM_NUM), parameters that
!    may be needed for Jacobi, Generalized-Laguerre, or Generalized Hermite rules.
!
!    Input, integer ( kind = 4 ) EXPON(DIM_NUM), the exponents.
!
!    Input, real ( kind = 8 ) WEIGHT(POINT_NUM), the quadrature weights.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the quadrature points.
!
!    Output, real ( kind = 8 ) QUAD_ERROR, the quadrature error.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real    ( kind = 8 ) alpha(dim_num)
  real    ( kind = 8 ) beta(dim_num)
  real    ( kind = 8 ) exact
  integer ( kind = 4 ) expon(dim_num)
  real    ( kind = 8 ) monomial_integral_mixed
  integer ( kind = 4 ) point_num
  real    ( kind = 8 ) quad
  real    ( kind = 8 ) quad_error
  integer ( kind = 4 ) rule(dim_num)
  real    ( kind = 8 ) value(point_num)
  real    ( kind = 8 ) weight(point_num)
  real    ( kind = 8 ) x(dim_num,point_num)
!
!  Get the exact value of the integral of the unscaled monomial.
!
  exact = monomial_integral_mixed ( dim_num, rule, alpha, beta, expon )
!
!  Evaluate the monomial at the quadrature points.
!
  call monomial_value ( dim_num, point_num, expon, x, value )
!
!  Compute the weighted sum and divide by the exact value.
!
  quad = dot_product ( weight, value )
!
!  Error:
!
  if ( exact == 0.0D+00 ) then
    quad_error = abs ( quad - exact )
  else
    quad_error = abs ( quad - exact ) / abs ( exact )
  end if

  return
end
subroutine monomial_value ( dim_num, point_num, expon, x, value )

!*****************************************************************************80
!
!! MONOMIAL_VALUE evaluates a monomial.
!
!  Discussion:
!
!    This routine evaluates a multidimensional monomial which is a product 
!    of 1D factors of the form x(dim)^expon(dim).
!
!    The exponents are nonnegative integers.  
!
!    Note that if the combination 0^0 is encountered, it should be treated
!    as 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points at which the
!    monomial is to be evaluated.
!
!    Input, integer ( kind = 4 ) EXPON(DIM_NUM), the exponents.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the point coordinates.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the value of the monomial.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) expon(dim_num)
  real    ( kind = 8 ) value(point_num)
  real    ( kind = 8 ) x(dim_num,point_num)

  value(1:point_num) = 1.0D+00

  do dim = 1, dim_num

    if ( 0 /= expon(dim) ) then
      value(1:point_num) = value(1:point_num) * x(dim,1:point_num)**expon(dim)
    end if

  end do

  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial of N.
!
!  Discussion:
!
!    factorial ( N ) = product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!    Output, real ( kind = 8 ) R8_FACTORIAL, the factorial of N.
!
  implicit none

  real    ( kind = 8 ) r8_factorial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  r8_factorial = 1.0D+00

  do i = 1, n
    r8_factorial = r8_factorial * real ( i, kind = 8 )
  end do

  return
end
function r8_factorial2 ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL2 computes the double factorial function N!!
!
!  Discussion:
!
!    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Example:
!
!     N    N!!
!
!     0     1
!     1     1
!     2     2
!     3     3
!     4     8
!     5    15
!     6    48
!     7   105
!     8   384
!     9   945
!    10  3840
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 September 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the double factorial 
!    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
!
!    Output, real ( kind = 8 ) R8_FACTORIAL2, the value of N!!.
!
  implicit none

  integer ( kind = 4 ) n
  real    ( kind = 8 ) r8_factorial2
  real    ( kind = 8 ) r8_n

  if ( n < 1 ) then
    r8_factorial2 = 1.0D+00
    return
  end if

  r8_n = real ( n, kind = 8 )
  r8_factorial2 = 1.0D+00

  do while ( 1.0D+00 < r8_n )
    r8_factorial2 = r8_factorial2 * r8_n
    r8_n = r8_n - 2.0D+00
  end do

  return
end
function r8_gamma ( x )

!*****************************************************************************80
!
!! R8_GAMMA evaluates Gamma(X) for a real argument.
!
!  Discussion:
!
!    This routine calculates the gamma function for a real argument X.
!
!    Computation is based on an algorithm outlined in reference 1.
!    The program uses rational functions that approximate the gamma
!    function to at least 20 significant decimal digits.  Coefficients
!    for the approximation over the interval (1,2) are unpublished.
!    Those for the approximation for 12 <= X are from reference 2.
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody,
!    An Overview of Software Development for Special Functions,
!    in Numerical Analysis Dundee, 1975,
!    edited by GA Watson,
!    Lecture Notes in Mathematics 506,
!    Springer, 1976.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_GAMMA, the value of the function.
!
  implicit none
!
!  Coefficients for minimax approximation over (12, INF).
!
  real    ( kind = 8 ), dimension ( 7 ) :: c = (/ &
   -1.910444077728D-03, &
    8.4171387781295D-04, &
   -5.952379913043012D-04, &
    7.93650793500350248D-04, &
   -2.777777777777681622553D-03, &
    8.333333333333333331554247D-02, &
    5.7083835261D-03 /)
  real    ( kind = 8 ), parameter :: eps = 2.22D-16
  real    ( kind = 8 ) fact
  real    ( kind = 8 ), parameter :: half = 0.5D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real    ( kind = 8 ), parameter :: one = 1.0D+00
  real    ( kind = 8 ), dimension ( 8 ) :: p = (/ &
    -1.71618513886549492533811D+00, &
     2.47656508055759199108314D+01, &
    -3.79804256470945635097577D+02, &
     6.29331155312818442661052D+02, &
     8.66966202790413211295064D+02, &
    -3.14512729688483675254357D+04, &
    -3.61444134186911729807069D+04, &
     6.64561438202405440627855D+04 /)
  logical parity
  real    ( kind = 8 ), parameter :: pi = 3.1415926535897932384626434D+00
  real    ( kind = 8 ), dimension ( 8 ) :: q = (/ &
    -3.08402300119738975254353D+01, &
     3.15350626979604161529144D+02, &
    -1.01515636749021914166146D+03, &
    -3.10777167157231109440444D+03, &
     2.25381184209801510330112D+04, &
     4.75584627752788110767815D+03, &
    -1.34659959864969306392456D+05, &
    -1.15132259675553483497211D+05 /)
  real    ( kind = 8 ) r8_gamma
  real    ( kind = 8 ) res
  real    ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real    ( kind = 8 ) sum
  real    ( kind = 8 ), parameter :: twelve = 12.0D+00
  real    ( kind = 8 ), parameter :: two = 2.0D+00
  real    ( kind = 8 ) x
  real    ( kind = 8 ), parameter :: xbig = 171.624D+00
  real    ( kind = 8 ) xden
  real    ( kind = 8 ), parameter :: xinf = 1.0D+30
  real    ( kind = 8 ), parameter :: xminin = 2.23D-308
  real    ( kind = 8 ) xnum
  real    ( kind = 8 ) y
  real    ( kind = 8 ) y1
  real    ( kind = 8 ) ysq
  real    ( kind = 8 ) z
  real    ( kind = 8 ), parameter :: zero = 0.0D+00

  parity = .false.
  fact = one
  n = 0
  y = x
!
!  Argument is negative.
!
  if ( y <= zero ) then

    y = - x
    y1 = aint ( y )
    res = y - y1

    if ( res /= zero ) then

      if ( y1 /= aint ( y1 * half ) * two ) then
        parity = .true.
      end if

      fact = - pi / sin ( pi * res )
      y = y + one

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Argument is positive.
!
  if ( y < eps ) then
!
!  Argument < EPS.
!
    if ( xminin <= y ) then
      res = one / y
    else
      res = xinf
      r8_gamma = res
      return
    end if

  else if ( y < twelve ) then

    y1 = y
!
!  0.0 < argument < 1.0.
!
    if ( y < one ) then

      z = y
      y = y + one
!
!  1.0 < argument < 12.0.
!  Reduce argument if necessary.
!
    else

      n = int ( y ) - 1
      y = y - real ( n, kind = 8 )
      z = y - one

    end if
!
!  Evaluate approximation for 1.0 < argument < 2.0.
!
    xnum = zero
    xden = one
    do i = 1, 8
      xnum = ( xnum + p(i) ) * z
      xden = xden * z + q(i)
    end do

    res = xnum / xden + one
!
!  Adjust result for case  0.0 < argument < 1.0.
!
    if ( y1 < y ) then

      res = res / y1
!
!  Adjust result for case 2.0 < argument < 12.0.
!
    else if ( y < y1 ) then

      do i = 1, n
        res = res * y
        y = y + one
      end do

    end if

  else
!
!  Evaluate for 12.0 <= argument.
!
    if ( y <= xbig ) then

      ysq = y * y
      sum = c(7)
      do i = 1, 6
        sum = sum / ysq + c(i)
      end do
      sum = sum / y - y + sqrtpi
      sum = sum + ( y - half ) * log ( y )
      res = exp ( sum )

    else

      res = xinf
      r8_gamma = res
      return

    end if

  end if
!
!  Final adjustments and return.
!
  if ( parity ) then
    res = - res
  end if

  if ( fact /= one ) then
    res = fact / res
  end if

  r8_gamma = res

  return
end
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a very large R8.
!
!  Discussion:
!
!    The value returned by this function is NOT required to be the
!    maximum representable R8.  This value varies from machine to machine,
!    from compiler to compiler, and may cause problems when being printed.
!    We simply want a "very large" but non-infinite number.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_HUGE, a "huge" value.
!
  implicit none

  real    ( kind = 8 ) r8_huge

  r8_huge = 1.0D+30

  return
end
subroutine r8_hyper_2f1 ( a, b, c, x, hf )

!*****************************************************************************80
!
!! R8_HYPER_2F1 evaluates the hypergeometric function 2F1(A,B,C,X).
!
!  Discussion:
!
!    A minor bug was corrected.  The HW variable, used in several places as
!    the "old" value of a quantity being iteratively improved, was not
!    being initialized.  JVB, 11 February 2008.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original F77 version by Shanjie Zhang, Jianming Jin.
!    F90 version by John Burkardt.
!
!    The F77 original version of this routine is copyrighted by 
!    Shanjie Zhang and Jianming Jin.  However, they give permission to 
!    incorporate this routine into user programs provided that the copyright 
!    is acknowledged.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, X, the arguments of the function.
!    C must not be equal to a nonpositive integer.
!    X < 1.
!
!    Output, real ( kind = 8 ) HF, the value of the function.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) a0
  real    ( kind = 8 ) aa
  real    ( kind = 8 ) b
  real    ( kind = 8 ) bb
  real    ( kind = 8 ) c
  real    ( kind = 8 ) c0
  real    ( kind = 8 ) c1
  real    ( kind = 8 ), parameter :: el = 0.5772156649015329D+00
  real    ( kind = 8 ) eps
  real    ( kind = 8 ) f0
  real    ( kind = 8 ) f1
  real    ( kind = 8 ) g0
  real    ( kind = 8 ) g1
  real    ( kind = 8 ) g2
  real    ( kind = 8 ) g3
  real    ( kind = 8 ) ga
  real    ( kind = 8 ) gabc
  real    ( kind = 8 ) gam
  real    ( kind = 8 ) gb
  real    ( kind = 8 ) gbm
  real    ( kind = 8 ) gc
  real    ( kind = 8 ) gca
  real    ( kind = 8 ) gcab
  real    ( kind = 8 ) gcb
  real    ( kind = 8 ) gm
  real    ( kind = 8 ) hf
  real    ( kind = 8 ) hw
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical l0
  logical l1
  logical l2
  logical l3
  logical l4
  logical l5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nm
  real    ( kind = 8 ) pa
  real    ( kind = 8 ) pb
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) r
  real    ( kind = 8 ) r0
  real    ( kind = 8 ) r1
  real    ( kind = 8 ) r8_gamma
  real    ( kind = 8 ) r8_psi
  real    ( kind = 8 ) rm
  real    ( kind = 8 ) rp
  real    ( kind = 8 ) sm
  real    ( kind = 8 ) sp
  real    ( kind = 8 ) sp0
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xx

  l0 = ( c == aint ( c ) ) .and. ( c < 0.0D+00 )
  l1 = ( 1.0D+00 - x < 1.0D-15 ) .and. ( c - a - b <= 0.0D+00 )
  l2 = ( a == aint ( a ) ) .and. ( a < 0.0D+00 )
  l3 = ( b == aint ( b ) ) .and. ( b < 0.0D+00 )
  l4 = ( c - a == aint ( c - a ) ) .and. ( c - a <= 0.0D+00 )
  l5 = ( c - b == aint ( c - b ) ) .and. ( c - b <= 0.0D+00 )

  if ( l0 .or. l1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_HYPER_2F1 - Fatal error!'
    write ( *, '(a)' ) '  The hypergeometric series is divergent.'
    return
  end if

  if ( 0.95D+00 < x ) then 
    eps = 1.0D-08
  else
    eps = 1.0D-15
  end if

  if ( x == 0.0D+00 .or. a == 0.0D+00 .or. b == 0.0D+00 ) then

    hf = 1.0D+00
    return

  else if ( 1.0D+00 - x == eps .and. 0.0D+00 < c - a - b ) then

    gc = r8_gamma ( c )
    gcab = r8_gamma ( c - a - b )
    gca = r8_gamma ( c - a )
    gcb = r8_gamma ( c - b )
    hf = gc * gcab / ( gca * gcb )
    return

  else if ( 1.0D+00 + x <= eps .and. abs ( c - a + b - 1.0D+00 ) <= eps ) then

    g0 = sqrt ( pi ) * 2.0D+00**( - a )
    g1 = r8_gamma ( c )
    g2 = r8_gamma ( 1.0D+00 + a / 2.0D+00 - b )
    g3 = r8_gamma ( 0.5D+00 + 0.5D+00 * a )
    hf = g0 * g1 / ( g2 * g3 )
    return

  else if ( l2 .or. l3 ) then

    if ( l2 ) then
      nm = int ( abs ( a ) )
    end if

    if ( l3 ) then
      nm = int ( abs ( b ) )
    end if

    hf = 1.0D+00
    r = 1.0D+00

    do k = 1, nm
      r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x
      hf = hf + r
    end do

    return

  else if ( l4 .or. l5 ) then

    if ( l4 ) then
      nm = int ( abs ( c - a ) )
    end if

    if ( l5 ) then
      nm = int ( abs ( c - b ) )
    end if

    hf = 1.0D+00
    r  = 1.0D+00
    do k = 1, nm
      r = r * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x
      hf = hf + r
    end do
    hf = ( 1.0D+00 - x )**( c - a - b ) * hf
    return

  end if

  aa = a
  bb = b
  xx = x
!
!  WARNING: ALTERATION OF INPUT ARGUMENTS A, B, and X, WHICH MIGHT BE CONSTANTS.
!
  if ( x < 0.0D+00 ) then
    x = x / ( x - 1.0D+00 )
    if ( a < c .and. b < a .and. 0.0D+00 < b ) then
      a = bb
      b = aa
    end if
    b = c - b
  end if

  if ( 0.75D+00 <= x ) then

    gm = 0.0D+00

    if ( abs ( c - a - b - aint ( c - a - b ) ) < 1.0D-15 ) then

      m = int ( c - a - b )
      ga = r8_gamma ( a )
      gb = r8_gamma ( b )
      gc = r8_gamma ( c )
      gam = r8_gamma ( a + m )
      gbm = r8_gamma ( b + m )

      pa = r8_psi ( a )
      pb = r8_psi ( b )

      if ( m /= 0 ) then
        gm = 1.0D+00
      end if

      do j = 1, abs ( m ) - 1
        gm = gm * j
      end do

      rm = 1.0D+00
      do j = 1, abs ( m )
        rm = rm * j
      end do

      f0 = 1.0D+00
      r0 = 1.0D+00
      r1 = 1.0D+00
      sp0 = 0.0D+00
      sp = 0.0D+00

      if ( 0 <= m ) then

        c0 = gm * gc / ( gam * gbm )
        c1 = - gc * ( x - 1.0D+00 )**m / ( ga * gb * rm )

        do k = 1, m - 1
          r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
            / ( k * ( k - m ) ) * ( 1.0D+00 - x )
          f0 = f0 + r0
        end do

        do k = 1, m
          sp0 = sp0 + 1.0D+00 / ( a + k - 1.0D+00 ) &
            + 1.0D+00 / ( b + k - 1.0D+00 ) - 1.0D+00 / real ( k, kind = 8 )
        end do

        f1 = pa + pb + sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
        hw = f1

        do k = 1, 250

          sp = sp + ( 1.0D+00 - a ) / ( k * ( a + k - 1.0D+00 ) ) &
            + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

          sm = 0.0D+00
          do j = 1, m
            sm = sm + ( 1.0D+00 - a ) &
              / ( ( j + k ) * ( a + j + k - 1.0D+00 ) ) &
              + 1.0D+00 / ( b + j + k - 1.0D+00 )
          end do

          rp = pa + pb + 2.0D+00 * el + sp + sm + log ( 1.0D+00 - x )

          r1 = r1 * ( a + m + k - 1.0D+00 ) * ( b + m + k - 1.0D+00 ) &
            / ( k * ( m + k ) ) * ( 1.0D+00 - x )

          f1 = f1 + r1 * rp

          if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
            exit
          end if

          hw = f1

        end do

        hf = f0 * c0 + f1 * c1

      else if ( m < 0 ) then

        m = - m
        c0 = gm * gc / ( ga * gb * ( 1.0D+00 - x )**m )
        c1 = - ( - 1 )**m * gc / ( gam * gbm * rm )

        do k = 1, m - 1
          r0 = r0 * ( a - m + k - 1.0D+00 ) * ( b - m + k - 1.0D+00 ) &
            / ( k * ( k - m ) ) * ( 1.0D+00 - x )
          f0 = f0 + r0
        end do

        do k = 1, m
          sp0 = sp0 + 1.0D+00 / real ( k, kind = 8 )
        end do

        f1 = pa + pb - sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
        hw = f1

        do k = 1, 250

          sp = sp + ( 1.0D+00 - a ) &
            / ( k * ( a + k - 1.0D+00 ) ) &
            + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

          sm = 0.0D+00
          do j = 1, m
            sm = sm + 1.0D+00 / real ( j + k, kind = 8 )
          end do

          rp = pa + pb + 2.0D+00 * el + sp - sm + log ( 1.0D+00 - x )

          r1 = r1 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
            / ( k * ( m + k ) ) * ( 1.0D+00 - x )

          f1 = f1 + r1 * rp

          if ( abs ( f1 - hw ) < abs ( f1 ) * eps ) then
            exit
          end if

          hw = f1

        end do

        hf = f0 * c0 + f1 * c1

      end if

    else

      ga = r8_gamma ( a )
      gb = r8_gamma ( b )
      gc = r8_gamma ( c )
      gca = r8_gamma ( c - a )
      gcb = r8_gamma ( c - b )
      gcab = r8_gamma ( c - a - b )
      gabc = r8_gamma ( a + b - c )
      c0 = gc * gcab / ( gca * gcb )
      c1 = gc * gabc / ( ga * gb ) * ( 1.0D+00 - x )**( c - a - b )
      hf = 0.0D+00
      hw = hf
      r0 = c0
      r1 = c1

      do k = 1, 250

        r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
          / ( k * ( a + b - c + k ) ) * ( 1.0D+00 - x )

        r1 = r1 * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) &
          / ( k * ( c - a - b + k ) ) * ( 1.0D+00 - x )

        hf = hf + r0 + r1

        if ( abs ( hf - hw ) < abs ( hf ) * eps ) then
          exit
        end if

        hw = hf

      end do

      hf = hf + c0 + c1

    end if

  else

    a0 = 1.0D+00

    if ( a < c .and. c < 2.0D+00 * a .and. b < c .and. c < 2.0D+00 * b ) then

      a0 = ( 1.0D+00 - x )**( c - a - b )
      a = c - a
      b = c - b

    end if

    hf = 1.0D+00
    hw = hf
    r = 1.0D+00

    do k = 1, 250

      r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) &
        / ( k * ( c + k - 1.0D+00 ) ) * x

      hf = hf + r

      if ( abs ( hf - hw ) <= abs ( hf ) * eps ) then
        exit
      end if

      hw = hf

    end do

    hf = a0 * hf

  end if

  if ( xx < 0.0D+00 ) then
    x = xx
    c0 = 1.0D+00 / ( 1.0D+00 - x )**aa
    hf = c0 * hf
  end if

  a = aa
  b = bb

  if ( 120 < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_HYPER_2F1 - Warning!'
    write ( *, '(a)' ) '  A large number of iterations were needed.'
    write ( *, '(a)' ) '  The accuracy of the results should be checked.'
  end if

  return
end
function r8_psi ( xx )

!*****************************************************************************80
!
!! R8_PSI evaluates the function Psi(X).
!
!  Discussion:
!
!    This routine evaluates the logarithmic derivative of the
!    Gamma function,
!
!      PSI(X) = d/dX ( GAMMA(X) ) / GAMMA(X)
!             = d/dX LN ( GAMMA(X) )
!
!    for real X, where either
!
!      - XMAX1 < X < - XMIN, and X is not a negative integer,
!
!    or
!
!      XMIN < X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by William Cody
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    William Cody, Anthony Strecok, Henry Thacher,
!    Chebyshev Approximations for the Psi Function,
!    Mathematics of Computation,
!    Volume 27, Number 121, January 1973, pages 123-127.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XX, the argument of the function.
!
!    Output, real ( kind = 8 ) R8_PSI, the value of the function.
!
  implicit none

  real    ( kind = 8 ) aug
  real    ( kind = 8 ) den
  real    ( kind = 8 ), parameter :: four = 4.0D+00
  real    ( kind = 8 ), parameter :: fourth = 0.25D+00
  real    ( kind = 8 ), parameter :: half = 0.5D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nq
  real    ( kind = 8 ), parameter :: one = 1.0D+00
  real    ( kind = 8 ), dimension ( 9 ) :: p1 = (/ &
   4.5104681245762934160D-03, &
   5.4932855833000385356D+00, &
   3.7646693175929276856D+02, &
   7.9525490849151998065D+03, &
   7.1451595818951933210D+04, &
   3.0655976301987365674D+05, &
   6.3606997788964458797D+05, &
   5.8041312783537569993D+05, &
   1.6585695029761022321D+05 /)
  real    ( kind = 8 ), dimension ( 7 ) :: p2 = (/ &
  -2.7103228277757834192D+00, &
  -1.5166271776896121383D+01, &
  -1.9784554148719218667D+01, &
  -8.8100958828312219821D+00, &
  -1.4479614616899842986D+00, &
  -7.3689600332394549911D-02, &
  -6.5135387732718171306D-21 /)
  real    ( kind = 8 ), parameter :: piov4 = 0.78539816339744830962D+00
  real    ( kind = 8 ), dimension ( 8 ) :: q1 = (/ &
   9.6141654774222358525D+01, &
   2.6287715790581193330D+03, &
   2.9862497022250277920D+04, &
   1.6206566091533671639D+05, &
   4.3487880712768329037D+05, &
   5.4256384537269993733D+05, &
   2.4242185002017985252D+05, &
   6.4155223783576225996D-08 /)
  real    ( kind = 8 ), dimension ( 6 ) :: q2 = (/ &
   4.4992760373789365846D+01, &
   2.0240955312679931159D+02, &
   2.4736979003315290057D+02, &
   1.0742543875702278326D+02, &
   1.7463965060678569906D+01, &
   8.8427520398873480342D-01 /)
  real    ( kind = 8 ) r8_psi
  real    ( kind = 8 ) sgn
  real    ( kind = 8 ), parameter :: three = 3.0D+00
  real    ( kind = 8 ) upper
  real    ( kind = 8 ) w
  real    ( kind = 8 ) x
  real    ( kind = 8 ), parameter :: x01 = 187.0D+00
  real    ( kind = 8 ), parameter :: x01d = 128.0D+00
  real    ( kind = 8 ), parameter :: x02 = 6.9464496836234126266D-04
  real    ( kind = 8 ), parameter :: xinf = 1.70D+38
  real    ( kind = 8 ), parameter :: xlarge = 2.04D+15
  real    ( kind = 8 ), parameter :: xmax1 = 3.60D+16
  real    ( kind = 8 ), parameter :: xmin1 = 5.89D-39
  real    ( kind = 8 ), parameter :: xsmall = 2.05D-09
  real    ( kind = 8 ) xx
  real    ( kind = 8 ) z
  real    ( kind = 8 ), parameter :: zero = 0.0D+00

  x = xx
  w = abs ( x )
  aug = zero
!
!  Check for valid arguments, then branch to appropriate algorithm.
!
  if ( xmax1 <= - x .or. w < xmin1 ) then

    if ( zero < x ) then
      r8_psi = - xinf
    else
      r8_psi = xinf
    end if

    return
  end if

  if ( x < half ) then
!
!  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
!  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
!
    if ( w <= xsmall ) then

      aug = - one / x
!
!  Argument reduction for cotangent.
!
    else

      if ( x < zero ) then
        sgn = piov4
      else
        sgn = - piov4
      end if

      w = w - real ( int ( w ), kind = 8 )
      nq = int ( w * four )
      w = four * ( w - real ( nq, kind = 8 ) * fourth )
!
!  W is now related to the fractional part of 4.0 * X.
!  Adjust argument to correspond to values in the first
!  quadrant and determine the sign.
!
      n = nq / 2

      if ( n + n /= nq ) then
        w = one - w
      end if

      z = piov4 * w

      if ( mod ( n, 2 ) /= 0 ) then
        sgn = - sgn
      end if
!
!  Determine the final value for  -pi * cotan(pi*x).
!
      n = ( nq + 1 ) / 2
      if ( mod ( n, 2 ) == 0 ) then
!
!  Check for singularity.
!
        if ( z == zero ) then

          if ( zero < x ) then
            r8_psi = -xinf
          else
            r8_psi = xinf
          end if

          return
        end if

        aug = sgn * ( four / tan ( z ) )

      else

        aug = sgn * ( four * tan ( z ) )

      end if

    end if

    x = one - x

  end if
!
!  0.5 <= X <= 3.0.
!
  if ( x <= three ) then

    den = x
    upper = p1(1) * x
    do i = 1, 7
      den = ( den + q1(i) ) * x
      upper = ( upper + p1(i+1) ) * x
    end do
    den = ( upper + p1(9) ) / ( den + q1(8) )
    x = ( x - x01 / x01d ) - x02
    r8_psi = den * x + aug
    return

  end if
!
!  3.0 < X.
!
  if ( x < xlarge ) then
    w = one / ( x * x )
    den = w
    upper = p2(1) * w
    do i = 1, 5
      den = ( den + q2(i) ) * w
      upper = ( upper + p2(i+1) ) * w
    end do
    aug = ( upper + p2(7) ) / ( den + q2(6) ) - half / x + aug
  end if

  r8_psi = aug + log ( x )

  return
end
subroutine r8mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_DATA_READ reads data from an R8MAT file.
!
!  Discussion:
!
!    The file may contain more than N points, but this routine will
!    return after reading N of them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2008
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
!    Output, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 )   m
  integer   ( kind = 4 )   n

  integer   ( kind = 4 )   ierror
  character ( len = * )    input_filename
  integer   ( kind = 4 )   input_status
  integer   ( kind = 4 )   input_unit
  integer   ( kind = 4 )   j
  character ( len = 255 )  line
  real      ( kind = 8 )   table(m,n)
  real      ( kind = 8 )   x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i8)' ) '  Number of values expected per line M = ', m
      write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

  end do

  close ( unit = input_unit )

  return
end
subroutine r8mat_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! R8MAT_HEADER_READ reads the header from an R8MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
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

  character ( len = * )  input_filename
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  call file_column_count ( input_filename, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
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
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
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
!    07 September 2004
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

  logical ch_eqi
  character c
  real    ( kind = 8 ) dval
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
  real    ( kind = 8 ) rbot
  real    ( kind = 8 ) rexp
  real    ( kind = 8 ) rtop
  character ( len = * ) s

  nchar = len_trim ( s )

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

    if ( nchar < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
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
!    07 September 2004
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
  integer ( kind = 4 ) lchar
  real    ( kind = 8 ) rvec(n)
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + lchar

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
