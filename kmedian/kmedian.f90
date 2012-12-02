program main

!*****************************************************************************80
!
!! MAIN is the main program for KMEDIAN.
!
!  Discussion:
!
!    KMEDIAN solves the weighted K median problem using relaxation.
!
!  Modified:
!
!    10 November 2002
!
!  Reference:
!
!    Gerard Cornuejols, Marshall Fisher, George Nemhauser,
!    Location of Bank Accounts to Optimize Float,
!    an Analytic Study of Exact and Approximate Algorithms,
!    Management Science,
!    Volume 23, Number 8, 1977, pages 789-810.
!
  implicit none

  real ambdk
  integer i
  real, allocatable, dimension(:,:) ::  dist
  character ( len = 80 ) file_dist
  character ( len = 80 ) file_main
  character ( len = 80 ) file_weight
  integer isw
  integer, allocatable, dimension(:) ::  itemp
  integer iter
  integer itera
  integer j
  integer, allocatable, dimension(:,:) ::  jsrt
  integer, allocatable, dimension(:) ::  jwarhs
  integer kopen
  character ( len = 80 ) :: lpt_file = 'kmedian.lpt'
  integer lpt_unit
  integer no
  integer, allocatable, dimension(:) ::  nv
  integer nvsq
  real rward
  real, allocatable, dimension(:,:) ::  shodi
  real, allocatable, dimension(:) ::  temp
  real tk
  real, allocatable, dimension(:) ::  u
  real w
  real w_best
  real, allocatable, dimension(:) :: weight

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'KMEDIAN'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Lagrangian relaxation method for'
  write ( *, '(a)' ) '  the K-median problem with weights.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Algorithm of Cornuejols, Fisher and Nemhauser.'

  call get_unit ( lpt_unit )

  open ( unit = lpt_unit, file = lpt_file, status = 'replace' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter the filename of your input.'

! read ( *, '(a)' ) file_main
  file_main = 'co04_main.txt'

  write ( *, '(a)' ) '  Reading main input file "' // trim ( file_main ) // '".'

  call main_read_size ( file_main, no )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of cities is ', no

  allocate ( dist(1:no,1:no) )
  allocate ( itemp(1:no) )
  allocate ( jsrt(1:no,1:no) )
  allocate ( jwarhs(1:no) )
  allocate ( nv(1:no) )
  allocate ( shodi(1:no,1:no) )
  allocate ( temp(1:no) )
  allocate ( u(1:no) )
  allocate ( weight(1:no) )

  call main_read_dist ( file_main, file_dist )

  write ( *, '(a)' ) '  Reading distance input file "' // trim ( file_dist ) // '".'

  call dist_read ( file_dist, no, dist )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Distance matrix:'
  write ( *, '(a)' ) ' '

  do i = 1, no
    write ( *, '(4f10.4)' ) dist(i,1:no)
  end do

  call main_read_weight ( file_main, file_weight )

  write ( *, '(a)' ) '  Reading weight input file "' // trim ( file_weight ) // '".'

  call weight_read ( file_weight, no, weight )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Weights'
  write ( *, '(a)' ) ' '

  do i = 1, no
    write ( *, '(f10.4)' ) weight(i)
  end do

  do i = 1, no
    shodi(i,1:no) = dist(i,1:no) * weight(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Weighted distance matrix:'
  write ( *, '(a)' ) ' '

  do i = 1, no
    write ( *, '(4f10.4)' ) shodi(i,1:no)
  end do

  do

    write ( *, '(a)' )
    write ( *, '(a)' ) '  Enter the number of facilities:'
    write ( *, '(a,i6)' ) '  This must be between 1 and ', no

    kopen = 4
!   read ( *, * ) kopen

    if ( 1 <= kopen .and. kopen <= no ) then
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEDIAN - Error!'
    write ( *, '(a)' ) '  That value was not acceptable!'

  end do
!
!  Get a greedy (nonoptimal) solution to the problem.
!
  call greedy ( kopen, no, shodi, jwarhs, rward )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Selected cities:'
  write ( *, '(a)' ) ' '

  do i = 1, kopen
    write ( *, * ) i, jwarhs(i)
  end do

  write ( *, * ) '  RWARD = ', rward

  stop
!
!  For each row I of SHODI,
!  * determine the rank of the elements in ascending order;
!  * store these ranks in JSRT(I,*);
!  * store the negative of the smallest element in U(I).
!
!  So I guess we expect all SHODI's to be nonnegative?
!
  do i = 1, no

    temp(1:no) = shodi(i,1:no)

    do j = 1, no
      itemp(j) = j
    end do

    call bubble_sort_a ( no, temp, itemp )

    jsrt(i,1:no) = itemp(1:no)

    u(i) = -shodi(i,itemp(1))

  end do
!
!  Solution of fixed Lagrangian multiplier.
!
  w_best = -1000000.0E+00
  iter = 0
  ambdk = 2.0E+00

  isw = 0
  itera = 0

  do

    call sol ( shodi, isw, itera, jsrt, jwarhs, kopen, lpt_unit, &
      no, rward, nv, nvsq, u, w, w_best )

    if ( isw /= 0 ) then
      exit
    end if

    if ( w_best <= w ) then

      w_best = w
      iter = 0

    else

      if ( 8 < iter ) then
        ambdk = ambdk / 2.0E+00
        iter = 0
      end if

    end if
!
!  Computation of a new multiplier vector for subgradient iteration.
!
    tk = ( ambdk * ( rward - w ) ) / nvsq
    iter = iter + 1

    u(1:no) = u(1:no) + tk * nv(1:no)

    write ( lpt_unit, '(a)' ) ' '
    write ( lpt_unit, '(a,g14.6)' ) '  Best heuristic solution ', rward
    write ( lpt_unit, '(a,g14.6)' ) '  Best Lagrangian bound   ', w_best

  end do

  close ( unit = lpt_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Your output is in file ' // trim ( lpt_file )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'KMEDIAN'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine bubble_sort_a ( n, a, ind )

!*****************************************************************************80
!
!! BUBBLE_SORT_A ascending sorts a real vector using bubble sort.
!
!  Discussion:
!
!    Bubble sort is simple to program, but inefficient.  It should not
!    be used for large arrays.
!
!    An associated index vector is adjusted to correspond to the
!    changes in the real vector.
!
!  Modified:
!
!    09 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, real A(N).
!    On input, an unsorted array.
!    On output, the array has been sorted.
!
!    Input/output, real IND(N).
!    On input, an integer array.
!    On output, the entries of this array have been shifted
!    in a way that corresponds to the changes to the real vector.
!
  implicit none

  integer n

  real a(n)
  integer i
  integer ind(n)
  integer j

  do i = 1, n-1
    do j = i+1, n
      if ( a(i) > a(j) ) then
        call r4_swap ( a(i), a(j) )
        call i4_swap ( ind(i), ind(j) )
      end if
    end do
  end do

  return
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
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
  integer itemp

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
!    Output, integer DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine dist_read ( file_name, n, dist )

!*****************************************************************************80
!
!! DIST_READ reads a distance matrix from a file.
!
!  Discussion:
!
!    The data is stored in a file, with each row of the distance
!    matrix written on one or more lines.
!
!    Blank lines and comment lines (beginning with '#') are ignored.
!
!    Individual data values should be separated by spaces or commas.
!
!  Modified:
!
!    08 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to read.
!
!    Input, integer N, the number of data items.
!
!    Output, real DIST(N,N), the distance matrix values.
!
  implicit none

  integer n

  real dist(n,n)
  character ( len = * ) file_name
  integer i
  integer ierror
  integer input
  integer ios
  integer j
  integer last
  integer length
  character ( len = 80 ) line
  integer line_num
  real value

  call get_unit ( input )

  open ( unit = input, file = file_name, status = 'old' )

  dist(1:n,1:n) = huge ( dist(1,1) )

  i = 1
  j = 0
  line_num = 0

  do
!
!  Have we read enough data?
!
    if ( i == n .and. j == n ) then
      exit
    end if
!
!  Have we read too much data?
!
    if ( i > n .or. j > n ) then
      exit
    end if
!
!  Read the next line from the file.
!
    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1
!
!  Skip blank lines and comment lines.
!
    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else
!
!  LAST points to the last character associated with the previous
!  data value read from the line.
!
      last = 0

      do
!
!  Try to read another value from the line.
!
        call s_to_r4 ( line(last+1:), value, ierror, length )
!
!  If we could not read a new value, it's time to read a new line.
!
        if ( ierror /= 0 ) then
          exit
        end if
!
!  Update the pointer.
!
        last = last + length
!
!  If we read a new value, where do we put it?
!
        j = j + 1

        if ( n < j ) then
          j = 1
          i = i + 1
          if ( n < i ) then
            exit
          end if
        end if

        dist(i,j) = value
!
!  If you reached the end of the row, it's time to read a new line.
!
        if ( j == n ) then
          exit
        end if

      end do

    end if

  end do

  close ( unit = input )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DIST_READ:'
  write ( *, '(a,i6)' ) '  Number of lines read was ', line_num

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
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 ) then

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
subroutine greedy ( kopen, no, shodi, jwarhs, rward )

!***********************************************************************
!
!! GREEDY finds a (nonoptimal) solution using a greedy algorithm.
!
!  Modified:
!
!    10 November 2002
!
!  Reference:
!
!    Gerard Cornuejols, Marshall Fisher, George Nemhauser,
!    Location of Bank Accounts to Optimize Float,
!    an Analytic Study of Exact and Approximate Algorithms,
!    Management Science,
!    Volume 23, Number 8, 1977, pages 789-810.
!
!  Parameters:
!
!    Input, integer KOPEN, the number of facilities.
!
!    Input, integer NO, the number of cities.
!
!    Input, real SHODI(NO,NO), the weighted distance matrix.
!
!    Output, integer JWARHS(NO), contains in entries 1 through KOPEN
!    the indices of the cities chosen for facilities.
!
!    Output, real RWARD, ?
!
  implicit none

  integer no

  real colm(no)
  real colto
  integer i
  integer is(no)
  integer j
  integer jwarhs(no)
  integer kopen
  integer maxidx
  integer nwar
  real rward
  real shodi(no,no)
  real uu(no)
  real yu
  real zva
!
!  Initialize UU(1:NO) to the largest weighted distance in each row.
!  Set ZVA to the total of the UU's.
!
  do i = 1, no
    uu(i) = maxval ( shodi(i,1:no) )
  end do

  zva = sum ( uu(1:no) )

  colto = 0.0E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GREEDY'
  write ( *, '(a,g14.6)' ) '  COLTO = ', colto

  nwar = 0
  is(1:no) = 0

  do

    colm(1:no) = 0.0E+00

    do j = 1, no

      if ( is(j) /= 1 ) then

        do i = 1, no
          if ( shodi(i,j) < uu(i) ) then
            colm(j) = colm(j) + uu(i) - shodi(i,j)
          end if
        end do

      end if

    end do
!
!  Find the largest entry of COLM associated with an unused city.
!
    yu = -100000.0E+00

    do j = 1, no
      if ( is(j) /= 1 ) then
        if ( yu < colm(j) ) then
          yu = colm(j)
          maxidx = j
        end if
      end if
    end do

    colto = colto + yu
    is(maxidx) = 1

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Adding city index ', maxidx
    write ( *, '(a,g14.6)' ) '  COLTO = ', colto

    nwar = nwar + 1
    jwarhs(nwar) = maxidx
!
!  Exit when we have assigned KOPEN sites.
!
    if ( nwar == kopen ) then
      exit
    end if
!
!  Is this backwards?
!
!  I think it says, if using the new city would be cheaper,
!  do it.  But I thought we would be saying, if using the new
!  city would be more expensive, then do it!
!
    do i = 1, no
      if ( shodi(i,maxidx) <= uu(i) ) then
        uu(i) = shodi(i,maxidx)
      end if
    end do

  end do

  rward = zva - colto

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GREEDY'
  write ( *, '(a,g14.6)' ) '  Maximum possible, ZVA = ', zva
  write ( *, '(a,g14.6)' ) '  Achieved, COLTO       = ', colto
  write ( *, '(a,g14.6)' ) '  Remaining RWARD       = ', rward

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
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
!    Input/output, integer I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer i
  integer j
  integer k

  k = i
  i = j
  j = k

  return
end
subroutine main_read_dist ( file_main, file_dist )

!*****************************************************************************80
!
!! MAIN_READ_DIST reads the name of the distance file from the main file.
!
!  Discussion:
!
!    FILE_DIST is the name of a file containing the city-to-city
!    distance matrix.
!
!    There MAY be a record in the main file of the form
!
!    "dist  key_dist.txt"
!
!  Modified:
!
!    07 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_MAIN, the name of the file to read.
!
!    Output, character ( len = * ) FILE_DIST, the name of the distance file,
!    or ' ' if no information was found.
!
  implicit none

  logical done
  character ( len = * ) file_dist
  character ( len = * ) file_main
  integer input
  integer ios
  character ( len = 80 ) line
  logical s_eqi
  character ( len = 80 ) word

  call get_unit ( input )

  open ( unit = input, file = file_main, status = 'old' )

  file_dist = ' '

  do
!
!  Read the next line from the file.
!
    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if
!
!  Skip blank lines and comment lines.
!
    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      done = .true.

      call word_next_read ( line, word, done )

      if ( done ) then
        cycle
      end if

      if ( .not. s_eqi ( word, 'dist' ) ) then
        cycle
      end if

      call word_next_read ( line, file_dist, done )

      exit

    end if

  end do

  close ( unit = input )

  return
end
subroutine main_read_size ( file_main, n )

!*****************************************************************************80
!
!! MAIN_READ_SIZE reads the problem size N from the main file.
!
!  Discussion:
!
!    The problem size is N, the number of cities.
!
!    There should always be a record in the main file of the form
!
!    "size  7"
!
!  Modified:
!
!    07 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_MAIN, the name of the file to read.
!
!    Output, integer N, the problem size.
!
  implicit none

  logical done
  character ( len = * ) file_main
  integer ierror
  integer input
  integer ios
  integer length
  character ( len = 80 ) line
  integer n
  logical s_eqi
  character ( len = 80 ) word

  call get_unit ( input )

  open ( unit = input, file = file_main, status = 'old' )

  n = 0

  do
!
!  Read the next line from the file.
!
    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if
!
!  Skip blank lines and comment lines.
!
    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      done = .true.

      call word_next_read ( line, word, done )

      if ( done ) then
        cycle
      end if

      if ( .not. s_eqi ( word, 'size' ) ) then
        cycle
      end if

      call word_next_read ( line, word, done )

      call s_to_i4 ( word, n, ierror, length )

      exit

    end if

  end do

  close ( unit = input )

  return
end
subroutine main_read_weight ( file_main, file_weight )

!*****************************************************************************80
!
!! MAIN_READ_WEIGHT reads the name of the weight file from the main file.
!
!  Discussion:
!
!    FILE_WEIGHT is the name of a file containing the weight
!    vector.
!
!    There MAY be a record in the main file of the form
!
!    "weight  key_weight.txt"
!
!  Modified:
!
!    09 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_MAIN, the name of the file to read.
!
!    Output, character ( len = * ) FILE_WEIGHT, the name of the weight file,
!    or ' ' if no information was found.
!
  implicit none

  logical done
  character ( len = * ) file_weight
  character ( len = * ) file_main
  integer input
  integer ios
  character ( len = 80 ) line
  logical s_eqi
  character ( len = 80 ) word

  call get_unit ( input )

  open ( unit = input, file = file_main, status = 'old' )

  file_weight = ' '

  do
!
!  Read the next line from the file.
!
    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if
!
!  Skip blank lines and comment lines.
!
    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      done = .true.

      call word_next_read ( line, word, done )

      if ( done ) then
        cycle
      end if

      if ( .not. s_eqi ( word, 'weight' ) ) then
        cycle
      end if

      call word_next_read ( line, file_weight, done )

      exit

    end if

  end do

  close ( unit = input )

  return
end
subroutine r4_swap ( x, y )

!*****************************************************************************80
!
!! R4_SWAP swaps two R4's.
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real x
  real y
  real z

  z = x
  x = y
  y = z

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Examples:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
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
  integer i
  integer len1
  integer len2
  integer lenc
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
subroutine sol ( shodi, isw, itera, jsrt, jwarhs, kopen, lpt_unit, &
  no, rward, nv, nvsq, u, w, w_best )

!***********************************************************************
!
!! SOL takes one step of the iteration.
!
!  Modified:
!
!    04 November 2002
!
!  Reference:
!
!    Gerard Cornuejols, Marshall Fisher, George Nemhauser,
!    Location of Bank Accounts to Optimize Float,
!    an Analytic Study of Exact and Approximate Algorithms,
!    Management Science,
!    Volume 23, Number 8, 1977, pages 789-810.
!
  implicit none

  integer no

  real compar
  integer i
  integer icjmax(no)
  integer indxj
  integer indx
  real shodi(no,no)
  integer isw
  integer itera
  integer itr
  integer j
  integer jj
  integer jjmax
  integer jk
  integer jsrt(no,no)
  integer jwarhs(no)
  integer k
  integer kopen
  integer ks(no)
  integer l
  integer lpt_unit
  integer mcomp
  integer nextj
  integer nreal
  integer nv(no)
  integer nvsq
  real rowj(no)
  real row_sum(no)
  real rward
  real sumind
  real u(no)
  real w
  real w_best
  real ww

  row_sum(1:no) = 0.0E+00

  do i = 1, no

    do itr = 1, no

      nextj = jsrt(i,itr)
      sumind = shodi(i,nextj) + u(i)

      if ( 0.0E+00 <= sumind ) then
        exit
      end if

      row_sum(nextj) = row_sum(nextj) + sumind

    end do

  end do

  rowj(1:no) = row_sum(1:no)

  do k = 1, kopen
    compar = 1000000.0E+00
    do j = 1, no
      if ( rowj(j) < compar ) then
        compar = rowj(j)
        indxj = j
      end if
    end do
    icjmax(k) = indxj
    rowj(indxj) = 1000000.0E+00
  end do

  ww = 0.0E+00
  do l = 1, kopen
    ww = ww + row_sum(icjmax(l))
  end do

  do i = 1, no
    ww = ww - u(i)
  end do

  w = ww

  if ( w_best - 0.001E+00 < w .and. w < w_best + 0.001E+00 ) then
    go to 72
  end if

  indx = 0
  go to 73

72    continue

  if ( kopen == no-1 ) then
    go to 74
  end if

  indx = indx + 1

73    continue

  if ( 20 <= indx ) then
    write ( lpt_unit, '(a)' ) ' '
    write ( lpt_unit, '(a)' ) 'SOL:'
    write ( lpt_unit, '(a)' ) '  The relaxation converges to a specific value.'
    isw = 3
    go to 27
  end if

74    continue

  if ( kopen /= 0 ) then

    nreal = 0
    do i = 1, no
      jjmax = 1000000
      do jj = 1, kopen
        jk = icjmax(jj)
        if ( shodi(i,jk) < jjmax ) then
          jjmax = shodi(i,jk)
        end if
      end do
      nreal = nreal + jjmax
    end do

    if ( nreal < rward ) then
      rward = nreal
      jwarhs(1:kopen) = icjmax(1:kopen)
    end if

  end if

  if ( w_best + 0.9E+00 <= rward ) then
    isw = 1
    go to 27
  end if

  do i = 1, no
    nv(i) = -1
    do j = 1, kopen
      if ( shodi(i,icjmax(j)) + u(i) <= 0.0E+00 ) then
        nv(i) = nv(i) + 1
      end if
    end do
  end do

  nvsq = sum ( nv(1:no)**2 )

  itera = itera + 1

  if ( nvsq == 0 ) then
    isw = 2
    go to 27
  end if

  if ( itera <= 300 ) then
    return
  end if

  isw = 4

  if ( w_best < w ) then
    w_best = w
  end if

  write ( lpt_unit, * ) 'The cost of the best solution found is:', rward
  go to 289

27 continue

  w_best = max ( w_best, w )

289   continue

  if ( isw /= 1 .and. isw /= 2 ) then
    write ( lpt_unit, '(a)' ) ' '
    write ( lpt_unit, '(a,g14.6)' ) '  Best feasible solution', rward
  else
    write ( lpt_unit, '(a)' ) ' '
    write ( lpt_unit, '(a,g14.6)' ) '  Cost of optimal solution is : ', rward
  end if

  write ( lpt_unit, '(a)' ) ' '
  write ( lpt_unit, 104 ) jwarhs(1:kopen)
104   format(' the stocks in the index fund are:'/20(1x,i4))
  write ( lpt_unit, '(a)' ) ' '
  write ( lpt_unit, '(a)' ) 'fund stock  tracks  target stocks'

  do i = 1, no
    mcomp = 10000000
    do j = 1, kopen
      if ( shodi(i,jwarhs(j)) <= mcomp ) then
        mcomp = shodi(i,jwarhs(j))
        ks(i) = jwarhs(j)
      end if
    end do
    write ( lpt_unit, '(a)' ) ' '
    write ( lpt_unit, 329 ) ks(i), i
329   format(i3,'            ----->    ',i5)
  end do

  if ( 300 < itera ) then
    write ( lpt_unit,*)' after 300 iterations lagrangian bound is:',w_best
  else
    write ( lpt_unit,*)' number of iterations:',itera
  end if

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
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
!    Output, integer IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LAST, the last character of S used to make IVAL.
!
  implicit none

  character c
  integer i
  integer ierror
  integer isgn
  integer istate
  integer ival
  integer last
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
!    Output, real R, the real value that was read from the string.
!
!    Output, integer IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character c
  logical ch_eqi
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer lchar
  integer nchar
  integer ndig
  real r
  real rbot
  real rexp
  real rtop
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
      else if ( ihave > 1 ) then
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
      else if ( ihave >= 6 .and. ihave <= 8 ) then
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
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
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
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
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
subroutine weight_read ( file_name, n, weight )

!*****************************************************************************80
!
!! WEIGHT_READ reads weights from a file.
!
!  Discussion:
!
!    The data is stored in a file, with each weight on one line.
!
!    Blank lines and comment lines (beginning with '#') are ignored.
!
!    Individual data values should be separated by spaces or commas.
!
!  Modified:
!
!    08 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to read.
!
!    Input, integer N, the number of data items.
!
!    Output, real WEIGHT(N), the weights.
!
  implicit none

  integer n
  integer, parameter :: ncol = 1

  character ( len = * ) file_name
  integer i
  integer ierror
  integer input
  integer ios
  integer j
  integer last
  integer length
  character ( len = 80 ) line
  integer line_num
  real value
  real weight(n)

  call get_unit ( input )

  open ( unit = input, file = file_name, status = 'old' )

  weight(1:n) = huge ( weight(1) )

  i = 1
  j = 0
  line_num = 0

  do
!
!  Have we read enough data?
!
    if ( i == n .and. j == ncol ) then
      exit
    end if
!
!  Have we read too much data?
!
    if ( n < i .or. ncol < j ) then
      exit
    end if
!
!  Read the next line from the file.
!
    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1
!
!  Skip blank lines and comment lines.
!
    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else
!
!  LAST points to the last character associated with the previous
!  data value read from the line.
!
      last = 0

      do
!
!  Try to read another value from the line.
!
        call s_to_r4 ( line(last+1:), value, ierror, length )
!
!  If we could not read a new value, it's time to read a new line.
!
        if ( ierror /= 0 ) then
          exit
        end if
!
!  Update the pointer.
!
        last = last + length
!
!  If we read a new value, where do we put it?
!
        j = j + 1

        if ( ncol < j ) then
          j = 1
          i = i + 1
          if ( n < i ) then
            exit
          end if
        end if

        weight(i) = value
!
!  If you reached the end of the row, it's time to read a new line.
!
!       if ( j == ncol ) then
!         exit
!       end if

      end do

    end if

  end do

  close ( unit = input )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WEIGHT_READ:'
  write ( *, '(a,i6)' ) '  Number of lines read was ', line_num

  return
end
subroutine word_next_read ( s, word, done )

!*****************************************************************************80
!
!! WORD_NEXT_READ "reads" words from a string, one at a time.
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
  integer ilo
  integer, save :: lenc = 0
  integer, save :: next = 1
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

    if ( ilo > lenc ) then
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
