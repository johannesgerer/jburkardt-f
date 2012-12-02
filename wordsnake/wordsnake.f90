program main

!*****************************************************************************80
!
!! MAIN is the main program for WORDSNAKE.
!
!  Discussion:
!
!    WORDSNAKE computes a good wordsnake from a file of words.
!
!  Usage:
!
!    wordsnake wordsnake.inp
!
!    where "wordsnake.inp" is a list of words, one word per line,
!    to be made into a wordsnake.
!
!  Best so far:
!
!    0    0    0  sea  invent
!    1    1    1  inven T errible
!    3    9   10  terri BLE mish
!    4   16   26  ble MISH apen
!    3    9   35  misha PEN ultimate
!    2    4   39  penultima TE nse
!    2    4   43  ten SE em
!    2    4   47  se EM erge
!    5   25   72  e MERGE r
!    3    9   81  mer GER iatric
!    4   16   97  geria TRIC ky
!    1    1   98  trick Y es
!    2    4  102  y ES sential
!    2    4  106  essenti AL ly
!    1    1  107  all Y et
!    2    4  111  y ET ernal
!    2    4  115  etern AL as
!    3    9  124  a LAS ting
!    5   25  149  la STING er
!    3    9  158  stin GER und
!    3    9  167  ger UND erdevelop
!    4   16  183  underdev ELOP ed
!    3    9  192  elo PED iatric
!    4   16  208  pedia TRIC e
!    3    9  217  tr ICE r
!    3    9  226  i CER tain
!    2    4  230  certa IN credible
!    3    9  239  incredi BLE nd
!    4   16  255  b LEND ing
!    3    9  264  lend ING rate
!    4   16  280  ing RATE s
!    3    9  289  ra TES silate
!    4   16  305  tessi LATE r
!    3    9  314  la TER restrial
!    5   25  339  terres TRIAL s
!    1    1  340  trial S udden
!    3    9  349  sud DEN ude
!    2    4  353  denu DE nse
!    2    4  357  den SE a
!
!    Total number of characters =  254
!    Reduced number of characters =  145
!    Score is  357
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 March 2003
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 39

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilen
  character ( len = 80 ) input_file_name  
  integer ( kind = 4 ) ipxfargc
  integer ( kind = 4 ) numarg
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) record_num
  integer ( kind = 4 ) score
  character ( len = 80 ), allocatable, dimension ( : ) :: word
!
!  Say hello.
!
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WORDSNAKE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Call WORDSNAKE to create a high scoring wordsnake'
  write ( *, '(a)' ) '  from a given set of words.'
!
!  Get the name of the input file.
!
  numarg = iargc ( )
!
!  Get the input file name.
!
  if ( numarg >= 1 ) then

    iarg = 1
    call getarg ( iarg, input_file_name )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'What is the name of the input file?'
    read ( *, '(a)' ) input_file_name
    if ( input_file_name == ' ' ) then
      stop
    end if

  end if
!
!  Count the number of records in the input file.
!
  call file_record_count ( input_file_name, record_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of records in the input file is ', &
    record_num
!
!  Allocate the word array.
!
  allocate ( word(1:record_num) )
!
!  Read the word array.
!
  call file_record_read ( input_file_name, record_num, word )
!
!  Print the word array.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The word list:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i2,2x,a)' ) i, trim ( word(i) )
  end do
!
!  Try to make the word snake.
!
  do i = 1, 3

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Trial number ', i

    call wordsnake ( n, word, perm )

    call wordsnake_print ( n, word, perm )

    call wordsnake_score ( n, word, perm, score )

    write ( *, '(a,i6)' ) '  The wordsnake score is ', score

  end do
!
!  Free memory.
!
  deallocate ( word )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WORDSNAKE'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine file_record_count ( file_in_name, record_num )

!*****************************************************************************80
!
!! FILE_RECORD_COUNT counts the number of records in a file.
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
!    11 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_IN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) RECORD_NUM, the number of records found.
!
  implicit none

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) comment_num
  character ( len = * ) file_in_name
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) record_num

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_RECORD_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( file_in_name )
    stop
  end if

  comment_num = 0
  record_num = 0
  line_num = 0
  bad_num = 0

  do

    read ( file_in_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = line_num
      exit
    end if

    line_num = line_num + 1

    if ( line(1:1) == '#' ) then
      comment_num = comment_num + 1
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      comment_num = comment_num + 1
      cycle
    end if

    record_num = record_num + 1

  end do

  close ( unit = file_in_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FILE_RECORD_COUNT:'
  write ( *, '(a,i6)' ) '  Number of lines:           ', line_num
  write ( *, '(a,i6)' ) '  Number of data records:    ', record_num
  write ( *, '(a,i6)' ) '  Number of comment records: ', comment_num

  return
end
subroutine file_record_read ( file_in_name, record_num, word )

!*****************************************************************************80
!
!! FILE_RECORD_READ reads the records in a file.
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
!    11 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_IN_NAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) RECORD_NUM, the number of records.
!
  implicit none

  integer ( kind = 4 ) record_num

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) comment_num
  character ( len = * ) file_in_name
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) record_num2
  character ( len = 80 ) word(record_num)

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_RECORD_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( file_in_name )
    stop
  end if

  record_num2 = 0

  do 

    read ( file_in_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = line_num
      exit
    end if

    line_num = line_num + 1

    if ( line(1:1) == '#' ) then
      comment_num = comment_num + 1
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      comment_num = comment_num + 1
      cycle
    end if

    record_num2 = record_num2 + 1

    if ( record_num2 <= record_num ) then
      word(record_num2) = line
    end if

  end do

  close ( unit = file_in_unit )

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Output, integer ( kind = 4 ) IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
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
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!  Example:
!
!        I     J     MOD  I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i6)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

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

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  real ( kind = 4 ) r
  real ( kind = 4 ) rhi
  real ( kind = 4 ) rlo
  logical, save :: seed = .false.

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
  rlo = real ( ilo, kind = 4 ) - 0.5E+00
  rhi = real ( ihi, kind = 4 ) + 0.5E+00
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
subroutine i4_swap3 ( i, j, k )

!*****************************************************************************80
!
!! I4_SWAP3 swaps three integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J, K.  On output, the values of I, J, and K
!    have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l

  l = i
  i = j
  j = k
  k = l

  return
end
subroutine i4_unswap3 ( i, j, k )

!*****************************************************************************80
!
!! I4_UNSWAP3 unswaps three integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J, K.  On output, the values of I, J, and K
!    have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l

  l = k
  k = j
  j = i
  i = l

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I4_WRAP
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, an integer value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds for the integer value.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) wide

  wide = ihi + 1 - ilo

  if ( wide == 0 ) then
    i4_wrap = ilo
  else
    i4_wrap = ilo + i4_modp ( ival-ilo, wide )
  end if

  return
end
subroutine i4vec_identity ( n, a )

!*****************************************************************************80
!
!! I4VEC_IDENTITY sets an I4VEC to the identity vector A(I)=I.
!
!  Modified:
!
!    09 November 2000
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
function lower ( s )

!*****************************************************************************80
!
!! LOWER returns a lowercase version of a string.
!
!  Discussion:
!
!    LOWER is a string function of undeclared length.  The length
!    of the argument returned is determined by the declaration of
!    LOWER in the calling routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string.
!
!    Output, character ( len = * ) LOWER, a lowercase copy of the string.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character ( len = * ) lower
  integer ( kind = 4 ) n
  character ( len = * ) s

  lower = s

  n = len_trim ( lower )

  do i = 1, n

    j = ichar ( lower(i:i) )

    if ( 65 <= j .and. j <= 90 ) then
      lower(i:i) = char ( j + 32 )
    end if

  end do

  return
end
subroutine overlap_table ( n, word, table )

!*****************************************************************************80
!
!! OVERLAP_TABLE computes a table of the overlap between pairs of words.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of words.
!
!    Input, character ( len = * ) WORD(N), a list of words.
!
!    Output, integer ( kind = 4 ) TABLE(N,N), a table containing, in TABLE(I,J),
!    the number of characters by which the end of WORD(I) and the
!    beginning of word J overlap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) table(n,n)
  character ( len = * ) word(n)
!
!  Construct the overlap table.
!
  do i = 1, n
    do j = 1, n

      call s_overlap ( word(i), word(j), table(i,j) )

    end do
  end do

  return
end
subroutine perm_random ( n, p, setup )

!*****************************************************************************80
!
!! PERM_RANDOM selects a random permutation of N objects.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!
!    Input/output, integer ( kind = 4 ) P(N), a permutation, in standard index form.
!
!    If SETUP is .TRUE., then the input value of P contains
!    the "current" labels of the objects.
!
!    Otherwise, P(I) is initialized to I.
!
!    On output, P(I) is the randomly permuted label of the I-th object.
!
!    Input, logical SETUP.
!
!    If SETUP is .TRUE. then the routine assumes the objects
!    are labeled 1, 2, ... N.
!
!    If SETUP is .FALSE., then the input values of P are used
!    as labels; that is, the I-th object is labeled P(I).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(n)
  logical setup

  if ( setup ) then
    call i4vec_identity ( n, p )
  end if

  do i = 1, n
    call i4_random ( i, n, j )
    call i4_swap ( p(i), p(j) )
  end do

  return
end
subroutine s_overlap ( s1, s2, overlap )

!*****************************************************************************80
!
!! S_OVERLAP determines the overlap between two strings.
!
!  Discussion:
!
!    To determine the overlap, write the first word followed immediately
!    by the second word.  Find the longest substring S which is both
!    a suffix of S1 and a prefix of S2.  The length of this substring
!    is the overlap.
!
!  Example:
!
!    S1              S2        OVERLAP
!
!    'timber'        'beret'   3
!    'timber'        'timber'  6
!    'beret'         'timber'  1
!    'beret'         'berets'  5
!    'beret'         'berth'   0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to be checked.
!
!    Output, integer ( kind = 4 ) OVERLAP, the length of the overlap.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) len3
  integer ( kind = 4 ) overlap
  character ( len = * ) s1
  character ( len = * ) s2

  overlap = 0

  len1 = len_trim ( s1 )
  len2 = len_trim ( s2 )
  len3 = min ( len1, len2 )

  do i = 1, len3
    if ( s1(len1+1-i:len1) == s2(1:i) ) then
      overlap = i
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
  character ( len = 10 ) time
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
function upper ( s )

!*****************************************************************************80
!
!! UPPER returns an uppercase version of a string.
!
!  Discussion:
!
!    UPPER is a string function of undeclared length.  The length
!    of the argument returned is determined by the declaration of
!    UPPER in the calling routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string.
!
!    Output, character ( len = * ) UPPER, an uppercase copy of the string.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  character ( len = * ) s
  character ( len = * ) upper

  upper = s

  n = len_trim ( upper )

  do i = 1, n

    j = ichar ( upper(i:i) )

    if ( 97 <= j .and. j <= 122 ) then
      upper(i:i) = char ( j - 32 )
    end if

  end do

  return
end
subroutine wordsnake ( n, word, perm )

!*****************************************************************************80
!
!! WORDSNAKE seeks a high scoring permutation of a set of words.
!
!  Discussion:
!
!    Note that MISHAPEN and TESSILATE, as given in the article,
!    are incorrectly spelled.
!
!    A wordsnake is formed from a list of words.  The words are rearranged
!    so that, where possible, the end of one word matches the beginning of
!    the next.  The highest scoring wordsnake is sought.  This is a hard
!    problem, and is similar, really, to the traveling salesman problem.
!
!  Scoring:
!
!    The score for a wordsnake is the sum of the squares of the number of
!    letters of overlap between pairs of successive words.  (The article
!    does not mention whether to count overlap for the possible wrap-around
!    from the last word to the first again.)
!
!    For instance, consider the wordsnake made up of:
!
!      "writes", "testate", "tater", "remote", "stew"
!
!    and consider the resulting wordsnake:
!
!      writestateremotestew
!
!    which is scored as:
!
!      Word pair     Overlap  Score
!      ------------  -------  -----
!      wri TES tate     3        +9
!      tes TATE r       4       +16
!      tate R emote     1        +1
!      remote stew      0        +0
!      ste W rites      1        +1
!
!      TOTAL SCORE               27
!      TOTAL CHARACTERS          20
!
!    Using this scoring system, and the set of 39 words given below, the 
!    article claims a score of 357 was achieved, with a wordsnake length of 
!    145 characters (that is, printing out the shared letters only once).
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
!  Reference:
!
!    Dennis Shasha,
!    Wordsnakes,
!    Dr Dobb's Journal,
!    July, 2000, pages 143-144.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of words.
!
!    Input, character ( len = * ) WORD(N), the words to be used.
!
!    Output, integer ( kind = 4 ) PERM(N), a permutation of the words to be used
!    in the wordsnake.
!
  implicit none

  integer ( kind = 4 ) n

  logical, parameter :: DEBUG = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_best(n)
  integer ( kind = 4 ) score
  integer ( kind = 4 ) score_best
  integer ( kind = 4 ) score_best_old
  integer ( kind = 4 ) table(n,n)
  character ( len = * ), dimension ( n ) :: word
!
!  Compute the overlap table.
!
  call overlap_table ( n, word, table )

  if ( DEBUG ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Overlap table:'
    write ( *, '(a)' ) ' '
    write ( *, '(6x,39i2)' ) ( i, i = 1, n )
    do i = 1, n
      write ( *, '(4x,40i2)' ) i, table(i,1:n)
   end do

  end if
!
!  Wordsnakes will be described by a permutation of our word list.
!  Our initial wordsnake is the identity permutation.
!
  call perm_random ( n, perm, .true. )

  if ( DEBUG ) then
    call wordsnake_print ( n, word, perm )
  end if

  call wordsnake_score ( n, word, perm, score )

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'Initial score is ', score
  end if

  perm_best(1:n) = perm(1:n)
  score_best = score
!
!  Do one naive "greedy" pass, where we follow word(I) by its optimal follower.
!
  call wordsnake_search_greedy ( n, word, perm, table )

  if ( DEBUG ) then
    call wordsnake_print ( n, word, perm )
  end if

  call wordsnake_score ( n, word, perm, score )

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'Greedy score is ', score
  end if

  if ( score > score_best ) then
    perm_best(1:n) = perm(1:n)
    score_best = score
  end if

  do

    score_best_old = score_best
!
!  Consider all possible inserts.
!
    call wordsnake_search_insert ( n, word, perm_best, score_best )
 
    if ( DEBUG ) then
      call wordsnake_print ( n, word, perm_best )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'INSERT score is ', score_best
    end if
!
!  Consider all possible pair swaps.
!
    call wordsnake_search_swap2 ( n, word, perm_best, score_best )
 
    if ( DEBUG ) then
      call wordsnake_print ( n, word, perm_best )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'SWAP2 score is ', score_best
    end if
!
!  Consider all possible trio swaps.
!
    call wordsnake_search_swap3 ( n, word, perm_best, score_best )
 
    if ( DEBUG ) then
      call wordsnake_print ( n, word, perm_best )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'SWAP3 score is ', score_best
    end if
!
!  Now consider all switches of the form -A-B-C-A- to -A-C-B-A-
!
    call wordsnake_search_transpose ( n, word, perm_best, score_best )

    if ( DEBUG ) then
      call wordsnake_print ( n, word, perm_best )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'TRANSPOSE score is ', score_best
    end if

    if ( score_best == score_best_old ) then
      exit
    end if

  end do

  return
end
subroutine wordsnake_print ( n, word, perm )

!*****************************************************************************80
!
!! WORDSNAKE_PRINT prints a wordsnake.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of words.
!
!    Input, character ( len = * ) WORD(N), a list of words.
!
!    Input, integer ( kind = 4 ) PERM(N), a permutation, the ordering of the
!    words for the wordsnake.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) inc
  character ( len = 25 ) lower
  integer ( kind = 4 ) na
  integer ( kind = 4 ) n_char
  integer ( kind = 4 ) n_overlap
  integer ( kind = 4 ) n_reduce
  integer ( kind = 4 ) overlap
  integer ( kind = 4 ) perm(n)
  character ( len = 25 ) s1
  character ( len = 25 ) s2
  character ( len = 25 ) sa
  character ( len = 25 ) sb
  character ( len = 25 ) sc
  integer ( kind = 4 ) score
  character, parameter :: space = ' '
  character ( len = 25 ) upper
  character ( len = * ) word(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WORDSNAKE'
  write ( *, '(a,i8)' ) '  Number of words = ', n
  write ( *, '(a)' ) ' '

  score = 0
  n_overlap = 0

  do i = 1, n

    s1 = word(perm(i))

    if ( i < n ) then
      s2 = word(perm(i+1))
    else
      s2 = word(perm(1))
    end if

    call s_overlap ( s1, s2, overlap )

    n_overlap = n_overlap + overlap

    inc = overlap**2
    score = score + inc 

    na = len_trim ( s1 )
    sa = lower ( s1(1:na-overlap) )
    sb = upper ( s2(1:overlap) )
    sc = lower ( s2(overlap+1:) )

    write ( *, '(3i8,2x,5a)' ) overlap, inc, score, &
      trim ( sa ), space, trim ( sb ), space, trim ( sc )

  end do

  n_char = 0
  do i = 1, n
    n_char = n_char + len_trim ( word(i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Total number of characters = ', n_char
  write ( *, '(a,i8)' ) '  Reduced number of characters = ', n_char - n_overlap

  return
end
subroutine wordsnake_score ( n, word, perm, score )

!*****************************************************************************80
!
!! WORDSNAKE_SCORE computes the score for a given wordsnake.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of words.
!
!    Input, character ( len = * ) WORD(N), a list of words.
!
!    Input, integer ( kind = 4 ) PERM(N), a permutation, the ordering of the
!    words for the wordsnake.
!
!    Output, integer ( kind = 4 ) SCORE, the overlap score for this wordsnake.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) overlap
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) score
  character ( len = * ) word(n)

  score = 0

  do i = 1, n

    j = perm(i)

    if ( i < n ) then
      ip1 = i + 1
    else
      ip1 = 1
    end if

    jp1 = perm(ip1)

    call s_overlap ( word(j), word(jp1), overlap )

    score = score + overlap**2

  end do

  return
end
subroutine wordsnake_search_greedy ( n, word, perm, table )

!*****************************************************************************80
!
!! WORDSNAKE_SEARCH_GREEDY constructs a wordsnake using a greedy algorithm.
!
!  Discussion:
!
!    The routine takes the first word, and puts it into the wordsnake.
!    Then it takes the next word that is available, and would produce the
!    best matching score with the previous word, and so on.
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
!    Input, integer ( kind = 4 ) N, the number of words.
!
!    Input, character ( len = * ) WORD(N), a list of words.
!
!    Input/output, integer ( kind = 4 ) PERM(N), a permutation, the ordering of the
!    words for the wordsnake.
!
!    Input, TABLE(N,N), a table of the overlap scores for each pair of words.
!
  implicit none

  integer ( kind = 4 ) :: n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) table(n,n)
  character ( len = * ), dimension ( n ) :: word

  do i = 1, n-1

    ip = perm(i)

    j = i+1
    jp = perm(i+1)

    do k = i+2, n
      kp = perm(k)

      if ( table(ip,kp) > table(ip,jp) ) then
        j = k
        jp = kp
      end if

    end do

    call i4_swap ( perm(i+1), perm(j) )

  end do

  return
end
subroutine wordsnake_search_insert ( n, word, perm_best, score_best )

!*****************************************************************************80
!
!! WORDSNAKE_SEARCH_INSERT tries to improve the score by inserting a word.
!
!  Discussion:
!
!    The operation considered here modifies the string by moving a single
!    word from one position to another.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of words.
!
!    Input, character ( len = * ) WORD(N), a list of words.
!
!    Input/output, integer ( kind = 4 ) PERM_BEST(N), a permutation, the ordering of the
!    words for the wordsnake.
!
!    Input/output, integer ( kind = 4 ) SCORE_BEST, the best score so far.
!
  implicit none

  integer ( kind = 4 ) n

  logical, parameter :: DEBUG = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) improve
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n_insert
  integer ( kind = 4 ) pi_new
  integer ( kind = 4 ) pi_old
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_best(n)
  integer ( kind = 4 ) score
  integer ( kind = 4 ) score_best
  character ( len = * ) word(n)

  n_insert = 0
  perm(1:n) = perm_best(1:n)
  improve = 0
!
!  For word I, which is now in position PI_OLD = PERM(I), consider inserting
!  it in position PI_NEW.
!
  do i = 1, n

    do pi_new = 1, n

      pi_old = perm(i)

      if ( pi_new /= pi_old ) then

        if ( pi_new < pi_old ) then
          do ii = 1, n
            if ( pi_new <= perm(ii) .and. perm(ii) < pi_old ) then
              perm(ii) = perm(ii) + 1
            else if ( perm(ii) == pi_old ) then
              perm(ii) = pi_new
            end if
          end do
        else if ( pi_new > pi_old ) then
          do ii = 1, n
            if ( pi_old < perm(ii) .and. perm(ii) <= pi_new ) then
              perm(ii) = perm(ii) - 1
            else if ( perm(ii) == pi_old ) then
              perm(ii) = pi_new
            end if
          end do
        end if

        call wordsnake_score ( n, word, perm, score )

        if ( score >= score_best ) then
          improve = improve + score - score_best
          score_best = score
          perm_best(1:n) = perm(1:n)
          n_insert = n_insert + 1
        else
          perm(1:n) = perm_best(1:n)
        end if

      end if
    end do
  end do

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Number of inserts = ', n_insert
    write ( *, '(a,i6)' ) '  Score improvement = ', improve
  end if

  return
end
subroutine wordsnake_search_swap2 ( n, word, perm_best, score_best )

!*****************************************************************************80
!
!! WORDSNAKE_SEARCH_SWAP2 tries to improve the score by swapping 2 words.
!
!  Discussion:
!
!    The operation essentially checks the possibility that the score would
!    be improved by the following 2-Swap operation:
!
!      Temp  <- Word1
!      Word1 <- Word2
!      Word2 <- Temp.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of words.
!
!    Input, character ( len = * ) WORD(N), a list of words.
!
!    Input/output, integer ( kind = 4 ) PERM_BEST(N), a permutation, the ordering of the
!    words for the wordsnake.
!
!    Input/output, integer ( kind = 4 ) SCORE_BEST, the best score so far.
!
  implicit none

  integer ( kind = 4 ) n

  logical, parameter :: DEBUG = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) improve
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n_swap
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_best(n)
  integer ( kind = 4 ) score
  integer ( kind = 4 ) score_best
  character ( len = * ) word(n)

  n_swap = 0
  perm(1:n) = perm_best(1:n)
  improve = 0
!
!  For word I, which is now in position PERM(I), consider swapping
!  it with the word in position PERM(J).
!
  do i = 1, n
    do j = 1, n
      if ( i /= j ) then
        call i4_swap ( perm(i), perm(j) )
        call wordsnake_score ( n, word, perm, score )
        if ( score >= score_best ) then
          improve = improve + score - score_best
          score_best = score
          perm_best(1:n) = perm(1:n)
          n_swap = n_swap + 1
        else
          call i4_swap ( perm(i), perm(j) )
        end if
      end if
    end do
  end do

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Number of swaps = ', n_swap
    write ( *, '(a,i6)' ) '  Score improvement = ', improve
  end if

  return
end
subroutine wordsnake_search_swap3 ( n, word, perm_best, score_best )

!*****************************************************************************80
!
!! WORDSNAKE_SEARCH_SWAP3 tries to improve the score by swapping 3 words.
!
!  Discussion:
!
!    The operation essentially checks the possibility that the score would
!    be improved by the following 3-Swap operation:
!
!      Temp  <- Word1
!      Word1 <- Word2
!      Word2 <- Word3
!      Word3 <- Temp.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of words.
!
!    Input, character ( len = * ) WORD(N), a list of words.
!
!    Input/output, integer ( kind = 4 ) PERM_BEST(N), a permutation, the ordering of the
!    words for the wordsnake.
!
!    Input/output, integer ( kind = 4 ) SCORE_BEST, the best score so far.
!
  implicit none

  integer ( kind = 4 ) n

  logical, parameter :: DEBUG = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) improve
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n_swap
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_best(n)
  integer ( kind = 4 ) score
  integer ( kind = 4 ) score_best
  character ( len = * ) word(n)

  n_swap = 0
  perm(1:n) = perm_best(1:n)
  improve = 0
!
!  For word I, which is now in position PERM(I), consider putting it
!  where PERM(J) is.
!
  do i = 1, n
    do j = 1, n
      if ( i /= j ) then
        do k = 1, n
          if ( k /= i .and. k /= j ) then
            call i4_swap3 ( perm(i), perm(j), perm(k) )
            call wordsnake_score ( n, word, perm, score )
            if ( score >= score_best ) then
              improve = improve + score - score_best
              score_best = score
              perm_best(1:n) = perm(1:n)
              n_swap = n_swap + 1
            else
              call i4_unswap3 ( perm(i), perm(j), perm(k) )
            end if
          end if
        end do
      end if
    end do
  end do

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Number of 3-swaps = ', n_swap
    write ( *, '(a,i6)' ) '  Score improvement = ', improve
  end if

  return
end
subroutine wordsnake_search_transpose ( n, word, perm_best, score_best )

!*****************************************************************************80
!
!! WORDSNAKE_SEARCH_TRANSPOSE tries to improve the score using transpositions.
!
!  Discussion:
!
!    In a transposition, the wordsnake is essentially cut in three places,
!    and two pieces are swapped.  We can think of the operation as 
!    replacing:
!
!      -A-B-C-A-
!
!    by
!
!      -A-C-B-A-
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
!    Input, integer ( kind = 4 ) N, the number of words.
!
!    Input, character ( len = * ) WORD(N), a list of words.
!
!    Input/output, integer ( kind = 4 ) PERM_BEST(N), a permutation, the ordering of the
!    words for the wordsnake.
!
!    Input/output, integer ( kind = 4 ) SCORE_BEST, the best score so far.
!
  implicit none

  integer ( kind = 4 ) n

  logical, parameter :: DEBUG = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) improve
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n_transpose
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm2(n)
  integer ( kind = 4 ) perm_best(n)
  integer ( kind = 4 ) score
  integer ( kind = 4 ) score_best
  character ( len = * ) word(n)

  n_transpose = 0
  perm(1:n) = perm_best(1:n)
  improve = 0
!
!  Consider I0:I1:I2:I3:I4:I5 and transposing to I0:I3:I4:I1:I2:I5.
!
  do i1 = 1, n-1
    do i2 = i1, n-1
      do i4 = i2+1, n

        i3 = i2+1

        perm2(1:n) = perm(1:n)
        perm2(i1:i4+i1-i3) = perm(i3:i4)
        perm2(i4+i1-i2:i4) = perm(i1:i2)

        call wordsnake_score ( n, word, perm2, score )

        if ( score >= score_best ) then
          improve = improve + score - score_best
          score_best = score
          perm_best(1:n) = perm2(1:n)
          n_transpose = n_transpose + 1
        end if

      end do
    end do
  end do

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Number of transpositions = ', n_transpose
    write ( *, '(a,i6)' ) '  Score improvement = ', improve
  end if

  return
end
