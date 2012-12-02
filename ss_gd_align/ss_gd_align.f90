subroutine a_index ( acid_num, acid_code, acid_index )

!*****************************************************************************80
!
!! A_INDEX sets up a reverse index for the amino acid codes.
!
!  Example:
!
!    Input:
!
!      ACID_CODE = 
!        'A', 'R', 'N', 'B', 'D', 'C', 'Q', 'Z', 'E', 'G',
!        'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W',
!        'Y', 'V', 'X'
!
!    Output:
!
!      ACID_INDEX = 
!        1,   4,   6,   5,   9,  16,  10,  11,  12,   0,
!       14,  13,  15,   3,   0,  17,   7,   2,  18,  19,
!        0,  22,  20,  23,  21,   8.
!
!  Discussion:
!
!    ACID_CODE allows us to discover the index of item 'R' only
!    by searching through all the entries of ACID_CODE until we 
!    encounter an 'R' or reach the end of the array.
!
!    ACID_INDEX allows us to discover the index of item 'R' by
!    converting 'R' to its numeric index of 17 and evaluating
!    ACID_INDEX(17), which is 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ACID_NUM, the number of entries in ACID_CODE.
!
!    Input, character ACID_CODE(ACID_NUM), a list of alphabetic
!    characters.  Normally, this list is upper case, and contains a
!    subset of the alphabetic characters 'A' through 'Z'.
!
!    Output, integer ACID_INDEX(26), indicates, for each alphabetic
!    character, the (last) index of ACID_CODE containing that character,
!    or 0 if the character does not occur in ACID_CODE.
!
  implicit none

  integer acid_num

  character a
  integer a_to_i4
  character acid_code(acid_num)
  integer acid_index(26)
  integer i
  integer j

  acid_index(1:26) = 0

  do i = 1, acid_num

    a = acid_code(i)
    j = a_to_i4 ( a )

    if ( j < 1 .or. j > 26 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'A_INDEX - Fatal error!'
      write ( *, '(a,i6)' ) '  Out-of-bounds acid index J = ', j
      write ( *, '(a)' ) '  Originating acid code character is ' // a
      write ( *, '(a,i6)' ) '  ASCII code = ', ichar ( a )
      stop
    end if

    acid_index(j) = i

  end do

  return
end
function a_to_i4 ( a )

!*****************************************************************************80
!
!! A_TO_I4 returns the index of an alphabetic character.
!
!  Example:
!
!    A  A_TO_I4
!
!    'a'   1
!    'b'   2
!    ...
!    'z'  26
!    'A'   1
!    'B'   2
!    ...   
!    'Z'  26
!    '$'  -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character A, a character.
!
!    Output, integer A_TO_I4, is the alphabetic index of the character,
!    between 1 and 26 if the character is alphabetic, and -1 otherwise.
!
  implicit none

  character a
  integer a_to_i4

  if ( lle ( 'A', a ) .and. lle ( a, 'Z' ) ) then
    a_to_i4 = 1 + ichar ( a ) - ichar ( 'A' )
  else if ( lle ( 'a', a ) .and. lle ( a, 'z' ) ) then
    a_to_i4 = 1 + ichar ( a ) - ichar ( 'a' )
  else
    a_to_i4 = -1
  end if

  return
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
  integer itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
subroutine chvec2_print ( m, a, n, b, title )

!*****************************************************************************80
!
!! CHVEC2_PRINT prints two vectors of characters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the length of the first sequence.
!
!    Input, character A(M), the first sequence.
!
!    Input, integer N, the length of the second sequence.
!
!    Input, character B(N), the second sequence.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer m
  integer n

  character a(m)
  character ai
  character b(n)
  character bi
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, max ( m, n )

    if ( i <= m ) then
      ai = a(i)
    else
      ai = ' '
    end if

    if ( i <= n ) then
      bi = b(i)
    else
      bi = ' '
    end if

    write ( *, '(i3,2x,a1,2x,a1)' ) i, ai, bi

  end do 

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
function i4_to_a ( i )

!*****************************************************************************80
!
!! I4_TO_A returns the I-th alphabetic character.
!
!  Example:
!
!    I  I4_TO_A
!
!   -8  ' '
!    0  ' '
!    1  'A'
!    2  'B'
!   ..
!   26  'Z'
!   27  'a'
!   52  'z'
!   53  ' '
!   99  ' '
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the index of the letter to be returned.
!    0 is a space;
!    1 through 26 requests 'A' through 'Z', (ASCII 65:90);
!    27 through 52 requests 'a' through 'z', (ASCII 97:122);
!
!    Output, character I4_TO_A, the requested alphabetic letter.
!
  implicit none

  integer, parameter :: cap_shift = 64
  integer i
  character i4_to_a
  integer, parameter :: low_shift = 96

  if ( i <= 0 ) then
    i4_to_a = ' '
  else if ( 1 <= i .and. i <= 26 ) then
    i4_to_a = char ( cap_shift + i )
  else if ( 27 <= i .and. i <= 52 ) then
    i4_to_a = char ( low_shift + i - 26 )
  else if ( i >= 53 ) then
    i4_to_a = ' '
  end if

  return
end
subroutine i4vec2_compare ( n, a1, a2, i, j, isgn )

!*****************************************************************************80
!
!! I4VEC2_COMPARE compares two I4VEC2's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data items.
!
!    Input, integer A1(N), A2(N), contain the two components of each item.
!
!    Input, integer I, J, the items to be compared.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item I > item J.
!
  implicit none

  integer n

  integer a1(n)
  integer a2(n)
  integer i
  integer isgn
  integer j

  isgn = 0

  if ( a1(i) < a1(j) ) then

    isgn = -1

  else if ( a1(i) == a1(j) ) then

    if ( a2(i) < a2(j) ) then
      isgn = -1
    else if ( a2(i) < a2(j) ) then
      isgn = 0
    else if ( a2(i) > a2(j) ) then
      isgn = +1
    end if

  else if ( a1(i) > a1(j) ) then

    isgn = +1

  end if

  return
end
subroutine i4vec2_print ( n, a, b, title )

!*****************************************************************************80
!
!! I4VEC2_PRINT prints a pair of I4VEC's, with an optional title.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), B(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer n

  integer a(n)
  integer b(n)
  integer i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) title
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,2i10)' ) i, a(i), b(i)
  end do

  return
end
subroutine i4vec2_sort_a ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC2_SORT_A ascending sorts an I4VEC2.
!
!  Discussion:
!
!    Each item to be sorted is a pair of integers (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items of data.
!
!    Input/output, integer A1(N), A2(N), the data to be sorted..
!
  implicit none

  integer n

  integer a1(n)
  integer a2(n)
  integer i
  integer indx
  integer isgn
  integer j
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

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( indx > 0 ) then

      call i4_swap ( a1(i), a1(j) )
      call i4_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4vec_reverse ( n, a )

!*****************************************************************************80
!
!! I4VEC_REVERSE reverses the elements of an I4VEC.
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N), the array to be reversed.
!
  implicit none

  integer n

  integer a(n)
  integer i

  do i = 1, n/2
    call i4_swap ( a(i), a(n+1-i) )
  end do

  return
end
subroutine pam120 ( acid_code, weight )

!*****************************************************************************80
!
!! PAM120 returns the PAM 120 substitution matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ACID_CODE(23), the nucleic acid codes.
!
!    Output, integer WEIGHT(23,23), the PAM120 substitution matrix.
!
  implicit none

  integer, parameter :: acid_num = 23

  character acid_code(acid_num)
  character, parameter, dimension ( acid_num ) :: acid_code2 = (/ &
    'A', 'R', 'N', 'B', 'D', 'C', 'Q', 'Z', 'E', 'G', &
    'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', &
    'Y', 'V', 'X' /)
  integer weight(acid_num,acid_num)
  integer, parameter, dimension ( acid_num, acid_num ) :: weight2 = &
     reshape ( (/&
     3,-3, 0, 0, 0,-3,-1,-1, 0, 1,-3,-1,-3,-2,-2,-4, 1, 1, 1,-7,-4, 0, 0, &
    -3, 6,-1,-2,-3,-4, 1,-1,-3,-4, 1,-2,-4, 2,-1,-4,-1,-1,-2, 1,-6,-3, 0, &
     0,-1, 4, 3, 2,-5, 0, 0, 1, 0, 2,-2,-4, 1,-3,-4,-2, 1, 0,-5,-2,-3, 0, &
     0,-2, 3, 4, 3,-6, 0, 1, 2, 0, 1,-3,-5, 0,-4,-6,-2, 0,-1,-7,-4,-3, 0, &
     0,-3, 2, 3, 5,-7, 1, 2, 3, 0, 0,-3,-5,-1,-4,-7,-2, 0,-1,-8,-5,-3, 0, &
    -3,-4,-5,-6,-7, 9,-7,-7,-7,-5,-4,-3,-7,-7,-6,-6,-3,-1,-3,-8,-1,-2, 0, &
    -1, 1, 0, 0, 1,-7, 6, 4, 2,-3, 3,-3,-2, 0,-1,-6, 0,-2,-2,-6,-5,-3, 0, &
    -1,-1, 0, 1, 2,-7, 4, 5, 3,-2, 1,-3,-3,-1,-3,-6,-1,-2,-2,-7,-5,-3, 0, &
     0,-3, 1, 2, 3,-7, 2, 3, 5,-1,-1,-3,-4,-1,-4,-6,-1,-1,-2,-8,-4,-3, 0, &
     1,-4, 0, 0, 0,-5,-3,-2,-1, 5,-4,-4,-5,-3,-4,-5,-2, 1,-1,-8,-6,-2, 0, &
    -3, 1, 2, 1, 0,-4, 3, 1,-1,-4, 7,-4,-3,-2,-4,-2,-1,-2,-3,-5,-1,-3, 0, &
    -1,-2,-2,-3,-3,-3,-3,-3,-3,-4,-4, 6, 1,-2, 1, 0,-3,-2, 0,-7,-2, 3, 0, &
    -3,-4,-4,-5,-5,-7,-2,-3,-4,-5,-3, 1, 5,-4, 3, 0,-3,-4,-3,-5,-3, 1, 0, &
    -2, 2, 1, 0,-1,-7, 0,-1,-1,-3,-2,-2,-4, 5, 0,-6,-2,-1,-1,-5,-6,-4, 0, &
    -2,-1,-3,-4,-4,-6,-1,-3,-4,-4,-4, 1, 3, 0, 8,-1,-3,-2,-1,-7,-4, 1, 0, &
    -4,-4,-4,-6,-7,-6,-6,-6,-6,-5,-2, 0, 0,-6,-1, 8,-5,-3,-4,-1, 4,-3, 0, &
     1,-1,-2,-2,-2,-3, 0,-1,-1,-2,-1,-3,-3,-2,-3,-5, 6, 1,-1,-7,-6,-2, 0, &
     1,-1, 1, 0, 0,-1,-2,-2,-1, 1,-2,-2,-4,-1,-2,-3, 1, 3, 2,-2,-3,-2, 0, &
     1,-2, 0,-1,-1,-3,-2,-2,-2,-1,-3, 0,-3,-1,-1,-4,-1, 2, 4,-6,-3, 0, 0, &
    -7, 1,-5,-7,-8,-8,-6,-7,-8,-8,-5,-7,-5,-5,-7,-1,-7,-2,-6,12,-1,-8, 0, &
    -4,-6,-2,-4,-5,-1,-5,-5,-4,-6,-1,-2,-3,-6,-4, 4,-6,-3,-3,-1, 8,-3, 0, &
     0,-3,-3,-3,-3,-2,-3,-3,-3,-2,-3, 3, 1,-4, 1,-3,-2,-2, 0,-8,-3, 5, 0, &
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  /), &
     (/ acid_num, acid_num /) )

  acid_code(1:acid_num) = acid_code2(1:acid_num)

  weight(1:acid_num,1:acid_num) = weight2(1:acid_num,1:acid_num)

  return
end
function pam120_score ( c1, c2 )

!*****************************************************************************80
!
!! PAM120_SCORE computes a single entry sequence/sequence matching score.
!
!  Discussion:
!
!    For each possible pair of characters C1 from sequence A and
!    C2 from sequence B, this routine returns a matching score,
!    which can be thought of as specifying how similar the two
!    characters are.  It is not necessary for the matching score
!    to be symmetric.  The character 'X' is commonly used to signify
!    a gap.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, two characters to be matched.
!    C1 is from sequence A, C2 from sequence B.
!
!    Output, real ( kind = 4 ) PAM120_SCORE, the score for matching 
!    the two characters.
!  
  implicit none

  integer, parameter :: acid_num = 23

  integer a_to_i4
  character acid_code(acid_num)
  integer, save, dimension ( 26 ) :: acid_index
  character c1
  character c2
  logical, parameter :: debug = .false.
  integer i
  integer j
  logical, save :: need_data = .true.
  real ( kind = 4 ) pam120_score
  integer, save, dimension ( acid_num, acid_num ) :: weight

  if ( need_data ) then

    call pam120 ( acid_code, weight )

    call a_index ( acid_num, acid_code, acid_index )

    need_data = .false.

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PAM120 substitution matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(3x,23(2x,a1))' ) acid_code(1:acid_num)
      write ( *, '(a)' ) ' '

      do i = 1, acid_num
        write ( *, '(a1,2x,23i3)' ) acid_code(i), weight(i,1:acid_num)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Acid index:'
      write ( *, '(a)' ) ' '
      write ( *, '(3x,26i3)' ) acid_index(1:26)

    end if

  end if

  pam120_score = 0.0E+00

  i = a_to_i4 ( c1 )
  j = a_to_i4 ( c2 )

  if ( i > 0 .and. j > 0 ) then

    i = acid_index ( i )
    j = acid_index ( j )

    if ( i > 0 .and. j > 0 ) then

      pam120_score = weight(i,j)

    end if

  end if

  return
end
subroutine pam200 ( acid_code, weight )

!*****************************************************************************80
!
!! PAM200 returns the PAM 200 substitution matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ACID_CODE(23), the nucleic acid codes.
!
!    Output, integer WEIGHT(23,23), the PAM200 substitution matrix.
!
  implicit none

  integer, parameter :: acid_num = 23

  character acid_code(acid_num)
  character, parameter, dimension ( acid_num ) :: acid_code2 = (/ &
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', &
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', &
    'B', 'Z', 'X' /)
  integer i
  integer j
  integer weight(acid_num,acid_num)
  integer, parameter, dimension ( acid_num, acid_num ) :: weight2 = &
    reshape ( (/ &
     2,-2, 0, 0,-2,-1, 0, 1,-1,-1,-2,-1,-1,-3, 1, 1, 1,-5,-3, 0, 1, 1, 0, &
    -2, 5, 0,-1,-3, 1,-1,-3, 1,-2,-2, 2, 0,-3, 0, 0,-1, 1,-4,-2, 0, 1, 0, &
     0, 0, 2, 2,-3, 0, 1, 0, 1,-1,-2, 1,-2,-2, 0, 1, 0,-4,-1,-2, 3, 2, 0, &
     0,-1, 2, 3,-4, 1, 3, 0, 0,-2,-3, 0,-2,-4,-1, 0, 0,-6,-3,-2, 4, 3, 0, &
    -2,-3,-3,-4, 8,-4,-4,-3,-3,-2,-5,-4,-4,-4,-2, 0,-2,-6, 0,-2,-3,-3, 0, &
    -1, 1, 0, 1,-4, 4, 2,-1, 2,-2,-1, 0,-1,-4, 0,-1,-1,-4,-3,-2, 2, 4, 0, &
     0,-1, 1, 3,-4, 2, 3, 0, 0,-2,-3, 0,-2,-4, 0, 0, 0,-6,-3,-2, 3, 4, 0, &
     1,-3, 0, 0,-3,-1, 0, 4,-2,-2,-3,-2,-3,-3,-1, 1, 0,-5,-4,-1, 1, 0, 0, &
    -1, 1, 1, 0,-3, 2, 0,-2, 5,-2,-2, 0,-2,-1, 0,-1,-1,-3, 0,-2, 2, 2, 0, &
    -1,-2,-1,-2,-2,-2,-2,-2,-2, 4, 2,-1, 2, 1,-2,-1, 0,-5,-1, 3,-1,-1, 0, &
    -2,-2,-2,-3,-5,-1,-3,-3,-2, 2, 4,-2, 3, 1,-2,-2,-1,-4,-1, 1,-2,-1, 0, &
    -1, 2, 1, 0,-4, 0, 0,-2, 0,-1,-2, 4, 1,-4,-1, 0, 0,-3,-4,-2, 1, 1, 0, &
    -1, 0,-2,-2,-4,-1,-2,-3,-2, 2, 3, 1, 5, 0,-2,-1, 0,-4,-2, 1,-1, 0, 0, &
    -3,-3,-2,-4,-4,-4,-4,-3,-1, 1, 1,-4, 0, 7,-4,-2,-2, 0, 5,-1,-2,-3, 0, &
     1, 0, 0,-1,-2, 0, 0,-1, 0,-2,-2,-1,-2,-4, 5, 1, 0,-5,-4,-1, 0, 1, 0, &
     1, 0, 1, 0, 0,-1, 0, 1,-1,-1,-2, 0,-1,-2, 1, 2, 1,-2,-2,-1, 1, 1, 0, &
     1,-1, 0, 0,-2,-1, 0, 0,-1, 0,-1, 0, 0,-2, 0, 1, 3,-4,-2, 0, 1, 0, 0, &
    -5, 1,-4,-6,-6,-4,-6,-5,-3,-5,-4,-3,-4, 0,-5,-2,-4,12, 0,-6,-3,-4, 0, &
    -3,-4,-1,-3, 0,-3,-3,-4, 0,-1,-1,-4,-2, 5,-4,-2,-2, 0, 7,-2,-1,-2, 0, &
     0,-2,-2,-2,-2,-2,-2,-1,-2, 3, 1,-2, 1,-1,-1,-1, 0,-6,-2, 4,-1,-1, 0, &
     1, 0, 3, 4,-3, 2, 3, 1, 2,-1,-2, 1,-1,-2, 0, 1, 1,-3,-1,-1, 4, 4, 0, &
     1, 1, 2, 3,-3, 4, 4, 0, 2,-1,-1, 1, 0,-3, 1, 1, 0,-4,-2,-1, 4, 5, 0, &
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  /), &
     (/ acid_num, acid_num /) )

  acid_code(1:acid_num) = acid_code2(1:acid_num)

  weight(1:acid_num,1:acid_num) = weight2(1:acid_num,1:acid_num)

  return
end
function pam200_score ( c1, c2 )

!*****************************************************************************80
!
!! PAM200_SCORE computes a single entry sequence/sequence matching score.
!
!  Discussion:
!
!    For each possible pair of characters C1 from sequence A and
!    C2 from sequence B, this routine returns a matching score,
!    which can be thought of as specifying how similar the two
!    characters are.  It is not necessary for the matching score
!    to be symmetric.  The character 'X' is commonly used to signify
!    a gap.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, two characters to be matched.
!    C1 is from sequence A, C2 from sequence B.
!
!    Output, real ( kind = 4 ) PAM200_SCORE, the score for matching the two characters.
!
  implicit none

  integer, parameter :: acid_num = 23

  integer a_to_i4
  character acid_code(acid_num)
  integer, save, dimension ( 26 ) :: acid_index
  character c1
  character c2
  logical, parameter :: debug = .false.
  integer i
  integer j
  logical, save :: need_data = .true.
  real ( kind = 4 ) pam200_score
  integer, save, dimension ( acid_num, acid_num ) :: weight

  if ( need_data ) then

    call pam200 ( acid_code, weight )
    call a_index ( acid_num, acid_code, acid_index )
    need_data = .false.

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PAM200 substitution matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(3x,23(2x,a1))' ) ( acid_code(j), j = 1, acid_num )
      write ( *, '(a)' ) ' '

      do i = 1, acid_num
        write ( *, '(a1,2x,23i3)' ) acid_code(i), weight(i,1:acid_num)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Acid index:'
      write ( *, '(a)' ) ' '
      write ( *, '(3x,26i3)' ) acid_index(1:26)

    end if

  end if

  pam200_score = 0.0E+00

  i = a_to_i4 ( c1 )
  j = a_to_i4 ( c2 )

  if ( i > 0 .and. j > 0 ) then

    i = acid_index ( i )
    j = acid_index ( j )

    if ( i > 0 .and. j > 0 ) then

      pam200_score = weight(i,j)

    end if

  end if

  return
end
subroutine r4vec2_sum_imax ( n, a, b, imax )

!*****************************************************************************80
!
!! R4VEC2_SUM_IMAX returns the index of the maximum sum of two R4VEC's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real ( kind = 4 ) A(N), B(N), two arrays whose sum is to be examined.
!
!    Output, integer IMAX, the index of the largest entry in A+B.
!
  implicit none

  integer n

  real ( kind = 4 ) a(n)
  real ( kind = 4 ) b(n)
  integer i
  integer imax
  real ( kind = 4 ) sum_max

  if ( n <= 0 ) then

    imax = 0

  else

    imax = 1
    sum_max = a(1) + b(1)

    do i = 2, n
      if ( a(i) + b(i) > sum_max ) then
        sum_max = a(i) + b(i)
        imax = i
      end if
    end do

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
subroutine s_to_chvec ( s, n, cvec )

!*****************************************************************************80
!
!! S_TO_CHVEC converts a string to a character vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string of characters.
!
!    Input/output, integer N.
!    if N is -1, extract characters from 1 to len(S);
!    if N is 0, extract characters up to the last nonblank;
!    if N is positive, extract characters from 1 to N.
!
!    On output, N is the number of characters successfully extracted.
!
!    Output, character CVEC(N), the characters extracted from S.
!
  implicit none

  character cvec(*)
  integer i
  integer n
  character ( len = * ) s

  if ( n <= - 1 ) then
    n = len ( s )
  else if ( n == 0 ) then
    n = len_trim ( s )
  else
    n = min ( n, len ( s ) )
  end if

  do i = 1, n
    cvec(i) = s(i:i)
  end do
 
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
!    Output, integer IVAL, the integer value read from the string.
!    If STRING is blank, then IVAL will be returned 0.
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
function simple_score ( c1, c2 )

!*****************************************************************************80
!
!! SIMPLE_SCORE computes a single entry sequence/sequence matching score.
!
!  Discussion:
!
!    This routine is a sample scoring function which returns a score of
!    -2 if either character is '-', representing a gap,
!    -1 if neither character is '-', but the characters are not identical,
!     0 if the characters are identical, and not '-'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, two characters to be matched.
!    C1 is from sequence A, C2 from sequence B.
!
!    Output, real ( kind = 4 ) SIMPLE_SCORE, the score for matching the two characters.
!
  implicit none

  character c1
  character c2
  real ( kind = 4 ) score
  real ( kind = 4 ) simple_score

  if ( c1 == '-' .or. c2 == '-' ) then
    score = -2.0E+00
  else if ( c1 /= c2 ) then
    score = -1.0E+00
  else
    score = 0.0E+00
  end if

  simple_score = score

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into linear order.
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
!    19 May 1999
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of items to be sorted.
!
!    Input/output, integer INDX, the main communication signal.
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
!      * set ISGN = -1 if I < J, ISGN = +1 if I > J;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ISGN, results of comparison of elements I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    ISGN => 0 means I is greater than or equal to J.
!
  implicit none

  integer i
  integer indx
  integer isgn
  integer j
  integer, save :: k = 0
  integer, save :: k1 = 0
  integer n
  integer, save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    n1 = n
    k = n / 2
    k1 = k
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i = i + 1
      end if

      j = k1
      k1 = i
      indx = - 1
      return

    end if

    if ( isgn > 0 ) then
      indx = 2
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        indx = 0
      else
        i = n1
        n1 = n1 - 1
        j = 1
        indx = 1
      end if

      return

    end if

    k = k - 1
    k1 = k
!
!  INDX > 0, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i = 2 * k1

    if ( i == n1 ) then
      j = k1
      k1 = i
      indx = - 1
      return
    else if ( i <= n1 ) then
      j = i + 1
      indx = - 2
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    indx = 0
  else
    i = n1
    n1 = n1 - 1
    j = 1
    indx = 1
  end if

  return
end
subroutine ss_gd_bpq ( m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

!*****************************************************************************80
!
!! SS_GD_BPQ does a global distance backward alignment path in quadratic space.
!
!  Discussion:
!
!    An effort has been made to handle the ambiguous case, where
!    more than one optimal path enters a cell.  In such a case,
!    the code tries to take a delete out if there was a delete in,
!    or an insert out if there was an insert in, since the optimal
!    score calculation includes a penalty for gap opening.
!
!    The routine is called "quadratic" because it uses an M by N array
!    to do the alignment.
!
!    The score table must have been computed before this routine is called.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!  Parameters:
!
!    Input, integer M, the number of entries in sequence A.
!
!    Input, integer M1, M2, the minimum and maximum rows of the computed score 
!    table.
!
!    Input, integer N, the number of entries in sequence B.
!
!    Input, integer N1, N2, the minimum and maximum columns of the computed
!    score table.
!
!    Input, integer LDS, the declared upper first dimension of S.
!    LDS must be at least M.
!
!    Input, integer T(0:LDS,0:N), the backward pointer table.
!
!    Output, integer NPATH, the number of points in the matching.
!
!    Input, integer PATHI(M+N+1), PATHJ(M+N+1), contains, in the first
!    NPATH entries, the indices of the aligned items.
!    The first entries are special marker values: 
!      PATHI(1) = 0 and PATHJ(1) = 0; 
!    A value of -1 for PATHI or PATHJ indicates a null partner.
!    Otherwise, A(PATHI(I)) is matched to B(PATHJ(I)).
!
  implicit none

  integer lds
  integer m
  integer n

  integer d_new
  integer d_old
  integer i
  integer i_new
  integer i_old
  integer ipath
  integer j
  integer j_new
  integer j_old
  integer m_new
  integer m_old
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  integer t(0:lds,0:n)
  integer tij
  integer tij_old

  npath = 0

  i = m1
  j = n1

  i_old = 0
  d_old = 0
  m_old = 1

  do while ( i <= m2 .and. j <= n2 )

    npath = npath + 1
    pathi(npath) = i
    pathj(npath) = j

    if ( i == m2 .and. j == n2 ) then
      exit
    end if

    tij = t(i,j)

    m_new = mod ( tij, 2 )
    tij = tij / 2
    i_new = mod ( tij, 2 )
    tij = tij / 2
    d_new = tij
!
!  Try to handle ambiguous cases.
!
    if ( i_old == 1 ) then

      if ( i_new == 1 ) then
        d_new = 0
        m_new = 0
      end if

    else if ( d_old == 1 ) then

      if ( d_new == 1 ) then
        i_new = 0
        m_new = 0
      end if

    end if

    if ( j < n2 .and. i_new == 1 ) then

      j = j + 1
      d_new = 0
      m_new = 0

    else if ( i < m2 .and. d_new == 1 ) then

      i = i + 1
      i_new = 0
      m_new = 0

    else if ( i < m2 .and. j < n2 .and. m_new == 1 ) then

      i = i + 1
      j = j + 1
      i_new = 0
      d_new = 0

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SS_GD_BPQ: Unexpected situation!'
      write ( *, '(a,2i6)' ) '  I, J = ', i, j
      write ( *, '(a,i6)' ) '  T(I,J) = ', t(i,j)
      stop

    end if
!
!  Copy the information.  Only one of these three values is now nonzero,
!  recording which direction we actually took.
!
    i_old = i_new
    d_old = d_new
    m_old = m_new

  end do
!
!  Mark DELETEs and INSERTs.
!
  i_new = -1
  j_new = -1

  do ipath = 1, npath

    i_old = i_new
    j_old = j_new

    i_new = pathi(ipath)
    j_new = pathj(ipath)

    if ( i_new == i_old ) then
      pathi(ipath) = -1
    else if ( j_new == j_old ) then
      pathj(ipath) = -1
    end if

  end do

  return
end
subroutine ss_gd_bsl ( a, b, m, m1, m2, n, n1, n2, ss_score, gap_extend, &
  s, e, f, t )

!*****************************************************************************80
!
!! SS_GD_BSL does a global distance backward alignment score in linear space.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!  Parameters:
!
!    Input, character A(M), B(N), two sequences to be aligned.
!
!    Input, integer M, the number of entries in sequence A.
!
!    Input, integer M1, M2, the minimum and maximum rows of the score 
!    matrix to compute.  0 <= M1 <= M2 <= M.
!
!    Input, integer N, the number of entries in sequence B.
!
!    Input, integer N1, N2, the minimum and maximum columns of the score 
!    matrix to compute.  0 <= N1 <= N2 <= N.
!
!    Input, external SS_SCORE, the name of a function of the form
!      function ss_score ( c1, c2 )
!    which returns a real value SS_SCORE for the matching of the characters
!    C1 from sequence A to C2 from sequence B.
!
!    Input, real ( kind = 4 ) GAP_EXTEND, the penalty for extending a gap.  A gap of 
!    length 7 incurs a penalty 7 * GAP_EXTEND.
!
!    Output, real ( kind = 4 ) S(0:N), the backward optimal score vector.
!    The maximum possible alignment score is in S(N1).
!
!    Output, real ( kind = 4 ) E(0:N), the backward final insertion score vector.
!
!    Output, real ( kind = 4 ) F(0:N), the backward final deletion score vector.
!
!    Output, integer T(0:N), the backward pointer vector.
!
  implicit none

  integer, parameter :: DUNNO = 0
  integer, parameter :: MATCH = 1
  integer, parameter :: INSERT = 2
  integer, parameter :: DELETE = 4

  integer m
  integer n

  character a(m)
  character b(n)
  real ( kind = 4 ) diag_new
  real ( kind = 4 ) diag_old
  real ( kind = 4 ) e(0:n)
  real ( kind = 4 ) f(0:n)
  real ( kind = 4 ) gap_extend
  integer i
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real ( kind = 4 ) s(0:n)
  real ( kind = 4 ), external :: ss_score
  integer t(0:n)
!
!  The last row, I = M2.
!
  e(n2) = 0.0
  f(n2) = 0.0
  s(n2) = 0.0
  t(n2) = DUNNO

  do j = n2-1, n1, -1
    e(j) = e(j+1) + gap_extend
    f(j) = f(j+1) + gap_extend
    s(j) = s(j+1) + gap_extend
    t(j) = INSERT
  end do
!
!  Upper rectangle.
!
  do i = m2-1, m1, -1

    diag_old = s(n2)

    e(n2) = e(n2) + gap_extend
    f(n2) = f(n2) + gap_extend
    s(n2) = s(n2) + gap_extend
    t(n2) = DELETE

    do j = n2-1, n1, -1
!
!  Insertion.
!
      e(j) = e(j+1) + gap_extend
      if ( s(j+1) > e(j+1) ) then
        e(j) = s(j+1) + gap_extend
      end if
!
!  Deletion.
!
      f(j) = f(j) + gap_extend
      if ( s(j) + gap_extend > f(j) ) then
        f(j) = s(j) + gap_extend
      end if
!
!  Best.
!
      diag_new = s(j)

      s(j) = diag_old + ss_score ( a(i+1), b(j+1) )
      t(j) = MATCH

      if ( e(j) == s(j) ) then
        t(j) = t(j) + INSERT
      else if ( e(j) > s(j) ) then
        s(j) = e(j)
        t(j) = INSERT
      end if

      if ( f(j) == s(j) ) then
        t(j) = t(j) + DELETE
      else if ( f(j) > s(j) ) then
        s(j) = f(j)
        t(j) = DELETE
      end if

      diag_old = diag_new

    end do

  end do

  return
end
subroutine ss_gd_bsq ( a, b, m, m1, m2, n, n1, n2, ss_score, gap_extend, &
  lds, s, e, f, t )

!*****************************************************************************80
!
!! SS_GD_BSQ does a global distance backward alignment score in quadratic space.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!  Parameters:
!
!    Input, character A(M), B(N), two sequences to be aligned.
!
!    Input, integer M, the number of entries in the sequence A.
!
!    Input, integer M1, M2, the minimum and maximum rows of the score 
!    matrix to compute.  0 <= M1 <= M2 <= M.
!
!    Input, integer N, the number of entries in the sequence B.
!
!    Input, integer N1, N2, the minimum and maximum columns of the score 
!    matrix to compute.  0 <= N1 <= N2 <= N.
!
!    Input, external SS_SCORE, the name of a function of the form
!      function ss_score ( c1, c2 )
!    which returns a real value SS_SCORE for the matching of the characters
!    C1 from sequence A to C2 from sequence B.
!
!    Input, real ( kind = 4 ) GAP_EXTEND, the penalty for extending a gap.  A gap of 
!    length 7 incurs a penalty 7 * GAP_EXTEND.
!
!    Input, integer LDS, the declared upper first dimension of S.
!    LDS must be at least M.
!
!    Output, real ( kind = 4 ) S(0:LDS,0:N), the M1:M by 0:N reverse score table.  
!    The maximum possible alignment score is in S(M1,0).
!
!    Output, real ( kind = 4 ) E(0:LDS,0:N), the backward final insertion table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
!    Output, real ( kind = 4 ) F(0:LDS,0:N), the backward final deletion table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
!    Output, integer T(0:LDS,0:N), the backward pointer table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
  implicit none

  integer, parameter :: DUNNO = 0
  integer, parameter :: MATCH = 1
  integer, parameter :: INSERT = 2
  integer, parameter :: DELETE = 4

  integer lds
  integer m
  integer n

  character a(m)
  character b(n)
  real ( kind = 4 ) e(0:lds,0:n)
  real ( kind = 4 ) f(0:lds,0:n)
  real ( kind = 4 ) gap_extend
  integer i
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real ( kind = 4 ) s(0:lds,0:n)
  real ( kind = 4 ), external :: ss_score
  integer t(0:lds,0:n)
!
!  Lower Right corner.
!
  e(m2,n2) = 0.0
  f(m2,n2) = 0.0
  s(m2,n2) = 0.0
  t(m2,n2) = DUNNO

  do j = n2-1, n1, -1
    e(m2,j) = e(m2,j+1) + gap_extend
    f(m2,j) = f(m2,j+1) + gap_extend
    s(m2,j) = s(m2,j+1) + gap_extend
    t(m2,j) = INSERT
  end do
!
!  Upper rectangle.
!
  do i = m2-1, m1, -1

    e(i,n2) = e(i+1,n2) + gap_extend
    f(i,n2) = f(i+1,n2) + gap_extend
    s(i,n2) = s(i+1,n2) + gap_extend
    t(i,n2) = DELETE

    do j = n2-1, n1, -1
!
!  Insertion.
!
      e(i,j) = e(i,j+1) + gap_extend

      if ( s(i,j+1) + gap_extend > e(i,j) ) then
        e(i,j) = s(i,j+1) + gap_extend
      end if
!
!  Deletion.
!
      f(i,j) = f(i+1,j) + gap_extend

      if ( s(i+1,j) + gap_extend > f(i,j) ) then
        f(i,j) = s(i+1,j) + gap_extend
      end if
!
!  Best.
!
      s(i,j) = s(i+1,j+1) + ss_score ( a(i+1), b(j+1) )
      t(i,j) = MATCH

      if ( e(i,j) == s(i,j) ) then
        t(i,j) = t(i,j) + INSERT
      else if ( e(i,j) > s(i,j) ) then
        s(i,j) = e(i,j)
        t(i,j) = INSERT
      end if

       if ( f(i,j) == s(i,j) ) then
        t(i,j) = t(i,j) + DELETE
       else if ( f(i,j) > s(i,j) ) then
        s(i,j) = f(i,j)
        t(i,j) = DELETE
      end if

    end do
  end do

  return
end
subroutine ss_gd_fpq ( m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

!*****************************************************************************80
!
!! SS_GD_FPQ does a global distance forward alignment path in quadratic space.
!
!  Discussion:
!
!    An effort has been made to handle the ambiguous case, where
!    more than one optimal path enters a cell.  In such a case,
!    the code tries to take a delete out if there was a delete in,
!    or an insert out if there was an insert in, since the optimal
!    score calculation includes a penalty for gap opening.
!
!    The routine is called "quadratic" because it uses an M by N array
!    to do the alignment.
!
!    The score table must have been computed before this routine is called.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!  Parameters:
!
!    Input, integer M, the number of entries in sequence A.
!
!    Input, integer M1, M2, the minimum and maximum rows of the computed score 
!    table.
!
!    Input, integer N, the number of entries in sequence B.
!
!    Input, integer N1, N2, the minimum and maximum columns of the computed
!    score table.
!
!    Input, integer LDS, the declared upper first dimension of S.
!    LDS must be at least M.
!
!    Input, integer T(0:LDS,0:N), the forward pointer table.
!
!    Output, integer NPATH, the number of points in the matching.
!
!    Input, integer PATHI(M+N+1), PATHJ(M+N+1), contains, in the first
!    NPATH entries, the indices of the aligned items.
!    The first entries are special marker values: 
!      PATHI(1) = 0 and PATHJ(1) = 0; 
!    A value of -1 for PATHI or PATHJ indicates a null partner.
!    Otherwise, A(PATHI(I)) is matched to B(PATHJ(I)).
!
  implicit none

  integer lds
  integer m
  integer n

  integer d_new
  integer d_old
  integer i
  integer i_new
  integer i_old
  integer ipath
  integer j
  integer j_new
  integer j_old
  integer m_new
  integer m_old
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  integer t(0:lds,0:n)
  integer tij
  integer tij_old

  npath = 0

  i = m2
  j = n2

  i_old = 0
  d_old = 0
  m_old = 1

  do while ( i >= m1 .and. j >= n1 )

    npath = npath + 1
    pathi(npath) = i
    pathj(npath) = j

    if ( i == m1 .and. j == n1 ) then
      exit
    end if

    tij = t(i,j)

    m_new = mod ( tij, 2 )
    tij = tij / 2
    i_new = mod ( tij, 2 )
    tij = tij / 2
    d_new = tij
!
!  Try to handle ambiguous cases.
!
    if ( i_old == 1 ) then

      if ( i_new == 1 ) then
        d_new = 0
        m_new = 0
      end if

    else if ( d_old == 1 ) then

      if ( d_new == 1 ) then
        i_new = 0
        m_new = 0
      end if

    end if

    if ( j > n1 .and. i_new == 1 ) then

      j = j - 1
      d_new = 0
      m_new = 0

    else if ( i > m1 .and. d_new == 1 ) then

      i = i - 1
      i_new = 0
      m_new = 0

    else if ( i > m1 .and. j > n1 .and. m_new == 1 ) then

      i = i - 1
      j = j - 1
      i_new = 0
      d_new = 0

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SS_GG_FPQ: Unexpected situation!'
      write ( *, '(a,2i6)' ) '  I, J = ', i, j
      write ( *, '(a,i6)' ) '  T(I,J) = ', t(i,j)
      stop

    end if
!
!  Copy the information.  Only one of these three values is now nonzero,
!  recording which direction we actually took.
!
    i_old = i_new
    d_old = d_new
    m_old = m_new

  end do
!
!  Put the path into proper order.
!
  call i4vec_reverse ( npath, pathi )
  call i4vec_reverse ( npath, pathj )
!
!  Mark DELETEs and INSERTs.
!
  i_new = -1
  j_new = -1

  do ipath = 1, npath

    i_old = i_new
    j_old = j_new

    i_new = pathi(ipath)
    j_new = pathj(ipath)

    if ( i_new == i_old ) then
      pathi(ipath) = -1
    else if ( j_new == j_old ) then
      pathj(ipath) = -1
    end if

  end do

  return
end
subroutine ss_gd_fsl ( a, b, m, m1, m2, n, n1, n2, ss_score, gap_extend, &
  s, e, f, t )

!*****************************************************************************80
!
!! SS_GD_FSL does a global distance forward alignment score in linear space.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!  Parameters:
!
!    Input, character A(M), B(N), two sequences to be aligned.
!
!    Input, integer M, the number of entries in sequence A.
!
!    Input, integer M1, M2, the minimum and maximum rows of the score 
!    matrix to compute.  0 <= M1 <= M2 <= M.
!
!    Input, integer N, the number of entries in sequence B.
!
!    Input, integer N1, N2, the minimum and maximum columns of the score 
!    matrix to compute.  0 <= N1 <= N2 <= N.
!
!    Input, external SS_SCORE, the name of a function of the form
!      function ss_score ( c1, c2 )
!    which returns a real value SS_SCORE for the matching of the characters
!    C1 from sequence A to C2 from sequence B.
!
!    Input, real ( kind = 4 ) GAP_EXTEND, the penalty for extending a gap.  A gap of 
!    length 7 incurs a penalty 7 * GAP_EXTEND.
!
!    Output, real ( kind = 4 ) S(0:N), the forward optimal score vector.  
!    The maximum possible alignment score is in S(N2).
!
!    Output, real ( kind = 4 ) E(0:N), the forward final insertion score vector.
!
!    Output, real ( kind = 4 ) F(0:N), the forward final deletion score vector.
!
!    Output, integer T(0:N), the forward pointer vector.
!
  implicit none

  integer, parameter :: DUNNO = 0
  integer, parameter :: MATCH = 1
  integer, parameter :: INSERT = 2
  integer, parameter :: DELETE = 4

  integer m
  integer n

  character a(m)
  character b(n)
  real ( kind = 4 ) diag_new
  real ( kind = 4 ) diag_old
  real ( kind = 4 ) e(0:n)
  real ( kind = 4 ) f(0:n)
  real ( kind = 4 ) gap_extend
  integer i
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real ( kind = 4 ) s(0:n)
  real ( kind = 4 ), external :: ss_score
  integer t(0:n)
!
!  The first row, I = M1.
!
  e(n1) = 0.0
  f(n1) = 0.0
  s(n1) = 0.0
  t(n1) = DUNNO

  do j = n1+1, n2
    e(j) = e(j-1) + gap_extend
    f(j) = f(j-1) + gap_extend
    s(j) = s(j-1) + gap_extend
    t(j) = INSERT
  end do
!
!  Subsequent rows.
!
  do i = m1+1, m2

    diag_old = s(n1)

    e(n1) = e(n1) + gap_extend
    f(n1) = f(n1) + gap_extend
    s(n1) = s(n1) + gap_extend
    t(n1) = DELETE

    do j = n1+1, n2
!
!  Insertion
!
      e(j) = e(j-1) + gap_extend

      if ( s(j-1) + gap_extend > e(j) ) then
        e(j) = s(j-1) + gap_extend
      end if
!
!  Deletion.
!
      f(j) = f(j) + gap_extend

      if ( s(j) + gap_extend > f(j) ) then
        f(j) = s(j) + gap_extend
      end if
!
!  Best.
!
      diag_new = s(j)

      s(j) = diag_old + ss_score ( a(i), b(j) )
      t(j) = MATCH

      if ( e(j) == s(j) ) then
        t(j) = t(j) + INSERT
      else if ( e(j) > s(j) ) then
        s(j) = e(j)
        t(j) = INSERT
      end if

      if ( f(j) == s(j) ) then
        t(j) = t(j) + DELETE
      else if ( f(j) > s(j) ) then
        s(j) = f(j)
        t(j) = DELETE
      end if

      diag_old = diag_new

    end do
  end do

  return
end
subroutine ss_gd_fsq ( a, b, m, m1, m2, n, n1, n2, ss_score, gap_extend, &
  lds, s, e, f, t )

!*****************************************************************************80
!
!! SS_GD_FSQ does a global distance forward alignment score in quadratic space.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!  Parameters:
!
!    Input, character A(M), B(N), two sequences to be aligned.
!
!    Input, integer M, the number of entries in sequence A.
!
!    Input, integer M1, M2, the minimum and maximum rows of the score 
!    matrix to compute.  0 <= M1 <= M2 <= M.
!
!    Input, integer N, the number of entries in sequence B.
!
!    Input, integer N1, N2, the minimum and maximum columns of the score 
!    matrix to compute.  0 <= N1 <= N2 <= N.
!
!    Input, external SS_SCORE, the name of a function of the form
!      function ss_score ( c1, c2 )
!    which returns a real value SS_SCORE for the matching of the characters
!    C1 from sequence A to C2 from sequence B.
!
!    Input, real ( kind = 4 ) GAP_EXTEND, the penalty for extending a gap.  A gap of 
!    length 7 incurs a penalty 7 * GAP_EXTEND.
!
!    Input, integer LDS, the declared upper first dimension of S.
!    LDS must be at least M.
!
!    Output, real ( kind = 4 ) S(0:LDS,0:N), the forward optimal score table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
!    Output, real ( kind = 4 ) E(0:LDS,0:N), the forward final insertion table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
!    Output, real ( kind = 4 ) F(0:LDS,0:N), the forward final deletion table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
!    Output, integer T(0:LDS,0:N), the forward pointer table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
  implicit none

  integer, parameter :: DUNNO = 0
  integer, parameter :: MATCH = 1
  integer, parameter :: INSERT = 2
  integer, parameter :: DELETE = 4

  integer lds
  integer m
  integer n

  character a(m)
  character b(n)
  real ( kind = 4 ) e(0:lds,0:n)
  real ( kind = 4 ) f(0:lds,0:n)
  real ( kind = 4 ) gap_extend
  integer i
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real ( kind = 4 ) s(0:lds,0:n)
  real ( kind = 4 ), external :: ss_score
  integer t(0:lds,0:n)
!
!  Upper Left corner.
!
  e(m1,n1) = 0.0
  f(m1,n1) = 0.0
  s(m1,n1) = 0.0
  t(m1,n1) = DUNNO
!
!  Upper Right row.
!
  do j = n1+1, n2
    e(m1,j) = e(m1,j-1) + gap_extend
    f(m1,j) = f(m1,j-1) + gap_extend
    s(m1,j) = s(m1,j-1) + gap_extend
    t(m1,j) = INSERT
  end do
!
!  Lower rectangle.
!
  do i = m1+1, m2

    e(i,n1) = e(i-1,n1) + gap_extend
    f(i,n1) = f(i-1,n1) + gap_extend
    s(i,n1) = s(i-1,n1) + gap_extend
    t(i,n1) = DELETE

    do j = n1+1, n2
!
!  Insertion.
!
      e(i,j) = e(i,j-1) + gap_extend

      if ( s(i,j-1) + gap_extend > e(i,j) ) then
        e(i,j) = s(i,j-1) + gap_extend
      end if
!
!  Deletion.
!
      f(i,j) = f(i-1,j) + gap_extend

      if ( s(i-1,j) + gap_extend > f(i,j) ) then
        f(i,j) = s(i-1,j) + gap_extend
      end if
!
!  Best.
!
      s(i,j) = s(i-1,j-1) + ss_score ( a(i), b(j) )
      t(i,j) = MATCH

      if ( e(i,j) == s(i,j) ) then
        t(i,j) = t(i,j) + INSERT
      else if ( e(i,j) > s(i,j) ) then
        s(i,j) = e(i,j)
        t(i,j) = INSERT
      end if

      if ( f(i,j) == s(i,j) ) then
        t(i,j) = t(i,j) + DELETE
      else if ( f(i,j) > s(i,j) ) then
        s(i,j) = f(i,j)
        t(i,j) = DELETE
      end if

    end do
  end do

  return
end
subroutine ss_gd_match_print ( a, b, m, n, npath, pathi, pathj, ss_score, &
  gap_extend )

!*****************************************************************************80
!
!! SS_GD_MATCH_PRINT prints a global distance alignment.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!  Parameters:
!
!    Input, character A(M), B(N), two sequences which have been aligned.
!
!    Input, integer M, the number of entries in sequence A.
!
!    Input, integer N, the number of entries in sequence B.
!
!    Input, integer NPATH, the number of points in the alignment.
!
!    Input, integer PATHI(M+N+1), PATHJ(M+N+1), contains, in the first
!    NPATH entries, the indices of the most recently matched items
!    from sequences A and B.  PATHI(1) = PATHJ(1) = 0.  PATHI(2)
!    and PATHJ(2) are the first elements of A and B to be aligned.
!
!    Input, external SS_SCORE, the name of a function of the form
!      function ss_score ( c1, c2 )
!    which returns a real value SS_SCORE for the matching of the character
!    C1 from sequence A to the character C2 from sequence B.
!
!    Input, real ( kind = 4 ) GAP_EXTEND, the penalty for extending a gap.  A gap of 
!    length 7 incurs a penalty 7 * GAP_EXTEND.
!
  implicit none

  integer m
  integer n

  character a(m)
  character b(n)
  character ( len = 3 ) c1
  character c2
  character c3
  character c4
  character ( len = 3 ) c5
  real ( kind = 4 ) gap_extend
  integer i
  real inc
  integer ipath
  integer j
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real ( kind = 4 ), external :: ss_score
  real ( kind = 4 ) sum2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Sequence/sequence matching,'
  write ( *, '(a)' ) 'Affine gap penalty:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' #  A     B    #      Increm     Score'
  write ( *, '(a)' ) ' '

  sum2 = 0.0
  ipath = 1

  write ( c1, '(i3)' ) 0
  c2 = ' '
  c3 = ' '
  c4 = ' '
  write ( c5, '(i3)' ) 0

  write ( *, '(a3,2x,a1,2x,a1,2x,a1,2x,a3,2x,10x,f10.2)' ) &
    c1, c2, c3, c4, c5, sum2

  do ipath = 2, npath

    i = pathi(ipath)
    j = pathj(ipath)

    if ( i == -1 ) then
      c1 = ' '
      c2 = '|'
      c3 = ' '
      c4 = b(j)
      write ( c5, '(i3)' ) j
      inc = + gap_extend
    else if ( j == -1 ) then
      write ( c1, '(i3)' ) i
      c2 = a(i)
      c3 = ' '
      c4 = '|'
      c5 = ' '
      inc = gap_extend
    else
      write ( c1, '(i3)' ) i
      c2 = a(i)
      if ( a(i) == b(j) ) then
        c3 = '='
      else
        c3 = '-'
      end if
      c4 = b(j)
      write ( c5, '(i3)' ) j
      inc = ss_score ( a(i), b(j) )
    end if

    sum2 = sum2 + inc

    write ( *, '(a3,2x,a1,2x,a1,2x,a1,2x,a3,2x,2f10.2)' ) &
      c1, c2, c3, c4, c5, inc, sum2

  end do

  return
end
subroutine ss_gd_rpl ( a, b, m, n, ss_score, gap_extend, npath, pathi, pathj )

!*****************************************************************************80
!
!! SS_GD_RPL does a global distance recursive alignment path in linear space.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!  Parameters:
!
!    Input, character A(M), B(N), two sequences to be aligned.
!
!    Input, integer M, the number of entries in sequence A.
!
!    Input, integer N, the number of entries in sequence B.
!
!    Input, external SS_SCORE, the name of a function of the form
!      function ss_score ( c1, c2 )
!    which returns a real value SS_SCORE for the matching of the characters
!    C1 from sequence A to C2 from sequence B.
!
!    Input, real ( kind = 4 ) GAP_EXTEND, the penalty for extending a gap.  A gap of 
!    length 7 incurs a penalty 7 * GAP_EXTEND.
!
!    Output, integer NPATH, the number of points in the matching.
!
!    Output, integer PATHI(M+N+1), PATHJ(M+N+1), contains, in the first
!    NPATH entries, the indices of the matched items from sequences A 
!    and B.
!
  implicit none

  integer m
  integer n

  character a(m)
  character b(n)
  real ( kind = 4 ) eb(0:n)
  real ( kind = 4 ) ef(0:n)
  real ( kind = 4 ) fb(0:n)
  real ( kind = 4 ) ff(0:n)
  real ( kind = 4 ) gap_extend
  integer i
  integer i1
  integer i2
  integer i3
  integer j
  integer j1
  integer j2
  integer j3
  integer maxstack
  integer npath
  integer nstack
  integer p1
  integer p2
  integer p3
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real ( kind = 4 ) sb(0:n)
  real ( kind = 4 ) sf(0:n)
  real ( kind = 4 ), external :: ss_score
  integer stack(2*(m+n),2)
  integer tb(0:n)
  integer tf(0:n)

  maxstack = 2 * ( m + n )
!
!  We know already where the path begins and ends.
!
  npath = 1
  pathi(npath) = 0
  pathj(npath) = 0

  npath = npath + 1
  pathi(npath) = m
  pathj(npath) = n
!
!  Each entry of the stack stores a pair of indices of PATH that
!  are endpoints of an unfinished interval.
!
  nstack = 1
  stack(nstack,1) = 1
  stack(nstack,2) = 2
!
!  If the stack has any jobs in it, work on one.
!
  do while ( nstack > 0 )

    p1 = stack(nstack,1)
    p3 = stack(nstack,2)
    nstack = nstack - 1

    i1 = pathi(p1)
    j1 = pathj(p1)

    i3 = pathi(p3)
    j3 = pathj(p3)
!
!  I2 is the row midway between I1 and I3.
!
    i2 = ( i1 + i3 ) / 2
!
!  I2 = I1, skip.
!
    if ( i2 == i1 ) then
!
!  I2 = I3, skip.
!
    else if ( i2 == i3 ) then
!
!  J1 = J3, skip.
!
    else if ( j1 == j3 ) then

      do i = i1 + 1, i2 - 1
        npath = npath + 1
        pathi(npath) = i
        pathj(npath) = j1
      end do
!
!  I1 < I2 < I3 and J1 < J3.
!  Compute J2, the optimal match in row I2.
!
    else

      call ss_gd_fsl ( a, b, m, i1, i2, n, j1, j3, ss_score, gap_extend, &
        sf, ef, ff, tf )

      call ss_gd_bsl ( a, b, m, i2, i3, n, j1, j3, ss_score, gap_extend, &
        sb, eb, fb, tb )

      call r4vec2_sum_imax ( j3+1-j1, sf(j1), sb(j1), j2 )

      j2 = j1 + j2 - 1

      npath = npath + 1
      pathi(npath) = i2
      pathj(npath) = j2

      p2 = npath
!
!  CONSIDER THE BLOCK (I1,J1) x (I2,J2).
!
!  ...a horizontal move.
!
      if ( i2 == i1 ) then

        do j = j1 + 1, j2 - 1
          npath = npath + 1
          pathi(npath) = i1
          pathj(npath) = j
        end do
!
!  ...a vertical move.
!
      else if ( j2 == j1 ) then

        do i = i1 + 1, i2 - 1
          npath = npath + 1
          pathi(npath) = i
          pathj(npath) = j1
        end do
!
!  ...a diagonal move.
!
      else if ( i2 == i1 + 1 .and. j2 == j1 + 1 ) then
!
!  ...a slide right and then down one.
!
      else if ( i2 == i1 + 1 .and. j2 > j1 + 1 ) then

        do j = j1 + 1, j2 - 1
          npath = npath + 1
          pathi(npath) = i1
          pathj(npath) = j
        end do
!
!  (I1,J1) x (I2,J2) still contains uncertainty, and has to go back
!  on the stack.
!
      else

        if ( nstack >= maxstack ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SS_GD_RPL - Fatal error!'
          write ( *, '(a)' ) '  NSTACK exceeds MAXSTACK.'
          return
        else
          nstack = nstack + 1
          stack(nstack,1) = p1
          stack(nstack,2) = p2
        end if

      end if
!
!  CONSIDER THE BLOCK (I2,J2) x (I3,J3).
!
!  ...a horizontal move.
!
     if ( i3 == i2 ) then

       do j = j2 + 1, j3 - 1
         npath = npath + 1
         pathi(npath) = i2
         pathj(npath) = j
       end do
!
!  ...a vertical move.
!
     else if ( j3 == j2  ) then

        do i = i2 + 1, i3 - 1
          npath = npath + 1
          pathi(npath) = i
          pathj(npath) = j2
        end do
!
!  ...a diagonal move.
!
     else if ( i3 == i2 + 1 .and. j3 == j2 + 1 ) then
!
!  ...a slide right and then down one.
!
      else if ( i3 == i2 + 1 .and. j3 > j2 + 1 ) then

        do j = j2 + 1, j3 - 1
          npath = npath + 1
          pathi(npath) = i2
          pathj(npath) = j
        end do
!
!  (I2,J2) x (I3,J3) still contains uncertainty, and goes back on the stack.
!
      else

        if ( nstack >= maxstack ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SS_GD_RPL - Fatal error!'
          write ( *, '(a)' ) '  NSTACK exceeds MAXSTACK.'
          return
        else
          nstack = nstack + 1
          stack(nstack,1) = p2
          stack(nstack,2) = p3
        end if

      end if

    end if

  end do
!
!  When the stack is empty, sort the path indices.
!
  call i4vec2_sort_a ( npath, pathi, pathj )
!
!  Now go through and mark gaps.
!
  do i = npath, 2, -1
    if ( pathi(i) == pathi(i-1) ) then
      pathi(i) = -1
    else if ( pathj(i) == pathj(i-1) ) then
      pathj(i) = -1
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
subroutine word_last_read ( line, word )

!*****************************************************************************80
!
!! WORD_LAST_READ returns the last word from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string containing words separated 
!    by spaces.
!
!    Output, character ( len = * ) WORD, the last word.
!
  implicit none

  integer first
  integer last
  character ( len = * ) line
  character ( len = * ) word

  last = len_trim ( line )

  if ( last <= 0 ) then
    word = ' '
    return
  end if

  first = last

  do while ( first > 1 )

    if ( line(first-1:first-1) /= ' ' ) then
      first = first - 1
    end if

  end do

  word = line(first:last)

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
