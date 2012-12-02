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
subroutine chvec_print ( m, a, title )

!*****************************************************************************80
!
!! CHVEC_PRINT prints a vector of characters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the length of the sequence.
!
!    Input, character A(M), the sequence.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer m

  character a(m)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, m
    write ( *, '(i3,2x,a1)' ) i, a(i)
  end do

  return
end
subroutine get_seed ( iseed )

!*****************************************************************************80
!
!! GET_SEED returns a seed for the random number generator.
!
!  Discussion:
!
!    The seed depends on the current time, and ought to be (slightly)
!    different every millisecond.  Once the seed is obtained, a random
!    number generator should be called a few times to further process
!    the seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ISEED, a pseudorandom seed value.
!
  implicit none

  integer, parameter :: I4_MAX = 2147483647

  integer iseed
  double precision temp
  integer values(8)

  character ( len = 10 ) time90
  character ( len = 8 ) today90
  character ( len = 5 ) zone

  call date_and_time ( today90, time90, zone, values )

  temp = 0.0

  temp = temp + dble ( values(2) - 1 ) / 11.0
  temp = temp + dble ( values(3) - 1 ) / 30.0
  temp = temp + dble ( values(5) ) / 23.0
  temp = temp + dble ( values(6) ) / 59.0
  temp = temp + dble ( values(7) ) / 59.0
  temp = temp + dble ( values(8) ) / 999.0
  temp = temp / 6.0

  if ( temp <= 0.0 ) then
    temp = 1.0 / 3.0
  else if ( temp >= 1.0 ) then
    temp = 2.0 / 3.0
  end if

  iseed = int ( dble ( I4_MAX ) * temp )
!
!  Never use a seed of 0 or I4_MAX.
!
  if ( iseed == 0 ) then
    iseed = 1
  end if

  if ( iseed == I4_MAX ) then
    iseed = I4_MAX - 1
  end if

  return
end
subroutine i4_random ( ilo, ihi, iseed, i)

!*****************************************************************************80
!
!! I4_RANDOM returns a random I4 in a given range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ILO, IHI, the minimum and maximum acceptable values.
!
!    Input/output, integer ISEED, a seed for the random number generator.
!
!    Output, integer I, the randomly chosen integer.
!
  implicit none

  integer i
  integer ihi
  integer ilo
  integer iseed
  real r
  real rhi
  real rlo
  real uniform_01_sample
!
!  Pick a random number in (0,1).
!
  r = uniform_01_sample ( iseed )
!
!  Set a real interval [RLO,RHI] which contains the integers [ILO,IHI],
!  each with a "neighborhood" of width 1.
!
  rlo = real ( ilo ) - 0.5
  rhi = real ( ihi ) + 0.5
!
!  Set I to the integer that is nearest the scaled value of R.
!
  i = nint ( ( 1.0 - r ) * rlo + r * rhi )
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
subroutine i4_to_amino_code ( i, c )

!*****************************************************************************80
!
!! I4_TO_AMINO_CODE converts an integer to an amino code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl Branden, John Tooze,
!    Introduction to Protein Structure,
!    Garland Publishing, 1991.
!
!  Parameters:
!
!    Input, integer I, the index of an amino acid, between 1 and 23.
!
!    Output, character C, the one letter code for an amino acid.
!
  implicit none

  integer, parameter :: n = 23

  character c
  character, dimension ( n ) :: c_table = (/ &
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', &
    'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', &
    'X', 'Y', 'Z' /)
  integer i

  if ( 1 <= i .and. i <= 23 ) then
    c = c_table(i)
  else
    c = '?'
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
!! I4VEC2_PRINT prints an I4VEC2
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
!! IVEC2_SORT_A ascending sorts an I4VEC2
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
!    27 June 2000
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
!! I4VEC_REVERSE reverses the elements of an I4VEC
!
!  Example:
!
!    Input:
!
!  N = 5,
!  A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!  A = ( 15, 14, 13, 12, 11 ).
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
subroutine mutate ( n_max, n, b, iseed )

!*****************************************************************************80
!
!! MUTATE applies a few mutations to a sequence.
!
!  Discussion:
!
!    This is just a draft routine, which applies some simple minded
!    mutations to a sequence.  There are lots of improvements that could
!    be made, but this was all I needed for now.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N_MAX, the dimension of B.
!
!    Input/output, integer N, the number of entries used in B.
!
!    Input/output, character B(N), the sequence.
!
!    Input/output, integer ISEED, a seed for the random number generator.
!
  implicit none

  integer n_max

  character b(n_max)
  character c
  integer i
  integer i_del
  integer i_ins
  integer i_mutate
  integer iseed
  integer j
  integer n
  integer n_del
  integer n_del_max
  integer n_ins
  integer n_ins_max
  integer n_mutate
  integer type
!
!  Pick the number of mutations to be applied.
!
  call i4_random ( 1, 3, iseed, n_mutate )

  i_mutate = 0

  do while ( i_mutate < n_mutate )

    call i4_random ( 1, 3, iseed, type )
!
!  Mutation of one item.
!
    if ( type == 1 ) then

      call i4_random ( 1, n, iseed, i )
      call i4_random ( 1, 23, iseed, j )

      call i4_to_amino_code ( j, c )

      b(i) = c

      i_mutate = i_mutate + 1
!
!  Deletion of several items.
!
    else if ( type == 2 ) then

      n_del_max = min ( n-1, 5 )

      if ( n_del_max > 0 ) then

        call i4_random ( 1, n_del_max, iseed, n_del )
        call i4_random ( 1, n+1-n_del, iseed, i_del )

        b(i_del:n-n_del) = b(i_del+n_del:n)

        n = n - n_del

        i_mutate = i_mutate + 1

      end if
!
!  Insertion of several items.
!
    else if ( type == 3 ) then

      n_ins_max = min ( n_max - n, 5 )

      if ( n_ins_max > 0 ) then

        call i4_random ( 1, n_ins_max, iseed, n_ins )

        call i4_random ( 0, n, iseed, i_ins )

        b(n+n_ins:i_ins+n_ins+1:-1) = b(n:i_ins+1:-1)

        do i = 1, n_ins
          call i4_random ( 1, 23, iseed, j )
          call i4_to_amino_code ( j, c )
          b(i_ins+i:i_ins+i) = c
        end do

        n = n + n_ins

        i_mutate = i_mutate + 1

      end if

    end if

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
!    Output, real PAM120_SCORE, the score for matching the two characters.
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
  real pam120_score
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
      write ( *, '(3x,23i3)' ) ( acid_index(j), j = 1, 26 )

    end if

  end if

  pam120_score = 0.0

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
!    Output, real PAM200_SCORE, the score for matching the two characters.
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
  real pam200_score
  integer, save, dimension ( acid_num, acid_num ) :: weight

  if ( need_data ) then

    call pam200 ( acid_code, weight )
    call a_index ( acid_num, acid_code, acid_index )
    need_data = .false.

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PAM200 substitution matrix.'
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

  pam200_score = 0.0

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
!    Input, real A(N), B(N), two arrays whose sum is to be examined.
!
!    Output, integer IMAX, the index of the largest entry in A+B.
!
  implicit none

  integer n

  real a(n)
  real b(n)
  integer i
  integer imax
  real sum_max

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
!    This routine is a sample scoring function which returns a score of:
!
!      -2 if either character is '-', representing a gap,
!      -1 if neither character is '-', but the characters are not identical,
!       0 if the characters are identical, and not '-'.
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
!    Output, real SIMPLE_SCORE, the score for matching the two characters.
!
  implicit none

  character c1
  character c2
  real score
  real simple_score

  if ( c1 == '-' .or. c2 == '-' ) then
    score = - 2.0
  else if ( c1 /= c2 ) then
    score = - 1.0
  else
    score = 0.0
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
subroutine ss_gg_bpq ( m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

!*****************************************************************************80
!
!! SS_GG_BPQ determines a global gap backward alignment path in quadratic space.
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
      write ( *, '(a)' ) 'SS_GG_BPQ: Unexpected situation!'
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
subroutine ss_gg_bsl ( a, b, m, m1, m2, n, n1, n2, ss_score, gap_open, &
  gap_extend, base, s, e, f, t )

!*****************************************************************************80
!
!! SS_GG_BSL determines a global gap backward alignment score in linear space.
!
!  Discussion:
!
!    The routine is called "linear" because it uses one N vector,
!    not an M by N array, to do the alignment.
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
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!    Michael Waterman,
!    Introduction to Computational Biology,
!    Chapman and Hall, 1995.
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
!    which returns a real value SS_SCORE for the matching of the character
!    C1 from sequence A to the character C2 from sequence B.
!
!    Input, real GAP_OPEN, GAP_EXTEND, the penalties for opening and
!    extending a gap.  A gap of length 7, for example, will result in a
!    penalty of GAP_OPEN + 7 * GAP_EXTEND.
!
!    Input, real BASE, an initial quantity to be added to certain penalties.
!
!    Output, real S(0:N), the backward optimal score vector.
!    The maximum possible alignment score is in S(N1).
!
!    Output, real E(0:N), the backward final insertion score vector.
!
!    Output, real F(0:N), the backward final deletion score vector.
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
  real base
  real diag_new
  real diag_old
  real e(0:n)
  real f(0:n)
  real gap_extend
  real gap_open
  integer i
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real s(0:n)
  real, external :: ss_score
  integer t(0:n)
!
!  The last row, I = M2.
!
  e(n2) = 0.0
  f(n2) = 0.0
  s(n2) = 0.0
  t(n2) = DUNNO

  if ( n2-1 >= n1 ) then
    e(n2-1) = e(n2) + gap_open + gap_extend
    f(n2-1) = f(n2) + 2.0 * gap_open + gap_extend
    s(n2-1) = s(n2) + gap_open + gap_extend
    t(n2-1) = INSERT
  end if

  do j = n2-2, n1, -1
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

    if ( i == m2-1 ) then
      e(n2) = base + gap_extend + gap_open
      f(n2) = base + gap_extend
      s(n2) = base + gap_extend
      t(n2) = DELETE
    else
      e(n2) = e(n2) + gap_extend
      f(n2) = f(n2) + gap_extend
      s(n2) = s(n2) + gap_extend
      t(n2) = DELETE
    end if

    do j = n2-1, n1, -1
!
!  Insertion.
!
      e(j) = e(j+1) + gap_extend
      if ( s(j+1) + gap_open + gap_extend > e(j) ) then
        e(j) = s(j+1) + gap_open + gap_extend
      end if
!
!  Deletion.
!
      f(j) = f(j) + gap_extend
      if ( s(j) + gap_open + gap_extend > f(j) ) then
        f(j) = s(j) + gap_open + gap_extend
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
subroutine ss_gg_bsq ( a, b, m, m1, m2, n, n1, n2, ss_score, gap_open, &
  gap_extend, base, lds, s, e, f, t )

!*****************************************************************************80
!
!! SS_GG_BSQ: global gap backward alignment score in quadratic space.
!
!  Discussion:
!
!    The routine can compute the full score table, or a sub-block.
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
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!    Michael Waterman,
!    Introduction to Computational Biology,
!    Chapman and Hall, 1995.
!
!  Parameters:
!
!    Input, character A(M), B(N), two sequences to be aligned.
!
!    Input, integer M, the number of entries in the sequence A.
!
!    Input, integer M1, M2, the lowest and highest table rows to compute.
!
!    Input, integer N, the number of entries in the sequence B.
!
!    Input, integer N1, N2, the lowest and highest table columns to compute.
!
!    Input, external SS_SCORE, the name of a function of the form
!      function ss_score ( c1, c2 )
!    which returns a real value SS_SCORE for the matching of the character
!    C1 from sequence A to the character C2 from sequence B.
!
!    Input, real GAP_OPEN, GAP_EXTEND, the penalties for opening and
!    extending a gap.  A gap of length 7, for example, will result in a
!    penalty of GAP_OPEN + 7 * GAP_EXTEND.
!
!    Input, real BASE, an initial quantity to be added to certain penalties.
!
!    Input, integer LDS, the declared upper first dimension of S.
!    LDS must be at least M.
!
!    Output, real S(0:LDS,0:N), the backward optimal score table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
!    Output, real E(0:LDS,0:N), the backward final insertion table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
!    Output, real F(0:LDS,0:N), the backward final deletion table.
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
  real base
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real gap_extend
  real gap_open
  integer i
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real s(0:lds,0:n)
  real, external :: ss_score
  integer t(0:lds,0:n)
!
!  Lower Right corner.
!
  e(m2,n2) = 0.0
  f(m2,n2) = 0.0
  s(m2,n2) = 0.0
  t(m2,n2) = DUNNO
!
!  Lower Left row.
!
  if ( n2-1 >= n1 ) then
    e(m2,n2-1) = e(m2,n2) + gap_open + gap_extend
    f(m2,n2-1) = f(m2,n2) + 2.0 * gap_open + gap_extend
    s(m2,n2-1) = s(m2,n2) + gap_open + gap_extend
    t(m2,n2-1) = INSERT
  end if

  do j = n2-2, n1, -1
    e(m2,j) = e(m2,j+1) + gap_extend
    f(m2,j) = f(m2,j+1) + gap_extend
    s(m2,j) = s(m2,j+1) + gap_extend
    t(m2,j) = INSERT
  end do
!
!  Upper rectangle.
!
  do i = m2-1, m1, -1

    if ( i == m2-1 ) then
      e(i,n2) = base + gap_extend + gap_open
      f(i,n2) = base + gap_extend
      s(i,n2) = base + gap_extend
      t(i,n2) = DELETE
    else
      e(i,n2) = e(i+1,n2) + gap_extend
      f(i,n2) = f(i+1,n2) + gap_extend
      s(i,n2) = s(i+1,n2) + gap_extend
      t(i,n2) = DELETE
    end if

    do j = n2-1, n1, -1
!
!  Insertion.
!
      e(i,j) = e(i,j+1) + gap_extend

      if ( s(i,j+1) + gap_open + gap_extend > e(i,j) ) then
        e(i,j) = s(i,j+1) + gap_open + gap_extend
      end if
!
!  Deletion.
!
      f(i,j) = f(i+1,j) + gap_extend

      if ( s(i+1,j) + gap_open + gap_extend > f(i,j) ) then
        f(i,j) = s(i+1,j) + gap_open + gap_extend
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
subroutine ss_gg_fpq ( m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

!*****************************************************************************80
!
!! SS_GG_FPQ determines a global gap forward alignment path in quadratic space.
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
subroutine ss_gg_fsl ( a, b, m, m1, m2, n, n1, n2, ss_score, gap_open, &
  gap_extend, base, s, e, f, t )

!*****************************************************************************80
!
!! SS_GG_FSL determines a global gap forward alignment score in linear space.
!
!  Discussion:
!
!    The routine is called "linear" because it uses one N vector,
!    not an M by N array, to do the alignment.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 June 2000
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
!    Michael Waterman,
!    Introduction to Computational Biology,
!    Chapman and Hall, 1995.
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
!    which returns a real value SS_SCORE for the matching of the character
!    C1 from sequence A to the character C2 from sequence B.
!
!    Input, real GAP_OPEN, GAP_EXTEND, the penalties for opening and
!    extending a gap.  A gap of length 7, for example, will result in a
!    penalty of GAP_OPEN + 7 * GAP_EXTEND.
!
!    Input, real BASE, an initial quantity to be added to certain penalties.
!
!    Output, real S(0:N), the forward optimal score vector.
!    The maximum possible alignment score is in S(N2).
!
!    Output, real E(0:N), the forward final insertion score vector.
!
!    Output, real F(0:N), the forward final deletion score vector.
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
  real base
  real diag_new
  real diag_old
  real e(0:n)
  real f(0:n)
  real gap_extend
  real gap_open
  integer i
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real s(0:n)
  real, external :: ss_score
  integer t(0:n)
!
!  The first row, I = M1.
!
  e(n1) = 0.0
  f(n1) = 0.0
  s(n1) = 0.0
  t(n1) = DUNNO

  if ( n1+1 <= n2 ) then
    e(n1+1) = e(n1) + gap_open + gap_extend
    f(n1+1) = f(n1) + 2.0 * gap_open + gap_extend
    s(n1+1) = s(n1) + gap_open + gap_extend
    t(n1+1) = INSERT
  end if

  do j = n1+2, n2
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

    if ( i == m1+1 ) then
      e(n1) = base + gap_extend + gap_open
      f(n1) = base + gap_extend
      s(n1) = base + gap_extend
      t(n1) = DELETE
    else
      e(n1) = e(n1) + gap_extend
      f(n1) = f(n1) + gap_extend
      s(n1) = s(n1) + gap_extend
      t(n1) = DELETE
    end if

    do j = n1+1, n2
!
!  Insertion
!
      e(j) = e(j-1) + gap_extend

      if ( s(j-1) + gap_open + gap_extend > e(j) ) then
        e(j) = s(j-1) + gap_open + gap_extend
      end if
!
!  Deletion.
!
      f(j) = f(j) + gap_extend

      if ( s(j) + gap_open + gap_extend > f(j) ) then
        f(j) = s(j) + gap_open + gap_extend
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
subroutine ss_gg_fsq ( a, b, m, m1, m2, n, n1, n2, ss_score, gap_open, &
  gap_extend, base, lds, s, e, f, t )

!*****************************************************************************80
!
!! SS_GG_FSQ determines a global gap forward alignment score in quadratic space.
!
!  Discussion:
!
!    The routine can compute the full score table, or a sub-block.
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
!  Reference:
!
!    Kun-Mao Chao, Ross Hardison, Webb Miller,
!    Recent Developments in Linear-Space Alignment Methods: A Survey,
!    Journal of Computational Biology,
!    Volume 1, Number 4, 1994, pages 271-291.
!
!    Michael Waterman,
!    Introduction to Computational Biology,
!    Chapman and Hall, 1995.
!
!  Parameters:
!
!    Input, character A(M), B(N), two sequences to be aligned.
!
!    Input, integer M, the number of entries in sequence A.
!
!    Input, integer M1, M2, the lowest and highest table rows to compute.
!
!    Input, integer N, the number of entries in sequence B.
!
!    Input, integer N1, N2, the lowest and highest table columns to compute.
!
!    Input, external SS_SCORE, the name of a function of the form
!      function ss_score ( c1, c2 )
!    which returns a real value SS_SCORE for the matching of the character
!    C1 from sequence A to the character C2 from sequence B.
!
!    Input, real GAP_OPEN, GAP_EXTEND, the penalties for opening and
!    extending a gap.  A gap of length 7, for example, will result in a
!    penalty of GAP_OPEN + 7 * GAP_EXTEND.
!
!    Input, real BASE, an initial quantity to be added to certain penalties.
!
!    Input, integer LDS, the declared upper first dimension of S, E, and F.
!    LDS must be at least M.
!
!    Output, real S(0:LDS,0:N), the forward optimal score table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
!    Output, real E(0:LDS,0:N), the forward final insertion table.
!    Entries in the M1:M2 by N1:N2 block have been computed.
!
!    Output, real F(0:LDS,0:N), the forward final deletion table.
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
  real base
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real gap_extend
  real gap_open
  integer i
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real s(0:lds,0:n)
  real, external :: ss_score
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
  if ( n1+1 <= n2 ) then
    e(m1,n1+1) = e(m1,n1) + gap_open + gap_extend
    f(m1,n1+1) = f(m1,n1) + 2.0 * gap_open + gap_extend
    s(m1,n1+1) = s(m1,n1) + gap_open + gap_extend
    t(m1,n1+1) = INSERT
  end if

  do j = n1+2, n2
    e(m1,j) = e(m1,j-1) + gap_extend
    f(m1,j) = f(m1,j-1) + gap_extend
    s(m1,j) = s(m1,j-1) + gap_extend
    t(m1,j) = INSERT
  end do
!
!  Lower rectangle.
!
  do i = m1+1, m2

    if ( i == m1+1 ) then
      e(i,n1) = base + gap_extend + gap_open
      f(i,n1) = base + gap_extend
      s(i,n1) = base + gap_extend
      t(i,n1) = DELETE
    else
      e(i,n1) = e(i-1,n1) + gap_extend
      f(i,n1) = f(i-1,n1) + gap_extend
      s(i,n1) = s(i-1,n1) + gap_extend
      t(i,n1) = DELETE
    end if

    do j = n1+1, n2
!
!  Insertion.
!
      e(i,j) = e(i,j-1) + gap_extend

      if ( s(i,j-1) + gap_open + gap_extend > e(i,j) ) then
        e(i,j) = s(i,j-1) + gap_open + gap_extend
      end if
!
!  Deletion.
!
      f(i,j) = f(i-1,j) + gap_extend

      if ( s(i-1,j) + gap_open + gap_extend > f(i,j) ) then
        f(i,j) = s(i-1,j) + gap_open + gap_extend
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
subroutine ss_gg_match_print ( a, b, m, n, npath, pathi, pathj, ss_score, &
  gap_open, gap_extend )

!*****************************************************************************80
!
!! SS_GG_MATCH_PRINT prints a global gap alignment.
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
!    Input, integer NPATH, the number of alignments.
!
!    Input, integer PATHI(M+N+1), PATHJ(M+N+1), contains, in the first
!    NPATH entries, the indices of the aligned items.
!    The first entries are special marker values:
!      PATHI(1) = 0 and PATHJ(1) = 0;
!    A value of -1 for PATHI or PATHJ indicates a null partner.
!    Otherwise, A(PATHI(I)) is matched to B(PATHJ(I)).
!
!    Input, external SS_SCORE, the name of a function of the form
!      function ss_score ( c1, c2 )
!    which returns a real value SS_SCORE for the matching of the character
!    C1 from sequence A to the character C2 from sequence B.
!
!    Input, real GAP_OPEN, GAP_EXTEND, the penalties for opening and
!    extending a gap.  A gap of length 7, for example, will result in a
!    penalty of GAP_OPEN + 7 * GAP_EXTEND.
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
  real gap_extend
  real gap_open
  integer i
  integer i_old
  real inc
  integer ipath
  integer j
  integer j_old
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: ss_score
  real sum2

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

  i_old = pathi(1)
  j_old = pathj(1)

  do ipath = 2, npath

    i = pathi(ipath)
    j = pathj(ipath)

    if ( i == -1 ) then

      c1 = ' '
      c2 = '|'
      c3 = ' '
      c4 = b(j)
      write ( c5, '(i3)' ) j

      if ( i_old /= -1 ) then
        inc = gap_open + gap_extend
      else
        inc = gap_extend
      end if

    else if ( j == -1 ) then

      write ( c1, '(i3)' ) i
      c2 = a(i)
      c3 = ' '
      c4 = '|'
      c5 = ' '

      if ( j_old /= -1 ) then
        inc = gap_open + gap_extend
      else
        inc = gap_extend
      end if

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

    i_old = i
    j_old = j

    sum2 = sum2 + inc

    write ( *, '(a3,2x,a1,2x,a1,2x,a1,2x,a3,2x,2f10.2)' ) &
      c1, c2, c3, c4, c5, inc, sum2

  end do

  return
end
subroutine ss_gg_match_score ( a, b, m, n, npath, pathi, pathj, ss_score, &
  gap_open, gap_extend, score )

!*****************************************************************************80
!
!! SS_GG_MATCH_SCORE scores a global gap alignment.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 June 2000
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
!    Input, integer NPATH, the number of alignments.
!
!    Input, integer PATHI(M+N+1), PATHJ(M+N+1), contains, in the first
!    NPATH entries, the indices of the aligned items.
!    The first entries are special marker values:
!      PATHI(1) = 0 and PATHJ(1) = 0;
!    A value of -1 for PATHI or PATHJ indicates a null partner.
!    Otherwise, A(PATHI(I)) is matched to B(PATHJ(I)).
!
!    Input, external SS_SCORE, the name of a function of the form
!      function ss_score ( c1, c2 )
!    which returns a real value SS_SCORE for the matching of the character
!    C1 from sequence A to the character C2 from sequence B.
!
!    Input, real GAP_OPEN, GAP_EXTEND, the penalties for opening and
!    extending a gap.  A gap of length 7, for example, will result in a
!    penalty of GAP_OPEN + 7 * GAP_EXTEND.
!
!    Output, real SCORE, the score for the given matching.
!
  implicit none

  integer m
  integer n

  character a(m)
  character b(n)
  real gap_extend
  real gap_open
  integer i
  integer i_old
  real inc
  integer ipath
  integer j
  integer j_old
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real score
  real, external :: ss_score
  real sum2

  sum2 = 0.0

  ipath = 1
  i_old = pathi(1)
  j_old = pathj(1)

  do ipath = 2, npath

    i = pathi(ipath)
    j = pathj(ipath)

    if ( i == -1 ) then

      if ( i_old /= -1 ) then
        inc = gap_open + gap_extend
      else
        inc = gap_extend
      end if

    else if ( j == -1 ) then

      if ( j_old /= -1 ) then
        inc = gap_open + gap_extend
      else
        inc = gap_extend
      end if

    else

      inc = ss_score ( a(i), b(j) )

    end if

    i_old = i
    j_old = j

    sum2 = sum2 + inc

  end do

  score = sum2

  return
end
subroutine ss_gg_rpl ( a, b, m, m1, m2, n, n1, n2, ss_score, gap_open, &
  gap_extend, npath, pathi, pathj )

!*****************************************************************************80
!
!! SS_GG_RPL determines a global gap recursive alignment path in linear space.
!
!  Discussion:
!
!    The routine is called "linear" because it uses a few vectors of
!    dimension N to determine the alignment.
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
!    Eugene Myers and Webb Miller,
!    Optimal Alignments in Linear Space,
!    CABIOS, volume 4, number 1, 1988, pages 11-17.
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
!    Input, integer M1, M2, the minimum and maximum entries of A to align.
!    For a full alignment, use M1 = 0, M2 = M.
!
!    Input, integer N, the number of entries in sequence B.
!
!    Input, integer N1, N2, the minimum and maximum entries of B to align.
!    For a full alignment, use N1 = 0, N2 = N.
!
!    Input, external SS_SCORE, the name of a function of the form
!      function ss_score ( c1, c2 )
!    which returns a real value SS_SCORE for the matching of the character
!    C1 from sequence A to the character C2 from sequence B.
!
!    Input, real GAP_OPEN, GAP_EXTEND, the penalties for opening and
!    extending a gap.  A gap of length 7, for example, will result in a
!    penalty of GAP_OPEN + 7 * GAP_EXTEND.
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

  integer, parameter :: DELETE = 4
  integer, parameter :: INSERT = 2
  integer, parameter :: MATCH = 1
  integer, parameter :: stack_max = 200

  integer m
  integer n

  character a(m)
  character b(n)
  real base
  real eb(0:n)
  real ef(0:n)
  real fb(0:n)
  real ff(0:n)
  real gap_open
  real gap_extend
  integer i
  integer i1sub
  integer i2sub
  integer i3sub
  integer j
  integer j1sub
  integer j1type
  integer j2sub
  integer j2type
  integer j3sub
  integer j3type
  integer ja
  integer jb
  integer m1
  integer m2
  integer n1
  integer n2
  integer nb1
  integer nb2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real sb(0:n)
  real sf(0:n)
  real, external :: ss_score
  integer stack1(6,stack_max)
  real stack2(2,stack_max)
  integer :: stack_num
  real t1
  real t2
  real t3
  integer tb(0:n)
  integer tf(0:n)
  real x
  real y
  real y1
  real y2
!
!  We begin with the problem of analyzing the entire M1:M2 by N1:N2 table.
!
  i1sub = m1
  i2sub = m2
  j1sub = n1
  j1type = MATCH
  j2sub = n2
  j2type = MATCH
  t1 = gap_open
  t2 = gap_open
!
!  Put the two endpoints on the path.
!
  npath = 1
  pathi(npath) = i1sub
  pathj(npath) = j1sub

  npath = 2
  pathi(npath) = i2sub
  pathj(npath) = j2sub
!
!  Put the initial problem on the stack.
!
  stack_num = 0

  call ss_gg_rpl_push ( i1sub, j1sub, j1type, t1, i2sub, j2sub, j2type, t2, &
    stack1, stack2, stack_num, stack_max )

  do while ( stack_num > 0 )
!
!  Pop the next problem off the stack.
!
    call ss_gg_rpl_pop ( i1sub, j1sub, j1type, t1, i3sub, j3sub, j3type, t3, &
      stack1, stack2, stack_num, stack_max )
!
!  Refuse to handle improperly described subregions.
!
    if ( i1sub > i3sub .or. j1sub > j3sub ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SS_GG_RPL - Fatal error!'
      write ( *, '(a)' ) '  The indices describing a subregion have the'
      write ( *, '(a)' ) '  wrong sense.'
      write ( *, '(a,i6)' ) '  I1 = ', i1sub
      write ( *, '(a,i6)' ) '  I3 = ', i3sub
      write ( *, '(a,i6)' ) '  J1 = ', j1sub
      write ( *, '(a,i6)' ) '  J3 = ', j3sub
      stop
!
!  Null regions require no processing.
!
    else if ( i1sub == i3sub .and. j1sub == j3sub ) then
!
!  A vertical strip is easy.
!
    else if ( j1sub == j3sub ) then

      do i = i1sub+1, i3sub-1
        npath = npath + 1
        pathi(npath) = i
        pathj(npath) = j1sub
      end do
!
!  A horizontal strip is easy.
!
    else if ( i1sub == i3sub ) then

      do j = j1sub+1, j3sub-1
        npath = npath + 1
        pathi(npath) = i1sub
        pathj(npath) = j
      end do
!
!  For the case where the uncertainty region is two units high, the path
!  can be described as going from J1SUB to JA in row I1SUB, and then
!  from JB in row I3SUB to J3SUB.
!
!  We need to know:
!
!    Does the incoming path to (I1SUB,J1SUB) represent
!      a match,
!      a deletion of A's,
!      an insertion of B's?
!
!    Does the outgoing path from (I3SUB,J3SUB) represent
!      a match,
!      a deletion of A's,
!      an insertion of B's.
!
!    When done, we have to record that the incoming and outgoing paths
!    to (I2SUB,J2SUB) represent:
!      a match,
!      a deletion of A's,
!      an insertion of B's.
!
    else if ( i3sub == i1sub + 1 ) then

      x = - huge ( x )
      ja = 0
      jb = 0
!
!  A: Cost of the path:
!
!    X X X ... X X . ... . . .
!    . . . ... . X X ... X X X
!
      do j = j1sub, j3sub

        nb1 = j - j1sub
        nb2 = j3sub - j

        y = 0.0

        if ( nb1 > 0 ) then
          if ( j1type /= INSERT ) then
            y = y + gap_open
          end if
          y = y + gap_extend * nb1
        end if

        y =  y + gap_open + gap_extend

        if ( nb2 > 0 ) then
          if ( j3type /= INSERT ) then
            y = y + gap_open
          end if
          y = y + gap_extend * nb2
        end if

        if ( y > x ) then
          x = y
          ja = j
          jb = j
        end if

      end do
!
!  B: Cost of the path:
!
!    X X X ... X . ... . . .
!    . . . ... . X ... X X X
!
      do j = j1sub, j3sub - 1

        nb1 = j - j1sub
        nb2 = j3sub - j - 1

        y = 0.0

        if ( nb1 > 0 ) then
          if ( j1type /= INSERT ) then
            y = y + gap_open
          end if
          y = y + gap_extend * nb1
        end if

        y = y + ss_score ( a(i3sub), b(j+1) )

        if ( nb2 > 0 ) then

          if ( j3type /= INSERT ) then
            y = y + gap_open
          end if

          y = y + gap_extend * nb2

        end if

        if ( y > x ) then
          x = y
          ja = j
          jb = j + 1
        end if

      end do
!
!  Now fill in the path.
!
      do j = j1sub + 1, ja
        npath = npath + 1
        pathi(npath) = i1sub
        pathj(npath) = j
      end do

      do j = jb, j3sub - 1
        npath = npath + 1
        pathi(npath) = i3sub
        pathj(npath) = j
      end do

    else

      i2sub = ( i1sub + i3sub ) / 2

      base = t1

      call ss_gg_fsl ( a, b, m, i1sub, i2sub, n, j1sub, j3sub, ss_score, &
        gap_open, gap_extend, base, sf, ef, ff, tf )

      base = t3

      call ss_gg_bsl ( a, b, m, i2sub, i3sub, n, j1sub, j3sub, ss_score, &
        gap_open, gap_extend, base, sb, eb, fb, tb )
!
!  Find J2SUB, the value of J between J1SUB and J3SUB that maximizes
!  SF(J)+SB(J) or FF(J)+FB(J)-GAP_OPEN.
!
      j = j1sub
      y1 = sf(j) + sb(j)

      x = y1
      j2sub = j
      j2type = MATCH

      y2 = ff(j) + fb(j) - gap_open

      if ( x < y2 ) then
        x = y2
        j2sub = j
        j2type = DELETE
      end if

      do j = j1sub+1, j3sub

        y1 = sf(j) + sb(j)

        if ( y1 > x ) then
          x = y1
          j2sub = j
          j2type = MATCH
        end if

        y2 = ff(j) + fb(j) - gap_open

        if ( y2 > x ) then
          x = y2
          j2sub = j
          j2type = DELETE
        end if

      end do

      npath = npath + 1
      pathi(npath) = i2sub
      pathj(npath) = j2sub

      if ( j2type == MATCH ) then

        t2 = gap_open

        call ss_gg_rpl_push ( i1sub, j1sub, j1type, t1, i2sub, j2sub, j2type, &
          t2, stack1, stack2, stack_num, stack_max )

        call ss_gg_rpl_push ( i2sub, j2sub, j2type, t2, i3sub, j3sub, j3type, &
          t3, stack1, stack2, stack_num, stack_max )

      else if ( j2type == DELETE ) then

        if ( ( i1sub < i2sub-1 .or. j1sub < j2sub ) .and. &
             ( i2sub+1 < i3sub .or. j2sub < j3sub ) ) then

          npath = npath + 1
          pathi(npath) = i2sub - 1
          pathj(npath) = j2sub

          npath = npath + 1
          pathi(npath) = i2sub + 1
          pathj(npath) = j2sub

          t2 = 0.0

          call ss_gg_rpl_push ( i1sub, j1sub, j1type, t1, i2sub-1, j2sub, &
            j2type, t2, stack1, stack2, stack_num, stack_max )

          call ss_gg_rpl_push ( i2sub+1, j2sub, j2type, t2, i3sub, j3sub, &
            j3type, t3, stack1, stack2, stack_num, stack_max )

        else if ( i2sub+1 < i3sub .or. j2sub < j3sub ) then

          npath = npath + 1
          pathi(npath) = i2sub + 1
          pathj(npath) = j2sub

          t2 = t1

          call ss_gg_rpl_push ( i2sub+1, j2sub, j2type, t2, i3sub, j3sub, &
            j3type, t3, stack1, stack2, stack_num, stack_max )

        else if ( i1sub < i2sub-1 .or. j1sub < j2sub ) then

          npath = npath + 1
          pathi(npath) = i2sub - 1
          pathj(npath) = j2sub

          t2 = t3

          call ss_gg_rpl_push ( i1sub, j1sub, j1type, t1, i2sub-1, j2sub, &
           j2type, t2, stack1, stack2, stack_num, stack_max )

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
subroutine ss_gg_rpl_pop ( i1sub, j1sub, j1type, t1, i2sub, j2sub, j2type, &
  t2, stack1, stack2, stack_num, stack_max )

!*****************************************************************************80
!
!! SS_GG_RPL_POP pops the data describing a subproblem off of the stack.
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
!  Parameters:
!
!    Output, integer I1SUB, J1SUB, J1TYPE, real T1, the row and column of the
!    first cell, the type of the match there, and the appropriate base.
!
!    Output, integer I2SUB, J2SUB, J2TYPE, real T2, the row and column of the
!    second cell, the type of the match there, and the appropriate base.
!
!    Input, integer STACK1(6,STACK_MAX), real STACK2(2,STACK_MAX), two
!    arrays in which stack data is stored.
!
!    Input/output, integer STACK_NUM, a pointer to the most recent item
!    on the stack.
!
!    Input, integer STACK_MAX, the maximum number of items in the stack.
!
  implicit none

  integer stack_max

  integer i1sub
  integer i2sub
  integer j1sub
  integer j1type
  integer j2sub
  integer j2type
  integer stack1(6,stack_max)
  real stack2(2,stack_max)
  integer stack_num
  real t1
  real t2

  if ( stack_num < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SS_GG_RPL_POP - Fatal error!'
    write ( *, '(a)' ) '  No more data to pop.'
    stop

  end if

  i1sub  = stack1(1,stack_num)
  i2sub  = stack1(2,stack_num)
  j1sub  = stack1(3,stack_num)
  j1type = stack1(4,stack_num)
  j2sub  = stack1(5,stack_num)
  j2type = stack1(6,stack_num)

  t1     = stack2(1,stack_num)
  t2     = stack2(2,stack_num)

  stack_num = stack_num - 1

  return
end
subroutine ss_gg_rpl_push ( i1sub, j1sub, j1type, t1, i2sub, j2sub, j2type, &
  t2, stack1, stack2, stack_num, stack_max )

!*****************************************************************************80
!
!! SS_GG_RPL_PUSH pushes the data describing a subproblem onto the stack.
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
!  Parameters:
!
!    Input, integer I1SUB, J1SUB, J1TYPE, real T1, the row and column of the
!    first cell, the type of the match there, and the appropriate base.
!
!    Input, integer I2SUB, J2SUB, J2TYPE, real T2, the row and column of the
!    second cell, the type of the match there, and the appropriate base.
!
!    Input/output, integer STACK1(6,STACK_MAX), real STACK2(2,STACK_MAX), two
!    arrays in which stack data is stored.
!
!    Input/output, integer STACK_NUM, a pointer to the most recent item
!    on the stack.
!
!    Input, integer STACK_MAX, the maximum number of items in the stack.
!
  implicit none

  integer stack_max

  integer i1sub
  integer i2sub
  integer j1sub
  integer j1type
  integer j2sub
  integer j2type
  integer stack1(6,stack_max)
  real stack2(2,stack_max)
  integer stack_num
  real t1
  real t2
!
!  You might be out of stack space.
!
  if ( stack_num >= stack_max ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SS_GG_RPL_PUSH - Fatal error!'
    write ( *, '(a)' ) '  No more room on the stack.'
    stop
  end if

  stack_num = stack_num + 1

  stack1(1,stack_num) = i1sub
  stack1(2,stack_num) = i2sub
  stack1(3,stack_num) = j1sub
  stack1(4,stack_num) = j1type
  stack1(5,stack_num) = j2sub
  stack1(6,stack_num) = j2type

  stack2(1,stack_num) = t1
  stack2(2,stack_num) = t2

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
function uniform_01_sample ( iseed )

!*****************************************************************************80
!
!! UNIFORM_01_SAMPLE is a portable random number generator.
!
!  Discussion:
!
!    ISEED = ISEED * (7**5) mod (2**31 - 1)
!    RANDOM = ISEED * / ( 2**31 - 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ISEED, the integer "seed" used to generate
!    the output random number, and updated in preparation for the
!    next one.  ISEED should not be zero.
!
!    Output, real UNIFORM_01_SAMPLE, a random value between 0 and 1.
!
!  Local Parameters:
!
!    IA = 7**5
!    IB = 2**15
!    IB16 = 2**16
!    IP = 2**31-1
!
  implicit none

  integer, parameter :: ia = 16807
  integer, parameter :: ib15 = 32768
  integer, parameter :: ib16 = 65536
  integer, parameter :: ip = 2147483647
  integer iprhi
  integer iseed
  integer ixhi
  integer k
  integer leftlo
  integer loxa
  real uniform_01_sample
!
!  Don't let ISEED be 0.
!
  if ( iseed == 0 ) then
    iseed = ip
  end if
!
!  Get the 15 high order bits of ISEED.
!
  ixhi = iseed / ib16
!
!  Get the 16 low bits of ISEED and form the low product.
!
  loxa = ( iseed - ixhi * ib16 ) * ia
!
!  Get the 15 high order bits of the low product.
!
  leftlo = loxa / ib16
!
!  Form the 31 highest bits of the full product.
!
  iprhi = ixhi * ia + leftlo
!
!  Get overflow past the 31st bit of full product.
!
  k = iprhi / ib15
!
!  Assemble all the parts and presubtract IP.  The parentheses are
!  essential.
!
  iseed = ( ( ( loxa - leftlo * ib16 ) - ip ) + ( iprhi - k * ib15 ) * ib16 ) &
    + k
!
!  Add IP back in if necessary.
!
  if ( iseed < 0 ) then
    iseed = iseed + ip
  end if
!
!  Multiply by 1 / (2**31-1).
!
  uniform_01_sample = real ( iseed ) * 4.656612875e-10

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
!  Discussion:
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
