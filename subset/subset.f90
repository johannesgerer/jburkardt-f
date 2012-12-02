subroutine asm_enum ( n, asm_num )

!*****************************************************************************80
!
!! ASM_ENUM returns the number of alternating sign matrices of a given order.
!
!  Discussion:
!
!    N     ASM_NUM
!
!    0       1
!    1       1
!    2       2
!    3       7
!    4      42
!    5     429
!    6    7436
!    7  218348
!
!    A direct formula is
!
!      ASM_NUM ( N ) = product ( 0 <= I <= N-1 ) ( 3 * I + 1 )! / ( N + I )!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrices.
!
!    Output, integer ( kind = 4 ) ASM_NUM, the number of alternating sign
!    matrices of order N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n+1)
  integer ( kind = 4 ) asm_num
  integer ( kind = 4 ) b(n+1)
  integer ( kind = 4 ) c(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nn

  asm_num = 0

  nn = n + 1

  if ( n + 1 <= 0 ) then
    return
  end if
!
!  Row 1
!
  a(1) = 1

  if ( n + 1 == 1 ) then
    asm_num = a(1)
    return
  end if
!
!  Row 2
!
  a(1) = 1

  if ( n + 1 == 2 ) then
    asm_num = a(1)
    return
  end if

  b(1) = 2
  c(1) = 2
  a(2) = 1
!
!  Row 3 and on.
!
  do nn = 3, n

    b(nn-1) = nn
    do i = nn-2, 2, -1
      b(i) = b(i) + b(i-1)
    end do
    b(1) = 2

    c(nn-1) = 2
    do i = nn-2, 2, -1
      c(i) = c(i) + c(i-1)
    end do
    c(1) = nn

    a(1) = sum ( a(1:nn-1) )
    do i = 2, nn
      a(i) = a(i-1) * c(i-1) / b(i-1)
    end do

  end do

  asm_num = sum ( a(1:n) )

  return
end
subroutine asm_triangle ( n, a )

!*****************************************************************************80
!
!! ASM_TRIANGLE returns a row of the alternating sign matrix triangle.
!
!  Discussion:
!
!    The first seven rows of the triangle are as follows:
!
!          1      2      3      4      5      6     7
!
!    0     1
!    1     1      1
!    2     2      3      2
!    3     7     14     14      7
!    4    42    105    135    105     42
!    5   429   1287   2002   2002   1287    429
!    6  7436  26026  47320  56784  47320  26026  7436
!
!    For a given N, the value of A(J) represents entry A(I,J) of
!    the triangular matrix, and gives the number of alternating sign matrices
!    of order N in which the (unique) 1 in row 1 occurs in column J.
!
!    Thus, of alternating sign matrices of order 3, there are
!    2 with a leading 1 in column 1:
!
!      1 0 0  1 0 0
!      0 1 0  0 0 1
!      0 0 1  0 1 0
!
!    3 with a leading 1 in column 2, and
!
!      0 1 0  0 1 0  0 1 0
!      1 0 0  0 0 1  1-1 1
!      0 0 1  1 0 0  0 1 0
!
!    2 with a leading 1 in column 3:
!
!      0 0 1  0 0 1
!      1 0 0  0 1 0
!      0 1 0  1 0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the desired row.
!
!    Output, integer ( kind = 4 ) A(N+1), the entries of the row.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n+1)
  integer ( kind = 4 ) b(n+1)
  integer ( kind = 4 ) c(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nn

  if ( n+1 <= 0 ) then
    return
  end if
!
!  Row 1
!
  a(1) = 1

  if ( n + 1 == 1 ) then
    return
  end if
!
!  Row 2
!
  nn = 2
  b(1) = 2
  c(1) = nn

  a(1) = sum ( a(1:nn-1) )
  do i = 2, nn
    a(i) = a(i-1) * c(i-1) / b(i-1)
  end do

  if ( n + 1 == 2 ) then
    return
  end if
!
!  Row 3 and on.
!
  do nn = 3, n + 1

    b(nn-1) = nn
    do i = nn-2, 2, -1
      b(i) = b(i) + b(i-1)
    end do
    b(1) = 2

    c(nn-1) = 2
    do i = nn - 2, 2, -1
      c(i) = c(i) + c(i-1)
    end do
    c(1) = nn

    a(1) = sum ( a(1:nn-1) )
    do i = 2, nn
      a(i) = a(i-1) * c(i-1) / b(i-1)
    end do

  end do

  return
end
subroutine bell ( n, b )

!*****************************************************************************80
!
!! BELL returns the Bell numbers from 0 to N.
!
!  Discussion:
!
!    The Bell number B(N) is defined as the number of partitions (of
!    any size) of a set of N distinguishable objects.
!
!    A partition of a set is a division of the objects of the set into
!    subsets.
!
!    The Bell number B(N) is the number of restricted growth functions
!    on N.
!
!    Note that the Stirling numbers of the second kind, S^m_n, count the
!    number of partitions of N objects into M classes, and so it is
!    true that
!
!      B(N) = S^1_N + S^2_N + ... + S^N_N.
!
!  Example:
!
!    There are 15 partitions of a set of 4 objects:
!
!      (1234), (123)(4), (124)(3), (12)(34), (12)(3)(4),
!      (134)(2), (13)(24), (13)(2)(4), (14)(23), (1)(234),
!      (1)(23)(4), (14)(2)(3), (1)(24)(3), (1)(2)(34), (1)(2)(3)(4)
!
!    and so B(4) = 15.
!
!  First values:
!
!     N         B(N)
!     0           1
!     1           1
!     2           2
!     3           5
!     4          15
!     5          52
!     6         203
!     7         877
!     8        4140
!     9       21147
!    10      115975
!
!  Recursion:
!
!    B(I) = sum ( 1 <= J <= I ) Binomial ( I-1, J-1 ) * B(I-J)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of Bell numbers desired.
!
!    Output, integer ( kind = 4 ) B(0:N), the Bell numbers from 0 to N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) b(0:n)
  integer ( kind = 4 ) combo
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) j

  b(0) = 1

  do i = 1, n
    b(i) = 0
    do j = 1, i
      combo = i4_choose ( i-1, j-1 )
      b(i) = b(i) + combo * b(i-j)
    end do
  end do

  return
end
subroutine bell_values ( n_data, n, c )

!*****************************************************************************80
!
!! BELL_VALUES returns some values of the Bell numbers for testing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the order of the Bell number.
!
!    Output, integer ( kind = 4 ) C, the value of the Bell number.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 11

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( nmax ) :: c_vec = (/ &
    1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
     0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine binary_vector_next ( n, bvec )

!*****************************************************************************80
!
!! BINARY_VECTOR_NEXT generates the next binary vector.
!
!  Discussion:
!
!    A binary vector is a vector whose entries are 0 or 1.
!
!    The user inputs an initial zero vector to start.  The program returns
!    the "next" vector.
!
!    The vectors are produced in the order:
!
!    ( 0, 0, 0, ..., 0 )
!    ( 1, 0, 0, ..., 0 )
!    ( 0, 1, 0, ..., 0 )
!    ( 1, 1, 0, ..., 0 )
!    ( 0, 0, 1, ..., 0 )
!    ( 1, 0, 1, ..., 0 )
!               ...
!    ( 1, 1, 1, ..., 1)
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
!    0 0 0  =>  1 0 0
!    1 0 0  =>  0 1 0
!    0 1 0  =>  1 1 0
!    1 1 0  =>  0 0 1
!    0 0 1  =>  1 0 1
!    1 0 1  =>  0 1 1
!    0 1 1  =>  1 1 1
!    1 1 1  =>  0 0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input/output, integer ( kind = 4 ) BVEC(N), on output, the successor
!    to the input vector.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) bvec(n)
  integer ( kind = 4 ) i

  do i = 1, n

    if ( bvec(i) == 1 ) then
      bvec(i) = 0
    else
      bvec(i) = 1
      exit
    end if

  end do

  return
end
subroutine bvec_add ( n, bvec1, bvec2, bvec3 )

!*****************************************************************************80
!
!! BVEC_ADD adds two (signed) binary vectors.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Example:
!
!    N = 5
!
!      BVEC1       +   BVEC2       =   BVEC3
!
!    ( 0 0 0 0 1 ) + ( 0 0 0 1 1 ) = ( 0 0 1 0 0 )
!
!              1   +           3   =           4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) BVEC1(N), BVEC2(N), the vectors to be added.
!
!    Output, integer ( kind = 4 ) BVEC3(N), the sum of the two input vectors.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 2
  integer ( kind = 4 ) bvec1(n)
  integer ( kind = 4 ) bvec2(n)
  integer ( kind = 4 ) bvec3(n)
  integer ( kind = 4 ) i
  logical overflow

  overflow = .false.

  bvec3(1:n) = bvec1(1:n) + bvec2(1:n)

  do i = n, 1, -1

    do while ( base <= bvec3(i) )

      bvec3(i) = bvec3(i) - base

      if ( 1 < i ) then
        bvec3(i-1) = bvec3(i-1) + 1
      else
        overflow = .true.
      end if

    end do

  end do

  return
end
subroutine bvec_and ( n, bvec1, bvec2, bvec3 )

!*****************************************************************************80
!
!! BVEC_AND computes the AND of two binary vectors.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) BVEC1(N), BVEC2(N), the binary vectors.
!
!    Input, integer ( kind = 4 ) BVEC3(N), the AND of the two vectors.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) bvec1(n)
  integer ( kind = 4 ) bvec2(n)
  integer ( kind = 4 ) bvec3(n)

  bvec3(1:n) = min ( bvec1(1:n), bvec2(1:n) )

  return
end
subroutine bvec_check ( n, bvec, ierror )

!*****************************************************************************80
!
!! BVEC_CHECK checks a binary vector.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!    The only check made is that the entries are all 0 or 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) BVEC(N), the vector to be checked.
!
!    Output, integer ( kind = 4 ) IERROR, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 2
  integer ( kind = 4 ) bvec(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror

  ierror = 0

  do i = 1, n
    if ( bvec(i) < 0 .or. base <= bvec(i) ) then
      ierror = i
      return
    end if
  end do

  return
end
subroutine bvec_complement2 ( n, bvec1, bvec2 )

!*****************************************************************************80
!
!! BVEC_COMPLEMENT2 computes the two's complement of a binary vector.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) BVEC1(N), the vector to be complemented.
!
!    Output, integer ( kind = 4 ) BVEC2(N), the two's complemented vector.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 2
  integer ( kind = 4 ) bvec1(n)
  integer ( kind = 4 ) bvec2(n)
  integer ( kind = 4 ) bvec3(n)
  integer ( kind = 4 ) bvec4(n)

  bvec3(1:n) = ( base - 1 ) - bvec1(1:n)

  bvec4(1:n-1) = 0
  bvec4(n) = 1

  call bvec_add ( n, bvec3, bvec4, bvec2 )

  return
end
subroutine bvec_mul ( n, bvec1, bvec2, bvec3 )

!*****************************************************************************80
!
!! BVEC_MUL computes the product of two binary vectors.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!    Since the user may want to make calls like
!
!      call bvec_mul ( n, bvec1, bvec1, bvec3 )
!    or even
!      call bvec_mul ( n, bvec1, bvec1, bvec1 )
!
!    we need to copy the arguments, work on them, and then copy out the result.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) BVEC1(N), BVEC2(N), the vectors to
!    be multiplied.
!
!    Output, integer ( kind = 4 ) BVEC3(N), the product of the two
!    input vectors.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 2
  integer ( kind = 4 ) carry
  integer ( kind = 4 ) bvec1(n)
  integer ( kind = 4 ) bvec2(n)
  integer ( kind = 4 ) bvec3(n)
  integer ( kind = 4 ) bveca(n)
  integer ( kind = 4 ) bvecb(n)
  integer ( kind = 4 ) bvecc(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) product_sign
!
!  Copy the input.
!
  bveca(1:n) = bvec1(1:n)
  bvecb(1:n) = bvec2(1:n)
!
!  Record the sign of the product.
!  Make the factors positive.
!
  product_sign = 1

  if ( bveca(1) /= 0 ) then
    product_sign = - product_sign
    call bvec_complement2 ( n, bveca, bveca )
  end if

  if ( bvecb(1) /= 0 ) then
    product_sign = - product_sign
    call bvec_complement2 ( n, bvecb, bvecb )
  end if

  bvecc(1:n) = 0
!
!  Multiply.
!
  do i = 2, n
    bvecc(2:n+2-i) = bvecc(2:n+2-i) + bveca(n+2-i) * bvecb(i:n)
  end do
!
!  Take care of carries.
!
  do i = n, 2, -1

    carry = bvecc(i) / base
    bvecc(i) = bvecc(i) - carry * base
!
!  Unlike the case of BVEC_ADD, we do NOT allow carries into
!  the first position when multiplying.
!
    if ( 2 < i ) then
      bvecc(i-1) = bvecc(i-1) + carry
    end if

  end do
!
!  Take care of the sign of the product.
!
  if ( product_sign < 0 ) then
    call bvec_complement2 ( n, bvecc, bvecc )
  end if
!
!  Copy the output.
!
  bvec3(1:n) = bvecc(1:n)

  return
end
subroutine bvec_next ( n, bvec )

!*****************************************************************************80
!
!! BVEC_NEXT generates the next BVEC.
!
!  Discussion:
!
!    A BVEC is a binary vector, an N vector whose entries are 0 or 1.
!
!    The vectors are produced in the order:
!
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
!    20 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input/output, integer ( kind = 4 ) BVEC(N), on output, the successor to the
!    input vector.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) bvec(n)
  integer ( kind = 4 ) i

  do i = n, 1, -1

    if ( bvec(i) == 0 ) then
      bvec(i) = 1
      return
    end if

    bvec(i) = 0

  end do

  return
end
subroutine bvec_not ( n, bvec1, bvec2 )

!*****************************************************************************80
!
!! BVEC_NOT "negates" or takes the 1's complement of a binary vector.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) BVEC1(N), the vector to be negated.
!
!    Output, integer ( kind = 4 ) BVEC2(N), the negated vector.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 2
  integer ( kind = 4 ) bvec1(n)
  integer ( kind = 4 ) bvec2(n)

  bvec2(1:n) = ( base - 1 ) - bvec1(1:n)

  return
end
subroutine bvec_or ( n, bvec1, bvec2, bvec3 )

!*****************************************************************************80
!
!! BVEC_OR computes the inclusive OR of two binary vectors.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) BVEC1(N), BVEC2(N), the binary vectors.
!
!    Input, integer ( kind = 4 ) BVEC3(N), the inclusive OR of the two vectors.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) bvec1(n)
  integer ( kind = 4 ) bvec2(n)
  integer ( kind = 4 ) bvec3(n)

  bvec3(1:n) = max ( bvec1(1:n), bvec2(1:n) )

  return
end
subroutine bvec_print ( n, bvec, title )

!*****************************************************************************80
!
!! BVEC_PRINT prints a BVEC, with an optional title.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) BVEC(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) bvec(n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do ilo = 1, n, 70
    ihi = min ( ilo + 70 - 1, n )
    write ( *, '(2x,70i1)' ) bvec(ilo:ihi)
  end do

  return
end
subroutine bvec_reverse ( n, bvec1, bvec2 )

!*****************************************************************************80
!
!! BVEC_REVERSE reverses a binary vector.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) BVEC1(N), the vector to be reversed.
!
!    Output, integer ( kind = 4 ) BVEC2(N), the reversed vector.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) bvec1(n)
  integer ( kind = 4 ) bvec2(n)

  bvec2(1:n) = bvec1(n:1:-1)

  return
end
subroutine bvec_sub ( n, bvec1, bvec2, bvec3 )

!*****************************************************************************80
!
!! BVEC_SUB subtracts two binary vectors.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Example:
!
!    N = 4
!
!    BVEC1         BVEC2         BVEC3
!    -------       -------       -------
!    0 1 0 0   -   0 0 0 1   =   0 0 1 1
!          4             1             3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) BVEC1(N), BVEC2(N), the vectors to
!    be subtracted.
!
!    Output, integer ( kind = 4 ) BVEC3(N), the value of BVEC1 - BVEC2.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) bvec1(n)
  integer ( kind = 4 ) bvec2(n)
  integer ( kind = 4 ) bvec3(n)
  integer ( kind = 4 ) bvec4(n)

  call bvec_complement2 ( n, bvec2, bvec4 )

  call bvec_add ( n, bvec1, bvec4, bvec3 )

  return
end
subroutine bvec_to_i4 ( n, bvec, i4 )

!*****************************************************************************80
!
!! BVEC_TO_I4 makes an integer from a (signed) binary vector.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Example:
!
!         BVEC   binary  I
!    ----------  -----  --
!    1  2  3  4
!    ----------
!    0  0  0  1       1  1
!    0  0  1  0      10  2
!    1  1  0  0    -100 -4
!    0  1  0  0     100  4
!    1  0  0  1    -111 -9
!    1  1  1  1      -0  0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input, integer ( kind = 4 ) BVEC(N), the binary representation.
!
!    Output, integer ( kind = 4 ) I4, the integer.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 2
  integer ( kind = 4 ) bvec(n)
  integer ( kind = 4 ) bvec2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_sign
  integer ( kind = 4 ) i4

  bvec2(1:n) = bvec(1:n)

  if ( bvec2(1) == base - 1 ) then
    i_sign = -1
    bvec2(1) = 0
    call bvec_complement2 ( n - 1, bvec2(2:n), bvec2(2:n) )
  else
    i_sign = 1
  end if

  i4 = 0
  do i = 2, n
    i4 = base * i4 + bvec2(i)
  end do

  i4 = i_sign * i4

  return
end
subroutine bvec_xor ( n, bvec1, bvec2, bvec3 )

!*****************************************************************************80
!
!! BVEC_XOR computes the exclusive OR of two binary vectors.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) BVEC1(N), BVEC2(N), the binary vectors
!    to be XOR'ed.
!
!    Input, integer ( kind = 4 ) BVEC3(N), the exclusive OR of the two vectors.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) bvec1(n)
  integer ( kind = 4 ) bvec2(n)
  integer ( kind = 4 ) bvec3(n)

  bvec3(1:n) = mod ( bvec1(1:n) + bvec2(1:n), 2 )

  return
end
subroutine catalan ( n, c )

!*****************************************************************************80
!
!! CATALAN computes the Catalan numbers, from C(0) to C(N).
!
!  Discussion:
!
!    The Catalan number C(N) counts:
!
!    1) the number of binary trees on N vertices;
!    2) the number of ordered trees on N+1 vertices;
!    3) the number of full binary trees on 2N+1 vertices;
!    4) the number of well formed sequences of 2N parentheses;
!    5) number of ways 2N ballots can be counted, in order,
!       with N positive and N negative, so that the running sum
!       is never negative;
!    6) the number of standard tables in a 2 by N rectangular Ferrers diagram;
!    7) the number of monotone functions from [1..N} to [1..N} which
!       satisfy f(i) <= i for all i;
!    8) the number of ways to triangulate a polygon with N+2 vertices.
!
!    When N = 3, here are the 5 well formed parentheses sets:
!
!      ()()()
!      ()(())
!      (()())
!      (())()
!      ((()))
!
!  Example:
!
!     0     1
!     1     1
!     2     2
!     3     5
!     4    14
!     5    42
!     6   132
!     7   429
!     8  1430
!     9  4862
!    10 16796
!
!  Formula:
!
!    C(N) = (2*N)! / ( (N+1) * (N!) * (N!) )
!         = 1 / (N+1) * COMB ( 2N, N )
!         = 1 / (2N+1) * COMB ( 2N+1, N+1).
!
!  Recursion:
!
!    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
!    C(N) = sum ( 1 <= I <= N-1 ) C(I) * C(N-I)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer, 1986,
!    ISBN: 0387963472,
!    LC: QA164.S79.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of Catalan numbers desired.
!
!    Output, integer ( kind = 4 ) C(0:N), the Catalan numbers from C(0) to C(N).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) c(0:n)
  integer ( kind = 4 ) i

  if ( n < 0 ) then
    return
  end if

  c(0) = 1
!
!  The extra parentheses ensure that the integer division is
!  done AFTER the integer multiplication.
!
  do i = 1, n
    c(i) = ( c(i-1) * 2 * ( 2 * i - 1 ) ) / ( i + 1 )
  end do

  return
end
subroutine catalan_row_next ( ido, n, irow )

!*****************************************************************************80
!
!! CATALAN_ROW_NEXT computes row N of Catalan's triangle.
!
!  Example:
!
!    I\J 0   1   2   3   4   5   6
!
!    0   1
!    1   1   1
!    2   1   2   2
!    3   1   3   5   5
!    4   1   4   9  14  14
!    5   1   5  14  28  42  42
!    6   1   6  20  48  90 132 132
!
!  Recursion:
!
!    C(0,0) = 1
!    C(I,0) = 1
!    C(I,J) = 0 for I < J
!    C(I,J) = C(I,J-1) + C(I-1,J)
!    C(I,I) is the I-th Catalan number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IDO, indicates whether this is a call for
!    the 'next' row of the triangle.
!    IDO = 0, this is a startup call.  Row N is desired, but
!    presumably this is a first call, or row N-1 was not computed
!    on the previous call.
!    IDO = 1, this is not the first call, and row N-1 was computed
!    on the previous call.  In this case, much work can be saved
!    by using the information from the previous values of IROW
!    to build the next values.
!
!    Input, integer ( kind = 4 ) N, the index of the row of the triangle
!    desired.
!
!    Input/output, integer ( kind = 4 ) IROW(0:N), the row of coefficients.
!    If IDO = 0, then IROW is not required to be set on input.
!    If IDO = 1, then IROW must be set on input to the value of
!    row N-1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) irow(0:n)
  integer ( kind = 4 ) j

  if ( n < 0 ) then
    return
  end if

  if ( ido == 0 ) then

    irow(0) = 1
    irow(1:n) = 0

    do i = 1, n

      irow(0) = 1

      do j = 1, i-1
        irow(j) = irow(j) + irow(j-1)
      end do

      irow(i) = irow(i-1)

    end do

  else

    irow(0) = 1

    do j = 1, n-1
      irow(j) = irow(j) + irow(j-1)
    end do

    if ( 1 <= n ) then
      irow(n) = irow(n-1)
    end if

  end if

  return
end
subroutine catalan_values ( n_data, n, c )

!*****************************************************************************80
!
!! CATALAN_VALUES returns some values of the Catalan numbers for testing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the order of the Catalan number.
!
!    Output, integer ( kind = 4 ) C, the value of the Catalan number.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 11

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( nmax ) :: c_vec = (/ &
    1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
     0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine cbt_traverse ( depth )

!*****************************************************************************80
!
!! CBT_TRAVERSE traverses a complete binary tree of given depth.
!
!  Discussion:
!
!    There will be 2^DEPTH terminal nodes of the complete binary tree.
!
!    This function traverses the tree, and prints out a binary code of 0's
!    and 1's each time it encounters a terminal node.  This results in a
!    printout of the binary digits from 0 to 2^DEPTH - 1.
!
!    The function is intended as a framework to be used to traverse a binary
!    tree.  Thus, in practice, a user would insert some action when a terminal
!    node is encountered.
!
!    Another use would occur when a combinatorial search is being made, for
!    example in a knapsack problem.  Each binary string then represents which
!    objects are to be included in the knapsack.  In that case, the traversal
!    could be speeded up by noticing cases where a nonterminal node has been
!    reached, but the knapsack is already full, in which case the only solution
!    uses none of the succeeding items, or overfull, in which case no solutions
!    exist that include this initial path segment.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DEPTH, the depth of the tree.
!
  implicit none

  integer ( kind = 4 ) depth

  integer ( kind = 4 ) b(depth)
  integer ( kind = 4 ) direction
  integer ( kind = 4 ) :: DOWNLEFT = 1
  integer ( kind = 4 ) k
  integer ( kind = 4 ) p
  integer ( kind = 4 ) :: UP = 3
  integer ( kind = 4 ) :: UPDOWNRIGHT = 2

  if ( depth < 1 ) then
    return
  end if

  b(1:depth) = 0
  p = 0
  direction = DOWNLEFT
  k = 0

  do
!
!  Try going in direction DOWNLEFT.
!
    if ( direction == DOWNLEFT ) then
      p = p + 1
      b(p) = 0
      if ( depth <= p ) then
        write ( *, '(2x,i4,2x,10i1)' ) k, b(1:depth)
        k = k + 1
        direction = UPDOWNRIGHT
      end if
    end if
!
!  Try going in direction UPDOWNRIGHT.
!
    if ( direction == UPDOWNRIGHT ) then
      b(p) = + 1
      if ( p < depth ) then
        direction = DOWNLEFT
      else
        write ( *, '(2x,i4,2x,10i1)' ) k, b(1:depth)
        k = k + 1
        direction = UP
      end if
    end if
!
!  Try going in direction UP.
!
    if ( direction == UP ) then
      p = p - 1
      if ( 1 <= p ) then
        if ( b(p) == 0 ) then
          direction = UPDOWNRIGHT
        end if
      else
        exit
      end if
    end if

  end do

  return
end
subroutine cfrac_to_rat ( n, a, p, q )

!*****************************************************************************80
!
!! CFRAC_TO_RAT converts a monic continued fraction to an ordinary fraction.
!
!  Discussion:
!
!    The routine is given the monic or "simple" continued fraction with
!    integer coefficients:
!
!      A(1) + 1 / ( A(2) + 1 / ( A(3) ... + 1 / A(N) ) )
!
!    and returns the N successive approximants P(I)/Q(I)
!    to the value of the rational number represented by the continued
!    fraction, with the value exactly equal to the final ratio P(N)/Q(N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi,
!    John Rice, Henry Thatcher, Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of continued fraction
!    coefficients.
!
!    Input, integer ( kind = 4 ) A(N), the continued fraction coefficients.
!
!    Output, integer ( kind = 4 ) P(N), Q(N), the N successive approximations
!    to the value of the continued fraction.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) q(n)

  do i = 1, n

    if ( i == 1 ) then
      p(i) = a(i) * 1 + 0
      q(i) = a(i) * 0 + 1
    else if ( i == 2 ) then
      p(i) = a(i) * p(i-1) + 1
      q(i) = a(i) * q(i-1) + 0
    else
      p(i) = a(i) * p(i-1) + p(i-2)
      q(i) = a(i) * q(i-1) + q(i-2)
    end if

  end do

  return
end
subroutine cfrac_to_rfrac ( m, g, h, p, q )

!*****************************************************************************80
!
!! CFRAC_TO_RFRAC converts polynomial fractions from continued to rational form.
!
!  Discussion:
!
!    The routine accepts a continued polynomial fraction:
!
!      G(1)     / ( H(1) +
!      G(2) * X / ( H(2) +
!      G(3) * X / ( H(3) + ...
!      G(M) * X / ( H(M) )...) ) )
!
!    and returns the equivalent rational polynomial fraction:
!
!      P(1) + P(2) * X + ... + P(L1) * X**(L1)
!      -------------------------------------------------------
!      Q(1) + Q(2) * X + ... + Q(L2) * X**(L2-1)
!
!    where
!
!      L1 = (M+1)/2
!      L2 = (M+2)/2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2004
!
!  Author:
!
!    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson,
!    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi,
!    John Rice, Henry Thatcher, Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of continued fraction
!    polynomial coefficients.
!
!    Input, real ( kind = 8 ) G(M), H(M), the continued polynomial
!    fraction coefficients.
!
!    Output, real ( kind = 8 ) P((M+1)/2), Q((M+2)/2), the rational
!    polynomial fraction coefficients.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,(m+2)/2)
  real ( kind = 8 ) g(m)
  real ( kind = 8 ) h(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p((m+1)/2)
  real ( kind = 8 ) q((m+2)/2)

  if ( m == 1 ) then
    p(1) = g(1)
    q(1) = h(1)
    return
  end if

  do i = 1, m
    do j = 1, (m+2)/2
      a(i,j) = 0.0D+00
    end do
  end do
!
!  Solve for P's.
!
  a(1,1) = g(1)
  a(2,1) = g(1) * h(2)

  do i = 3, m
    a(i,1) = h(i) * a(i-1,1)
    do j = 2, (i+1)/2
      a(i,j) = h(i) * a(i-1,j) + g(i) * a(i-2,j-1)
    end do
  end do

  do j = 1, (m+1)/2
    p(j) = a(m,j)
  end do
!
!  Solve for Q's.
!
  a(1,1) = h(1)
  a(2,1) = h(1) * h(2)
  a(2,2) = g(2)

  do i = 3, m
    a(i,1) = h(i) * a(i-1,1)
    do j = 2, (i+2) / 2
      a(i,j) = h(i) * a(i-1,j) + g(i) * a(i-2,j-1)
    end do
  end do

  do j = 1, (m+2)/2
    q(j) = a(m,j)
  end do

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

  character              c
  integer   ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
subroutine change_greedy ( total, coin_num, coin_value, change_num, change )

!*****************************************************************************80
!
!! CHANGE_GREEDY makes change for a given total using the biggest coins first.
!
!  Discussion:
!
!    The algorithm is simply to use as many of the largest coin first,
!    then the next largest, and so on.
!
!    It is assumed that there is always a coin of value 1.  The
!    algorithm will otherwise fail!
!
!  Example:
!
!    Total = 17
!    COIN_NUM = 3
!    COIN_VALUE = (/ 1, 5, 10 /)
!
!
!    #  CHANGE              COIN_VALUE(CHANGE)
!
!    4  3 2 1 1             10 5 1 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TOTAL, the total for which change is to be
!    made.
!
!    Input, integer ( kind = 4 ) COIN_NUM, the number of types of coins.
!
!    Input, integer ( kind = 4 ) COIN_VALUE(COIN_NUM), the value of each coin.
!    The values should be in ascending order, and if they are not,
!    they will be sorted.
!
!    Output, integer ( kind = 4 ) CHANGE_NUM, the number of coins given in
!    change.
!
!    Output, integer ( kind = 4 ) CHANGE(TOTAL), the indices of the coins
!    will be in entries 1 through CHANGE_NUM.
!
  implicit none

  integer ( kind = 4 ) coin_num
  integer ( kind = 4 ) total

  integer ( kind = 4 ) change(total)
  integer ( kind = 4 ) change_num
  integer ( kind = 4 ) coin_value(coin_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) total_copy

  change_num = 0
!
!  Find the largest coin smaller than the total.
!
  j = coin_num

  do while ( 0 < j )
    if ( coin_value(j) <= total ) then
      exit
    end if
    j = j - 1
  end do

  if ( j <= 0 ) then
    return
  end if
!
!  Subtract the current coin from the total.
!  Once that coin is too big, use the next coin.
!
  total_copy = total

  do while ( 0 < total_copy )

    if ( coin_value(j) <= total_copy ) then

      total_copy = total_copy - coin_value(j)
      change_num = change_num + 1
      change(change_num) = j

    else

      j = j - 1
      if ( j <= 0 ) then
        exit
      end if

    end if

  end do

  return
end
subroutine change_next ( total, coin_num, coin_value, change_num, change, &
  done  )

!*****************************************************************************80
!
!! CHANGE_NEXT computes the next set of change for a given sum.
!
!  Example:
!
!    Total = 17
!    COIN_NUM = 3
!    COIN_VALUE = (/ 1, 5, 10 /)
!
!
!        #  CHANGE              COIN_VALUE(CHANGE)
!
!    1   4  3 2 1 1             10 5 1 1
!    2   8  3 1 1 1 1 1 1 1     10 1 1 1 1 1 1 1
!    3   5  2 2 2 1 1            5 5 5 1 1
!    4   9  2 2 1 1 1 1 1 1 1    5 5 1 1 1 1 1 1 1
!    5  13  2 1 1 1 1 1 1 1 1 1  5 1 1 1 1 1 1 1 1 1
!           1 1 1                1 1 1
!    6  17  1 1 1 1 1 1 1 1 1 1  1 1 1 1 1 1 1 1 1 1 1
!           1 1 1 1 1 1 1        1 1 1 1 1 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TOTAL, the total for which change is to be
!    made.
!
!    Input, integer ( kind = 4 ) COIN_NUM, the number of types of coins.
!
!    Input, integer ( kind = 4 ) COIN_VALUE(COIN_NUM), the value of each coin.
!    The values must be in ascending order.
!
!    Input/output, integer ( kind = 4 ) CHANGE_NUM, the number of coins given
!    in change for this form of the change.
!
!    Input/output, integer ( kind = 4 ) CHANGE(CHANGE_NUM), the indices of the
!    coins.  The user must dimension this array to have dimension TOTAL!
!
!    Input/output, logical DONE.  The user sets DONE = TRUE on
!    first call to tell the routine this is the beginning of a computation.
!    The program resets DONE to FALSE and it stays that way until
!    the last possible change combination is made, at which point the
!    program sets DONE to TRUE again.
!
  implicit none

  integer ( kind = 4 ) coin_num
  integer ( kind = 4 ) total

  integer ( kind = 4 ) change(total)
  integer ( kind = 4 ) change_num
  integer ( kind = 4 ) change_num2
  integer ( kind = 4 ) coin_num2
  integer ( kind = 4 ) coin_value(coin_num)
  logical done
  integer ( kind = 4 ) i
  logical i4vec_ascends
  integer ( kind = 4 ) last
  integer ( kind = 4 ) total2

  if ( done ) then
!
!  Make sure the coin values are sorted.
!
    if ( .not. i4vec_ascends ( coin_num, coin_value ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CHANGE_NEXT - Fatal error!'
      write ( *, '(a)' ) '  The COIN_VALUE array is not in ascending order.'
      stop
    end if
!
!  Start with the greedy change.
!
    call change_greedy ( total, coin_num, coin_value, change_num, change )
!
!  In a few cases, like change for 4 cents, we're done after the first call.
!
    if ( change_num == total ) then
      done = .true.
    else
      done = .false.
    end if

  else
!
!  Find the last location in the input change which is NOT a penny.
!
    last = 0

    do i = change_num, 1, -1
      if ( change(i) /= 1 ) then
        last = i
        exit
      end if
    end do
!
!  If that location is still 0, an error was made.
!
    if ( last == 0 ) then
      done = .true.
      return
    end if
!
!  Sum the entries from that point to the end.
!
    total2 = sum ( coin_value ( change(last:change_num) ) )
!
!  Make greedy change for the partial sum using coins smaller than that one.
!
    coin_num2 = change(last) - 1

    call change_greedy ( total2, coin_num2, coin_value, change_num2, &
      change(last:total) )

    change_num = ( last - 1 ) + change_num2

  end if

  return
end
subroutine chinese_check ( n, m, ierror )

!*****************************************************************************80
!
!! CHINESE_CHECK checks the Chinese remainder moduluses.
!
!  Discussion:
!
!    For a Chinese remainder representation, the moduluses M(I) must
!    be positive and pairwise prime.  Also, in case this is not obvious,
!    no more than one of the moduluses may be 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of moduluses.
!
!    Input, integer ( kind = 4 ) M(N), the moduluses.  These should be positive
!    and pairwise prime.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error was detected.
!    nonzero, an error was detected.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  logical i4vec_pairwise_prime
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m(n)

  ierror = 0
!
!  Do not allow nonpositive entries.
!
  if ( any ( m(1:n) <= 0 ) ) then
    ierror = 1
    return
  end if
!
!  Allow one entry to be 1, but not two entries.
!
  do i = 1, n
    do j = i+1, n
      if ( m(i) == 1 .and. m(j) == 1 ) then
        ierror = 2
        return
      end if
    end do
  end do
!
!  Now check pairwise primeness.
!
  if ( .not. i4vec_pairwise_prime ( n, m ) ) then
    ierror = 3
    return
  end if

  return
end
subroutine chinese_to_i4 ( n, m, r, j )

!*****************************************************************************80
!
!! CHINESE_TO_I4 converts a set of Chinese remainders to an equivalent integer.
!
!  Discussion:
!
!    Given a set of N pairwise prime, positive moduluses M(I), and
!    a corresponding set of remainders R(I), this routine finds an
!    integer J such that, for all I,
!
!      J = R(I) mod M(I)
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
!    Input, integer ( kind = 4 ) N, the number of moduluses.
!
!    Input, integer ( kind = 4 ) M(N), the moduluses.  These should be
!    positive and pairwise prime.
!
!    Input, integer ( kind = 4 ) R(N), the Chinese remainder representation
!    of the integer.
!
!    Output, integer ( kind = 4 ) J, the corresponding integer.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b(n)
  integer ( kind = 4 ) big_m
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m(n)
  integer ( kind = 4 ) r(n)

  call chinese_check ( n, m, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHINESE_TO_I4 - Fatal error!'
    write ( *, '(a)' ) '  The moduluses are not legal.'
    stop
  end if
!
!  Set BIG_M.
!
  big_m = product ( m(1:n) )
!
!  Solve BIG_M / M(I) * B(I) = 1, mod M(I)
!
  do i = 1, n
    a = big_m / m(i)
    c = 1
    call congruence ( a, m(i), c, ierror, b(i) )
  end do
!
!  Set J = sum ( 1 <= I <= N ) ( R(I) * B(I) * BIG_M / M(I) ) mod M
!
  j = 0
  do i = 1, n
    j = mod ( j + r(i) * b(i) * ( big_m / m(i) ), big_m )
  end do

  return
end
subroutine comb_next ( n, k, a, done )

!*****************************************************************************80
!
!! COMB_NEXT computes combinations of K things out of N.
!
!  Discussion:
!
!    The combinations are computed one at a time, in lexicographical order.
!
!    10 April 1009: Thanks to "edA-qa mort-ora-y" for supplying a
!    correction to this code!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Charles Mifsud,
!    Algorithm 154:
!    Combination in Lexicographic Order,
!    Communications of the ACM,
!    March 1963.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the total number of things.
!
!    Input, integer ( kind = 4 ) K, the number of things in each combination.
!
!    Input/output, integer ( kind = 4 ) A(K), contains the list of elements in
!    the current combination.
!
!    Input/output, logical DONE.  On first call, set DONE to TRUE,
!    and thereafter, its input value should be the output value from
!    the previous call.  The output value will normally be FALSE,
!    indicating that there are further combinations that can be
!    returned.  When DONE is returned TRUE, the sequence is exhausted.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  logical done
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n

  if ( done ) then

    if ( k <= 0 ) then
      return
    end if

    call i4vec_indicator ( k, a )

    done = .false.

  else

    if ( a(k) < n ) then
      a(k) = a(k) + 1
      return
    end if

    do i = k, 2, -1

      if ( a(i-1) < n-k+i-1 ) then

        a(i-1) = a(i-1) + 1

        do j = i, k
          a(j) = a(i-1) + j - ( i-1 )
        end do

        return

      end if

    end do

    done = .true.

  end if

  return
end
subroutine comb_row ( ido, n, irow )

!*****************************************************************************80
!
!! COMB_ROW computes row N of Pascal's triangle.
!
!  Discussion:
!
!    Row N contains the combinatorial coefficients
!
!      C(N,0), C(N,1), C(N,2), ... C(N,N)
!
!    The sum of the elements of row N is equal to 2**N.
!
!    The formula is
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  First terms:
!
!     N K:0  1   2   3   4   5   6   7  8  9 10
!
!     0   1
!     1   1  1
!     2   1  2   1
!     3   1  3   3   1
!     4   1  4   6   4   1
!     5   1  5  10  10   5   1
!     6   1  6  15  20  15   6   1
!     7   1  7  21  35  35  21   7   1
!     8   1  8  28  56  70  56  28   8  1
!     9   1  9  36  84 126 126  84  36  9  1
!    10   1 10  45 120 210 252 210 120 45 10  1
!
!  Recursion:
!
!    C(N,K) = C(N-1,K-1)+C(N-1,K)
!
!  Special values:
!
!    C(N,0) = C(N,N) = 1
!    C(N,1) = C(N,N-1) = N
!    C(N,N-2) = sum ( 1 <= I <= N ) N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IDO, indicates whether this is a call for
!    the 'next' row of the triangle.
!
!    0 means this is a startup call.  Row N is desired, but
!    presumably this is a first call, or row N-1 was not computed
!    on the previous call.
!
!    1 means this is not the first call, and row N-1 was computed
!    on the previous call.  In this case, much work can be saved
!    by using the information from the previous values of IROW
!    to build the next values.
!
!    Input, integer ( kind = 4 ) N, the row of the triangle desired.  The
!    triangle begins with row 0.
!
!    Output, integer ( kind = 4 ) IROW(N+1), the row of coefficients.
!    IROW(I) = C(N,I-1).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) irow(n+1)
  integer ( kind = 4 ) j

  if ( n < 0 ) then
    return
  end if

  if ( ido == 1 ) then

    do i = n, 2, -1
      irow(i) = irow(i) + irow(i-1)
    end do

    irow(n+1) = 1

  else

    irow(1) = 1
    irow(2:n+1) = 0

    do j = 1, n
      do i = j+1, 2, -1
        irow(i) = irow(i) + irow(i-1)
      end do
    end do

  end if

  return
end
subroutine comb_unrank ( m, n, rank, a )

!*****************************************************************************80
!
!! COMB_UNRANK returns the RANK-th combination of N things out of M.
!
!  Discussion:
!
!    Going from a rank to a thing is called "unranking".
!
!    The combinations are ordered lexically.
!
!    Lexical order can be illustrated for the general case of N and M as
!    follows:
!
!    1:       1,     2,     3,     ..., N-2, N-1, N
!    2:       1,     2,     3,     ..., N-2, N-1, N+1
!    3:       1,     2,     3,     ..., N-2, N-1, N+2
!    ...
!    M-N+1:   1,     2,     3,     ..., N-2, N-1, M
!    M-N+2:   1,     2,     3,     ..., N-2, N,   N+1
!    M-N+3:   1,     2,     3,     ..., N-2, N,   N+2
!    ...
!    LAST-2:  M-N,   M-N+1, M-N+3, ..., M-2, M-1, M
!    LAST-1:  M-N,   M-N+2, M-N+3, ..., M-2, M-1, M
!    LAST:    M-N+1, M-N+2, M-N+3, ..., M-2, M-1, M
!
!    There are a total of M!/(N!*(M-N)!) combinations of M
!    things taken N at a time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Bill Buckles, Matthew Lybanon,
!    Algorithm 515:
!    Generation of a Vector from the Lexicographical Index,
!    ACM Transactions on Mathematical Software,
!    Volume 3, Number 2, pages 180-182, June 1977.
!
!    David Crouse,
!    Remark on Algorithm 515,
!    ACM Transactions on Mathematical Software,
!    Volume 33, Number 2, Article 15, June 2007.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the size of the set.
!
!    Input, integer ( kind = 4 ) N, the number of things in the combination.
!    N must be greater than 0, and no greater than M.
!
!    Input, integer ( kind = 4 ) RANK, the lexicographical rank of the
!    combination sought.  RANK must be at least 1, and no greater
!    than M!/(N!*(M-N)!).
!
!    Output, integer ( kind = 4 ) A(N), the combination.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) rank

  if ( n < 0 ) then
    return
  else if ( n == 0 ) then
    a(1) = rank
    return
  end if
!
!  Initialize the lower bound index.
!
  k = 0
!
!  Select elements in ascending order.
!
  do i = 1, n - 1
!
!  Set the lower bound element number for next element value.
!
    a(i) = 0

    if ( 1 < i ) then
      a(i) = a(i-1)
    end if
!
!  Check each element value.
!
    do

      a(i) = a(i) + 1
      j = i4_choose ( m-a(i), n-i )
      k = k + j

      if ( rank <= k ) then
        exit
      end if

    end do

    k = k - j

  end do

  a(n) = a(n-1) + rank - k

  return
end
subroutine comp_enum ( n, k, number )

!*****************************************************************************80
!
!! COMP_ENUM returns the number of compositions of the integer N into K parts.
!
!  Discussion:
!
!    A composition of the integer N into K parts is an ordered sequence
!    of K nonnegative integers which sum to N.  The compositions (1,2,1)
!    and (1,1,2) are considered to be distinct.
!
!    The 28 compositions of 6 into three parts are:
!
!      6 0 0,  5 1 0,  5 0 1,  4 2 0,  4 1 1,  4 0 2,
!      3 3 0,  3 2 1,  3 1 2,  3 0 3,  2 4 0,  2 3 1,
!      2 2 2,  2 1 3,  2 0 4,  1 5 0,  1 4 1,  1 3 2,
!      1 2 3,  1 1 4,  1 0 5,  0 6 0,  0 5 1,  0 4 2,
!      0 3 3,  0 2 4,  0 1 5,  0 0 6.
!
!    The formula for the number of compositions of N into K parts is
!
!      Number = ( N + K - 1 )! / ( N! * ( K - 1 )! )
!
!    Describe the composition using N '1's and K-1 dividing lines '|'.
!    The number of distinct permutations of these symbols is the number
!    of compositions.  This is equal to the number of permutations of
!    N+K-1 things, with N identical of one kind and K-1 identical of another.
!
!    Thus, for the above example, we have:
!
!      Number = ( 6 + 3 - 1 )! / ( 6! * (3-1)! ) = 28
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
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
!    Output, integer ( kind = 4 ) NUMBER, the number of compositions of N
!    into K parts.
!
  implicit none

  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) number

  number = i4_choose ( n + k - 1, n )

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
!    Input/output, integer ( kind = 4 ) H, T, two internal parameters needed
!    for the computation.  The user should allocate space for these in the
!    calling program, include them in the calling sequence, but never alter
!    them!
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
!
!  If the first entry A(1) is positive, then set H to zero,
!  so that when we increment H, it points to A(1); we will decrement A(1) by 1
!  and increment A(2).
!
    if ( 1 < t ) then
      h = 0
    end if
!
!  Otherwise, A(1) is 0.  Then by H + 1 is the entry we incremented last time.
!  Set H = H + 1, zero A(H), adding all but one of its value to A(1),
!  and incrementing A(H+1) by 1.
!
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
subroutine comp_random ( n, k, seed, a )

!*****************************************************************************80
!
!! COMP_RANDOM selects a random composition of the integer N into K parts.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2003
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
!    Input, integer ( kind = 4 ) N, the integer to be decomposed.
!
!    Input, integer ( kind = 4 ) K, the number of parts in the composition.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, integer ( kind = 4 ) A(K), the parts of the composition.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed

  call ksub_random ( n + k - 1, k - 1, seed, a )

  a(k) = n + k
  l = 0

  do i = 1, k
    m = a(i)
    a(i) = a(i) - l - 1
    l = m
  end do

  return
end
subroutine compnz_enum ( n, k, number )

!*****************************************************************************80
!
!! COMPNZ_ENUM returns the number of nonzero compositions of the N into K parts.
!
!  Discussion:
!
!    A composition of the integer N into K nonzero parts is an ordered sequence
!    of K positive integers which sum to N.  The compositions (1,2,1)
!    and (1,1,2) are considered to be distinct.
!
!    The routine computes one composition on each call until there are no more.
!    For instance, one composition of 6 into 3 parts is 3+2+1, another would
!    be 4+1+1 but 5+1+0 is not allowed since it includes a zero part.
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
!    The 10 compositions of 6 into three nonzero parts are:
!
!      4 1 1,  3 2 1,  3 1 2,  2 3 1,  2 2 2,  2 1 3,
!      1 4 1,  1 3 2,  1 2 3,  1 1 4.
!
!    The formula for the number of compositions of N into K nonzero
!    parts is
!
!      Number = ( N - 1 )! / ( ( N - K )! * ( K - 1 )! )
!
!    (Describe the composition using N-K '1's and K-1 dividing lines '|'.
!    The number of distinct permutations of these symbols is the number
!    of compositions into nonzero parts.  This is equal to the number of
!    permutations of  N-1 things, with N-K identical of one kind
!    and K-1 identical of another.)
!
!    Thus, for the above example, we have:
!
!      Number = ( 6 - 1 )! / ( ( 6 - 3 )! * ( 3 - 1 )! ) = 10
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
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
!    Output, integer ( kind = 4 ) NUMBER, the number of compositions of N into
!    K nonzero parts.
!
  implicit none

  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) number

  number = i4_choose ( n - 1, n - k )

  return
end
subroutine compnz_next ( n, k, a, more )

!*****************************************************************************80
!
!! COMPNZ_NEXT computes the compositions of the integer N into K nonzero parts.
!
!  Discussion:
!
!    A composition of the integer N into K nonzero parts is an ordered sequence
!    of K positive integers which sum to N.  The compositions (1,2,1)
!    and (1,1,2) are considered to be distinct.
!
!    The routine computes one composition on each call until there are no more.
!    For instance, one composition of 6 into 3 parts is 3+2+1, another would
!    be 4+1+1 but 5+1+0 is not allowed since it includes a zero part.
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
!    The 10 compositions of 6 into three nonzero parts are:
!
!      4 1 1,  3 2 1,  3 1 2,  2 3 1,  2 2 2,  2 1 3,
!      1 4 1,  1 3 2,  1 2 3,  1 1 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 December 2005
!
!  Author:
!
!    John Burkardt
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
!    K must be less than or equal to N.
!
!    Input/output, integer ( kind = 4 ) A(K), the parts of the composition.
!
!    Input/output, logical MORE, set by the user to start the computation,
!    and by the routine to terminate it.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ), save :: h = 0
  logical more
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: t = 0
!
!  We use the trick of computing ordinary compositions of (N-K)
!  into K parts, and adding 1 to each part.
!
  if ( n < k ) then
    more = .false.
    a(1:k) = -1
    return
  end if
!
!  The first computation.
!
  if ( .not. more ) then

    t = n - k
    h = 0
    a(1) = n - k
    a(2:k) = 0
!
!  The next computation.
!
  else

    a(1:k) = a(1:k) - 1

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
  more = ( a(k) /= ( n - k ) )

  a(1:k) = a(1:k) + 1

  return
end
subroutine compnz_random ( n, k, seed, a )

!*****************************************************************************80
!
!! COMPNZ_RANDOM selects a random composition of N into K nonzero parts.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 December 2005
!
!  Author:
!
!    John Burkardt
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
!    Input, integer ( kind = 4 ) N, the integer to be decomposed.
!
!    Input, integer ( kind = 4 ) K, the number of parts in the composition.
!    K must be no greater than N.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, integer ( kind = 4 ) A(K), the parts of the composition.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed

  if ( n < k ) then
    a(1:k) = -1
    return
  end if

  call ksub_random ( n - 1, k - 1, seed, a )

  a(k) = n
  l = 0

  do i = 1, k
    m = a(i)
    a(i) = a(i) - l - 1
    l = m
  end do

  a(1:k) = a(1:k) + 1

  return
end
subroutine congruence ( a, b, c, ierror, x )

!*****************************************************************************80
!
!! CONGRUENCE solves a congruence of the form A * X = C ( mod B ).
!
!  Discussion:
!
!    A, B and C are given integers.  The equation is solvable if and only
!    if the greatest common divisor of A and B also divides C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein, editor,
!    CRC Concise Encylopedia of Mathematics,
!    CRC Press, 1998, page 446.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, C, the coefficients of the Diophantine
!    equation.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error, X was computed.
!    1, A = B = 0, C is nonzero.
!    2, A = 0, B and C nonzero, but C is not a multiple of B.
!    3, A nonzero, B zero, C nonzero, but C is not a multiple of A.
!    4, A, B, C nonzero, but GCD of A and B does not divide C.
!    5, algorithm ran out of internal space.
!
!    Output, integer ( kind = 4 ) X, the solution of the Diophantine equation.
!    X will be between 0 and B-1.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 100

  integer ( kind = 4 ) a
  integer ( kind = 4 ) a_copy
  integer ( kind = 4 ) a_mag
  integer ( kind = 4 ) a_sign
  integer ( kind = 4 ) b
  integer ( kind = 4 ) b_copy
  integer ( kind = 4 ) b_mag
  integer ( kind = 4 ) b_sign
  integer ( kind = 4 ) c
  integer ( kind = 4 ) c_copy
  integer ( kind = 4 ) g
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) q(nmax)
  logical swap
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
  integer ( kind = 4 ) z
!
!  Defaults for output parameters.
!
  ierror = 0
  x = 0
  y = 0
!
!  Special cases.
!
  if ( a == 0 .and. b == 0 .and. c == 0 ) then
    x = 0
    return
  else if ( a == 0 .and. b == 0 .and. c /= 0 ) then
    ierror = 1
    x = 0
    return
  else if ( a == 0 .and. b /= 0 .and. c == 0 ) then
    x = 0
    return
  else if ( a == 0 .and. b /= 0 .and. c /= 0 ) then
    x = 0
    if ( mod ( c, b ) /= 0 ) then
      ierror = 2
    end if
    return
  else if ( a /= 0 .and. b == 0 .and. c == 0 ) then
    x = 0
    return
  else if ( a /= 0 .and. b == 0 .and. c /= 0 ) then
    x = c / a
    if ( mod ( c, a ) /= 0 ) then
      ierror = 3
    end if
    return
  else if ( a /= 0 .and. b /= 0 .and. c == 0 ) then
    g = i4_gcd ( a, b )
    x = b / g
    return
  end if
!
!  Now handle the "general" case: A, B and C are nonzero.
!
!  Step 1: Compute the GCD of A and B, which must also divide C.
!
  g = i4_gcd ( a, b )

  if ( mod ( c, g ) /= 0 ) then
    ierror = 4
    return
  end if

  a_copy = a / g
  b_copy = b / g
  c_copy = c / g
!
!  Step 2: Split A and B into sign and magnitude.
!
  a_mag = abs ( a_copy )
  a_sign = sign ( 1, a_copy )
  b_mag = abs ( b_copy )
  b_sign = sign ( 1, b_copy )
!
!  Another special case, A_MAG = 1 or B_MAG = 1.
!
  if ( a_mag == 1 ) then
    x = a_sign * c_copy
    return
  else if ( b_mag == 1 ) then
    x = 0
    return
  end if
!
!  Step 3: Produce the Euclidean remainder sequence.
!
  if ( b_mag <= a_mag ) then

    swap = .false.
    q(1) = a_mag
    q(2) = b_mag

  else

    swap = .true.
    q(1) = b_mag
    q(2) = a_mag

  end if

  n = 3

  do

    q(n) = mod ( q(n-2), q(n-1) )

    if ( q(n) == 1 ) then
      exit
    end if

    n = n + 1

    if ( nmax < n ) then
      ierror = 5
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CONGRUENCE - Fatal error!'
      write ( *, '(a)' ) '  Exceeded number of iterations.'
      stop
    end if

  end do
!
!  Step 4: Now go backwards to solve X * A_MAG + Y * B_MAG = 1.
!
  y = 0
  do k = n, 2, -1
    x = y
    y = ( 1 - x * q(k-1) ) / q(k)
  end do
!
!  Step 5: Undo the swapping.
!
  if ( swap ) then
    z = x
    x = y
    y = z
  end if
!
!  Step 6: Now apply signs to X and Y so that X * A + Y * B = 1.
!
  x = x * a_sign
!
!  Step 7: Multiply by C, so that X * A + Y * B = C.
!
  x = x * c_copy
!
!  Step 8: Now force 0 <= X < B.
!
  x = mod ( x, b )
!
!  Force positivity.
!
  if ( x < 0 ) then
    x = x + b
  end if

  return
end
subroutine count_pose_random ( seed, blocks, goal )

!*****************************************************************************80
!
!! COUNT_POSE_RANDOM poses a problem for the game "The Count is Good"
!
!  Discussion:
!
!    The French television show "The Count is Good" has a game that goes
!    as follows:
!
!      A number is chosen at random between 100 and 999.  This is the GOAL.
!
!      Six numbers are randomly chosen from the set 1, 2, 3, 4, 5, 6, 7, 8,
!      9, 10, 25, 50, 75, 100.  These numbers are the BLOCKS.
!
!      The player must construct a formula, using some or all of the blocks,
!      (but not more than once), and the operations of addition, subtraction,
!      multiplication and division.  Parentheses should be used to remove
!      all ambiguity.  However, it is forbidden to use subtraction in a
!      way that produces a negative result, and all division must come out
!      exactly, with no remainder.
!
!    This routine poses a sample problem from the show.  The point is,
!    to determine how to write a program that can solve such a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Raymond Seroul,
!    Programming for Mathematicians,
!    Springer Verlag, 2000, pages 355-357.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, integer ( kind = 4 ) BLOCKS(6), the six numbers available for
!    the formula.
!
!    Output, integer ( kind = 4 ) GOAL, the goal number.
!
  implicit none

  integer ( kind = 4 ) blocks(6)
  integer ( kind = 4 ) goal
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) ind(6)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), dimension ( 14 ) :: stuff = &
    (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 25, 50, 75, 100 /)

  goal = i4_uniform ( 100, 999, seed )

  call ksub_random ( 14, 6, seed, ind )

  blocks(1:6) = stuff(ind(1:6))

  return
end
subroutine debruijn ( m, n, string )

!*****************************************************************************80
!
!! DEBRUIJN constructs a de Bruijn sequence.
!
!  Discussion:
!
!    Suppose we have an alphabet of M letters, and we are interested in
!    all possible strings of length N.  If M = 2 and N = 3, then we are
!    interested in the M**N strings:
!
!      000
!      001
!      010
!      011
!      100
!      101
!      110
!      111
!
!    Now, instead of making a list like this, we prefer, if possible, to
!    write a string of letters, such that every consecutive sequence of
!    N letters is one of the strings, and every string occurs once, if
!    we allow wraparound.
!
!    For the above example, a suitable sequence would be the 8 characters:
!
!      00011101(00...
!
!    where we have suggested the wraparound feature by repeating the first
!    two characters at the end.
!
!    Such a sequence is called a de Bruijn sequence.  It can easily be
!    constructed by considering a directed graph, whose nodes are all
!    M**(N-1) strings of length N-1.  A node I has a directed edge to
!    node J (labeled with character K) if the string at node J can
!    be constructed by beheading the string at node I and adding character K.
!
!    In this setting, a de Bruijn sequence is simply an Eulerian circuit
!    of the directed graph, with the edge labels being the entries of the
!    sequence.  In general, there are many distinct de Bruijn sequences
!    for the same parameter M and N.  This program will only find one
!    of them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of letters in the alphabet.
!
!    Input, integer ( kind = 4 ) N, the number of letters in a codeword.
!
!    Output, integer ( kind = 4 ) STRING(M**N), a deBruijn string.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) inode(m**n)
  integer ( kind = 4 ) ivec(n-1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(m**n)
  integer ( kind = 4 ) jvec(n-1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) knode(m**n)
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) string(m**n)
  logical success
  integer ( kind = 4 ) trail(m**n)
!
!  Construct the adjacency information.
!
  nnode = m**(n-1)
  nedge = m**n

  iedge = 0

  do i = 1, nnode

    call index_unrank0 ( n - 1, m, i, ivec )

    do k = 1, m
      jvec(1:n-2) = ivec(2:n-1)
      jvec(n-1) = k
      call index_rank0 ( n - 1, m, jvec, j )
      iedge = iedge + 1
      inode(iedge) = i
      jnode(iedge) = j
      knode(iedge) = k
    end do

  end do
!
!  Determine a circuit.
!
  call digraph_arc_euler ( nnode, nedge, inode, jnode, success, trail )
!
!  The string is constructed from the labels of the edges in the circuit.
!
  string(1:nedge) = knode(trail(1:nedge))

  return
end
subroutine dec_add ( mantissa1, exponent1, mantissa2, exponent2, &
  dec_digit, mantissa, exponent )

!*****************************************************************************80
!
!! DEC_ADD adds two decimal quantities.
!
!  Discussion:
!
!    A decimal value is represented by MANTISSA * 10**EXPONENT.
!
!    The routine computes
!
!      MANTISSA * 10**EXPONENT =
!        MANTISSA1 * 10**EXPONENT1
!      + MANTISSA2 * 10**EXPONENT2
!
!    while trying to avoid integer overflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MANTISSA1, EXPONENT1, the first number to
!    be added.
!
!    Input, integer ( kind = 4 ) MANTISSA2, EXPONENT2, the second number to
!    be added.
!
!    Input, integer ( kind = 4 ) DEC_DIGIT, the number of decimal digits.
!
!    Output, integer ( kind = 4 ) MANTISSA, EXPONENT, the sum.
!
  implicit none

  integer ( kind = 4 ) dec_digit
  integer ( kind = 4 ) exponent
  integer ( kind = 4 ) exponent1
  integer ( kind = 4 ) exponent2
  integer ( kind = 4 ) mantissa
  integer ( kind = 4 ) mantissa1
  integer ( kind = 4 ) mantissa2
  integer ( kind = 4 ) mantissa3
  integer ( kind = 4 ) mantissa4

  if ( mantissa1 == 0 ) then
    mantissa = mantissa2
    exponent = exponent2
    return
  else if ( mantissa2 == 0 ) then
    mantissa = mantissa1
    exponent = exponent1
    return
  else if ( exponent1 == exponent2 ) then
    mantissa = mantissa1 + mantissa2
    exponent = exponent1
    call dec_round ( mantissa, exponent, dec_digit, mantissa, exponent )
    return
  end if
!
!  Line up the exponents.
!
  mantissa3 = mantissa1
  mantissa4 = mantissa2

  if ( exponent1 < exponent2 ) then
    mantissa4 = mantissa4 * 10**( exponent2 - exponent1 )
  else
    mantissa3 = mantissa3 * 10**( exponent1 - exponent2 )
  end if
!
!  Add the coefficients.
!
  mantissa = mantissa3 + mantissa4
  exponent = min ( exponent1, exponent2 )
!
!  Clean up the result.
!
  call dec_round ( mantissa, exponent, dec_digit, mantissa, exponent )

  return
end
subroutine dec_div ( mantissa1, exponent1, mantissa2, exponent2, &
  dec_digit, mantissa, exponent, ierror )

!*****************************************************************************80
!
!! DEC_DIV divides two decimal values.
!
!  Discussion:
!
!    A decimal value is represented by MANTISSA * 10**EXPONENT.
!
!    The routine computes
!
!      MANTISSA * 10**EXPONENT
!      = ( MANTISSA1 * 10**EXPONENT1 ) / ( MANTISSA2 * 10**EXPONENT2 )
!      = ( MANTISSA1 / MANTISSA2 ) * 10**( EXPONENT1 - EXPONENT2 )
!
!    while avoiding integer overflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MANTISSA1, EXPONENT1, the numerator.
!
!    Input, integer ( kind = 4 ) MANTISSA2, EXPONENT2, the denominator.
!
!    Input, integer ( kind = 4 ) DEC_DIGIT, the number of decimal digits.
!
!    Output, integer ( kind = 4 ) MANTISSA, EXPONENT, the result.
!
!    Output, integer ( kind = 4 ) IERROR.
!    0, no error occurred.
!    1, an error occurred.
!
  implicit none

  integer ( kind = 4 ) dec_digit
  real ( kind = 8 ) dval
  integer ( kind = 4 ) exponent
  integer ( kind = 4 ) exponent1
  integer ( kind = 4 ) exponent2
  integer ( kind = 4 ) exponent3
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) mantissa
  integer ( kind = 4 ) mantissa1
  integer ( kind = 4 ) mantissa2
  integer ( kind = 4 ) mantissa3
!
!  First special case, top fraction is 0.
!
  if ( mantissa1 == 0 ) then
    mantissa = 0
    exponent = 0
    return
  end if
!
!  First error, bottom of fraction is 0.
!
  if ( mantissa2 == 0 ) then
    ierror = 1
    mantissa = 0
    exponent = 0
    return
  end if
!
!  Second special case, result is 1.
!
  if ( mantissa1 == mantissa2 .and. exponent1 == exponent2 ) then
    mantissa = 1
    exponent = 0
    return
  end if
!
!  Third special case, result is power of 10.
!
  if ( mantissa1 == mantissa2 ) then
    mantissa = 1
    exponent = exponent1 - exponent2
    return
  end if
!
!  Fourth special case: MANTISSA1/MANTISSA2 is exact.
!
  if ( ( mantissa1 / mantissa2 ) * mantissa2 == mantissa1 ) then
    mantissa = mantissa1 / mantissa2
    exponent = exponent1 - exponent2
    return
  end if
!
!  General case.
!
  dval = real ( mantissa1, kind = 8 ) / real ( mantissa2, kind = 8 )

  call r8_to_dec ( dval, dec_digit, mantissa3, exponent3 )

  mantissa = mantissa3
  exponent = exponent3 + exponent1 - exponent2

  return
end
subroutine dec_mul ( mantissa1, exponent1, mantissa2, exponent2, &
  dec_digit, mantissa, exponent )

!*****************************************************************************80
!
!! DEC_MUL multiplies two decimals.
!
!  Discussion:
!
!    A decimal value is represented by MANTISSA * 10**EXPONENT.
!
!    The routine computes
!
!      MANTISSA * 10**EXPONENT
!      = ( MANTISSA1 * 10**EXPONENT1) * (MANTISSA2 * 10**EXPONENT2)
!      = ( MANTISSA1 * MANTISSA2 ) * 10**( EXPONENT1 + EXPONENT2 )
!
!    while avoiding integer overflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MANTISSA1, EXPONENT1, the first multiplier.
!
!    Input, integer ( kind = 4 ) MANTISSA2, EXPONENT2, the second multiplier.
!
!    Input, integer ( kind = 4 ) DEC_DIGIT, the number of decimal digits.
!
!    Output, integer ( kind = 4 ) MANTISSA, EXPONENT, the product.
!
  implicit none

  integer ( kind = 4 ) dec_digit
  real ( kind = 8 ) dval
  integer ( kind = 4 ) exponent
  integer ( kind = 4 ) exponent1
  integer ( kind = 4 ) exponent2
  integer ( kind = 4 ) exponent3
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) mantissa
  integer ( kind = 4 ) mantissa1
  integer ( kind = 4 ) mantissa2
  integer ( kind = 4 ) mantissa3
  real ( kind = 8 ) temp

  i_max = i4_huge ( )
!
!  The result is zero if either MANTISSA1 or MANTISSA2 is zero.
!
  if ( mantissa1 == 0 .or. mantissa2 == 0 ) then
    mantissa = 0
    exponent = 0
    return
  end if
!
!  The result is simple if either MANTISSA1 or MANTISSA2 is one.
!
  if ( abs ( mantissa1 ) == 1 .or. abs ( mantissa2 ) == 1 ) then
    mantissa = mantissa1 * mantissa2
    exponent = exponent1 + exponent2
    return
  end if

  temp = log ( real ( abs ( mantissa1 ), kind = 8 ) ) &
       + log ( real ( abs ( mantissa2 ), kind = 8 ) )

  if ( temp < log ( real ( i_max, kind = 8 ) ) ) then

    mantissa = mantissa1 * mantissa2
    exponent = exponent1 + exponent2

  else

    dval = real ( mantissa1, kind = 8 ) * real ( mantissa2, kind = 8 )

    call r8_to_dec ( dval, dec_digit, mantissa3, exponent3 )

    mantissa = mantissa3
    exponent = exponent3 + ( exponent1 + exponent2 )

  end if

  call dec_round ( mantissa, exponent, dec_digit, mantissa, exponent )

  return
end
subroutine dec_round ( mantissa1, exponent1, dec_digit, mantissa2, exponent2 )

!*****************************************************************************80
!
!! DEC_ROUND rounds a decimal fraction to a given number of digits.
!
!  Discussion:
!
!    A decimal value is represented by MANTISSA * 10**EXPONENT.
!
!    The routine takes an arbitrary decimal fraction makes sure that
!    MANTISSA has no more than DEC_DIGIT digits.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MANTISSA1, EXPONENT1, the coefficient and
!    exponent of a decimal fraction to be rounded.
!
!    Input, integer ( kind = 4 ) DEC_DIGIT, the number of decimal digits.
!
!    Output, integer ( kind = 4 ) MANTISSA2, EXPONENT2, the rounded coefficient
!    and  exponent of a decimal fraction.  MANTISSA2 has no more than
!    DEC_DIGIT decimal digits.
!
  implicit none

  integer ( kind = 4 ) dec_digit
  integer ( kind = 4 ) exponent1
  integer ( kind = 4 ) exponent2
  integer ( kind = 4 ) mantissa1
  integer ( kind = 4 ) mantissa2

  mantissa2 = mantissa1
  exponent2 = exponent1

  if ( mantissa2 == 0 ) then
    exponent2 = 0
    return
  end if

  do while ( 10**dec_digit <= abs ( mantissa2 ) )
    mantissa2 = nint ( real ( mantissa2, kind = 8 ) / 10.0D+00 )
    exponent2 = exponent2 + 1
  end do
!
!  Absorb trailing 0's into the exponent.
!
  do while ( ( mantissa2 / 10 ) * 10 == mantissa2 )
    mantissa2 = mantissa2 / 10
    exponent2 = exponent2 + 1
  end do

  return
end
subroutine dec_to_r8 ( mantissa, exponent, r )

!*****************************************************************************80
!
!! DEC_TO_R8 converts a decimal to an R8.
!
!  Discussion:
!
!    A decimal value is represented by MANTISSA * 10**EXPONENT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MANTISSA, EXPONENT, the coefficient and
!    exponent of the decimal value.
!
!    Output, real ( kind = 8 ) R, the equivalent real value.
!
  implicit none

  integer ( kind = 4 ) exponent
  integer ( kind = 4 ) mantissa
  real ( kind = 8 ) r

  r = mantissa * 10.0D+00**exponent

  return
end
subroutine dec_to_rat ( mantissa, exponent, rat_top, rat_bot )

!*****************************************************************************80
!
!! DEC_TO_RAT converts a decimal to a rational representation.
!
!  Discussion:
!
!    A decimal value is represented by MANTISSA * 10**EXPONENT.
!
!    A rational value is represented by RAT_TOP / RAT_BOT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MANTISSA, EXPONENT, the decimal number.
!
!    Output, integer ( kind = 4 ) RAT_TOP, RAT_BOT, the rational value.
!
  implicit none

  integer ( kind = 4 ) gcd
  integer ( kind = 4 ) exponent
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) mantissa
  integer ( kind = 4 ) rat_bot
  integer ( kind = 4 ) rat_top

  if ( exponent == 0 ) then
    rat_top = mantissa
    rat_bot = 1
  else if ( 0 < exponent ) then
    rat_top = mantissa * 10**exponent
    rat_bot = 1
  else
    rat_top = mantissa
    rat_bot = 10**( - exponent )
    gcd = i4_gcd ( rat_top, rat_bot )
    rat_top = rat_top / gcd
    rat_bot = rat_bot / gcd
  end if

  return
end
subroutine dec_to_s ( mantissa, exponent, s )

!*****************************************************************************80
!
!! DEC_TO_S returns a string representation of a decimal.
!
!  Discussion:
!
!    A decimal value is represented by MANTISSA * 10**EXPONENT.
!
!  Example:
!
!    MANTISSA EXPONENT   S
!    ----     ----       ------
!       0        0       0
!      21        3       21000
!      -3        0       -3
!     147       -2       14.7
!      16       -5       0.00016
!      34       30       Inf
!     123      -21       0.0000000000000000012
!      34      -30       0.0
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
!    Input, integer ( kind = 4 ) MANTISSA, EXPONENT, integers which represent
!    the decimal.
!
!    Output, character ( len = * ) S, the representation of the value.
!    The string is 'Inf' or '0.0' if the value was too large
!    or small to represent with a fixed point format.
!
  implicit none

  character ( len = 22 ) chrrep
  integer   ( kind = 4 ) exponent
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) iget1
  integer   ( kind = 4 ) iget2
  integer   ( kind = 4 ) iput1
  integer   ( kind = 4 ) iput2
  integer   ( kind = 4 ) mantissa
  integer   ( kind = 4 ) maxdigit
  integer   ( kind = 4 ) ndigit
  integer   ( kind = 4 ) nleft
  character ( len = * ) s

  s = ' '

  if ( mantissa == 0 ) then
    s = '0'
    return
  end if

  maxdigit = len ( s )
!
!  Store a representation of MANTISSA in CHRREP.
!
  write ( chrrep, '(i22)' ) mantissa
  call s_blank_delete ( chrrep )
  ndigit = len_trim ( chrrep )
!
!  Overflow if EXPONENT is positive, and MAXDIGIT < NDIGIT + EXPONENT.
!
  if ( 0 < exponent ) then
    if ( maxdigit < ndigit + exponent ) then
      s = 'Inf'
      return
    end if
  end if
!
!  Underflow if EXPONENT is negative, and MAXDIGIT < 3 + NDIGIT - EXPONENT.
!
  if ( exponent < 0 ) then
    if ( 0 < mantissa ) then
      if ( maxdigit < 3 - ndigit - exponent ) then
        s = '0.0'
        return
      end if
    else
      if ( maxdigit < 5 - ndigit - exponent ) then
        s = '0.0'
        return
      end if
    end if
  end if
!
!  If EXPONENT is nonnegative, insert trailing zeros.
!
  if ( 0 <= exponent ) then

    s(1:ndigit) = chrrep(1:ndigit)

    do i = ndigit + 1, ndigit + exponent
      s(i:i) = '0'
    end do

  else if ( exponent < 0 ) then

    iput2 = 0
    iget2 = 0
!
!  Sign.
!
    if ( mantissa < 0 ) then
      iput1 = 1
      iput2 = 1
      iget2 = 1
      s(iput1:iput2) = '-'
      ndigit = ndigit - 1
    end if
!
!  Digits of the integral part.
!
    if ( 0 < ndigit + exponent ) then
      iput1 = iput2 + 1
      iput2 = iput1 + ndigit + exponent -1
      iget1 = iget2 + 1
      iget2 = iget1 + ndigit + exponent - 1
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
    do i = 1, - exponent - ndigit
      iput1 = iput2 + 1
      iput2 = iput1
      s(iput1:iput2) = '0'
    end do

    nleft = min ( -exponent, ndigit )
    nleft = min ( nleft, maxdigit - iput2 )
    iput1 = iput2 + 1
    iput2 = iput1 + nleft - 1
    iget1 = iget2 + 1
    iget2 = iget1 + nleft - 1
    s(iput1:iput2) = chrrep(iget1:iget2)

  end if

  return
end
function dec_width ( mantissa, exponent )

!*****************************************************************************80
!
!! DEC_WIDTH returns the "width" of a decimal number.
!
!  Discussion:
!
!    A decimal value is represented by MANTISSA * 10**EXPONENT.
!
!    The "width" of a decimal number is the number of characters
!    required to print it.
!
!  Example:
!
!    Mantissa  Exponent Width  Representation:
!
!         523      -1       4           5.23
!         134       2       5       13400
!           0      10       1           0
!      123456     -10      12           0.0000123456
!      123456      -3       7         123.456
!      123456       0       6      123456
!      123456       3       9   123456000
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MANTISSA, EXPONENT, the decimal number.
!
!    Output, integer ( kind = 4 ) DEC_WIDTH, the "width" of the decimal number.
!
  implicit none

  integer ( kind = 4 ) dec_width
  integer ( kind = 4 ) exponent
  integer ( kind = 4 ) mantissa
  integer ( kind = 4 ) mantissa_abs
  integer ( kind = 4 ) ten_pow
  integer ( kind = 4 ) value

  value = 1
  ten_pow = 10

  if ( mantissa == 0 ) then
    dec_width = value
    return
  end if

  mantissa_abs = abs ( mantissa )

  do while ( ten_pow <= mantissa_abs )
    value = value + 1
    ten_pow = ten_pow * 10
  end do

  if ( 0 < exponent ) then
    value = value + exponent
  else if ( exponent < 0 ) then
    value = max ( value, 1-exponent )
!
!  An internal decimal point adds one position.
!
    if ( 0 < value ) then
      value = value + 1
!
!  A leading "0." adds two positions.
!
    else
      value = 2 - value
    end if
  end if

  if ( mantissa < 0 ) then
    value = value + 1
  end if

  dec_width = value

  return
end
subroutine decmat_det ( n, atop, abot, dec_digit, dtop, dbot, ierror )

!*****************************************************************************80
!
!! DECMAT_DET finds the determinant of an N by N matrix of decimal entries.
!
!  Discussion:
!
!    The brute force method is used.  The routine should only be used for
!    small matrices, since this calculation requires the summation of N!
!    products of N numbers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of A.
!
!    Input, integer ( kind = 4 ) ATOP(N,N), ABOT(N,N), the decimal
!    representation of the matrix.
!
!    Output, integer ( kind = 4 ) DTOP, DBOT, the decimal determinant of
!    the matrix.
!
!    Output, integer ( kind = 4 ) IERROR.
!    0, no error occurred.
!    1, an error occurred (probably overflow).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) dec_digit
  logical even
  integer ( kind = 4 ) i
  integer ( kind = 4 ) abot(n,n)
  integer ( kind = 4 ) atop(n,n)
  integer ( kind = 4 ) iarray(n)
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ibot1
  integer ( kind = 4 ) ibot2
  integer ( kind = 4 ) dbot
  integer ( kind = 4 ) dtop
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) itop1
  integer ( kind = 4 ) itop2
  logical more

  ierror = 0
  more = .false.
  dtop = 0
  dbot = 1
!
!  Compute the next permutation.
!
  do

    call perm_next ( n, iarray, more, even )
!
!  The sign of this term depends on the sign of the permutation.
!
    if ( even ) then
      itop = 1
    else
      itop = -1
    end if
!
!  Choose one item from each row, as specified by the permutation,
!  and multiply them
!
    ibot = 0

    do i = 1, n

      itop1 = itop
      ibot1 = ibot
      itop2 = atop(i,iarray(i))
      ibot2 = abot(i,iarray(i))

      call dec_mul ( itop1, ibot1, itop2, ibot2, dec_digit, itop, ibot )

    end do
!
!  Add this term to the total.
!
    itop1 = itop
    ibot1 = ibot

    call dec_add ( itop1, ibot1, dtop, dbot, dec_digit, itop, ibot )

    dtop = itop
    dbot = ibot

    if ( .not. more ) then
      exit
    end if

  end do

  return
end
subroutine decmat_print ( m, n, a, b, title )

!*****************************************************************************80
!
!! DECMAT_PRINT prints out decimal vectors and matrices.
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the matrix.
!
!    Input, integer ( kind = 4 ) A(M,N), B(M,N), the decimal matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) a(m,n)
  integer   ( kind = 4 ) b(m,n)
  character ( len = 22 ) chrtmp
  character ( len = 10 ) chrtmp2
  character ( len = 10 ) chrtmp3
  character ( len = 40 ) format2
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) imax
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) jmax
  integer   ( kind = 4 ) jmin
  integer   ( kind = 4 ) khi
  integer   ( kind = 4 ) klo
  integer   ( kind = 4 ) kmax
  integer   ( kind = 4 ) lenc
  integer   ( kind = 4 ), parameter :: ncolum = 80
  integer   ( kind = 4 ) npline
  character ( len = 100 ) output
  character ( len = * ) title
!
!  Figure out how wide we must make each column.
!
  imax = 0
  jmax = 0

  do i = 1, m
    do j = 1, n

      call dec_to_s ( a(i,j), b(i,j), chrtmp )
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

  do jmin = 1, n, npline

    jmax = min ( jmin + npline - 1, n )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(a)' ) ' '

    if ( 1 < jmin .or. jmax < n ) then
      write ( output, * ) 'Columns ', jmin, ' to ', jmax
      call s_blanks_delete ( output )
      write ( *, '(a)' ) trim ( output )
      write ( *, '(a)' ) ' '
    end if

    do i = 1, m

      output = ' '

      do j = jmin, jmax
        klo = 4 + ( j - jmin ) * kmax + 1
        khi = 4 + ( j - jmin ) * kmax + kmax
        call dec_to_s ( a(i,j), b(i,j), chrtmp )
        output(klo:khi) = adjustr ( chrtmp(1:kmax) )
      end do

      write ( *, '(a)' ) trim ( output )

    end do

  end do

  return
end
subroutine derange_back_candidate ( n, maxstack, a, k, nstack, stack, ncan )

!*****************************************************************************80
!
!! DERANGE_BACK_CANDIDATE finds values for the K-th entry of a derangement.
!
!  Discussion:
!
!    A derangement of N objects is a permutation which leaves no object
!    unchanged.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the derangement.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the maximum stack length.
!
!    Input, integer ( kind = 4 ) A(N).  The first K-1 entries of A
!    record the currently set values of the derangement.
!
!    Input, integer ( kind = 4 ) K, the entry of the derangement for which
!    candidates are to be found.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the length of the stack.
!
!    Input/output, integer ( kind = 4 ) STACK(MAXSTACK).  On output, we have
!    added the candidates for entry K to the end of the stack.
!
!    Input/output, integer ( kind = 4 ) NCAN(N), the number of candidates
!    for each level.
!
  implicit none

  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) ican
  integer ( kind = 4 ) ifree(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ncan(n)
  integer ( kind = 4 ) nfree
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) stack(maxstack)
!
!  Consider all the integers from 1 through N that have not been used yet.
!
  nfree = n - k + 1

  call perm_free ( k - 1, a, nfree, ifree )
!
!  Everything but K is a legitimate candidate for the K-th entry.
!
  ncan(k) = 0

  do ican = 1, nfree

    if ( ifree(ican) /= k ) then

      ncan(k) = ncan(k) + 1
      nstack = nstack + 1

      if ( maxstack < nstack ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DERANGE_BACK_CANDIDATE - Fatal error!'
        write ( *, '(a,i8)' ) '  Exceeding stacksize limit of ', maxstack
        stop
      end if

      stack(nstack) = ifree(ican)

    end if

  end do

  return
end
subroutine derange_back_next ( n, a, more )

!*****************************************************************************80
!
!! DERANGE_BACK_NEXT returns the next derangement of N items.
!
!  Discussion:
!
!    A derangement of N objects is a permutation which leaves no object
!    unchanged.
!
!    This routine uses backtracking.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be deranged.
!    N should be at least 2.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On the first call, the input value of A is not important.
!    On return with MORE = TRUE, A contains the next derangement.
!    On subsequent input, A should not be changed.
!
!    Input/output, logical MORE.
!    On first call, set MORE to FALSE, and do not alter it after.
!    On return, MORE is TRUE if another derangement is being treturned in A,
!    and FALSE if no more can be found.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ), save :: indx = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: maxstack = 0
  logical more
  integer ( kind = 4 ), save, allocatable, dimension (:) :: ncan
  integer ( kind = 4 ), save :: nstack = 0
  integer ( kind = 4 ), save, allocatable, dimension (:) :: stack

  if ( .not. more ) then

    if ( n < 2 ) then
      more = .false.
      return
    end if

    indx = 0
    k = 0
    maxstack = ( n * ( n + 1 ) ) / 2
    nstack = 0

    if ( allocated ( stack ) ) then
      deallocate ( stack )
    end if

    allocate ( stack(1:maxstack) )
    stack(1:maxstack) = 0

    if ( allocated ( ncan ) ) then
      deallocate ( ncan )
    end if

    allocate ( ncan(1:n) )
    ncan(1:n) = 0

    more = .true.

  end if

  do

    call i4vec_backtrack ( n, maxstack, stack, a, indx, k, nstack, ncan )

    if ( indx == 1 ) then

      exit

    else if ( indx == 2 ) then

      call derange_back_candidate ( n, maxstack, a, k, nstack, stack, ncan )

    else

      more = .false.
      exit

    end if

  end do

  return
end
subroutine derange_check ( n, a, deranged )

!*****************************************************************************80
!
!! DERANGE_CHECK determines whether a permutation is a derangement.
!
!  Discussion:
!
!    A derangement of N objects is a permutation which leaves no object
!    unchanged.
!
!    A derangement of the integers 1 through N is a permutation of the
!    integers such that the first value is not 1, the second is not 2,
!    and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects permuted.
!
!    Input, integer ( kind = 4 ) A(N), a permutation of the integers 1
!    through N.
!
!    Output, logical DERANGED, is TRUE if A is a derangement, and
!    FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  logical deranged
  integer ( kind = 4 ) i

  do i = 1, n
    if ( a(i) == i ) then
      deranged = .false.
      return
    end if
  end do

  deranged = .true.

  return
end
function derange_enum ( n )

!*****************************************************************************80
!
!! DERANGE_ENUM returns the number of derangements of N objects.
!
!  Discussion:
!
!    A derangement of N objects is a permutation which leaves no object
!    unchanged.
!
!    A derangement of N objects is a permutation with no fixed
!    points.  If we symbolize the permutation operation by "P",
!    then for a derangment, P(I) is never equal to I.
!
!    The number of derangements of N objects is sometimes called
!    the subfactorial function, or the derangement number D(N).
!
!    D(N) is the number of ways of placing N non-attacking rooks on
!    an N by N chessboard with one diagonal deleted.
!
!    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
!
!    The number of permutations with exactly K items in the right
!    place is COMB(N,K) * D(N-K).
!
!  Recursion:
!
!      D(0) = 1
!      D(1) = 0
!      D(2) = 1
!      D(N) = (N-1) * ( D(N-1) + D(N-2) )
!
!    or
!
!      D(0) = 1
!      D(1) = 0
!      D(N) = N * D(N-1) + (-1)**N
!
!  Formula:
!
!    D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
!
!    Based on the inclusion/exclusion law.
!
!    D(N) = nint ( N! / E )
!
!  First values:
!
!     N         D(N)
!     0           1
!     1           0
!     2           1
!     3           2
!     4           9
!     5          44
!     6         265
!     7        1854
!     8       14833
!     9      133496
!    10     1334961
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!
!    Output, integer ( kind = 4 ) DERANGE_ENUM, the number of derangements
!    of N objects.
!
  implicit none

  integer ( kind = 4 ) derange_enum
  integer ( kind = 4 ) dn
  integer ( kind = 4 ) dnm1
  integer ( kind = 4 ) dnm2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  if ( n < 0 ) then

    dn = 0

  else if ( n == 0 ) then

    dn = 1

  else if ( n == 1 ) then

    dn = 0

  else if ( n == 2 ) then

    dn = 1

  else

    dnm1 = 0
    dn = 1

    do i = 3, n
      dnm2 = dnm1
      dnm1 = dn
      dn = ( i - 1 ) * ( dnm1 + dnm2 )
    end do

  end if

  derange_enum = dn

  return
end
subroutine derange_enum2 ( n, d )

!*****************************************************************************80
!
!! DERANGE_ENUM2 returns the number of derangements of 0 through N objects.
!
!  Discussion:
!
!    A derangement of N objects is a permutation which leaves no object
!    unchanged.
!
!    A derangement of N objects is a permutation with no fixed
!    points.  If we symbolize the permutation operation by "P",
!    then for a derangment, P(I) is never equal to I.
!
!    The number of derangements of N objects is sometimes called
!    the subfactorial function, or the derangement number D(N).
!
!    D(N) is the number of ways of placing N non-attacking rooks on
!    an N by N chessboard with one diagonal deleted.
!
!    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
!
!    The number of permutations with exactly K items in the right
!    place is COMB(N,K) * D(N-K).
!
!  Recursion:
!
!      D(0) = 1
!      D(1) = 0
!      D(2) = 1
!      D(N) = (N-1) * ( D(N-1) + D(N-2) )
!
!    or
!
!      D(0) = 1
!      D(1) = 0
!      D(N) = N * D(N-1) + (-1)**N
!
!  Formula:
!
!    D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
!
!    Based on the inclusion/exclusion law.
!
!    D(N) = nint ( N! / E )
!
!  First values:
!
!     N         D(N)
!     0           1
!     1           0
!     2           1
!     3           2
!     4           9
!     5          44
!     6         265
!     7        1854
!     8       14833
!     9      133496
!    10     1334961
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the maximum number of objects to be
!    permuted.
!
!    Output, integer ( kind = 4 ) D(0:N); D(I) is the number of derangements of
!    I objects.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) d(0:n)
  integer ( kind = 4 ) i

  d(0) = 1
  d(1) = 0

  do i = 2, n
    d(i) = ( i - 1 ) * ( d(i-1) + d(i-2) )
  end do

  return
end
function derange_enum3 ( n )

!*****************************************************************************80
!
!! DERANGE_ENUM3 returns the number of derangements of 0 through N objects.
!
!  Discussion:
!
!    A derangement of N objects is a permutation which leaves no object
!    unchanged.
!
!    A derangement of N objects is a permutation with no fixed
!    points.  If we symbolize the permutation operation by "P",
!    then for a derangment, P(I) is never equal to I.
!
!    The number of derangements of N objects is sometimes called
!    the subfactorial function, or the derangement number D(N).
!
!    D(N) is the number of ways of placing N non-attacking rooks on
!    an N by N chessboard with one diagonal deleted.
!
!    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
!
!    The number of permutations with exactly K items in the right
!    place is COMB(N,K) * D(N-K).
!
!  Recursion:
!
!      D(0) = 1
!      D(1) = 0
!      D(2) = 1
!      D(N) = (N-1) * ( D(N-1) + D(N-2) )
!
!    or
!
!      D(0) = 1
!      D(1) = 0
!      D(N) = N * D(N-1) + (-1)**N
!
!  Formula:
!
!    D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
!
!    Based on the inclusion/exclusion law.
!
!    D(N) = nint ( N! / E )
!
!  First values:
!
!     N         D(N)
!     0           1
!     1           0
!     2           1
!     3           2
!     4           9
!     5          44
!     6         265
!     7        1854
!     8       14833
!     9      133496
!    10     1334961
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the maximum number of objects to be
!    permuted.
!
!    Output, integer ( kind = 4 ) DERANGE_ENUM3, the number of derangements
!    of N objects.
!
  implicit none

  integer ( kind = 4 ) derange_enum3
  real ( kind = 8 ), parameter :: e = 2.718281828459045D+00
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_factorial

  if ( n < 0 ) then
    derange_enum3 = -1
  else if ( n == 0 ) then
    derange_enum3 = 1
  else if ( n == 1 ) then
    derange_enum3 = 0
  else
    derange_enum3 = nint ( r8_factorial ( n ) / e )
  end if

  return
end
subroutine derange_weed_next ( n, a, more )

!*****************************************************************************80
!
!! DERANGE_WEED_NEXT computes all derangements of N objects, one at a time.
!
!  Discussion:
!
!    A derangement of N objects is a permutation which leaves no object
!    unchanged.
!
!    A derangement of N objects is a permutation with no fixed
!    points.  If we symbolize the permutation operation by "P",
!    then for a derangment, P(I) is never equal to I.
!
!    The number of derangements of N objects is sometimes called
!    the subfactorial function, or the derangement number D(N).
!
!    This routine simply generates all permutations, one at a time,
!    and weeds out those that are not derangements.
!
!  Example:
!
!    Here are the derangements when N = 4:
!
!    2143  3142  4123
!    2341  3412  4312
!    2413  3421  4321
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
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On first call, the input contents of A are unimportant.  But
!    on the second and later calls, the input value of A should be
!    the output value returned on the previous call.
!    On output, A contains the next derangement.
!
!    Input/output, logical MORE.
!    Set MORE = FALSE before the first call.
!    MORE will be reset to TRUE and a derangement will be returned.
!    Each new call produces a new derangement until MORE is returned FALSE.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  logical deranged
  integer ( kind = 4 ) derange_enum
  integer ( kind = 4 ), save :: maxder = 0
  logical more
  integer ( kind = 4 ), save :: numder = 0
!
!  Initialization on call with MORE = FALSE.
!
  if ( .not. more ) then

    maxder = derange_enum ( n )
    numder = 0

  end if
!
!  Watch out for cases where there are no derangements.
!
  if ( maxder == 0 ) then
    more = .false.
    return
  end if
!
!  Get the next permutation.
!
  do

    call perm_lex_next ( n, a, more )
!
!  See if it is a derangment.
!
    call derange_check ( n, a, deranged )

    if ( deranged ) then
      exit
    end if

  end do

  numder = numder + 1

  if ( maxder <= numder ) then
    more = .false.
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
subroutine digraph_arc_euler ( nnode, nedge, inode, jnode, success, trail )

!*****************************************************************************80
!
!! DIGRAPH_ARC_EULER returns an Euler circuit in a digraph.
!
!  Discussion:
!
!    An Euler circuit of a digraph is a path which starts and ends at
!    the same node and uses each directed edge exactly once.  A digraph is
!    eulerian if it has an Euler circuit.  The problem is to decide whether
!    a given digraph is eulerian and to find an Euler circuit if the
!    answer is affirmative.
!
!
!    A digraph has an Euler circuit if and only if the number of incoming
!    edges is equal to the number of outgoing edges at each node.
!
!    This characterization gives a straightforward procedure to decide whether
!    a digraph is eulerian.  Furthermore, an Euler circuit in an eulerian
!    digraph G of NEDGE edges can be determined by the following method:
!
!      STEP 1: Choose any node U as the starting node, and traverse any edge
!      ( U, V ) incident to node U, and than traverse any unused edge incident
!      to node U.  Repeat this process of traversing unused edges until the
!      starting node U is reached.  Let P be the resulting walk consisting of
!      all used edges.  If all edges of G are in P, than stop.
!
!      STEP 2: Choose any unused edge ( X,  Y) in G such that X is
!      in P and Y is not in P.  Use node X as the starting node and
!      find another walk Q using all unused edges as in step 1.
!
!      STEP 3: Walk P and walk Q share a common node X, they can be merged
!      to form a walk R by starting at any node S of P and to traverse P
!      until node X is reached; than, detour and traverse all edges of Q
!      until node X is reached and continue to traverse the edges of P until
!      the starting node S is reached.  Set P = R.
!
!      STEP 4: Repeat steps 2 and 3 until all edges are used.
!
!    The running time of the algorithm is O ( NEDGE ).
!
!    The digraph is assumed to be connected.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2000
!
!  Author:
!
!    Original FORTRAN77 version by Hang Tong Lau.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); the I-th edge
!    starts at node INODE(I) and ends at node JNODE(I).
!
!    Output, logical SUCCESS, is TRUE if an Euler circuit was found,
!    and FALSE otherwise.
!
!    Output, integer ( kind = 4 ) TRAIL(NEDGE).  TRAIL(I) is the edge number
!    of the I-th edge in the Euler circuit.
!
  implicit none

  integer ( kind = 4 ) nedge

  logical candid(nedge)
  integer ( kind = 4 ) endnod(nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) istak
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) len
  integer ( kind = 4 ) lensol
  integer ( kind = 4 ) lenstk
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) stack(2*nedge)
  logical success
  integer ( kind = 4 ) trail(nedge)
!
!  Check if the digraph is eulerian.
!
  trail(1:nedge) = 0
  endnod(1:nedge) = 0

  do i = 1, nedge
    j = inode(i)
    trail(j) = trail(j) + 1
    j = jnode(i)
    endnod(j) = endnod(j) + 1
  end do

  do i = 1, nnode
    if ( trail(i) /= endnod(i) ) then
      success = .false.
      return
    end if
  end do
!
!  The digraph is eulerian; find an Euler circuit.
!
  success = .true.
  lensol = 1
  lenstk = 0
!
!  Find the next edge.
!
  do

    if ( lensol == 1 ) then

      endnod(1) = inode(1)
      stack(1) = 1
      stack(2) = 1
      lenstk = 2

    else

      l = lensol - 1

      if ( lensol /= 2 ) then
        endnod(l) = inode(trail(l)) + jnode(trail(l)) - endnod(l-1)
      end if

      k = endnod(l)

      do i = 1, nedge
        candid(i) = ( k == jnode(i) )
      end do

      do i = 1, l
        candid(trail(i)) = .false.
      end do

      len = lenstk

      do i = 1, nedge

        if ( candid(i) ) then
          len = len + 1
          stack(len) = i
        end if

      end do

      stack(len+1) = len - lenstk
      lenstk = len + 1

    end if

    do

      istak = stack(lenstk)
      lenstk = lenstk - 1

      if ( istak /= 0 ) then
        exit
      end if

      lensol = lensol - 1

      if ( lensol == 0 ) then
        call i4vec_reverse ( nedge, trail )
        return
      end if

    end do

    trail(lensol) = stack(lenstk)
    stack(lenstk) = istak - 1

    if ( lensol == nedge ) then
      exit
    end if

    lensol = lensol + 1

  end do

  call i4vec_reverse ( nedge, trail )

  return
end
subroutine digraph_arc_print ( nedge, inode, jnode, title )

!*****************************************************************************80
!
!! DIGRAPH_ARC_PRINT prints out a digraph from an edge list.
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
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the beginning
!    and end nodes of the edges.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) jnode(nedge)
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nedge
    write ( *, '(i8,4x,2i8)' ) i, inode(i), jnode(i)
  end do

  return
end
subroutine diophantine ( a, b, c, ierror, x, y )

!*****************************************************************************80
!
!! DIOPHANTINE solves a Diophantine equation A * X + B * Y = C.
!
!  Discussion:
!
!    Given integers A, B and C, produce X and Y so that
!
!      A * X + B * Y = C.
!
!    In general, the equation is solvable if and only if the
!    greatest common divisor of A and B also divides C.
!
!    A solution (X,Y) of the Diophantine equation also gives the solution
!    X to the congruence equation:
!
!      A * X = C mod ( B ).
!
!    Generally, if there is one nontrivial solution, there are an infinite
!    number of solutions to a Diophantine problem.
!    If (X0,Y0) is a solution, then so is ( X0+T*B/D, Y0-T*A/D ) where
!    T is any integer, and D is the greatest common divisor of A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein, editor,
!    CRC Concise Encylopedia of Mathematics,
!    CRC Press, 1998, page 446.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, C, the coefficients of the Diophantine
!    equation.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error, X and Y were computed.
!    1, A = B = 0, C is nonzero.
!    2, A = 0, B and C nonzero, but C is not a multiple of B.
!    3, A nonzero, B zero, C nonzero, but C is not a multiple of A.
!    4, A, B, C nonzero, but GCD of A and B does not divide C.
!    5, the algorithm ran out of internal space.
!
!    Output, integer ( kind = 4 ) X, Y, the solution of the Diophantine
!    equation.  Note that the algorithm will attempt to return a solution with
!    smallest Euclidean norm.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 100

  integer ( kind = 4 ) a
  integer ( kind = 4 ) a_copy
  integer ( kind = 4 ) a_mag
  integer ( kind = 4 ) a_sign
  integer ( kind = 4 ) b
  integer ( kind = 4 ) b_copy
  integer ( kind = 4 ) b_mag
  integer ( kind = 4 ) b_sign
  integer ( kind = 4 ) c
  integer ( kind = 4 ) c_copy
  integer ( kind = 4 ) g
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) q(nmax)
  logical swap
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
!
!  Defaults for output parameters.
!
  ierror = 0
  x = 0
  y = 0
!
!  Special cases.
!
  if ( a == 0 .and. b == 0 .and. c == 0 ) then
    x = 0
    y = 0
    return
  else if ( a == 0 .and. b == 0 .and. c /= 0 ) then
    ierror = 1
    x = 0
    y = 0
    return
  else if ( a == 0 .and. b /= 0 .and. c == 0 ) then
    x = 0
    y = 0
    return
  else if ( a == 0 .and. b /= 0 .and. c /= 0 ) then
    x = 0
    y = c / b
    if ( mod ( c, b ) /= 0 ) then
      ierror = 2
    end if
    return
  else if ( a /= 0 .and. b == 0 .and. c == 0 ) then
    x = 0
    y = 0
    return
  else if ( a /= 0 .and. b == 0 .and. c /= 0 ) then
    x = c / a
    y = 0
    if ( mod ( c, a ) /= 0 ) then
      ierror = 3
    end if
    return
  else if ( a /= 0 .and. b /= 0 .and. c == 0 ) then
    g = i4_gcd ( a, b )
    x = b / g
    y = -a / g
    return
  end if
!
!  Now handle the "general" case: A, B and C are nonzero.
!
!  Step 1: Compute the GCD of A and B, which must also divide C.
!
  g = i4_gcd ( a, b )

  if ( mod ( c, g ) /= 0 ) then
    ierror = 4
    return
  end if

  a_copy = a / g
  b_copy = b / g
  c_copy = c / g
!
!  Step 2: Split A and B into sign and magnitude.
!
  a_mag = abs ( a_copy )
  a_sign = sign ( 1, a_copy )
  b_mag = abs ( b_copy )
  b_sign = sign ( 1, b_copy )
!
!  Another special case, A_MAG = 1 or B_MAG = 1.
!
  if ( a_mag == 1 ) then
    x = a_sign * c_copy
    y = 0
    return
  else if ( b_mag == 1 ) then
    x = 0
    y = b_sign * c_copy
    return
  end if
!
!  Step 3: Produce the Euclidean remainder sequence.
!
  if ( b_mag <= a_mag ) then

    swap = .false.
    q(1) = a_mag
    q(2) = b_mag

  else

    swap = .true.
    q(1) = b_mag
    q(2) = a_mag

  end if

  n = 3

  do

    q(n) = mod ( q(n-2), q(n-1) )

    if ( q(n) == 1 ) then
      exit
    end if

    n = n + 1

    if ( nmax < n ) then
      ierror = 5
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIOPHANTINE - Fatal error!'
      write ( *, '(a)' ) '  Exceeded number of iterations.'
      stop
    end if

  end do
!
!  Step 4: Now go backwards to solve X * A_MAG + Y * B_MAG = 1.
!
  y = 0
  do k = n, 2, -1
    x = y
    y = ( 1 - x * q(k-1) ) / q(k)
  end do
!
!  Step 5: Undo the swapping.
!
  if ( swap ) then
    call i4_swap ( x, y )
  end if
!
!  Step 6: Now apply signs to X and Y so that X * A + Y * B = 1.
!
  x = x * a_sign
  y = y * b_sign
!
!  Step 7: Multiply by C, so that X * A + Y * B = C.
!
  x = x * c_copy
  y = y * c_copy
!
!  Step 8: Given a solution (X,Y), try to find the solution of
!  minimal magnitude.
!
  call diophantine_solution_minimize ( a_copy, b_copy, x, y )

  return
end
subroutine diophantine_solution_minimize ( a, b, x, y )

!*****************************************************************************80
!
!! DIOPHANTINE_SOLUTION_MINIMIZE: minimal solution of a Diophantine equation.
!
!  Discussion:
!
!    Given a solution (X,Y) of a Diophantine equation:
!
!      A * X + B * Y = C.
!
!    then there are an infinite family of solutions of the form
!
!      ( X(i), Y(i) ) = ( X + i * B, Y - i * A )
!
!    An integral solution of minimal Euclidean norm can be found by
!    tentatively moving along the vectors (B,-A) and (-B,A) one step
!    at a time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein, editor,
!    CRC Concise Encylopedia of Mathematics,
!    CRC Press, 1998, page 446.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the coefficients of the Diophantine
!    equation.  A and B are assumed to be relatively prime.
!
!    Input/output, integer ( kind = 4 ) X, Y, on input, a solution of the
!    Diophantine equation.  On output, a solution of minimal Euclidean norm.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) norm
  real ( kind = 8 ) norm_new
  real ( kind = 8 ) t
  integer ( kind = 4 ) x
  integer ( kind = 4 ) xnew
  integer ( kind = 4 ) y
  integer ( kind = 4 ) ynew
!
!  Compute the minimum for T real, and then look nearby.
!
  t = ( - real ( b, kind = 8 ) * real ( x, kind = 8 ) &
    + real ( a, kind = 8 ) * real ( y, kind = 8 ) ) &
    / ( real ( a, kind = 8 )**2 + real ( b, kind = 8 )**2 )

  x = x + nint ( t ) * b
  y = y - nint ( t ) * a
!
!  Now look nearby.
!
  norm = ( real ( x, kind = 8 ) )**2 + ( real ( y, kind = 8 ) )**2

  do

    xnew = x + b
    ynew = y - a

    norm_new = ( real ( xnew, kind = 8 ) )**2 + ( real ( ynew, kind = 8 ) )**2

    if ( norm <= norm_new ) then
      exit
    end if

    x = xnew
    y = ynew
    norm = norm_new

  end do

  do

    xnew = x - b
    ynew = y + a

    norm_new = ( real ( xnew, kind = 8 ) )**2 + ( real ( ynew, kind = 8 ) )**2

    if ( norm <= norm_new ) then
      exit
    end if

    x = xnew
    y = ynew
    norm = norm_new

  end do

  return
end
subroutine dvec_add ( n, dvec1, dvec2, dvec3 )

!*****************************************************************************80
!
!! DVEC_ADD adds two (signed) DVEC's.
!
!  Discussion:
!
!    A DVEC is an integer vector of decimal digits, intended to
!    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
!    is the coefficient of 10**(N-2), and DVEC(N) contains sign
!    information.  It is 0 if the number is positive, and 9 if
!    the number is negative.
!
!  Example:
!
!    N = 4
!
!      DVEC1     +   DVEC2     =   DVEC3
!
!    ( 0 0 1 7 ) + ( 0 1 0 4 ) = ( 0 0 1 2 1 )
!
!          17    +       104   =         121
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) DVEC1(N), DVEC2(N), the vectors to be added.
!
!    Output, integer ( kind = 4 ) DVEC3(N), the sum of the two input vectors.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 10
  integer ( kind = 4 ) dvec1(n)
  integer ( kind = 4 ) dvec2(n)
  integer ( kind = 4 ) dvec3(n)
  integer ( kind = 4 ) i
  logical overflow

  overflow = .false.

  dvec3(1:n) = dvec1(1:n) + dvec2(1:n)

  do i = 1, n

    do while ( base <= dvec3(i) )

      dvec3(i) = dvec3(i) - base

      if ( i < n ) then
        dvec3(i+1) = dvec3(i+1) + 1
      else
        overflow = .true.
      end if

    end do

  end do

  return
end
subroutine dvec_complementx ( n, dvec1, dvec2 )

!*****************************************************************************80
!
!! DVEC_COMPLEMENTX computes the ten's complement of a DVEC.
!
!  Discussion:
!
!    A DVEC is an integer vector of decimal digits, intended to
!    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
!    is the coefficient of 10**(N-2), and DVEC(N) contains sign
!    information.  It is 0 if the number is positive, and 9 if
!    the number is negative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) DVEC1(N), the vector to be complemented.
!
!    Output, integer ( kind = 4 ) DVEC2(N), the complemented vector.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 10
  integer ( kind = 4 ) dvec1(n)
  integer ( kind = 4 ) dvec2(n)
  integer ( kind = 4 ) dvec3(n)
  integer ( kind = 4 ) dvec4(n)

  dvec3(1:n) = ( base - 1 ) - dvec1(1:n)

  dvec4(1) = 1
  dvec4(2:n) = 0

  call dvec_add ( n, dvec3, dvec4, dvec2 )

  return
end
subroutine dvec_mul ( n, dvec1, dvec2, dvec3 )

!*****************************************************************************80
!
!! DVEC_MUL computes the product of two DVEC's.
!
!  Discussion:
!
!    A DVEC is an integer vector of decimal digits, intended to
!    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
!    is the coefficient of 10**(N-2), and DVEC(N) contains sign
!    information.  It is 0 if the number is positive, and 9 if
!    the number is negative.
!
!    Since the user may want to make calls like
!
!      call dvec_mul ( n, dvec1, dvec1, dvec3 )
!    or even
!      call dvec_mul ( n, dvec1, dvec1, dvec1 )
!
!    we need to copy the arguments, work on them, and then copy out the result.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) DVEC1(N), DVEC2(N), the vectors to be
!    multiplied.
!
!    Output, integer ( kind = 4 ) DVEC3(N), the product of the two input
!    vectors.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 10
  integer ( kind = 4 ) carry
  integer ( kind = 4 ) dvec1(n)
  integer ( kind = 4 ) dvec2(n)
  integer ( kind = 4 ) dvec3(n)
  integer ( kind = 4 ) dveca(n)
  integer ( kind = 4 ) dvecb(n)
  integer ( kind = 4 ) dvecc(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) product_sign
!
!  Copy the input.
!
  dveca(1:n) = dvec1(1:n)
  dvecb(1:n) = dvec2(1:n)
!
!  Record the sign of the product.
!  Make the factors positive.
!
  product_sign = 1

  if ( dveca(n) /= 0 ) then
    product_sign = - product_sign
    call dvec_complementx ( n, dveca, dveca )
  end if

  if ( dvecb(n) /= 0 ) then
    product_sign = - product_sign
    call dvec_complementx ( n, dvecb, dvecb )
  end if

  dvecc(1:n) = 0
!
!  Multiply.
!
  do i = 1, n-1
    dvecc(i:n-1) = dvecc(i:n-1) + dveca(i) * dvecb(1:n-i)
  end do
!
!  Take care of carries.
!  Unlike the DVEC_ADD routine, we do NOT allow carries into the
!  N-th position.
!
  do i = 1, n-1

    carry = dvecc(i) / base
    dvecc(i) = dvecc(i) - carry * base

    if ( i < n - 1 ) then
      dvecc(i+1) = dvecc(i+1) + carry
    end if

  end do
!
!  Take care of the sign of the product.
!
  if ( product_sign < 0 ) then
    call dvec_complementx ( n, dvecc, dvecc )
  end if
!
!  Copy the output.
!
  dvec3(1:n) = dvecc(1:n)

  return
end
subroutine dvec_print ( n, dvec, title )

!*****************************************************************************80
!
!! DVEC_PRINT prints a DVEC, with an optional title.
!
!  Discussion:
!
!    A DVEC is an integer vector of decimal digits, intended to
!    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
!    is the coefficient of 10**(N-2), and DVEC(N) contains sign
!    information.  It is 0 if the number is positive, and 9 if
!    the number is negative.
!
!    The vector is printed "backwards", that is, the first entry
!    printed is DVEC(N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) DVEC(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) dvec(n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do ihi = n, 1, -80
    ilo = max ( ihi - 79, 1 )
    write ( *, '(2x,80i1)' ) dvec(ihi:ilo:-1)
  end do

  return
end
subroutine dvec_sub ( n, dvec1, dvec2, dvec3 )

!*****************************************************************************80
!
!! DVEC_SUB subtracts two DVEC's.
!
!  Discussion:
!
!    A DVEC is an integer vector of decimal digits, intended to
!    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
!    is the coefficient of 10**(N-2), and DVEC(N) contains sign
!    information.  It is 0 if the number is positive, and 9 if
!    the number is negative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) DVEC1(N), DVEC2(N), the vectors to be
!    subtracted.
!
!    Output, integer ( kind = 4 ) DVEC3(N), the value of DVEC1 - DVEC2.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) dvec1(n)
  integer ( kind = 4 ) dvec2(n)
  integer ( kind = 4 ) dvec3(n)

  dvec3(1:n) = dvec2(1:n)
  call dvec_complementx ( n, dvec3, dvec3 )

  call dvec_add ( n, dvec1, dvec3, dvec3 )

  return
end
subroutine dvec_to_i4 ( n, dvec, i4 )

!*****************************************************************************80
!
!! DVEC_TO_I4 makes an integer from a (signed) DVEC.
!
!  Discussion:
!
!    A DVEC is an integer vector of decimal digits, intended to
!    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
!    is the coefficient of 10**(N-2), and DVEC(N) contains sign
!    information.  It is 0 if the number is positive, and 9 if
!    the number is negative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input, integer ( kind = 4 ) DVEC(N), the decimal vector.
!
!    Output, integer ( kind = 4 ) I4, the integer.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 10
  integer ( kind = 4 ) dvec(n)
  integer ( kind = 4 ) dvec2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_sign
  integer ( kind = 4 ) i4

  dvec2(1:n) = dvec(1:n)

  i_sign = 1

  if ( dvec2(n) == base - 1 ) then
    i_sign = -1
    call dvec_complementx ( n - 1, dvec2, dvec2 )
  end if

  i4 = 0
  do i = n-1, 1, -1
    i4 = base * i4 + dvec2(i)
  end do

  i4 = i_sign * i4

  return
end
subroutine equiv_next ( n, npart, jarray, iarray, more )

!*****************************************************************************80
!
!! EQUIV_NEXT computes the partitions of a set one at a time.
!
!  Discussion:
!
!    A partition of a set assigns each element to exactly one subset.
!
!    The number of partitions of a set of size N is the Bell number B(N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 July 2004
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
!    Input, integer ( kind = 4 ) N, the number of elements in the set to
!    be partitioned.
!
!    Input/output, integer ( kind = 4 ) NPART, the number of subsets in
!    the partition.
!
!    Input/output, integer ( kind = 4 ) JARRAY(N), the number of elements
!    in each subset of the partition.
!
!    Input/output, integer ( kind = 4 ) IARRAY(N), the subset to which
!    each element belongs.
!
!    Input/output, logical MORE.  Set MORE = FALSE before first call.
!    It is reset and held at TRUE as long as
!    the partition returned is not the last one.
!    When MORE is returned FALSE, all the partitions
!    have been computed and returned.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iarray(n)
  integer ( kind = 4 ) jarray(n)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) npart

  if ( .not. more ) then

    npart = 1
    iarray(1:n) = 1
    jarray(1) = n

  else

    m = n

    do while ( jarray(iarray(m)) == 1 )
      iarray(m) = 1
      m = m - 1
    end do

    l = iarray(m)
    npart = npart + m - n
    jarray(1) = jarray(1) + n - m

    if ( l == npart ) then
      npart = npart + 1
      jarray(npart) = 0
    end if

    iarray(m) = l + 1
    jarray(l) = jarray(l) - 1
    jarray(l+1) = jarray(l+1) + 1

  end if

  more = npart /= n

  return
end
subroutine equiv_next2 ( n, a, done )

!*****************************************************************************80
!
!! EQUIV_NEXT2 computes, one at a time, the partitions of a set.
!
!  Discussion:
!
!    A partition of a set assigns each element to exactly one subset.
!
!    The number of partitions of a set of size N is the Bell number B(N).
!
!    The entries of A are the partition subset to which each
!    element of the original set belongs.  If there are NPART distinct
!    parts of the partition, then each entry of A will be a
!    number between 1 and NPART.  Every number from 1 to NPART will
!    occur somewhere in the list.  If the entries of A are
!    examined in order, then each time a new partition subset occurs,
!    it will be the next unused integer.
!
!    For instance, for N = 4, the program will describe the set
!    where each element is in a separate subset as 1, 2, 3, 4,
!    even though such a partition might also be described as
!    4, 3, 2, 1 or even 1, 5, 8, 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the set.
!
!    Input/output, integer ( kind = 4 ) A(N), contains the information
!    defining the current partition.  The user should not alter
!    A between calls.  Except for the very first
!    call, the routine uses the previous output value of A to compute
!    the next value.
!
!    Input/output, logical DONE.  Before the very first call, the
!    user should set DONE to TRUE, which prompts the program
!    to initialize its data, and return the first partition.
!    Thereafter, the user should call again, for the next
!    partition, and so on, until the routine returns with DONE
!    equal to TRUE, at which point there are no more partitions
!    to compute.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  logical done
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jmax

  if ( done ) then

    done = .false.

    a(1:n) = 1

  else
!
!  Find the last element J that can be increased by 1.
!  This is the element that is not equal to its maximum possible value,
!  which is the maximum value of all preceding elements +1.
!
    jmax = a(1)
    imax = 1

    do j = 2, n

      if ( jmax < a(j) ) then
        jmax = a(j)
      else
        imax = j
      end if

    end do
!
!  If no element can be increased by 1, we are done.
!
    if ( imax == 1 ) then
      done = .true.
      return
    end if
!
!  Increase the value of the IMAX-th element by 1, set its successors to 1.
!
    done = .false.
    a(imax) = a(imax) + 1
    a(imax+1:n) = 1

  end if

  return
end
subroutine equiv_print ( n, a, title )

!*****************************************************************************80
!
!! EQUIV_PRINT prints a partition of a set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, number of elements in set to be partitioned.
!
!    Input, integer ( kind = 4 ) A(N), defines the partition or set of
!    equivalence classes.  Element I belongs to subset A(I).
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) karray(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) s
  integer ( kind = 4 ) s_max
  integer ( kind = 4 ) s_min
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Set  Size'

  s_min = minval ( a(1:n) )
  s_max = maxval ( a(1:n) )

  do s = s_min, s_max

    k = 0

    do j = 1, n

      if ( a(j) == s ) then
        k = k + 1
        karray(k) = j
      end if

    end do

    if ( 0 < k ) then
      write ( *, '(i8,i8,a,(10i4))' ) s, k, ' :: ', karray(1:k)
    end if

  end do

  return
end
subroutine equiv_print2 ( n, s, title )

!*****************************************************************************80
!
!! EQUIV_PRINT2 prints a partition of a set.
!
!  Discussion:
!
!    The partition is printed using the parenthesis format.
!
!    For example, here are the partitions of a set of 4 elements:
!
! ÊÊÊÊÊ(1,2,3,4)
! ÊÊÊÊÊ(1,2,3)(4)
! ÊÊÊÊÊ(1,2,4)(3)
! ÊÊÊÊÊ(1,2)(3,4)
! ÊÊÊÊÊ(1,2)(3)(4)
! ÊÊÊÊÊ(1,3,4)(2)
! ÊÊÊÊÊ(1,3)(2,4)
! ÊÊÊÊÊ(1,3)(2)(4)
! ÊÊÊÊÊ(1,4)(2,3)
! ÊÊÊÊÊ(1)(2,3,4)
! ÊÊÊÊÊ(1)(2,3)(4)
! ÊÊÊÊÊ(1,4)(2)(3)
! ÊÊÊÊÊ(1)(2,4)(3)
! ÊÊÊÊÊ(1)(2)(3,4)
! ÊÊÊÊÊ(1)(2)(3)(4)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the set.
!
!    Input, integer ( kind = 4 ) S(N), defines the partition.
!    Element I belongs to subset S(I).
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s(n)
  integer ( kind = 4 ) s_max
  integer ( kind = 4 ) s_min
  integer ( kind = 4 ) size
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  s_min = minval ( s(1:n) )
  s_max = maxval ( s(1:n) )

  do j = 1, s_max
    write ( *, '(a)', ADVANCE = 'NO' ) '('
    size = 0
    do i = 1, n
      if ( s(i) == j ) then
        if ( 0 < size ) then
          write ( *, '(a)', ADVANCE = 'NO' ) ','
        end if
        write ( *, '(i2)', ADVANCE = 'NO' ) i
        size = size + 1
      end if
    end do
    write ( *, '(a)', ADVANCE = 'NO' ) ')'

  end do

  write ( *, '(a)', ADVANCE = 'YES' )

  return
end
subroutine equiv_random ( n, seed, npart, a, b )

!*****************************************************************************80
!
!! EQUIV_RANDOM selects a random partition of a set.
!
!  Discussion:
!
!    The user does not control the number of parts in the partition.
!
!    The equivalence classes are numbered in no particular order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 2002
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
!    Input, integer ( kind = 4 ) N, the number of elements in the set to
!    be partitioned.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, integer ( kind = 4 ) NPART, the number of classes or parts in the
!    partition.  NPART will be between 1 and N.
!
!    Output, integer ( kind = 4 ) A(N), indicates the class to which each
!    element is assigned.
!
!    Output, real ( kind = 8 ) B(N).  B(K) = C(K)/(K!), where
!    C(K) = number of partitions of a set of K objects.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) npart
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) sum1
  real ( kind = 8 ) z

  b(1) = 1.0D+00

  do l = 1, n-1

    sum1 = 1.0D+00 / real ( l, kind = 8 )
    do k = 1, l-1
      sum1 = ( sum1 + b(k) ) / real ( l - k, kind = 8 )
    end do

    b(l+1) = ( sum1 + b(l) ) / real ( l + 1, kind = 8 )

  end do

  m = n
  npart = 0

  do

    z = r8_uniform_01 ( seed )
    z = real ( m, kind = 8 ) * b(m) * z
    k = 0
    npart = npart + 1

    do while ( 0.0D+00 <= z )

      a(m) = npart
      m = m - 1

      if ( m == 0 ) then
        exit
      end if

      z = z - b(m)
      k = k + 1
      z = z * k

    end do

    if ( m == 0 ) then
      exit
    end if

  end do
!
!  Randomly permute the assignments.
!
  call perm_random2 ( n, seed, a )

  return
end
subroutine euler ( n, ieuler )

!*****************************************************************************80
!
!! EULER returns the N-th row of Euler's triangle.
!
!  Discussion:
!
!    E(N,K) counts the number of permutations of the N digits that have
!    exactly K "ascents", that is, K places where the Ith digit is
!    less than the (I+1)th digit.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the row of Euler's triangle desired.
!
!    Output, integer ( kind = 4 ) IEULER(0:N), the N-th row of Euler's
!    triangle, IEULER(K) contains the value of E(N,K).  Note
!    that IEULER(0) should be 1 and IEULER(N) should be 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ieuler(0:n)
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) k

  ieuler(0) = 1

  if ( 0 < n ) then

    ieuler(1) = 0

    do irow = 2, n

      ieuler(irow) = 0

      do k = irow-1, 1, -1

        ieuler(k) = ( k + 1 ) * ieuler(k) + ( irow - k ) * ieuler(k-1)

      end do

      ieuler(0) = 1

    end do

  end if

  return
end
function frobenius_number_order2 ( c1, c2 )

!*****************************************************************************80
!
!! FROBENIUS_NUMBER_ORDER2 returns the Frobenius number for order 2.
!
!  Discussion:
!
!    The Frobenius number of order N is the solution of the Frobenius
!    coin sum problem for N coin denominations.
!
!    The Frobenius coin sum problem assumes the existence of
!    N coin denominations, and asks for the largest value that cannot
!    be formed by any combination of coins of these denominations.
!
!    The coin denominations are assumed to be distinct positive integers.
!
!    For general N, this problem is fairly difficult to handle.
!
!    For N = 2, it is known that:
!
!    * if C1 and C2 are not relatively prime, then
!      there are infinitely large values that cannot be formed.
!
!    * otherwise, the largest value that cannot be formed is
!      C1 * C2 - C1 - C2, and that exactly half the values between
!      1 and C1 * C2 - C1 - C2 + 1 cannot be represented.
!
!    As a simple example, if C1 = 2 and C2 = 7, then the largest
!    unrepresentable value is 5, and there are (5+1)/2 = 3
!    unrepresentable values, namely 1, 3, and 5.
!
!    For a general N, and a set of coin denominations C1, C2, ..., CN,
!    the Frobenius number F(N, C(1:N) ) is defined as the largest value
!    B for which the equation
!
!      C1*X1 + C2*X2 + ... + CN*XN = B
!
!    has no nonnegative integer solution X(1:N).
!
!    In the Mathematica Package "NumberTheory", the Frobenius number
!    can be determined by
!
!    <<NumberTheory`Frobenius`
!    FrobeniusF[ {C1,...,CN} ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Sylvester,
!    Question 7382,
!    Mathematical Questions with their Solutions,
!    Educational Times,
!    Volume 41, page 21, 1884.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) C1, C2, the coin denominations. C1 and C2
!    should be positive and relatively prime.
!
!    Output, integer ( kind = 4 ) FROBENIUS_NUMBER_ORDER2, the Frobenius
!    number of (C1,C2).
!
  implicit none

  integer ( kind = 4 ) c1
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) frobenius_number_order2
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) value

  if ( c1 <= 0 ) then
    value = i4_huge ( )
  else if ( c2 <= 0 ) then
    value = i4_huge ( )
  else if ( i4_gcd ( c1, c2 ) /= 1 ) then
    value = i4_huge ( )
  else
    value = c1 * c2 - c1 - c2
  end if

  frobenius_number_order2 = value

  return
end
subroutine frobenius_number_order2_values ( n_data, c1, c2, f )

!*****************************************************************************80
!
!! FROBENIUS_NUMBER_ORDER2_VALUES: values of the order 2 Frobenius number.
!
!  Discussion:
!
!    The Frobenius number of order N is the solution of the Frobenius
!    coin sum problem for N coin denominations.
!
!    The Frobenius coin sum problem assumes the existence of
!    N coin denominations, and asks for the largest value that cannot
!    be formed by any combination of coins of these denominations.
!
!    The coin denominations are assumed to be distinct positive integers.
!
!    For general N, this problem is fairly difficult to handle.
!
!    For N = 2, it is known that:
!
!    * if C1 and C2 are not relatively prime, then
!      there are infinitely large values that cannot be formed.
!
!    * otherwise, the largest value that cannot be formed is
!      C1 * C2 - C1 - C2, and that exactly half the values between
!      1 and C1 * C2 - C1 - C2 + 1 cannot be represented.
!
!    As a simple example, if C1 = 2 and C2 = 7, then the largest
!    unrepresentable value is 5, and there are (5+1)/2 = 3
!    unrepresentable values, namely 1, 3, and 5.
!
!    For a general N, and a set of coin denominations C1, C2, ..., CN,
!    the Frobenius number F(N, C(1:N) ) is defined as the largest value
!    B for which the equation
!
!      C1*X1 + C2*X2 + ... + CN*XN = B
!
!    has no nonnegative integer solution X(1:N).
!
!    In the Mathematica Package "NumberTheory", the Frobenius number
!    can be determined by
!
!    <<NumberTheory`Frobenius`
!    FrobeniusF[ {C1,...,CN} ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    James Sylvester,
!    Question 7382,
!    Mathematical Questions with their Solutions,
!    Educational Times,
!    Volume 41, page 21, 1884.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) C1, integer C2, the parameters of the
!    function.
!
!    Output, integer ( kind = 4 ) F, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 6

  integer ( kind = 4 ) c1
  integer ( kind = 4 ), save, dimension ( n_max ) :: c1_vec = (/ &
     2, &
     3, &
     4, &
     5, &
    12, &
    99 /)
  integer ( kind = 4 ) c2
  integer ( kind = 4 ), save, dimension ( n_max ) :: c2_vec = (/ &
      5, &
     17, &
     19, &
     13, &
     11, &
    100 /)
  integer ( kind = 4 ) f
  integer ( kind = 4 ), save, dimension ( n_max ) :: f_vec = (/ &
    3, &
    31, &
    23, &
    47, &
    109, &
    9701 /)
  integer ( kind = 4 ) n_data

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    c1 = 0
    c2 = 0
    f = 0
  else
    c1 = c1_vec(n_data)
    c2 = c2_vec(n_data)
    f = f_vec(n_data)
  end if

  return
end
subroutine gamma_log_values ( n_data, x, fx )

!*****************************************************************************80
!
!! GAMMA_LOG_VALUES returns some values of the Log Gamma function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Log[Gamma[x]]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Wolfram Media / Cambridge University Press, 1999.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1524063822430784D+01, &
     0.7966778177017837D+00, &
     0.3982338580692348D+00, &
     0.1520596783998375D+00, &
     0.0000000000000000D+00, &
    -0.4987244125983972D-01, &
    -0.8537409000331584D-01, &
    -0.1081748095078604D+00, &
    -0.1196129141723712D+00, &
    -0.1207822376352452D+00, &
    -0.1125917656967557D+00, &
    -0.9580769740706586D-01, &
    -0.7108387291437216D-01, &
    -0.3898427592308333D-01, &
    0.00000000000000000D+00, &
    0.69314718055994530D+00, &
    0.17917594692280550D+01, &
    0.12801827480081469D+02, &
    0.39339884187199494D+02, &
    0.71257038967168009D+02 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.20D+00, &
     0.40D+00, &
     0.60D+00, &
     0.80D+00, &
     1.00D+00, &
     1.10D+00, &
     1.20D+00, &
     1.30D+00, &
     1.40D+00, &
     1.50D+00, &
     1.60D+00, &
     1.70D+00, &
     1.80D+00, &
     1.90D+00, &
     2.00D+00, &
     3.00D+00, &
     4.00D+00, &
    10.00D+00, &
    20.00D+00, &
    30.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine get_seed ( seed )

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
!    17 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) SEED, a pseudorandom seed value.
!
  implicit none

  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) seed
  real ( kind = 8 ) temp
  character ( len = 10 ) time
  character ( len = 8 ) today
  integer ( kind = 4 ) values(8)
  character ( len = 5 ) zone

  call date_and_time ( today, time, zone, values )

  temp = 0.0D+00

  temp = temp + real ( values(2) - 1, kind = 8 ) / 11.0D+00
  temp = temp + real ( values(3) - 1, kind = 8 ) / 30.0D+00
  temp = temp + real ( values(5),     kind = 8 ) / 23.0D+00
  temp = temp + real ( values(6),     kind = 8 ) / 59.0D+00
  temp = temp + real ( values(7),     kind = 8 ) / 59.0D+00
  temp = temp + real ( values(8),     kind = 8 ) / 999.0D+00
  temp = temp / 6.0D+00
!
!  Force 0 < TEMP <= 1.
!
  do while ( temp <= 0.0D+00 )
    temp = temp + 1.0D+00
  end do

  do while ( 1.0D+00 < temp )
    temp = temp - 1.0D+00
  end do

  seed = int ( real ( i4_huge ( ), kind = 8 ) * temp )
!
!  Never use a seed of 0 or maximum integer.
!
  if ( seed == 0 ) then
    seed = 1
  end if

  if ( seed == i4_huge ( ) ) then
    seed = seed - 1
  end if

  return
end
subroutine gray_next ( n, change )

!*****************************************************************************80
!
!! GRAY_NEXT generates the next Gray code by switching one item at a time.
!
!  Discussion:
!
!    On the first call only, the user must set CHANGE = -N.
!    This initializes the routine to the Gray code for N zeroes.
!
!    Each time it is called thereafter, it returns in CHANGE the index
!    of the item to be switched in the Gray code.  The sign of CHANGE
!    indicates whether the item is to be added or subtracted (or
!    whether the corresponding bit should become 1 or 0).  When
!    CHANGE is equal to N+1 on output, all the Gray codes have been
!    generated.
!
!    The routine has internal memory that is set up on call with
!    CHANGE = -N, and released on final return.
!
!  Example:
!
!    N  CHANGE         Subset in/out   Binary Number
!                      Interpretation  Interpretation
!                       1 2 4 8
!   --  ---------      --------------  --------------
!
!    4   -4 / 0         0 0 0 0         0
!
!        +1             1 0 0 0         1
!          +2           1 1 0 0         3
!        -1             0 1 0 0         2
!            +3         0 1 1 0         6
!        +1             1 1 1 0         7
!          -2           1 0 1 0         5
!        -1             0 0 1 0         4
!              +4       0 0 1 1        12
!        +1             1 0 1 1        13
!          +2           1 1 1 1        15
!        -1             0 1 1 1        14
!            -3         0 1 0 1        10
!        +1             1 1 0 1        11
!          -2           1 0 0 1         9
!        -1             0 0 0 1         8
!              -4       0 0 0 0         0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 April 2003
!
!  Author:
!
!    John Burkardt
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
!    Input, integer ( kind = 4 ) N, the order of the total set from which
!    subsets will be drawn.
!
!    Input/output, integer ( kind = 4 ) CHANGE.  This is an input item only
!    on the first call for a particular sequence of Gray codes,
!    at which time it must be set to -N.  This corresponds to
!    all items being excluded, or all bits being 0, in the Gray code.
!    On output, CHANGE indicates which of the N items must be "changed",
!    and the sign indicates whether the item is to be added or removed
!    (or the bit is to become 1 or 0).  On return from the
!    first call, CHANGE is set to 0, indicating that the first set
!    is the empty set.
!
  implicit none

  integer ( kind = 4 ), save, allocatable, dimension ( : ) :: a
  integer ( kind = 4 ) change
  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n_save = 0

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAY_NEXT - Fatal error!'
    write ( *, '(a)' ) '  Input value of N <= 0.'
    stop
  end if

  if ( change == -n ) then

    if ( allocated ( a ) ) then
      deallocate ( a )
    end if

    allocate ( a(1:n) )
    a(1:n) = 0

    n_save = n
    k = 1
    change = 0

    return
  end if

  if ( n /= n_save ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAY_NEXT - Fatal error!'
    write ( *, '(a)' ) '  Input value of N has changed from definition value.'
    stop
  end if
!
!  First determine WHICH item is to be changed.
!
  if ( mod ( k, 2 ) == 1 ) then

    change = 1

  else

    do i = 1, n_save
      if ( a(i) == 1 ) then
        change = i + 1
        exit
      end if
    end do

  end if
!
!  Take care of the terminal case CHANGE = N + 1.
!
  if ( change == n + 1 ) then
    change = n
  end if
!
!  Now determine HOW the item is to be changed.
!
  if ( a(change) == 0 ) then
    a(change) = 1
  else if ( a(change) == 1 ) then
    a(change) = 0
    change = - change
  end if
!
!  Update the counter.
!
  k = k + 1
!
!  If the output CHANGE = -N_SAVE, then we're done.
!
  if ( change == - n_save ) then

    deallocate ( a )
    n_save = 0
    k = 0

  end if

  return
end
subroutine gray_rank ( gray, rank )

!*****************************************************************************80
!
!! GRAY_RANK ranks a Gray code.
!
!  Discussion:
!
!    Given the number GRAY, its ranking is the order in which it would be
!    visited in the Gray code ordering.  The Gray code ordering begins
!
!    Rank  Gray  Gray
!          (Dec) (Bin)
!
!       0     0  0000
!       1     1  0001
!       2     3  0011
!       3     2  0010
!       4     6  0110
!       5     7  0111
!       6     5  0101
!       7     4  0100
!       8    12  0110
!       etc
!
!   This routine is given a Gray code, and has to return the rank.
!
!  Example:
!
!    Gray  Gray  Rank
!    (Dec) (Bin)
!
!     0       0     0
!     1       1     1
!     2      10     3
!     3      11     2
!     4     100     7
!     5     101     6
!     6     110     4
!     7     111     5
!     8    1000    15
!     9    1001    14
!    10    1010    12
!    11    1011    13
!    12    1100     8
!    13    1101     9
!    14    1110    11
!    15    1111    10
!    16   10000    31
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) GRAY, the Gray code to be ranked.
!
!    Output, integer ( kind = 4 ) RANK, the rank of GRAY, and the integer
!    whose Gray code is GRAY.
!
  implicit none

  integer ( kind = 4 ) gray
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: nbits = 32
  integer ( kind = 4 ) rank

  j = 0

  if ( btest ( gray, nbits-1 ) ) then
    j = ibset ( j, nbits-1 )
  end if

  do i = nbits-2, 0, -1

    if ( btest ( j, i+1 ) .and. ( .not. btest ( gray, i ) ) ) then
      j = ibset ( j, i )
    end if

    if ( ( .not. btest ( j, i+1 ) ) .and. btest ( gray, i ) ) then
      j = ibset ( j, i )
    end if

  end do

  rank = j

  return
end
subroutine gray_rank2 ( gray, rank )

!*****************************************************************************80
!
!! GRAY_RANK2 ranks a Gray code.
!
!  Discussion:
!
!    In contrast to GRAY_RANK, this routine is entirely arithmetical,
!    and does not require access to bit testing and setting routines.
!
!    Given the number GRAY, its ranking is the order in which it would be
!    visited in the Gray code ordering.  The Gray code ordering begins
!
!    Rank  Gray  Gray
!          (Dec) (Bin)
!
!       0     0  0000
!       1     1  0001
!       2     3  0011
!       3     2  0010
!       4     6  0110
!       5     7  0111
!       6     5  0101
!       7     4  0100
!       8    12  0110
!       etc
!
!   This routine is given a Gray code, and has to return the rank.
!
!  Example:
!
!    Gray  Gray  Rank
!    (Dec) (Bin)
!
!     0       0     0
!     1       1     1
!     2      10     3
!     3      11     2
!     4     100     7
!     5     101     6
!     6     110     4
!     7     111     5
!     8    1000    15
!     9    1001    14
!    10    1010    12
!    11    1011    13
!    12    1100     8
!    13    1101     9
!    14    1110    11
!    15    1111    10
!    16   10000    31
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) GRAY, the Gray code to be ranked.
!
!    Output, integer ( kind = 4 ) RANK, the rank of GRAY, and the integer
!    whose Gray code is GRAY.
!
  implicit none

  integer ( kind = 4 ) gray
  integer ( kind = 4 ) gray_copy
  integer ( kind = 4 ) k
  logical last
  logical next
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) two_k

  gray_copy = gray

  if ( gray_copy < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAY_RANK2 - Fatal error!'
    write ( *, '(a)' ) '  Input value of GRAY < 0.'
    stop
  end if

  if ( gray_copy == 0 ) then
    rank = 0
    return
  end if
!
!  Find TWO_K, the largest power of 2 less than or equal to GRAY.
!
  k = 0
  two_k = 1
  do while ( 2 * two_k <= gray_copy )
    two_k = two_k * 2
    k = k + 1
  end do

  rank = two_k
  last = .true.
  gray_copy = gray_copy - two_k

  do while ( 0 < k )

    two_k = two_k / 2
    k = k - 1

    next = ( two_k <= gray_copy .and. gray_copy < two_k * 2 )

    if ( next ) then
      gray_copy = gray_copy - two_k
    end if

    if ( next .neqv. last ) then
      rank = rank + two_k
      last = .true.
    else
      last = .false.
    end if

  end do

  return
end
subroutine gray_unrank ( rank, gray )

!*****************************************************************************80
!
!! GRAY_UNRANK unranks a Gray code.
!
!  Discussion:
!
!    The binary values of the Gray codes of successive integers differ in
!    just one bit.
!
!    The sequence of Gray codes for 0 to (2**N)-1 can be interpreted as a
!    Hamiltonian cycle on a graph of the cube in N dimensions.
!
!  Example:
!
!    Rank  Gray  Gray
!          (Dec) (Bin)
!
!     0     0       0
!     1     1       1
!     2     3      11
!     3     2      10
!     4     6     110
!     5     7     111
!     6     5     101
!     7     4     100
!     8    12    1100
!     9    14    1001
!    10    12    1010
!    11    13    1011
!    12     8    1100
!    13     9    1101
!    14    11    1110
!    15    10    1111
!    16    31   10000
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the integer whose Gray code is desired.
!
!    Output, integer ( kind = 4 ) GRAY, the Gray code of the given rank.
!
  implicit none

  integer ( kind = 4 ) gray
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: nbits = 32
  integer ( kind = 4 ) rank

  j = 0

  if ( btest ( rank, nbits-1 ) ) then
    j = ibset ( j, nbits-1 )
  end if

  do i = nbits-2, 0, -1

    if ( btest ( rank, i+1 ) .neqv. btest ( rank, i ) ) then
      j = ibset ( j, i )
    end if

  end do

  gray = j

  return
end
subroutine gray_unrank2 ( rank, gray )

!*****************************************************************************80
!
!! GRAY_UNRANK2 unranks a Gray code.
!
!  Discussion:
!
!    In contrast to GRAY_UNRANK, this routine is entirely arithmetical,
!    and does not require access to bit testing and setting routines.
!
!    The binary values of the Gray codes of successive integers differ in
!    just one bit.
!
!    The sequence of Gray codes for 0 to (2**N)-1 can be interpreted as a
!    Hamiltonian cycle on a graph of the cube in N dimensions.
!
!  Example:
!
!    Rank  Gray  Gray
!          (Dec) (Bin)
!
!     0     0       0
!     1     1       1
!     2     3      11
!     3     2      10
!     4     6     110
!     5     7     111
!     6     5     101
!     7     4     100
!     8    12    1100
!     9    14    1001
!    10    12    1010
!    11    13    1011
!    12     8    1100
!    13     9    1101
!    14    11    1110
!    15    10    1111
!    16    31   10000
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the integer whose Gray code is desired.
!
!    Output, integer ( kind = 4 ) GRAY, the Gray code of the given rank.
!
  implicit none

  integer ( kind = 4 ) gray
  integer ( kind = 4 ) k
  logical last
  logical next
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_copy
  integer ( kind = 4 ) two_k

  if ( rank <= 0 ) then
    gray = 0
    return
  end if

  rank_copy = rank
  k = 0
  two_k = 1
  do while ( 2 * two_k <= rank_copy )
    two_k = two_k * 2
    k = k + 1
  end do

  gray = two_k
  rank_copy = rank_copy - two_k
  next = .true.

  do while ( 0 < k )

    two_k = two_k / 2
    k = k - 1

    last = next
    next = ( two_k <= rank_copy .and. rank_copy <= two_k * 2 )

    if ( next .neqv. last ) then
      gray = gray + two_k
    end if

    if ( next ) then
      rank_copy = rank_copy - two_k
    end if

  end do

  return
end
function i4_bclr ( i4, pos )

!*****************************************************************************80
!
!! I4_BCLR returns a copy of an I4 in which the POS-th bit is set to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Military Standard 1753,
!    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
!    9 November 1978.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I4, the integer to be tested.
!
!    Input, integer ( kind = 4 ) POS, the bit position, between 0 and 31.
!
!    Output, integer ( kind = 4 ) I4_BCLR, a copy of I4, but with the POS-th bit
!    set to 0.
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_bclr
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) pos
  integer ( kind = 4 ) sub
  integer ( kind = 4 ) value

  value = i4

  if ( pos < 0 ) then

  else if ( pos < 31 ) then

    sub = 1

    if ( 0 <= i4 ) then
      j = i4
    else
      j = ( i4_huge ( ) + i4 ) + 1
    end if

    do k = 1, pos
      j = j / 2
      sub = sub * 2
    end do

    if ( mod ( j, 2 ) == 1 ) then
      value = i4 - sub
    end if

  else if ( pos == 31 ) then

    if ( i4 < 0 ) then
      value = ( i4_huge ( ) + i4 ) + 1
    end if

  else if ( 31 < pos ) then

    value = i4

  end if

  i4_bclr = value

  return
end
function i4_bset ( i4, pos )

!*****************************************************************************80
!
!! I4_BSET returns a copy of an I4 in which the POS-th bit is set to 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Military Standard 1753,
!    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
!    9 November 1978.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I4, the integer to be tested.
!
!    Input, integer ( kind = 4 ) POS, the bit position, between 0 and 31.
!
!    Output, integer ( kind = 4 ) I4_BSET, a copy of I4, but with the POS-th bit
!    set to 1.
!
  implicit none

  integer ( kind = 4 ) add
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_bset
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) pos
  integer ( kind = 4 ) value

  value = i4

  if ( pos < 0 ) then

  else if ( pos < 31 ) then

    add = 1

    if ( 0 <= i4 ) then
      j = i4
    else
      j = ( i4_huge ( ) + i4 ) + 1
    end if

    do k = 1, pos
      j = j / 2
      add = add * 2
    end do

    if ( mod ( j, 2 ) == 0 ) then
      value = i4 + add
    end if

  else if ( pos == 31 ) then

    if ( 0 < i4 ) then
      value = - ( i4_huge ( ) - i4 ) - 1
    end if

  else if ( 31 < pos ) then

    value = i4

  end if

  i4_bset = value

  return
end
function i4_btest ( i4, pos )

!*****************************************************************************80
!
!! I4_BTEST returns TRUE if the POS-th bit of an I4 is 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Military Standard 1753,
!    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
!    9 November 1978.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I4, the integer to be tested.
!
!    Input, integer ( kind = 4 ) POS, the bit position, between 0 and 31.
!
!    Output, logical I4_BTEST, is TRUE if the POS-th bit of I4 is 1.
!
  implicit none

  integer ( kind = 4 ) i4
  logical i4_btest
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) pos

  if ( pos < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_BTEST - Fatal error!'
    write ( *, '(a)' ) '  POS < 0.'
    stop

  else if ( pos < 31 ) then

    if ( 0 <= i4 ) then
      j = i4
    else
      j = ( i4_huge ( ) + i4 ) + 1
    end if

    do k = 1, pos
      j = j / 2
    end do

    if ( mod ( j, 2 ) == 0 ) then
      i4_btest = .false.
    else
      i4_btest = .true.
    end if

  else if ( pos == 31 ) then

    if ( i4 < 0 ) then
      i4_btest = .true.
    else
      i4_btest = .false.
    end if

  else if ( 31 < pos ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_BTEST - Fatal error!'
    write ( *, '(a)' ) '  31 < POS.'
    stop

  end if

  return
end
function i4_choose ( n, k )

!*****************************************************************************80
!
!! I4_CHOOSE computes the binomial coefficient C(N,K).
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, K, are the values of N and K.
!
!    Output, integer ( kind = 4 ) I4_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0

  else if ( mn == 0 ) then

    value = 1

  else

    mx = max ( k, n - k )
    value = mx + 1

    do i = 2, mn
      value = ( value * ( mx + i ) ) / i
    end do

  end if

  i4_choose = value

  return
end
subroutine i4_factor ( n, factor_max, factor_num, factor, power, nleft )

!*****************************************************************************80
!
!! I4_FACTOR factors an I4 into prime factors.
!
!  Discussion:
!
!    The factorization has the form
!
!      N = NLEFT * Product ( 1 <= I <= FACTOR_NUM ) FACTOR(I)**POWER(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be factored.  N may be
!    positive, negative, or 0.
!
!    Input, integer ( kind = 4 ) FACTOR_MAX, the maximum number of prime
!    factors for which storage has been allocated.
!
!    Output, integer ( kind = 4 ) FACTOR_NUM, the number of prime factors
!    of N discovered by the routine.
!
!    Output, integer ( kind = 4 ) FACTOR(FACTOR_MAX), the prime factors of N.
!
!    Output, integer ( kind = 4 ) POWER(FACTOR_MAX).  POWER(I) is the power of
!    the FACTOR(I) in the representation of N.
!
!    Output, integer ( kind = 4 ) NLEFT, the factor of N that the routine
!    could not divide out.  If NLEFT is 1, then N has been completely factored.
!    Otherwise, NLEFT represents factors of N involving large primes.
!
  implicit none

  integer ( kind = 4 ) factor_max

  integer ( kind = 4 ) factor(factor_max)
  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nleft
  integer ( kind = 4 ) p
  integer ( kind = 4 ) power(factor_max)
  integer ( kind = 4 ) prime
  integer ( kind = 4 ) prime_max

  factor_num = 0

  factor(1:factor_max) = 0
  power(1:factor_max) = 0

  nleft = n

  if ( n == 0 ) then
    return
  end if

  if ( abs ( n ) == 1 ) then
    factor_num = 1
    factor(1) = 1
    power(1) = 1
    return
  end if
!
!  Find out how many primes we stored.
!
  prime_max = prime ( -1 )
!
!  Try dividing the remainder by each prime.
!
  do i = 1, prime_max

    p = prime ( i )

    if ( mod ( abs ( nleft ), p ) == 0 ) then

      if ( factor_num < factor_max ) then

        factor_num = factor_num + 1
        factor(factor_num) = p

        do

          power(factor_num) = power(factor_num) + 1
          nleft = nleft / p

          if ( mod ( abs ( nleft ), p ) /= 0 ) then
            exit
          end if

        end do

        if ( abs ( nleft ) == 1 ) then
          exit
        end if

        if ( factor_max <= factor_num ) then
          exit
        end if

      end if

    end if

  end do

  return
end
function i4_factorial ( n )

!*****************************************************************************80
!
!! I4_FACTORIAL computes the factorial of N.
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
!    26 June 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!    0 <= N <= 13 is required.
!
!    Output, integer ( kind = 4 ) I4_FACTORIAL, the factorial of N.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_factorial
  integer ( kind = 4 ) n

  i4_factorial = 1

  if ( 13 < n ) then
    i4_factorial = - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_FACTORIAL - Fatal error!'
    write ( *, '(a)' ) '  I4_FACTORIAL(N) cannot be computed as an integer'
    write ( *, '(a)' ) '  for 13 < N.'
    write ( *, '(a,i8)' ) '  Input value N = ', n
    stop
  end if

  do i = 1, n
    i4_factorial = i4_factorial * i
  end do

  return
end
function i4_fall ( x, n )

!*****************************************************************************80
!
!! I4_FALL computes the falling factorial function [X]_N.
!
!  Discussion:
!
!    Note that the number of "injections" or 1-to-1 mappings from
!    a set of N elements to a set of M elements is [M]_N.
!
!    The number of permutations of N objects out of M is [M]_N.
!
!    Moreover, the Stirling numbers of the first kind can be used
!    to convert a falling factorial into a polynomial, as follows:
!
!      [X]_N = S^0_N + S^1_N * X + S^2_N * X^2 + ... + S^N_N X^N.
!
!    The formula used is:
!
!      [X]_N = X * ( X - 1 ) * ( X - 2 ) * ... * ( X - N + 1 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the falling factorial
!    function.
!
!    Input, integer ( kind = 4 ) N, the order of the falling factorial function.
!    If N = 0, FALL = 1, if N = 1, FALL = X.  Note that if N is
!    negative, a "rising" factorial will be computed.
!
!    Output, integer ( kind = 4 ) I4_FALL, the value of the falling
!    factorial function.
!
  implicit none

  integer ( kind = 4 ) arg
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_fall
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value
  integer ( kind = 4 ) x

  value = 1

  arg = x

  if ( 0 < n ) then

    do i = 1, n
      value = value * arg
      arg = arg - 1
    end do

  else if ( n < 0 ) then

    do i = -1, n, -1
      value = value * arg
      arg = arg + 1
    end do

  end if

  i4_fall = value

  return
end
function i4_gcd ( i, j )

!*****************************************************************************80
!
!! I4_GCD finds the greatest common divisor of two I4's.
!
!  Discussion:
!
!    Only the absolute values of I and J are
!    considered, so that the result is always nonnegative.
!
!    If I or J is 0, I4_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
!
!    If I and J have no common factor, I4_GCD is returned as 1.
!
!    Otherwise, using the Euclidean algorithm, I4_GCD is the
!    largest common factor of I and J.
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
!    Output, integer ( kind = 4 ) I4_GCD, the greatest common divisor of
!    I and J.
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
function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) I4_HUGE, a "huge" integer.
!
  implicit none

  integer ( kind = 4 ) i4_huge

  i4_huge = 2147483647

  return
end
function i4_log_10 ( i )

!*****************************************************************************80
!
!! I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
!
!  Discussion:
!
!    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
!
!  Example:
!
!        I  I4_LOG_10
!    -----  --------
!        0    0
!        1    0
!        2    0
!        9    0
!       10    1
!       11    1
!       99    1
!      100    2
!      101    2
!      999    2
!     1000    3
!     1001    3
!     9999    3
!    10000    4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose logarithm base 10 is
!    desired.
!
!    Output, integer ( kind = 4 ) I4_LOG_10, the integer part of the logarithm
!    base 10 of the absolute value of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_abs
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ) ten_pow

  if ( i == 0 ) then

    i4_log_10 = 0

  else

    i4_log_10 = 0
    ten_pow = 10

    i_abs = abs ( i )

    do while ( ten_pow <= i_abs )
      i4_log_10 = i4_log_10 + 1
      ten_pow = ten_pow * 10
    end do

  end if

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    The computation may be summarized as:
!
!      Let
!        NREM = I4_MODP ( I, J )
!      and
!        NMULT = ( I - NREM ) / J
!      then
!        I = J * NMULT + NREM
!
!      where NREM is always nonnegative.
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
!        I         J     MOD  I4_MODP   I4_MODP Factorization
!
!      107        50       7       7    107 =  2 *  50 + 7
!      107       -50       7       7    107 = -2 * -50 + 7
!     -107        50      -7      43   -107 = -3 *  50 + 43
!     -107       -50      -7      43   -107 =  3 * -50 + 43
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
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder
!    when I is divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

  return
end
subroutine i4_moebius ( n, mu )

!*****************************************************************************80
!
!! I4_MOEBIUS returns the value of MU(N), the Moebius function of N.
!
!  Discussion:
!
!    MU(N) is defined as follows:
!
!      MU(N) = 1 if N = 1;
!              0 if N is divisible by the square of a prime;
!              (-1)**K, if N is the product of K distinct primes.
!
!    As special cases, MU(N) is -1 if N is a prime, and MU(N) is 0
!    if N is a square, cube, etc.
!
!    The Moebius function MU(D) is related to Euler's totient
!    function PHI(N):
!
!      PHI(N) = sum ( D divides N ) MU(D) * ( N / D ).
!
!  First values:
!
!     N  MU(N)
!
!     1    1
!     2   -1
!     3   -1
!     4    0
!     5   -1
!     6    1
!     7   -1
!     8    0
!     9    0
!    10    1
!    11   -1
!    12    0
!    13   -1
!    14    1
!    15    1
!    16    0
!    17   -1
!    18    0
!    19   -1
!    20    0
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
!    Input, integer ( kind = 4 ) N, the value to be analyzed.
!
!    Output, integer ( kind = 4 ) MU, the value of MU(N).
!    If N is less than or equal to 0, MU will be returned as -2.
!    If there was not enough internal space for factoring, MU
!    is returned as -3.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxfactor = 20

  integer ( kind = 4 ) exponent(maxfactor)
  integer ( kind = 4 ) factor(maxfactor)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nfactor
  integer ( kind = 4 ) nleft

  if ( n <= 0 ) then
    mu = -2
    return
  end if

  if ( n == 1 ) then
    mu = 1
    return
  end if
!
!  Factor N.
!
  call i4_factor ( n, maxfactor, nfactor, factor, exponent, nleft )

  if ( nleft /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MOEBIUS - Fatal error!'
    write ( *, '(a)' ) '  Not enough factorization space.'
    mu = -3
    stop
  end if

  mu = 1

  do i = 1, nfactor

    mu = -mu

    if ( 1 < exponent(i) ) then
      mu = 0
      return
    end if

  end do

  return
end
subroutine i4_partition_conj ( n, iarray1, mult1, npart1, iarray2, mult2, &
  npart2 )

!*****************************************************************************80
!
!! I4_PARTITION_CONJ computes the conjugate of a partition.
!
!  Discussion:
!
!    A partition of an integer N is a representation of the integer
!    as the sum of nonzero positive integers.  The order of the summands
!    does not matter.  Thus, the number 5 has the following partitions
!    and no more:
!
!    5 = 5
!      = 4 + 1
!      = 3 + 2
!      = 3 + 1 + 1
!      = 2 + 2 + 1
!      = 2 + 1 + 1 + 1
!      = 1 + 1 + 1 + 1 + 1
!
!    so the number of partitions of 5 is 7.
!
!    The conjugate of a partition P1 of N is another partition P2 of N
!    obtained in the following way:
!
!      The first element of P2 is the number of parts of P1 greater than
!      or equal to 1.
!
!      The K-th element of P2 is the number of parts of P1 greater than
!      or equal to K.
!
!    Clearly, P2 will have no more than N elements; it may be surprising
!    to find that P2 is guaranteed to be a partition of N.  However, if
!    we symbolize the initial partition P1 by rows of X's, then we can
!    see that P2 is simply produced by grouping by columns:
!
!        6 3 2 2 1
!      5 X X X X X
!      4 X X X X
!      2 X X
!      1 X
!      1 X
!      1 X
!
!  Example:
!
!    14 = 5 + 4 + 2 + 1 + 1 + 1
!
!    The conjugate partition is:
!
!    14 = 6 + 3 + 2 + 2 + 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!
!    Input, integer ( kind = 4 ) IARRAY1(NPART1).  IARRAY1 contains the parts of
!    the partition.  The value of N is represented by
!      sum ( 1 <= I <= NPART1 ) MULT1(I) * IARRAY1(I).
!
!    Input, integer ( kind = 4 ) MULT1(NPART1).  MULT1 counts the multiplicity
!    of the parts of the partition.  MULT1(I) is the multiplicity
!    of the part IARRAY1(I), for I = 1 to NPART1.
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!
!    Input, integer ( kind = 4 ) NPART1, the number of "parts" in the partition.
!
!    Output, integer ( kind = 4 ) IARRAY2(N).  IARRAY contains the parts of
!    the conjugate partition in entries 1 through NPART2.
!
!    Output, integer ( kind = 4 ) MULT2(N).  MULT2 counts the multiplicity of
!    the parts of the conjugate partition in entries 1 through NPART2.
!
!    Output, integer ( kind = 4 ) NPART2, the number of "parts" in the
!    conjugate partition.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) npart1

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarray1(npart1)
  integer ( kind = 4 ) iarray2(n)
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) mult1(npart1)
  integer ( kind = 4 ) mult2(n)
  integer ( kind = 4 ) npart2

  iarray2(1:n) = 0
  mult2(1:n) = 0
  npart2 = 0

  itest = 0

  do

    itest = itest + 1

    itemp = 0

    do i = 1, npart1
      if ( itest <= iarray1(i) ) then
        itemp = itemp + mult1(i)
      end if
    end do

    if ( itemp <= 0 ) then
      exit
    end if

    if ( 0 < npart2 ) then
      if ( itemp == iarray2(npart2) ) then
        mult2(npart2) = mult2(npart2) + 1
      else
        npart2 = npart2 + 1
        iarray2(npart2) = itemp
        mult2(npart2) = 1
      end if
    else
      npart2 = npart2 + 1
      iarray2(npart2) = itemp
      mult2(npart2) = 1
    end if

  end do

  return
end
subroutine i4_partition_count ( n, p )

!*****************************************************************************80
!
!! I4_PARTITION_COUNT computes the number of partitions of an I4.
!
!  Discussion:
!
!    A partition of an integer N is a representation of the integer
!    as the sum of nonzero positive integers.  The order of the summands
!    does not matter.  Thus, the number 5 has the following partitions
!    and no more:
!
!    5 = 5
!      = 4 + 1
!      = 3 + 2
!      = 3 + 1 + 1
!      = 2 + 2 + 1
!      = 2 + 1 + 1 + 1
!      = 1 + 1 + 1 + 1 + 1
!
!    so the number of partitions of 5 is 7.
!
!    Partition numbers are difficult to compute.  This routine uses
!    Euler's method, which observes that:
!
!      P(0) = 1
!      P(N) =   P(N-1)  + P(N-2)
!             - P(N-5)  - P(N-7)
!             + P(N-12) + P(N-15)
!             - ...
!
!      where the numbers 1, 2, 5, 7, ... to be subtracted from N in the
!      indices are the successive pentagonal numbers, (with both positive
!      and negative indices) with the summation stopping when a negative
!      index is reached.
!
!  First values:
!
!    N   P
!
!    0   1
!    1   1
!    2   2
!    3   3
!    4   5
!    5   7
!    6  11
!    7  15
!    8  22
!    9  30
!   10  42
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
!    John Conway, Richard Guy,
!    The Book of Numbers,
!    Springer, 1998,
!    ISBN: 038797993X.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the highest partition
!    number desired.
!
!    Output, integer ( kind = 4 ) P(0:N), the partition numbers.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(0:n)
  integer ( kind = 4 ) pj
  integer ( kind = 4 ) sgn

  p(0) = 1

  do i = 1, n

    p(i) = 0

    j = 0
    sgn = 1

    do

      j = j + 1
      call pent_enum ( j, pj )

      if ( i < pj ) then
        exit
      end if

      p(i) = p(i) + sgn * p(i-pj)
      sgn = -sgn

    end do

    j = 0
    sgn = 1

    do

      j = j - 1
      call pent_enum ( j, pj )

      if ( i < pj ) then
        exit
      end if

      p(i) = p(i) + sgn * p(i-pj)
      sgn = -sgn

    end do

  end do

  return
end
subroutine i4_partition_count2 ( n, p )

!*****************************************************************************80
!
!! I4_PARTITION_COUNT2 computes the number of partitions of an I4.
!
!  Discussion:
!
!    A partition of an integer N is a representation of the integer
!    as the sum of nonzero positive integers.  The order of the summands
!    does not matter.  Thus, the number 5 has the following partitions
!    and no more:
!
!    5 = 5
!      = 4 + 1
!      = 3 + 2
!      = 3 + 1 + 1
!      = 2 + 2 + 1
!      = 2 + 1 + 1 + 1
!      = 1 + 1 + 1 + 1 + 1
!
!    so the number of partitions of 5 is 7.
!
!  First values:
!
!    N   P
!
!    0   1
!    1   1
!    2   2
!    3   3
!    4   5
!    5   7
!    6  11
!    7  15
!    8  22
!    9  30
!   10  42
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2004
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
!    Input, integer ( kind = 4 ) N, the largest integer to be considered.
!
!    Output, integer ( kind = 4 ) P(0:N), the partition numbers.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(0:n)
  integer ( kind = 4 ) s
  integer ( kind = 4 ) t
  integer ( kind = 4 ) total

  if ( n < 0 ) then
    return
  end if

  p(0) = 1

  if ( n < 1 ) then
    return
  end if

  p(1) = 1

  do i = 2, n

    total = 0

    do t = 1, i

      s = 0
      j = i

      do

        j = j - t

        if ( j < 0 ) then
          exit
        end if

        s = s + p(j)

      end do

      total = total + s * t

    end do

    p(i) = total / i

  end do

  return
end
subroutine i4_partition_count_values ( n_data, n, c )

!*****************************************************************************80
!
!! I4_PARTITION_COUNT_VALUES returns some values of the integer partition count.
!
!  Discussion:
!
!    A partition of an integer N is a representation of the integer
!    as the sum of nonzero positive integers.  The order of the summands
!    does not matter.  Thus, the number 5 has the following partitions
!    and no more:
!
!    5 = 5
!      = 4 + 1
!      = 3 + 2
!      = 3 + 1 + 1
!      = 2 + 2 + 1
!      = 2 + 1 + 1 + 1
!      = 1 + 1 + 1 + 1 + 1
!
!    so the number of partitions of 5 is 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.
!    On input, if N_DATA is 0, the first test data is returned, and N_DATA
!    is set to 1.  On each subsequent call, the input value of N_DATA is
!    incremented and that test data item is returned, if available.  When
!    there is no more test data, N_DATA is set to 0.
!
!    Output, integer ( kind = 4 ) N, the integer.
!
!    Output, integer ( kind = 4 ) C, the number of partitions of the integer.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 21

  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( nmax ) :: c_vec = (/ &
      1, &
      1,   2,   3,   5,   7,  11,  15,  22,  30,  42, &
     56,  77, 101, 135, 176, 231, 297, 385, 490, 627 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( nmax ) :: n_vec = (/ &
     0,  &
     1,  2,  3,  4,  5,  6,  7,  8,  9, 10, &
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( nmax < n_data ) then
    n_data = 0
    n = 0
    c = 0
  else
    n = n_vec(n_data)
    c = c_vec(n_data)
  end if

  return
end
subroutine i4_partition_next ( n, npart, a, mult, done )

!*****************************************************************************80
!
!! I4_PARTITION_NEXT generates the partitions of an I4, one at a time.
!
!  Discussion:
!
!    A partition of an integer N is a representation of the integer
!    as the sum of nonzero positive integers.  The order of the summands
!    does not matter.  Thus, the number 5 has the following partitions
!    and no more:
!
!    5 = 5
!      = 4 + 1
!      = 3 + 2
!      = 3 + 1 + 1
!      = 2 + 2 + 1
!      = 2 + 1 + 1 + 1
!      = 1 + 1 + 1 + 1 + 1
!
!    so the number of partitions of 5 is 7.
!
!    The number of partitions of N is:
!
!      1     1
!      2     2
!      3     3
!      4     5
!      5     7
!      6    11
!      7    15
!      8    22
!      9    30
!     10    42
!     11    56
!     12    77
!     13   101
!     14   135
!     15   176
!     16   231
!     17   297
!     18   385
!     19   490
!     20   627
!     21   792
!     22  1002
!     23  1255
!     24  1575
!     25  1958
!     26  2436
!     27  3010
!     28  3718
!     29  4565
!     30  5604
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
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
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!
!    Input/output, integer ( kind = 4 ) NPART, the number of "parts" in
!    the partition.
!
!    Input/output, integer ( kind = 4 ) A(N).  A contains the parts of
!    the partition.  The value of N is represented by
!      N = sum ( 1 <= I <= NPART ) MULT(I) * A(I).
!
!    Input/output, integer ( kind = 4 ) MULT(N).  MULT counts the multiplicity
!    of the parts of the partition.  MULT(I) is the multiplicity
!    of the part A(I), for I = 1 to NPART.
!
!    Input/output, logical DONE.
!    On first call, the user should set DONE to TRUE to signal
!    that the program should initialize data.
!    On each return, the programs sets DONE to FALSE if it
!    has another partition to return.  If the program returns
!    with DONE TRUE, then there are no more partitions.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  logical done
  integer ( kind = 4 ) is
  integer ( kind = 4 ) iu
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) mult(n)
  integer ( kind = 4 ) npart

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_PARTITION_NEXT - Fatal error!'
    write ( *, '(a)' ) '  N must be positive.'
    write ( *, '(a,i8)' ) '  The input value of N was ', n
    stop
  end if

  if ( done ) then

    a(1) = n
    a(2:n) = 0

    mult(1) = 1
    mult(2:n) = 0

    npart = 1
    done = .false.

  else

    if ( 1 < a(npart) .or. 1 < npart ) then

      done = .false.

      if ( a(npart) == 1 ) then
        is = a(npart-1) + mult(npart)
        k = npart - 1
      else
        is = a(npart)
        k = npart
      end if

      iw = a(k) - 1
      iu = is / iw
      iv = mod ( is, iw )
      mult(k) = mult(k) - 1

      if ( mult(k) == 0 ) then
        k1 = k
      else
        k1 = k + 1
      end if

      mult(k1) = iu
      a(k1) = iw

      if ( iv == 0 ) then
        npart = k1
      else
        mult(k1+1) = 1
        a(k1+1) = iv
        npart = k1 + 1
      end if

    else
      done = .true.
    end if

  end if

  return
end
subroutine i4_partition_next2 ( n, npart, a, mult, more )

!*****************************************************************************80
!
!! I4_PARTITION_NEXT2 computes the partitions of the integer N one at a time.
!
!  Discussion:
!
!    A partition of an integer N is a representation of the integer
!    as the sum of nonzero positive integers.  The order of the summands
!    does not matter.  Thus, the number 5 has the following partitions
!    and no more:
!
!    5 = 5
!      = 4 + 1
!      = 3 + 2
!      = 3 + 1 + 1
!      = 2 + 2 + 1
!      = 2 + 1 + 1 + 1
!      = 1 + 1 + 1 + 1 + 1
!
!    so the number of partitions of 5 is 7.
!
!    Unlike compositions, order is not important in a partition.  Thus
!    the sequences 3+2+1 and 1+2+3 represent distinct compositions, but
!    not distinct partitions.  Also 0 is never returned as one of the
!    elements of the partition.
!
!  Example:
!
!    Sample partitions of 6 include:
!
!      6 = 4+1+1 = 3+2+1 = 2+2+2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer whose partitions are desired.
!
!    Input/output, integer ( kind = 4 ) NPART, the number of distinct, nonzero
!    parts in the output partition.
!
!    Input/output, integer ( kind = 4 ) A(N).  A(I) is the I-th distinct part
!    of the partition, for I = 1, NPART.  Note that if a certain number
!    shows up several times in the partition, it is listed only
!    once in A, and its multiplicity is counted in MULT.
!
!    Input/output, integer ( kind = 4 ) MULT(N).  MULT(I) is the multiplicity
!    of A(I) in the partition, for I = 1, NPART; that is, the number of repeated
!    times that A(I) is used in the partition.
!
!    Input/output, logical MORE.  Set MORE = FALSE on first call.  It
!    will be reset TRUE on return with the first partition.
!    Keep calling for more partitions until MORE
!    is returned FALSE.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) iff
  integer ( kind = 4 ) is
  integer ( kind = 4 ) isum
  logical more
  integer ( kind = 4 ) mult(n)
  integer ( kind = 4 ), save :: nlast = 0
  integer ( kind = 4 ) npart
!
!  On the first call, set NLAST to 0.
!
  if ( .not. more ) then
    nlast = 0
  end if

  if ( n /= nlast .or. ( .not. more ) ) then
    nlast = n
    npart = 1
    a(npart) = n
    mult(npart) = 1
    more = mult(npart) /= n
    return
  end if

  isum = 1

  if ( a(npart) <= 1 ) then
    isum = mult(npart) + 1
    npart = npart - 1
  end if

  iff = a(npart) - 1

  if ( mult(npart) /= 1 ) then
    mult(npart) = mult(npart) - 1
    npart = npart + 1
  end if

  a(npart) = iff
  mult(npart) = 1 + ( isum / iff )
  is = mod ( isum, iff )

  if ( 0 < is ) then
    npart = npart + 1
    a(npart) = is
    mult(npart) = 1
  end if
!
!  There are more partitions, as long as we haven't just computed
!  the last one, which is N copies of 1.
!
  more = mult(npart) /= n

  return
end
subroutine i4_partition_print ( n, npart, a, mult )

!*****************************************************************************80
!
!! I4_PARTITION_PRINT prints a partition of an I4.
!
!  Discussion:
!
!    A partition of an integer N is a representation of the integer as
!    the sum of nonzero integers:
!
!      N = A1 + A2 + A3 + ...
!
!    It is standard practice to gather together all the values that
!    are equal, and replace them in the sum by a single term, multiplied
!    by its "multiplicity":
!
!      N = M1 * A1 + M2 * A2 + ... + M(NPART) * A(NPART)
!
!    In this representation, every A is a unique positive number, and
!    no M is zero (or negative).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!
!    Input, integer ( kind = 4 ) NPART, the number of "parts" in the partition.
!
!    Input, integer ( kind = 4 ) A(NPART), the parts of the partition.
!
!    Input, integer ( kind = 4 ) MULT(NPART), the multiplicies of the parts.
!
  implicit none

  integer ( kind = 4 ) npart

  integer   ( kind = 4 ) a(npart)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) mult(npart)
  integer   ( kind = 4 ) n
  character ( len = 80 ) string

  string(1:2) = '  '
  call i4_to_s_left ( n, string(3:) )
  j = len_trim ( string )
  string (j+1:j+3) = ' = '
  j = j + 3

  do i = 1, npart

    if ( 1 < i ) then
      string(j+1:j+3) = ' + '
      j = j + 3
    end if

    call i4_to_s_left ( mult(i), string(j+1:) )
    j = len_trim ( string )

    string(j+1:j+3) = ' * '
    j = j + 3

    call i4_to_s_left ( a(i), string(j+1:) )
    j = len_trim ( string )

  end do

  write ( *, '(a)' ) trim ( string )

  return
end
subroutine i4_partition_random ( n, table, seed, a, mult, npart )

!*****************************************************************************80
!
!! I4_PARTITION_RANDOM selects a random partition of the integer N.
!
!  Discussion:
!
!    A partition of an integer N is a representation of the integer
!    as the sum of nonzero positive integers.  The order of the summands
!    does not matter.  Thus, the number 5 has the following partitions
!    and no more:
!
!    5 = 5
!      = 4 + 1
!      = 3 + 2
!      = 3 + 1 + 1
!      = 2 + 2 + 1
!      = 2 + 1 + 1 + 1
!      = 1 + 1 + 1 + 1 + 1
!
!    so the number of partitions of 5 is 7.
!
!    Note that some elements of the partition may be 0.  The partition is
!    returned as (MULT(I),I), with NPART nonzero entries in MULT, and
!
!      N = sum ( 1 <= I <= N ) MULT(I) * I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2004
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
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!
!    Input, integer ( kind = 4 ) TABLE(N), the number of partitions of each
!    integer from 1 to N.  This table may be computed by I4_PARTITION_COUNT2.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, integer ( kind = 4 ) A(N), contains in A(1:NPART) the parts of the
!    partition.
!
!    Output, integer ( kind = 4 ) MULT(N), contains in MULT(1:NPART) the
!    multiplicity of the parts.
!
!    Output, integer ( kind = 4 ) NPART, the number of parts in the partition
!    chosen, that is, the number of integers I with nonzero multiplicity
!    MULT(I).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) id
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mult(n)
  integer ( kind = 4 ) npart
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) table(n)
  real ( kind = 8 ) z

  m = n
  npart = 0
  mult(1:n) = 0

  do while ( 0 < m )

    z = r8_uniform_01 ( seed )
    z = m * table(m) * z
    id = 1
    i1 = m
    j = 0

    do

      j = j + 1
      i1 = i1 - id

      if ( i1 < 0 ) then
        id = id + 1
        i1 = m
        j = 0
        cycle
      end if

      if ( i1 == 0 ) then
        z = z - id
        if ( 0.0D+00 < z ) then
          id = id + 1
          i1 = m
          j = 0
          cycle
        else
          exit
        end if
      end if

      if ( 0 < i1 ) then
        z = z - id * table(i1)
        if ( z <= 0.0D+00 ) then
          exit
        end if
      end if

    end do

    mult(id) = mult(id) + j
    npart = npart + j
    m = i1

  end do
!
!  Reformulate the partition in the standard form.
!  NPART is the number of distinct parts.
!
  npart = 0

  do i = 1, n
    if ( mult(i) /= 0 ) then
      npart = npart + 1
      a(npart) = i
      mult(npart) = mult(i)
    end if
  end do

  mult(npart+1:n) = 0

  return
end
subroutine i4_partitions_next ( s, m )

!*****************************************************************************80
!
!! I4_PARTITIONS_NEXT: next partition into S parts.
!
!  Discussion:
!
!    This function generates, one at a time, entries from the list of
!    nondecreasing partitions of the integers into S or fewer parts.
!
!    The list is ordered first by the integer that is partitioned
!    (the sum of the entries), and second by decreasing lexical order
!    in the partition vectors.
!
!    The first value returned is the only such partition of 0.
!
!    Next comes the only partition of 1.
!
!    There follow two partitions of 2, and so on.
!
!    Typical use of this function begins with an initialization call,
!    and then repeated calls in which the output from the previous call
!    is used as input to the next call:
!
!    m = [ 0, 0, 0 ];
!
!    while ( condition )
!      m = i4_partitions_next ( s, m );
!    end
!
!  Example:
!
!    S = 3
!
!    P  D    M
!    _  _  _____
!    1  0  0 0 0
!    2  1  1 0 0
!    3  2  2 0 0
!    4  2  1 1 0
!    5  3  3 0 0
!    6  3  2 1 0
!    7  3  1 1 1
!    8  4  4 0 0
!    9  4  3 1 0
!   10  4  2 2 0
!   11  4  2 1 1
!   12  5  5 0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2010
!
!  Author:
!
!    Original MATLAB version by Alan Genz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) S, the number of items in the partition.
!
!    Input/output, integer ( kind = 4 ) M(S).  On input, the current partition.
!    On first call, this should be a nondecreasing partition.  Thereafter, it
!    should be the output partition from the previous call.  On output, the
!    next partition.
!
  implicit none

  integer ( kind = 4 ) s

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m(s)
  integer ( kind = 4 ) msum

  msum = m(1)

  do i = 2, s

    msum = msum + m(i)

    if ( m(1) <= m(i) + 1 ) then
      m(i) = 0
    else
      m(1) = msum - ( i - 1 ) * ( m(i) + 1 )
      m(2:i) = m(i) + 1
      return
    end if

  end do
!
!  If we failed to find a suitable index I, put
!  the entire sum into M(1), increment by 1, and
!  prepare to partition the next integer.
!
  m(1) = msum + 1

  return
end
function i4_rise ( x, n )

!*****************************************************************************80
!
!! I4_RISE computes the rising factorial function [X]^N.
!
!  Discussion:
!
!    [X]^N = X * ( X + 1 ) * ( X + 2 ) * ... * ( X + N - 1 ).
!
!    Note that the number of ways of arranging N objects in M ordered
!    boxes is [M]^N.  (Here, the ordering of the objects in each box matters).
!    Thus, 2 objects in 2 boxes have the following 6 possible arrangements:
!
!      -|12, 1|2, 12|-, -|21, 2|1, 21|-.
!
!    Moreover, the number of non-decreasing maps from a set of
!    N to a set of M ordered elements is [M]^N / N!.  Thus the set of
!    nondecreasing maps from (1,2,3) to (a,b,c,d) is the 20 elements:
!
!      aaa, abb, acc, add, aab, abc, acd, aac, abd, aad
!      bbb, bcc, bdd, bbc, bcd, bbd, ccc, cdd, ccd, ddd.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the argument of the rising factorial
!    function.
!
!    Input, integer ( kind = 4 ) N, the order of the rising factorial function.
!    If N = 0, RISE = 1, if N = 1, RISE = X.  Note that if N is
!    negative, a "falling" factorial will be computed.
!
!    Output, integer ( kind = 4 ) I4_RISE, the value of the rising factorial
!    function.
!
  implicit none

  integer ( kind = 4 ) arg
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_rise
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value
  integer ( kind = 4 ) x

  value = 1

  arg = x

  if ( 0 < n ) then

    do i = 1, n
      value = value * arg
      arg = arg + 1
    end do

  else if ( n < 0 ) then

    do i = -1, n, -1
      value = value * arg
      arg = arg - 1
    end do

  end if

  i4_rise = value

  return
end
subroutine i4_sqrt ( n, q, r )

!*****************************************************************************80
!
!! I4_SQRT finds the integer square root of N by solving N = Q*Q + R.
!
!  Discussion:
!
!    The integer square root of N is an integer Q such that
!    Q*Q <= N but N < (Q+1)*(Q+1).
!
!    A simpler calculation would be something like
!
!      Q = int ( sqrt ( real ( N ) ) )
!
!    but this calculation has the virtue of using only integer arithmetic.
!
!    To avoid the tedium of worrying about negative arguments, the routine
!    automatically considers the absolute value of the argument.
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
!   John Burkardt
!
!  Reference:
!
!    Mark Herkommer,
!    Number Theory, A Programmer's Guide,
!    McGraw Hill, 1999,
!    ISBN: 0-07-913074-7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number whose integer square root is
!    desired.  Actually, only the absolute value of N is considered.
!
!    Output, integer ( kind = 4 ) Q, R, the integer square root, and positive
!    remainder, of N.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_abs
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r

  n_abs = abs ( n )

  q = n_abs

  if ( 0 < n_abs ) then

    do while ( ( n_abs / q ) < q )
      q = ( q + ( n_abs / q ) ) / 2
    end do

  end if

  r = n_abs - q * q

  return
end
subroutine i4_sqrt_cf ( n, max_term, n_term, b )

!*****************************************************************************80
!
!! I4_SQRT_CF: continued fraction representation of a square root of an integer.
!
!  Discussion:
!
!    The continued fraction representation of the square root of an integer
!    has the form
!
!      [ B0, (B1, B2, B3, ..., BM), ... ]
!
!    where
!
!      B0 = int ( sqrt ( real ( N ) ) )
!      BM = 2 * B0
!      the sequence ( B1, B2, B3, ..., BM ) repeats in the representation.
!      the value M is termed the period of the representation.
!
!  Example:
!
!     N  Period  Continued Fraction
!
!     2       1  [ 1, 2, 2, 2, ... ]
!     3       2  [ 1, 1, 2, 1, 2, 1, 2... ]
!     4       0  [ 2 ]
!     5       1  [ 2, 4, 4, 4, ... ]
!     6       2  [ 2, 2, 4, 2, 4, 2, 4, ... ]
!     7       4  [ 2, 1, 1, 1, 4, 1, 1, 4, 1, 1, 4... ]
!     8       2  [ 2, 1, 4, 1, 4, 1, 4, 1, 4, ... ]
!     9       0  [ 3 ]
!    10       1  [ 3, 6, 6, 6, ... ]
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
!   John Burkardt
!
!  Reference:
!
!    Mark Herkommer,
!    Number Theory, A Programmer's Guide,
!    McGraw Hill, 1999,
!    ISBN: 0-07-913074-7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number whose continued fraction square
!    root is desired.
!
!    Input, integer ( kind = 4 ) MAX_TERM, the maximum number of terms that may
!    be computed.
!
!    Output, integer ( kind = 4 ) N_TERM, the number of terms computed beyond
!    the 0 term.  The routine should stop if it detects that the period
!    has been reached.
!
!    Output, integer ( kind = 4 ) B(0:MAX_TERM), contains the continued fraction
!    coefficients in entries B(0:N_TERM).
!
  implicit none

  integer ( kind = 4 ) max_term

  integer ( kind = 4 ) b(0:max_term)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_term
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s

  n_term = 0

  call i4_sqrt ( n, s, r )
  b(0) = s

  if ( 0 < r ) then

    p = 0
    q = 1

    do

      p = b(n_term) * q - p
      q = ( n - p * p ) / q

      if ( max_term <= n_term ) then
        return
      end if

      n_term = n_term + 1
      b(n_term) = ( p + s ) / q

      if ( q == 1 ) then
        exit
      end if

    end do

  end if

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
subroutine i4_to_bvec ( i4, n, bvec )

!*****************************************************************************80
!
!! I4_TO_BVEC makes a signed binary vector from an I4.
!
!  Discussion:
!
!    A BVEC is a vector of binary digits representing an integer.
!
!    BVEC(1) is 0 for positive values and 1 for negative values, which
!    are stored in 2's complement form.
!
!    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
!    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!    Negative values have a two's complement operation applied.
!
!    To guarantee that there will be enough space for any
!    value of I, it would be necessary to set N = 32.
!
!  Example:
!
!    I4       BVEC         binary
!    --  ----------------  ------
!     1  1  0  0  0  0  1      1
!     2  0  0  0  0  1  0     10
!     3  0  0  0  0  1  1     11
!     4  0  0  0  1  0  0    100
!     9  0  0  1  0  0  1   1001
!    -9  1  1  0  1  1  1  -1001 = 110111 (2's complement)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I4, an integer to be represented.
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Output, integer ( kind = 4 ) BVEC(N), the signed binary representation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 2
  integer ( kind = 4 ) bvec(n)
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_copy
  integer ( kind = 4 ) j

  i4_copy = abs ( i4 )

  do j = n, 2, - 1

    bvec(j) = mod ( i4_copy, base )

    i4_copy = i4_copy / base

  end do

  bvec(1) = 0

  if ( i4 < 0 ) then
    call bvec_complement2 ( n, bvec, bvec )
  end if

  return
end
subroutine i4_to_chinese ( j, n, m, r )

!*****************************************************************************80
!
!! I4_TO_CHINESE converts an I4 to its Chinese remainder form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) J, the integer to be converted.
!
!    Input, integer ( kind = 4 ) N, the number of moduluses.
!
!    Input, integer ( kind = 4 ) M(N), the moduluses.  These should be positive
!    and pairwise prime.
!
!    Output, integer ( kind = 4 ) R(N), the Chinese remainder representation
!    of the integer.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m(n)
  integer ( kind = 4 ) r(n)

  call chinese_check ( n, m, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_CHINESE - Fatal error!'
    write ( *, '(a)' ) '  The moduluses are not legal.'
    stop
  end if

  do i = 1, n
    r(i) = i4_modp ( j, m(i) )
  end do

  return
end
subroutine i4_to_dvec ( i4, n, dvec )

!*****************************************************************************80
!
!! I4_TO_DVEC makes a signed DVEC from an I4.
!
!  Discussion:
!
!    A DVEC is an integer vector of decimal digits, intended to
!    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
!    is the coefficient of 10**(N-2), and DVEC(N) contains sign
!    information.  It is 0 if the number is positive, and 9 if
!    the number is negative.
!
!    Negative values have a ten's complement operation applied.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I4, an integer to be represented.
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Output, integer ( kind = 4 ) DVEC(N), the signed decimal representation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 10
  integer ( kind = 4 ) dvec(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_copy

  i4_copy = abs ( i4 )

  do i = 1, n-1

    dvec(i) = mod ( i4_copy, base )

    i4_copy = i4_copy / base

  end do

  dvec(n) = 0

  if ( i4 < 0 ) then
    call dvec_complementx ( n, dvec, dvec )
  end if

  return
end
subroutine i4_to_i4poly ( intval, base, degree_max, degree, a )

!*****************************************************************************80
!
!! I4_TO_I4POLY converts an I4 to an I4POLY in a given base.
!
!  Example:
!
!    INTVAL  BASE  Degree     A (in reverse order!)
!
!         1     2       0     1
!         6     2       2     1  1  0
!        23     2       5     1  0  1  1  1
!        23     3       3     2  1  2
!        23     4       3     1  1  3
!        23     5       2     4  3
!        23     6       2     3  5
!        23    23       1     1  0
!        23    24       0    23
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.
!
!    Input, integer ( kind = 4 ) BASE, the base, which should be greater than 1.
!
!    Input, integer ( kind = 4 ) DEGREE_MAX, the maximum degree.
!
!    Output, integer ( kind = 4 ) DEGREE, the degree of the polynomial.
!
!    Output, integer ( kind = 4 ) A(0:DEGREE_MAX), contains the coefficients
!    of the polynomial expansion of INTVAL in base BASE.
!
  implicit none

  integer ( kind = 4 ) degree_max

  integer ( kind = 4 ) a(0:degree_max)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) degree
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) j

  a(0:degree_max) = 0

  j = abs ( intval )

  degree = 0

  a(degree) = mod ( j, base )

  j = j - a(degree)
  j = j / base

  do while ( 0 < j )

    degree = degree + 1

    if ( degree < degree_max ) then
      a(degree) = mod ( j, base )
    end if

    j = j - a(degree)
    j = j / base

  end do

  if ( intval < 0 ) then
    a(0:degree_max) = -a(0:degree_max)
  end if

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

    ival = -ival
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
subroutine i4_to_van_der_corput ( seed, base, r )

!*****************************************************************************80
!
!! I4_TO_VAN_DER_CORPUT computes an element of a van der Corput sequence.
!
!  Discussion:
!
!    The van der Corput sequence is often used to generate a "subrandom"
!    sequence of points which have a better covering property
!    than pseudorandom points.
!
!    The van der Corput sequence generates a sequence of points in [0,1]
!    which (theoretically) never repeats.  Except for SEED = 0, the
!    elements of the van der Corput sequence are strictly between 0 and 1.
!
!    The van der Corput sequence writes an integer in a given base B,
!    and then its digits are "reflected" about the decimal point.
!    This maps the numbers from 1 to N into a set of numbers in [0,1],
!    which are especially nicely distributed if N is one less
!    than a power of the base.
!
!    Hammersley suggested generating a set of N nicely distributed
!    points in two dimensions by setting the first component of the
!    Ith point to I/N, and the second to the van der Corput
!    value of I in base 2.
!
!    Halton suggested that in many cases, you might not know the number
!    of points you were generating, so Hammersley's formulation was
!    not ideal.  Instead, he suggested that to generate a nicely
!    distributed sequence of points in M dimensions, you simply
!    choose the first M primes, P(1:M), and then for the J-th component of
!    the I-th point in the sequence, you compute the van der Corput
!    value of I in base P(J).
!
!    Thus, to generate a Halton sequence in a 2 dimensional space,
!    it is typical practice to generate a pair of van der Corput sequences,
!    the first with prime base 2, the second with prime base 3.
!    Similarly, by using the first K primes, a suitable sequence
!    in K-dimensional space can be generated.
!
!    The generation is quite simple.  Given an integer SEED, the expansion
!    of SEED in base BASE is generated.  Then, essentially, the result R
!    is generated by writing a decimal point followed by the digits of
!    the expansion of SEED, in reverse order.  This decimal value is actually
!    still in base BASE, so it must be properly interpreted to generate
!    a usable value.
!
!  Example:
!
!    BASE = 2
!
!    SEED     SEED      van der Corput
!    decimal  binary    binary   decimal
!    -------  ------    ------   -------
!        0  =     0  =>  .0     = 0.0
!        1  =     1  =>  .1     = 0.5
!        2  =    10  =>  .01    = 0.25
!        3  =    11  =>  .11    = 0.75
!        4  =   100  =>  .001   = 0.125
!        5  =   101  =>  .101   = 0.625
!        6  =   110  =>  .011   = 0.375
!        7  =   111  =>  .111   = 0.875
!        8  =  1000  =>  .0001  = 0.0625
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Halton,
!    On the efficiency of certain quasi-random sequences of points
!    in evaluating multi-dimensional integrals,
!    Numerische Mathematik,
!    Volume 2, pages 84-90, 1960.
!
!    John Hammersley,
!    Monte Carlo methods for solving multivariable problems,
!    Proceedings of the New York Academy of Science,
!    Volume 86, pages 844-874, 1960.
!
!    Johannes van der Corput,
!    Verteilungsfunktionen I & II,
!    Proceedings of the Koninklijke Nederlandsche Akademie
!    van Wetenschappen,
!    Volume 38, 1935, pages 813-820, pages 1058-1066.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED, the seed or index of the desired element.
!    SEED should be nonnegative.
!    SEED = 0 is allowed, and returns R = 0.
!
!    Input, integer ( kind = 4 ) BASE, the van der Corput base, which is
!    typically a prime number.  BASE must be greater than 1.
!
!    Output, real ( kind = 8 ) R, the SEED-th element of the van der
!    Corput sequence for base BASE.
!
  implicit none

  integer ( kind = 4 ) base
  real ( kind = 8 ) base_inv
  integer ( kind = 4 ) digit
  real ( kind = 8 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed2

  if ( base <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_VAN_DER_CORPUT - Fatal error!'
    write ( *, '(a)' ) '  The input base BASE is <= 1!'
    write ( *, '(a,i8)' ) '  BASE = ', base
    stop
  end if

  if ( seed < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_VAN_DER_CORPUT - Fatal error!'
    write ( *, '(a)' ) '  The input base SEED is < 0!'
    write ( *, '(a,i8)' ) '  SEED = ', seed
    stop
  end if

  seed2 = seed

  r = 0.0D+00

  base_inv = 1.0D+00 / real ( base, kind = 8 )

  do while ( seed2 /= 0 )
    digit = mod ( seed2, base )
    r = r + real ( digit, kind = 8 ) * base_inv
    base_inv = base_inv / real ( base, kind = 8 )
    seed2 = seed2 / base
  end do

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
!    12 November 2006
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
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
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
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
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
    seed = seed + huge ( seed )
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
subroutine i4mat_01_rowcolsum ( m, n, r, c, a, ierror )

!*****************************************************************************80
!
!! I4MAT_01_ROWCOLSUM creates a 0/1 I4MAT with given row and column sums.
!
!  Discussion:
!
!    Given an M vector R and N vector C, there may exist one or more
!    M by N matrices with entries that are 0 or 1, whose row sums are R
!    and column sums are C.
!
!    For convenience, this routine requires that the entries of R and C
!    be given in nonincreasing order.
!
!    There are several requirements on R and C.  The simple requirements
!    are that the entries of R and C must be nonnegative, that the entries
!    of R must each be no greater than N, and those of C no greater than M,
!    and that the sum of the entries of R must equal the sum of the entries
!    of C.
!
!    The final technical requirement is that if we form R*, the conjugate
!    partition of R, then C is majorized by R*, that is, that every partial
!    sum from 1 to K of the entries of C is no bigger than the sum of the same
!    entries of R*, for every K from 1 to N.
!
!    Given these conditions on R and C, there is at least one 0/1 matrix
!    with the given row and column sums.
!
!    The conjugate partition of R is constructed as follows:
!      R*(1) is the number of entries of R that are 1 or greater.
!      R*(2) is the number of entries of R that are 2 or greater.
!      ...
!      R*(N) is the number of entries of R that are N (can't be greater).
!
!  Example:
!
!    M = N = 5
!    R = ( 3, 2, 2, 1, 1 )
!    C = ( 2, 2, 2, 2, 1 )
!
!    A =
!      1 0 1 0 1
!      1 0 0 1 0
!      0 1 0 1 0
!      0 1 0 0 0
!      0 0 1 0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack vanLint, Richard Wilson,
!    A Course in Combinatorics,
!    Cambridge, 1992,
!    ISBN: 0-521-42260-4.
!
!    James Sandeson,
!    Testing Ecological Patterns,
!    American Scientist,
!    Volume 88, July-August 2000, pages 332-339.
!
!    Ian Saunders,
!    Algorithm AS 205,
!    Enumeration of R x C Tables with Repeated Row Totals,
!    Applied Statistics,
!    Volume 33, Number 3, pages 340-352, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in the
!    array.
!
!    Input, integer ( kind = 4 ) R(M), C(N), the row and column sums desired
!    for the array.  Both vectors must be arranged in descending order.
!    The elements of R must be between 0 and N.
!    The elements of C must be between 0 and M.
!
!    Output, integer ( kind = 4 ) A(M,N), the M by N matrix with the given
!    row and column sums.
!    Each entry of A is 0 or 1.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error was encountered, and the array was computed.
!    1, R and C do not have the same total.
!    2, R is not monotone decreasing, or has illegal entries.
!    3, C is not monotone decreasing, or has illegal entries.
!    4, R and C are not a possible set of row and column sums.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) c_sum
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  logical i4vec_descends
  integer ( kind = 4 ) i4vec_maxloc_last
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) r(m)
  integer ( kind = 4 ) r_conj(n)
  integer ( kind = 4 ) r_sum
  integer ( kind = 4 ) r2(m)

  a(1:m,1:n) = 0
!
!  Check conditions.
!
  ierror = 0

  if ( sum ( r(1:m) ) /= sum ( c(1:n) ) ) then
    ierror = 1
    return
  end if

  if ( .not. i4vec_descends ( m, r ) ) then
    ierror = 2
    return
  end if

  if ( n < r(1) .or. r(m) < 0 ) then
    ierror = 2
    return
  end if

  if ( .not. i4vec_descends ( n, c ) ) then
    ierror = 3
    return
  end if

  if ( m < c(1) .or. c(n) < 0 ) then
    ierror = 3
    return
  end if
!
!  Compute the conjugate of R.
!
  r_conj(1:n) = 0

  do i = 1, m
    do j = 1, r(i)
      r_conj(j) = r_conj(j) + 1
    end do
  end do
!
!  C must be majorized by R_CONJ.
!
  r_sum = 0
  c_sum = 0
  do i = 1, n
    r_sum = r_sum + r_conj(i)
    c_sum = c_sum + c(i)
    if ( r_sum < c_sum ) then
      ierror = 4
      return
    end if
  end do
!
!  We need a temporary copy of R that we can decrement.
!
  r2(1:m) = r(1:m)

  do j = n, 1, -1

    i = i4vec_maxloc_last ( m, r2 )

    do k = 1, c(j)
!
!  By adding 1 rather than setting A(I,J) to 1, we were able to spot
!  an error where the index was "sticking".
!
      a(i,j) = a(i,j) + 1

      r2(i) = r2(i) - 1
!
!  There's a special case you have to watch out for.
!  If I was 1, and when you decrement R2(1), I is going to be 1 again,
!  and you're staying in the same column, that's not good.
!
      if ( 1 < i ) then
        i = i - 1
      else
        i = i4vec_maxloc_last ( m, r2 )
        if ( i == 1 .and. k < c(j) ) then
          i = 1 + i4vec_maxloc_last ( m-1, r2(2:m) )
        end if
      end if

    end do

  end do

  return
end
subroutine i4mat_01_rowcolsum2 ( m, n, r, c, a, ierror )

!*****************************************************************************80
!
!! I4MAT_01_ROWCOLSUM2 creates a 0/1 I4MAT with given row and column sums.
!
!  Discussion:
!
!    This routine uses network flow optimization to compute the results.
!
!    Given an M vector R and N vector C, there may exist one or more
!    M by N matrices with entries that are 0 or 1, whose row sums are R
!    and column sums are C.
!
!    For convenience, this routine requires that the entries of R and C
!    be given in nonincreasing order.
!
!    There are several requirements on R and C.  The simple requirements
!    are that the entries of R and C must be nonnegative, that the entries
!    of R must each no greater than N, and those of C no greater than M,
!    and that the sum of the entries of R must equal the sum of the
!    entries of C.
!
!    The final technical requirement is that if we form R*, the conjugate
!    partition of R, then C is majorized by R*, that is, that every partial
!    sum from 1 to K of the entries of C is no bigger than the sum of the same
!    entries of R*, for every K from 1 to N.
!
!    Given these conditions on R and C, there is at least one 0/1 matrix
!    with the given row and column sums.
!
!    The conjugate partition of R is constructed as follows:
!      R*(1) is the number of entries of R that are 1 or greater.
!      R*(2) is the number of entries of R that are 2 or greater.
!      ...
!      R*(N) is the number of entries of R that are N (can't be greater).
!
!  Example:
!
!    M = N = 5
!    R = ( 3, 2, 2, 1, 1 )
!    C = ( 2, 2, 2, 2, 1 )
!
!    A =
!      1 0 1 0 1
!      1 0 0 1 0
!      0 1 0 1 0
!      0 1 0 0 0
!      0 0 1 0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2003
!
!  Author:
!
!    John Burkardt
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
!    Jack van Lint, Richard Wilson,
!    A Course in Combinatorics,
!    Oxford, 1992, pages 148-156.
!
!    James Sandeson,
!    Testing Ecological Patterns,
!    American Scientist,
!    Volume 88, July-August 2000, pages 332-339.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the array.  These values do not have to be equal.
!
!    Input, integer ( kind = 4 ) R(M), C(N), the row and column sums desired
!    for the array.  Both vectors must be arranged in descending order.
!    The elements of R must be between 0 and N.
!    The elements of C must be between 0 and M.
!    One of the conditions for a solution to exist is that the sum of the
!    elements in R equal the sum of the elements in C.
!
!    Output, integer ( kind = 4 ) A(M,N), the matrix with the given row and
!    column sums.  Each entry of A is 0 or 1.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error was encountered, and the array was computed.
!    1, R and C are not consistent.  A partial solution may be constructed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) capflo(2,2*(m+m*n+n))
  integer ( kind = 4 ) cut(m+n+2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iendpt(2,2*(m+m*n+n))
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) node_flow(m+n+2)
  integer ( kind = 4 ) r(m)
  integer ( kind = 4 ) sink
  integer ( kind = 4 ) source

  ierror = 0
!
!  There are M + N + 2 nodes.  The last two are the special source and sink.
!
  source = m + n + 1
  sink = m + n + 2
  nnode = m + n + 2
!
!  The source is connected to each of the R nodes.
!
  k = 0

  do i = 1, m

    k = k + 1
    iendpt(1,k) = source
    iendpt(2,k) = i
    capflo(1,k) = r(i)
    capflo(2,k) = 0

    k = k + 1
    iendpt(1,k) = i
    iendpt(2,k) = source
    capflo(1,k) = r(i)
    capflo(2,k) = 0

  end do
!
!  Every R node is connected to every C node, with capacity 1.
!
  do i = 1, m
    do j = 1, n

      k = k + 1
      iendpt(1,k) = i
      iendpt(2,k) = j+m
      capflo(1,k) = 1
      capflo(2,k) = 0

      k = k + 1
      iendpt(1,k) = j+m
      iendpt(2,k) = i
      capflo(1,k) = 1
      capflo(2,k) = 0

    end do
  end do
!
!  Every C node is connected to the sink.
!
  do j = 1, n

    k = k + 1
    iendpt(1,k) = j+m
    iendpt(2,k) = sink
    capflo(1,k) = c(j)
    capflo(2,k) = 0

    k = k + 1
    iendpt(1,k) = sink
    iendpt(2,k) = j+m
    capflo(1,k) = c(j)
    capflo(2,k) = 0

  end do
!
!  Determine the maximum flow on the network.
!
  nedge = k

  call network_flow_max ( nnode, nedge, iendpt, capflo, source, sink, &
    cut, node_flow )
!
!  We have a perfect solution if, and only if, the edges leading from the
!  source, and the edges leading to the sink, are all saturated.
!
  do k = 1, nedge

    i = iendpt(1,k)
    j = iendpt(2,k) - m

    if ( i <= m .and. 1 <= j .and. j <= n ) then
      if ( capflo(2,k) /= 0 .and. capflo(2,k) /= 1 ) then
        ierror = 1
      end if
    end if

    if ( iendpt(1,k) == source ) then
      if ( capflo(1,k) /= capflo(2,k) ) then
        ierror = 1
      end if
    end if

    if ( iendpt(2,k) == sink ) then
      if ( capflo(1,k) /= capflo(2,k) ) then
        ierror = 1
      end if
    end if

  end do
!
!  If we have a solution, then A(I,J) = the flow on the edge from
!  R node I to C node J.
!
  a(1:m,1:n) = 0

  do k = 1, nedge

    i = iendpt(1,k)
    j = iendpt(2,k) - m

    if ( i <= m .and. 1 <= j .and. j <= n ) then
      a(i,j) = capflo(2,k)
    end if

  end do

  return
end
subroutine i4mat_perm ( n, a, p )

!*****************************************************************************80
!
!! I4MAT_PERM permutes the rows and columns of a square I4MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2000
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, integer ( kind = 4 ) A(N,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) is the new
!    number of row and column I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) p(n)

  call perm_cycle ( n, 1, p, is, nc )

  do i = 1, n

    i1 = -p(i)

    if ( 0 < i1 ) then

      lc = 0

      do

        i1 = p(i1)
        lc = lc + 1

        if ( i1 <= 0 ) then
          exit
        end if

      end do

      i1 = i

      do j = 1, n

        if ( p(j) <= 0 ) then

          j2 = j
          k = lc

          do

            j1 = j2
            it = a(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( p(j1) )

              call i4_swap ( a(i1,j1), it )

              if ( j1 /= j2 ) then
                cycle
              end if

              k = k - 1

              if ( i1 == i ) then
                exit
              end if

            end do

            j2 = abs ( p(j2) )

            if ( k == 0 ) then
              exit
            end if

          end do

        end if

      end do

    end if

  end do
!
!  Restore the positive signs of the data.
!
  p(1:n) = abs ( p(1:n) )

  return
end
subroutine i4mat_perm2 ( m, n, a, p, q )

!*****************************************************************************80
!
!! I4MAT_PERM2 permutes the rows and columns of a rectangular I4MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 October 1999
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
!    Input, integer ( kind = 4 ) M, the number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) N, the number of columns in the matrix.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer ( kind = 4 ) P(M), the row permutation.  P(I) is the new
!    number of row I.
!
!    Input, integer ( kind = 4 ) Q(N), the column permutation.  Q(I) is the
!    new number of column I.  Note that this routine allows you to pass a
!    single array as both P and Q.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) q(n)

  call perm_cycle ( m, 1, p, is, nc )

  if ( 0 < q(1) ) then
    call perm_cycle ( n, 1, q, is, nc )
  end if

  do i = 1, m

    i1 = -p(i)

    if ( 0 < i1 ) then

      lc = 0

      do

        i1 = p(i1)
        lc = lc + 1

        if ( i1 <= 0 ) then
          exit
        end if

      end do

      i1 = i

      do j = 1, n

        if ( q(j) <= 0 ) then

          j2 = j
          k = lc

          do

            j1 = j2
            it = a(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( q(j1) )

              call i4_swap ( a(i1,j1), it )

              if ( j1 /= j2 ) then
                cycle
              end if

              k = k - 1

              if ( i1 == i ) then
                exit
              end if

            end do

            j2 = abs ( q(j2) )

            if ( k == 0 ) then
              exit
            end if

          end do

        end if

      end do

    end if

  end do
!
!  Restore the positive signs of the data.
!
  p(1:m) = abs ( p(1:m) )

  if ( q(1) <= 0 ) then
    q(1:n) = abs ( q(1:n) )
  end if

  return
end
subroutine i4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_PRINT prints an I4MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2004
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
!    Input, integer ( kind = 4 ) A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = * ) title

  call i4mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_PRINT_SOME prints some of an I4MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 10
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = 7 ) ctemp(incx)
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
      write ( ctemp(j2), '(i7)') j
    end do

    write ( *, '(''  Col   '',10a7)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(i7)' ) a(i,j)

      end do

      write ( *, '(i8,2x,10a7)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine i4mat_u1_inverse ( n, a, b )

!*****************************************************************************80
!
!! I4MAT_U1_INVERSE inverts a unit upper triangular I4MAT.
!
!  Discussion:
!
!    A unit upper triangular matrix is a matrix with only 1's on the main
!    diagonal, and only 0's below the main diagonal.  Above the main
!    diagonal, the entries may be assigned any value.
!
!    It may be surprising to note that the inverse of an integer unit upper
!    triangular matrix is also an integer unit upper triangular matrix.
!
!    Note that this routine can invert a matrix in place, that is, with no
!    extra storage.  If the matrix is stored in A, then the call
!
!      call i4mat_u1_inverse ( n, a, a )
!
!    will result in A being overwritten by its inverse, which can
!    save storage if the original value of A is not needed later.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 1999
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
!    Input, integer ( kind = 4 ) N, the number of rows and columns in
!    the matrix.
!
!    Input, integer ( kind = 4 ) A(N,N), the unit upper triangular matrix
!    to be inverted.
!
!    Output, integer ( kind = 4 ) B(N,N), the inverse matrix.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) b(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isum
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  do j = n, 1, -1

    do i = n, 1, -1

      if ( i == j ) then
        isum = 1
      else
        isum = 0
      end if

      do k = i + 1, j
        isum = isum - a(i,k) * b(k,j)
      end do

      b(i,j) = isum

    end do
  end do

  return
end
subroutine i4poly ( n, a, x0, iopt, value )

!*****************************************************************************80
!
!! I4POLY performs operations on I4POLY's in power or factorial form.
!
!  Discussion:
!
!    The power sum form of a polynomial is
!
!      P(X) = A1 + A2*X + A3*X^2 + ... + (AN+1)*X^N
!
!    The Taylor expansion at C has the form
!
!      P(X) = A1 + A2*(X-C) + A3*(X-C)^2 + ... + (AN+1)*(X-C)^N
!
!    The factorial form of a polynomial is
!
!      P(X) = A1 + A2*X + A3*(X)*(X-1) + A4*(X)*(X-1)*(X-2)+...
!        + (AN+1)*(X)*(X-1)*...*(X-N+1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
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
!    Input, integer ( kind = 4 ) N, the number of coefficients in the
!    polynomial (in other words, the polynomial degree + 1)
!
!    Input/output, integer ( kind = 4 ) A(N), the coefficients of the
!    polynomial.  Depending on the option chosen, these coefficients may be
!    overwritten by those of a different form of the polynomial.
!
!    Input, integer ( kind = 4 ) X0, for IOPT = -1, 0, or positive, the value
!    of the argument at which the polynomial is to be evaluated, or the
!    Taylor expansion is to be carried out.
!
!    Input, integer ( kind = 4 ) IOPT, a flag describing which algorithm is to
!    be carried out:
!    -3: Reverse Stirling.  Input the coefficients of the polynomial in
!    factorial form, output them in power sum form.
!    -2: Stirling.  Input the coefficients in power sum form, output them
!    in factorial form.
!    -1: Evaluate a polynomial which has been input in factorial form.
!    0:  Evaluate a polynomial input in power sum form.
!    1 or more:  Given the coefficients of a polynomial in
!    power sum form, compute the first IOPT coefficients of
!    the polynomial in Taylor expansion form.
!
!    Output, integer ( kind = 4 ) VALUE, for IOPT = -1 or 0, the value of the
!    polynomial at the point X0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iopt
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) value
  integer ( kind = 4 ) w
  integer ( kind = 4 ) x0
  integer ( kind = 4 ) z

  n1 = min ( n, iopt )
  n1 = max ( 1, n1 )

  if ( iopt < -1 ) then
    n1 = n
  end if

  eps = mod ( max ( -iopt, 0 ), 2 )

  w = -n * eps

  if ( -2 < iopt ) then
    w = w + x0
  end if

  do m = 1, n1

    value = 0
    z = w

    do i = m, n
      z = z + eps
      value = a(n+m-i) + z * value
      if ( iopt /= 0 .and. iopt /= -1 ) then
        a(n+m-i) = value
      end if
    end do

    if ( iopt < 0 ) then
      w = w + 1
    end if

  end do

  return
end
subroutine i4poly_cyclo ( n, phi )

!*****************************************************************************80
!
!! I4POLY_CYCLO computes a cyclotomic polynomial.
!
!  Discussion:
!
!    For 1 <= N, let
!
!      I = SQRT ( - 1 )
!      L = EXP ( 2 * PI * I / N )
!
!    Then the N-th cyclotomic polynomial is defined by
!
!      PHI(N;X) = Product ( 1 <= K <= N and GCD(K,N) = 1 ) ( X - L^K )
!
!    We can use the Moebius MU function to write
!
!      PHI(N;X) = Product ( mod ( D, N ) = 0 ) ( X**D - 1 )^MU(N/D)
!
!    There is a sort of inversion formula:
!
!      X^N - 1 = Product ( mod ( D, N ) = 0 ) PHI(D;X)
!
!  Example:
!
!     N  PHI
!
!     0  1
!     1  X - 1
!     2  X + 1
!     3  X^2 + X + 1
!     4  X^2 + 1
!     5  X^4 + X^3 + X^2 + X + 1
!     6  X^2 - X + 1
!     7  X^6 + X^5 + X^4 + X^3 + X^2 + X + 1
!     8  X^4 + 1
!     9  X^6 + X^3 + 1
!    10  X^4 - X^3 + X^2 - X + 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Raymond Seroul,
!    Programming for Mathematicians,
!    Springer, 2000,
!    ISBN: 3-540-66422-X.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the cyclotomic polynomial
!    desired.
!
!    Output, integer ( kind = 4 ) PHI(0:N), the N-th cyclotomic polynomial.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: max_poly = 100

  integer ( kind = 4 ) d
  integer ( kind = 4 ) den(0:max_poly)
  integer ( kind = 4 ) den_n
  integer ( kind = 4 ) factor(0:n)
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) num(0:max_poly)
  integer ( kind = 4 ) num_n
  integer ( kind = 4 ) phi(0:n)
  integer ( kind = 4 ) rem(0:n)

  num(0) = 1
  num(1:max_poly) = 0
  num_n = 0

  den(0) = 1
  den(1:max_poly) = 0
  den_n = 0

  phi(0:n) = 0

  do d = 1, n
!
!  For each divisor D of N, ...
!
    if ( mod ( n, d ) == 0 ) then

      call i4_moebius ( n / d, mu )
!
!  ...multiply the numerator or denominator by (X^D-1).
!
      factor(0) = -1
      factor(1:d-1) = 0
      factor(d) = 1

      if ( mu == + 1 ) then

        if ( max_poly < num_n + d ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I4POLY_CYCLO - Fatal error!'
          write ( *, '(a)' ) '  Numerator polynomial degree too high.'
          stop
        end if

        call i4poly_mul ( num_n, num, d, factor, num )

        num_n = num_n + d

      else if ( mu == -1 ) then

        if ( max_poly < den_n + d ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I4POLY_CYCLO - Fatal error!'
          write ( *, '(a)' ) '  Denominator polynomial degree too high.'
          stop
        end if

        call i4poly_mul ( den_n, den, d, factor, den )

        den_n = den_n + d

      end if

    end if

  end do
!
!  PHI = NUM / DEN
!
  call i4poly_div ( num_n, num, den_n, den, nq, phi, nr, rem )

  return
end
subroutine i4poly_degree ( na, a, degree )

!*****************************************************************************80
!
!! I4POLY_DEGREE returns the degree of an I4POLY.
!
!  Discussion:
!
!    The degree of a polynomial is the index of the highest power
!    of X with a nonzero coefficient.
!
!    The degree of a constant polynomial is 0.  The degree of the
!    zero polynomial is debatable, but this routine returns the
!    degree as 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, the dimension of A.
!
!    Input, integer ( kind = 4 ) A(0:NA), the coefficients of the polynomials.
!
!    Output, integer ( kind = 4 ) DEGREE, the degree of A.
!
  implicit none

  integer ( kind = 4 ) na

  integer ( kind = 4 ) a(0:na)
  integer ( kind = 4 ) degree

  degree = na

  do while ( 0 < degree )

    if ( a(degree) /= 0 ) then
      return
    end if

    degree = degree - 1

  end do

  return
end
subroutine i4poly_div ( na, a, nb, b, nq, q, nr, r )

!*****************************************************************************80
!
!! I4POLY_DIV computes the quotient and remainder of two I4POLY's.
!
!  Discussion:
!
!    Normally, the quotient and remainder would have rational coefficients.
!    This routine assumes that the special case applies that the quotient
!    and remainder are known beforehand to be integral.
!
!    The polynomials are assumed to be stored in power sum form.
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, the degree of polynomial A.
!
!    Input, integer ( kind = 4 ) A(0:NA), the coefficients of the polynomial
!    to be divided.
!
!    Input, integer ( kind = 4 ) NB, the degree of polynomial B.
!
!    Input, integer ( kind = 4 ) B(0:NB), the coefficients of the divisor
!    polynomial.
!
!    Output, integer ( kind = 4 ) NQ, the degree of polynomial Q.
!    If the divisor polynomial is zero, NQ is returned as -1.
!
!    Output, integer ( kind = 4 ) Q(0:NA-NB), contains the quotient of A/B.
!    If A and B have full degree, Q should be dimensioned Q(0:NA-NB).
!    In any case, Q(0:NA) should be enough.
!
!    Output, integer ( kind = 4 ) NR, the degree of polynomial R.
!    If the divisor polynomial is zero, NR is returned as -1.
!
!    Output, integer ( kind = 4 ) R(0:NB-1), contains the remainder of A/B.
!    If B has full degree, R should be dimensioned R(0:NB-1).
!    Otherwise, R will actually require less space.
!
  implicit none

  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb

  integer ( kind = 4 ) a(0:na)
  integer ( kind = 4 ) a2(0:na)
  integer ( kind = 4 ) b(0:nb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) na2
  integer ( kind = 4 ) nb2
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) q(0:*)
  integer ( kind = 4 ) r(0:*)

  call i4poly_degree ( na, a, na2 )

  call i4poly_degree ( nb, b, nb2 )

  if ( b(nb2) == 0 ) then
    nq = -1
    nr = -1
    return
  end if

  a2(0:na2) = a(0:na2)

  nq = na2 - nb2
  nr = nb2 - 1

  do i = nq, 0, -1
    q(i) = a2(i+nb2) / b(nb2)
    a2(i+nb2) = 0
    a2(i:i+nb2-1) = a2(i:i+nb2-1) - q(i) * b(0:nb2-1)
  end do

  r(0:nr) = a2(0:nr)

  return
end
subroutine i4poly_mul ( na, a, nb, b, c )

!*****************************************************************************80
!
!! I4POLY_MUL computes the product of two I4POLY's.
!
!  Discussion:
!
!    The polynomials are in power sum form.
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
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
!    Input, integer ( kind = 4 ) NA, the degree of polynomial A.
!
!    Input, integer ( kind = 4 ) A(0:NA), the coefficients of the first
!    polynomial factor.
!
!    Input, integer ( kind = 4 ) NB, the degree of polynomial B.
!
!    Input, integer ( kind = 4 ) B(0:NB), the coefficients of the
!    second polynomial factor.
!
!    Output, integer ( kind = 4 ) C(0:NA+NB), the coefficients of A * B.
!
  implicit none

  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb

  integer ( kind = 4 ) a(0:na)
  integer ( kind = 4 ) b(0:nb)
  integer ( kind = 4 ) c(0:na+nb)
  integer ( kind = 4 ) d(0:na+nb)
  integer ( kind = 4 ) i

  d(0:na+nb) = 0

  do i = 0, na
    d(i:i+nb) = d(i:i+nb) + a(i) * b(0:nb)
  end do

  c(0:na+nb) = d(0:na+nb)

  return
end
subroutine i4poly_print ( n, a, title )

!*****************************************************************************80
!
!! I4POLY_PRINT prints an I4POLY.
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the polynomial of A.
!
!    Input, integer ( kind = 4 ) A(0:N), the polynomial coefficients.
!    A(0) is the constant term and
!    A(N) is the coefficient of X**N.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) mag
  integer ( kind = 4 ) n2
  character plus_minus
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  call i4poly_degree ( n, a, n2 )

  if ( a(n2) < 0 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( a(n2) )

  if ( 2 <= n2 ) then
    write ( *, '( ''  p(x) = '', a1, i8, '' * x ^ '', i3 )' ) &
      plus_minus, mag, n2
  else if ( n2 == 1 ) then
    write ( *, '( ''  p(x) = '', a1, i8, '' * x'' )' ) plus_minus, mag
  else if ( n2 == 0 ) then
    write ( *, '( ''  p(x) = '', a1, i8 )' ) plus_minus, mag
  end if

  do i = n2-1, 0, -1

    if ( a(i) < 0.0D+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( a(i) )

    if ( mag /= 0 ) then

      if ( 2 <= i ) then
        write ( *, ' ( ''         '', a1, i8, '' * x ^ '', i3 )' ) &
          plus_minus, mag, i
      else if ( i == 1 ) then
        write ( *, ' ( ''         '', a1, i8, '' * x'' )' ) plus_minus, mag
      else if ( i == 0 ) then
        write ( *, ' ( ''         '', a1, i8 )' ) plus_minus, mag
      end if
    end if

  end do

  return
end
subroutine i4poly_to_i4 ( n, a, x, value )

!*****************************************************************************80
!
!! I4POLY_TO_I4 evaluates an I4POLY.
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the degree of the polynomial.
!
!    Input, integer ( kind = 4 ) A(0:N), the polynomial coefficients.
!    A(0) is the constant term and
!    A(N) is the coefficient of X**N.
!
!    Input, integer ( kind = 4 ) X, the point at which the polynomial is
!    to be evaluated.
!
!    Output, integer ( kind = 4 ) VALUE, the value of the polynomial.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) value
  integer ( kind = 4 ) x

  value = 0

  do i = n, 0, -1
    value = value * x + a(i)
  end do

  return
end
function i4vec_ascends ( n, x )

!*****************************************************************************80
!
!! I4VEC_ASCENDS determines if an I4VEC is (weakly) ascending.
!
!  Discussion:
!
!    The sequence is not required to be strictly ascending.
!
!  Example:
!
!    X = ( -8, 1, 2, 3, 7, 7, 9 )
!
!    I4VEC_ASCENDS = TRUE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, integer ( kind = 4 ) X(N), the array to be examined.
!
!    Output, logical I4VEC_ASCENDS, is TRUE if the entries of X ascend.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical i4vec_ascends
  integer ( kind = 4 ) x(n)

  i4vec_ascends = .false.

  do i = 1, n-1
    if ( x(i+1) < x(i) ) then
      return
    end if
  end do

  i4vec_ascends = .true.

  return
end
subroutine i4vec_backtrack ( n, maxstack, stack, x, indx, k, nstack, ncan )

!*****************************************************************************80
!
!! I4VEC_BACKTRACK supervises a backtrack search for an I4VEC.
!
!  Discussion:
!
!    The routine tries to construct an integer vector one index at a time,
!    using possible candidates as supplied by the user.
!
!    At any time, the partially constructed vector may be discovered to be
!    unsatisfactory, but the routine records information about where the
!    last arbitrary choice was made, so that the search can be
!    carried out efficiently, rather than starting out all over again.
!
!    First, call the routine with INDX = 0 so it can initialize itself.
!
!    Now, on each return from the routine, if INDX is:
!      1, you've just been handed a complete candidate vector;
!         Admire it, analyze it, do what you like.
!      2, please determine suitable candidates for position X(K).
!         Return the number of candidates in NCAN(K), adding each
!         candidate to the end of STACK, and increasing NSTACK.
!      3, you're done.  Stop calling the routine;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 July 2004
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
!    Input, integer ( kind = 4 ) N, the number of positions to be filled in
!    the vector.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the maximum length of the stack.
!
!    Input, integer ( kind = 4 ) STACK(MAXSTACK), a list of all current
!    candidates for all positions 1 through K.
!
!    Input/output, integer ( kind = 4 ) X(N), the partially filled in
!    candidate vector.
!
!    Input/output, integer ( kind = 4 ) INDX, a communication flag.
!    On input,
!      0, to begin a backtracking search.
!      2, the requested candidates for position K have been added to
!      STACK, and NCAN(K) was updated.
!    On output:
!      1, a complete output vector has been determined and returned in X(1:N);
!      2, candidates are needed for position X(K);
!      3, no more possible vectors exist.
!
!    Input/output, integer ( kind = 4 ) K, the index in X that we are trying
!    to fill.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!
!    Input/output, integer ( kind = 4 ) NCAN(N), lists the current number
!    of candidates for all positions 1 through K.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) indx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ncan(n)
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) stack(maxstack)
  integer ( kind = 4 ) x(n)
!
!  If this is the first call, request a candidate for position 1.
!
  if ( indx == 0 ) then
    k = 1
    nstack = 0
    indx = 2
    return
  end if
!
!  Examine the stack.
!
  do
!
!  If there are candidates for position K, take the first available
!  one off the stack, and increment K.
!
!  This may cause K to reach the desired value of N, in which case
!  we need to signal the user that a complete set of candidates
!  is being returned.
!
    if ( 0 < ncan(k) ) then

      x(k) = stack(nstack)
      nstack = nstack - 1

      ncan(k) = ncan(k) - 1

      if ( k /= n ) then
        k = k + 1
        indx = 2
      else
        indx = 1
      end if

      exit
!
!  If there are no candidates for position K, then decrement K.
!  If K is still positive, repeat the examination of the stack.
!
    else

      k = k - 1

      if ( k <= 0 ) then
        indx = 3
        exit
      end if

    end if

  end do

  return
end
function i4vec_descends ( n, x )

!*****************************************************************************80
!
!! I4VEC_DESCENDS determines if an I4VEC is decreasing.
!
!  Example:
!
!    X = ( 9, 7, 7, 3, 2, 1, -8 )
!
!    I4VEC_DESCENDS = TRUE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, integer ( kind = 4 ) X(N), the array to be examined.
!
!    Output, logical I4VEC_DESCEND, is TRUE if the entries of X descend.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical i4vec_descends
  integer ( kind = 4 ) x(n)

  i4vec_descends = .false.

  do i = 1, n-1
    if ( x(i) < x(i+1) ) then
      return
    end if
  end do

  i4vec_descends = .true.

  return
end
subroutine i4vec_frac ( n, a, k, afrac )

!*****************************************************************************80
!
!! I4VEC_FRAC searches for the K-th smallest element in an I4VEC.
!
!  Discussion:
!
!    Hoare's algorithm is used.
!
!  Modified:
!
!    17 July 2000
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, integer ( kind = 4 ) A(N), array to search.  On output,
!    the elements of A have been somewhat rearranged.
!
!    Input, integer ( kind = 4 ) K, the fractile to be sought.  If K = 1, the
!    minimum entry is sought.  If K = N, the maximum is sought.
!    Other values of K search for the entry which is K-th in size.
!    K must be at least 1, and no greater than N.
!
!    Output, integer ( kind = 4 ) AFRAC, the value of the K-th fractile of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) afrac
  integer ( kind = 4 ) iryt
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) left

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal nonpositive value of N = ', n
    stop
  end if

  if ( k <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal nonpositive value of K = ', k
    stop
  end if

  if ( n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal N < K, K = ', k
    stop
  end if

  left = 1
  iryt = n

  do

    if ( iryt <= left ) then
      afrac = a(k)
      exit
    end if

    ix = a(k)
    i = left
    j = iryt

    do

      if ( j < i ) then

        if ( j < k ) then
          left = i
        end if

        if ( k < i ) then
          iryt = j
        end if

        exit

      end if
!
!  Find I so that IX <= A(I).
!
      do while ( a(i) < ix )
        i = i + 1
      end do
!
!  Find J so that A(J) <= IX.
!
      do while ( ix < a(j) )
        j = j - 1
      end do

      if ( i <= j ) then
        call i4_swap ( a(i), a(j) )
        i = i + 1
        j = j - 1
      end if

    end do

  end do

  return
end
subroutine i4vec_heap_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_D reorders an I4VEC into an descending heap.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    A descending heap is an array A with the property that, for every index J,
!    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n/2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(m) < a(m+1) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) <= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
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
function i4vec_index ( n, a, aval )

!*****************************************************************************80
!
!! I4VEC_INDEX returns the location of the first occurrence of a given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be searched.
!
!    Input, integer ( kind = 4 ) AVAL, the value to be indexed.
!
!    Output, integer ( kind = 4 ) I4VEC_INDEX, the first location in A
!    which has the value AVAL, or -1 if no such index exists.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4vec_index

  do i = 1, n
    if ( a(i) == aval ) then
      i4vec_index = i
      return
    end if
  end do

  i4vec_index = -1

  return
end
function i4vec_maxloc_last ( n, x )

!*****************************************************************************80
!
!! I4VEC_MAXLOC_LAST returns the index of the last maximal I4VEC entry.
!
!  Example:
!
!    X = ( 5, 1, 2, 5, 0, 5, 3 )
!
!    I4VEC_MAXLOC_LAST = 6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, integer ( kind = 4 ) X(N), the array to be examined.
!
!    Output, integer ( kind = 4 ) I4VEC_MAXLOC_LAST, the index of the
!    last element of X of maximal value.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4vec_maxloc_last
  integer ( kind = 4 ) i4vec_maxval_last
  integer ( kind = 4 ) x(n)

  i4vec_maxloc_last = 0

  do i = 1, n
    if ( i == 1 ) then
      i4vec_maxloc_last = 1
      i4vec_maxval_last = x(1)
    else if ( i4vec_maxval_last <= x(i) ) then
      i4vec_maxloc_last = i
      i4vec_maxval_last = x(i)
    end if
  end do

  return
end
function i4vec_pairwise_prime ( n, a )

!*****************************************************************************80
!
!! I4VEC_PAIRWISE_PRIME checks whether an I4VEC is pairwise prime.
!
!  Discussion:
!
!    Two positive integers I and J are pairwise prime if they have no common
!    factor greater than 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values to check.
!
!    Input, integer ( kind = 4 ) A(N), the vector of integers.
!
!    Output, logical I4VEC_PAIRWISE_PRIME, is TRUE if the vector of integers
!    is pairwise prime.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_gcd
  logical i4vec_pairwise_prime
  integer ( kind = 4 ) j

  i4vec_pairwise_prime = .false.

  do i = 1, n
    do j = i+1, n
      if ( i4_gcd ( a(i), a(j) ) /= 1 ) then
        return
      end if
    end do
  end do

  i4vec_pairwise_prime = .true.

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC, with an optional title.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(i8,i10)' ) i, a(i)
  end do

  return
end
subroutine i4vec_print_some ( n, a, max_print, title )

!*****************************************************************************80
!
!! I4VEC_PRINT_SOME prints "some" of an I4VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(i8,2x,i10)' ) i, a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print-2
      write ( *, '(i8,2x,i10)' ) i, a(i)
    end do
    write ( *, '(a)' ) '......  ..............'
    i = n
    write ( *, '(i8,2x,i10)' ) i, a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i8,2x,i10)' ) i, a(i)
    end do
    i = max_print
    write ( *, '(i8,2x,i10,2x,a)' ) i, a(i), '...more entries...'

  end if

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
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be reversed.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n/2
    call i4_swap ( a(i), a(n+1-i) )
  end do

  return
end
subroutine i4vec_sort_bubble_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_BUBBLE_A ascending sorts an I4VEC using bubble sort.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n-1
    do j = i+1, n
      if ( a(j) < a(i) ) then
        call i4_swap ( a(i), a(j) )
      end if
    end do
  end do

  return
end
subroutine i4vec_sort_heap_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) n1

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into descending heap form.
!
  call i4vec_heap_d ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
!  Move it to position A(N).
!
  call i4_swap ( a(1), a(n) )
!
!  Consider the diminished heap of size N1.
!
  do n1 = n - 1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call i4vec_heap_d ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    call i4_swap ( a(1), a(n1) )

  end do

  return
end
subroutine i4vec_sort_heap_index_d ( n, a, indx )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an I4VEC.
!
!  Discussion:
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(I)), I = 1 to N is sorted,
!
!    after which A(I), I = 1 to N is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), contains the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  call i4vec_indicator ( n, indx )

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j+1)) < a(indx(j)) ) then
          j = j + 1
        end if
      end if

      if ( a(indx(j)) < aval ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine i4vec_transpose_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
!
!  Example:
!
!    A = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 /)
!    TITLE = 'My vector:  '
!
!    My vector:      1    2    3    4    5
!                    6    7    8    9   10
!                   11
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 July 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = 11 ) string
  character ( len = * ) title
  integer ( kind = 4 ) title_len

  if ( 0 < len_trim ( title ) ) then

    title_len = len ( title )

    write ( string, '(a,i3,a)' ) '(', title_len, 'x,5i12)'

    do ilo = 1, n, 5
      ihi = min ( ilo + 5 - 1, n )
      if ( ilo == 1 ) then
        write ( *, '(a, 5i12)' ) title, a(ilo:ihi)
      else
        write ( *, string      )        a(ilo:ihi)
      end if
    end do

  else

    do ilo = 1, n, 5
      ihi = min ( ilo + 5 - 1, n )
      write ( *, '(5i12)' ) a(ilo:ihi)
    end do

  end if

  return
end
subroutine i4vec_uniform ( n, a, b, seed, x )

!*****************************************************************************80
!
!! I4VEC_UNIFORM returns a scaled pseudorandom I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer ( kind = 4 ) values.
!
!    The pseudorandom numbers should be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) X(N), a vector of numbers between A and B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value
  integer ( kind = 4 ) x(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + huge ( seed )
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

    x(i) = value

  end do

  return
end
subroutine index_box2_next_2d ( n1, n2, ic, jc, i, j, more )

!*****************************************************************************80
!
!! INDEX_BOX2_NEXT_2D produces index vectors on the surface of a box in 2D.
!
!  Discussion:
!
!    The box is has center at (IC,JC), and has half-widths N1 and N2.
!    The index vectors are exactly those which are between (IC-N1,JC-N1) and
!    (IC+N1,JC+N2) with the property that at least one of I and J
!    is an "extreme" value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the half-widths of the box, that is,
!    the maximum distance allowed between (IC,JC) and (I,J).
!
!    Input, integer ( kind = 4 ) IC, JC, the central cell of the box.
!
!    Input/output, integer ( kind = 4 ) I, J.  On input, the previous index set.
!    On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I and J are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  logical more
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  if ( .not. more ) then
    more = .true.
    i = ic - n1
    j = jc - n2
    return
  end if

  if ( i == ic + n1 .and. j == jc + n2 ) then
    more = .false.
    return
  end if
!
!  Increment J.
!
  j = j + 1
!
!  Check J.
!
  if ( jc + n2 < j ) then
    j = jc - n2
    i = i + 1
  else if ( j < jc + n2 .and. ( i == ic - n1 .or. i == ic + n1 ) ) then
    return
  else
    j = jc + n2
    return
  end if

  return
end
subroutine index_box2_next_3d ( n1, n2, n3, ic, jc, kc, i, j, k, more )

!*****************************************************************************80
!
!! INDEX_BOX2_NEXT_3D produces index vectors on the surface of a box in 3D.
!
!  Discussion:
!
!    The box has a central cell of (IC,JC,KC), with a half widths of
!    (N1,N2,N3).  The index vectors are exactly those between
!    (IC-N1,JC-N2,KC-N3) and (IC+N1,JC+N2,KC+N3) with the property that
!    at least one of I, J, and K is an "extreme" value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, N3, the "half widths" of the box,
!    that is, the maximum distances from the central cell allowed for
!    I, J and K.
!
!    Input, integer ( kind = 4 ) IC, JC, KC, the central cell of the box.
!
!    Input/output, integer ( kind = 4 ) I, J, K.  On input, the previous
!    index set.  On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I, J, and K are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  logical more
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  if ( .not. more ) then
    more = .true.
    i = ic - n1
    j = jc - n2
    k = kc - n3
    return
  end if

  if ( i == ic + n1 .and. j == jc + n2 .and. k == kc + n3 ) then
    more = .false.
    return
  end if
!
!  Increment K.
!
  k = k + 1
!
!  Check K.
!
  if ( kc + n3 < k ) then
    k = kc - n3
    j = j + 1
  else if ( k < kc + n3 .and. &
    ( i == ic - n1 .or. i == ic + n1 .or. &
      j == jc - n2 .or. j == jc + n2 ) ) then
    return
  else
    k = kc + n3
    return
  end if
!
!  Check J.
!
  if ( jc + n2 < j ) then
    j = jc - n2
    i = i + 1
  else if ( j < jc + n2 .and. &
    ( i == ic - n1 .or. i == ic + n1 .or. &
      k == kc - n3 .or. k == kc + n3 ) ) then
    return
  else
    j = jc + n2
    return
  end if

  return
end
subroutine index_box_next_2d ( n1, n2, i, j, more )

!*****************************************************************************80
!
!! INDEX_BOX_NEXT_2D produces index vectors on the surface of a box in 2D.
!
!  Discussion:
!
!    The index vectors are exactly those which are between (1,1) and
!    (N1,N2) with the property that at least one of I and J
!    is an "extreme" value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the "dimensions" of the box, that is,
!    the maximum values allowed for I and J.  The minimum values are
!    assumed to be 1.
!
!    Input/output, integer ( kind = 4 ) I, J.  On input, the previous index set.
!    On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I and J are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  logical more
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  if ( .not. more ) then
    more = .true.
    i = 1
    j = 1
    return
  end if

  if ( i == n1 .and. j == n2 ) then
    more = .false.
    return
  end if
!
!  Increment J.
!
  j = j + 1
!
!  Check J.
!
  if ( n2 < j ) then
    j = 1
    i = i + 1
  else if ( j < n2 .and. ( i == 1 .or. i == n1 ) ) then
    return
  else
    j = n2
    return
  end if

  return
end
subroutine index_box_next_3d ( n1, n2, n3, i, j, k, more )

!*****************************************************************************80
!
!! INDEX_BOX_NEXT_3D produces index vectors on the surface of a box in 3D.
!
!  Discussion:
!
!    The index vectors are exactly those which are between (1,1,1) and
!    (N1,N2,N3) with the property that at least one of I, J, and K
!    is an "extreme" value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, N3, the "dimensions" of the box,
!    that is, the maximum values allowed for I, J and K.  The minimum values
!    are assumed to be 1.
!
!    Input/output, integer ( kind = 4 ) I, J, K.  On input, the previous
!    index set.  On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I, J, and K are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical more
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  if ( .not. more ) then
    more = .true.
    i = 1
    j = 1
    k = 1
    return
  end if

  if ( i == n1 .and. j == n2 .and. k == n3 ) then
    more = .false.
    return
  end if
!
!  Increment K.
!
  k = k + 1
!
!  Check K.
!
  if ( n3 < k ) then
    k = 1
    j = j + 1
  else if ( k < n3 .and. &
    ( i == 1 .or. i == n1 .or. j == 1 .or. j == n2 ) ) then
    return
  else
    k = n3
    return
  end if
!
!  Check J.
!
  if ( n2 < j ) then
    j = 1
    i = i + 1
  else if ( j < n2 .and. &
    ( i == 1 .or. i == n1 .or. k == 1 .or. k == n3 ) ) then
    return
  else
    j = n2
    return
  end if

  return
end
subroutine index_next0 ( n, hi, a, more )

!*****************************************************************************80
!
!! INDEX_NEXT0 generates all index vectors within given upper limits.
!
!  Discussion:
!
!    The index vectors are generated in such a way that the reversed
!    sequences are produced in lexicographic order.
!
!  Example:
!
!    N = 3,
!    HI = 3
!
!    1   2   3
!    ---------
!    1   1   1
!    2   1   1
!    3   1   1
!    1   2   1
!    2   2   1
!    3   2   1
!    1   3   1
!    2   3   1
!    3   3   1
!    1   1   2
!    2   1   2
!    3   1   2
!    1   2   2
!    2   2   2
!    3   2   2
!    1   3   2
!    2   3   2
!    3   3   2
!    1   1   3
!    2   1   3
!    3   1   3
!    1   2   3
!    2   2   3
!    3   2   3
!    1   3   3
!    2   3   3
!    3   3   3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, integer ( kind = 4 ) HI, the upper limit for the array indices.
!    The lower limit is implicitly 1 and HI must be at least 1.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On startup calls, with MORE = FALSE, the input value of A
!    doesn't matter, because the routine initializes it.
!    On calls with MORE = TRUE, the input value of A must be
!    the output value of A from the previous call.  (In other words,
!    just leave it alone!).
!    On output, A contains the successor set of indices to the input
!    value.
!
!    Input/output, logical MORE.  Set this variable FALSE before
!    the first call.  Normally, MORE will be returned TRUE but
!    once all the vectors have been generated, MORE will be
!    reset to FALSE and you should stop calling the program.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inc
  logical more

  if ( .not. more ) then

    a(1:n) = 1

    if ( hi < 1 ) then
      more = .false.
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INDEX_NEXT0 - Fatal error!'
      write ( *, '(a,i8)' ) '  HI is ', hi
      write ( *, '(a)' ) '  but HI must be at least 1.'
      stop
    end if

  else

    inc = 1

    do while ( hi <= a(inc) )
      a(inc) = 1
      inc = inc + 1
    end do

    a(inc) = a(inc) + 1

  end if
!
!  See if there are more entries to compute.
!
  more = .false.

  do i = 1, n
    if ( a(i) < hi ) then
      more = .true.
    end if
  end do

  return
end
subroutine index_next1 ( n, hi, a, more )

!*****************************************************************************80
!
!! INDEX_NEXT1 generates all index vectors within given upper limits.
!
!  Discussion:
!
!    The index vectors are generated in such a way that the reversed
!    sequences are produced in lexicographic order.
!
!  Example:
!
!    N = 3,
!    HI(1) = 4, HI(2) = 2, HI(3) = 3
!
!    1   2   3
!    ---------
!    1   1   1
!    2   1   1
!    3   1   1
!    4   1   1
!    1   2   1
!    2   2   1
!    3   2   1
!    4   2   1
!    1   1   2
!    2   1   2
!    3   1   2
!    4   1   2
!    1   2   2
!    2   2   2
!    3   2   2
!    4   2   2
!    1   1   3
!    2   1   3
!    3   1   3
!    4   1   3
!    1   2   3
!    2   2   3
!    3   2   3
!    4   2   3
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
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, integer ( kind = 4 ) HI(N), the upper limits for the array indices.
!    The lower limit is implicitly 1, and each HI(I) should be at least 1.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On startup calls, with MORE = FALSE, the input value of A
!    doesn't matter, because the routine initializes it.
!    On calls with MORE = TRUE, the input value of A must be
!    the output value of A from the previous call.  (In other words,
!    just leave it alone!).
!    On output, A contains the successor set of indices to the input
!    value.
!
!    Input/output, logical MORE.  Set this variable FALSE before
!    the first call.  Normally, MORE will be returned TRUE but
!    once all the vectors have been generated, MORE will be
!    reset FALSE and you should stop calling the program.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) hi(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inc
  logical more

  if ( .not. more ) then

    a(1:n) = 1

    do i = 1, n
      if ( hi(i) < 1 ) then
        more = .false.
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INDEX_NEXT1 - Fatal error!'
        write ( *, '(a,i8,a,i8)' ) '  Entry ', i, ' of HI is ', hi(i)
        write ( *, '(a)' ) '  but all entries must be at least 1.'
        stop
      end if
    end do

  else

    inc = 1

    do while ( hi(inc) <= a(inc) )
      a(inc) = 1
      inc = inc + 1
    end do

    a(inc) = a(inc) + 1

  end if
!
!  See if there are more entries to compute.
!
  more = .false.

  do i = 1, n
    if ( a(i) < hi(i) ) then
      more = .true.
    end if
  end do

  return
end
subroutine index_next2 ( n, lo, hi, a, more )

!*****************************************************************************80
!
!! INDEX_NEXT2 generates all index vectors within given lower and upper limits.
!
!  Example:
!
!    N = 3,
!    LO(1) = 1, LO(2) = 10, LO(3) = 4
!    HI(1) = 2, HI(2) = 11, HI(3) = 6
!
!    1   2   3
!    ---------
!    1  10   4
!    2  10   4
!    1  11   4
!    2  11   4
!    1  10   5
!    2  10   5
!    1  11   5
!    2  11   5
!    1  10   6
!    2  10   6
!    1  11   6
!    2  11   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.  The rank of
!    the object being indexed.
!
!    Input, integer ( kind = 4 ) LO(N), HI(N), the lower and upper limits
!    for the array indices.  LO(I) should be less than or equal to HI(I),
!    for each I.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On startup calls, with MORE = FALSE, the input value of A
!    doesn't matter, because the routine initializes it.
!    On calls with MORE = TRUE, the input value of A must be
!    the output value of A from the previous call.  (In other words,
!    just leave it alone!).
!    On output, A contains the successor set of indices to the input
!    value.
!
!    Input/output, logical MORE.  Set this variable FALSE before
!    the first call.  Normally, MORE will be returned TRUE but
!    once all the vectors have been generated, MORE will be
!    reset FALSE and you should stop calling the program.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) hi(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lo(n)
  logical more

  if ( .not. more ) then

    a(1:n) = lo(1:n)

    do i = 1, n
      if ( hi(i) < lo(i) ) then
        more = .false.
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INDEX_NEXT2 - Fatal error!'
        write ( *, '(a,i8,a,i8)' ) '  Entry ', i, ' of HI is ', hi(i)
        write ( *, '(a,i8,a,i8)' ) '  Entry ', i, ' of LO is ', lo(i)
        write ( *, '(a)' ) '  but LO(I) <= HI(I) is required.'
        stop
      end if
    end do

  else

    inc = 1

    do while ( hi(inc) <= a(inc) )
      a(inc) = lo(inc)
      inc = inc + 1
    end do

    a(inc) = a(inc) + 1

  end if
!
!  See if there are more entries to compute.
!
  more = .false.

  do i = 1, n
    if ( a(i) < hi(i) ) then
      more = .true.
    end if
  end do

  return
end
subroutine index_rank0 ( n, hi, a, rank )

!*****************************************************************************80
!
!! INDEX_RANK0 ranks an index vector within given upper limits.
!
!  Example:
!
!    N = 3,
!    HI = 3
!    A = ( 3, 1, 2 )
!
!    RANK = 12
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, integer ( kind = 4 ) HI, the upper limit for the array indices.
!    The lower limit is implicitly 1, and HI should be at least 1.
!
!    Input, integer ( kind = 4 ) A(N), the index vector to be ranked.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the index vector, or -1 if A
!    is not a legal index.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) range
  integer ( kind = 4 ) rank

  rank = -1
  do i = 1, n
    if ( a(i) < 1 .or. hi < a(i) ) then
      return
    end if
  end do

  rank = 0
  do i = n, 1, -1
    rank = hi * rank + a(i)
  end do

  rank = 1
  range = 1
  do i = 1, n
    rank = rank + ( a(i) - 1 ) * range
    range = range * hi
  end do

  return
end
subroutine index_rank1 ( n, hi, a, rank )

!*****************************************************************************80
!
!! INDEX_RANK1 ranks an index vector within given upper limits.
!
!  Example:
!
!    N = 3,
!    HI(1) = 4, HI(2) = 2, HI(3) = 3
!    A = ( 4, 1, 2 )
!
!    RANK = 12
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, integer ( kind = 4 ) HI(N), the upper limits for the array indices.
!    The lower limit is implicitly 1, and each HI(I) should be at least 1.
!
!    Input, integer ( kind = 4 ) A(N), the index to be ranked.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the index vector, or -1 if A
!    is not a legal index.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) hi(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) range
  integer ( kind = 4 ) rank

  rank = -1
  do i = 1, n
    if ( a(i) < 1 .or. hi(i) < a(i) ) then
      return
    end if
  end do

  rank = 0
  do i = n, 1, -1
    rank = hi(i) * rank + a(i)
  end do

  rank = 1
  range = 1
  do i = 1, n
    rank = rank + ( a(i) - 1 ) * range
    range = range * hi(i)
  end do

  return
end
subroutine index_rank2 ( n, lo, hi, a, rank )

!*****************************************************************************80
!
!! INDEX_RANK2 ranks an index vector within given lower and upper limits.
!
!  Example:
!
!    N = 3,
!    LO(1) = 1, LO(2) = 10, LO(3) = 4
!    HI(1) = 2, HI(2) = 11, HI(3) = 6
!    A = ( 1, 11, 5 )
!
!    RANK = 7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, integer ( kind = 4 ) LO(N), HI(N), the lower and upper limits
!    for the array indices.  LO(I) should be less than or equal to HI(I),
!    for each I.
!
!    Input, integer ( kind = 4 ) A(N), the index vector to be ranked.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the index vector,
!    or -1 if A is not a legal index vector.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) hi(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lo(n)
  integer ( kind = 4 ) range
  integer ( kind = 4 ) rank

  do i = 1, n
    if ( a(i) < lo(i) .or. hi(i) < a(i) ) then
      rank = -1
      return
    end if
  end do

  rank = 1
  range = 1
  do i = 1, n
    rank = rank + ( a(i) - lo(i) ) * range
    range = range * ( hi(i) + 1 - lo(i) )
  end do

  return
end
subroutine index_unrank0 ( n, hi, rank, a )

!*****************************************************************************80
!
!! INDEX_UNRANK0 unranks an index vector within given upper limits.
!
!  Example:
!
!    N = 3,
!    HI = 3
!    RANK = 12
!
!    A = ( 3, 1, 2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, integer ( kind = 4 ) HI, the upper limit for the array indices.
!    The lower limit is implicitly 1, and HI should be at least 1.
!
!    Input, integer ( kind = 4 ) RANK, the rank of the desired index vector.
!
!    Output, integer ( kind = 4 ) A(N), the index vector of the given rank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) range
  integer ( kind = 4 ) rank

  a(1:n) = 0
!
!  The rank might be too small.
!
  if ( rank < 1 ) then
    return
  end if

  range = hi**n
!
!  The rank might be too large.
!
  if ( range < rank ) then
    return
  end if

  k = rank - 1
  do i = n, 1, -1
    range = range / hi
    j = k / range
    a(i) = j + 1
    k = k - j * range
  end do

  return
end
subroutine index_unrank1 ( n, hi, rank, a )

!*****************************************************************************80
!
!! INDEX_UNRANK1 unranks an index vector within given upper limits.
!
!  Example:
!
!    N = 3,
!    HI(1) = 4, HI(2) = 2, HI(3) = 3
!    RANK = 11
!
!    A = ( 3, 1, 2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, integer ( kind = 4 ) HI(N), the upper limits for the array indices.
!    The lower limit is implicitly 1, and each HI(I) should be at least 1.
!
!    Input, integer ( kind = 4 ) RANK, the rank of the desired index vector.
!
!    Output, integer ( kind = 4 ) A(N), the index vector of the given rank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) hi(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) range
  integer ( kind = 4 ) rank

  a(1:n) = 0
!
!  The rank might be too small.
!
  if ( rank < 1 ) then
    return
  end if

  range = product ( hi )
!
!  The rank might be too large.
!
  if ( range < rank ) then
    return
  end if

  k = rank - 1
  do i = n, 1, -1
    range = range / hi(i)
    j = k / range
    a(i) = j + 1
    k = k - j * range
  end do

  return
end
subroutine index_unrank2 ( n, lo, hi, rank, a )

!*****************************************************************************80
!
!! INDEX_UNRANK2 unranks an index vector within given lower and upper limits.
!
!  Example:
!
!    N = 3,
!    LO(1) = 1, LO(2) = 10, LO(3) = 4
!    HI(1) = 2, HI(2) = 11, HI(3) = 6
!    RANK = 7
!
!    A = ( 1, 11, 5 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in A.
!
!    Input, integer ( kind = 4 ) LO(N), HI(N), the lower and upper limits
!    for the array indices.  It should be the case that LO(I) <= HI(I)
!    for each I.
!
!    Input, integer ( kind = 4 ) RANK, the rank of the desired index.
!
!    Output, integer ( kind = 4 ) A(N), the index vector of the given rank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) hi(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lo(n)
  integer ( kind = 4 ) range
  integer ( kind = 4 ) rank

  a(1:n) = 0
!
!  The rank might be too small.
!
  if ( rank < 1 ) then
    return
  end if

  range = 1
  do i = 1, n
    range = range * ( hi(i) + 1 - lo(i) )
  end do
!
!  The rank might be too large.
!
  if ( range < rank ) then
    return
  end if

  k = rank - 1
  do i = n, 1, -1
    range = range / ( hi(i) + 1 - lo(i) )
    j = k / range
    a(i) = j + lo(i)
    k = k - j * range
  end do

  return
end
subroutine ins_perm ( n, ins, p )

!*****************************************************************************80
!
!! INS_PERM computes a permutation from its inversion sequence.
!
!  Discussion:
!
!    For a given permutation P acting on objects 1 through N, the
!    inversion sequence INS is defined as:
!
!      INS(1) = 0
!      INS(I) = number of values J < I for which P(I) < P(J).
!
!  Example:
!
!    Input:
!
!      ( 0, 0, 2, 1, 3 )
!
!    Output:
!
!      ( 3, 5, 1, 4, 2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer, 1986,
!    ISBN: 0387963472,
!    LC: QA164.S79.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input, integer ( kind = 4 ) INS(N), the inversion sequence of a
!    permutation.  It must be the case that 0 <= INS(I) < I for I = 1 to N.
!
!    Output, integer ( kind = 4 ) P(N), the permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ins(n)
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(n)

  call i4vec_indicator ( n, p )

  do i = n, 2, -1

    itemp = p(i-ins(i))

    do j = i-ins(i), i-1
      p(j) = p(j+1)
    end do

    p(i) = itemp

  end do

  return
end
subroutine inverse_mod_n ( b, n, y )

!*****************************************************************************80
!
!! INVERSE_MOD_N computes the inverse of B mod N.
!
!  Discussion:
!
!    If
!
!      Y = inverse_mod_n ( B, N )
!
!    then
!
!      mod ( B * Y, N ) = 1
!
!    The value Y will exist if and only if B and N are relatively prime.
!
!  Examples:
!
!    B  N  Y
!
!    1  2  1
!
!    1  3  1
!    2  3  2
!
!    1  4  1
!    2  4  0
!    3  4  3
!
!    1  5  1
!    2  5  3
!    3  5  2
!    4  5  4
!
!    1  6  1
!    2  6  0
!    3  6  0
!    4  6  0
!    5  6  5
!
!    1  7  1
!    2  7  4
!    3  7  5
!    4  7  2
!    5  7  3
!    6  7  6
!
!    1  8  1
!    2  8  0
!    3  8  3
!    4  8  0
!    5  8  5
!    6  8  0
!    7  8  7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) B, the number whose inverse mod N is desired.
!    B should be positive.  Normally, B < N, but this is not required.
!
!    Input, integer ( kind = 4 ) N, the number with respect to which the
!    modulus is computed.  N should be positive.
!
!    Output, integer ( kind = 4 ) Y, the inverse of B mod N, or 0 if there
!    is not inverse for B mode N.  1 <= Y < N if the inverse exists.
!
  implicit none

  integer ( kind = 4 ) b
  integer ( kind = 4 ) b0
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r
  integer ( kind = 4 ) t
  integer ( kind = 4 ) t0
  integer ( kind = 4 ) temp
  integer ( kind = 4 ) y

  n0 = n
  b0 = b
  t0 = 0
  t = 1

  q = n / b
  r = n - q * b

  do while ( 0 < r )

    temp = t0 - q * t

    if ( 0 <= temp ) then
      temp = mod ( temp, n )
    end if

    if ( temp < 0 ) then
      temp = n - mod ( - temp, n )
    end if

    t0 = t
    t = temp
    n0 = b0
    b0 = r
    q = n0 / b0
    r = n0 - q * b0

  end do

  if ( b0 /= 1 ) then
    y = 0
    return
  end if

  y = mod ( t, n )

  return
end
subroutine involute_enum ( n, s )

!*****************************************************************************80
!
!! INVOLUTE_ENUM enumerates the involutions of N objects.
!
!  Discussion:
!
!    An involution is a permutation consisting only of fixed points and
!    pairwise transpositions.
!
!    An involution is its own inverse permutation.
!
!  Recursion:
!
!    S(0) = 1
!    S(1) = 1
!    S(N) = S(N-1) + (N-1) * S(N-2)
!
!  First values:
!
!     N         S(N)
!     0           1
!     1           1
!     2           2
!     3           4
!     4          10
!     5          26
!     6          76
!     7         232
!     8         764
!     9        2620
!    10        9496
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!
!    Output, integer ( kind = 4 ) S(0:N), the number of involutions of
!    0, 1, 2, ... N objects.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) s(0:n)

  if ( n < 0 ) then
    return
  end if

  s(0) = 1

  if ( n <= 0 ) then
    return
  end if

  s(1) = 1

  do i = 2, n
    s(i) = s(i-1) + ( i - 1 ) * s(i-2)
  end do

  return
end
subroutine jfrac_to_rfrac ( m, r, s, p, q )

!*****************************************************************************80
!
!! JFRAC_TO_RFRAC converts a J-fraction into a rational polynomial fraction.
!
!  Discussion:
!
!    The routine accepts a J-fraction:
!
!        R(1) / ( X + S(1)
!      + R(2) / ( X + S(2)
!      + R(3) / ...
!      + R(M) / ( X + S(M) )... ))
!
!    and returns the equivalent rational polynomial fraction:
!
!      P(1) + P(2) * X + ... + P(M) * X**(M-1)
!      -------------------------------------------------------
!      Q(1) + Q(2) * X + ... + Q(M) * X**(M-1) + Q(M+1) * X**M
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2000
!
!  Author:
!
!    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson,
!    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi,
!    John Rice, Henry Thatcher, Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, defines the number of P, R, and S
!    coefficients, and is one less than the number of Q
!    coefficients.
!
!    Input, real ( kind = 8 ) R(M), S(M), the coefficients defining
!    the J-fraction.
!
!    Output, real ( kind = 8 ) P(M), Q(M+1), the coefficients defining
!    the rational polynomial fraction.  The algorithm used normalizes
!    the coefficients so that Q(M+1) = 1.0.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m,m)
  real ( kind = 8 ) b(m,m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) p(m)
  real ( kind = 8 ) q(m+1)
  real ( kind = 8 ) r(m)
  real ( kind = 8 ) s(m)

  a(1,1) = r(1)
  b(1,1) = s(1)

  if ( 1 < m ) then

    do k = 2, m
      a(k,k) = r(1)
      b(k,k) = b(k-1,k-1) + s(k)
    end do

    a(1,2) = r(1) * s(2)
    b(1,2) = r(2) + s(1) * s(2)

    do k = 3, m
      a(1,k) = s(k) * a(1,k-1) + r(k) * a(1,k-2)
      a(k-1,k) = a(k-2,k-1) + s(k) * r(1)
      b(1,k) = s(k) * b(1,k-1) + r(k) * b(1,k-2)
      b(k-1,k) = b(k-2,k-1) + s(k) * b(k-1,k-1) + r(k)
    end do

    do k = 4, m
      do i = 2, k-2
        a(i,k) = a(i-1,k-1) + s(k) * a(i,k-1) + r(k) * a(i,k-2)
        b(i,k) = b(i-1,k-1) + s(k) * b(i,k-1) + r(k) * b(i,k-2)
      end do
    end do

  end if

  p(1:m) = a(1:m,m)

  q(1:m) = b(1:m,m)
  q(m+1) = 1.0D+00

  return
end
subroutine josephus ( n, m, k, x )

!*****************************************************************************80
!
!! JOSEPHUS returns the position X of the K-th man to be executed.
!
!  Discussion:
!
!    The classic Josephus problem concerns a circle of 41 men.
!    Every third man is killed and removed from the circle.  Counting
!    and executing continues until all are dead.  Where was the last
!    survivor sitting?
!
!    Note that the first person killed was sitting in the third position.
!    Moreover, when we get down to 2 people, and we need to count the
!    "third" one, we just do the obvious thing, which is to keep counting
!    around the circle until our count is completed.
!
!    The process may be regarded as generating a permutation of
!    the integers from 1 to N.  The permutation would be the execution
!    list, that is, the list of the executed men, by position number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    WW Rouse Ball,
!    Mathematical Recreations and Essays,
!    Macmillan, 1962, pages 32-36.
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 1, Fundamental Algorithms,
!    Addison Wesley, 1968, pages 158-159.
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 3, Sorting and Searching,
!    Addison Wesley, 1968,
!    ISBN: 0201896850,
!    LC: QA76.6.K64.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of men.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) M, the counting index.
!    M must not be zero.  Ordinarily, M is positive, and no greater than N.
!
!    Input, integer ( kind = 4 ) K, the index of the executed man of interest.
!    K must be between 1 and N.
!
!    Output, integer ( kind = 4 ) X, the position of the K-th man.
!    X will be between 1 and N.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JOSEPHUS - Fatal error!'
    write ( *, '(a)' ) '  N <= 0.'
    stop
  end if

  if ( m == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JOSEPHUS - Fatal error!'
    write ( *, '(a)' ) '  M = 0.'
    stop
  end if

  if ( k <= 0 .or. n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JOSEPHUS - Fatal error!'
    write ( *, '(a)' ) '  J <= 0 or N < K.'
    stop
  end if
!
!  In case M is bigger than N, or negative, get the
!  equivalent positive value between 1 and N.
!  You can skip this operation if 1 <= M <= N.
!
  m2 = i4_modp ( m, n )

  x = k * m2

  do while ( n < x )
    x = ( m2 * ( x - n ) - 1 ) / ( m2 - 1 )
  end do

  return
end
subroutine ksub_next ( n, k, a, more )

!*****************************************************************************80
!
!! KSUB_NEXT generates the subsets of size K from a set of size N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
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
!    Input, integer ( kind = 4 ) N, the size of the set from which subsets
!    are drawn.
!
!    Input, integer ( kind = 4 ) K, the desired size of the subsets.  K must
!    be between 0 and N.
!
!    Input/output, integer ( kind = 4 ) A(K).  A(I) is the I-th element of the
!    subset.  Thus A(I) will be an integer between 1 and N.
!    Note that the routine will return the values in A
!    in sorted order: 1 <= A(1) < A(2) < ... < A(K) <= N
!
!    Input/output, logical MORE.  Set MORE = FALSE before first call
!    for a new sequence of subsets.  It then is set and remains
!    TRUE as long as the subset computed on this call is not the
!    final one.  When the final subset is computed, MORE is set to
!    FALSE as a signal that the computation is done.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: m = 0
  integer ( kind = 4 ), save :: m2 = 0
  logical more
  integer ( kind = 4 ) n

  if ( k < 0 .or. n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_NEXT - Fatal error!'
    write ( *, '(a,i8)' ) 'N = ', n
    write ( *, '(a,i8)' ) 'K = ', k
    write ( *, '(a)' ) 'but 0 <= K <= N is required!'
    stop
  end if

  if ( .not. more ) then
    m2 = 0
    m = k
  else
    if ( m2 < n - m ) then
      m = 0
    end if
    m = m + 1
    m2 = a(k+1-m)
  end if

  do j = 1, m
    a(k+j-m) = m2 + j
  end do

  more = ( a(1) /= (n-k+1) )

  return
end
subroutine ksub_next2 ( n, k, a, in, iout )

!*****************************************************************************80
!
!! KSUB_NEXT2 generates the subsets of size K from a set of size N.
!
!  Discussion:
!
!    This routine uses the revolving door method.  It has no "memory".
!    It simply calculates the successor of the input set,
!    and will start from the beginning after the last set.
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
!    Input, integer ( kind = 4 ) N, the size of the set from which subsets
!    are drawn.  N must be positive.
!
!    Input, integer ( kind = 4 ) K, the size of the desired subset.  K must be
!    between 0 and N.
!
!    Input/output, integer ( kind = 4 ) A(K).  On input, the user must
!    supply a subset of size K in A.  That is, A must
!    contain K unique numbers, in order, between 1 and N.  On
!    output, A(I) is the I-th element of the output subset.
!    The output array is also in sorted order.
!
!    Output, integer ( kind = 4 ) IN, the element of the output subset which
!    was not in the input set.  Each new subset differs from the
!    last one by adding one element and deleting another.
!
!    Output, integer ( kind = 4 ) IOUT, the element of the input subset which
!    is not in the output subset.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) in
  integer ( kind = 4 ) iout
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_NEXT2 - Fatal error!'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a)' ) '  but 0 < N is required!'
    stop
  end if

  if ( k < 0 .or. n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_NEXT2 - Fatal error!'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K <= N is required!'
    stop
  end if

  j = 0

  do

    if ( 0 < j .or. mod ( k, 2 ) == 0 ) then

      j = j + 1

      if ( k < j ) then
        a(k) = k
        in = k
        iout = n
        return
      end if

      if ( a(j) /= j ) then

        iout = a(j)
        in = iout - 1
        a(j) = in

        if ( j /= 1 ) then
          in = j - 1
          a(j-1) = in
        end if

        return

      end if

    end if

    j = j + 1
    m = n

    if ( j < k ) then
      m = a(j+1) - 1
    end if

    if ( m /= a(j) ) then
      exit
    end if

  end do

  in = a(j) + 1
  a(j) = in
  iout = in - 1

  if ( j /= 1 ) then
    a(j-1) = iout
    iout = j - 1
  end if

  return
end
subroutine ksub_next3 ( n, k, a, more, in, iout )

!*****************************************************************************80
!
!! KSUB_NEXT3 generates the subsets of size K from a set of size N.
!
!  Discussion:
!
!    The routine uses the revolving door method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 2001
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
!    Input, integer ( kind = 4 ) N, the size of the set from which subsets
!    are drawn.  N must be positive.
!
!    Input, integer ( kind = 4 ) K, the size of the desired subsets.  K must be
!    between 0 and N.
!
!    Input/output, integer ( kind = 4 ) A(K).  A(I) is the I-th element of the
!    output subset.  The elements of A are sorted.
!
!    Input/output, logical MORE.  On first call, set MORE = FALSE
!    to signal the beginning.  MORE will be set to TRUE, and on
!    each call, the routine will return another K-subset.
!    Finally, when the last subset has been returned,
!    MORE will be set FALSE and you may stop calling.
!
!    Output, integer ( kind = 4 ) IN, the element of the output subset which
!    was not in the input set.  Each new subset differs from the
!    last one by adding one element and deleting another.  IN is not
!    defined the first time that the routine returns, and is
!    set to zero.
!
!    Output, integer ( kind = 4 ) IOUT, the element of the input subset which is
!    not in the output subset.  IOUT is not defined the first time
!    the routine returns, and is set to zero.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) in
  integer ( kind = 4 ) iout
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) n

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_NEXT3 - Fatal error!'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a)' ) '  but 0 < N is required!'
    stop
  end if

  if ( k < 0 .or. n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_NEXT3 - Fatal error!'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K <= N is required!'
    stop
  end if

  if ( .not. more ) then
    in = 0
    iout = 0
    call i4vec_indicator ( k, a )
    more = ( k /= n )
    return
  end if

  j = 0

  do

    if ( 0 < j .or. mod ( k, 2 ) == 0 ) then

      j = j + 1

      if ( a(j) /= j ) then

        iout = a(j)
        in = iout - 1
        a(j) = in

        if ( j /= 1 ) then
          in = j - 1
          a(j-1) = in
        end if

        if ( k /= 1 ) then
          more = ( a(k-1) == k-1 )
        end if

        more = ( .not. more ) .or. ( a(k) /= n )

        return

      end if

    end if

    j = j + 1
    m = n

    if ( j < k ) then
      m = a(j+1) - 1
    end if

    if ( m /= a(j) ) then
      exit
    end if

  end do

  in = a(j) + 1
  a(j) = in
  iout = in - 1

  if ( j /= 1 ) then
    a(j-1) = iout
    iout = j - 1
  end if

  if ( k /= 1 ) then
    more = ( a(k-1) == k-1 )
  end if

  more = ( .not. more ) .or. ( a(k) /= n )

  return
end
subroutine ksub_next4 ( n, k, a, done )

!*****************************************************************************80
!
!! KSUB_NEXT4 generates the subsets of size K from a set of size N.
!
!  Discussion:
!
!    The subsets are generated one at a time.
!
!    The routine should be used by setting DONE to TRUE, and then calling
!    repeatedly.  Each call returns with DONE equal to FALSE, the array
!    A contains information defining a new subset.  When DONE returns
!    equal to TRUE, there are no more subsets.
!
!    There are ( N*(N-1)*...*(N+K-1)) / ( K*(K-1)*...*2*1) such subsets.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the entire set.
!
!    Input, integer ( kind = 4 ) K, the size of the desired subset.  K must be
!    between 0 and N.
!
!    Input/output, integer ( kind = 4 ) A(K), contains information about
!    the subsets.  On the first call with DONE = TRUE, the input contents
!    of A don't matter.  Thereafter, the input value of A
!    should be the same as the output value of the previous call.
!    In other words, leave the array alone!
!    On output, as long as DONE is returned FALSE, A contains
!    information defining a subset of K elements of a set of N elements.
!    In other words, A will contain K distinct numbers (in order)
!    between 1 and N.
!
!    Input/output, logical DONE.
!    On the first call, DONE is an input quantity with a value
!    of TRUE which tells the program to initialize data and
!    return the first subset.
!    On return, DONE is an output quantity that is TRUE as long as
!    the routine is returning another subset, and FALSE when
!    there are no more.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  logical done
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jsave
  integer ( kind = 4 ) n

  if ( k < 0 .or. n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_NEXT4 - Fatal error!'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K <= N is required!'
    stop
  end if
!
!  First call:
!
  if ( done ) then

    call i4vec_indicator ( k, a )

    if ( 0 < n ) then
      done = .false.
    else
      done = .true.
    end if
!
!  Next call.
!
  else

    if ( a(1) < n-k+1 ) then

      done = .false.

      jsave = k

      do j = 1, k-1

        if ( a(j) + 1 < a(j+1) ) then
          jsave = j
          exit
        end if

      end do

      call i4vec_indicator ( jsave - 1, a )
      a(jsave) = a(jsave) + 1

    else

      done = .true.

    end if

  end if

  return
end
subroutine ksub_random ( n, k, seed, a )

!*****************************************************************************80
!
!! KSUB_RANDOM selects a random subset of size K from a set of size N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2003
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
!    Input, integer ( kind = 4 ) N, the size of the set from which subsets
!    are drawn.
!
!    Input, integer ( kind = 4 ) K, number of elements in desired subsets.
!    K must be between 0 and N.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, integer ( kind = 4 ) A(K).  A(I) is the I-th element of the
!    output set.  The elements of A are in order.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) ids
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) is
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed

  if ( k < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K is required!'
    stop
  else if ( n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  K <= N is required!'
    stop
  end if

  if ( k == 0 ) then
    return
  end if

  do i = 1, k
    a(i) = ( ( i - 1 ) * n ) / k
  end do

  do i = 1, k

    do

      ix = i4_uniform ( 1, n, seed )

      l = 1 + ( ix * k - 1 ) / n

      if ( a(l) < ix ) then
        exit
      end if

    end do

    a(l) = a(l) + 1

  end do

  ip = 0
  is = k

  do i = 1, k

    m = a(i)
    a(i) = 0

    if ( m /= ( ( i - 1 ) * n ) / k ) then
      ip = ip + 1
      a(ip) = m
    end if

  end do

  ihi = ip

  do i = 1, ihi
    ip = ihi + 1 - i
    l = 1 + ( a(ip) * k - 1 ) / n
    ids = a(ip) - ( ( l - 1 ) * n ) / k
    a(ip) = 0
    a(is) = l
    is = is - ids
  end do

  do ll = 1, k

    l = k + 1 - ll

    if ( a(l) /= 0 ) then
      ir = l
      m0 = 1 + ( ( a(l) - 1 ) * n ) / k
      m = ( a(l) * n ) / k - m0 + 1
    end if

    ix = i4_uniform ( m0, m0 + m - 1, seed )

    i = l + 1

    do while ( i <= ir )

      if ( ix < a(i) ) then
        exit
      end if

      ix = ix + 1
      a(i-1) = a(i)
      i = i + 1

    end do

    a(i-1) = ix
    m = m - 1

  end do

  return
end
subroutine ksub_random2 ( n, k, seed, a )

!*****************************************************************************80
!
!! KSUB_RANDOM2 selects a random subset of size K from a set of size N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2003
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
!    Input, integer ( kind = 4 ) N, the size of the set.
!
!    Input, integer ( kind = 4 ) K, the size of the subset, between 0 and N.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) A(K), the indices of the selected elements.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) available
  integer ( kind = 4 ) candidate
  integer ( kind = 4 ) have
  integer ( kind = 4 ) n
  integer ( kind = 4 ) need
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( k < 0 .or. n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM2 - Fatal error!'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K <= N is required!'
    stop
  end if

  if ( k == 0 ) then
    return
  end if

  need = k
  have = 0

  available = n
  candidate = 0

  do

    candidate = candidate + 1

    r = r8_uniform_01 ( seed )

    if ( real ( available, kind = 8 ) * r <= real ( need, kind = 8 ) ) then

      need = need - 1
      have = have + 1
      a(have) = candidate

      if ( need <= 0 ) then
        exit
      end if

    end if

    available = available - 1

  end do

  return
end
subroutine ksub_random3 ( n, k, seed, a )

!*****************************************************************************80
!
!! KSUB_RANDOM3 selects a random subset of size K from a set of size N.
!
!  Discussion:
!
!    This routine uses Floyd's algorithm.
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
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the set from which subsets
!    are drawn.
!
!    Input, integer ( kind = 4 ) K, number of elements in desired subsets.
!    K must be between 0 and N.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) A(N).  I is an element of the subset
!    if A(I) = 1, and I is not an element if A(I)=0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed

  if ( k < 0 .or. n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM3 - Fatal error!'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K <= N is required!'
    stop
  end if

  a(1:n) = 0

  if ( k == 0 ) then
    return
  end if

  do i = n - k + 1, n

    j = i4_uniform ( 1, i, seed )

    if ( a(j) == 0 ) then
      a(j) = 1
    else
      a(i) = 1
    end if

  end do

  return
end
subroutine ksub_random4 ( n, k, seed, a )

!*****************************************************************************80
!
!! KSUB_RANDOM4 selects a random subset of size K from a set of size N.
!
!  Discussion:
!
!    This routine is somewhat impractical for the given problem, but
!    it is included for comparison, because it is an interesting
!    approach that is superior for certain applications.
!
!    The approach is mainly interesting because it is "incremental";
!    it proceeds by considering every element of the set, and does not
!    need to know how many elements there are.
!
!    This makes this approach ideal for certain cases, such as the
!    need to pick 5 lines at random from a text file of unknown length,
!    or to choose 6 people who call a certain telephone number on a
!    given day.  Using this technique, it is possible to make the
!    selection so that, whenever the input stops, a valid uniformly
!    random subset has been chosen.
!
!    Obviously, if the number of items is known in advance, and
!    it is easy to extract K items directly, there is no need for
!    this approach, and it is less efficient since, among other costs,
!    it has to generate a random number for each item, and make an
!    acceptance/rejection test.
!
!    This routine is based on "8.6: Picking a Random Line from a File",
!    in the Perl Cookbook.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Tom Christiansen, Nathan Torkington,
!    Perl Cookbook,
!    OReilly, 1999.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the set from which subsets
!    are drawn.
!
!    Input, integer ( kind = 4 ) K, number of elements in desired subsets.
!    K must be between 0 and N.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) A(K), contains the indices of the
!    selected items.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  integer ( kind = 4 ) next
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  next = 0
!
!  Here, we use a DO WHILE to suggest that the algorithm
!  proceeds to the next item, without knowing how many items
!  there are in total.
!
!  Note that this is really the only place where N occurs,
!  so other termination criteria could be used, and we really
!  don't need to know the value of N!
!
  do while ( next < n )

    next = next + 1

    if ( next <= k ) then

      i = next
      a(i) = next

    else

      r = r8_uniform_01 ( seed )

      if ( r * real ( next, kind = 8 ) <= real ( k, kind = 8 ) ) then
        i = i4_uniform ( 1, k, seed )
        a(i) = next
      end if

    end if

  end do

  return
end
subroutine ksub_random5 ( n, k, seed, a )

!*****************************************************************************80
!
!! KSUB_RANDOM5 selects a random subset of size K from a set of size N.
!
!  Discussion:
!
!    Consider the set A(1:N) = 1, 2, 3, ... N.  
!    Choose a random index I1 between 1 and N, and swap items A(1) and A(I1).
!    Choose a random index I2 between 2 and N, and swap items A(2) and A(I2).
!    repeat K times.
!    A(1:K) is your random K-subset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the set from which subsets
!    are drawn.
!
!    Input, integer ( kind = 4 ) K, number of elements in desired subsets.
!    1 <= K <= N.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) A(K), the indices of the randomly
!    chosen elements.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) t
!
!  Let B index the set.
!
  do i = 1, n
    b(i) = i
  end do
!
!  Choose item 1 from N things,
!  choose item 2 from N-1 things,
!  choose item K from N-K+1 things.
!
  do i = 1, k

    j = i4_uniform ( i, n, seed )

    t    = b(i)
    b(i) = b(j)
    b(j) = t

  end do
!
!  Copy the first K elements.
!
  a(1:k) = b(1:k)
!
!  Put the elements in ascending order.
!
  call i4vec_sort_heap_a ( k, a )

  return
end
subroutine ksub_rank ( k, a, rank )

!*****************************************************************************80
!
!! KSUB_RANK computes the rank of a K subset of an N set.
!
!  Discussion:
!
!    The routine accepts an array representing a subset of size K from a set
!    of size N, and returns the rank (or order) of that subset.
!
!    It uses the same ranking that KSUB_NEXT2 uses to generate all the subsets
!    one at a time.
!
!    Note the value of N is not input, and is not, in fact,
!    needed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2003
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
!    Input, integer ( kind = 4 ) K, the number of elements in the subset.
!
!    Input, integer ( kind = 4 ) A(K), contains K distinct numbers between
!    1 and N, in order.
!
!    Output, integer ( kind = 4 ) RANK, the rank of this subset.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iprod
  integer ( kind = 4 ) j
  integer ( kind = 4 ) rank

  rank = 0

  do i = 1, k

    iprod = 1

    do j = i+1, a(i)-1
      iprod = iprod * j
    end do

    do j = 1, a(i)-i-1
      iprod = iprod / j
    end do

    if ( a(i) == 1 ) then
      iprod = 0
    end if

    rank = rank + iprod

  end do

  rank = rank + 1

  return
end
subroutine ksub_unrank ( k, rank, a )

!*****************************************************************************80
!
!! KSUB_UNRANK returns the subset of a given rank.
!
!  Discussion:
!
!    The routine is given a rank and returns the corresponding subset of K
!    elements of a set of N elements.
!
!    It uses the same ranking that KSUB_NEXT2 uses to generate all the subsets
!    one at a time.
!
!    Note that the value of N itself is not input, nor is it needed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2004
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
!    Input, integer ( kind = 4 ) K, the number of elements in the subset.
!
!    Input, integer ( kind = 4 ) RANK, the rank of the desired subset.
!    There are ( N*(N-1)*...*(N+K-1)) / ( K*(K-1)*...*2*1) such
!    subsets, so RANK must be between 1 and that value.
!
!    Output, integer ( kind = 4 ) A(K), K distinct integers in order between
!    1 and N, which define the subset.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iprod
  integer ( kind = 4 ) jrank
  integer ( kind = 4 ) rank

  jrank = rank - 1

  do i = k, 1, -1

    ip = i - 1
    iprod = 1

    do

      ip = ip + 1

      if ( ip /= i ) then
        iprod = ( ip * iprod ) / ( ip - i )
      end if

      if ( jrank < iprod ) then
        exit
      end if

    end do

    if ( ip /= i ) then
      iprod = ( ( ip - i ) * iprod ) / ip
    end if

    jrank = jrank - iprod
    a(i) = ip

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
!    In the following discussion, we will let '0' stand for FALSE and
!    '1' for TRUE.
!
!    The logical vectors have the order
!
!      (0,0,...,0),
!      (0,0,...,1),
!      ...
!      (1,1,...,1)
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
!    31 May 2008
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
  logical lvec(n)

  do i = n, 1, -1

    if ( .not. lvec(i) ) then
      lvec(i) = .true.
      return
    end if

    lvec(i) = .false.

  end do

  return
end
subroutine matrix_product_opt ( n, rank, cost, order )

!*****************************************************************************80
!
!! MATRIX_PRODUCT_OPT determines the optimal cost of a matrix product.
!
!  Discussion:
!
!    The cost of multiplying an LxM matrix by an M by N matrix is
!    assessed as L*M*N.
!
!    Any particular order of multiplying a set of N matrices is equivalent
!    to parenthesizing an expression of N objects.
!
!    The actual number of ways of parenthesizing an expression
!    of N objects is C(N), the N-th Catalan number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Sedgewick,
!    Algorithms in C,
!    Addison-Wesley, 1990,
!    ISBN: 0-201-51425-7,
!    LC: QA76.73.C15S43.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of matrices to be multiplied.
!
!    Input, integer ( kind = 4 ) RANK(N+1), the rank information for the
!    matrices.  Matrix I has RANK(I) rows and RANK(I+1) columns.
!
!    Output, integer ( kind = 4 ) COST, the cost of the multiplication if the
!    optimal order is used.
!
!    Output, integer ( kind = 4 ) ORDER(N-1), indicates the order in which the
!    N-1 multiplications are to be carried out.  ORDER(1) is the first
!    multiplication to do, and so on.
!
  implicit none

  integer ( kind = 4 ), parameter :: stack_max = 100
  integer ( kind = 4 ) n

  integer ( kind = 4 ) best(n,n)
  integer ( kind = 4 ) cost
  integer ( kind = 4 ) cost2(n,n)
  integer ( kind = 4 ) cost3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) order(n-1)
  integer ( kind = 4 ) rank(n+1)
  integer ( kind = 4 ) stack(stack_max)
  integer ( kind = 4 ) stack_num
  integer ( kind = 4 ) step
!
!  Initialize the cost matrix.
!
  do i = 1, n

    cost2(i,1:i) = 0
    cost2(i,i+1:n) = i4_huge ( )

  end do
!
!  Initialize the BEST matrix.
!
  best(1:n,1:n) = 0
!
!  Compute the cost and best matrices.
!
  do j = 1, n-1
    do i = 1, n-j
      do k = i+1, i+j
        cost3 = cost2(i,k-1) + cost2(k,i+j) + rank(i) * rank(k) * rank(i+j+1)
        if ( cost3 < cost2(i,i+j) ) then
          cost2(i,i+j) = cost3
          best(i,i+j) = k
        end if
      end do
    end do
  end do
!
!  Pick off the optimal cost.
!
  cost = cost2(1,n)
!
!  Backtrack to determine the optimal order.
!
  stack_num = 0

  i1 = 1
  i2 = n

  if ( i1 + 1 < i2 ) then
    stack_num = stack_num + 1
    stack(stack_num) = i1
    stack_num = stack_num + 1
    stack(stack_num) = i2
  end if

  step = n - 1
!
!  Take an item off the stack.
!
  do while ( 0 < stack_num )

    i3 = stack(stack_num)
    stack_num = stack_num - 1
    i1 = stack(stack_num)
    stack_num = stack_num - 1

    i2 = best(i1,i3)

    order(step) = i2 - 1
    step = step - 1
!
!  The left chunk is matrices (I1...I2-1)
!
    if ( i1 == i2 - 1 ) then

    else if ( i1 + 1 == i2 - 1 ) then
      order(step) = i2 - 2
      step = step - 1
    else
      stack_num = stack_num + 1
      stack(stack_num) = i1
      stack_num = stack_num + 1
      stack(stack_num) = i2 - 1
    end if
!
!  The right chunk is matrices (I2...I3)
!
    if ( i2 == i3 ) then

    else if ( i2 + 1 == i3 ) then
      order(step) = i2
      step = step - 1
    else
      stack_num = stack_num + 1
      stack(stack_num) = i2
      stack_num = stack_num + 1
      stack(stack_num) = i3
    end if

  end do

  return
end
subroutine moebius_matrix ( n, a, mu )

!*****************************************************************************80
!
!! MOEBIUS_MATRIX finds the Moebius matrix from a covering relation.
!
!  Discussion:
!
!    This routine can be called with A and MU being the same matrix.
!    The routine will correctly compute the Moebius matrix, which
!    will, in this case, overwrite the input matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 2004
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
!    Input, integer ( kind = 4 ) N, number of elements in the partially ordered
!    set.
!
!    Input, integer ( kind = 4 ) A(N,N).  A(I,J) = 1 if I is covered by J,
!    0 otherwise.
!
!    Output, integer ( kind = 4 ) MU(N,N), the Moebius matrix as computed
!    by the routine.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mu(n,n)
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) q(n)
!
!  Compute a reordering P of the elements of the partially ordered matrix.
!
  call triang ( n, a, p )
!
!  Copy the matrix.
!
  mu(1:n,1:n) = a(1:n,1:n)
!
!  Apply the reordering to MU.
!
  call i4mat_perm2 ( n, n, mu, p, p )
!
!  Negate the (strict) upper triangular elements of MU.
!
  do i = 1, n-1
    mu(i,i+1:n) = -mu(i,i+1:n)
  end do
!
!  Compute the inverse of MU.
!
  call i4mat_u1_inverse ( n, mu, mu )
!
!  All nonzero elements are reset to 1.
!
  do i = 1, n
    do j = i, n
      if ( mu(i,j) /= 0 ) then
        mu(i,j) = 1
      end if
    end do
  end do
!
!  Invert the matrix again.
!
  call i4mat_u1_inverse ( n, mu, mu )
!
!  Compute the inverse permutation.
!
  do i = 1, n
    q(p(i)) = i
  end do
!
!  Unpermute the rows and columns of MU.
!
  call i4mat_perm2 ( n, n, mu, q, q )

  return
end
subroutine monomial_count ( degree_max, dim, total )

!*****************************************************************************80
!
!! MONOMIAL_COUNT counts the number of monomials up to a given degree.
!
!  Discussion:
!
!    In 3D, there are 10 monomials of degree 3 or less:
!
!    Degree  Count  List
!    ------  -----  ----
!         0      1  1
!         1      3  x y z
!         2      6  xx xy xz yy yz zz
!         3     10  xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
!
!    Total      20
!
!    The formula is
!
!      COUNTS(DEGREE,DIM) = (DIM-1+DEGREE)! / (DIM-1)! / DEGREE!
!
!      TOTAL              = (DIM  +DEGREE)! / (DIM)!   / DEGREE!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DEGREE_MAX, the maximum degree.
!
!    Input, integer ( kind = 4 ) DIM, the spatial dimension.
!
!    Output, integer ( kind = 4 ) TOTAL, the total number of monomials
!    of degrees 0 through DEGREE_MAX.
!
  implicit none

  integer ( kind = 4 ) bot
  integer ( kind = 4 ) degree_max
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) top
  integer ( kind = 4 ) total

  total = 1

  if ( degree_max < dim ) then

    top = dim + 1
    do bot = 1, degree_max
      total = ( total * top ) / bot
      top = top + 1
    end do

  else

    top = degree_max + 1
    do bot = 1, dim
      total = ( total * top ) / bot
      top = top + 1
    end do

  end if

  return
end
subroutine monomial_counts ( degree_max, dim, counts )

!*****************************************************************************80
!
!! MONOMIAL_COUNTS counts the number of monomials up to a given degree.
!
!  Discussion:
!
!    In 3D, there are 10 monomials of degree 3 or less:
!
!    Degree  Count  List
!    ------  -----  ----
!         0      1  1
!         1      3  x y z
!         2      6  xx xy xz yy yz zz
!         3     10  xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
!
!    Total      20
!
!    The formula is
!
!      COUNTS(DEGREE,DIM) = (DIM-1+DEGREE)! / (DIM-1)! / DEGREE!
!
!      TOTAL              = (DIM  +DEGREE)! / (DIM)!   / DEGREE!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DEGREE_MAX, the maximum degree.
!
!    Input, integer ( kind = 4 ) DIM, the spatial dimension.
!
!    Output, integer ( kind = 4 ) COUNTS(0:DEGREE_MAX), the number of
!    monomials of each degree.
!
  implicit none

  integer ( kind = 4 ) degree_max

  integer ( kind = 4 ) counts(0:degree_max)
  integer ( kind = 4 ) degree
  integer ( kind = 4 ) dim

  degree = 0
  counts(degree) = 1

  do degree = 1, degree_max
    counts(degree) = ( counts(degree-1) * ( dim - 1 + degree ) ) / degree
  end do

  return
end
subroutine morse_thue ( i, s )

!*****************************************************************************80
!
!! MORSE_THUE generates a Morse_Thue number.
!
!  Discussion:
!
!    The Morse_Thue sequence can be defined in a number of ways.
!
!    A) Start with the string containing the single letter '0'; then
!       repeatedly apply the replacement rules '0' -> '01' and
!       '1' -> '10' to the letters of the string.  The Morse_Thue sequence
!       is the resulting letter sequence.
!
!    B) Starting with the string containing the single letter '0',
!       repeatedly append the binary complement of the string to itself.
!       Thus, '0' becomes '0' + '1' or '01', then '01' becomes
!       '01' + '10', which becomes '0110' + '1001', and so on.
!
!    C) Starting with I = 0, the I-th Morse-Thue number is determined
!       by taking the binary representation of I, adding the digits,
!       and computing the remainder modulo 2.
!
!  Example:
!
!     I  binary   S
!    --  ------  --
!     0       0   0
!     1       1   1
!     2      10   1
!     3      11   0
!     4     100   1
!     5     101   0
!     6     110   0
!     7     111   1
!     8    1000   1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the Morse-Thue number.
!    Normally, I is 0 or greater, but any value is allowed.
!
!    Output, integer ( kind = 4 ) S, the Morse-Thue number of index I.
!
  implicit none

  integer ( kind = 4 ), parameter :: nbits = 32

  integer ( kind = 4 ) b(nbits)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_copy
  integer ( kind = 4 ) s

  i_copy = abs ( i )
!
!  Expand I into binary form.
!
  call ui4_to_ubvec ( i_copy, nbits, b )
!
!  Sum the 1's in the binary representation.
!
  s = sum ( b(1:nbits) )
!
!  Take the value modulo 2.
!
  s = mod ( s, 2 )

  return
end
subroutine multinomial_coef1 ( nfactor, factor, ncomb )

!*****************************************************************************80
!
!! MULTINOMIAL_COEF1 computes a multinomial coefficient.
!
!  Discussion:
!
!    The multinomial coefficient is a generalization of the binomial
!    coefficient.  It may be interpreted as the number of combinations of
!    N objects, where FACTOR(1) objects are indistinguishable of type 1,
!    ... and FACTOR(NFACTOR) are indistinguishable of type NFACTOR,
!    and N is the sum of FACTOR(1) through FACTOR(NFACTOR).
!
!    NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
!
!    The logarithm of the Gamma function is used, to avoid overflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NFACTOR, the number of factors.
!
!    Input, integer ( kind = 4 ) FACTOR(NFACTOR), contains the factors.
!    0 <= FACTOR(I)
!
!    Output, integer ( kind = 4 ) NCOMB, the value of the multinomial
!    coefficient.
!
  implicit none

  integer ( kind = 4 ) nfactor

  real ( kind = 8 ) arg
  real ( kind = 8 ) fack
  real ( kind = 8 ) facn
  integer ( kind = 4 ) factor(nfactor)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncomb
  real ( kind = 8 ) r8_gamma_log
!
!  Each factor must be nonnegative.
!
  do i = 1, nfactor

    if ( factor(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MULTINOMIAL_COEF1 - Fatal error!'
      write ( *, '(a,i8,a,i8)' ) '  Factor ', i, ' = ', factor(i)
      write ( *, '(a)' ) '  But this value must be nonnegative.'
      stop
    end if

  end do
!
!  The factors sum to N.
!
  n = sum ( factor(1:nfactor) )

  arg = real ( n + 1, kind = 8 )
  facn = r8_gamma_log ( arg )

  do i = 1, nfactor

    arg = real ( factor(i) + 1, kind = 8 )
    fack = r8_gamma_log ( arg )
    facn = facn - fack

  end do

  ncomb = nint ( exp ( facn ) )

  return
end
subroutine multinomial_coef2 ( nfactor, factor, ncomb )

!*****************************************************************************80
!
!! MULTINOMIAL_COEF2 computes a multinomial coefficient.
!
!  Discussion:
!
!    The multinomial coefficient is a generalization of the binomial
!    coefficient.  It may be interpreted as the number of combinations of
!    N objects, where FACTOR(1) objects are indistinguishable of type 1,
!    ... and FACTOR(NFACTOR) are indistinguishable of type NFACTOR,
!    and N is the sum of FACTOR(1) through FACTOR(NFACTOR).
!
!    NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
!
!    A direct method is used, which should be exact.  However, there
!    is a possibility of intermediate overflow of the result.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NFACTOR, the number of factors.
!
!    Input, integer ( kind = 4 ) FACTOR(NFACTOR), contains the factors.
!    0 <= FACTOR(I)
!
!    Output, integer ( kind = 4 ) NCOMB, the multinomial coefficient.
!
  implicit none

  integer ( kind = 4 ) nfactor

  integer ( kind = 4 ) factor(nfactor)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ncomb
!
!  Each factor must be nonnegative.
!
  do i = 1, nfactor

    if ( factor(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MULTINOMIAL_COEF2 - Fatal error!'
      write ( *, '(a,i8,a,i8)' ) '  Factor ', i, ' = ', factor(i)
      write ( *, '(a)' ) '  But this value must be nonnegative.'
      stop
    end if

  end do

  ncomb = 1
  k = 0

  do i = 1, nfactor

    do j = 1, factor(i)
      k = k + 1
      ncomb = ( ncomb * k ) / j
    end do

  end do

  return
end
subroutine multiperm_enum ( n, k, counts, number )

!*****************************************************************************80
!
!! MULTIPERM_ENUM enumerates multipermutations.
!
!  Discussion:
!
!    A multipermutation is a permutation of objects, some of which are
!    identical.
!
!    While there are 6 permutations of the distinct objects A,B,C, there
!    are only 3 multipermutations of the objects A,B,B.
!
!    In general, there are N! permutations of N distinct objects, but
!    there are N! / ( ( M1! ) ( M2! ) ... ( MK! ) ) multipermutations
!    of N objects, in the case where the N objects consist of K
!    types, with M1 examples of type 1, M2 examples of type 2 and so on,
!    and for which objects of the same type are indistinguishable.
!
!  Example:
!
!    Input:
!
!      N = 5, K = 3, COUNTS = (/ 1, 2, 2 /)
!
!    Output:
!
!      Number = 30
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the multipermutation.
!
!    Input, integer ( kind = 4 ) K, the number of types of items.
!    1 <= K.  Ordinarily, K <= N, but we allow any positive K, because
!    we also allow entries in COUNTS to be 0.
!
!    Input, integer ( kind = 4 ) COUNTS(K), the number of items of each type.
!    0 <= COUNTS(1:K) <= N and sum ( COUNTS(1:K) ) = N.
!
!    Output, integer ( kind = 4 ) NUMBER, the number of multipermutations.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) counts(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) number
  integer ( kind = 4 ) top

  if ( n < 0 ) then
    number = -1
    return
  end if

  if ( n == 0 ) then
    number = 1
    return
  end if

  if ( k < 1 ) then
    number = -1
    return
  end if

  if ( any ( counts(1:k) < 0 ) ) then
    number = -1
    return
  end if

  if ( sum ( counts(1:k) ) /= n ) then
    number = -1
    return
  end if
!
!  Ready for computation.
!  By design, the integer division should never have a remainder.
!
  top = 0
  number = 1

  do i = 1, k

    do j = 1, counts ( i )
      top = top + 1
      number = ( number * top ) / j
    end do

  end do

  return
end
subroutine multiperm_next ( n, a, more )

!*****************************************************************************80
!
!! MULTIPERM_NEXT returns the next multipermutation.
!
!  Discussion:
!
!    A multipermutation is a permutation of objects, some of which are
!    identical.
!
!    While there are 6 permutations of the distinct objects A,B,C, there
!    are only 3 multipermutations of the objects A,B,B.
!
!    In general, there are N! permutations of N distinct objects, but
!    there are N! / ( ( M1! ) ( M2! ) ... ( MK! ) ) multipermutations
!    of N objects, in the case where the N objects consist of K
!    types, with M1 examples of type 1, M2 examples of type 2 and so on,
!    and for which objects of the same type are indistinguishable.
!
!    To begin the computation, the user must set up the first multipermutation.
!    To compute ALL possible multipermutations, this first permutation should
!    list the values in ascending order.
!
!    The routine will compute, one by one, the next multipermutation,
!    in lexicographical order.  On the call after computing the last
!    multipermutation, the routine will return MORE = FALSE (and will
!    reset the multipermutation to the FIRST one again.)
!
!  Example:
!
!    1  1 2 2 3 3
!    2  1 2 3 2 3
!    3  1 2 3 3 2
!    4  1 3 2 2 3
!    5  1 3 2 3 2
!    6  1 3 3 2 2
!    7  2 1 2 3 3
!    8  2 1 3 2 3
!    ...
!   30  3 3 2 2 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the multipermutation.
!
!    Input/output, integer ( kind = 4 ) A(N); on input, the current
!    multipermutation.  On output, the next multipermutation.
!
!    Output, logical MORE, is TRUE if the next multipermutation
!    was computed, or FALSE if no further multipermutations could
!    be computed.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) temp
!
!  Step 1:
!  Find M, the last location in A for which A(M) < A(M+1).
!
  m = 0
  do i = 1, n-1
    if ( a(i) < a(i+1) ) then
      m = i
    end if
  end do
!
!  Step 2:
!  If no M was found, we've run out of multipermutations.
!
  if ( m == 0 ) then
    more = .false.
    call i4vec_sort_heap_a ( n, a )
    return
  else
    more = .true.
  end if
!
!  Step 3:
!  Ascending sort A(M+1:N).
!
  if ( m + 1 < n ) then
    call i4vec_sort_heap_a ( n-m, a(m+1:n) )
  end if
!
!  Step 4:
!  Locate the first larger value after A(M).
!
  i = 1
  do while ( a(m+i) <= a(m) )
    i = i + 1
  end do
!
!  Step 5:
!  Interchange A(M) and the next larger value.
!
  temp = a(m)
  a(m) = a(m+i)
  a(m+i) = temp

  return
end
subroutine network_flow_max ( nnode, nedge, iendpt, icpflo, source, sink, &
  cut, node_flow )

!*****************************************************************************80
!
!! NETWORK_FLOW_MAX finds the maximal flow and a minimal cut in a network.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2003
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
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input/output, integer ( kind = 4 ) IENDPT(2,NEDGE), the edges of the
!    network, defined as pairs of nodes.  Each edge should be listed TWICE,
!    the second time in reverse order.  On output, the edges have
!    been reordered, and so the columns of IENDPT have been rearranged.
!
!    Input/output, integer ( kind = 4 ) ICPFLO(2,NEDGE).
!    On input, ICPFLO(1,I) is the capacity of edge I.  On output,
!    ICPFLO(2,I) is the flow on edge I and ICPFLO(1,I) has
!    been rearranged to match the reordering of IENDPT.
!
!    Input, integer ( kind = 4 ) SOURCE, the designated source node.
!
!    Input, integer ( kind = 4 ) SINK, the designated sink node.
!
!    Output, integer ( kind = 4 ) CUT(NNODE).  CUT(I) = 1 if node I is in the
!    minimal cut set, otherwise 0.
!
!    Output, integer ( kind = 4 ) NODE_FLOW(NNODE), the flow through each node.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) cut(nnode)
  integer ( kind = 4 ) del
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarray(nnode)
  integer ( kind = 4 ) icpflo(2,nedge)
  integer ( kind = 4 ) :: ien1 = 0
  integer ( kind = 4 ) :: ien2 = 0
  integer ( kind = 4 ) iendpt(2,nedge)
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iparm
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) iread
  integer ( kind = 4 ) irite
  integer ( kind = 4 ) is
  integer ( kind = 4 ) isort
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kz
  integer ( kind = 4 ) lst
  integer ( kind = 4 ) m
  integer ( kind = 4 ) node_flow(nnode)
  integer ( kind = 4 ) sink
  integer ( kind = 4 ) source
  integer ( kind = 4 ) work1(nnode)
  integer ( kind = 4 ) work2(nnode)
!
!  Initialization.
!
  iarray(1:nnode) = 0
  del = 0
  icpflo(2,1:nedge) = 0

  do i = 1, nedge

    ip = iendpt(1,i)

    if ( ip == source ) then
      del = del + icpflo(1,i)
    end if

    iarray(ip) = iarray(ip) + 1

  end do

  node_flow(source) = del
  is = 1

  do i = 1, nnode
    it = iarray(i)
    iarray(i) = is
    work1(i) = is
    is = is + it
  end do

  isort = 0
!
!  Sorting.
!
10 continue

  indx = 0

50 continue

  do

    call sort_heap_external ( nedge, indx, ien1, ien2, is )

    if ( indx < 0 ) then

      is = iendpt(1,ien1) - iendpt(1,ien2)

      if ( is == 0 ) then
        is = iendpt(2,ien1) - iendpt(2,ien2)
      end if

    else if ( 0 < indx ) then

      do ir = 1, 2
        call i4_swap ( iendpt(ir,ien1), iendpt(ir,ien2) )
        call i4_swap ( icpflo(ir,ien1), icpflo(ir,ien2) )
      end do

    else

      if ( 0 < isort ) then
        return
      end if

      do i = 1, nedge
        iq = iendpt(2,i)
        iendpt(1,i) = work1(iq)
        work1(iq) = work1(iq) + 1
      end do

      go to 100

    end if

  end do

80 continue

  iendpt(1,iendpt(1,ien1)) = ien2
  iendpt(1,iendpt(1,ien2)) = ien1

  do ir = 1, 2
    call i4_swap ( iendpt(ir,ien1), iendpt(ir,ien2) )
    call i4_swap ( icpflo(ir,ien1), icpflo(ir,ien2) )
  end do

  if ( indx < 0 ) then
    work2(iq) = ien2
    go to 280
  end if

  if ( indx == 0 ) then
    go to 170
  end if

  go to 50

100   continue

  indx = 0

  do i = 1, nnode

    if ( i /= source ) then
      node_flow(i) = 0
    end if

    work2(i) = nedge + 1

    if ( i < nnode ) then
      work2(i) = iarray(i+1)
    end if

    cut(i) = 0

  end do

  iread = 0
  irite = 1
  work1(1) = source
  cut(source) = -1

120   continue

  iread = iread + 1

  if ( iread <= irite ) then

    ip = work1(iread)
    lst = work2(ip) - 1
    i = iarray(ip) - 1

    do

      i = i + 1

      if ( lst < i ) then
        go to 120
      end if

      iq = iendpt(2,i)
      del = icpflo(1,i) - icpflo(2,i)

      if ( cut(iq) == 0 .and. del /= 0 ) then

        if ( iq /= sink ) then
          irite = irite + 1
          work1(irite) = iq
        end if

        cut(iq) = -1

      end if

    end do

  end if

  if ( cut(sink) == 0 ) then

    cut(1:nnode) = -cut(1:nnode)

    do i = 1, nedge
      ip = iendpt(2,iendpt(1,i))
      if ( icpflo(2,i) < 0 ) then
        node_flow(ip) = node_flow(ip) - icpflo(2,i)
      end if
      iendpt(1,i) = ip
    end do

    node_flow(source) = node_flow(sink)
    isort = 1
    go to 10

  end if

  cut(sink) = 1

160   continue

  iread = iread - 1

  if ( iread == 0 ) then
    go to 180
  end if

  ip = work1(iread)
  ien1 = iarray(ip) - 1
  ien2 = work2(ip) - 1

170   continue

  if ( ien1 /= ien2 ) then

    iq = iendpt(2,ien2)

    if ( cut(iq) <= 0 .or. icpflo(1,ien2) == icpflo(2,ien2) ) then
      ien2 = ien2 - 1
      go to 170
    end if

    iendpt(2,ien2) = -iq
    icpflo(1,ien2) = icpflo(1,ien2) - icpflo(2,ien2)
    icpflo(2,ien2) = 0
    ien1 = ien1 + 1

    if ( ien1 < ien2 ) then
      go to 80
    end if

  end if

  if ( iarray(ip) <= ien1 ) then
    cut(ip) = ien1
  end if

  go to 160

180   continue

  kz = 0

  do ir = 1, irite
    if ( 0 < cut(work1(ir)) ) then
      kz = kz + 1
      work1(kz) = work1(ir)
    end if
  end do

  indx = -1
  m = 1

200   continue

  ip = work1(m)

  if ( 0 < node_flow(ip) ) then
    go to 250
  end if

210   continue

  m = m + 1

  if ( m <= kz ) then
    go to 200
  end if

  iparm = 0

220   continue

  m = m - 1

  if ( m == 1 ) then

    do i = 1, nedge

      iq = -iendpt(2,i)

      if ( 0 <= iq ) then

        iendpt(2,i) = iq
        j = iendpt(1,i)
        icpflo(1,i) = icpflo(1,i) - icpflo(2,j)

        del = icpflo(2,i) - icpflo(2,j)
        icpflo(2,i) = del
        icpflo(2,j) = -del

      end if

    end do

    go to 100

  end if

  ip = work1(m)

  if ( node_flow(ip) < 0 ) then
    go to 220
  end if

  if ( node_flow(ip) == 0 ) then

    lst = nedge + 1

    if ( ip < nnode ) then
      lst = iarray(ip+1)
    end if

    i = work2(ip)
    work2(ip) = lst

    do

      if ( i == lst ) then
        go to 220
      end if

      j = iendpt(1,i)
      del = icpflo(2,j)
      icpflo(2,j) = 0
      icpflo(1,j) = icpflo(1,j) - del
      icpflo(2,i) = icpflo(2,i) - del
      i = i + 1

    end do

  end if

  if ( cut(ip) < iarray(ip) ) then
    go to 300
  end if

250   continue

  i = cut(ip) + 1

260  continue

  do

    i = i - 1

    if ( i < iarray(ip) ) then
      go to 290
    end if

    iq = -iendpt(2,i)

    if ( 0 <= node_flow(iq) ) then
      exit
    end if

  end do

  del = min ( icpflo(1,i) - icpflo(2,i), node_flow(ip) )
  icpflo(2,i) = icpflo(2,i) + del
  node_flow(ip) = node_flow(ip) - del
  node_flow(iq) = node_flow(iq) + del
  iparm = 1
  ien1 = iendpt(1,i)
  ien2 = work2(iq) - 1

  if ( ien1 < ien2 ) then
    go to 80
  end if

  if ( ien1 == ien2 ) then
    work2(iq) = ien2
  end if

280   continue

  if ( 0 < node_flow(ip) ) then
    go to 260
  end if

  if ( icpflo(1,i) == icpflo(2,i) ) then
    i = i - 1
  end if

290   continue

  cut(ip) = i

  if ( iparm /= 0 ) then
    go to 210
  end if

300   continue

  i = work2(ip)

  do

    j = iendpt(1,i)
    del = min ( icpflo(2,j), node_flow(ip) )
    icpflo(2,j) = icpflo(2,j) - del
    node_flow(ip) = node_flow(ip) - del
    iq = iendpt(2,i)
    node_flow(iq) = node_flow(iq) + del
    i = i + 1

    if ( node_flow(ip) <= 0 ) then
      exit
    end if

  end do

  node_flow(ip) = -1
  go to 220

end
subroutine nim_sum ( i, j, k )

!*****************************************************************************80
!
!! NIM_SUM computes the Nim sum of two integers.
!
!  Discussion:
!
!    If K is the Nim sum of I and J, then each bit of K is the exclusive
!    OR of the corresponding bits of I and J.
!
!  Example:
!
!     I     J     K     I base 2    J base 2    K base 2
!   ----  ----  ----  ----------  ----------  ----------
!      0     0     0           0           0           0
!      1     0     1           1           0           1
!      1     1     0           1           1           0
!      2     7     5          10         111         101
!     11    28    23        1011       11100       10111
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, the integers to be Nim-summed.
!
!    Output, integer ( kind = 4 ) K, the Nim sum of I and J.
!
  implicit none

  integer ( kind = 4 ), parameter :: nbits = 32

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ivec(nbits)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jvec(nbits)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kvec(nbits)

  call ui4_to_ubvec ( i, nbits, ivec )
  call ui4_to_ubvec ( j, nbits, jvec )
  call bvec_xor ( nbits, ivec, jvec, kvec )
  call ubvec_to_ui4 ( nbits, kvec, k )

  return
end
subroutine padovan ( n, p )

!*****************************************************************************80
!
!! PADOVAN returns the first N values of the Padovan sequence.
!
!  Discussion:
!
!    The Padovan sequence has the initial values:
!
!      P(0) = 1
!      P(1) = 1
!      P(2) = 1
!
!    and subsequent entries are generated by the recurrence
!
!      P(I+1) = P(I-1) + P(I-2)
!
!  Example:
!
!    0   1
!    1   1
!    2   1
!    3   2
!    4   2
!    5   3
!    6   4
!    7   5
!    8   7
!    9   9
!   10  12
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
!    Ian Stewart,
!    "A Neglected Number",
!    Scientific American,
!    Volume 274, pages 102-102, June 1996.
!
!    Ian Stewart,
!    Math Hysteria,
!    Oxford, 2004,
!    ISBN: 0198613369.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of terms.
!
!    Output, integer ( kind = 4 ) P(N), terms 0 though N-1 of the sequence.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(n)

  if ( n < 1 ) then
    return
  end if

  p(1) = 1

  if ( n < 2 ) then
    return
  end if

  p(2) = 1

  if ( n < 3 ) then
    return
  end if

  p(3) = 1

  do i = 4, n
    p(i) = p(i-2) + p(i-3)
  end do

  return
end
subroutine pell_basic ( d, x0, y0 )

!*****************************************************************************80
!
!! PELL_BASIC returns the fundamental solution for Pell's basic equation.
!
!  Discussion:
!
!    Pell's equation has the form:
!
!      X*X - D * Y*Y = 1
!
!    where D is a given non-square integer, and X and Y may be assumed
!    to be positive integers.
!
!  Example:
!
!     D   X0   Y0
!
!     2    3    2
!     3    2    1
!     5    9    4
!     6    5    2
!     7    8    3
!     8    3    1
!    10   19    6
!    11   10    3
!    12    7    2
!    13  649  180
!    14   15    4
!    15    4    1
!    17   33    8
!    18   17    4
!    19  170   39
!    20    9    2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2000
!
!  Author:
!
!   John Burkardt
!
!  Reference:
!
!    Mark Herkommer,
!    Number Theory, A Programmer's Guide,
!    McGraw Hill, 1999,
!    ISBN: 0-07-913074-7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) D, the coefficient in Pell's equation.  D
!    should be positive, and not a perfect square.
!
!    Output, integer ( kind = 4 ) X0, Y0, the fundamental or 0'th solution.
!    If X0 = Y0 = 0, then the calculation was canceled because of an error.
!    Both X0 and Y0 will be nonnegative.
!
  implicit none

  integer ( kind = 4 ), parameter :: max_term = 100

  integer ( kind = 4 ) b(0:max_term)
  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n_term
  integer ( kind = 4 ) p
  integer ( kind = 4 ) pm1
  integer ( kind = 4 ) pm2
  integer ( kind = 4 ) q
  integer ( kind = 4 ) qm1
  integer ( kind = 4 ) qm2
  integer ( kind = 4 ) r
  integer ( kind = 4 ) x0
  integer ( kind = 4 ) y0
!
!  If these values are returned, an error has occurred.
!
  x0 = 0
  y0 = 0
!
!  Check D.
!
  if ( d <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PELL_BASIC - Fatal error!'
    write ( *, '(a)' ) '  Pell coefficient D <= 0.'
    stop
  end if

  call i4_sqrt ( d, q, r )

  if ( r == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PELL_BASIC - Fatal error!'
    write ( *, '(a)' ) '  Pell coefficient is a perfect square.'
    stop
  end if
!
!  Find the continued fraction representation of sqrt ( D ).
!
  call i4_sqrt_cf ( d, max_term, n_term, b )
!
!  If necessary, go for two periods.
!
  if ( mod ( n_term, 2 ) == 1 ) then

    do i = n_term + 1, 2*n_term
      b(i) = b(i-n_term)
    end do

    n_term = 2 * n_term

  end if
!
!  Evaluate the continued fraction using the forward recursion algorithm.
!
  pm2 = 0
  pm1 = 1
  qm2 = 1
  qm1 = 0

  do i = 0, n_term-1
    p = b(i) * pm1 + pm2
    q = b(i) * qm1 + qm2
    pm2 = pm1
    pm1 = p
    qm2 = qm1
    qm1 = q
  end do
!
!  Get the fundamental solution.
!
  x0 = p
  y0 = q

  return
end
subroutine pell_next ( d, x0, y0, xn, yn, xnp1, ynp1 )

!*****************************************************************************80
!
!! PELL_NEXT returns the next solution of Pell's equation.
!
!  Discussion:
!
!    Pell's equation has the form:
!
!      X*X - D * Y*Y = 1
!
!    where D is a given non-square integer, and X and Y may be assumed
!    to be positive integers.
!
!    To compute X0, Y0, call PELL_BASIC.
!    To compute X1, Y1, call this routine, with XN and YN set to X0 and Y0.
!    To compute further solutions, call again with X0, Y0 and the previous
!    solution.
!
!  Example:
!
!    ------INPUT--------  --OUTPUT--
!
!    D  X0  Y0   XN   YN  XNP1  YNP1
!
!    2   3   2    3    2    17    12
!    2   3   2   17   12    99    70
!    2   3   2   99   70   577   408
!    2   3   2  577  408  3363  2378
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 May 2000
!
!  Author:
!
!   John Burkardt
!
!  Reference:
!
!    Mark Herkommer,
!    Number Theory, A Programmer's Guide,
!    McGraw Hill, 1999,
!    ISBN: 0-07-913074-7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) D, the coefficient in Pell's equation.
!
!    Input, integer ( kind = 4 ) X0, Y0, the fundamental or 0'th solution.
!
!    Input, integer ( kind = 4 ) XN, YN, the N-th solution.
!
!    Output, integer ( kind = 4 ) XNP1, YNP1, the N+1-th solution.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) x0
  integer ( kind = 4 ) xn
  integer ( kind = 4 ) xnp1
  integer ( kind = 4 ) y0
  integer ( kind = 4 ) yn
  integer ( kind = 4 ) ynp1

  xnp1 = x0 * xn + d * y0 * yn
  ynp1 = x0 * yn +     y0 * xn

  return
end
subroutine pent_enum ( n, p )

!*****************************************************************************80
!
!! PENT_ENUM computes the N-th pentagonal number.
!
!  Discussion:
!
!    The pentagonal number P(N) counts the number of dots in a figure of
!    N nested pentagons.  The pentagonal numbers are defined for both
!    positive and negative N.
!
!    The pentagonal numbers are also useful in determining the
!    number of partitions of an integer.
!
!    P(N) = ( N * ( 3 * N - 1 ) ) / 2
!
!  First values:
!
!     N    P
!
!    -5   40
!    -4   26
!    -3   15
!    -2    7
!    -1    2
!     0    0
!     1    1
!     2    5
!     3   12
!     4   22
!     5   35
!     6   51
!     7   70
!     8   92
!     9  117
!    10  145
!
!  Modified:
!
!    22 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the pentagonal number desired.
!
!    Output, integer ( kind = 4 ) P, the value of the N-th pentagonal number.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) p

  p = ( n * ( 3 * n - 1 ) ) / 2

  return
end
subroutine perm_ascend ( n, a, length, sub )

!*****************************************************************************80
!
!! PERM_ASCEND computes the longest ascending subsequence of a permutation.
!
!  Discussion:
!
!    Although this routine is intended to be applied to a permutation,
!    it will work just as well for an arbitrary vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the permutation.
!
!    Input, integer ( kind = 4 ) A(N), the permutation to be examined.
!
!    Output, integer ( kind = 4 ) LENGTH, the length of the longest
!    increasing subsequence.
!
!    Output, integer ( kind = 4 ) SUB(N), contains in entries 1 through LENGTH
!    a longest increasing subsequence of A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) length
  integer ( kind = 4 ) sub(n)
  integer ( kind = 4 ) top(n)
  integer ( kind = 4 ) top_prev(n)

  top(1:n) = 0
  top_prev(1:n) = 0
  sub(1:n) = 0

  if ( n <= 0 ) then
    length = 0
    return
  end if

  length = 0

  do i = 1, n

    k = 0

    do j = 1, length
      if ( a(i) <= a(top(j)) ) then
        k = j
        exit
      end if
    end do

    if ( k == 0 ) then
      length = length + 1
      k = length
    end if

    top(k) = i

    if ( 1 < k ) then
      top_prev(i) = top(k-1)
    else
      top_prev(i) = 0
    end if

  end do

  j = top(length)
  sub(length) = a(j)

  do i = length-1, 1, -1
    j = top_prev(j)
    sub(i) = a(j)
  end do

  return
end
subroutine perm_break_count ( n, p, break_count )

!*****************************************************************************80
!
!! PERM_BREAK_COUNT counts the number of "breaks" in a permutation.
!
!  Discussion:
!
!    We begin with a permutation of order N.  We prepend an element
!    labeled "0" and append an element labeled "N+1".  There are now
!    N+1 pairs of neighbors.  A "break" is a pair of neighbors whose
!    value differs by more than 1.
!
!    The identity permutation has a break count of 0.  The maximum
!    break count is N+1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the permutation.
!
!    Input, integer ( kind = 4 ) P(N), a permutation, in standard index form.
!
!    Output, integer ( kind = 4 ) BREAK_COUNT, the number of breaks in
!    the permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) break_count
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) p(n)

  break_count = 0
!
!  Make sure the permutation is a legal one.
!
  call perm_check ( n, p, ierror )

  if ( p(1) /= 1 ) then
    break_count = break_count + 1
  end if

  do i = 1, n-1
    if ( abs ( p(i+1) - p(i) ) /= 1 ) then
      break_count = break_count + 1
    end if
  end do

  if ( p(n) /= n ) then
    break_count = break_count + 1
  end if

  return
end
subroutine perm_canon_to_cycle ( n, p1, p2 )

!*****************************************************************************80
!
!! PERM_CANON_TO_CYCLE converts a permutation from canonical to cycle form.
!
!  Example:
!
!    Input:
!
!      4 5 2 1 6 3
!
!    Output:
!
!      -4 5 -2 -1 6 3,
!      indicating the cycle structure
!      ( 4, 5 ) ( 2 ) ( 1, 6, 3 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 1, Fundamental Algorithms,
!    Addison Wesley, 1968, page 176.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects permuted.
!
!    Input, integer ( kind = 4 ) P1(N), the permutation, in canonical form.
!
!    Output, integer ( kind = 4 ) P2(N), the permutation, in cycle form.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) p1(n)
  integer ( kind = 4 ) p2(n)
  integer ( kind = 4 ) pmin

  p2(1:n) = p1(1:n)

  pmin = p2(1) + 1

  do i = 1, n

    if ( p2(i) < pmin ) then
      pmin = p2(i)
      p2(i) = -p2(i)
    end if

  end do

  return
end
subroutine perm_check ( n, p, ierror )

!*****************************************************************************80
!
!! PERM_CHECK checks that a vector represents a permutation.
!
!  Discussion:
!
!    The routine verifies that each of the integers from 1
!    to N occurs among the N entries of the permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, integer ( kind = 4 ) P(N), the permutation, in standard index form.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, the array does represent a permutation.
!    nonzero, the array does not represent a permutation.  The smallest
!    missing value is equal to IERROR.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifind
  integer ( kind = 4 ) iseek
  integer ( kind = 4 ) p(n)

  ierror = 0

  do iseek = 1, n

    ierror = iseek

    do ifind = 1, n
      if ( p(ifind) == iseek ) then
        ierror = 0
        exit
      end if
    end do

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
      write ( *, '(a)' ) '  The input array does not represent'
      write ( *, '(a)' ) '  a proper permutation.  In particular, the'
      write ( *, '(a,i8)' ) '  array is missing the value ', ierror
      stop
    end if

  end do

  return
end
subroutine perm_cycle ( n, iopt, p, isgn, ncycle )

!*****************************************************************************80
!
!! PERM_CYCLE analyzes a permutation.
!
!  Discussion:
!
!    The routine will count cycles, find the sign of a permutation,
!    and tag a permutation.
!
!  Example:
!
!    Input:
!
!      N = 9
!      IOPT = 1
!      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
!
!    Output:
!
!      NCYCLE = 3
!      ISGN = +1
!      P = -2, 3, 9, -6, -7, 8, 5, 4, 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
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
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input, integer ( kind = 4 ) IOPT, requests tagging.
!    0, the permutation will not be tagged.
!    1, the permutation will be tagged.
!
!    Input/output, integer ( kind = 4 ) P(N).  On input, P describes a
!    permutation, in the sense that entry I is to be moved to P(I).
!    If IOPT = 0, then P will not be changed by this routine.
!    If IOPT = 1, then on output, P will be "tagged".  That is,
!    one element of every cycle in P will be negated.  In this way,
!    a user can traverse a cycle by starting at any entry I1 of P
!    which is negative, moving to I2 = ABS(P(I1)), then to
!    P(I2), and so on, until returning to I1.
!
!    Output, integer ( kind = 4 ) ISGN, the "sign" of the permutation, which is
!    +1 if the permutation is even, -1 if odd.  Every permutation
!    may be produced by a certain number of pairwise switches.
!    If the number of switches is even, the permutation itself is
!    called even.
!
!    Output, integer ( kind = 4 ) NCYCLE, the number of cycles in the
!    permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iopt
  integer ( kind = 4 ) is
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) ncycle
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, ierror )

  is = 1
  ncycle = n

  do i = 1, n

    i1 = p(i)

    do while ( i < i1 )
      ncycle = ncycle - 1
      i2 = p(i1)
      p(i1) = -i2
      i1 = i2
    end do

    if ( iopt /= 0 ) then
      is = - sign ( 1, p(i) )
    end if

    p(i) = sign ( p(i), is )

  end do

  isgn = 1 - 2 * mod ( n - ncycle, 2 )

  return
end
subroutine perm_cycle_to_canon ( n, p1, p2 )

!*****************************************************************************80
!
!! PERM_CYCLE_TO_CANON converts a permutation from cycle to canonical form.
!
!  Discussion:
!
!    The procedure is to "rotate" the elements of each cycle so that
!    the smallest element is first:
!
!      ( 1, 6, 3 ) ( 4, 5 ) ( 2 )
!
!    and then to sort the cycles in decreasing order of their first
!    (and lowest) element:
!
!      ( 4, 5 ) ( 2 ) ( 1, 6, 3 )
!
!    and then to drop the parentheses:
!
!      4 5 2 1 6 3
!
!  Example:
!
!    Input:
!
!      -6 3 1 -5, 4 -2,
!      indicating the cycle structure
!      ( 6, 3, 1 ) ( 5, 4 ) ( 2 )
!
!    Output:
!
!      4 5 2 1 6 3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 1, Fundamental Algorithms,
!    Addison Wesley, 1968, pages 176.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects permuted.
!
!    Input, integer ( kind = 4 ) P1(N), the permutation, in cycle form.
!
!    Output, integer ( kind = 4 ) P2(N), the permutation, in canonical form.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) hi(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lo(n)
  integer ( kind = 4 ) ncycle
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nhi
  integer ( kind = 4 ) nlo
  integer ( kind = 4 ) nmin
  integer ( kind = 4 ) p1(n)
  integer ( kind = 4 ) p2(n)
  integer ( kind = 4 ) pmin(n)
  integer ( kind = 4 ) ptemp(n)

  p2(1:n) = p1(1:n)
!
!  Work on the next cycle.
!
  nlo = 1
  ncycle = 0

  do while ( nlo <= n )
!
!  Identify NHI, the last index in this cycle.
!
    ncycle = ncycle + 1

    nhi = nlo

    do while ( nhi < n )
      if ( p2(nhi+1) < 0 ) then
        exit
      end if
      nhi = nhi + 1
    end do
!
!  Identify the smallest value in this cycle.
!
    p2(nlo) = -p2(nlo)
    pmin(ncycle) = p2(nlo)
    nmin = nlo

    do i = nlo+1, nhi
      if ( p2(i) < pmin(ncycle) ) then
        pmin(ncycle) = p2(i)
        nmin = i
      end if
    end do
!
!  Rotate the cycle so A_MIN occurs first.
!
    ptemp(nlo+nhi+1-nmin:nhi) = p2(nlo:nmin-1)
    ptemp(nlo:nlo+nhi-nmin) = p2(nmin:nhi)

    lo(ncycle) = nlo
    hi(ncycle) = nhi
!
!  Prepare to operate on the next cycle.
!
    nlo = nhi + 1

  end do
!
!  Compute a sorting index for the cycle minima.
!
  call i4vec_sort_heap_index_d ( ncycle, pmin, indx )
!
!  Copy the cycles out of the temporary array in sorted order.
!
  j = 0
  do i = 1, ncycle
    next = indx(i)
    nlo = lo(next)
    nhi = hi(next)
    do k = nlo, nhi
      j = j + 1
      p2(j) = ptemp(k)
    end do
  end do

  return
end
subroutine perm_cycle_to_index ( n, p1, p2 )

!*****************************************************************************80
!
!! PERM_CYCLE_TO_INDEX converts a permutation from cycle to standard index form.
!
!  Example:
!
!    Input:
!
!      N = 9
!      P1 = -1, 2, 3, 9, -4, 6, 8, -5, 7
!
!    Output:
!
!      P2 = 2, 3, 9, 6, 7, 8, 5, 4, 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
!
!  Author:
!
!    John Burkardt
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
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input, integer ( kind = 4 ) P1(N), the permutation, in cycle form.
!
!    Output, integer ( kind = 4 ) P2(N), the permutation, in standard index
!    form.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
  integer ( kind = 4 ) p1(n)
  integer ( kind = 4 ) p2(n)

  do j = 1, n

    k1 = p1(j)

    if ( k1 < 0 ) then
      k1 = -k1
      k3 = k1
    end if

    if ( j + 1 <= n ) then
      k2 = p1(j+1)
      if ( k2 < 0 ) then
        k2 = k3
      end if
    else
      k2 = k3
    end if

    p2(k1) = k2

  end do

  return
end
subroutine perm_distance ( n, a, b, k )

!*****************************************************************************80
!
!! PERM_DISTANCE computes the Ulam metric distance of two permutations.
!
!  Discussion:
!
!    If we let N be the order of the permutations A and B, and L(P) be
!    the length of the longest ascending subsequence of a permutation P,
!    then the Ulam metric distance between A and B is
!
!      N - L ( A * inverse ( B ) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the permutation.
!
!    Input, integer ( kind = 4 ) A(N), B(N), the permutations to be examined.
!
!    Output, integer ( kind = 4 ) K, the Ulam metric distance between A and B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b(n)
  integer ( kind = 4 ) binv(n)
  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) length
  integer ( kind = 4 ) sub(n)

  binv(1:n) = b(1:n)

  call perm_inverse ( n, binv )

  call perm_mul ( n, a, binv, c )

  call perm_ascend ( n, c, length, sub )

  k = n - length

  return
end
subroutine perm_fixed_enum ( n, m, fnm )

!*****************************************************************************80
!
!! PERM_FIXED_ENUM enumerates the permutations of N objects with M fixed.
!
!  Discussion:
!
!    A permutation of N objects with M fixed is a permutation in which
!    exactly M of the objects retain their original positions.  If
!    M = 0, the permutation is a "derangement".  If M = N, the
!    permutation is the identity.
!
!    F(N,M) = ( N! / M! ) * ( 1 - 1/1! + 1/2! - 1/3! ... 1/(N-M)! )
!           = COMB(N,M) * D(N-M)
!
!    where D(N-M) is the number of derangements of N-M objects.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!    N should be at least 1.
!
!    Input, integer ( kind = 4 ) M, the number of objects that retain their
!    position.  M should be between 0 and N.
!
!    Output, integer ( kind = 4 ) FNM, the number of derangements of N objects
!    in which M objects retain their positions.
!
  implicit none

  integer ( kind = 4 ) derange_enum
  integer ( kind = 4 ) fnm
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  if ( n <= 0 ) then

    fnm = 1

  else if ( m < 0 ) then

    fnm = 0

  else if ( n < m ) then

    fnm = 0

  else if ( m == n ) then

    fnm = 1

  else if ( n == 1 ) then

    if ( m == 1 ) then
      fnm = 1
    else
      fnm = 0
    end if

  else

    fnm = i4_choose ( n, m ) * derange_enum ( n - m )

  end if

  return
end
subroutine perm_free ( npart, ipart, nfree, ifree )

!*****************************************************************************80
!
!! PERM_FREE reports the unused items in a partial permutation.
!
!  Discussion:
!
!    It is assumed that the N objects being permuted are the integers
!    from 1 to N, and that IPART contains a "partial" permutation, that
!    is, the NPART entries of IPART represent the beginning of a
!    permutation of all N items.
!
!    The routine returns in IFREE the items that have not been used yet.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPART, the number of entries in IPART.
!    NPART may be 0.
!
!    Input, integer ( kind = 4 ) IPART(NPART), the partial permutation, which
!    should contain, at most once, some of the integers between 1 and
!    NPART+NFREE.
!
!    Input, integer ( kind = 4 ) NFREE, the number of integers that have not
!    been used in IPART.  This is simply N - NPART.  NFREE may be zero.
!
!    Output, integer ( kind = 4 ) IFREE(NFREE), the integers between 1 and
!    NPART+NFREE that were not used in IPART.
!
  implicit none

  integer ( kind = 4 ) nfree
  integer ( kind = 4 ) npart

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree(nfree)
  integer ( kind = 4 ) ipart(npart)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) match
  integer ( kind = 4 ) n

  n = npart + nfree

  if ( npart < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
    write ( *, '(a)' ) '  NPART < 0.'
    write ( *, '(a,i8)' ) '  NPART = ', npart
    stop

  else if ( npart == 0 ) then

    call i4vec_indicator ( n, ifree )

  else if ( nfree < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
    write ( *, '(a)' ) '  NFREE < 0.'
    write ( *, '(a,i8)' ) '  NFREE = ', nfree
    stop

  else if ( nfree == 0 ) then

    return

  else

    k = 0

    do i = 1, n

      match = 0

      do j = 1, npart
        if ( ipart(j) == i ) then
          match = j
          exit
        end if
      end do

      if ( match == 0 ) then

        k = k + 1

        if ( nfree < k ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
          write ( *, '(a)' ) '  The partial permutation is illegal.'
          write ( *, '(a)' ) '  It should contain, at most once, some of'
          write ( *, '(a,i8)' ) '  the integers between 1 and ', n
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Our error is that NFREE < K,'
          write ( *, '(a)' ) '  We have TOO MANY missing values.'
          write ( *, '(a,i8)' ) '  Value of NFREE = ', nfree
          write ( *, '(a,i8)' ) '  Value of K =     ', k
          call i4vec_print ( npart, ipart, '  Partial permutation:' )
          stop
        end if

        ifree(k) = i

      end if

    end do

  end if

  return
end
subroutine perm_index_to_cycle ( n, p1, p2 )

!*****************************************************************************80
!
!! PERM_INDEX_TO_CYCLE converts a permutation from standard index to cycle form.
!
!  Example:
!
!    Input:
!
!      N = 9
!      P1 = 2, 3, 9, 6, 7, 8, 5, 4, 1
!
!    Output:
!
!      P2 = -1, 2, 3, 9, -4, 6, 8, -5, 7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
!
!  Author:
!
!    John Burkardt
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
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input, integer ( kind = 4 ) P1(N), the permutation, in standard index form.
!
!    Output, integer ( kind = 4 ) P2(N), the permutation, in cycle form.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) p1(n)
  integer ( kind = 4 ) p2(n)

  i = 0
  j = 1

  do while ( j <= n )

    if ( p1(j) < 0 ) then

      j = j + 1

    else

      k = j

      i = i + 1
      p2(i) = -k

      do while ( p1(k) /= j )
        i = i + 1
        p2(i) = p1(k)
        p1(k) = -p1(k)
        k = abs ( p1(k) )
      end do

      p1(k) = -p1(k)

    end if

  end do

  p1(1:n) = abs ( p1(1:n) )

  return
end
subroutine perm_ins ( n, p, ins )

!*****************************************************************************80
!
!! PERM_INS computes the inversion sequence of a permutation.
!
!  Discussion:
!
!    For a given permutation P acting on objects 1 through N, the inversion
!    sequence INS is defined as:
!
!      INS(1) = 0
!      INS(I) = number of values J < I for which P(I) < P(J).
!
!    The original permutation can be recovered from the inversion sequence.
!
!  Example:
!
!    Input:
!
!      ( 3, 5, 1, 4, 2 )
!
!    Output:
!
!      ( 0, 0, 2, 1, 3 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer, 1986,
!    ISBN: 0387963472,
!    LC: QA164.S79.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input, integer ( kind = 4 ) P(N), the permutation, in standard index form.
!    The I-th item has been mapped to P(I).
!
!    Output, integer ( kind = 4 ) INS(N), the inversion sequence of the
!    permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ins(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, ierror )

  ins(1:n) = 0

  do i = 1, n
    do j = 1, i-1
      if ( p(i) < p(j) ) then
        ins(i) = ins(i) + 1
      end if
    end do
  end do

  return
end
subroutine perm_inverse ( n, p )

!*****************************************************************************80
!
!! PERM_INVERSE inverts a permutation "in place".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2000
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) P(N), the permutation, in standard
!    index form.  On output, P describes the inverse permutation
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) is
  integer ( kind = 4 ) p(n)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    stop
  end if

  call perm_check ( n, p, ierror )

  is = 1

  do i = 1, n

    i1 = p(i)

    do while ( i < i1 )
      i2 = p(i1)
      p(i1) = -i2
      i1 = i2
    end do

    is = - sign ( 1, p(i) )
    p(i) = sign ( p(i), is )

  end do

  do i = 1, n

    i1 = -p(i)

    if ( 0 <= i1 ) then

      i0 = i

      do

        i2 = p(i1)
        p(i1) = i0

        if ( i2 < 0 ) then
          exit
        end if

        i0 = i1
        i1 = i2

      end do

    end if

  end do

  return
end
subroutine perm_inverse2 ( n, p )

!*****************************************************************************80
!
!! PERM_INVERSE2 inverts a permutation "in place".
!
!  Discussion:
!
!    The routine needs no extra vector storage in order to compute the
!    inverse of a permutation.
!
!    This feature might be useful if the permutation is large.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects in the permutation.
!
!    Input/output, integer ( kind = 4 ) P(N), the permutation, in standard
!    index form.  On output, the inverse permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, ierror )

  do ii = 1, n

    m = n + 1 - ii
    i = p(m)

    if ( i < 0 ) then

      p(m) = -i

    else if ( i /= m ) then

      k = m

      do

        j = p(i)
        p(i) = -k

        if ( j == m ) then
          p(m) = i
          exit
        end if

        k = i
        i = j

      end do

    end if

  end do

  return
end
subroutine perm_inverse3 ( n, perm, perm_inv )

!*****************************************************************************80
!
!! PERM_INVERSE3 produces the inverse of a given permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items permuted.
!
!    Input, integer ( kind = 4 ) PERM(N), a permutation.
!
!    Output, integer ( kind = 4 ) PERM_INV(N), the inverse permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm_inv(n)

  do i = 1, n
    perm_inv(perm(i)) = i
  end do

  return
end
subroutine perm_lex_next ( n, p, more )

!*****************************************************************************80
!
!! PERM_LEX_NEXT generates permutations in lexical order, one at a time.
!
!  Example:
!
!    N = 3
!
!    1   1 2 3
!    2   1 3 2
!    3   2 1 3
!    4   2 3 1
!    5   3 1 2
!    6   3 2 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 September 1998
!
!  Reference:
!
!    Mok-Kong Shen,
!    Algorithm 202: Generation of Permutations in Lexicographical Order,
!    Communications of the ACM,
!    Volume 6, September 1963, page 517.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements being permuted.
!
!    Input/output, integer ( kind = 4 ) P(N); on first call with MORE = FALSE,
!    this value is not used.  Otherwise, the input value is the previous
!    permutation.  The output value is the next permutation.
!
!    Input/output, logical MORE.
!    On the first call, set MORE = FALSE, to request initialization.
!    On return, if MORE is TRUE, another permutation has been
!    computed and returned, while if MORE is FALSE, there are no more
!    permutations.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical more
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) t
  integer ( kind = 4 ) u
  integer ( kind = 4 ) w
!
!  Initialization.
!
  if ( .not. more ) then

    call i4vec_indicator ( n, p )
    more = .true.

  else

    if ( n <= 1 ) then
      more = .false.
      return
    end if

    w = n

    do while ( p(w) < p(w-1) )

      if ( w == 2 ) then
        more = .false.
        return
      end if

      w = w - 1
    end do

    u = p(w-1)

    do j = n, w, -1

      if ( u < p(j) ) then

        p(w-1) = p(j)
        p(j) = u

        do k = 0, ( n - w - 1 ) / 2
          t      = p(n-k)
          p(n-k) = p(w+k)
          p(w+k) = t
        end do

        return

      end if

    end do

  end if

  return
end
subroutine perm_mul ( n, p1, p2, p3 )

!*****************************************************************************80
!
!! PERM_MUL "multiplies" two permutations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the permutations.
!
!    Input, integer ( kind = 4 ) P1(N), P2(N), the permutations, in standard
!    index form.
!
!    Output, integer ( kind = 4 ) P3(N), the product permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) p1(n)
  integer ( kind = 4 ) p2(n)
  integer ( kind = 4 ) p3(n)

  call perm_check ( n, p1, ierror )

  call perm_check ( n, p2, ierror )

  p3(1:n) = p2(p1(1:n))

  return
end
subroutine perm_next ( n, p, more, even )

!*****************************************************************************80
!
!! PERM_NEXT computes all of the permutations of N objects, one at a time.
!
!  Discussion:
!
!    The routine is initialized by calling with MORE = TRUE, in which case
!    it returns the identity permutation.
!
!    If the routine is called with MORE = FALSE, then the successor of the
!    input permutation is computed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 2001
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
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) P(N), the permutation, in standard
!    index form.  On the first call, the input value is unimportant.
!    On subsequent calls, the input value should be the same
!    as the output value from the previous call.  In other words, the
!    user should just leave P alone.
!    On output, contains the "next" permutation.
!
!    Input/output, logical MORE.
!    Set MORE = FALSE before the first call.
!    MORE will be reset to TRUE and a permutation will be returned.
!    Each new call produces a new permutation until
!    MORE is returned FALSE.
!
!    Input/output, logical EVEN.
!    The input value of EVEN should simply be its output value from the
!    previous call; (the input value on the first call doesn't matter.)
!    On output, EVEN is TRUE if the output permutation is even, that is,
!    involves an even number of transpositions.
!
  implicit none

  integer ( kind = 4 ) n

  logical even
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) id
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) p(n)

  if ( .not. more ) then

    call i4vec_indicator ( n, p )
    more = .true.
    even = .true.

    if ( n == 1 ) then
      more = .false.
      return
    end if

    if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, 2 ) ) then
      return
    end if

    do i = 1, n-3
      if ( p(i+1) /= p(i)+1 ) then
        return
      end if
    end do

    more = .false.

  else

    if ( n == 1 ) then
      p(1) = 0
      more = .false.
      return
    end if

    if ( even ) then

      ia = p(1)
      p(1) = p(2)
      p(2) = ia
      even = .false.

      if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, 2 ) ) then
        return
      end if

      do i = 1, n-3
        if ( p(i+1) /= p(i)+1 ) then
          return
        end if
      end do

      more = .false.
      return

    else

      more = .false.

      is = 0

      do i1 = 2, n

        ia = p(i1)
        i = i1 - 1
        id = 0

        do j = 1, i
          if ( ia < p(j) ) then
            id = id + 1
          end if
        end do

        is = id + is
        if ( id /= i * mod ( is, 2 ) ) then
          more = .true.
          exit
        end if

      end do

      if ( .not. more ) then
        p(1) = 0
        return
      end if

    end if

    m = mod ( is + 1, 2 ) * ( n + 1 )

    do j = 1, i

      if ( sign ( 1, p(j)-ia ) /= sign ( 1, p(j)-m ) ) then
        m = p(j)
        l = j
      end if

    end do

    p(l) = ia
    p(i1) = m
    even = .true.

  end if

  return
end
subroutine perm_next2 ( n, p, done, invers )

!*****************************************************************************80
!
!! PERM_NEXT2 generates all the permutations of N objects.
!
!  Discussion:
!
!    The routine generates the permutations one at a time.  It uses a
!    particular ordering of permutations, generating them from the first
!    (which is the identity permutation) to the N!-th.  The same ordering
!    is used by the routines PERM_RANK and PERM_UNRANK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 July 2004
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer, 1986,
!    ISBN: 0387963472,
!    LC: QA164.S79.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the set to
!    be permuted.
!
!    Input/output, integer ( kind = 4 ) P(N), the permutation, in standard
!    index form.
!
!    Input/output, logical DONE.  The user should set the input value of
!    DONE only once, before the first call to compute the permutations.
!    The user should set DONE to TRUE, which signals the routine
!    that it is to initialize itself.
!    Thereafter, the routine will set DONE to FALSE and will
!    compute a new permutation on each call.
!    However, when there are no more permutations to compute, the
!    routine will not return a new permutation, but instead will
!    return DONE with the value TRUE.  At this point, all the
!    permutations have been computed.
!
!    Output, integer ( kind = 4 ) INVERS(N), the inverse permutation of P.
!
!  Local Parameters:
!
!    Local, integer ACTIVE(N), DIR(N).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), save, allocatable, dimension ( : ) :: active
  integer ( kind = 4 ), save, allocatable, dimension ( : ) :: dir
  logical done
  integer ( kind = 4 ) i
  integer ( kind = 4 ) invers(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nactiv
  integer ( kind = 4 ) p(n)
!
!  An input value of TRUE for DONE is assumed to mean a new
!  computation is beginning.
!
  if ( done ) then

    call i4vec_indicator ( n, p )
    invers(1:n) = p(1:n)

    if ( allocated ( active ) ) then
      deallocate ( active )
    end if

    if ( allocated ( dir ) ) then
      deallocate ( dir )
    end if

    allocate ( active(1:n) )
    allocate ( dir(1:n) )

    dir(1:n) = -1

    active(1) = 0
    active(2:n) = 1
!
!  Set the DONE flag to FALSE, signifying there are more permutations
!  to come.  Except, of course, that we must take care of the trivial case!
!
    if ( 1 < n ) then
      done = .false.
    else
      done = .true.
    end if
!
!  Otherwise, assume we are in a continuing computation
!
  else

    nactiv = 0

    do i = 1, n
      if ( active(i) /= 0 ) then
        nactiv = i
      end if
    end do

    if ( nactiv <= 0 ) then

      done = .true.

    else

      j = invers(nactiv)

      p(j) = p(j+dir(nactiv))
      p(j+dir(nactiv)) = nactiv

      invers(nactiv) = invers(nactiv) + dir(nactiv)
      invers(p(j)) = j

      if ( j + 2 * dir(nactiv) < 1 .or. n < j + 2 * dir(nactiv) ) then
        dir(nactiv) = -dir(nactiv)
        active(nactiv) = 0
      else if ( nactiv < p(j+2*dir(nactiv)) ) then
        dir(nactiv) = -dir(nactiv)
        active(nactiv) = 0
      end if

      active(nactiv+1:n) = 1

    end if

  end if

  if ( done ) then
    deallocate ( active )
    deallocate ( dir )
  end if

  return
end
subroutine perm_next3 ( n, p, more )

!*****************************************************************************80
!
!! PERM_NEXT3 computes all of the permutations of N objects, one at a time.
!
!  Discussion:
!
!    The routine is initialized by calling with MORE = TRUE, in which case
!    it returns the identity permutation.
!
!    If the routine is called with MORE = FALSE, then the successor of the
!    input permutation is computed.
!
!    Trotter's algorithm is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2003
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Hale Trotter,
!    Algorithm 115:
!    PERM,
!    Communications of the Association for Computing Machinery,
!    Volume 5, 1962, pages 434-435.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) P(N), the permutation, in standard
!    index form.  If MORE is TRUE, then P is assumed to contain the
!    "previous" permutation, and on P(I) is the value
!    of the I-th object under the next permutation.
!    Otherwise, P will be set to the "first" permutation.
!
!    Input/output, logical MORE.
!    Set MORE = FALSE before first calling this routine.
!    MORE will be reset to TRUE and a permutation will be returned.
!    Each new call produces a new permutation until MORE is returned FALSE.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i4_factorial
  integer ( kind = 4 ) m2
  logical more
  integer ( kind = 4 ) n2
  integer ( kind = 4 ), save :: nfact = 0
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) q
  integer ( kind = 4 ), save :: rank = 0
  integer ( kind = 4 ) s
  integer ( kind = 4 ) t

  if ( .not. more ) then

    call i4vec_indicator ( n, p )
    more = .true.
    rank = 1

    nfact = i4_factorial ( n )

  else

    n2 = n
    m2 = rank
    s = n

    do

      q = mod ( m2, n2 )
      t = mod ( m2, 2 * n2 )

      if ( q /= 0 ) then
        exit
      end if

      if ( t == 0 ) then
        s = s - 1
      end if

      m2 = m2 / n2
      n2 = n2 - 1

    end do

    if ( q == t ) then
      s = s - q
    else
      s = s + q - n2
    end if

    call i4_swap ( p(s), p(s+1) )

    rank = rank + 1

    if ( rank == nfact ) then
      more = .false.
    end if

  end if

  return
end
subroutine perm_print ( n, p, title )

!*****************************************************************************80
!
!! PERM_PRINT prints a permutation.
!
!  Example:
!
!    Input:
!
!      P = 7 2 4 1 5 3 6
!
!    Printed output:
!
!      "This is the permutation:"
!
!      1 2 3 4 5 6 7
!      7 2 4 1 5 3 6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects permuted.
!
!    Input, integer ( kind = 4 ) P(N), the permutation, in standard index form.
!
!    Input, character ( len = * ) TITLE, a title.
!    If no title is supplied, then only the permutation is printed.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ), parameter :: inc = 20
  integer ( kind = 4 ) p(n)
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )

    do ilo = 1, n, inc
      ihi = min ( n, ilo + inc - 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(2x,20i4)' ) ( i, i = ilo, ihi )
      write ( *, '(2x,20i4)' ) p(ilo:ihi)
    end do

  else

    do ilo = 1, n, inc
      ihi = min ( n, ilo + inc - 1 )
      write ( *, '(2x,20i4)' ) p(ilo:ihi)
    end do

  end if

  return
end
subroutine perm_random ( n, seed, p )

!*****************************************************************************80
!
!! PERM_RANDOM selects a random permutation of N objects.
!
!  Discussion:
!
!    The routine assumes the objects are labeled 1, 2, ... N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 2002
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
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, integer ( kind = 4 ) P(N), a permutation of ( 1, 2, ..., N ),
!    in standard index form.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  call i4vec_indicator ( n, p )

  do i = 1, n - 1
    j = i4_uniform ( i, n, seed )
    call i4_swap ( p(i), p(j) )
  end do

  return
end
subroutine perm_random2 ( n, seed, p )

!*****************************************************************************80
!
!! PERM_RANDOM2 selects a random permutation of N objects.
!
!  Discussion:
!
!    The input values of P are used as labels; that is, the I-th object
!    is labeled P(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 2002
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
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Input/output, integer ( kind = 4 ) P(N), on input, a list of labels.
!    On output, the list has been permuted randomly.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  do i = 1, n - 1
    j = i4_uniform ( i, n, seed )
    call i4_swap ( p(i), p(j) )
  end do

  return
end
subroutine perm_random3 ( n, seed, p )

!*****************************************************************************80
!
!! PERM_RANDOM3 selects a random permutation of N elements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    Original FORTRAN77 version by James Filliben.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Karla Hoffman, Douglas Shier,
!    Algorithm 564:
!    A Test Problem Generator for Discrete Linear L1 Approximation Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 615-617.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of the array.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, integer ( kind = 4 ) P(N), a permutation, in standard index form.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) iadd
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_RANDOM3 - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of N  = ', n
    write ( *, '(a)' ) '  N must be at least 1!'
    stop
  end if

  if ( n == 1 ) then
    p(1) = 1
    return
  end if

  call i4vec_indicator ( n, p )

  do i = 1, n

    iadd = i4_uniform ( 1, n, seed )

    j = i + iadd

    if ( n < j ) then
      j = j - n
    end if

    if ( i /= j ) then
      call i4_swap ( p(j), p(i) )
    end if

  end do

  return
end
subroutine perm_rank ( n, p, rank )

!*****************************************************************************80
!
!! PERM_RANK computes the rank of a given permutation.
!
!  Discussion:
!
!    This is the same as asking for the step at which PERM_NEXT2
!    would compute the permutation.  The value of the rank will be
!    between 1 and N!.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    Original FORTRAN77 version by Dennis Stanton, Dennis White.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer, 1986,
!    ISBN: 0387963472,
!    LC: QA164.S79.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the set that
!    is permuted by P.
!
!    Input, integer ( kind = 4 ) P(N), a permutation, in standard index form.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the permutation.  This
!    gives the order of the given permutation in the set of all
!    the permutations on N elements.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) count
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) invers(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rem
!
!  Make sure the permutation is a legal one.
!
  call perm_check ( n, p, ierror )
!
!  Compute the inverse permutation.
!
  invers(1:n) = p(1:n)

  call perm_inverse2 ( n, invers )

  rank = 0

  do i = 1, n

    count = 0

    do j = 1, invers(i)
      if ( p(j) < i ) then
        count = count + 1
      end if
    end do

    if ( mod ( rank, 2 ) == 1 ) then
      rem = count
    else
      rem = i - 1 - count
    end if

    rank = i * rank + rem

  end do

  rank = rank + 1

  return
end
subroutine perm_sign ( n, p, p_sign )

!*****************************************************************************80
!
!! PERM_SIGN returns the sign of a permutation.
!
!  Discussion:
!
!    A permutation can always be replaced by a sequence of pairwise
!    transpositions.  A given permutation can be represented by
!    many different such transposition sequences, but the number of
!    such transpositions will always be odd or always be even.
!    If the number of transpositions is even or odd, the permutation is
!    said to be even or odd.
!
!  Example:
!
!    Input:
!
!      N = 9
!      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
!
!    Output:
!
!      P_SIGN = +1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 2000
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
!    Input, integer ( kind = 4 ) N, the number of objects permuted.
!
!    Input, integer ( kind = 4 ) P(N), a permutation, in standard index form.
!
!    Output, integer ( kind = 4 ) P_SIGN, the "sign" of the permutation.
!    +1, the permutation is even,
!    -1, the permutation is odd.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) i4vec_index
  integer ( kind = 4 ) j
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) p_sign
  integer ( kind = 4 ) q(n)

  call perm_check ( n, p, ierror )
!
!  Make a temporary copy of the permutation.
!
  q(1:n) = p(1:n)
!
!  Start with P_SIGN indicating an even permutation.
!  Restore each element of the permutation to its correct position,
!  updating P_SIGN as you go.
!
  p_sign = 1

  do i = 1, n - 1

    j = i4vec_index ( n, q, i )

    if ( j /= i ) then
      call i4_swap ( q(i), q(j) )
      p_sign = -p_sign
    end if

  end do

  return
end
subroutine perm_to_equiv ( n, p, npart, jarray, iarray )

!*****************************************************************************80
!
!! PERM_TO_EQUIV computes the partition induced by a permutation.
!
!  Example:
!
!    Input:
!
!      N = 9
!      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
!
!    Output:
!
!      NPART = 3
!      JARRAY = 4, 3, 2
!      IARRAY = 1, 1, 1, 2  3  2  3  2, 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input, integer ( kind = 4 ) P(N), a permutation, in standard index form.
!
!    Output, integer ( kind = 4 ) NPART, number of subsets in the partition.
!
!    Output, integer ( kind = 4 ) JARRAY(N).  JARRAY(I) is the number of
!    elements in the I-th subset of the partition.
!
!    Output, integer ( kind = 4 ) IARRAY(N).  IARRAY(I) is the class to which
!    element I belongs.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) iarray(n)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jarray(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, ierror )
!
!  Initialize.
!
  iarray(1:n) = 0
  jarray(1:n) = 0

  npart = 0
!
!  Search for the next item J which has not been assigned a subset/orbit.
!
  do j = 1, n

    if ( iarray(j) /= 0 ) then
      cycle
    end if
!
!  Begin a new subset/orbit.
!
    npart = npart + 1
    k = j
!
!  Add the item to the subset/orbit.
!
    do

      jarray(npart) = jarray(npart) + 1
      iarray(k) = npart
!
!  Apply the permutation.  If the permuted object isn't already in the
!  subset/orbit, add it.
!
      k = p(k)

      if ( iarray(k) /= 0 ) then
        exit
      end if

    end do

  end do

  return
end
subroutine perm_to_ytb ( n, p, lambda, a )

!*****************************************************************************80
!
!! PERM_TO_YTB converts a permutation to a Young table.
!
!  Discussion:
!
!    The mapping is not invertible.  In most cases, several permutations
!    correspond to the same table.
!
!  Example:
!
!    N = 7
!    P = 7 2 4 1 5 3 6
!
!    YTB =
!      1 2 3 6
!      4 5
!      7
!
!    LAMBDA = 4 2 1 0 0 0 0
!
!    A = 1 1 1 2 2 1 3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 April 2001
!
!  Author:
!
!    John Burkardt
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
!    Input, integer ( kind = 4 ) N, the integer to be partitioned.
!
!    Input, integer ( kind = 4 ) P(N), a permutation, in standard index form.
!
!    Output, integer ( kind = 4 ) LAMBDA(N).  LAMBDA(I) is the length of
!    the I-th row.
!
!    Output, integer ( kind = 4 ) A(N).  A(I) is the row containing I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  logical another
  integer ( kind = 4 ) compare
  integer ( kind = 4 ) lambda(n)
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) put_index
  integer ( kind = 4 ) put_row
  integer ( kind = 4 ) put_value
!
!  Initialize.
!
  lambda(1:n) = 0
  a(1:n) = 0
!
!  Now insert each item of the permutation.
!
  do put_index = 1, n

    put_value = p(put_index)
    put_row = 1

    do

      another = .false.

      do compare = put_value + 1, n

        if ( a(compare) == put_row ) then
          another = .true.
          a(put_value) = put_row
          a(compare) = 0
          put_value = compare
          put_row = put_row + 1
          exit
        end if

      end do

      if ( .not. another ) then
        exit
      end if

    end do

    a(put_value) = put_row
    lambda(put_row) = lambda(put_row) + 1

  end do

  return
end
subroutine perm_unrank ( n, rank, p )

!*****************************************************************************80
!
!! PERM_UNRANK "unranks" a permutation.
!
!  Discussion:
!
!    That is, given a rank, it computes the corresponding permutation.
!    This is the same as asking for the permutation which PERM_NEXT2
!    would compute at the RANK-th step.
!
!    The value of the rank should be between 1 and N!.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2004
!
!  Author:
!
!    Original FORTRAN77 version by Dennis Stanton, Dennis White.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer, 1986,
!    ISBN: 0387963472,
!    LC: QA164.S79.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the set.
!
!    Input, integer ( kind = 4 ) RANK, the desired rank of the permutation.
!    This gives the order of the given permutation in the set of all
!    the permutations on N elements, using the ordering of PERM_NEXT2.
!
!    Output, integer ( kind = 4 ) P(N), the permutation, in standard index form.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) iprev
  integer ( kind = 4 ) irem
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jdir
  integer ( kind = 4 ) jrank
  integer ( kind = 4 ) nfact
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) rank

  p(1:n) = 0

  nfact = 1

  do i = 1, n
    nfact = nfact * i
  end do

  if ( rank < 1 .or. nfact < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value for RANK.'
    write ( *, '(a,i8)' ) '  RANK must be between 1 and ', nfact
    write ( *, '(a,i8)' ) '  but the input value is ', rank
    stop
  end if

  jrank = rank - 1

  do i = 1, n

    iprev = n + 1 - i
    irem = mod ( jrank, iprev )
    jrank = jrank / iprev

    if ( mod ( jrank, 2 ) == 1 ) then
      j = 0
      jdir = 1
    else
      j = n + 1
      jdir = -1
    end if

    icount = 0

    do

      j = j + jdir

      if ( p(j) == 0 ) then
        icount = icount + 1
      end if

      if ( irem < icount ) then
        exit
      end if

    end do

    p(j) = iprev

  end do

  return
end
subroutine perrin ( n, p )

!*****************************************************************************80
!
!! PERRIN returns the first N values of the Perrin sequence.
!
!  Discussion:
!
!    The Perrin sequence has the initial values:
!
!      P(0) = 3
!      P(1) = 0
!      P(2) = 2
!
!    and subsequent entries are generated by the recurrence
!
!      P(I+1) = P(I-1) + P(I-2)
!
!    Note that if N is a prime, then N must evenly divide P(N).
!
!  Example:
!
!    0   3
!    1   0
!    2   2
!    3   3
!    4   2
!    5   5
!    6   5
!    7   7
!    8  10
!    9  12
!   10  17
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
!    Ian Stewart,
!    "A Neglected Number",
!    Scientific American, Volume 274, pages 102-102, June 1996.
!
!    Ian Stewart,
!    Math Hysteria,
!    Oxford, 2004.
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of terms.
!
!    Output, integer ( kind = 4 ) P(N), the terms 0 through N-1 of the sequence.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(n)

  if ( n < 1 ) then
    return
  end if

  p(1) = 3

  if ( n < 2 ) then
    return
  end if

  p(2) = 0

  if ( n < 3 ) then
    return
  end if

  p(3) = 2

  do i = 4, n
    p(i) = p(i-2) + p(i-3)
  end do

  return
end
subroutine pord_check ( n, a, ierror )

!*****************************************************************************80
!
!! PORD_CHECK checks a matrix representing a partial ordering.
!
!  Discussion:
!
!    The array A is supposed to represent a partial ordering of
!    the elements of a set of N objects.
!
!    For distinct indices I and J, the value of A(I,J) is:
!
!      1, if I << J
!      0, otherwise ( I and J may be unrelated, or perhaps J << I).
!
!    Diagonal elements of A are ignored.
!
!    This routine checks that the values of A do represent
!    a partial ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in the set.
!
!    Input, integer ( kind = 4 ) A(N,N), the partial ordering.
!    1 if I is less than J in the partial ordering,
!    0 otherwise.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors detected.  A is a partial ordering.
!    1, N <= 0.
!    2, 0 < A(I,J) and 0 < A(J,I) for some I and J.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j

  ierror = 0

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PORD_CHECK - Fatal error!'
    write ( *, '(a,i8)' ) '  N must be positive, but N = ', n
    ierror = 1
    stop
  end if

  do i = 1, n
    do j = i+1, n

      if ( 0 < a(i,j) ) then
        if ( 0 < a(j,i) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PORD_CHECK - Fatal error!'
          write ( *, '(a,i8)' ) '  For indices I = ', i
          write ( *, '(a,i8)' ) '  and J = ', j
          write ( *, '(a,i8)' ) '  A(I,J) = ', a(i,j)
          write ( *, '(a,i8)' ) '  A(J,I) = ', a(j,i)
          ierror = 2
          stop
        end if
      end if

    end do
  end do

  return
end
subroutine power_mod ( a, n, m, x )

!*****************************************************************************80
!
!! POWER_MOD computes mod ( A^N, M ).
!
!  Discussion:
!
!    Some programming tricks are used to speed up the computation, and to
!    allow computations in which the value A**N is much too large to
!    store in an integer word.
!
!    First, for efficiency, the power A**N is computed by determining
!    the binary expansion of N, then computing A, A^2, A^4, and so on
!    by repeated squaring, and multiplying only those factors that
!    contribute to A^N.
!
!    Secondly, the intermediate products are immediately "mod'ed", which
!    keeps them small.
!
!    For instance, to compute mod ( A^13, 11 ), we essentially compute
!
!       13 = 1 + 4 + 8
!
!       A^13 = A * A^4 * A^8
!
!       mod ( A^13, 11 ) = mod ( A, 11 ) * mod ( A^4, 11 ) * mod ( A^8, 11 ).
!
!    Fermat's little theorem says that if P is prime, and A is not divisible
!    by P, then ( A^(P-1) - 1 ) is divisible by P.
!
!  Example:
!
!     A  N  M  X
!
!    13  0 31  1
!    13  1 31 13
!    13  2 31 14
!    13  3 31 27
!    13  4 31 10
!    13  5 31  6
!    13  6 31 16
!    13  7 31 22
!    13  8 31  7
!    13  9 31 29
!    13 10 31  5
!    13 11 31  3
!    13 12 31  8
!    13 13 31 11
!    13 14 31 19
!    13 15 31 30
!    13 16 31 18
!    13 17 31 17
!    13 18 31  4
!    13 19 31 21
!    13 20 31 25
!    13 21 31 15
!    13 22 31  9
!    13 23 31 24
!    13 24 31  2
!    13 25 31 26
!    13 26 31 28
!    13 27 31 23
!    13 28 31 20
!    13 29 31 12
!    13 30 31  1
!    13 31 31 13
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, the base of the expression to be tested.
!    0 <= A is required.
!
!    Input, integer ( kind = 4 ) N, the power to which the base is raised.
!    0 <= N is required.
!
!    Input, integer ( kind = 4 ) M, the divisor against which the expression
!    is tested.  0 < M is required.
!
!    Output, integer ( kind = 4 ) X, the remainder when A^N is divided by M.
!    If any input quantity is unacceptable, then the nonsensical value
!    X = -1 is returned.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 8 ) a_square2
  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  integer ( kind = 8 ) m2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncopy
  integer ( kind = 4 ) x
  integer ( kind = 8 ) x2

  if ( a < 0 ) then
    x = -1
    return
  end if

  if ( m <= 0 ) then
    x = -1
    return
  end if

  if ( n < 0 ) then
    x = -1
    return
  end if
!
!  A_SQUARE2 contains the successive squares of A.
!
  a_square2 = int ( a, kind = 8 )
  m2 = int ( m, kind = 8 )
  x2 = int ( 1, kind = 8 )

  ncopy = n

  do while ( 0 < ncopy )

    d = mod ( ncopy, 2 )

    if ( d == 1 ) then
      x2 = mod ( x2 * a_square2, m2 )
    end if

    a_square2 = mod ( a_square2 * a_square2, m2 )
    ncopy = ( ncopy - d ) / 2

  end do
!
!  Ensure that X is nonnegative.
!
  do while ( x2 < 0 )
    x2 = x2 + m
  end do

  x = int ( x2, kind = 4 )

  return
end
subroutine power_series1 ( n, alpha, a, b )

!*****************************************************************************80
!
!! POWER_SERIES1 computes the power series for G(Z) = (1+F(Z))**ALPHA.
!
!  Discussion:
!
!    The power series for F(Z) is given.
!
!    The form of the power series are:
!
!      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
!
!      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2003
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
!    Input, integer ( kind = 4 ) N, the number of terms in the power series.
!
!    Input, real ( kind = 8 ) ALPHA, the exponent of 1+F(Z) in the
!    definition of G(Z).
!
!    Input, real ( kind = 8 ) A(N), the power series coefficients for F(Z).
!
!    Output, real ( kind = 8 ) B(N), the power series coefficients for G(Z).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) v

  do j = 1, n

    v = 0.0D+00

    do i = 1, j-1
      v = v + b(i) * a(j-i) * ( alpha * ( j - i ) - i )
    end do

    b(j) = ( alpha * a(j) + v / real ( j, kind = 8 ) )

  end do

  return
end
subroutine power_series2 ( n, a, b )

!*****************************************************************************80
!
!! POWER_SERIES2 computes the power series for G(Z) = exp(F(Z)) - 1.
!
!  Discussion:
!
!    The power series for F(Z) is given.
!
!    The power series have the form:
!
!      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
!
!      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2003
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
!    Input, integer ( kind = 4 ) N, the number of terms in the power series.
!
!    Input, real ( kind = 8 ) A(N), the power series coefficients for F(Z).
!
!    Output, real ( kind = 8 ) B(N), the power series coefficients for G(Z).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) v

  do j = 1, n

    v = 0.0D+00

    do i = 1, j-1
      v = v + b(i) * a(j-i) * real ( j - i, kind = 8 )
    end do

    b(j) = a(j) + v / real ( j, kind = 8 )

  end do

  return
end
subroutine power_series3 ( n, a, b, c )

!*****************************************************************************80
!
!! POWER_SERIES3 computes the power series for H(Z) = G(F(Z)).
!
!  Discussion:
!
!    The power series for G and H are given.
!
!    We assume that
!
!      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
!      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
!      H(Z) = C1*Z + C2*Z**2 + C3*Z**3 + ... + CN*Z**N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2003
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
!    Input, integer ( kind = 4 ) N, the number of terms in the power series.
!
!    Input, real ( kind = 8 ) A(N), the power series for F.
!
!    Input, real ( kind = 8 ) B(N), the power series for G.
!
!    Output, real ( kind = 8 ) C(N), the power series for H.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) r
  real ( kind = 8 ) v
  real ( kind = 8 ) work(n)

  work(1:n) = b(1) * a(1:n)
!
!  Search for IQ, the index of the first nonzero entry in A.
!
  iq = 0

  do i = 1, n

    if ( a(i) /= 0.0D+00 ) then
      iq = i
      exit
    end if

  end do

  if ( iq /= 0 ) then

    m = 1

    do

      m = m + 1

      if ( n < m * iq ) then
        exit
      end if

      if ( b(m) == 0.0D+00 ) then
        cycle
      end if

      r = b(m) * a(iq)**m
      work(m*iq) = work(m*iq) + r

      do j = 1, n-m*iq

        v = 0.0D+00
        do i = 1, j-1
          v = v + c(i) * a(j-i+iq) * real ( m * ( j - i ) - i, kind = 8 )
        end do

        c(j) = ( real ( m, kind = 8 ) * a(j) + v &
          / real ( j, kind = 8 ) ) / a(iq)

      end do

      do i = 1, n-m*iq
        work(i+m*iq) = work(i+m*iq) + c(i) * r
      end do

    end do

  end if

  c(1:n) = work(1:n)

  return
end
subroutine power_series4 ( n, a, b, c )

!*****************************************************************************80
!
!! POWER_SERIES4 computes the power series for H(Z) = G ( 1/F(Z) ).
!
!  Discussion:
!
!    The routine is given the power series for the functions F and G.
!
!    We assume that
!
!      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
!      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
!      H(Z) = C1*Z + C2*Z**2 + C3*Z**3 + ... + CN*Z**N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2003
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
!    Input, integer ( kind = 4 ) N, the number of terms in the power series.
!
!    Input, real ( kind = 8 ) A(N), the power series for F.
!    A(1) may not be 0.0.
!
!    Input, real ( kind = 8 ) B(N), the power series for G.
!
!    Output, real ( kind = 8 ) C(N), the power series for H.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) work(n)

  if ( a(1) == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POWER_SERIES4 - Fatal error!'
    write ( *, '(a)' ) '  A(1) is zero.'
    stop
  end if

  t = 1.0D+00

  do i = 1, n
    t = t / a(1)
    c(i) = b(i) * t
    work(i) = a(i) * t
  end do

  do k = 2, n
    s = -work(k)
    do i = k, n
      do j = i, n
        c(j) = c(j) + s * c(j+1-k)
        work(j) = work(j) + s * work(j+1-k)
      end do
    end do
  end do

  return
end
function prime ( n )

!*****************************************************************************80
!
!! PRIME returns any of the first PRIME_MAX prime numbers.
!
!  Discussion:
!
!    PRIME_MAX is 1600, and the largest prime stored is 13499.
!
!    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, pages 95-98.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the desired prime number.
!    In general, is should be true that 0 <= N <= PRIME_MAX.
!    N = -1 returns PRIME_MAX, the index of the largest prime available.
!    N = 0 is legal, returning PRIME = 1.
!
!    Output, integer ( kind = 4 ) PRIME, the N-th prime.  If N is out of range,
!    PRIME is returned as -1.
!
  implicit none

  integer ( kind = 4 ), parameter :: prime_max = 1600

  integer ( kind = 4 ), save :: call = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save, dimension ( prime_max ) :: npvec
  integer ( kind = 4 ) prime

  if ( call == 0 ) then

    call = 1

    npvec(1:100) = (/ &
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29, &
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71, &
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113, &
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173, &
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229, &
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281, &
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349, &
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409, &
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463, &
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541 /)

    npvec(101:200) = (/ &
      547,  557,  563,  569,  571,  577,  587,  593,  599,  601, &
      607,  613,  617,  619,  631,  641,  643,  647,  653,  659, &
      661,  673,  677,  683,  691,  701,  709,  719,  727,  733, &
      739,  743,  751,  757,  761,  769,  773,  787,  797,  809, &
      811,  821,  823,  827,  829,  839,  853,  857,  859,  863, &
      877,  881,  883,  887,  907,  911,  919,  929,  937,  941, &
      947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013, &
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, &
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, &
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 /)

    npvec(201:300) = (/ &
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, &
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, &
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, &
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, &
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, &
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, &
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, &
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, &
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987 /)

    npvec(301:400) = (/ &
     1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, &
     2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, &
     2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, &
     2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, &
     2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, &
     2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, &
     2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, &
     2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, &
     2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, &
     2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741 /)

    npvec(401:500) = (/ &
     2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, &
     2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, &
     2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, &
     3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, &
     3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, &
     3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, &
     3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, &
     3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, &
     3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, &
     3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571 /)

    npvec(501:600) = (/ &
     3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, &
     3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, &
     3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, &
     3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, &
     3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, &
     4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, &
     4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, &
     4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, &
     4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, &
     4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409 /)

    npvec(601:700) = (/ &
     4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, &
     4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, &
     4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, &
     4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, &
     4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, &
     4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, &
     4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, &
     5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, &
     5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, &
     5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279 /)

    npvec(701:800) = (/ &
     5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, &
     5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, &
     5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, &
     5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, &
     5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, &
     5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, &
     5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, &
     5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, &
     5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, &
     6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133 /)

    npvec(801:900) = (/ &
     6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, &
     6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, &
     6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, &
     6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, &
     6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, &
     6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, &
     6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, &
     6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, &
     6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, &
     6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997 /)

    npvec(901:1000) = (/ &
     7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, &
     7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, &
     7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, &
     7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, &
     7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, &
     7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, &
     7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, &
     7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, &
     7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, &
     7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /)

    npvec(1001:1100) = (/ &
     7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, &
     8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, &
     8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, &
     8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, &
     8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, &
     8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, &
     8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, &
     8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, &
     8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, &
     8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831 /)

    npvec(1101:1200) = (/ &
     8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, &
     8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, &
     9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, &
     9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, &
     9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, &
     9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, &
     9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, &
     9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, &
     9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, &
     9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733 /)

    npvec(1201:1300) = (/ &
     9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, &
     9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, &
     9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007, &
    10009,10037,10039,10061,10067,10069,10079,10091,10093,10099, &
    10103,10111,10133,10139,10141,10151,10159,10163,10169,10177, &
    10181,10193,10211,10223,10243,10247,10253,10259,10267,10271, &
    10273,10289,10301,10303,10313,10321,10331,10333,10337,10343, &
    10357,10369,10391,10399,10427,10429,10433,10453,10457,10459, &
    10463,10477,10487,10499,10501,10513,10529,10531,10559,10567, &
    10589,10597,10601,10607,10613,10627,10631,10639,10651,10657 /)

    npvec(1301:1400) = (/ &
    10663,10667,10687,10691,10709,10711,10723,10729,10733,10739, &
    10753,10771,10781,10789,10799,10831,10837,10847,10853,10859, &
    10861,10867,10883,10889,10891,10903,10909,10937,10939,10949, &
    10957,10973,10979,10987,10993,11003,11027,11047,11057,11059, &
    11069,11071,11083,11087,11093,11113,11117,11119,11131,11149, &
    11159,11161,11171,11173,11177,11197,11213,11239,11243,11251, &
    11257,11261,11273,11279,11287,11299,11311,11317,11321,11329, &
    11351,11353,11369,11383,11393,11399,11411,11423,11437,11443, &
    11447,11467,11471,11483,11489,11491,11497,11503,11519,11527, &
    11549,11551,11579,11587,11593,11597,11617,11621,11633,11657 /)

    npvec(1401:1500) = (/ &
    11677,11681,11689,11699,11701,11717,11719,11731,11743,11777, &
    11779,11783,11789,11801,11807,11813,11821,11827,11831,11833, &
    11839,11863,11867,11887,11897,11903,11909,11923,11927,11933, &
    11939,11941,11953,11959,11969,11971,11981,11987,12007,12011, &
    12037,12041,12043,12049,12071,12073,12097,12101,12107,12109, &
    12113,12119,12143,12149,12157,12161,12163,12197,12203,12211, &
    12227,12239,12241,12251,12253,12263,12269,12277,12281,12289, &
    12301,12323,12329,12343,12347,12373,12377,12379,12391,12401, &
    12409,12413,12421,12433,12437,12451,12457,12473,12479,12487, &
    12491,12497,12503,12511,12517,12527,12539,12541,12547,12553 /)

   npvec(1501:1600) = (/ &
    12569,12577,12583,12589,12601,12611,12613,12619,12637,12641, &
    12647,12653,12659,12671,12689,12697,12703,12713,12721,12739, &
    12743,12757,12763,12781,12791,12799,12809,12821,12823,12829, &
    12841,12853,12889,12893,12899,12907,12911,12917,12919,12923, &
    12941,12953,12959,12967,12973,12979,12983,13001,13003,13007, &
    13009,13033,13037,13043,13049,13063,13093,13099,13103,13109, &
    13121,13127,13147,13151,13159,13163,13171,13177,13183,13187, &
    13217,13219,13229,13241,13249,13259,13267,13291,13297,13309, &
    13313,13327,13331,13337,13339,13367,13381,13397,13399,13411, &
    13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 /)

  end if

  if ( n == -1 ) then
    prime = prime_max
  else if ( n == 0 ) then
    prime = 1
  else if ( n <= prime_max ) then
    prime = npvec(n)
  else
    prime = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRIME - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal prime index N = ', n
    write ( *, '(a,i8)' ) '  N should be between 1 and PRIME_MAX =', prime_max
    stop
  end if

  return
end
subroutine pythag_triple_next ( i, j, a, b, c )

!*****************************************************************************80
!
!! PYTHAG_TRIPLE_NEXT computes the next Pythagorean triple.
!
!  Example:
!
!     I       J       A       B       C    A^2+B^2 = C^2
!
!     2       1       3       4       5      25
!     3       2       5      12      13     169
!     4       1      15       8      17     289
!     4       3       7      24      25     625
!     5       2      21      20      29     841
!     5       4       9      40      41    1681
!     6       1      35      12      37    1369
!     6       3      27      36      45    2025
!     6       5      11      60      61    3721
!     7       2      45      28      53    2809
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J, the generators.
!    On first call, set I = J = 0.  On repeated calls, leave I and J
!    at their output values from the previous call.
!
!    Output, integer ( kind = 4 ) A, B, C, the next Pythagorean triple.
!    A, B, and C are positive integers which have no common factors,
!    and A*A + B*B = C*C.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
!
!  I starts at 2, and when it increases, increases by 1 and resets J;
!
!  When I is reset, J starts out at 2 if I is odd, or 1 if I is even;
!  Then I is held fixed and J increases by 2, as long as it remains less than I.
!
  if ( i == 0 .and. j == 0 ) then
    i = 2
    j = 1
  else if ( j + 2 < i ) then
    j = j + 2
  else
    i = i + 1
    j = mod ( i, 2 ) + 1
  end if

  a = i * i - j * j
  b = 2 * i * j
  c = i * i + j * j

  return
end
function r8_agm ( a, b )

!*****************************************************************************80
!
!! R8_AGM finds the arithmetic-geometric mean of two numbers.
!
!  Discussion:
!
!    The AGM of (A,B) is produced by the following iteration:
!
!      (A,B) -> ( (A+B)/2, SQRT(A*B) ).
!
!    The sequence of successive values of (A,B) quickly converge to the AGM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the numbers whose AGM is desired.
!    A and B should both be non-negative.
!
!    Output, real ( kind = 8 ) R8_AGM, the AGM of the two numbers.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) r8_agm
  real ( kind = 8 ) tol

  if ( a < 0.0D+00 ) then
    r8_agm = -1.0D+00
    return
  end if

  if ( b < 0.0D+00 ) then
    r8_agm = -1.0D+00
    return
  end if

  if ( a == 0.0D+00 .or. b == 0.0D+00 ) then
    r8_agm = 0.0D+00
    return
  end if

  if ( a == b ) then
    r8_agm = a
    return
  end if

  tol = epsilon ( tol ) * ( a + b + 1.0D+00 )

  a1 = a
  b1 = b

  do

    a2 = ( a1 + b1 ) / 2.0D+00
    b2 = sqrt ( a1 * b1 )

    if ( abs ( a2 - b2 ) < tol ) then
      r8_agm = ( a2 + b2 ) / 2.0D+00
      exit
    end if

    a1 = a2
    b1 = b2

  end do

  return
end
function r8_choose ( n, k )

!*****************************************************************************80
!
!! R8_CHOOSE computes the combinatorial coefficient C(N,K).
!
!  Discussion:
!
!    Real arithmetic is used, and C(N,K) is computed directly, via
!    Gamma functions, rather than recursively.
!
!    C(N,K) is the number of distinct combinations of K objects
!    chosen from a set of N distinct objects.  A combination is
!    like a set, in that order does not matter.
!
!  Example:
!
!    The number of combinations of 2 things chosen from 5 is 10.
!
!    C(5,2) = ( 5 * 4 * 3 * 2 * 1 ) / ( ( 3 * 2 * 1 ) * ( 2 * 1 ) ) = 10.
!
!    The actual combinations may be represented as:
!
!      (1,2), (1,3), (1,4), (1,5), (2,3),
!      (2,4), (2,5), (3,4), (3,5), (4,5).
!
!  Formula:
!
!    C(N,K) = N! / ( (N-K)! * K! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the value of N.
!
!    Input, integer ( kind = 4 ) K, the value of K.
!
!    Output, real ( kind = 8 ) R8_CHOOSE, the value of C(N,K)
!
  implicit none

  real ( kind = 8 ) arg
  real ( kind = 8 ) fack
  real ( kind = 8 ) facn
  real ( kind = 8 ) facnmk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) value

  if ( n < 0 ) then

    value = 0.0D+00

  else if ( k == 0 ) then

    value = 1.0D+00

  else if ( k == 1 ) then

    value = real ( n, kind = 8 )

  else if ( 1 < k .and. k < n-1 ) then

    arg = real ( n + 1, kind = 8 )
    facn = r8_gamma_log ( arg )

    arg = real ( k + 1, kind = 8 )
    fack = r8_gamma_log ( arg )

    arg = real ( n - k + 1, kind = 8 )
    facnmk = r8_gamma_log ( arg )

    value = anint ( exp ( facn - fack - facnmk ) )

  else if ( k == n-1 ) then

    value = real ( n, kind = 8 )

  else if ( k == n ) then

    value = 1.0D+00

  else

    value = 0.0D+00

  end if

  r8_choose = value

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

  real ( kind = 8 ) r8_factorial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  r8_factorial = 1.0D+00

  do i = 1, n
    r8_factorial = r8_factorial * real ( i, kind = 8 )
  end do

  return
end
function r8_fall ( x, n )

!*****************************************************************************80
!
!! R8_FALL computes the falling factorial function [X]_N.
!
!  Discussion:
!
!    Note that the number of "injections" or 1-to-1 mappings from
!    a set of N elements to a set of M elements is [M]_N.
!
!    The number of permutations of N objects out of M is [M]_N.
!
!    Moreover, the Stirling numbers of the first kind can be used
!    to convert a falling factorial into a polynomial, as follows:
!
!      [X]_N = S^0_N + S^1_N * X + S^2_N * X^2 + ... + S^N_N X^N.
!
!    The formula used is:
!
!      [X]_N = X * ( X - 1 ) * ( X - 2 ) * ... * ( X - N + 1 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the falling factorial function.
!
!    Input, integer ( kind = 4 ) N, the order of the falling factorial function.
!    If N = 0, FALL = 1, if N = 1, FALL = X.  Note that if N is
!    negative, a "rising" factorial will be computed.
!
!    Output, real ( kind = 8 ) R8_FALL, the value of the falling
!    factorial function.
!
  implicit none

  real ( kind = 8 ) arg
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_fall
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  value = 1.0D+00

  arg = x

  if ( 0 < n ) then

    do i = 1, n
      value = value * arg
      arg = arg - 1.0D+00
    end do

  else if ( n < 0 ) then

    do i = -1, n, -1
      value = value * arg
      arg = arg + 1.0D+00
    end do

  end if

  r8_fall = value

  return
end
function r8_gamma_log ( x )

!*****************************************************************************80
!
!! R8_GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
!
!  Discussion:
!
!    The program uses rational functions that theoretically approximate
!    LOG(GAMMA(X)) to at least 18 significant decimal digits.  The
!    approximation for 12 < X is from Hart et al, while approximations
!    for X < 12.0 are similar to those in Cody and Hillstrom, but are
!    unpublished.  The accuracy achieved depends on the arithmetic system,
!    the compiler, intrinsic functions, and proper selection of the
!    machine-dependent constants.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 1999
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Kenneth Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
!    Mathematics of Computation,
!    Volume 21, 1967, pages 198-203.
!
!    Kenneth Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi,
!    John Rice, Henry Thatcher, Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Gamma function.
!    X must be positive.
!
!    Output, real ( kind = 8 ) GAMMA_LOG, the logarithm of the Gamma
!    function of X.  If X <= 0.0, or if overflow would occur, the program
!    returns the value XINF, the largest representable floating point number.
!
!  Local Parameters:
!
!    BETA   - radix for the floating-point representation.
!
!    MAXEXP - the smallest positive power of BETA that overflows.
!
!    XBIG   - largest argument for which LN(GAMMA(X)) is representable
!           in the machine, i.e., the solution to the equation
!             LN(GAMMA(XBIG)) = BETA**MAXEXP.
!
!    FRTBIG - Rough estimate of the fourth root of XBIG
!
!
!  Approximate values for some important machines are:
!
!                            BETA      MAXEXP         XBIG
!
!  CRAY-1        (S.P.)        2        8191       9.62D+2461
!  Cyber 180/855
!    under NOS   (S.P.)        2        1070       1.72D+319
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)        2         128       4.08D+36
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)        2        1024       2.55D+305
!  IBM 3033      (D.P.)       16          63       4.29D+73
!  VAX D-Format  (D.P.)        2         127       2.05D+36
!  VAX G-Format  (D.P.)        2        1023       1.28D+305
!
!
!                           FRTBIG
!
!  CRAY-1        (S.P.)   3.13D+615
!  Cyber 180/855
!    under NOS   (S.P.)   6.44D+79
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   1.42D+9
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)   2.25D+76
!  IBM 3033      (D.P.)   2.56D+18
!  VAX D-Format  (D.P.)   1.20D+9
!  VAX G-Format  (D.P.)   1.89D+76
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728D-03, &
     8.4171387781295D-04, &
    -5.952379913043012D-04, &
     7.93650793500350248D-04, &
    -2.777777777777681622553D-03, &
     8.333333333333333331554247D-02, &
     5.7083835261D-03 /)
  real ( kind = 8 ) corr
  real ( kind = 8 ), parameter :: d1 =  -5.772156649015328605195174D-01
  real ( kind = 8 ), parameter :: d2 =   4.227843350984671393993777D-01
  real ( kind = 8 ), parameter :: d4 =   1.791759469228055000094023D+00
  real ( kind = 8 ), parameter :: frtbig = 1.42D+09
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888D+00, &
    2.018112620856775083915565D+02, &
    2.290838373831346393026739D+03, &
    1.131967205903380828685045D+04, &
    2.855724635671635335736389D+04, &
    3.848496228443793359990269D+04, &
    2.637748787624195437963534D+04, &
    7.225813979700288197698961D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064D+00, &
    5.424138599891070494101986D+02, &
    1.550693864978364947665077D+04, &
    1.847932904445632425417223D+05, &
    1.088204769468828767498470D+06, &
    3.338152967987029735917223D+06, &
    5.106661678927352456275255D+06, &
    3.074109054850539556250927D+06 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062D+04, &
    2.426813369486704502836312D+06, &
    1.214755574045093227939592D+08, &
    2.663432449630976949898078D+09, &
    2.940378956634553899906876D+010, &
    1.702665737765398868392998D+011, &
    4.926125793377430887588120D+011, &
    5.606251856223951465078242D+011 /)
  real ( kind = 8 ), parameter :: pnt68 = 0.6796875D+00
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036D+01, &
    1.113332393857199323513008D+03, &
    7.738757056935398733233834D+03, &
    2.763987074403340708898585D+04, &
    5.499310206226157329794414D+04, &
    6.161122180066002127833352D+04, &
    3.635127591501940507276287D+04, &
    8.785536302431013170870835D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942D+02, &
    7.765049321445005871323047D+03, &
    1.331903827966074194402448D+05, &
    1.136705821321969608938755D+06, &
    5.267964117437946917577538D+06, &
    1.346701454311101692290052D+07, &
    1.782736530353274213975932D+07, &
    9.533095591844353613395747D+06 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843D+03, &
    6.393885654300092398984238D+05, &
    4.135599930241388052042842D+07, &
    1.120872109616147941376570D+09, &
    1.488613728678813811542398D+010, &
    1.016803586272438228077304D+011, &
    3.417476345507377132798597D+011, &
    4.463158187419713286462081D+011 /)
  real ( kind = 8 ) r8_gamma_log
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 4.08D+36
  real ( kind = 8 ) xden
  real ( kind = 8 ) xm1
  real ( kind = 8 ) xm2
  real ( kind = 8 ) xm4
  real ( kind = 8 ) xnum
  real ( kind = 8 ) xsq
!
!  Return immediately if the argument is out of range.
!
  if ( x <= 0.0D+00 .or. xbig < x ) then
    r8_gamma_log = huge ( r8_gamma_log )
    return
  end if

  if ( x <= epsilon ( x ) ) then

    res = -log ( x )

  else if ( x <= 1.5D+00 ) then

    if ( x < pnt68 ) then
      corr = -log ( x )
      xm1 = x
    else
      corr = 0.0D+00
      xm1 = ( x - 0.5D+00 ) - 0.5D+00
    end if

    if ( x <= 0.5D+00 .or. pnt68 <= x ) then

      xden = 1.0D+00
      xnum = 0.0D+00

      do i = 1, 8
        xnum = xnum * xm1 + p1(i)
        xden = xden * xm1 + q1(i)
      end do

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

    else

      xm2 = ( x - 0.5D+00 ) - 0.5D+00
      xden = 1.0D+00
      xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

    end if

  else if ( x <= 4.0D+00 ) then

    xm2 = x - 2.0D+00
    xden = 1.0D+00
    xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm2 + p2(i)
      xden = xden * xm2 + q2(i)
    end do

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

  else if ( x <= 12.0D+00 ) then

    xm4 = x - 4.0D+00
    xden = -1.0D+00
    xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm4 + p4(i)
      xden = xden * xm4 + q4(i)
    end do

    res = d4 + xm4 * ( xnum / xden )

  else

    res = 0.0D+00

    if ( x <= frtbig ) then

      res = c(7)
      xsq = x * x

      do i = 1, 6
        res = res / xsq + c(i)
      end do

    end if

    res = res / x
    corr = log ( x )
    res = res + sqrtpi - 0.5D+00 * corr
    res = res + x * ( corr - 1.0D+00 )

  end if

  r8_gamma_log = res

  return
end
function r8_is_int ( r )

!*****************************************************************************80
!
!! R8_IS_INT determines if a real number represents an integer value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the number to be checked.
!
!    Output, logical R8_IS_INT, is TRUE if R is an integer value.
!
  implicit none

  integer ( kind = 4 ) i4_huge
  real ( kind = 8 ) r
  logical r8_is_int

  if ( real ( i4_huge ( ), kind = 8 ) < r ) then
    r8_is_int = .false.
  else if ( r < - real ( i4_huge ( ), kind = 8 ) ) then
    r8_is_int = .false.
  else if ( r == real ( int ( r ), kind = 8 ) ) then
    r8_is_int = .true.
  else
    r8_is_int = .false.
  end if

  return
end
function r8_pi ( )

!*****************************************************************************80
!
!! R8_PI returns the value of pi.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_PI, the value of pi.
!
  implicit none

  real ( kind = 8 ) r8_pi

  r8_pi = 3.141592653589793D+00

  return
end
function r8_rise ( x, n )

!*****************************************************************************80
!
!! R8_RISE computes the rising factorial function [X]^N.
!
!  Discussion:
!
!    [X]^N = X * ( X + 1 ) * ( X + 2 ) * ... * ( X + N - 1 ).
!
!    Note that the number of ways of arranging N objects in M ordered
!    boxes is [M]^N.  (Here, the ordering of the objects in each box matters).
!    Thus, 2 objects in 2 boxes have the following 6 possible arrangements:
!
!      -|12, 1|2, 12|-, -|21, 2|1, 21|-.
!
!    Moreover, the number of non-decreasing maps from a set of
!    N to a set of M ordered elements is [M]^N / N!.  Thus the set of
!    nondecreasing maps from (1,2,3) to (a,b,c,d) is the 20 elements:
!
!      aaa, abb, acc, add, aab, abc, acd, aac, abd, aad
!      bbb, bcc, bdd, bbc, bcd, bbd, ccc, cdd, ccd, ddd.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the rising factorial function.
!
!    Input, integer ( kind = 4 ) N, the order of the rising factorial function.
!    If N = 0, RISE = 1, if N = 1, RISE = X.  Note that if N is
!    negative, a "falling" factorial will be computed.
!
!    Output, real ( kind = 8 ) R8_RISE, the value of the rising factorial
!    function.
!
  implicit none

  real ( kind = 8 ) arg
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_rise
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  value = 1.0D+00

  arg = x

  if ( 0 < n ) then

    do i = 1, n
      value = value * arg
      arg = arg + 1.0D+00
    end do

  else if ( n < 0 ) then

    do i = -1, n, -1
      value = value * arg
      arg = arg - 1.0D+00
    end do

  end if

  r8_rise = value

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP switches two real values.
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
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
subroutine r8_to_cfrac ( r, n, a, p, q )

!*****************************************************************************80
!
!! R8_TO_CFRAC converts a real value to a continued fraction.
!
!  Discussion:
!
!    The routine is given a real number R.  It computes a sequence of
!    continued fraction approximations to R, returning the results as
!    simple fractions of the form P(I) / Q(I).
!
!  Example:
!
!    X = 2 * PI
!    N = 7
!
!    A = [ *, 6,  3,  1,  1,   7,   2,    146,      3 ]
!    P = [ 1, 6, 19, 25, 44, 333, 710, 103993, 312689 ]
!    Q = [ 0, 1,  3,  4,  7,  53, 113,  16551,  49766 ]
!
!    (This ignores roundoff error, which will cause later terms to differ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Norman Richert,
!    Strang's Strange Figures,
!    American Mathematical Monthly,
!    Volume 99, Number 2, February 1992, pages 101-107.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the real value.
!
!    Input, integer ( kind = 4 ) N, the number of convergents to compute.
!
!    Output, integer ( kind = 4 ) A(0:N), the partial quotients.
!
!    Output, integer ( kind = 4 ) P(-1:N), Q(-1:N), the numerators and
!    denominators of the continued fraction approximations.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(-1:n)
  integer ( kind = 4 ) q(-1:n)
  real ( kind = 8 ) r
  real ( kind = 8 ) r_copy
  real ( kind = 8 ) x(0:n)

  if ( r == 0.0D+00 ) then
    a(0:n) = 0
    p(-1:n) = 0
    q(-1:n) = 1
    return
  end if

  r_copy = abs ( r )

  p(-1) = 1
  q(-1) = 0

  p(0) = int ( r_copy )
  q(0) = 1
  x(0) = r_copy
  a(0) = int ( x(0) )

  do i = 1, n
    x(i) = 1.0D+00 / ( x(i-1) - real ( a(i-1), kind = 8 ) )
    a(i) = int ( x(i) )
    p(i) = a(i) * p(i-1) + p(i-2)
    q(i) = a(i) * q(i-1) + q(i-2)
  end do

  if ( r < 0.0D+00 ) then
    p(-1:n) = -p(-1:n)
  end if

  return
end
subroutine r8_to_dec ( dval, dec_digit, mantissa, exponent )

!*****************************************************************************80
!
!! R8_TO_DEC converts a real quantity to a decimal representation.
!
!  Discussion:
!
!    Given the real ( kind = 8 ) value DVAL, the routine computes integers
!    MANTISSA and EXPONENT so that it is approximately true that:
!
!      DVAL = MANTISSA * 10 ** EXPONENT
!
!    In particular, only DEC_DIGIT digits of DVAL are used in constructing the
!    representation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2004
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
!    Input, integer ( kind = 4 ) DEC_DIGIT, the number of decimal digits to use.
!
!    Output, integer ( kind = 4 ) MANTISSA, EXPONENT, the approximate decimal
!    representation of DVAL.
!
  implicit none

  integer ( kind = 4 ) dec_digit
  real ( kind = 8 ) dval
  integer ( kind = 4 ) exponent
  integer ( kind = 4 ) mantissa
  real ( kind = 8 ) mantissa_double
  real ( kind = 8 ) ten1
  real ( kind = 8 ) ten2
!
!  Special cases.
!
  if ( dval == 0.0D+00 ) then
    mantissa = 0
    exponent = 0
    return
  end if
!
!  Factor DVAL = MANTISSA_DOUBLE * 10**EXPONENT
!
  mantissa_double = dval
  exponent = 0
!
!  Now normalize so that
!  10**(DEC_DIGIT-1) <= ABS(MANTISSA_DOUBLE) < 10**(DEC_DIGIT)
!
  ten1 = 10.0D+00**( dec_digit - 1 )
  ten2 = 10.0D+00**dec_digit

  do while ( abs ( mantissa_double ) < ten1 )
    mantissa_double = mantissa_double * 10.0D+00
    exponent = exponent - 1
  end do

  do while ( ten2 <= abs ( mantissa_double ) )
    mantissa_double = mantissa_double / 10.0D+00
    exponent = exponent + 1
  end do
!
!  MANTISSA is the integer part of MANTISSA_DOUBLE, rounded.
!
  mantissa = nint ( mantissa_double )
!
!  Now divide out any factors of ten from MANTISSA.
!
  if ( mantissa /= 0 ) then
    do while ( 10 * ( mantissa / 10 ) == mantissa )
      mantissa = mantissa / 10
      exponent = exponent + 1
    end do
  end if

  return
end
subroutine r8_to_rat ( r, ndig, iatop, iabot )

!*****************************************************************************80
!
!! R8_TO_RAT converts a real value to a rational value.
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
!    15 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the real value to be converted.
!
!    Input, integer ( kind = 4 ) NDIG, the number of decimal digits used.
!
!    Output, integer ( kind = 4 ) IATOP, IABOT, the numerator and denominator
!    of the rational value that approximates the real number.
!
  implicit none

  real ( kind = 8 ) factor
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) iabot
  integer ( kind = 4 ) iatop
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) jfac
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) r

  factor = 10.0D+00**ndig

  if ( 0 < ndig ) then
    ifac = 10**ndig
    jfac = 1
  else
    ifac = 1
    jfac = 10**(-ndig)
  end if

  itop = nint ( r * factor ) * jfac
  ibot = ifac
!
!  Factor out the greatest common factor.
!
  itemp = i4_gcd ( itop, ibot )

  iatop = itop / itemp
  iabot = ibot / itemp

  return
end
subroutine r8_to_s_left ( rval, s )

!*****************************************************************************80
!
!! R8_TO_S_LEFT represents a real using 14 left_justified characters.
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
!    Input, real ( kind = 8 ) RVAL, a real number.
!
!    Output, character ( len = * ) S, a left-justified character variable
!    containing the representation of RVAL.
!
  implicit none

  character ( len = 14 ) chrtmp
  integer ( kind = 4 ) i
  real ( kind = 8 ) rval
  character ( len = * ) s
!
!  We can't seem to write directly into the string because of compiler
!  quibbles.
!
  if ( real ( int ( rval ), kind = 8 ) == rval .and. &
       abs ( rval ) < 1.0D+13 ) then

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
function r8_uniform ( a, b, seed )

!*****************************************************************************80
!
!! R8_UNIFORM returns a scaled pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM, a number strictly between A and B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge ( )
  end if

  r8_uniform = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
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
!    Pierre LEcuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge ( )
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8mat_det ( n, a, det )

!*****************************************************************************80
!
!! R8MAT_DET finds the determinant of an N by N R8MAT.
!
!  Discussion:
!
!    A brute force calculation is made.
!
!    This routine should only be used for small matrices, since this
!    calculation requires the summation of N! products of N numbers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of A.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) det
  logical even
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarray(n)
  logical more
  real ( kind = 8 ) term

  more = .false.
  det = 0.0D+00

  do

    call perm_next ( n, iarray, more, even )

    if ( even ) then
      term = 1.0D+00
    else
      term = -1.0D+00
    end if

    do i = 1, n
      term = term * a(i,iarray(i))
    end do

    det = det + term

    if ( .not. more ) then
      exit
    end if

  end do

  return
end
subroutine r8mat_perm ( n, a, p )

!*****************************************************************************80
!
!! R8MAT_PERM permutes the rows and columns of a square R8MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2002
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer ( kind = 4 ) P(N), a permutation to be applied to the rows
!    and columns.  P(I) is the new number of row and column I.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) is
  real ( kind = 8 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) p(n)

  call perm_cycle ( n, 1, p, is, nc )

  do i = 1, n

    i1 = -p(i)

    if ( 0 < i1 ) then

      lc = 0

      do

        i1 = p(i1)
        lc = lc + 1

        if ( i1 <= 0 ) then
          exit
        end if

      end do

      i1 = i

      do j = 1, n

        if ( p(j) <= 0 ) then

          j2 = j
          k = lc

          do

            j1 = j2
            it = a(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( p(j1) )

              call r8_swap ( a(i1,j1), it )

              if ( j1 /= j2 ) then
                cycle
              end if

              k = k - 1

              if ( i1 == i ) then
                exit
              end if

            end do

            j2 = abs ( p(j2) )

            if ( k == 0 ) then
              exit
            end if

          end do

        end if

      end do

    end if

  end do
!
!  Restore the positive signs of the data.
!
  p(1:n) = abs ( p(1:n) )

  return
end
subroutine r8mat_perm2 ( m, n, a, p, q )

!*****************************************************************************80
!
!! R8MAT_PERM2 permutes rows and columns of a rectangular R8MAT, in place.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 October 1999
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
!    Input, integer ( kind = 4 ) M, number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) N, number of columns in the matrix.
!
!    Input/output, real ( kind = 8 ) A(M,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer ( kind = 4 ) P(M), the row permutation.  P(I) is the
!    new number of row I.
!
!    Input, integer ( kind = 4 ) Q(N), the column permutation.  Q(I) is the new
!    number of column I.  Note that the routine allows you to pass a single
!    array as both P and Q.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) q(n)
  real ( kind = 8 ) t

  call perm_cycle ( m, 1, p, is, nc )

  if ( 0 < q(1) ) then
    call perm_cycle ( n, 1, q, is, nc )
  end if

  do i = 1, m

    i1 = -p(i)

    if ( 0 < i1 ) then

      lc = 0

      do

        i1 = p(i1)
        lc = lc + 1

        if ( i1 <= 0 ) then
          exit
        end if

      end do

      i1 = i

      do j = 1, n

        if ( q(j) <= 0 ) then

          j2 = j
          k = lc

          do

            j1 = j2
            t = a(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( q(j1) )

              call r8_swap ( a(i1,j1), t )

              if ( j1 /= j2 ) then
                cycle
              end if

              k = k - 1

              if ( i1 == i ) then
                exit
              end if

            end do

            j2 = abs ( q(j2) )

            if ( k == 0 ) then
              exit
            end if

          end do

        end if

      end do

    end if

  end do

  p(1:m) = abs ( p(1:m) )

  if ( q(1) <= 0 ) then

    q(1:n) = abs ( q(1:n) )

  end if

  return
end
subroutine r8mat_permanent ( n, a, perm )

!*****************************************************************************80
!
!! R8MAT_PERMANENT computes the permanent of an R8MAT.
!
!  Discussion:
!
!    The permanent function is similar to the determinant.  Recall that
!    the determinant of a matrix may be defined as the sum of all the
!    products:
!
!      S * A(1,I(1)) * A(2,I(2)) * ... * A(N,I(N))
!
!    where I is any permutation of the columns of the matrix, and S is the
!    sign of the permutation.  By contrast, the permanent function is
!    the (unsigned) sum of all the products
!
!      A(1,I(1)) * A(2,I(2)) * ... * A(N,I(N))
!
!    where I is any permutation of the columns of the matrix.  The only
!    difference is that there is no permutation sign multiplying each summand.
!
!    Symbolically, then, the determinant of a 2 by 2 matrix
!
!      a b
!      c d
!
!    is a*d-b*c, whereas the permanent of the same matrix is a*d+b*c.
!
!
!    The permanent is invariant under row and column permutations.
!    If a row or column of the matrix is multiplied by S, then the
!      permanent is likewise multiplied by S.
!    If the matrix is square, then the permanent is unaffected by
!      transposing the matrix.
!    Unlike the determinant, however, the permanent does change if
!      one row is added to another, and it is not true that the
!      permanent of the product is the product of the permanents.
!
!
!    Note that if A is a matrix of all 1's and 0's, then the permanent
!    of A counts exactly which permutations hit exactly 1's in the matrix.
!    This fact can be exploited for various combinatorial purposes.
!
!    For instance, setting the diagonal of A to 0, and the offdiagonals
!    to 1, the permanent of A counts the number of derangements of N
!    objects.
!
!    Setting the diagonal of A to 0, and ensuring that the offdiagonal
!    entries are symmetric, then A is the adjacency matrix of a graph,
!    and its permanent counts the number of perfect matchings.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2003
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the value of the matrix.
!
!    Output, real ( kind = 8 ) PERM, the value of the permanent of the matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iadd
  integer ( kind = 4 ) iwork(n)
  logical more
  integer ( kind = 4 ) ncard
  real ( kind = 8 ) p
  real ( kind = 8 ) perm
  real ( kind = 8 ) sgn
  real ( kind = 8 ) work(n)
  real ( kind = 8 ) z

  more = .false.

  do i = 1, n
    work(i) = a(i,n) - 0.5D+00 * sum ( a(i,1:n) )
  end do

  p = 0.0D+00
  sgn = -1.0D+00

  do

    sgn = - sgn
    call subset_gray_next ( n - 1, iwork, more, ncard, iadd )

    if ( ncard /= 0 ) then
      z = real ( 2 * iwork(iadd) - 1, kind = 8 )
      work(1:n) = work(1:n) + z * a(1:n,iadd)
    end if

    p = p + sgn * product ( work )

    if ( .not. more ) then
      exit
    end if

  end do

  perm = p * real ( 4 * mod ( n, 2 ) - 2, kind = 8 )

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 May 2004
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2003
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
  logical r8_is_int
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
      write ( ctemp(j2), '(i7,7x)') j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( r8_is_int ( a(i,j) ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8poly ( n, a, x0, iopt, val )

!*****************************************************************************80
!
!! R8POLY performs operations on real polynomials in power or factorial form.
!
!  Discussion:
!
!    The power sum form of a polynomial is
!
!      P(X) = A1 + A2 * X + A3 * X**2 + ... + (AN+1) * X**N
!
!    The Taylor expansion at C has the form
!
!      P(X) = A1 + A2 * (X-C) + A3 * (X-C)**2+... + (AN+1) * (X-C)**N
!
!    The factorial form of a polynomial is
!
!      P(X) = A1 + A2 * X + A3 * (X) * (X-1) + A4 * (X) * (X-1) * (X-2) + ...
!        + (AN+1) * (X) * (X-1) *...* (X-N+1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
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
!    Input, integer ( kind = 4 ) N, the number of coefficients in the polynomial
!    (in other words, the polynomial degree + 1)
!
!    Input/output, real ( kind = 8 ) A(N), the coefficients of the polynomial.
!    Depending on the option chosen, these coefficients may be overwritten
!    by those of a different form of the polynomial.
!
!    Input, real ( kind = 8 ) X0, for IOPT = -1, 0, or positive, the value of
!    the argument at which the polynomial is to be evaluated, or the
!    Taylor expansion is to be carried out.
!
!    Input, integer ( kind = 4 ) IOPT, a flag describing which algorithm is to
!    be carried out:
!
!    -3: Reverse Stirling.  Input the coefficients of the polynomial in
!    factorial form, output them in power sum form.
!
!    -2: Stirling.  Input the coefficients in power sum
!    form, output them in factorial form.
!
!    -1: Evaluate a polynomial which has been input
!    in factorial form.
!
!    0:  Evaluate a polynomial input in power sum form.
!
!    1 or more:  Given the coefficients of a polynomial in
!    power sum form, compute the first IOPT coefficients of
!    the polynomial in Taylor expansion form.
!
!    Output, real ( kind = 8 ) VAL, for IOPT = -1 or 0, the value of the
!    polynomial at the point X0.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iopt
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n1
  real ( kind = 8 ) val
  real ( kind = 8 ) w
  real ( kind = 8 ) x0
  real ( kind = 8 ) z

  n1 = min ( n, iopt )
  n1 = max ( 1, n1 )

  if ( iopt < -1 ) then
    n1 = n
  end if

  eps = real ( mod ( max ( -iopt, 0 ), 2 ), kind = 8 )

  w = - real ( n, kind = 8 ) * eps

  if ( -2 < iopt ) then
    w = w + x0
  end if

  do m = 1, n1

    val = 0.0D+00
    z = w

    do i = m, n
      z = z + eps
      val = a(n+m-i) + z * val
      if ( iopt /= 0 .and. iopt /= -1 ) then
        a(n+m-i) = val
      end if
    end do

    if ( iopt < 0 ) then
      w = w + 1.0D+00
    end if

  end do

  return
end
subroutine r8poly_degree ( na, a, degree )

!*****************************************************************************80
!
!! R8POLY_DEGREE returns the degree of a polynomial in power sum form.
!
!  Discussion:
!
!    The power sum form of a polynomial is:
!
!      p(x) = a(0) + a(1) * x + ... + a(n-1) * x**(n-1) + a(n) * x**(n)
!
!    The degree of a polynomial is the index of the highest power
!    of X with a nonzero coefficient.
!
!    The degree of a constant polynomial is 0.  The degree of the
!    zero polynomial is debatable, but this routine returns the
!    degree as 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:NA), the coefficients of the polynomials.
!
!    Output, integer ( kind = 4 ) DEGREE, the degree of A.
!
  implicit none

  integer ( kind = 4 ) na

  real ( kind = 8 ) a(0:na)
  integer ( kind = 4 ) degree

  degree = na

  do while ( 0 < degree )

    if ( a(degree) /= 0.0D+00 ) then
      return
    end if

    degree = degree - 1

  end do

  return
end
subroutine r8poly_div ( na, a, nb, b, nq, q, nr, r )

!*****************************************************************************80
!
!! R8POLY_DIV computes the quotient and remainder of two polynomials.
!
!  Discussion:
!
!    The polynomials are assumed to be stored in power sum form.
!
!    The power sum form of a polynomial is:
!
!      p(x) = a(0) + a(1) * x + ... + a(n-1) * x**(n-1) + a(n) * x**(n)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:NA), the coefficients of the polynomial
!    to be divided.
!
!    Input, integer ( kind = 4 ) NB, the dimension of B.
!
!    Input, real ( kind = 8 ) B(0:NB), the coefficients of the divisor
!    polynomial.
!
!    Output, integer ( kind = 4 ) NQ, the degree of Q.
!    If the divisor polynomial is zero, NQ is returned as -1.
!
!    Output, real ( kind = 8 ) Q(0:NA-NB), contains the quotient of A/B.
!    If A and B have full degree, Q should be dimensioned Q(0:NA-NB).
!    In any case, Q(0:NA) should be enough.
!
!    Output, integer ( kind = 4 ) NR, the degree of R.
!    If the divisor polynomial is zero, NR is returned as -1.
!
!    Output, real ( kind = 8 ) R(0:NB-1), contains the remainder of A/B.
!    If B has full degree, R should be dimensioned R(0:NB-1).
!    Otherwise, R will actually require less space.
!
  implicit none

  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb

  real ( kind = 8 ) a(0:na)
  real ( kind = 8 ) a2(0:na)
  real ( kind = 8 ) b(0:nb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) na2
  integer ( kind = 4 ) nb2
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nr
  real ( kind = 8 ) q(0:*)
  real ( kind = 8 ) r(0:*)

  call r8poly_degree ( na, a, na2 )
  call r8poly_degree ( nb, b, nb2 )

  if ( b(nb2) == 0.0D+00 ) then
    nq = -1
    nr = -1
    return
  end if

  a2(0:na) = a(0:na)

  nq = na2 - nb2
  nr = nb2 - 1

  do i = nq, 0, -1
    q(i) = a2(i+nb2) / b(nb2)
    a2(i+nb2) = 0.0D+00
    a2(i:i+nb2-1) = a2(i:i+nb2-1) - q(i) * b(0:nb2-1)
  end do

  r(0:nr) = a2(0:nr)

  return
end
subroutine r8poly_f2p ( n, a )

!*****************************************************************************80
!
!! R8POLY_F2P converts a real polynomial from factorial form to power sum form.
!
!  Discussion:
!
!    The (falling) factorial form is
!
!      p(x) =   a(1)
!             + a(2) * x
!             + a(3) * x*(x-1)
!             ...
!             + a(n) * x*(x-1)*...*(x-(n-2))
!
!    The power sum form is
!
!      p(x) = a(1) + a(2)*x + a(3)*x**2 + ... + a(n)*x**(n-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of A.
!
!    Input/output, real ( kind = 8 ) A(N), on input, the polynomial
!    coefficients in factorial form.  On output, the polynomial
!    coefficients in power sum form.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) val
  real ( kind = 8 ) w
  real ( kind = 8 ) z

  w = - real ( n, kind = 8 )

  do m = 1, n

    val = 0.0D+00
    z = w

    do i = m, n
      z = z + 1.0D+00
      val = a(n+m-i) + z * val
      a(n+m-i) = val
    end do

    w = w + 1.0D+00

  end do

  return
end
subroutine r8poly_fval ( n, a, x, val )

!*****************************************************************************80
!
!! R8POLY_FVAL evaluates a real polynomial in factorial form.
!
!  Discussion:
!
!    The (falling) factorial form of a polynomial is:
!
!      p(x) = a(1)
!           + a(2)  *x
!           + a(3)  *x*(x-1)
!           +...
!           + a(n-1)*x*(x-1)*(x-2)...*(x-(n-3))
!           + a(n)  *x*(x-1)*(x-2)...*(x-(n-3))*(x-(n-2))
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of A.
!
!    Input, real ( kind = 8 ) A(N), the coefficients of the polynomial.
!    A(1) is the constant term.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomial is
!    to be evaluated.
!
!    Output, real ( kind = 8 ) VAL, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) val
  real ( kind = 8 ) x

  val = 0.0D+00
  do i = 1, n
    val = a(n+1-i) + ( x - n + i ) * val
  end do

  return
end
subroutine r8poly_mul ( na, a, nb, b, c )

!*****************************************************************************80
!
!! R8POLY_MUL computes the product of two real polynomials A and B.
!
!  Discussion:
!
!    The polynomials are in power sum form.
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1) * x + ... + a(n-1) * x**(n-1) + a(n) * x**(n)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:NA), the coefficients of the first
!    polynomial factor.
!
!    Input, integer ( kind = 4 ) NB, the dimension of B.
!
!    Input, real ( kind = 8 ) B(0:NB), the coefficients of the second
!    polynomial factor.
!
!    Output, real ( kind = 8 ) C(0:NA+NB), the coefficients of A * B.
!
  implicit none

  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb

  real ( kind = 8 ) a(0:na)
  real ( kind = 8 ) b(0:nb)
  real ( kind = 8 ) c(0:na+nb)
  real ( kind = 8 ) d(0:na+nb)
  integer ( kind = 4 ) i

  d(0:na+nb) = 0.0D+00

  do i = 0, na
    d(i:i+nb) = d(i:i+nb) + a(i) * b(0:nb)
  end do

  c(0:na+nb) = d(0:na+nb)

  return
end
subroutine r8poly_n2p ( n, a, xarray )

!*****************************************************************************80
!
!! R8POLY_N2P converts a real polynomial from Newton form to power sum form.
!
!  Discussion:
!
!    This is done by shifting all the Newton abscissas to zero.
!
!    Actually, what happens is that the abscissas of the Newton form
!    are all shifted to zero, which means that A is the power sum
!    polynomial description and A, XARRAY is the Newton polynomial
!    description.  It is only because all the abscissas are shifted to
!    zero that A can be used as both a power sum and Newton polynomial
!    coefficient array.
!
!    The Newton form of a polynomial is described by an array of N coefficients
!    A and N abscissas X:
!
!      p(x) =   a(1)
!             + a(2) * (x-x(1))
!             + a(3) * (x-x(1)) * (x-x(2))
!             ...
!             + a(n) * (x-x(1)) * (x-x(2)) * ... * (x-x(n-1))
!
!    X(N) does not occur explicitly in the formula for the evaluation of p(x),
!    although it is used in deriving the coefficients A.
!
!    The power sum form of a polynomial is:
!
!      p(x) = a(1) + a(2)*x + ... + a(n-1)*x**(n-2) + a(n)*x**(n-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of A.
!
!    Input/output, real ( kind = 8 ) A(N).  On input, the coefficients
!    of the polynomial in Newton form, and on output, the coefficients
!    in power sum form.
!
!    Input/output, real ( kind = 8 ) XARRAY(N).  On input, the abscissas of
!    the Newton form of the polynomial.  On output, these values
!    have all been set to zero.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) xarray(n)
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  do i = 1, n
    call r8poly_nx ( n, a, xarray, zero )
  end do

  return
end
subroutine r8poly_nval ( n, a, xarray, x, val )

!*****************************************************************************80
!
!! R8POLY_NVAL evaluates a real polynomial in Newton form.
!
!  Discussion:
!
!    The Newton form of a polynomial is;
!
!      p(x) = a(1)
!           + a(2)  *(x-x1)
!           + a(3)  *(x-x1)*(x-x2)
!           +...
!           + a(n-1)*(x-x1)*(x-x2)*(x-x3)...*(x-x(n-2))
!           + a(n)  *(x-x1)*(x-x2)*(x-x3)...*(x-x(n-2))*(x-x(n-1))
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of A.
!
!    Input, real ( kind = 8 ) A(N), the coefficients of the polynomial.
!    A(1) is the constant term.
!
!    Input, real ( kind = 8 ) XARRAY(N-1), the N-1 points X which are part
!    of the definition of the polynomial.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomial
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) VAL, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) val
  real ( kind = 8 ) x
  real ( kind = 8 ) xarray(n-1)

  val = a(n)
  do i = n-1, 1, -1
    val = a(i) + ( x - xarray(i) ) * val
  end do

  return
end
subroutine r8poly_nx ( n, a, xarray, x )

!*****************************************************************************80
!
!! R8POLY_NX replaces one of the base points in a polynomial in Newton form.
!
!  Discussion:
!
!    The Newton form of a polynomial is described by an array of N coefficients
!    A and N abscissas X:
!
!      p(x) =   a(1)
!             + a(2) * (x-x(1))
!             + a(3) * (x-x(1)) * (x-x(2))
!             ...
!             + a(n) * (x-x(1)) * (x-x(2)) * ... * (x-x(n-1))
!
!    X(N) does not occur explicitly in the formula for the evaluation of p(x),
!    although it is used in deriving the coefficients A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of A.
!
!    Input/output, real ( kind = 8 ) A(N), the polynomial coefficients
!    of the Newton form.
!
!    Input/output, real ( kind = 8 ) XARRAY(N), the set of abscissas that
!    are part of the Newton form of the polynomial.  On output,
!    the abscissas have been shifted up one index, so that
!    the first location now holds X, and the original contents
!    of XARRAY(N) have been discarded.
!
!    Input, real ( kind = 8 ) X, the new point to be shifted into XARRAY.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x
  real ( kind = 8 ) xarray(n)

  do i = n-1, 1, -1
    a(i) = a(i) + ( x - xarray(i) ) * a(i+1)
  end do

  do i = n, 2, -1
    xarray(i) = xarray(i-1)
  end do

  xarray(1) = x

  return
end
subroutine r8poly_p2f ( n, a )

!*****************************************************************************80
!
!! R8POLY_P2F converts a real polynomial from power sum form to factorial form.
!
!  Discussion:
!
!    The power sum form is
!
!      p(x) = a(1) + a(2) * x + a(3) * x**2 + ... + a(n) * x**(n-1)
!
!    The (falling) factorial form is
!
!      p(x) =   a(1)
!             + a(2) * x
!             + a(3) * x * (x-1)
!             ...
!             + a(n) * x * (x-1) *...* (x-(n-2))
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of A.
!
!    Input/output, real ( kind = 8 ) A(N), on input, the polynomial
!    coefficients in the power sum form, on output, the polynomial
!    coefficients in factorial form.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) val

  do m = 1, n
    val = 0.0D+00
    do i = m, n
      val = a(n+m-i) + real ( m - 1, kind = 8 ) * val
      a(n+m-i) = val
    end do
  end do

  return
end
subroutine r8poly_p2n ( n, a, xarray )

!*****************************************************************************80
!
!! R8POLY_P2N converts a real polynomial from power sum form to Newton form.
!
!  Discussion:
!
!    This is done by shifting all the Newton abscissas from zero.
!
!    The power sum form of a polynomial is:
!
!      p(x) = a(1) + a(2) * x + ... + a(n-1) * x**(n-2) + a(n) * x**(n-1)
!
!    The Newton form of a polynomial is described by an array of N coefficients
!    A and N abscissas X:
!
!      p(x) =   a(1)
!             + a(2) * (x-x(1))
!             + a(3) * (x-x(1)) * (x-x(2))
!             ...
!             + a(n) * (x-x(1)) * (x-x(2)) * ... * (x-x(n-1))
!
!    X(N) does not occur explicitly in the formula for the evaluation of p(x),
!    although it is used in deriving the coefficients A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of A.
!
!    Input/output, real ( kind = 8 ) A(N).  On input, the coefficients
!    of the polynomial in power sum form, and on output, the
!    coefficients in Newton form.
!
!    Input/output, real ( kind = 8 ) XARRAY(N).  On input, the desired
!    abscissas of the Newton form of the polynomial.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) xarray(n)
  real ( kind = 8 ) work(n)

  work(1:n) = 0.0D+00

  do i = n, 1, -1
    call r8poly_nx ( n, a, work, xarray(i) )
  end do

  return
end
subroutine r8poly_p2t ( n, a, x )

!*****************************************************************************80
!
!! R8POLY_P2T converts a real polynomial from power sum form to Taylor form.
!
!  Discussion:
!
!    The power sum form is
!
!      p(x) = a(1) + a(2)*x + a(3)*x**2 + ... + a(n)*x**(n-1)
!
!    The Taylor form of a polynomial based at X0 is
!
!      p(x) =   a(1)
!             + a(2) * (x-x0)
!             + a(3) * (x-x0)**2
!             ...
!             + a(n) * (x-x0)**(n-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of A.
!
!    Input/output, real ( kind = 8 ) A(N), on input, the coefficients in
!    power sum form, and on output, the coefficients in Taylor form.
!
!    Input, real ( kind = 8 ) X, the point at which the Taylor form of the
!    polynomial is to be based.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) val
  real ( kind = 8 ) x

  do m = 1, n
    val = 0.0D+00
    do i = m, n
      val = a(n+m-i) + x * val
      a(n+m-i) = val
    end do
  end do

  return
end
subroutine r8poly_power ( na, a, p, b )

!*****************************************************************************80
!
!! R8POLY_POWER computes a positive integer power of a polynomial.
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NA, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:NA), the polynomial to be raised to the power.
!
!    Input, integer ( kind = 4 ) P, the nonnegative power to which A is raised.
!
!    Output, real ( kind = 8 ) B(0:P*NA), the power of the polynomial.
!
  implicit none

  integer ( kind = 4 ) na
  integer ( kind = 4 ) p

  real ( kind = 8 ) a(0:na)
  real ( kind = 8 ) b(0:p*na)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nonzer
!
!  Zero out B.
!
  b(0:p*na) = 0.0D+00
!
!  Search for the first nonzero element in A.
!
  nonzer = 0

  do i = 0, na
    if ( a(i) /= 0.0D+00 ) then
      nonzer = i
      exit
    end if
  end do

  if ( nonzer == 0 ) then
    return
  end if

  b(0) = a(nonzer)**p

  do i = 1, p*(na-nonzer)

    if ( i + nonzer <= na ) then
      b(i) = real ( i * p, kind = 8 ) * b(0) * a(i+nonzer)
    else
      b(i) = 0.0D+00
    end if

    do j = 1, i-1

      if ( j+nonzer <= na ) then
        b(i) = b(i) - real ( i - j, kind = 8 ) * a(j+nonzer) * b(i-j)
      end if

      if ( i-j+nonzer <= na ) then
        b(i) = b(i) + real ( i - j, kind = 8 ) * real ( p, kind = 8 ) &
          * b(j) * a(i-j+nonzer)
      end if

    end do

    b(i) = b(i) / ( real ( i, kind = 8 ) * a(nonzer) )

  end do
!
!  Shift B up.
!
  do i = p*nonzer, p*na
    b(i) = b(i-p*nonzer)
  end do

  do i = 0, p * nonzer-1
    b(i) = 0.0D+00
  end do

  return
end
subroutine r8poly_print ( n, a, title )

!*****************************************************************************80
!
!! R8POLY_PRINT prints out a polynomial.
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1) * x + ... + a(n-1) * x**(n-1) + a(n) * x**(n)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:N), the polynomial coefficients.
!    A(0) is the constant term and
!    A(N) is the coefficient of X**N.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) mag
  integer ( kind = 4 ) n2
  character plus_minus
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  call r8poly_degree ( n, a, n2 )

  if ( a(n2) < 0.0D+00 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( a(n2) )

  if ( 2 <= n2 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x ^ '', i3 )' ) &
      plus_minus, mag, n2
  else if ( n2 == 1 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x'' )' ) plus_minus, mag
  else if ( n2 == 0 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6 )' ) plus_minus, mag
  end if

  do i = n2-1, 0, -1

    if ( a(i) < 0.0D+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( a(i) )

    if ( mag /= 0.0D+00 ) then

      if ( 2 <= i ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x ^ '', i3 )' ) &
          plus_minus, mag, i
      else if ( i == 1 ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x'' )' ) plus_minus, mag
      else if ( i == 0 ) then
        write ( *, ' ( ''         '', a1, g14.6 )' ) plus_minus, mag
      end if
    end if

  end do

  return
end
subroutine r8poly_pval ( n, a, x, val )

!*****************************************************************************80
!
!! R8POLY_PVAL evaluates a real polynomial in power sum form.
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1) * x + ... + a(n-1) * x**(n-1) + a(n) * x**(n)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of A.
!
!    Input, real ( kind = 8 ) A(0:N), the coefficients of the polynomial.
!    A(0) is the constant term.
!
!    Input, real ( kind = 8 ) X, the point at which the polynomial
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) VAL, the value of the polynomial at X.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) val
  real ( kind = 8 ) x

  val = 0.0D+00
  do i = n, 0, -1
    val = val * x + a(i)
  end do

  return
end
subroutine r8poly_t2p ( n, a, x )

!*****************************************************************************80
!
!! R8POLY_T2P converts a real polynomial from Taylor form to power sum form.
!
!  Discussion:
!
!    The Taylor form of a polynomial based at X0 is
!
!      p(x) =   a(1)
!             + a(2) * (x-x0)
!             + a(3) * (x-x0)**2
!             ...
!             + a(n) * (x-x0)**(n-1)
!
!    The power sum form is
!
!      p(x) = a(1) + a(2)*x + a(3)*x**2 + ... + a(n)*x**(n-1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of A.
!
!    Input/output, real ( kind = 8 ) A(N).  On input, the coefficients
!    in Taylor form, and on output, the coefficients in power sum form.
!
!    Input, real ( kind = 8 ) X, the point at which the Taylor form
!    polynomial is based.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x

  do i = n, 1, -1
    do j = i, n-1
      a(j) = a(j) - a(j+1) * x
    end do
  end do

  return
end
subroutine r8vec_backtrack ( n, maxstack, stack, x, indx, k, nstack, ncan )

!*****************************************************************************80
!
!! R8VEC_BACKTRACK supervises a backtrack search for an R8VEC.
!
!  Discussion:
!
!    The routine tries to construct a real vector one index at a time,
!    using possible candidates as supplied by the user.
!
!    At any time, the partially constructed vector may be discovered to be
!    unsatisfactory, but the routine records information about where the
!    last arbitrary choice was made, so that the search can be
!    carried out efficiently, rather than starting out all over again.
!
!    First, call the routine with INDX = 0 so it can initialize itself.
!
!    Now, on each return from the routine, if INDX is:
!      1, you've just been handed a complete candidate vector;
!         Admire it, analyze it, do what you like.
!      2, please determine suitable candidates for position X(K).
!         Return the number of candidates in NCAN(K), adding each
!         candidate to the end of STACK, and increasing NSTACK.
!      3, you're done.  Stop calling the routine;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 July 2004
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
!    Input, integer ( kind = 4 ) N, the number of positions to be filled in
!    the vector.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the maximum length of the stack.
!
!    Input, real ( kind = 8 ) STACK(MAXSTACK), a list of all current
!    candidates for all positions 1 through K.
!
!    Input/output, real ( kind = 8 ) X(N), the partially filled in
!    candidate vector.
!
!    Input/output, integer ( kind = 4 ) INDX, a communication flag.
!    On input,
!      0, to begin a backtracking search.
!      2, the requested candidates for position K have been added to
!      STACK, and NCAN(K) was updated.
!    On output:
!      1, a complete output vector has been determined and returned in X(1:N);
!      2, candidates are needed for position X(K);
!      3, no more possible vectors exist.
!
!    Input/output, integer ( kind = 4 ) K, the index in X that we are trying
!    to fill.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!
!    Input/output, integer ( kind = 4 ) NCAN(N), lists the current number of
!    candidates for all positions 1 through K.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) indx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ncan(n)
  integer ( kind = 4 ) nstack
  real ( kind = 8 ) stack(maxstack)
  real ( kind = 8 ) x(n)
!
!  If this is the first call, request a candidate for position 1.
!
  if ( indx == 0 ) then
    k = 1
    nstack = 0
    indx = 2
    return
  end if
!
!  Examine the stack.
!
  do
!
!  If there are candidates for position K, take the first available
!  one off the stack, and increment K.
!
!  This may cause K to reach the desired value of N, in which case
!  we need to signal the user that a complete set of candidates
!  is being returned.
!
    if ( 0 < ncan(k) ) then

      x(k) = stack(nstack)
      nstack = nstack - 1

      ncan(k) = ncan(k) - 1

      if ( k /= n ) then
        k = k + 1
        indx = 2
      else
        indx = 1
      end if

      exit
!
!  If there are no candidates for position K, then decrement K.
!  If K is still positive, repeat the examination of the stack.
!
    else

      k = k - 1

      if ( k <= 0 ) then
        indx = 3
        exit
      end if

    end if

  end do

  return
end
subroutine r8vec_frac ( n, a, k, afrac )

!*****************************************************************************80
!
!! R8VEC_FRAC searches for the K-th smallest entry in an R8VEC.
!
!  Discussion:
!
!    Hoare's algorithm is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, A is the array to search.
!    On output, the elements of A have been somewhat rearranged.
!
!    Input, integer ( kind = 4 ) K, the fractile to be sought.  If K = 1, the
!    minimum entry is sought.  If K = N, the maximum is sought.  Other values
!    of K search for the entry which is K-th in size.  K must be at
!    least 1, and no greater than N.
!
!    Output, real ( kind = 8 ) AFRAC, the value of the K-th fractile of A.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) afrac
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iryt
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) left
  real ( kind = 8 ) x

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_FRAC - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal nonpositive value of N = ', n
    stop
  end if

  if ( k <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_FRAC - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal nonpositive value of K = ', k
    stop
  end if

  if ( n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_FRAC - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal N < K, K = ', k
    stop
  end if

  left = 1
  iryt = n

  do

    if ( iryt <= left ) then
      afrac = a(k)
      exit
    end if

    x = a(k)
    i = left
    j = iryt

    do

      if ( j < i ) then
        if ( j < k ) then
          left = i
        end if
        if ( k < i ) then
          iryt = j
        end if
        exit
      end if
!
!  Find I so that X <= A(I)
!
      do while ( a(i) < x )
        i = i + 1
      end do
!
!  Find J so that A(J) <= X
!
      do while ( x < a(j) )
        j = j - 1
      end do

      if ( i <= j ) then
        call r8_swap ( a(i), a(j) )
        i = i + 1
        j = j - 1
      end if

    end do

  end do

  return
end
subroutine r8vec_indicator ( n, a )

!*****************************************************************************80
!
!! R8VEC_INDICATOR sets an R8VEC to the indicator vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, real ( kind = 8 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = real ( i, kind = 8 )
  end do

  return
end
subroutine r8vec_mirror_next ( n, a, done )

!*****************************************************************************80
!
!! R8VEC_MIRROR_NEXT steps through all sign variations of an R8VEC.
!
!  Discussion:
!
!    In normal use, the user would set every element of A to be positive.
!    The routine will take the input value of A, and output a copy in
!    which the signs of one or more entries have been changed.  Repeatedly
!    calling the routine with the output from the previous call will generate
!    every distinct "variation" of A; that is, all possible sign variations.
!
!    When the output variable DONE is TRUE (or equal to 1), then the
!    output value of A_NEW is the last in the series.
!
!    Note that A may have some zero values.  The routine will essentially
!    ignore such entries; more exactly, it will not stupidly assume that -0
!    is a proper "variation" of 0!
!
!    Also, it is possible to call this routine with the signs of A set
!    in any way you like.  The routine will operate properly, but it
!    will nonethess terminate when it reaches the value of A in which
!    every nonzero entry has negative sign.
!
!
!    More efficient algorithms using the Gray code seem to require internal
!    memory in the routine, which is not one of MATLAB's strong points,
!    or the passing back and forth of a "memory array", or the use of
!    global variables, or unnatural demands on the user.  This form of
!    the routine is about as clean as I can make it.
!
!  Example:
!
!      Input         Output
!    ---------    --------------
!    A            A_NEW     DONE
!    ---------    --------  ----
!     1  2  3     -1  2  3  false
!    -1  2  3      1 -2  3  false
!     1 -2  3     -1 -2  3  false
!    -1 -2  3      1  2 -3  false
!     1  2 -3     -1  2 -3  false
!    -1  2 -3      1 -2 -3  false
!     1 -2 -3     -1 -2 -3  false
!    -1 -2 -3      1  2  3  true
!
!     1  0  3     -1  0  3  false
!    -1  0  3      1  0 -3  false
!     1  0 -3     -1  0 -3  false
!    -1  0 -3      1  0  3  true
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 May 2004
!
!  Author:
!
!    John Burkardt
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
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, real ( kind = 8 ) A(N), a vector of real numbers.
!    On output, some signs have been changed.
!
!    Output, logical DONE, is TRUE if the input vector A was the last element
!    in the series (every entry was nonpositive); the output vector is reset
!    so that all entries are nonnegative, but presumably the ride is over!
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  logical done
  integer ( kind = 4 ) i
  integer ( kind = 4 ) positive
!
!  Seek the first strictly positive entry of A.
!
  positive = 0
  do i = 1, n
    if ( 0.0D+00 < a(i) ) then
      positive = i
      exit
    end if
  end do
!
!  If there is no strictly positive entry of A, there is no successor.
!
  if ( positive == 0 ) then
    a(1:n) = -a(1:n)
    done = .true.
    return
  end if
!
!  Otherwise, negate A up to the positive entry.
!
  a(1:positive) = -a(1:positive)
  done = .false.

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(i8,g14.6)' ) i, a(i)
  end do

  return
end
subroutine r8vec_uniform ( n, a, b, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM returns a scaled pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer variable.
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
!    Input, integer ( kind = 4 ) M, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer variable.
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
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine random_initialize ( seed )

!*****************************************************************************80
!
!! RANDOM_INITIALIZE initializes the FORTRAN 90 random number seed.
!
!  Discussion:
!
!    If you don't initialize the random number generator, its behavior
!    is not specified.  If you initialize it simply by:
!
!      call random_seed ( )
!
!    its behavior is not specified.  On the DEC ALPHA, if that's all you
!    do, the same random number sequence is returned.  In order to actually
!    try to scramble up the random number generator a bit, this routine
!    goes through the tedious process of getting the size of the random
!    number seed, making up values based on the current time, and setting
!    the random number seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED.
!    If SEED is zero on input, then you're asking this routine to come up
!    with a seed value, which is returned as output.
!    If SEED is nonzero on input, then you're asking this routine to
!    use the input value of SEED to initialize the random number generator.
!
  implicit none

  integer ( kind = 4 ) count
  integer ( kind = 4 ) count_max
  integer ( kind = 4 ) count_rate
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), allocatable :: seed_vector(:)
  integer ( kind = 4 ) seed_size
  real ( kind = 8 ) t
!
!  Initialize the random number seed.
!
  call random_seed ( )
!
!  Determine the size of the random number seed.
!
  call random_seed ( size = seed_size )
!
!  Allocate a seed of the right size.
!
  allocate ( seed_vector(seed_size) )

  if ( seed /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_INITIALIZE'
    write ( *, '(a,i12)' ) '  Initialize RANDOM_NUMBER with user SEED = ', seed

  else

    call system_clock ( count, count_rate, count_max )

    seed = count

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_INITIALIZE'
    write ( *, '(a,i12)' ) &
      '  Initialize RANDOM_NUMBER with arbitrary SEED = ', seed

  end if
!
!  Now set the seed.
!
  seed_vector(1:seed_size) = seed

  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Call the random number routine a bunch of times.
!
  do i = 1, 100
    call random_number ( harvest = t )
  end do

  return
end
subroutine rat_add ( itop1, ibot1, itop2, ibot2, itop, ibot, ierror )

!*****************************************************************************80
!
!! RAT_ADD adds two rational values.
!
!  Discussion:
!
!    The routine computes
!
!      ITOP/IBOT = ITOP1/IBOT1 + ITOP2/IBOT2
!
!    while trying to avoid integer overflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ITOP1, IBOT1, the first value to add.
!
!    Input, integer ( kind = 4 ) ITOP2, IBOT2, the second value to add.
!
!    Output, integer ( kind = 4 ) ITOP, IBOT, the sum.
!
!    Output, integer ( kind = 4 ) IERROR.
!    0, no error occurred.
!    1, an error occurred.  The addition of the two values
!    requires a numerator or denominator larger than the
!    maximum legal integer.
!
  implicit none

  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ibot1
  integer ( kind = 4 ) ibot2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) itop1
  integer ( kind = 4 ) itop2
  integer ( kind = 4 ) jbot1
  integer ( kind = 4 ) jbot2
  integer ( kind = 4 ) jbot3
  integer ( kind = 4 ) jtop1
  integer ( kind = 4 ) jtop2

  i_max = i4_huge ( )

  ierror = 0

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
  if ( real ( i_max, kind = 8 ) &
    < abs ( real ( jtop1, kind = 8 ) * real ( jbot2, kind = 8 ) ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_ADD - Fatal error!'
    write ( *, '(a)' ) '  Overflow of top of rational sum.'
    itop = 0
    stop
  end if

  jtop1 = jtop1 * jbot2

  if ( real ( i_max, kind = 8 ) &
    < abs ( real ( jtop2, kind = 8 ) * real ( jbot1, kind = 8 ) ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_ADD - Fatal error!'
    write ( *, '(a)' ) '  Overflow of top of rational sum.'
    itop = 0
    stop
  end if

  jtop2 = jtop2 * jbot1

  if ( real ( i_max, kind = 8 ) &
    < abs ( real ( jtop1, kind = 8 ) + real ( jtop2, kind = 8 ) ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_ADD - Fatal error!'
    write ( *, '(a)' ) '  Overflow of top of rational sum.'
    itop = 0
    stop
  end if

  itop = jtop1 + jtop2
!
!  Check the bottom for overflow.
!
  if ( real ( i_max, kind = 8 ) < &
    abs ( real ( jbot1, kind = 8 ) * real ( jbot2, kind = 8 ) &
    * real ( jbot3, kind = 8 ) ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_ADD - Fatal error!'
    write ( *, '(a)' ) '  Overflow of bottom of rational sum.'
    ibot = 1
    stop
  end if

  ibot = jbot1 * jbot2 * jbot3
!
!  Put the fraction in lowest terms.
!
  itemp = i4_gcd ( itop, ibot )
  itop = itop / itemp
  ibot = ibot / itemp
!
!  The bottom should be positive.
!
  if ( ibot < 0 ) then
    ibot = -ibot
    itop = -itop
  end if

  return
end
subroutine rat_div ( itop1, ibot1, itop2, ibot2, itop, ibot, ierror )

!*****************************************************************************80
!
!! RAT_DIV divides one rational value by another.
!
!  Discussion:
!
!    The routine computes
!
!      ITOP / IBOT = ( ITOP1 / IBOT1 ) / ( ITOP2 / IBOT2 ).
!
!    while avoiding integer overflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2003
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
!    Output, integer ( kind = 4 ) IERROR.
!    0, no error occurred.
!    1, an error occurred.  One of the quantities IBOT1, IBOT2,
!    or ITOP2 is zero, or the result of the division
!    requires a numerator or denominator larger than the
!    maximum legal integer.
!
  implicit none

  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) i_max
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

  ierror = 0

  i_max = i4_huge ( )

  if ( ibot1 == 0 .or. itop2 == 0 .or. ibot2 == 0 ) then
    ierror = 1
    return
  end if

  if ( itop1 == 0 ) then
    itop = 0
    ibot = 1
    return
  end if
!
!  Make copies of the input arguments, since we will change them.
!  Implicitly invert the divisor fraction here.  The rest of
!  the code will be a multiply operation.
!
  jbot1 = ibot1
  jbot2 = itop2
  jtop1 = itop1
  jtop2 = ibot2
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
!  The fraction (ITOP1*IBOT2)/(IBOT1*ITOP2) is in lowest terms.
!
!  Check the top for overflow.
!
  if ( real ( i_max, kind = 8 ) &
    < abs ( real ( jtop1, kind = 8 ) * real ( jtop2, kind = 8 ) ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_DIV - Fatal error!'
    write ( *, '(a)' ) '  Overflow of top of rational fraction.'
    itop = 0
    stop
  end if

  itop = jtop1 * jtop2
!
!  Check the bottom IBOT1*ITOP2 for overflow.
!
  if ( real ( i_max, kind = 8 ) &
    < abs ( real ( jbot1, kind = 8 ) * real ( jbot2, kind = 8 ) ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_DIV - Fatal error!'
    write ( *, '(a)' ) '  Overflow of bottom of rational fraction.'
    ibot = 1
    stop
  end if

  ibot = jbot1 * jbot2
!
!  The bottom should be positive.
!
  if ( ibot < 0 ) then
    ibot = -ibot
    itop = -itop
  end if
!
!  The fraction is ITOP/IBOT with no loss of accuracy.
!
  return
end
subroutine rat_farey ( n, max_frac, num_frac, a, b )

!*****************************************************************************80
!
!! RAT_FAREY computes the N-th row of the Farey fraction table.
!
!  Example:
!
!    N = 5
!
!    NUM_FRAC = 11
!    A =  0  1  1  1  2  1  3  2  3  4  1
!    B =  1  5  4  3  5  2  5  3  4  5  1
!
!  Discussion:
!
!    In this form of the Farey fraction table, fractions in row N lie between
!    0 and 1, are in lowest terms, and have a denominator that is no greater
!    than N.  Row N is computed directly, and does not require the computation
!    of previous rows.
!
!    The data satisfy the relationship:
!
!      A(K+1) * B(K) - A(K) * B(K+1) = 1
!
!    The number of items in the N-th row is roughly N**2 / PI**2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 1, Fundamental Algorithms,
!    Addison Wesley, 1968, page 157.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the desired row number.  N must be positive.
!
!    Input, integer ( kind = 4 ) MAX_FRAC, the maximum number of fractions to
!    compute.
!
!    Output, integer ( kind = 4 ) NUM_FRAC, the number of fractions computed.
!
!    Output, integer ( kind = 4 ) A(MAX_FRAC), B(MAX_FRAC), contains the
!    NUM_FRAC numerators and denominators of the N-th row of the Farey
!    fraction table.
!
  implicit none

  integer ( kind = 4 ) max_frac

  integer ( kind = 4 ) a(max_frac)
  integer ( kind = 4 ) b(max_frac)
  integer ( kind = 4 ) c
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) num_frac

  if ( n <= 0 ) then
    num_frac = 0
    return
  end if

  if ( max_frac <= 0 ) then
    num_frac = 0
    return
  end if

  k = 1
  a(k) = 0
  b(k) = 1

  if ( max_frac <= 1 ) then
    num_frac = k
    return
  end if

  k = 2
  a(k) = 1
  b(k) = n

  do while ( k < max_frac )

    if ( a(k) == 1 .and. b(k) == 1 ) then
      exit
    end if

    k = k + 1
    c = ( b(k-2) + n ) / b(k-1)
    a(k) = c * a(k-1) - a(k-2)
    b(k) = c * b(k-1) - b(k-2)

  end do

  num_frac = k

  return
end
subroutine rat_farey2 ( n, a, b )

!*****************************************************************************80
!
!! RAT_FAREY2 computes the next row of the Farey fraction table.
!
!  Example:
!
!    Input:
!
!      N = 3
!      A =  0  1  1  2  1
!      B =  1  3  2  3  1
!
!    Output:
!
!      A =  0  1  1  2  1  3  2  3  1
!      B =  1  4  3  5  2  5  3  4  1
!
!  Discussion:
!
!    In this form of the Farey fraction table, fractions in row N lie between
!    0 and 1, and are in lowest terms.  For every adjacent pair of input
!    fractions, A1/B1 and A2/B2, the mediant (A1+A2)/(B1+B2) is computed
!    and inserted between them.
!
!    The number of items in the N-th row is 1+2**(N-1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the input row number.  N must be
!    nonnegative.  If N is zero, then the input is ignored, and the entries of
!    row 1 are computed directly.
!
!    Input/output, integer ( kind = 4 ) A(1+2**N), B(1+2**N).
!    On input, entries 1 through 1+2**(N-1) contain the entries of row N.
!    On output, entries 1 through 1+2**N contain the entries of row N+1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(1+2**n)
  integer ( kind = 4 ) b(1+2**n)
  integer ( kind = 4 ) i

  if ( n == 0 ) then
    a(1) = 0
    b(1) = 1
    a(2) = 1
    b(2) = 1
    return
  end if
!
!  Shift the current data.
!
  do i = 1+2**(n-1), 1, -1
    a(2*i-1) = a(i)
    b(2*i-1) = b(i)
  end do
!
!  Compute the mediants.
!
  do i = 2, 2**n, 2
    a(i) = a(i-1) + a(i+1)
    b(i) = b(i-1) + b(i+1)
  end do

  return
end
subroutine rat_mul ( itop1, ibot1, itop2, ibot2, itop, ibot, ierror )

!*****************************************************************************80
!
!! RAT_MUL multiplies two fractions.
!
!  Discussion:
!
!    The routine computes
!
!      ITOP / IBOT = ( ITOP1 / IBOT1 ) * ( ITOP2 / IBOT2 ).
!
!    while avoiding integer overflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ITOP1, IBOT1, the first factor.
!
!    Input, integer ( kind = 4 ) ITOP2, IBOT2, the second factor.
!
!    Output, integer ( kind = 4 ) ITOP, IBOT, the product.
!
!    Output, integer ( kind = 4 ) IERROR.
!    0, no error occurred.
!    1, an error occurred.  The multiplication of the two values
!    requires a numerator or denominator larger than the
!    maximum legal integer.
!
  implicit none

  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) i_max
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

  ierror = 0

  i_max = i4_huge ( )

  if ( itop1 == 0 .or. itop2 == 0 ) then
    itop = 0
    ibot = 1
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
  if ( real ( i_max, kind = 8 ) &
    < abs ( real ( jtop1, kind = 8 ) * real ( jtop2, kind = 8 ) ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_MUL - Fatal error!'
    write ( *, '(a)' ) '  Overflow of top of rational product.'
    itop = 0
    stop
  end if

  itop = jtop1 * jtop2
!
!  Check the bottom IBOT1*IBOT2 for overflow.
!
  if ( real ( i_max, kind = 8 ) &
    < abs ( real ( jbot1, kind = 8 ) * real ( jbot2, kind = 8 ) ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_MUL - Fatal error!'
    write ( *, '(a)' ) '  Overflow of bottom of rational product.'
    ibot = 1
    stop
  end if

  ibot = jbot1 * jbot2
!
!  The bottom should be positive.
!
  if ( ibot < 0 ) then
    ibot = -ibot
    itop = -itop
  end if
!
!  The fraction is ITOP/IBOT with no loss of accuracy.
!
  return
end
subroutine rat_normalize ( a, b )

!*****************************************************************************80
!
!! RAT_NORMALIZE normalizes a rational number.
!
!  Discussion:
!
!    If A = B = 0, return.
!
!    If A = 0 (and B nonzero) set B => 1 and return.
!
!    If A nonzero, and B = 0, then A => 1 and return.
!
!    If A = B, then set A => 1, B => 1 and return.
!
!    If B < 0, then A => -A, B => -B.
!
!    If 1 < C = GCD(|A|,|B|), A => A/C, B => B/C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) A, B, the rational number.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i4_gcd
!
!  Cases where one or both is 0.
!
  if ( a == 0 .and. b == 0 ) then
    return
  else if ( a == 0 .and. b /= 0 ) then
    b = 1
    return
  else if ( a /= 0 .and. b == 0 ) then
    a = 1
    return
  end if

  if ( a == b ) then
    a = 1
    b = 1
    return
  end if

  if ( b < 0 ) then
    a = -a
    b = -b
  end if

  c = i4_gcd ( abs ( a ), abs ( b ) )

  if ( 1 < c ) then
    a = a / c
    b = b / c
  end if

  return
end
subroutine rat_sum_formula ( n, a, b )

!*****************************************************************************80
!
!! RAT_SUM_FORMULA computes the formulas for sums of powers of integers.
!
!  Example:
!
!    N = 6
!
!        1    2    3    4    5    6    7
!    -----------------------------------
!    0 | 1    0    0    0    0    0    0
!      |
!    1 | 1    1    0    0    0    0    0
!      | 2    2
!      |
!    2 | 1    1    1    0    0    0    0
!      | 3    2    6
!      |
!    3 | 1    1    1    0    0    0    0
!      | 4    2    4
!      |
!    4 | 1    1    1    0   -1    0    0
!      | 5    2    3        30
!      |
!    5 | 1    1    5    0   -1    0    0
!      | 6    2   12        12
!      |
!    6 | 1    1    1    0   -1    0    1
!      | 7    2    2         6        42
!
!    The interpretation of row 2, for instance, is:
!
!      sum ( 1 <= I <= N ) I**2 = 1/3 N**3 + 1/2 N**2 + 1/6 N
!
!    This suggests that a more sensible way to display the table
!    is to reverse the order of the entries in the row, so that
!    the entry in column J is the coeficient of N**J, which is
!    not the case now.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Owens,
!    Sums of Powers of Integers,
!    Mathematics Magazine,
!    Volume 65, Number 1, February 1992, pages 38-40.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows of coefficients
!    to compute.
!
!    Output, integer ( kind = 4 ) A(0:N,N+1), B(0:N,N+1), the numerator and
!    denominator of the coefficients.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(0:n,1:n+1)
  integer ( kind = 4 ) asum
  integer ( kind = 4 ) b(0:n,1:n+1)
  integer ( kind = 4 ) bsum
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j

  a(0,1) = 1
  b(0,1) = 1
  a(0,2:n+1) = 0
  b(0,2:n+1) = 1

  do i = 1, n

    asum = 0
    bsum = 0
!
!  Subdiagonal entries are multiples of entries above them.
!
    do j = 1, i

      call rat_mul ( a(i-1,j), b(i-1,j), i, i+2-j, a(i,j), b(i,j), ierror )

      call rat_add ( asum, bsum, a(i,j), b(i,j), asum, bsum, ierror )

    end do
!
!  Diagonal entry is 1 - sum of previous entries in row.
!
    asum = -asum
    call rat_add ( 1, 1, asum, bsum, a(i,i+1), b(i,i+1), ierror )
!
!  Superdiagonal entries are zero.
!
    a(i,i+2:n+1) = 0
    b(i,i+2:n+1) = 1

  end do

  return
end
subroutine rat_to_cfrac ( ip, iq, m, n, a, ierror )

!******************************************************************************
!
!! RAT_TO_CFRAC converts a rational value to a continued fraction.
!
!  Discussion:
!
!    The routine is given a rational number represented by IP/IQ, and
!    computes the monic or "simple" continued fraction representation
!    with integer coefficients of the number:
!
!      A(1) + 1/ (A(2) + 1/ (A(3) + ... + 1/A(N) ...))
!
!    The user must dimension A to a value M which is "large enough".
!    The actual number of terms needed in the continued fraction
!    representation cannot be known beforehand.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson,
!    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi,
!    John Rice, Henry Thatcher, Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IP, IQ, the numerator and denominator of the
!    rational value whose continued fraction representation is
!    desired.
!
!    Input, integer ( kind = 4 ) M, the dimension of A.  If M is not great
!    enough, the algorithm may run out of space.
!
!    Output, integer ( kind = 4 ) N, the actual number of entries used in A.
!
!    Output, integer ( kind = 4 ) A(M), contains the continued fraction
!    representation of the number.
!
!    Output, integer ( kind = 4 ) IERROR, error indicator.  0 if no error,
!    1 if there was an error, namely, M is not large enough.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) a(m)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) jq
  integer ( kind = 4 ) n

  jp = ip
  jq = iq

  n = 0

  do

    n = n + 1

    if ( m < n ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RAT_TO_CFRAC - Fatal error!'
      write ( *, '(a)' ) '  M < N.'
      write ( *, '(a)' ) '  M = ', m
      write ( *, '(a)' ) '  N = ', n
      ierror = 1
      stop
    end if

    a(n) = jp / jq
    jp = mod ( jp, jq )

    if ( jp == 0 ) then
      return
    end if

    n = n + 1

    if ( m < n ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RAT_TO_CFRAC - Fatal error!'
      write ( *, '(a)' ) '  M < N.'
      write ( *, '(a)' ) '  M = ', m
      write ( *, '(a)' ) '  N = ', n
      ierror = 1
      stop
    end if

    a(n) = jq / jp
    jq = mod ( jq, jp )

    if ( jq == 0 ) then
      exit
    end if

  end do

  return
end
subroutine rat_to_dec ( rat_top, rat_bot, mantissa, exponent )

!*****************************************************************************80
!
!! RAT_TO_DEC converts a rational to a decimal representation.
!
!  Discussion:
!
!    A rational value is represented by RAT_TOP / RAT_BOT.
!
!    A decimal value is represented as MANTISSA * 10**EXPONENT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RAT_TOP, RAT_BOT, the rational value.
!
!    Output, integer ( kind = 4 ) MANTISSA, EXPONENT, the decimal number.
!
  implicit none

  integer ( kind = 4 ) exponent
  integer ( kind = 4 ) gcd
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) i4_huge
  integer ( kind = 4 ) mantissa
  real ( kind = 8 ) r
  real ( kind = 8 ) r_max
  integer ( kind = 4 ) rat_bot
  integer ( kind = 4 ) rat_bot2
  integer ( kind = 4 ) rat_top
  integer ( kind = 4 ) rat_top2
  integer ( kind = 4 ) s

  if ( rat_top == 0 ) then
    mantissa = 0
    exponent = 0
    return
  end if

  gcd = i4_gcd ( rat_top, rat_bot )
  rat_top2 = rat_top / gcd
  rat_bot2 = rat_bot / gcd

  if ( rat_bot2 < 0 ) then
    rat_top2 = -rat_top2
    rat_bot2 = -rat_bot2
  end if

  if ( rat_bot2 == 1 ) then
    mantissa = rat_top2
    exponent = 0
    return
  end if

  exponent = 0

  do while ( mod ( rat_bot2, 10 ) == 0 )
    exponent = exponent - 1
    rat_bot2 = rat_bot2 / 10
  end do

  do while ( mod ( rat_top2, 10 ) == 0 )
    exponent = exponent + 1
    rat_top2 = rat_top2 / 10
  end do

  r = real ( rat_top2, kind = 8 ) / real ( rat_bot2, kind = 8 )

  if ( r < 0.0D+00 ) then
    s = -1
    r = -r
  else
    s = 1
  end if

  r_max = real ( i4_huge ( ), kind = 8 ) / 10.0D+00

  do while ( r /= real ( int ( r ), kind = 8 ) .and. r < r_max )
    r = r * 10.0D+00
    exponent = exponent - 1
  end do

  mantissa = s * int ( r, kind = 8 )

  return
end
subroutine rat_to_r8 ( a, b, r )

!*****************************************************************************80
!
!! RAT_TO_R8 converts rational values to real values.
!
!  Example:
!
!    A    B    R
!   --   --    ---
!    1    2    0.5
!    7    5    1.4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the rational quantity to be converted.
!
!    Output, real ( kind = 8 ) R, the value of the rational quantity
!    as a real number.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) r

  if ( b == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_TO_R8 - Warning!'
    write ( *, '(a)' ) '  The input fraction to be converted had a'
    write ( *, '(a)' ) '  zero denominator.'
    r = 0.0D+00
  else
    r = real ( a, kind = 8 ) / real ( b, kind = 8 )
  end if

  return
end
subroutine rat_to_s_left ( a, b, s )

!*****************************************************************************80
!
!! RAT_TO_S_LEFT returns a left-justified representation of A/B.
!
!  Discussion:
!
!    If the ratio is negative, a minus sign precedes A.
!    A slash separates A and B.
!
!    Note that if A is nonzero and B is 0, S will
!    be returned as "Inf" or "-Inf" (Infinity), and if both
!    A and B are zero, S will be returned as "NaN"
!    (Not-a-Number).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the numerator and denominator.
!
!    Output, character ( len = * ) S, a left-justified string
!    containing the representation of A/B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  character ( len = * ) s
  character ( len = 25 ) s2
!
!  Take care of simple cases right away.
!
  if ( a == 0 ) then

    if ( b /= 0 ) then
      s2 = '0'
    else
      s2= 'NaN'
    end if

  else if ( b == 0 ) then

    if ( 0 < a ) then
      s2 = 'Inf'
    else
      s2 = '-Inf'
    end if
!
!  Make copies of A and B.
!
  else

    if ( b == 1 ) then
      write ( s2, '(i12)' ) a
    else
      write ( s2, '(i12, ''/'', i12)' ) a, b
    end if

    call s_blank_delete ( s2 )

  end if

  s = s2

  return
end
function rat_width ( a, b )

!*****************************************************************************80
!
!! RAT_WIDTH returns the "width" of a rational number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the rational number.
!
!    Output, integer ( kind = 4 ) RAT_WIDTH, the "width" of the rational number.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) abs_a
  integer ( kind = 4 ) abs_b
  integer ( kind = 4 ) b
  integer ( kind = 4 ) rat_width
  integer ( kind = 4 ) ten_pow
  integer ( kind = 4 ) value

  value = 1
  ten_pow = 10

  if ( a == 0 ) then
    rat_width = 1
    return
  end if

  abs_a = abs ( a )

  do while ( ten_pow <= abs_a )
    value = value + 1
    ten_pow = ten_pow * 10
  end do
!
!  If the fraction is negative, a minus sign will be prepended to the
!  numerator.
!
  if ( a * b < 0 ) then
    value = value + 1
    ten_pow = ten_pow * 10
  end if

  abs_b = abs ( b )

  do while ( ten_pow <= abs_b )
    value = value + 1
    ten_pow = ten_pow * 10
  end do

  rat_width = value

  return
end
subroutine ratmat_det ( n, iatop, iabot, idtop, idbot, ierror )

!*****************************************************************************80
!
!! RATMAT_DET finds the determinant of an N by N matrix of rational entries.
!
!  Discussion:
!
!    The brute force method is used.
!
!    This routine should only be used for small matrices, since this
!    calculation requires the summation of N! products of N numbers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of A.
!
!    Input, integer ( kind = 4 ) IATOP(N,N), IABOT(N,N), the numerators
!    and denominators of the entries of the matrix.
!
!    Output, integer ( kind = 4 ) IDTOP, IDBOT, the determinant of the matrix,
!    expressed as IDTOP/IDBOT.
!
!    Output, integer ( kind = 4 ) IERROR.
!    0, the determinant was computed.
!    1, an overflow error occurred, and the determinant was not
!    computed.
!
  implicit none

  integer ( kind = 4 ) n

  logical even
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iabot(n,n)
  integer ( kind = 4 ) iatop(n,n)
  integer ( kind = 4 ) iarray(n)
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ibot1
  integer ( kind = 4 ) ibot2
  integer ( kind = 4 ) idbot
  integer ( kind = 4 ) idtop
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) itop1
  integer ( kind = 4 ) itop2
  logical more

  ierror = 0

  more = .false.
  idtop = 0
  idbot = 1

  do

    call perm_next ( n, iarray, more, even )

    if ( even ) then
      itop = 1
    else
      itop = -1
    end if

    ibot = 1

    do i = 1, n

      itop1 = itop
      ibot1 = ibot
      itop2 = iatop(i,iarray(i))
      ibot2 = iabot(i,iarray(i))

      call rat_mul ( itop1, ibot1, itop2, ibot2, itop, ibot, ierror )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RATMAT_DET - Fatal error!'
        write ( *, '(a)' ) '  An overflow occurred.'
        write ( *, '(a)' ) '  The determinant calculation cannot be done'
        write ( *, '(a)' ) '  for this matrix.'
        idtop = 0
        idbot = 1
        stop
      end if

    end do

    itop1 = itop
    ibot1 = ibot

    itop2 = idtop
    ibot2 = idbot

    call rat_add ( itop1, ibot1, itop2, ibot2, itop, ibot, ierror )

    if ( ierror == 0 ) then
      idtop = itop
      idbot = ibot
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RATMAT_DET - Fatal error!'
      write ( *, '(a)' ) '  An overflow occurred.'
      write ( *, '(a)' ) '  The determinant calculation cannot be done'
      write ( *, '(a)' ) '  for this matrix.'
      idtop = 0
      idbot = 1
      stop
    end if

    if ( .not. more ) then
      exit
    end if

  end do
!
!  The bottom should be positive.
!
  if ( idbot < 0 ) then
    idbot = -idbot
    idtop = -idtop
  end if

  return
end
subroutine ratmat_print ( m, n, a, b, title )

!*****************************************************************************80
!
!! RATMAT_PRINT prints out rational vectors or matrices.
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the matrix.
!
!    Input, integer ( kind = 4 ) A(M,N), B(M,N), the current rational or
!    decimal matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) b(m,n)
  character ( len = 10 ) chrtmp2
  character ( len = 10 ) chrtmp3
  character ( len = 40 ) format1
  character ( len = 40 ) format2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ione
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) jmin
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ), parameter :: ncolum = 80
  integer ( kind = 4 ) none
  integer ( kind = 4 ) npline
  character ( len = 100 ) output
  character ( len = * ) title
!
!  Figure out how many rationals we can get in NCOLUM columns.
!
  kmax = 3

  do i = 1, m
    do j = 1, n

      itemp = abs ( a(i,j) )

      do while ( 10**(kmax-2) <= itemp )
        kmax = kmax + 1
      end do

      itemp = abs ( b(i,j) )

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

  do jmin = 1, n, npline

    jmax = min ( jmin + npline - 1, n )

    write ( *, '(a)' ) ' '

    if ( jmin == 1 ) then
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
    end if

    if ( 1 < jmin .or. jmax < n ) then
      write ( output, * ) 'Columns ', jmin, ' to ', jmax
      call s_blanks_delete ( output )
      write ( *, '(a)' ) trim ( output )
      write ( *, '(a)' ) ' '
    end if

    do i = 1, m

      write ( *, format1 ) a(i,jmin:jmax)
      write ( output, format1 ) b(i,jmin:jmax)
!
!  Delete each denominator that is 1.  If all are 1, don't
!  even print out the line.
!
      none = 0

      do j = jmin, jmax

        if ( b(i,j) == 1 ) then
          ione = ( j - jmin + 1 ) * kmax
          output(ione:ione) = ' '
        else
          none = 1
        end if

      end do

      write ( *, '(a)' ) trim ( output )

      if ( jmax == n .and. i == m ) then
      else
        write ( *, '(a)' ) ' '
      end if

    end do

  end do

  return
end
subroutine regro_next ( n, v, vmax, done )

!*****************************************************************************80
!
!! REGRO_NEXT computes restricted growth functions one at a time.
!
!  Discussion:
!
!    A restricted growth function on N is a vector (V(1), ..., V(N) )
!    of values V(I) between 1 and N, satisfying the requirements:
!      V(1) = 1;
!      V(I) <= 1 + max ( V(1), V(2), ..., V(I-1) ).
!
!    The number of restricted growth functions on N is equal to
!    the Bell number B(N).
!
!    There is a bijection between restricted growth functions on N
!    and set partitions of N.
!
!  Example:
!
!    The 15 restricted growth functions for N = 4 are:
!
!    (1111), (1112), (1121), (1122), (1123),
!    (1211), (1212), (1213), (1221), (1222),
!    (1223), (1231), (1232), (1233), (1234).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer, 1986,
!    ISBN: 0387963472,
!    LC: QA164.S79.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components in the restricted
!    growth function.
!
!    Input/output, integer ( kind = 4 ) V(N).  The user need not set this
!    quantity before the initial call, and should not alter it between
!    successive calls.  On each return from the routine, with DONE = FALSE,
!    V will contain the componentwise values of the next restricted
!    growth function.
!
!    Input/output, integer ( kind = 4 ) VMAX(N).  The user need not set this
!    quantity before the initial call, and should not alter it between calls.
!    VMAX(I) records the largest value that component V(I) could take,
!    given the values of components 1 through I-1.
!
!    Input/output, logical DONE.
!    On first call, set DONE to TRUE, and then do not alter it.
!    On output, DONE will be FALSE if the routine has computed another
!    restricted growth function, or TRUE if all the restricted
!    growth functions have been returned.
!
  implicit none

  integer ( kind = 4 ) n

  logical done
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) v(n)
  integer ( kind = 4 ) vmax(n)
!
!  First call:
!
  if ( done ) then

    v(1:n) = 1

    vmax(1) = 1
    vmax(2:n) = 2

    done = .false.
!
!  Later calls.
!
  else

    j = n

    do

      if ( j == 1 ) then
        done = .true.
        return
      end if

      if ( v(j) /= vmax(j) ) then
        exit
      end if

      j = j - 1

    end do

    v(j) = v(j) + 1

    do i = j+1, n

      v(i) = 1

      if ( v(j) == vmax(j) ) then
        vmax(i) = vmax(j) + 1
      else
        vmax(i) = vmax(j)
      end if

    end do

  end if

  return
end
subroutine rfrac_to_cfrac ( m, p, q, t, ierror )

!*****************************************************************************80
!
!! RFRAC_TO_CFRAC converts rational polynomial fractions to continued fractions.
!
!  Discussion:
!
!    That is, it accepts
!
!      P(1) + P(2) * X + ... + P(M) * X**(M-1)
!      -------------------------------------------------------
!      Q(1) + Q(2) * X + ... + Q(M) * X**(M-1) + Q(M+1) * X**M
!
!    and returns the equivalent continued fraction:
!
!      1 / ( T(1) + X / ( T(2) + X / (...T(2*M-1) + X / ( T(2*M) ... )))
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2000
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi,
!    John Rice, Henry Thatcher, Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Author:
!
!    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson,
!    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, defines the number of P coefficients,
!    and is one less than the number of Q coefficients, and one
!    half the number of T coefficients.
!
!    Input, real ( kind = 8 ) P(M), Q(M+1), the coefficients defining
!    the rational polynomial fraction.
!
!    Output, real ( kind = 8 ) T(2*M), the coefficients defining the
!    continued fraction.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error;
!    nonzero, the algorithm broke down at some point with a zero divisor.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m+1,2*m+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) k
  real ( kind = 8 ) p(m)
  real ( kind = 8 ) q(m+1)
  real ( kind = 8 ) t(2*m)
  real ( kind = 8 ) ta

  ierror = 0

  a(1:m+1,1) = q(1:m+1)
  a(1:m,  2) = p(1:m)

  t(1) = a(1,1) / a(1,2)
  ta = a(m+1,1)

  do i = 1, m
    a(m-i+1,2*i+1) = ta
  end do

  do k = 1, 2*m-2

    do i = 1, (2*m-k)/2
      a(i,k+2) = a(i+1,k) - t(k) * a(i+1,k+1)
    end do

    if ( a(1,k+2) == 0.0D+00 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RFRAC_TO_CFRAC - Fatal error!'
      write ( *, '(a,i8)' ) '  A(1,K+2) is zero for K = ', k
      stop
    end if

    t(k+1) = a(1,k+1) / a(1,k+2)

  end do

  t(2*m) = a(1,2*m) / a(1,2*m+1)

  return
end
subroutine rfrac_to_jfrac ( m, p, q, r, s )

!*****************************************************************************80
!
!! RFRAC_TO_JFRAC converts a rational polynomial fraction to a J fraction.
!
!  Discussion:
!
!    The routine accepts
!
!    P(1) + P(2) * X + ... + P(M) * X**(M-1)
!    -------------------------------------------------------
!    Q(1) + Q(2) * X + ... + Q(M) * X**(M-1) + Q(M+1) * X**M
!
!    and returns the equivalent J-fraction:
!
!    R(1) / ( X + S(1) +
!    R(2) / ( X + S(2) +
!    R(3) / ...        +
!    R(M) / ( X + S(M) )... ))
!
!    Thanks to Henry Amuasi for noticing and correcting an error in a
!    previous formulation of this routine, 02 October 2010.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2010
!
!  Author:
!
!    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson,
!    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi,
!    John Rice, Henry Thatcher, Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, defines the number of P, R, and S
!    coefficients, and is one less than the number of Q coefficients.
!    1 <= M.
!
!    Input, real ( kind = 8 ) P(M), Q(M+1), the coefficients defining
!    the rational polynomial fraction.
!
!    Output, real ( kind = 8 ) R(M), S(M), the coefficients defining the
!    J-fraction.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m+1,m+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) p(m)
  real ( kind = 8 ) q(m+1)
  real ( kind = 8 ) r(m)
  real ( kind = 8 ) s(m)

  if ( m < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RFRAC_TO_JFRAC - Fatal error!'
    write ( *, '(a)' ) '  Input M < 1.'
    stop
  end if

  a(1:m+1,1) = q(1:m+1)
  a(1:m,  2) = p(1:m)

  if ( 1 < m ) then

    r(1) = a(m,2) / a(m+1,1)
    s(1) = ( r(1) * a(m,1) - a(m-1,2) ) / a(m,2)

    do k = 1, m - 2

      a(1,k+2) = r(k) * a(1,k) - s(k) * a(1,k+1)
      a(2:m-k,k+2) = r(k) * a(2:m-k,k) - a(1:m-k-1,k+1) - s(k) * a(2:m-k,k+1)

      if ( a(m-k,k+2) == 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RFRAC_TO_JFRAC - Fatal error!'
        write ( *, '(a,i8)' ) '  A(M-K,K+2) = 0 for K=', k
        stop
      end if

      r(k+1) = a(m-k,k+2) / a(m-k+1,k+1)
      s(k+1) = ( r(k+1) * a(m-k,k+1) - a(m-k-1,k+2) ) / a(m-k,k+2)

    end do

    a(1,m+1) = r(m-1) * a(1,m-1) - s(m-1) * a(1,m)

  end if

  r(m) = a(1,m+1) / a(2,m)
  s(m) = a(1,m) / a(2,m)

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
subroutine s_blanks_delete ( s )

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
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character newchr
  character oldchr
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  j = 0
  newchr = ' '

  do i = 1, len ( s )

    oldchr = newchr
    newchr = s(i:i)

    if ( newchr == TAB ) then
      newchr = ' '
    end if

    s(i:i) = ' '

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
subroutine schroeder ( n, s )

!*****************************************************************************80
!
!! SCHROEDER generates the Schroeder numbers.
!
!  Discussion:
!
!    The Schroeder number S(N) counts the number of ways to insert
!    parentheses into an expression of N items, with two or more items within
!    a parenthesis.
!
!    Note that the Catalan number C(N) counts the number of ways
!    to legally arrange a set of N left and N right parentheses.
!
!  Example:
!
!    N = 4
!
!    1234
!    12(34)
!    1(234)
!    1(2(34))
!    1(23)4
!    1((23)4)
!    (123)4
!    (12)34
!    (12)(34)
!    (1(23))4
!    ((12)3)4
!
!  First Values:
!
!           1
!           1
!           3
!          11
!          45
!         197
!         903
!        4279
!       20793
!      103049
!      518859
!     2646723
!    13648869
!    71039373
!
!  Formula:
!
!    S(N) = ( P(N)(3.0) - 3 P(N-1)(3.0) ) / ( 4 * ( N - 1 ) )
!    where P(N)(X) is the N-th Legendre polynomial.
!
!  Recursion:
!
!    S(1) = 1
!    S(2) = 1
!    S(N) = ( ( 6 * N - 9 ) * S(N-1) - ( N - 3 ) * S(N-2) ) / N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    RP Stanley,
!    Hipparchus, Plutarch, Schroeder, and Hough,
!    American Mathematical Monthly,
!    Volume 104, Number 4, 1997, pages 344-350.
!
!    Laurent Habsieger, Maxim Kazarian, Sergei Lando,
!    On the Second Number of Plutarch,
!    American Mathematical Monthly,
!    May 1998, page 446.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of Schroeder numbers desired.
!
!    Output, integer ( kind = 4 ) S(N), the Schroeder numbers.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) s(n)

  if ( n <= 0 ) then
    return
  end if

  s(1) = 1

  if ( n <= 1 ) then
    return
  end if

  s(2) = 1

  if ( n <= 2 ) then
    return
  end if

  do i = 3, n
    s(i) = ( ( 6 * i - 9 ) * s(i-1) - ( i - 3 ) * s(i-2) ) / i
  end do

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
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
!    05 February 2004
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
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
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
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements I
!    and J. (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

  return
end
subroutine subcomp_next ( n, k, a, more, h, t )

!*****************************************************************************80
!
!! SUBCOMP_NEXT computes the next subcomposition of N into K parts.
!
!  Discussion:
!
!    A composition of the integer N into K parts is an ordered sequence
!    of K nonnegative integers which sum to a value of N.
!
!    A subcomposition of the integer N into K parts is a composition
!    of M into K parts, where 0 <= M <= N.
!
!    A subcomposition of the integer N into K parts is also a lattice
!    point in the simplex whose vertices are the origin, and the K direction
!    vectors N*E(I) for I = 1 to K.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer whose subcompositions
!    are desired.
!
!    Input, integer ( kind = 4 ) K, the number of parts in the subcomposition.
!
!    Input/output, integer ( kind = 4 ) A(K), the parts of the subcomposition.
!
!    Input/output, logical MORE, set by the user to start the computation,
!    and by the routine to terminate it.
!
!    Input/output, integer H, T, two internal parameters needed for the
!    computation.  The user should allocate space for these in the calling
!    program, include them in the calling sequence, but never alter them!
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) h
  logical more
  logical, save :: more2 = .false.
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n2 = 0
  integer ( kind = 4 ) t
!
!  The first computation.
!
  if ( .not. more ) then

    n2 = 0
    a(1:k) = 0
    more2 = .false.
    h = 0
    t = 0

    more = .true.
!
!  Do the next element at the current value of N.
!
  else if ( more2 ) then

    call comp_next ( n2, k, a, more2, h, t )

  else

    more2 = .false.
    n2 = n2 + 1

    call comp_next ( n2, k, a, more2, h, t )

  end if
!
!  Termination occurs if MORE2 = FALSE and N2 = N.
!
  if ( .not. more2 .and. n2 == n ) then
    more = .false.
  end if

  return
end
subroutine subcompnz_next ( n, k, a, more )

!*****************************************************************************80
!
!! SUBCOMPNZ_NEXT computes the next subcomposition of N into K nonzero parts.
!
!  Discussion:
!
!    A composition of the integer N into K nonzero parts is an ordered sequence
!    of K positive integers which sum to a value of N.
!
!    A subcomposition of the integer N into K nonzero parts is a composition
!    of M into K nonzero parts, where 0 < M <= N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 December 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer whose subcompositions are
!    desired.
!
!    Input, integer ( kind = 4 ) K, the number of parts in the subcomposition.
!    K must be no greater than N.
!
!    Input/output, integer ( kind = 4 ) A(K), the parts of the subcomposition.
!
!    Input/output, logical MORE, set by the user to start the computation,
!    and by the routine to terminate it.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  logical more
  logical, save :: more2 = .false.
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n2 = 0

  if ( n < k ) then
    more = .false.
    a(1:k) = -1
    return
  end if
!
!  The first computation.
!
  if ( .not. more ) then

    more = .true.

    a(1:k) = 1
    n2 = k
    more2 = .false.
!
!  Do the next element at the current value of N.
!
  else if ( more2 ) then

    call compnz_next ( n2, k, a, more2 )

  else

    more2 = .false.
    n2 = n2 + 1

    call compnz_next ( n2, k, a, more2 )

  end if
!
!  Termination occurs if MORE2 = FALSE and N2 = N.
!
  if ( .not. more2 .and. n2 == n ) then
    more = .false.
  end if

  return
end
subroutine subcompnz2_next ( n_lo, n_hi, k, a, more )

!*****************************************************************************80
!
!! SUBCOMPNZ2_NEXT computes the next subcomposition of N into K nonzero parts.
!
!  Discussion:
!
!    A composition of the integer N into K nonzero parts is an ordered sequence
!    of K positive integers which sum to a value of N.
!
!    A subcomposition of the integer N into K nonzero parts is a composition
!    of M into K nonzero parts, where 0 < M <= N.
!
!    This routine computes all compositions of K into nonzero parts which sum
!    to values between N_LO and N_HI.
!
!    The routine SUBCOMPNZ_NEXT can be regarded as a special case where
!    N_LO = K.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_LO, N_HI, the range of values N for which
!    compositions are desired.
!
!    Input, integer ( kind = 4 ) K, the number of parts in the subcomposition.
!    K must be no greater than N_HI.
!
!    Input/output, integer ( kind = 4 ) A(K), the parts of the subcomposition.
!
!    Input/output, logical MORE, set by the user to start the computation,
!    and by the routine to terminate it.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  logical more
  logical, save :: more2 = .false.
  integer ( kind = 4 ) n_hi
  integer ( kind = 4 ) n_lo
  integer ( kind = 4 ), save :: n2 = 0

  if ( n_hi < k ) then
    more = .false.
    a(1:k) = -1
    return
  end if

  if ( n_hi < n_lo ) then
    more = .false.
    a(1:k) = -1
    return
  end if
!
!  The first computation.
!
  if ( .not. more ) then

    more = .true.

    n2 = max ( k, n_lo )
    more2 = .false.

    call compnz_next ( n2, k, a, more2 )
!
!  Do the next element at the current value of N.
!
  else if ( more2 ) then

    call compnz_next ( n2, k, a, more2 )

  else

    n2 = n2 + 1

    call compnz_next ( n2, k, a, more2 )

  end if
!
!  Termination occurs if MORE2 = FALSE and N2 = N_HI.
!
  if ( .not. more2 .and. n2 == n_hi ) then
    more = .false.
  end if

  return
end
subroutine subset_by_size_next ( n, a, size, more )

!*****************************************************************************80
!
!! SUBSET_BY_SIZE_NEXT returns all subsets of an N set, in order of size.
!
!  Example:
!
!    N = 4:
!
!    1 2 3 4
!    1 2 3
!    1 2 4
!    1 3 4
!    1 3
!    1 4
!    2 3
!    1
!    2
!    3
!    (the empty set)
!
!  Discussion:
!
!    The subsets are returned in decreasing order of size, with the
!    empty set last.
!
!    For a given size K, the K subsets are returned in lexicographic order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the set.
!
!    Input/output, integer ( kind = 4 ) A(N).  The entries A(1:SIZE) contain
!    the elements of the subset.  The elements are given in ascending
!    order.
!
!    Input/output, integer ( kind = 4 ) SIZE, the number of elements in the
!    subset.
!
!    Input/output, logical MORE.  Set MORE = FALSE before first call
!    for a new sequence of subsets.  It then is set and remains
!    TRUE as long as the subset computed on this call is not the
!    final one.  When the final subset is computed, MORE is set to
!    FALSE as a signal that the computation is done.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  logical more
  logical, save :: more2 = .false.
  integer ( kind = 4 ) size

  if ( .not. more ) then
    more = .true.
    more2 = .false.
    size = n
  else if ( .not. more2 ) then
    size = size - 1
  end if
!
!  Compute the next subset of size SIZE.
!
  if ( 0 < size ) then
    call ksub_next ( n, size, a, more2 )
  else if ( size == 0 ) then
    more = .false.
  end if

  return
end
subroutine subset_gray_next ( n, a, more, ncard, iadd )

!*****************************************************************************80
!
!! SUBSET_GRAY_NEXT generates all subsets of a set of order N, one at a time.
!
!  Discussion:
!
!    It generates the subsets one at a time, by adding or subtracting
!    exactly one element on each step.
!
!    This uses a Gray code ordering of the subsets.
!
!    The user should set MORE = FALSE and the value of N before
!    the first call.  On return, the user may examine A which contains
!    the definition of the new subset, and must check MORE, because
!    as soon as it is FALSE on return, all the subsets have been
!    generated and the user probably should cease calling.
!
!    The first set returned is the empty set.
!
!  Example:
!
!    N = 4
!
!    0 0 0 0
!    1 0 0 0
!    1 1 0 0
!    0 1 0 0
!    0 1 1 0
!    1 1 1 0
!    1 0 1 0
!    0 0 1 0
!    0 0 1 1
!    1 0 1 1
!    1 1 1 1
!    0 1 1 1
!    0 1 0 1
!    1 1 0 1
!    1 0 0 1
!    0 0 0 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2003
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
!    Input, integer ( kind = 4 ) N, the order of the total set from which
!    subsets will be drawn.
!
!    Input/output, integer ( kind = 4 ) A(N).  On each return, the Gray code
!    for the newly generated subset.  A(I) = 0 if element I is in the subset,
!    1 otherwise.
!
!    Input/output, logical MORE.  Set this variable FALSE before
!    the first call.  Normally, MORE will be returned TRUE but once
!    all the subsets have been generated, MORE will be
!    reset FALSE on return and you should stop calling the program.
!
!    Input/output, integer ( kind = 4 ) NCARD, the cardinality of the set
!    returned, which may be any value between 0 (the empty set) and N (the
!    whole set).
!
!    Output, integer ( kind = 4 ) IADD, the element which was added or removed
!    to the previous subset to generate the current one.  Exception:
!    the empty set is returned on the first call, and IADD is set to 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) iadd
  logical more
  integer ( kind = 4 ) ncard
!
!  The first set returned is the empty set.
!
  if ( .not. more ) then

    a(1:n) = 0

    iadd = 0
    ncard = 0
    more = .true.

  else

    iadd = 1

    if ( mod ( ncard, 2 ) /= 0 ) then

      do

        iadd = iadd + 1
        if ( a(iadd-1) /= 0 ) then
          exit
        end if

      end do

    end if

    a(iadd) = 1 - a(iadd)
    ncard = ncard + 2 * a(iadd) - 1
!
!  The last set returned is the singleton A(N).
!
    if ( ncard == a(n) ) then
      more = .false.
    end if

  end if

  return
end
subroutine subset_gray_rank ( n, a, rank )

!*****************************************************************************80
!
!! SUBSET_GRAY_RANK ranks a subset of an N set, using the Gray code ordering.
!
!  Example:
!
!    N = 4
!
!       A       Rank
!    -------   -----
!
!    0 0 0 0       1
!    0 0 0 1       2
!    0 0 1 1       3
!    0 0 1 0       4
!    0 1 1 0       5
!    0 1 1 1       6
!    0 1 0 1       7
!    0 1 0 0       8
!    1 1 0 0       9
!    1 1 0 1      10
!    1 1 1 1      11
!    1 1 1 0      12
!    1 0 1 0      13
!    1 0 1 1      14
!    1 0 0 1      15
!    1 0 0 0      16
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the total set from which
!    subsets will be drawn.
!
!    Input, integer ( kind = 4 ) A(N); A(I) is 1 if element I is in the set,
!    and 0 otherwise.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the subset in the Gray
!    code ordering.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) gray
  integer ( kind = 4 ) rank

  call ubvec_to_ui4 ( n, a, gray )

  call gray_rank ( gray, rank )

  rank = rank + 1

  return
end
subroutine subset_gray_unrank ( rank, n, a )

!*****************************************************************************80
!
!! SUBSET_GRAY_UNRANK produces a subset of an N set of the given Gray code rank.
!
!  Example:
!
!    N = 4
!
!     Rank     A
!    -----  -------
!
!        1  0 0 0 0
!        2  0 0 0 1
!        3  0 0 1 1
!        4  0 0 1 0
!        5  0 1 1 0
!        6  0 1 1 1
!        7  0 1 0 1
!        8  0 1 0 0
!        9  1 1 0 0
!       10  1 1 0 1
!       11  1 1 1 1
!       12  1 1 1 0
!       13  1 0 1 0
!       14  1 0 1 1
!       15  1 0 0 1
!       16  1 0 0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the subset in the Gray
!    code ordering.
!
!    Input, integer ( kind = 4 ) N, the order of the total set from which
!    subsets will be drawn.
!
!    Output, integer ( kind = 4 ) A(N); A(I) is 1 if element I is in the set,
!    and 0 otherwise.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) gray
  integer ( kind = 4 ) rank

  call gray_unrank ( rank - 1, gray )

  call ui4_to_ubvec ( gray, n, a )

  return
end
subroutine subset_lex_next ( n, jmp, ndim, k, a )

!*****************************************************************************80
!
!! SUBSET_LEX_NEXT generates the subsets of a set of N elements, one at a time.
!
!  Discussion:
!
!    The subsets are generated in lexicographical order.
!
!    The routine can also be forced to generate only those subsets whose
!    size is no greater than some user-specified maximum.
!
!  Example:
!
!    N = 5, JMP = ( K == 3 )
!
!    1
!    1 2
!    1 2 3
!    1 2 4
!    1 2 5
!    1 3
!    1 3 4
!    1 3 5
!    1 4
!    1 4 5
!    1 5
!    2
!    2 3
!    2 3 4
!    2 3 5
!    2 4
!    2 4 5
!    2 5
!    3
!    3 4
!    3 4 5
!    3 5
!    4
!    4 5
!    5
!    empty set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2004
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
!    Input, integer ( kind = 4 ) N, the order of the main set from which subsets
!    are chosen.
!
!    Input, logical JMP.  In the simplest case, set JMP = FALSE for
!    a normal computation.  But to jump over supersets of the input set,
!    set JMP = TRUE.  Setting JMP = ( K == 3 ) before every new call
!    will, for example, force all the subsets returned
!    to have cardinality 3 or less.
!
!    Input, integer ( kind = 4 ) NDIM, the allowed storage for A.  If NDIM < N,
!    JMP must be used to avoid creation of a subset too large to store in A.
!
!    Input/output, integer ( kind = 4 ) K.  On first call, the user must set
!    K = 0 as a startup signal to the program.  Thereafter, the routine returns
!    the size of the computed subset in K.  On the last return,
!    the empty set is returned and K is 0, which is a signal to
!    the user that the computation is complete.
!
!    Input/output, integer ( kind = 4 ) A(NDIM).  A(I) is the I-th element of
!    the subset, listed in increasing order, with 0's in entries
!    beyond entry K.
!
  implicit none

  integer ( kind = 4 ) ndim

  integer ( kind = 4 ) a(ndim)
  integer ( kind = 4 ) is
  logical jmp
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  if ( k <= 0 ) then

    if ( jmp ) then
      return
    end if

    is = 0
    k = 1
    a(1) = 1

  else if ( a(k) /= n ) then

    is = a(k)

    if ( .not. jmp ) then
      k = k + 1
    end if

    a(k) = is + 1

  else

    k = k - 1

    if ( k /= 0 ) then
      a(k) = a(k) + 1
    end if

  end if

  return
end
subroutine subset_random ( n, seed, a )

!*****************************************************************************80
!
!! SUBSET_RANDOM selects a random subset of an N-set.
!
!  Example:
!
!    N = 4
!
!    0 0 1 1
!    0 1 0 1
!    1 1 0 1
!    0 0 1 0
!    0 0 0 1
!    1 1 0 0
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
!    Input, integer ( kind = 4 ) N, the size of the full set.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, integer ( kind = 4 ) A(N).  A vector to hold the information about
!    the set chosen.  On return, if A(I) = 1, then
!    I is in the random subset, otherwise, A(I) = 0
!    and I is not in the random subset.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) seed

  do i = 1, n
    a(i) = i4_uniform ( 0, 1, seed )
  end do

  return
end
subroutine subtriangle_next ( n, more, i1, j1, i2, j2, i3, j3 )

!*****************************************************************************80
!
!! SUBTRIANGLE_NEXT computes the next subtriangle of a triangle.
!
!  Discussion:
!
!    The three sides of a triangle have been subdivided into N segments,
!    inducing a natural subdivision of the triangle into N*N subtriangles.
!    It is desired to consider each subtriangle, one at a time, in some
!    definite order.  This routine can produce information defining each
!    of the subtriangles, one after another.
!
!    The subtriangles are described in terms of the integer coordinates
!    (I,J) of their vertices.  These coordinates both range from 0 to N,
!    with the additional restriction that I + J <= N.
!
!    The vertices of each triangle are listed in counterclockwise order.
!
!  Example:
!
!    N = 4
!
!    4  *
!       |\
!       16\
!    3  *--*
!       |14|\
!       13\15\
!    2  *--*--*
!       |\9|11|\
!       |8\10\12\
!    1  *--*--*--*
!       |\2|\4|\6|\
!       |1\|3\|5\|7\
!   0   *--*--*--*--*
!
!       0  1  2  3  4
!
!    Rank  I1 J1  I2 J2  I3 J3
!    ----  -----  -----  -----
!       1   0  0   1  0   0  1
!       2   1  1   0  1   1  0
!       3   1  0   2  0   1  1
!       4   2  1   1  1   2  0
!       5   2  0   3  0   2  1
!       6   3  1   1  1   3  0
!       7   3  0   4  0   3  1
!       8   0  1   1  1   0  2
!       9   1  2   0  2   1  1
!      10   1  1   2  1   1  2
!      11   2  2   1  2   2  1
!      12   2  1   3  1   2  2
!      13   0  2   1  2   0  3
!      14   1  3   0  3   1  2
!      15   1  2   2  2   1  3
!      16   0  3   1  3   0  4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, indicates the number of subdivisions of each
!    side of the original triangle.
!
!    Input/output, logical MORE.
!    On first call, set MORE to FALSE.  Thereafter, the output value of MORE
!    will be TRUE if there are more subtriangles that can be generated by
!    further calls.  However, if MORE is returned as FALSE, the accompanying
!    subtriangle information refers to the last subtriangle that can be
!    generated.
!
!    Input/output, integer ( kind = 4 ) I1, J1, I2, J2, I3, J3, the indices of
!    the vertices of the subtriangle.
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  logical more
  integer ( kind = 4 ) n

  if ( n <= 0 ) then
    more = .false.
    return
  end if

  if ( .not. more ) then

    i1 = 0
    j1 = 0
    i2 = 1
    j2 = 0
    i3 = 0
    j3 = 1

    if ( n == 1 ) then
      more = .false.
    else
      more = .true.
    end if
!
!  We last generated a triangle like:
!
!    2---1
!     \  |
!      \ |
!       \|
!        3
!
  else if ( i2 < i3 ) then

    i1 = i3
    j1 = j3
    i2 = i1 + 1
    j2 = j1
    i3 = i1
    j3 = j1 + 1
!
!  We last generated a triangle like
!
!    3
!    |\
!    | \
!    |  \
!    1---2
!
  else if ( i1 + 1 + j1 + 1 <= n ) then

    i1 = i1 + 1
    j1 = j1 + 1
    i2 = i1 - 1
    j2 = j1
    i3 = i1
    j3 = j1 - 1
!
!  We must be at the end of a row.
!
  else

    i1 = 0
    j1 = j1 + 1
    i2 = i1 + 1
    j2 = j1
    i3 = i1
    j3 = j1 + 1

    if ( n <= j1 + 1 ) then
      more = .false.
    end if

  end if

  return
end
subroutine thue_binary_next ( n, thue )

!*****************************************************************************80
!
!! THUE_BINARY_NEXT returns the next element in a binary Thue sequence.
!
!  Discussion:
!
!    Thue demonstrated that arbitrarily long sequences of 0's and
!    1's could be generated which had the "cubefree" property.  In
!    other words, for a given string S, there was no substring W
!    such that S contained "WWW".  In fact, a stronger result holds:
!    if "a" is the first letter of W, it is never the case that S
!    contains the substring "WWa".
!
!    In this example, the digits allowed are binary, that is, just
!    "0" and "1".  The replacement rules are:
!
!    "0" -> "01"
!    "1" -> "10"
!
!    This routine produces the next binary Thue sequence in a given series.
!    However, the input sequence must be a Thue sequence in order for
!    us to guarantee that the output sequence will also have the
!    cubic nonrepetition property.
!
!    Also, enough space must be set aside in THUE to hold the
!    output sequence.  This will always be twice the input
!    value of N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 November 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N.  On input, the length of the input
!    sequence.  On output, the length of the output sequence.
!
!    Input/output, integer ( kind = 4 ) THUE(N).  On input, the initial Thue
!    sequence, and on output, the result of applying the substitution rules
!    once.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n_out
  integer ( kind = 4 ) thue(*)
  integer ( kind = 4 ) thue_out(2*n)

  n_out = 0

  do i = 1, n

    if ( thue(i) == 0 ) then
      n_out = n_out + 1
      thue_out(n_out) = 0
      n_out = n_out + 1
      thue_out(n_out) = 1
    else if ( thue(i) == 1 ) then
      n_out = n_out + 1
      thue_out(n_out) = 1
      n_out = n_out + 1
      thue_out(n_out) = 0
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'THUE_BINARY_NEXT - Fatal error!'
      write ( *, '(a)' ) '  The input sequence contains a non-binary digit'
      write ( *, '(a,i8,a,i8)' ) '  THUE(', i, ') = ', thue(i)
      stop
    end if

  end do

  n = n_out
  thue(1:n) = thue_out(1:n)

  return
end
subroutine thue_ternary_next ( n, thue )

!*****************************************************************************80
!
!! THUE_TERNARY_NEXT returns the next element in a ternary Thue sequence.
!
!  Discussion:
!
!    Thue was interested in showing that there were arbitrarily long
!    sequences of digits which never displayed a pair of contiguous
!    repetitions of any length.  That is, there was no occurrence of
!    "00" or "1010" or "121121", anywhere in the string.  This makes
!    the string "squarefree".
!
!    To do this, he demonstrated a way to start with a single digit,
!    and to repeatedly apply a series of transformation rules to each
!    digit of the sequence, deriving nonrepeating sequences of ever
!    greater length.
!
!    In this example, the digits allowed are ternary, that is, just
!    "0", "1" and "2".  The replacement rules are:
!
!    "0" -> "12"
!    "1" -> "102"
!    "2" -> "0"
!
!    This routine produces the next Thue sequence in a given series.
!    However, the input sequence must be a Thue sequence in order for
!    us to guarantee that the output sequence will also have the
!    nonrepetition property.
!
!    Also, enough space must be set aside in THUE to hold the
!    output sequence.  This will never be more than 3 times the input
!    value of N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Brian Hayes,
!    Third Base,
!    American Scientist,
!    Volume 89, Number 6, pages 490-494, November-December 2001.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N.  On input, the length of the input
!    sequence.  On output, the length of the output sequence.
!
!    Input/output, integer ( kind = 4 ) THUE(N).  On input, the initial Thue
!    sequence, and on output, the result of applying the substitution rules
!    once.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n_out
  integer ( kind = 4 ) thue(*)
  integer ( kind = 4 ) thue_out(3*n)

  n_out = 0

  do i = 1, n

    if ( thue(i) == 0 ) then
      n_out = n_out + 1
      thue_out(n_out) = 1
      n_out = n_out + 1
      thue_out(n_out) = 2
    else if ( thue(i) == 1 ) then
      n_out = n_out + 1
      thue_out(n_out) = 1
      n_out = n_out + 1
      thue_out(n_out) = 0
      n_out = n_out + 1
      thue_out(n_out) = 2
    else if ( thue(i) == 2 ) then
      n_out = n_out + 1
      thue_out(n_out) = 0
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'THUE_TERNARY_NEXT - Fatal error!'
      write ( *, '(a)' ) '  The input sequence contains a non-ternary digit'
      write ( *, '(a,i8,a,i8)' ) '  THUE(', i, ') = ', thue(i)
      stop
    end if

  end do

  n = n_out
  thue(1:n) = thue_out(1:n)

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
subroutine timestring ( string )

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = '31 May 2001   9:45:54.872 AM'
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
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
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
  character ( len = *  ) string
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

  write ( string, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine triang ( n, zeta, p )

!*****************************************************************************80
!
!! TRIANG renumbers elements in accordance with a partial ordering.
!
!  Discussion:
!
!    TRIANG is given a partially ordered set.  The partial ordering
!    is defined by a matrix ZETA, where element I is partially less than
!    or equal to element J if and only if ZETA(I,J) = 1.
!
!    TRIANG renumbers the elements with a permutation P so that if
!    element I is partially less than element J in the partial ordering,
!    then P(I) < P(J) in the usual, numerical ordering.
!
!    In other words, the elements are relabeled so that their labels
!    reflect their ordering.  This is equivalent to relabeling the
!    matrix so that, on unscrambling it, the matrix would be upper
!    triangular.
!
!    Calling I4MAT_PERM or R8MAT_PERM with P used for both the row
!    and column permutations applied to matrix ZETA will result in
!    an upper triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 May 2003
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
!    Input, integer ( kind = 4 ) N, the number of elements in the set.
!
!    Input, integer ( kind = 4 ) ZETA(N,N), describes the partial ordering.
!    ZETA(I,J) =:
!      0, for diagonal elements (I = J), or
!         for unrelated elements, or
!         if J << I.
!      1, if I << J.
!
!    Output, integer ( kind = 4 ) P(N), a permutation of the elements that
!    reflects their partial ordering.  P(I) is the new label of element I, with
!    the property that if ZETA(I,J) = 1, that is, I << J,
!    then P(I) < P(J) (in the usual ordering).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) it
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) zeta(n,n)
!
!  Make sure ZETA represents a partially ordered set.  In other words,
!  if ZETA(I,J) = 1, then ZETA(J,I) must NOT be 1.
!
  call pord_check ( n, zeta, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANG - Fatal error!'
    write ( *, '(a)' ) '  The matrix ZETA does not represent a'
    write ( *, '(a)' ) '  partial ordering.'
    stop
  end if

  m = 0
  l = 0
  p(1:n) = 0
!
!  Find the next value of M for which P(M) is 0.
!
  do

    m = m + 1

    if ( p(m) == 0 ) then
      exit
    end if

    if ( m == n ) then
      return
    end if

  end do

  it = m + 1
  ir = m + 1

  do

    if ( ir <= n ) then

      if ( p(ir) == 0 .and. zeta(ir,m) /= 0 ) then
        p(ir) = m
        m = ir
        ir = it
      else
        ir = ir + 1
      end if

    else

      l = l + 1
      iq = p(m)
      p(m) = l

      if ( iq /= 0 ) then

        ir = m + 1
        m = iq

      else if ( m == n ) then

        exit

      else

        do

          m = m + 1

          if ( p(m) == 0 ) then
            exit
          end if

          if ( m == n ) then
            return
          end if

        end do

        it = m + 1
        ir = m + 1

      end if

    end if

  end do

  return
end
subroutine tuple_next ( m1, m2, n, rank, x )

!*****************************************************************************80
!
!! TUPLE_NEXT computes the next element of a tuple space.
!
!  Discussion:
!
!    The elements are N vectors.  Each entry is constrained to lie
!    between M1 and M2.  The elements are produced one at a time.
!    The first element is
!      (M1,M1,...,M1),
!    the second element is
!      (M1,M1,...,M1+1),
!    and the last element is
!      (M2,M2,...,M2)
!    Intermediate elements are produced in lexicographic order.
!
!  Example:
!
!    N = 2, M1 = 1, M2 = 3
!
!    INPUT        OUTPUT
!    -------      -------
!    Rank  X      Rank   X
!    ----  ---    -----  ---
!    0     * *    1      1 1
!    1     1 1    2      1 2
!    2     1 2    3      1 3
!    3     1 3    4      2 1
!    4     2 1    5      2 2
!    5     2 2    6      2 3
!    6     2 3    7      3 1
!    7     3 1    8      3 2
!    8     3 2    9      3 3
!    9     3 3    0      0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M1, M2, the minimum and maximum entries.
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input/output, integer ( kind = 4 ) RANK, counts the elements.
!    On first call, set RANK to 0.  Thereafter, the output value of RANK
!    will indicate the order of the element returned.  When there are no
!    more elements, RANK will be returned as 0.
!
!    Input/output, integer ( kind = 4 ) X(N), on input the previous tuple.
!    On output, the next tuple.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)

  if ( m2 < m1 ) then
    rank = 0
    return
  end if

  if ( rank <= 0 ) then

    x(1:n) = m1
    rank = 1

  else

    rank = rank + 1
    i = n

    do

      if ( x(i) < m2 ) then
        x(i) = x(i) + 1
        exit
      end if

      x(i) = m1

      if ( i == 1 ) then
        rank = 0
        x(1:n) = m1
        exit
      end if

      i = i - 1

    end do

  end if

  return
end
subroutine tuple_next_fast ( m, n, rank, x )

!*****************************************************************************80
!
!! TUPLE_NEXT_FAST computes the next element of a tuple space, "fast".
!
!  Discussion:
!
!    The elements are N vectors.  Each entry is constrained to lie
!    between 1 and M.  The elements are produced one at a time.
!    The first element is
!      (1,1,...,1)
!    and the last element is
!      (M,M,...,M)
!    Intermediate elements are produced in lexicographic order.
!
!    This code was written as a possibly faster version of TUPLE_NEXT.
!
!  Example:
!
!    N = 2,
!    M = 3
!
!    INPUT        OUTPUT
!    -------      -------
!    Rank          X
!    ----          ----
!   -1            -1 -1
!
!    0             1  1
!    1             1  2
!    2             1  3
!    3             2  1
!    4             2  2
!    5             2  3
!    6             3  1
!    7             3  2
!    8             3  3
!    9             1  1
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the maximum entry in any component.
!    M must be greater than 0.
!
!    Input, integer ( kind = 4 ) N, the number of components.
!    N must be greater than 0.
!
!    Input, integer ( kind = 4 ) RANK, indicates the rank of the tuple.
!    Typically, 0 <= RANK < N**M.  Values of RANK greater than
!    N**M are legal and meaningful; they are equivalent to the
!    corresponding value mod (N**M).  If RANK < 0, this indicates
!    that this is the first call for the given values of (M,N).
!    Initialization is done, and X is set to a dummy value.
!
!    Output, integer ( kind = 4 ) X(N), the next tuple, or a dummy value if
!    initialization has just been done.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), save, allocatable, dimension ( : ) :: base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)

  if ( rank < 0 ) then

    if ( m <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TUPLE_NEXT_FAST - Fatal error!'
      write ( *, '(a)' ) '  The value M <= 0 is not allowed.'
      write ( *, '(a,i8)' ) '  M = ', m
      stop
    end if

    if ( n <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TUPLE_NEXT_FAST - Fatal error!'
      write ( *, '(a)' ) '  The value N <= 0 is not allowed.'
      write ( *, '(a,i8)' ) '  N = ', n
      stop
    end if

    if ( allocated ( base ) ) then
      deallocate ( base )
    end if
    allocate ( base(1:n) )

    base(n) = 1
    do i = n-1, 1, -1
      base(i) = base(i+1) * m
    end do

    x(1:n) = -1

  else

    x(1:n) = mod ( rank / base(1:n), m ) + 1

  end if

  return
end
subroutine tuple_next_ge ( m, n, rank, x )

!*****************************************************************************80
!
!! TUPLE_NEXT_GE computes the next "nondecreasing" element of a tuple space.
!
!  Discussion:
!
!    The elements are N vectors.  Each element is constrained to lie
!    between 1 and M, and to have components that are nondecreasing.
!    That is, for an element X, and any positive RANK,
!      X(I) <= X(I+RANK)
!
!    The elements are produced one at a time.
!    The first element is
!      (1,1,...,1)
!    and the last element is
!      (M,M,...,M)
!    Intermediate elements are produced in lexicographic order.
!
!  Example:
!
!    N = 3, M = 3
!
!    RANK  X
!    ----  -----
!       1  1 1 1
!       2  1 1 2
!       3  1 1 3
!       4  1 2 2
!       5  1 2 3
!       6  1 3 3
!       7  2 2 2
!       8  2 2 3
!       9  2 3 3
!      10  3 3 3
!       0  0 0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the maximum entry.
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input/output, integer ( kind = 4 ) RANK, counts the elements.
!    On first call, set RANK to 0.  Thereafter, RANK will indicate the
!    order of the element returned.  When there are no more elements,
!    RANK will be returned as 0.
!
!    Input/output, integer ( kind = 4 ) X(N), on input the previous tuple
!    (except on the first call, when the input value of X is not needed.)
!    On output, the next tuple.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)

  if ( m < 1 ) then
    return
  end if

  if ( rank <= 0 ) then
    x(1:n) = 1
    rank = 1
    return
  end if

  do i = n, 1, -1

    if ( x(i) < m ) then
      x(i) = x(i) + 1
      x(i+1:n) = x(i)
      rank = rank + 1
      return
    end if

  end do

  rank = 0
  x(1:n) = 0

  return
end
subroutine tuple_next2 ( n, xmin, xmax, rank, x )

!*****************************************************************************80
!
!! TUPLE_NEXT2 computes the next element of an integer tuple space.
!
!  Discussion:
!
!    The elements X are N vectors.
!
!    Each entry X(I) is constrained to lie between XMIN(I) and XMAX(I).
!
!    The elements are produced one at a time.
!
!    The first element is
!      (XMIN(1), XMIN(2), ..., XMIN(N)),
!    the second is (probably)
!      (XMIN(1), XMIN(2), ..., XMIN(N)+1),
!    and the last element is
!      (XMAX(1), XMAX(2), ..., XMAX(N))
!
!    Intermediate elements are produced in a lexicographic order, with
!    the first index more important than the last, and the ordering of
!    values at a fixed index implicitly defined by the sign of
!    XMAX(I) - XMIN(I).
!
!  Example:
!
!    N = 2,
!    XMIN = (/ 1, 10 /)
!    XMAX = (/ 3,  8 /)
!
!    RANK    X
!    ----  -----
!      1   1 10
!      2   1  9
!      3   1  8
!      4   2 10
!      5   2  9
!      6   2  8
!      7   3 10
!      8   3  9
!      9   3  8
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input, integer ( kind = 4 ) XMIN(N), XMAX(N), the "minimum" and "maximum"
!    entry values.  These values are minimum and maximum only in the sense of
!    the lexicographic ordering.  In fact, XMIN(I) may be less than, equal to,
!    or greater than XMAX(I).
!
!    Input/output, integer ( kind = 4 ) RANK, the rank of the item.  On first
!    call, set RANK to 0 to start up the sequence.  On return, if RANK is zero,
!    there are no more items in the sequence.
!
!    Input/output, integer ( kind = 4 ) X(N), on input the previous tuple.
!    On output, the next tuple.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) xmin(n)
  integer ( kind = 4 ) xmax(n)

  if ( rank < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of RANK = ', rank
    stop
  end if

  if ( product ( 1 + abs ( xmax(1:n) - xmin(1:n) ) ) < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of RANK = ', rank
    stop
  end if

  if ( rank == 0 ) then
    x(1:n) = xmin(1:n)
    rank = 1
    return
  end if

  rank = rank + 1
  i = n

  do

    if ( x(i) /= xmax(i) ) then
      x(i) = x(i) + sign ( 1, xmax(i) - xmin(i) )
      exit
    end if

    x(i) = xmin(i)

    if ( i == 1 ) then
      rank = 0
      exit
    end if

    i = i - 1

  end do

  return
end
subroutine ubvec_add ( n, bvec1, bvec2, bvec3 )

!*****************************************************************************80
!
!! UBVEC_ADD adds two unsigned binary vectors.
!
!  Discussion:
!
!    A UBVEC is a vector of binary digits representing an unsigned integer.
!
!    UBVEC(N) contains the units digit, UBVEC(N-1)
!    the coefficient of 2, UBVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Example:
!
!    N = 4
!
!     UBVEC1       +  UBVEC2       =  UBVEC3
!
!    ( 0 0 0 1 )   + ( 0 0 1 1 )   = ( 0 1 0 0 )
!
!      1           +   3           =   4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the vectors.
!
!    Input, integer ( kind = 4 ) BVEC1(N), BVEC2(N), the vectors to be added.
!
!    Output, integer ( kind = 4 ) BVEC3(N), the sum of the two input vectors.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 2
  integer ( kind = 4 ) bvec1(n)
  integer ( kind = 4 ) bvec2(n)
  integer ( kind = 4 ) bvec3(n)
  integer ( kind = 4 ) i
  logical overflow

  overflow = .false.

  bvec3(1:n) = bvec1(1:n) + bvec2(1:n)

  do i = n, 1, - 1
    do while ( base <= bvec3(i) )
      bvec3(i) = bvec3(i) - base
      if ( 1 < i ) then
        bvec3(i-1) = bvec3(i-1) + 1
      else
        overflow = .true.
      end if
    end do
  end do

  return
end
subroutine ubvec_to_ui4 ( n, bvec, ui4 )

!*****************************************************************************80
!
!! UBVEC_TO_UI4 makes an unsigned integer from an unsigned binary vector.
!
!  Discussion:
!
!    A UBVEC is a vector of binary digits representing an unsigned integer.
!
!    UBVEC(N) contains the units digit, UBVEC(N-1)
!    the coefficient of 2, UBVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!  Example:
!
!    N = 4
!
!         BVEC   binary UI4
!    ----------  -----  --
!    1  2  3  4
!    ----------
!    0  0  0  1       1  1
!    0  0  1  0      10  2
!    0  0  1  1      11  3
!    0  1  0  0     100  4
!    1  0  0  1    1001  9
!    1  1  1  1    1111 15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input, integer ( kind = 4 ) BVEC(N), the binary representation.
!
!    Output, integer ( kind = 4 ) UI4, the integer.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 2
  integer ( kind = 4 ) bvec(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ui4

  ui4 = 0
  do i = 1, n
    ui4 = base * ui4 + bvec(i)
  end do

  return
end
subroutine ui4_to_ubvec ( ui4, n, bvec )

!*****************************************************************************80
!
!! UI4_TO_UBVEC makes an unsigned binary vector from an unsigned integer.
!
!  Discussion:
!
!    A UBVEC is a vector of binary digits representing an unsigned integer.
!
!    UBVEC(N) contains the units digit, UBVEC(N-1)
!    the coefficient of 2, UBVEC(N-2) the coefficient of 4 and so on,
!    so that printing the digits in order gives the binary form of the number.
!
!    To guarantee that there will be enough space for any
!    value of I, it would be necessary to set N = 32.
!
!  Example:
!
!     I       BVEC         binary
!    --  ----------------  ------
!     1  1  0  0  0  0  1       1
!     2  0  1  0  0  1  0      10
!     3  1  1  0  0  1  1      11
!     4  0  0  0  1  0  0     100
!     9  0  0  1  0  0  1    1001
!    57  1  1  1  0  1  1  110111
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) UI4, an integer to be represented.
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Output, integer ( kind = 4 ) BVEC(N), the binary representation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), parameter :: base = 2
  integer ( kind = 4 ) bvec(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ui4
  integer ( kind = 4 ) ui4_copy

  ui4_copy = ui4

  do i = n, 1, -1

    bvec(i) = mod ( ui4_copy, base )

    ui4_copy = ui4_copy / base

  end do

  return
end
subroutine vec_colex_next ( dim_num, base, a, more )

!*****************************************************************************80
!
!! VEC_COLEX_NEXT generates vectors in colex order.
!
!  Discussion:
!
!    The vectors are produced in colexical order, starting with
!    (0,0,...,0),
!    (1,0,...,0),
!    ...
!    (BASE-1,BASE-1,...,BASE-1).
!
!  Example:
!
!    DIM_NUM = 2,
!    BASE = 3
!
!    0   0
!    1   0
!    2   0
!    0   1
!    1   1
!    2   1
!    0   2
!    1   2
!    2   2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) BASE, the base to be used.  BASE = 2 will
!    give vectors of 0's and 1's, for instance.
!
!    Input/output, integer ( kind = 4 ) A(DIM_NUM), the next vector.
!
!    Input/output, logical MORE.  Set this variable FALSE before
!    the first call.  On return, MORE is TRUE if another vector has
!    been computed.  If MORE is returned FALSE, ignore the output
!    vector and stop calling the routine.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) a(dim_num)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  logical more

  if ( .not. more ) then

    a(1:dim_num) = 0
    more = .true.

  else

    do i = 1, dim_num

      a(i) = a(i) + 1

      if ( a(i) < base ) then
        return
      end if

      a(i) = 0

    end do

    more = .false.

  end if

  return
end
subroutine vec_colex_next2 ( dim_num, base, a, more )

!*****************************************************************************80
!
!! VEC_COLEX_NEXT2 generates vectors in colex order.
!
!  Discussion:
!
!    The vectors are produced in colexical order, starting with
!
!    (0,        0,        ...,0),
!    (1,        0,        ...,0),
!     ...
!    (BASE(1)-1,0,        ...,0)
!
!    (0,        1,        ...,0)
!    (1,        1,        ...,0)
!    ...
!    (BASE(1)-1,1,        ...,0)
!
!    (0,        2,        ...,0)
!    (1,        2,        ...,0)
!    ...
!    (BASE(1)-1,BASE(2)-1,...,BASE(DIM_NUM)-1).
!
!  Example:
!
!    DIM_NUM = 2,
!    BASE = ( 3, 3 )
!
!    0   0
!    1   0
!    2   0
!    0   1
!    1   1
!    2   1
!    0   2
!    1   2
!    2   2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the bases to be used in each
!    dimension.  In dimension I, entries will range from 0 to BASE(I)-1.
!
!    Input/output, integer ( kind = 4 ) A(DIM_NUM).  On each return, A
!    will contain entries in the range 0 to N-1.
!
!    Input/output, logical MORE.  Set this variable FALSE before
!    the first call.  On return, MORE is TRUE if another vector has
!    been computed.  If MORE is returned FALSE, ignore the output
!    vector and stop calling the routine.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) a(dim_num)
  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) i
  logical more

  if ( .not. more ) then

    a(1:dim_num) = 0
    more = .true.

  else

    do i = 1, dim_num

      a(i) = a(i) + 1

      if ( a(i) < base(i) ) then
        return
      end if

      a(i) = 0

    end do

    more = .false.

  end if

  return
end
subroutine vec_colex_next3 ( dim_num, base, a, more )

!*****************************************************************************80
!
!! VEC_COLEX_NEXT3 generates vectors in colex order.
!
!  Discussion:
!
!    The vectors are produced in colexical order, starting with
!
!    (1,        1,        ...,1),
!    (2,        1,        ...,1),
!     ...
!    (BASE(1),  1,        ...,1)
!
!    (1,        2,        ...,1)
!    (2,        2,        ...,1)
!    ...
!    (BASE(1),  2,        ...,0)
!
!    (1,        3,        ...,1)
!    (2,        3,        ...,1)
!    ...
!    (BASE(1),  BASE(2),  ...,BASE(DIM_NUM)).
!
!  Example:
!
!    DIM_NUM = 2,
!    BASE = ( 3, 3 )
!
!    1   1
!    2   1
!    3   1
!    1   2
!    2   2
!    3   2
!    1   3
!    2   3
!    3   3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the bases to be used in each
!    dimension.  In dimension I, entries will range from 1 to BASE(I).
!
!    Input/output, integer ( kind = 4 ) A(DIM_NUM).
!
!    Input/output, logical MORE.  Set this variable FALSE before
!    the first call.  On return, MORE is TRUE if another vector has
!    been computed.  If MORE is returned FALSE, ignore the output
!    vector and stop calling the routine.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) a(dim_num)
  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) i
  logical more

  if ( .not. more ) then

    a(1:dim_num) = 1
    more = .true.

  else

    do i = 1, dim_num

      a(i) = a(i) + 1

      if ( a(i) <= base(i) ) then
        return
      end if

      a(i) = 1

    end do

    more = .false.

  end if

  return
end
subroutine vec_gray_rank ( n, base, a, rank )

!*****************************************************************************80
!
!! VEC_GRAY_RANK computes the rank of a product space element.
!
!  Discussion:
!
!    The rank applies only to the elements as produced by the routine
!    VEC_GRAY_NEXT
!
!  Example:
!
!    N = 2, BASE = (/ 2, 3 /), A = ( 1, 2 ),
!
!    RANK = 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 May 2007
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer, 1986,
!    ISBN: 0387963472,
!    LC: QA164.S79.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input, integer ( kind = 4 ) BASE(N), contains the number of degrees of
!    freedom of each component.  The output values of A will
!    satisfy 0 <= A(I) < BASE(I).
!
!    Input, integer ( kind = 4 ) A(N), the product space element, with the
!    property that 0 <= A(I) < BASE(I) for each entry I.
!
!    Output, integer ( kind = 4 ) RANK, the rank, or order, of the element in
!    the list of all elements.  The rank count begins at 1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) base(n)
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank

  rank = 0

  do i = 1, n

    if ( mod ( rank, 2 ) == 1 ) then
      c = base(i) - a(i) - 1
    else
      c = a(i)
    end if

    rank = base(i) * rank + c

  end do

  rank = rank + 1

  return
end
subroutine vec_gray_unrank ( n, base, rank, a )

!*****************************************************************************80
!
!! VEC_GRAY_UNRANK computes the product space element of a given rank.
!
!  Discussion:
!
!    The rank applies only to the elements as produced by the routine
!    VEC_GRAY_NEXT.
!
!  Example:
!
!    N = 2, BASE = ( 2, 3 ), RANK = 4.
!
!    A = ( 1, 2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 May 2007
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer, 1986,
!    ISBN: 0387963472,
!    LC: QA164.S79.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input, integer ( kind = 4 ) BASE(N), contains the number of degrees of
!    freedom of each component.  The output values of A will
!    satisfy 0 <= A(I) < BASE(I).
!
!    Input, integer ( kind = 4 ) RANK, the desired rank, or order, of the
!    element in the list of all elements.  The rank count begins at 1 and
!    extends to RANK_MAX = Product ( 1 <= I <= N ) BASE(I).
!
!    Output, integer ( kind = 4 ) A(N), the product space element of the
!    given rank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) base(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) s

  s = rank - 1

  do i = n, 1, -1

    a(i) = mod ( s, base(i) )
    s = s / base(i)

    if ( mod ( s, 2 ) == 1 ) then
      a(i) = base(i) - a(i) - 1
    end if

  end do

  return
end
subroutine vec_gray_next ( n, base, a, done, change )

!*****************************************************************************80
!
!! VEC_GRAY_NEXT computes the elements of a product space.
!
!  Discussion:
!
!    The elements are produced one at a time.
!
!    This routine handles the case where the number of degrees of freedom may
!    differ from one component to the next.
!
!    A method similar to the Gray code is used, so that successive
!    elements returned by this routine differ by only a single element.
!
!    The routine uses internal static memory.
!
!  Example:
!
!    N = 2, BASE = ( 2, 3 ), DONE = TRUE
!
!     A    DONE  CHANGE
!    ---  -----  ------
!    0 0  FALSE    1
!    0 1  FALSE    2
!    0 2  FALSE    2
!    1 2  FALSE    1
!    1 1  FALSE    2
!    1 0  FALSE    2
!    1 0   TRUE   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer, 1986,
!    ISBN: 0387963472,
!    LC: QA164.S79.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input, integer ( kind = 4 ) BASE(N), contains the number of degrees of
!    freedom of each component.  The output values of A will
!    satisfy 0 <= A(I) < BASE(I).
!
!    Input/output, integer ( kind = 4 ) A(N).  On the first call, the input
!    value of A doesn't matter.  Thereafter, it should be the same as
!    its output value from the previous call.  On output, if DONE
!    is FALSE, then A contains the next element of the space.
!
!    Input/output, logical DONE.  On the first call, the user must
!    set DONE to TRUE.  This signals the program to initialize data.
!    On every return, if DONE is FALSE, the program has computed
!    another entry, which is contained in A.  If DONE is TRUE,
!    then there are no more entries, and the program should not be
!    called for any more.
!
!    Output, integer ( kind = 4 ) CHANGE, is set to the index of the element
!    whose value was changed.  On return from the first call, CHANGE
!    is 1, even though all the elements have been "changed".  On
!    return with DONE equal to TRUE, CHANGE is -1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ), save, allocatable, dimension ( : ) :: active
  integer ( kind = 4 ) base(n)
  integer ( kind = 4 ) change
  integer ( kind = 4 ), save, allocatable, dimension ( : ) :: dir
  logical done
  integer ( kind = 4 ) i
!
!  The user is calling for the first time.
!
  if ( done ) then

    done = .false.
    a(1:n) = 0

    if ( allocated ( active ) ) then
      deallocate ( active )
    end if

    if ( allocated ( dir ) ) then
      deallocate ( dir )
    end if

    allocate ( active(1:n) )
    allocate ( dir(1:n) )

    dir(1:n) = 1
    active(1:n) = 1

    do i = 1, n

      if ( base(i) < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'VEC_GRAY_NEXT - Warning!'
        write ( *, '(a,i8)' ) '  For index I = ',i
        write ( *, '(a,i8)' ) '  the nonpositive value of BASE(I) = ', base(i)
        write ( *, '(a)' ) '  which was reset to 1!'
        base(i) = 1
        active(i) = 0
      else if ( base(i) == 1 ) then
        active(i) = 0
      end if

    end do

    change = 1

    return

  end if
!
!  Seek the maximum active index.
!
  change = -1

  do i = n, 1, -1
    if ( active(i) == 1 ) then
      change = i
      exit
    end if
  end do
!
!  If there are NO active indices, we have generated all vectors.
!
  if ( change == -1 ) then
    done = .true.
    deallocate ( active )
    deallocate ( dir )
    return
  end if
!
!  Increment the element with maximum active index.
!
  a(change) = a(change) + dir(change)
!
!  If we attained a minimum or maximum value, reverse the direction
!  vector, and deactivate the index.
!
  if ( a(change) == 0 .or. a(change) == base(change) - 1 ) then
    dir(change) = -dir(change)
    active(change) = 0
  end if
!
!  Activate all subsequent indices.
!
  do i = change + 1, n
    if ( 1 < base(i) ) then
      active(i) = 1
    end if
  end do

  return
end
subroutine vec_lex_next ( dim_num, base, a, more )

!*****************************************************************************80
!
!! VEC_LEX_NEXT generates vectors in lex order.
!
!  Discussion:
!
!    The vectors are produced in lexical order, starting with
!    (0,0,...,0),
!    (0,0,...,1),
!    ...
!    (BASE-1,BASE-1,...,BASE-1).
!
!  Example:
!
!    DIM_NUM = 2,
!    BASE = 3
!
!    0   0
!    0   1
!    0   2
!    1   0
!    1   1
!    1   2
!    2   0
!    2   1
!    2   2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 May 2007
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the size of the vectors to be used.
!
!    Input, integer ( kind = 4 ) BASE, the base to be used.  BASE = 2 will
!    give vectors of 0's and 1's, for instance.
!
!    Input/output, integer ( kind = 4 ) A(DIM_NUM), the next vector.
!
!    Input/output, logical MORE.  Set this variable FALSE before
!    the first call.  On return, MORE is TRUE if another vector has
!    been computed.  If MORE is returned FALSE, ignore the output
!    vector and stop calling the routine.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) a(dim_num)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  logical more

  if ( .not. more ) then

    a(1:dim_num) = 0
    more = .true.

  else

    do i = dim_num, 1, -1

      a(i) = a(i) + 1

      if ( a(i) < base ) then
        return
      end if

      a(i) = 0

    end do

    more = .false.

  end if

  return
end
subroutine vec_random ( n, base, seed, a )

!*****************************************************************************80
!
!! VEC_RANDOM selects a random N-vector of integers modulo a given base.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the vector to be generated.
!
!    Input, integer ( kind = 4 ) BASE, the base to be used.
!
!    Input/output, integer ( kind = 4 ) SEED, a random number seed.
!
!    Output, integer ( kind = 4 ) A(N), a list of N random values between
!    0 and BASE-1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) seed

  do i = 1, n
    a(i) = i4_uniform ( 0, base-1, seed )
  end do

  return
end
subroutine vector_constrained_next ( n, x_min, x_max, x, constraint, more )

!*****************************************************************************80
!
!! VECTOR_CONSTRAINED_NEXT returns the "next" constrained vector.
!
!  Discussion:
!
!    We consider all vectors of dimension N satisfying:
!
!      X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
!
!    We are only interested in the subset of these vectors which
!    satisfy the following constraint:
!
!      sum ( 1 <= I <= N ) ( ( X(I) - 1 ) / X_MAX(I) ) <= 1
!
!    We can carry out this check using integer arithmetic if we
!    multiply through by P = product ( X_MAX(1:N) ):
!
!      sum ( 1 <= I <= N ) ( ( X(I) - 1 ) * ( P / X_MAX(I) ) ) <= P.
!
!    This routine returns, one at a time, and in right-to-left
!    lexicographic order, exactly those vectors which satisfy
!    the constraint.
!
!  Example:
!
!    N = 3
!    X_MIN:   2   2   1
!    X_MAX:   4   5   3
!
!    P = 60
!
!    #  X(1)  X(2)  X(3)  CONSTRAINT
!
!    1    2     2     1       27
!    2    3     2     1       42
!    3    4     2     1       57
!    4    2     3     1       39
!    5    3     3     1       54
!    6    2     4     1       51
!    7    2     2     2       47
!    8    2     3     2       59
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components in the vector.
!
!    Input, integer ( kind = 4 ) X_MIN(N), X_MAX(N), the minimum and maximum
!    values allowed in each component.
!
!    Input/output, integer ( kind = 4 ) X(N).  On first call, with
!    MORE = FALSE), the input value of X is not important.  On subsequent calls,
!    the input value of X should be the output value from the previous call.
!    On output, (with MORE = TRUE), the value of X will be the "next"
!    vector in the reverse lexicographical list of vectors that satisfy
!    the condition.  However, on output with MORE = FALSE, the vector
!    X is meaningless, because there are no more vectors in the list.
!
!    Output, integer ( kind = 4 ) CONSTRAINT, the constraint value for X.
!    Valid vectors X will have a CONSTRAINT value between product(X_MIN(1:N))
!    (automatically) and product(X_MAX(1:N)) (because we skip over vectors
!    with a constraint larger than this value).
!
!    Input/output, logical MORE.  On input, if the user has set MORE
!    FALSE, the user is requesting the initiation of a new sequence
!    of values.  If MORE is TRUE, then the user is requesting "more"
!    values in the current sequence.  On output, if MORE is TRUE,
!    then another value was found and returned in X, but if MORE is
!    FALSE, then there are no more values in the sequence, and X is
!    NOT the next value.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) constraint
  integer ( kind = 4 ) i
  logical more
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x_max(n)
  integer ( kind = 4 ) x_min(n)
  integer ( kind = 4 ), save :: x_prod

  if ( .not. more ) then

    x(1:n) = x_min(1:n)

    x_prod = product ( x_max(1:n) )

    constraint = sum ( ( x(1:n) - 1 ) * ( x_prod / x_max(1:n) ) )

    if ( x_prod < constraint ) then
      more = .false.
    else
      more = .true.
    end if

    return

  else

    i = 1

    do

      if ( x(i) < x_max(i) ) then

        x(i) = x(i) + 1

        constraint = sum ( ( x(1:n) - 1 ) * ( x_prod / x_max(1:n) ) )

        if ( constraint <= x_prod ) then
          exit
        end if

      end if

      x(i) = x_min(i)

      i = i + 1

      if ( n < i ) then
        more = .false.
        exit
      end if

    end do

  end if

  return
end
subroutine vector_constrained_next2 ( n, x_min, x_max, x, constraint, more )

!*****************************************************************************80
!
!! VECTOR_CONSTRAINED_NEXT2 returns the "next" constrained vector.
!
!  Discussion:
!
!    We consider all vectors of dimension N satisfying
!
!      X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
!
!    We are only interested in the subset of these vectors which
!    satisfy the following constraint:
!
!      sum ( 1 <= I <= N ) ( X(I) / X_MAX(I) ) <= 1
!
!    We can carry out this check using integer arithmetic if we
!    multiply through by P = product ( X_MAX(1:N) ):
!
!      sum ( 1 <= I <= N ) ( X(I) * ( P / X_MAX(I) ) ) <= P.
!
!    This routine returns, one at a time, and in right-to-left
!    lexicographic order, exactly those vectors which satisfy
!    the constraint.
!
!  Example:
!
!    N = 3
!    X_MIN:   1   1   1
!    X_MAX:   5   6   4
!
!    P = 120
!
!    #  X(1)  X(2)  X(3)  CONSTRAINT
!
!    1    1     1     1       74
!    2    2     1     1       98
!    3    1     2     1       94
!    4    2     2     1      119
!    5    1     3     1      114
!    6    1     1     2      104
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components in the vector.
!
!    Input, integer ( kind = 4 ) X_MIN(N), X_MAX(N), the minimum and maximum
!    values allowed in each component.
!
!    Input/output, integer ( kind = 4 ) X(N).  On first call, with
!    MORE = FALSE, the input value of X is not important.  On subsequent calls,
!    the input value of X should be the output value from the previous call.
!    On output, (with MORE = TRUE), the value of X will be the "next"
!    vector in the reverse lexicographical list of vectors that satisfy
!    the condition.  However, on output with MORE = FALSE, the vector
!    X is meaningless, because there are no more vectors in the list.
!
!    Output, integer ( kind = 4 ) CONSTRAINT, the constraint value for X.
!    Valid vectors X will have a CONSTRAINT value between product(X_MIN(1:N))
!    (automatically) and product(X_MAX(1:N)) (because we skip over vectors
!    with a constraint larger than this value).
!
!    Input/output, logical MORE.  On input, if the user has set MORE
!    FALSE, the user is requesting the initiation of a new sequence
!    of values.  If MORE is TRUE, then the user is requesting "more"
!    values in the current sequence.  On output, if MORE is TRUE,
!    then another value was found and returned in X, but if MORE is
!    FALSE, then there are no more values in the sequence, and X is
!    NOT the next value.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) constraint
  integer ( kind = 4 ) i
  logical more
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x_max(n)
  integer ( kind = 4 ) x_min(n)
  integer ( kind = 4 ), save :: x_prod

  if ( .not. more ) then

    x(1:n) = x_min(1:n)

    x_prod = product ( x_max(1:n) )

    constraint = sum ( x(1:n) * ( x_prod / x_max(1:n) ) )

    if ( x_prod < constraint ) then
      more = .false.
    else
      more = .true.
    end if

    return

  else

    i = 1

    do

      if ( x(i) < x_max(i) ) then

        x(i) = x(i) + 1

        constraint = sum ( x(1:n) * ( x_prod / x_max(1:n) ) )

        if ( constraint <= x_prod ) then
          exit
        end if

      end if

      x(i) = x_min(i)

      i = i + 1

      if ( n < i ) then
        more = .false.
        exit
      end if

    end do

  end if

  return
end
subroutine vector_constrained_next3 ( n, x_min, x_max, x, constraint, more )

!*****************************************************************************80
!
!! VECTOR_CONSTRAINED_NEXT3 returns the "next" constrained vector.
!
!  Discussion:
!
!    This routine addresses the same problem as VECTOR_CONSTRAINED_NEXT2,
!    and differs only in that real arithmetic is used, rather than
!    integer arithmetic.  Integer arithmetic allows us to do an exact
!    calculation, but we run into overflow problems in simple cases
!    where N is 10 and the X_MAX entries are of order 10, for instance.
!
!    We consider all vectors X of dimension N satisfying
!
!      X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
!
!    We are only interested in the subset of these vectors which
!    satisfy the following constraint:
!
!      sum ( 1 <= I <= N ) ( X(I) / X_MAX(I) ) <= 1
!
!    This routine returns, one at a time, and in right-to-left
!    lexicographic order, exactly those vectors which satisfy
!    the constraint.
!
!  Example:
!
!    N = 3
!    X_MIN:   1   1   1
!    X_MAX:   5   6   4
!
!    P = 120
!
!    #  X(1)  X(2)  X(3)  CONSTRAINT
!
!    1    1     1     1       0.62
!    2    2     1     1       0.82
!    3    1     2     1       0.78
!    4    2     2     1       0.98
!    5    1     3     1       0.95
!    6    1     1     2       0.87
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 April 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components in the vector.
!
!    Input, integer ( kind = 4 ) X_MIN(N), X_MAX(N), the minimum and maximum
!    values allowed in each component.
!
!    Input/output, integer ( kind = 4 ) X(N).  On first call, with
!    MORE = FALSE, the input value of X is not important.  On subsequent calls,
!    the input value of X should be the output value from the previous call.
!    On output, (with MORE = TRUE), the value of X will be the "next"
!    vector in the reverse lexicographical list of vectors that satisfy
!    the condition.  However, on output with MORE = FALSE, the vector
!    X is meaningless, because there are no more vectors in the list.
!
!    Output, real ( kind = 8 ) CONSTRAINT, the constraint value for X.
!    Valid vectors X will have a CONSTRAINT value between
!      product(X_MIN(1:N)) / product(X_MAX(1:N))
!    and 1.0.
!
!    Input/output, logical MORE.  On input, if the user has set MORE
!    FALSE, the user is requesting the initiation of a new sequence
!    of values.  If MORE is TRUE, then the user is requesting "more"
!    values in the current sequence.  On output, if MORE is TRUE,
!    then another value was found and returned in X, but if MORE is
!    FALSE, then there are no more values in the sequence, and X is
!    NOT the next value.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) constraint
  integer ( kind = 4 ) i
  logical more
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x_max(n)
  integer ( kind = 4 ) x_min(n)

  if ( .not. more ) then

    x(1:n) = x_min(1:n)

    constraint = sum ( real ( x(1:n), kind = 8 ) &
                     / real ( x_max(1:n), kind = 8 ) )

    if ( 1.0D+00 < constraint ) then
      more = .false.
    else
      more = .true.
    end if

    return

  else

    i = 1

    do

      if ( x(i) < x_max(i) ) then

        x(i) = x(i) + 1

        constraint = sum ( real ( x(1:n), kind = 8 ) &
                         / real ( x_max(1:n), kind = 8 ) )

        if ( constraint <= 1.0D+00 ) then
          exit
        end if

      end if

      x(i) = x_min(i)

      i = i + 1

      if ( n < i ) then
        more = .false.
        exit
      end if

    end do

  end if

  return
end
subroutine vector_constrained_next4 ( n, alpha, x_min, x_max, x, q, more )

!*****************************************************************************80
!
!! VECTOR_CONSTRAINED_NEXT4 returns the "next" constrained vector.
!
!  Discussion:
!
!    This routine is similar to VECTOR_CONSTRAINED_NEXT2 and
!    VECTOR_CONSTRAINED_NEXT3.
!
!    We consider all vectors X of dimension N satisfying
!
!      X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
!
!    We are only interested in the subset of these vectors which
!    satisfy the following constraint:
!
!      sum ( 1 <= I <= N ) ALPHA(I) * X(I) <= Q
!
!    This routine returns, one at a time, and in right-to-left
!    lexicographic order, exactly those vectors which satisfy
!    the constraint.
!
!  Example:
!
!    N = 3
!    ALPHA    4.0  3.0  5.0
!    Q       20.0
!    X_MIN:   1   1   1
!    X_MAX:   5   6   4
!
!    #  X(1)  X(2)  X(3)     Total
!
!    1    1     1     1       12.0
!    2    2     1     1       16.0
!    3    3     1     1       20.0
!    4    1     2     1       15.0
!    5    2     2     1       19.0
!    6    1     3     1       18.0
!    7    1     1     2       17.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components in the vector.
!
!    Input, real ( kind = 8 ) ALPHA(N), the coefficient vector.
!
!    Input, integer ( kind = 4 ) X_MIN(N), X_MAX(N), the minimum and maximum
!    values allowed in each component.
!
!    Input/output, integer ( kind = 4 ) X(N).  On first call, with
!    MORE = FALSE, the input value of X is not important.  On subsequent calls,
!    the input value of X should be the output value from the previous call.
!    On output, (with MORE = TRUE), the value of X will be the "next"
!    vector in the reverse lexicographical list of vectors that satisfy
!    the condition.  However, on output with MORE = FALSE, the vector
!    X is meaningless, because there are no more vectors in the list.
!
!    Input, real ( kind = 8 ) Q, the limit on the sum.
!
!    Input/output, logical MORE.  On input, if the user has set MORE
!    FALSE, the user is requesting the initiation of a new sequence
!    of values.  If MORE is TRUE, then the user is requesting "more"
!    values in the current sequence.  On output, if MORE is TRUE,
!    then another value was found and returned in X, but if MORE is
!    FALSE, then there are no more values in the sequence, and X is
!    NOT the next value.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha(n)
  integer ( kind = 4 ) i
  logical more
  real ( kind = 8 ) q
  real ( kind = 8 ) total
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x_max(n)
  integer ( kind = 4 ) x_min(n)

  if ( .not. more ) then

    x(1:n) = x_min(1:n)

    total = dot_product ( alpha(1:n), real ( x(1:n), kind = 8 ) )

    if ( q < total ) then
      more = .false.
    else
      more = .true.
    end if

    return

  else

    i = 1

    do

      if ( x(i) < x_max(i) ) then

        x(i) = x(i) + 1

        total = dot_product ( alpha(1:n), real ( x(1:n), kind = 8 ) )

        if ( total <= q ) then
          exit
        end if

      end if

      x(i) = x_min(i)

      i = i + 1

      if ( n < i ) then
        more = .false.
        exit
      end if

    end do

  end if

  return
end
subroutine vector_constrained_next5 ( n, x, sum_min, sum_max, more )

!*****************************************************************************80
!
!! VECTOR_CONSTRAINED_NEXT5 returns the "next" constrained vector.
!
!  Discussion:
!
!    We consider all positive integer vectors X dimension N satisfying:
!
!      SUM_MIN <= X(1:N) <= SUM_MAX.
!
!    This routine returns, one at a time, and in right-to-left
!    lexicographic order, exactly those vectors which satisfy
!    the constraint.
!
!  Example:
!
!    N = 3
!    SUM_MIN = 5
!    SUM_MAX = 6
!
!    #  X(1)  X(2)  X(3)     SUM
!
!    1    3     1     1        5
!    2    2     2     1        5
!    3    2     1     2        5
!    4    1     3     1        5
!    5    1     2     2        5
!    6    1     1     3        5
!
!    7    4     1     1        6
!    8    3     2     1        6
!    9    3     1     2        6
!   10    2     3     1        6
!   11    2     2     2        6
!   12    2     1     3        6
!   13    1     4     1        6
!   14    1     3     2        6
!   15    1     2     3        6
!   16    1     1     4        6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components in the vector.
!
!    Input, integer ( kind = 4 ) SUM_MIN, SUM_MAX, the minimum and maximum sums.
!
!    Input/output, integer ( kind = 4 ) X(N).  On first call, with
!    MORE = FALSE, the input value of X is not important.  On subsequent calls,
!    the input value of X should be the output value from the previous call.
!    On output, (with MORE = TRUE), the value of X will be the "next"
!    vector in the reverse lexicographical list of vectors that satisfy
!    the condition.  However, on output with MORE = FALSE, the vector
!    X is meaningless, because there are no more vectors in the list.
!
!    Input/output, logical MORE.  On input, if the user has set MORE
!    FALSE, the user is requesting the initiation of a new sequence
!    of values.  If MORE is TRUE, then the user is requesting "more"
!    values in the current sequence.  On output, if MORE is TRUE,
!    then another value was found and returned in X, but if MORE is
!    FALSE, then there are no more values in the sequence, and X is
!    NOT the next value.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), save :: base = 0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  logical more
  integer ( kind = 4 ) sum_max
  integer ( kind = 4 ) sum_min
  integer ( kind = 4 ) x(n)
!
!  Initialization.
!
  if ( .not. more ) then

    if ( sum_max < n ) then
      more = .false.
      return
    end if

    if ( sum_max < sum_min ) then
      more = .false.
      return
    end if

    more = .true.

    base = max ( sum_min, n )

    x(1) = base - n + 1
    x(2:n) = 1

    return
!
!  Next element.
!
  else
!
!  Search from the right, seeking an index I < N for which 1 < X(I).
!
    do i = n-1, 1, -1
!
!  If you find such an I, decrease X(I) by 1, and add that to X(I+1).
!
      if ( 1 < x(i) ) then

        x(i)   = x(i)   - 1
        x(i+1) = x(i+1) + 1
!
!  Now grab all the "excess" 1's from the entries to the right of X(I+1).
!
        do j = i+2, n
          if ( 1 < x(j) ) then
            x(i+1) = x(i+1) + x(j) - 1
            x(j) = 1
          end if
        end do

        return

      end if

    end do
!
!  The current vector is (1,1,1,...,BASE-N+1).
!  If BASE < SUM_MAX, then increase BASE by 1, and start the new series.
!
    if ( base < sum_max ) then
      base = base + 1
      x(1) = base - n + 1
      x(2:n) = 1
      return
    end if
!
!  We returned the last legal vector on the previous call.
!  The calculation is done.
!
    more = .false.

  end if

  return
end
subroutine vector_constrained_next6 ( n, alpha, x_min, x_max, x, q_min, &
  q_max, more )

!*****************************************************************************80
!
!! VECTOR_CONSTRAINED_NEXT6 returns the "next" constrained vector.
!
!  Discussion:
!
!    We consider vectors X of dimension N satisfying:
!
!      X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
!
!    We are only interested in the subset of these vectors which
!    satisfy the following constraint:
!
!      Q_MIN <= sum ( 1 <= I <= N ) ALPHA(I) * X(I) <= Q_MAX
!
!    This routine returns, one at a time, and in right-to-left
!    lexicographic order, exactly those vectors which satisfy
!    the constraint.
!
!  Example:
!
!    N = 3
!    ALPHA    4.0  3.0  5.0
!    Q_MIN   16.0
!    Q_MAX   20.0
!    X_MIN:   1   1   1
!    X_MAX:   5   6   4
!
!    #  X(1)  X(2)  X(3)     Total
!
!    1    2     1     1       16.0
!    2    3     1     1       20.0
!    3    2     2     1       19.0
!    4    1     3     1       18.0
!    5    1     1     2       17.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components in the vector.
!
!    Input, real ( kind = 8 ) ALPHA(N), the coefficient vector.
!
!    Input, integer ( kind = 4 ) X_MIN(N), X_MAX(N), the minimum and maximum
!    values allowed in each component.
!
!    Input/output, integer ( kind = 4 ) X(N).  On first call, with
!    MORE = FALSE, the input value of X is not important.  On subsequent calls,
!    the input value of X should be the output value from the previous call.
!    On output, (with MORE = TRUE), the value of X will be the "next"
!    vector in the reverse lexicographical list of vectors that satisfy
!    the condition.  However, on output with MORE = FALSE, the vector
!    X is meaningless, because there are no more vectors in the list.
!
!    Input, real ( kind = 8 ) Q_MIN, Q_MAX, the lower and upper
!    limits on the sum.
!
!    Input/output, logical MORE.  On input, if the user has set MORE
!    FALSE, the user is requesting the initiation of a new sequence
!    of values.  If MORE is TRUE, then the user is requesting "more"
!    values in the current sequence.  On output, if MORE is TRUE,
!    then another value was found and returned in X, but if MORE is
!    FALSE, then there are no more values in the sequence, and X is
!    NOT the next value.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha(n)
  integer ( kind = 4 ) i
  logical more
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  real ( kind = 8 ) total
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x_max(n)
  integer ( kind = 4 ) x_min(n)

  if ( .not. more ) then

    more = .true.

    x(1:n) = x_min(1:n)

    total = dot_product ( alpha(1:n), real ( x(1:n), kind = 8 ) )

    if ( q_min <= total .and. total <= q_max ) then
      return
    end if

  end if

  do

    i = n

    do

      if ( x(i) < x_max(i) ) then
        exit
      end if

      if ( i <= 1 ) then
        more = .false.
        return
      end if

      i = i - 1

    end do

    x(i) = x(i) + 1
    x(i+1:n) = x_min(i+1:n)

    total = dot_product ( alpha(1:n), real ( x(1:n), kind = 8 ) )

    if ( q_min <= total .and. total <= q_max ) then
      exit
    end if

  end do

  return
end
subroutine vector_constrained_next7 ( n, level_weight, x_max, x, q_min, q_max, &
  more )

!*****************************************************************************80
!
!! VECTOR_CONSTRAINED_NEXT7 returns the "next" constrained vector.
!
!  Discussion:
!
!    We consider vectors X of dimension N satisfying:
!
!      0 <= X(1:N) <= X_MAX(1:N).
!
!    and the following constraint:
!
!      Q_MIN < sum ( 1 <= I <= N ) LEVEL_WEIGHT(I) * X(I) <= Q_MAX
!
!    This routine returns, one at a time, and in right-to-left
!    lexicographic order, exactly those vectors which satisfy
!    the constraint.
!
!  Example:
!
!    N = 3
!    LEVEL_WEIGHT    4.0  3.0  5.0
!    Q_MIN   16.0
!    Q_MAX   20.0
!    X_MAX:   5   6   4
!
!    #  X(1)  X(2)  X(3)     Total
!
!    1    3     1     1       20.0
!    2    2     2     1       19.0
!    3    1     3     1       18.0
!    4    1     1     2       17.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components in the vector.
!
!    Input, real ( kind = 8 ) LEVEL_WEIGHT(N), the coefficient vector.
!
!    Input, integer ( kind = 4 ) X_MAX(N), the maximum
!    values allowed in each component.
!
!    Input/output, integer ( kind = 4 ) X(N).  On first call, with
!    MORE = FALSE, the input value of X is not important.  On subsequent calls,
!    the input value of X should be the output value from the previous call.
!    On output, (with MORE = TRUE), the value of X will be the "next"
!    vector in the reverse lexicographical list of vectors that satisfy
!    the condition.  However, on output with MORE = FALSE, the vector
!    X is meaningless, because there are no more vectors in the list.
!
!    Input, real ( kind = 8 ) Q_MIN, Q_MAX, the lower and upper
!    limits on the sum.
!
!    Input/output, logical MORE.  On input, if the user has set MORE
!    FALSE, the user is requesting the initiation of a new sequence
!    of values.  If MORE is TRUE, then the user is requesting "more"
!    values in the current sequence.  On output, if MORE is TRUE,
!    then another value was found and returned in X, but if MORE is
!    FALSE, then there are no more values in the sequence, and X is
!    NOT the next value.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) level_weight(n)
  logical more
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  real ( kind = 8 ) total
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x_max(n)

  if ( .not. more ) then

    more = .true.

    x(1:n) = 0

    total = dot_product ( level_weight(1:n), real ( x(1:n), kind = 8 ) )

    if ( q_min < total .and. total <= q_max ) then
      return
    end if

  end if

  do

    i = n

    do

      if ( x(i) < x_max(i) ) then
        exit
      end if

      if ( i <= 1 ) then
        more = .false.
        return
      end if

      i = i - 1

    end do

    x(i) = x(i) + 1
    x(i+1:n) = 0

    total = dot_product ( level_weight(1:n), real ( x(1:n), kind = 8 ) )

    if ( q_min < total .and. total <= q_max ) then
      exit
    end if

  end do

  return
end
subroutine vector_next ( n, x_min, x_max, x, more )

!*****************************************************************************80
!
!! VECTOR_NEXT returns the "next" integer vector between two ranges.
!
!  Discussion:
!
!    We consider all integer vectors of dimension N satisfying:
!
!      X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
!
!    This routine returns, one at a time, and in right-to-left
!    lexicographic order, all these vectors.
!
!  Example:
!
!    N = 3
!    X_MIN:   2   2   0
!    X_MAX:   4   3   1
!
!    #  X(1)  X(2)  X(3)
!
!    1    2     2     0
!    2    3     2     0
!    3    4     2     0
!    4    2     3     0
!    5    3     3     0
!    6    4     3     0
!    7    2     2     1
!    8    3     2     1
!    9    4     2     1
!   10    2     3     1
!   11    3     3     1
!   12    4     3     1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components in the vector.
!
!    Input, integer ( kind = 4 ) X_MIN(N), X_MAX(N), the minimum and maximum
!    values allowed in each component.
!
!    Input/output, integer ( kind = 4 ) X(N).  On first call, with
!    MORE = FALSE, the input value of X is not important.  On subsequent calls,
!    the input value of X should be the output value from the previous call.
!    On output, with MORE = TRUE, the value of X will be the "next"
!    vector in the reverse lexicographical list of vectors.  However, on
!    output with MORE = FALSE, the vector X is meaningless, because there
!    are no more vectors in the list.
!
!    Input/output, logical MORE.  On input, if the user has set MORE
!    FALSE, the user is requesting the initiation of a new sequence
!    of values.  If MORE is TRUE, then the user is requesting "more"
!    values in the current sequence.  On output, if MORE is TRUE,
!    then another value was found and returned in X, but if MORE is
!    FALSE, then there are no more values in the sequence, and X is
!    NOT the next value.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical more
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x_max(n)
  integer ( kind = 4 ) x_min(n)

  if ( .not. more ) then

    x(1:n) = x_min(1:n)
    more = .true.

  else

    i = 1

    do

      if ( x(i) < x_max(i) ) then
        x(i) = x(i) + 1
        exit
      end if

      x(i) = x_min(i)

      i = i + 1

      if ( n < i ) then
        more = .false.
        exit
      end if

    end do

  end if

  return
end
subroutine ytb_enum ( n, ytb_num )

!*****************************************************************************80
!
!! YTB_ENUM enumerates the Young tables of size N.
!
!  Discussion:
!
!    If A(N) is the number of Young table of size N, then A(1) = 1,
!    A(2) = 2, and
!
!    A(N) = A(N-1) + (N-1) * A(N-2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer which is to be partitioned.
!
!    Output, integer ( kind = 4 ) YTB_NUM, the number of Young tables of N.
!
  implicit none

  integer ( kind = 4 ) a1
  integer ( kind = 4 ) a2
  integer ( kind = 4 ) a3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ytb_num

  if ( n <= 0 ) then
    ytb_num = 0
  else if ( n == 1 ) then
    ytb_num = 1
  else if ( n == 2 ) then
    ytb_num = 2
  else
    a2 = 1
    a3 = 2
    do i = 3, n
      a1 = a2
      a2 = a3
      a3 = a2 + ( i - 1 ) * a1
    end do
    ytb_num = a3
  end if

  return
end
subroutine ytb_next ( n, lambda, a, more )

!*****************************************************************************80
!
!! YTB_NEXT computes the next Young table for a given shape.
!
!  Discussion:
!
!    When the routine is called with MORE = FALSE (the first time), and
!    if LAMBDA on this call has M parts, with M < N, then the user
!    must also make sure that LAMBDA(M+1) = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 June 2004
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
!    Input, integer ( kind = 4 ) N, the integer which is to be partitioned.
!
!    Input/output, integer ( kind = 4 ) LAMBDA(N), contains a partition of N.
!    The elements of LAMBDA are nonnegative integers that sum to N.
!    On the first call, with MORE = FALSE, the user sets LAMBDA.
!    After the first call, the input value of LAMBDA is not important.
!    On output, the value of LAMBDA is the partition corresponding
!    to the Young table.
!
!    Input/output, integer ( kind = 4 ) A(N).  On the first call, with
!    MORE = FALSE, no value of A needs to be set.  After the first call, the
!    input value of A should be its output value from the previous call.
!    The output value of A is the next Young table.  A(I) is the
!    row containing I in the output table.
!
!    Input/output, logical MORE.  Set MORE to FALSE before the first call.
!    It is reset to TRUE as the program returns a new table
!    on each call, until the last table is computed, when
!    the program also sets MORE = FALSE.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lambda(n)
  integer ( kind = 4 ) isave
  logical more

  it = n

  if ( more ) then

    lambda(1) = 1
    lambda(2:n) = 0

    isave = 0

    do i = 2, n

      lambda(a(i)) = lambda(a(i)) + 1

      if ( a(i) < a(i-1) ) then
        isave = i
        exit
      end if

    end do

    if ( isave == 0 ) then
      more = .false.
      return
    end if

    it = lambda(1+a(isave))

    do i = n, 1, -1

      if ( lambda(i) == it ) then
        a(isave) = i
        lambda(i) = lambda(i) - 1
        it = isave - 1
        exit
      end if

    end do

  end if

  k = 1
  ir = 1

  do

    if ( n < ir ) then
      exit
    end if

    if ( lambda(ir) /= 0 ) then
      a(k) = ir
      lambda(ir) = lambda(ir) - 1
      k = k + 1
      ir = ir + 1
      cycle
    end if

    if ( it < k ) then
      exit
    end if

    ir = 1

  end do

  if ( n == 1 ) then
    more = .false.
    return
  end if

  do j = 2, n
    if ( a(j) < a(j-1) ) then
      more = .true.
      return
    end if
  end do

  more = .false.

  return
end
subroutine ytb_print ( n, a, title )

!*****************************************************************************80
!
!! YTB_PRINT prints a Young table.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer that is partitioned.
!
!    Input, integer ( kind = 4 ) A(N), describes the Young table.
!    A(I) is the row of the table on which I occurs.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jarray(n)
  integer ( kind = 4 ) row_i
  integer ( kind = 4 ) row_length
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  row_i = 0

  do

    row_i = row_i + 1

    row_length = 0

    do j = 1, n

      if ( a(j) == row_i ) then
        row_length = row_length + 1
        jarray(row_length) = j
      end if

    end do

    if ( row_length <= 0 ) then
      exit
    end if

    write ( *, '(20i4)' ) jarray(1:row_length)

  end do

  return
end
subroutine ytb_random ( n, lambda, seed, a )

!*****************************************************************************80
!
!! YTB_RANDOM selects a random Young table of a given shape.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2003
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
!    Input, integer ( kind = 4 ) N, the integer which has been partitioned.
!
!    Input, integer ( kind = 4 ) LAMBDA(N), the partition of N.
!    N = sum ( 1 <= I <= N ) LAMBDA(I).
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, integer ( kind = 4 ) A(N), the vector describing the Young table.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) ih
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lambda(n)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) seed

  a(1:n) = 0

  i = 0
  k = 0

  do

    i = i + 1
    do j = 1, lambda(i)
      a(j) = a(j) + 1
      k = k + 1
    end do

    if ( n <= k ) then
      exit
    end if

  end do

  do m = 1, n

    do

      i = i4_uniform ( 1, a(1), seed )
      j = i4_uniform ( 1, lambda(1), seed )

      if ( i <= a(j) .and. j <= lambda(i) ) then
        exit
      end if

    end do

    do

      ih = a(j) + lambda(i) - i - j

      if ( ih == 0 ) then
        exit
      end if

      k = i4_uniform ( 1, ih, seed )

      if ( k <= lambda(i)-j ) then
        j = j + k
      else
        i = k - lambda(i) + i + j
      end if

    end do

    lambda(i) = lambda(i) - 1
    a(j) = a(j) - 1
    a(n+1-m) = i

  end do

  do i = 1, n
    lambda(a(i)) = lambda(a(i)) + 1
  end do

  return
end
