function index0 ( i_min, i, i_max )

!*****************************************************************************80
!
!! INDEX0 indexes a 1D vector using a zero base.
!
!  Discussion:
!
!    Index       Element
!    ---------   --------
!    0           I_MIN
!    INDEX0      I
!   (INDEX_MAX)  I_MAX
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_MIN, I, I_MAX, for the first index,
!    the minimum, the index, and the maximum.
!
!    Output, integer ( kind = 4 ) INDEX0, the index of element I.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ), parameter :: index_min = 0
  integer ( kind = 4 ) index0

  index0 = index_min + ( i - i_min )

  return
end
function index01 ( i_min, i, i_max, j_min, j, j_max )

!*****************************************************************************80
!
!! INDEX01 indexes a 2D array by columns, with a zero base.
!
!  Discussion:
!
!    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
!    and increasing the row index first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_MIN, I, I_MAX, for row indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) J_MIN, J, J_MAX, for column indices,
!    the minimum, the index, and the maximum.
!
!    Output, integer ( kind = 4 ) INDEX01, the index of element (I,J).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ), parameter :: index_min = 0
  integer ( kind = 4 ) index01
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_max
  integer ( kind = 4 ) j_min

  index01 = &
    index_min &
    + (         i - i_min ) &
    + ( i_max + 1 - i_min ) * ( j - j_min )

  return
end
function index012 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max )

!*****************************************************************************80
!
!! INDEX012 indexes a 3D array by columns with zero base.
!
!  Discussion:
!
!    Entries of the array are indexed starting at entry (I_MIN,J_MIN,K_MIN), 
!    and increasing the row index first, then the column index.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_MIN, I, I_MAX, for row indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) J_MIN, J, J_MAX, for column indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) K_MIN, K, K_MAX, for plane indices,
!    the minimum, the index, and the maximum.
!
!    Output, integer ( kind = 4 ) INDEX012, the index of element (I,J,K).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ), parameter :: index_min = 0
  integer ( kind = 4 ) index012
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_max
  integer ( kind = 4 ) j_min
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) k_min

  index012 = &
      index_min &
    + (         i - i_min ) &
    + ( i_max + 1 - i_min ) * (         j - j_min ) *  &
    + ( i_max + 1 - i_min ) * ( j_max + 1 - j_min ) * ( k - k_min )

  return
end
function index0123 ( i1_min, i1, i1_max, i2_min, i2, i2_max, i3_min, i3, &
  i3_max, i4_min, i4, i4_max )

!*****************************************************************************80
!
!! INDEX0123 indexes a 4D array by columns, with a zero base.
!
!  Discussion:
!
!    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
!    and increasing the initial index first, then the second, third and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I1_MIN, I1, I1_MAX, for index 1,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I2_MIN, I2, I2_MAX, for index 2,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I3_MIN, I3, I3_MAX, for index 3,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I4_MIN, I4, I4_MAX, for index 4,
!    the minimum, the index, and the maximum.
!
!    Output, integer ( kind = 4 ) INDEX0123, the index of (I1,I2,I3,I4).
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i1_max
  integer ( kind = 4 ) i1_min
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2_max
  integer ( kind = 4 ) i2_min
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i3_max
  integer ( kind = 4 ) i3_min
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_max
  integer ( kind = 4 ) i4_min
  integer ( kind = 4 ), parameter :: index_min = 0
  integer ( kind = 4 ) index0123

  index0123 = &
      index_min &
    + (         i1 - i1_min ) &
    + ( i1_max + 1 - i1_min ) * (         i2 - i2_min ) &
    + ( i1_max + 1 - i1_min ) * ( i2_max + 1 - i2_min ) &
    * (         i3 - i3_min ) &
    + ( i1_max + 1 - i1_min ) * ( i2_max + 1 - i2_min ) &
    * ( i3_max + 1 - i3_min ) * (         i4 - i4_min )

  return
end
function index0n ( n, i_min, i, i_max )

!*****************************************************************************80
!
!! INDEX0N indexes an N-dimensional array by columns, with zero base.
!
!  Discussion:
!
!    Entries of the array are indexed starting at entry 
!      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
!    and increasing the first index up to I_MAX(1), 
!    then the second and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of indices.
!
!    Input, integer ( kind = 4 ) I_MIN(N), the minimum indices.
!
!    Input, integer ( kind = 4 ) I(N), the indices.
!
!    Input, integer ( kind = 4 ) I_MAX(N), for maximum indices.
!
!    Output, integer ( kind = 4 ) INDEX0N, the index of element I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i(n)
  integer ( kind = 4 ) i_max(n)
  integer ( kind = 4 ) i_min(n)
  integer ( kind = 4 ), parameter :: index_min = 0
  integer ( kind = 4 ) index0n
  integer ( kind = 4 ) j
  integer ( kind = 4 ) value

  value = ( i(n) - i_min(n) )

  do j = n - 1, 1, - 1
    value = value * ( i_max(j) + 1 - i_min(j) ) + ( i(j) - i_min(j) )
  end do
  value = value + index_min

  index0n = value

  return
end
function index1 ( i_min, i, i_max )

!*****************************************************************************80
!
!! INDEX1 indexes a 1D vector using a unit base.
!
!  Discussion:
!
!    Index       Element
!    ---------   --------
!    1           I_MIN
!    INDEX1      I
!   (INDEX_MAX)  I_MAX
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_MIN, I, I_MAX, for the first index,
!    the minimum, the index, and the maximum.
!
!    Output, integer ( kind = 4 ) INDEX1, the index of element I.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ), parameter :: index_min = 1
  integer ( kind = 4 ) index1

  index1 = index_min + ( i - i_min )

  return
end
function index10 ( i_min, i, i_max, j_min, j, j_max )

!*****************************************************************************80
!
!! INDEX10 indexes a 2D array by rows, with a zero base.
!
!  Discussion:
!
!    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
!    and increasing the column index first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_MIN, I, I_MAX, for row indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) J_MIN, J, J_MAX, for column indices,
!    the minimum, the index, and the maximum.
!
!    Output, integer ( kind = 4 ) INDEX10, the index of element (I,J).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ), parameter :: index_min = 0
  integer ( kind = 4 ) index10
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_max
  integer ( kind = 4 ) j_min

  index10 = index_min &
             +                         ( j - j_min ) &
             + ( i - i_min ) * ( j_max + 1 - j_min )

  return
end
function index12 ( i_min, i, i_max, j_min, j, j_max )

!*****************************************************************************80
!
!! INDEX12 indexes a 2D array by columns, with a unit base.
!
!  Discussion:
!
!    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
!    and increasing the row index first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_MIN, I, I_MAX, for row indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) J_MIN, J, J_MAX, for column indices,
!    the minimum, the index, and the maximum.
!
!    Output, integer ( kind = 4 ) INDEX12, the index of element (I,J).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ), parameter :: index_min = 1
  integer ( kind = 4 ) index12
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_max
  integer ( kind = 4 ) j_min

  index12 = &
    index_min &
    + (         i - i_min ) &
    + ( i_max + 1 - i_min ) * ( j - j_min )

  return
end
function index123 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max )

!*****************************************************************************80
!
!! INDEX123 indexes a 3D array by columns with unit base.
!
!  Discussion:
!
!    Entries of the array are indexed starting at entry (I_MIN,J_MIN,K_MIN), 
!    and increasing the row index first, then the column index.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_MIN, I, I_MAX, for row indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) J_MIN, J, J_MAX, for column indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) K_MIN, K, K_MAX, for plane indices,
!    the minimum, the index, and the maximum.
!
!    Output, integer ( kind = 4 ) INDEX123, the index of element (I,J,K).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ), parameter :: index_min = 1
  integer ( kind = 4 ) index123
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_max
  integer ( kind = 4 ) j_min
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) k_min

  index123 = &
      index_min &
    + (         i - i_min ) &
    + ( i_max + 1 - i_min ) * (         j - j_min ) *  &
    + ( i_max + 1 - i_min ) * ( j_max + 1 - j_min ) * ( k - k_min )

  return
end
function index1234 ( i1_min, i1, i1_max, i2_min, i2, i2_max, i3_min, i3, &
  i3_max, i4_min, i4, i4_max )

!*****************************************************************************80
!
!! INDEX1234 indexes a 4D array by columns, with a unit base.
!
!  Discussion:
!
!    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
!    and increasing the initial index first, then the second, third and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I1_MIN, I1, I1_MAX, for index 1,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I2_MIN, I2, I2_MAX, for index 2,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I3_MIN, I3, I3_MAX, for index 3,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I4_MIN, I4, I4_MAX, for index 4,
!    the minimum, the index, and the maximum.
!
!    Output, integer ( kind = 4 ) INDEX1234, the index of (I1,I2,I3,I4).
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i1_max
  integer ( kind = 4 ) i1_min
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2_max
  integer ( kind = 4 ) i2_min
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i3_max
  integer ( kind = 4 ) i3_min
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_max
  integer ( kind = 4 ) i4_min
  integer ( kind = 4 ), parameter :: index_min = 1
  integer ( kind = 4 ) index1234

  index1234 = &
      index_min &
    + (         i1 - i1_min ) &
    + ( i1_max + 1 - i1_min ) * (         i2 - i2_min ) &
    + ( i1_max + 1 - i1_min ) * ( i2_max + 1 - i2_min ) &
    * (         i3 - i3_min ) &
    + ( i1_max + 1 - i1_min ) * ( i2_max + 1 - i2_min ) &
    * ( i3_max + 1 - i3_min ) * (         i4 - i4_min )

  return
end
function index1n ( n, i_min, i, i_max )

!*****************************************************************************80
!
!! INDEX1N indexes an N-dimensional array by columns, with unit base.
!
!  Discussion:
!
!    Entries of the array are indexed starting at entry 
!      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
!    and increasing the first index up to I_MAX(1), 
!    then the second and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of indices.
!
!    Input, integer ( kind = 4 ) I_MIN(N), the minimum indices.
!
!    Input, integer ( kind = 4 ) I(N), the indices.
!
!    Input, integer ( kind = 4 ) I_MAX(N), for maximum indices.
!
!    Output, integer ( kind = 4 ) INDEX1N, the index of element I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i(n)
  integer ( kind = 4 ) i_max(n)
  integer ( kind = 4 ) i_min(n)
  integer ( kind = 4 ), parameter :: index_min = 1
  integer ( kind = 4 ) index1n
  integer ( kind = 4 ) j
  integer ( kind = 4 ) value

  value = ( i(n) - i_min(n) )

  do j = n - 1, 1, - 1
    value = value * ( i_max(j) + 1 - i_min(j) ) + ( i(j) - i_min(j) )
  end do
  value = value + index_min

  index1n = value

  return
end
function index21 ( i_min, i, i_max, j_min, j, j_max )

!*****************************************************************************80
!
!! INDEX21 indexes a 2D array by rows, with a unit base.
!
!  Discussion:
!
!    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
!    and increasing the column index first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_MIN, I, I_MAX, for row indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) J_MIN, J, J_MAX, for column indices,
!    the minimum, the index, and the maximum.
!
!    Output, integer ( kind = 4 ) INDEX21, the index of element (I,J).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ), parameter :: index_min = 1
  integer ( kind = 4 ) index21
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_max
  integer ( kind = 4 ) j_min

  index21 = index_min &
             +                         ( j - j_min ) &
             + ( i - i_min ) * ( j_max + 1 - j_min )

  return
end
function index210 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max )

!*****************************************************************************80
!
!! INDEX210 indexes a 3D array by rows, with zero base.
!
!  Discussion:
!
!    When we say "by rows", we really just mean that entries of the array are 
!    indexed starting at entry (I_MIN,J_MIN,K_MIN), and the increasing the LAST
!    index first, then the next-to-the-last, and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_MIN, I, I_MAX, for row indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) J_MIN, J, J_MAX, for column indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) K_MIN, K, K_MAX, for plane indices,
!    the minimum, the index, and the maximum.
!
!    Output, integer ( kind = 4 ) INDEX210, the index of element (I,J,K).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ), parameter :: index_min = 0
  integer ( kind = 4 ) index210
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_max
  integer ( kind = 4 ) j_min
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) k_min

  index210 = &
      index_min &
    +                                                 ( k - k_min ) &
    +                         ( j - j_min ) * ( k_max + 1 - k_min ) &
    + ( i - i_min ) * ( j_max + 1 - j_min ) * ( k_max + 1 - k_min )

  return
end
function index321 ( i_min, i, i_max, j_min, j, j_max, k_min, k, k_max )

!*****************************************************************************80
!
!! INDEX321 indexes a 3D array by rows, with zero base.
!
!  Discussion:
!
!    When we say "by rows", we really just mean that entries of the array are 
!    indexed starting at entry (I_MIN,J_MIN,K_MIN), and the increasing the LAST
!    index first, then the next-to-the-last, and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I_MIN, I, I_MAX, for row indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) J_MIN, J, J_MAX, for column indices,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) K_MIN, K, K_MAX, for plane indices,
!    the minimum, the index, and the maximum.
!
!    Output, integer ( kind = 4 ) INDEX321, the index of element (I,J,K).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_max
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ), parameter :: index_min = 1
  integer ( kind = 4 ) index321
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_max
  integer ( kind = 4 ) j_min
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  integer ( kind = 4 ) k_min

  index321 = &
      index_min &
    +                                                 ( k - k_min ) &
    +                         ( j - j_min ) * ( k_max + 1 - k_min ) &
    + ( i - i_min ) * ( j_max + 1 - j_min ) * ( k_max + 1 - k_min )

  return
end
function index3210 ( i1_min, i1, i1_max, i2_min, i2, i2_max, i3_min, i3, &
  i3_max, i4_min, i4, i4_max )

!*****************************************************************************80
!
!! INDEX3210 indexes a 4D array by rows, with zero base.
!
!  Discussion:
!
!    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
!    and increasing the last index, then the next to last, and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I1_MIN, I1, I1_MAX, for index 1,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I2_MIN, I2, I2_MAX, for index 2,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I3_MIN, I3, I3_MAX, for index 3,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I4_MIN, I4, I4_MAX, for index 4,
!    the minimum, the index, and the maximum.
!
!    Output, integer ( kind = 4 ) INDEX3210, the index of (I1,I2,I3,I4).
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i1_max
  integer ( kind = 4 ) i1_min
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2_max
  integer ( kind = 4 ) i2_min
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i3_max
  integer ( kind = 4 ) i3_min
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_max
  integer ( kind = 4 ) i4_min
  integer ( kind = 4 ), parameter :: index_min = 0
  integer ( kind = 4 ) index3210

  index3210 = &
    index_min &
    +         ( i4 - i4_min ) &
    +                                                     ( i3 - i3_min ) &
    * ( i4_max + 1 - i4_min ) &
    +                           ( i2 - i2_min ) * ( i3_max + 1 - i3_min ) &
    * ( i4_max + 1 - i4_min ) &
    + ( i1 - i1_min ) * ( i2_max + 1 - i2_min ) * ( i3_max + 1 - i3_min ) &
    * ( i4_max + 1 - i4_min )

  return
end
function index4321 ( i1_min, i1, i1_max, i2_min, i2, i2_max, i3_min, i3, &
  i3_max, i4_min, i4, i4_max )

!*****************************************************************************80
!
!! INDEX4321 indexes a 4D array by rows, with unit base.
!
!  Discussion:
!
!    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
!    and increasing the last index, then the next to last, and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I1_MIN, I1, I1_MAX, for index 1,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I2_MIN, I2, I2_MAX, for index 2,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I3_MIN, I3, I3_MAX, for index 3,
!    the minimum, the index, and the maximum.
!
!    Input, integer ( kind = 4 ) I4_MIN, I4, I4_MAX, for index 4,
!    the minimum, the index, and the maximum.
!
!    Output, integer ( kind = 4 ) INDEX4321, the index of (I1,I2,I3,I4).
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i1_max
  integer ( kind = 4 ) i1_min
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2_max
  integer ( kind = 4 ) i2_min
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i3_max
  integer ( kind = 4 ) i3_min
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_max
  integer ( kind = 4 ) i4_min
  integer ( kind = 4 ), parameter :: index_min = 1
  integer ( kind = 4 ) index4321

  index4321 = &
    index_min &
    +         ( i4 - i4_min ) &
    +                                                     ( i3 - i3_min ) &
    * ( i4_max + 1 - i4_min ) &
    +                           ( i2 - i2_min ) * ( i3_max + 1 - i3_min ) &
    * ( i4_max + 1 - i4_min ) &
    + ( i1 - i1_min ) * ( i2_max + 1 - i2_min ) * ( i3_max + 1 - i3_min ) &
    * ( i4_max + 1 - i4_min )

  return
end
function indexn0 ( n, i_min, i, i_max )

!*****************************************************************************80
!
!! INDEXN0 indexes an N-dimensional array by rows, with zero base.
!
!  Discussion:
!
!    Entries of the array are indexed starting at entry 
!      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
!    and increasing the last index up to I_MAX(N), 
!    then the next-to-last and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of indices.
!
!    Input, integer ( kind = 4 ) I_MIN(N), the minimum indices.
!
!    Input, integer ( kind = 4 ) I(N), the indices.
!
!    Input, integer ( kind = 4 ) I_MAX(N), for maximum indices.
!
!    Output, integer ( kind = 4 ) INDEXN0, the index of element I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i(n)
  integer ( kind = 4 ) i_max(n)
  integer ( kind = 4 ) i_min(n)
  integer ( kind = 4 ), parameter :: index_min = 0
  integer ( kind = 4 ) indexn0
  integer ( kind = 4 ) j
  integer ( kind = 4 ) value

  value = ( i(1) - i_min(1) )

  do j = 2, n
    value = value * ( i_max(j) + 1 - i_min(j) ) + ( i(j) - i_min(j) )
  end do
  value = value + index_min

  indexn0 = value

  return
end
function indexn1 ( n, i_min, i, i_max )

!*****************************************************************************80
!
!! INDEXN1 indexes an N-dimensional array by rows, with unit base.
!
!  Discussion:
!
!    Entries of the array are indexed starting at entry 
!      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
!    and increasing the last index up to I_MAX(N), 
!    then the next-to-last and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of indices.
!
!    Input, integer ( kind = 4 ) I_MIN(N), the minimum indices.
!
!    Input, integer ( kind = 4 ) I(N), the indices.
!
!    Input, integer ( kind = 4 ) I_MAX(N), for maximum indices.
!
!    Output, integer ( kind = 4 ) INDEXN1, the index of element I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i(n)
  integer ( kind = 4 ) i_max(n)
  integer ( kind = 4 ) i_min(n)
  integer ( kind = 4 ), parameter :: index_min = 1
  integer ( kind = 4 ) indexn1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) value

  value = ( i(1) - i_min(1) )

  do j = 2, n
    value = value * ( i_max(j) + 1 - i_min(j) ) + ( i(j) - i_min(j) )
  end do
  value = value + index_min

  indexn1 = value

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

