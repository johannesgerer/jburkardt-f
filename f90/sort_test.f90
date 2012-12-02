program main

!*****************************************************************************80
!
!! MAIN is the main program for SORT_TEST.
!
!  Discussion:
!
!    SORT_TEST tests out a sorting routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 10

  real, dimension(n) :: a
  real ahi
  real alo
  integer i

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SORT_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate the use of R4VEC_SORT_BUBBLE'
  write ( *, '(a)' ) '  to sort a real array using bubble sort.'

  do i = 1, 2

    alo = 10.0E+00
    ahi = 25.0E+00

    call r4vec_random ( alo, ahi, n, a )

    call r4vec_print ( n, a, '  Unsorted array:' )

    call r4vec_sort_bubble ( n, a )

    call r4vec_print ( n, a, '  Sorted array:' )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SORT_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine r4vec_random ( alo, ahi, n, a )

!*****************************************************************************80
!
!! R4VEC_RANDOM sets a real vector to random values.
!
!  Discussion:
!
!    The random values are chosen by a call to the intrinsic routine
!    RANDOM_NUMBER.  Hence, the user may set the seed by a call to
!    RANDOM_SEED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ALO, AHI, the range of values desired in the vector.
!
!    Input, integer N, the dimension of the vector.
!
!    Output, real A(N), random values between ALO and AHI.
!
  implicit none

  integer, intent ( in ) :: n

  real, intent ( in ) :: ahi
  real, intent ( in ) :: alo
  real, dimension ( n ), intent ( out ) :: a

  call random_number ( a )

  a = alo + a * ( ahi - alo )

end
subroutine r4vec_sort_bubble ( n, a )

!*****************************************************************************80
!
!! R4VEC_SORT_BUBBLE sorts a real vector using bubble sort.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vector.
!
!    Input/output, real A(N), the vector to be sorted.
!
  implicit none

  integer, intent ( in ) :: n

  integer :: i
  integer :: j
  real, dimension ( n ), intent ( inout ) :: a

  do i = 1, n

    do j = i+1, n

      if ( a(j) < a(i) ) then
        call r4_swap ( a(i), a(j) )
      end if

    end do

  end do

  return
end
subroutine r4_swap ( a1, a2 )

!*****************************************************************************80
!
!! R4_SWAP swaps two real values
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real A1, A2, two real values to swap.
!
  implicit none

  real, intent ( inout ) :: a1
  real, intent ( inout ) :: a2
  real :: a3

  a3 = a1
  a1 = a2
  a2 = a3

  return
end
subroutine r4vec_print ( n, a, title )

!*****************************************************************************80
!
!! R4VEC_PRINT prints a real vector.
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
!    Input, integer N, the number of components of the vector.
!
!    Input, real A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer n

  real a(n)
  integer i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, a(i)
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

  character ( len = 8  ) ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
