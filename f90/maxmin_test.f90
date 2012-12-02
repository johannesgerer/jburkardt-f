program main

!*****************************************************************************80
!
!! MAIN is the main program for MAXMIN_TEST.
!
!  Discussion:
!
!    MAXMIN_TEST tests the F90 MAXVAL, MAXLOC, MINVAL, MINLOC statements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2003
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MAXMIN_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate the use of:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    MAXVAL ( V ) returns maximum value of a vector V;'
  write ( *, '(a)' ) '    MAXVAL ( A ) returns a vector of maximum values;'
  write ( *, '(a)' ) '    MAXLOC ( V ) returns the index of the maximum value;'
  write ( *, '(a)' ) '    MAXLOC ( A ) returns a vector of maximum value indices.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  and the corresponding minimum operators:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    MINVAL ( V );'
  write ( *, '(a)' ) '    MINVAL ( A );'
  write ( *, '(a)' ) '    MINLOC ( V );'
  write ( *, '(a)' ) '    MINLOC ( A ).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Warning: the MAXLOC and MINLOC routines have a very'
  write ( *, '(a)' ) '  fussy and nonintuitive interface, designed to make'
  write ( *, '(a)' ) '  it equally impossible to work with vectors as with'
  write ( *, '(a)' ) '  multi-index arrays.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MAXMIN_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tries out MAXVAL, MAXLOC, MINVAL, MINLOC on a vector.
!
!  Discussion:
!
!    Note that, even though MAXLOC is returning a scalar value,
!    V_MAXLOC MUST be declared a vector quantity, not a scalar,
!    and the assignment must have the form:
!
!      v_maxloc = maxloc ( v )
!
!    but NOT
!
!      v_maxloc(1) = maxloc ( v )
!
!    and NOT
!
!      v_maxloc(1:1) = maxloc ( v )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2003
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 10

  integer i
  real, dimension ( n ) :: v = (/ &
    5.0E+00, 9.0E+00, 2.0E+00, 0.0E+00, 4.0E+00, &
    5.0E+00, 9.0E+00, 0.0E+00, 8.0E+00, 5.0E+00 /)
  integer v_maxloc(1)
  real v_maxval
  integer v_minloc(1)
  real v_minval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For a real vector V, use'
  write ( *, '(a)' ) '  v_maxval = maxval ( v )'
  write ( *, '(a)' ) '  v_maxloc = maxloc ( v )'
  write ( *, '(a)' ) '  v_minval = minval ( v )'
  write ( *, '(a)' ) '  v_minloc = minloc ( v )'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I   V(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, v(i)
  end do

  v_maxval =  maxval ( v )
  v_maxloc =  maxloc ( v )
  v_minval =  minval ( v )
  v_minloc =  minloc ( v )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  v_maxval = ', v_maxval
  write ( *, '(a,i6)' )    '  v_maxloc = ', v_maxloc(1)
  write ( *, '(a,g14.6)' ) '  v_minval = ', v_minval
  write ( *, '(a,i6)' )    '  v_minloc = ', v_minloc(1)

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tries out MAXVAL, MAXLOC, MINVAL, MINLOC on a 5 x 3 array.
!
!  Discussion:
!
!    The strange protocol of these routines in the vector case
!    starts to make a little more sense in the double-dimension
!    array case!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2003
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 5
  integer, parameter :: n = 3

  real, dimension ( m, n ) :: a = reshape ( (/ &
    11.0E+00, 21.0E+00, 31.0E+00, 41.0E+00, 51.0E+00, &
    12.0E+00, 22.0E+00, 32.0E+00, 42.0E+00, 52.0E+00, &
    13.0E+00, 23.0E+00, 33.0E+00, 43.0E+00, 53.0E+00 /), (/ 5, 3 /) )
  integer a_maxloc(2)
  real a_maxval
  integer a_minloc(2)
  real a_minval
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For a real 2D array A, use'
  write ( *, '(a)' ) '  a_maxval = maxval ( a )'
  write ( *, '(a)' ) '  a_maxloc = maxloc ( a )'
  write ( *, '(a)' ) '  a_minval = minval ( a )'
  write ( *, '(a)' ) '  a_minloc = minloc ( a )'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I   A(I,1:N)'
  write ( *, '(a)' ) ' '
  do i = 1, m
    write ( *, '(i6,3g14.6)' ) i, a(i,1:n)
  end do

  a_maxval =  maxval ( a )
  a_maxloc =  maxloc ( a )
  a_minval =  minval ( a )
  a_minloc =  minloc ( a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  a_maxval = ', a_maxval
  write ( *, '(a,2i6)'   ) '  a_maxloc = ', a_maxloc
  write ( *, '(a,g14.6)' ) '  a_minval = ', a_minval
  write ( *, '(a,2i6)'   ) '  a_minloc = ', a_minloc

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
