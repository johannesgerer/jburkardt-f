program main

!*****************************************************************************80
!
!! MAIN is the main program for BOUNDER.
!
!  Discussion:
!
!    BOUNDER is a program that generates an out-of-bounds array access.
!
!  Modified:
!
!    27 January 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BOUNDER'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Generate an out-of-bounds array access.'

  call test01
  call test02

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BOUNDER'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
subroutine test01

!*****************************************************************************80
!
!! TEST01 is the CORRECT code.
!
  implicit none

  integer, parameter :: n = 10

  integer a(n)
  integer b(n)
  integer i
  integer j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  The CORRECT code.'

  do i = 1, n
    a(i) = i + 1
  end do
  a(n) = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  A(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2i6)' ) i, a(i)
  end do

  do i = 1, n
    b(i) = - 1000 - i
  end do

  do i = 1, n
    j = a(i)
    b(j) = j + 1
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  B(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2i6)' ) i, b(i)
  end do

  return
end
subroutine test02

!*****************************************************************************80
!
!! TEST02 is the INCORRECT code.
!
  implicit none

  integer, parameter :: n = 10

  integer a(n)
  integer b(n)
  integer i
  integer j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  The INCORRECT code.'

  do i = 1, n
    a(i) = i + 1
  end do
  a(n) = 100000

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  A(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2i6)' ) i, a(i)
  end do

  do i = 1, n
    b(i) = - 1000 - i
  end do

  do i = 1, n
    j = a(i)
    b(j) = j + 1
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  B(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2i6)' ) i, b(i)
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
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

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
