function d1mach ( i )

!*****************************************************************************80
!
!! D1MACH returns double precision real machine constants.
!
!  Discussion:
!
!    Assuming that the internal representation of a double precision real
!    number is in base B, with T the number of base-B digits in the mantissa,
!    and EMIN the smallest possible exponent and EMAX the largest possible 
!    exponent, then
!
!      D1MACH(1) = B^(EMIN-1), the smallest positive magnitude.
!      D1MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.
!      D1MACH(3) = B^(-T), the smallest relative spacing.
!      D1MACH(4) = B^(1-T), the largest relative spacing.
!      D1MACH(5) = log10(B).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer,
!    Algorithm 528:
!    Framework for a Portable Library,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, page 176-188.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, chooses the parameter to be returned.
!    1 <= I <= 5.
!
!    Output, real ( kind = 8 ) D1MACH, the value of the chosen parameter.
!
  implicit none

  real ( kind = 8 ) d1mach
  integer ( kind = 4 ) i

  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'D1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    d1mach = 0.0D+00
    stop
  else if ( i == 1 ) then
    d1mach = 4.450147717014403D-308
  else if ( i == 2 ) then
    d1mach = 8.988465674311579D+307
  else if ( i == 3 ) then
    d1mach = 1.110223024625157D-016
  else if ( i == 4 ) then
    d1mach = 2.220446049250313D-016
  else if ( i == 5 ) then
    d1mach = 0.301029995663981D+000
  else if ( 5 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'D1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    d1mach = 0.0D+00
    stop
  end if

  return
end
function i1mach ( i )

!*****************************************************************************80
!
!! I1MACH returns integer machine constants.
!
!  Discussion:
!
!    Input/output unit numbers.
!
!      I1MACH(1) = the standard input unit.
!      I1MACH(2) = the standard output unit.
!      I1MACH(3) = the standard punch unit.
!      I1MACH(4) = the standard error message unit.
!
!    Words.
!
!      I1MACH(5) = the number of bits per integer storage unit.
!      I1MACH(6) = the number of characters per integer storage unit.
!
!    Integers.
!
!    Assume integers are represented in the S digit base A form:
!
!      Sign * (X(S-1)*A^(S-1) + ... + X(1)*A + X(0))
!
!    where 0 <= X(1:S-1) < A.
!
!      I1MACH(7) = A, the base.
!      I1MACH(8) = S, the number of base A digits.
!      I1MACH(9) = A**S-1, the largest integer.
!
!    Floating point numbers
!
!    Assume floating point numbers are represented in the T digit 
!    base B form:
!
!      Sign * (B^E) * ((X(1)/B) + ... + (X(T)/B^T) )
!
!    where 0 <= X(I) < B for I=1 to T, 0 < X(1) and EMIN <= E <= EMAX.
!
!      I1MACH(10) = B, the base.
!
!    Single precision
!
!      I1MACH(11) = T, the number of base B digits.
!      I1MACH(12) = EMIN, the smallest exponent E.
!      I1MACH(13) = EMAX, the largest exponent E.
!
!    Double precision
!
!      I1MACH(14) = T, the number of base B digits.
!      I1MACH(15) = EMIN, the smallest exponent E.
!      I1MACH(16) = EMAX, the largest exponent E.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer,
!    Algorithm 528,
!    Framework for a Portable Library,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, page 176-188.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, chooses the parameter to be returned.
!    1 <= I <= 16.
!
!    Output, integer ( kind = 4 ) I1MACH, the value of the chosen parameter.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1mach

  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
    write ( *, '(a,i12)' ) '  I = ', i
    i1mach = 0
    stop
  else if ( i == 1 ) then
    i1mach = 5
  else if ( i == 2 ) then
    i1mach = 6
  else if ( i == 3 ) then
    i1mach = 7
  else if ( i == 4 ) then
    i1mach = 6
  else if ( i == 5 ) then
    i1mach = 32
  else if ( i == 6 ) then
    i1mach = 4
  else if ( i == 7 ) then
    i1mach = 2
  else if ( i == 8 ) then
    i1mach = 31
  else if ( i == 9 ) then
    i1mach = 2147483647
  else if ( i == 10 ) then
    i1mach = 2
  else if ( i == 11 ) then
    i1mach = 24
  else if ( i == 12 ) then
    i1mach = -125
  else if ( i == 13 ) then
    i1mach = 128
  else if ( i == 14 ) then
    i1mach = 53
  else if ( i == 15 ) then
    i1mach = -1021
  else if ( i == 16 ) then
    i1mach = 1024
  else if ( 16 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
    write ( *, '(a,i12)' ) '  I = ', i
    i1mach = 0
    stop
  end if

  return
end
function r1mach ( i )

!*****************************************************************************80
!
!! R1MACH returns single precision real machine constants.
!
!  Discussion:
!
!    Assume that single precision real numbers are stored with a mantissa 
!    of T digits in base B, with an exponent whose value must lie 
!    between EMIN and EMAX.  Then for values of I between 1 and 5, 
!    R1MACH will return the following values:
!
!      R1MACH(1) = B^(EMIN-1), the smallest positive magnitude.
!      R1MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.
!      R1MACH(3) = B^(-T), the smallest relative spacing.
!      R1MACH(4) = B^(1-T), the largest relative spacing.
!      R1MACH(5) = log10(B)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer,
!    Algorithm 528,
!    Framework for a Portable Library,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, page 176-188.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, chooses the parameter to be returned.
!    1 <= I <= 5.
!
!    Output, real ( kind = 4 ) R1MACH, the value of the chosen parameter.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 4 ) r1mach

  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    r1mach = 0.0E+00
    stop
  else if ( i == 1 ) then
    r1mach = 1.1754944E-38
  else if ( i == 2 ) then
    r1mach = 3.4028235E+38
  else if ( i == 3 ) then
    r1mach = 5.9604645E-08
  else if ( i == 4 ) then
    r1mach = 1.1920929E-07
  else if ( i == 5 ) then
    r1mach = 0.3010300E+00
  else if ( 5 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R1MACH - Fatal error!'
    write ( *, '(a)' ) '  The input argument I is out of bounds.'
    write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
    write ( *, '(a,i12)' ) '  I = ', i
    r1mach = 0.0E+00
    stop
  end if
 
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
