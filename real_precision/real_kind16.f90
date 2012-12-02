program main

!*****************************************************************************80
!
!! MAIN is the main program for REAL_KIND16.
!
!  Discussion:
!
!    REAL_KIND16 demonstrates how to set up real numbers of various kinds.
!
!    The FORTRAN90 standard attempted to merge the real and double
!    precision data types into a single real data type.  Single precision
!    and double precision numbers were now simply two "kinds" of
!    real number.  This allowed for the possibility of defining and
!    using real numbers of even higher precisions - if the compiler
!    writers chose to provide them.
!
!    The mechanism for doing this was to assign each "kind" a numeric
!    label.  You don't need to specify the kind, but then you only
!    get the default real variable, presumably what we'd think of
!    as single precision:
!
!      real x
!
!    To specify the kind, the FORM of the declaration is easy:
!
!      real ( kind = label ) x
!
!    However, the FORTRAN90 standard does not specify what the numeric
!    value of the label should be.  In fact, most compiler writers
!    have chosen label to reflect the number of bytes associated with
!    the data type.  THIS IS NOT GUARANTEED TO BE THE CASE ON THE
!    COMPILER YOU USE.  But if it is, it makes it easy, because
!    single precision can be explicitly requested by the statement
!
!      real ( kind = 4 ) x
!
!    double precision is requested by
!
!      real ( kind = 8 ) x
!
!    and quad precision (if it happens to be provided) is gotten by
!
!      real ( kind = 16 ) x
!
!
!    In case you don't know the appropriate labels, then
!    What is "guaranteed" to work is something like
!
!      integer, parameter :: kind_val = selected_real_kind ( p, r )
!
!    which asks that KIND_VAL be set to the appropriate KIND value
!    that is able to represent real numbers with at least P decimal
!    digits ("precision") and a decimal exponent range of R.
!
!
!    If you want, you can specify just the number of decimal
!    digits, using the command:
!
!      integer, parameter :: kind_val = selected_real_kind ( p )
!    or
!      integer, parameter :: kind_val = selected_real_kind ( P = p )
!
!    In these examples, "p" is to be replaced by a number, but "P"
!    is literally literally "P".
!
!
!    If you want, you can specify just the exponent range, using
!    a command like:
!
!      integer, parameter :: kind_val = selected_real_kind ( R = r )
!
!    In this example, "r" is to be replaced by a number, but "R"
!    is literally literally "R".
!
!
!    If there is no available real data type that can
!    deliver the precision you're asking for, you're out of luck.
!
!    However, assuming there is such a data type, you can now refer
!    to it using by using label value stored in KIND_VAL in your
!    declarations:
!
!      real ( kind_val ) :: r
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'REAL_KIND16'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Attempt to declare higher precision'
  write ( *, '(a)' ) '  real data using the NONSTANDARD statements:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  REAL ( KIND = 4 ) for single precision,'
  write ( *, '(a)' ) '  REAL ( KIND = 8 ) for double precision,'
  write ( *, '(a)' ) '  REAL ( KIND = 16 ) for quadruple precision.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'REAL_KIND16'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 examines the default real numbers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real r1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use reals of the default type.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HUGE returns the largest value of the given type.'
  write ( *, '(a)' ) '  TINY returns the smallest value of the given type.'
  write ( *, '(a)' ) '  EPSILON returns the precision of a real type.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '    HUGE(R1) =    ', huge ( r1 )
  write ( *, '(a,g14.6)' ) '    TINY(R1) =    ', tiny ( r1 )
  write ( *, '(a,g14.6)' ) '    EPSILON(R1) = ', epsilon ( r1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DIGITS counts the significant binary digits.'
  write ( *, '(a)' ) '  RANGE provides the decimal exponent range.'
  write ( *, '(a)' ) '  PRECISION provides the decimal precision.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '    DIGITS(R1) =    ', digits ( r1 )
  write ( *, '(a,i12)' ) '    RANGE(R1) =     ', range ( r1 )
  write ( *, '(a,i12)' ) '    PRECISION(R1) = ', precision ( r1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RADIX provides the base of the model.'
  write ( *, '(a)' ) '  MAXEXPONENT returns the maximum exponent of a variable.'
  write ( *, '(a)' ) '  MINEXPONENT returns the minimum exponent of a variable.'
  write ( *, '(a)' ) '  (These are exponents of the RADIX)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '    RADIX(R1) =       ', radix ( r1 )
  write ( *, '(a,i12)' ) '    MAXEXPONENT(R1) = ', maxexponent ( r1 )
  write ( *, '(a,i12)' ) '    MINEXPONENT(R1) = ', minexponent ( r1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  KIND returns the "kind" of a variable.'
  write ( *, '(a,i12)' ) '    KIND(R1) = ', kind ( r1 )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 shows what you get with KIND = 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) r1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Declare "real ( kind = 4 ) r1".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HUGE returns the largest value of the given type.'
  write ( *, '(a)' ) '  TINY returns the smallest value of the given type.'
  write ( *, '(a)' ) '  EPSILON returns the precision of a real type.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '    HUGE(R1) =    ', huge ( r1 )
  write ( *, '(a,g14.6)' ) '    TINY(R1) =    ', tiny ( r1 )
  write ( *, '(a,g14.6)' ) '    EPSILON(R1) = ', epsilon ( r1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DIGITS counts the significant binary digits.'
  write ( *, '(a)' ) '  RANGE provides the decimal exponent range.'
  write ( *, '(a)' ) '  PRECISION provides the decimal precision.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '    DIGITS(R1) =    ', digits ( r1 )
  write ( *, '(a,i12)' ) '    RANGE(R1) =     ', range ( r1 )
  write ( *, '(a,i12)' ) '    PRECISION(R1) = ', precision ( r1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RADIX provides the base of the model.'
  write ( *, '(a)' ) '  MAXEXPONENT returns the maximum exponent of a variable.'
  write ( *, '(a)' ) '  MINEXPONENT returns the minimum exponent of a variable.'
  write ( *, '(a)' ) '  (These are exponents of the RADIX)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '    RADIX(R1) =       ', radix ( r1 )
  write ( *, '(a,i12)' ) '    MAXEXPONENT(R1) = ', maxexponent ( r1 )
  write ( *, '(a,i12)' ) '    MINEXPONENT(R1) = ', minexponent ( r1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  KIND returns the "kind" of a variable.'
  write ( *, '(a,i12)' ) '    KIND(R1) = ', kind ( r1 )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 shows what you get with KIND = 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) r1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Declare "real ( kind = 8 ) r1".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HUGE returns the largest value of the given type.'
  write ( *, '(a)' ) '  TINY returns the smallest value of the given type.'
  write ( *, '(a)' ) '  EPSILON returns the precision of a real type.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '    HUGE(R1) =    ', huge ( r1 )
  write ( *, '(a,g14.6)' ) '    TINY(R1) =    ', tiny ( r1 )
  write ( *, '(a,g14.6)' ) '    EPSILON(R1) = ', epsilon ( r1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DIGITS counts the significant binary digits.'
  write ( *, '(a)' ) '  RANGE provides the decimal exponent range.'
  write ( *, '(a)' ) '  PRECISION provides the decimal precision.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '    DIGITS(R1) =    ', digits ( r1 )
  write ( *, '(a,i12)' ) '    RANGE(R1) =     ', range ( r1 )
  write ( *, '(a,i12)' ) '    PRECISION(R1) = ', precision ( r1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RADIX provides the base of the model.'
  write ( *, '(a)' ) '  MAXEXPONENT returns the maximum exponent of a variable.'
  write ( *, '(a)' ) '  MINEXPONENT returns the minimum exponent of a variable.'
  write ( *, '(a)' ) '  (These are exponents of the RADIX)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '    RADIX(R1) =       ', radix ( r1 )
  write ( *, '(a,i12)' ) '    MAXEXPONENT(R1) = ', maxexponent ( r1 )
  write ( *, '(a,i12)' ) '    MINEXPONENT(R1) = ', minexponent ( r1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  KIND returns the "kind" of a variable.'
  write ( *, '(a,i12)' ) '    KIND(R1) = ', kind ( r1 )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 shows what you get with KIND = 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 16 ) r1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Declare "real ( kind = 16 ) r1".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  HUGE returns the largest value of the given type.'
  write ( *, '(a)' ) '  TINY returns the smallest value of the given type.'
  write ( *, '(a)' ) '  EPSILON returns the precision of a real type.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '    HUGE(R1) =    ', huge ( r1 )
  write ( *, '(a,g14.6)' ) '    TINY(R1) =    ', tiny ( r1 )
  write ( *, '(a,g14.6)' ) '    EPSILON(R1) = ', epsilon ( r1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DIGITS counts the significant binary digits.'
  write ( *, '(a)' ) '  RANGE provides the decimal exponent range.'
  write ( *, '(a)' ) '  PRECISION provides the decimal precision.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '    DIGITS(R1) =    ', digits ( r1 )
  write ( *, '(a,i12)' ) '    RANGE(R1) =     ', range ( r1 )
  write ( *, '(a,i12)' ) '    PRECISION(R1) = ', precision ( r1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RADIX provides the base of the model.'
  write ( *, '(a)' ) '  MAXEXPONENT returns the maximum exponent of a variable.'
  write ( *, '(a)' ) '  MINEXPONENT returns the minimum exponent of a variable.'
  write ( *, '(a)' ) '  (These are exponents of the RADIX)'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '    RADIX(R1) =       ', radix ( r1 )
  write ( *, '(a,i12)' ) '    MAXEXPONENT(R1) = ', maxexponent ( r1 )
  write ( *, '(a,i12)' ) '    MINEXPONENT(R1) = ', minexponent ( r1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  KIND returns the "kind" of a variable.'
  write ( *, '(a,i12)' ) '    KIND(R1) = ', kind ( r1 )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 uses SELECTED_REAL_KIND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: k1 = selected_real_kind ( 5, 10 )
  real ( kind = k1 ) r1
  integer, parameter :: k2 = selected_real_kind ( 6, 70 )
  real ( kind = k2 ) r2
  integer, parameter :: k3 = selected_real_kind ( P = 14 )
  real ( kind = k3 ) r3
  integer, parameter :: k4 = selected_real_kind ( R = 4000 )
  real ( kind = k4 ) r4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Use SELECTED_REAL_KIND'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  k1 = selected_real_kind ( 5, 10 )'
  write ( *, '(a)' ) '  real ( kind = k1 ) r1'
  write ( *, '(a,i12)' ) '  KIND(R1) = ', kind ( r1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  k2 = selected_real_kind ( 6, 70 )'
  write ( *, '(a)' ) '  real ( kind = k2 ) r2'
  write ( *, '(a,i12)' ) '  KIND(R2) = ', kind ( r2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  k3 = selected_real_kind ( P = 14 )'
  write ( *, '(a)' ) '  real ( kind = k3 ) r3'
  write ( *, '(a,i12)' ) '  KIND(R3) = ', kind ( r3 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  k4 = selected_real_kind ( R = 4000 )'
  write ( *, '(a)' ) '  real ( kind = k4 ) r4'
  write ( *, '(a,i12)' ) '  KIND(R4) = ', kind ( r4 )
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
