program main

!*****************************************************************************80
!
!! MAIN is the main program for MACHAR_PRB.
!
!  Discussion:
!
!    MACHAR_PRB tests the MACHAR library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MACHAR_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the MACHAR library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MACHAR_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests R4_MACHAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) eps
  real ( kind = 4 ) epsneg
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  R4_MACHAR computes single'
  write ( *, '(a)' ) '  precision machine constants.'

  call r4_machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, &
    minexp, maxexp, eps, epsneg, xmin, xmax )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  IBETA is the internal base for machine arithmetic.'
  write ( *, '(a,i8)' ) '    IBETA =  ', ibeta
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  IT is the number of digits, base IBETA, in the'
  write ( *, '(a)' ) '  floating point significand.'
  write ( *, '(a,i8)' ) '    IT =     ', it
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  IRND reports on floating point addition rounding:'
  write ( *, '(a)' ) '  0, for chopping;'
  write ( *, '(a)' ) '  1, for non-IEEE rounding;'
  write ( *, '(a)' ) '  2, for IEEE rounding;'
  write ( *, '(a)' ) '  3, for chopping with partial underflow;'
  write ( *, '(a)' ) '  4, for non-IEEE rounding with partial underflow.'
  write ( *, '(a)' ) '  5, for IEEE rounding with partial underflow.'
  write ( *, '(a,i8)' ) '    IRND =   ', irnd
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NGRD is the number of guard digits for floating point'
  write ( *, '(a)' ) '  multiplication with truncating arithmetic.'
  write ( *, '(a,i8)' ) '    NGRD =   ', ngrd
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  MACHEP is the largest negative integer such that'
  write ( *, '(a)' ) '  1.0 < 1.0 + BETA^MACHEP.'
  write ( *, '(a,i8)' ) '    MACHEP = ', machep
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NEGEPS is the largest negative integer such that'
  write ( *, '(a)' ) '  1.0 - BETA^NEGEPS < 1.0:'
  write ( *, '(a,i8)' ) '    NEGEP =  ', negep
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  IEXP is the number of bits reserved for the exponent'
  write ( *, '(a)' ) '  of a floating point number:'
  write ( *, '(a,i8)' ) '    IEXP =   ', iexp
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  MINEXP is the most negative power of BETA such that'
  write ( *, '(a)' ) '  BETA^MINEXP is positive and normalized.'
  write ( *, '(a,i8)' ) '    MINEXP = ', minexp
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  MAXEXP is the smallest positive power of BETA that'
  write ( *, '(a)' ) '  overflows:'
  write ( *, '(a,i8)' ) '    MAXEXP = ', maxexp
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EPS is a small positive floating point number'
  write ( *, '(a)' ) '  such that 1.0 < 1.0 + EPS.'
  write ( *, '(a,e26.16)' ) '    EPS    = ', eps
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EPSNEG is a small positive floating point number'
  write ( *, '(a)' ) '  such that 1.0 - EPSNEG < 1.0.'
  write ( *, '(a,e26.16)' ) '    EPSNEG = ', epsneg
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XMIN is the smallest positive normalized floating'
  write ( *, '(a)' ) '  point power of the radix:'
  write ( *, '(a,e26.16)' ) '    XMIN =   ', xmin
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XMAX is the largest finite floating point number:'
  write ( *, '(a,e26.16)' ) '    XMAX   = ', xmax

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat floating point data using * format:'
  write ( *, '(a)' ) ' '
  write ( *, * ) '    EPS    = ', eps
  write ( *, * ) '    EPSNEG = ', epsneg
  write ( *, * ) '    XMIN   = ', xmin
  write ( *, * ) '    XMAX   = ', xmax

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests R8_MACHAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  integer ( kind = 8 ) ibeta
  integer ( kind = 8 ) iexp
  integer ( kind = 8 ) irnd
  integer ( kind = 8 ) it
  integer ( kind = 8 ) machep
  integer ( kind = 8 ) maxexp
  integer ( kind = 8 ) minexp
  integer ( kind = 8 ) negep
  integer ( kind = 8 ) ngrd
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  R8_MACHAR computes double'
  write ( *, '(a)' ) '  precision machine constants.'

  call r8_machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, &
    minexp, maxexp, eps, epsneg, xmin, xmax )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  IBETA is the internal base for machine arithmetic.'
  write ( *, '(a,i8)' ) '    IBETA =  ', ibeta
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  IT is the number of digits, base IBETA, in the'
  write ( *, '(a)' ) '  floating point significand.'
  write ( *, '(a,i8)' ) '    IT =     ', it
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  IRND reports on floating point addition rounding:'
  write ( *, '(a)' ) '  0, for chopping;'
  write ( *, '(a)' ) '  1, for non-IEEE rounding;'
  write ( *, '(a)' ) '  2, for IEEE rounding;'
  write ( *, '(a)' ) '  3, for chopping with partial underflow;'
  write ( *, '(a)' ) '  4, for non-IEEE rounding with partial underflow.'
  write ( *, '(a)' ) '  5, for IEEE rounding with partial underflow.'
  write ( *, '(a,i8)' ) '    IRND =   ', irnd
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NGRD is the number of guard digits for floating point'
  write ( *, '(a)' ) '  multiplication with truncating arithmetic.'
  write ( *, '(a,i8)' ) '    NGRD =   ', ngrd
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  MACHEP is the largest negative integer such that'
  write ( *, '(a)' ) '  1.0 < 1.0 + BETA^MACHEP.'
  write ( *, '(a,i8)' ) '    MACHEP = ', machep
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NEGEPS is the largest negative integer such that'
  write ( *, '(a)' ) '  1.0 - BETA^NEGEPS < 1.0:'
  write ( *, '(a,i8)' ) '    NEGEP =  ', negep
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  IEXP is the number of bits reserved for the exponent'
  write ( *, '(a)' ) '  of a floating point number:'
  write ( *, '(a,i8)' ) '    IEXP =   ', iexp
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  MINEXP is the most negative power of BETA such that'
  write ( *, '(a)' ) '  BETA^MINEXP is positive and normalized.'
  write ( *, '(a,i8)' ) '    MINEXP = ', minexp
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  MAXEXP is the smallest positive power of BETA that'
  write ( *, '(a)' ) '  overflows:'
  write ( *, '(a,i8)' ) '    MAXEXP = ', maxexp
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EPS is a small positive floating point number'
  write ( *, '(a)' ) '  such that 1.0 < 1.0 + EPS.'
  write ( *, '(a,e26.16)' ) '    EPS    = ', eps
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EPSNEG is a small positive floating point number'
  write ( *, '(a)' ) '  such that 1.0 - EPSNEG < 1.0.'
  write ( *, '(a,e26.16)' ) '    EPSNEG = ', epsneg
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XMIN is the smallest positive normalized floating'
  write ( *, '(a)' ) '  point power of the radix:'
  write ( *, '(a,e26.16)' ) '    XMIN =   ', xmin
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  XMAX is the largest finite floating point number:'
  write ( *, '(a,e26.16)' ) '    XMAX   = ', xmax
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat floating point data using * format:'
  write ( *, '(a)' ) ' '
  write ( *, * ) '    EPS    = ', eps
  write ( *, * ) '    EPSNEG = ', epsneg
  write ( *, * ) '    XMIN   = ', xmin
  write ( *, * ) '    XMAX   = ', xmax

  return
end
