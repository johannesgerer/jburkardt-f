program main

!*****************************************************************************80
!
!! MAIN is the main program for MACHINE_PRB.
!
!  Discussion:
!
!    MACHINE_PRB runs the MACHINE tests.
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
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MACHINE_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the MACHINE library.'

  call d1mach_prb ( )
  call i1mach_prb ( )
  call r1mach_prb ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MACHINE_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine d1mach_prb ( )

!*****************************************************************************80
!
!! D1MACH_PRB reports the constants returned by D1MACH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  double precision d1mach

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'D1MACH_PRB'
  write ( *, '(a)' ) '  D1MACH reports the value of constants associated'
  write ( *, '(a)' ) '  with real double precision computer arithmetic.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Assume that double precision numbers are stored '
  write ( *, '(a)' ) '  with a mantissa of T digits in base B, with an '
  write ( *, '(a)' ) '  exponent whose value must lie between EMIN and EMAX.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For input arguments of 1 <= I <= 5,'
  write ( *, '(a)' ) '  D1MACH will return the following values:'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.'
  write ( *, * ) d1mach(1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  D1MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude.'
  write ( *, * ) d1mach(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  D1MACH(3) = B**(-T), the smallest relative spacing.'
  write ( *, * ) d1mach(3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  D1MACH(4) = B**(1-T), the largest relative spacing.'
  write ( *, * ) d1mach(4)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  D1MACH(5) = log10(B).'
  write ( *, * ) d1mach(5)

  return
end
subroutine i1mach_prb ( )

!*****************************************************************************80
!
!! I1MACH_PRB reports the constants returned by I1MACH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i1mach

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I1MACH_PRB'
  write ( *, '(a)' ) '  I1MACH reports the value of constants associated'
  write ( *, '(a)' ) '  with integer computer arithmetic.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Numbers associated with input/output units:'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I1MACH(1) = the standard input unit.'
  write ( *,     * ) i1mach(1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I1MACH(2) = the standard output unit.'
  write ( *,     * ) i1mach(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I1MACH(3) = the standard punch unit.'
  write ( *,     * ) i1mach(3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I1MACH(4) = the standard error message unit.'
  write ( *,     * ) i1mach(4)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Numbers associated with words:'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I1MACH(5) = the number of bits per integer.'
  write ( *,     * ) i1mach(5)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I1MACH(6) = the number of characters per integer.'
  write ( *,     * ) i1mach(6)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Numbers associated with integer values:'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Assume integers are represented in the S digit '
  write ( *, '(a)' ) '  base A form:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Sign * (X(S-1)*A**(S-1) + ... + X(1)*A + X(0))'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where the digits X satisfy 0 <= X(1:S-1) < A.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I1MACH(7) = A, the base.'
  write ( *,     * ) i1mach(7)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I1MACH(8) = S, the number of base A digits.'
  write ( *,     * ) i1mach(8)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I1MACH(9) = A**S-1, the largest integer.'
  write ( *,     * ) i1mach(9)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Numbers associated with floating point values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Assume floating point numbers are represented '
  write ( *, '(a)' ) '  in the T digit base B form:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Sign * (B**E) * ((X(1)/B) + ... + (X(T)/B**T) )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    0 <= X(1:T) < B,'
  write ( *, '(a)' ) '    0 < X(1) (unless the value being represented is 0),'
  write ( *, '(a)' ) '    EMIN <= E <= EMAX.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I1MACH(10) = B, the base.'
  write ( *,     * ) i1mach(10)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Numbers associated with single precision values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I1MACH(11) = T, the number of base B digits.'
  write ( *,     * ) i1mach(11)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I1MACH(12) = EMIN, the smallest exponent E.'
  write ( *,     * ) i1mach(12)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I1MACH(13) = EMAX, the largest exponent E.'
  write ( *,     * ) i1mach(13)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Numbers associated with double precision values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I1MACH(14) = T, the number of base B digits.'
  write ( *,     * ) i1mach(14)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I1MACH(15) = EMIN, the smallest exponent E.'
  write ( *,     * ) i1mach(15)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I1MACH(16) = EMAX, the largest exponent E.'
  write ( *,     * ) i1mach(16)

  return
end
subroutine r1mach_prb ( )

!*****************************************************************************80
!
!! R1MACH_PRB reports the constants returned by R1MACH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real r1mach

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R1MACH_PRB'
  write ( *, '(a)' ) '  R1MACH reports the value of constants associated'
  write ( *, '(a)' ) '  with real single precision computer arithmetic.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Assume that single precision numbers are stored '
  write ( *, '(a)' ) '  with a mantissa of T digits in base B, with an '
  write ( *, '(a)' ) '  exponent whose value must lie between EMIN and EMAX.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For input arguments of 1 <= I <= 5,'
  write ( *, '(a)' ) '  R1MACH will return the following values:'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.'
  write ( *, * ) r1mach(1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  R1MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude.'
  write ( *, * ) r1mach(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  R1MACH(3) = B**(-T), the smallest relative spacing.'
  write ( *, * ) r1mach(3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  R1MACH(4) = B**(1-T), the largest relative spacing.'
  write ( *, * ) r1mach(4)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R1MACH(5) = log10(B).'
  write ( *, * ) r1mach(5)

  return
end
