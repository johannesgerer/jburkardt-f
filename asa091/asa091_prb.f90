program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA091_PRB.
!
!  Modified:
!
!    24 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA091_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA091 library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA091_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 makes a single simple calculation with PPCHI2.
!
!  Modified:
!
!    24 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) alngam
  real ( kind = 8 ) g
  integer ( kind = 4 ) ifault
  real ( kind = 8 ) p
  real ( kind = 8 ) ppchi2
  real ( kind = 8 ) v
  real ( kind = 8 ) value
  real ( kind = 8 ), parameter :: value_correct = 0.4D+00

  p = 0.017523D+00
  v = 4.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Perform a simple sample calculation using'
  write ( *, '(a)' ) '  PPCHI2 to invert the Chi-Squared CDF.'

  g = alngam ( v / 2.0D+00, ifault )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g24.16)' ) '  P =                  ', p
  write ( *, '(a,g24.16)' ) '  V =                  ', v
  write ( *, '(a,g24.16)' ) '  G Log(Gamma(V/2)) =  ', g

  value = ppchi2 ( p, v, g, ifault )

  write ( *, '(a,g24.16)' ) '  VALUE =              ', value
  write ( *, '(a,g24.16)' ) '  VALUE (correct) =    ', value_correct

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'  ) '  Error flag IFAULT = ', ifault

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 compare PPCHI2 against tabulated values.
!
!  Modified:
!
!    24 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) alngam
  integer ( kind = 4 ) a
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) g
  real ( kind = 8 ) ppchi2
  integer ( kind = 4 ) ifault
  real ( kind = 8 ) v
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  Compare tabulated values of the Chi-Squared '
  write ( *, '(a)' ) '  Cumulative Density Function against values'
  write ( *, '(a)' ) '  computed by PPCHI2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
  '         N        CDF       X                        ', &
  ' X2                    DIFF'
  write ( *, '(a,a)' ) &
  '                           (tabulated)               ', &
  '(PPCHI2)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call chi_square_cdf_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    v = dble ( a )

    g = alngam ( v / 2.0D+00, ifault )

    x2 = ppchi2 ( fx, v, g, ifault )

    write ( *, '(2x,i8,2x,f10.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    a, fx, x, x2, abs ( x - x2 )

  end do

  return
end
