program main

!*****************************************************************************80
!
!! MAIN is the main program for BETA_NC_PRB.
!
!  Discussion:
!
!    BETA_NC_PRB calls the BETA_NC routines.
!
!  Modified:
!
!    12 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BETA_NC_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BETA_NC library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BETA_NC_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BETA_NONCENTRAL_CDF against tabulated values.
!
!  Modified:
!
!    12 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) error_max
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  error_max = 1.0D-10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Compare tabulated values of the noncentral'
  write ( *, '(a)' ) '  incomplete Beta Function against values'
  write ( *, '(a)' ) '  computed by BETA_NONCENTRAL_CDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      A        B     LAMBDA        X       ', &
    ' CDF                         CDF                    DIFF'
  write ( *, '(a,a)' ) &
    '                                           ', &
    '(tabulated)       (BETA_NC)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_noncentral_cdf_values ( n_data, a, b, lambda, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call beta_noncentral_cdf ( a, b, lambda, x, error_max, fx2 )

    write ( *, &
    '(2x,f7.1,2x,f7.1,2x,f7.1,2x,f10.4, &
    2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    a, b, lambda, x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
