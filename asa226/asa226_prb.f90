program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA226_PRB.
!
!  Discussion:
!
!    ASA226_PRB calls the ASA226 routines.
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
  write ( *, '(a)' ) 'ASA226_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA226 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA226_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BETANC against tabulated values.
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

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) betanc
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ifault
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Compare tabulated values of the noncentral'
  write ( *, '(a)' ) '  incomplete Beta Function against values'
  write ( *, '(a)' ) '  computed by BETANC.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      A        B     LAMBDA        X       ', &
    ' CDF               CDF          DIFF'
  write ( *, '(a,a)' ) &
    '                                           ', &
    '(tabulated)       (BETANC)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_noncentral_cdf_values ( n_data, a, b, lambda, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = betanc ( x, a, b, lambda, ifault )

    write ( *, &
    '(2x,f7.1,2x,f7.1,2x,f7.1,2x,f10.4, &
    2x,g14.6,2x,g14.6,2x,g10.4)' ) &
    a, b, lambda, x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
