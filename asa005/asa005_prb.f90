program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA005_PRB.
!
!  Discussion:
!
!    ASA005_PRB calls the ASA005 routines.
!
!  Modified:
!
!    26 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA005_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA005 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA005_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates the use of PRNCST.
!
!  Modified:
!
!    26 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) df
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ifault
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) prncst
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  PRNCST computes the noncentral Student T '
  write ( *, '(a)' ) '  Cumulative Density Function.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '      X   LAMBDA  DF     ', &
    ' CDF                       CDF                     DIFF'
  write ( *, '(a,a)' ) '                         ', &
    ' Tabulated                 PRNCST'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call student_noncentral_cdf_values ( n_data, df, lambda, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = prncst ( x, df, lambda, ifault )

    write ( *, &
    '(2x,f6.2,2x,f6.2,2x,i2,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, lambda, df, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
