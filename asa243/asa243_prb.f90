program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA243_PRB.
!
!  Discussion:
!
!    ASA243_PRB calls the ASA243 routines.
!
!  Modified:
!
!    25 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA243_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA243 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA243_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates the use of TNC.
!
!  Modified:
!
!    25 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) delta
  integer ( kind = 4 ) df
  real ( kind = 8 ) df_real
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) tnc
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  TNC computes the noncentral Student T '
  write ( *, '(a)' ) '  Cumulative Density Function.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '        X         LAMBDA        DF     ', &
    ' CDF             CDF           DIFF'
  write ( *, '(a,a)' ) '                                       ', &
    ' Tabulated       PRNCST'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call student_noncentral_cdf_values ( n_data, df, delta, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    df_real = real ( df, kind = 8 )

    fx2 = tnc ( x, df_real, delta, ifault )

    write ( *, &
    '(2x,f10.4,2x,f10.4,2x,i8,2x,g14.6,2x,g14.6,2x,g10.4)' ) &
    x, delta, df, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
