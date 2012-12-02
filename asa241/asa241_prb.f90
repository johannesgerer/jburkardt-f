program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA241_PRB.
!
!  Discussion:
!
!    ASA241_PRB tests ASA241.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA241_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA241 library.'

  call test01 ( )
  call test02 ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA241_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests R4_NORMAL_01_CDF_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 4 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_normal_01_cdf_inverse
  real ( kind = 8 ) x
  real ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Let FX = NormalCDF ( X ).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NORMAL_01_CDF_VALUES returns some values of ( X, FX ).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R4_NORMAL_01_CDF_INVERSE takes the value of FX, and'
  write ( *, '(a)' ) '  computes an estimate X2, of the corresponding input,'
  write ( *, '(a)' ) '  argument, accurate to about 7 decimal places.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      FX                         X                         X2' &
    // '                     DIFF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_cdf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = real ( fx, kind = 4 )
    x2 = r4_normal_01_cdf_inverse ( fx2 )

    write ( *, '(2x,g24.16,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    fx, x, x2, abs ( x - x2 )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests R8_NORMAL_01_CDF_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) r8_normal_01_cdf_inverse
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  Let FX = NormalCDF ( X ).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NORMAL_01_CDF_VALUES returns some values of ( X, FX ).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8_NORMAL_01_CDF_INVERSE takes the value of FX, and'
  write ( *, '(a)' ) '  computes an estimate X2, of the corresponding input, '
  write ( *, '(a)' ) '  argument, accurate to about 16 decimal places.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) &
    '      FX                         X                         X2' &
    // '                     DIFF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_cdf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = r8_normal_01_cdf_inverse ( fx )

    write ( *, '(2x,g24.16,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    fx, x, x2, abs ( x - x2 )

  end do

  return
end
