program main

!*****************************************************************************80
!
!! TOMS462_PRB tests BIVNOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS462_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS462 library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS462_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BIVNOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bivnor
  real ( kind = 8 ) cdf
  real ( kind = 8 ) expect
  real ( kind = 8 ) r
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Compare BIVNOR with some simple data'
  write ( *, '(a)' ) '  with 3 digit accuracy.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '       X         Y          R          P               P'
  write ( *, '(a)' ) &
    '                                      (Tabulated)     (BIVNOR)'
  write ( *, '(a)' ) ' '

  x =  0.8D+00
  y = -1.5D+00
  r =  -0.9D+00
  expect = 0.148D+00
  cdf = bivnor ( x, y, r )
  write ( *, '(2x,f9.5,2x,f9.5,2x,f9.5,2x,g14.6,2x,g14.6)' ) &
    x, y, r, expect, cdf

stop
  x =  0.6D+00
  y = -1.4D+00
  r =  -0.7D+00
  expect = 0.208D+00
  cdf = bivnor ( x, y, r )
  write ( *, '(2x,f9.5,2x,f9.5,2x,f9.5,2x,g14.6,2x,g14.6)' ) &
    x, y, r, expect, cdf

  x =  0.2D+00
  y = -1.0D+00
  r =  -0.5D+00
  expect = 0.304D+00
  cdf = bivnor ( x, y, r )
  write ( *, '(2x,f9.5,2x,f9.5,2x,f9.5,2x,g14.6,2x,g14.6)' ) &
    x, y, r, expect, cdf

  x = -1.2D+00
  y =  0.1D+00
  r =   0.0D+00
  expect = 0.407D+00
  cdf = bivnor ( x, y, r )
  write ( *, '(2x,f9.5,2x,f9.5,2x,f9.5,2x,g14.6,2x,g14.6)' ) &
    x, y, r, expect, cdf

  x = -1.2D+00
  y = -0.1D+00
  r =   0.3D+00
  expect = 0.501D+00
  cdf = bivnor ( x, y, r )
  write ( *, '(2x,f9.5,2x,f9.5,2x,f9.5,2x,g14.6,2x,g14.6)' ) &
    x, y, r, expect, cdf

  x = -0.4D+00
  y = -0.9D+00
  r =   0.6D+00
  expect = 0.601D+00
  cdf = bivnor ( x, y, r )
  write ( *, '(2x,f9.5,2x,f9.5,2x,f9.5,2x,g14.6,2x,g14.6)' ) &
    x, y, r, expect, cdf

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests BIVNOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bivnor
  real ( kind = 8 ) fxy1
  real ( kind = 8 ) fxy2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) r
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Compare BIVNOR with some tabulated data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a,a)' ) &
    '      X          Y          ', &
    'R           P                         P', &
    '                      DIFF'
  write ( *, '(a,a)' ) &
    '                                ', &
    '       (Tabulated)               (BIVNOR)' 
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bivariate_normal_cdf_values ( n_data, x, y, r, fxy1 )

    if ( n_data == 0 ) then
      exit
    end if
!
!  BIVNOR computes the "tail" of the probability, and we want the
!  initial part!
!
    x = - x
    y = - y
    fxy2 = bivnor ( x, y, r )
    x = - x
    y = - y

    write ( *, &
      '(2x,f8.4,2x,f8.4,2x,f8.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
       x, y, r, fxy1, fxy2, abs ( fxy1 - fxy2 )

  end do

  return
end
