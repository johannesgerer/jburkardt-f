program main

!*****************************************************************************80
!
!! MAIN is the main program for OWENS_PRB.
!
!  Discussion:
!
!    OWENS_PRB calls the OWENS routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'OWENS_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the OWENS library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'OWENS_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates the use of T.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) h
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  T computes the Owen T function.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '          H               A           ', &
  'T                         T                       DIFF'
  write ( *, '(a,a)' ) '                                     ', &
  '(Tabulated)               (TFN)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call owen_values ( n_data, h, a, t1 )

    if ( n_data == 0 ) then
      exit
    end if

    t2 = t ( h, a )

    write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    h, a, t1, t2, abs ( t1 - t2 )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 demonstrates the use of BIVNOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 April 2012
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
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  BIVNOR computes the bivariate normal probability.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The value R is the correlation between X and Y.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
  '          X               Y               ', &
  'R           P                         P                       DIFF'
  write ( *, '(a,a)' ) &
  '                                          ', &
  '           (Tabulated)               (BIVNOR)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bivariate_normal_cdf_values ( n_data, x, y, r, fxy1 )

    if ( n_data == 0 ) then
      exit
    end if
!
!  BIVNOR computes the "tail" of the probability, and we want the
!  initial part!  To get that value, negate X and Y.
!
    fxy2 = bivnor (  - x, - y, r )

    write ( *, '(2x,f14.6,2x,f14.6,2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, y, r, fxy1, fxy2, abs ( fxy1 - fxy2 )

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 demonstrates the use of ZNORM1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ) znorm1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  ZNORM1 computes the normal CDF starting at 0.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
  '          X           P                         P                       DIFF'
  write ( *, '(a,a)' ) &
  '                     (Tabulated)               (ZNORM1)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_cdf_values ( n_data, x, fx1 )

    fx1 = fx1 - 0.5D+00

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = znorm1 ( x )

    write ( *, '(2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx1, fx2, abs ( fx1 - fx2 )

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 demonstrates the use of ZNORM2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 November 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx1
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ) znorm2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04:'
  write ( *, '(a)' ) '  ZNORM2 computes the complementary normal CDF.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
  '          X           P                         P                       DIFF'
  write ( *, '(a,a)' ) &
  '                     (Tabulated)               (ZNORM2)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_cdf_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx1 = 1.0D+00 - fx1
    fx2 = znorm2 ( x )

    write ( *, '(2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx1, fx2, abs ( fx1 - fx2 )

  end do

  return
end
