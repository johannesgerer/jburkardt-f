program main

!*****************************************************************************80
!
!! MAIN is the main program for TETRAHEDRON_MONTE_CARLO_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TETRAHEDRON_MONTE_CARLO_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TETRAHEDRON_MONTE_CARLO library.'
!
!  Try each sampler on the unit tetrahedron, integrating quadratics.
!
  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Try each sampler on a general tetrahedron, integrating a selection of functions.
!
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TETRAHEDRON_MONTE_CARLO_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses TETRAHEDRON_SAMPLE_01 with an increasing number of points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 6

  integer ( kind = 4 ) p_num
  real ( kind = 8 ) result(f_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: t(3,4) = reshape ( (/ &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, 0.0D+00 /), (/ 3, 4 /) )
  external tetrahedron_integrand_03
  external tetrahedron_unit_sample_01

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Sample using TETRAHEDRON_UNIT_SAMPLE_01'
  write ( *, '(a)' ) '  Integrate TETRAHEDRON_UNIT_INTEGRAND_03'
  write ( *, '(a)' ) '  Integration region is the unit tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use an increasing number of points P_NUM.'
  write ( *, '(a)' ) '  Note that the sample routine is a bad sampler.'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P_NUM      X^2             X*Y             X*Z' // &
    '             Y^2             Y*Z             Z^2'
  write ( *, '(a)' ) ' '

  p_num = 1

  do while ( p_num <= 65536 )

    call tetrahedron_monte_carlo ( t, p_num, f_num, tetrahedron_unit_sample_01, &
      tetrahedron_integrand_03, seed, result )

    write ( *, '(2x,i8,6(2x,g14.6))' ) p_num, result(1:f_num)

    p_num = 2 * p_num

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 uses TETRAHEDRON_SAMPLE_02 with an increasing number of points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 6

  integer ( kind = 4 ) p_num
  real ( kind = 8 ) result(f_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: t(3,4) = reshape ( (/ &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, 0.0D+00 /), (/ 3, 4 /) )
  external             tetrahedron_integrand_03
  external             tetrahedron_unit_sample_02

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Sample using TETRAHEDRON_UNIT_SAMPLE_02'
  write ( *, '(a)' ) '  Integrate TETRAHEDRON_UNIT_INTEGRAND_03'
  write ( *, '(a)' ) '  Integration region is the unit tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use an increasing number of points P_NUM.'
  write ( *, '(a)' ) '  Note that the sample routine is a good sampler.'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P_NUM      X^2             X*Y             X*Z' // &
    '             Y^2             Y*Z             Z^2'
  write ( *, '(a)' ) ' '

  p_num = 1

  do while ( p_num <= 65536 )

    call tetrahedron_monte_carlo ( t, p_num, f_num, tetrahedron_unit_sample_02, &
      tetrahedron_integrand_03, seed, result )

    write ( *, '(2x,i8,6(2x,g14.6))' ) p_num, result(1:f_num)

    p_num = 2 * p_num

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 uses TETRAHEDRON_SAMPLE_03 with an increasing number of points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 6

  integer ( kind = 4 ) p_num
  real ( kind = 8 ) result(f_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: t(3,4) = reshape ( (/ &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, 0.0D+00 /), (/ 3, 4 /) )
  external             tetrahedron_integrand_03
  external             tetrahedron_unit_sample_03

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Sample using TETRAHEDRON_UNIT_SAMPLE_03'
  write ( *, '(a)' ) '  Integrate TETRAHEDRON_UNIT_INTEGRAND_03'
  write ( *, '(a)' ) '  Integration region is the unit tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use an increasing number of points P_NUM.'
  write ( *, '(a)' ) '  Note that the sample routine is a good sampler.'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P_NUM      X^2             X*Y             X*Z' // &
    '             Y^2             Y*Z             Z^2'
  write ( *, '(a)' ) ' '

  p_num = 1

  do while ( p_num <= 65536 )

    call tetrahedron_monte_carlo ( t, p_num, f_num, tetrahedron_unit_sample_03, &
      tetrahedron_integrand_03, seed, result )

    write ( *, '(2x,i8,6(2x,g14.6))' ) p_num, result(1:f_num)

    p_num = 2 * p_num

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 uses TETRAHEDRON_SAMPLE_04 with an increasing number of points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 6

  integer ( kind = 4 ) p_num
  real ( kind = 8 ) result(f_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: t(3,4) = reshape ( (/ &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, 0.0D+00 /), (/ 3, 4 /) )
  external             tetrahedron_integrand_03
  external             tetrahedron_unit_sample_04

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Sample using TETRAHEDRON_UNIT_SAMPLE_04'
  write ( *, '(a)' ) '  Integrate TETRAHEDRON_UNIT_INTEGRAND_03'
  write ( *, '(a)' ) '  Integration region is the unit tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use an increasing number of points P_NUM.'
  write ( *, '(a)' ) '  Note that the sample routine is a good sampler.'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P_NUM      X^2             X*Y             X*Z' // &
    '             Y^2             Y*Z             Z^2'
  write ( *, '(a)' ) ' '

  p_num = 1

  do while ( p_num <= 65536 )

    call tetrahedron_monte_carlo ( t, p_num, f_num, tetrahedron_unit_sample_04, &
      tetrahedron_integrand_03, seed, result )

    write ( *, '(2x,i8,6(2x,g14.6))' ) p_num, result(1:f_num)

    p_num = 2 * p_num

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 uses TETRAHEDRON_SAMPLE_01 on a general tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 6

  integer ( kind = 4 ) p_num
  real ( kind = 8 ) result(f_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: t(3,4) = reshape ( (/ &
    1.0D+00, 2.0D+00, 3.0D+00, &
    4.0D+00, 1.0D+00, 2.0D+00, &
    2.0D+00, 4.0D+00, 4.0D+00, &
    3.0D+00, 2.0D+00, 5.0D+00 /), (/ 3, 4 /) )
  external             tetrahedron_integrand_user
  external             tetrahedron_unit_sample_01

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Sample using TETRAHEDRON_UNIT_SAMPLE_01'
  write ( *, '(a)' ) '  Integrate TETRAHEDRON_UNIT_INTEGRAND_USER'
  write ( *, '(a)' ) '  Integration region is over a general tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use an increasing number of points P_NUM.'
  write ( *, '(a)' ) '  Note that the sample routine is a bad sampler.'

  seed = 123456789

  call r8mat_transpose_print ( 3, 4, t, '  Tetrahedron vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P_NUM'
  write ( *, '(a)' ) ' '

  p_num = 1

  do while ( p_num <= 65536 )

    call tetrahedron_monte_carlo ( t, p_num, f_num, tetrahedron_unit_sample_01, &
      tetrahedron_integrand_user, seed, result )

    write ( *, '(2x,i8,6(2x,f12.6))' ) p_num, result(1:f_num)

    p_num = 2 * p_num

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 uses TETRAHEDRON_SAMPLE_02 on a general tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 6

  integer ( kind = 4 ) p_num
  real ( kind = 8 ) result(f_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: t(3,4) = reshape ( (/ &
    1.0D+00, 2.0D+00, 3.0D+00, &
    4.0D+00, 1.0D+00, 2.0D+00, &
    2.0D+00, 4.0D+00, 4.0D+00, &
    3.0D+00, 2.0D+00, 5.0D+00 /), (/ 3, 4 /) )
  external             tetrahedron_integrand_user
  external             tetrahedron_unit_sample_02

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Sample using TETRAHEDRON_UNIT_SAMPLE_02'
  write ( *, '(a)' ) '  Integrate TETRAHEDRON_UNIT_INTEGRAND_USER'
  write ( *, '(a)' ) '  Integration region is over a general tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use an increasing number of points P_NUM.'
  write ( *, '(a)' ) '  Note that the sample routine is a good sampler.'

  seed = 123456789

  call r8mat_transpose_print ( 3, 4, t, '  Tetrahedron vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P_NUM'
  write ( *, '(a)' ) ' '

  p_num = 1

  do while ( p_num <= 65536 )

    call tetrahedron_monte_carlo ( t, p_num, f_num, tetrahedron_unit_sample_02, &
      tetrahedron_integrand_user, seed, result )

    write ( *, '(2x,i8,6(2x,f12.6))' ) p_num, result(1:f_num)

    p_num = 2 * p_num

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 uses TETRAHEDRON_SAMPLE_03 on a general tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 6

  integer ( kind = 4 ) p_num
  real ( kind = 8 ) result(f_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: t(3,4) = reshape ( (/ &
    1.0D+00, 2.0D+00, 3.0D+00, &
    4.0D+00, 1.0D+00, 2.0D+00, &
    2.0D+00, 4.0D+00, 4.0D+00, &
    3.0D+00, 2.0D+00, 5.0D+00 /), (/ 3, 4 /) )
  external             tetrahedron_integrand_user
  external             tetrahedron_unit_sample_03

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  Sample using TETRAHEDRON_UNIT_SAMPLE_03'
  write ( *, '(a)' ) '  Integrate TETRAHEDRON_UNIT_INTEGRAND_USER'
  write ( *, '(a)' ) '  Integration region is over a general tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use an increasing number of points P_NUM.'
  write ( *, '(a)' ) '  Note that the sample routine is a good sampler.'

  seed = 123456789

  call r8mat_transpose_print ( 3, 4, t, '  Tetrahedron vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P_NUM'
  write ( *, '(a)' ) ' '

  p_num = 1

  do while ( p_num <= 65536 )

    call tetrahedron_monte_carlo ( t, p_num, f_num, tetrahedron_unit_sample_03, &
      tetrahedron_integrand_user, seed, result )

    write ( *, '(2x,i8,6(2x,f12.6))' ) p_num, result(1:f_num)

    p_num = 2 * p_num

  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 uses TETRAHEDRON_SAMPLE_04 on a general tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 6

  integer ( kind = 4 ) p_num
  real ( kind = 8 ) result(f_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: t(3,4) = reshape ( (/ &
    1.0D+00, 2.0D+00, 3.0D+00, &
    4.0D+00, 1.0D+00, 2.0D+00, &
    2.0D+00, 4.0D+00, 4.0D+00, &
    3.0D+00, 2.0D+00, 5.0D+00 /), (/ 3, 4 /) )
  external             tetrahedron_integrand_user
  external             tetrahedron_unit_sample_04

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  Sample using TETRAHEDRON_UNIT_SAMPLE_04'
  write ( *, '(a)' ) '  Integrate TETRAHEDRON_UNIT_INTEGRAND_USER'
  write ( *, '(a)' ) '  Integration region is over a general tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use an increasing number of points P_NUM.'
  write ( *, '(a)' ) '  Note that the sample routine is a good sampler.'

  seed = 123456789

  call r8mat_transpose_print ( 3, 4, t, '  Tetrahedron vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P_NUM'
  write ( *, '(a)' ) ' '

  p_num = 1

  do while ( p_num <= 65536 )

    call tetrahedron_monte_carlo ( t, p_num, f_num, tetrahedron_unit_sample_04, &
      tetrahedron_integrand_user, seed, result )

    write ( *, '(2x,i8,6(2x,f12.6))' ) p_num, result(1:f_num)

    p_num = 2 * p_num

  end do

  return
end
subroutine tetrahedron_integrand_user ( p_num, p, f_num, fp )

!*****************************************************************************80
!
!! TETRAHEDRON_INTEGRAND_USER evaluates 6 integrand functions defined by the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P_NUM, the number of points.
!
!    Input, real ( kind = 8 ) P(3,P_NUM), the evaluation points.
!
!    Input, integer ( kind = 4 ) F_NUM, the number of integrands.
!
!    Output, real ( kind = 8 ) FP(F_NUM,P_NUM), the integrand values.
!
  implicit none

  integer ( kind = 4 ) f_num
  integer ( kind = 4 ) p_num

  real ( kind = 8 ) fp(f_num,p_num)
  real ( kind = 8 ) p(3,p_num)

  fp(1,1:p_num) = 1.0D+00
  fp(2,1:p_num) = p(1,1:p_num)
  fp(3,1:p_num) =                   p(2,1:p_num)**2
  fp(4,1:p_num) =                                     p(3,1:p_num)**3
  fp(5,1:p_num) = p(1,1:p_num)    * p(2,1:p_num)    * p(3,1:p_num)**2
  fp(6,1:p_num) = p(1,1:p_num)**2 * p(2,1:p_num)**2 * p(3,1:p_num)**2

  return
end

