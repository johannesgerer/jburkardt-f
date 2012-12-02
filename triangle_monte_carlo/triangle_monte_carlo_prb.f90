program main

!*****************************************************************************80
!
!! MAIN is the main program for TRIANGLE_MONTE_CARLO_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_MONTE_CARLO_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TRIANGLE_MONTE_CARLO library.'
!
!  Try each sampler on the unit triangle, integrating X^2, X*Y, Y^2.
!
  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Try each sampler on a general triangle, integrating a selection of functions.
!
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_MONTE_CARLO_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses TRIANGLE_SAMPLE_01 with an increasing number of points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 3

  integer ( kind = 4 ) p_num
  real ( kind = 8 ) result(f_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: t(2,3) = reshape ( (/ &
    1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00 /), (/ 2, 3 /) )
  external triangle_integrand_03
  external triangle_unit_sample_01

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Sample using TRIANGLE_UNIT_SAMPLE_01'
  write ( *, '(a)' ) '  Integrate TRIANGLE_UNIT_INTEGRAND_03'
  write ( *, '(a)' ) '  Integration region is the unit triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use an increasing number of points P_NUM.'
  write ( *, '(a)' ) '  Note that the sample routine is a "bad" sampler.'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P_NUM      X^2             X*Y             Y^2'
  write ( *, '(a)' ) ' '

  p_num = 1

  do while ( p_num <= 65536 )

    call triangle_monte_carlo ( t, p_num, f_num, triangle_unit_sample_01, &
      triangle_integrand_03, seed, result )

    write ( *, '(2x,i8,3(2x,g14.6))' ) p_num, result(1:f_num)

    p_num = 2 * p_num

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 uses TRIANGLE_SAMPLE_02 with an increasing number of points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 3

  integer ( kind = 4 ) p_num
  real ( kind = 8 ) result(f_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: t(2,3) = reshape ( (/ &
    1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00 /), (/ 2, 3 /) )
  external triangle_integrand_03
  external triangle_unit_sample_02

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Sample using TRIANGLE_UNIT_SAMPLE_02'
  write ( *, '(a)' ) '  Integrate TRIANGLE_UNIT_INTEGRAND_03'
  write ( *, '(a)' ) '  Integration region is the unit triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use an increasing number of points P_NUM.'
  write ( *, '(a)' ) '  Note that the sample routine is a good" sampler.'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P_NUM      X^2             X*Y             Y^2'
  write ( *, '(a)' ) ' '

  p_num = 1

  do while ( p_num <= 65536 )

    call triangle_monte_carlo ( t, p_num, f_num, triangle_unit_sample_02, &
      triangle_integrand_03, seed, result )

    write ( *, '(2x,i8,3(2x,g14.6))' ) p_num, result(1:f_num)

    p_num = 2 * p_num

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 uses TRIANGLE_SAMPLE_03 with an increasing number of points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 3

  integer ( kind = 4 ) p_num
  real ( kind = 8 ) result(f_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: t(2,3) = reshape ( (/ &
    1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00 /), (/ 2, 3 /) )
  external triangle_integrand_03
  external triangle_unit_sample_03

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Sample using TRIANGLE_UNIT_SAMPLE_03'
  write ( *, '(a)' ) '  Integrate TRIANGLE_UNIT_INTEGRAND_03'
  write ( *, '(a)' ) '  Integration region is the unit triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use an increasing number of points P_NUM.'
  write ( *, '(a)' ) '  Note that the sample routine is a good" sampler.'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P_NUM      X^2             X*Y             Y^2'
  write ( *, '(a)' ) ' '

  p_num = 1

  do while ( p_num <= 65536 )

    call triangle_monte_carlo ( t, p_num, f_num, triangle_unit_sample_03, &
      triangle_integrand_03, seed, result )

    write ( *, '(2x,i8,3(2x,g14.6))' ) p_num, result(1:f_num)

    p_num = 2 * p_num

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 uses TRIANGLE_SAMPLE_04 with an increasing number of points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 3

  integer ( kind = 4 ) p_num
  real ( kind = 8 ) result(f_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: t(2,3) = reshape ( (/ &
    1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00 /), (/ 2, 3 /) )
  external triangle_integrand_03
  external triangle_unit_sample_04

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Sample using TRIANGLE_UNIT_SAMPLE_04'
  write ( *, '(a)' ) '  Integrate TRIANGLE_UNIT_INTEGRAND_03'
  write ( *, '(a)' ) '  Integration region is the unit triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use an increasing number of points P_NUM.'
  write ( *, '(a)' ) '  Note that the sample routine is a good" sampler.'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P_NUM      X^2             X*Y             Y^2'
  write ( *, '(a)' ) ' '

  p_num = 1

  do while ( p_num <= 65536 )

    call triangle_monte_carlo ( t, p_num, f_num, triangle_unit_sample_04, &
      triangle_integrand_03, seed, result )

    write ( *, '(2x,i8,3(2x,g14.6))' ) p_num, result(1:f_num)

    p_num = 2 * p_num

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 uses TRIANGLE_SAMPLE_01 on a general triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 8

  integer ( kind = 4 ) p_num
  real ( kind = 8 ) result(f_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: t(2,3) = reshape ( (/ &
    4.0D+00, 1.0D+00, &
    8.0D+00, 3.0D+00, &
    0.0D+00, 9.0D+00 /), (/ 2, 3 /) )
  external triangle_integrand_user
  external triangle_unit_sample_01

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Sample using TRIANGLE_UNIT_SAMPLE_01'
  write ( *, '(a)' ) '  Integrate TRIANGLE_UNIT_INTEGRAND_USER'
  write ( *, '(a)' ) '  Integration region is over a general triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use an increasing number of points P_NUM.'
  write ( *, '(a)' ) '  Note that the sample routine is a "bad" sampler.'

  seed = 123456789

  call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P_NUM'
  write ( *, '(a)' ) ' '

  p_num = 1

  do while ( p_num <= 65536 )

    call triangle_monte_carlo ( t, p_num, f_num, triangle_unit_sample_01, &
      triangle_integrand_user, seed, result )

    write ( *, '(2x,i8,8(2x,f12.6))' ) p_num, result(1:f_num)

    p_num = 2 * p_num

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 uses TRIANGLE_SAMPLE_02 on a general triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 8

  integer ( kind = 4 ) p_num
  real ( kind = 8 ) result(f_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: t(2,3) = reshape ( (/ &
    4.0D+00, 1.0D+00, &
    8.0D+00, 3.0D+00, &
    0.0D+00, 9.0D+00 /), (/ 2, 3 /) )
  external triangle_integrand_user
  external triangle_unit_sample_02

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Sample using TRIANGLE_UNIT_SAMPLE_02'
  write ( *, '(a)' ) '  Integrate TRIANGLE_UNIT_INTEGRAND_USER'
  write ( *, '(a)' ) '  Integration region is over a general triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use an increasing number of points P_NUM.'
  write ( *, '(a)' ) '  Note that the sample routine is a "good" sampler.'

  seed = 123456789

  call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P_NUM'
  write ( *, '(a)' ) ' '

  p_num = 1

  do while ( p_num <= 65536 )

    call triangle_monte_carlo ( t, p_num, f_num, triangle_unit_sample_02, &
      triangle_integrand_user, seed, result )

    write ( *, '(2x,i8,8(2x,f12.6))' ) p_num, result(1:f_num)

    p_num = 2 * p_num

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 uses TRIANGLE_SAMPLE_03 on a general triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 8

  integer ( kind = 4 ) p_num
  real ( kind = 8 ) result(f_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: t(2,3) = reshape ( (/ &
    4.0D+00, 1.0D+00, &
    8.0D+00, 3.0D+00, &
    0.0D+00, 9.0D+00 /), (/ 2, 3 /) )
  external triangle_integrand_user
  external triangle_unit_sample_03

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  Sample using TRIANGLE_UNIT_SAMPLE_03'
  write ( *, '(a)' ) '  Integrate TRIANGLE_UNIT_INTEGRAND_USER'
  write ( *, '(a)' ) '  Integration region is over a general triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use an increasing number of points P_NUM.'
  write ( *, '(a)' ) '  Note that the sample routine is a "good" sampler.'

  seed = 123456789

  call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P_NUM'
  write ( *, '(a)' ) ' '

  p_num = 1

  do while ( p_num <= 65536 )

    call triangle_monte_carlo ( t, p_num, f_num, triangle_unit_sample_03, &
      triangle_integrand_user, seed, result )

    write ( *, '(2x,i8,8(2x,f12.6))' ) p_num, result(1:f_num)

    p_num = 2 * p_num

  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 uses TRIANGLE_SAMPLE_04 on a general triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: f_num = 8

  integer ( kind = 4 ) p_num
  real ( kind = 8 ) result(f_num)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) :: t(2,3) = reshape ( (/ &
    4.0D+00, 1.0D+00, &
    8.0D+00, 3.0D+00, &
    0.0D+00, 9.0D+00 /), (/ 2, 3 /) )
  external triangle_integrand_user
  external triangle_unit_sample_04

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  Sample using TRIANGLE_UNIT_SAMPLE_04'
  write ( *, '(a)' ) '  Integrate TRIANGLE_UNIT_INTEGRAND_USER'
  write ( *, '(a)' ) '  Integration region is over a general triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use an increasing number of points P_NUM.'
  write ( *, '(a)' ) '  Note that the sample routine is a "good" sampler.'

  seed = 123456789

  call r8mat_transpose_print ( 2, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P_NUM'
  write ( *, '(a)' ) ' '

  p_num = 1

  do while ( p_num <= 65536 )

    call triangle_monte_carlo ( t, p_num, f_num, triangle_unit_sample_04, &
      triangle_integrand_user, seed, result )

    write ( *, '(2x,i8,8(2x,f12.6))' ) p_num, result(1:f_num)

    p_num = 2 * p_num

  end do

  return
end
subroutine triangle_integrand_user ( p_num, p, f_num, fp )

!*****************************************************************************80
!
!! TRIANGLE_INTEGRAND_USER evaluates 8 integrand functions defined by the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P_NUM, the number of points.
!
!    Input, real ( kind = 8 ) P(2,P_NUM), the evaluation points.
!
!    Input, integer ( kind = 4 ) F_NUM, the number of integrands.
!
!    Output, real ( kind = 8 ) FP(F_NUM,P_NUM), the integrand values.
!
  implicit none

  integer ( kind = 4 ) f_num
  integer ( kind = 4 ) p_num

  real ( kind = 8 ) fp(f_num,p_num)
  real ( kind = 8 ) p(2,p_num)

  fp(1,1:p_num) = 1.0D+00
  fp(2,1:p_num) = p(1,1:p_num)
  fp(3,1:p_num) =                   p(2,1:p_num)
  fp(4,1:p_num) = p(1,1:p_num)**2
  fp(5,1:p_num) = p(1,1:p_num)    * p(2,1:p_num)
  fp(6,1:p_num) =                   p(2,1:p_num)**2
  fp(7,1:p_num) = p(1,1:p_num)**2 * p(2,1:p_num)
  fp(8,1:p_num) = p(1,1:p_num)**2 * p(2,1:p_num)**2

  return
end

