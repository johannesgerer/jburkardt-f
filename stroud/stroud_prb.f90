program main

!*****************************************************************************80
!
!! MAIN is the main program for STROUD_PRB.
!
!  Discussion:
!
!    STROUD_PRB tests the routines in the STROUD library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STROUD_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the STROUD library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test045 ( )
  call test05 ( )
  call test052 ( )
  call test054 ( )
  call test07 ( )
  call test08 ( )
  call test085 ( )
  call test09 ( )

  call test10 ( )
  call test11 ( )
  call test12 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
  call test16 ( )
  call test163 ( )
  call test165 ( )
  call test167 ( )
  call test17 ( )
  call test18 ( )
  call test19 ( )

  call test20 ( )
  call test205 ( )
  call test207 ( )
  call test2075 ( )
  call test208 ( )
  call test21 ( )
  call test215 ( )
  call test22 ( )
  call test23 ( )
  call test24 ( )
  call test25 ( )
  call test255 ( )
  call test26 ( )
  call test27 ( )
  call test28 ( )
  call test29 ( )

  call test30 ( )
  call test31 ( )
  call test32 ( )
  call test322 ( )
  call test324 ( )
  call test326 ( )
  call test33 ( )
  call test335 ( )
  call test34 ( )
  call test345 ( )
  call test35 ( )
  call test36 ( )
  call test37 ( )
  call test38 ( )
  call test39 ( )

  call test40 ( )
  call test41 ( )
  call test42 ( )
  call test425 ( )
  call test43 ( )
  call test44 ( )
  call test45 ( )
  call test46 ( )
  call test47 ( )
  call test48 ( )
  call test49 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STROUD_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BALL_F1_ND, BALL_F3_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 3

  real ( kind = 8 ) ball_volume_nd
  real ( kind = 8 ), dimension ( n_max ) :: center = (/ &
    1.0D+00, -1.0D+00, 2.0D+00 /)
  real ( kind = 8 ), external :: function_nd
  character ( len = 7 ) function_nd_name
  integer ( kind = 4 ) function_nd_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) num
  real ( kind = 8 ) result1
  real ( kind = 8 ) result2
  real ( kind = 8 ) r

  r = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For integrals in a ball in ND:'
  write ( *, '(a)' ) '  BALL_F1_ND approximates the integral;'
  write ( *, '(a)' ) '  BALL_F3_ND approximates the integral.'
  write ( *, '(a)' ) ' '

  do n = 2, n_max

    do i = 1, n
      center(i) = real ( i, kind = 8 )
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Spatial dimension N = ', n
    write ( *, '(a)' ) '  Ball center:'
    write ( *, '(2x,3f10.4)' ) center(1:n)
    write ( *, '(a,g14.6)') '  Ball radius = ', r
    write ( *, '(a,g14.6)' ) '  Ball volume = ', ball_volume_nd ( n, r )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    Rule:      F1          F3'
    write ( *, '(a)' ) '    F(X)'
    write ( *, '(a)' ) ' '

    num = function_nd_num ( )

    do i = 1, num

      call function_nd_set ( 'SET', i )

      call ball_f1_nd ( function_nd, n, center, r, result1 )
      call ball_f3_nd ( function_nd, n, center, r, result2 )
      write ( *, '(2x,a7,2f14.8)' ) function_nd_name(i), result1, result2

    end do

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests BALL_MONOMIAL_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) ball_monomial_nd
  real ( kind = 8 ) ball_volume_nd
  real ( kind = 8 ), dimension ( dim_num ) :: center = (/ &
    0.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), external :: mono_000_3d
  real ( kind = 8 ), external :: mono_111_3d
  real ( kind = 8 ), external :: mono_202_3d
  real ( kind = 8 ), external :: mono_422_3d
  integer ( kind = 4 ) p(dim_num)
  real ( kind = 8 ) result1
  real ( kind = 8 ) result2
  real ( kind = 8 ) :: r = 2.0D+00
  character ( len = 10 ) string
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  For the integral of a monomial in a ball in ND:'
  write ( *, '(a)' ) '  BALL_MONOMIAL_ND approximates the integral.'
  write ( *, '(a)' ) '  BALL_F1_ND, which can handle general integrands,'
  write ( *, '(a)' ) '    will be used for comparison.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension N = ', dim_num
  write ( *, '(a,g14.6)' ) '  Ball radius = ', r
  write ( *, '(a,g14.6)' ) '  Ball volume = ', &
    ball_volume_nd ( dim_num, r )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Rule:     MONOMIAL    F1'
  write ( *, '(a)' ) '    F(X)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    if ( test == 1 ) then
      string = '1'
      p = (/ 0, 0, 0 /)
      call ball_f1_nd ( mono_000_3d, dim_num, center, r, result2 )
    else if ( test == 2 ) then
      string = 'xyz'
      p = (/ 1, 1, 1 /)
      call ball_f1_nd ( mono_111_3d, dim_num, center, r, result2 )
    else if ( test == 3 ) then
      string = 'x^2z^2'
      p = (/ 2, 0, 2 /)
      call ball_f1_nd ( mono_202_3d, dim_num, center, r, result2 )
    else if ( test == 4 ) then
      string = 'x^4y^2z^2'
      p = (/ 4, 2, 2 /)
      call ball_f1_nd ( mono_422_3d, dim_num, center, r, result2 )
    end if

    result1 = ball_monomial_nd ( dim_num, p, r )

    write ( *, '(2x,a10,2f14.8)' ) string, result1, result2

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests BALL_UNIT_**_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) ball_unit_volume_nd
  character ( len = 7 ) function_3d_name
  integer ( kind = 4 ) function_3d_num
  real ( kind = 8 ), external :: function_3d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) num
  real ( kind = 8 ) result1
  real ( kind = 8 ) result2
  real ( kind = 8 ) result3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For integrals in the unit ball in 3D:'
  write ( *, '(a)' ) '  BALL_UNIT_07_3D uses a formula of degree 7;'
  write ( *, '(a)' ) '  BALL_UNIT_14_3D uses a formula of degree 14;'
  write ( *, '(a)' ) '  BALL_UNIT_15_3D uses a formula of degree 15.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Unit ball volume = ', ball_unit_volume_nd ( 3 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Rule:      #7             #14           #15'
  write ( *, '(a)' ) '    F(X)'
  write ( *, '(a)' ) ' '
 
  num = function_3d_num ( )

  do i = 1, num

    call function_3d_set ( 'SET', i )

    call ball_unit_07_3d ( function_3d, result1 )
    call ball_unit_14_3d ( function_3d, result2 )
    call ball_unit_15_3d ( function_3d, result3 )

    write ( *, '(2x,a7,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      function_3d_name(i), result1, result2, result3

  end do
 
  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests BALL_UNIT_F1_ND, BALL_UNIT_F3_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 3

  real ( kind = 8 ) ball_unit_volume_nd
  real ( kind = 8 ), external :: function_nd
  character ( len = 7 ) function_nd_name
  integer ( kind = 4 ) function_nd_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) num
  real ( kind = 8 ) result1
  real ( kind = 8 ) result2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For integrals inside the unit ball in ND:'
  write ( *, '(a)' ) '  BALL_UNIT_F1_ND approximates the integral;'
  write ( *, '(a)' ) '  BALL_UNIT_F3_ND approximates the integral.'
  write ( *, '(a)' ) ' '
 
  do n = 2, n_max
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Spatial dimension N = ', n
    write ( *, '(a,g14.6)' ) '  Unit ball volume = ', ball_unit_volume_nd ( n )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    Rule:      F1          F3'
    write ( *, '(a)' ) '    F(X)'
    write ( *, '(a)' ) ' '
 
    num = function_nd_num ( )

    do i = 1, num

      call function_nd_set ( 'SET', i )

      call ball_unit_f1_nd ( function_nd, n, result1 )
      call ball_unit_f3_nd ( function_nd, n, result2 )
      write ( *, '(2x,a7,2f14.8)' ) function_nd_name(i), result1, result2
 
    end do
 
  end do
 
  return
end
subroutine test045 ( )

!*****************************************************************************80
!
!! TEST045 tests BALL_UNIT_VOLUME_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) ball_unit_volume_3d
  real ( kind = 8 ) ball_unit_volume_nd
  integer ( kind = 4 ), parameter :: dim_num = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST045'
  write ( *, '(a)' ) '  In 3 dimensions:'
  write ( *, '(a)' ) '  BALL_UNIT_VOLUME_3D gets the volume of the unit ball.'
  write ( *, '(a)' ) '  BALL_UNIT_VOLUME_ND will be called for comparison.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    N    Volume    Method'
  write ( *, '(a)' ) ' '

  write ( *, '(2x,i3,f14.8,2x,a)' ) dim_num, &
    ball_unit_volume_3d (   ), 'BALL_UNIT_VOLUME_3D'
  write ( *, '(2x,i3,f14.8,2x,a)' ) dim_num, &
    ball_unit_volume_nd ( dim_num ), 'BALL_UNIT_VOLUME_ND'

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests BALL_UNIT_VOLUME_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none
 
  real ( kind = 8 ) ball_unit_volume_nd
  integer ( kind = 4 ) dim_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  BALL_UNIT_VOLUME_ND computes the volume '
  write ( *, '(a)' ) '    of the unit ball in ND.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    N    Volume'
  write ( *, '(a)' ) ' '

  do dim_num = 2, 10
    write ( *, '(2x,i3,f14.8)' ) dim_num, &
      ball_unit_volume_nd ( dim_num )
  end do

  return
end
subroutine test052 ( )

!*****************************************************************************80
!
!! TEST052 tests BALL_VOLUME_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) ball_volume_3d
  real ( kind = 8 ) ball_volume_nd
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: n = 3
  real ( kind = 8 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST052'
  write ( *, '(a)' ) '  In 3 dimensions:'
  write ( *, '(a)' ) '  BALL_VOLUME_3D computes the volume of a unit ball.'
  write ( *, '(a)' ) '  BALL_VOLUME_ND will be called for comparison.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    N    R      Volume    Method'
  write ( *, '(a)' ) ' '

  r = 1.0D+00

  do i = 1, 3

    write ( *, '(2x,i3,2f14.8,2x,a)' ) &
      n, r, ball_volume_3d (    r ), 'BALL_VOLUME_3D'

    write ( *, '(2x,i3,2f14.8,2x,a)' ) &
      n, r, ball_volume_nd ( n, r ), 'BALL_VOLUME_ND'

    r = r * 2.0D+00

  end do

  return
end
subroutine test054 ( )

!*****************************************************************************80
!
!! TEST054 tests BALL_VOLUME_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) ball_volume_nd
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST054'
  write ( *, '(a)' ) '  BALL_UNIT_VOLUME_ND computes the volume of '
  write ( *, '(a)' ) '    the unit ball in N dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    N        R      Volume'
  write ( *, '(a)' ) ' '

  do n = 2, 10
    r = 0.5D+00
    do i = 1, 3
      write ( *, '(2x,i3,2f14.8)' ) n, r, ball_volume_nd ( n, r )
      r = r * 2.0D+00
    end do
  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests CIRCLE_ANNULUS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) area
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: center_test = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) circle_annulus_area_2d
  real ( kind = 8 ), external :: function_2d
  character ( len = 7 ) function_2d_name
  integer ( kind = 4 ) function_2d_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) num
  integer ( kind = 4 ) nr
  real ( kind = 8 ) radius1
  real ( kind = 8 ), dimension ( test_num ) :: radius1_test = (/ &
    0.0D+00, 1.0D+00 /)
  real ( kind = 8 ) radius2
  real ( kind = 8 ), dimension ( test_num ) :: radius2_test = (/ &
    1.0D+00, 2.0D+00 /)
  real ( kind = 8 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  CIRCLE_ANNULUS estimates integrals '
  write ( *, '(a)' ) '    in a circular annulus.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '        F       CENTER         Radius1   Radius2   NR  Result'
  write ( *, '(a)' ) ' '

  do i = 1, test_num

    center(1:dim_num) = center_test(1:dim_num,i)

    radius1 = radius1_test(i)
    radius2 = radius2_test(i)

    area = circle_annulus_area_2d ( radius1, radius2 )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a7,4f10.6,2x,f10.6)' ) &
      '   Area', center(1:dim_num), radius1, radius2, area

    num = function_2d_num ( )

    do j = 1, num

      call function_2d_set ( 'SET', j )

      do nr = 1, 4
        call circle_annulus ( function_2d, center, radius1, radius2, nr, result )
        write ( *, '(2x,a7,4f10.6,i2,f10.6)' ) &
          function_2d_name(j), center(1:dim_num), radius1, radius2, nr, result
      end do

    end do

  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests CIRCLE_ANNULUS, CIRCLE_RT_SET, CIRCLE_RT_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) area
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: center_test = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) circle_annulus_area_2d
  real ( kind = 8 ), external :: function_2d
  character ( len = 7 ) function_2d_name
  integer ( kind = 4 ) function_2d_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) num
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nr2
  integer ( kind = 4 ) nt
  real ( kind = 8 ) ra(5)
  real ( kind = 8 ) radius1
  real ( kind = 8 ), dimension ( test_num ) :: radius1_test = (/ &
    0.0D+00, 1.0D+00, 1.0D+00 /)
  real ( kind = 8 ) radius2
  real ( kind = 8 ), dimension ( test_num ) :: radius2_test = (/ &
    1.0D+00, 2.0D+00, 3.0D+00 /)
  real ( kind = 8 ) result1
  real ( kind = 8 ) result2
  real ( kind = 8 ) result3
  real ( kind = 8 ) rw(5)
  integer ( kind = 4 ) rule
  real ( kind = 8 ) ta(20)
  real ( kind = 8 ) tw(20)
  real ( kind = 8 ) zw

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  CIRCLE_ANNULUS estimates integrals in a '
  write ( *, '(a)' ) '    circular annulus.'
  write ( *, '(a)' ) '  CIRCLE_RT_SET sets up a rule for a circle;'
  write ( *, '(a)' ) '  CIRCLE_RT_SUM applies the rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RESULT1 = CIRCLE_ANNULUS result.'
  write ( *, '(a)' ) '  RESULT2 = Difference of two CIRCLE_RT_SUM results.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        F      CENTER       Radius1   Radius2   Result1 Result2'
  write ( *, '(a)' ) ' '

  do i = 1, test_num

    center(1:dim_num) = center_test(1:dim_num,i)

    radius1 = radius1_test(i)
    radius2 = radius2_test(i)

    area = circle_annulus_area_2d ( radius1, radius2 )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a7,5f11.5)' ) &
      '   Area', center(1:dim_num), radius1, radius2, area

    rule = 9
    call circle_rt_size ( rule, nr2, nt, nc )
    call circle_rt_set ( rule, nr2, nt, nc, ra, rw, ta, tw, zw )

    num = function_2d_num ( )

    do j = 1, num

      call function_2d_set ( 'SET', j )

      nr = 5
      call circle_annulus ( function_2d, center, radius1, radius2, nr, result1 )

      call circle_rt_sum ( function_2d, center, radius1, nr2, ra, rw, nt, ta, tw, &
        zw, result2 )

      call circle_rt_sum ( function_2d, center, radius2, nr2, ra, rw, nt, ta, tw, &
        zw, result3 )

      write ( *, '(2x,a7,6f11.5)' ) function_2d_name(j), center(1:dim_num), &
        radius1, radius2, result1, result3 - result2

    end do

  end do

  return
end
subroutine test085 ( )

!*****************************************************************************80
!
!! TEST085 tests CIRCLE_ANNULUS_AREA_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) area
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: center_test = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    3.0D+00, 4.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) circle_annulus_area_2d
  integer ( kind = 4 ) i
  real ( kind = 8 ) radius1
  real ( kind = 8 ), dimension ( test_num ) :: radius1_test = (/ &
    0.0D+00, 1.0D+00, 1.0D+00 /)
  real ( kind = 8 ) radius2
  real ( kind = 8 ), dimension ( test_num ) :: radius2_test = (/ &
    1.0D+00, 2.0D+00, 3.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST085'
  write ( *, '(a)' ) '  CIRCLE_ANNULUS_AREA_2D computes the area of a '
  write ( *, '(a)' ) '  circular annulus.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      CENTER       Radius1   Radius2   Area'
  write ( *, '(a)' ) ' '

  do i = 1, test_num

    center(1:dim_num) = center_test(1:dim_num,i)

    radius1 = radius1_test(i)
    radius2 = radius2_test(i)

    area = circle_annulus_area_2d ( radius1, radius2 )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,5f11.5)' ) center(1:dim_num), radius1, radius2, area

  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests CIRCLE_ANNULUS_SECTOR, CIRCLE_RT_SET, CIRCLE_RT_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) as1
  real ( kind = 8 ) as2
  real ( kind = 8 ) as3
  real ( kind = 8 ) as4
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ), external :: function_2d
  character ( len = 7 ) function_2d_name
  integer ( kind = 4 ) function_2d_num
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) num
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nr2
  integer ( kind = 4 ) nt
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) ra(5)
  real ( kind = 8 ) radius
  real ( kind = 8 ) radius1a
  real ( kind = 8 ) radius2a
  real ( kind = 8 ) radius1b
  real ( kind = 8 ) radius2b
  real ( kind = 8 ) radius1c
  real ( kind = 8 ) radius2c
  real ( kind = 8 ) radius1d
  real ( kind = 8 ) radius2d
  real ( kind = 8 ) result1
  real ( kind = 8 ) result2
  integer ( kind = 4 ) rule
  real ( kind = 8 ) rw(5)
  real ( kind = 8 ) ta(20)
  real ( kind = 8 ) theta1a
  real ( kind = 8 ) theta2a
  real ( kind = 8 ) theta1b
  real ( kind = 8 ) theta2b
  real ( kind = 8 ) theta1c
  real ( kind = 8 ) theta2c
  real ( kind = 8 ) theta1d
  real ( kind = 8 ) theta2d
  real ( kind = 8 ) tw(20)
  real ( kind = 8 ) zw

  nr = 5

  rule = 9
  call circle_rt_size ( rule, nr2, nt, nc )
  call circle_rt_set ( rule, nr2, nt, nc, ra, rw, ta, tw, zw )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  CIRCLE_ANNULUS_SECTOR estimates an integral in a '
  write ( *, '(a)' ) '  circular annulus sector.'
  write ( *, '(a)' ) '  CIRCLE_RT_SET sets an integration rule in a circle.'
  write ( *, '(a)' ) '  CIRCLE_RT_SUM uses an integration rule in a circle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  To test CIRCLE_ANNULUS_SECTOR, we estimate an integral'
  write ( *, '(a)' ) '  over 4 annular sectors that make up the unit circle, '
  write ( *, '(a)' ) '  and add to get RESULT1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We will also estimate the integral over the unit circle'
  write ( *, '(a)' ) '  using CIRCLE_RT_SET and CIRCLE_RT_SUM to get RESULT2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We will then compare RESULT1 and RESULT2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  CIRCLE_ANNULUS_SECTOR computations will use NR = ',nr
  write ( *, '(a,i8)' ) '  CIRCLE_RT_SET/CIRCLE_RT_SUM will use rule ', rule
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "RESULT1" is the sum of Annulus Sector calculations.'
  write ( *, '(a)' ) '  "RESULT2" is for CIRCLE_RT_SET/CIRCLE_RT_SUM.'
  write ( *, '(a)' ) ' '

  center(1:dim_num) = (/ 0.0D+00, 0.0D+00 /)
  radius = 1.0D+00

  radius1a = 0.0D+00
  radius2a = 0.25D+00
  theta1a = 0.0D+00
  theta2a = 0.5D+00 * pi

  radius1b = 0.0D+00
  radius2b = 0.25D+00
  theta1b = 0.5D+00 * pi
  theta2b = 2.0D+00 * pi

  radius1c = 0.25D+00
  radius2c = 1.0D+00
  theta1c = 0.0D+00
  theta2c = 0.25D+00 * pi

  radius1d = 0.25D+00
  radius2d = 1.0D+00
  theta1d = 0.25D+00 * pi
  theta2d = 2.0D+00 * pi

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       F  Result1  Result2'
  write ( *, '(a)' ) ' '

  num = function_2d_num ( )

  do j = 1, num

    call function_2d_set ( 'SET', j )

    call circle_annulus_sector ( function_2d, center, radius1a, radius2a, theta1a, &
      theta2a, nr, as1 )

    call circle_annulus_sector ( function_2d, center, radius1b, radius2b, theta1b, &
      theta2b, nr, as2 )

    call circle_annulus_sector ( function_2d, center, radius1c, radius2c, theta1c, &
      theta2c, nr, as3 )

    call circle_annulus_sector ( function_2d, center, radius1d, radius2d, theta1d, &
      theta2d, nr, as4 )

    result1 = as1 + as2 + as3 + as4

    call circle_rt_sum ( function_2d, center, radius, nr2, ra, rw, nt, ta, tw, zw, &
      result2 )

    write ( *, '(2x,a7,2g14.6)' ) function_2d_name(j), result1, result2

  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests CIRCLE_CUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ), dimension ( dim_num ) :: center = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), external :: function_2d
  character ( len = 7 ) function_2d_name
  integer ( kind = 4 ) function_2d_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) num
  integer ( kind = 4 ) order
  real ( kind = 8 ) :: r = 3.0D+00
  real ( kind = 8 ) result(4)
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  CIRCLE_CUM approximates an integral over a circle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  We use radius R = ', r
  write ( *, '(a)' ) '  and center:'
  write ( *, '(a,2g14.6)' ) '  CENTER = ', center(1:dim_num)
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Order:      2             4              8     ' // &
    '       16'
  write ( *, '(a)' ) '  F(X)'
  write ( *, '(a)' ) ' '
 
  num = function_2d_num ( )

  do i = 1, num

    call function_2d_set ( 'SET', i )

    do j = 1, 4

      order = 2**j

      call circle_cum ( function_2d, center, r, order, result(j) )

    end do

    write ( *, '(2x,a7,5f14.8)' ) function_2d_name(i), result(1:4)
 
  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests LENS_HALF_AREA_2D, CIRCLE_SECTOR_AREA_2D, CIRCLE_TRIANGLE_AREA_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) area1
  real ( kind = 8 ) area2
  real ( kind = 8 ) area3
  real ( kind = 8 ) lens_half_area_2d
  real ( kind = 8 ) circle_sector_area_2d
  real ( kind = 8 ) circle_triangle_area_2d
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2

  r = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  LENS_HALF_AREA_2D computes the area of a'
  write ( *, '(a)' ) '  circular half lens, defined by joining the endpoints'
  write ( *, '(a)' ) '  of a circular arc.'
  write ( *, '(a)' ) '  CIRCLE_SECTOR_AREA_2D computes the area of a'
  write ( *, '(a)' ) '  circular sector, defined by joining the endpoints'
  write ( *, '(a)' ) '  of a circular arc to the center.'
  write ( *, '(a)' ) '  CIRCLE_TRIANGLE_AREA_2D computes the signed area of a'
  write ( *, '(a)' ) '  triangle, defined by joining the endpoints'
  write ( *, '(a)' ) '  of a circular arc and the center.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      R       Theta1 Theta2        ' // &
    'Sector       Triangle     Half Lens'
  write ( *, '(a)' ) ' '

  do i = 0, 12

    theta1 = 0.0D+00
    theta2 = real ( i, kind = 8 ) * 2.0D+00 * pi / 12.0D+00

    area1 = circle_sector_area_2d ( r, theta1, theta2 )

    area2 = circle_triangle_area_2d ( r, theta1, theta2 )

    area3 = lens_half_area_2d ( r, theta1, theta2 )

    write ( *, '(2x,f9.3,f9.3,4f14.8)' ) r, theta1, theta2, area1, area2, area3

  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests LENS_HALF_AREA_2D, LENS_HALF_H_AREA_2D, LENS_HALF_W_AREA_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) area1
  real ( kind = 8 ) area2
  real ( kind = 8 ) area3
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  real ( kind = 8 ) lens_half_area_2d
  real ( kind = 8 ) lens_half_h_area_2d
  real ( kind = 8 ) lens_half_w_area_2d
  real ( kind = 8 ) r
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  real ( kind = 8 ) w

  r = 50.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  For the area of a circular half lens,'
  write ( *, '(a)' ) '  LENS_HALF_AREA_2D uses two angles;'
  write ( *, '(a)' ) '  LENS_HALF_H_AREA_2D works from the height;'
  write ( *, '(a)' ) '  LENS_HALF_W_AREA_2D works from the width.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The circle has radius R = ', r
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  THETA1 THETA2  H     W  Area(THETA) Area(H)  Area(W)'
  write ( *, '(a)' ) ' '

  do i = 0, 12

    theta1 = 0.0D+00
    theta2 = real ( i, kind = 8 ) * 2.0D+00 * pi / 12.0D+00
    w = 2.0D+00 * r * sin ( 0.5D+00 * ( theta2 - theta1 ) ) 
    h = r * ( 1.0D+00 - cos ( 0.5D+00 * ( theta2 - theta1 ) ) )

    area1 = lens_half_area_2d ( r, theta1, theta2 )

    area2 = lens_half_h_area_2d ( r, h )

    area3 = lens_half_w_area_2d ( r, w )

    write ( *, '(2x,4f6.2,3f10.4)' ) theta1, theta2, h, w, area1, area2, area3

  end do

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests CIRCLE_SECTOR, CIRCLE_SECTOR_AREA_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: nrlo = 1
  integer ( kind = 4 ), parameter :: nrhi = 5
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) area
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: center_test = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) circle_sector_area_2d
  real ( kind = 8 ), external :: function_2d
  character ( len = 7 ) function_2d_name
  integer ( kind = 4 ) function_2d_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) num
  integer ( kind = 4 ) nr
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) radius
  real ( kind = 8 ), dimension ( test_num ) :: radius_test = (/ &
    1.0D+00, 2.0D+00, 4.0D+00, 8.0D+00 /)
  real ( kind = 8 ) result(nrlo:nrhi)
  real ( kind = 8 ) theta1
  real ( kind = 8 ), dimension ( test_num ) :: theta1_test = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) theta2
  real ( kind = 8 ), dimension ( test_num ) :: theta2_test = (/ &
    2.0D+00, 1.0D+00, 0.5D+00, 0.25D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  CIRCLE_SECTOR_AREA_2D computes the area '
  write ( *, '(a)' ) '  of a circular sector.'
  write ( *, '(a)' ) '  CIRCLE_SECTOR estimates an integral '
  write ( *, '(a)' ) '  in a circular sector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The user can specify NR, the number of radial values'
  write ( *, '(a)' ) '  used to approximated the integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, computations will use values of NR'
  write ( *, '(a,i8)' ) '  from ', nrlo
  write ( *, '(a,i8)' ) '  to   ', nrhi
  write ( *, '(a)' ) ' '

  do i = 1, test_num

    center(1:dim_num) = center_test(1:dim_num,i)

    radius = radius_test(i)
    theta1 = theta1_test(i) * pi
    theta2 = theta2_test(i) * pi

    area = circle_sector_area_2d ( radius, theta1, theta2 )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '       CENTER      RADIUS  THETA1  THETA2  Area'
    write ( *, '(a)' ) ' '
    write ( *, '(6f8.4)' ) center(1:dim_num), radius, theta1, theta2, area
    write ( *, '(a)' ) ' '
    write ( *, '(a7,14(6x,i2,6x))' ) '   F   ', ( nr, nr = nrlo, nrhi )
    write ( *, '(a)' ) ' '

    num = function_2d_num ( )

    do j = 1, num

      call function_2d_set ( 'SET', j )

      do nr = nrlo, nrhi
        call circle_sector ( function_2d, center, radius, theta1, theta2, nr, &
          result(nr) )
      end do

      write ( *, '(2x,a7,5g14.6)' ) function_2d_name(j), result(nrlo:nrhi)

    end do

  end do

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests CIRCLE_SECTOR, CIRCLE_RT_SET, CIRCLE_RT_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) area1
  real ( kind = 8 ) area2
  real ( kind = 8 ) area3
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: center_test = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) circle_area_2d
  real ( kind = 8 ) circle_sector_area_2d
  real ( kind = 8 ), external :: function_2d
  character ( len = 7 ) function_2d_name
  integer ( kind = 4 ) function_2d_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) num
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nr2
  integer ( kind = 4 ) nt
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) ra(5)
  real ( kind = 8 ) radius
  real ( kind = 8 ), dimension ( test_num ) :: radius_test = (/ &
    1.0D+00, 2.0D+00, 4.0D+00, 8.0D+00 /)
  real ( kind = 8 ) result1
  real ( kind = 8 ) result2
  real ( kind = 8 ) resulta
  real ( kind = 8 ) resultb
  integer ( kind = 4 ) rule
  real ( kind = 8 ) rw(5)
  real ( kind = 8 ) ta(20)
  real ( kind = 8 ) theta1
  real ( kind = 8 ), dimension ( test_num ) :: theta1_test = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) theta2
  real ( kind = 8 ), dimension ( test_num ) :: theta2_test = (/ &
    2.0D+00, 1.0D+00, 0.5D+00, 0.25D+00 /)
  real ( kind = 8 ) theta3
  real ( kind = 8 ) tw(20)
  real ( kind = 8 ) zw

  nr = 5

  rule = 9
  call circle_rt_size ( rule, nr2, nt, nc )
  call circle_rt_set ( rule, nr2, nt, nc, ra, rw, ta, tw, zw )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  CIRCLE_SECTOR estimates integrals in a circular sector.'
  write ( *, '(a)' ) '  CIRCLE_RT_SET sets an integration rule in a circle.'
  write ( *, '(a)' ) '  CIRCLE_RT_SUM uses an integration rule in a circle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  To test CIRCLE_SECTOR, we estimate an integral over'
  write ( *, '(a)' ) '  a sector, and over its complement and add the results'
  write ( *, '(a)' ) '  to get RESULT1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We also estimate the integral over the whole circle'
  write ( *, '(a)' ) '  using CIRCLE_RT_SET and CIRCLE_RT_SUM to get RESULT2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We will then compare RESULT1 and RESULT2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  CIRCLE_SECTOR computations will use NR = ',nr
  write ( *, '(a,i8)' ) '  CIRCLE_RT_SET/CIRCLE_RT_SUM will use rule ', rule
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "Sector1" and "Sector2" are the CIRCLE_SECTOR'
  write ( *, '(a)' ) '  computations'
  write ( *, '(a)' ) '  for the sector and its complement.'
  write ( *, '(a)' ) '  "Sum" is the sum of Sector1 and Sector2.'
  write ( *, '(a)' ) '  "Circle" is the computation for '
  write ( *, '(a)' ) '  CIRCLE_RT_SET + CIRCLE_RT_SUM.'
  write ( *, '(a)' ) ' '

  do i = 1, test_num

    center(1:dim_num) = center_test(1:dim_num,i)

    radius = radius_test(i)

    theta1 = theta1_test(i) * pi
    theta2 = theta2_test(i) * pi
    theta3 = theta2 + 2.0D+00 * pi - ( theta2 - theta1 )

    area1 = circle_sector_area_2d ( radius, theta1, theta2 )
    area2 = circle_sector_area_2d ( radius, theta2, theta3 )
    area3 = circle_area_2d ( radius )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '      CENTER       RADIUS   THETA1   THETA2   Area1   Area2  Circle'
    write ( *, '(a)' ) ' '
    write ( *, '(8f9.4)' ) center(1:dim_num), radius, theta1, theta2, area1, area2, area3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '       F   Sector1       Sector2         Sum         Circle'
    write ( *, '(a)' ) ' '

    num = function_2d_num ( )

    do j = 1, num

      call function_2d_set ( 'SET', j )

      call circle_sector ( function_2d, center, radius, theta1, theta2, nr, &
        resulta )

      call circle_sector ( function_2d, center, radius, theta2, theta3, nr, &
        resultb )

      result1 = resulta + resultb

      call circle_rt_sum ( function_2d, center, radius, nr2, ra, rw, nt, ta, &
        tw, zw, result2 )

      write ( *, '(2x,a7,4g14.6)' ) function_2d_name(j), resulta, resultb, &
        result1, result2

    end do

  end do

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests CIRCLE_RT_SET and CIRCLE_RT_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: rule_max = 9

  real ( kind = 8 ), dimension ( dim_num ) :: center = (/ 1.0D+00, 1.0D+00 /)
  real ( kind = 8 ), external :: function_2d
  character ( len = 7 ) function_2d_name
  integer ( kind = 4 ) function_2d_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) num
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nt
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: r = 1.0D+00
  real ( kind = 8 ), allocatable, dimension ( : ) :: ra
  real ( kind = 8 ) result(rule_max)
  real ( kind = 8 ), allocatable, dimension ( : ) :: rw
  integer ( kind = 4 ) rule
  real ( kind = 8 ), allocatable, dimension ( : ) :: ta
  real ( kind = 8 ), allocatable, dimension ( : ) :: tw
  real ( kind = 8 ) zw

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  For R, Theta product rules on the unit circle,'
  write ( *, '(a)' ) '  CIRCLE_RT_SET sets a rule.'
  write ( *, '(a)' ) '  CIRCLE_RT_SUM uses the rule in an arbitrary circle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  We use a radius ', r
  write ( *, '(a)' ) '  and center:'
  write ( *, '(a,2g14.6)' ) '  CENTER = ', center(1:dim_num)
  write ( *, '(a)' ) ' '

  do ilo = 1, rule_max, 5

    ihi = min ( ilo +  4, rule_max )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a7,7x,5(i7,7x))' ) 'Rule:  ', ( rule, rule = ilo, ihi )
    write ( *, '(2x,a)' ) 'Function '

    num = function_2d_num ( )

    do i = 1, num

      call function_2d_set ( 'SET', i )

      do rule = ilo, ihi

        call circle_rt_size ( rule, nr, nt, nc )

        allocate ( ra(1:nr) )
        allocate ( rw(1:nr) )
        allocate ( ta(1:nt) )
        allocate ( tw(1:nt) )

        call circle_rt_set ( rule, nr, nt, nc, ra, rw, ta, tw, zw )

        call circle_rt_sum ( function_2d, center, r, nr, ra, rw, nt, ta, tw, zw, &
          result(rule) )

        deallocate ( ra )
        deallocate ( rw )
        deallocate ( ta )
        deallocate ( tw )

      end do

      write ( *, '(2x,a7,5f14.8)' ) function_2d_name(i), result(ilo:ihi)

    end do

  end do

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests CIRCLE_XY_SET and CIRCLE_XY_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: rule_max = 13

  real ( kind = 8 ), dimension ( dim_num ) :: center = (/ 1.0D+00, 1.0D+00 /)
  real ( kind = 8 ), external :: function_2d
  character ( len = 7 ) function_2d_name
  integer ( kind = 4 ) function_2d_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) num
  integer ( kind = 4 ) order
  real ( kind = 8 ) r
  real ( kind = 8 ) result(rule_max)
  integer ( kind = 4 ) rule
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight
  real ( kind = 8 ), allocatable, dimension ( : ) :: xtab
  real ( kind = 8 ), allocatable, dimension ( : ) :: ytab

  r = 1.0D+00
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  CIRCLE_XY_SET sets a quadrature rule '
  write ( *, '(a)' ) '  for the unit circle.'
  write ( *, '(a)' ) '  CIRCLE_XY_SUM evaluates the quadrature rule'
  write ( *, '(a)' ) '  in an arbitrary circle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  We use a radius ', r
  write ( *, '(a)' ) '  and center:'
  write ( *, '(a,2g14.6)' ) '  CENTER = ', center(1:dim_num)
  write ( *, '(a)' ) ' '

  do ilo = 1, rule_max, 5

    ihi = min ( ilo +  4, rule_max )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a7,7x,5(i7,7x))' ) 'Rule:  ', ( rule, rule = ilo, ihi )
    write ( *, '(2x,a)' ) 'Function '

    num = function_2d_num ( )

    do i = 1, num

      call function_2d_set ( 'SET', i )

      do rule = ilo, ihi

        call circle_xy_size ( rule, order )

        allocate ( xtab(1:order) )
        allocate ( ytab(1:order) )
        allocate ( weight(1:order) )

        call circle_xy_set ( rule, order, xtab, ytab, weight )

        call circle_xy_sum ( function_2d, center, r, order, xtab, ytab, weight, &
          result(rule) )

        deallocate ( weight )
        deallocate ( xtab )
        deallocate ( ytab )

      end do

      write ( *, '(2x,a7,5f14.8)' ) function_2d_name(i), result(ilo:ihi)

    end do

  end do

  return
end
subroutine test163 ( )

!*****************************************************************************80
!
!! TEST163 tests the rules for CN with Gegenbauer weight on monomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  real ( kind = 8 ) alpha
  real ( kind = 8 ), dimension ( test_num ) :: alpha_test = (/ &
    - 0.5D+00, 0.0D+00, 0.5D+00, 1.0D+00, 1.5D+00 /)
  integer ( kind = 4 ), allocatable :: expon(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST163'
  write ( *, '(a)' ) '  Demonstrate the use of quadrature rules for the region'
  write ( *, '(a)' ) '  CN_GEG, that is, the hypercube [-1,+1]^N, with the'
  write ( *, '(a)' ) '  weight W(ALPHA;X) = product ( 1 <= I <= N )'
  write ( *, '(a)' ) '  7771-X(I)^2)^ALPHA'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use the formulas to integrate various monomials of'
  write ( *, '(a)' ) '  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)'
  write ( *, '(a)' ) '  and compare to the exact integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The precision of each formula is known, and we only use'
  write ( *, '(a)' ) '  a formula if its precision indicates it should be able to'
  write ( *, '(a)' ) '  produce an exact result.'

  do n = 1, 6

    allocate ( expon(1:n) )

    do test = 1, test_num

      alpha = alpha_test(test)

      expon(1:n) = 0
      call cn_geg_test ( n, alpha, expon )

    end do

    do test = 1, test_num

      alpha = alpha_test(test)

      expon(1:n) = 0
      expon(n) = 1
      call cn_geg_test ( n, alpha, expon )

    end do

    if ( 2 <= n ) then

      do test = 1, test_num

        alpha = alpha_test(test)

        expon(1:n) = 0
        expon(1) = 1
        expon(2) = 1
        call cn_geg_test ( n, alpha, expon )

      end do

    end if

    do test = 1, test_num

      alpha = alpha_test(test)

      expon(1:n) = 0
      expon(1) = 2
      call cn_geg_test ( n, alpha, expon )

    end do

    deallocate ( expon )

  end do

  return
end
subroutine cn_geg_test ( n, alpha, expon )

!*****************************************************************************80
!
!! CN_GEG_TEST tests the rules for CN with Gegenbauer weight on a monomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) c1
  integer ( kind = 4 ) d
  real ( kind = 8 ) delta0
  real ( kind = 8 ) err
  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon(n)
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option
  integer ( kind = 4 ) p
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) quad
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:,:)

  write ( *, '(a)'      ) ' '
  write ( *, '(a,i4)'   ) '  N = ', n
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,10i4)' ) '  EXPON = ', expon(1:n)
  d = sum ( expon(1:n) )
  write ( *, '(a,i4)'   ) '  Degree = ', d
  write ( *, '(a)'      ) ' '

  call cn_geg_monomial_integral ( n, alpha, expon, exact )

  p = 0

  if ( d <= p ) then
    call cn_geg_00_1_size ( n, alpha, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_geg_00_1 ( n, alpha, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_GEG_00_1:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 1

  if ( d <= p ) then
    call cn_geg_01_1_size ( n, alpha, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_geg_01_1 ( n, alpha, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_GEG_01_1:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 2

  if ( d <= p ) then

    call cn_geg_02_xiu_size ( n, alpha, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_geg_02_xiu ( n, alpha, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_GEG_02_XIU: ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call gw_02_xiu_size ( n, o )
    gamma0 = 1.0D+00
    delta0 = 0.0D+00
    c1 = 1.0D+00 / ( 2.0D+00 * alpha + 3.0D+00 )
    volume_1d = sqrt ( pi ) * r8_gamma ( alpha + 1.0D+00 ) &
      / r8_gamma ( alpha + 1.5D+00 )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  GW_02_XIU:     ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

  end if

  p = 3

  if ( d <= p ) then
    call cn_geg_03_xiu_size ( n, alpha, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_geg_03_xiu ( n, alpha, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_GEG_03_XIU: ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  write ( *, '(a,2x,6x,2x,g14.6)' ) '  EXACT:         ', exact

  return
end
subroutine test165 ( )

!*****************************************************************************80
!
!! TEST165 tests the rules for CN with Jacobi weight on monomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) alpha
  real ( kind = 8 ), dimension ( test_num ) :: alpha_test = (/ &
    0.0D+00, 1.0D+00, 0.0D+00, 0.5D+00 /)
  real ( kind = 8 ) beta
  real ( kind = 8 ), dimension ( test_num ) :: beta_test = (/ &
    0.0D+00, 0.0D+00, 2.0D+00, 1.5D+00 /)
  integer ( kind = 4 ), allocatable :: expon(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST165'
  write ( *, '(a)' ) '  Demonstrate the use of quadrature rules for the region'
  write ( *, '(a)' ) '  CN_JAC, that is, the hypercube [-1,+1]^N, with the'
  write ( *, '(a)' ) '  weight W(ALPHA,BETA;X) = product ( 1 <= I <= N )'
  write ( *, '(a)' ) '    (1-X(I))^ALPHA (1+X(I))^BETA'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use the formulas to integrate various monomials of'
  write ( *, '(a)' ) '  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)'
  write ( *, '(a)' ) '  and compare to the exact integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The precision of each formula is known, and we only use'
  write ( *, '(a)' ) '  a formula if its precision indicates it should be able to'
  write ( *, '(a)' ) '  produce an exact result.'

  do n = 1, 6

    allocate ( expon(1:n) )

    do test = 1, test_num

      alpha = alpha_test(test)
      beta  = beta_test(test)

      expon(1:n) = 0
      call cn_jac_test ( n, alpha, beta, expon )

    end do

    do test = 1, test_num

      alpha = alpha_test(test)
      beta  = beta_test(test)

      expon(1:n) = 0
      expon(n) = 1
      call cn_jac_test ( n, alpha, beta, expon )

    end do

    if ( 2 <= n ) then

      do test = 1, test_num

        alpha = alpha_test(test)
        beta  = beta_test(test)

        expon(1:n) = 0
        expon(1) = 1
        expon(2) = 1
        call cn_jac_test ( n, alpha, beta, expon )

      end do

    end if

    do test = 1, test_num

      alpha = alpha_test(test)
      beta  = beta_test(test)

      expon(1:n) = 0
      expon(1) = 2
      call cn_jac_test ( n, alpha, beta, expon )

    end do

    deallocate ( expon )

  end do

  return
end
subroutine cn_jac_test ( n, alpha, beta, expon )

!*****************************************************************************80
!
!! CN_JAC_TEST tests the rules for CN with Jacobi weight on a monomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) c1
  integer ( kind = 4 ) d
  real ( kind = 8 ) delta0
  real ( kind = 8 ) err
  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon(n)
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option
  integer ( kind = 4 ) p
  real ( kind = 8 ) quad
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:,:)

  write ( *, '(a)'      ) ' '
  write ( *, '(a,i4)'   ) '  N = ', n
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,g14.6)' ) '  BETA =  ', beta
  write ( *, '(a,10i4)' ) '  EXPON = ', expon(1:n)
  d = sum ( expon(1:n) )
  write ( *, '(a,i4)'   ) '  Degree = ', d
  write ( *, '(a)'      ) ' '

  call cn_jac_monomial_integral ( n, alpha, beta, expon, exact )

  p = 0

  if ( d <= p ) then
    call cn_jac_00_1_size ( n, alpha, beta, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_jac_00_1 ( n, alpha, beta, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_JAC_00_1:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 1

  if ( d <= p ) then
    call cn_jac_01_1_size ( n, alpha, beta, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_jac_01_1 ( n, alpha, beta, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_JAC_01_1:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 2

  if ( d <= p ) then

    call cn_jac_02_xiu_size ( n, alpha, beta, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_jac_02_xiu ( n, alpha, beta, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_JAC_02_XIU: ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call gw_02_xiu_size ( n, o )
    gamma0 = ( alpha + beta + 2.0D+00 ) / 2.0D+00
    delta0 = ( alpha - beta ) / 2.0D+00
    c1 = 2.0D+00 * ( alpha + 1.0D+00 ) * ( beta + 1.0D+00 ) &
      / ( alpha + beta + 3.0D+00 ) / ( alpha + beta + 2.0D+00 )
    volume_1d = 2.0D+00 ** ( alpha + beta + 1.0D+00 ) &
      * r8_gamma ( alpha + 1.0D+00 ) * r8_gamma ( beta + 1.0D+00 ) &
      / ( alpha + beta + 1.0D+00 ) / r8_gamma ( alpha + beta + 1.0D+00 )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  GW_02_XIU:     ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

  end if

  write ( *, '(a,2x,6x,2x,g14.6)' ) '  EXACT:         ', exact

  return
end
subroutine test167 ( )

!*****************************************************************************80
!
!! TEST167 tests the rules for CN with Legendre weight on monomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable :: expon(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST167'
  write ( *, '(a)' ) '  Demonstrate the use of quadrature rules for the region'
  write ( *, '(a)' ) '  CN_LEG, that is, the hypercube [-1,+1]^N, with the'
  write ( *, '(a)' ) '  Legendre weight W(X) = 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use the formulas to integrate various monomials of'
  write ( *, '(a)' ) '  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)'
  write ( *, '(a)' ) '  and compare to the exact integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The precision of each formula is known, and we only use'
  write ( *, '(a)' ) '  a formula if its precision indicates it should be able to'
  write ( *, '(a)' ) '  produce an exact result.'

  do n = 1, 6

    allocate ( expon(1:n) )

    expon(1:n) = 0
    call cn_leg_test ( n, expon )

    expon(1:n) = 0
    expon(n) = 1
    call cn_leg_test ( n, expon )

    if ( 2 <= n ) then

      expon(1:n) = 0
      expon(1) = 1
      expon(2) = 1
      call cn_leg_test ( n, expon )

    end if

    expon(1:n) = 0
    expon(1) = 2
    call cn_leg_test ( n, expon )

    expon(1:n) = 0
    expon(1) = 3
    call cn_leg_test ( n, expon )

    expon(1:n) = 0
    expon(n) = 4
    call cn_leg_test ( n, expon )

    if ( 2 <= n ) then
      expon(1:n) = 0
      expon(1) = 3
      expon(2) = 2
      call cn_leg_test ( n, expon )
    end if

    deallocate ( expon )

  end do

  return
end
subroutine cn_leg_test ( n, expon )

!*****************************************************************************80
!
!! CN_LEG_TEST tests the rules for CN with Legendre weight on a monomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c1
  integer ( kind = 4 ) d
  real ( kind = 8 ) delta0
  real ( kind = 8 ) err
  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon(n)
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option
  integer ( kind = 4 ) p
  real ( kind = 8 ) quad
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:,:)

  write ( *, '(a)'      ) ' '
  write ( *, '(a,i4)'   ) '  N = ', n
  write ( *, '(a,10i4)' ) '  EXPON = ', expon(1:n)
  d = sum ( expon(1:n) )
  write ( *, '(a,i4)'   ) '  Degree = ', d
  write ( *, '(a)'      ) ' '

  call cn_leg_monomial_integral ( n, expon, exact )

  p = 1

  if ( d <= p ) then
    call cn_leg_01_1_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_leg_01_1 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_LEG_01_1:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 2

  if ( d <= p ) then

    call cn_leg_02_xiu_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_leg_02_xiu ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_LEG_02_XIU: ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call gw_02_xiu_size ( n, o )
    gamma0 = 1.0D+00
    delta0 = 0.0D+00
    c1 = 1.0D+00 / 3.0D+00
    volume_1d = 2.0D+00
    allocate ( x(n,o) )
    allocate ( w(o) )
    call gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  GW_02_XIU:     ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

  end if

  p = 3

  if ( d <= p ) then

    call cn_leg_03_1_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_leg_03_1 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_LEG_03_1:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call cn_leg_03_xiu_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call cn_leg_03_xiu ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_LEG_03_XIU: ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

  end if

  p = 5

  if ( d <= p ) then

    if ( 4 <= n .and. n <= 6 ) then
      call cn_leg_05_1_size ( n, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      option = 1
      call cn_leg_05_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_LEG_05_1(1):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

    if ( 4 <= n .and. n <= 5 ) then
      call cn_leg_05_1_size ( n, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      option = 2
      call cn_leg_05_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_LEG_05_1(2):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

    if ( 2 <= n ) then
      call cn_leg_05_2_size ( n, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call cn_leg_05_2 ( n, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  CN_LEG_05_2:   ', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

  end if

  write ( *, '(a,2x,6x,2x,g14.6)' ) '  EXACT:         ', exact

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests CONE_UNIT_3D, CONE_VOLUME_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) cone_volume_3d
  character ( len = 7 ) function_3d_name
  integer ( kind = 4 ) function_3d_num
  real ( kind = 8 ), external :: function_3d
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) num
  real ( kind = 8 ), parameter :: r = 1.0D+00
  real ( kind = 8 ) result

  h = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  CONE_UNIT_3D approximates integrals in a unit cone.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Volume = ', cone_volume_3d ( r, h )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    F(X)    CONE_3D'
  write ( *, '(a)' ) ' '

  num = function_3d_num ( )

  do i = 1, num

    call function_3d_set ( 'SET', i )

    call cone_unit_3d ( function_3d, result )

    write ( *, '(2x,a7,f14.8)' ) function_3d_name(i), result

  end do

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests CUBE_SHELL_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 4
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) cube_shell_volume_nd
  real ( kind = 8 ), external :: function_nd
  character ( len = 7 ) function_nd_name
  integer ( kind = 4 ) function_nd_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) num
  real ( kind = 8 ) r1
  real ( kind = 8 ), dimension ( test_num ) :: r1_test = (/ &
    0.0D+00, 1.0D+00 /)
  real ( kind = 8 ) :: r2
  real ( kind = 8 ), dimension ( test_num ) :: r2_test = (/ &
    1.0D+00, 2.0D+00 /)
  real ( kind = 8 ) result
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  CUBE_SHELL_ND approximates integrals in a'
  write ( *, '(a)' ) '  cubical shell in ND.'

  do test = 1, test_num

    r1 = r1_test(test)
    r2 = r2_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Inner radius = ', r1
    write ( *, '(a,g14.6)' ) '  Outer radius = ', r2
    write ( *, '(a)' ) ' '
 
    do n = 2, n_max
 
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Spatial dimension N = ', n
      write ( *, '(a,g14.6)' ) '  Volume = ', cube_shell_volume_nd ( n, r1, r2 )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    F(X)      CUBE_SHELL_ND'
      write ( *, '(a)' ) ' '
 
      num = function_nd_num ( )

      do i = 1, num

        call function_nd_set ( 'SET', i )

        call cube_shell_nd ( function_nd, n, r1, r2, result )

        write ( *, '(2x,a7,f14.8)' ) function_nd_name(i), result

      end do

    end do
 
  end do
 
  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests CUBE_UNIT_3D, QMULT_3D, RECTANGLE_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a1
  real ( kind = 8 ) a(3)
  real ( kind = 8 ) b1
  real ( kind = 8 ) b(3)
  real ( kind = 8 ), external :: fl18
  real ( kind = 8 ), external :: fl28
  character ( len = 7 ) function_3d_name
  integer ( kind = 4 ) function_3d_num
  real ( kind = 8 ), external :: fu18
  real ( kind = 8 ), external :: fu28
  real ( kind = 8 ), external :: function_3d
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ) num
  real ( kind = 8 ) qmult_3d
  real ( kind = 8 ) result1
  real ( kind = 8 ) result2
  real ( kind = 8 ) result3

  a1 = -1.0D+00
  b1 = +1.0D+00

  a(1) = -1.0D+00
  a(2) = -1.0D+00
  a(3) = -1.0D+00
  b(1) = 1.0D+00
  b(2) = 1.0D+00
  b(3) = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  CUBE_UNIT_3D approximates integrals '
  write ( *, '(a)' ) '  in the unit cube in 3D.'
  write ( *, '(a)' ) '  QMULT_3D approximates triple integrals.'
  write ( *, '(a)' ) '  RECTANGLE_3D approximates integrals '
  write ( *, '(a)' ) '  in a rectangular block.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    F(X)    CUBE_UNIT_3D  QMULT_3D      RECTANGLE_3D'
  write ( *, '(a)' ) ' '

  num = function_3d_num ( )

  do i = 1, num

    call function_3d_set ( 'SET', i )

    call cube_unit_3d ( function_3d, result1 )
    result2 = qmult_3d ( function_3d, a1, b1, fu18, fl18, fu28, fl28 )
    call rectangle_3d ( function_3d, a, b, result3 )

    write ( *, '(2x,a7,3f14.8)' ) function_3d_name(i), result1, result2, result3

  end do

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests CUBE_UNIT_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_k = 10
  integer ( kind = 4 ), parameter :: max_test = 2

  real ( kind = 8 ), external :: function_nd
  character ( len = 7 ) function_nd_name
  integer ( kind = 4 ) function_nd_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) k
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  integer ( kind = 4 ), parameter, dimension ( max_test ) :: k_test = (/ 10, 5 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter, dimension ( max_test ) :: n_test = (/ 2, 3 /)
  integer ( kind = 4 ) num
  real ( kind = 8 ) qa(max_k)
  real ( kind = 8 ) qb(max_k)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  CUBE_UNIT_ND approximates integrals inside '
  write ( *, '(a)' ) '  the unit cube in ND.'

  do i_test = 1, max_test
 
    n = n_test(i_test)
    k = k_test(i_test)
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Spatial dimension N = ', n
    write ( *, '(a,i8)' ) '  Value of K = ', k
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    F(X)    CUBE_UNIT_ND'
    write ( *, '(a)' ) ' '

    num = function_nd_num ( )

    do i = 1, num

      call function_nd_set ( 'SET', i )

      call cube_unit_nd ( function_nd, qa, qb, n, k )

      do klo = 1, k, 5
        khi = min ( klo + 4, k )
        if ( klo == 1 ) then
          write ( *, '(2x,a7,5f14.8)' ) function_nd_name(i), qa(klo:khi)
        else
          write ( *, '(2x,7x,5f14.8)' )           qa(klo:khi)
        end if
      end do

      do klo = 1, k, 5
        khi = min ( klo + 4, k )
        write ( *, '(2x,7x,5f14.8)' )           qb(klo:khi)
      end do
                                           
    end do

  end do

  return
end
subroutine test205 ( )

!*****************************************************************************80
!
!! TEST205 tests ELLIPSE_AREA_2D, ELLIPSE_CIRCUMFERENCE_2D, ELLIPSE_ECCENTRICITY_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) e
  real ( kind = 8 ) ellipse_area_2d
  real ( kind = 8 ) ellipse_circumference_2d
  real ( kind = 8 ) ellipse_eccentricity_2d
  integer ( kind = 4 ) i
  real ( kind = 8 ) p
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST205'
  write ( *, '(a)' ) '  ELLIPSE_AREA_2D returns the area of an ellipse.'
  write ( *, '(a)' ) '  ELLIPSE_ECCENTRICITY_2D returns the '
  write ( *, '(a)' ) '  eccentricity of an ellipse.'
  write ( *, '(a)' ) '  ELLIPSE_CIRCUMFERENCE_2D returns the '
  write ( *, '(a)' ) '  circumference of an ellipse.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        R1        R2        E         Circum    Area'
  write ( *, '(a)' ) ' '

  do i = 1, 5

    if ( i == 1 ) then
      r1 = 25.0D+00
      r2 = 20.0D+00
    else
      r1 = r8_uniform_01 ( seed )
      r2 = r8_uniform_01 ( seed )
    end if

    e = ellipse_eccentricity_2d ( r1, r2 )
    p = ellipse_circumference_2d ( r1, r2 )    
    area = ellipse_area_2d ( r1, r2 )

    write ( *, '(2x,5f10.4)' ) r1, r2, e, p, area

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (For the first example, '
  write ( *, '(a)' ) '  the eccentricity should be 0.6,'
  write ( *, '(a)' ) '  the circumference should be about 141.8).'

  return
end
subroutine test207 ( )

!*****************************************************************************80
!
!! TEST207 tests the Stroud EN_R2 rules on monomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable :: alpha(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST207'
  write ( *, '(a)' ) '  Demonstrate the use of Stroud rules for the region'
  write ( *, '(a)' ) '  EN_R2, that is, all of N-dimensional space, with the'
  write ( *, '(a)' ) '  weight function W(X) = exp ( - X1^2 - X2^2 ... -XN^2 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use the formulas to integrate various monomials of'
  write ( *, '(a)' ) '  the form X1^ALPHA1 * X2^ALPHA2 * ... XN^ALPHAN'
  write ( *, '(a)' ) '  and compare to the exact integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The precision of each formula is known, and we only use'
  write ( *, '(a)' ) '  a formula if its precision indicates it should be able to'
  write ( *, '(a)' ) '  produce an exact result.'

  do n = 1, 7

    allocate ( alpha(1:n) )

    alpha(1:n) = 0
    call en_r2_test ( n, alpha )

    alpha(1:n) = 0
    alpha(1) = 2
    call en_r2_test ( n, alpha )

    alpha(1:n) = 0
    alpha(2) = 4
    call en_r2_test ( n, alpha )

    alpha(1:n) = 0
    i = mod ( 3 - 1, n ) + 1
    alpha(i) = 6
    call en_r2_test ( n, alpha )

    alpha(1:n) = 0
    alpha(1) = 2
    alpha(2) = 4
    call en_r2_test ( n, alpha )

    alpha(1:n) = 0
    i = mod ( 4 - 1, n ) + 1
    alpha(i) = 8
    call en_r2_test ( n, alpha )

    alpha(1:n) = 0
    i = mod ( 5 - 1, n ) + 1
    alpha(i) = 10
    call en_r2_test ( n, alpha )

    do i = 1, n
      alpha(i) = i
    end do
    call en_r2_test ( n, alpha )

    alpha(1:n) = 2
    call en_r2_test ( n, alpha )

    deallocate ( alpha )

  end do

  return
end
subroutine en_r2_test ( n, expon )

!*****************************************************************************80
!
!! EN_R2_TEST tests the Stroud EN_R2 rules on a monomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c1
  integer ( kind = 4 ) d
  real ( kind = 8 ) delta0
  real ( kind = 8 ) err
  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon(n)
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option
  integer ( kind = 4 ) p
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) quad
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:,:)

  write ( *, '(a)'      ) ' '
  write ( *, '(a,i4)'   ) '  N = ', n
  write ( *, '(a,10i4)' ) '  EXPON = ', expon(1:n)
  d = sum ( expon(1:n) )
  write ( *, '(a,i4)'   ) '  Degree = ', d
  write ( *, '(a)'      ) ' '

  call en_r2_monomial_integral ( n, expon, exact )

  p = 1

  if ( d <= p ) then

    call en_r2_01_1_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call en_r2_01_1 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_01_1:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 2

  if ( d <= p ) then

    call en_r2_02_xiu_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call en_r2_02_xiu ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_02_XIU: ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call gw_02_xiu_size ( n, o )
    gamma0 = 2.0D+00
    delta0 = 0.0D+00
    c1 = 1.0D+00
    volume_1d = sqrt ( pi )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  GW_02_XIU:    ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

  end if

  p = 3

  if ( d <= p ) then

    call en_r2_03_1_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call en_r2_03_1 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_03_1:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call en_r2_03_2_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call en_r2_03_2 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_03_2:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call en_r2_03_xiu_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call en_r2_03_xiu ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_03_XIU: ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

  end if

  p = 5

  if ( d <= p ) then

    if ( 2 <= n .and. n <= 7 ) then

      option = 1
      call en_r2_05_1_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_r2_05_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_05_1(1):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )

      if ( n == 3 .or. n == 5 .or. n == 6 ) then
        option = 2
        call en_r2_05_1_size ( n, option, o )
        allocate ( x(n,o) )
        allocate ( w(o) )
        call en_r2_05_1 ( n, option, o, x, w )
        allocate ( v(1:o) )
        call monomial_value ( n, o, x, expon, v )
        quad = dot_product ( w(1:o), v(1:o) )
        err = abs ( quad - exact )
        write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_05_1(2):', o, quad, err
        deallocate ( v )
        deallocate ( w )
        deallocate ( x )
      end if

    end if

    call en_r2_05_2_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call en_r2_05_2 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_05_2:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    if ( 3 <= n ) then
      call en_r2_05_3_size ( n, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_r2_05_3 ( n, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_05_3:   ', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

    call en_r2_05_4_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call en_r2_05_4 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_05_4:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call en_r2_05_5_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call en_r2_05_5 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_05_5:   ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    if ( 5 <= n ) then
      call en_r2_05_6_size ( n, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_r2_05_6 ( n, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_05_6:   ', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

  end if

  p = 7

  if ( d <= p ) then

    if ( n == 3 .or. n == 4 .or. n == 6 .or. n == 7 ) then
      option = 1
      call en_r2_07_1_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_r2_07_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_07_1(1):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

    if ( n == 3 .or. n == 4 ) then
      option = 2
      call en_r2_07_1_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_r2_07_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_07_1(2):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

    if ( 3 <= n ) then
      call en_r2_07_2_size ( n, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_r2_07_2 ( n, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_07_2:   ', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

    if ( 3 <= n .and. n <= 6 ) then
      option = 1
      call en_r2_07_3_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_r2_07_3 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_07_3(1):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )

      option = 2
      call en_r2_07_3_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_r2_07_3 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_07_3(2):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

  end if

  p = 9

  if ( d <= p ) then

    if ( 3 <= n .and. n <= 6 ) then
      option = 1
      call en_r2_09_1_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_r2_09_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_09_1(1):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )

      option = 2
      call en_r2_09_1_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_r2_09_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_09_1(2):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

  end if

  p = 11

  if ( d <= p ) then

    if ( 3 <= n .and. n <= 5 ) then
      option = 1
      call en_r2_11_1_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_r2_11_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_11_1(1):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )

      option = 2
      call en_r2_11_1_size ( n, option, o )
      allocate ( x(n,o) )
      allocate ( w(o) )
      call en_r2_11_1 ( n, option, o, x, w )
      allocate ( v(1:o) )
      call monomial_value ( n, o, x, expon, v )
      quad = dot_product ( w(1:o), v(1:o) )
      err = abs ( quad - exact )
      write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EN_R2_11_1(2):', o, quad, err
      deallocate ( v )
      deallocate ( w )
      deallocate ( x )
    end if

  end if

  write ( *, '(a,2x,6x,2x,g14.6)' ) '  EXACT:        ', exact

  return
end
subroutine test2075 ( )

!*****************************************************************************80
!
!! TEST2075 tests the rules for EPN with GLG weight on monomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 5

  real ( kind = 8 ) alpha
  real ( kind = 8 ), dimension ( test_num ) :: alpha_test = (/ &
    - 0.5D+00, 0.0D+00, 0.5D+00, 1.0D+00, 2.0D+00 /)
  integer ( kind = 4 ), allocatable :: expon(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2075'
  write ( *, '(a)' ) '  Demonstrate the use of quadrature rules for the region'
  write ( *, '(a)' ) '  EPN_GLG, that is, the positive half space [0,+oo)^N, with the'
  write ( *, '(a)' ) '  weight W(ALPHA;X) = product ( 1 <= I <= N ) X(I)^ALPHA exp ( -X(I) )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use the formulas to integrate various monomials of'
  write ( *, '(a)' ) '  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)'
  write ( *, '(a)' ) '  and compare to the exact integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The precision of each formula is known, and we only use'
  write ( *, '(a)' ) '  a formula if its precision indicates it should be able to'
  write ( *, '(a)' ) '  produce an exact result.'

  do n = 1, 6

    allocate ( expon(1:n) )

    expon(1:n) = 0
    do test = 1, test_num
      alpha = alpha_test(test)
      call epn_glg_test ( n, expon, alpha )
    end do

    expon(1:n) = 0
    expon(n) = 1
    do test = 1, test_num
      alpha = alpha_test(test)
      call epn_glg_test ( n, expon, alpha )
    end do

    if ( 2 <= n ) then
      expon(1:n) = 0
      expon(1) = 1
      expon(2) = 1
      do test = 1, test_num
        alpha = alpha_test(test)
        call epn_glg_test ( n, expon, alpha )
      end do
    end if

    expon(1:n) = 0
    expon(1) = 2
    do test = 1, test_num
      alpha = alpha_test(test)
      call epn_glg_test ( n, expon, alpha )
    end do

    deallocate ( expon )

  end do

  return
end
subroutine epn_glg_test ( n, expon, alpha )

!*****************************************************************************80
!
!! EPN_GLG_TEST tests the rules for EPN with GLG weight on a monomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) c1
  integer ( kind = 4 ) d
  real ( kind = 8 ) delta0
  real ( kind = 8 ) err
  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon(n)
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option
  integer ( kind = 4 ) p
  real ( kind = 8 ) quad
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:,:)

  write ( *, '(a)'      ) ' '
  write ( *, '(a,i4)'   ) '  N = ', n
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a,10i4)' ) '  EXPON = ', expon(1:n)
  d = sum ( expon(1:n) )
  write ( *, '(a,i4)'   ) '  Degree = ', d
  write ( *, '(a)'      ) ' '

  call epn_glg_monomial_integral ( n, expon, alpha, exact )

  p = 0

  if ( d <= p ) then
    call epn_glg_00_1_size ( n, alpha, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call epn_glg_00_1 ( n, alpha, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EPN_GLG_00_1:  ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 1

  if ( d <= p ) then
    call epn_glg_01_1_size ( n, alpha, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call epn_glg_01_1 ( n, alpha, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EPN_GLG_01_1:  ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 2

  if ( d <= p ) then

    call epn_glg_02_xiu_size ( n, alpha, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call epn_glg_02_xiu ( n, alpha, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EPN_GLG_02_XIU:', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call gw_02_xiu_size ( n, o )
    gamma0 = - 1.0D+00
    delta0 = alpha + 1.0D+00
    c1 = - alpha - 1.0D+00
    volume_1d = r8_gamma ( 1.0D+00 + alpha )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  GW_02_XIU:     ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

  end if

  write ( *, '(a,2x,6x,2x,g14.6)' ) '  EXACT:         ', exact

  return
end
subroutine test208 ( )

!*****************************************************************************80
!
!! TEST208 tests the rules for EPN with Laguerre weight on monomials.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable :: expon(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST208'
  write ( *, '(a)' ) '  Demonstrate the use of quadrature rules for the region'
  write ( *, '(a)' ) '  EPN_LAG, that is, the positive half space [0,+oo)^N, with the'
  write ( *, '(a)' ) '  weight W(X) = product ( 1 <= I <= N ) exp ( -X(I) )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use the formulas to integrate various monomials of'
  write ( *, '(a)' ) '  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)'
  write ( *, '(a)' ) '  and compare to the exact integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The precision of each formula is known, and we only use'
  write ( *, '(a)' ) '  a formula if its precision indicates it should be able to'
  write ( *, '(a)' ) '  produce an exact result.'

  do n = 1, 6

    allocate ( expon(1:n) )

    expon(1:n) = 0
    call epn_lag_test ( n, expon )

    expon(1:n) = 0
    expon(n) = 1
    call epn_lag_test ( n, expon )

    if ( 2 <= n ) then
      expon(1:n) = 0
      expon(1) = 1
      expon(2) = 1
      call epn_lag_test ( n, expon )
    end if

    expon(1:n) = 0
    expon(1) = 2
    call epn_lag_test ( n, expon )

    deallocate ( expon )

  end do

  return
end
subroutine epn_lag_test ( n, expon )

!*****************************************************************************80
!
!! EPN_LAG_TEST tests the rules for EPN with Laguerre weight on a monomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c1
  integer ( kind = 4 ) d
  real ( kind = 8 ) delta0
  real ( kind = 8 ) err
  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon(n)
  real ( kind = 8 ) gamma0
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option
  integer ( kind = 4 ) p
  real ( kind = 8 ) quad
  real ( kind = 8 ), allocatable :: v(:)
  real ( kind = 8 ) volume_1d
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:,:)

  write ( *, '(a)'      ) ' '
  write ( *, '(a,i4)'   ) '  N = ', n
  write ( *, '(a,10i4)' ) '  EXPON = ', expon(1:n)
  d = sum ( expon(1:n) )
  write ( *, '(a,i4)'   ) '  Degree = ', d
  write ( *, '(a)'      ) ' '

  call epn_lag_monomial_integral ( n, expon, exact )

  p = 0

  if ( d <= p ) then
    call epn_lag_00_1_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call epn_lag_00_1 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EPN_LAG_00_1:  ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 1

  if ( d <= p ) then
    call epn_lag_01_1_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call epn_lag_01_1 ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EPN_LAG_01_1:  ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )
  end if

  p = 2

  if ( d <= p ) then

    call epn_lag_02_xiu_size ( n, o )
    allocate ( x(n,o) )
    allocate ( w(o) )
    call epn_lag_02_xiu ( n, o, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  EPN_LAG_02_XIU:', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

    call gw_02_xiu_size ( n, o )
    gamma0 = - 1.0D+00
    delta0 = 1.0D+00
    c1 = - 1.0D+00
    volume_1d = 1.0D+00
    allocate ( x(n,o) )
    allocate ( w(o) )
    call gw_02_xiu ( n, o, gamma0, delta0, c1, volume_1d, x, w )
    allocate ( v(1:o) )
    call monomial_value ( n, o, x, expon, v )
    quad = dot_product ( w(1:o), v(1:o) )
    err = abs ( quad - exact )
    write ( *, '(a,2x,i6,2x,g14.6,2x,g14.6)' ) '  GW_02_XIU:     ', o, quad, err
    deallocate ( v )
    deallocate ( w )
    deallocate ( x )

  end if

  write ( *, '(a,2x,6x,2x,g14.6)' ) '  EXACT:         ', exact

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests HEXAGON_UNIT_SET and HEXAGON_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: rule_max = 4

  real ( kind = 8 ), dimension ( dim_num ) :: center = (/ &
    0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), external :: function_2d
  character ( len = 7 ) function_2d_name
  integer ( kind = 4 ) function_2d_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) num
  integer ( kind = 4 ) order
  real ( kind = 8 ) rad
  real ( kind = 8 ) result(rule_max)
  integer ( kind = 4 ) rule
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight
  real ( kind = 8 ), allocatable, dimension ( : ) :: xtab
  real ( kind = 8 ), allocatable, dimension ( : ) :: ytab

  rad = 2.0D+00
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  HEXAGON_UNIT_SET sets a quadrature rule for the '
  write ( *, '(a)' ) '  unit hexagon.'
  write ( *, '(a)' ) '  HEXAGON_SUM evaluates the quadrature rule'
  write ( *, '(a)' ) '  in an arbitrary hexagon.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  We use a radius ', rad
  write ( *, '(a)' ) '  and center:'
  write ( *, '(a,2g14.6)' ) '  CENTER = ', center(1:dim_num)
  write ( *, '(a)' ) ' '

  do ilo = 1, rule_max, 5

    ihi = min ( ilo + 4, rule_max )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a,5i6)' ) 'Rule:    ', ( rule, rule = ilo, ihi )
    write ( *, '(2x,a)' ) 'Function '
    write ( *, '(a)' ) ' '

    num = function_2d_num ( )

    do i = 1, num

      call function_2d_set ( 'SET', i )

      do rule = ilo, ihi

        call hexagon_unit_size ( rule, order )

        allocate ( xtab(1:order) )
        allocate ( ytab(1:order) )
        allocate ( weight(1:order) )

        call hexagon_unit_set ( rule, order, xtab, ytab, weight )

        call hexagon_sum ( function_2d, center, rad, order, xtab, ytab, weight, &
          result(rule) )

        deallocate ( xtab )
        deallocate ( ytab )
        deallocate ( weight )

      end do

      write ( *, '(2x,a7,5f14.8)' ) function_2d_name(i), result(ilo:ihi)

    end do

  end do

  return
end
subroutine test215 ( )

!*****************************************************************************80
!
!! TEST215 tests LENS_HALF_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) area
  real ( kind = 8 ), dimension ( dim_num ) :: center = (/ &
    0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), external :: f_1_2d
  real ( kind = 8 ), external :: f_x_2d
  real ( kind = 8 ), external :: f_r_2d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: i_max = 8
  real ( kind = 8 ) lens_half_2d
  real ( kind = 8 ) lens_half_area_2d
  integer ( kind = 4 ) order
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  real ( kind = 8 ) value

  r = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST215'
  write ( *, '(a)' ) '  LENS_HALF_2D approximates an integral within a'
  write ( *, '(a)' ) '  circular half lens, defined by joining the endpoints'
  write ( *, '(a)' ) '  of a circular arc.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrate F(X,Y) = 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      R            Theta1      Theta2        ' // &
    'Area        Order Integral'
  write ( *, '(a)' ) ' '

  do i = 0, i_max

    theta1 = 0.0D+00
    theta2 = real ( i, kind = 8 ) * 2.0D+00 * pi &
      / real ( i_max, kind = 8 )

    area = lens_half_area_2d ( r, theta1, theta2 )

    write ( *, '(a)' ) ' '

    do order = 2, 16, 2
      value = lens_half_2d ( f_1_2d, center, r, theta1, theta2, order )
      write ( *, '(4g14.6,i8,g14.6)' ) r, theta1, theta2, area, order, value
    end do

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrate F(X,Y) = X'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      R            Theta1      Theta2        ' // &
    'Area        Order Integral'
  write ( *, '(a)' ) ' '

  do i = 0, i_max

    theta1 = 0.0D+00
    theta2 = real ( i, kind = 8 ) * 2.0D+00 * pi &
      / real ( i_max, kind = 8 )

    area = lens_half_area_2d ( r, theta1, theta2 )

    write ( *, '(a)' ) ' '

    do order = 2, 16, 2
      value = lens_half_2d ( f_x_2d, center, r, theta1, theta2, order )
      write ( *, '(4g14.6,i8,g14.6)' ) r, theta1, theta2, area, order, value
    end do

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrate F(X,Y) = R'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      R            Theta1      Theta2        ' // &
    'Area        Order Integral'
  write ( *, '(a)' ) ' '

  do i = 0, i_max

    theta1 = 0.0D+00
    theta2 = real ( i, kind = 8 ) * 2.0D+00 * pi &
      / real ( i_max, kind = 8 )

    area = lens_half_area_2d ( r, theta1, theta2 )

    write ( *, '(a)' ) ' '

    do order = 2, 16, 2
      value = lens_half_2d ( f_r_2d, center, r, theta1, theta2, order )
      write ( *, '(4g14.6,i8,g14.6)' ) r, theta1, theta2, area, order, value
    end do

  end do

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests OCTAHEDRON_UNIT_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 3

  real ( kind = 8 ), external :: function_nd
  character ( len = 7 ) function_nd_name
  integer ( kind = 4 ) function_nd_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) num
  real ( kind = 8 ) result(3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  OCTAHEDRON_UNIT_ND approximates integrals in a unit' 
  write ( *, '(a)' ) '  octahedron in N dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    F(X)    N = 1    N = 2   N = 3 '
  write ( *, '(a)' ) ' '
 
  num = function_nd_num ( )

  do i = 1, num

    call function_nd_set ( 'SET', i )

    do n = 1, n_max
      call octahedron_unit_nd ( function_nd, n, result(n) )
    end do

    write ( *, '(2x,a7,3f14.8)' ) function_nd_name(i), result(1:n_max)

  end do

  return
end
subroutine test23 ( )

!*****************************************************************************80
!
!! TEST23 tests PARALLELIPIPED_VOLUME_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) parallelipiped_volume_nd
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: v
  real ( kind = 8 ) volume

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  PARALLELIPIPED_VOLUME_ND computes the volume of a'
  write ( *, '(a)' ) '  parallelipiped in N dimensions.'
  write ( *, '(a)' ) ' '

  do n = 2, 4

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Spatial dimension N = ', n

    allocate ( v(1:n,1:n+1) )
!
!  Set the values of the parallelipiped.
!
    call setsim ( n, v )

    call r8mat_print ( n, n + 1, v, '  Parallelipiped vertices:' )

    volume = parallelipiped_volume_nd ( n, v )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a,g14.6)' ) '  Volume is ', volume

    deallocate ( v )

  end do

  return
end
subroutine test24 ( )

!*****************************************************************************80
!
!! TEST24 tests POLYGON_**_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: npts = 4

  character ( len = 4 ) name
  real ( kind = 8 ) result
  real ( kind = 8 ), dimension ( npts) :: x = (/ &
    0.0D+00, 1.0D+00, 1.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( npts) :: y = (/ &
    0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  For a polygon in 2D:' 
  write ( *, '(a)' ) '  POLYGON_1_2D integrates 1'
  write ( *, '(a)' ) '  POLYGON_X_2D integrates X'
  write ( *, '(a)' ) '  POLYGON_Y_2D integrates Y'
  write ( *, '(a)' ) '  POLYGON_XX_2D integrates X^2'
  write ( *, '(a)' ) '  POLYGON_XY_2D integrates X*Y'
  write ( *, '(a)' ) '  POLYGON_YY_2D integrates Y^2'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  F(X,Y)    Integral'
  write ( *, '(a)' ) ' '
 
  name = '1'
  call polygon_1_2d ( npts, x, y, result )
  write ( *, '(2x,a4,4x,g14.6)' ) name, result
 
  name = 'X'
  call polygon_x_2d ( npts, x, y, result )
  write ( *, '(2x,a4,4x,g14.6)' ) name, result
 
  name = 'Y'
  call polygon_y_2d ( npts, x, y, result )
  write ( *, '(2x,a4,4x,g14.6)' ) name, result
 
  name = 'X^2'
  call polygon_xx_2d ( npts, x, y, result )
  write ( *, '(2x,a4,4x,g14.6)' ) name, result
 
  name = 'XY'
  call polygon_xy_2d ( npts, x, y, result )
  write ( *, '(2x,a4,4x,g14.6)' ) name, result
 
  name = 'Y^2'
  call polygon_yy_2d ( npts, x, y, result )
  write ( *, '(2x,a4,4x,g14.6)' ) name, result
 
  return
end
subroutine test25 ( )

!*****************************************************************************80
!
!! TEST25 tests PYRAMID_UNIT_O**_3D, PYRAMID_VOLUME_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: rule_num = 10

  character ( len = 7 ), allocatable, dimension ( : ) :: column_label
  character ( len = 7 ) function_3d_name
  integer ( kind = 4 ) function_3d_num
  real ( kind = 8 ), external :: function_3d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) num
  real ( kind = 8 ) pyramid_unit_volume_3d
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: result
  integer ( kind = 4 ), dimension ( rule_num ) :: order = (/ &
    1, 5, 6, 8, 8, 9, 13, 18, 27, 48 /)

  num = function_3d_num ( )

  allocate ( column_label(1:num) )
  allocate ( result(1:rule_num,1:num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) '  For the unit pyramid, we approximate integrals with:'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O01_3D, a 1 point rule.'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O05_3D, a 5 point rule.'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O06_3D, a 6 point rule.'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O08_3D, an 8 point rule.'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O08b_3D, an 8 point rule.'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O09_3D, a 9 point rule.'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O13_3D, a 13 point rule.'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O18_3D, a 18 point rule.'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O27_3D, a 27 point rule.'
  write ( *, '(a)' ) '  PYRAMID_UNIT_O48_3D, a 48 point rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PYRAMID_UNIT_VOLUME_3D computes the volume of a unit pyramid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Volume = ', pyramid_unit_volume_3d ( )
  write ( *, '(a)' ) ' '
    
  do j = 1, num

    call function_3d_set ( 'SET', j )

    column_label(j) = function_3d_name(j)

    call pyramid_unit_o01_3d ( function_3d, result(1,j) )
    call pyramid_unit_o05_3d ( function_3d, result(2,j) )
    call pyramid_unit_o06_3d ( function_3d, result(3,j) )
    call pyramid_unit_o08_3d ( function_3d, result(4,j) )
    call pyramid_unit_o08b_3d ( function_3d, result(5,j) )
    call pyramid_unit_o09_3d ( function_3d, result(6,j) )
    call pyramid_unit_o13_3d ( function_3d, result(7,j) )
    call pyramid_unit_o18_3d ( function_3d, result(8,j) )
    call pyramid_unit_o27_3d ( function_3d, result(9,j) )
    call pyramid_unit_o48_3d ( function_3d, result(10,j) )

  end do

  do jlo = 1, num, 5
    jhi = min ( jlo + 4, num )
    write ( *, '(a)' ) ' '
    write ( *, '(2x,a5,3x,a7,7x,a7,7x,a7,7x,a7,7x,a7)' ) &
      'Order', column_label(jlo:jhi)
    write ( *, '(a)' ) ' '
    do i = 1, rule_num
      write ( *, '(2x,i5,5g14.6)' ) order(i), result(i,jlo:jhi)
    end do
  end do

  deallocate ( column_label)
  deallocate ( result )

  return
end
subroutine test255 ( )

!*****************************************************************************80
!
!! TEST255 tests PYRAMID_UNIT_MONOMIAL_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) alpha
  integer ( kind = 4 ) beta
  integer ( kind = 4 ) :: degree_max = 4
  integer ( kind = 4 ) gamma
  real ( kind = 8 ) pyramid_unit_monomial_3d
  real ( kind = 8 ) pyramid_unit_volume_3d
  real ( kind = 8 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST255'
  write ( *, '(a)' ) '  For the unit pyramid,'
  write ( *, '(a)' ) '  PYRAMID_UNIT_MONOMIAL_3D returns the exact value of the'
  write ( *, '(a)' ) ' integral of X^ALPHA Y^BETA Z^GAMMA'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Volume = ', pyramid_unit_volume_3d ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     ALPHA      BETA     GAMMA      INTEGRAL'
  write ( *, '(a)' ) ' '

  do alpha = 0, degree_max
    do beta = 0, degree_max - alpha
      do gamma = 0, degree_max - alpha - beta
        value = pyramid_unit_monomial_3d ( alpha, beta, gamma )
        write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) alpha, beta, gamma, value
      end do
    end do
  end do

  return
end
subroutine test26 ( )

!*****************************************************************************80
!
!! TEST26 tests QMULT_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  character ( len = 7 ) function_1d_name
  real ( kind = 8 ), external :: function_1d
  integer ( kind = 4 ) function_1d_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) num
  real ( kind = 8 ) qmult_1d
  real ( kind = 8 ) result

  a = -1.0D+00
  b = 1.0D+00
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) '  QMULT_1D approximates an integral on a'
  write ( *, '(a)' ) '  one-dimensional interval.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use the interval:'
  write ( *, '(a,g14.6)' ) '  A = ', a
  write ( *, '(a,g14.6)' ) '  B = ', b
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    F(X)     QMULT_1D'
  write ( *, '(a)' ) ' '
 
  num = function_1d_num ( )

  do i = 1, num

    call function_1d_set ( 'SET', i )

    result = qmult_1d ( function_1d, a, b )
    write ( *, '(2x,a7,f14.8)' ) function_1d_name ( i ), result
 
  end do

  return
end
subroutine test27 ( )

!*****************************************************************************80
!
!! TEST27 tests SIMPLEX_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), external :: function_nd
  character ( len = 7 ) function_nd_name
  integer ( kind = 4 ) function_nd_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) num
  real ( kind = 8 ) result
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: v

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)' ) '  SIMPLEX_ND approximates integrals inside an' 
  write ( *, '(a)' ) '  arbitrary simplex in ND.'
  write ( *, '(a)' ) ' '
 
  do n = 2, 4
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Spatial dimension N = ', n

    allocate ( v(1:n,1:n+1) )
!
!  Restore values of simplex.
!
    call setsim ( n, v )

    call r8mat_print ( n, n + 1, v, '  Simplex vertices:' )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  F(X)    SIMPLEX_ND'
    write ( *, '(a)' ) ' '

    num = function_nd_num ( )

    do i = 1, num

      call function_nd_set ( 'SET', i )

      call simplex_nd ( function_nd, n, v, result )
      write ( *, '(2x,a7,f14.8)' ) function_nd_name(i), result
 
      call setsim ( n, v )

    end do

    deallocate ( v )

  end do
 
  return
end
subroutine test28 ( )

!*****************************************************************************80
!
!! TEST28 tests SIMPLEX_VOLUME_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) simplex_volume_nd
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: v
  real ( kind = 8 ) volume

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST28'
  write ( *, '(a)' ) '  SIMPLEX_VOLUME_ND computes the volume of a simplex'
  write ( *, '(a)' ) '  in N dimensions.'
  write ( *, '(a)' ) ' '

  do n = 2, 4

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Spatial dimension N = ', n

    allocate ( v(1:n,1:n+1) )
!
!  Set the values of the simplex.
!
    call setsim ( n, v )

    call r8mat_print ( n, n + 1, v, '  Simplex vertices:' )

    volume = simplex_volume_nd ( n, v )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Volume is ', volume

    deallocate ( v )

  end do

  return
end
subroutine test29 ( )

!*****************************************************************************80
!
!! TEST29 tests SIMPLEX_UNIT_**_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), external :: function_nd
  character ( len = 7 ) function_nd_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) result1
  real ( kind = 8 ) result2
  real ( kind = 8 ) result3
  real ( kind = 8 ) result4
  real ( kind = 8 ) simplex_unit_volume_nd
  real ( kind = 8 ) volume

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST29'
  write ( *, '(a)' ) '  For integrals in the unit simplex in ND,'
  write ( *, '(a)' ) '  SIMPLEX_UNIT_01_ND uses a formula of degree 1.'
  write ( *, '(a)' ) '  SIMPLEX_UNIT_03_ND uses a formula of degree 3.'
  write ( *, '(a)' ) '  SIMPLEX_UNIT_05_ND uses a formula of degree 5.'
  write ( *, '(a)' ) '  SIMPLEX_UNIT_05_2_ND uses a formula of degree 5.'

  do i = 1, 6

    call function_nd_set ( 'SET', i )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Check the integral of ', function_nd_name ( i )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  N     Volume         #1              #3' // &
      '              #5              #5.2'
    write ( *, '(a)' ) ' '

    do n = 2, 16

      call simplex_unit_01_nd ( function_nd, n, result1 )
      call simplex_unit_03_nd ( function_nd, n, result2 )
      call simplex_unit_05_nd ( function_nd, n, result3 )
      call simplex_unit_05_2_nd ( function_nd, n, result4 )

      volume = simplex_unit_volume_nd ( n )
      write ( *, '(2x,i2,2x,g13.5,2x,g13.5,2x,g13.5,2x,g13.5,2x,g13.5)' ) &
        n, volume, result1, result2, result3, result4

    end do

  end do

  return
end
subroutine test30 ( )

!*****************************************************************************80
!
!! TEST30 tests SPHERE_UNIT_**_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 7 ) function_3d_name
  integer ( kind = 4 ) function_3d_num
  real ( kind = 8 ), external :: function_3d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) num
  real ( kind = 8 ) result1
  real ( kind = 8 ) result2
  real ( kind = 8 ) result3
  real ( kind = 8 ) result4
  real ( kind = 8 ) sphere_unit_area_nd

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST30'
  write ( *, '(a)' ) '  For integrals on the unit sphere in 3D:'
  write ( *, '(a)' ) '  SPHERE_UNIT_07_3D uses a formula of degree 7.'
  write ( *, '(a)' ) '  SPHERE_UNIT_11_3D uses a formula of degree 11.'
  write ( *, '(a)' ) '  SPHERE_UNIT_14_3D uses a formula of degree 14.'
  write ( *, '(a)' ) '  SPHERE_UNIT_15_3D uses a formula of degree 15.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Unit sphere area = ', sphere_unit_area_nd ( 3 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    F(X)    S3S07        S3S11         S3S14' // &
    '         S3S15      '
  write ( *, '(a)' ) ' '
 
  num = function_3d_num ( )

  do i = 1, num

    call function_3d_set ( 'SET', i )

    call sphere_unit_07_3d ( function_3d, result1 ) 
    call sphere_unit_11_3d ( function_3d, result2 )
    call sphere_unit_14_3d ( function_3d, result3 )
    call sphere_unit_15_3d ( function_3d, result4 )
 
    write ( *, '(2x,a7,4f14.8)' ) &
      function_3d_name(i), result1, result2, result3, result4

  end do
 
  return
end
subroutine test31 ( )

!*****************************************************************************80
!
!! TEST31 tests SPHERE_UNIT_**_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), external :: function_nd
  character ( len = 7 ) function_nd_name
  integer ( kind = 4 ) function_nd_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) num
  real ( kind = 8 ) result1
  real ( kind = 8 ) result2
  real ( kind = 8 ) result3
  real ( kind = 8 ) result4
  real ( kind = 8 ) result5
  real ( kind = 8 ) result6
  real ( kind = 8 ) sphere_unit_area_nd

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST31'
  write ( *, '(a)' ) '  For integrals on the unit sphere in ND:'
  write ( *, '(a)' ) '  SPHERE_UNIT_03_ND uses a formula of degree 3;'
  write ( *, '(a)' ) '  SPHERE_UNIT_04_ND uses a formula of degree 4;'
  write ( *, '(a)' ) '  SPHERE_UNIT_05_ND uses a formula of degree 5.'
  write ( *, '(a)' ) '  SPHERE_UNIT_07_1_ND uses a formula of degree 7.'
  write ( *, '(a)' ) '  SPHERE_UNIT_07_2_ND uses a formula of degree 7.'
  write ( *, '(a)' ) '  SPHERE_UNIT_11_ND uses a formula of degree 11.'
  write ( *, '(a)' ) ' '

  do n = 3, 10

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Spatial dimension N = ', n
    write ( *, '(a,g14.6)' ) '  Unit sphere area = ', sphere_unit_area_nd ( n )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    Rule:     #3            #4            #5'          
    write ( *, '(a)' ) '              #7.1          #7.2          #11'
    write ( *, '(a)' ) '    Function'
    write ( *, '(a)' ) ' '

    num = function_nd_num ( )

    do i = 1, num

      call function_nd_set ( 'SET', i )

      call sphere_unit_03_nd ( function_nd, n, result1 )
      call sphere_unit_04_nd ( function_nd, n, result2 )
      call sphere_unit_05_nd ( function_nd, n, result3 )
      call sphere_unit_07_1_nd ( function_nd, n, result4 )
      call sphere_unit_07_2_nd ( function_nd, n, result5 )
      call sphere_unit_11_nd ( function_nd, n, result6 )

      write ( *, '(2x,a7,3f14.8)' ) function_nd_name(i), result1, result2, result3
      write ( *, '(2x,7x,3f14.8)' )           result4, result5, result6

    end do

  end do

  return
end
subroutine test32 ( )

!*****************************************************************************80
!
!! TEST32 tests SPHERE_05_ND, SPHERE_07_1_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 5

  real ( kind = 8 ) center(n_max)
  real ( kind = 8 ), external :: function_nd
  character ( len = 7 ) function_nd_name
  integer ( kind = 4 ) function_nd_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) num
  real ( kind = 8 ) r
  real ( kind = 8 ) result1
  real ( kind = 8 ) result2
  real ( kind = 8 ) sphere_area_nd

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST32'
  write ( *, '(a)' ) '  For integrals on a sphere in ND:'
  write ( *, '(a)' ) '  SPHERE_05_ND uses a formula of degree 5.'
  write ( *, '(a)' ) '  SPHERE_07_1_ND uses a formula of degree 7.'
  write ( *, '(a)' ) ' '
  r = 2.0D+00
  center(1:n_max) = 1.0D+00

  do n = 2, 4

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Spatial dimension N = ', n
    write ( *, '(a)' ) '  Sphere center = '
    write ( *, '(5g14.6)' ) center(1:n)
    write ( *, '(a,g14.6)' ) '  Sphere radius = ', r
    write ( *, '(a,g14.6)' ) '  Sphere area = ', sphere_area_nd ( n, r )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    Rule:     #5           #7.1'
    write ( *, '(a)' ) '    Function'
    write ( *, '(a)' ) ' '

    num = function_nd_num ( )

    do i = 1, num

      call function_nd_set ( 'SET', i )

      call sphere_05_nd ( function_nd, n, center, r, result1 )
      call sphere_07_1_nd ( function_nd, n, center, r, result2 )

      write ( *, '(2x,a7,2f14.8)' ) function_nd_name(i), result1, result2

    end do

  end do

  return
end
subroutine test322 ( )

!*****************************************************************************80
!
!! TEST322 tests SPHERE_CAP_AREA_3D, SPHERE_CAP_AREA_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) area1
  real ( kind = 8 ) area2
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: test_num = 12
  real ( kind = 8 ) :: r = 1.0D+00
  real ( kind =8 ) sphere_area_3d

  center(1:3) = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST322'
  write ( *, '(a)' ) '  SPHERE_CAP_AREA_3D computes the volume of a'
  write ( *, '(a)' ) '  3D spherical cap, defined by a plane that cuts the'
  write ( *, '(a)' ) '  sphere to a thickness of H units.'
  write ( *, '(a)' ) '  SPHERE_CAP_AREA_ND computes the volume of an'
  write ( *, '(a)' ) '  ND spherical cap, defined by a plane that cuts the'
  write ( *, '(a)' ) '  sphere to a thickness of H units.'

  area1 = sphere_area_3d ( r )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Area of the total sphere in 3D = ', area1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        R           H           Cap         Cap'
  write ( *, '(a)' ) '                                area_3d     area_nd'
  write ( *, '(a)' ) ' '

  do i = 0, test_num+1

    h = 2.0D+00 * r * real ( i, kind = 8 ) / real ( test_num, kind = 8 )

    call sphere_cap_area_3d ( r, h, area1 )

    call sphere_cap_area_nd ( dim_num, r, h, area2 )

    write ( *, '(2x,5f12.6)' ) r, h, area1, area2

  end do

  return
end
subroutine test324 ( )

!*****************************************************************************80
!
!! TEST324 tests SPHERE_CAP_VOLUME_2D, SPHERE_CAP_VOLUME_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: test_num = 12
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) :: r = 1.0D+00
  real ( kind = 8 ) volume1
  real ( kind = 8 ) volume2

  center(1:2) = (/ 0.0D+00, 0.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST324'
  write ( *, '(a)' ) '  SPHERE_CAP_VOLUME_2D computes the volume (area) of a'
  write ( *, '(a)' ) '  spherical cap, defined by a plane that cuts the'
  write ( *, '(a)' ) '  sphere to a thickness of H units.'
  write ( *, '(a)' ) '  SPHERE_CAP_VOLUME_ND does the same operation,'
  write ( *, '(a)' ) '  but in N dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Using a radius R = ', r

  call sphere_volume_2d ( r, volume1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Volume of the total sphere in 2D = ', volume1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        H           Cap        Cap'
  write ( *, '(a)' ) '                    vol_2d     vol_nd'
  write ( *, '(a)' ) ' '

  do i = 0, test_num+1

    h = 2.0D+00 * r * real ( i, kind = 8 ) / real ( test_num, kind = 8 )

    call sphere_cap_volume_2d ( r, h, volume1 )

    call sphere_cap_volume_nd ( dim_num, r, h, volume2 )

    write ( *, '(2x,3f12.6)' ) h, volume1, volume2

  end do

  return
end
subroutine test326 ( )

!*****************************************************************************80
!
!! TEST326 tests SPHERE_CAP_VOLUME_3D, SPHERE_CAP_VOLUME_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: test_num = 12
  real ( kind = 8 ) :: r = 1.0D+00
  real ( kind = 8 ) volume1
  real ( kind = 8 ) volume2

  center(1:3) = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST326'
  write ( *, '(a)' ) '  SPHERE_CAP_VOLUME_3D computes the volume of a'
  write ( *, '(a)' ) '  spherical cap, defined by a plane that cuts the'
  write ( *, '(a)' ) '  sphere to a thickness of H units.'
  write ( *, '(a)' ) '  SPHERE_CAP_VOLUME_ND does the same operation,'
  write ( *, '(a)' ) '  but in N dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Using a radius R = ', r

  call sphere_volume_3d ( r, volume1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Volume of the total sphere in 3D = ', volume1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        H           Cap        Cap'
  write ( *, '(a)' ) '                    volume_3d  volume_nd'
  write ( *, '(a)' ) ' '

  do i = 0, test_num+1

    h = 2.0D+00 * r * real ( i, kind = 8 ) / real ( test_num, kind = 8 )

    call sphere_cap_volume_3d ( r, h, volume1 )

    call sphere_cap_volume_nd ( dim_num, r, h, volume2 )

    write ( *, '(2x,3f12.6)' ) h, volume1, volume2

  end do

  return
end
subroutine test33 ( )

!*****************************************************************************80
!
!! TEST33 tests SPHERE_CAP_AREA_ND, SPHERE_CAP_VOLUME_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: n = 12
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) sphere_area_nd
  real ( kind = 8 ) volume

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST33'
  write ( *, '(a)' ) '  For a sphere in ND:'
  write ( *, '(a)' ) '  SPHERE_CAP_AREA_ND computes the area '
  write ( *, '(a)' ) '  of a spherical cap.'
  write ( *, '(a)' ) '  SPHERE_CAP_VOLUME_ND computes the volume '
  write ( *, '(a)' ) '  of a spherical cap.'
  write ( *, '(a)' ) ' '

  r = 1.0D+00

  do dim_num = 2, 5

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Spatial dimension N = ', dim_num
    write ( *, '(a,g14.6)' ) '  Radius =       ', r
    write ( *, '(a,g14.6)' ) '  Area =         ', sphere_area_nd ( dim_num, r )
    call sphere_volume_nd ( dim_num, r, volume )
    write ( *, '(a,g14.6)' ) '  Volume =       ', volume

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '                 Sphere         Sphere' 
    write ( *, '(a)' ) '                 cap            cap'
    write ( *, '(a)' ) '      H          area           volume'
    write ( *, '(a)' ) ' '

    do i = 0, n + 1
      h = real ( 2 * i, kind = 8 ) * r / real ( n, kind = 8 )
      call sphere_cap_area_nd ( dim_num, r, h, area )
      call sphere_cap_volume_nd ( dim_num, r, h, volume )
      write ( *, '(2x,f8.4,2x,g14.6,2x,g14.6)' ) h, area, volume
    end do

  end do

  return
end
subroutine test335 ( )

!*****************************************************************************80
!
!! TEST335 tests SPHERE_SHELL_03_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 3

  real ( kind = 8 ) center(n_max)
  real ( kind = 8 ), external :: function_nd
  character ( len = 7 ) function_nd_name
  integer ( kind = 4 ) function_nd_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) num
  real ( kind = 8 ) result1
  real ( kind = 8 ) result3
  real ( kind = 8 ) result4
  real ( kind = 8 ) result5
  real ( kind = 8 ) result6
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) sphere_shell_volume_nd

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST335'
  write ( *, '(a)' ) '  For integrals inside a spherical shell in ND:'
  write ( *, '(a)' ) '  SPHERE_SHELL_03_ND approximates the integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We compare these results with those computed by'
  write ( *, '(a)' ) '  from the difference of two ball integrals:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BALL_F1_ND approximates the integral;'
  write ( *, '(a)' ) '  BALL_F3_ND approximates the integral.'
  write ( *, '(a)' ) ' '

  do j = 1, 2

    if ( j == 1 ) then
      r1 = 0.0D+00
      r2 = 1.0D+00
      center(1:n_max) = 0.0D+00
    else
      r1 = 2.0D+00
      r2 = 3.0D+00
      center(1:n_max) = (/ 1.0D+00, -1.0D+00, 2.0D+00 /)
    end if

    do n = 2, n_max

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Spatial dimension N = ', n
      write ( *, '(a)' ) '  Sphere center:'
      write ( *, '(2x,3f10.4)' ) center(1:n)
      write ( *, '(a,g14.6)' ) '  Inner sphere radius = ', r1
      write ( *, '(a,g14.6)' ) '  Outer sphere radius = ', r2
      write ( *, '(a,g14.6)' ) '  Spherical shell volume = ', &
        sphere_shell_volume_nd ( n, r1, r2 )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    Rule:      #3       F1(R2)-F1(R1)  ' // &      
        'F3(R2)-F3(R1)'
      write ( *, '(a)' ) '    F(X)'
      write ( *, '(a)' ) ' '

      num = function_nd_num ( )

      do i = 1, num

        call function_nd_set ( 'SET', i )

        call sphere_shell_03_nd ( function_nd, n, center, r1, r2, result1 )

        call ball_f1_nd ( function_nd, n, center, r1, result3 )
        call ball_f1_nd ( function_nd, n, center, r2, result4 )
      
        call ball_f3_nd ( function_nd, n, center, r1, result5 )
        call ball_f3_nd ( function_nd, n, center, r2, result6 )

        write ( *, '(2x,a7,4f14.8)' ) function_nd_name(i), result1, &
          result4-result3, result6-result5

      end do

    end do

  end do

  return
end
subroutine test34 ( )

!*****************************************************************************80
!
!! TEST34 tests SPHERE_UNIT_AREA_ND, SPHERE_UNIT_AREA_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ) area2
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sphere_unit_area_nd

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST34:'
  write ( *, '(a)' ) '  SPHERE_UNIT_AREA_ND evaluates the area of the unit'
  write ( *, '(a)' ) '  sphere in N dimensions.'
  write ( *, '(a)' ) '  SPHERE_UNIT_AREA_VALUES returns some test values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     dim_num    Exact          Computed'
  write ( *, '(a)' ) '             Area           Area'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sphere_unit_area_values ( n_data, dim_num, area )

    if ( n_data == 0 ) then
      exit
    end if

    area2 = sphere_unit_area_nd ( dim_num )

    write ( *, '(2x,i8,2x,f10.6,2x,f10.6)' ) dim_num, area, area2

  end do

  return
end
subroutine test345 ( )

!*****************************************************************************80
!
!! TEST345 tests SPHERE_UNIT_VOLUME_ND, SPHERE_UNIT_VOLUME_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sphere_unit_volume_nd
  real ( kind = 8 ) volume
  real ( kind = 8 ) volume2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST345:'
  write ( *, '(a)' ) '  SPHERE_UNIT_VOLUME_ND evaluates the area of the unit'
  write ( *, '(a)' ) '  sphere in N dimensions.'
  write ( *, '(a)' ) '  SPHERE_UNIT_VOLUME_VALUES returns some test values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     dim_num    Exact          Computed'
  write ( *, '(a)' ) '             Volume         Volume'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sphere_unit_volume_values ( n_data, dim_num, volume )

    if ( n_data == 0 ) then
      exit
    end if

    volume2 = sphere_unit_volume_nd ( dim_num )

    write ( *, '(2x,i8,2x,f10.6,2x,f10.6)' ) dim_num, volume, volume2

  end do

  return
end
subroutine test35 ( )

!*****************************************************************************80
!
!! TEST35 tests SQUARE_UNIT_SET, RECTANGLE_SUB_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), external :: function_2d
  character ( len = 7 ) function_2d_name
  integer ( kind = 4 ) function_2d_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) num
  integer ( kind = 4 ) order
  integer ( kind = 4 ) nsub(2)
  real ( kind = 8 ) result
  integer ( kind = 4 ) rule
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight
  real ( kind = 8 ), allocatable, dimension ( : ) :: xtab
  real ( kind = 8 ) xval(2)
  real ( kind = 8 ), allocatable, dimension ( : ) :: ytab
  real ( kind = 8 ) yval(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST35'
  write ( *, '(a)' ) '  SQUARE_UNIT_SET sets up a quadrature rule '
  write ( *, '(a)' ) '  on a unit square.'
  write ( *, '(a)' ) '  RECTANGLE_SUB_2D applies it to subrectangles of an'
  write ( *, '(a)' ) '  arbitrary rectangle.'
  write ( *, '(a)' ) ' '
!
!  Set the location of the square.
!
  xval(1) = 1.0D+00
  yval(1) = 2.0D+00

  xval(2) = 3.0D+00
  yval(2) = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The corners of the rectangle are:'
  write ( *, '(a)' ) ' '
  write ( *, '(2g14.6)' ) xval(1), yval(1)
  write ( *, '(2g14.6)' ) xval(2), yval(2)
!
!  Get the quadrature abscissas and weights for a unit square.
!
  rule = 2
  call square_unit_size ( rule, order )

  allocate ( xtab(1:order) )
  allocate ( ytab(1:order) )
  allocate ( weight(1:order) )

  call square_unit_set ( rule, order, xtab, ytab, weight )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Using unit square integration rule number ', rule
  write ( *, '(a,i8)' ) '  Order of rule is ', order
!
!  Set the function.
!
  num = function_2d_num ( )

  do i = 1, num

    call function_2d_set ( 'SET', i )
!
!  Try an increasing number of subdivisions.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    Function  Subdivisions  Integral'
    write ( *, '(a)' ) ' '

    do j = 1, 5

      nsub(1) = j
      nsub(2) = 2 * j

      call rectangle_sub_2d ( function_2d, xval, yval, nsub, order, xtab, ytab, &
        weight, result )

      write ( *, '(2x,a7,2i4, f14.8)' ) function_2d_name(i), nsub(1), nsub(2), result

    end do

  end do

  deallocate ( weight )
  deallocate ( xtab )
  deallocate ( ytab )

  return
end
subroutine test36 ( )

!*****************************************************************************80
!
!! TEST36 tests SQUARE_UNIT_SET and SQUARE_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: rule_max = 6

  real ( kind = 8 ), dimension ( dim_num ) :: center = (/ 2.0D+00, 2.0D+00 /)
  real ( kind = 8 ), external :: function_2d
  character ( len = 7 ) function_2d_name
  integer ( kind = 4 ) function_2d_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) num
  integer ( kind = 4 ) order
  real ( kind = 8 ) r
  real ( kind = 8 ) result(rule_max)
  integer ( kind = 4 ) rule
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight
  real ( kind = 8 ), allocatable, dimension ( : ) :: xtab
  real ( kind = 8 ), allocatable, dimension ( : ) :: ytab

  r = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST36'
  write ( *, '(a)' ) '  SQUARE_UNIT_SET sets up quadrature on the unit square;'
  write ( *, '(a)' ) '  SQUARE_SUM carries it out on an arbitrary square.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Square center:'
  write ( *, '(a,2g14.6)' ) '  CENTER = ', center(1:dim_num)
  write ( *, '(a,g14.6)' ) '  Square radius is ', r

  do ilo = 1, rule_max, 5

    ihi = min ( ilo + 4, rule_max )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a,5(i6,7x))' ) 'Rule:  ', ( rule, rule = ilo, ihi )
    write ( *, '(a)' ) '  Function '
    write ( *, '(a)' ) ' '

    num = function_2d_num ( )

    do i = 1, num

      call function_2d_set ( 'SET', i )

      do rule = ilo, ihi

        call square_unit_size ( rule, order )

        allocate ( xtab(1:order) )
        allocate ( ytab(1:order) )
        allocate ( weight(1:order) )

        call square_unit_set ( rule, order, xtab, ytab, weight )

        call square_sum ( function_2d, center, r, order, xtab, ytab, weight, &
          result(rule) )

        deallocate ( weight )
        deallocate ( xtab )
        deallocate ( ytab )

      end do

      write ( *, '(2x,a7,5f13.7)' ) function_2d_name(i), result(ilo:ihi)

    end do

  end do
 
  return
end
subroutine test37 ( )

!*****************************************************************************80
!
!! TEST37 tests SQUARE_UNIT_SET and SQUARE_UNIT_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: rule_max = 6

  real ( kind = 8 ), external :: function_2d
  character ( len = 7 ) function_2d_name
  integer ( kind = 4 ) function_2d_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) num
  integer ( kind = 4 ) order
  real ( kind = 8 ) result(rule_max)
  integer ( kind = 4 ) rule
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight
  real ( kind = 8 ), allocatable, dimension ( : ) :: xtab
  real ( kind = 8 ), allocatable, dimension ( : ) :: ytab

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST37'
  write ( *, '(a)' ) '  SQUARE_UNIT_SET sets up quadrature on the unit square;'
  write ( *, '(a)' ) '  SQUARE_UNIT_SUM carries it out on the unit square.'
  write ( *, '(a)' ) ' '
 
  do ilo = 1, rule_max, 5

    ihi = min ( ilo + 4, rule_max )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a,5i6)' ) 'Rule:    ', ( rule, rule = ilo, ihi )
    write ( *, '(2x,a)' ) 'Function '
    write ( *, '(a)' ) ' '

    num = function_2d_num ( )

    do i = 1, num

      call function_2d_set ( 'SET', i )

      do rule = ilo, ihi

        call square_unit_size ( rule, order )

        allocate ( xtab(1:order) )
        allocate ( ytab(1:order) )
        allocate ( weight(1:order) )

        call square_unit_set ( rule, order, xtab, ytab, weight )

        call square_unit_sum ( function_2d, order, xtab, ytab, weight, &
          result(rule) )

        deallocate ( weight )
        deallocate ( xtab )
        deallocate ( ytab )

      end do

      write ( *, '(2x,a7,5f13.7)' ) function_2d_name(i), result(ilo:ihi)

    end do

  end do
 
  return
end
subroutine test38 ( )

!*****************************************************************************80
!
!! TEST38 tests TETRA_07, TETRA_TPRODUCT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 9

  character ( len = 7 ) function_3d_name
  integer ( kind = 4 ) function_3d_num
  real ( kind = 8 ), external :: function_3d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) num
  integer ( kind = 4 ) order
  real ( kind = 8 ) result2
  real ( kind = 8 ) result3(order_max)
  real ( kind = 8 ) tetra_unit_volume
  real ( kind = 8 ) tetra_volume
  real ( kind = 8 ), dimension ( 4 ) :: x = (/ &
    1.0D+00, 4.0D+00, 1.0D+00, 1.0D+00 /)
  real ( kind = 8 ), dimension ( 4 ) :: y = (/ &
    2.0D+00, 2.0D+00, 3.0D+00, 2.0D+00 /)
  real ( kind = 8 ), dimension ( 4 ) :: z = (/ &
    6.0D+00, 6.0D+00, 6.0D+00, 8.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST38'
  write ( *, '(a)' ) '  For integrals inside an arbitrary tetrahedron:'
  write ( *, '(a)' ) '  TETRA_07 uses a formula of degree 7;'
  write ( *, '(a)' ) '  TETRA_TPRODUCT uses a triangular product formula '
  write ( *, '(a)' ) '  of varying degree.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tetrahedron vertices:'
  write ( *, '(a)' ) ' '
  do i = 1, 4
    write ( *, '(2x,3f4.0)' ) x(i), y(i), z(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Tetrahedron unit volume = ', tetra_unit_volume ( )
  write ( *, '(a,g14.6)' ) '  Tetrahedron Volume = ', tetra_volume ( x, y, z )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  F(X)    TETRA_07'
  write ( *, '(a)' ) '          TETRA_TPRODUCT(1:4)'
  write ( *, '(a)' ) '          TETRA_TPRODUCT(5:8)'
  write ( *, '(a)' ) '          TETRA_TPRODUCT(9)'
  write ( *, '(a)' ) ' '

  num = function_3d_num ( )

  do i = 1, num

    call function_3d_set ( 'SET', i )

    call tetra_07 ( function_3d, x, y, z, result2 )

    do order = 1, order_max
      call tetra_tproduct ( function_3d, order, x, y, z, result3(order) )
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a7,4f16.10)' ) function_3d_name(i), result2
    write ( *, '(2x,7x,4f16.10)' )                  result3(1:4)
    write ( *, '(2x,7x,4f16.10)' )                  result3(5:8)
    write ( *, '(2x,7x,4f16.10)' )                  result3(9)

  end do

  return
end
subroutine test39 ( )

!*****************************************************************************80
!
!! TEST39 tests TETRA_UNIT_SET and TETRA_UNIT_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: rule_max = 8

  character ( len = 7 ) function_3d_name
  integer ( kind = 4 ) function_3d_num
  real ( kind = 8 ), external :: function_3d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) num
  integer ( kind = 4 ) order
  real ( kind = 8 ) result(rule_max)
  integer ( kind = 4 ) rule
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight
  real ( kind = 8 ), allocatable, dimension ( : ) :: xtab
  real ( kind = 8 ), allocatable, dimension ( : ) :: ytab
  real ( kind = 8 ), allocatable, dimension ( : ) :: ztab

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST39'
  write ( *, '(a)' ) '  TETRA_UNIT_SET sets quadrature rules '
  write ( *, '(a)' ) '  for the unit tetrahedron;'
  write ( *, '(a)' ) '  TETRA_UNIT_SUM applies them to the unit tetrahedron.'
  write ( *, '(a)' ) ' '

  do ilo = 1, rule_max, 5

    ihi = min ( ilo +  4, rule_max )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a,5i6)' ) 'Rule:   ', ( rule, rule = ilo, ihi )
    write ( *, '(2x,a)' ) 'Function'
    write ( *, '(a)' ) ' '

    num = function_3d_num ( )

    do i = 1, num

      call function_3d_set ( 'SET', i )

      do rule = ilo, ihi

        call tetra_unit_size ( rule, order )

        allocate ( weight(1:order) )
        allocate ( xtab(1:order) )
        allocate ( ytab(1:order) )
        allocate ( ztab(1:order) )

        call tetra_unit_set ( rule, order, xtab, ytab, ztab, weight )
 
        call tetra_unit_sum ( function_3d, order, xtab, ytab, ztab, weight, &
          result(rule) )

        deallocate ( weight )
        deallocate ( xtab )
        deallocate ( ytab )
        deallocate ( ztab )

      end do

      write ( *, '(2x,a7,5f14.6)' ) function_3d_name(i), result(ilo:ihi)

    end do

  end do

  return
end
subroutine test40 ( )

!*****************************************************************************80
!
!! TEST40 tests TETRA_UNIT_SET and TETRA_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: rule_max = 8

  character ( len = 7 ) function_3d_name
  integer ( kind = 4 ) function_3d_num
  real ( kind = 8 ), external :: function_3d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) num
  integer ( kind = 4 ) order
  real ( kind = 8 ) result(rule_max)
  integer ( kind = 4 ) rule
  real ( kind = 8 ) value
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight
  real ( kind = 8 ), dimension ( 4 ) :: x = (/ &
    1.0D+00, 4.0D+00, 1.0D+00, 1.0D+00 /)
  real ( kind = 8 ), allocatable, dimension ( : ) :: xtab
  real ( kind = 8 ), dimension ( 4 ) :: y = (/ &
    2.0D+00, 2.0D+00, 3.0D+00, 2.0D+00 /)
  real ( kind = 8 ), allocatable, dimension ( : ) :: ytab
  real ( kind = 8 ), dimension ( 4 ) :: z = (/ &
    6.0D+00, 6.0D+00, 6.0D+00, 8.0D+00 /)
  real ( kind = 8 ), allocatable, dimension ( : ) :: ztab

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST40'
  write ( *, '(a)' ) '  TETRA_UNIT_SET sets quadrature rules '
  write ( *, '(a)' ) '  for the unit tetrahedron;'
  write ( *, '(a)' ) '  TETRA_SUM applies them to an arbitrary tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tetrahedron vertices:'
  write ( *, '(a)' ) ' '
  do i = 1, 4
    write ( *, '(2x,3f6.2)' ) x(i), y(i), z(i)
  end do

  do ilo = 1, rule_max, 5

    ihi = min ( ilo +  4, rule_max )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a,5(i7,7x))' ) 'Rule:    ', ( rule, rule = ilo, ihi )
    write ( *, '(2x,a)' ) 'Function '
    write ( *, '(a)' ) ' '
!
!  Believe it or not, if I remove these unused debugging statements,
!  the program will fail under the G95 compiler.  Somehow, somewhere,
!  there must be a memory error, but I have not figured out where it is.
!
    num = function_3d_num ( )

    do i = 1, num

      call function_3d_set ( 'SET', i )

      do rule = ilo, ihi

        call tetra_unit_size ( rule, order )

        allocate ( weight(1:order) )
        allocate ( xtab(1:order) )
        allocate ( ytab(1:order) )
        allocate ( ztab(1:order) )

        call tetra_unit_set ( rule, order, xtab, ytab, ztab, weight )

        call tetra_sum ( function_3d, x, y, z, order, xtab, ytab, ztab, weight, &
          value )

        result(rule) = value

        deallocate ( weight )
        deallocate ( xtab )
        deallocate ( ytab )
        deallocate ( ztab )

      end do

      write ( *, '(2x,a7,5g14.7)' ) function_3d_name(i), result(ilo:ihi)

    end do

  end do

  return
end
subroutine test41 ( )

!*****************************************************************************80
!
!! TEST41 tests TRIANGLE_UNIT_SET, TRIANGLE_SUB.
!
!  Discussion:
!
!    Break up the triangle into NSUB*NSUB equal subtriangles.  Approximate 
!    the integral over the triangle by the sum of the integrals over each
!    subtriangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), external :: function_2d
  character ( len = 7 ) function_2d_name
  integer ( kind = 4 ) function_2d_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) num
  integer ( kind = 4 ) order
  integer ( kind = 4 ) nsub
  real ( kind = 8 ) result
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) triangle_unit_size
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight
  real ( kind = 8 ), allocatable, dimension ( : ) :: xtab
  real ( kind = 8 ), dimension ( 3 ) :: xval = (/ &
    0.0D+00, 0.0D+00, 1.0D+00 /)
  real ( kind = 8 ), allocatable, dimension ( : ) :: ytab
  real ( kind = 8 ), dimension ( 3 ) :: yval = (/ &
    0.0D+00, 1.0D+00, 0.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST41'
  write ( *, '(a)' ) '  TRIANGLE_UNIT_SET sets up a quadrature rule '
  write ( *, '(a)' ) '  on a triangle.'
  write ( *, '(a)' ) '  TRIANGLE_SUB applies it to subtriangles of an'
  write ( *, '(a)' ) '  arbitrary triangle.'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Triangle vertices:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,2g14.6)' ) xval(1), yval(1)
  write ( *, '(2x,2g14.6)' ) xval(2), yval(2)
  write ( *, '(2x,2g14.6)' ) xval(3), yval(3)
!
!  Get the quadrature abscissas and weights for a unit triangle.
!
  rule = 3
  order = triangle_unit_size ( rule )

  allocate ( xtab(1:order) )
  allocate ( ytab(1:order) )
  allocate ( weight(1:order) )

  call triangle_unit_set ( rule, order, xtab, ytab, weight )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Using unit triangle quadrature rule ', rule
  write ( *, '(a,i8)' ) '  Rule order = ', order
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Function Nsub  Result'
  write ( *, '(a)' ) ' '
!
!  Set the function.
!
  num = function_2d_num ( )

  do i = 1, num

    call function_2d_set ( 'SET', i )
!
!  Try an increasing number of subdivisions.
!
    do nsub = 1, 5

      call triangle_sub ( function_2d, xval, yval, nsub, order, xtab, ytab,  &
        weight, result )

      write ( *, '(2x,a7,i4, f14.8)' ) function_2d_name(i), nsub, result
 
    end do
  
  end do

  deallocate ( xtab )
  deallocate ( ytab )
  deallocate ( weight )

  return
end
subroutine test42 ( )

!*****************************************************************************80
!
!! TEST42 tests TRIANGLE_UNIT_SET and TRIANGLE_UNIT_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: rule_max = 20

  real ( kind = 8 ), external :: function_2d
  character ( len = 7 ) function_2d_name
  integer ( kind = 4 ) function_2d_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) num
  integer ( kind = 4 ) order
  real ( kind = 8 ) result(rule_max)
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) triangle_unit_size
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight
  real ( kind = 8 ), allocatable, dimension ( : ) :: xtab
  real ( kind = 8 ), allocatable, dimension ( : ) :: ytab

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST42'
  write ( *, '(a)' ) '  TRIANGLE_UNIT_SET sets up a quadrature '
  write ( *, '(a)' ) '  in the unit triangle,'
  write ( *, '(a)' ) '  TRIANGLE_UNIT_SUM applies it.'
  write ( *, '(a)' ) ' '

  do ilo = 1, rule_max, 5

    ihi = min ( ilo + 4, rule_max )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a,5(i7,7x))' ) 'Rule:    ', ( rule, rule = ilo, ihi )
    write ( *, '(2x,a)' ) 'Function '
    write ( *, '(a)' ) ' '

    num = function_2d_num ( )

    do i = 1, num

      call function_2d_set ( 'SET', i )

      do rule = ilo, ihi

        order = triangle_unit_size ( rule )

        allocate ( xtab(1:order) )
        allocate ( ytab(1:order) )
        allocate ( weight(1:order) )

        call triangle_unit_set ( rule, order, xtab, ytab, weight )
 
        call triangle_unit_sum ( function_2d, order, xtab, ytab, weight, &
          result(rule) )

        deallocate ( xtab )
        deallocate ( ytab )
        deallocate ( weight )

      end do

      write ( *, '(2x,a7,5f14.8)' ) function_2d_name(i), result(ilo:ihi)

    end do

  end do

  return
end
subroutine test425 ( )

!*****************************************************************************80
!
!! TEST425 tests TRIANGLE_UNIT_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 20008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) coef
  real ( kind = 8 ) err
  real ( kind = 8 ) exact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  real ( kind = 8 ) quad
  integer ( kind = 4 ) rule
  integer ( kind = 4 ), parameter :: rule_max = 20
  integer ( kind = 4 ) triangle_unit_size
  real ( kind = 8 ) value
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight
  real ( kind = 8 ), allocatable, dimension ( : ) :: xtab
  real ( kind = 8 ), allocatable, dimension ( : ) :: ytab

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST425'
  write ( *, '(a)' ) '  TRIANGLE_UNIT_SET sets up a quadrature '
  write ( *, '(a)' ) '  in the unit triangle,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimate integral of X^A * Y^B.'

  do a = 0, 10

    do b = 0, 10 - a

      coef = real ( a + b + 2, kind = 8 ) * real ( a + b + 1, kind = 8 )
      do i = 1, b
        coef = coef * real ( a + i, kind = 8 ) / real ( i, kind = 8 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a,i8)' ) '  A = ', a, '  B = ', b
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Rule       QUAD           ERROR'
      write ( *, '(a)' ) ' '

      do rule = 1, rule_max

        order = triangle_unit_size ( rule )

        allocate ( xtab(1:order) )
        allocate ( ytab(1:order) )
        allocate ( weight(1:order) )

        call triangle_unit_set ( rule, order, xtab, ytab, weight )
 
        quad = 0.0D+00

        do i = 1, order

          if ( a == 0 .and. b == 0 ) then
            value = coef
          else if ( a == 0 .and. b /= 0 ) then
            value = coef              * ytab(i)**b
          else if ( a /= 0 .and. b == 0 ) then
            value = coef * xtab(i)**a 
          else if ( a /= 0 .and. b /= 0 ) then
            value = coef * xtab(i)**a * ytab(i)**b
          end if

          quad = quad + 0.5D+00 * weight(i) * value

        end do

        exact = 1.0D+00
        err = abs ( exact - quad )

        write ( *, '(2x,i4,2x,g14.6,2x,f14.8)' ) rule, quad, err

        deallocate ( xtab )
        deallocate ( ytab )
        deallocate ( weight )

      end do

    end do

  end do

  return
end
subroutine test43 ( )

!*****************************************************************************80
!
!! TEST43 tests TRIANGLE_UNIT_PRODUCT_SET and TRIANGLE_UNIT_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: rule_max = 8

  real ( kind = 8 ), external :: function_2d
  character ( len = 7 ) function_2d_name
  integer ( kind = 4 ) function_2d_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) num
  integer ( kind = 4 ) order
  real ( kind = 8 ) result(rule_max)
  integer ( kind = 4 ) rule
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight
  real ( kind = 8 ), allocatable, dimension ( : ) :: xtab
  real ( kind = 8 ), allocatable, dimension ( : ) :: ytab

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST43'
  write ( *, '(a)' ) '  TRIANGLE_UNIT_PRODUCT_SET sets up a product quadrature'
  write ( *, '(a)' ) '  rule in the unit triangle,'
  write ( *, '(a)' ) '  TRIANGLE_UNIT_SUM applies it.'
  write ( *, '(a)' ) ' '

  do ilo = 1, rule_max, 5

    ihi = min ( ilo +  4, rule_max )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a,5i6)' ) 'Rule Order: ', ( rule, rule = ilo, ihi )
    write ( *, '(2x,a)' ) 'Function '
    write ( *, '(a)' ) ' '

    num = function_2d_num ( )

    do i = 1, num

      call function_2d_set ( 'SET', i )

      do rule = ilo, ihi

        call triangle_unit_product_size ( rule, order )

        allocate ( xtab(1:order) )
        allocate ( ytab(1:order) )
        allocate ( weight(1:order) )

        call triangle_unit_product_set ( rule, order, xtab, ytab, weight )

        call triangle_unit_sum ( function_2d, order, xtab, ytab, weight, &
          result(rule) )

        deallocate ( xtab )
        deallocate ( ytab )
        deallocate ( weight )

      end do

      write ( *, '(2x,a7,5f14.8)' ) function_2d_name(i), result(ilo:ihi)

    end do

  end do

  return
end
subroutine test44 ( )

!*****************************************************************************80
!
!! TEST44 tests TRIANGLE_UNIT_SET and TRIANGLE_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: rule_max = 20

  real ( kind = 8 ), external :: function_2d
  character ( len = 7 ) function_2d_name
  integer ( kind = 4 ) function_2d_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) num
  integer ( kind = 4 ) order
  real ( kind = 8 ) result(rule_max)
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) triangle_unit_size
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight
  real ( kind = 8 ), allocatable, dimension ( : ) :: xtab
  real ( kind = 8 ), dimension ( 3 ) :: xval = (/ &
    1.0D+00, 3.0D+00, 1.0D+00 /)
  real ( kind = 8 ), allocatable, dimension ( : ) :: ytab
  real ( kind = 8 ), dimension ( 3 ) :: yval = (/ &
    1.0D+00, 1.0D+00, 4.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST44'
  write ( *, '(a)' ) '  TRIANGLE_UNIT_SET sets up quadrature '
  write ( *, '(a)' ) '  in the unit triangle,'
  write ( *, '(a)' ) '  TRIANGLE_SUM applies it to an arbitrary triangle.'
  write ( *, '(a)' ) ' '

  do ilo = 1, rule_max, 5

    ihi = min ( ilo + 4, rule_max )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a,5i6)' ) 'Rule:    ', ( rule, rule = ilo, ihi )
    write ( *, '(2x,a)' ) 'Function '
    write ( *, '(a)' ) ' '

    num = function_2d_num ( )

    do i = 1, num

      call function_2d_set ( 'SET', i )

      do rule = ilo, ihi

        order = triangle_unit_size ( rule )

        allocate ( xtab(1:order) )
        allocate ( ytab(1:order) )
        allocate ( weight(1:order) )

        call triangle_unit_set ( rule, order, xtab, ytab, weight )
 
        call triangle_sum ( function_2d, xval, yval, order, xtab, ytab, weight, &
          result(rule) )

        deallocate ( xtab )
        deallocate ( ytab )
        deallocate ( weight )

      end do

      write ( *, '(2x,a7,5f14.8)' ) function_2d_name(i), result(ilo:ihi)

    end do

  end do

  return
end
subroutine test45 ( )

!*****************************************************************************80
!
!! TEST45 tests TORUS_1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 7 ) function_3d_name
  integer ( kind = 4 ) function_3d_num
  real ( kind = 8 ), external :: function_3d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) num
  real ( kind = 8 ) result(5)
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) torus_area_3d

  r1 = 0.5D+00
  r2 = 1.0D+00
  n = 10
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST45'
  write ( *, '(a)' ) '  TORUS_1 approximates integrals on a torus.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The degree N will be varied.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Inner radius = ', r1
  write ( *, '(a,g14.6)' ) '  Outer radius = ', r2
  write ( *, '(a,g14.6)' ) '  Area = ', torus_area_3d ( r1, r2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,5i5)' ) '    F(X)  ', ( 2**j, j = 0, 8, 2 )
  write ( *, '(a)' ) ' '
 
  num = function_3d_num ( )

  do i = 1, num

    call function_3d_set ( 'SET', i )

    do j = 1, 5

      j2 = 2 * ( j - 1 )
      n = 2**j2
      call torus_1 ( function_3d, r1, r2, n, result(j) )

    end do

    write ( *, '(2x,a7,5f14.8)' ) function_3d_name(i), result(1:5)

  end do
 
  return
end
subroutine test46 ( )

!*****************************************************************************80
!
!! TEST46 tests TORUS_5S2, TORUS_6S2 and TORUS_14S.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 7 ) function_3d_name
  integer ( kind = 4 ) function_3d_num
  real ( kind = 8 ), external :: function_3d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) num
  real ( kind = 8 ) result1
  real ( kind = 8 ) result2
  real ( kind = 8 ) result3
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) torus_volume_3d

  r1 = 0.5D+00
  r2 = 1.0D+00
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST46'
  write ( *, '(a)' ) '  For the interior of a torus,'
  write ( *, '(a)' ) '  TORUS_5S2,'
  write ( *, '(a)' ) '  TORUS_6S2, and'
  write ( *, '(a)' ) '  TORUS_5S2 approximate integrals.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Inner radius = ', r1
  write ( *, '(a,g14.6)' ) '  Outer radius = ', r2
  write ( *, '(a,g14.6)' ) '  Volume = ', torus_volume_3d ( r1, r2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Rule:        #5S2          #6S2          #14S'
  write ( *, '(a)' ) '    F(X)'
  write ( *, '(a)' ) ' '
 
  num = function_3d_num ( )

  do i = 1, num

    call function_3d_set ( 'SET', i )

    call torus_5s2 ( function_3d, r1, r2, result1 )
    call torus_6s2 ( function_3d, r1, r2, result2 )
    call torus_14s ( function_3d, r1, r2, result3 )

    write ( *, '(2x,a7,3f14.8)' ) function_3d_name(i), result1, result2, result3

  end do
 
  return
end
subroutine test47 ( )

!*****************************************************************************80
!
!! TEST47 tests TORUS_SQUARE_5C2 and TORUS_SQUARE_14C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 7 ) function_3d_name
  integer ( kind = 4 ) function_3d_num
  real ( kind = 8 ), external :: function_3d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) num
  real ( kind = 8 ) result1
  real ( kind = 8 ) result2
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) torus_square_volume_3d

  r1 = 1.0D+00
  r2 = 0.125D+00
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST47'
  write ( *, '(a)' ) '  For integrals inside a torus with square cross-section:'
  write ( *, '(a)' ) '  TORUS_SQUARE_5C2 approximates the integral;'
  write ( *, '(a)' ) '  TORUS_SQUARE_14C approximates the integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Inner radius = ', r1
  write ( *, '(a,g14.6)' ) '  Outer radius = ', r2
  write ( *, '(a,g14.6)' ) '  Volume = ', torus_square_volume_3d ( r1, r2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    F(X)    5C2           14C'
  write ( *, '(a)' ) ' '
 
  num = function_3d_num ( )

  do i = 1, num

    call function_3d_set ( 'SET', i )

    call torus_square_5c2 ( function_3d, r1, r2, result1 )
    call torus_square_14c ( function_3d, r1, r2, result2 )

    write ( *, '(2x,a7,2f14.8)' ) function_3d_name(i), result1, result2

  end do
 
  return
end
subroutine test48 ( )

!*****************************************************************************80
!
!! TEST48 tests TVEC_EVEN, TVEC_EVEN2 and TVEC_EVEN3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) nt
  real ( kind = 8 ), allocatable, dimension ( : ) :: t

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST48'
  write ( *, '(a)' ) '  For evenly spaced angles between 0 and 2*PI:'
  write ( *, '(a)' ) '  TVEC_EVEN'
  write ( *, '(a)' ) '  TVEC_EVEN2'
  write ( *, '(a)' ) '  TVEC_EVEN3'

  nt = 4
  allocate ( t(1:nt) )
  call tvec_even ( nt, t )
  call r8vec_print ( nt, t, '  TVEC_EVEN:' )
  deallocate ( t )

  nt = 4
  allocate ( t(1:nt) )
  call tvec_even2 ( nt, t )
  call r8vec_print ( nt, t, '  TVEC_EVEN2:' )
  deallocate ( t )

  nt = 4
  allocate ( t(1:nt) )
  call tvec_even3 ( nt, t )
  call r8vec_print ( nt, t, '  TVEC_EVEN3:' )
  deallocate ( t )

  return
end
subroutine test49 ( )

!*****************************************************************************80
!
!! TEST49 tests TVEC_EVEN_BRACKET, TVEC_EVEN_BRACKET2 and TVEC_EVEN_BRACKET3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) nt
  real ( kind = 8 ), allocatable, dimension ( : ) :: t
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST49'
  write ( *, '(a)' ) '  For evenly spaced angles between THETA1 and THETA2:'
  write ( *, '(a)' ) '  TVEC_EVEN_BRACKET'
  write ( *, '(a)' ) '  TVEC_EVEN_BRACKET2.'
  write ( *, '(a)' ) '  TVEC_EVEN_BRACKET3.'

  theta1 = 30.0D+00
  theta2 = 90.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '    THETA1 = ', theta1
  write ( *, '(a,g14.6)' ) '    THETA2 = ', theta2

  nt = 4
  allocate ( t(1:nt) )
  call tvec_even_bracket ( nt, theta1, theta2, t )
  call r8vec_print ( nt, t, '  TVEC_EVEN_BRACKET' )
  deallocate ( t )

  nt = 5
  allocate ( t(1:nt) )
  call tvec_even_bracket2 ( nt, theta1, theta2, t )
  call r8vec_print ( nt, t, '  TVEC_EVEN_BRACKET2' )
  deallocate ( t )

  nt = 3
  allocate ( t(1:nt) )
  call tvec_even_bracket3 ( nt, theta1, theta2, t )
  call r8vec_print ( nt, t, '  TVEC_EVEN_BRACKET3' )
  deallocate ( t )

  return
end
function function_1d ( x )

!*****************************************************************************80
!
!! FUNCTION_1D evaluates a function F(X) of one variable.
!
!  Discussion:
!
!    The actual form of the function can be determined by calling FUNCSET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value of the variable.
!
!    Output, real ( kind = 8 ) FUNCTION_1D, the value of the function.
!
  implicit none

  real ( kind = 8 ) function_1d
  integer ( kind = 4 ) ifunc
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  call function_1d_set ( 'GET', ifunc )

  if ( ifunc == 1 ) then
    value = 1.0D+00
  else if ( ifunc == 2 ) then
    value = x
  else if ( ifunc == 3 ) then
    value = x**2
  else if ( ifunc == 4 ) then
    value = x**3
  else if ( ifunc == 5 ) then
    value = x**4
  else if ( ifunc == 6 ) then
    value = x**5
  else if ( ifunc == 7 ) then
    value = x**6
  else if ( ifunc == 8 ) then
    value = abs ( x )
  else if ( ifunc == 9 ) then
    value = sin ( x )
  else if ( ifunc == 10 ) then
    value = exp ( x )
  else if ( ifunc == 11 ) then
    value = 1.0D+00 / ( 1.0D+00 + abs ( x ) )
  else if ( ifunc == 12 ) then
    value = sqrt ( abs ( x ) )
  else
    value = 0.0D+00
  end if

  function_1d = value

  return
end
function function_1d_name ( ifunc )

!*****************************************************************************80
!
!! FUNCTION_1D_NAME returns the name of the current 1D function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IFUNC, the index of the function.
!
!    Output, character ( len = 7 ) FUNCTION_1D_NAME, the name of the 1D function.
!
  implicit none

  character ( len = 7 ) function_1d_name
  integer ( kind = 4 ) ifunc

  if ( ifunc == 1 ) then
    function_1d_name = '      1'
  else if ( ifunc == 2 ) then
    function_1d_name = '      X'
  else if ( ifunc == 3 ) then
    function_1d_name = '    X^2'
  else if ( ifunc == 4 ) then
    function_1d_name = '    X^3'
  else if ( ifunc == 5 ) then
    function_1d_name = '    X^4'
  else if ( ifunc == 6 ) then
    function_1d_name = '    X^5'
  else if ( ifunc == 7 ) then
    function_1d_name = '    X^6'
  else if ( ifunc == 8 ) then
    function_1d_name = '      R'
  else if ( ifunc == 9 ) then
    function_1d_name = ' SIN(X)'
  else if ( ifunc == 10 ) then
    function_1d_name = ' EXP(X)'
  else if ( ifunc == 11 ) then
    function_1d_name = '1/(1+R)'
  else if ( ifunc == 12 ) then
    function_1d_name = 'SQRT(R)'
  else
    function_1d_name = '???????'
  end if

  return
end
function function_1d_num ( )

!*****************************************************************************80
!
!! FUNCTION_1D_NUM returns the number of 1D functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FUNCTION_1D_NUM, the number of 1D functions.
!
  implicit none

  integer ( kind = 4 ) function_1d_num

  function_1d_num = 12

  return
end
subroutine function_1d_set ( action, i )

!*****************************************************************************80
!
!! FUNCTION_1D_SET sets or reports the index of the current 1D function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION:
!    'GET', please return the index of the current function.
!    'SET', please set the function according to the input index.
!
!    Input/output, integer ( kind = 4 ) I.
!    If ACTION = 'SET', then I is input, and is the index of the desired
!    function.
!    If ACTION = 'GET', then I is output, and is the index of the current
!    function.
!
  implicit none

  character ( len = * ) action
  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: ival = 0

  if ( action == 'SET' ) then
    ival = i
  else if ( action == 'GET' ) then
    i = ival
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FUNCTION_1D_SET - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop
  end if

  return
end
function function_2d ( x, y )

!*****************************************************************************80
!
!! FUNCTION_2D evaluates a function F(X,Y) of two variables.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the value of the variables.
!
!    Output, real ( kind = 8 ) FUNCTION_2D, the value of the function.
!
  implicit none

  real ( kind = 8 ) function_2d
  integer ( kind = 4 ) ifunc
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  call function_2d_set ( 'GET', ifunc )

  if ( ifunc == 1 ) then
    value = 1.0D+00
  else if ( ifunc == 2 ) then
    value = x
  else if ( ifunc == 3 ) then
    value = x * x
  else if ( ifunc == 4 ) then
    value = x * x * x
  else if ( ifunc == 5 ) then
    value = x * x * x * x
  else if ( ifunc == 6 ) then
    value = x**5
  else if ( ifunc == 7 ) then
    value = x**6
  else if ( ifunc == 8 ) then
    value = sqrt ( x * x + y * y )
  else if ( ifunc == 9 ) then
    value = sin ( x )
  else if ( ifunc == 10 ) then
    value = exp ( x )
  else if ( ifunc == 11 ) then
    value = 1.0D+00 / ( 1.0D+00 + sqrt ( x * x + y * y ) )
  else if ( ifunc == 12 ) then
    value = sqrt ( sqrt ( x * x + y * y ) )
  else
    value = 0.0D+00
  end if

  function_2d = value

  return
end
function function_2d_name ( ifunc )

!*****************************************************************************80
!
!! FUNCTION_2D_NAME returns the name of the current 2D function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IFUNC, the index of the function.
!
!    Output, character ( len = 7 ) FUNCTION_2D_NAME, the name of the 2D function.
!
  implicit none

  character ( len = 7 ) function_2d_name
  integer ( kind = 4 ) ifunc

  if ( ifunc == 1 ) then
    function_2d_name = '      1'
  else if ( ifunc == 2 ) then
    function_2d_name = '      X'
  else if ( ifunc == 3 ) then
    function_2d_name = '    X^2'
  else if ( ifunc == 4 ) then
    function_2d_name = '    X^3'
  else if ( ifunc == 5 ) then
    function_2d_name = '    X^4'
  else if ( ifunc == 6 ) then
    function_2d_name = '    X^5'
  else if ( ifunc == 7 ) then
    function_2d_name = '    X^6'
  else if ( ifunc == 8 ) then
    function_2d_name = '      R'
  else if ( ifunc == 9 ) then
    function_2d_name = ' SIN(X)'
  else if ( ifunc == 10 ) then
    function_2d_name = ' EXP(X)'
  else if ( ifunc == 11 ) then
    function_2d_name = '1/(1+R)'
  else if ( ifunc == 12 ) then
    function_2d_name = 'SQRT(R)'
  else
    function_2d_name = '???????'
  end if

  return
end
function function_2d_num ( )

!*****************************************************************************80
!
!! FUNCTION_2D_NUM returns the number of 2D functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FUNCTION_2D_NUM, the number of 2D functions.
!
  implicit none

  integer ( kind = 4 ) function_2d_num

  function_2d_num = 12

  return
end
subroutine function_2d_set ( action, i )

!*****************************************************************************80
!
!! FUNCTION_2D_SET sets or reports the index of the current 2D function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION:
!    'GET', please return the index of the current function.
!    'SET', please set the function according to the input index.
!
!    Input/output, integer ( kind = 4 ) I.
!    If ACTION = 'SET', then I is input, and is the index of the desired
!    function.
!    If ACTION = 'GET', then I is output, and is the index of the current
!    function.
!
  implicit none

  character ( len = * ) action
  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: ival = 0

  if ( action == 'SET' ) then
    ival = i
  else if ( action == 'GET' ) then
    i = ival
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FUNCTION_2D_SET - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop
  end if

  return
end
function function_3d ( x, y, z )

!*****************************************************************************80
!
!! FUNCTION_3D evaluates a the current 3D function.
!
!  Discussion:
!
!    The actual form of the function can be determined by calling FUNCTION_3D_SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, Z, the value of the variables.
!
!    Output, real ( kind = 8 ) FUNCTION_3D, the value of the function.
!
  implicit none

  real ( kind = 8 ) function_3d
  integer ( kind = 4 ) ifunc
  real ( kind = 8 ) value
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  call function_3d_set ( 'GET', ifunc )

  if ( ifunc == 1 ) then
    value = 1.0D+00

  else if ( ifunc == 2 ) then
    value = x
  else if ( ifunc == 3 ) then
    value = y
  else if ( ifunc == 4 ) then
    value = z

  else if ( ifunc == 5 ) then
    value = x * x
  else if ( ifunc == 6 ) then
    value = x * y
  else if ( ifunc == 7 ) then
    value = x * z
  else if ( ifunc == 8 ) then
    value = y * y
  else if ( ifunc == 9 ) then
    value = y * z
  else if ( ifunc == 10 ) then
    value = z * z

  else if ( ifunc == 11 ) then
    value = x * x * x
  else if ( ifunc == 12 ) then
    value = x * y * z
  else if ( ifunc == 13 ) then
    value = z * z * z

  else if ( ifunc == 14 ) then
    value = x * x * x * x
  else if ( ifunc == 15 ) then
    value = x * x * z * z
  else if ( ifunc == 16 ) then
    value = z * z * z * z

  else if ( ifunc == 17 ) then
    value = x * x * x * x * x
  else if ( ifunc == 18 ) then
    value = x**6
  else if ( ifunc == 19 ) then
    value = sqrt ( x * x + y * y + z * z )
  else if ( ifunc == 20 ) then
    value = sin ( x )
  else if ( ifunc == 21 ) then
    value = exp ( x )
  else if ( ifunc == 22 ) then
    value = 1.0D+00 / sqrt ( 1.0D+00 + x * x + y * y + z * z )
  else if ( ifunc == 23 ) then
    value = sqrt ( sqrt ( x * x + y * y + z * z ) )
  else
    value = 0.0D+00
  end if

  function_3d = value

  return
end
function function_3d_name ( ifunc )

!*****************************************************************************80
!
!! FUNCTION_3D_NAME returns the name of the current 3D function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IFUNC, the index of the function.
!
!    Output, character ( len = 7 ) FUNCTION_3D_NAME, the name of the function.
!
  implicit none

  character ( len = 7 ) function_3d_name
  integer ( kind = 4 ) ifunc

  if ( ifunc == 1 ) then
    function_3d_name = '      1'

  else if ( ifunc == 2 ) then
    function_3d_name = '      X'
  else if ( ifunc == 3 ) then
    function_3d_name = '      Y'
  else if ( ifunc == 4 ) then
    function_3d_name = '      Z'

  else if ( ifunc == 5 ) then
    function_3d_name = '    X*X'
  else if ( ifunc == 6 ) then
    function_3d_name = '    X*Y'
  else if ( ifunc == 7 ) then
    function_3d_name = '    X*Z'
  else if ( ifunc == 8 ) then
    function_3d_name = '    Y*Y'
  else if ( ifunc == 9 ) then
    function_3d_name = '    Y*Z'
  else if ( ifunc == 10 ) then
    function_3d_name = '    Z*Z'

  else if ( ifunc == 11 ) then
    function_3d_name = '    X^3'
  else if ( ifunc == 12 ) then
    function_3d_name = '  X*Y*Z'
  else if ( ifunc == 13 ) then
    function_3d_name = '  Z*Z*Z'

  else if ( ifunc == 14 ) then
    function_3d_name = '    X^4'
  else if ( ifunc == 15 ) then
    function_3d_name = 'X^2 Z^2'
  else if ( ifunc == 16 ) then
    function_3d_name = '    Z^4'

  else if ( ifunc == 17 ) then
    function_3d_name = '    X^5'
  else if ( ifunc == 18 ) then
    function_3d_name = '    X^6'
  else if ( ifunc == 19 ) then
    function_3d_name = '      R'
  else if ( ifunc == 20 ) then
    function_3d_name = ' SIN(X)'
  else if ( ifunc == 21 ) then
    function_3d_name = ' EXP(X)'
  else if ( ifunc == 22 ) then
    function_3d_name = '1/(1+R)'
  else if ( ifunc == 23 ) then
    function_3d_name = 'SQRT(R)'
  else
    function_3d_name = '???????'
  end if

  return
end
function function_3d_num ( )

!*****************************************************************************80
!
!! FUNCTION_3D_NUM returns the number of 3D functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FUNCTION_3D_NUM, the number of 3D functions.
!
  implicit none

  integer ( kind = 4 ) function_3d_num

  function_3d_num = 23

  return
end
subroutine function_3d_set ( action, i )

!*****************************************************************************80
!
!! FUNCTION_3D_SET sets or reports the index of the current 3D function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION:
!    'GET', please return the index of the current function.
!    'SET', please set the function according to the input index.
!
!    Input/output, integer ( kind = 4 ) I.
!    If ACTION = 'SET', then I is input, and is the index of the desired
!    function.
!    If ACTION = 'GET', then I is output, and is the index of the current
!    function.
!
  implicit none

  character ( len = * ) action
  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: ival = 0

  if ( action == 'SET' ) then
    ival = i
  else if ( action == 'GET' ) then
    i = ival
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FUNCTION_3D_SET - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop
  end if

  return
end
function function_nd ( n, x )

!*****************************************************************************80
!
!! FUNCTION_ND evaluates a function of N variables.
!
!  Discussion:
!
!    The actual form of the function can be determined by calling FUNCSET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) X(N), the value of the variables.
!
!    Output, real ( kind = 8 ) FUNCTION_ND, the value of the function.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) function_nd
  integer ( kind = 4 ) ifunc
  real ( kind = 8 ) value
  real ( kind = 8 ) x(n)

  call function_nd_set ( 'GET', ifunc )

  if ( ifunc == 1 ) then
    value = 1.0D+00
  else if ( ifunc == 2 ) then
    value = x(1)
  else if ( ifunc == 3 ) then
    value = x(1)**2
  else if ( ifunc == 4 ) then
    value = x(1)**3
  else if ( ifunc == 5 ) then
    value = x(1)**4
  else if ( ifunc == 6 ) then
    value = x(1)**5
  else if ( ifunc == 7 ) then
    value = x(1)**6
  else if ( ifunc == 8 ) then
    value = sqrt ( sum ( x(1:n)**2 ) )
  else if ( ifunc == 9 ) then
    value = sin ( x(1) )
  else if ( ifunc == 10 ) then
    value = exp ( x(1) )
  else if ( ifunc == 11 ) then
    value = 1.0D+00 / ( 1.0D+00 + sqrt ( sum ( x(1:n)**2 ) ) )
  else if ( ifunc == 12 ) then
    value = sqrt ( sqrt ( sum ( x(1:n)**2 ) ) )
  else
    value = 0.0D+00
  end if

  function_nd = value

  return
end
function function_nd_name ( ifunc )

!*****************************************************************************80
!
!! FUNCTION_ND_NAME returns the name of the current ND function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IFUNC, the index of the function.
!
!    Output, character ( len = 7 ) FUNCTION_ND_NAME, the name of the function.
!
  implicit none

  character ( len = 7 ) function_nd_name
  integer ( kind = 4 ) ifunc

  if ( ifunc == 1 ) then
    function_nd_name = '      1'
  else if ( ifunc == 2 ) then
    function_nd_name = '      X'
  else if ( ifunc == 3 ) then
    function_nd_name = '    X^2'
  else if ( ifunc == 4 ) then
    function_nd_name = '    X^3'
  else if ( ifunc == 5 ) then
    function_nd_name = '    X^4'
  else if ( ifunc == 6 ) then
    function_nd_name = '    X^5'
  else if ( ifunc == 7 ) then
    function_nd_name = '    X^6'
  else if ( ifunc == 8 ) then
    function_nd_name = '      R'
  else if ( ifunc == 9 ) then
    function_nd_name = ' SIN(X)'
  else if ( ifunc == 10 ) then
    function_nd_name = ' EXP(X)'
  else if ( ifunc == 11 ) then
    function_nd_name = '1/(1+R)'
  else if ( ifunc == 12 ) then
    function_nd_name = 'SQRT(R)'
  else
    function_nd_name = '???????'
  end if

  return
end
function function_nd_num ( )

!*****************************************************************************80
!
!! FUNCTION_ND_NUM returns the number of ND functions.
!
!  Modified:
!
!    06 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FUNCTION_ND_NUM, the number of ND functions.
!
  implicit none

  integer ( kind = 4 ) function_nd_num

  function_nd_num = 12

  return
end
subroutine function_nd_set ( action, i )

!*****************************************************************************80
!
!! FUNCTION_ND_SET sets or reports the index of the current ND function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION:
!    'GET', please return the index of the current function.
!    'SET', please set the function according to the input index.
!
!    Input/output, integer ( kind = 4 ) I.
!    If ACTION = 'SET', then I is input, and is the index of the desired
!    function.
!    If ACTION = 'GET', then I is output, and is the index of the current
!    function.
!
  implicit none

  character ( len = * ) action
  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: ival = 0

  if ( action == 'SET' ) then
    ival = i
  else if ( action == 'GET' ) then
    i = ival
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FUNCSET - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized action = "' // trim ( action ) // '".'
    stop
  end if

  return
end
function fup7 ( x )

!*****************************************************************************80
!
!! FUP7 is the upper limit function for QMULT_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fup7
  real ( kind = 8 ) x

  fup7 = 1.0D+00

  return
end
function flo7 ( x )

!*****************************************************************************80
!
!! FLO7 is the lower limit function for QMULT_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) flo7
  real ( kind = 8 ) x

  flo7 = -1.0D+00

  return
end
function fu18 ( x )

!*****************************************************************************80
!
!! FU18 is the upper limit of integration for x.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fu18
  real ( kind = 8 ) x

  fu18 = 1.0D+00

  return
end
function fl18 ( x )

!*****************************************************************************80
!
!! FL18 is the lower limit of integration for x.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fl18
  real ( kind = 8 ) x

  fl18 = -1.0D+00

  return
end
function fu28 ( x, y )

!*****************************************************************************80
!
!! FU28 computes the upper limit of integration for (x,y).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fu28
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  fu28 = 1.0D+00

  return
end
function fl28 ( x, y )

!*****************************************************************************80
!
!! FL28 computes the lower limit of integration for (x,y).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fl28
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  fl28 = -1.0D+00

  return
end
function mono_000_3d ( n, x )

!*****************************************************************************80
!
!! MONO_000_3D evaluates X**0 Y**0 Z**0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) mono_000_3d
  real ( kind = 8 ) x(n)

  mono_000_3d = 1.0D+00

  return
end
function mono_111_3d ( n, x )

!*****************************************************************************80
!
!! MONO_111_3D evaluates X^1 Y^1 Z^1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) mono_111_3d
  real ( kind = 8 ) x(n)

  mono_111_3d = x(1) * x(2) * x(3)

  return
end
function mono_202_3d ( n, x )

!*****************************************************************************80
!
!! MONO_202_3D evaluates X^2 Y^0 Z^2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) mono_202_3d
  real ( kind = 8 ) x(n)

  mono_202_3d = x(1)**2 * x(3)**2

  return
end
function mono_422_3d ( n, x )

!*****************************************************************************80
!
!! MONO_422_3D evaluates X^4 Y^2 Z^2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) mono_422_3d
  real ( kind = 8 ) x(n)

  mono_422_3d = x(1)**4 * x(2)**2 * x(3)**2

  return
end
function f_1_2d ( x, y )

!*****************************************************************************80
!
!! F_1_2D evaluates the function 1 in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the arguments.
!
!    Output, real ( kind = 8 ) F_1_2D, the value of the function.
!
  implicit none

  real ( kind = 8 ) f_1_2d
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  f_1_2d = 1.0D+00

  return
end
function f_x_2d ( x, y )

!*****************************************************************************80
!
!! F_X_2D evaluates the function X in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the arguments.
!
!    Output, real ( kind = 8 ) F_X_2D, the value of the function.
!
  implicit none

  real ( kind = 8 ) f_x_2d
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  f_x_2d = x

  return
end
function f_r_2d ( x, y )

!*****************************************************************************80
!
!! F_R_2D evaluates the function sqrt ( X^2 + Y^2 ) in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the arguments.
!
!    Output, real ( kind = 8 ) F_R_2D, the value of the function.
!
  implicit none

  real ( kind = 8 ) f_r_2d
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  f_r_2d = sqrt ( x * x + y * y )

  return
end
