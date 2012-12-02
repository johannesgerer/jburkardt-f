program main

!*****************************************************************************80
!
!! MAIN is the main program for GEOMETRY_PRB.
!
!  Discussion:
!
!    GEOMETRY_PRB tests routines from the GEOMETRY library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 November 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEOMETRY_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the GEOMETRY library.'

  call test000 ( )
  call test0005 ( )
  call test001 ( )
  call test002 ( )
  call test0023 ( )
  call test0025 ( )
  call test003 ( )
  call test0032 ( )
  call test0035 ( )
  call test004 ( )
  call test0045 ( )
  call test005 ( )
  call test006 ( )
  call test007 ( )
  call test0075 ( )
  call test008 ( )
  call test0085 ( )
  call test009 ( )

  call test010 ( )
  call test011 ( )
  call test012 ( )
  call test0125 ( )
  call test0126 ( )
  call test0127 ( )
  call test013 ( )
  call test014 ( )
  call test015 ( )
  call test0155 ( )
  call test0156 ( )
  call test016 ( )
  call test0165 ( )
  call test017 ( )
  call test018 ( )
  call test0183 ( )
  call test0185 ( )
  call test019 ( )

  call test020 ( )
  call test0201 ( )
  call test02015 ( )
  call test0202 ( )
  call test0203 ( )
  call test02035 ( )
  call test0204 ( )
  call test0205 ( )
  call test021 ( )
  call test023 ( )
  call test0232 ( )
  call test0234 ( )
  call test0235 ( )
  call test0236 ( )
  call test0238 ( )
  call test024 ( )
  call test0243 ( )
  call test0245 ( )
  call test025 ( )
  call test0255 ( )
  call test026 ( )
  call test027 ( )
  call test028 ( )
  call test029 ( )

  call test030 ( )
  call test031 ( )
  call test0315 ( )
  call test032 ( )
  call test0321 ( )
  call test0322 ( )
  call test0323 ( )
  call test0325 ( )
  call test0327 ( )
  call test033 ( )
  call test0335 ( )
  call test0336 ( )
  call test0337 ( )
  call test034 ( )
  call test0345 ( )
  call test0346 ( )
  call test035 ( )
  call test038 ( )
  call test0385 ( )
  call test03855 ( )
  call test0386 ( )
  call test039 ( )

  call test040 ( )
  call test041 ( )
  call test0415 ( )
  call test0416 ( )
  call test0418 ( )
  call test042 ( )
  call test043 ( )
  call test044 ( )
  call test045 ( )
  call test046 ( )
  call test047 ( )
  call test0475 ( )
  call test0477 ( )
  call test0478 ( )
  call test048 ( )
  call test0485 ( )
  call test049 ( )
  call test0493 ( )
  call test0495 ( )

  call test050 ( )
  call test051 ( )
  call test052 ( )
  call test053 ( )
  call test054 ( )
  call test055 ( )
  call test056 ( )
  call test057 ( )
  call test058 ( )
  call test059 ( )

  call test060 ( )
  call test061 ( )
  call test0615 ( )
  call test0616 ( )
  call test0617 ( )
  call test062 ( )
  call test063 ( )
  call test064 ( )
  call test065 ( )
  call test066 ( )
  call test067 ( )
  call test068 ( )
  call test0685 ( )

  call test0755 ( )
  call test0757 ( )
  call test076 ( )
  call test0765 ( )
  call test078 ( )
  call test0782 ( )
  call test0784 ( )
  call test0786 ( )
  call test079 ( )

  call test080 ( )
  call test0801 ( )
  call test0803 ( )
  call test0805 ( )
  call test0807 ( )
  call test081 ( )
  call test082 ( )
  call test0825 ( )
  call test083 ( )
  call test084 ( )
  call test0844 ( )
  call test0845 ( )
  call test0846 ( )
  call test085 ( )

  call test170 ( )
  call test171 ( )
  call test1712 ( )
  call test1715 ( )
  call test172 ( )
  call test173 ( )
  call test174 ( )
  call test1745 ( )
  call test1746 ( )
  call test175 ( )
  call test176 ( )
  call test177 ( )
  call test178 ( )
  call test1787 ( )
  call test1893 ( )
  call test036 ( )
  call test0364 ( )
  call test0365 ( )
  call test0366 ( )
  call test0367 ( )
  call test0368 ( )
  call test037 ( )
  call test1788 ( )
  call test1789 ( )
  call test179 ( )

  call test180 ( )
  call test1804 ( )
  call test1805 ( )
  call test181 ( )
  call test182 ( )
  call test183 ( )
  call test1835 ( )
  call test1836 ( )
  call test188 ( )
  call test189 ( )
  call test1895 ( )

  call test190 ( )
  call test191 ( )
  call test192 ( )
  call test193 ( )
  call test194 ( )
  call test195 ( )
  call test1955 ( )
  call test196 ( )
  call test197 ( )
  call test198 ( )
  call test199 ( )

  call test200 ( )
  call test201 ( )
  call test202 ( )
  call test203 ( )
  call test2031 ( )
  call test2032 ( )
  call test20321 ( )
  call test20322 ( )
  call test203224 ( )
  call test203225 ( )
  call test20323 ( )
  call test203232 ( )
  call test203233 ( )
  call test203234 ( )
  call test203235 ( )
  call test20324 ( )
  call test20325 ( )
  call test2033 ( )
  call test204 ( )
  call test205 ( )
  call test206 ( )
  call test20605 ( )
  call test2061 ( )
  call test2062 ( )
  call test209 ( )
  call test20655 ( )
  call test2066 ( )
  call test2094 ( )
  call test2101 ( )
  call test21011 ( )
  call test2067 ( )
  call test21015 ( )
  call test2068 ( )
  call test2069 ( )
  call test207 ( )
  call test2075 ( )
  call test208 ( )

  call test2102 ( )
  call test2070 ( )
  call test20701 ( )
  call test2104 ( )
  call test2105 ( )
  call test211 ( )
  call test2103 ( )
  call test2071 ( )
  call test20715 ( )
  call test2095 ( )
  call test2072 ( )
  call test2115 ( )
  call test212 ( )
  call test213 ( )
  call test219 ( )

  call test220 ( )
  call test221 ( )
  call test222 ( )
  call test2225 ( )
  call test223 ( )
  call test224 ( )
  call test2245 ( )
  call test225 ( )
  call test226 ( )
  call test227 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEOMETRY_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test000 ( )

!*****************************************************************************80
!
!! TEST000 tests RANDOM_SEED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST000'
  write ( *, '(a)' ) '  Call RANDOM_SEED to initialize the random number'
  write ( *, '(a)' ) '  generator.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Input SEED = ', seed

  call random_seed ( seed )

  return
end
subroutine test0005 ( )

!*****************************************************************************80
!
!! TEST0005 tests ANGLE_BOX_2D.
!
!  Discussion:
!
!    Test 1:
!
!      y = 0
!      y = 2x-6
!
!    Test 2:
!
!      y = 0
!      y = 2x-6
!
!    Test 3:
!
!      By setting P1 = P2, we are asking to extend the line
!      y = 2x-6
!      from P3 to P2 through to the other side.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) dist
  real ( kind = 8 ), dimension ( test_num ) :: dist_test = (/ &
    1.0D+00, 1.0D+00, 1.0D+00 /)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p1_test = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, &
    3.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p2_test = reshape ( (/ &
    3.0D+00, 0.0D+00, &
    3.0D+00, 0.0D+00, &
    3.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) p3(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p3_test = reshape ( (/ &
    4.0D+00,  2.0D+00, &
    2.0D+00, -2.0D+00, &
    2.0D+00, -2.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) p4(dim_num)
  real ( kind = 8 ) p5(dim_num)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0005'
  write ( *, '(a)' ) '  ANGLE_BOX_2D'
  write ( *, '(a)' ) '  Compute points P4 and P5, normal to '
  write ( *, '(a)' ) '  line through P1 and P2, and'
  write ( *, '(a)' ) '  line through P2 and P3, '
  write ( *, '(a)' ) '  and DIST units from P2.'

  do test = 1, test_num

    p1(1:dim_num) = p1_test(1:dim_num,test)
    p2(1:dim_num) = p2_test(1:dim_num,test)
    p3(1:dim_num) = p3_test(1:dim_num,test)
    dist = dist_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  DIST = ', dist
    write ( *, '(a,2g14.6)' ) '  P1:', p1(1:dim_num)
    write ( *, '(a,2g14.6)' ) '  P2:', p2(1:dim_num)
    write ( *, '(a,2g14.6)' ) '  P3:', p3(1:dim_num)
 
    call angle_box_2d ( dist, p1, p2, p3, p4, p5 )
 
    write ( *, '(a,2g14.6)' ) '  P4:', p4(1:dim_num)
    write ( *, '(a,2g14.6)' ) '  P5:', p5(1:dim_num)

  end do

  return
end
subroutine test001 ( )

!*****************************************************************************80
!
!! TEST001 tests ANGLE_CONTAINS_POINT_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 6

  integer ( kind = 4 ) angle
  integer ( kind = 4 ), parameter :: angle_num = 12
  real ( kind = 8 ) angle_rad
  logical inside
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p1_test = reshape ( (/ &
    1.0D+00,  0.0D+00, &
    1.0D+00,  0.0D+00, &
    1.0D+00, -1.0D+00, &
    1.0D+00,  0.0D+00, &
    1.0D+00,  0.0D+00, &
    1.0D+00,  0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p2_test = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) p3(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p3_test = reshape ( (/ &
     1.0D+00,  1.0D+00, &
     0.0D+00,  1.0D+00, &
     0.0D+00,  1.0D+00, &
    -1.0D+00,  0.0D+00, &
     0.0D+00, -1.0D+00, &
     1.0D+00, -0.01D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST001'
  write ( *, '(a)' ) '  ANGLE_CONTAINS_POINT_2D sees if a point'
  write ( *, '(a)' ) '  lies within an angle.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p1(1:dim_num) = p1_test(1:dim_num,test)
    p2(1:dim_num) = p2_test(1:dim_num,test)
    p3(1:dim_num) = p3_test(1:dim_num,test)

    call r8vec_print ( dim_num, p1, '  Vertex P1' )
    call r8vec_print ( dim_num, p2, '  Vertex P2' )
    call r8vec_print ( dim_num, p3, '  Vertex P3' )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '       X            Y       Inside?'
    write ( *, '(a)' ) ' '

    do angle = 0, angle_num

      angle_rad = real ( angle, kind = 8 ) * 2.0D+00 * pi &
        / real ( angle_num, kind = 8 )

      p(1:2) = (/ cos ( angle_rad ), sin ( angle_rad ) /)

      call angle_contains_point_2d ( p1, p2, p3, p, inside )

      write ( *, '(2x,2g14.6,2x,l1)' ) p(1:2), inside

    end do

  end do
 
  return
end
subroutine test002 ( )

!*****************************************************************************80
!
!! TEST002 tests ANGLE_DEG_2D and ANGLE_RAD_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) angle_deg_2d
  integer ( kind = 4 ), parameter :: angle_num = 12
  real ( kind = 8 ) angle_rad_nd
  real ( kind = 8 ) degrees_to_radians
  integer ( kind = 4 ) i
  real ( kind = 8 ) radians_to_degrees
  real ( kind = 8 ) temp1
  real ( kind = 8 ) temp2
  real ( kind = 8 ) temp3
  real ( kind = 8 ) thetad
  real ( kind = 8 ) thetar
  real ( kind = 8 ), dimension(dim_num) :: v1 = (/ 1.0D+00, 0.0D+00 /)
  real ( kind = 8 ) v2(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: v3 = (/ 0.0D+00, 0.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST002'
  write ( *, '(a)' ) '  ANGLE_DEG_2D computes an angle;'
  write ( *, '(a)' ) '  ANGLE_RAD_ND computes an angle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X  Y  Theta  ATAN2(y, x), ' // &
    'ANGLE_RAD_ND, ANGLE_DEG_2D'
  write ( *, '(a)' ) ' '

  do i = 0, angle_num

    thetad = real ( i, kind = 8 ) * 360.0D+00 / real ( angle_num, kind = 8 )
    thetar = degrees_to_radians ( thetad )

    v2(1) = cos ( thetar )
    v2(2) = sin ( thetar )

    temp1 = radians_to_degrees ( atan2 ( v2(2), v2(1) ) )

    temp2 = angle_rad_nd ( dim_num, v1, v2 )

    temp3 = angle_deg_2d ( v1, v3, v2 )

    write ( *, '(2x,7f10.3)') v2(1:2), thetad, temp1, temp2, temp3
 
  end do
 
  return
end
subroutine test0023 ( )

!*****************************************************************************80
!
!! TEST0023 tests ANGLE_HALF_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) angle_deg
  real ( kind = 8 ) cos_deg
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) p3(dim_num)
  real ( kind = 8 ) p4(dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ) sin_deg

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0023'
  write ( *, '(a)' ) '  ANGLE_HALF_2D computes the half angle between two rays;'
  write ( *, '(a)' ) '  The angle is defined by the points (P1,P2,P3)'
  write ( *, '(a)' ) '  or by the rays P2-->P3, P2-->P1.'

  p2(1:dim_num) = (/ 5.0D+00, 3.0D+00 /)

  angle_deg = 75.0D+00
  r = 3.0D+00
  p1(1:dim_num) = p2(1:dim_num) &
    + r * (/ cos_deg ( angle_deg ), sin_deg ( angle_deg ) /)

  angle_deg = 15.0D+00
  r = 2.0D+00
  p3(1:dim_num) = p2(1:dim_num) &
    + r * (/ cos_deg ( angle_deg ), sin_deg ( angle_deg ) /)

  call r8vec_print ( dim_num, p1, '  Point P1:' )
  call r8vec_print ( dim_num, p2, '  Point P2:' )
  call r8vec_print ( dim_num, p3, '  Point P3:' )

  call angle_half_2d ( p1, p2, p3, p4 )

  call r8vec_print ( dim_num, p4, &
    '  End point of unit ray from P2, defining half angle, P4:' )

  angle_deg = 45.0D+00
  r = 1.0D+00
  p4(1:dim_num) = p2(1:dim_num) &
    + r * (/ cos_deg ( angle_deg ), sin_deg ( angle_deg ) /)

  call r8vec_print ( dim_num, p4, '  Expected value of P4:' )

  return
end
subroutine test0025 ( )

!*****************************************************************************80
!
!! TEST0025 tests ANGLE_RAD_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 6

  real ( kind = 8 ) angle_rad
  real ( kind = 8 ) angle_rad_2d
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p1_test = reshape ( (/ &
    1.0D+00,  0.0D+00, &
    1.0D+00,  0.0D+00, &
    1.0D+00, -1.0D+00, &
    1.0D+00,  0.0D+00, &
    1.0D+00,  0.0D+00, &
    1.0D+00,  0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p2_test = reshape ( (/ &
    0.0D+00,  0.0D+00, &
    0.0D+00,  0.0D+00, &
    0.0D+00,  0.0D+00, &
    0.0D+00,  0.0D+00, &
    0.0D+00,  0.0D+00, &
    0.0D+00,  0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) p3(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p3_test = reshape ( (/ &
    1.0D+00,  1.0D+00, &
    0.0D+00,  1.0D+00, &
    0.0D+00,  1.0D+00, &
   -1.0D+00,  0.0D+00, &
    0.0D+00, -1.0D+00, &
    1.0D+00, -0.01D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0025'
  write ( *, '(a)' ) '  ANGLE_RAD_2D computes the angle between two rays;'

  do test = 1, test_num

    p1(1:dim_num) = p1_test(1:dim_num,test)
    p2(1:dim_num) = p2_test(1:dim_num,test)
    p3(1:dim_num) = p3_test(1:dim_num,test)
 
    angle_rad = angle_rad_2d ( p1, p2, p3 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Angle = ', angle_rad

  end do
 
  return
end
subroutine test003 ( )

!*****************************************************************************80
!
!! TEST003 tests ANGLE_RAD_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) angle_rad_3d
  real ( kind = 8 ) radians_to_degrees
  real ( kind = 8 ), dimension ( dim_num ) :: p1
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p1_test = reshape ( (/ &
    1.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 2.0D+00, 3.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num ) :: p2 = (/ &
    0.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p3 = (/ &
    0.0D+00, 0.0D+00, 1.0D+00 /)
  real ( kind = 8 ) temp1
  real ( kind = 8 ) temp2
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST003'
  write ( *, '(a)' ) '  ANGLE_RAD_3D computes an angle;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  P, ANGLE_RAD_3D, (Degrees)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p1(1:dim_num) = p1_test(1:dim_num,test)

    temp1 = angle_rad_3d ( p1, p2, p3 )
    temp2 = radians_to_degrees ( temp1 )

    write ( *, '(2x,6g12.4)') p1(1:dim_num), temp1, temp2
 
  end do
 
  return
end
subroutine test0032 ( )

!*****************************************************************************80
!
!! TEST0032 tests ANGLE_TURN_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 13

  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: p2 = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p3 = (/ 1.0D+00, 0.0D+00 /)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) test
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta_degrees
  real ( kind = 8 ) turn

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0032'
  write ( *, '(a)' ) '  ANGLE_TURN_2D computes the turning angle '
  write ( *, '(a)' ) '  defined by the line segments [P1,P2] and [P2,P3].'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our three points are:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    P1 = (C,S)'
  write ( *, '(a)' ) '    P2 = (0,0)'
  write ( *, '(a)' ) '    P3 = (1,0)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  C = cosine ( theta ), S = sine ( theta ).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test  Theta    Turn'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    theta = 2.0D+00 * pi * real ( test     - 1, kind = 8 ) &
                         / real ( test_num - 1, kind = 8 )

    theta_degrees = 360.0D+00 * real ( test     - 1, kind = 8 ) &
                         /      real ( test_num - 1, kind = 8 )

    p1(1:dim_num) = (/ cos ( theta ), sin ( theta ) /)

    call angle_turn_2d ( p1, p2, p3, turn )

    write ( *, '(2x,i4,2x,f5.0,2x,g14.6)' ) test, theta_degrees, turn

  end do

  return
end
subroutine test0035 ( )

!*****************************************************************************80
!
!! TEST0035 tests ANNULUS_SECTOR_CENTROID_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 December 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) centroid(dim_num)
  real ( kind = 8 ) degrees_to_radians
  real ( kind = 8 ), dimension(dim_num) :: pc = (/ 5.0D+00, 3.0D+00 /)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) :: r1 = 2.0D+00
  real ( kind = 8 ) :: r2 = 3.0D+00
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2

  theta1 = degrees_to_radians ( 30.0D+00 )
  theta2 = degrees_to_radians ( 60.0D+00 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0035'
  write ( *, '(a)' ) '  ANNULUS_SECTOR_CENTROID_2D computes the centroid of a'
  write ( *, '(a)' ) '  circular annulus.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  The circle has center        ', pc(1:2)
  write ( *, '(a,g14.6)'  ) '  The inner radius is R1 =     ', r1
  write ( *, '(a,g14.6)'  ) '  The outer radius is R2 =     ', r2
  write ( *, '(a,g14.6)'  ) '  The first angle is THETA1 =  ', theta1
  write ( *, '(a,g14.6)'  ) '  The second angle is THETA2 = ', theta2

  call annulus_sector_centroid_2d ( pc, r1, r2, theta1, theta2, centroid )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2f14.8)' ) '  Centroid: ', centroid(1:2)

  return
end
subroutine test004 ( )

!*****************************************************************************80
!
!! TEST004 tests R8_ACOS;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 9

  real ( kind = 8 ) r8_acos
  real ( kind = 8 ) radians_to_degrees
  real ( kind = 8 ) temp1
  real ( kind = 8 ) temp2
  integer ( kind = 4 ) test
  real ( kind = 8 ) x
  real ( kind = 8 ), dimension ( test_num ) :: x_test = (/ &
    5.0D+00, 1.2D+00, 1.0D+00, 0.9D+00, 0.5D+00, &
    0.0D+00, -0.9D+00, -1.0D+00, -1.01D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST004'
  write ( *, '(a)' ) '  R8_ACOS computes an angle with a given cosine;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X  R8_ACOS(X)  (Degrees)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    x = x_test(test)

    temp1 = r8_acos ( x )
    temp2 = radians_to_degrees ( temp1 )

    write ( *, '(2x,6g12.4)') x, temp1, temp2
 
  end do
 
  return
end
subroutine test0045 ( )

!*****************************************************************************80
!
!! TEST0045 tests R8_ASIN;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 9

  real ( kind = 8 ) r8_asin
  real ( kind = 8 ) radians_to_degrees
  real ( kind = 8 ) temp1
  real ( kind = 8 ) temp2
  integer ( kind = 4 ) test
  real ( kind = 8 ) x
  real ( kind = 8 ), dimension ( test_num ) :: x_test = (/ &
    5.0D+00, 1.2D+00, 1.0D+00, 0.9D+00, 0.5D+00, &
    0.0D+00, -0.9D+00, -1.0D+00, -1.01D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0045'
  write ( *, '(a)' ) '  R8_ASIN computes an angle with a given sine;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X  R8_ASIN(X)  (Degrees)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    x = x_test(test)

    temp1 = r8_asin ( x )
    temp2 = radians_to_degrees ( temp1 )

    write ( *, '(2x,6g12.4)') x, temp1, temp2
 
  end do
 
  return
end
subroutine test005 ( )

!*****************************************************************************80
!
!! TEST005 tests R8_ATAN;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 8

  real ( kind = 8 ) r8_atan
  real ( kind = 8 ) radians_to_degrees
  real ( kind = 8 ) temp1
  real ( kind = 8 ) temp2
  real ( kind = 8 ) temp3
  integer ( kind = 4 ) test
  real ( kind = 8 ) x
  real ( kind = 8 ), dimension ( test_num ) :: x_test = (/ &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
    -1.0D+00, -1.0D+00, 0.0D+00 /)
  real ( kind = 8 ) y
  real ( kind = 8 ), dimension ( test_num ) :: y_test = (/ &
    0.0D+00, 1.0D+00, 2.0D+00, 0.0D+00, -1.0D+00, &
    -1.0D+00, -1.0D+00, -1.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005'
  write ( *, '(a)' ) '  R8_ATAN computes an angle with a given tangent.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X  Y  ATAN(Y/X)  ATAN2(Y,X)  R8_ATAN(Y,X)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    x = x_test(test)
    y = y_test(test)

    if ( x /= 0.0D+00 ) then
      temp1 = atan ( y / x )
    else
      temp1 = huge ( y )
    end if

    temp2 = atan2 ( y, x )
    temp3 = r8_atan ( y, x )
    
    write ( *, '(2x,6g12.4)') x, y, temp1, temp2, temp3
 
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat, but display answers in degrees.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    x = x_test(test)
    y = y_test(test)

    if ( x /= 0.0D+00 ) then
      temp1 = radians_to_degrees ( atan ( y / x ) )
    else
      temp1 = huge ( y )
    end if
    
    temp2 = radians_to_degrees ( atan2 ( y, x ) )
    temp3 = radians_to_degrees ( r8_atan ( y, x ) )
    
    write ( *, '(2x,6g12.4)') x, y, temp1, temp2, temp3
 
  end do
 
  return
end
subroutine test006 ( )

!*****************************************************************************80
!
!! TEST006 tests BALL_UNIT_SAMPLE_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) average(dim_num)
  real ( kind = 8 ) average_r
  real ( kind = 8 ) average_theta
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: sample_num = 1000
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_atan
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST006'
  write ( *, '(a)' ) '  For the unit ball in 2 dimensions (the disk):'
  write ( *, '(a)' ) '  BALL_UNIT_SAMPLE_2D samples;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A few sample values:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call ball_unit_sample_2d ( seed, x )
    write ( *, '(2x,2f8.4)' ) x(1:dim_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of sample points = ', sample_num

  average(1:dim_num) = 0.0D+00

  do i = 1, sample_num
    call ball_unit_sample_2d ( seed, x )
    average(1:dim_num) = average(1:dim_num) + x(1:dim_num)
  end do

  average(1:dim_num) = average(1:dim_num) / real ( sample_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now average the points, which should get a value'
  write ( *, '(a)' ) '  close to zero, and closer as sample_num increases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2f8.4)' ) '  Average:        ', average(1:dim_num)

  average_r = 0.0D+00

  do i = 1, sample_num
    call ball_unit_sample_2d ( seed, x )
    average_r = average_r + sqrt ( sum ( x(1:dim_num)**2 ) )
  end do

  average_r = average_r / real ( sample_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now average the distance of the points from' 
  write ( *, '(a,f8.4)' ) '  the center, which should be 1/sqrt(2) = ', &
    1.0D+00 / sqrt ( 2.0D+00 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,2f8.4)' ) '  Average:        ', average_r

  average_theta = 0.0D+00

  do i = 1, sample_num
    call ball_unit_sample_2d ( seed, x )
    theta = r8_atan ( x(2), x(1) )
    average_theta = average_theta + theta
  end do

  average_theta = average_theta / real ( sample_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now average the angle THETA,'
  write ( *, '(a)' ) '  which should be PI.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2f8.4)' ) '  Average:        ', average_theta

  return
end
subroutine test007 ( )

!*****************************************************************************80
!
!! TEST007 tests BALL_UNIT_SAMPLE_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) average(dim_num)
  real ( kind = 8 ) average_phi
  real ( kind = 8 ) average_r
  real ( kind = 8 ) average_theta
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: sample_num = 1000
  real ( kind = 8 ) phi
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_atan
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST007'
  write ( *, '(a)' ) '  For the unit ball in 3 dimensions:'
  write ( *, '(a)' ) '  BALL_UNIT_SAMPLE_3D samples;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A few sample values:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call ball_unit_sample_3d ( seed, x )
    write ( *, '(2x,3f8.4)' ) x(1:dim_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of sample points = ', sample_num

  average(1:dim_num) = 0.0D+00

  do i = 1, sample_num
    call ball_unit_sample_3d ( seed, x )
    average(1:dim_num) = average(1:dim_num) + x(1:dim_num)
  end do

  average(1:dim_num) = average(1:dim_num) / real ( sample_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now average the points, which should get a value'
  write ( *, '(a)' ) '  close to zero, and closer as sample_num increases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,3f8.4)' ) '  Average:        ', average(1:dim_num)

  seed = 123456789

  average_r = 0.0D+00

  do i = 1, sample_num
    call ball_unit_sample_3d ( seed, x )
    r = sqrt ( sum ( x(1:dim_num)**2 ) )
    average_r = average_r + r
  end do

  average_r = average_r / real ( sample_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now average the distance of the points from'
  write ( *, '(a,f8.4)' ) '  the center, which should be the '
  write ( *, '(a,f8.4)' ) '  1/2**(1/dim_num) = ', &
    0.5D+00**( 1.0D+00 / real ( dim_num, kind = 8 ) )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4)' ) '  Average:        ', average_r

  average_theta = 0.0D+00

  do i = 1, sample_num
    call ball_unit_sample_3d ( seed, x )
    theta = r8_atan ( x(2), x(1) )
    average_theta = average_theta + theta
  end do

  average_theta = average_theta / real ( sample_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now average the angle THETA,'
  write ( *, '(a)' ) '  which should be PI.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4)' ) '  Average:        ', average_theta

  average_phi = 0.0D+00

  do i = 1, sample_num
    call ball_unit_sample_3d ( seed, x )
    r = sqrt ( sum ( x(1:dim_num)**2 ) )
    phi = acos ( x(3) / r )
    average_phi = average_phi + phi
  end do

  average_phi = average_phi / real ( sample_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now average the angle PHI,'
  write ( *, '(a)' ) '  which should be PI/2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4)' ) '  Average:        ', average_phi

  return
end
subroutine test0075 ( )

!*****************************************************************************80
!
!! TEST0075 tests BALL_UNIT_SAMPLE_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) average(dim_num)
  real ( kind = 8 ) average_phi
  real ( kind = 8 ) average_r
  real ( kind = 8 ) average_theta
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: sample_num = 1000
  real ( kind = 8 ) phi
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_atan
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0075'
  write ( *, '(a)' ) '  For the unit ball in N dimensions:'
  write ( *, '(a)' ) '  BALL_UNIT_SAMPLE_ND samples;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A few sample values:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call ball_unit_sample_nd ( dim_num, seed, x )
    write ( *, '(2x,3f8.4)' ) x(1:dim_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of sample points = ', sample_num

  average(1:dim_num) = 0.0D+00

  do i = 1, sample_num
    call ball_unit_sample_nd ( dim_num, seed, x )
    average(1:dim_num) = average(1:dim_num) + x(1:dim_num)
  end do

  average(1:dim_num) = average(1:dim_num) / real ( sample_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now average the points, which should get a value'
  write ( *, '(a)' ) '  close to zero, and closer as N increases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,3f8.4)' ) '  Average:        ', average(1:dim_num)

  seed = 123456789

  average_r = 0.0D+00

  do i = 1, sample_num
    call ball_unit_sample_nd ( dim_num, seed, x )
    r = sqrt ( sum ( x(1:dim_num)**2 ) )
    average_r = average_r + r
  end do

  average_r = average_r / real ( sample_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now average the distance of the points from'
  write ( *, '(a,f8.4)' ) '  the center, which should be the '
  write ( *, '(a,f8.4)' ) '  1/2**(1/dim_num) = ', &
    0.5D+00**( 1.0D+00 / real ( dim_num, kind = 8 ) )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4)' ) '  Average:        ', average_r

  average_theta = 0.0D+00

  do i = 1, sample_num
    call ball_unit_sample_nd ( dim_num, seed, x )
    theta = r8_atan ( x(2), x(1) )
    average_theta = average_theta + theta
  end do

  average_theta = average_theta / real ( sample_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now average the angle THETA,'
  write ( *, '(a)' ) '  which should be PI.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4)' ) '  Average:        ', average_theta

  average_phi = 0.0D+00

  do i = 1, sample_num
    call ball_unit_sample_nd ( dim_num, seed, x )
    r = sqrt ( sum ( x(1:dim_num)**2 ) )
    phi = acos ( x(3) / r )
    average_phi = average_phi + phi
  end do

  average_phi = average_phi / real ( sample_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now average the angle PHI,'
  write ( *, '(a)' ) '  which should be PI/2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4)' ) '  Average:        ', average_phi

  return
end
subroutine test008 ( )

!*****************************************************************************80
!
!! TEST008 tests BASIS_MAP_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a(3,3)
  real ( kind = 8 ) c(3,3)
  integer ( kind = 4 ) ierror
  real ( kind = 8 ), dimension ( 3, 3 ) :: u = reshape ( (/ &
    1.0D+00, 2.0D+00, 3.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00, &
    1.0D+00, 0.0D+00, 2.0D+00 /), (/ 3, 3 /) )
  real ( kind = 8 ), dimension ( 3, 3 ) :: v = reshape ( (/ &
    14.0D+00, 4.0D+00, 4.0D+00, &
     3.0D+00, 1.0D+00, 0.0D+00, &
     7.0D+00, 3.0D+00, 2.0D+00 /), (/ 3, 3 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST008'
  write ( *, '(a)' ) '  BASIS_MAP_3D computes the linear transform A'
  write ( *, '(a)' ) '  which maps vectors U1, U2 and U3 to vectors'
  write ( *, '(a)' ) '  V1, V2 and V3.'

  call r8mat_print ( 3, 3, u, '  The matrix U' )

  call r8mat_print ( 3, 3, v, '  The matrix V' )

  call basis_map_3d ( u, v, a, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The matrix [ U1 | U2 | U3 ] was singular.'
    write ( *, '(a)' ) '  No transformation was computed.'
    return
  end if

  call r8mat_print ( 3, 3, a, '  The transformation matrix' )

  c(1:3,1:3) = matmul ( a(1:3,1:3), u(1:3,1:3) )

  call r8mat_print ( 3, 3, c, &
    '  The product matrix A * [ U1 | U2 | U3 ]' )

  return
end
subroutine test0085 ( )

!*****************************************************************************80
!
!! TEST0085 tests BOX_01_CONTAINS_POINT_2D
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

  integer ( kind = 4 ), parameter :: n = 46
  integer ( kind = 4 ), parameter :: dim_num = 2

  logical box_01_contains_point_2d
  character dot(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), parameter :: xhi =  1.2D+00
  real ( kind = 8 ), parameter :: xlo = -0.3D+00
  real ( kind = 8 ), parameter :: yhi =  1.4D+00
  real ( kind = 8 ), parameter :: ylo = -0.1D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0085'
  write ( *, '(a)' ) '  BOX_01_CONTAINS_POINT_2D reports if the unit box'
  write ( *, '(a)' ) '  contains a point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We will call the function repeatedly, and draw'
  write ( *, '(a)' ) '  a sketch of the unit square.'
  write ( *, '(a)' ) ' '

  do i = 1, n

    p(2) = ( real ( n - i,     kind = 8 ) * yhi   &
           + real (     i - 1, kind = 8 ) * ylo ) &
           / real ( n     - 1, kind = 8 )

    do j = 1, n

      p(1) = ( real ( n - j,     kind = 8 ) * xlo   &
             + real (     j - 1, kind = 8 ) * xhi ) &
             / real ( n     - 1, kind = 8 )

      if ( box_01_contains_point_2d ( p ) ) then
        dot(j) = '*'
      else
        dot(j) = '-'
      end if

    end do
    write ( *, '(2x,46a1)' ) dot(1:n)
  end do

  return
end
subroutine test0087 ( )

!*****************************************************************************80
!
!! TEST0087 tests BOX_CONTAINS_POINT_2D
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 46
  integer ( kind = 4 ), parameter :: dim_num = 2

  logical box_contains_point_2d
  character dot(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: p1 = (/ -0.1D+00, 0.3D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p2 = (/  1.1D+00, 0.9D+00 /)
  real ( kind = 8 ), parameter :: xhi =  1.2D+00
  real ( kind = 8 ), parameter :: xlo = -0.3D+00
  real ( kind = 8 ), parameter :: yhi =  1.4D+00
  real ( kind = 8 ), parameter :: ylo = -0.1D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0087'
  write ( *, '(a)' ) '  BOX_CONTAINS_POINT_2D reports if a box'
  write ( *, '(a)' ) '  contains a point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We will call the function repeatedly, and draw'
  write ( *, '(a)' ) '  a sketch of the box.'
  write ( *, '(a)' ) ' '

  do i = 1, n

    p(2) = ( real ( n - i,     kind = 8 ) * yhi   &
           + real (     i - 1, kind = 8 ) * ylo ) &
           / real ( n     - 1, kind = 8 )

    do j = 1, n

      p(1) = ( real ( n - j,     kind = 8 ) * xlo   &
             + real (     j - 1, kind = 8 ) * xhi ) &
             / real ( n     - 1, kind = 8 )

      if ( box_contains_point_2d ( p1, p2, p ) ) then
        dot(j) = '*'
      else
        dot(j) = '-'
      end if

    end do
    write ( *, '(2x,45a1)' ) dot(1:n)
  end do

  return
end
subroutine test009 ( )

!*****************************************************************************80
!
!! TEST009 tests BOX_SEGMENT_CLIP_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 5

  integer ( kind = 4 ) ival
  real ( kind = 8 ), dimension(dim_num) :: p1 = (/ -10.0D+00, 10.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p2 = (/  10.0D+00, 20.0D+00 /)
  real ( kind = 8 ) pa(dim_num)
  real ( kind = 8 ) pb(dim_num)
  real ( kind = 8 ) qa(dim_num)
  real ( kind = 8 ) qb(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p1_test = reshape ( (/ &
     1.0D+00,  2.0D+00, &
    -3.0D+00, 12.0D+00, &
   -20.0D+00, 20.0D+00, &
   -20.0D+00, 40.0D+00, &
    10.0D+00, 40.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p2_test = reshape ( (/ &
     8.0D+00, 16.0D+00, &
     5.0D+00, 12.0D+00, &
     7.0D+00, 20.0D+00, &
     0.0D+00,  0.0D+00, &
    20.0D+00, 30.0D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST009'
  write ( *, '(a)' ) '  BOX_SEGMENT_CLIP_2D clips a line with respect to a box.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The lower left box corner is:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,4f8.4)' ) p1(1:dim_num)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The upper right box corner is:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,4f8.4)' ) p2(1:dim_num)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We list the points PA and PB, and then'
  write ( *, '(a)' ) '  the clipped values.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    pa(1:2) = p1_test(1:2,test)
    pb(1:2) = p2_test(1:2,test)

    qa(1:2) = pa(1:2)
    qb(1:2) = pb(1:2)
    
    call box_segment_clip_2d ( p1, p2, qa, qb, ival )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,4f8.4)' ) pa(1:dim_num), pb(1:dim_num)
    if ( ival == -1 ) then
      write ( *, '(a)' ) '  Line is outside the box.'
    else if ( ival == 0 ) then
      write ( *, '(a)' ) '  Line is inside the box.'
    else if ( ival == 1 ) then
      write ( *, '(2x,2f8.4)' ) qa(1:dim_num)
    else if ( ival == 2 ) then
      write ( *, '(2x,16x,2f8.4)' )         qb(1:dim_num)
    else if ( ival == 3 ) then
      write ( *, '(2x,4f8.4)' ) qa(1:dim_num), qb(1:dim_num)
    end if

  end do

  return
end
subroutine test010 ( )

!*****************************************************************************80
!
!! TEST010 tests BOX_RAY_INT_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ), dimension(dim_num) :: p1 = (/  0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p2 = (/  5.0D+00, 3.0D+00 /)
  real ( kind = 8 ) pa(1:dim_num)
  real ( kind = 8 ), dimension(1:dim_num,test_num) :: pa_test = reshape ( (/ &
    3.0D+00, 1.0D+00, &
    4.0D+00, 1.0D+00, &
    3.0D+00, 1.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) pb(1:dim_num)
  real ( kind = 8 ), dimension(1:dim_num,test_num) :: pb_test = reshape ( (/ &
    5.0D+00, 5.0D+00, &
    3.0D+00, 1.0D+00, &
    4.0D+00, 2.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension(1:dim_num,test_num) :: pc_test = reshape ( (/ &
    4.0D+00, 3.0D+00, &
    0.0D+00, 1.0D+00, &
    5.0D+00, 3.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) pint(1:dim_num)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST010'
  write ( *, '(a)' ) '  For a box with coordinate line sides in 2D,'
  write ( *, '(a)' ) '  BOX_RAY_INT_2D computes the intersection of'
  write ( *, '(a)' ) '  a shape and a ray whose origin is within'
  write ( *, '(a)' ) '  the shape.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Lower left box corner:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,2g14.6)' ) p1(1:dim_num)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Upper right box corner:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,2g14.6)' ) p2(1:dim_num)
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    pa(1:2) = pa_test(1:2,test)
    pb(1:2) = pb_test(1:2,test)

    call box_ray_int_2d ( p1, p2, pa, pb, pint )

    write ( *, '(a)' ) ' '
    write ( *, '(a,2f12.4)' ) '  Origin:       ', pa(1:2)
    write ( *, '(a,2f12.4)' ) '  Point 2:      ', pb(1:2)
    write ( *, '(a,2f12.4)' ) '  Intersection: ', pint(1:2)
    write ( *, '(a,2f12.4)' ) '  Correct:      ', pc_test(1:2,test)

  end do

  return
end
subroutine test011 ( )

!*****************************************************************************80
!
!! TEST011 tests CIRCLE_DIA2IMP_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ) theta

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST011'
  write ( *, '(a)' ) '  CIRCLE_DIA2IMP_2D converts a diameter to an'
  write ( *, '(a)' ) '  implicit circle in 2D.'
 
  theta = 2.0D+00

  p1(1:dim_num) = 2.0D+00 + 5.0D+00 * (/ cos ( theta ),  sin ( theta ) /)
  p2(1:dim_num) = 2.0D+00 - 5.0D+00 * (/ cos ( theta ),  sin ( theta ) /)

  call r8vec_print ( dim_num, p1, '  P1:' )
  call r8vec_print ( dim_num, p2, '  P2:' )

  call circle_dia2imp_2d ( p1, p2, r, pc )

  call circle_imp_print_2d ( r, pc, '  The implicit circle:' )
  
  return
end
subroutine test012 ( )

!*****************************************************************************80
!
!! TEST012 tests CIRCLE_LUNE_AREA_2D, CIRCLE_SECTOR_AREA_2D,  CIRCLE_TRIANGLE_AREA_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) area1
  real ( kind = 8 ) area2
  real ( kind = 8 ) area3
  real ( kind = 8 ), dimension(dim_num) :: pc = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) :: r = 1.0D+00
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 12
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST012'
  write ( *, '(a)' ) '  CIRCLE_LUNE_AREA_2D computes the area of a'
  write ( *, '(a)' ) '  circular lune, defined by joining the endpoints'
  write ( *, '(a)' ) '  of a circular arc.'
  write ( *, '(a)' ) '  CIRCLE_SECTOR_AREA_2D computes the area of a'
  write ( *, '(a)' ) '  circular sector, defined by joining the endpoints'
  write ( *, '(a)' ) '  of a circular arc to the center.'
  write ( *, '(a)' ) '  CIRCLE_TRIANGLE_AREA_2D computes the signed area of a'
  write ( *, '(a)' ) '  triangle, defined by joining the endpoints'
  write ( *, '(a)' ) '  of a circular arc and the center.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      R            Theta1      Theta2        Sector       Triangle     Lune'
  write ( *, '(a)' ) ' '

  do test = 0, test_num

    theta1 = 0.0D+00
    theta2 = real ( test, kind = 8 ) * 2.0D+00 * pi &
           / real ( test_num, kind = 8 )

    call circle_sector_area_2d ( r, pc, theta1, theta2, area1 )

    call circle_triangle_area_2d ( r, pc, theta1, theta2, area2 )

    call circle_lune_area_2d ( r, pc, theta1, theta2, area3 )

    write ( *, '(2x,6f12.6)' ) r, theta1, theta2, area1, area2, area3

  end do

  return
end
subroutine test0125 ( )

!*****************************************************************************80
!
!! TEST0125 tests CIRCLE_LUNE_AREA_2D and SPHERE_CAP_VOLUME_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) area
  real ( kind = 8 ) h
  real ( kind = 8 ) haver_sine
  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_asin
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 12
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  real ( kind = 8 ) volume1
  real ( kind = 8 ) volume2

  pc(1:2) = (/ 0.0D+00, 0.0D+00 /)
  r = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0125'
  write ( *, '(a)' ) '  CIRCLE_LUNE_AREA_2D computes the area of a'
  write ( *, '(a)' ) '  circular lune, defined by joining the endpoints'
  write ( *, '(a)' ) '  of a circular arc (THETA1,THETA2).'
  write ( *, '(a)' ) '  SPHERE_CAP_VOLUME_2D computes the volume (area) of a'
  write ( *, '(a)' ) '  spherical cap, defined by a plane that cuts the'
  write ( *, '(a)' ) '  sphere to a thickness of H units.'
  write ( *, '(a)' ) '  SPHERE_CAP_VOLUME_ND does the same operation,'
  write ( *, '(a)' ) '  but in N dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The two routines should get the same results'
  write ( *, '(a)' ) '  if THETA1, THETA2 and H correspond.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Using a radius R = ', r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '        Theta1      Theta2      H           Lune        Cap        Cap'
  write ( *, '(a)' ) &
    '                                            area        vol_3d     vol_nd'
  write ( *, '(a)' ) ' '

  do test = 0, test_num

    h = 2.0D+00 * r * real ( test, kind = 8 ) / real ( test_num, kind = 8 )

    haver_sine = sqrt ( r * r - ( r - h )**2 )

    if ( h <= r ) then
      theta2 = r8_asin ( haver_sine / r )
    else
      theta2 = ( pi - r8_asin ( haver_sine / r ) )
    end if

    theta1 = -theta2

    call circle_lune_area_2d ( r, pc, theta1, theta2, area )

    call sphere_cap_volume_2d ( r, h, volume1 )

    call sphere_cap_volume_nd ( dim_num, r, h, volume2 )

    write ( *, '(2x,6f12.6)' ) theta1, theta2, h, area, volume1, volume2

  end do

  return
end
subroutine test0126 ( )

!*****************************************************************************80
!
!! TEST0126 tests SPHERE_CAP_VOLUME_3D and SPHERE_CAP_VOLUME_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) h
  real ( kind = 8 ) :: r = 1.0D+00
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 12
  real ( kind = 8 ) volume1
  real ( kind = 8 ) volume2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0126'
  write ( *, '(a)' ) '  SPHERE_CAP_VOLUME_3D computes the volume of a'
  write ( *, '(a)' ) '  spherical cap, defined by a plane that cuts the'
  write ( *, '(a)' ) '  sphere to a thickness of H units.'
  write ( *, '(a)' ) '  SPHERE_CAP_VOLUME_ND does the same operation,'
  write ( *, '(a)' ) '  but in N dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Using a radius R = ', r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        H           Cap        Cap'
  write ( *, '(a)' ) '                    volume_3d  volume_nd'
  write ( *, '(a)' ) ' '

  do test = 0, test_num

    h = 2.0D+00 * r * real ( test, kind = 8 ) / real ( test_num, kind = 8 )

    call sphere_cap_volume_3d ( r, h, volume1 )

    call sphere_cap_volume_nd ( dim_num, r, h, volume2 )

    write ( *, '(2x,3f12.6)' ) h, volume1, volume2

  end do

  return
end
subroutine test0127 ( )

!*****************************************************************************80
!
!! TEST0127 tests SPHERE_CAP_AREA_3D and SPHERE_CAP_AREA_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) area1
  real ( kind = 8 ) area2
  real ( kind = 8 ) h
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 12
  real ( kind = 8 ) r

  r = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0127'
  write ( *, '(a)' ) '  SPHERE_CAP_AREA_3D computes the volume of a'
  write ( *, '(a)' ) '  3D spherical cap, defined by a plane that cuts the'
  write ( *, '(a)' ) '  sphere to a thickness of H units.'
  write ( *, '(a)' ) '  SPHERE_CAP_AREA_ND computes the volume of an'
  write ( *, '(a)' ) '  ND spherical cap, defined by a plane that cuts the'
  write ( *, '(a)' ) '  sphere to a thickness of H units.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        R           H           Cap         Cap'
  write ( *, '(a)' ) '                                area_3d     area_nd'
  write ( *, '(a)' ) ' '

  do test = 0, test_num

    h = 2.0D+00 * r * real ( test, kind = 8 ) / real ( test_num, kind = 8 )

    call sphere_cap_area_3d ( r, h, area1 )

    call sphere_cap_area_nd ( dim_num, r, h, area2 )

    write ( *, '(2x,5f12.6)' ) r, h, area1, area2

  end do

  return
end
subroutine test013 ( )

!*****************************************************************************80
!
!! TEST013 tests CIRCLE_LUNE_CENTROID_2D and CIRCLE_SECTOR_CENTROID_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) centroid1(dim_num)
  real ( kind = 8 ) centroid2(dim_num)
  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 12
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2

  r = 2.0D+00
  pc(1:2) = (/ 5.0D+00, 3.0D+00 /)
  theta1 = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST013'
  write ( *, '(a)' ) '  CIRCLE_LUNE_CENTROID_2D computes the centroid of a'
  write ( *, '(a)' ) '  circular lune, defined by joining the endpoints'
  write ( *, '(a)' ) '  of a circular arc.'
  write ( *, '(a)' ) '  CIRCLE_SECTOR_CENTROID_2D computes the centroid of a'
  write ( *, '(a)' ) '  circular sector, defined by joining the endpoints'
  write ( *, '(a)' ) '  of a circular arc to the center.'

  call circle_imp_print_2d ( r, pc, '  The implicit circle:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first angle of our lune and sector is always 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '                         Lune                       Sector'
  write ( *, '(a)' ) &
    '  THETA2           X             Y             X             Y'
  write ( *, '(a)' ) ' '

  do test = 0, test_num

    theta2 = real ( test, kind = 8 ) * 2.0D+00 * pi &
      / real ( test_num, kind = 8 )

    call circle_lune_centroid_2d ( r, pc, theta1, theta2, centroid1 )

    call circle_sector_centroid_2d ( r, pc, theta1, theta2, centroid2 )

    write ( *, '(2x,5f14.8)' ) theta2, centroid1(1:2), centroid2(1:2)

  end do

  return
end
subroutine test014 ( )

!*****************************************************************************80
!
!! TEST014 tests CIRCLE_EXP_CONTAINS_POINT_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) inside
  character ( len = 60 ) message(-1:7)
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) p3(dim_num)

  message(-1) = 'The point is inside the circle.'
  message(0) = 'The point is on the circle.'
  message(1) = 'The point is outside the circle'
  message(2) = 'Colinear data, the point is on the line.'
  message(3) = 'Colinear data, the point is not on the line.'
  message(4) = 'Two equal data points, the point is on the line.'
  message(5) = 'Two equal data points, the point is not on the line.'
  message(6) = 'All data points equal, the point is equal.'
  message(7) = 'All data points equal, the point is not equal.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST014'
  write ( *, '(a)' ) '  CIRCLE_EXP_CONTAINS_POINT_2D determines if a'
  write ( *, '(a)' ) '  point lies inside a circle.'
!
!  This point is inside.
!
  p1(1:dim_num) = (/  4.0D+00, 2.0D+00 /)
  p2(1:dim_num) = (/  1.0D+00, 5.0D+00 /)
  p3(1:dim_num) = (/ -2.0D+00, 2.0D+00 /)
  p(1:dim_num)  = (/  2.0D+00, 3.0D+00 /)

  call r8vec_print ( dim_num, p1, '  P1:' )
  call r8vec_print ( dim_num, p2, '  P2:' )
  call r8vec_print ( dim_num, p3, '  P3:' )
  call r8vec_print ( dim_num, p, '  P:' )

  call circle_exp_contains_point_2d ( p1, p2, p3, p, inside )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  INSIDE = ', inside
  write ( *, '(2x,a)' ) message(inside)
!
!  This point is actually right on the circle.
!
  p1(1:dim_num) = (/  4.0D+00,  2.0D+00 /)
  p2(1:dim_num) = (/  1.0D+00,  5.0D+00 /)
  p3(1:dim_num) = (/ -2.0D+00,  2.0D+00 /)
  p(1:dim_num)  = (/  1.0D+00, -1.0D+00 /)

  call r8vec_print ( dim_num, p1, '  P1:' )
  call r8vec_print ( dim_num, p2, '  P2:' )
  call r8vec_print ( dim_num, p3, '  P3:' )
  call r8vec_print ( dim_num, p, '  P:' )

  call circle_exp_contains_point_2d ( p1, p2, p3, p, inside )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  INSIDE = ', inside
  write ( *, '(2x,a)' ) message(inside)
!
!  This point is outside.
!
  p1(1:dim_num) = (/  4.0D+00, 2.0D+00 /)
  p2(1:dim_num) = (/  1.0D+00, 5.0D+00 /)
  p3(1:dim_num) = (/ -2.0D+00, 2.0D+00 /)
  p(1:dim_num)  = (/  4.0D+00, 6.0D+00 /)

  call r8vec_print ( dim_num, p1, '  P1:' )
  call r8vec_print ( dim_num, p2, '  P2:' )
  call r8vec_print ( dim_num, p3, '  P3:' )
  call r8vec_print ( dim_num, p, '  P:' )

  call circle_exp_contains_point_2d ( p1, p2, p3, p, inside )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  INSIDE = ', inside
  write ( *, '(2x,a)' ) message(inside)

  return
end
subroutine test015 ( )

!*****************************************************************************80
!
!! TEST015 tests CIRCLE_EXP2IMP_2D and TRIANGLE_CIRCUMCIRCLE_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) p3(dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ) t(dim_num,3)
  integer ( kind = 4 ) test
  real ( kind = 8 ), dimension(dim_num,3,test_num) :: test_t = reshape ( (/ &
    4.0D+00, 2.0D+00, &
    1.0D+00, 5.0D+00, &
   -2.0D+00, 2.0D+00, &
    4.0D+00, 2.0D+00, &
    5.0D+00, 4.0D+00, &
    6.0D+00, 6.0D+00, &
    4.0D+00, 2.0D+00, &
    1.0D+00, 5.0D+00, &
    4.0D+00, 2.0D+00 /), (/ dim_num, 3, test_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST015'
  write ( *, '(a)' ) '  CIRCLE_EXP2IMP_2D computes the radius and '
  write ( *, '(a)' ) '  center of the circle through three points.'
  write ( *, '(a)' ) '  TRIANGLE_CIRCUMCIRCLE_2D computes the radius and '
  write ( *, '(a)' ) '  center of the circle through the vertices of'
  write ( *, '(a)' ) '  a triangle.'

  do test = 1, test_num

    p1(1:dim_num) = test_t(1:dim_num,1,test)
    p2(1:dim_num) = test_t(1:dim_num,2,test)
    p3(1:dim_num) = test_t(1:dim_num,3,test)

    call r8vec_print ( dim_num, p1, '  P1:' )
    call r8vec_print ( dim_num, p2, '  P2:' )
    call r8vec_print ( dim_num, p3, '  P3:' )

    call circle_exp2imp_2d ( p1, p2, p3, r, pc )

    call circle_imp_print_2d ( r, pc, '  The implicit circle:' )

    t(1:dim_num,1:3) = test_t(1:dim_num,1:3,test)

    call triangle_circumcircle_2d ( t, r, pc )

    call circle_imp_print_2d ( r, pc, '  The triangle''s circumcircle:' )

  end do

  return
end
subroutine test0155 ( )

!*****************************************************************************80
!
!! TEST0155 tests CIRCLE_EXP2IMP_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 13

  real ( kind = 8 ) curvature
  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) p3(dim_num)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta_degrees
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0155'
  write ( *, '(a)' ) '  CIRCLE_EXP2IMP_2D computes the radius and '
  write ( *, '(a)' ) '  center of the circle through three points.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We can use this routine to compute, for three'
  write ( *, '(a)' ) '  points in space, the circle incident to those'
  write ( *, '(a)' ) '  points, and hence the radius of that circle,'
  write ( *, '(a)' ) '  and hence the "curvature" of those points.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our three points are:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    (0,0)'
  write ( *, '(a)' ) '    (1,0)'
  write ( *, '(a)' ) '    (C,S)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  C = cosine ( theta), S = sine ( theta ).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test  Theta  Curvature'
  write ( *, '(a)' ) ' '

  p1(1:dim_num) = (/ 0.0D+00, 0.0D+00 /)
  p2(1:dim_num) = (/ 1.0D+00, 0.0D+00 /)

  do test = 1, test_num

    theta = 2.0D+00 * pi * real ( test     - 1, kind = 8 ) &
                         / real ( test_num - 1, kind = 8 )

    theta_degrees = 360.0D+00 * real ( test     - 1, kind = 8 ) &
                         /      real ( test_num - 1, kind = 8 )

    p3(1:dim_num) = (/ cos ( theta ), sin ( theta ) /)

    call circle_exp2imp_2d ( p1, p2, p3, r, pc )

    if ( 0.0D+00 < r ) then
      curvature = 1.0D+00 / r
    else
      curvature = 0.0D+00
    end if

    write ( *, '(2x,i4,2x,f5.0,2x,g14.6)' ) test, theta_degrees, curvature

  end do

  return
end
subroutine test0156 ( )

!*****************************************************************************80
!
!! TEST0156 tests CIRCLE_EXP2IMP_2D and CIRCLE_IMP2EXP_2D.
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

  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) p3(dim_num)
  real ( kind = 8 ) pc1(dim_num)
  real ( kind = 8 ) pc2(dim_num)
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0156'
  write ( *, '(a)' ) '  CIRCLE_EXP2IMP_2D converts an explicit circle'
  write ( *, '(a)' ) '  to an implicit circle.'
  write ( *, '(a)' ) '  CIRCLE_IMP2EXP_2D converts an implicit circle'
  write ( *, '(a)' ) '  to an explicit circle.'

  pc1(1) = 10.0D+00
  pc1(2) =  5.0D+00
  r1 = 3.0D+00

  call circle_imp_print_2d ( r1, pc1, '  The implicit circle:' )

  call circle_imp2exp_2d ( r1, pc1, p1, p2, p3 )

  call r8vec_print ( dim_num, p1, '  P1:' )
  call r8vec_print ( dim_num, p2, '  P2:' )
  call r8vec_print ( dim_num, p3, '  P3:' )

  call circle_exp2imp_2d ( p1, p2, p3, r2, pc2 )

  call circle_imp_print_2d ( r2, pc2, '  The recovered implicit circle:' )

  return
end
subroutine test016 ( )

!*****************************************************************************80
!
!! TEST016 tests CIRCLE_IMP_POINTS_2D and POLYGON_AREA_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: p
  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) :: r = 2.0D+00
  real ( kind = 8 ) result

  pc(1:2) =  (/ 5.0D+00, -2.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST016'
  write ( *, '(a)' ) '  CIRCLE_IMP_POINTS_2D gets points on a circle;'
  write ( *, '(a)' ) '  POLYGON_AREA_2D finds the area of a polygon.'

  call circle_imp_print_2d ( r, pc, '  The implicit circle:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The area = ', pi * r * r

  n = 8
  allocate ( p(1:dim_num,1:n) )

  call circle_imp_points_2d ( r, pc, n, p )

  call r8mat_transpose_print ( dim_num, n, p, '  Sample results:' )

  deallocate ( p )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For any N, the sampled points define a polygon'
  write ( *, '(a)' ) '  whose area approximates the circle area.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  N      Area'
  write ( *, '(a)' ) ' '
 
  do n = 3, 24

    allocate ( p(1:dim_num,1:n) )

    call circle_imp_points_2d ( r, pc, n, p )
    call polygon_area_2d ( n, p, result )
    write ( *, '(2x,i8,2x,g14.6)' ) n, result

    deallocate ( p )

  end do
 
  return
end
subroutine test0165 ( )

!*****************************************************************************80
!
!! TEST0165 tests CIRCLE_IMP_POINTS_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: n = 12

  real ( kind = 8 ), dimension ( dim_num ) :: nc = (/ &
    1.0D+00, 1.0D+00, 1.0D+00 /)
  real ( kind = 8 ), dimension(dim_num,n) :: p
  real ( kind = 8 ), dimension ( dim_num ) :: pc = (/ &
    5.0D+00, -2.0D+00, 1.0D+00 /)
  real ( kind = 8 ) :: r = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0165'
  write ( *, '(a)' ) '  CIRCLE_IMP_POINTS_3D gets points on a circle in 3D;'

  call circle_imp_print_3d ( r, pc, nc, '  The implicit circle:' )

  call circle_imp_points_3d ( r, pc, nc, n, p )

  call r8mat_transpose_print ( dim_num, n, p, '  Points on the circle:' )
 
  return
end
subroutine test017 ( )

!*****************************************************************************80
!
!! TEST017 tests CIRCLE_IMP_POINTS_ARC_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 13
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) p(dim_num,n)
  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) :: r = 2.0D+00
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2

  pc(1:2) =  (/ 5.0D+00, -2.0D+00 /)
  theta1 = pi / 2.0D+00
  theta2 = 3.0D+00 * pi / 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST017'
  write ( *, '(a)' ) '  CIRCLE_IMP_POINTS_ARC_2D returns points on a'
  write ( *, '(a)' ) '  circular arc.'

  call circle_imp_print_2d ( r, pc, '  The implicit circle:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  The arc extends from THETA1 = ', theta1
  write ( *, '(a,2g14.6)' ) '  to THETA2 = ', theta2

  call circle_imp_points_arc_2d ( r, pc, theta1, theta2, n, p )

  call r8mat_transpose_print ( dim_num, n, p, '  Sample results:' )

  return
end
subroutine test018 ( )

!*****************************************************************************80
!
!! TEST018 tests CIRCLE_IMP_POINT_DIST_2D and CIRCLES_IMP_INT_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 5

  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  integer ( kind = 4 ) int_num
  real ( kind = 8 ), dimension(dim_num) :: pc1 = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) pc2(dim_num)
  real ( kind = 8 ), parameter, dimension (dim_num,test_num) :: &
    pc2_test = reshape ( (/ &
    5.0D+00,       5.0D+00, &
    7.0710678D+00, 7.0710678D+00, &
    4.0D+00,       0.0D+00, &
    6.0D+00,       0.0D+00, &
    0.0D+00,       0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) pint(dim_num,2)
  real ( kind = 8 ), parameter :: r1 = 5.0D+00
  real ( kind = 8 ) r2
  real ( kind = 8 ), parameter, dimension ( test_num ) :: r2_test = &
    (/ 0.5D+00, 5.0D+00, 3.0D+00, 3.0D+00, 5.0D+00 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST018'
  write ( *, '(a)' ) '  CIRCLE_IMP_POINT_DIST_2D finds the'
  write ( *, '(a)' ) '  distance from a point to a circle.'
  write ( *, '(a)' ) '  CIRCLES_IMP_INT_2D determines the intersections of'
  write ( *, '(a)' ) '  two circles in 2D.'

  call circle_imp_print_2d ( r1, pc1, '  The first circle:' )

  do test = 1, test_num

    r2 = r2_test(test)
    pc2(1:2) = pc2_test(1:2,test)

    call circle_imp_print_2d ( r2, pc2, '  The second circle:' )

    call circles_imp_int_2d ( r1, pc1, r2, pc2, int_num, pint )

    if ( int_num == 0 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The circles do not intersect.'

    else if ( int_num == 1 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The circles intersect at one point:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        P          Dist 1  Dist 2'
      write ( *, '(a)' ) ' '
      call circle_imp_point_dist_2d ( r1, pc1, pint(1:2,1), d1 )
      call circle_imp_point_dist_2d ( r2, pc2, pint(1:2,1), d2 )
      write ( *, '(2x,4f8.4)' ) pint(1:2,1), d1, d2

    else if ( int_num == 2 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The circles intersect at two points:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        P          Dist 1  Dist 2'
      write ( *, '(a)' ) ' '
      call circle_imp_point_dist_2d ( r1, pc1, pint(1:2,1), d1 )
      call circle_imp_point_dist_2d ( r2, pc2, pint(1:2,1), d2 )
      write ( *, '(2x,4f8.4)' ) pint(1:2,1), d1, d2
      call circle_imp_point_dist_2d ( r1, pc1, pint(1:2,2), d1 )
      call circle_imp_point_dist_2d ( r2, pc2, pint(1:2,2), d2 )
      write ( *, '(2x,4f8.4)' ) pint(1:2,2), d1, d2

    else if ( int_num == 3 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The circles coincide (infinite intersection).'

    end if

  end do

  return
end
subroutine test0183 ( )

!*****************************************************************************80
!
!! TEST0183 tests CIRCLE_LLR2IMP_2D.
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

  real ( kind = 8 ) :: d1
  real ( kind = 8 ) :: d2
  real ( kind = 8 ) :: d3
  real ( kind = 8 ) :: d4
  real ( kind = 8 ) :: p_hi =  10.0D+00
  real ( kind = 8 ) :: p_lo = -10.0D+00
  real ( kind = 8 ) pc(dim_num,4)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) q1(dim_num)
  real ( kind = 8 ) q2(dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ) r_hi
  real ( kind = 8 ) r_lo
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0183'
  write ( *, '(a)' ) '  CIRCLE_LLR2IMP_2D is given:'
  write ( *, '(a)' ) '  a line through P1 and P2,'
  write ( *, '(a)' ) '  a line through Q1 and Q2,'
  write ( *, '(a)' ) '  and a radius R,'
  write ( *, '(a)' ) '  and determines the centers C of 4 circles'
  write ( *, '(a)' ) '  of the given radius, tangent to both lines.'

  seed = 123456789

  do test = 1, test_num

    call r8vec_uniform_ab ( dim_num, p_lo, p_hi, seed, p1 )
    call r8vec_uniform_ab ( dim_num, p_lo, p_hi, seed, p2 )
    call r8vec_uniform_ab ( dim_num, p_lo, p_hi, seed, q1 )
    call r8vec_uniform_ab ( dim_num, p_lo, p_hi, seed, q2 )

    r_lo = 1.0D+00
    r_hi = 5.0D+00
    r = r8_uniform ( r_lo, r_hi, seed )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Radius R = ', r

    write ( *, '(a,2g14.6)' ) '  Point #P1: ', p1(1:dim_num)
    write ( *, '(a,2g14.6)' ) '  Point #P2: ', p2(1:dim_num)
    write ( *, '(a,2g14.6)' ) '  Point #Q1: ', q1(1:dim_num)
    write ( *, '(a,2g14.6)' ) '  Point #Q2: ', q2(1:dim_num)

    call circle_llr2imp_2d ( p1, p2, q1, q2, r, pc )

    write ( *, '(a,2g14.6)' ) '  Center #1: ', pc(1:dim_num,1)
    write ( *, '(a,2g14.6)' ) '  Center #2: ', pc(1:dim_num,2)
    write ( *, '(a,2g14.6)' ) '  Center #3: ', pc(1:dim_num,3)
    write ( *, '(a,2g14.6)' ) '  Center #4: ', pc(1:dim_num,4)

    call line_exp_point_dist_2d ( p1, p2, pc(1:dim_num,1), d1 )
    call line_exp_point_dist_2d ( p1, p2, pc(1:dim_num,2), d2 )
    call line_exp_point_dist_2d ( p1, p2, pc(1:dim_num,3), d3 )
    call line_exp_point_dist_2d ( p1, p2, pc(1:dim_num,4), d4 )

    write ( *, '(2x,4g14.6)' ) d1, d2, d3, d4

    call line_exp_point_dist_2d ( q1, q2, pc(1:dim_num,1), d1 )
    call line_exp_point_dist_2d ( q1, q2, pc(1:dim_num,2), d2 )
    call line_exp_point_dist_2d ( q1, q2, pc(1:dim_num,3), d3 )
    call line_exp_point_dist_2d ( q1, q2, pc(1:dim_num,4), d4 )

    write ( *, '(2x,4g14.6)' ) d1, d2, d3, d4

  end do

  return
end
subroutine test0185 ( )

!*****************************************************************************80
!
!! TEST0185 tests CIRCLE_PPPR2IMP_3D.
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

  real ( kind = 8 ) :: d11
  real ( kind = 8 ) :: d12
  real ( kind = 8 ) :: d21
  real ( kind = 8 ) :: d22
  real ( kind = 8 ) :: p_hi =  10.0D+00
  real ( kind = 8 ) :: p_lo = -10.0D+00
  real ( kind = 8 ) normal(dim_num)
  real ( kind = 8 ) pc(dim_num,2)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) p3(dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ) r_hi
  real ( kind = 8 ) r_lo
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0185'
  write ( *, '(a)' ) '  CIRCLE_PPPR2IMP_3D is given 3D points P1, P2, P3,'
  write ( *, '(a)' ) '  and a radius R,'
  write ( *, '(a)' ) '  and determines the centers C of two circles'
  write ( *, '(a)' ) '  of the given radius, passing through P1 and P2'
  write ( *, '(a)' ) '  and lying in the plane of P1, P2 and P3.'

  seed = 123456789

  do test = 1, test_num

    call r8vec_uniform_ab ( dim_num, p_lo, p_hi, seed, p1 )
    call r8vec_uniform_ab ( dim_num, p_lo, p_hi, seed, p2 )
    call r8vec_uniform_ab ( dim_num, p_lo, p_hi, seed, p3 )

    r_lo = sqrt ( sum ( ( p1(1:dim_num) - p2(1:dim_num) )**2 ) )
    r_hi = r_lo + 5.0D+00
    r = r8_uniform ( r_lo, r_hi, seed )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Radius R = ', r

    write ( *, '(a,3g14.6)' ) '  Point #1: ', p1(1:dim_num)
    write ( *, '(a,3g14.6)' ) '  Point #2: ', p2(1:dim_num)
    write ( *, '(a,3g14.6)' ) '  Point #3: ', p3(1:dim_num)

    call circle_pppr2imp_3d ( p1, p2, p3, r, pc, normal )

    write ( *, '(a,3g14.6)' ) '  Center #1: ', pc(1:dim_num,1)
    write ( *, '(a,3g14.6)' ) '  Center #2: ', pc(1:dim_num,2)
!
!  Check that the points are the right distance from the center.
!
    d11 = sqrt ( sum ( ( p1(1:dim_num) - pc(1:dim_num,1) )**2 ) )
    d21 = sqrt ( sum ( ( p2(1:dim_num) - pc(1:dim_num,1) )**2 ) )
    d12 = sqrt ( sum ( ( p1(1:dim_num) - pc(1:dim_num,2) )**2 ) )
    d22 = sqrt ( sum ( ( p2(1:dim_num) - pc(1:dim_num,2) )**2 ) )

    write ( *, '(2x,4g14.6)' ) d11, d21, d12, d22
!
!  Check that the radial vector to the point is perpendicular to NORMAL.
!
    d11 = dot_product ( normal(1:dim_num), p1(1:dim_num) - pc(1:dim_num,1) )
    d21 = dot_product ( normal(1:dim_num), p2(1:dim_num) - pc(1:dim_num,1) )
    d12 = dot_product ( normal(1:dim_num), p1(1:dim_num) - pc(1:dim_num,2) )
    d22 = dot_product ( normal(1:dim_num), p2(1:dim_num) - pc(1:dim_num,2) )

    write ( *, '(2x,4g14.6)' ) d11, d21, d12, d22

  end do

  return
end
subroutine test019 ( )

!*****************************************************************************80
!
!! TEST019 tests CIRCLE_PPR2IMP_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 November 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) :: d11
  real ( kind = 8 ) :: d12
  real ( kind = 8 ) :: d21
  real ( kind = 8 ) :: d22
  real ( kind = 8 ) :: p_hi =  10.0D+00
  real ( kind = 8 ) :: p_lo = -10.0D+00
  real ( kind = 8 ) pc(dim_num,2)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ) r_hi
  real ( kind = 8 ) r_lo
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST019'
  write ( *, '(a)' ) '  CIRCLE_PPR2IMP_2D is given 2D points P1 and P2,'
  write ( *, '(a)' ) '  and a radius R,'
  write ( *, '(a)' ) '  and determines the centers C of two circles'
  write ( *, '(a)' ) '  of the given radius, passing through P1 and P2.'

  seed = 123456789

  do test = 1, test_num

    call r8vec_uniform_ab ( dim_num, p_lo, p_hi, seed, p1 )
    call r8vec_uniform_ab ( dim_num, p_lo, p_hi, seed, p2 )

    r_lo = sqrt ( sum ( ( p1(1:dim_num) - p2(1:dim_num) )**2 ) )
    r_hi = r_lo + 5.0D+00
    r = r8_uniform ( r_lo, r_hi, seed );

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Radius R = ', r

    write ( *, '(a,2g14.6)' ) '  Point #1: ', p1(1:dim_num)
    write ( *, '(a,2g14.6)' ) '  Point #2: ', p2(1:dim_num)

    call circle_ppr2imp_2d ( p1, p2, r, pc )

    write ( *, '(a,2g14.6)' ) '  Center #1: ', pc(1:dim_num,1)
    write ( *, '(a,2g14.6)' ) '  Center #2: ', pc(1:dim_num,2)

    d11 = sqrt ( sum ( ( p1(1:dim_num) - pc(1:dim_num,1) )**2 ) )
    d21 = sqrt ( sum ( ( p2(1:dim_num) - pc(1:dim_num,1) )**2 ) )
    d12 = sqrt ( sum ( ( p1(1:dim_num) - pc(1:dim_num,2) )**2 ) )
    d22 = sqrt ( sum ( ( p2(1:dim_num) - pc(1:dim_num,2) )**2 ) )

    write ( *, '(2x,4g14.6)' ) d11, d21, d12, d22

  end do

  return
end
subroutine test020 ( )

!*****************************************************************************80
!
!! TEST020 tests CUBE_SIZE_3D and CUBE_SHAPE_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ), allocatable, dimension ( : ) ::  face_order
  integer ( kind = 4 ) face_order_max
  integer ( kind = 4 ), allocatable, dimension ( :, : ) ::  face_point
  integer ( kind = 4 ) point_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_coord

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST020'
  write ( *, '(a)' ) '  For the cube,'
  write ( *, '(a)' ) '  CUBE_SIZE_3D returns dimension information;'
  write ( *, '(a)' ) '  CUBE_SHAPE_3D returns face and order information.'
  write ( *, '(a)' ) '  SHAPE_PRINT_3D prints this information.'
!
!  Get the data sizes.
!
  call cube_size_3d ( point_num, edge_num, face_num, face_order_max )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of vertices = ', point_num
  write ( *, '(a,i8)' ) '  Number of edges =    ', edge_num
  write ( *, '(a,i8)' ) '  Number of faces =    ', face_num
  write ( *, '(a,i8)' ) '  Maximum face order = ', face_order_max
!
!  Make room for the data.
!
  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )
  allocate ( point_coord(1:3,1:point_num) )
!
!  Get the data.
!
  call cube_shape_3d ( point_num, face_num, face_order_max, point_coord, &
    face_order, face_point )
!
!  Print the data.
!
  call shape_print_3d ( point_num, face_num, face_order_max, &
    point_coord, face_order, face_point )

  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )

  return
end
subroutine test0201 ( )

!*****************************************************************************80
!
!! TEST0201 tests CYLINDER_POINT_DIST_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 6

  real ( kind = 8 ) :: dist
  real ( kind = 8 ), dimension ( test_num ) :: dist_test = (/ &
      3.0D+00, 0.5D+00, 5.0D+00, 8.0D+00, 1.0D+00, 0.25D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p_test = reshape ( (/ &
      4.0D+00,    0.5D+00,  0.0D+00, &
     -0.5D+00,   -1.0D+00,  0.0D+00, &
      4.0D+00,    6.0D+00,  0.0D+00, &
      0.75D+00, -10.0D+00,  0.0D+00, &
      0.0D+00,    0.0D+00,  0.0D+00, &
      0.25D+00,   1.75D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num ) :: p1 = (/ &
    0.0D+00, -2.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p2 = (/ &
    0.0D+00,  2.0D+00, 0.0D+00 /)
  real ( kind = 8 ) :: r = 1.0D+00
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0201'
  write ( *, '(a)' ) '  CYLINDER_POINT_DIST_3D computes the distance'
  write ( *, '(a)' ) '  to a cylinder.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)'  ) '  Radius R = ', r
  write ( *, '(a,3g14.6)' ) '  Center of bottom disk = ', p1(1:dim_num)
  write ( *, '(a,3g14.6)' ) '  Center of top disk =    ', p2(1:dim_num)

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,3g14.6)' ) '  P =             ', p(1:dim_num)

    call cylinder_point_dist_3d ( p1, p2, r, p, dist )

    write ( *, '(a,g14.6)' ) '  DIST (computed) = ', dist
    write ( *, '(a,g14.6)' ) '  DIST (exact) =    ', dist_test(test)

  end do

  return
end
subroutine test02015 ( )

!*****************************************************************************80
!
!! TEST02015 tests CYLINDER_POINT_DIST_SIGNED_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 6

  real ( kind = 8 ) :: dist
  real ( kind = 8 ), dimension ( test_num ) :: dist_test = (/ &
      3.0D+00, -0.5D+00, 5.0D+00, 8.0D+00, -1.0D+00, -0.25D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p_test = reshape ( (/ &
      4.0D+00,    0.5D+00,  0.0D+00, &
     -0.5D+00,   -1.0D+00,  0.0D+00, &
      4.0D+00,    6.0D+00,  0.0D+00, &
      0.75D+00, -10.0D+00,  0.0D+00, &
      0.0D+00,    0.0D+00,  0.0D+00, &
      0.25D+00,   1.75D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num ) :: p1 = (/ &
    0.0D+00, -2.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p2 = (/ &
    0.0D+00,  2.0D+00, 0.0D+00 /)
  real ( kind = 8 ) :: r = 1.0D+00
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02015'
  write ( *, '(a)' ) '  CYLINDER_POINT_DIST_SIGNED_3D computes the signed'
  write ( *, '(a)' ) '  distance to a cylinder.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)'  ) '  Radius R = ', r
  write ( *, '(a,3g14.6)' ) '  Center of bottom disk = ', p1(1:dim_num)
  write ( *, '(a,3g14.6)' ) '  Center of top disk =    ', p2(1:dim_num)

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,3g14.6)' ) '  P =             ', p(1:dim_num)

    call cylinder_point_dist_signed_3d ( p1, p2, r, p, dist )

    write ( *, '(a,g14.6)' ) '  Signed distance (computed) = ', dist
    write ( *, '(a,g14.6)' ) '  Signed distance (exact) =    ', dist_test(test)

  end do

  return
end
subroutine test0202 ( )

!*****************************************************************************80
!
!! TEST0202 tests CYLINDER_POINT_INSIDE_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 6

  logical :: inside
  logical, dimension ( test_num ) :: inside_test = (/ &
      .false., .true., .false., .false., .true., .true. /)
  real ( kind = 8 ), dimension ( dim_num ) :: p
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p_test = reshape ( (/ &
      4.0D+00,    0.5D+00,  0.0D+00, &
     -0.5D+00,   -1.0D+00,  0.0D+00, &
      4.0D+00,    6.0D+00,  0.0D+00, &
      0.75D+00, -10.0D+00,  0.0D+00, &
      0.0D+00,    0.0D+00,  0.0D+00, &
      0.25D+00,   1.75D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num ) :: p1 = (/ &
    0.0D+00, -2.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p2 = (/ &
    0.0D+00,  2.0D+00, 0.0D+00 /)
  real ( kind = 8 ) :: r = 1.0D+00
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0202'
  write ( *, '(a)' ) '  CYLINDER_POINT_INSIDE_3D determines if a point is'
  write ( *, '(a)' ) '  inside a cylinder.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)'  ) '  Radius R = ', r
  write ( *, '(a,3g14.6)' ) '  Center of bottom disk = ', p1(1:dim_num)
  write ( *, '(a,3g14.6)' ) '  Center of top disk =    ', p2(1:dim_num)

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,3g14.6)' ) '  P =             ', p(1:dim_num)

    call cylinder_point_inside_3d ( p1, p2, r, p, inside )

    write ( *, '(a,l1)' ) '  INSIDE (computed) = ', inside
    write ( *, '(a,l1)' ) '  INSIDE (exact) =    ', inside_test(test)

  end do

  return
end
subroutine test0203 ( )

!*****************************************************************************80
!
!! TEST0203 tests CYLINDER_POINT_NEAR_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 6

  real ( kind = 8 ), dimension ( dim_num ) :: p
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p_test = reshape ( (/ &
      4.0D+00,    0.5D+00,  0.0D+00, &
     -0.5D+00,   -1.0D+00,  0.0D+00, &
      4.0D+00,    6.0D+00,  0.0D+00, &
      0.75D+00, -10.0D+00,  0.0D+00, &
      0.0D+00,    0.0D+00,  0.0D+00, &
      0.25D+00,   1.75D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num ) :: p1 = (/ &
    0.0D+00, -2.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p2 = (/ &
    0.0D+00,  2.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: pn
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: pn_test = reshape ( (/ &
      1.0D+00,   0.5D+00,  0.0D+00, &
     -1.0D+00,  -1.0D+00,  0.0D+00, &
      1.0D+00,   2.0D+00,  0.0D+00, &
      0.75D+00, -2.0D+00,  0.0D+00, &
      1.0D+00,   0.0D+00,  0.0D+00, &
      0.25D+00,  2.0D+00,  0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) :: r = 1.0D+00
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0203'
  write ( *, '(a)' ) '  CYLINDER_POINT_NEAR_3D computes the nearest point'
  write ( *, '(a)' ) '  on a cylinder.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)'  ) '  Radius R = ', r
  write ( *, '(a,3g14.6)' ) '  Center of bottom disk = ', p1(1:dim_num)
  write ( *, '(a,3g14.6)' ) '  Center of top disk =    ', p2(1:dim_num)

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,3g14.6)' ) '  P =             ', p(1:dim_num)

    call cylinder_point_near_3d ( p1, p2, r, p, pn )

    write ( *, '(a,3g14.6)' ) '  PN (computed) = ', pn(1:dim_num)
    write ( *, '(a,3g14.6)' ) '  PN (exact) =    ', pn_test(1:dim_num,test)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (Note that case 5 is ambiguous.  The set of nearest'
  write ( *, '(a)' ) '  points forms a circle, any of which will do.)'

  return
end
subroutine test02035 ( )

!*****************************************************************************80
!
!! TEST02035 tests CYLINDER_SAMPLE_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: n = 20

  real ( kind = 8 ), dimension ( dim_num, n ) :: p
  real ( kind = 8 ), dimension ( dim_num ) :: p1 = (/ &
    0.0D+00, -2.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p2 = (/ &
    0.0D+00,  2.0D+00, 0.0D+00 /)
  real ( kind = 8 ) :: r = 1.0D+00
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02035'
  write ( *, '(a)' ) '  CYLINDER_SAMPLE_3D samples points in a cylinder.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)'  ) '  Radius R = ', r
  write ( *, '(a,3g14.6)' ) '  Center of bottom disk = ', p1(1:dim_num)
  write ( *, '(a,3g14.6)' ) '  Center of top disk =    ', p2(1:dim_num)

  call cylinder_sample_3d ( p1, p2, r, n, seed, p )

  call r8mat_transpose_print ( dim_num, n, p, '  Sample points:' )

  return
end
subroutine test0204 ( )

!*****************************************************************************80
!
!! TEST0204 tests CYLINDER_VOLUME_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) r8_pi
  real ( kind = 8 ), dimension ( dim_num ) :: p1 = (/ &
    1.0D+00, 2.0D+00, 3.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p2 = (/ &
    5.0D+00, 6.0D+00, 5.0D+00 /)
  real ( kind = 8 ) :: r = 5.0D+00
  real ( kind = 8 ) volume

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0204'
  write ( *, '(a)' ) '  CYLINDER_VOLUME_3D computes the volume of a cylinder.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)'  ) '  Radius R = ', r
  write ( *, '(a,3g14.6)' ) '  Center of bottom disk = ', p1(1:dim_num)
  write ( *, '(a,3g14.6)' ) '  Center of top disk =    ', p2(1:dim_num)

  call cylinder_volume_3d ( p1, p2, r, volume )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Volume (computed) = ', volume
  write ( *, '(a,g14.6)' ) '  Volume (exact)    = ', r8_pi ( ) * 150.0D+00

  return
end
subroutine test0205 ( )

!*****************************************************************************80
!
!! TEST0205 tests DEGREES_TO_RADIANS and RADIANS_TO_DEGREES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) angle_deg
  real ( kind = 8 ) angle_deg2
  real ( kind = 8 ) angle_rad
  real ( kind = 8 ) degrees_to_radians
  integer ( kind = 4 ) i
  real ( kind = 8 ) radians_to_degrees

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0205'
  write ( *, '(a)' ) '  DEGREES_TO_RADIANS converts an angle from degrees'
  write ( *, '(a)' ) '  to radians;'
  write ( *, '(a)' ) '  RADIANS_TO_DEGREES converts an angle from radians'
  write ( *, '(a)' ) '  to degrees;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Degrees     Radians     Degrees'
  write ( *, '(a)' ) ' '

  do i = -2, 14
    angle_deg = real ( 30 * i, kind = 8 )
    angle_rad = degrees_to_radians ( angle_deg )
    angle_deg2 = radians_to_degrees ( angle_rad )
    write ( *, '(2x,f10.6,2x,f10.6,2x,f10.6)' ) angle_deg, angle_rad, angle_deg2
  end do

  return
end
subroutine test021 ( )

!*****************************************************************************80
!
!! TEST021 tests DIRECTION_PERT_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 3

  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) sigma(test_num)
  integer ( kind = 4 ) test
  real ( kind = 8 ) vbase(dim_num)
  real ( kind = 8 ) vran(dim_num)

  vbase(1:dim_num) = (/ 1.0D+00,  0.0D+00, 0.0D+00 /)
  sigma(1:test_num) = (/ 0.99D+00, 0.5D+00, 0.1D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST021'
  write ( *, '(a)' ) '  DIRECTION_PERT_3D perturbs a direction vector.'

  call r8vec_print ( dim_num, vbase, '  The base vector:' )

  do test = 1, test_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Using Sigma = ', sigma(test)
    write ( *, '(a)' ) ' '

    do i = 1, 20
      call direction_pert_3d ( sigma(test), vbase, seed, vran )
      write ( *, '(2x,3f8.4)' ) vran(1:dim_num)
    end do

  end do

  return
end
subroutine test023 ( )

!*****************************************************************************80
!
!! TEST023 tests R8VEC_UNIFORM_UNIT;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 4
  integer ( kind = 4 ), parameter :: test_num = 10

  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  real ( kind = 8 ) vran(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST023'
  write ( *, '(a)' ) '  R8VEC_UNIFORM_UNIT picks a random direction vector.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    call r8vec_uniform_unit ( dim_num, seed, vran )
    write ( *, '(2x,4f8.4)' ) vran(1:dim_num)
  end do

  return
end
subroutine test0232 ( )

!*****************************************************************************80
!
!! TEST0232 tests DISK_POINT_DIST_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 5

  real ( kind = 8 ), dimension(dim_num) :: axis = (/ &
    0.0D+00, 1.0D+00, 1.0D+00 /)
  real ( kind = 8 ) dist
  real ( kind = 8 ), dimension(test_num) :: dist_test = (/ &
    2.0D+00, 0.0D+00, 0.0D+00, 8.0D+00, 10.0D+00 /)
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: pc = (/ &
    0.0D+00, 1.4142135D+00, 1.4142135D+00 /)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = &
  reshape ( (/ &
     0.0D+00,  0.0D+00,        0.0D+00, &
     0.0D+00,  0.70710677D+00, 2.1213202D+00, &
     2.0D+00,  1.4142135D+00,  1.4142135D+00, &
    10.0D+00,  1.4142135D+00,  1.4142135D+00, &
    10.0D+00,  5.6568542D+00,  5.6568542D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) :: r = 2.0D+00
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0232'
  write ( *, '(a)' ) '  DISK_POINT_DIST_3D finds the distance from'
  write ( *, '(a)' ) '  a disk to a point in 3D.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Disk radius = ', r
  call r8vec_print ( dim_num, pc, '  Disk center: ' )
  call r8vec_print ( dim_num, axis, '  Disk axis: ' )

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call r8vec_print ( dim_num, p, '  Point: ' )

    call disk_point_dist_3d ( pc, r, axis, p, dist )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6,a,g14.6)' ) '  Distance = ', dist, &
      '  Expected = ', dist_test(test)

  end do

  return
end
subroutine test0234 ( )

!*****************************************************************************80
!
!! TEST0234 tests R8MAT_SOLVE_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2005
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: n = 2

  real ( kind = 8 ), dimension (n,n) :: a
  real ( kind = 8 ), dimension ( n ) :: b
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5
  real ( kind = 8 ), dimension ( n ) :: x
  real ( kind = 8 ), dimension ( n ) :: x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0234'
  write ( *, '(a)' ) '  R8MAT_SOLVE_2D solves 2D linear systems.'

  seed = 123456789

  do test = 1, test_num

    call r8mat_uniform_01 ( n, n, seed, a )
    call r8vec_uniform_01 ( n, seed, x )
    b(1:n) = matmul ( a(1:n,1:n), x(1:n) )

    call r8mat_solve_2d ( a, b, det, x2 )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Solution / Computed:'
    write ( *, '(a)' ) ' '

    do i = 1, n
      write ( *, '(2x,g14.6,2x,g14.6)' ) x(i), x2(i)
    end do

  end do

  return
end
subroutine test0235 ( )

!*****************************************************************************80
!
!! TEST0235 tests DMS_TO_RADIANS and RADIANS_TO_DMS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) angle_deg
  integer ( kind = 4 ) angle_min
  real ( kind = 8 ) angle_rad
  real ( kind = 8 ) angle_rad2
  integer ( kind = 4 ) angle_sec
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0235'
  write ( *, '(a)' ) '  DMS_TO_RADIANS converts an angle from '
  write ( *, '(a)' ) '  degrees/minutes/seconds to radians;'
  write ( *, '(a)' ) '  RADIANS_TO_DEGREES converts an angle from radians'
  write ( *, '(a)' ) '  to degrees/minutes/seconds;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Radians     DMS     Radians'
  write ( *, '(a)' ) ' '

  do i = -2, 15
    angle_rad = pi * real ( i, kind = 8 ) / 7.0D+00
    call radians_to_dms ( angle_rad, angle_deg, angle_min, angle_sec )
    call dms_to_radians ( angle_deg, angle_min, angle_sec, angle_rad2 )
    write ( *, '(2x,f10.6,2x,i4,2x,i3,2x,i3,2x,f10.6)' ) &
      angle_rad, angle_deg, angle_min, angle_sec, angle_rad2
  end do

  return
end
subroutine test0236 ( )

!*****************************************************************************80
!
!! TEST0236 tests DODEC_SIZE_3D and DODEC_SHAPE_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ), allocatable, dimension ( : ) ::  face_order
  integer ( kind = 4 ) face_order_max
  integer ( kind = 4 ), allocatable, dimension ( :, : ) ::  face_point
  integer ( kind = 4 ) point_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_coord

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0236'
  write ( *, '(a)' ) '  For the dodecahedron,'
  write ( *, '(a)' ) '  DODEC_SIZE_3D returns dimension information;'
  write ( *, '(a)' ) '  DODEC_SHAPE_3D returns face and order information.'
  write ( *, '(a)' ) '  SHAPE_PRINT_3D prints this information.'

  call dodec_size_3d ( point_num, edge_num, face_num, face_order_max )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of vertices = ', point_num
  write ( *, '(a,i8)' ) '  Number of edges =    ', edge_num
  write ( *, '(a,i8)' ) '  Number of faces =    ', face_num
  write ( *, '(a,i8)' ) '  Maximum face order = ', face_order_max

  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )
  allocate ( point_coord(1:3,1:point_num) )

  call dodec_shape_3d ( point_num, face_num, face_order_max, point_coord, &
    face_order, face_point )

  call shape_print_3d ( point_num, face_num, face_order_max, &
    point_coord, face_order, face_point )

  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )

  return
end
subroutine test0238 ( )

!*****************************************************************************80
!
!! TEST0238 tests DUAL_SIZE_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) edge_num1
  integer ( kind = 4 ) edge_num2
  integer ( kind = 4 ), allocatable, dimension ( :, : ) ::  edge_point1
  integer ( kind = 4 ) face_num1
  integer ( kind = 4 ) face_num2
  integer ( kind = 4 ), allocatable, dimension ( : ) ::  face_order1
  integer ( kind = 4 ) face_order_max1
  integer ( kind = 4 ) face_order_max2
  integer ( kind = 4 ), allocatable, dimension ( :, : ) ::  face_point1
  integer ( kind = 4 ) point_num1
  integer ( kind = 4 ) point_num2
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_coord1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0238'
  write ( *, '(a)' ) '  DUAL_SIZE_3D finds the "sizes" of the dual of a'
  write ( *, '(a)' ) '  polyhedron;'
!
!  Get the CUBE shape.
!
  call cube_size_3d ( point_num1, edge_num1, face_num1, face_order_max1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The cube:'
  write ( *, '(a,i8)' ) '    Number of vertices: ', point_num1
  write ( *, '(a,i8)' ) '    Number of edges:    ', edge_num1
  write ( *, '(a,i8)' ) '    Number of faces:    ', face_num1
  write ( *, '(a,i8)' ) '    Maximum face order: ', face_order_max1

  allocate ( face_order1(1:face_num1) )
  allocate ( face_point1(1:face_order_max1,1:face_num1) )
  allocate ( point_coord1(1:dim_num,1:point_num1) )

  call cube_shape_3d ( point_num1, face_num1, face_order_max1, point_coord1, &
    face_order1, face_point1 )

  call dual_size_3d ( point_num1, edge_num1, face_num1, face_order_max1, &
    point_coord1, face_order1, face_point1, point_num2, edge_num2, &
    face_num2, face_order_max2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The dual of the cube:'
  write ( *, '(a,i8)' ) '    Number of vertices: ', point_num2
  write ( *, '(a,i8)' ) '    Number of edges:    ', edge_num2
  write ( *, '(a,i8)' ) '    Number of faces:    ', face_num2
  write ( *, '(a,i8)' ) '    Maximum face order: ', face_order_max2

  deallocate ( face_order1 )
  deallocate ( face_point1 )
  deallocate ( point_coord1 )
!
!  Get the DODECAHEDRON shape.
!
  call dodec_size_3d ( point_num1, edge_num1, face_num1, face_order_max1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The dodecahedron:'
  write ( *, '(a,i8)' ) '    Number of vertices: ', point_num1
  write ( *, '(a,i8)' ) '    Number of edges:    ', edge_num1
  write ( *, '(a,i8)' ) '    Number of faces:    ', face_num1
  write ( *, '(a,i8)' ) '    Maximum face order: ', face_order_max1

  allocate ( face_order1(1:face_num1) )
  allocate ( face_point1(1:face_order_max1,1:face_num1) )
  allocate ( point_coord1(1:dim_num,1:point_num1) )

  call dodec_shape_3d ( point_num1, face_num1, face_order_max1, &
    point_coord1, face_order1, face_point1 )

  call dual_size_3d ( point_num1, edge_num1, face_num1, face_order_max1, &
    point_coord1, face_order1, face_point1, point_num2, edge_num2, &
    face_num2, face_order_max2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The dual of the dodecahedron:'
  write ( *, '(a,i8)' ) '    Number of vertices: ', point_num2
  write ( *, '(a,i8)' ) '    Number of edges:    ', edge_num2
  write ( *, '(a,i8)' ) '    Number of faces:    ', face_num2
  write ( *, '(a,i8)' ) '    Maximum face order: ', face_order_max2

  deallocate ( face_order1 )
  deallocate ( face_point1 )
  deallocate ( point_coord1 )
!
!  Get the ICOSAHEDRON shape.
!
  call icos_size ( point_num1, edge_num1, face_num1, face_order_max1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The icosahedron:'
  write ( *, '(a,i8)' ) '    Number of vertices: ', point_num1
  write ( *, '(a,i8)' ) '    Number of edges:    ', edge_num1
  write ( *, '(a,i8)' ) '    Number of faces:    ', face_num1
  write ( *, '(a,i8)' ) '    Maximum face order: ', face_order_max1

  allocate ( edge_point1(1:2,1:edge_num1) )
  allocate ( face_order1(1:face_num1) )
  allocate ( face_point1(1:face_order_max1,1:face_num1) )
  allocate ( point_coord1(1:dim_num,1:point_num1) )

  call icos_shape ( point_num1, edge_num1, face_num1, face_order_max1, &
    point_coord1, edge_point1, face_order1, face_point1 )

  call dual_size_3d ( point_num1, edge_num1, face_num1, face_order_max1, &
    point_coord1, face_order1, face_point1, point_num2, edge_num2, &
    face_num2, face_order_max2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The dual of the icosahedron:'
  write ( *, '(a,i8)' ) '    Number of vertices: ', point_num2
  write ( *, '(a,i8)' ) '    Number of edges:    ', edge_num2
  write ( *, '(a,i8)' ) '    Number of faces:    ', face_num2
  write ( *, '(a,i8)' ) '    Maximum face order: ', face_order_max2

  deallocate ( edge_point1 )
  deallocate ( face_order1 )
  deallocate ( face_point1 )
  deallocate ( point_coord1 )
!
!  Get the OCTAHEDRON shape.
!
  call octahedron_size_3d ( point_num1, edge_num1, face_num1, face_order_max1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The octahedron:'
  write ( *, '(a,i8)' ) '    Number of vertices: ', point_num1
  write ( *, '(a,i8)' ) '    Number of edges:    ', edge_num1
  write ( *, '(a,i8)' ) '    Number of faces:    ', face_num1
  write ( *, '(a,i8)' ) '    Maximum face order: ', face_order_max1

  allocate ( face_order1(1:face_num1) )
  allocate ( face_point1(1:face_order_max1,1:face_num1) )
  allocate ( point_coord1(1:dim_num,1:point_num1) )

  call octahedron_shape_3d ( point_num1, face_num1, face_order_max1, &
    point_coord1, face_order1, face_point1 )

  call dual_size_3d ( point_num1, edge_num1, face_num1, face_order_max1, &
    point_coord1, face_order1, face_point1, point_num2, edge_num2, &
    face_num2, face_order_max2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The dual of the octahedron:'
  write ( *, '(a,i8)' ) '    Number of vertices: ', point_num2
  write ( *, '(a,i8)' ) '    Number of edges:    ', edge_num2
  write ( *, '(a,i8)' ) '    Number of faces:    ', face_num2
  write ( *, '(a,i8)' ) '    Maximum face order: ', face_order_max2

  deallocate ( face_order1 )
  deallocate ( face_point1 )
  deallocate ( point_coord1 )
!
!  Get the SOCCER BALL shape.
!
  call soccer_size_3d ( point_num1, edge_num1, face_num1, face_order_max1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The soccer ball:'
  write ( *, '(a,i8)' ) '    Number of vertices: ', point_num1
  write ( *, '(a,i8)' ) '    Number of edges:    ', edge_num1
  write ( *, '(a,i8)' ) '    Number of faces:    ', face_num1
  write ( *, '(a,i8)' ) '    Maximum face order: ', face_order_max1

  allocate ( face_order1(1:face_num1) )
  allocate ( face_point1(1:face_order_max1,1:face_num1) )
  allocate ( point_coord1(1:dim_num,1:point_num1) )

  call soccer_shape_3d ( point_num1, face_num1, face_order_max1, &
    point_coord1, face_order1, face_point1 )

  call dual_size_3d ( point_num1, edge_num1, face_num1, face_order_max1, & 
    point_coord1, face_order1, face_point1, point_num2, edge_num2, face_num2, &
    face_order_max2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The dual of the "soccer ball":'
  write ( *, '(a,i8)' ) '    Number of vertices: ', point_num2
  write ( *, '(a,i8)' ) '    Number of edges:    ', edge_num2
  write ( *, '(a,i8)' ) '    Number of faces:    ', face_num2
  write ( *, '(a,i8)' ) '    Maximum face order: ', face_order_max2

  deallocate ( face_order1 )
  deallocate ( face_point1 )
  deallocate ( point_coord1 )
!
!  Get the TETRAHEDRON shape.
!
  call tetrahedron_size_3d ( point_num1, edge_num1, face_num1, face_order_max1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The tetrahedron:'
  write ( *, '(a,i8)' ) '    Number of vertices: ', point_num1
  write ( *, '(a,i8)' ) '    Number of edges:    ', edge_num1
  write ( *, '(a,i8)' ) '    Number of faces:    ', face_num1
  write ( *, '(a,i8)' ) '    Maximum face order: ', face_order_max1

  allocate ( face_order1(1:face_num1) )
  allocate ( face_point1(1:face_order_max1,1:face_num1) )
  allocate ( point_coord1(1:dim_num,1:point_num1) )

  call tetrahedron_shape_3d ( point_num1, face_num1, face_order_max1, &
    point_coord1, face_order1, face_point1 )

  call dual_size_3d ( point_num1, edge_num1, face_num1, face_order_max1, &
    point_coord1, face_order1, face_point1, point_num2, edge_num2, &
    face_num2, face_order_max2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The dual of the tetrahedron:'
  write ( *, '(a,i8)' ) '    Number of vertices: ', point_num2
  write ( *, '(a,i8)' ) '    Number of edges:    ', edge_num2
  write ( *, '(a,i8)' ) '    Number of faces:    ', face_num2
  write ( *, '(a,i8)' ) '    Maximum face order: ', face_order_max2

  deallocate ( face_order1 )
  deallocate ( face_point1 )
  deallocate ( point_coord1 )

  return
end
subroutine test024 ( )

!*****************************************************************************80
!
!! TEST024 tests DUAL_SHAPE_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) edge_num1
  integer ( kind = 4 ) edge_num2
  integer ( kind = 4 ) face_num1
  integer ( kind = 4 ) face_num2
  integer ( kind = 4 ), allocatable, dimension ( : ) ::  face_order1
  integer ( kind = 4 ), allocatable, dimension ( : ) ::  face_order2
  integer ( kind = 4 ) face_order_max1
  integer ( kind = 4 ) face_order_max2
  integer ( kind = 4 ), allocatable, dimension ( :, : ) ::  face_point1
  integer ( kind = 4 ), allocatable, dimension ( :, : ) ::  face_point2
  integer ( kind = 4 ) point_num1
  integer ( kind = 4 ) point_num2
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_coord1
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_coord2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST024'
  write ( *, '(a)' ) '  DUAL_SHAPE_3D finds the dual of a polyhedron.'
!
!  Get the dodecahedron shape.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The dodecahedron:'

  call dodec_size_3d ( point_num1, edge_num1, face_num1, face_order_max1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of vertices = ', point_num1
  write ( *, '(a,i8)' ) '  Number of edges =    ', edge_num1
  write ( *, '(a,i8)' ) '  Number of faces =    ', face_num1
  write ( *, '(a,i8)' ) '  Maximum face order = ', face_order_max1

  allocate ( face_order1(1:face_num1) )
  allocate ( face_point1(1:face_order_max1,1:face_num1) )
  allocate ( point_coord1(1:dim_num,1:point_num1) )

  call dodec_shape_3d ( point_num1, face_num1, face_order_max1, point_coord1, &
    face_order1, face_point1 )
!
!  Get the dual.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The dual of the dodecahedron:'

  call dual_size_3d ( point_num1, edge_num1, face_num1, face_order_max1, &
    point_coord1, face_order1, face_point1, point_num2, edge_num2, &
    face_num2, face_order_max2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of vertices = ', point_num2
  write ( *, '(a,i8)' ) '  Number of edges =    ', edge_num2
  write ( *, '(a,i8)' ) '  Number of faces =    ', face_num2
  write ( *, '(a,i8)' ) '  Maximum face order = ', face_order_max2

  allocate ( face_order2(1:face_num2) )
  allocate ( face_point2(1:face_order_max2,1:face_num2) )
  allocate ( point_coord2(1:dim_num,1:point_num2) )

  call dual_shape_3d ( point_num1, face_num1, face_order_max1, point_coord1, &
    face_order1, face_point1, point_num2, face_num2, face_order_max2, &
    point_coord2, face_order2, face_point2 )

  call shape_print_3d ( point_num2, face_num2, face_order_max2, &
    point_coord2, face_order2, face_point2 )

  deallocate ( face_order1 )
  deallocate ( face_order2 )
  deallocate ( face_point1 )
  deallocate ( face_point2 )
  deallocate ( point_coord1 )
  deallocate ( point_coord2 )

  return
end
subroutine test0243 ( )

!*****************************************************************************80
!
!! TEST0243 tests R8VEC_ANY_NORMAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 10
  integer ( kind = 4 ), parameter :: test_num = 5

  real ( kind = 8 ) r8vec_norm
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  real ( kind = 8 ) v1(dim_num)
  real ( kind = 8 ) v1_length
  real ( kind = 8 ) v1v2_dot
  real ( kind = 8 ) v2(dim_num)
  real ( kind = 8 ) v2_length

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0243'
  write ( *, '(a)' ) '  R8VEC_ANY_NORMAL computes a vector V2 that is normal'
  write ( *, '(a)' ) '  to a given vector V1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Test    ||V1||      ||V2||        V1.V2'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do test = 1, test_num

    call r8vec_uniform_01 ( dim_num, seed, v1 )
    v1_length = r8vec_norm ( dim_num, v1 )
    call r8vec_any_normal ( dim_num, v1, v2 )
    v2_length = r8vec_norm ( dim_num, v2 )
    v1v2_dot = dot_product ( v1(1:dim_num), v2(1:dim_num) )
    write ( *, '(2x,i8,2x,f10.6,2x,f10.6,2x,f10.6)' ) &
      test, v1_length, v2_length, v1v2_dot
  end do

  return
end
subroutine test0245 ( )

!*****************************************************************************80
!
!! TEST0245 tests R8VEC_NORMAL_01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 1000

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n_max)
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_mean
  real ( kind = 8 ) x_min
  real ( kind = 8 ) x_var

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0245'
  write ( *, '(a)' ) '  R8VEC_NORMAL_01 computes a vector of normally'
  write ( *, '(a)' ) '  distributed random numbers.'
  write ( *, '(a,i12)' ) '  Using initial random number seed = ', seed
!
!  Test 1:
!  Simply call 5 times for 1 value, and print.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test #1: Call 5 times, 1 value each time.'
  write ( *, '(a)' ) ' '

  n = 1
  do i = 1, 5
    call r8vec_normal_01 ( n, seed, x )
    write ( *, '(2x,i8,g14.6)' ) i, x(1)
  end do
!
!  Test 2:
!  Restore the random number seed, and repeat.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test #2: Restore the random number seed.'
  write ( *, '(a)' ) '  Call 5 times, 1 value each time.'
  write ( *, '(a)' ) '  The results should be identical.'
  write ( *, '(a)' ) ' '

  n = -1
  call r8vec_normal_01 ( n, seed, x )

  seed = 123456789

  n = 1
  do i = 1, 5
    call r8vec_normal_01 ( n, seed, x )
    write ( *, '(2x,i8,g14.6)' ) i, x(1)
  end do
!
!  Test 3:
!  Restore the random number seed, compute all 5 values at once.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test #3: Restore the random number seed.'
  write ( *, '(a)' ) '  Call 1 time for 5 values.'
  write ( *, '(a)' ) '  The results should be identical.'
  write ( *, '(a)' ) ' '

  n = -1
  call r8vec_normal_01 ( n, seed, x )

  seed = 123456789

  n = 5
  call r8vec_normal_01 ( n, seed, x )

  do i = 1, n
    write ( *, '(2x,i8,g14.6)' ) i, x(i)
  end do
!
!  Test 4:
!  Restore the random number seed, compute all 5 values at once.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test #4: Restore the random number seed.'
  write ( *, '(a)' ) '  Call for 2, 1, and 2 values.'
  write ( *, '(a)' ) '  The results should be identical.'
  write ( *, '(a)' ) ' '

  n = -1
  call r8vec_normal_01 ( n, seed, x )

  seed = 123456789

  n = 2
  call r8vec_normal_01 ( n, seed, x )

  do i = 1, n
    write ( *, '(2x,i8,g14.6)' ) i, x(i)
  end do

  n = 1
  call r8vec_normal_01 ( n, seed, x )

  do i = 1, n
    write ( *, '(2x,i8,g14.6)' ) i, x(i)
  end do

  n = 2
  call r8vec_normal_01 ( n, seed, x )

  do i = 1, n
    write ( *, '(2x,i8,g14.6)' ) i, x(i)
  end do
!
!  Test 5:
!  Determine the minimum, maximum, mean and variance.
!
  n = n_max
  call r8vec_normal_01 ( n, seed, x )
  x_min = minval ( x(1:n) )
  x_max = maxval ( x(1:n) )
  x_mean = sum ( x(1:n) ) / real ( n, kind = 8 )
  x_var = sum ( ( x(1:n) - x_mean )**2 ) / real ( n - 1, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test #5:'
  write ( *, '(a,i12)' ) '  Number of samples was ', n
  write ( *, '(a,g14.6)' ) '  Minimum value was ', x_min
  write ( *, '(a,g14.6)' ) '  Maximum value was ', x_max
  write ( *, '(a,g14.6)' ) '  Average value was ', x_mean
  write ( *, '(a,g14.6)' ) '  Variance was      ', x_var
  write ( *, '(a,g14.6)' ) '  Expected average  ', 0.0D+00
  write ( *, '(a,g14.6)' ) '  Expected variance ', 1.0D+00

  return
end
subroutine test025 ( )

!*****************************************************************************80
!
!! TEST025 tests ELLIPSE_POINT_DIST_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) dist
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: n = 10
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) :: r1 = 3.0D+00
  real ( kind = 8 ) :: r2 = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST025:'
  write ( *, '(a)' ) '  ELLIPSE_POINT_DIST_2D is given a point P, and'
  write ( *, '(a)' ) '  finds the distance to an ellipse in 2D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The ellipse is (X/R1)^2 + (Y/R2)^2 = 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f14.6)' ) '  R1 = ', r1
  write ( *, '(a,f14.6)' ) '  R2 = ', r2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           P            DIST'
  write ( *, '(a)' ) ' '

  do i = -3, n + 3

    p(1) = ( real ( n - i, kind = 8 ) * 0.0D+00   &
           + real (     i, kind = 8 ) * 4.0D+00 ) &
           / real ( n,     kind = 8 )

    p(2) = ( real ( n - i, kind = 8 ) * 3.0D+00   &
           + real (     i, kind = 8 ) * 0.0D+00 ) &
           / real ( n,     kind = 8 )

    call ellipse_point_dist_2d ( r1, r2, p, dist )

    write ( *, '(2x,2f8.4,2x,f8.4)' ) p(1:dim_num), dist

  end do

  return
end
subroutine test0255 ( )

!*****************************************************************************80
!
!! TEST0255 tests ELLIPSE_POINT_NEAR_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: n = 10
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) pn(dim_num)
  real ( kind = 8 ) :: r1 = 3.0D+00
  real ( kind = 8 ) :: r2 = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0255:'
  write ( *, '(a)' ) '  ELLIPSE_POINT_NEAR_2D is given a point P, and'
  write ( *, '(a)' ) '  finds the nearest point PN on an ellipse in 2D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The ellipse is (X/R1)^2 + (Y/R2)^2 = 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a,f14.6)' ) '  R1 = ', r1
  write ( *, '(a,f14.6)' ) '  R2 = ', r2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           P                PN'
  write ( *, '(a)' ) ' '

  do i = -3, n + 3

    p(1) = ( real ( n - i, kind = 8 ) * 0.0D+00   &
           + real (     i, kind = 8 ) * 4.0D+00 ) &
           / real ( n,     kind = 8 )

    p(2) = ( real ( n - i, kind = 8 ) * 3.0D+00   &
           + real (     i, kind = 8 ) * 0.0D+00 ) &
           / real ( n,     kind = 8 )

    call ellipse_point_near_2d ( r1, r2, p, pn )

    write ( *, '(2x,2f8.4,2x,2f8.4)' ) p(1:dim_num), pn(1:dim_num)

  end do

  return
end
subroutine test026 ( )

!*****************************************************************************80
!
!! TEST026 tests ELLIPSE_POINTS_2D and ELLIPSE_AREA_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n_max = 24

  real ( kind = 8 ) area
  integer ( kind = 4 ) n
  real ( kind = 8 ), dimension ( dim_num ) :: pc = (/ 5.0D+00, -2.0D+00 /)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) psi
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) result
  real ( kind = 8 ) v(dim_num,n_max)

  r1 = 3.0D+00
  r2 = 1.0D+00
  psi = pi / 6.0D+00
  n = 16

  call ellipse_area_2d ( r1, r2, area )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST026'
  write ( *, '(a)' ) '  ELLIPSE_POINTS_2D returns points on an ellipse;'
  write ( *, '(a)' ) '  ELLIPSE_AREA_2D returns the area of an ellipse;'
  write ( *, '(a)' ) '  POLYGON_AREA_2D finds the area of a polygon.'

  call r8vec_print ( dim_num, pc, '  Ellipse center:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  radii R1 = ', r1, ' R2 = ', r2
  write ( *, '(a,g14.6)' ) '  and angle PSI = ', psi
  write ( *, '(a,g14.6)' ) '  and area = ', area

  call ellipse_points_2d ( pc, r1, r2, psi, n, v )

  call r8mat_transpose_print ( dim_num, n, v, '  Sample points:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For any N, the sampled points define a polygon'
  write ( *, '(a)' ) '  whose area approximates the ellipse area.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  N      Area'
  write ( *, '(a)' ) ' '
 
  do n = 3, n_max
    call ellipse_points_2d ( pc, r1, r2, psi, n, v )
    call polygon_area_2d ( n, v, result )
    write ( *, '(2x,i8,g14.6)' ) n, result
  end do
 
  return
end
subroutine test027 ( )

!*****************************************************************************80
!
!! TEST027 tests ELLIPSE_POINTS_ARC_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 13
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) p(dim_num,n)
  real ( kind = 8 ), dimension(dim_num) :: pc = (/ 5.0D+00, -2.0D+00 /)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) psi
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2

  r1 = 3.0D+00
  r2 = 1.0D+00
  psi = pi / 6.0D+00
  theta1 = pi / 2.0D+00
  theta2 = 2.0D+00 * pi

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST027'
  write ( *, '(a)' ) '  ELLIPSE_POINTS_ARC_2D returns points on an'
  write ( *, '(a)' ) '  elliptical arc.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  The ellipse has center ', pc(1), pc(2)
  write ( *, '(a,g14.6,a,g14.6)' ) '  radii R1 = ', r1, ' R2 = ', r2
  write ( *, '(a,g14.6)' ) '  and angle PSI = ', psi
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The arc extends from THETA1 = ', theta1
  write ( *, '(a,g14.6)') '  to THETA2 = ', theta2

  call ellipse_points_arc_2d ( pc, r1, r2, psi, theta1, theta2, n, p )

  call r8mat_transpose_print ( dim_num, n, p, '  Sample points:' )
 
  return
end
subroutine test028 ( )

!*****************************************************************************80
!
!! TEST028 tests HALFPLANE_CONTAINS_POINT_2D
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  logical, dimension ( test_num ) :: expected = (/ &
    .true., .false., .true., .false. /)
  logical halfplane_contains_point_2d
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p = reshape ( (/ &
    1.0D+00,   1.0D+00, &
    1.0D+00,  -1.0D+00, &
   -1.0D+00,   1.0D+00, &
    2.0D+00, 200.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p1 = reshape ( (/ &
    0.0D+00,   0.0D+00, &
    0.0D+00,   0.0D+00, &
   -5.0D+00,  -5.0D+00, &
    3.0D+00, 150.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p2 = reshape ( (/ &
    2.0D+00,   0.0D+00, &
    2.0D+00,   0.0D+00, &
   10.0D+00,  10.0D+00, &
    1.0D+00,  50.0D+00 /), (/ dim_num, test_num /) )
  logical temp
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST028'
  write ( *, '(a)' ) '  HALFPLANE_CONTAINS_POINT_2D determines whether a'
  write ( *, '(a)' ) '  halfplane bounded by PA:PB contains the'
  write ( *, '(a)' ) '  point P.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    temp = halfplane_contains_point_2d ( p1(1:dim_num,test), &
      p2(1:dim_num,test), p(1:dim_num,test) )

    write ( *, '(a)' ) ' '
    write ( *, '(a,2g14.6)' ) '  P1 = ', p1(1:dim_num,test)
    write ( *, '(a,2g14.6)' ) '  P2 = ', p2(1:dim_num,test)
    write ( *, '(a,2g14.6)' ) '  P =  ', p(1:dim_num,test)
    write ( *, '(a,l1,a,l1)' ) '  Contains? = ', temp, &
      '  Correct = ', expected(test)

  end do

  return
end
subroutine test029 ( )

!*****************************************************************************80
!
!! TEST029 tests HALFSPACE_IMP_TRIANGLE_INT_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 6

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  integer ( kind = 4 ) int_num
  real ( kind = 8 ) pint(dim_num,4)
  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ), dimension(dim_num,3,test_num) :: t_test = reshape ( (/ &
        0.0D+00,  0.0D+00,  0.0D+00, &
        0.0D+00, -1.0D+00,  0.0D+00, &
        0.0D+00,  0.0D+00, -2.0D+00, &
       -6.0D+00,  0.0D+00,  0.0D+00, &
        0.0D+00, -1.0D+00,  0.0D+00, &
        0.0D+00,  0.0D+00, -2.0D+00, &
        0.0D+00,  0.0D+00,  0.0D+00, &
        0.0D+00,  3.0D+00,  0.0D+00, &
        0.0D+00,  0.0D+00,  2.0D+00, &
       -6.0D+00,  0.0D+00,  0.0D+00, &
        0.0D+00,  4.0D+00,  0.0D+00, &
        0.0D+00,  0.0D+00,  3.0D+00, &
       -8.0D+00,  0.0D+00,  0.0D+00, &
        0.0D+00, -1.0D+00,  0.0D+00, &
        0.0D+00,  0.0D+00, -2.0D+00, &
        0.0D+00,  0.0D+00,  0.0D+00, &
        0.0D+00,  4.0D+00,  0.0D+00, &
        0.0D+00,  0.0D+00,  4.0D+00 /), (/ dim_num, 3, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST029'
  write ( *, '(a)' ) '  HALFSPACE_IMP_TRIANGLE_INT_3D finds'
  write ( *, '(a)' ) '  intersection points of an implicit'
  write ( *, '(a)' ) '  halfspace and a triangle.'

  a =   1.0D+00
  b = - 2.0D+00
  c = - 3.0D+00
  d =   6.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The implicitly defined bounding plane'
  write ( *, '(a)' ) '  has the form: A*X + B*Y + C*Z + D = 0.'
  write ( *, '(a,4g14.6)' ) '  A,B,C,D = ', a, b, c, d

  do test = 1, test_num

    t(1:dim_num,1:3) = t_test(1:dim_num,1:3,test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Case ', test
    write ( *, '(a)' ) ' '

    call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices' )

    call halfspace_imp_triangle_int_3d ( a, b, c, d, t, int_num, pint )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of intersection points is ', int_num
    write ( *, '(a)' ) ' '

    call r8mat_transpose_print ( dim_num, int_num, pint, '  Intersections:' )

  end do

  return
end
subroutine test030 ( )

!*****************************************************************************80
!
!! TEST030 tests HALFSPACE_NORMAL_TRIANGLE_INT_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 6

  integer ( kind = 4 ) int_num
  real ( kind = 8 ), dimension(dim_num) :: normal = (/ &
    2.0D+00, -4.0D+00, -6.0D+00 /)
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) pint(dim_num,4)
  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ), dimension(dim_num,3,test_num) :: t_test = reshape ( (/ &
        0.0D+00,  0.0D+00,  0.0D+00, &
        0.0D+00, -1.0D+00,  0.0D+00, &
        0.0D+00,  0.0D+00, -2.0D+00, &
       -6.0D+00,  0.0D+00,  0.0D+00, &
        0.0D+00, -1.0D+00,  0.0D+00, &
        0.0D+00,  0.0D+00, -2.0D+00, &
        0.0D+00,  0.0D+00,  0.0D+00, &
        0.0D+00,  3.0D+00,  0.0D+00, &
        0.0D+00,  0.0D+00,  2.0D+00, &
       -6.0D+00,  0.0D+00,  0.0D+00, &
        0.0D+00,  4.0D+00,  0.0D+00, &
        0.0D+00,  0.0D+00,  3.0D+00, &
       -8.0D+00,  0.0D+00,  0.0D+00, &
        0.0D+00, -1.0D+00,  0.0D+00, &
        0.0D+00,  0.0D+00, -2.0D+00, &
        0.0D+00,  0.0D+00,  0.0D+00, &
        0.0D+00,  4.0D+00,  0.0D+00, &
        0.0D+00,  0.0D+00,  4.0D+00 /), (/ dim_num, 3, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST030'
  write ( *, '(a)' ) '  HALFSPACE_NORMAL_TRIANGLE_INT_3D finds'
  write ( *, '(a)' ) '  intersection points of a normal form'
  write ( *, '(a)' ) '  halfspace and a triangle.'

  p(1:dim_num) = (/ -6.0D+00,  0.0D+00, 0.0D+00 /)

  call r8vec_print ( dim_num, p, '  Plane point P:' )
  call r8vec_print ( dim_num, normal, '  Plane normal:' )

  do test = 1, test_num

    t(1:dim_num,1:3) = t_test(1:dim_num,1:3,test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Case ', test
    write ( *, '(a)' ) ' '
    call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

    call halfspace_normal_triangle_int_3d ( p, normal, t, int_num, pint )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of intersection points is ', int_num
    write ( *, '(a)' ) ' '

    call r8mat_transpose_print ( dim_num, int_num, pint, '  Intersections:' )

  end do

  return
end
subroutine test031 ( )

!*****************************************************************************80
!
!! TEST031 tests HAVERSINE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) d
  real ( kind = 8 ) haversine
  real ( kind = 8 ) hx
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) radians_to_degrees
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 12
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST031'
  write ( *, '(a)' ) '  HAVERSINE computes the haversine of an angle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Degrees  Radians  Haversine'
  write ( *, '(a)' ) ' '

  do test = 0, test_num
    x = real ( test, kind = 8 ) * 2.0D+00 * pi / real ( test_num, kind = 8 )
    d = radians_to_degrees ( x )
    hx = haversine ( x )
    write ( *, '(2x,2f8.4,g14.6)' ) d, x, hx
  end do

  return
end
subroutine test0315 ( )

!*****************************************************************************80
!
!! TEST0315 tests HEXAGON_CONTAINS_POINT_2D
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 40
  integer ( kind = 4 ), parameter :: dim_num = 2

  character dot(n)
  real ( kind = 8 ), dimension(dim_num,6) :: h = reshape ( (/ &
    0.2D+00, 0.4D+00, &
    0.4D+00, 0.2D+00, &
    0.8D+00, 0.0D+00, &
    1.0D+00, 0.6D+00, &
    0.4D+00, 1.0D+00, &
    0.2D+00, 0.8D+00 /), (/ dim_num, 6 /) )
  logical hexagon_contains_point_2d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0315'
  write ( *, '(a)' ) '  HEXAGON_CONTAINS_POINT_2D reports if a hexagon'
  write ( *, '(a)' ) '  contains a point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We will call the function repeatedly, and draw'
  write ( *, '(a)' ) '  a sketch of an irregular hexagon in the unit square.'
  write ( *, '(a)' ) ' '

  do i = 1, n

    p(2) = real ( n - i, kind = 8 ) &
         / real ( n - 1, kind = 8 )

    do j = 1, n

      p(1) = real ( j - 1, kind = 8 ) &
           / real ( n - 1, kind = 8 )

      if ( hexagon_contains_point_2d ( h, p ) ) then
        dot(j) = '*'
      else
        dot(j) = '-'
      end if

    end do
    write ( *, '(2x,40a1)' ) dot(1:n)
  end do

  return
end
subroutine test032 ( )

!*****************************************************************************80
!
!! TEST032 tests HEXAGON_SHAPE_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  real ( kind = 8 ) p(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST032'
  write ( *, '(a)' ) '  HEXAGON_SHAPE_2D: points on a unit hexagon.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Angle    X    Y '
  write ( *, '(a)' ) ' '

  do i = -10, 370, 10
    angle = real ( i, kind = 8 )
    call hexagon_shape_2d ( angle, p )
    write ( *, '(2x,3g14.6)' ) angle, p(1:dim_num)
  end do

  return
end
subroutine test0321 ( )

!*****************************************************************************80
!
!! TEST0321 tests HEXAGON_VERTICES_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) p(dim_num,6)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0321'
  write ( *, '(a)' ) '  HEXAGON_VERTICES_2D: the vertices of the unit hexagon.'

  call hexagon_vertices_2d ( p )

  call r8mat_transpose_print ( dim_num, 6, p, '  Vertices:' )

  return
end
subroutine test0322 ( )

!*****************************************************************************80
!
!! TEST0322 tests I4COL_FIND_ITEM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: test_num = 3

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) item
  integer ( kind = 4 ), dimension ( test_num ) :: item_test = (/ &
    34, 12, 90 /)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0322'
  write ( *, '(a)' ) '  I4COL_FIND_ITEM finds the first occurrence of'
  write ( *, '(a)' ) '  an item in an integer array of columns.'
 
  do i = 1, m
    do j = 1, n
      a(i,j) = 10 * i + j
    end do
  end do

  call i4mat_print ( m, n, a, '  The matrix of columns:' )

  do test = 1, test_num
 
    item = item_test(test)

    call i4col_find_item ( m, n, a, item, row, col )

    write ( *, '(a,i8,a,i8,a,i8)' ) '  Item ', item, '  occurs in row ', &
      row, ' and column ', col

  end do

  return
end
subroutine test0323 ( )

!*****************************************************************************80
!
!! TEST0323 tests I4COL_FIND_PAIR_WRAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: test_num = 5

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) item1
  integer ( kind = 4 ), dimension ( test_num ) :: item1_test = (/ &
    22, 32, 22, 54, 54 /)
  integer ( kind = 4 ) item2
  integer ( kind = 4 ), dimension ( test_num ) :: item2_test = (/ &
    32, 22, 23, 14, 11 /)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0323'
  write ( *, '(a)' ) '  I4COL_FIND_PAIR_WRAP finds the first occurrence of'
  write ( *, '(a)' ) '  a pair of item in an integer array of columns.'
  write ( *, '(a)' ) '  Items in the array are ordered by column, and'
  write ( *, '(a)' ) '  wraparound is allowed.'
 
  do i = 1, m
    do j = 1, n
      a(i,j) = 10 * i + j
    end do
  end do

  call i4mat_print ( m, n, a, '  The matrix of columns:' )

  do test = 1, test_num
 
    item1 = item1_test(test)
    item2 = item2_test(test)

    call i4col_find_pair_wrap ( m, n, a, item1, item2, row, col )

    write ( *, '(a,i8,a,i8,a,i8,a,i8)' ) '  Item ', item1, &
      ' followed by item ', item2, ' occurs in row ', &
      row, ' and column ', col

  end do

  return
end
subroutine test0325 ( )

!*****************************************************************************80
!
!! TEST0325 tests ICOS_SIZE and ICOS_SHAPE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) ::  edge_point
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ), allocatable, dimension ( : ) ::  face_order
  integer ( kind = 4 ) face_order_max
  integer ( kind = 4 ), allocatable, dimension ( :, : ) ::  face_point
  integer ( kind = 4 ) point_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_coord

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0325'
  write ( *, '(a)' ) '  For the icosahedron,'
  write ( *, '(a)' ) '  ICOS_SIZE returns dimension information;'
  write ( *, '(a)' ) '  ICOS_SHAPE returns face and order information.'
  write ( *, '(a)' ) '  SHAPE_PRINT_3D prints this information.'

  call icos_size ( point_num, edge_num, face_num, face_order_max )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of vertices = ', point_num
  write ( *, '(a,i8)' ) '  Number of edges =    ', edge_num
  write ( *, '(a,i8)' ) '  Number of faces =    ', face_num
  write ( *, '(a,i8)' ) '  Maximum face order = ', face_order_max

  allocate ( edge_point(1:2,1:edge_num) )
  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )
  allocate ( point_coord(1:3,1:point_num) )

  call icos_shape ( point_num, edge_num, face_num, face_order_max, &
    point_coord, edge_point, face_order, face_point )

  call shape_print_3d ( point_num, face_num, face_order_max, &
    point_coord, face_order, face_point )

  deallocate ( edge_point )
  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )

  return
end
subroutine test0327 ( )

!*****************************************************************************80
!
!! TEST0327 tests LINE_EXP_NORMAL_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) normal(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0327'
  write ( *, '(a)' ) '  LINE_EXP_NORMAL_2D determines a unit normal vector'
  write ( *, '(a)' ) '  to a given explicit line.'

  p1(1:dim_num) = (/ 1.0D+00, 3.0D+00 /)
  p2(1:dim_num) = (/ 4.0D+00, 0.0D+00 /)

  call r8vec_print ( dim_num, p1, '  Point 1: ' )
  call r8vec_print ( dim_num, p2, '  Point 2: ' )

  call line_exp_normal_2d ( p1, p2, normal )

  call r8vec_print ( dim_num, normal, '  Normal vector N:' )

  return
end
subroutine test033 ( )

!*****************************************************************************80
!
!! TEST033 tests LINE_EXP_PERP_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  logical              flag
  real ( kind = 8 ), dimension(dim_num) :: p1 = (/ 1.0D+00, 3.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p2 = (/ 4.0D+00, 0.0D+00 /)
  real ( kind = 8 ) p3(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p3_test = reshape ( (/ &
    0.0D+00,  0.0D+00, &
    5.0D+00, -1.0D+00, &
    5.0D+00,  3.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) p4(dim_num)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST033'
  write ( *, '(a)' ) '  LINE_EXP_PERP_2D is given an explicit line (P1,P2),'
  write ( *, '(a)' ) '  and another point P3.  It then finds a point'
  write ( *, '(a)' ) '  P4 on (P1,P2) so that (P1,P2) is perpendicular'
  write ( *, '(a)' ) '  to (P3,P4).'

  call r8vec_print ( dim_num, p1, '  Point P1:' )
  call r8vec_print ( dim_num, p2, '  Point P2:' )

  do test = 1, test_num

    p3(1:dim_num) = p3_test(1:dim_num,test)

    call r8vec_print ( dim_num, p3, '  Point P3:' )

    call line_exp_perp_2d ( p1, p2, p3, p4, flag )

    call r8vec_print ( dim_num, p4, '  Point P4:' )

  end do

  return
end
subroutine test0335 ( )

!*****************************************************************************80
!
!! TEST0335 tests LINE_EXP_POINT_DIST_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) dist
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: p1 = (/ 1.0D+00, 3.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p2 = (/ 4.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = &
  reshape ( (/ &
    0.0D+00,  0.0D+00, &
    5.0D+00, -1.0D+00, &
    5.0D+00,  3.0D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0335'
  write ( *, '(a)' ) '  LINE_EXP_POINT_DIST_2D finds the distance from'
  write ( *, '(a)' ) '  an explicit line to a point in 2D.'

  call r8vec_print ( dim_num, p1, '  Point 1: ' )
  call r8vec_print ( dim_num, p2, '  Point 2: ' )

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call r8vec_print ( dim_num, p, '  Point: ' )

    call line_exp_point_dist_2d ( p1, p2, p, dist )

    write ( *, '(a,g14.6)' ) '  Distance = ', dist

  end do

  return
end
subroutine test0336 ( )

!*****************************************************************************80
!
!! TEST0336 tests LINE_EXP_POINT_DIST_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) dist
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: p1 = (/ 1.0D+00, 3.0D+00, 2.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p2 = (/ 4.0D+00, 0.0D+00, 1.0D+00 /)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = reshape ( (/ &
    0.0D+00,  0.0D+00, 2.0D+00, &
    5.0D+00, -1.0D+00, 1.0D+00, &
    5.0D+00,  3.0D+00, 3.0D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0336'
  write ( *, '(a)' ) '  LINE_EXP_POINT_DIST_3D finds the distance'
  write ( *, '(a)' ) '  from an explicit line to a point in 3D.'

  call r8vec_print ( dim_num, p1, '  Point 1: ' )
  call r8vec_print ( dim_num, p2, '  Point 2: ' )

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call r8vec_print ( dim_num, p, '  Point: ' )

    call line_exp_point_dist_3d ( p1, p2, p, dist )

    write ( *, '(a,g14.6)' ) '  Distance = ', dist

  end do

  return
end
subroutine test0337 ( )

!*****************************************************************************80
!
!! TEST0337 tests LINE_EXP_POINT_DIST_SIGNED_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) dist
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: p1 = (/ 1.0D+00, 3.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p2 = (/ 4.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = reshape ( (/ &
    0.0D+00,  0.0D+00, &
    5.0D+00, -1.0D+00, &
    5.0D+00,  3.0D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0337'
  write ( *, '(a)' ) '  LINE_EXP_POINT_DIST_SIGNED_2D finds the signed'
  write ( *, '(a)' ) '  distance to a point from an explicit line.'

  call r8vec_print ( dim_num, p1, '  Point 1: ' )
  call r8vec_print ( dim_num, p2, '  Point 2: ' )

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call r8vec_print ( dim_num, p, '  Point: ' )

    call line_exp_point_dist_signed_2d ( p1, p2, p, dist )

    write ( *, '(a,g14.6)' ) '  Signed distance = ', dist

  end do

  return
end
subroutine test034 ( )

!*****************************************************************************80
!
!! TEST034 tests LINE_EXP_POINT_NEAR_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) dist
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: p1 = (/ 1.0D+00, 3.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p2 = (/ 4.0D+00, 0.0D+00 /)
  real ( kind = 8 ) pn(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = reshape ( (/ &
    0.0D+00,  0.0D+00, &
    5.0D+00, -1.0D+00, &
    5.0D+00,  3.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) t
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST034'
  write ( *, '(a)' ) '  LINE_EXP_POINT_NEAR_2D finds the point on'
  write ( *, '(a)' ) '  a line nearest in point in 2D.'

  call r8vec_print ( dim_num, p1, '  The point P1:' )
  call r8vec_print ( dim_num, p2, '  The point P2:' )

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call r8vec_print ( dim_num, p, '  The point P:' )

    call line_exp_point_near_2d ( p1, p2, p, pn, dist, t )

    call r8vec_print ( dim_num, pn, '  Nearest point PN:' )

    write ( *, '(a,g14.6)' ) '  Distance = ', dist
    write ( *, '(a,g14.6)' ) '  Relative line position T = ', t

  end do

  return
end
subroutine test0345 ( )

!*****************************************************************************80
!
!! TEST0345 tests LINE_EXP2IMP_2D and LINE_IMP2EXP_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) :: a1 = 1.0D+00
  real ( kind = 8 ) a2
  real ( kind = 8 ) :: b1 = 2.0D+00
  real ( kind = 8 ) b2
  real ( kind = 8 ) :: c1 = 3.0D+00
  real ( kind = 8 ) c2
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0345'
  write ( *, '(a)' ) '  LINE_EXP2IMP_2D converts explicit to implicit lines.'
  write ( *, '(a)' ) '  LINE_IMP2EXP_2D converts implicit to explicit lines.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,3f8.4)' ) '  Implicit line A, B, C = ', a1, b1, c1

  call line_imp2exp_2d ( a1, b1, c1, p1, p2 )

  call r8vec_print ( dim_num, p1, '  The point P1:' )
  call r8vec_print ( dim_num, p2, '  The point P2:' )

  call line_exp2imp_2d ( p1, p2, a2, b2, c2 )

  write ( *, '(a,3f8.4)' ) '  Recovered implicit line A, B, C = ', a2, b2, c2

  return
end
subroutine test0346 ( )

!*****************************************************************************80
!
!! TEST0346 tests LINE_EXP2PAR_2D and LINE_PAR2EXP_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) :: f1 = 1.0D+00
  real ( kind = 8 ) f2
  real ( kind = 8 ) :: g1 = 2.0D+00
  real ( kind = 8 ) g2
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) :: x1 = 3.0D+00
  real ( kind = 8 ) x2
  real ( kind = 8 ) :: y1 = 4.0D+00
  real ( kind = 8 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0346'
  write ( *, '(a)' ) '  LINE_EXP2PAR_2D converts explicit to parametric lines.'
  write ( *, '(a)' ) '  LINE_PAR2EXP_2D converts parametric to explicit lines.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,3f8.4)' ) '  Parametric line:'
  write ( *, '(a,2f8.4)' ) '    F,  G =  ', f1, g1
  write ( *, '(a,2f8.4)' ) '    X0, Y0 = ', x1, y1

  call line_par2exp_2d ( f1, g1, x1, y1, p1, p2 )

  call r8vec_print ( dim_num, p1, '  The point P1:' )
  call r8vec_print ( dim_num, p2, '  The point P2:' )

  call line_exp2par_2d ( p1, p2, f2, g2, x2, y2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,3f8.4)' ) '  Recovered parametric line:'
  write ( *, '(a,2f8.4)' ) '    F,  G =  ', f2, g2
  write ( *, '(a,2f8.4)' ) '    X0, Y0 = ', x2, y2

  return
end
subroutine test035 ( )

!*****************************************************************************80
!
!! TEST035 tests LINE_IMP_POINT_DIST_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) a
  real ( kind = 8 ), parameter, dimension ( test_num ) :: a_test = (/ &
    2.0D+00, 2.0D+00, 2.0D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), parameter, dimension ( test_num ) :: b_test = (/ &
    5.0D+00, 5.0D+00, 5.0D+00 /)
  real ( kind = 8 ) c
  real ( kind = 8 ), parameter, dimension ( test_num ) :: c_test = (/ &
    3.0D+00, 3.0D+00, 3.0D+00 /)
  real ( kind = 8 ) dist
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), parameter, dimension ( dim_num, test_num ) :: &
    p_test = reshape ( (/ &
    0.0D+00, 6.0D+00, &
    0.0D+00, 5.0D+00, &
    0.0D+00, 4.0D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST035'
  write ( *, '(a)' ) '  LINE_IMP_POINT_DIST_2D finds the distance from'
  write ( *, '(a)' ) '  a point P to a line A * X + B * Y + C = 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   X       Y       A       B       C       DIST'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    a = a_test(test)
    b = b_test(test)
    c = c_test(test)
    p(1:dim_num) = p_test(1:dim_num,test)

    call line_imp_point_dist_2d ( a, b, c, p, dist )

    write ( *, '(2x,6f8.4)' ) p(1:dim_num), a, b, c, dist

  end do

  return
end
subroutine test038 ( )

!*****************************************************************************80
!
!! TEST038 tests LINES_EXP_ANGLE_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 July 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) angle
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p1_test = &
    reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 2.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p2_test = &
    reshape ( (/ &
    1.0D+00, 2.0D+00, 0.0D+00, &
    1.0D+00, 2.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) q1(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: q1_test = &
    reshape ( (/ &
    0.0D+00, 3.0D+00,  3.0D+00, &
    1.0D+00, 2.0D+00, -1.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) q2(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: q2_test = &
    reshape ( (/ &
    3.0D+00, 0.0D+00, 3.0D+00, &
    1.0D+00, 2.0D+00, 3.0D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST038'
  write ( *, '(a)' ) '  LINES_EXP_ANGLE_3D finds the angle between'
  write ( *, '(a)' ) '  two explicit lines in 3D;'

  do test = 1, test_num

    p1(1:dim_num) = p1_test(1:dim_num,test)
    p2(1:dim_num) = p2_test(1:dim_num,test)
    q1(1:dim_num) = q1_test(1:dim_num,test)
    q2(1:dim_num) = q2_test(1:dim_num,test)

    call lines_exp_angle_3d ( p1, p2, q1, q2, angle )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Angle between lines is ', angle

  end do

  return
end
subroutine test0385 ( )

!*****************************************************************************80
!
!! TEST0385 tests LINES_EXP_DIST_3D and LINES_EXP_DIST_3D_2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) dist
  real ( kind = 8 ) dist2
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p1_test = &
    reshape ( (/ &
    0.0D+00,  0.0D+00, 0.0D+00, &
    4.0D+00, -3.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p2_test = &
    reshape ( (/ &
    1.0D+00, 2.0D+00, 0.0D+00, &
   -8.0D+00, 6.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) q1(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: q1_test = &
    reshape ( (/ &
    0.0D+00, 3.0D+00,  3.0D+00, &
    3.0D+00, 4.0D+00, -1.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) q2(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: q2_test = &
    reshape ( (/ &
    3.0D+00, 0.0D+00, 3.0D+00, &
    3.0D+00, 4.0D+00, 3.0D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0385'
  write ( *, '(a)' ) '  LINES_EXP_DIST_3D finds the distance between'
  write ( *, '(a)' ) '  two explicit lines in 3D.'
  write ( *, '(a)' ) '  LINES_EXP_DIST_3D_2 finds the distance between'
  write ( *, '(a)' ) '  two explicit lines in 3D.'

  do test = 1, test_num

    p1(1:dim_num) = p1_test(1:dim_num,test)
    p2(1:dim_num) = p2_test(1:dim_num,test)
    q1(1:dim_num) = q1_test(1:dim_num,test)
    q2(1:dim_num) = q2_test(1:dim_num,test)

    call lines_exp_dist_3d ( p1, p2, q1, q2, dist )
    call lines_exp_dist_3d_2 ( p1, p2, q1, q2, dist2 )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,3g14.6)' ) '  P1:', p1(1:dim_num)
    write ( *, '(a,3g14.6)' ) '  P2:', p2(1:dim_num)
    write ( *, '(a,3g14.6)' ) '  Q1:', q1(1:dim_num)
    write ( *, '(a,3g14.6)' ) '  Q2:', q2(1:dim_num)
    write ( *, '(a,g14.6)' ) '  LINES_EXP_DIST_3D =   ', dist
    write ( *, '(a,g14.6)' ) '  LINES_EXP_DIST_3D_2 = ', dist2

  end do

  return
end
subroutine test03855 ( )

!*****************************************************************************80
!
!! TEST03855 tests LINES_EXP_NEAR_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p1_test = &
    reshape ( (/ &
    0.0D+00,  0.0D+00, 0.0D+00, &
    4.0D+00, -3.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p2_test = &
    reshape ( (/ &
    1.0D+00, 2.0D+00, 0.0D+00, &
   -8.0D+00, 6.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) pn(dim_num)
  real ( kind = 8 ) q1(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: q1_test = &
    reshape ( (/ &
    0.0D+00, 3.0D+00,  3.0D+00, &
    3.0D+00, 4.0D+00, -1.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) q2(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: q2_test = &
    reshape ( (/ &
    3.0D+00, 0.0D+00, 3.0D+00, &
    3.0D+00, 4.0D+00, 3.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) qn(dim_num)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03855'
  write ( *, '(a)' ) '  LINES_EXP_NEAR_3D finds nearest points on'
  write ( *, '(a)' ) '  two explicit lines in 3D.'

  do test = 1, test_num

    p1(1:dim_num) = p1_test(1:dim_num,test)
    p2(1:dim_num) = p2_test(1:dim_num,test)
    q1(1:dim_num) = q1_test(1:dim_num,test)
    q2(1:dim_num) = q2_test(1:dim_num,test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,3g14.6)' ) '  P1:', p1(1:dim_num)
    write ( *, '(a,3g14.6)' ) '  P2:', p2(1:dim_num)
    write ( *, '(a,3g14.6)' ) '  Q1:', q1(1:dim_num)
    write ( *, '(a,3g14.6)' ) '  Q2:', q2(1:dim_num)

    call lines_exp_near_3d ( p1, p2, q1, q2, pn, qn )

    write ( *, '(a)' ) ' '
    write ( *, '(a,3g14.6)' ) '  PN:', pn(1:dim_num)
    write ( *, '(a,3g14.6)' ) '  QN:', qn(1:dim_num)

  end do

  return
end
subroutine test0386 ( )

!*****************************************************************************80
!
!! TEST0386 tests LINES_EXP_EQUAL_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 6

  logical equal
  logical lines_exp_equal_2d
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p1_test = &
    reshape ( (/ &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p2_test = &
    reshape ( (/ &
    1.0D+00, 2.0D+00, &
    1.0D+00, 2.0D+00, &
    1.0D+00, 2.0D+00, &
    1.0D+00, 2.0D+00, &
    1.0D+00, 2.0D+00, &
    1.0D+00, 2.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) q1(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: q1_test = &
    reshape ( (/ &
    0.0D+00,  0.0D+00, &
    1.0D+00,  2.0D+00, &
    0.0D+00,  0.0D+00, &
    7.0D+00, 14.0D+00, &
    1.0D+00,  2.0D+00, &
    0.0D+00, 10.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) q2(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: q2_test = &
    reshape ( (/ &
    1.0D+00,  2.0D+00, &
    0.0D+00,  0.0D+00, &
    2.0D+00,  4.0D+00, &
    5.5D+00, 11.0D+00, &
    3.0D+00,  5.0D+00, &
    1.0D+00, 12.0D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0386'
  write ( *, '(a)' ) '  LINES_EXP_EQUAL_2D tries to determine if two '
  write ( *, '(a)' ) '  explicit lines in 2D are equal.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p1(1:dim_num) = p1_test(1:dim_num,test)
    p2(1:dim_num) = p2_test(1:dim_num,test)
    q1(1:dim_num) = q1_test(1:dim_num,test)
    q2(1:dim_num) = q2_test(1:dim_num,test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,2g14.6)' ) '  P1', p1(1:dim_num)
    write ( *, '(a,2g14.6)' ) '  P2', p2(1:dim_num)
    write ( *, '(a)' ) ' '
    write ( *, '(a,2g14.6)' ) '  Q1', q1(1:dim_num)
    write ( *, '(a,2g14.6)' ) '  Q2', q2(1:dim_num)
 
    equal = lines_exp_equal_2d ( p1, p2, q1, q2 )
 
    if ( equal ) then
      write ( *, '(a)' ) '  The lines are equal.'
    else
      write ( *, '(a)' ) '  The lines are distinct.'
    end if

  end do
 
  return
end
subroutine test039 ( )

!*****************************************************************************80
!
!! TEST039 tests LINES_EXP_INT_2D.
!
!  Discussion:
!
!    Test #1:
!
!      x + 2y -  4 = 0
!      x -  y -  1 = 0
!
!    Test #2: 
!
!      x + 2y -  4 = 0
!     2x + 4y -  1 = 0
!
!    Test #3:
!
!      x + 2y -  4 = 0
!    -3x - 6y + 12 = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 July 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  integer ( kind = 4 ) ival
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p1_test = &
    reshape ( (/ &
    0.0D+00, 0.0D+00, &
    0.0D+00, 2.0D+00, &
    0.0D+00, 2.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p2_test = &
    reshape ( (/ &
    4.0D+00, 0.0D+00, &
    4.0D+00, 0.0D+00, &
    4.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) q1(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: q1_test = &
    reshape ( (/ &
    0.0D+00, -1.0D+00, &
    0.0D+00,  0.25D+00, &
    0.0D+00,  2.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) q2(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: q2_test = &
    reshape ( (/ &
    1.0D+00,  0.0D+00, &
    0.5D+00,  0.0D+00, &
    4.0D+00,  0.0D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST039'
  write ( *, '(a)' ) '  LINES_EXP_INT_2D finds intersections of '
  write ( *, '(a)' ) '  two explicit lines in 2D.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p1(1:dim_num) = p1_test(1:dim_num,test)
    p2(1:dim_num) = p2_test(1:dim_num,test)
    q1(1:dim_num) = q1_test(1:dim_num,test)
    q2(1:dim_num) = q2_test(1:dim_num,test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,2g14.6)' ) '  P1', p1(1:dim_num)
    write ( *, '(a,2g14.6)' ) '  P2', p2(1:dim_num)
    write ( *, '(a)' ) ' '
    write ( *, '(a,2g14.6)' ) '  Q1', q1(1:dim_num)
    write ( *, '(a,2g14.6)' ) '  Q2', q2(1:dim_num)
 
    call lines_exp_int_2d ( p1, p2, q1, q2, ival, p )
 
    if ( ival == 1 ) then
      write ( *, '(a,2g14.6)' ) '  Intersection at ', p(1:dim_num)
    else if ( ival == 0 ) then
      write ( *, '(a)' ) '  Lines are parallel, no intersection.'
    else if ( ival == 2 ) then
      write ( *, '(a)' ) '  Lines are coincident.'
    else
      write ( *, '(a,i8)' ) '  Unknown return value of IVAL = ', ival
    end if

  end do
 
  return
end
subroutine test040 ( )

!*****************************************************************************80
!
!! TEST040 tests LINES_IMP_ANGLE_2D.
!
!  Discussion:
!
!    Test 1:
!
!      x + 2y - 4 = 0
!      x - y - 1 = 0
!
!    Test 2:
!
!      x + 2y - 4 = 0
!     2x + 4y - 1 = 0
!
!    Test 3:
!
!      x + 2y - 4 = 0
!    -3x - 6y +12 = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ), dimension ( test_num ) :: a2_test = (/ &
    1.0D+00, 2.0D+00, -3.0D+00 /)
  real ( kind = 8 ) angle
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ), dimension ( test_num ) :: b2_test = (/ &
    -1.0D+00, 4.0D+00, -6.0D+00 /)
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ), dimension ( test_num ) :: c2_test = (/ &
    -1.0D+00, -1.0D+00, 12.0D+00 /)
  real ( kind = 8 ) radians_to_degrees
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST040'
  write ( *, '(a)' ) '  LINES_IMP_ANGLE_2D finds the angle between'
  write ( *, '(a)' ) '  two lines written in implicit form.'
  write ( *, '(a)' ) ' '

  a1 =  1.0D+00
  b1 =  2.0D+00
  c1 = -4.0D+00

  do test = 1, test_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,3g14.6)' ) '  Line 1 coefficients:', a1, b1, c1

    a2 =  a2_test(test)
    b2 =  b2_test(test)
    c2 =  c2_test(test)

    write ( *, '(a,3g14.6)' ) '  Line 2 coefficients:', a2, b2, c2

    call lines_imp_angle_2d ( a1, b1, c1, a2, b2, c2, angle )

    angle = radians_to_degrees ( angle )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Angle between lines is ', angle

  end do

  return
end
subroutine test041 ( )

!*****************************************************************************80
!
!! TEST041 tests LINES_IMP_DIST_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) a1
  real ( kind = 8 ), parameter, dimension ( test_num ) :: a1_test = &
    (/  4.0D+00,  2.0D+00, 1.0D+00 /)
  real ( kind = 8 ) a2
  real ( kind = 8 ), parameter, dimension ( test_num ) :: a2_test = &
    (/  4.0D+00,  4.0D+00, 2.0D+00 /)
  real ( kind = 8 ) b1
  real ( kind = 8 ), parameter, dimension ( test_num ) :: b1_test = &
    (/ -1.0D+00, -1.0D+00, 2.0D+00 /)
  real ( kind = 8 ) b2
  real ( kind = 8 ), parameter, dimension ( test_num ) :: b2_test = &
    (/ -1.0D+00, -2.0D+00, 3.0D+00 /)
  real ( kind = 8 ) c1
  real ( kind = 8 ), parameter, dimension ( test_num ) :: c1_test = &
    (/  3.0D+00,  0.0D+00, 2.0D+00 /)
  real ( kind = 8 ) c2 
  real ( kind = 8 ), parameter, dimension ( test_num ) :: c2_test = &
    (/ 12.0D+00,  6.0D+00, 1.0D+00 /)
  real ( kind = 8 ) dist
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST041'
  write ( *, '(a)' ) '  LINES_IMP_DIST_3D finds the distance between'
  write ( *, '(a)' ) '  two implicit lines in 2D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   A1      B1      C1      A2      B2      C2   DIST'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    a1 = a1_test(test)
    b1 = b1_test(test)
    c1 = c1_test(test)
    a2 = a2_test(test)
    b2 = b2_test(test)
    c2 = c2_test(test)

    call lines_imp_dist_2d ( a1, b1, c1, a2, b2, c2, dist )

    write ( *, '(2x,7f8.4)' ) a1, b1, c1, a2, b2, c2, dist

  end do

  return
end
subroutine test0415 ( )

!*****************************************************************************80
!
!! TEST0415 tests LINES_IMP_INT_2D.
!
!  Discussion:
!
!    Test 1:
!
!      x + 2y - 4 = 0
!      x - y - 1 = 0
!
!    Test 2:
!
!      x + 2y - 4 = 0
!     2x + 4y - 1 = 0
!
!    Test 3:
!
!      x + 2y - 4 = 0
!    -3x - 6y +12 = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ), dimension ( test_num ) :: a2_test = (/ &
    1.0D+00, 2.0D+00, -3.0D+00 /)
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ), dimension ( test_num ) :: b2_test = (/ &
    -1.0D+00, 4.0D+00, -6.0D+00 /)
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ), dimension ( test_num ) :: c2_test = (/ &
    -1.0D+00, -1.0D+00, 12.0D+00 /)
  integer ( kind = 4 ) ival
  real ( kind = 8 ) p(dim_num)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0415'
  write ( *, '(a)' ) '  LINES_IMP_INT_2D finds the intersection of'
  write ( *, '(a)' ) '  two lines written in implicit form.'
  write ( *, '(a)' ) ' '

  a1 =  1.0D+00
  b1 =  2.0D+00
  c1 = -4.0D+00

  do test = 1, test_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,3g14.6)' ) '  Line 1 coefficients:', a1, b1, c1

    a2 = a2_test(test)
    b2 = b2_test(test)
    c2 = c2_test(test)

    write ( *, '(a,3g14.6)' ) '  Line 2 coefficients:', a2, b2, c2

    call lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, p )

    if ( ival == 1 ) then
      write ( *, '(a,2g14.6)' ) '  Intersection at ', p(1:dim_num)
    else if ( ival == 0 ) then
      write ( *, '(a)' ) '  Lines are parallel, no intersection.'
    else if ( ival == 2 ) then
      write ( *, '(a)' ) '  Lines are coincident.'
    else
      write ( *, '(a,i8)' ) '  Unknown return value of ival = ', ival
    end if

  end do
 
  return
end
subroutine test0416 ( )

!*****************************************************************************80
!
!! TEST0416 tests LINES_PAR_INT_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) g1
  real ( kind = 8 ) g2
  real ( kind = 8 ) pint(dim_num)
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0416'
  write ( *, '(a)' ) '  LINES_PAR_INT_2D finds the intersection of'
  write ( *, '(a)' ) '  two lines written in parametric form.'
  write ( *, '(a)' ) ' '
!
!  x - 2y = -1
!
  x1 =  0.0D+00
  y1 =  1.0D+00
  f1 =  2.0D+00
  g1 =  1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,4g14.6)' ) '  Line 1 parameters:', x1, y1, f1, g1
!
!  x + y - 8 = 0
!
  x2 = 10.0D+00
  y2 = -2.0D+00
  f2 =  1.0D+00
  g2 =  1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,4g14.6)' ) '  Line 2 parameters:', x2, y2, f2, g2

  call lines_par_int_2d ( f1, g1, x1, y1, f2, g2, x2, y2, t1, t2, pint )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Line 1 evaluated at T1:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '    T1 =   ', t1
  write ( *, '(a,g14.6)' ) '    X(T1)= ', x1 + f1 * t1
  write ( *, '(a,g14.6)' ) '    Y(T1)= ', y1 + g1 * t1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Line 2 evaluated at T2:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '    T2 =   ', t2
  write ( *, '(a,g14.6)' ) '    X(T2)= ', x2 + f2 * t2
  write ( *, '(a,g14.6)' ) '    Y(T2)= ', y2 + g2 * t2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reported intersection PINT:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  ', pint(1)
  write ( *, '(a,g14.6)' ) '  ', pint(2)

  return
end
subroutine test0418 ( )

!*****************************************************************************80
!
!! TEST0418 tests SEGMENTS_CURVATURE_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 13

  real ( kind = 8 ) curvature
  real ( kind = 8 ), dimension(dim_num) :: p1 = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p2 = (/ 1.0D+00, 0.0D+00 /)
  real ( kind = 8 ) p3(dim_num)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta_degrees
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0418'
  write ( *, '(a)' ) '  SEGMENTS_CURVATURE_2D computes the local curvature '
  write ( *, '(a)' ) '  defined by the line segments [P1,P2] and [P2,P3].'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Our three points are:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    P1 = (0,0)'
  write ( *, '(a)' ) '    P2 = (1,0)'
  write ( *, '(a)' ) '    P3 = (C,S)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  C = cosine ( theta), S = sine ( theta ).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test  Theta  Curvature'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    theta = 2.0D+00 * pi * real ( test     - 1, kind = 8 ) &
                         / real ( test_num - 1, kind = 8 )

    theta_degrees = 360.0D+00 * real ( test     - 1, kind = 8 ) &
                         /      real ( test_num - 1, kind = 8 )

    p3(1:dim_num) = (/ cos ( theta ), sin ( theta ) /)

    call segments_curvature_2d ( p1, p2, p3, curvature )

    write ( *, '(2x,i4,2x,f5.0,2x,g14.6)' ) test, theta_degrees, curvature

  end do

  return
end
subroutine test042 ( )

!*****************************************************************************80
!
!! TEST042 tests SEGMENTS_DIST_2D.
!
!  Discussion:
!
!    Case 1, parallel, not coincident.
!    Case 2, parallel, coincident, overlapping.
!    Case 3, parallel, coincident, disjoint.
!    Case 4, nonparallel, intersecting.
!    Case 5, nonparallel, disjoint.
!    Case 6 and 7, should be same, because simply a translation by 50;
!    Case 8 and 9, answers should be same.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 9

  real ( kind = 8 ) dist
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p1_test = &
    reshape ( (/ &
    2.0D+00,   3.0D+00, &
    2.0D+00,   3.0D+00, &
    2.0D+00,   3.0D+00, &
    2.0D+00,   3.0D+00, &
    2.0D+00,   3.0D+00, &
   57.0D+00,  53.0D+00, &
    7.0D+00,   3.0D+00, &
    0.0D+00,   0.0D+00, &
  -10.0D+00, -10.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p2_test = &
    reshape ( (/ &
    8.0D+00,   6.0D+00, &
    8.0D+00,   6.0D+00, &
    8.0D+00,   6.0D+00, &
    8.0D+00,   6.0D+00, &
    8.0D+00,   6.0D+00, &
   58.0D+00,  53.0D+00, &
    8.0D+00,   3.0D+00, &
  100.0D+00, 100.0D+00, &
  100.0D+00, 100.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) q1(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: q1_test = &
    reshape ( (/ &
    8.0D+00,  3.0D+00, &
    4.0D+00,  4.0D+00, &
   14.0D+00,  9.0D+00, &
    0.0D+00,  8.0D+00, &
    7.0D+00,  3.0D+00, &
   65.0D+00, 45.0D+00, &
   15.0D+00, -5.0D+00, &
   50.0D+00,  0.0D+00, &
   50.0D+00,  0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) q2(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: q2_test = &
    reshape ( (/ &
   14.0D+00,  6.0D+00, &
   14.0D+00,  9.0D+00, &
   16.0D+00, 10.0D+00, &
    5.0D+00,  3.0D+00, &
    9.0D+00, -1.0D+00, &
   57.0D+00, 53.0D+00, &
    7.0D+00,  3.0D+00, &
   60.0D+00,  0.0D+00, &
   60.0D+00,  0.0D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST042'
  write ( *, '(a)' ) '  SEGMENTS_DIST_2D computes the distance between'
  write ( *, '(a)' ) '  line segments in 2D.'

  do test = 1, test_num

    p1(1:dim_num) = p1_test(1:dim_num,test)
    p2(1:dim_num) = p2_test(1:dim_num,test)
    q1(1:dim_num) = q1_test(1:dim_num,test)
    q2(1:dim_num) = q2_test(1:dim_num,test)

    call segments_dist_2d ( p1, p2, q1, q2, dist )

    if ( test == 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Same slope, different intercepts.'
    else if ( test == 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Same slope, same intercepts, overlapping.'
      write ( *, '(a)' ) '  Distance should be 0.'
    else if ( test == 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Same slope, same intercepts, disjoint.'
      write ( *, '(a)' ) '  Distance should be sqrt(45)=6.7082038'
    else if ( test == 4 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Different slopes, intersecting.'
      write ( *, '(a)' ) '  Distance should be 0.'
    else if ( test == 5 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Different slopes, not intersecting.'
    else if ( test == 6 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Simple problem.'
      write ( *, '(a)' ) '  Distance should be 0'
    else if ( test == 7 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Same data, translated by 50.'
      write ( *, '(a)' ) '  Distance should be 0'
    else if ( test == 8 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Diagonal and horizontal.'
      write ( *, '(a)' ) '  Distance should be sqrt(2500/2)=35.355339'
    else if ( test == 9 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Same data, except first segment extended.'
      write ( *, '(a)' ) '  Distance should be sqrt(2500/2)=35.355339'
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,2g14.6)' ) '  P1 = ', p1(1:dim_num)
    write ( *, '(a,2g14.6)' ) '  P2 = ', p2(1:dim_num)
    write ( *, '(a,2g14.6)' ) '  Q1 = ', q1(1:dim_num)
    write ( *, '(a,2g14.6)' ) '  Q2 = ', q2(1:dim_num)
    write ( *, '(a,g14.6)' ) '  Distance([P1,P2],[Q1,Q2]) = ', dist

  end do

  return
end
subroutine test043 ( )

!*****************************************************************************80
!
!! TEST043 tests SEGMENTS_DIST_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) dist
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) q1(dim_num)
  real ( kind = 8 ) q2(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST043'
  write ( *, '(a)' ) '  SEGMENTS_DIST_3D computes the distance between'
  write ( *, '(a)' ) '  line segments in 3D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Case   Computed    True'
  write ( *, '(a)' ) ' '
!
!  Case 1, parallel, not coincident.
!
!  LS1: (2,3,0) + t * (2,1,0) for t = 0 to 3.
!  LS2: (11,6,4) + t * (2,1,0) for t = 0 to 3.
!  Distance is 5.
!
  p1(1:dim_num) = (/  2.0D+00, 3.0D+00, 0.0D+00 /)
  p2(1:dim_num) = (/  8.0D+00, 6.0D+00, 0.0D+00 /)
  q1(1:dim_num) = (/ 11.0D+00, 6.0D+00, 4.0D+00 /)
  q2(1:dim_num) = (/ 17.0D+00, 9.0D+00, 4.0D+00 /)

  call segments_dist_3d ( p1, p2, q1, q2, dist )

  write ( *, '(2x,i8,2g14.6)' ) 1, dist, 5.0D+00
!
!  Case 2, parallel, coincident, overlapping.
!
!  (1,2,3) + t * ( 1,-1,2)
!  LS1: t = 0 to t = 3.
!  Distance is 0.
!
  p1(1:dim_num) = (/  1.0D+00,  2.0D+00,  3.0D+00 /)
  p2(1:dim_num) = (/  4.0D+00, -1.0D+00,  9.0D+00 /)
  q1(1:dim_num) = (/  3.0D+00,  0.0D+00,  7.0D+00 /)
  q2(1:dim_num) = (/  6.0D+00, -3.0D+00, 13.0D+00 /)

  call segments_dist_3d ( p1, p2, q1, q2, dist )

  write ( *, '(2x,i8,2g14.6)' ) 2, dist, 0.0D+00
!
!  Case 3, parallel, coincident, disjoint.
!
!  LS1: (3,4,5) + t * ( 2,2,1) for 0 <= t <= 2.
!  LS2: (3,4,5) + t * ( 2,2,1) for 3 <= t <= 5.
!  Distance = 3.
!
  p1(1:dim_num) = (/  3.0D+00,  4.0D+00,  5.0D+00 /)
  p2(1:dim_num) = (/  7.0D+00,  8.0D+00,  7.0D+00 /)
  q1(1:dim_num) = (/  9.0D+00, 10.0D+00,  8.0D+00 /)
  q2(1:dim_num) = (/ 13.0D+00, 14.0D+00, 10.0D+00 /)

  call segments_dist_3d ( p1, p2, q1, q2, dist )

  write ( *, '(2x,i8,2g14.6)' ) 3, dist, 3.0D+00
!
!  Case 4, nonparallel, could intersect, and does intersect.
!
!  L1: (1,1,1) + t * (0,1,2)
!  L2: (0,2,3) + t * (1,0,0)
!  intersect at (1,2,3)
!  Distance is 0.
!
  p1(1:dim_num) = (/  1.0D+00,  1.0D+00,  1.0D+00 /)
  p2(1:dim_num) = (/  1.0D+00,  4.0D+00,  7.0D+00 /)
  q1(1:dim_num) = (/  0.0D+00,  2.0D+00,  3.0D+00 /)
  q2(1:dim_num) = (/  5.0D+00,  2.0D+00,  3.0D+00 /)

  call segments_dist_3d ( p1, p2, q1, q2, dist )

  write ( *, '(2x,i8,2g14.6)' ) 4, dist, 0.0D+00
!
!  Case 5, nonparallel, could intersect, and does not intersect.
!
!  L1: (1,1,1) + t * (0,1,2)
!  L2: (0,2,3) + t * (1,0,0)
!  lines intersect at (1,2,3), line segments do not intersect.
!  Distance is 1.0D+00
!
  p1(1:dim_num) = (/  1.0D+00,  1.0D+00,  1.0D+00 /)
  p2(1:dim_num) = (/  1.0D+00,  4.0D+00,  7.0D+00 /)
  q1(1:dim_num) = (/  0.0D+00,  2.0D+00,  3.0D+00 /)
  q2(1:dim_num) = (/ -5.0D+00,  2.0D+00,  3.0D+00 /)

  call segments_dist_3d ( p1, p2, q1, q2, dist )

  write ( *, '(2x,i8,2g14.6)' ) 5, dist, 1.0D+00
!
!  Case 6, nonparallel, can not intersect, "end-to-end".
!
!  L1: (2,2,1) + t * (0,1,2)  0 <= t <= 5
!  L2: (0,0,0) + t * (-1,-1,-1) 0 <= t <= 5
!  Distance is 3.
!
  p1(1:dim_num) = (/  2.0D+00,  2.0D+00,  1.0D+00 /)
  p2(1:dim_num) = (/  2.0D+00,  7.0D+00,  11.0D+00 /)
  q1(1:dim_num) = (/  0.0D+00,  0.0D+00,  0.0D+00 /)
  q2(1:dim_num) = (/ -5.0D+00, -5.0D+00, -5.0D+00 /)

  call segments_dist_3d ( p1, p2, q1, q2, dist )

  write ( *, '(2x,i8,2g14.6)' ) 6, dist, 3.0D+00
!
!  Case 7, nonparallel, can not intersect, "end-to-mid".
!
!  L1: (1,1,1) + t * (0,1,2) 0 <= t <= 5
!  L2: (0,4,7) + t * (-1,0,0) 0 <= t <= 5
!  Distance is 1.
!
  p1(1:dim_num) = (/  1.0D+00,  1.0D+00,  1.0D+00 /)
  p2(1:dim_num) = (/  1.0D+00,  6.0D+00, 11.0D+00 /)
  q1(1:dim_num) = (/  0.0D+00,  4.0D+00,  7.0D+00 /)
  q2(1:dim_num) = (/ -5.0D+00,  4.0D+00,  7.0D+00 /)

  call segments_dist_3d ( p1, p2, q1, q2, dist )

  write ( *, '(2x,i8,2g14.6)' ) 7, dist, 1.0D+00
!
!  Case 8, nonparallel, can not intersect, "mid-to-mid".
!
!  L1: (0,5,10) + t * (1,-1,0) 0 <= t <= 5
!  L2: (0,0,0) + t * (1,1,0) 0 <= t <= 6
!  Distance = 10.
!
  p1(1:dim_num) = (/  0.0D+00,  5.0D+00, 10.0D+00 /)
  p2(1:dim_num) = (/  5.0D+00,  0.0D+00, 10.0D+00 /)
  q1(1:dim_num) = (/  0.0D+00,  0.0D+00,  0.0D+00 /)
  q2(1:dim_num) = (/  6.0D+00,  6.0D+00,  0.0D+00 /)

  call segments_dist_3d ( p1, p2, q1, q2, dist )

  write ( *, '(2x,i8,2g14.6)' ) 8, dist, 10.0D+00
!
!  Case 9, nonparallel, can not intersect, "mid-to-end".
!
!  L1: (-2,0,0) + t * (1,0,0) 0 <= t <= 12
!  L2: (-2,8,1) + t * (9,-4,-1) 0 <= t <= 1
!  Distance = 4.
!
  p1(1:dim_num) = (/ -2.0D+00,  0.0D+00,  0.0D+00 /)
  p2(1:dim_num) = (/ 10.0D+00,  0.0D+00,  0.0D+00 /)
  q1(1:dim_num) = (/ -2.0D+00,  8.0D+00,  1.0D+00 /)
  q2(1:dim_num) = (/  7.0D+00,  4.0D+00,  0.0D+00 /)

  call segments_dist_3d ( p1, p2, q1, q2, dist )

  write ( *, '(2x,i8,2g14.6)' ) 9, dist, 4.0D+00

  return
end
subroutine test044 ( )

!*****************************************************************************80
!
!! TEST044 tests SEGMENTS_INT_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 7

  real ( kind = 8 ) dist
  integer ( kind = 4 ) test
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) q1
  real ( kind = 8 ), dimension ( test_num ) :: q1_test = &
    (/ -1.0D+00, 3.0D+00, 1.0D+00, 0.5D+00, 0.25D+00, 0.5D+00, 2.0D+00 /)
  real ( kind = 8 ) q2
  real ( kind = 8 ), dimension ( test_num ) :: q2_test = &
    (/  1.0D+00, 2.0D+00, 2.0D+00, -3.0D+00, 0.50D+00, 0.5D+00, 2.0D+00 /)
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST044'
  write ( *, '(a)' ) '  SEGMENTS_INT_1D determines the intersection [R1,R2]'
  write ( *, '(a)' ) '  of line segments [P1,P2] and [Q1,Q2] in 1D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DIST is negative for overlap,'
  write ( *, '(a)' ) '  0 for point intersection,'
  write ( *, '(a)' ) '  positive if there is no overlap.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Test      P1        P2        Q1        Q2        R1        R2        DIST'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p1 = -1.0D+00
    p2 = 1.0D+00
    q1 = q1_test(test)
    q2 = q2_test(test)

    call segments_int_1d ( p1, p2, q1, q2, dist, r1, r2 )

    write ( *, '(2x,i4,3(2x,f8.4,2x,f8.4),2x,f8.4)' ) &
      test, p1, p2, q1, q2, r1, r2, dist

  end do

  return
end
subroutine test045 ( )

!*****************************************************************************80
!
!! TEST045 tests SEGMENTS_INT_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) flag
  real ( kind = 8 ), dimension ( dim_num ) :: p1 = (/ -1.0D+00, 3.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p2 = (/  1.0D+00, 1.0D+00 /)
  real ( kind = 8 ) q1(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: q1_test = reshape ( (/ &
    -1.0D+00,  1.0D+00, &
     3.0D+00, -1.0D+00, &
     0.0D+00,  0.0D+00, &
     1.0D+00,  2.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) q2(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: q2_test = reshape ( (/ &
    1.0D+00, -1.0D+00, &
    2.0D+00,  0.0D+00, &
    0.0D+00,  9.0D+00, &
    3.0D+00,  2.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) r(dim_num)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST045'
  write ( *, '(a)' ) '  SEGMENTS_INT_2D searches for an intersection of two'
  write ( *, '(a)' ) '  line segments in 2D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  All tests use the same line segment 1:'
  write ( *, '(a,2g14.6)' ) '  P1 = ', p1(1:dim_num)
  write ( *, '(a,2g14.6)' ) '  P2 = ', p2(1:dim_num)

  do test = 1, test_num

    q1(1:dim_num) = q1_test(1:dim_num,test)
    q2(1:dim_num) = q2_test(1:dim_num,test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,2g14.6)' ) '  Q1 = ', q1(1:dim_num)
    write ( *, '(a,2g14.6)' ) '  Q2 = ', q2(1:dim_num)

    call segments_int_2d ( p1, p2, q1, q2, flag, r )

    if ( flag == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The line segments do not intersect.'
    else if ( flag == 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The line segments intersect at:'
      write ( *, '(2x,2g14.6)' ) r(1:dim_num)
    end if

  end do

  return
end
subroutine test046 ( )

!*****************************************************************************80
!
!! TEST046 tests MINABS.
!
!  Discussion:
!
!    Case 1: the three points lie on a straight line.
!    (XMIN=9,YMIN=2).
!
!    Case 2: the three points straddle a minimum.
!    (XMIN=7, YMIN=2).
!
!    Case 3: the three points straddle a maximum.
!    (XMIN=2, YMIN=5).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  integer ( kind = 4 ) test
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) ymin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST046'
  write ( *, '(a)' ) '  MINABS finds the minimum of a function'
  write ( *, '(a)' ) '  F(X) = a * ABS ( X ) + B'
  write ( *, '(a)' ) '  within an interval, given three data points.'

  do test = 1, test_num

    if ( test == 1 ) then

      x1 = 14.0D+00
      y1 = 7.0D+00
      x2 = 9.0D+00
      y2 = 2.0D+00
      x3 = 12.0D+00
      y3 = 5.0D+00

    else if ( test == 2 ) then

      x1 = 3.0D+00
      y1 = 6.0D+00
      x2 = 12.0D+00
      y2 = 7.0D+00
      x3 = 9.0D+00
      y3 = 4.0D+00

    else if ( test == 3 ) then

      x1 = 11.0D+00
      y1 = 6.0D+00
      x2 = 6.0D+00
      y2 = 9.0D+00
      x3 = 2.0D+00
      y3 = 5.0D+00

    end if

    call minabs ( x1, y1, x2, y2, x3, y3, xmin, ymin )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The points lie on a straight line.'
    write ( *, '(a,2g14.6)' ) '  XMIN, YMIN = ', xmin, ymin

  end do

  return
end
subroutine test047 ( )

!*****************************************************************************80
!
!! TEST047 tests MINQUAD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) ymin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST047'
  write ( *, '(a)' ) '  MINQUAD finds the minimum of a function'
  write ( *, '(a)' ) '  F(X) = A * X * X + B * X + C'
  write ( *, '(a)' ) '  within an interval, given three data points.'
!
!  Case 1: a minimum is in the interval.
!  y = ( x - 1 )**2 + 4
!
  x1 = 0.0D+00
  y1 = ( x1 - 1.0D+00 )**2 + 4.0D+00

  x2 = 2.0D+00
  y2 = ( x2 - 1.0D+00 )**2 + 4.0D+00

  x3 = 3.0D+00
  y3 = ( x3 - 1.0D+00 )**2 + 4.0D+00

  call minquad ( x1, y1, x2, y2, x3, y3, xmin, ymin )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The minimum lies in the interval.'
  write ( *, '(a,2g14.6)' ) '  X1,   Y1   = ', x1, y1
  write ( *, '(a,2g14.6)' ) '  X2,   Y2   = ', x2, y2
  write ( *, '(a,2g14.6)' ) '  X3,   Y3   = ', x3, y3
  write ( *, '(a,2g14.6)' ) '  XMIN, YMIN = ', xmin, ymin
!
!  Case 2: the minimum is to the left of the interval.
!  y = ( x - 1 )**2 + 4
!
  x1 = 2.0D+00
  y1 = ( x1 - 1.0D+00 )**2 + 4.0D+00

  x2 = 4.0D+00
  y2 = ( x2 - 1.0D+00 )**2 + 4.0D+00

  x3 = 5.0D+00
  y3 = ( x3 - 1.0D+00 )**2 + 4.0D+00

  call minquad ( x1, y1, x2, y2, x3, y3, xmin, ymin )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The minimum is to the left of the interval'
  write ( *, '(a,2g14.6)' ) '  X1,   Y1   = ', x1, y1
  write ( *, '(a,2g14.6)' ) '  X2,   Y2   = ', x2, y2
  write ( *, '(a,2g14.6)' ) '  X3,   Y3   = ', x3, y3
  write ( *, '(a,2g14.6)' ) '  XMIN, YMIN = ', xmin, ymin
!
!  Case 3: the function is flat.
!
  x1 = 11.0D+00
  y1 = 6.0D+00

  x2 = 6.0D+00
  y2 = 6.0D+00

  x3 = 2.0D+00
  y3 = 6.0D+00

  call minquad ( x1, y1, x2, y2, x3, y3, xmin, ymin )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The function is flat.'
  write ( *, '(a,2g14.6)' ) '  X1,   Y1   = ', x1, y1
  write ( *, '(a,2g14.6)' ) '  X2,   Y2   = ', x2, y2
  write ( *, '(a,2g14.6)' ) '  X3,   Y3   = ', x3, y3
  write ( *, '(a,2g14.6)' ) '  XMIN, YMIN = ', xmin, ymin
!
!  Case 4: the function has a maximum.
!  y = - ( x - 1 )**2 + 4
!
  x1 = 0.0D+00
  y1 = - ( x1 - 1.0D+00 )**2 + 4.0D+00

  x2 = 2.0D+00
  y2 = - ( x2 - 1.0D+00 )**2 + 4.0D+00

  x3 = 3.0D+00
  y3 = - ( x3 - 1.0D+00 )**2 + 4.0D+00

  call minquad ( x1, y1, x2, y2, x3, y3, xmin, ymin )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The function has a maximum.'
  write ( *, '(a,2g14.6)' ) '  X1,   Y1   = ', x1, y1
  write ( *, '(a,2g14.6)' ) '  X2,   Y2   = ', x2, y2
  write ( *, '(a,2g14.6)' ) '  X3,   Y3   = ', x3, y3
  write ( *, '(a,2g14.6)' ) '  XMIN, YMIN = ', xmin, ymin

  return
end
subroutine test0475 ( )

!*****************************************************************************80
!
!! TEST0475 tests OCTAHEDRON_SIZE_3D and OCTAHEDRON_SHAPE_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) area
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ), allocatable, dimension ( : ) ::  face_order
  integer ( kind = 4 ) face_order_max
  integer ( kind = 4 ), allocatable, dimension ( :, : ) ::  face_point
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) normal(dim_num)
  integer ( kind = 4 ) point_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_coord
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: v
  real ( kind = 8 ) vave(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0475'
  write ( *, '(a)' ) '  For the octahedron:'
  write ( *, '(a)' ) '  OCTAHEDRON_SIZE_3D returns dimension information;'
  write ( *, '(a)' ) '  OCTAHEDRON_SHAPE_3D returns face and order information.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We will use this information to compute the'
  write ( *, '(a)' ) '  areas and centers of each face.'

  call octahedron_size_3d ( point_num, edge_num, face_num, face_order_max )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of vertices = ', point_num
  write ( *, '(a,i8)' ) '  Number of edges =    ', edge_num
  write ( *, '(a,i8)' ) '  Number of faces =    ', face_num
  write ( *, '(a,i8)' ) '  Maximum face order = ', face_order_max

  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )
  allocate ( point_coord(1:dim_num,1:point_num) )
  allocate ( v(1:3,1:face_order_max) )

  call octahedron_shape_3d ( point_num, face_num, face_order_max, point_coord, &
    face_order, face_point )
!
!  Compute the area of each face.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Face  Order  Area'
  write ( *, '(a)' ) ' '

  do i = 1, face_num

    do j = 1, face_order(i)
      k = face_point(j,i)
      v(1:3,j) = point_coord(1:dim_num,k)
    end do

    call polygon_area_3d ( face_order(i), v, area, normal )

    write ( *, '(2x,i8,i7,f8.4)' ) i, face_order(i), area

  end do
!
!  Find the center of each face.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Face  Center'
  write ( *, '(a)' ) ' '

  do i = 1, face_num

    vave(1:dim_num) = 0.0D+00

    do j = 1, face_order(i)
      k = face_point(j,i)
      vave(1:dim_num) = vave(1:dim_num) + point_coord(1:dim_num,k)
    end do

    vave(1:dim_num) = vave(1:dim_num) / real ( face_order(i), kind = 8 )

    write ( *, '(2x,i8,3f8.4)' ) i, vave(1:dim_num)

  end do

  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )
  deallocate ( v )

  return
end
subroutine test0477 ( )

!*****************************************************************************80
!
!! TEST0477 tests PARALLELOGRAM_AREA_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 May 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ) :: p(2,4) = reshape  ( (/ &
    2.0D+00, 7.0D+00, &
    5.0D+00, 7.0D+00, &
    6.0D+00, 9.0D+00, &
    3.0D+00, 9.0D+00 &
    /), (/ 2, 4 /) ) 

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0477'
  write ( *, '(a)' ) '  PARALLELOGRAM_AREA_2D finds the area of a'
  write ( *, '(a)' ) '  parallelogram in 2D.'

  call r8mat_transpose_print ( 2, 4, p, '  Vertices:' )

  call parallelogram_area_2d ( p, area )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  AREA = ', area

  return
end
subroutine test0478 ( )

!*****************************************************************************80
!
!! TEST0478 tests PARALLELOGRAM_AREA_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 May 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ) :: p(3,4) = reshape  ( (/ &
    1.0D+00,       2.0D+00,       3.0D+00, &
    2.4142137D+00, 3.4142137D+00, 3.0D+00, &
    1.7071068D+00, 2.7071068D+00, 4.0D+00, &
    0.2928931D+00, 0.2928931D+00, 4.0D+00  &
  /), (/ 3, 4 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0478'
  write ( *, '(a)' ) '  PARALLELOGRAM_AREA_3D finds the area of a'
  write ( *, '(a)' ) '  parallelogram in 3D.'

  call r8mat_transpose_print ( 3, 4, p, '  Vertices:' )

  call parallelogram_area_3d ( p, area )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  AREA = ', area

  return
end
subroutine test048 ( )

!*****************************************************************************80
!
!! TEST048 tests PARALLELOGRAM_CONTAINS_POINT_2D.
!
!  Discussion:
!
!    The four points are In, Out, Out, and Out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  logical inside
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = &
    reshape ( (/ &
    1.0D+00,  0.5D+00, &
    2.0D+00,  0.0D+00, &
    0.5D+00, -0.1D+00, &
    0.1D+00,  0.5D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension(dim_num) :: p1 = (/ 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p2 = (/ 1.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p3 = (/ 1.0D+00, 1.0D+00 /)
  logical parallelogram_contains_point_2d
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST048'
  write ( *, '(a)' ) '  PARALLELOGRAM_CONTAINS_POINT_2D determines if a point '
  write ( *, '(a)' ) '  is within a parallelogram in 2D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Parallelogram defined by P2-P1, P3-P1:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  P1 = ', p1(1:dim_num)
  write ( *, '(a,2g14.6)' ) '  P2 = ', p2(1:dim_num)
  write ( *, '(a,2g14.6)' ) '  P3 = ', p3(1:dim_num)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       P     Inside?'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    inside = parallelogram_contains_point_2d ( p1, p2, p3, p )

    write ( *, '(2x,2g14.6,2x,l1)' ) p(1:dim_num), inside

  end do

  return
end
subroutine test0485 ( )

!*****************************************************************************80
!
!! TEST0485 tests PARALLELOGRAM_CONTAINS_POINT_2D
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 51
  integer ( kind = 4 ), parameter :: dim_num = 2

  logical parallelogram_contains_point_2d
  character dot(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: p1 = (/ 0.2D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p2 = (/ 0.4D+00, 0.6D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p3 = (/ 0.6D+00, 0.4D+00 /)
  real ( kind = 8 ), parameter :: xhi =  1.0D+00
  real ( kind = 8 ), parameter :: xlo =  0.0D+00
  real ( kind = 8 ), parameter :: yhi =  1.0D+00
  real ( kind = 8 ), parameter :: ylo =  0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0485'
  write ( *, '(a)' ) '  PARALLELOGRAM_CONTAINS_POINT_2D reports if a'
  write ( *, '(a)' ) '  parallelogram contains a point.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We will call the function repeatedly, and draw'
  write ( *, '(a)' ) '  a sketch of the unit square.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Parallelogram defined by P2-P1, P3-P1:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  P1 = ', p1(1:dim_num)
  write ( *, '(a,2g14.6)' ) '  P2 = ', p2(1:dim_num)
  write ( *, '(a,2g14.6)' ) '  P3 = ', p3(1:dim_num)
  write ( *, '(a)' ) ' '

  do i = 1, n

    p(2) = ( real ( n - i,     kind = 8 ) * yhi   &
           + real (     i - 1, kind = 8 ) * ylo ) &
           / real ( n     - 1, kind = 8 )

    do j = 1, n

      p(1) = ( real ( n - j,     kind = 8 ) * xlo   &
             + real (     j - 1, kind = 8 ) * xhi ) &
             / real ( n     - 1, kind = 8 )

      if ( parallelogram_contains_point_2d ( p1, p2, p3, p ) ) then
        dot(j) = '*'
      else
        dot(j) = '-'
      end if

    end do
    write ( *, '(2x,51a1)' ) dot(1:n)
  end do

  return
end
subroutine test049 ( )

!*****************************************************************************80
!
!! TEST049 tests PARALLELOGRAM_CONTAINS_POINT_3D.
!
!  Discussion:
!
!    The points are In, Out, Out, Out, Out
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 5

  logical inside
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = &
    reshape ( (/ &
    1.0D+00,  1.0D+00,  0.5D+00, &
    3.0D+00,  3.0D+00,  0.0D+00, &
    0.5D+00,  0.5D+00, -0.1D+00, &
    0.1D+00,  0.1D+00,  0.5D+00, &
    1.5D+00,  1.6D+00,  0.5D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension(dim_num) :: p1 = (/ &
    0.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p2 = (/ &
    2.0D+00, 2.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p3 = (/ &
    1.0D+00, 1.0D+00, 1.0D+00 /)
  logical parallelogram_contains_point_3d
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST049'
  write ( *, '(a)' ) '  PARALLELOGRAM_CONTAINS_POINT_3D determines if a point '
  write ( *, '(a)' ) '  is within a parallelogram in 3D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Parallelogram defined by P2-P1, P3-P1:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,3g14.6)' ) '  P1 = ', p1(1:dim_num)
  write ( *, '(a,3g14.6)' ) '  P2 = ', p2(1:dim_num)
  write ( *, '(a,3g14.6)' ) '  P3 = ', p3(1:dim_num)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           P           Inside?'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    inside = parallelogram_contains_point_3d ( p1, p2, p3, p )

    write ( *, '(2x,3g14.6,2x,l1)' ) p(1:dim_num), inside

  end do

  return
end
subroutine test0493 ( )

!*****************************************************************************80
!
!! TEST0493 tests PARABOLA_EX and PARABOLA_EX2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) ymin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0493'
  write ( *, '(a)' ) '  PARABOLA_EX finds the extreme value of a parabola'
  write ( *, '(a)' ) '  determined by three points.'
  write ( *, '(a)' ) '  PARABOLA_EX2 finds the extreme value of a parabola'
  write ( *, '(a)' ) '  determined by three points.'

  a =  2.0D+00
  b = -4.0D+00
  c = 10.0D+00

  x1 = 1.0D+00
  y1 = a * x1 * x1 + b * x1 + c
  x2 = 2.0D+00
  y2 = a * x2 * x2 + b * x2 + c
  x3 = 3.0D+00
  y3 = a * x3 * x3 + b * x3 + c

  write ( *, '(a)' ) ' '
  write ( *, '(a,3g14.6)' ) '  Parabolic coefficients (A,B,C) = ', a, b, c
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, Y data'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  X1, Y1 = ', x1, y1
  write ( *, '(a,2g14.6)' ) '  X2, Y2 = ', x2, y2
  write ( *, '(a,2g14.6)' ) '  X3, Y3 = ', x3, y3 

  a = 0.0D+00
  b = 0.0D+00
  c = 0.0D+00

  call parabola_ex ( x1, y1, x2, y2, x3, y3, xmin, ymin, ierror )

  if ( ierror == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,2g14.6)' ) '  PARABOLA_EX returns (XMIN,YMIN) = ', xmin, ymin
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  PARABOLA_EX returns error code ', ierror
  end if

  call parabola_ex2 ( x1, y1, x2, y2, x3, y3, xmin, ymin, a, b, c, ierror )

  if ( ierror == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,2g14.6)' ) '  PARABOLA_EX2 returns (XMIN,YMIN) = ', xmin, ymin
    write ( *, '(a,3g14.6)' ) '  and (A,B,C) = ', a, b, c
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  PARABOLA_EX2 returns error code ', ierror
  end if

  return
end
subroutine test0495 ( )

!*****************************************************************************80
!
!! TEST0495 tests PARALLELEPIPED_POINT_DIST_3D.
!
!  Discussion:
!
!    The points tested are:
!
!    1: Center of box.
!    2: The middle of a face.
!    3: The middle of an edge.
!    4: A corner.
!    5: Close to a face.
!    6: Close to an edge.
!    7: Close to a corner.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 7

  real ( kind = 8 ) dist
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = &
    reshape ( (/ &
    1.0D+00,  4.0D+00,  0.5D+00, &
    1.0D+00,  0.0D+00,  0.5D+00, &
    0.0D+00,  4.0D+00,  1.0D+00, &
    2.0D+00,  8.0D+00,  1.0D+00, &
   -0.5D+00,  4.0D+00,  0.5D+00, &
    1.0D+00, -1.0D+00, -1.0D+00, &
    3.0D+00,  9.0D+00,  2.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num ) :: p1 = (/ &
    0.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p2 = (/ &
    2.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p3 = (/ &
    0.0D+00, 8.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p4 = (/ &
    0.0D+00, 0.0D+00, 1.0D+00 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0495'
  write ( *, '(a)' ) '  PARALLELEPIPED_POINT_DIST_3D computes the distance '
  write ( *, '(a)' ) '  from a point to a box (parallelipiped) in 3D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The 4 box corners that are specified:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,3f10.2)' ) p1(1:dim_num)
  write ( *, '(2x,3f10.2)' ) p2(1:dim_num)
  write ( *, '(2x,3f10.2)' ) p3(1:dim_num)
  write ( *, '(2x,3f10.2)' ) p4(1:dim_num)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Test          P                  Distance to box'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call parallelepiped_point_dist_3d ( p1, p2, p3, p4, p, dist )

    write ( *, '(2x,i3,2x,3f10.2,2x,g14.6)' ) test, p(1:dim_num), dist

  end do

  return
end
subroutine test050 ( )

!*****************************************************************************80
!
!! TEST050 tests PLANE_EXP_NORMAL_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) normal(dim_num)
  real ( kind = 8 ), dimension ( dim_num) :: p1 = (/ &
    -10.56D+00, -10.56D+00, 78.09D+00 /)
  real ( kind = 8 ), dimension ( dim_num) :: p2 = (/ &
     44.66D+00, -65.77D+00,  0.00D+00 /)
  real ( kind = 8 ), dimension ( dim_num) :: p3 = (/ &
     44.66D+00,  44.66D+00,  0.00D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST050'
  write ( *, '(a)' ) '  PLANE_EXP_NORMAL_3D finds the normal '
  write ( *, '(a)' ) '  to a plane.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,3g14.6)' ) '  P1: ', p1(1:dim_num)
  write ( *, '(a,3g14.6)' ) '  P2: ', p2(1:dim_num)
  write ( *, '(a,3g14.6)' ) '  P3: ', p3(1:dim_num)

  call plane_exp_normal_3d ( p1, p2, p3, normal )

  call r8vec_print ( dim_num, normal, '  The normal vector:' )

  return
end
subroutine test051 ( )

!*****************************************************************************80
!
!! TEST051 tests PLANE_EXP2IMP_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ), dimension ( dim_num ) :: p1
  real ( kind = 8 ), dimension ( dim_num,test_num ) :: p1_test = reshape ( (/ &
     -1.0D+00, 0.0D+00, -1.0D+00, &
    -16.0D+00,  2.0D+00, 4.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num ) :: p2
  real ( kind = 8 ), dimension ( dim_num,test_num ) :: p2_test = reshape ( (/ &
    -4.0D+00, 0.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num ) :: p3
  real ( kind = 8 ), dimension ( dim_num,test_num ) :: p3_test = reshape ( (/ &
    -20.0D+00, 2.0D+00,  4.0D+00, &
      4.0D+00, -2.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST051'
  write ( *, '(a)' ) '  PLANE_EXP2IMP_3D puts a plane defined by '
  write ( *, '(a)' ) '  3 points into A*X+B*Y+C*Z+D = 0 form.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 1: Correct answer is a multiple of 1, 2, 3, 4.'
  write ( *, '(a)' ) '  Test 2: Correct answer is a multiple of 1, 2, 3, 0.'

  do test = 1, test_num

    p1(1:dim_num) = p1_test(1:dim_num,test)
    p2(1:dim_num) = p2_test(1:dim_num,test)
    p3(1:dim_num) = p3_test(1:dim_num,test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,3g14.6)' ) '  P1: ', p1(1:dim_num)
    write ( *, '(a,3g14.6)' ) '  P2: ', p2(1:dim_num)
    write ( *, '(a,3g14.6)' ) '  P3: ', p3(1:dim_num)

    call plane_exp2imp_3d ( p1, p2, p3, a, b, c, d )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (A,B,C,D)= '
    write ( *, '(2x,4g14.6)' ) a, b, c, d
    write ( *, '(a)' ) ' '

  end do
 
  return
end
subroutine test052 ( )

!*****************************************************************************80
!
!! TEST052 tests PLANE_EXP2NORMAL_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) normal(dim_num)
  real ( kind = 8 ), dimension ( dim_num ) :: p1
  real ( kind = 8 ), dimension ( dim_num,test_num ) :: p1_test = reshape ( (/ &
     -1.0D+00, 0.0D+00, -1.0D+00, &
    -16.0D+00,  2.0D+00, 4.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num ) :: p2
  real ( kind = 8 ), dimension ( dim_num,test_num ) :: p2_test = reshape ( (/ &
    -4.0D+00, 0.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num ) :: p3
  real ( kind = 8 ), dimension ( dim_num,test_num ) :: p3_test = reshape ( (/ &
    -20.0D+00, 2.0D+00,  4.0D+00, &
      4.0D+00, -2.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) pp(dim_num)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST052'
  write ( *, '(a)' ) '  PLANE_EXP2NORMAL_3D puts a plane defined by '
  write ( *, '(a)' ) '  3 points into point, normal form.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p1(1:dim_num) = p1_test(1:dim_num,test)
    p2(1:dim_num) = p2_test(1:dim_num,test)
    p3(1:dim_num) = p3_test(1:dim_num,test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,3g14.6)' ) '  P1: ', p1(1:dim_num)
    write ( *, '(a,3g14.6)' ) '  P2: ', p2(1:dim_num)
    write ( *, '(a,3g14.6)' ) '  P3: ', p3(1:dim_num)

    call plane_exp2normal_3d ( p1, p2, p3, pp, normal )
 
    call r8vec_print ( dim_num, pp, '  The point PP:' )

    call r8vec_print ( dim_num, normal, '  Normal vector:' )

  end do
 
  return
end
subroutine test053 ( )

!*****************************************************************************80
!
!! TEST053 tests PLANE_EXP_PROJECT_3D.
!
!  Discussion:
!
!    1: Projection   is ( 0, 0.5, 0.5 ), IVIS is 3.
!    2: Projection   is ( 4, 5, -8 ), IVIS is 2.
!    3: Projection   is ( 0.33, 0.33, 0.33), IVIS is 1.
!    4: "Projection" is ( 0, 0, 0 ), IVIS is 0.
!    5: Projection   is ( 1, 0, 0 ), IVIS is -1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 5

  integer ( kind = 4 ) ivis(test_num)
  real ( kind = 8 ), dimension ( dim_num ) :: p1 = (/ &
    1.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p2 = (/ &
    0.0D+00, 1.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p3 = (/ &
    0.0D+00, 0.0D+00, 1.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: pf = (/ &
    0.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension(dim_num,test_num) :: po = reshape ( (/ &
     0.00D+00,  2.00D+00,  2.00D+00, &
     4.00D+00,  5.00D+00, -8.00D+00, &
     0.25D+00,  0.25D+00,  0.25D+00, &
     5.00D+00, -2.00D+00, -3.00D+00, &
    -2.00D+00,  0.00D+00,  0.00D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) pp(dim_num,test_num)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST053'
  write ( *, '(a)' ) '  PLANE_EXP_PROJECT_3D projects a point through'
  write ( *, '(a)' ) '  a focus point into a plane.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      PO      PP      IVIS'
  write ( *, '(a)' ) ' '
 
  call plane_exp_project_3d ( p1, p2, p3, pf, test_num, po, pp, ivis )
 
  do test = 1, test_num
    write ( *, '(2x,6g12.4,i4)' ) po(1:dim_num,test), &
      pp(1:dim_num,test), ivis(test_num)
  end do
 
  return
end
subroutine test054 ( )

!*****************************************************************************80
!
!! TEST054 tests PLANE_IMP2EXP_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) p3(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST054'
  write ( *, '(a)' ) '  PLANE_IMP2EXP_3D converts a plane in implicit'
  write ( *, '(a)' ) '  (A,B,C,D) form to explicit form.'

  a = 1.0D+00
  b = -2.0D+00
  c = -3.0D+00
  d = 6.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    (A,B,C,D) = '
  write ( *, '(2x,4g14.6)' ) a, b, c, d

  call plane_imp2exp_3d ( a, b, c, d, p1, p2, p3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,3g14.6)' ) '  P1: ', p1(1:dim_num)
  write ( *, '(a,3g14.6)' ) '  P2: ', p2(1:dim_num)
  write ( *, '(a,3g14.6)' ) '  P3: ', p3(1:dim_num)

  return
end
subroutine test055 ( )

!*****************************************************************************80
!
!! TEST055 tests PLANE_IMP2NORMAL_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) normal(dim_num)
  real ( kind = 8 ) pp(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST055'
  write ( *, '(a)' ) '  PLANE_IMP2NORMAL_3D converts a plane in implicit'
  write ( *, '(a)' ) '  (A,B,C,D) form to point, normal form.'

  a = 1.0D+00
  b = -2.0D+00
  c = -3.0D+00
  d = 6.0D+00
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,4g14.6)' ) '    (A,B,C,D) = ', a, b, c, d

  call plane_imp2normal_3d ( a, b, c, d, pp, normal )

  call r8vec_print ( dim_num, pp, '  The point PP:' )

  call r8vec_print ( dim_num, normal, '  Normal vector:' )

  return
end
subroutine test056 ( )

!*****************************************************************************80
!
!! TEST056 tests PLANE_IMP_LINE_PAR_INT_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  logical intersect
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) x0
  real ( kind = 8 ) y0
  real ( kind = 8 ) z0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST056'
  write ( *, '(a)' ) '  PLANE_IMP_LINE_PAR_INT_3D finds the '
  write ( *, '(a)' ) '  intersection of an implicit plane and'
  write ( *, '(a)' ) '  a parametric line, in 3D.'

  a = 1.0D+00
  b = -2.0D+00
  c = -3.0D+00
  d = 6.0D+00
 
  f = 2.0D+00
  g = 1.0D+00
  h = 5.0D+00
  x0 = 3.0D+00
  y0 = 0.0D+00
  z0 = -7.0D+00
 
  call plane_imp_line_par_int_3d ( a, b, c, d, x0, y0, z0, f, g, h, &
    intersect, p )
 
  if ( intersect ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The plane and line intersect at '
    write ( *, '(2x,3g14.6)' ) p(1:dim_num)
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The plane and the line do not intersect.'
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Expected answer:'
  write ( *, '(a)' ) '    The plane and line intersect at '
  write ( *, '(a)' ) '    7, 2, 3.'
 
  return
end
subroutine test057 ( )

!*****************************************************************************80
!
!! TEST057 tests PLANE_IMP_SEGMENT_NEAR_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) dist
  real ( kind = 8 ) p(1:dim_num)
  real ( kind = 8 ) p1(1:dim_num)
  real ( kind = 8 ) p2(1:dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p2_test = reshape ( (/ &
    9.0D+00, 3.0D+00,  8.0D+00, &
    5.0D+00, 1.0D+00, -2.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) pn(1:dim_num)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST057'
  write ( *, '(a)' ) '  PLANE_IMP_SEGMENT_NEAR_3D finds the point'
  write ( *, '(a)' ) '  on a line segment nearest a plane.'
 
  do test = 1, test_num
 
    p1(1:dim_num) = (/ 3.0D+00, 0.0D+00, -7.0D+00 /) 
    p2(1:dim_num) = p2_test(1:dim_num,test)

    a = 1.0D+00
    b = -2.0D+00
    c = -3.0D+00
    d = 6.0D+00
 
    call plane_imp_segment_near_3d ( p1, p2, a, b, c, d, dist, p, pn )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The distance between the plane and the'
    write ( *, '(a,g14.6)' ) '  line segment is ', dist
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  A nearest point on the line segment is '
    write ( *, '(2x,3g14.6)' ) pn(1:dim_num)
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  A nearest point on the plane is '
    write ( *, '(2x,3g14.6)' ) p(1:dim_num)
 
  end do
 
  return
end
subroutine test058 ( )

!*****************************************************************************80
!
!! TEST058 tests PLANE_IMP_POINT_DIST_3D and PLANE_IMP_POINT_DIST_SIGNED_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_signed
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p_test = &
    reshape ( (/ &
    -12.0D+00, 14.0D+00,  0.0D+00, &
      7.0D+00,  8.0D+00,  9.0D+00, &
      1.0D+00,  2.0D+00, 10.0D+00, &
      0.0D+00,  0.0D+00, 12.0D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test
!
!  This is the plane Z = 10.
!
  a =    0.0D+00
  b =    0.0D+00
  c =    1.0D+00
  d = - 10.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST058'
  write ( *, '(a)' ) '  PLANE_IMP_POINT_DIST_3D computes the distance'
  write ( *, '(a)' ) '  between an implicit plane and a point in 3D;'
  write ( *, '(a)' ) '  PLANE_IMP_POINT_DIST_SIGNED 3D computes the '
  write ( *, '(a)' ) '  signed distance between an implicit plane '
  write ( *, '(a)' ) '  and a point in 3D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For all tests, we use the implicit plane with'
  write ( *, '(a,4g14.6)' ) '  (A,B,C,D) = ', a, b, c, d
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (X,Y,Z)  DISTANCE   SIGNED_DISTANCE'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call plane_imp_point_dist_3d ( a, b, c, d, p, dist )

    call plane_imp_point_dist_signed_3d ( a, b, c, d, p, dist_signed )

    write ( *, '(2x,5g14.6)' ) p(1:dim_num), dist, dist_signed

  end do

  return
end
subroutine test059 ( )

!*****************************************************************************80
!
!! TEST059 tests PLANE_IMP_TRIANGLE_NEAR_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) dist
  integer ( kind = 4 ) near_num
  real ( kind = 8 ) pn(dim_num,6)
  real ( kind = 8 ), dimension(dim_num,3) :: t
  real ( kind = 8 ), dimension(dim_num,3,test_num) :: t_test = reshape ( (/ &
     3.0D+00,  0.0D+00, -7.0D+00, &
    13.0D+00, -4.0D+00, -1.0D+00, &
     5.0D+00,  1.0D+00, -2.0D+00, &
     3.0D+00,  0.0D+00, -7.0D+00, &
    13.0D+00, -4.0D+00, -1.0D+00, &
     9.0D+00,  3.0D+00,  8.0D+00 /), (/ dim_num, 3, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST059'
  write ( *, '(a)' ) '  PLANE_IMP_TRIANGLE_NEAR_3D finds the nearest'
  write ( *, '(a)' ) '  points on an implicit plane and a triangle.'
 
  a = 1.0D+00
  b = -2.0D+00
  c = -3.0D+00
  d = 6.0D+00
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Implicit plane: A*X + B*Y + C*Z + D = 0.'
  write ( *, '(a)' ) '  A,B,C,D = '
  write ( *, '(2x,4g14.6)' ) a, b, c, d
 
  do test = 1, test_num
 
    t(1:dim_num,1:3) = t_test(1:dim_num,1:3,test)
 
    call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

    call plane_imp_triangle_near_3d ( t, a, b, c, d, dist, near_num, pn )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Triangle to plane distance is ', dist

    call r8mat_transpose_print ( dim_num, near_num, pn, '  Nearest points:' )
 
  end do
 
  return
end
subroutine test060 ( )

!*****************************************************************************80
!
!! TEST060 tests PLANE_IMP_TRIANGLE_INT_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  integer ( kind = 4 ) int_num
  real ( kind = 8 ) pint(dim_num,3)
  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ), dimension(dim_num,3,test_num) :: t_test = reshape ( (/ &
         3.0D+00,  0.0D+00, -7.0D+00, &
        13.0D+00, -4.0D+00, -1.0D+00, &
         5.0D+00,  1.0D+00, -2.0D+00, &
         3.0D+00,  0.0D+00, -7.0D+00, &
        13.0D+00, -4.0D+00, -1.0D+00, &
         9.0D+00,  3.0D+00,  8.0D+00, &
        -6.0D+00,  0.0D+00,  0.0D+00, &
         0.0D+00,  3.0D+00,  0.0D+00, &
         0.0D+00,  0.0D+00,  2.0D+00, &
        -4.0D+00,  1.0D+00,  0.0D+00, &
         0.0D+00,  6.0D+00, -2.0D+00, &
         0.0D+00,  0.0D+00,  1.0D+00 /), (/ dim_num, 3, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST060'
  write ( *, '(a)' ) '  PLANE_IMP_TRIANGLE_INT_3D finds the'
  write ( *, '(a)' ) '  intersection points of an implicit plane'
  write ( *, '(a)' ) '  and a triangle.'
 
  a = 1.0D+00
  b = -2.0D+00
  c = -3.0D+00
  d = 6.0D+00
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The implicit plane: A*X + B*Y + C*Z + D = 0.'
  write ( *, '(a,4g14.6)' ) '  A,B,C,D = ', a, b, c, d

  do test = 1, test_num
 
    t(1:dim_num,1:3) = t_test(1:dim_num,1:3,test)
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Case ', test
    write ( *, '(a)' ) ' '

    call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

    call plane_imp_triangle_int_3d ( a, b, c, d, t, int_num, pint )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of intersection points is ', int_num
    write ( *, '(a)' ) ' '

    call r8mat_transpose_print ( dim_num, int_num, pint, &
      '  Intersection points:' )
 
  end do
 
  return
end
subroutine test061 ( )

!*****************************************************************************80
!
!! TEST061 tests PLANE_NORMAL_BASIS_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 5

  real ( kind = 8 ) b(dim_num,dim_num)
  real ( kind = 8 ), dimension(dim_num) :: normal
  real ( kind = 8 ), dimension(dim_num) :: pp
  real ( kind = 8 ), dimension(dim_num) :: pq
  real ( kind = 8 ), dimension(dim_num) :: pr
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST061'
  write ( *, '(a)' ) '  PLANE_NORMAL_BASIS_3D, given a plane in'
  write ( *, '(a)' ) '  point, normal form (P,N), finds two unit'
  write ( *, '(a)' ) '  vectors Q and R that "lie" in the plane'
  write ( *, '(a)' ) '  and are mutually orthogonal.'

  do test = 1, test_num

    call r8vec_uniform_01 ( dim_num, seed, normal )
    call r8vec_uniform_01 ( dim_num, seed, pp )

    call plane_normal_basis_3d ( pp, normal, pq, pr )

    if ( test == 1 ) then

      call r8vec_print ( dim_num, pp, '  Point PP:' )
      call r8vec_print ( dim_num, normal, '  Normal vector N:' )
      call r8vec_print ( dim_num, pq, '  Vector PQ:' )
      call r8vec_print ( dim_num, pr, '  Vector PR:' )

    end if

    b(1,1) = dot_product ( normal(1:dim_num), normal(1:dim_num) )
    b(1,2) = dot_product ( normal(1:dim_num), pq(1:dim_num) )
    b(1,3) = dot_product ( normal(1:dim_num), pr(1:dim_num) )
    b(2,1) = dot_product ( pq(1:dim_num), normal(1:dim_num) )
    b(2,2) = dot_product ( pq(1:dim_num), pq(1:dim_num) )
    b(2,3) = dot_product ( pq(1:dim_num), pr(1:dim_num) )
    b(3,1) = dot_product ( pr(1:dim_num), normal(1:dim_num) )
    b(3,2) = dot_product ( pr(1:dim_num), pq(1:dim_num) )
    b(3,3) = dot_product ( pr(1:dim_num), pr(1:dim_num) )

    call r8mat_print ( dim_num, dim_num, b, '  Dot product matrix:' )

  end do

  return
end
subroutine test0615 ( )

!*****************************************************************************80
!
!! TEST0615 tests PLANE_NORMAL_LINE_EXP_INT_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) ival
  real ( kind = 8 ), dimension(dim_num) :: normal = (/ &
    1.0D+00, -2.0D+00, -3.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p1 = (/ &
    3.0D+00, 0.0D+00, -7.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p2 = (/ &
    5.0D+00, 1.0D+00, -2.0D+00 /)
  real ( kind = 8 ) pint(dim_num)
  real ( kind = 8 ), dimension ( dim_num ) :: pp = (/ &
    -1.0D+00, +1.0D+00, +1.0D+00 /)
  real ( kind = 8 ) temp

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0615'
  write ( *, '(a)' ) '  PLANE_NORMAL_LINE_EXP_INT_3D finds the '
  write ( *, '(a)' ) '  intersection of a normal plane and'
  write ( *, '(a)' ) '  an explicit line, in 3D.'

  temp = sqrt ( sum ( normal(1:dim_num)**2 ) )
  normal(1:dim_num) = normal(1:dim_num) / temp

  call r8vec_print ( dim_num, pp, '  Plane point PP:' )
  call r8vec_print ( dim_num, normal, '  Plane Normal:' )

  call r8vec_print ( dim_num, p1, '  Line point P1:' )
  call r8vec_print ( dim_num, p2, '  Line point P2:' )

  call plane_normal_line_exp_int_3d ( pp, normal, p1, p2, ival, pint )
 
  write ( *, '(a)' ) ' '

  if ( ival == 0 ) then
    write ( *, '(a)' ) '  The plane and line do not intersect.'
  else if ( ival == 1 ) then
    write ( *, '(a)' ) '  The plane and line intersect at '
    write ( *, '(2x,3g14.6)' ) pint(1:dim_num)
  else if ( ival == 2 ) then
    write ( *, '(a)' ) '  The plane and line are coincident.'
    write ( *, '(a)' ) '  One of the infinitely many points of intersection:'
    write ( *, '(2x,3g14.6)' ) pint(1:dim_num)
  else
    write ( *, '(a)' ) '  The plane and the line do not intersect.'
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Expected answer:'
  write ( *, '(a)' ) '    The plane and line intersect at '
  write ( *, '(a)' ) '    7, 2, 3.'
 
  return
end
subroutine test0616 ( )

!*****************************************************************************80
!
!! TEST0616 tests PLANE_NORMAL_QR_TO_XYZ and PLANE_NORMAL_XYZ_TO_QR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) dif
  integer ( kind = 4 ) j
  real ( kind = 8 ) normal(m)
  real ( kind = 8 ) pp(m)
  real ( kind = 8 ) pq(m)
  real ( kind = 8 ) pr(m)
  real ( kind = 8 ) qr1(m-1,n)
  real ( kind = 8 ) qr2(m-1,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) xyz(m,n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0616'
  write ( *, '(a)' ) '  For a normal plane, with point PP and NORMAL vector,'
  write ( *, '(a)' ) '  and in-plane basis vectors PQ and PR,'
  write ( *, '(a)' ) '  PLANE_NORMAL_QR_TO_XYZ converts QR to XYZ coordinates;'
  write ( *, '(a)' ) '  PLANE_NORMAL_XYZ_TO_QR converts XYZ to QR coordinates.'
!
!  Choose PP and NORMAL at random.
!
  call r8vec_uniform_01 ( m, seed, pp )

  call r8vec_uniform_01 ( m, seed, normal )
!
!  Compute in-plane basis vectors PQ and PR.
!
  call plane_normal_basis_3d ( pp, normal, pq, pr )
!
!  Choose random Q, R coordinates.
!
  call r8mat_uniform_01 ( m - 1, n, seed, qr1 )
  call r8mat_transpose_print ( m - 1, n, qr1, '  QR1' )
!
!  Convert to XYZ.
!
  call plane_normal_qr_to_xyz ( pp, normal, pq, pr, n, qr1, xyz )
  call r8mat_transpose_print ( m, n, xyz, '  XYZ' )
!
!  Convert XYZ to QR.
!
  call plane_normal_xyz_to_qr ( pp, normal, pq, pr, n, xyz, qr2 )
  call r8mat_transpose_print ( m - 1, n, qr2, '  QR2' )

  dif = 0.0D+00
  do j = 1, n
    dif = max ( dif, sqrt ( sum ( ( qr1(1:m-1,j) - qr2(1:m-1,j) )**2 ) ) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Maximum difference was ', dif

  return
end
subroutine test0617 ( )

!*****************************************************************************80
!
!! TEST0617 tests PLANE_NORMAL_TETRAHEDRON_INTERSECT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) int_num
  real ( kind = 8 ) normal(3)
  real ( kind = 8 ) pint(3,4)
  real ( kind = 8 ) pp(3)
  real ( kind = 8 ), dimension ( 3, 4 ) :: t = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00  &
    /), (/ 3, 4 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0617'
  write ( *, '(a)' ) '  PLANE_NORMAL_TETRAHEDRON_INTERSECT determines'
  write ( *, '(a)' ) '  the intersection of a plane and tetrahedron.'

  do k = 1, 2

    if ( k == 1 ) then
      normal(1:3) = (/ 0.0D+00, 0.0D+00, 1.0D+00 /)
    else
      normal(1:3) = (/ 1.0D+00, 1.0D+00, 0.0D+00 /) / sqrt ( 2.0D+00 )
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Plane normal vector number ', k
    write ( *, '(a)' ) ' '
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) normal(1:3)

    do i = 0, 6

      pp(1:3) = normal(1:3) * real ( i, kind = 8 ) / 5.0D+00

      call plane_normal_tetrahedron_intersect ( pp, normal, t, int_num, pint )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Point on plane:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) pp(1:3)
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Number of intersection points = ', int_num
      write ( *, '(a)' ) ' '
      do j = 1, int_num
        write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) j, pint(1:3,j)
      end do

    end do

  end do

  return
end
subroutine test062 ( )

!*****************************************************************************80
!
!! TEST062 tests PLANE_NORMAL_TRIANGLE_INT_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) int_num
  real ( kind = 8 ), dimension(dim_num) :: normal = (/ &
    1.0D+00, -2.0D+00, -3.0D+00 /)
  real ( kind = 8 ) pint(dim_num,3)
  real ( kind = 8 ), dimension(dim_num) :: pp = (/ &
    0.0D+00, 0.0D+00, 2.0D+00 /)
  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ), dimension(dim_num,3,test_num) :: t_test = reshape ( (/ &
         3.0D+00,  0.0D+00, -7.0D+00, &
        13.0D+00, -4.0D+00, -1.0D+00, &
         5.0D+00,  1.0D+00, -2.0D+00, &
         3.0D+00,  0.0D+00, -7.0D+00, &
        13.0D+00, -4.0D+00, -1.0D+00, &
         9.0D+00,  3.0D+00,  8.0D+00, &
        -6.0D+00,  0.0D+00,  0.0D+00, &
         0.0D+00,  3.0D+00,  0.0D+00, &
         0.0D+00,  0.0D+00,  2.0D+00, &
        -4.0D+00,  1.0D+00,  0.0D+00, &
         0.0D+00,  6.0D+00, -2.0D+00, &
         0.0D+00,  0.0D+00,  1.0D+00 /), (/ dim_num, 3, test_num /) )

  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST062'
  write ( *, '(a)' ) '  PLANE_NORMAL_TRIANGLE_INT_3D finds the'
  write ( *, '(a)' ) '  intersection points of a normal form plane'
  write ( *, '(a)' ) '  and a triangle.'
 
  call r8vec_print ( dim_num, pp, '  The point PP:' )
  call r8vec_print ( dim_num, normal, '  The normal vector N:' )

  do test = 1, test_num
 
    t(1:dim_num,1:3) = t_test(1:dim_num,1:3,test)
 
    call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

    call plane_normal_triangle_int_3d ( pp, normal, t, int_num, pint )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of intersection points is ', int_num
    write ( *, '(a)' ) ' '

    call r8mat_transpose_print ( dim_num, int_num, pint, &
      '  Intersection points:' )
 
  end do
 
  return
end
subroutine test063 ( )

!*****************************************************************************80
!
!! TEST063 tests PLANE_NORMAL2EXP_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ), dimension(dim_num) :: normal = (/ &
    -0.2672612D+00, -0.5345225D+00, -0.8017837D+00 /)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) p3(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: pp = (/ &
    -1.0D+00,   0.0D+00, -1.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST063'
  write ( *, '(a)' ) '  PLANE_NORMAL2EXP_3D puts a plane defined by '
  write ( *, '(a)' ) '  point, normal form into explicit form.'

  call r8vec_print ( dim_num, pp, '  The point PP:' )

  call r8vec_print ( dim_num, normal, '  Normal vector:' )
 
  call plane_normal2exp_3d ( pp, normal, p1, p2, p3 )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,3g14.6)' ) '  P1: ', p1(1:dim_num)
  write ( *, '(a,3g14.6)' ) '  P2: ', p2(1:dim_num)
  write ( *, '(a,3g14.6)' ) '  P3: ', p3(1:dim_num)

  return
end
subroutine test064 ( )
!
!*****************************************************************************80
!
!! TEST064 tests PLANE_NORMAL2IMP_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ), dimension(dim_num) :: normal
  real ( kind = 8 ), dimension(dim_num,test_num) :: normal_test = reshape ( (/ &
    -0.2672612D+00, -0.5345225D+00, -0.8017837D+00, &
    -0.2672612D+00, -0.5345225D+00, -0.8017837D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension(dim_num) :: pp
  real ( kind = 8 ), dimension(dim_num,test_num) :: pp_test = reshape ( (/ &
     -1.0D+00,   0.0D+00, -1.0D+00, &
    -16.0D+00, 2.0D+00, 4.0D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST064'
  write ( *, '(a)' ) '  PLANE_NORMAL2IMP_3D puts a plane defined by '
  write ( *, '(a)' ) '  point, normal form into implicit ABCD form.'
 
  do test = 1, test_num

    pp(1:dim_num) = pp_test(1:dim_num,test)
    normal(1:dim_num) = normal_test(1:dim_num,test)

    call r8vec_print ( dim_num, pp, '  The point PP:' )

    call r8vec_print ( dim_num, normal, '  Normal vector:' )

    call plane_normal2imp_3d ( pp, normal, a, b, c, d )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Output:'
    write ( *, '(a)' ) ' '
    write ( *, '(a,4g14.6)' ) '    (A,B,C,D)= ', a, b, c, d

  end do

  return
end
subroutine test065 ( )

!*****************************************************************************80
!
!! TEST065 tests POINTS_CENTROID_2D.
!
!  Diagram:
!
!    !....3&11....
!    !............
!    !............
!    X..9.........
!    !.....5......
!    !...........6
!    !.4.2...10...
!    !.....8...12.
!    V............
!    !..7.........
!    !......1.....
!    !............
!    !............
!    !----V----X--
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) centroid_index
  real ( kind = 8 ), dimension ( dim_num, n ) :: p = reshape ( (/ &
     7.0D+00,  3.0D+00, &
     4.0D+00,  7.0D+00, &
     5.0D+00, 13.0D+00, &
     2.0D+00,  7.0D+00, &
     6.0D+00,  9.0D+00, &
    12.0D+00,  8.0D+00, &
     3.0D+00,  4.0D+00, &
     6.0D+00,  6.0D+00, &
     3.0D+00, 10.0D+00, &
     8.0D+00,  7.0D+00, &
     5.0D+00, 13.0D+00, &
    10.0D+00,  6.0D+00 /), (/ dim_num, n /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST065'
  write ( *, '(a)' ) '  POINTS_CENTROID_2D computes the centroid of a'
  write ( *, '(a)' ) '  discrete set of points.'

  call r8mat_transpose_print ( dim_num, n, p, '  The points:' )

  call points_centroid_2d ( n, p, centroid_index )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The centroid is point #:', centroid_index

  return
end
subroutine test066 ( )

!*****************************************************************************80
!
!! TEST066 tests POINTS_COLIN_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) colin
  real ( kind = 8 ), dimension(dim_num) :: p1
  real ( kind = 8 ), dimension(dim_num,test_num) :: p1_test = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension(dim_num) :: p2
  real ( kind = 8 ), dimension(dim_num,test_num) :: p2_test = reshape ( (/ &
    10.0D+00, 10.0D+00, &
     0.0D+00, 1.0D+00, &
     1.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension(dim_num) :: p3
  real ( kind = 8 ), dimension(dim_num,test_num) :: p3_test = reshape ( (/ &
      5.0D+00, 4.99D+00, &
    100.0D+00, 0.0D+00, &
      0.5D+00, 0.86602539D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST066'
  write ( *, '(a)' ) '  POINTS_COLIN_2D estimates the colinearity'
  write ( *, '(a)' ) '  of three points.'

  do test = 1, test_num

    p1(1:dim_num) = p1_test(1:dim_num,test)
    p2(1:dim_num) = p2_test(1:dim_num,test)
    p3(1:dim_num) = p3_test(1:dim_num,test)

    write ( *, '(a)' ) ' '

    if ( test == 1 ) then
      write ( *, '(a)' ) '  Points almost on a line: Expect tiny COLIN.'
    else if ( test == 2 ) then
      write ( *, '(a)' ) '  Two points close, one far: Expect tiny COLIN.'
    else if ( test == 3 ) then
      write ( *, '(a)' ) '  Points on an equilateral triangle: Expect COLIN = 1.'
    end if

    call points_colin_2d ( p1, p2, p3, colin )

    write ( *, '(a)' ) ' '
    write ( *, '(a,3g14.6)' ) '  P1: ', p1(1:dim_num)
    write ( *, '(a,3g14.6)' ) '  P2: ', p2(1:dim_num)
    write ( *, '(a,3g14.6)' ) '  P3: ', p3(1:dim_num)
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Colinearity index = ', colin

  end do

  return
end
subroutine test067 ( )

!*****************************************************************************80
!
!! TEST067 tests POINTS_HULL_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: node_num = 7

  integer ( kind = 4 ) hull_num
  real ( kind = 8 ) node_xy(2,node_num)
  integer ( kind = 4 ) hull(node_num)
  real ( kind = 8 ) hull_xy(2,node_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST067'
  write ( *, '(a)' ) '  POINTS_HULL_2D computes the convex hull'
  write ( *, '(a)' ) '  of a set of N 2D points using an algorithm'
  write ( *, '(a)' ) '  that is order NlogH.'
  write ( *, '(a)' ) '  (H is the number of points on the convex hull.)'

  node_xy(1:2,1:node_num) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 2.0D+00, &
    2.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 2.0D+00, &
    1.0D+00, 3.0D+00, &
    2.0D+00, 2.0D+00 /), (/ 2, node_num /) )

  call r8mat_transpose_print ( 2, node_num, node_xy, &
    '  Coordinates of the points:' )

  call points_hull_2d ( node_num, node_xy, hull_num, hull )

  hull_xy(1:2,1:hull_num) = node_xy(1:2,hull(1:hull_num))

  call r8mat_transpose_print ( 2, hull_num, hull_xy, &
    '  Coordinates of the convex hull:' )

  return
end
subroutine test068 ( )

!*****************************************************************************80
!
!! TEST068 tests the SPHERE_DISTANCE routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: test_num = 6

  real      ( kind = 8 ) dist1
  real      ( kind = 8 ) dist2
  real      ( kind = 8 ) dist3
  character ( len = 18 ), dimension ( test_num ) :: name = (/ &
    'Atlanta, Georgia  ', &
    'North Pole        ', &
    'South Pole        ', &
    'Timbuktu          ', &
    'San Antonio, Texas', &
    'Savannah, Georgia ' /)
  integer   ( kind = 4 ), dimension ( test_num ) :: lat_d =  (/ 33, 90, -90, 16, 29, 32 /)
  integer   ( kind = 4 ), dimension ( test_num ) :: lat_m =  (/ 11,  0,   0, 49, 25,  5 /)
  integer   ( kind = 4 ), dimension ( test_num ) :: long_d = (/ 82,  0,   0,  3, 98, 81 /)
  integer   ( kind = 4 ), dimension ( test_num ) :: long_m = (/ 34,  0,   0,  0, 30,  6 /)
  real      ( kind = 8 ) lat1
  real      ( kind = 8 ) lat2
  real      ( kind = 8 ) long1
  real      ( kind = 8 ) long2
  real      ( kind = 8 ), parameter :: radius = 3957.0D+00
  integer   ( kind = 4 ) test1
  integer   ( kind = 4 ) test2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST068'
  write ( *, '(a)' ) '  SPHERE_DISTANCE1, SPHERE_DISTANCE2 and SPHERE_DISTANCE3'
  write ( *, '(a)' ) '  measure the distance between two points on a sphere.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  All tests uses RADIUS = ', radius
  write ( *, '(a)' ) '  which is the radius of the earth in miles.'
  write ( *, '(a)' ) ' '

  do test1 = 1, test_num-1

    call dms_to_radians ( lat_d(test1), lat_m(test1), 0, lat1 )
    call dms_to_radians ( long_d(test1), long_m(test1), 0, long1 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,a)' ) '  Distance from ', name(test1)

    do test2 = test1+1, test_num

      call dms_to_radians ( lat_d(test2), lat_m(test2), 0, lat2 )
      call dms_to_radians ( long_d(test2), long_m(test2), 0, long2 )

      call sphere_distance1 ( lat1, long1, lat2, long2, radius, dist1 )
      call sphere_distance2 ( lat1, long1, lat2, long2, radius, dist2 )
      call sphere_distance3 ( lat1, long1, lat2, long2, radius, dist3 )

      write ( *, '(a,a,3g14.6)' ) '             to ', &
        name(test2), dist1, dist2, dist3

    end do

  end do

  return
end
subroutine test0685 ( )

!*****************************************************************************80
!
!! TEST0685 tests POLAR_TO_XY and XY_TO_POLAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10
  real ( kind = 8 ) xy1(2)
  real ( kind = 8 ) xy2(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0685'
  write ( *, '(a)' ) '  POLAR_TO_XY converts (R,Theta) to (X,Y);'
  write ( *, '(a)' ) '  XY_TO_POLAR converts (X,Y) to (R,Theta).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         X           Y     ===>  R           T   =>      X           Y'
  write ( *, '(a)' ) ' '

  b = -1.0D+00
  c = +1.0D+00
  seed = 123456789

  do test = 1, test_num

    xy1(1) = r8_uniform ( b, c, seed )
    xy1(2) = r8_uniform ( b, c, seed )

    call xy_to_polar ( xy1, r, t )
    call polar_to_xy ( r, t, xy2 )

    write ( *, '(2x,6f12.5)' ) xy1(1:2), r, t, xy2(1:2)

  end do

  return
end
subroutine test0755 ( )

!*****************************************************************************80
!
!! TEST0755 tests POLYGON_*_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) result
  real ( kind = 8 ), dimension(dim_num,n) :: v = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00 /), (/ dim_num, n /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0755'
  write ( *, '(a)' ) '  For a polygon in 2D:'
  write ( *, '(a)' ) '  POLYGON_1_2D integrates 1'
  write ( *, '(a)' ) '  POLYGON_X_2D integrates X'
  write ( *, '(a)' ) '  POLYGON_Y_2D integrates Y'
  write ( *, '(a)' ) '  POLYGON_XX_2D integrates X*X'
  write ( *, '(a)' ) '  POLYGON_XY_2D integrates X*Y'
  write ( *, '(a)' ) '  POLYGON_YY_2D integrates Y*Y'

  call r8mat_transpose_print ( dim_num, n, v, '  The polygon vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  F(X,Y)    Integral'
  write ( *, '(a)' ) ' '

  call polygon_1_2d ( n, v, result )
  write ( *, '(2x,a4,4x,g14.6)' ) '  1', result

  call polygon_x_2d ( n, v, result )
  write ( *, '(2x,a4,4x,g14.6)' ) '  X', result

  call polygon_y_2d ( n, v, result )
  write ( *, '(2x,a4,4x,g14.6)' ) '  Y', result

  call polygon_xx_2d ( n, v, result )
  write ( *, '(2x,a4,4x,g14.6)' ) 'X*X', result

  call polygon_xy_2d ( n, v, result )
  write ( *, '(2x,a4,4x,g14.6)' ) 'X*Y', result

  call polygon_yy_2d ( n, v, result )
  write ( *, '(2x,a4,4x,g14.6)' ) 'Y*Y', result

  return
end
subroutine test0757 ( )

!*****************************************************************************80
!
!! TEST0757 tests POLYGON_ANGLES_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: n = 6

  real ( kind = 8 ) angle(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) radians_to_degrees
  real ( kind = 8 ), dimension (dim_num,n) :: v = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    2.0D+00, 1.0D+00, &
    3.0D+00, 0.0D+00, &
    3.0D+00, 2.0D+00, &
    1.0D+00, 2.0D+00 /), (/ dim_num, n /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0757'
  write ( *, '(a)' ) '  For a polygon in 2D:'
  write ( *, '(a)' ) '  POLYGON_ANGLES_2D computes the angles.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of polygonal vertices = ', n

  call r8mat_transpose_print ( dim_num, n, v, '  The polygon vertices:' )

  call polygon_angles_2d ( n, v, angle )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Polygonal angles in degrees:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,g14.6)' ) i, radians_to_degrees ( angle(i) )
  end do

  return
end
subroutine test076 ( )

!*****************************************************************************80
!
!! TEST076 tests POLYGON_AREA_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) area
  real ( kind = 8 ) area_exact
  real ( kind = 8 ), dimension ( test_num ) :: area_exact_test = (/ &
    2.0D+00, 6.0D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ), dimension ( test_num ) :: n_test = (/ 4, 8 /)
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: v

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST076'
  write ( *, '(a)' ) '  For a polygon in 2D:'
  write ( *, '(a)' ) '  POLYGON_AREA_2D computes the area.'

  do test = 1, test_num

    n = n_test(test)
    area_exact = area_exact_test(test)

    allocate ( v(1:dim_num,1:n) )

    if ( test == 1 ) then

      v(1:dim_num,1:n) = reshape ( (/ &
        1.0D+00, 0.0D+00, &
        2.0D+00, 1.0D+00, &
        1.0D+00, 2.0D+00, &
        0.0D+00, 1.0D+00 /), (/ dim_num, n /) )

    else if ( test == 2 ) then

      v(1:dim_num,1:n) = reshape ( (/ &
        0.0D+00, 0.0D+00, &
        3.0D+00, 0.0D+00, &
        3.0D+00, 3.0D+00, &
        2.0D+00, 3.0D+00, &
        2.0D+00, 1.0D+00, &
        1.0D+00, 1.0D+00, &
        1.0D+00, 2.0D+00, &
        0.0D+00, 2.0D+00 /), (/ dim_num, n /) )

    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of polygonal vertices = ', n

    call r8mat_transpose_print ( dim_num, n, v, '  The polygon vertices:' )

    call polygon_area_2d ( n, v, area )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Exact area is        ', area_exact
    write ( *, '(a,g14.6)' ) '  The computed area is ', area
 
    deallocate ( v )

  end do

  return
end
subroutine test0765 ( )

!*****************************************************************************80
!
!! TEST0765 tests POLYGON_AREA_2D_2;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) area
  real ( kind = 8 ) area_exact
  real ( kind = 8 ), dimension ( test_num ) :: area_exact_test = (/ &
    2.0D+00, 6.0D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ), dimension ( test_num ) :: n_test = (/ 4, 8 /)
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: v

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0765'
  write ( *, '(a)' ) '  For a polygon in 2D:'
  write ( *, '(a)' ) '  POLYGON_AREA_2D_2 computes the area.'

  do test = 1, test_num

    n = n_test(test)
    area_exact = area_exact_test(test)

    allocate ( v(1:dim_num,1:n) )

    if ( test == 1 ) then

      v(1:dim_num,1:n) = reshape ( (/ &
        1.0D+00, 0.0D+00, &
        2.0D+00, 1.0D+00, &
        1.0D+00, 2.0D+00, &
        0.0D+00, 1.0D+00 /), (/ dim_num, n /) )

    else if ( test == 2 ) then

      v(1:dim_num,1:n) = reshape ( (/ &
        0.0D+00, 0.0D+00, &
        3.0D+00, 0.0D+00, &
        3.0D+00, 3.0D+00, &
        2.0D+00, 3.0D+00, &
        2.0D+00, 1.0D+00, &
        1.0D+00, 1.0D+00, &
        1.0D+00, 2.0D+00, &
        0.0D+00, 2.0D+00 /), (/ dim_num, n /) )

    end if

    call r8mat_transpose_print ( dim_num, n, v, '  The polygon vertices:' )

    call polygon_area_2d_2 ( n, v, area )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Exact area is        ', area_exact
    write ( *, '(a,g14.6)' ) '  The computed area is ', area
 
    deallocate ( v )

  end do

  return
end
subroutine test078 ( )

!*****************************************************************************80
!
!! TEST078 tests POLYGON_AREA_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) area
  real ( kind = 8 ) :: area_exact
  real ( kind = 8 ), dimension ( test_num ) :: area_exact_test = (/ &
    2.4494898D+00, 6.0D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ), dimension ( test_num ) :: n_test = (/ 4, 8 /)
  real ( kind = 8 ), dimension ( dim_num ) ::  normal
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension (:,:) :: v

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST078'
  write ( *, '(a)' ) '  For a polygon in 3D:'
  write ( *, '(a)' ) '  POLYGON_AREA_3D computes the area;'

  do test = 1, test_num

    area_exact = area_exact_test(test)
    n = n_test(test)

    allocate ( v(1:dim_num,1:n) )

    if ( test == 1 ) then

      v = reshape ( (/ &
        1.0D+00, 0.0D+00, 0.0D+00, &
        2.0D+00, 1.0D+00, 1.0D+00, &
        1.0D+00, 2.0D+00, 1.0D+00, &
        0.0D+00, 1.0D+00, 0.0D+00 /), (/ dim_num, n /) )

    else if ( test == 2 ) then

      v = reshape ( (/ &
        0.00000D+00,  0.00000D+00,  0.00000D+00, &    
        2.62679D+00,  1.26009D+00, -0.715657D+00, &    
        1.48153D+00,  3.97300D+00, -0.142512D+00, &    
        0.605932D+00, 3.55297D+00,  0.960401D-01, &
        1.36944D+00,  1.74437D+00, -0.286056D+00, &
        0.493842D+00, 1.32433D+00, -0.475041D-01, &
        0.112090D+00, 2.22864D+00,  0.143544D+00, &    
       -0.763505D+00, 1.80861D+00,  0.382097D+00 /), (/ dim_num, n /) )

    end if

    call r8mat_transpose_print ( dim_num, n, v, '  The polygon vertices:' )

    call polygon_area_3d ( n, v, area, normal )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Exact area is        ', area_exact
    write ( *, '(a,g14.6)' ) '  The computed area is ', area

    deallocate ( v )
 
  end do

  return
end
subroutine test0782 ( )

!*****************************************************************************80
!
!! TEST0782 tests POLYGON_AREA_3D_2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) area
  real ( kind = 8 ) :: area_exact
  real ( kind = 8 ), dimension ( test_num ) :: area_exact_test = (/ &
    2.4494898D+00, 6.0D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ), dimension ( test_num ) :: n_test = (/ 4, 8 /)
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension (:,:) :: v

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0782'
  write ( *, '(a)' ) '  For a polygon in 3D:'
  write ( *, '(a)' ) '  POLYGON_AREA_3D_2 computes the area;'

  do test = 1, test_num

    area_exact = area_exact_test(test)
    n = n_test(test)

    allocate ( v(1:dim_num,1:n) )

    if ( test == 1 ) then

      v = reshape ( (/ &
        1.0D+00, 0.0D+00, 0.0D+00, &
        2.0D+00, 1.0D+00, 1.0D+00, &
        1.0D+00, 2.0D+00, 1.0D+00, &
        0.0D+00, 1.0D+00, 0.0D+00 /), (/ dim_num, n /) )

    else if ( test == 2 ) then

      v = reshape ( (/ &
        0.00000D+00,  0.00000D+00,  0.00000D+00, &    
        2.62679D+00,  1.26009D+00, -0.715657D+00, &    
        1.48153D+00,  3.97300D+00, -0.142512D+00, &    
        0.605932D+00, 3.55297D+00,  0.960401D-01, &
        1.36944D+00,  1.74437D+00, -0.286056D+00, &
        0.493842D+00, 1.32433D+00, -0.475041D-01, &
        0.112090D+00, 2.22864D+00,  0.143544D+00, &    
       -0.763505D+00, 1.80861D+00,  0.382097D+00 /), (/ dim_num, n /) )

    end if

    call r8mat_transpose_print ( dim_num, n, v, '  The polygon vertices:' )

    call polygon_area_3d_2 ( n, v, area )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Exact area is        ', area_exact
    write ( *, '(a,g14.6)' ) '  The computed area is ', area

    deallocate ( v )
 
  end do

  return
end
subroutine test0784 ( )

!*****************************************************************************80
!
!! TEST0784 tests POLYGON_CENTROID_2D and POLYGON_CENTROID_2D_2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) centroid(dim_num)
  real ( kind = 8 ), dimension (dim_num,n) :: v = reshape ( (/ &
    1.0D+00, 0.0D+00, &
    2.0D+00, 1.0D+00, &
    1.0D+00, 2.0D+00, &
    0.0D+00, 1.0D+00 /), (/ dim_num, n /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0784'
  write ( *, '(a)' ) '  For a polygon in 2D:'
  write ( *, '(a)' ) '  POLYGON_CENTROID_2D computes the centroid.'
  write ( *, '(a)' ) '  POLYGON_CENTROID_2D_2 computes the centroid.'

  call r8mat_transpose_print ( dim_num, n, v, '  The polygon vertices:' )

  call polygon_centroid_2d ( n, v, centroid )

  call r8vec_print ( dim_num, centroid, '  POLYGON_CENTROID_2D:' )

  call polygon_centroid_2d_2 ( n, v, centroid )

  call r8vec_print ( dim_num, centroid, '  POLYGON_CENTROID_2D_2:' )

  return
end
subroutine test0786 ( )

!*****************************************************************************80
!
!! TEST0786 tests POLYGON_CENTROID_3D.
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

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) centroid(dim_num)
  real ( kind = 8 ), dimension (dim_num,n) :: v = reshape ( (/ &
    1.0D+00, 0.0D+00, 0.0D+00, &
    2.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 2.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00 /), (/ dim_num, n /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0786'
  write ( *, '(a)' ) '  For a polygon in 3D:'
  write ( *, '(a)' ) '  POLYGON_CENTROID_3D computes the centroid.'

  call r8mat_transpose_print ( dim_num, n, v, '  The polygon vertices:' )

  call polygon_centroid_3d ( n, v, centroid )

  call r8vec_print ( dim_num, centroid, '  The centroid:' )
 
  return
end
subroutine test079 ( )

!*****************************************************************************80
!
!! TEST079 tests POLYGON_CONTAINS_POINT_2D and POLYGON_CONTAINS_POINT_2D_2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) i
  logical inside1
  logical inside2
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension (dim_num,test_num) :: p_test = reshape ( (/ &
    1.0D+00,  1.0D+00, &
    3.0D+00,  4.0D+00, &
    0.0D+00,  2.0D+00, &
    0.5D+00, -0.25D+00 /), (/dim_num, test_num /) )
  integer ( kind = 4 ) test
  real ( kind = 8 ), dimension(dim_num,n) :: v = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    2.0D+00, 1.0D+00, &
    1.0D+00, 2.0D+00, &
    0.0D+00, 2.0D+00 /), (/ dim_num, n /) )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST079'
  write ( *, '(a)' ) '  POLYGON_CONTAINS_POINT_2D determines if '
  write ( *, '(a)' ) '  a point is in a polygon.'
  write ( *, '(a)' ) '  POLYGON_CONTAINS_POINT_2D_2 determines if'
  write ( *, '(a)' ) '  a point is in a polygon.'

  call r8mat_transpose_print ( dim_num, n, v, '  The polygon vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          P          In1  In2'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
 
    p(1:dim_num) = p_test(1:dim_num,test)
 
    call polygon_contains_point_2d ( n, v, p, inside1 )

    call polygon_contains_point_2d_2 ( n, v, p, inside2 )

    write ( *, '(2x,2g14.6,4x,l1,4x,l1)' ) p(1:dim_num), inside1, inside2

  end do
 
  return
end
subroutine test080 ( )

!*****************************************************************************80
!
!! TEST080 tests POLYGON_DIAMETER_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) diameter
  real ( kind = 8 ) :: diameter_exact = 2.0D+00
  real ( kind = 8 ), dimension ( dim_num, n ) :: v = reshape ( (/ &
    1.0D+00, 0.0D+00, &
    2.0D+00, 1.0D+00, &
    1.0D+00, 2.0D+00, &
    0.0D+00, 1.0D+00 /), (/ dim_num, n /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST080'
  write ( *, '(a)' ) '  For a polygon in 2D:'
  write ( *, '(a)' ) '  POLYGON_DIAMETER_2D computes the diameter;'

  call r8mat_transpose_print ( dim_num, n, v, '  The polygon vertices:' )

  call polygon_diameter_2d ( n, v, diameter )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Diameter ( computed ) ', diameter
  write ( *, '(a,g14.6)' ) '  Diameter ( exact )    ', diameter_exact
 
  return
end
subroutine test0801 ( )

!*****************************************************************************80
!
!! TEST0801 tests POLYGON_EXPAND_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) h
  real ( kind = 8 ), dimension(dim_num,n) :: v = reshape ( (/ &
    1.0D+00, 1.0D+00, &
    5.0D+00, 1.0D+00, &
    2.0D+00, 4.0D+00, &
    1.0D+00, 3.0D+00 /), (/ dim_num, n /) )
  real ( kind = 8 ), dimension(dim_num,n) :: w

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0801'
  write ( *, '(a)' ) '  For a polygon in 2D:'
  write ( *, '(a)' ) '  POLYGON_EXPAND_2D "expands" it by an amount H.'

  h = 0.5D+00

  call r8mat_transpose_print ( dim_num, n, v, '  The polygon vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The expansion amount H = ', h

  call polygon_expand_2d ( n, v, h, w )

  call r8mat_transpose_print ( dim_num, n, w, '  The expanded polygon:' )

  return
end
subroutine test0803 ( )

!*****************************************************************************80
!
!! TEST0803 tests POLYGON_INRAD_DATA_2D, POLYGON_OUTRAD_DATA_2D, POLYGON_SIDE_DATA_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) area
  integer ( kind = 4 ) n
  real ( kind = 8 ) radin
  real ( kind = 8 ) radout
  real ( kind = 8 ) side

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0803'
  write ( *, '(a)' ) '  For a REGULAR polygon in 2D:'
  write ( *, '(a)' ) '  the inradius, outradius and side are related.'
  write ( *, '(a)' ) '  POLYGON_INRAD_DATA_2D uses the inradius;'
  write ( *, '(a)' ) '  POLYGON_OUTRAD_DATA_2D uses the inradius;'
  write ( *, '(a)' ) '  POLYGON_SIDE_DATA_2D uses the inradius;'

  do n = 3, 5

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of polygonal sides = ', n

    side = 1.0D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Assuming SIDE = ', side
    call polygon_side_data_2d ( n, side, area, radin, radout )
    write ( *, '(a,g14.6)' ) '    AREA =   ', area
    write ( *, '(a,g14.6)' ) '    RADIN =  ', radin
    write ( *, '(a,g14.6)' ) '    RADOUT = ', radout
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Assuming RADIN = ', radin
    call polygon_inrad_data_2d ( n, radin, area, radout, side )
    write ( *, '(a,g14.6)' ) '    AREA =   ', area
    write ( *, '(a,g14.6)' ) '    RADOUT = ', radout
    write ( *, '(a,g14.6)' ) '    SIDE =   ', side
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Assuming RADOUT = ', radout
    call polygon_outrad_data_2d ( n, radout, area, radin, side )
    write ( *, '(a,g14.6)' ) '    AREA =   ', area
    write ( *, '(a,g14.6)' ) '    RADIN =  ', radin
    write ( *, '(a,g14.6)' ) '    SIDE =   ', side

  end do

  return
end
subroutine test0805 ( )

!*****************************************************************************80
!
!! TEST0805 tests POLYGON_IS_CONVEX_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 10
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  character ( len = 80 ) message(-1:2)
  integer ( kind = 4 ) n
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) polygon_is_convex_2d
  integer ( kind = 4 ) result
  real ( kind = 8 ) v(dim_num,n_max)

  message(-1) = 'The polygon is not convex.'
  message( 0) = 'The polygon is degenerate and convex.'
  message( 1) = 'The polygon is convex and counterclockwise.'
  message( 2) = 'The polygon is convex and clockwise.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0805'
  write ( *, '(a)' ) '  POLYGON_IS_CONVEX_2D determines if a polygon'
  write ( *, '(a)' ) '  is convex.'
!
!  Shape 1: a point
!
  n = 1

  v(1:dim_num,1:n) = reshape ( (/ &
    0.0D+00, 0.0D+00 /), (/ dim_num, n /) )

  call r8mat_transpose_print ( dim_num, n, v, '  Shape #1, a point:' )

  result = polygon_is_convex_2d ( n, v )

  write ( *, '(2x,a)' ) message(result)
!
!  Shape 2: a line
!
  n = 2

  v(1:dim_num,1:n) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 2.0D+00 /), (/ dim_num, n /) )

  call r8mat_transpose_print ( dim_num, n, v, '  Shape #2, a line:' )

  result = polygon_is_convex_2d ( n, v )

  write ( *, '(2x,a)' ) message(result)
!
!  Shape 3: a flat "triangle."
!
  n = 3
  v(1:dim_num,1:n) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    2.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00 /), (/ dim_num, n /) )

  call r8mat_transpose_print ( dim_num, n, v, '  Shape #3, a flat triangle:' )

  result = polygon_is_convex_2d ( n, v )

  write ( *, '(2x,a)' ) message(result)
!
!  Shape 4: a triangle.
!
  n = 3
  v(1:dim_num,1:n) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    0.0D+00, 2.0D+00 /), (/ dim_num, n /) )

  call r8mat_transpose_print ( dim_num, n, v, '  Shape #4, a CCW triangle:' )

  result = polygon_is_convex_2d ( n, v )

  write ( *, '(2x,a)' ) message(result)
!
!  Shape 5: CW triangle.
!
  n = 3
  v(1:dim_num,1:n) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    0.0D+00, 2.0D+00, &
    1.0D+00, 0.0D+00 /), (/ dim_num, n /) )

  call r8mat_transpose_print ( dim_num, n, v, '  Shape #5, a CW triangle:' )

  result = polygon_is_convex_2d ( n, v )

  write ( *, '(2x,a)' ) message(result)
!
!  Shape 6: polygon with an interior angle of more than 90.
!
  n = 4
  v(1:dim_num,1:n) = reshape ( (/ &
    1.0D+00, 0.0D+00, &
    2.0D+00, 0.0D+00, &
    3.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00 /), (/ dim_num, n /) )

  call r8mat_transpose_print ( dim_num, n, v, '  #6, Polygon with large angle:' )

  result = polygon_is_convex_2d ( n, v )

  write ( *, '(2x,a)' ) message(result)
!
!  Shape 7: polygon with an interior angle of more than 180.
!
  n = 5
  v(1:dim_num,1:n) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    0.5D+00, 0.5D+00, &
    1.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00 /), (/ dim_num, n /) )

  call r8mat_transpose_print ( dim_num, n, v, '  #7, Polygon with huge angle:' )

  result = polygon_is_convex_2d ( n, v )

  write ( *, '(2x,a)' ) message(result)
!
!  Shape 8: star
!
  n = 5

  do i = 1, n
    angle = real ( i - 1, kind = 8 ) * 4.0D+00 * pi / real ( n, kind = 8 )
    v(1,i) = cos ( angle )
    v(2,i) = sin ( angle )
  end do

  call r8mat_transpose_print ( dim_num, n, v, '  #8, a star:' )

  result = polygon_is_convex_2d ( n, v )

  write ( *, '(2x,a)' ) message(result)
!
!  Shape 9: regular hexagon
!
  n = 6

  do i = 1, n
    angle = real ( i - 1, kind = 8 ) * 2.0D+00 * pi / real ( n, kind = 8 )
    v(1,i) = cos ( angle )
    v(2,i) = sin ( angle )
  end do

  call r8mat_transpose_print ( dim_num, n, v, '  #9, regular hexagon:' )

  result = polygon_is_convex_2d ( n, v )

  write ( *, '(2x,a)' ) message(result)
!
!  Shape 10: double triangle
!
  n = 6
  v(1:dim_num,1:n) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    2.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    2.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00 /), (/ dim_num, n /) )

  call r8mat_transpose_print ( dim_num, n, v, '  #10, double hexagon:' )

  result = polygon_is_convex_2d ( n, v )

  write ( *, '(2x,a)' ) message(result)
!
!  Shape 11: "square knot"
!
  n = 8
  v(1:dim_num,1:n) = reshape ( (/ &
    1.0D+00, 0.0D+00, &
    3.0D+00, 0.0D+00, &
    3.0D+00, 3.0D+00, &
    0.0D+00, 3.0D+00, &
    0.0D+00, 1.0D+00, &
    2.0D+00, 1.0D+00, &
    2.0D+00, 2.0D+00, &
    1.0D+00, 2.0D+00 /), (/ dim_num, n /) )

  call r8mat_transpose_print ( dim_num, n, v, '  #11, square knot:' )

  result = polygon_is_convex_2d ( n, v )

  write ( *, '(2x,a)' ) message(result)

  return
end
subroutine test0807 ( )

!*****************************************************************************80
!
!! TEST0807 tests POLYGON_SOLID_ANGLE_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) n
  real ( kind = 8 ), dimension ( dim_num ) :: p
  real ( kind = 8 ) solid_angle
  integer ( kind = 4 ) test
  real ( kind = 8 ), allocatable, dimension (:,:) :: v

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0807'
  write ( *, '(a)' ) '  For a polygon in 3D:'
  write ( *, '(a)' ) '  POLYGON_SOLID_ANGLE_3D computes the solid angle'
  write ( *, '(a)' ) '  subtended by a planar polygon as viewed from'
  write ( *, '(a)' ) '  a point P.'

  do test = 1, test_num
!
!  One eighth of sphere surface, on the unit sphere surface.
!
    if ( test == 1 ) then

      n = 3

      allocate ( v(1:dim_num,1:n) )

      v = reshape ( (/ &
        1.0D+00, 0.0D+00, 0.0D+00, &
        0.0D+00, 1.0D+00, 0.0D+00, &
        0.0D+00, 0.0D+00, 1.0D+00 /), (/ dim_num, n /) )

      p(1:3) = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
!
!  Reverse order of vertices.
!
    else if ( test == 2 ) then

      n = 3

      allocate ( v(1:dim_num,1:n) )

      v = reshape ( (/ &
        1.0D+00, 0.0D+00, 0.0D+00, &
        0.0D+00, 0.0D+00, 1.0D+00, &
        0.0D+00, 1.0D+00, 0.0D+00 /), (/ dim_num, n /) )

      p(1:3) = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
!
!  One eighth of sphere surface, on the unit sphere surface, 
!  translated by (1,2,3).
!
    else if ( test == 3 ) then

      n = 3

      allocate ( v(1:dim_num,1:n) )

      v = reshape ( (/ &
        2.0D+00, 2.0D+00, 3.0D+00, &
        1.0D+00, 3.0D+00, 3.0D+00, &
        1.0D+00, 2.0D+00, 4.0D+00 /), (/ dim_num, n /) )

      p(1:3) = (/ 1.0D+00, 2.0D+00, 3.0D+00 /)
!
!  One eighth of sphere surface, but on sphere of radius 2.
!
    else if ( test == 4 ) then

      n = 3

      allocate ( v(1:dim_num,1:n) )

      v = reshape ( (/ &
        2.0D+00, 0.0D+00, 0.0D+00, &
        0.0D+00, 2.0D+00, 0.0D+00, &
        0.0D+00, 0.0D+00, 2.0D+00 /), (/ dim_num, n /) )

      p(1:3) = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)

    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  TEST # ', test
    write ( *, '(a)' ) ' '

    call r8vec_print ( dim_num, p, '  The viewing point P:' )

    call r8mat_transpose_print ( dim_num, n, v, '  The polygon vertices V:' )

    call polygon_solid_angle_3d ( n, v, p, solid_angle )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Solid angle subtended: ', solid_angle

    deallocate ( v )
 
  end do

  return
end
subroutine test081 ( )

!*****************************************************************************80
!
!! TEST081 tests POLYHEDRON_AREA_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 3
  integer ( kind = 4 ), parameter :: face_num = 4
  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: node_num = 4

  real ( kind = 8 ) area
  real ( kind = 8 ), parameter :: area_exact = 2.366025D+00
  real ( kind = 8 ), dimension ( dim_num, node_num ) :: coord = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00 /), (/ dim_num, node_num /) )
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), dimension(face_num,order_max) :: node = reshape ( (/ &
    3, 1, 1, 2, &
    2, 2, 4, 3, &
    1, 4, 3, 4 /), (/ face_num,order_max /) )
  integer ( kind = 4 ), dimension ( face_num ) :: order = (/ 3, 3, 3, 3 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST081'
  write ( *, '(a)' ) '  For a polyhedron in 3D:'
  write ( *, '(a)' ) '  POLYHEDRON_AREA_3D computes surface area;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of faces is ', face_num

  call i4vec_print ( face_num, order, '  Order of each face:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Nodes per face:'
  write ( *, '(a)' ) ' '
  do i = 1, face_num
    write ( *, '(5i8)' ) i, ( node(i,j), j = 1, order(i) )
  end do

  call r8mat_transpose_print ( dim_num, node_num, coord, '  Polyhedron nodes' )

  call polyhedron_area_3d ( coord, order_max, face_num, node, node_num, &
    order, area )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Surface area = ', area
  write ( *, '(a,g14.6)' ) '  Exact area =   ', area_exact

  return
end
subroutine test082 ( )

!*****************************************************************************80
!
!! TEST082 tests POLYHEDRON_CENTROID_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 3
  integer ( kind = 4 ), parameter :: face_num = 4
  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: node_num = 4

  real ( kind = 8 ) centroid(dim_num)
  real ( kind = 8 ), dimension ( dim_num ) :: centroid_exact = (/ &
    0.25D+00, 0.25D+00, 0.25D+00 /)
  real ( kind = 8 ), dimension ( dim_num, node_num ) :: coord = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00 /), (/ dim_num, node_num /) )
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), dimension(face_num,order_max) :: node = reshape ( (/ &
    3, 1, 1, 2, &
    2, 2, 4, 3, &
    1, 4, 3, 4 /), (/ face_num,order_max /) )
  integer ( kind = 4 ), dimension ( face_num ) :: order = (/ 3, 3, 3, 3 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST082'
  write ( *, '(a)' ) '  For a polyhedron in 3D:'
  write ( *, '(a)' ) '  POLYHEDRON_CENTROID_3D computes the centroid;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of faces is ', face_num

  call i4vec_print ( face_num, order, '  Order of each face:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Nodes per face:'
  write ( *, '(a)' ) ' '
  do i = 1, face_num
    write ( *, '(5i8)' ) i, ( node(i,j), j = 1, order(i) )
  end do

  call r8mat_transpose_print ( dim_num, node_num, coord, '  Polyhedron nodes:' )

  call polyhedron_centroid_3d ( coord, order_max, face_num, node, node_num, &
    order, centroid )

  call r8vec_print ( dim_num, centroid, '  Computed centroid:' )
  call r8vec_print ( dim_num, centroid_exact, '  Exact centroid:' )

  return
end
subroutine test0825 ( )

!*****************************************************************************80
!
!! TEST0825 tests POLYHEDRON_CONTAINS_POINT_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: face_num = 4
  integer ( kind = 4 ), parameter :: face_order_max = 3
  integer ( kind = 4 ), parameter :: node_num = 4
  integer ( kind = 4 ), parameter :: test_num = 100

  real ( kind = 8 ) area
  real ( kind = 8 ) c(dim_num+1)
  integer ( kind = 4 ), dimension(face_num) :: face_order = (/ 3, 3, 3, 3 /)
  integer ( kind = 4 ), dimension (face_order_max,face_num) :: face_point = reshape ( (/ &
    1, 2, 4, &
    1, 3, 2, &
    1, 4, 3, &
    2, 3, 4 /), (/ face_order_max, face_num /) )
  logical inside1
  logical inside2
  real ( kind = 8 ) p(dim_num)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  real ( kind = 8 ), dimension ( dim_num, node_num ) :: v = reshape ( (/ & 
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00 /), (/ dim_num, node_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0825'
  write ( *, '(a)' ) '  POLYHEDRON_CONTAINS_POINT_3D determines if a point'
  write ( *, '(a)' ) '  is inside a polyhedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We test this routine by using a tetrahedron as '
  write ( *, '(a)' ) '  the polyhedron.'
  write ( *, '(a)' ) '  For this shape, an independent check can be made,'
  write ( *, '(a)' ) '  using barycentric coordinates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We label these checks IN1 and IN2, and '
  write ( *, '(a)' ) '  we expect them to agree.'

  call r8mat_transpose_print ( dim_num, node_num, v, '  The vertices:' )

  call i4vec_print ( face_num, face_order, '  The face orders:' )

  call i4mat_transpose_print ( face_order_max, face_num, face_point, &
    '  The nodes making each face:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X           Y           Z      IN1 IN2'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do test = 1, test_num
  
    call r8vec_uniform_01 ( dim_num, seed, p )

    call polyhedron_contains_point_3d ( node_num, face_num, &
      face_order_max, v, face_order, face_point, p, inside1 )

    call tetrahedron_barycentric_3d ( v, p, c )

    inside2 =  ( 0.0D+00 <= c(1) ) .and. ( c(1) <= 1.0D+00 ) .and. &
               ( 0.0D+00 <= c(2) ) .and. ( c(2) <= 1.0D+00 ) .and. &
               ( 0.0D+00 <= c(3) ) .and. ( c(3) <= 1.0D+00 ) .and. &
               ( 0.0D+00 <= c(4) ) .and. ( c(4) <= 1.0D+00 ) .and. &
               ( c(1) + c(2) + c(3) + c(4) <= 1.0D+00 )

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,l1,2x,l1)' ) &
      p(1:3), inside1, inside2

    if ( inside1 .neqv. inside2 ) then
      write ( *, '(a)' ) '??? Disagreement!  Barycentric coordinates:'
      write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) c(1:4)
    end if

  end do

  return
end
subroutine test083 ( )

!*****************************************************************************80
!
!! TEST083 tests POLYHEDRON_VOLUME_3D and POLYHEDRON_VOLUME_3D_2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order_max = 3
  integer ( kind = 4 ), parameter :: face_num = 4
  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: node_num = 4

  real ( kind = 8 ), dimension ( dim_num, node_num ) :: coord = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00 /), (/ dim_num, node_num /) )
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), dimension(face_num,order_max) :: node = reshape ( (/ &
    3, 1, 1, 2, &
    2, 2, 4, 3, &
    1, 4, 3, 4 /), (/ face_num,order_max /) )
  integer ( kind = 4 ), dimension ( face_num ) :: order = (/ 3, 3, 3, 3 /)
  real ( kind = 8 ) :: volume_exact = 1.0D+00 / 6.0D+00
  real ( kind = 8 ) volume1
  real ( kind = 8 ) volume2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST083'
  write ( *, '(a)' ) '  For a polyhedron in 3D:'
  write ( *, '(a)' ) '  POLYHEDRON_VOLUME_3D computes volume.'
  write ( *, '(a)' ) '  POLYHEDRON_VOLUME_3D_2 computes volume.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of faces is ', face_num

  call i4vec_print ( face_num, order, '  Order of each face:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Nodes per face:'
  write ( *, '(a)' ) ' '
  do i = 1, face_num
    write ( *, '(5i8)' ) i, ( node(i,j), j = 1, order(i) )
  end do

  call r8mat_transpose_print ( dim_num, node_num, coord, '  Polyhedron nodes' )

  call polyhedron_volume_3d ( coord, order_max, face_num, node, node_num, &
    order, volume1 )

  call polyhedron_volume_3d_2 ( coord, order_max, face_num, node, node_num, &
    order, volume2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Volume ( method 1 ) = ', volume1
  write ( *, '(a,g14.6)' ) '  Volume ( method 2 ) = ', volume2
  write ( *, '(a,g14.6)' ) '  Volume ( exact ) =    ', volume_exact

  return
end
subroutine test084 ( )

!*****************************************************************************80
!
!! TEST084 tests POLYLINE_ARCLENGTH_ND and POLYLINE_INDEX_POINT_ND;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) i
  real ( kind = 8 ), dimension(dim_num,n) :: p = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    2.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00 /), (/ dim_num, n /) )
  real ( kind = 8 ) pt(dim_num)
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) t

  t = 2.0D+00
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST084'
  write ( *, '(a)' ) '  POLYLINE_INDEX_POINT_ND finds a point on a '
  write ( *, '(a)' ) '  polyline with given arclength.'
  write ( *, '(a)' ) '  POLYLINE_ARCLENGTH_ND computes the arclength '
  write ( *, '(a)' ) '  of the polyline, and its nodes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The line we examine is defined by these points:'
!
!  The call to POLYLINE_ARCLENGTH_ND is just to help us believe 
!  the final result.
!
  call polyline_arclength_nd ( dim_num, n, p, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           P              Arclength(X,Y)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,3g14.6)' ) p(1:dim_num,i), s(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  We search for the point with coordinate ', t
 
  call polyline_index_point_nd ( dim_num, n, p, t, pt )

  call r8vec_print ( dim_num, pt, '  The computed point:' ) 

  return
end
subroutine test0844 ( )

!*****************************************************************************80
!
!! TEST0844 tests POLYLINE_POINTS_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: nk = 4
  integer ( kind = 4 ), parameter :: nt = 13

  real ( kind = 8 ), dimension(dim_num,nk) :: pk = reshape ( (/ &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    1.0D+00, 2.0D+00 /), (/ dim_num, nk /) )
  real ( kind = 8 ) pt(dim_num,nt)
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0844'
  write ( *, '(a)' ) '  POLYLINE_POINTS_ND computes points on a polyline.'

  call r8mat_transpose_print ( dim_num, nk, pk, '  The defining points:' )

  call polyline_points_nd ( dim_num, nk, pk, nt, pt )
 
  call r8mat_transpose_print ( dim_num, nt, pt, '  The computed points:' )

  return
end
subroutine test0845 ( )

!*****************************************************************************80
!
!! TEST0845 tests POLYLOOP_ARCLENGTH_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  real ( kind = 8 ), dimension(dim_num,n) :: p = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    2.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00 /), (/ dim_num, n /) )
  real ( kind = 8 ) s(n+1)
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0845'
  write ( *, '(a)' ) '  POLYLOOP_ARCLENGTH_ND computes the arclength '
  write ( *, '(a)' ) '  of the nodes of a polyloop.'

  call polyloop_arclength_nd ( dim_num, n, p, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           P            Arclength(P)'
  write ( *, '(a)' ) ' '

  do j = 1, n + 1
    j2 = i4_wrap ( j, 1, n )
    write ( *, '(2x,3g14.6)' ) p(1:dim_num,j2), s(j)
  end do

  return
end
subroutine test0846 ( )

!*****************************************************************************80
!
!! TEST0846 tests POLYLOOP_POINTS_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: nk = 4
  integer ( kind = 4 ), parameter :: nt = 12

  real ( kind = 8 ), dimension(dim_num,nk) :: pk = reshape ( (/ &
    0.0D+00, 2.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    1.0D+00, 2.0D+00 /), (/ dim_num, nk /) )
  real ( kind = 8 ) pt(dim_num,nt)
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0846'
  write ( *, '(a)' ) '  POLYLOOP_POINTS_ND computes points on a polyloop.'

  call r8mat_transpose_print ( dim_num, nk, pk, '  The defining points:' )

  call polyloop_points_nd ( dim_num, nk, pk, nt, pt )
 
  call r8mat_transpose_print ( dim_num, nt, pt, '  The computed points:' )

  return
end
subroutine test085 ( )

!*****************************************************************************80
!
!! TEST085 tests PLANE_EXP_PRO3.
!
!  Discussion:
!
!    Projection is ( -1, 1, 1 ).
!    Projection is ( 4, 5, -8 ).
!    Projection is ( 0.33, 0.33, 0.33).
!    Projection is ( 5.33, -1.66, -2.66 ).
!    Projection is ( -1, 1, 1 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 5

  integer ( kind = 4 ) i
  real ( kind = 8 ), dimension ( dim_num ) :: p1 = (/ &
    1.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p2 = (/ &
    0.0D+00, 1.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p3 = (/ &
    0.0D+00, 0.0D+00, 1.0D+00 /)
  real ( kind = 8 ), dimension (dim_num,test_num) :: po = reshape ( (/ &
    0.0D+00,  2.0D+00,  2.0D+00, &
    4.0D+00,  5.0D+00, -8.0D+00, &
    0.25D+00, 0.25D+00, 0.25D+00, &
    5.0D+00, -2.0D+00, -3.0D+00, &
   -2.0D+00,  0.0D+00,  0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) pp(dim_num,test_num)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST085'
  write ( *, '(a)' ) '  PLANE_EXP_PRO3 projects an object point '
  write ( *, '(a)' ) '  orthographically into a plane.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '            PO                PP'
  write ( *, '(a)' ) ' '

  call plane_exp_pro3 ( p1, p2, p3, test_num, po, pp )
 
  do test = 1, test_num
    write ( *, '(2x,6g12.4)' ) po(1:dim_num,test), pp(1:dim_num,test)
  end do
 
  return
end
subroutine test170 ( )

!*****************************************************************************80
!
!! TEST170 tests PROVEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 2

  real ( kind = 8 ), dimension(m,n) :: base = reshape ( (/ &
    4.0D+00, 3.0D+00, 2.0D+00, 1.0D+00, &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00 /), (/ m, n /) )
  real ( kind = 8 ), dimension ( m ) :: vecm = (/ &
    1.0D+00, 1.0D+00, 1.0D+00, 2.0D+00 /)
  real ( kind = 8 ) vecn(n)
  real ( kind = 8 ) vecnm(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST170'
  write ( *, '(a)' ) '  PROVEC projects a vector onto a subspace.'
 
  call r8mat_transpose_print ( m, n, base, '  Base vectors' )
  
  call r8vec_print ( m, vecm, '  Vector to be projected:' )
 
  call provec ( m, n, base, vecm, vecn, vecnm )
 
  call r8vec_print ( n, vecn, '  Projected vector in BASE coordinates:' )

  call r8vec_print ( m, vecnm, '  Projected vector in original coordinates:' )
 
  return
end
subroutine test171 ( )

!*****************************************************************************80
!
!! TEST171 tests QUAD_AREA_2D, QUAD_AREA2_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 May 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ), dimension(2,4) :: q = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, 4 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST171'
  write ( *, '(a)' ) '  For a quadrilateral in 2D:'
  write ( *, '(a)' ) '  QUAD_AREA_2D finds the area;'
  write ( *, '(a)' ) '  QUAD_AREA2_2D finds the area;'

  call r8mat_transpose_print ( 2, 4, q, '  The vertices:' )

  call quad_area_2d ( q, area )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  QUAD_AREA_2D area is  ', area

  call quad_area2_2d ( q, area )

  write ( *, '(a,g14.6)' ) '  QUAD_AREA2_2D area is ', area
 
  return
end
subroutine test1712 ( )

!*****************************************************************************80
!
!! TEST1712 tests QUAD_AREA_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 May 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ) area1
  real ( kind = 8 ) area2
  real ( kind = 8 ), dimension(3,4) :: q = reshape ( (/ &
    2.0D+00, 2.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, &
    3.0D+00, 3.0D+00, 1.0D+00  &
    /), (/ 3, 4 /) )
  real ( kind = 8 ) t(3,3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1712'
  write ( *, '(a)' ) '  For a quadrilateral in 3D:'
  write ( *, '(a)' ) '  QUAD_AREA_3D finds the area.'

  call r8mat_transpose_print ( 3, 4, q, '  The vertices:' )

  call quad_area_3d ( q, area )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  QUAD_AREA_3D area is     ', area

  t(1:3,1:3) = q(1:3,1:3)
  call triangle_area_3d ( t, area1 )
  t(1:3,1:2) = q(1:3,3:4)
  t(1:3,  3) = q(1:3,1  )
  call triangle_area_3d ( t, area2 )
  write ( *, '(a,g14.6)' ) '  Sum of TRIANGLE_AREA_3D: ', area1 + area2

  return
end
subroutine test1715 ( )

!*****************************************************************************80
!
!! TEST1715 tests QUAD_CONTAINS_POINT_2D, QUAD_POINT_DIST_2D, QUAD_POINT_DIST_SIGNED_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 7

  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_signed
  logical inside
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = reshape ( (/ &
     0.25D+00,   0.25D+00, &
     0.75D+00,   0.25D+00, &
     1.00D+00,   1.00D+00, &
    11.00D+00,   0.50D+00, &
     0.00D+00,   0.50D+00, &
     0.50D+00, -10.00D+00, &
     2.00D+00,   2.00D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension(dim_num,4) :: q = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 1.0D+00 /), (/ dim_num, 4 /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1715'
  write ( *, '(a)' ) '  For a quadrilateral in 2D:'
  write ( *, '(a)' ) '  QUAD_AREA_2D finds the area;'
  write ( *, '(a)' ) '  QUAD_CONTAINS_POINT_2D tells if a point is inside;'
  write ( *, '(a)' ) '  QUAD_POINT_DIST_2D computes the distance.'
  write ( *, '(a)' ) '  QUAD_POINT_DIST_SIGNED_2D computes signed distance.'

  call r8mat_transpose_print ( dim_num, 4, q, '  The vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        P        Contains  Dist    Dist'
  write ( *, '(a)' ) '                          Signed  Unsigned'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call quad_contains_point_2d ( q, p, inside )

    call quad_point_dist_signed_2d ( q, p, dist_signed )

    call quad_point_dist_2d ( q, p, dist )

    write ( *, '(2x,2g14.6,2x,l1,2x,2f12.4)' ) &
      p(1:dim_num), inside, dist_signed, dist

  end do
 
  return
end
subroutine test172 ( )

!*****************************************************************************80
!
!! TEST172 tests QUAT_CONJ, QUAT_INV, QUAT_MUL, and QUAT_NORM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 4

  real ( kind = 8 ), dimension ( dim_num ) :: q1 = (/ &
    2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00 /)
  real ( kind = 8 ) q2(dim_num)
  real ( kind = 8 ) q3(dim_num)
  real ( kind = 8 ) quat_norm

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST172'
  write ( *, '(a)' ) '  QUAT_CONJ conjugates a quaternion;'
  write ( *, '(a)' ) '  QUAT_INV inverts a quaternion;'
  write ( *, '(a)' ) '  QUAT_MUL multiplies quaternions.'
  write ( *, '(a)' ) '  QUAT_NORM computes the norm.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,4g14.6)' ) '  Q1 =            ', q1(1:dim_num)

  write ( *, '(a,g14.6)' ) '  Norm ( Q1 ) = ', quat_norm ( q1 )

  call quat_conj ( q1, q2 )

  write ( *, '(a,4g14.6)' ) '  Q2 = conj(Q1) = ', q2(1:dim_num)

  call quat_mul ( q1, q2, q3 )

  write ( *, '(a,4g14.6)' ) '  Q3 = Q1*Q2 =    ', q3(1:dim_num)

  call quat_inv ( q1, q2 )

  write ( *, '(a,4g14.6)' ) '  Q2 = inv(Q1) =  ', q2(1:dim_num)

  call quat_mul ( q1, q2, q3 )

  write ( *, '(a,4g14.6)' ) '  Q3 = Q1*Q2 =    ', q3(1:dim_num)

  return
end
subroutine test173 ( )

!*****************************************************************************80
!
!! TEST173 tests RADEC_DISTANCE_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 6

  real ( kind = 8 ) dec1
  real ( kind = 8 ) dec2
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = reshape ( (/ &
     1.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00,  1.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00,  1.0D+00, &
     1.0D+00,  1.0D+00,  1.0D+00, &
     5.0D+00, -2.0D+00, -1.0D+00, &
    -2.0D+00, -2.0D+00, -2.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) ra1
  real ( kind = 8 ) ra2
  real ( kind = 8 ) radians_to_degrees
  integer ( kind = 4 ) test1
  integer ( kind = 4 ) test2
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta_deg

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST173'
  write ( *, '(a)' ) '  RADEC_DISTANCE_3D computes the angular separation'
  write ( *, '(a)' ) '  between two points on a sphere described in terms of'
  write ( *, '(a)' ) '  right ascension and declination.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '     RA1       DEC1      RA2       DEC2    Radians  Degrees'
  write ( *, '(a)' ) ' '

  do test1 = 1, test_num

    p1(1:dim_num) = p_test(1:dim_num,test1)

    call xyz_to_radec ( p1, ra1, dec1 )

    do test2 = test1+1, test_num

      p2(1:dim_num) = p_test(1:dim_num,test2)

      call xyz_to_radec ( p2, ra2, dec2 )
      call radec_distance_3d ( ra1, dec1, ra2, dec2, theta )
      theta_deg = radians_to_degrees ( theta )
      write ( *, '(2x,6f10.4)' ) ra1, dec1, ra2, dec2, theta, theta_deg

    end do

  end do

  return
end
subroutine test174 ( )

!*****************************************************************************80
!
!! TEST174 tests RADEC_TO_XYZ and XYZ_TO_RADEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 6

  real ( kind = 8 ) dec
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = reshape ( (/ &
     1.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00,  1.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00,  1.0D+00, &
     1.0D+00,  1.0D+00,  1.0D+00, &
     5.0D+00, -2.0D+00, -1.0D+00, &
    -2.0D+00, -2.0D+00, -2.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) ra
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST174'
  write ( *, '(a)' ) '  RADEC_TO_XYZ converts XYZ to RADEC coordinates.'
  write ( *, '(a)' ) '  XYZ_TO_RADEC converts RADEC to XYZ coordinates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          P1          RA     DEC           P2'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p1(1:dim_num) = p_test(1:dim_num,test)

    call xyz_to_radec ( p1, ra, dec )
    call radec_to_xyz ( ra, dec, p2 )

    write ( *, '(2x,8f7.3)' ) p1(1:dim_num), ra, dec, p2(1:dim_num)

  end do

  return
end
subroutine test1745 ( )

!*****************************************************************************80
!
!! TEST1745 tests R8MAT_SOLVE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: rhs_num = 2

  real ( kind = 8 ), dimension (n,n+rhs_num) :: a = reshape ( (/ &
     1.0D+00,  4.0D+00,  7.0D+00, &
     2.0D+00,  5.0D+00,  8.0D+00, &
     3.0D+00,  6.0D+00,  0.0D+00, &
    14.0D+00, 32.0D+00, 23.0D+00, &
     7.0D+00, 16.0D+00,  7.0D+00 /), (/ n, n+rhs_num /) )
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1745'
  write ( *, '(a)' ) '  R8MAT_SOLVE solves linear systems.'
  write ( *, '(a)' ) ' '
!
!  Print out the matrix to be inverted.
!
  call r8mat_print ( n, n+rhs_num, a, '  The linear system:' )
!
!  Solve the systems.
!
  call r8mat_solve ( n, rhs_num, a, info )
 
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The input matrix was singular.'
    write ( *, '(a)' ) '  The solutions could not be computed.'
    write ( *, '(a)' ) ' '
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The computed solutions:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,2g14.6)' ) a(i,n+1:n+rhs_num)
  end do
 
  return
end
subroutine test1746 ( )

!*****************************************************************************80
!
!! TEST1746 tests R8MAT_INVERSE_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ), dimension ( n, n ) :: a = reshape ( (/ &
    3.0D+00, 2.0D+00, 0.0D+00, &
    2.0D+00, 2.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00 /), (/ n, n /) )
  real ( kind = 8 ) b(n,n)
  real ( kind = 8 ) det

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1746'
  write ( *, '(a)' ) '  R8MAT_INVERSE_3D inverts a 3 by 3 matrix.'

  call r8mat_print ( n, n, a, '  Matrix A:' )

  call r8mat_inverse_3d ( a, b, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Determinant of A is ', det

  call r8mat_print ( n, n, b, '  Inverse matrix B:' )

  return
end
subroutine test175 ( )

!*****************************************************************************80
!
!! TEST175 tests ROTATION_AXIS_VECTOR_3D, ROTATION_MAT_VECTOR_3D, ROTATION_QUAT_VECTOR_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) a(dim_num,dim_num)
  real ( kind = 8 ) :: angle = 1.159804D+00
  real ( kind = 8 ), dimension(dim_num) :: axis = (/ &
    0.2361737D+00, -0.8814124D+00, -0.4090649D+00 /)
  real ( kind = 8 ) degrees_to_radians
  real ( kind = 8 ) q(4)
  real ( kind = 8 ), dimension(dim_num) :: v = (/ &
    1.0D+00, 4.0D+00, 10.0D+00 /)
  real ( kind = 8 ) w(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST175'
  write ( *, '(a)' ) '  ROTATION_AXIS_VECTOR_3D applies an axis'
  write ( *, '(a)' ) '  rotation to a vector;'
  write ( *, '(a)' ) '  ROTATION_MAT_VECTOR_3D applies a matrix'
  write ( *, '(a)' ) '  rotation to a vector.'
  write ( *, '(a)' ) '  ROTATION_QUAT_VECTOR_3D applies a quaternion'
  write ( *, '(a)' ) '  rotation to a vector.'

  call r8vec_print ( dim_num, v, '  The vector:' )

  call r8vec_print ( dim_num, axis, '  The rotation axis:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The rotation angle is ', angle

  call rotation_axis_vector_3d ( axis, angle, v, w )

  call r8vec_print ( dim_num, w, '  The rotated vector:' )

  call rotation_axis2mat_3d ( axis, angle, a )

  call r8mat_print ( dim_num, dim_num, a, '  The rotation matrix:' )

  call rotation_mat_vector_3d ( a, v, w )

  call r8vec_print ( dim_num, w, '  The rotated vector:' )

  call rotation_axis2quat_3d ( axis, angle, q )

  call r8vec_print ( 4, q, '  The rotation quaternion:' )

  call rotation_quat_vector_3d ( q, v, w )

  call r8vec_print ( dim_num, w, '  The rotated vector:' )
!
!  Another test of ROTATION_AXIS_VECTOR_3D with an axis vector
!  that does not have unit length.
!
  v(1:3) = (/ 1.0D+00, 1.0D+00, 1.0D+00 /)

  call r8vec_print ( dim_num, v, '  The vector:' )

  axis(1:3) = (/ 0.0D+00, 0.0D+00, 2.0D+00 /)

  call r8vec_print ( dim_num, axis, '  The rotation axis:' )

  angle = 90.0D+00
  angle = degrees_to_radians ( angle )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The rotation angle is ', angle

  call rotation_axis_vector_3d ( axis, angle, v, w )

  call r8vec_print ( dim_num, w, '  The rotated vector:' )

  return
end
subroutine test176 ( )

!*****************************************************************************80
!
!! TEST176 tests ROTATION_AXIS2MAT_3D and ROTATION_MAT2AXIS_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ), dimension(dim_num,dim_num) :: a = reshape ( (/ &
    0.43301269D+00, -0.5D+00,        0.75D+00, &
    0.25D+00,        0.86602539D+00, 0.43301269D+00, &
   -0.86602539D+00,  0.0D+00,        0.5D+00 /), (/ dim_num, dim_num /) )
  real ( kind = 8 ) angle
  real ( kind = 8 ) axis(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST176'
  write ( *, '(a)' ) '  ROTATION_MAT2AXIS_3D computes a rotation axis'
  write ( *, '(a)' ) '  and angle from a rotation matrix.'
  write ( *, '(a)' ) '  ROTATION_AXIS2MAT_3D computes a rotation matrix'
  write ( *, '(a)' ) '  from a rotation axis and angle.'

  call r8mat_print ( dim_num, dim_num, a, '  The rotation matrix:' )

  call rotation_mat2axis_3d ( a, axis, angle )

  call r8vec_print ( dim_num, axis, '  The rotation axis:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The rotation angle is ', angle

  call rotation_axis2mat_3d ( axis, angle, a )

  call r8mat_print ( dim_num, dim_num, a, '  The rotation matrix:' )

  return
end
subroutine test177 ( )

!*****************************************************************************80
!
!! TEST177 tests ROTATION_AXIS2QUAT_3D and ROTATION_QUAT2AXIS_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) :: angle = 1.159804D+00
  real ( kind = 8 ), dimension(dim_num) :: axis = (/ &
    0.2361737D+00, -0.8814124D+00, -0.4090649D+00 /)
  real ( kind = 8 ) q(4)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST177'
  write ( *, '(a)' ) '  ROTATION_QUAT2AXIS_3D computes a rotation axis'
  write ( *, '(a)' ) '  and angle from a rotation quaternion.'
  write ( *, '(a)' ) '  ROTATION_AXIS2QUAT_3D computes a rotation'
  write ( *, '(a)' ) '  quaternion from a rotation axis and angle.'

  call r8vec_print ( dim_num, axis, '  Rotation axis:' )

  write ( *, '(a,g14.6)' ) '  Rotation angle is ', angle

  call rotation_axis2quat_3d ( axis, angle, q )

  call r8vec_print ( 4, q, '  The rotation quaternion:' )

  call rotation_quat2axis_3d ( q, axis, angle )

  call r8vec_print ( dim_num, axis, '  The rotation axis:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The rotation angle is ', angle

  return
end
subroutine test178 ( )

!*****************************************************************************80
!
!! TEST178 tests ROTATION_MAT2QUAT_3D and ROTATION_QUAT2MAT_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ), dimension(dim_num,dim_num) :: a = reshape ( (/ &
    0.43301269D+00, -0.5D+00,        0.75D+00, &
    0.25D+00,        0.86602539D+00, 0.43301269D+00, &
   -0.86602539D+00,  0.0D+00,        0.5D+00 /), (/ dim_num, dim_num /) )
  real ( kind = 8 ) q(4)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST178'
  write ( *, '(a)' ) '  ROTATION_MAT2QUAT_3D computes a rotation'
  write ( *, '(a)' ) '  quaternion from a rotation matrix.'
  write ( *, '(a)' ) '  ROTATION_QUAT2MAT_3D computes a rotation matrix'
  write ( *, '(a)' ) '  from a rotation quaternion.'

  call r8mat_print ( dim_num, dim_num, a, '  The rotation matrix:' )

  call rotation_mat2quat_3d ( a, q )

  call r8vec_print ( 4, q, '  The rotation quaternion:' )

  call rotation_quat2mat_3d ( q, a )

  call r8mat_print ( dim_num, dim_num, a, '  The rotation matrix:' )

  return
end
subroutine test1787 ( )

!*****************************************************************************80
!
!! TEST1787 tests R8GE_FA and R8GE_SL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) alu(n,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) job
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1787'
  write ( *, '(a)' ) '  DGE_FA factors a general linear system,'
  write ( *, '(a)' ) '  DGE_SL solves a factored system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  seed = 123456789
  call r8mat_uniform_01 ( n, n, seed, a )
  call r8mat_print ( n, n, a, '  Matrix A.' )
!
!  Set the desired solution.
!
  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do
!
!  Compute the corresponding right hand side.
!
  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )
!
!  Make a copy of the matrix.
!
  alu(1:n,1:n) = a(1:n,1:n)
!
!  Factor the matrix.
!
  call r8ge_fa ( n, alu, pivot, info )

  call r8mat_print ( n, n, alu, '  Factored ALU matrix.' )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Fatal error!'
    write ( *, '(a)' ) '  DGE_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call r8ge_sl ( n, alu, pivot, b, job )

  call r8vec_print ( n, b, '  Solution: (Should be 1, 2, 3,...)' )
!
!  Set another the desired solution.
!
  x(1:n) = 1.0D+00
!
!  Compute the corresponding right hand side.
!
  b(1:n) = matmul ( a(1:n,1:n), x(1:n) )
!
!  Solve the system
!
  job = 0
  call r8ge_sl ( n, alu, pivot, b, job )

  call r8vec_print ( n, b, '  Solution: (Should be 1, 1, 1,...)' )
!
!  Set the desired solution.
!
  do i = 1, n
    x(i) = real ( i, kind = 8 )
  end do
!
!  Compute the corresponding right hand side.
!
  b(1:n) = matmul ( transpose ( a(1:n,1:n) ), x(1:n) )
!
!  Solve the system.
!
  job = 1
  call r8ge_sl ( n, alu, pivot, b, job )

  call r8vec_print ( n, b, &
    '  Solution of transposed system: (Should be 1, 2, 3,...)' )

  return
end
subroutine test1893 ( )

!*****************************************************************************80
!
!! TEST1893 tests RTP_TO_XYZ and XYZ_TO_RTP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) :: a = -2.0D+00
  real ( kind = 8 ) :: b =  3.0D+00
  real ( kind = 8 ) phi
  real ( kind = 8 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5
  real ( kind = 8 ) theta
  real ( kind = 8 ) xyz1(3)
  real ( kind = 8 ) xyz2(3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1893'
  write ( *, '(a)' ) '  RTP_TO_XYZ converts XYZ to (R,Theta,Phi) coordinates.'
  write ( *, '(a)' ) '  XYZ_TO_RTP converts (R,Theta,Phi) to XYZ coordinates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      X1     Y1     Z1     R    THETA    PHI    X2     Y2     Z2'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do test = 1, test_num

    call r8vec_uniform_ab ( 3, a, b, seed, xyz1 )

    call xyz_to_rtp ( xyz1, r, theta, phi )
    call rtp_to_xyz ( r, theta, phi, xyz2 )

    write ( *, '(2x,9f7.3)' ) xyz1(1:3), r, theta, phi, xyz2(1:3)

  end do

  return
end
subroutine test036 ( )

!*****************************************************************************80
!
!! TEST036 tests SEGMENT_CONTAINS_POINT_1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) p
  real ( kind = 8 ), dimension ( test_num ) :: p_test = (/ &
    3.0D+00,   7.5D+00, 20.0D+00,  5.0D+00 /)
  real ( kind = 8 ) p1
  real ( kind = 8 ), dimension ( test_num ) :: p1_test = (/ &
    2.0D+00,  10.0D+00,  8.0D+00, 88.0D+00 /)
  real ( kind = 8 ) p2
  real ( kind = 8 ), dimension ( test_num ) :: p2_test = (/ &
    6.0D+00, -10.0D+00, 10.0D+00, 88.0D+00 /)
  real ( kind = 8 ) t
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST036'
  write ( *, '(a)' ) '  SEGMENT_CONTAINS_POINT_1D determines if a point'
  write ( *, '(a)' ) '  lies within a line segment in 1D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       P1     P       T'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p1 = p1_test(test)
    p2 = p2_test(test)
    p = p_test(test)

    call segment_contains_point_1d ( p1, p2, p, t )
    write ( *, '(2x,3f7.2,g14.6)' ) p1, p2, p, t

  end do

  return
end
subroutine test0364 ( )

!*****************************************************************************80
!
!! TEST0364 tests SEGMENT_POINT_COORDS_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 6

  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = reshape ( (/ &
    3.0D+00, 1.0D+00, &
    4.0D+00, 1.0D+00, &
    100.0D+00, 1.0D+00, &
      5.0D+00, 100.0D+00, &
      7.0D+00, -1.0D+00, &
      0.0D+00, 5.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) s
  real ( kind = 8 ), dimension ( test_num ) :: s_test = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 99.0D+00, 2.0D+00, 4.0D+00 /)
  real ( kind = 8 ) t
  real ( kind = 8 ), dimension ( test_num ) :: t_test = (/ &
    0.0D+00, 0.25D+00, 24.250D+00, 0.50D+00, 1.0D+00, -0.75D+00 /)
  integer ( kind = 4 ) test

  p1(1:dim_num) = (/ 3.0D+00, 1.0D+00 /)
  p2(1:dim_num) = (/ 7.0D+00, 1.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0364'
  write ( *, '(a)' ) '  SEGMENT_POINT_COORDS_2D computes coordinates'
  write ( *, '(a)' ) '  (S,T) for a point relative to a line segment in 2D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2f9.4)' ) '  P1 =   ', p1(1:dim_num)
  write ( *, '(a,2f9.4)' ) '  P2 =   ', p2(1:dim_num)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '            P1              S       T'
  write ( *, '(a)' ) '   -----------------     ------   ------'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call segment_point_coords_2d ( p1, p2, p, s, t )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,2f9.4,2x,2f9.4,2x,a)' ) p(1:dim_num), s, t, 'Computed'
    write ( *, '(2x,18x,2x,2f9.4,2x,a)' ) s_test(test), t_test(test), 'Expected'

  end do

  return
end
subroutine test0365 ( )

!*****************************************************************************80
!
!! TEST0365 tests SEGMENT_POINT_DIST_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) dist
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0365'
  write ( *, '(a)' ) '  SEGMENT_POINT_DIST_2D computes the distance'
  write ( *, '(a)' ) '  between a line segment and point in 2D.'

  do test = 1, test_num

    call r8vec_uniform_01 ( dim_num, seed, p1 )
    call r8vec_uniform_01 ( dim_num, seed, p2 )
    call r8vec_uniform_01 ( dim_num, seed, p )

    call segment_point_dist_2d ( p1, p2, p, dist )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2)' )    '  TEST = ', test
    write ( *, '(a,2f9.4)' ) '  P1 =   ', p1(1:dim_num)
    write ( *, '(a,2f9.4)' ) '  P2 =   ', p2(1:dim_num)
    write ( *, '(a,2f9.4)' ) '  P =    ', p(1:dim_num)
    write ( *, '(a, f9.4)' ) '  DIST = ', dist

  end do

  return
end
subroutine test0366 ( )

!*****************************************************************************80
!
!! TEST0366 tests SEGMENT_POINT_DIST_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) dist
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0366'
  write ( *, '(a)' ) '  SEGMENT_POINT_DIST_3D computes the distance'
  write ( *, '(a)' ) '  between a line segment and point in 3D.'

  do test = 1, test_num

    call r8vec_uniform_01 ( dim_num, seed, p1 )
    call r8vec_uniform_01 ( dim_num, seed, p2 )
    call r8vec_uniform_01 ( dim_num, seed, p )

    call segment_point_dist_3d ( p1, p2, p, dist )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2)' )    '  TEST = ', test
    write ( *, '(a,3f9.4)' ) '  P1 =   ', p1(1:dim_num)
    write ( *, '(a,3f9.4)' ) '  P2 =   ', p2(1:dim_num)
    write ( *, '(a,3f9.4)' ) '  P =    ', p(1:dim_num)
    write ( *, '(a, f9.4)' ) '  DIST = ', dist

  end do

  return
end
subroutine test0367 ( )

!*****************************************************************************80
!
!! TEST0367 tests SEGMENT_POINT_NEAR_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) dist
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) pn(dim_num)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) t
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0367'
  write ( *, '(a)' ) '  SEGMENT_POINT_NEAR_2D computes the nearest point'
  write ( *, '(a)' ) '  on a line segment to a point in 2D.'

  do test = 1, test_num

    call r8vec_uniform_01 ( dim_num, seed, p1 )
    call r8vec_uniform_01 ( dim_num, seed, p2 )
    call r8vec_uniform_01 ( dim_num, seed, p )

    call segment_point_near_2d ( p1, p2, p, pn, dist, t )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2)' )    '  TEST = ', test
    write ( *, '(a,2f9.4)' ) '  P1 =   ', p1(1:dim_num)
    write ( *, '(a,2f9.4)' ) '  P2 =   ', p2(1:dim_num)
    write ( *, '(a,2f9.4)' ) '  P  =   ', p(1:dim_num)
    write ( *, '(a,2f9.4)' ) '  PN =   ', pn(1:dim_num)
    write ( *, '(a, f9.4)' ) '  DIST = ', dist
    write ( *, '(a, f9.4)' ) '  T =    ', t

  end do

  return
end
subroutine test0368 ( )

!*****************************************************************************80
!
!! TEST0368 tests SEGMENT_POINT_NEAR_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) dist
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) pn(dim_num)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) t
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0368'
  write ( *, '(a)' ) '  SEGMENT_POINT_NEAR_3D computes the nearest point'
  write ( *, '(a)' ) '  on a line segment to a point in 3D.'

  do test = 1, test_num

    call r8vec_uniform_01 ( dim_num, seed, p1 )
    call r8vec_uniform_01 ( dim_num, seed, p2 )
    call r8vec_uniform_01 ( dim_num, seed, p )

    call segment_point_near_3d ( p1, p2, p, pn, dist, t )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2)' )    '  TEST = ', test
    write ( *, '(a,3f9.4)' ) '  P1 =   ', p1(1:dim_num)
    write ( *, '(a,3f9.4)' ) '  P2 =   ', p2(1:dim_num)
    write ( *, '(a,3f9.4)' ) '  P  =   ', p(1:dim_num)
    write ( *, '(a,3f9.4)' ) '  PN =   ', pn(1:dim_num)
    write ( *, '(a, f9.4)' ) '  DIST = ', dist
    write ( *, '(a, f9.4)' ) '  T =    ', t

  end do

  return
end
subroutine test037 ( )

!*****************************************************************************80
!
!! TEST037 tests SEGMENT_POINT_NEAR_3D.
!
!  Discussion:
!
!    Case 1, point is nearest end of segment.
!
!      LS: (2,3,0) + t * (2,1,0) for t = 0 to 3.
!      P (11,6,4)
!      Distance is 5.
!
!    Case 2, point is nearest interior point of segment.
!
!      LS: (2,3,0) + t * (2,1,0) for t = 0 to 3.
!      P (4,4,1)
!      Distance is 1.
!
!    Case 3, point is on the line.
!
!      LS: (2,3,0) + t * (2,1,0) for t = 0 to 3.
!      P (6,5,0)
!      Distance is 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) dist
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = reshape ( (/ &
    11.0D+00, 6.0D+00, 4.0D+00, &
     4.0D+00, 4.0D+00, 1.0D+00, &
     6.0D+00, 5.0D+00, 0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension(dim_num) :: p1 = (/ 2.0D+00, 3.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p2 = (/ 8.0D+00, 6.0D+00, 0.0D+00 /)
  real ( kind = 8 ) pn(dim_num)
  real ( kind = 8 ) t
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST037'
  write ( *, '(a)' ) '  SEGMENT_POINT_NEAR_3D computes the nearest'
  write ( *, '(a)' ) '  point on a line segment, to a given point,'
  write ( *, '(a)' ) '  in 3D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test  T   Distance   PN.'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call segment_point_near_3d ( p1, p2, p, pn, dist, t )
 
    write ( *, '(2x,i2,5f9.4)' ) test, t, dist, pn(1:dim_num)

  end do

  return
end
subroutine test1788 ( )

!*****************************************************************************80
!
!! TEST17888 tests SIMPLEX_LATTICE_LAYER_POINT_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) layer
  logical              more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) :: n_test(test_num) = (/ 1, 2, 3, 4 /)
  integer ( kind = 4 ) test
  integer ( kind = 4 ), allocatable :: v(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1788'
  write ( *, '(a)' ) '  SIMPLEX_LATTICE_LAYER_POINT_NEXT returns the next'
  write ( *, '(a)' ) '  point in an N-dimensional simplex lattice layer defined by:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    C(N+1) - 1 <= X(1)/C(1) + X(2)/C(2) + ... + X(N)/C(N) <= C(N+1).'

  do test = 1, test_num

    n = n_test(test)

    allocate ( c(1:n+1) )
    allocate ( v(1:n) )

    do i = 1, n
      c(i) = i + 1
    end do
    v(1:n) = 0

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  N = ', n
    write ( *, '(a)', ADVANCE = 'NO' ) '  C =       '
    do i = 1, n
      write ( *, '(2x,i4)', ADVANCE = 'NO' ) c(i)
    end do
    write ( *, '(a)', ADVANCE = 'YES' ) 
    write ( *, '(a)' ) ' '

    do layer = 0, 2

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Layer ', layer
      write ( *, '(a)' ) ' '

      c(n+1) = layer
      more = .false.
      i = 0

      do
        call simplex_lattice_layer_point_next ( n, c, v, more )
        if ( .not. more ) then
          write ( *, '(a)' ) '  No more.'
          exit
        end if
        i = i + 1
        write ( *, '(2x,i4,6x,10(2x,i4))' ) i, v(1:n)
      end do

    end do

    deallocate ( c )
    deallocate ( v )

  end do

  return
end
subroutine test1789 ( )

!*****************************************************************************80
!
!! TEST1789 tests SIMPLEX_LATTICE_POINT_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ), allocatable :: c(:)
  integer ( kind = 4 ) i
  logical              more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) :: n_test(test_num) = (/ 1, 2, 3, 4 /)
  integer ( kind = 4 ) test
  integer ( kind = 4 ), allocatable :: v(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1789'
  write ( *, '(a)' ) '  SIMPLEX_LATTICE_POINT_NEXT returns the next lattice'
  write ( *, '(a)' ) '  point in an N-dimensional simplex defined by:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    0 <= X(1)/C(1) + X(2)/C(2) + ... + X(N)/C(N) <= C(N+1).'

  do test = 1, test_num

    n = n_test(test)

    allocate ( c(1:n+1) )
    allocate ( v(1:n) )

    do i = 1, n + 1
      c(i) = n + 2 - i
    end do
    v(1:n) = 0
    more = .false.

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  N = ', n
    write ( *, '(a)', ADVANCE = 'NO' ) '  C =       '
    do i = 1, n + 1
      write ( *, '(2x,i4)', ADVANCE = 'NO' ) c(i)
    end do
    write ( *, '(a)', ADVANCE = 'YES' ) 
    write ( *, '(a)' ) ' '

    i = 0

    do
      call simplex_lattice_point_next ( n, c, v, more )
      if ( .not. more ) then
        write ( *, '(a)' ) '  No more.'
        exit
      end if
      i = i + 1
      write ( *, '(2x,i4,6x,10(2x,i4))' ) i, v(1:n)
    end do

    deallocate ( c )
    deallocate ( v )

  end do

  return
end
subroutine test179 ( )

!*****************************************************************************80
!
!! TEST179 tests SOCCER_SIZE_3D and SOCCER_SHAPE_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) area
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ), allocatable, dimension ( : ) ::  face_order
  integer ( kind = 4 ) face_order_max
  integer ( kind = 4 ), allocatable, dimension ( :, : ) ::  face_point
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) normal(dim_num)
  integer ( kind = 4 ) point_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_coord
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: v
  real ( kind = 8 ) vave(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST179'
  write ( *, '(a)' ) '  For the truncated icosahedron, or soccer ball,'
  write ( *, '(a)' ) '  SOCCER_SIZE_3D returns dimension information;'
  write ( *, '(a)' ) '  SOCCER_SHAPE_3D returns face and order information.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We will use this information to compute the'
  write ( *, '(a)' ) '  areas and centers of each face.'

  call soccer_size_3d ( point_num, edge_num, face_num, face_order_max )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of vertices = ', point_num
  write ( *, '(a,i8)' ) '  Number of edges =    ', edge_num
  write ( *, '(a,i8)' ) '  Number of faces =    ', face_num
  write ( *, '(a,i8)' ) '  Maximum face order = ', face_order_max

  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )
  allocate ( point_coord(1:dim_num,1:point_num) )
  allocate ( v(1:dim_num,1:face_order_max) )

  call soccer_shape_3d ( point_num, face_num, face_order_max, point_coord, &
    face_order, face_point )
!
!  Compute the area of each face.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Face  Order  Area'
  write ( *, '(a)' ) ' '

  do i = 1, face_num

    do j = 1, face_order(i)
      k = face_point(j,i)
      v(1:dim_num,j) = point_coord(1:dim_num,k)
    end do

    call polygon_area_3d ( face_order(i), v, area, normal )

    write ( *, '(2x,i8,i7,f8.4)' ) i, face_order(i), area

  end do
!
!  Find the center of each face.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Face  Center'
  write ( *, '(a)' ) ' '

  do i = 1, face_num

    vave(1:dim_num) = 0.0D+00

    do j = 1, face_order(i)
      k = face_point(j,i)
      vave(1:dim_num) = vave(1:dim_num) + point_coord(1:dim_num,k)
    end do

    vave(1:dim_num) = vave(1:dim_num) / real ( face_order(i), kind = 8 )

    write ( *, '(2x,i8,3f8.4)' ) i, vave(1:dim_num)

  end do

  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )
  deallocate ( v )

  return
end
subroutine test180 ( )

!*****************************************************************************80
!
!! TEST180 tests SORT_HEAP_EXTERNAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST180'
  write ( *, '(a)' ) '  SORT_HEAP_EXTERNAL sorts objects externally.'
  write ( *, '(a)' ) ' '

  indx = 0
  i = 0
  j = 0
  isgn = 0

  call i4vec_uniform ( n, 1, n, seed, a )
 
  call i4vec_print ( n, a, '  Unsorted array' )
  
  do

    call sort_heap_external ( n, indx, i, j, isgn )
 
    if ( indx < 0 ) then
      if ( a(i) <= a(j) ) then
        isgn = -1
      else
        isgn = +1
      end if
    else if ( 0 < indx ) then
      call i4_swap ( a(i), a(j) )
    else
      exit
    end if

  end do

  call i4vec_print ( n, a, '  Sorted array' )
 
  return
end
subroutine test1804 ( )

!*****************************************************************************80
!
!! TEST1804 tests SIMPLEX_UNIT_LATTICE_POINT_NUM_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 11

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), dimension ( test_num ) :: dim_num_test = (/ &
    2, 2, 2, 2, 3, 3, 3, 3, 4, 5, 6 /)
  integer ( kind = 4 ) t
  integer ( kind = 4 ), dimension ( test_num ) :: t_test = (/ &
    1, 2, 3, 4, 1, 2, 3, 10, 3, 3, 3 /)
  integer ( kind = 4 ) test
  integer ( kind = 4 ) volume

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1804'
  write ( *, '(a)' ) '  For an N-dimensional unit simplex'
  write ( *, '(a)' ) '    0 <= X(1:N),'
  write ( *, '(a)' ) '    sum X(1:N) <= T'
  write ( *, '(a)' ) '  where T is an integer,'
  write ( *, '(a)' ) '  SIMPLEX_UNIT_LATTICE_POINT_NUM_ND computes the lattice volume,'
  write ( *, '(a)' ) '  that is, the number of lattice points it contains.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N         T    Volume'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    dim_num = dim_num_test(test)
    t = t_test(test)
    call simplex_unit_lattice_point_num_nd ( dim_num, t, volume )

    write ( *, '(2x,i8,2x,i8,2x,i8)' ) dim_num, t, volume

  end do

  return
end
subroutine test1805 ( )

!*****************************************************************************80
!
!! TEST1805 tests SIMPLEX_VOLUME_ND and TETRAHEDRON_VOLUME_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ), dimension(dim_num,4) :: tetra = reshape ( (/&
     0.000000D+00,  0.942809D+00, -0.333333D+00, &
    -0.816496D+00, -0.816496D+00, -0.333333D+00, &
     0.816496D+00, -0.816496D+00, -0.333333D+00, &
     0.000000D+00,  0.000000D+00,  1.000000D+00 /), (/ dim_num, 4 /) )
  real ( kind = 8 ) volume

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1805'
  write ( *, '(a)' ) '  For an N-dimensional simplex,'
  write ( *, '(a)' ) '  SIMPLEX_VOLUME_ND computes the volume.'
  write ( *, '(a)' ) '  Here, we check the routine by comparing it'
  write ( *, '(a)' ) '  with TETRAHEDRON_VOLUME_3D.'

  call r8mat_transpose_print ( dim_num, 4, tetra, '  Simplex vertices:' )

  call tetrahedron_volume_3d ( tetra, volume )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Volume computed by TETRAHEDRON_VOLUME_3D:'
  write ( *, '(2x,g14.6)' ) volume

  call simplex_volume_nd ( dim_num, tetra, volume )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Volume computed by SIMPLEX_VOLUME_ND:'
  write ( *, '(2x,g14.6)' ) volume

  return
end
subroutine test181 ( )

!*****************************************************************************80
!
!! TEST181 tests SPHERE_DIA2IMP_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: p1 = (/ &
    -1.0D+00, -1.0D+00, 4.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p2 = (/ &
     5.0D+00,  7.0D+00, 4.0D+00 /)
  real ( kind = 8 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST181'
  write ( *, '(a)' ) '  SPHERE_DIA2IMP_3D converts a sphere from'
  write ( *, '(a)' ) '  diameter to implicit form.'

  call r8vec_print ( dim_num, p1, '  Point P1:' )
  call r8vec_print ( dim_num, p2, '  Point P2:' )

  call sphere_dia2imp_3d ( p1, p2, r, pc )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '    Radius: ', r
 
  call r8vec_print ( dim_num, pc, '  The center:' )

  return
end
subroutine test182 ( )

!*****************************************************************************80
!
!! TEST182 tests SPHERE_EXP_CONTAINS_POINT_3D and SPHERE_IMP_CONTAINS_POINT_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) i
  logical inside
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = reshape ( (/ &
    1.0D+00, 2.0D+00, 3.0D+00, &
    7.0D+00, 2.0D+00, 3.0D+00, &
    1.0D+00, 5.0D+00, 3.0D+00, &
    2.5D+00, 3.5D+00, 4.5D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num ) :: pc = (/ &
    1.0D+00, 2.0D+00, 3.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p1 = (/ &
    4.0D+00, 2.0D+00, 3.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p2 = (/ &
    1.0D+00, 5.0D+00, 3.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p3 = (/ &
    1.0D+00, 2.0D+00, 6.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p4 = (/ &
   -2.0D+00, 2.0D+00, 3.0D+00 /)
  real ( kind = 8 ) :: r = 3.0D+00
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST182'
  write ( *, '(a)' ) '  SPHERE_EXP_CONTAINS_POINT_3D determines if a'
  write ( *, '(a)' ) '  point is within an explicit sphere;'
  write ( *, '(a)' ) '  SPHERE_IMP_CONTAINS_POINT_3D determines if a'
  write ( *, '(a)' ) '  point is within an implicit sphere;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SPHERE_EXP_CONTAINS_POINT_3D:'
  write ( *, '(a)' ) '    Inside, P'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call sphere_exp_contains_point_3d ( p1, p2, p3, p4, p, inside )

    write ( *, '(2x,l1,3g14.6)' ) inside, p(1:dim_num)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SPHERE_IMP_CONTAINS_POINT_3D:'
  write ( *, '(a)' ) '    Inside, P'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call sphere_imp_contains_point_3d ( r, pc, p, inside )

    write ( *, '(2x,l1,3g14.6)' ) inside, p(1:dim_num)

  end do

  return
end
subroutine test183 ( )

!*****************************************************************************80
!
!! TEST183 tests SPHERE_EXP_POINT_NEAR_3D and SPHERE_IMP_POINT_NEAR_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = reshape ( (/ &
    1.0D+00, 2.0D+00, 3.0D+00, &
    7.0D+00, 2.0D+00, 3.0D+00, &
    1.0D+00, 5.0D+00, 3.0D+00, &
    2.5D+00, 3.5D+00, 4.5D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num ) :: p1 = (/ &
    4.0D+00, 2.0D+00, 3.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p2 = (/ &
    1.0D+00, 5.0D+00, 3.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p3 = (/ &
    1.0D+00, 2.0D+00, 6.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p4 = (/ &
   -2.0D+00, 2.0D+00, 3.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: pc = (/ &
    1.0D+00, 2.0D+00, 3.0D+00 /)
  real ( kind = 8 ) pn(dim_num)
  real ( kind = 8 ) :: r = 3.0D+00
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST183'
  write ( *, '(a)' ) '  SPHERE_EXP_POINT_NEAR_3D determines if a'
  write ( *, '(a)' ) '  point is within an explicit sphere;'
  write ( *, '(a)' ) '  SPHERE_IMP_POINT_NEAR_3D determines if a'
  write ( *, '(a)' ) '  point is within an implicit sphere;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Sphere radius ', r

  call r8vec_print ( dim_num, pc, '  Sphere center:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SPHERE_EXP_POINT_NEAR_3D:'
  write ( *, '(a)' ) '            P               PN'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call sphere_exp_point_near_3d ( p1, p2, p3, p4, p, pn )

    write ( *, '(2x,6f10.4)' ) p(1:dim_num), pn(1:dim_num)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SPHERE_IMP_POINT_NEAR_3D:'
  write ( *, '(a)' ) '         P           PN'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call sphere_imp_point_near_3d ( r, pc, p, pn )

    write ( *, '(2x,6f10.4)' ) p(1:dim_num), pn(1:dim_num)

  end do

  return
end
subroutine test1835 ( )

!*****************************************************************************80
!
!! TEST1835 tests SPHERE_EXP2IMP_3D and SPHERE_IMP2EXP_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ), dimension ( dim_num ) :: pc = (/ &
    1.0D+00, 2.0D+00, 3.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p1 = (/ &
    4.0D+00, 2.0D+00, 3.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p2 = (/ &
    1.0D+00, 5.0D+00, 3.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p3 = (/ &
    1.0D+00, 2.0D+00, 6.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: p4 = (/ &
   -2.0D+00, 2.0D+00, 3.0D+00 /)
  real ( kind = 8 ) :: r = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1835'
  write ( *, '(a)' ) '  SPHERE_EXP2IMP_3D: explicit sphere => implicit form;'
  write ( *, '(a)' ) '  SPHERE_IMP2EXP_3D: implicit sphere => explicit form.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initial form of explicit sphere:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,3g14.6)' ) p1(1:dim_num)
  write ( *, '(2x,3g14.6)' ) p2(1:dim_num)
  write ( *, '(2x,3g14.6)' ) p3(1:dim_num)
  write ( *, '(2x,3g14.6)' ) p4(1:dim_num)

  call sphere_exp2imp_3d ( p1, p2, p3, p4, r, pc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed form of implicit sphere:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Imputed radius = ', r

  call r8vec_print ( dim_num, pc, '  Imputed center' )

  call sphere_imp2exp_3d ( r, pc, p1, p2, p3, p4 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed form of explicit sphere:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,3g14.6)' ) p1(1:dim_num)
  write ( *, '(2x,3g14.6)' ) p2(1:dim_num)
  write ( *, '(2x,3g14.6)' ) p3(1:dim_num)
  write ( *, '(2x,3g14.6)' ) p4(1:dim_num)

  return
end
subroutine test1836 ( )

!*****************************************************************************80
!
!! TEST1836 tests SPHERE_EXP2IMP_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ), dimension (n,n+1) :: p = reshape ( (/ &
    4.0D+00, 2.0D+00, 3.0D+00, &
    1.0D+00, 5.0D+00, 3.0D+00, &
    1.0D+00, 2.0D+00, 6.0D+00, &
   -2.0D+00, 2.0D+00, 3.0D+00 /), (/ n, n + 1 /) )
  real ( kind = 8 ) pc(n)
  real ( kind = 8 ), dimension ( n ) :: pc_true = (/ &
    1.0D+00, 2.0D+00, 3.0D+00 /)
  real ( kind = 8 ) r
  real ( kind = 8 ) :: r_true = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1836'
  write ( *, '(a)' ) '  SPHERE_EXP2IMP_ND: explicit sphere => implicit form;'

  call r8mat_transpose_print ( n, n + 1, p, '  Initial form of explicit sphere:' )

  call sphere_exp2imp_nd ( n, p, r, pc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed form of implicit sphere:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Imputed radius = ', r
  write ( *, '(a,g14.6)' ) '  True radius =    ', r_true

  call r8vec_print ( n, pc, '  Imputed center' )

  call r8vec_print ( n, pc_true, '  True center' )

  return
end
subroutine test188 ( )

!*****************************************************************************80
!
!! TEST188 tests SPHERE_IMP_POINT_PROJECT_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ), dimension ( dim_num,test_num) :: p_test = reshape ( (/ &
    2.0D+00, 0.0D+00,  0.0D+00, &
    0.0D+00, 4.0D+00,  0.0D+00, &
    2.0D+00, 4.0D+00, 10.0D+00, &
    3.0D+00, 5.0D+00,  0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: pc = (/ &
    2.0D+00, 4.0D+00, 0.0D+00 /)
  real ( kind = 8 ) :: r = 2.0D+00
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST188'
  write ( *, '(a)' ) '  SPHERE_IMP_POINT_PROJECT_3D projects a 3D point'
  write ( *, '(a)' ) '  onto a sphere.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        P1       projection P2'
  write ( *, '(a)' ) ' '
    
  do test = 1, test_num

    p1(1:dim_num) = p_test(1:dim_num,test)

    call sphere_imp_point_project_3d ( r, pc, p1, p2 )

    write ( *, '(6g12.4)' ) p1(1:dim_num), p2(1:dim_num)

  end do
 
  return
end
subroutine test189 ( )

!*****************************************************************************80
!
!! TEST189 tests SPHERE_IMP_AREA_ND and SPHERE_IMP_VOLUME_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) area
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ), parameter :: r = 1.0D+00
  real ( kind = 8 ) volume

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST189'
  write ( *, '(a)' ) '  For a implicit sphere in N dimensions:'
  write ( *, '(a)' ) '  SPHERE_IMP_AREA_ND computes the area;'
  write ( *, '(a)' ) '  SPHERE_IMP_VOLUME_ND computes the volume.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  We use a radius of R = ', r
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    DIM_NUM Area    Volume'
  write ( *, '(a)' ) ' '

  do dim_num = 2, 10
    call sphere_imp_area_nd ( dim_num, r, area )
    call sphere_imp_volume_nd ( dim_num, r, volume )
    write ( *, '(2x,i3,2g14.6)' ) dim_num, area, volume
  end do

  return
end
subroutine test1895 ( )

!*****************************************************************************80
!
!! TEST1895 tests SPHERE_UNIT_AREA_ND and SPHERE_UNIT_AREA_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
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
  write ( *, '(a)' ) 'TEST1895:'
  write ( *, '(a)' ) '  SPHERE_UNIT_AREA_ND evaluates the area of the unit'
  write ( *, '(a)' ) '  sphere in N dimensions.'
  write ( *, '(a)' ) '  SPHERE_UNIT_AREA_VALUES returns some test values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     DIM_NUM    Exact          Computed'
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
subroutine test190 ( )

!*****************************************************************************80
!
!! TEST190 tests SPHERE_UNIT_SAMPLE_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) average(dim_num)
  real ( kind = 8 ) dot_average
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: sample_num = 1000
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) v(dim_num)
  real ( kind = 8 ) x(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST190'
  write ( *, '(a)' ) '  For the unit sphere in 2 dimensions (the circle):'
  write ( *, '(a)' ) '  SPHERE_UNIT_SAMPLE_2D samples;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A few sample values:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call sphere_unit_sample_2d ( seed, x )
    write ( *, '(2x,2f8.4)' ) x(1:dim_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of sample points = ', sample_num

  average(1:dim_num) = 0.0D+00

  do i = 1, sample_num
    call sphere_unit_sample_2d ( seed, x )
    average(1:dim_num) = average(1:dim_num) + x(1:dim_num)
  end do

  average(1:dim_num) = average(1:dim_num) / real ( sample_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now average the points, which should get a value'
  write ( *, '(a)' ) '  close to zero, and closer as sample_num increases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2f8.4)' ) '  Average:        ', average(1:dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now choose a random direction, sample the same'
  write ( *, '(a)' ) '  number of points, and compute the dot product with'
  write ( *, '(a)' ) '  the direction.'
  write ( *, '(a)' ) '  Take the absolute value of each dot product '
  write ( *, '(a)' ) '  and sum and average.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We expect a value near 2 / PI = 0.6366...'

  do j = 1, 5

    call sphere_unit_sample_2d ( seed, v )

    dot_average = 0.0D+00

    do i = 1, sample_num
      call sphere_unit_sample_2d ( seed, x )
      dot_average = dot_average &
        + abs ( dot_product ( x(1:dim_num), v(1:dim_num) ) )
    end do

    dot_average = dot_average / real ( sample_num, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,2f8.4)' ) '  V:                ', v(1:dim_num)
    write ( *, '(a, f8.4)' ) '  Average |(XdotV)| ', dot_average

  end do

  return
end
subroutine test191 ( )

!*****************************************************************************80
!
!! TEST191 tests SPHERE_UNIT_SAMPLE_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) average(dim_num)
  real ( kind = 8 ) dot_average
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: sample_num = 1000
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) v(dim_num)
  real ( kind = 8 ) x(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST191'
  write ( *, '(a)' ) '  For the unit sphere in 3 dimensions:'
  write ( *, '(a)' ) '  SPHERE_UNIT_SAMPLE_3D samples;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A few sample values:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call sphere_unit_sample_3d ( seed, x )
    write ( *, '(2x,3f8.4)' ) x(1:dim_num)
  end do

  average(1:dim_num) = 0.0D+00

  do i = 1, sample_num
    call sphere_unit_sample_3d ( seed, x )
    average(1:dim_num) = average(1:dim_num) + x(1:dim_num)
  end do

  average(1:dim_num) = average(1:dim_num) / real ( sample_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now average the points, which should get a value'
  write ( *, '(a)' ) '  close to zero, and closer as sample_num increases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,3f8.4)' ) '  Average:        ', average(1:dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now choose a random direction, sample the same'
  write ( *, '(a)' ) '  number of points, and compute the dot product with'
  write ( *, '(a)' ) '  the direction.'
  write ( *, '(a)' ) '  Take the absolute value of each dot product '
  write ( *, '(a)' ) '  and sum and average.'

  do j = 1, 5

    call sphere_unit_sample_3d ( seed, v )

    dot_average = 0.0D+00

    do i = 1, sample_num
      call sphere_unit_sample_3d ( seed, x )
      dot_average = dot_average &
        + abs ( dot_product ( x(1:dim_num), v(1:dim_num) ) )
    end do

    dot_average = dot_average / real ( sample_num, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,3f8.4)' ) '  V:                ', v(1:dim_num)
    write ( *, '(a, f8.4)' ) '  Average |(XdotV)| ', dot_average

  end do

  return
end
subroutine test192 ( )

!*****************************************************************************80
!
!! TEST192 tests SPHERE_UNIT_SAMPLE_3D_2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) average(dim_num)
  real ( kind = 8 ) dot_average
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: sample_num = 1000
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) v(dim_num)
  real ( kind = 8 ) x(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST192'
  write ( *, '(a)' ) '  For the unit sphere in 3 dimensions:'
  write ( *, '(a)' ) '  SPHERE_UNIT_SAMPLE_3D_2 samples;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Warning: SPHERE_UNIT_SAMPLE_3D_2 is NOT a good code!'
  write ( *, '(a)' ) '  I only implemented it for comparison.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A few sample values:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call sphere_unit_sample_3d_2 ( seed, x )
    write ( *, '(2x,3f8.4)' ) x(1:dim_num)
  end do

  average(1:dim_num) = 0.0D+00

  do i = 1, sample_num
    call sphere_unit_sample_3d_2 ( seed, x )
    average(1:dim_num) = average(1:dim_num) + x(1:dim_num)
  end do

  average(1:dim_num) = average(1:dim_num) / real ( sample_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now average the points, which should get a value'
  write ( *, '(a)' ) '  close to zero, and closer as sample_num increases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,3f8.4)' ) '  Average:        ', average(1:dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now choose a random direction, sample the same'
  write ( *, '(a)' ) '  number of points, and compute the dot product with'
  write ( *, '(a)' ) '  the direction.'
  write ( *, '(a)' ) '  Take the absolute value of each dot product '
  write ( *, '(a)' ) '  and sum and average.'

  do j = 1, 5

    call sphere_unit_sample_3d_2 ( seed, v )

    dot_average = 0.0D+00

    do i = 1, sample_num
      call sphere_unit_sample_3d_2 ( seed, x )
      dot_average = dot_average &
        + abs ( dot_product ( x(1:dim_num), v(1:dim_num) ) )
    end do

    dot_average = dot_average / real ( sample_num, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,3f8.4)' ) '  V:                ', v(1:dim_num)
    write ( *, '(a, f8.4)' ) '  Average |(XdotV)| ', dot_average

  end do

  return
end
subroutine test193 ( )

!*****************************************************************************80
!
!! TEST193 tests SPHERE_UNIT_SAMPLE_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) average(dim_num)
  real ( kind = 8 ) dot_average
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: sample_num = 1000
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) v(dim_num)
  real ( kind = 8 ) x(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST193'
  write ( *, '(a)' ) '  For the unit sphere in N dimensions:'
  write ( *, '(a)' ) '  SPHERE_UNIT_SAMPLE_ND samples;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A few sample values:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call sphere_unit_sample_nd ( dim_num, seed, x )
    write ( *, '(2x,3f8.4)' ) x(1:dim_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of sample points = ', sample_num

  average(1:dim_num) = 0.0D+00

  do i = 1, sample_num
    call sphere_unit_sample_nd ( dim_num, seed, x )
    average(1:dim_num) = average(1:dim_num) + x(1:dim_num)
  end do

  average(1:dim_num) = average(1:dim_num) / real ( sample_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now average the points, which should get a value'
  write ( *, '(a)' ) '  close to zero, and closer as N increases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,3f8.4)' ) '  Average:        ', average(1:dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now choose a random direction, sample the same'
  write ( *, '(a)' ) '  number of points, and compute the dot product with'
  write ( *, '(a)' ) '  the direction.'
  write ( *, '(a)' ) '  Take the absolute value of each dot product '
  write ( *, '(a)' ) '  and sum and average.'

  do j = 1, 5

    call sphere_unit_sample_nd ( dim_num, seed, v )

    dot_average = 0.0D+00

    do i = 1, sample_num
      call sphere_unit_sample_nd ( dim_num, seed, x )
      dot_average = dot_average &
        + abs ( dot_product ( x(1:dim_num), v(1:dim_num) ) )
    end do

    dot_average = dot_average / real ( sample_num, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,3f8.4)' ) '  V:                ', v(1:dim_num)
    write ( *, '(a, f8.4)' ) '  Average |(XdotV)| ', dot_average

  end do

  return
end
subroutine test194 ( )

!*****************************************************************************80
!
!! TEST194 tests SPHERE_UNIT_SAMPLE_ND_2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) average(dim_num)
  real ( kind = 8 ) dot_average
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: sample_num = 1000
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) v(dim_num)
  real ( kind = 8 ) x(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST194'
  write ( *, '(a)' ) '  For the unit sphere in N dimensions:'
  write ( *, '(a)' ) '  SPHERE_UNIT_SAMPLE_ND_2 samples;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A few sample values:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call sphere_unit_sample_nd_2 ( dim_num, seed, x )
    write ( *, '(2x,3f8.4)' ) x(1:dim_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of sample points = ', sample_num

  average(1:dim_num) = 0.0D+00

  do i = 1, sample_num
    call sphere_unit_sample_nd_2 ( dim_num, seed, x )
    average(1:dim_num) = average(1:dim_num) + x(1:dim_num)
  end do

  average(1:dim_num) = average(1:dim_num) / real ( sample_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now average the points, which should get a value'
  write ( *, '(a)' ) '  close to zero, and closer as sample_num increases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,3f8.4)' ) '  Average:        ', average(1:dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now choose a random direction, sample the same'
  write ( *, '(a)' ) '  number of points, and compute the dot product with'
  write ( *, '(a)' ) '  the direction.'
  write ( *, '(a)' ) '  Take the absolute value of each dot product '
  write ( *, '(a)' ) '  and sum and average.'

  do j = 1, 5

    call sphere_unit_sample_nd_2 ( dim_num, seed, v )

    dot_average = 0.0D+00

    do i = 1, sample_num
      call sphere_unit_sample_nd_2 ( dim_num, seed, x )
      dot_average = dot_average &
        + abs ( dot_product ( x(1:dim_num), v(1:dim_num) ) )
    end do

    dot_average = dot_average / real ( sample_num, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,3f8.4)' ) '  V:                ', v(1:dim_num)
    write ( *, '(a, f8.4)' ) '  Average |(XdotV)| ', dot_average

  end do

  return
end
subroutine test195 ( )

!*****************************************************************************80
!
!! TEST195 tests SPHERE_UNIT_SAMPLE_ND_3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) average(dim_num)
  real ( kind = 8 ) dot_average
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: sample_num = 1000
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) v(dim_num)
  real ( kind = 8 ) x(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST195'
  write ( *, '(a)' ) '  For the unit sphere in N dimensions:'
  write ( *, '(a)' ) '  SPHERE_UNIT_SAMPLE_ND_3 samples;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A few sample values:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call sphere_unit_sample_nd_3 ( dim_num, seed, x )
    write ( *, '(2x,3f8.4)' ) x(1:dim_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension = ', dim_num
  write ( *, '(a,i8)' ) '  Number of sample points = ', sample_num

  average(1:dim_num) = 0.0D+00

  do i = 1, sample_num
    call sphere_unit_sample_nd_3 ( dim_num, seed, x )
    average(1:dim_num) = average(1:dim_num) + x(1:dim_num)
  end do

  average(1:dim_num) = average(1:dim_num) / real ( sample_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now average the points, which should get a value'
  write ( *, '(a)' ) '  close to zero, and closer as sample_num increases.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,3f8.4)' ) '  Average:        ', average(1:dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now choose a random direction, sample the same'
  write ( *, '(a)' ) '  number of points, and compute the dot product with'
  write ( *, '(a)' ) '  the direction.'
  write ( *, '(a)' ) '  Take the absolute value of each dot product '
  write ( *, '(a)' ) '  and sum and average.'

  do j = 1, 5

    call sphere_unit_sample_nd_3 ( dim_num, seed, v )

    dot_average = 0.0D+00

    do i = 1, sample_num
      call sphere_unit_sample_nd_3 ( dim_num, seed, x )
      dot_average = dot_average &
        + abs ( dot_product ( x(1:dim_num), v(1:dim_num) ) )
    end do

    dot_average = dot_average / real ( sample_num, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,3f8.4)' ) '  V:                ', v(1:dim_num)
    write ( *, '(a, f8.4)' ) '  Average |(XdotV)| ', dot_average

  end do

  return
end
subroutine test1955 ( )

!*****************************************************************************80
!
!! TEST1955 tests SPHERE_UNIT_VOLUME_ND and SPHERE_UNIT_VOLUME_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
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
  write ( *, '(a)' ) 'TEST1955:'
  write ( *, '(a)' ) '  SPHERE_UNIT_VOLUME_ND evaluates the area of the unit'
  write ( *, '(a)' ) '  sphere in N dimensions.'
  write ( *, '(a)' ) '  SPHERE_UNIT_VOLUME_VALUES returns some test values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     DIM_NUM    Exact          Computed'
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
subroutine test196 ( )

!*****************************************************************************80
!
!! TEST196 tests SHAPE_POINT_DIST_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: side_num = 4
  integer ( kind = 4 ), parameter :: test_num = 9

  real ( kind = 8 ) dist
  real ( kind = 8 ) p(1:dim_num)
  real ( kind = 8 ),dimension ( dim_num ) :: p1 = (/ 5.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: pc = (/ 3.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension (dim_num,test_num) :: p_test = reshape ( (/ &
     3.0D+00,  0.0D+00, &
     5.0D+00,  0.0D+00, &
     4.0D+00,  0.0D+00, &
    10.0D+00,  0.0D+00, &
     8.0D+00,  5.0D+00, &
     6.0D+00,  6.0D+00, &
     1.0D+00,  2.0D+00, &
     2.5D+00, -0.5D+00, &
     4.0D+00, -1.0D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST196'
  write ( *, '(a)' ) '  For a shape in 2D,'
  write ( *, '(a)' ) '  SHAPE_POINT_DIST_2D computes the distance'
  write ( *, '(a)' ) '  to a point;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Number of sides:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i8)' ) side_num

  call r8vec_print ( dim_num, pc, '  Center of square:' )

  call r8vec_print ( dim_num, p1, '  Square vertex #1' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TEST       X            Y            DIST'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call shape_point_dist_2d ( pc, p1, side_num, p, dist ) 

    write ( *, '(2x,i8,3g14.6)' ) test, p(1:dim_num), dist

  end do
 
  return
end
subroutine test197 ( )

!*****************************************************************************80
!
!! TEST197 tests SHAPE_POINT_DIST_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: side_num = 6
  integer ( kind = 4 ), parameter :: test_num = 8

  real ( kind = 8 ) dist
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: p1 = (/ 5.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension (dim_num) :: pc = (/ 3.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension (dim_num,test_num) :: p_test = reshape ( (/ &
     3.0D+00, 0.0D+00, &
     5.0D+00, 0.0D+00, &
     4.0D+00, 0.0D+00, &
    10.0D+00, 0.0D+00, &
     4.0D+00, 1.7320508D+00, &
     5.0D+00, 3.4641016D+00,&
     3.0D+00, 1.7320508D+00, &
     3.0D+00, 0.86602539D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST197'
  write ( *, '(a)' ) '  For a shape in 2D,'
  write ( *, '(a)' ) '  SHAPE_POINT_DIST_2D computes the distance'
  write ( *, '(a)' ) '  to a point;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Number of sides:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i8)' ) side_num

  call r8vec_print ( dim_num, pc, '  Center of hexagon:' )

  call r8vec_print ( dim_num, p1, '  Hexagon vertex #1' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TEST         X		Y	     DIST'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call shape_point_dist_2d ( pc, p1, side_num, p, dist ) 

    write ( *, '(2x,i8,3g14.6)' ) test, p(1:dim_num), dist

  end do
 
  return
end
subroutine test198 ( )

!*****************************************************************************80
!
!! TEST198 tests SHAPE_POINT_NEAR_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: side_num = 6
  integer ( kind = 4 ), parameter :: test_num = 8

  real ( kind = 8 ) dist
  real ( kind = 8 ) p(1:dim_num)
  real ( kind = 8 ), dimension(1:dim_num,test_num) :: p_test = &
    reshape ( (/ &
     3.0D+00, 0.0D+00, &
     5.0D+00, 0.0D+00, &
     4.0D+00, 0.0D+00, &
    10.0D+00, 0.0D+00, &
     4.0D+00, 1.7320508D+00, &
     5.0D+00, 3.4641016D+00, &
     3.0D+00, 1.7320508D+00, &
     3.0D+00, 0.86602539D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num ) :: p1 = (/ &
    5.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num ) :: pc = (/ &
    3.0D+00, 0.0D+00 /)
  real ( kind = 8 ) pn(1:dim_num)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST198'
  write ( *, '(a)' ) '  For a shape in 2D,'
  write ( *, '(a)' ) '  SHAPE_POINT_NEAR_2D computes the nearest'
  write ( *, '(a)' ) '  point to a point;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Number of sides:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i8)' ) side_num

  call r8vec_print ( dim_num, pc, '  Hexagon center:' )

  call r8vec_print ( dim_num, p1, '  Hexagon vertex #1' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TEST       X            Y            ' // &
    '  PN     Dist'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call shape_point_near_2d ( pc, p1, side_num, p, pn, dist ) 

    write ( *, '(2x,i8,5f12.4)' ) test, p(1:dim_num), pn(1:dim_num), dist

  end do
 
  return
end
subroutine test199 ( )

!*****************************************************************************80
!
!! TEST199 tests SHAPE_RAY_INT_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: side_num = 6
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ), dimension ( dim_num ) :: p1 = (/ 5.0D+00, 0.0D+00 /)
  real ( kind = 8 ) pa(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: pa_test = reshape ( (/ &
    3.0D+00,  0.0D+00, &
    3.0D+00,  0.0D+00, &
    3.0D+00, -1.0D+00, &
    3.0D+00, -1.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) pb(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: pb_test = reshape ( (/ &
    4.0D+00,  0.0D+00, &
    3.0D+00,  1.0D+00, &
    3.0D+00,  1.0D+00, &
    7.0D+00,  5.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num ) :: pc = (/ 3.0D+00, 0.0D+00 /)
  real ( kind = 8 ) pint(dim_num)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST199'
  write ( *, '(a)' ) '  For a shape in 2D,'
  write ( *, '(a)' ) '  SHAPE_RAY_INT_2D computes the intersection of'
  write ( *, '(a)' ) '  a shape and a ray whose origin is within'
  write ( *, '(a)' ) '  the shape.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Number of sides:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i8)' ) side_num

  call r8vec_print ( dim_num, pc, '  Hexagon center:' )

  call r8vec_print ( dim_num, p1, '  Hexagon vertex #1' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TEST       XA          YA          XB' // &
    '          YB          XI          YI'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    pa(1:dim_num) = pa_test(1:dim_num,test)
    pb(1:dim_num) = pb_test(1:dim_num,test)

    call shape_ray_int_2d ( pc, p1, side_num, pa, pb, pint ) 

    write ( *, '(2x,i8,6f12.4)' ) &
      test, pa(1:dim_num), pb(1:dim_num), pint(1:dim_num)

  end do
 
  return
end
subroutine test200 ( )

!*****************************************************************************80
!
!! TEST200 tests SPHERE_TRIANGLE_SIDES_TO_ANGLES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) as
  real ( kind = 8 ) b
  real ( kind = 8 ) bs
  real ( kind = 8 ) c
  real ( kind = 8 ) cs
  real ( kind = 8 ) degrees_to_radians
  real ( kind = 8 ), parameter :: r = 10.0D+00
  real ( kind = 8 ) radians_to_degrees

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST200'
  write ( *, '(a)' ) '  SPHERE_TRIANGLE_SIDES_TO_ANGLES takes the sides of a'
  write ( *, '(a)' ) '  spherical triangle and determines the angles.'

  as = 121.0D+00 + ( 15.4D+00 / 60.0D+00 )
  bs = 104.0D+00 + ( 54.7D+00 / 60.0D+00 )
  cs =  65.0D+00 + ( 42.5D+00 / 60.0D+00 )

  as = degrees_to_radians ( as )
  bs = degrees_to_radians ( bs )
  cs = degrees_to_radians ( cs )

  as = r * as
  bs = r * bs 
  cs = r * cs
!
!  Get the spherical angles.
!
  call sphere_triangle_sides_to_angles ( r, as, bs, cs, a, b, c )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  A       = ', a,  ' (radians)'
  a = radians_to_degrees ( a )
  write ( *, '(a,f8.4,a)' ) '          = ', a,  ' ( degrees )'
  a = 117.0D+00 + ( 58.0D+00 / 60.0D+00 )
  write ( *, '(a,f8.4,a)' ) '  Correct = ', a, ' (degrees)'

  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  B       = ', b,  ' (radians)'
  b = radians_to_degrees ( b )
  write ( *, '(a,f8.4,a)' ) '          = ', b,  ' ( degrees )'
  b = 93.0D+00 + ( 13.8D+00 / 60.0D+00 )
  write ( *, '(a,f8.4,a)' ) '  Correct = ', b, ' (degrees)'

  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.4,a)' ) '  C       = ', c,  ' (radians)'
  c = radians_to_degrees ( c )
  write ( *, '(a,f8.4,a)' ) '          = ', c,  ' ( degrees )'
  c = 70.0D+00 + ( 20.6D+00 / 60.0D+00 )
  write ( *, '(a,f8.4,a)' ) '  Correct = ', c, ' (degrees)'

  return
end
subroutine test201 ( )

!*****************************************************************************80
!
!! TEST201 tests STRING_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: vec_num = 15

  integer ( kind = 4 ) i
  integer ( kind = 4 ) jstrng
  integer ( kind = 4 ) order(vec_num)
  real ( kind = 8 ), dimension ( dim_num, vec_num ) :: p1 = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    3.0D+00, 4.0D+00, &
    2.0D+00, 2.0D+00, &
    3.0D+00, 2.0D+00, &
    2.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, &
    0.0D+00, 5.0D+00, &
    1.0D+00, 2.0D+00, &
    3.0D+00, 2.0D+00, &
    0.0D+00, 0.0D+00, &
    5.0D+00, 5.0D+00, &
    3.0D+00, 3.0D+00, &
    2.0D+00, 4.0D+00, &
    7.0D+00, 4.0D+00, &
    1.0D+00, 0.0D+00 /), (/ dim_num, vec_num /) )
  real ( kind = 8 ), dimension ( dim_num, vec_num ) :: p2 = reshape ( (/ &
    1.0D+00, 1.0D+00, &
    2.0D+00, 4.0D+00, &
    1.0D+00, 3.0D+00, &
    2.0D+00, 3.0D+00, &
    2.0D+00, 2.0D+00, &
    1.0D+00, 2.0D+00, &
    1.0D+00, 6.0D+00, &
    1.0D+00, 3.0D+00, &
    3.0D+00, 3.0D+00, &
    1.0D+00, 0.0D+00, &
    6.0D+00, 6.0D+00, &
    3.0D+00, 4.0D+00, &
    2.0D+00, 3.0D+00, &
    5.0D+00, 5.0D+00, &
    2.0D+00, 1.0D+00 /), (/ dim_num, vec_num /) )
  integer ( kind = 4 ) string(vec_num)
  integer ( kind = 4 ) string_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST201'
  write ( *, '(a)' ) '  STRING_2D takes a set of line segments, and'
  write ( *, '(a)' ) '  "strings" them together.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I     P1      P2'
  write ( *, '(a)' ) ' '
  do i = 1, vec_num
    write ( *, '(2x,i8,4g14.6)' ) i, p1(1:2,i), p2(1:2,i)
  end do
 
  call string_2d ( vec_num, p1, p2, string_num, order, string )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Found ', string_num, ' groups of segments.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  STRING  ORDER      P1      P2'
  write ( *, '(a)' ) ' '

  jstrng = 1

  do i = 1, vec_num
 
    if ( jstrng < string(i) ) then
      write ( *, '(a)' ) ' '
      jstrng = jstrng + 1
    end if
 
    write ( *, '(2x,i3,1x,i3,4f10.4)' ) string(i), order(i), p1(1:2,i), &
      p2(1:2,i)

  end do
 
  return
end
subroutine test202 ( )

!*****************************************************************************80
!
!! TEST202 tests SUPER_ELLIPSE_POINTS_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 24
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ), dimension ( dim_num ) :: pc = (/ &
    5.0D+00, -2.0D+00 /)
  real ( kind = 8 ) expo
  real ( kind = 8 ) p(dim_num,n)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) psi
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2

  r1 = 3.0D+00
  r2 = 1.0D+00
  expo = 1.5D+00
  psi = pi / 6.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST202'
  write ( *, '(a)' ) &
    '  SUPER_ELLIPSE_POINTS_2D returns points on a super ellipse;'

  call r8vec_print ( dim_num, pc, '  Superellipse center:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  radii R1 = ', r1, ' R2 = ', r2
  write ( *, '(a,g14.6)' ) '  exponent EXPO = ', expo
  write ( *, '(a,g14.6)' ) '  and angle PSI = ', psi

  call super_ellipse_points_2d ( pc, r1, r2, expo, psi, n, p )

  call r8mat_transpose_print ( dim_num, n, p, '  Sample points:' )
 
  return
end
subroutine test203 ( )

!*****************************************************************************80
!
!! TEST203 tests TETRAHEDRON_CENTROID_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) centroid(dim_num)
  real ( kind = 8 ), dimension (dim_num,4) :: tetra = reshape ( (/&
     0.000000D+00,  0.942809D+00, -0.333333D+00, &
    -0.816496D+00, -0.816496D+00, -0.333333D+00, &
     0.816496D+00, -0.816496D+00, -0.333333D+00, &
     0.000000D+00,  0.000000D+00,  1.000000D+00 /), (/ dim_num, 4 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST203'
  write ( *, '(a)' ) '  For a tetrahedron in 3D,'
  write ( *, '(a)' ) '  TETRAHEDRON_CENTROID_3D computes the centroid;'

  call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

  call tetrahedron_centroid_3d ( tetra, centroid )

  call r8vec_print ( dim_num, centroid, '  Centroid:' )

  return
end
subroutine test2031 ( )

!*****************************************************************************80
!
!! TEST2031 tests TETRAHEDRON_CONTAINS_POINT_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) c(4)
  real ( kind = 8 ), dimension(4,test_num) :: c_test = reshape ( (/ &
     0.0D+00, 0.1D+00,  0.2D+00, 0.7D+00, &
    -1.3D+00, 2.0D+00,  0.2D+00, 0.1D+00, &
     0.8D+00, 0.6D+00, -0.5D+00, 0.1D+00 /), (/ 4, test_num /) )
  logical inside
  real ( kind = 8 ) p(dim_num)
  integer ( kind = 4 ) test
  real ( kind = 8 ), dimension(dim_num,4) :: tetra = reshape ( (/&
     0.000000D+00,  0.942809D+00, -0.333333D+00, &
    -0.816496D+00, -0.816496D+00, -0.333333D+00, &
     0.816496D+00, -0.816496D+00, -0.333333D+00, &
     0.000000D+00,  0.000000D+00,  1.000000D+00 /), (/ dim_num, 4 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2031'
  write ( *, '(a)' ) '  For a tetrahedron in 3D,'
  write ( *, '(a)' ) '  TETRAHEDRON_CONTAINS_POINT_3D finds if a point '
  write ( *, '(a)' ) '  is inside;'

  call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  P     Inside_Tetra?'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    c(1:4) = c_test(1:4,test)

    p(1:dim_num) = matmul ( tetra(1:dim_num,1:4), c(1:4) )

    call tetrahedron_contains_point_3d ( tetra, p, inside )

    write ( *, '(2x,3g14.6,2x,l1)' ) p(1:dim_num), inside

  end do

  return
end
subroutine test2032 ( )

!*****************************************************************************80
!
!! TEST2032 tests TETRAHEDRON_CIRCUMSPHERE_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ), dimension(dim_num,4) :: tetra = reshape ( (/&
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.816496580927726D+00 /), &
    (/ dim_num, 4 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2032'
  write ( *, '(a)' ) '  For a tetrahedron in 3D,'
  write ( *, '(a)' ) '  TETRAHEDRON_CIRCUMSPHERE_3D computes the circumsphere;'

  call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

  call tetrahedron_circumsphere_3d ( tetra, r, pc )

  call r8vec_print ( dim_num, pc, '  Circumsphere center:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Circumsphere radius is ', r
 
  return
end
subroutine test20321 ( )

!*****************************************************************************80
!
!! TEST20321 tests TETRAHEDRON_EDGE_LENGTH_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) edge_length(6)
  real ( kind = 8 ), dimension(dim_num,4) :: tetra = reshape ( (/&
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.816496580927726D+00 /), &
     (/ dim_num, 4 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20321'
  write ( *, '(a)' ) '  For a tetrahedron in 3D,'
  write ( *, '(a)' ) '  TETRAHEDRON_EDGE_LENGTH_3D computes the edge lengths;'

  call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

  call tetrahedron_edge_length_3d ( tetra, edge_length )

  call r8vec_print ( 6, edge_length, '  Edge lengths:' )
 
  return
end
subroutine test20322 ( )

!*****************************************************************************80
!
!! TEST20322 tests TETRAHEDRON_INSPHERE_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ), dimension(dim_num,4) :: tetra = reshape ( (/&
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.816496580927726D+00 /), &
        (/ dim_num, 4 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20322'
  write ( *, '(a)' ) '  For a tetrahedron in 3D,'
  write ( *, '(a)' ) '  TETRAHEDRON_INSPHERE_3D computes the insphere;'

  call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

  call tetrahedron_insphere_3d ( tetra, r, pc )

  call r8vec_print ( dim_num, pc, '  Insphere center:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Insphere radius is ', r
 
  return
end
subroutine test203224 ( )

!*****************************************************************************80
!
!! TEST203224 tests TETRAHEDRON_LATTICE_LAYER_POINT_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ) c(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) layer
  logical              more
  integer ( kind = 4 ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST203224'
  write ( *, '(a)' ) '  TETRAHEDRON_LATTICE_LAYER_POINT_NEXT returns the next'
  write ( *, '(a)' ) '  point in a tetrahedron lattice layer defined by:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    C(4) - 1 < X(1)/C(1) + X(2)/C(2) +X(3)/C(3) <= C(4).'

  c(1) = 2
  c(2) = 3
  c(3) = 4
  v(1:n) = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  N = ', n
  write ( *, '(a)', ADVANCE = 'NO' ) '  C =       '
  do i = 1, n
    write ( *, '(2x,i4)', ADVANCE = 'NO' ) c(i)
  end do
  write ( *, '(a)', ADVANCE = 'YES' ) 

  do layer = 0, 2

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Layer ', layer
    write ( *, '(a)' ) ' '

    c(4) = layer
    more = .false.
    i = 0

    do
      call tetrahedron_lattice_layer_point_next ( c, v, more )
      if ( .not. more ) then
        write ( *, '(a)' ) '  No more.'
        exit
      end if
      i = i + 1
      write ( *, '(2x,i4,6x,10(2x,i4))' ) i, v(1:n)

    end do

  end do

  return
end
subroutine test203225 ( )

!*****************************************************************************80
!
!! TEST203225 tests TETRAHEDRON_LATTICE_POINT_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ) c(n+1)
  integer ( kind = 4 ) i
  logical              more
  integer ( kind = 4 ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST203225'
  write ( *, '(a)' ) '  TETRAHEDRON_LATTICE_POINT_NEXT returns the next lattice'
  write ( *, '(a)' ) '  point in a tetrahedron defined by:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    0 <= X(1)/C(1) + X(2)/C(2) + X(3)/C(3) <= C(4).'

  do i = 1, n + 1
    c(i) = n + 2 - i
  end do
  v(1:n) = 0
  more = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  N = ', n
  write ( *, '(a)', ADVANCE = 'NO' ) '  C =       '
  do i = 1, n + 1
    write ( *, '(2x,i4)', ADVANCE = 'NO' ) c(i)
  end do
  write ( *, '(a)', ADVANCE = 'YES' ) 
  write ( *, '(a)' ) ' '

  i = 0

  do
    call tetrahedron_lattice_point_next ( c, v, more )
    if ( .not. more ) then
      write ( *, '(a)' ) '  No more.'
      exit
    end if
    i = i + 1
    write ( *, '(2x,i4,6x,10(2x,i4))' ) i, v(1:n)
  end do

  return
end
subroutine test20323 ( )

!*****************************************************************************80
!
!! TEST20323 tests TETRAHEDRON_QUALITY1_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) quality
  real ( kind = 8 ), dimension(dim_num,4) :: tetra
  real ( kind = 8 ), dimension(dim_num,4,test_num) :: tetra_test = reshape ( (/&
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.816496580927726D+00, &
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.408248290463863D+00 /), &
        (/ dim_num, 4, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20323'
  write ( *, '(a)' ) '  For a tetrahedron in 3D,'
  write ( *, '(a)' ) '  TETRAHEDRON_QUALITY1_3D computes quality measure #1;'

  do test = 1, test_num

    tetra(1:dim_num,1:4) = tetra_test(1:dim_num,1:4,test)

    call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

    call tetrahedron_quality1_3d ( tetra, quality )

    write ( *, '(a,g14.6)' ) '  Tetrahedron quality is ', quality

  end do
 
  return
end
subroutine test203232 ( )

!*****************************************************************************80
!
!! TEST203232 tests TETRAHEDRON_QUALITY2_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) quality2
  real ( kind = 8 ), dimension(dim_num,4) :: tetra
  real ( kind = 8 ), dimension(dim_num,4,test_num) :: tetra_test = reshape ( (/&
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.816496580927726D+00, &
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.408248290463863D+00 /), &
        (/ dim_num, 4, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST203232'
  write ( *, '(a)' ) '  For a tetrahedron in 3D,'
  write ( *, '(a)' ) '  TETRAHEDRON_QUALITY_3D computes quality measure #2;'

  do test = 1, test_num

    tetra(1:dim_num,1:4) = tetra_test(1:dim_num,1:4,test)

    call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

    call tetrahedron_quality2_3d ( tetra, quality2 )

    write ( *, '(a,g14.6)' ) '  Tetrahedron quality is ', quality2

  end do

  return
end
subroutine test203233 ( )

!*****************************************************************************80
!
!! TEST203233 tests TETRAHEDRON_QUALITY3_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) quality3
  real ( kind = 8 ), dimension(dim_num,4) :: tetra
  real ( kind = 8 ), dimension(dim_num,4,test_num) :: tetra_test = reshape ( (/&
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.816496580927726D+00, &
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.408248290463863D+00 /), &
        (/ dim_num, 4, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST203233'
  write ( *, '(a)' ) '  For a tetrahedron in 3D,'
  write ( *, '(a)' ) '  TETRAHEDRON_QUALITY3_3D computes quality measure #3;'

  do test = 1, test_num

    tetra(1:dim_num,1:4) = tetra_test(1:dim_num,1:4,test)

    call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

    call tetrahedron_quality3_3d ( tetra, quality3 )

    write ( *, '(a,g14.6)' ) '  Tetrahedron quality is ', quality3

  end do

  return
end
subroutine test203234 ( )

!*****************************************************************************80
!
!! TEST203234 tests TETRAHEDRON_QUALITY4_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 2

  real ( kind = 8 ) quality4
  real ( kind = 8 ), dimension(dim_num,4) :: tetra
  real ( kind = 8 ), dimension(dim_num,4,test_num) :: tetra_test = reshape ( (/&
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.816496580927726D+00, &
     0.577350269189626D+00,  0.0D+00, 0.0D+00, &
    -0.288675134594813D+00,  0.5D+00, 0.0D+00, &
    -0.288675134594813D+00, -0.5D+00, 0.0D+00, &
     0.0D+00,                0.0D+00, 0.408248290463863D+00 /), &
        (/ dim_num, 4, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST203234'
  write ( *, '(a)' ) '  For a tetrahedron in 3D,'
  write ( *, '(a)' ) '  TETRAHEDRON_QUALITY4_3D computes quality measure #4;'

  do test = 1, test_num

    tetra(1:dim_num,1:4) = tetra_test(1:dim_num,1:4,test)

    call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices:' )

    call tetrahedron_quality4_3d ( tetra, quality4 )

    write ( *, '(a,g14.6)' ) '  Tetrahedron quality is ', quality4

  end do

  return
end
subroutine test203235 ( )

!*****************************************************************************80
!
!! TEST203235 tests TETRAHEDRON_RHOMBIC_SIZE_3D and TETRAHEDRON_RHOMBIC_SHAPE_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ), allocatable, dimension ( : ) ::  face_order
  integer ( kind = 4 ) face_order_max
  integer ( kind = 4 ), allocatable, dimension ( :, : ) ::  face_point
  integer ( kind = 4 ) point_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_coord

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST203235'
  write ( *, '(a)' ) '  For the cube,'
  write ( *, '(a)' ) '  TETRAHEDRON_RHOMBIC_SIZE_3D returns dimension information;'
  write ( *, '(a)' ) '  TETRAHEDRON_RHOMBIC_SHAPE_3D returns face and order information.'
  write ( *, '(a)' ) '  SHAPE_PRINT_3D prints this information.'
!
!  Get the data sizes.
!
  call tetrahedron_rhombic_size_3d ( point_num, edge_num, face_num, &
    face_order_max )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of vertices = ', point_num
  write ( *, '(a,i8)' ) '  Number of edges =    ', edge_num
  write ( *, '(a,i8)' ) '  Number of faces =    ', face_num
  write ( *, '(a,i8)' ) '  Maximum face order = ', face_order_max
!
!  Make room for the data.
!
  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )
  allocate ( point_coord(1:3,1:point_num) )
!
!  Get the data.
!
  call tetrahedron_rhombic_shape_3d ( point_num, face_num, face_order_max, &
    point_coord, face_order, face_point )
!
!  Print the data.
!
  call shape_print_3d ( point_num, face_num, face_order_max, &
    point_coord, face_order, face_point )

  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )

  return
end
subroutine test20324 ( )

!*****************************************************************************80
!
!! TEST20324 tests TETRAHEDRON_SAMPLE_3D and TETRAHEDRON_BARYCENTRIC_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 10

  real ( kind = 8 ) p(dim_num,test_num)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ), dimension(dim_num,4) :: t = reshape ( (/ &
     1.0D+00, 4.0D+00, 3.0D+00, &
     2.0D+00, 4.0D+00, 3.0D+00, &
     1.0D+00, 6.0D+00, 3.0D+00, &
     1.0D+00, 4.0D+00, 4.0D+00 /), (/ dim_num, 4 /) )
  integer ( kind = 4 ) test
  real ( kind = 8 ) xsi(dim_num+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20324'
  write ( *, '(a)' ) '  TETRAHEDRON_SAMPLE_3D samples a tetrahedron.'
  write ( *, '(a)' ) '  TETRAHEDRON_BARYCENTRIC_3D converts Cartesian to'
  write ( *, '(a)' ) '  barycentric coordinates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We are computing the barycentric coordinates just to'
  write ( *, '(a)' ) '  verify that the points are inside the tetrahedron.'

  call r8mat_transpose_print ( dim_num, 4, t, '  Tetrahedron vertices' )

  call tetrahedron_sample_3d ( t, test_num, seed, p )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      P                           Barycentric:'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    call tetrahedron_barycentric_3d ( t, p(1:3,test), xsi )
    write ( *, '(2x,3f8.4,4x,4f8.4)' ) p(1:dim_num,test), xsi(1:dim_num+1)
  end do

  return
end
subroutine test20325 ( )

!*****************************************************************************80
!
!! TEST20325 tests TETRAHEDRON_SIZE_3D and TETRAHEDRON_SHAPE_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) area
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) face
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ), allocatable, dimension ( : ) ::  face_order
  integer ( kind = 4 ) face_order_max
  integer ( kind = 4 ), allocatable, dimension ( :, : ) ::  face_point
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) normal(dim_num)
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point_coord
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: v
  real ( kind = 8 ) vave(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20325'
  write ( *, '(a)' ) '  For the tetrahedron,'
  write ( *, '(a)' ) '  TETRAHEDRON_SIZE_3D returns dimension information;'
  write ( *, '(a)' ) '  TETRAHEDRON_SHAPE_3D returns face and order info.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We will use this information to compute the'
  write ( *, '(a)' ) '  areas and centers of each face.'

  call tetrahedron_size_3d ( point_num, edge_num, face_num, face_order_max )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of vertices = ', point_num
  write ( *, '(a,i8)' ) '  Number of edges =    ', edge_num
  write ( *, '(a,i8)' ) '  Number of faces =    ', face_num
  write ( *, '(a,i8)' ) '  Maximum face order = ', face_order_max

  allocate ( face_order(1:face_num) )
  allocate ( face_point(1:face_order_max,1:face_num) )
  allocate ( point_coord(1:dim_num,1:point_num) )
  allocate ( v(1:dim_num,1:face_order_max) )

  call tetrahedron_shape_3d ( point_num, face_num, face_order_max, &
    point_coord, face_order, face_point )

  call shape_print_3d ( point_num, face_num, face_order_max, &
    point_coord, face_order, face_point )
!
!  Compute the area of each face.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Face  Order  Area'
  write ( *, '(a)' ) ' '

  do face = 1, face_num

    do j = 1, face_order(face)
      point = face_point(j,face)
      v(1:dim_num,j) = point_coord(1:dim_num,point)
    end do

    call polygon_area_3d ( face_order(face), v, area, normal )

    write ( *, '(2x,i8,i7,f8.4)' ) face, face_order(face), area

  end do
!
!  Find the center of each face.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Face  Center'
  write ( *, '(a)' ) ' '

  do i = 1, face_num

    vave(1:dim_num) = 0.0D+00

    do j = 1, face_order(i)
      k = face_point(j,i)
      vave(1:dim_num) = vave(1:dim_num) + point_coord(1:dim_num,k)
    end do

    vave(1:dim_num) = vave(1:dim_num) / real ( face_order(i), kind = 8 )

    write ( *, '(2x,i8,3f8.4)' ) i, vave(1:dim_num)

  end do

  deallocate ( face_order )
  deallocate ( face_point )
  deallocate ( point_coord )
  deallocate ( v )

  return
end
subroutine test2033 ( )

!*****************************************************************************80
!
!! TEST2033 tests TETRAHEDRON_VOLUME_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ), dimension(dim_num,4) :: tetra = reshape ( (/&
     0.000000D+00,  0.942809D+00, -0.333333D+00, &
    -0.816496D+00, -0.816496D+00, -0.333333D+00, &
     0.816496D+00, -0.816496D+00, -0.333333D+00, &
     0.000000D+00,  0.000000D+00,  1.000000D+00 /), (/ dim_num, 4 /) )
  real ( kind = 8 ) volume

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2033'
  write ( *, '(a)' ) '  For a tetrahedron in 3D,'
  write ( *, '(a)' ) '  TETRAHEDRON_VOLUME_3D computes the volume;'

  call r8mat_transpose_print ( dim_num, 4, tetra, '  Tetrahedron vertices' )

  call tetrahedron_volume_3d ( tetra, volume )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Volume = ', volume

  return
end
subroutine test204 ( )

!*****************************************************************************80
!
!! TEST204 tests TMAT_INIT, TMAT_ROT_AXIS, TMAT_ROT_VECTOR, TMAT_SCALE, TMAT_SHEAR, TMAT_TRANS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) a(4,4)
  real ( kind = 8 ) angle
  real ( kind = 8 ) axis(dim_num)
  character axis1
  character ( len = 2 ) axis2
  real ( kind = 8 ) b(4,4)
  real ( kind = 8 ) s
  real ( kind = 8 ) v(3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST204'
  write ( *, '(a)' ) '  TMAT geometric transformation matrix routines:'
  write ( *, '(a)' ) '  TMAT_INIT initializes,'
  write ( *, '(a)' ) '  TMAT_ROT_AXIS for rotation about an axis,'
  write ( *, '(a)' ) '  TMAT_ROT_VECTOR for rotation about a vector,'
  write ( *, '(a)' ) '  TMAT_SCALE for scaling,'
  write ( *, '(a)' ) '  TMAT_SHEAR for shear,'
  write ( *, '(a)' ) '  TMAT_TRANS for translation'
!
!  Initialization.
!
  call tmat_init ( a )

  call r8mat_print ( 4, 4, a, '  Initial transformation matrix:' )
!
!  Rotation about an axis.
!
  angle = 30.0D+00
  axis1 = 'x'
  call tmat_rot_axis ( a, angle, axis1, b )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Transformation matrix for'
  write ( *, '(a)' ) '  rotation about ' // axis1
  write ( *, '(a,g14.6)' ) '  by ' , angle

  call r8mat_print ( 4, 4, b, ' ' )
!
!  Rotation about a vector.
!
  angle = 30.0D+00
  axis(1:dim_num) = (/ 1.0D+00, 2.0D+00, 3.0D+00 /)
  call tmat_rot_vector ( a, angle, axis, b )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Transformation matrix for'
  write ( *, '(a,3g14.6)' ) '  rotation about ', axis(1:dim_num)
  write ( *, '(a,g14.6)' ) '  of ', angle

  call r8mat_print ( 4, 4, b, ' ' )
!
!  Scaling.
!
  v(1:3) = (/ 2.0D+00, 0.5D+00, 10.0D+00 /)
  call tmat_scale ( a, v, b )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Transformation matrix for'
  write ( *, '(a,3g14.6)' ) '  scaling by ', v(1:3)

  call r8mat_print ( 4, 4, b, ' ' )
!
!  Shear.
!
  axis2 = 'xy'
  s = 0.5D+00
  call tmat_shear ( a, axis2, s, b )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Transformation matrix for'
  write ( *, '(2x,a)' ) axis2
  write ( *, '(a,g14.6)' ) '  shear coefficient of ', s

  call r8mat_print ( 4, 4, b, ' ' )
!
!  Translation.
!
  v(1:3) = (/ 1.0D+00, 2.0D+00, 3.0D+00 /)
  call tmat_trans ( a, v, b )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Transformation matrix for'
  write ( *, '(a,3g14.6)' ) '  translation by ', v(1:3)

  call r8mat_print ( 4, 4, b, ' ' )

  return
end
subroutine test205 ( )

!*****************************************************************************80
!
!! TEST205 tests TMAT_MXP2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) a(4,4)
  real ( kind = 8 ) angle
  real ( kind = 8 ) axis(dim_num)
  character axis1
  character ( len = 2 ) axis2
  real ( kind = 8 ) b(4,4)
  real ( kind = 8 ), dimension ( dim_num, n ) :: point = reshape ( (/ &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00 /), (/ dim_num, n /) )
  real ( kind = 8 ) point2(dim_num,n)
  real ( kind = 8 ) s
  real ( kind = 8 ) v(1:3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST205'
  write ( *, '(a)' ) '  TMAT_MXP2 applies a geometric transformation'
  write ( *, '(a)' ) '  matrix to a set of points.'

  call r8mat_transpose_print ( 3, n, point, '  Points:' )
!
!  Initialization of transformation matrix.
!
  call tmat_init ( a )

  call r8mat_print ( 4, 4, a, '  Initial transformation matrix:' )
!
!  Rotation about an axis.
!
  angle = 30.0D+00
  axis1 = 'x'
  call tmat_rot_axis ( a, angle, axis1, b )

  call tmat_mxp2 ( b, n, point, point2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rotation about ' // axis1
  write ( *, '(a,g14.6)' ) '  by ' , angle

  call r8mat_transpose_print ( 3, n, point2, ' ' )
!
!  Rotation about a vector.
!
  angle = 30.0D+00
  axis(1:3) = (/ 1.0D+00, 2.0D+00, 3.0D+00 /)

  call tmat_rot_vector ( a, angle, axis, b )

  call tmat_mxp2 ( b, n, point, point2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,3g14.6)' ) '  Rotation about ', axis(1:3)
  write ( *, '(a,g14.6)' ) '  of ', angle

  call r8mat_transpose_print ( 3, n, point2, ' ' )
!
!  Scaling.
!
  v(1:3) = (/ 2.0D+00, 0.5D+00, 10.0D+00 /)
  call tmat_scale ( a, v, b )

  call tmat_mxp2 ( b, n, point, point2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,3g14.6)' ) '  Scaling by ', v(1:3)

  call r8mat_transpose_print ( 3, n, point2, ' ' )
!
!  Shear.
!
  axis2 = 'xy'
  s = 0.5D+00
  call tmat_shear ( a, axis2, s, b )

  call tmat_mxp2 ( b, n, point, point2 )

  write ( *, '(a)' ) ' '
  write ( *, '(2x,a)' ) axis2
  write ( *, '(a,g14.6)' ) ' shear coefficient of ', s

  call r8mat_transpose_print ( 3, n, point2, ' ' )
!
!  Translation.
!
  v(1:3) = (/ 1.0D+00, 2.0D+00, 3.0D+00 /)
  call tmat_trans ( a, v, b )

  call tmat_mxp2 ( b, n, point, point2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,3g14.6)' ) '  Translation by ', v(1:3)

  call r8mat_transpose_print ( 3, n, point2, ' ' )

  return
end
subroutine test206 ( )

!*****************************************************************************80
!
!! TEST206 tests TRIANGLE_ANGLES_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) angle(3)
  integer ( kind = 4 ) i
  real ( kind = 8 ) radians_to_degrees
  real ( kind = 8 ), dimension ( dim_num, 3 ) :: t = reshape ( (/ &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00 /), (/ dim_num, 3 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST206'
  write ( *, '(a)' ) '  For a triangle in 2D,'
  write ( *, '(a)' ) '  TRIANGLE_ANGLES_2D computes the angles;'

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

  call triangle_angles_2d ( t, angle )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Radians      Degrees'
  write ( *, '(a)' ) ' '
  do i = 1, 3
    write ( *, '(2x,g14.6,2x,g14.6)' ) angle(i), radians_to_degrees ( angle(i) )
  end do

  return
end
subroutine test20605 ( )

!*****************************************************************************80
!
!! TEST20605 tests TRIANGLE_ANGLES_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) angle(3)
  integer ( kind = 4 ) i
  real ( kind = 8 ) radians_to_degrees
  real ( kind = 8 ), dimension ( dim_num, 3 ) :: t = reshape ( (/ &
    1.0D+00,       2.0D+00,       3.0D+00, &
    2.4142137D+00, 3.4142137D+00, 3.0D+00, &
    1.7071068D+00, 2.7071068D+00, 4.0D+00 /), (/ dim_num, 3 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20605'
  write ( *, '(a)' ) '  For a triangle in 3D:'
  write ( *, '(a)' ) '  TRIANGLE_ANGLES_3D computes the angles;'
  write ( *, '(a)' ) ' '

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices: ' )

  call triangle_angles_3d ( t, angle )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Radians      Degrees'
  write ( *, '(a)' ) ' '

  do i = 1, 3
    write ( *, '(2x,g14.6,2x,g14.6)' ) angle(i), radians_to_degrees ( angle(i) )
  end do

  return
end
subroutine test2061 ( )

!*****************************************************************************80
!
!! TEST2061 tests TRIANGLE_AREA_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) area
  real ( kind = 8 ), dimension ( dim_num, 3 ) :: t = reshape ( (/ &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00 /), (/ dim_num, 3 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2061'
  write ( *, '(a)' ) '  For a triangle in 2D,'
  write ( *, '(a)' ) '  TRIANGLE_AREA_2D computes the area;'

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

  call triangle_area_2d ( t, area )

  write ( *, '(a,g14.6)' ) '  Triangle area is ', area

  return
end
subroutine test2062 ( )

!*****************************************************************************80
!
!! TEST2062 tests TRIANGLE_AREA_HERON;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) area
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  real ( kind = 8 ) s(3)
  real ( kind = 8 ), dimension ( dim_num, 3 ) :: t = reshape ( (/ &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00 /), (/ dim_num, 3 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2062'
  write ( *, '(a)' ) '  For a triangle in any dimension,'
  write ( *, '(a)' ) '  TRIANGLE_AREA_HERON computes the area;'

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

  do j = 1, 3

    s(j) = 0.0D+00

    jp1 = mod ( j, 3 ) + 1

    do i = 1, dim_num
      s(j) = s(j) + ( t(i,j) - t(i,jp1) )**2
    end do

    s(j) = sqrt ( s(j) )

  end do

  call r8vec_print ( 3, s, '  Side lengths:' )

  call triangle_area_heron ( s, area )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The area is ', area

  return
end
subroutine test209 ( )

!*****************************************************************************80
!
!! TEST209 tests TRIANGLE_AREA_3D, TRIANGLE_AREA_3D_2, TRIANGLE_AREA_3D_3;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) area
  real ( kind = 8 ), dimension ( dim_num, 3 ) :: t = reshape ( (/ &
    1.0D+00,       2.0D+00,       3.0D+00, &
    2.4142137D+00, 3.4142137D+00, 3.0D+00, &
    1.7071068D+00, 2.7071068D+00, 4.0D+00 /), (/ dim_num, 3 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST209'
  write ( *, '(a)' ) '  For a triangle in 3D:'
  write ( *, '(a)' ) '  TRIANGLE_AREA_3D   computes the area;'
  write ( *, '(a)' ) '  TRIANGLE_AREA_3D_2 computes the area;'
  write ( *, '(a)' ) '  TRIANGLE_AREA_3D_3 computes the area;'

  call r8mat_print ( dim_num, 3, t, '  Triangle (vertices are columns)' )

  call triangle_area_3d ( t, area )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Area #1 ', area

  call triangle_area_3d_2 ( t, area )

  write ( *, '(a,g14.6)' ) '  Area #2 ', area

  call triangle_area_3d_3 ( t, area )

  write ( *, '(a,g14.6)' ) '  Area #3 ', area

  return
end
subroutine test20655 ( )

!*****************************************************************************80
!
!! TEST20655 tests TRIANGLE_BARYCENTRIC_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 7

  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p_test = reshape ( (/ &
     0.25D+00,   0.25D+00, &
     0.75D+00,   0.25D+00, &
     1.00D+00,   1.00D+00, &
    11.00D+00,   0.50D+00, &
     0.00D+00,   1.00D+00, &
     0.50D+00, -10.00D+00, &
     0.60D+00,   0.60D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num, 3 ) :: t = reshape ( (/ &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00 /), (/ dim_num, 3 /) )
  integer ( kind = 4 ) test
  real ( kind = 8 ) xsi(dim_num+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20655'
  write ( *, '(a)' ) '  For a triangle in 2D,'
  write ( *, '(a)' ) '  TRIANGLE_BARYCENTRIC_2D converts XY coordinates'
  write ( *, '(a)' ) '  to barycentric XSI coordinates;'

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          P       XSI'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call triangle_barycentric_2d ( t, p, xsi )

    write ( *, '(2x,2f8.3,2x,3f8.3)' ) p(1:dim_num), xsi(1:dim_num+1)

  end do
 
  return
end
subroutine test2066 ( )

!*****************************************************************************80
!
!! TEST2066 tests TRIANGLE_CENTROID_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) centroid(dim_num)
  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ), dimension(dim_num,3,test_num) :: t_test = reshape ( (/ &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
         0.0D+00,  1.0D+00, &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
         0.5D+00,  0.86602539D+00, &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
         0.5D+00, 10.0D+00, &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
        10.0D+00,  2.0D+00 /), (/ dim_num, 3, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2066'
  write ( *, '(a)' ) '  For a triangle in 2D:'
  write ( *, '(a)' ) '  TRIANGLE_CENTROID_2D computes the centroid.'

  do test = 1, test_num

    t(1:dim_num,1:3) = t_test(1:dim_num,1:3,test)

    call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

    call triangle_centroid_2d ( t, centroid )

    call r8vec_print ( dim_num, centroid, '  Centroid:' )
 
  end do

  return
end
subroutine test2094 ( )

!*****************************************************************************80
!
!! TEST2094 tests TRIANGLE_CENTROID_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) centroid(dim_num)
  real ( kind = 8 ), dimension (dim_num,3) :: t = reshape ( (/ &
    1.0D+00,       2.0D+00,       3.0D+00, &
    2.4142137D+00, 3.4142137D+00, 3.0D+00, &
    1.7071068D+00, 2.7071068D+00, 4.0D+00 /), (/ dim_num, 3 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2094'
  write ( *, '(a)' ) '  For a triangle in 3D:'
  write ( *, '(a)' ) '  TRIANGLE_CENTROID_3D computes the centroid.'

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

  call triangle_centroid_3d ( t, centroid )

  call r8vec_print ( dim_num, centroid, '  Centroid:' )

  return
end
subroutine test2101 ( )

!*****************************************************************************80
!
!! TEST2101 tests TRIANGLE_CIRCUMCENTER_2D and others.
!
!  Discussion:
!
!    The functions tested include
!    * TRIANGLE_CIRCUMCENTER_2D;
!    * TRIANGLE_CIRCUMCENTER_2D_2;
!    * TRIANGLE_CIRCUMCENTER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) pc(m)
  real ( kind = 8 ) t(m,3)
  real ( kind = 8 ), dimension(m,3,test_num) :: t_test = reshape ( (/ &
         10.0D+00,  5.0D+00, &
         11.0D+00,  5.0D+00, &
         10.0D+00,  6.0D+00, &
         10.0D+00,  5.0D+00, &
         11.0D+00,  5.0D+00, &
         10.5D+00,  5.86602539D+00, &
         10.0D+00,  5.0D+00, &
         11.0D+00,  5.0D+00, &
         10.5D+00, 15.0D+00, &
         10.0D+00,  5.0D+00, &
         11.0D+00,  5.0D+00, &
        20.0D+00,   7.0D+00 /), (/ m, 3, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2101'
  write ( *, '(a)' ) &
    '  For a triangle in 2D, the circumenter can be computed by:'
  write ( *, '(a)' ) '  TRIANGLE_CIRCUMCENTER_2D;'
  write ( *, '(a)' ) '  TRIANGLE_CIRCUMCENTER_2D_2;'
  write ( *, '(a)' ) '  TRIANGLE_CIRCUMCENTER (any dimension);'

  do test = 1, test_num

    t(1:m,1:3) = t_test(1:m,1:3,test)

    call r8mat_transpose_print ( m, 3, t, '  Triangle vertices:' )

    call triangle_circumcenter_2d ( t, pc )

    call r8vec_print ( m, pc, '  Circumcenter by TRIANGLE_CIRCUMCENTER_2D:' )

    call triangle_circumcenter_2d_2 ( t, pc )

    call r8vec_print ( m, pc, '  Circumcenter by TRIANGLE_CIRCUMCENTER_2D_2:' )

    call triangle_circumcenter ( m, t, pc )

    call r8vec_print ( m, pc, '  Circumcenter by TRIANGLE_CIRCUMCENTER:' )

  end do

  return
end
subroutine test21011 ( )

!*****************************************************************************80
!
!! TEST21011 tests TRIANGLE_CIRCUMCENTER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m1 = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ), allocatable :: a12(:,:)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m2
  real ( kind = 8 ), allocatable :: o1(:)
  real ( kind = 8 ), allocatable :: o2(:)
  real ( kind = 8 ) pc1(m1)
  real ( kind = 8 ), allocatable :: pc2(:)
  real ( kind = 8 ) r8vec_norm_affine
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t1(m1,3)
  real ( kind = 8 ), allocatable :: t2(:,:)
  real ( kind = 8 ), dimension(m1,3,test_num) :: t_test = reshape ( (/ &
         10.0D+00,  5.0D+00, &
         11.0D+00,  5.0D+00, &
         10.0D+00,  6.0D+00, &
         10.0D+00,  5.0D+00, &
         11.0D+00,  5.0D+00, &
         10.5D+00,  5.86602539D+00, &
         10.0D+00,  5.0D+00, &
         11.0D+00,  5.0D+00, &
         10.5D+00, 15.0D+00, &
         10.0D+00,  5.0D+00, &
         11.0D+00,  5.0D+00, &
        20.0D+00,   7.0D+00 /), (/ m1, 3, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21011'
  write ( *, '(a)' ) &
    '  For a triangle in M dimensions, the circumenter can be computed by:'
  write ( *, '(a)' ) '  TRIANGLE_CIRCUMCENTER;'

!
!  Vary the dimension.
!
  do m2 = 2, 5

    seed = 123456789

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  M2 = ', m2

    allocate ( a12(1:m2,1:m1) )
    allocate ( o1(1:m1) )
    allocate ( o2(1:m2) )
    allocate ( pc2(1:m2) )
    allocate ( t2(1:m2,1:3) )
!
!  Randomly choose a mapping P2 = O2 + A12 * ( P1 - O1 )
!
    call r8mat_uniform_01 ( m2, m1, seed, a12 )
    call r8vec_uniform_01 ( m1, seed, o1 )
    call r8vec_uniform_01 ( m2, seed, o2 )
!
!  Map each M1-dimensional triangle into M2 space.
!
    do test = 1, test_num

      t1(1:m1,1:3) = t_test(1:m1,1:3,test)

      do j = 1, 3
        t1(1:m1,j) = t1(1:m1,j) - o1(1:m1)
      end do

      t2(1:m2,1:3) = matmul ( a12(1:m2,1:m1), t1(1:m1,1:3) )

      do j = 1, 3
        t2(1:m2,j) = t2(1:m2,j) + o2(1:m2)
      end do

      call triangle_circumcenter ( m2, t2, pc2 )

      call r8vec_print ( m2, pc2, '  Circumcenter by TRIANGLE_CIRCUMCENTER:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Distances from circumcenter to vertices:'
      write ( *, '(a)' ) ' '
      do j = 1, 3
        write ( *, '(2x,g14.6)' ) r8vec_norm_affine ( m2, pc2,      t2(1:m2,j) )
      end do

    end do

    deallocate ( a12 )
    deallocate ( o1 )
    deallocate ( o2 )
    deallocate ( pc2 )
    deallocate ( t2 )

  end do

  return
end
subroutine test2067 ( )

!*****************************************************************************80
!
!! TEST2067 tests TRIANGLE_CIRCUMCIRCLE_2D and TRIANGLE_CIRCUMCIRCLE_2D_2;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ), dimension(dim_num,3,test_num) :: t_test = reshape ( (/ &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
         0.0D+00,  1.0D+00, &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
         0.5D+00,  0.86602539D+00, &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
         0.5D+00, 10.0D+00, &
         0.0D+00,  0.0D+00, &
         1.0D+00,  0.0D+00, &
        10.0D+00,  2.0D+00 /), (/ dim_num, 3, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2067'
  write ( *, '(a)' ) '  For a triangle in 2D:'
  write ( *, '(a)' ) '  TRIANGLE_CIRCUMCIRCLE_2D computes the circumcenter.'
  write ( *, '(a)' ) '  TRIANGLE_CIRCUMCIRCLE_2D_2 computes the circumcenter.'

  do test = 1, test_num

    t(1:dim_num,1:3) = t_test(1:dim_num,1:3,test)

    call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

    call triangle_circumcircle_2d ( t, r, pc )

    call r8vec_print ( dim_num, pc, '  Circumcenter' )
    write ( *, '(a,g14.6)' ) '  Circumradius: ', r

    call triangle_circumcircle_2d_2 ( t, r, pc )

    call r8vec_print ( dim_num, pc, '  Circumcenter2' )
    write ( *, '(a,g14.6)' ) '  Circumradius2: ', r

  end do

  return
end
subroutine test21015 ( )

!*****************************************************************************80
!
!! TEST21015 tests TRIANGLE_CIRCUMRADIUS_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 July 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) r
  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ), dimension(dim_num,3,test_num) :: t_test = reshape ( (/ &
        0.0D+00, 0.0D+00, &
        1.0D+00, 0.0D+00, &
        0.0D+00, 1.0D+00, &
        0.0D+00, 0.0D+00, &
        1.0D+00, 0.0D+00, &
        0.5D+00, 0.86602539D+00, &
        0.0D+00,  0.0D+00, &
        1.0D+00,  0.0D+00, &
        0.5D+00, 10.0D+00, &
         0.0D+00, 0.0D+00, &
         1.0D+00, 0.0D+00, &
        10.0D+00, 2.0D+00 /), (/ dim_num, 3, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21015'
  write ( *, '(a)' ) '  For a triangle in 2D:'
  write ( *, '(a)' ) '  TRIANGLE_CIRCUMRADIUS_2D computes the circumradius.'

  do test = 1, test_num

    t(1:dim_num,1:3) = t_test(1:dim_num,1:3,test)

    call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

    call triangle_circumradius_2d ( t, r )

    write ( *, '(a,g14.6)' ) '  Circumradius: ', r

  end do

  return
end
subroutine test2068 ( )

!*****************************************************************************80
!
!! TEST2068 tests TRIANGLE_CONTAINS_LINE_EXP_3D.
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

  logical inside
  real ( kind = 8 ), dimension(dim_num) :: p1 = (/ &
    3.0D+00, 0.0D+00, -7.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: p2 = (/ &
    5.0D+00, 1.0D+00, -2.0D+00 /)
  real ( kind = 8 ) pint(dim_num)
  real ( kind = 8 ), dimension(dim_num,3) :: t = reshape ( (/  &
    8.0D+00, 4.0D+00, 2.0D+00, &
    9.0D+00, 0.0D+00, 5.0D+00, &
    2.0D+00, 1.0D+00, 2.0D+00 /), (/ dim_num, 3 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2068'
  write ( *, '(a)' ) '  TRIANGLE_CONTAINS_LINE_EXP_3D determines whether '
  write ( *, '(a)' ) '  a triangle "contains" an explicit line in 3D.'

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

  call r8vec_print ( dim_num, p1, '  Line point P1:' )
  call r8vec_print ( dim_num, p2, '  Line point P2:' )

  call triangle_contains_line_exp_3d ( t, p1, p2, inside, pint )
 
  if ( inside ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The triangle contains the line.'
    call r8vec_print ( dim_num, pint, '  Intersection point:' )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The triangle does not contain the line.'
    call r8vec_print ( dim_num, pint, '  The intersection point:' )
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Expected answer:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    The triangle contains the line, and'
  write ( *, '(a)' ) '    the intersection point is at:'
  write ( *, '(a)' ) '      ( 7, 2, 3 ).'
 
  return
end
subroutine test2069 ( )

!*****************************************************************************80
!
!! TEST2069 tests TRIANGLE_CONTAINS_LINE_PAR_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  logical inside
  real ( kind = 8 ) norm
  real ( kind = 8 ), dimension(dim_num) :: p0 = (/ &
    3.0D+00, 0.0D+00, -7.0D+00 /)
  real ( kind = 8 ), dimension(dim_num) :: pd = (/ &
    2.0D+00, 1.0D+00, 5.0D+00 /)
  real ( kind = 8 ) pint(dim_num)
  real ( kind = 8 ), dimension(dim_num,3) :: t = reshape ( (/  &
    8.0D+00, 4.0D+00, 2.0D+00, &
    9.0D+00, 0.0D+00, 5.0D+00, &
    2.0D+00, 1.0D+00, 2.0D+00 /), (/ dim_num, 3 /) )
  logical triangle_contains_line_exp_3d

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2069'
  write ( *, '(a)' ) '  TRIANGLE_CONTAINS_LINE_PAR_3D determines whether '
  write ( *, '(a)' ) '  a triangle "contains" a parametric line in 3D.'

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

  norm = sqrt ( sum ( pd(1:dim_num)**2 ) )
! pd(1:dim_num) = pd(1:dim_num) / norm

  call r8vec_print ( dim_num, p0, '  Parametric base point P0:' )
  call r8vec_print ( dim_num, pd, '  Parametric direction PD:' )

  call triangle_contains_line_par_3d ( t, p0, pd, inside, pint )
 
  if ( inside ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The triangle contains the line.'
    call r8vec_print ( dim_num, pint, '  Intersection point:' )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The triangle does not contain the line.'
    call r8vec_print ( dim_num, pint, '  The intersection point:' )
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Expected answer:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    The triangle contains the line, and'
  write ( *, '(a)' ) '    the intersection point is at:'
  write ( *, '(a)' ) '      ( 7, 2, 3 ).'
 
  return
end
subroutine test207 ( )

!*****************************************************************************80
!
!! TEST207 tests TRIANGLE_CONTAINS_POINT_2D_*.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 7

  logical inside1
  logical inside2
  logical inside3
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p_test = reshape ( (/ &
     0.25D+00,   0.25D+00, &
     0.75D+00,   0.25D+00, &
     1.00D+00,   1.00D+00, &
    11.00D+00,   0.50D+00, &
     0.00D+00,   1.00D+00, &
     0.50D+00, -10.00D+00, &
     0.60D+00,   0.60D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num, 3 ) :: t = reshape ( (/ &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00 /), (/ dim_num, 3 /) )
  real ( kind = 8 ), dimension ( dim_num, 3 ) :: t2
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST207'
  write ( *, '(a)' ) '  For a triangle in 2D,'
  write ( *, '(a)' ) '  TRIANGLE_CONTAINS_POINT_2D_1 reports if a point '
  write ( *, '(a)' ) '  is inside a triangle (and doesn''t care about'
  write ( *, '(a)' ) '  the ordering of the vertices);'
  write ( *, '(a)' ) '  TRIANGLE_CONTAINS_POINT_2D_2 reports if a point '
  write ( *, '(a)' ) '  is inside a triangle (and DOES care about'
  write ( *, '(a)' ) '  the ordering of the vertices);'
  write ( *, '(a)' ) '  TRIANGLE_CONTAINS_POINT_2D_3 reports if a point '
  write ( *, '(a)' ) '  is inside a triangle (and doesn''t care about'
  write ( *, '(a)' ) '  the ordering of the vertices);'

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X       Y     In1  In2  In3'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)
 
    call triangle_contains_point_2d_1 ( t, p, inside1 )
    call triangle_contains_point_2d_2 ( t, p, inside2 )
    call triangle_contains_point_2d_3 ( t, p, inside3 )

    write ( *, '(2x,2f8.3,5x,l1,4x,l1,4x,l1)' ) &
      p(1:dim_num), inside1, inside2, inside3

  end do
!
!  Make a copy of the triangle with vertices in reverse order.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat the test, but reverse the triangle vertex'
  write ( *, '(a)' ) '  ordering.'
 
  do j = 1, 3
    t2(1:2,j) = t(1:2,4-j)
  end do

  call r8mat_transpose_print ( dim_num, 3, t2, &
    '  Triangle vertices (reversed):' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X       Y     In1  In2  In3'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)
 
    call triangle_contains_point_2d_1 ( t2, p, inside1 )
    call triangle_contains_point_2d_2 ( t2, p, inside2 )
    call triangle_contains_point_2d_3 ( t2, p, inside3 )

    write ( *, '(2x,2f8.3,5x,l1,4x,l1,4x,l1)' ) &
      p(1:dim_num), inside1, inside2, inside3

  end do
 
  return
end
subroutine test2075 ( )

!*****************************************************************************80
!
!! TEST2075 tests TRIANGLE_DIAMETER_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 July 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) diameter
  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ), dimension(dim_num,3,test_num) :: t_test = reshape ( (/ &
     4.0D+00, 2.0D+00, &
     1.0D+00, 5.0D+00, &
    -2.0D+00, 2.0D+00, &
     4.0D+00, 2.0D+00, &
     5.0D+00, 4.0D+00, &
     6.0D+00, 6.0D+00, &
     4.0D+00, 2.0D+00, &
     1.0D+00, 5.0D+00, &
     4.0D+00, 2.0D+00 /), (/ dim_num, 3, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2075'
  write ( *, '(a)' ) '  TRIANGLE_DIAMETER_2D computes the diameter of '
  write ( *, '(a)' ) '  the SMALLEST circle around the triangle.'

  do test = 1, test_num

    t(1:dim_num,1:3) = t_test(1:dim_num,1:3,test)

    call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

    call triangle_diameter_2d ( t, diameter )

    write ( *, '(a,g14.6)' ) '  Diameter =       ', diameter

  end do

  return
end
subroutine test208 ( )

!*****************************************************************************80
!
!! TEST208 tests TRIANGLE_GRIDPOINTS_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: grid_max = 50

  real ( kind = 8 ) g(dim_num,grid_max)
  integer ( kind = 4 ) grid_num
  integer ( kind = 4 ) sub_num
  real ( kind = 8 ), dimension ( dim_num, 3 ) :: t = reshape ( (/ &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00 /), (/ dim_num, 3 /) )

  sub_num = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST208'
  write ( *, '(a)' ) '  For a triangle in 2D,'
  write ( *, '(a)' ) '  TRIANGLE_GRIDPOINTS_2D produces a set of'
  write ( *, '(a)' ) '  gridpoints in or on the triangle.'

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

  call triangle_gridpoints_2d ( t, sub_num, grid_max, grid_num, g )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of grid points is ', grid_num

  call r8mat_print ( dim_num, grid_num, g, '  Grid points: ' )

  return
end
subroutine test2102 ( )

!*****************************************************************************80
!
!! TEST2102 tests TRIANGLE_INCENTER_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ), dimension(dim_num,3,test_num) :: t_test = reshape ( (/ &
        0.0D+00, 0.0D+00, &
        1.0D+00, 0.0D+00, &
        0.0D+00, 1.0D+00, &
        0.0D+00, 0.0D+00, &
        1.0D+00, 0.0D+00, &
        0.5D+00, 0.86602539D+00, &
        0.0D+00,  0.0D+00, &
        1.0D+00,  0.0D+00, &
        0.5D+00, 10.0D+00, &
         0.0D+00, 0.0D+00, &
         1.0D+00, 0.0D+00, &
        10.0D+00, 2.0D+00 /), (/ dim_num, 3, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2102'
  write ( *, '(a)' ) '  For a triangle in 2D:'
  write ( *, '(a)' ) '  TRIANGLE_INCENTER_2D computes the incenter.'

  do test = 1, test_num

    t(1:dim_num,1:3) = t_test(1:dim_num,1:3,test)

    call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

    call triangle_incenter_2d ( t, pc )

    call r8vec_print ( dim_num, pc, '  Incenter' )

  end do

  return
end
subroutine test2070 ( )

!*****************************************************************************80
!
!! TEST2070 tests TRIANGLE_INCIRCLE_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ), dimension ( dim_num, 3 ) :: t = reshape ( (/ &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00 /), (/ dim_num, 3 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2070'
  write ( *, '(a)' ) '  For a triangle in 2D,'
  write ( *, '(a)' ) '  TRIANGLE_INCIRCLE_2D computes the incircle.'

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

  call triangle_incircle_2d ( t, r, pc )

  call r8vec_print ( dim_num, pc, '  Incenter' )

  write ( *, '(a,g14.6)' ) '  Incircle radius is ', r

  return
end
subroutine test20701 ( )

!*****************************************************************************80
!
!! TEST20701 tests TRIANGLE_INRADIUS_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) r
  real ( kind = 8 ), dimension ( dim_num, 3 ) :: t = reshape ( (/ &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00 /), (/ dim_num, 3 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20701'
  write ( *, '(a)' ) '  For a triangle in 2D,'
  write ( *, '(a)' ) '  TRIANGLE_INRADIUS_2D computes the inradius.'

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

  call triangle_inradius_2d ( t, r )

  write ( *, '(a,g14.6)' ) '  Incircle radius is ', r

  return
end
subroutine test2104 ( )

!*****************************************************************************80
!
!! TEST2104 tests TRIANGLE_LATTICE_LAYER_POINT_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2

  integer ( kind = 4 ) c(n+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) layer
  logical              more
  integer ( kind = 4 ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2104'
  write ( *, '(a)' ) '  TRIANGLE_LATTICE_LAYER_POINT_NEXT returns the next'
  write ( *, '(a)' ) '  point in a triangle lattice layer defined by:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    C(3) - 1 < X(1)/C(1) + X(2)/C(2) <= C(3).'

  c(1) = 2
  c(2) = 3
  v(1:n) = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  N = ', n
  write ( *, '(a)', ADVANCE = 'NO' ) '  C =       '
  do i = 1, n
    write ( *, '(2x,i4)', ADVANCE = 'NO' ) c(i)
  end do
  write ( *, '(a)', ADVANCE = 'YES' ) 

  do layer = 0, 4

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  Layer ', layer
    write ( *, '(a)' ) ' '

    c(3) = layer
    more = .false.
    i = 0

    do
      call triangle_lattice_layer_point_next ( c, v, more )
      if ( .not. more ) then
        write ( *, '(a)' ) '  No more.'
        exit
      end if
      i = i + 1
      write ( *, '(2x,i4,6x,10(2x,i4))' ) i, v(1:n)

    end do

  end do

  return
end
subroutine test2105 ( )

!*****************************************************************************80
!
!! TEST2105 tests TRIANGLE_LATTICE_POINT_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2

  integer ( kind = 4 ) c(n+1)
  integer ( kind = 4 ) i
  logical              more
  integer ( kind = 4 ) v(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2105'
  write ( *, '(a)' ) '  TRIANGLE_LATTICE_POINT_NEXT returns the next lattice'
  write ( *, '(a)' ) '  point in a triangle defined by:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    0 <= X(1)/C(1) + X(2)/C(2) <= C(3).'

  do i = 1, n + 1
    c(i) = n + 2 - i
  end do
  v(1:n) = 0
  more = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  N = ', n
  write ( *, '(a)', ADVANCE = 'NO' ) '  C =       '
  do i = 1, n + 1
    write ( *, '(2x,i4)', ADVANCE = 'NO' ) c(i)
  end do
  write ( *, '(a)', ADVANCE = 'YES' ) 
  write ( *, '(a)' ) ' '

  i = 0

  do
    call triangle_lattice_point_next ( c, v, more )
    if ( .not. more ) then
      write ( *, '(a)' ) '  No more.'
      exit
    end if
    i = i + 1
    write ( *, '(2x,i4,6x,10(2x,i4))' ) i, v(1:n)
  end do

  return
end
subroutine test211 ( )

!*****************************************************************************80
!
!! TEST211 tests TRIANGLE_ORIENTATION_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  integer ( kind = 4 ) triangle_orientation_2d
  integer ( kind = 4 ) i
  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ), dimension(dim_num,3,test_num) :: t_test = reshape ( (/ &
     4.0D+00,  2.0D+00, &
     1.0D+00,  5.0D+00, &
    -2.0D+00,  2.0D+00, &
     1.0D+00,  5.0D+00, &
     4.0D+00,  2.0D+00, &
     1.0D+00, -1.0D+00, &
     1.0D+00,  5.0D+00, &
     2.0D+00,  7.0D+00, &
     3.0D+00,  9.0D+00, &
     1.0D+00,  5.0D+00, &
     4.0D+00,  2.0D+00, &
     1.0D+00,  5.0D+00 /), (/ dim_num, 3, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST211'
  write ( *, '(a)' ) '  TRIANGLE_ORIENTATION_2D determines orientation'
  write ( *, '(a)' ) '  of a triangle.'

  do test = 1, test_num

    t(1:dim_num,1:3) = t_test(1:dim_num,1:3,test)

    i = triangle_orientation_2d ( t )

    call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

    if ( i == 0 ) then
      write ( *, '(a)' ) '  The points are counterclockwise.'
    else if ( i == 1 ) then
      write ( *, '(a)' ) '  The points are clockwise.'
    else if ( i == 2 ) then
      write ( *, '(a)' ) '  The points are colinear.'
    else if ( i == 3 ) then
      write ( *, '(a)' ) '  The points are not distinct.'
    else
      write ( *, '(a)' ) '  The return value makes no sense.'
    end if

  end do

  return
end
subroutine test2103 ( )

!*****************************************************************************80
!
!! TEST2103 tests TRIANGLE_ORTHOCENTER_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  logical              flag
  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ), dimension(dim_num,3,test_num) :: t_test = reshape ( (/ &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.0D+00,  1.0D+00, &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.5D+00,  0.86602539D+00, &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.5D+00, 10.0D+00, &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
    10.0D+00,  2.0D+00 /), (/ dim_num, 3, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2103'
  write ( *, '(a)' ) '  For a triangle in 2D:'
  write ( *, '(a)' ) '  TRIANGLE_ORTHOCENTER_2D computes the orthocenter.'

  do test = 1, test_num

    t(1:dim_num,1:3) = t_test(1:dim_num,1:3,test)

    call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

    call triangle_orthocenter_2d ( t, pc, flag )

    call r8vec_print ( dim_num, pc, '  Orthocenter' )

  end do

  return
end
subroutine test2071 ( )

!*****************************************************************************80
!
!! TEST2071 tests TRIANGLE_POINT_DIST_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 7

  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_signed
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = reshape ( (/ &
     0.25D+00,   0.25D+00, &
     0.75D+00,   0.25D+00, &
     1.00D+00,   1.00D+00, &
    11.00D+00,   0.50D+00, &
     0.00D+00,   1.00D+00, &
     0.50D+00, -10.00D+00, &
     0.60D+00,   0.60D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num, 3 ) :: t = reshape ( (/ &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00 /), (/ dim_num, 3 /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2071'
  write ( *, '(a)' ) '  For a triangle in 2D,'
  write ( *, '(a)' ) '  TRIANGLE_POINT_DIST_2D computes the distance'
  write ( *, '(a)' ) '  to a point;'

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           P            DIST'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call triangle_point_dist_2d ( t, p, dist ) 

    write ( *, '(2x,2f8.3,2x,f8.3)' ) p(1:dim_num), dist

  end do
 
  return
end
subroutine test20715 ( )

!*****************************************************************************80
!
!! TEST20715 tests TRIANGLE_POINT_DIST_SIGNED_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 7

  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_signed
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = reshape ( (/ &
     0.25D+00,   0.25D+00, &
     0.75D+00,   0.25D+00, &
     1.00D+00,   1.00D+00, &
    11.00D+00,   0.50D+00, &
     0.00D+00,   1.00D+00, &
     0.50D+00, -10.00D+00, &
     0.60D+00,   0.60D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ), dimension ( dim_num, 3 ) :: t = reshape ( (/ &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00 /), (/ dim_num, 3 /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20715'
  write ( *, '(a)' ) '  For a triangle in 2D,'
  write ( *, '(a)' ) '  TRIANGLE_POINT_DIST_SIGNED_2D computes signed'
  write ( *, '(a)' ) '  distance to a point;'

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           P       DIST_SIGNED'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call triangle_point_dist_signed_2d ( t, p, dist_signed )

    write ( *, '(2x,2f8.3,2x,f8.3)' ) p(1:dim_num), dist_signed

  end do
 
  return
end
subroutine test2095 ( )

!*****************************************************************************80
!
!! TEST2095 tests TRIANGLE_POINT_DIST_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) dist
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = reshape ( (/ &
    1.0D+00,       2.0D+00,       3.0D+00, &
    1.3535534D+00, 2.3535534D+00, 3.0D+00, &
    0.0D+00,       0.0D+00,       0.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) , dimension ( dim_num, 3 ) :: t = reshape ( (/ &
    1.0D+00,       2.0D+00,       3.0D+00, &
    2.4142137D+00, 3.4142137D+00, 3.0D+00, &
    1.7071068D+00, 2.7071068D+00, 4.0D+00 /), (/ dim_num, 3 /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2095'
  write ( *, '(a)' ) '  For a triangle in 3D:'
  write ( *, '(a)' ) '  TRIANGLE_POINT_DIST_3D computes the distance'
  write ( *, '(a)' ) '  to a point;'

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                   P                          DIST'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call triangle_point_dist_3d ( t, p, dist ) 

    write ( *, '(2x,3g12.4,2x,g14.6)' ) p(1:dim_num), dist

  end do

  return
end
subroutine test2072 ( )

!*****************************************************************************80
!
!! TEST2072 tests TRIANGLE_POINT_NEAR_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 7

  real ( kind = 8 ) dist
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: p_test = reshape ( (/ &
     0.25D+00,   0.25D+00, &
     0.75D+00,   0.25D+00, &
     1.00D+00,   1.00D+00, &
    11.00D+00,   0.50D+00, &
     0.00D+00,   1.00D+00, &
     0.50D+00, -10.00D+00, &
     0.60D+00,   0.60D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) pn(dim_num)
  real ( kind = 8 ), dimension ( dim_num, 3 ) :: t = reshape ( (/ &
    0.0D+00, 1.0D+00, &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00 /), (/ dim_num, 3 /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2072'
  write ( *, '(a)' ) '  For a triangle in 2D,'
  write ( *, '(a)' ) '  TRIANGLE_POINT_NEAR_2D computes the nearest'
  write ( *, '(a)' ) '  point to a point.'

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           P                PN'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p(1:dim_num) = p_test(1:dim_num,test)

    call triangle_point_near_2d ( t, p, pn, dist ) 

    write ( *, '(2x,2f8.3,2x,2f8.3)' ) p(1:dim_num), pn(1:dim_num)

  end do
 
  return
end
subroutine test2115 ( )

!*****************************************************************************80
!
!! TEST2115 tests TRIANGLE_QUALITY_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) quality
  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ), dimension (dim_num,3,test_num) :: t_test = reshape ( (/ &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.0D+00,  1.0D+00, &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.5D+00,  0.86602539D+00, &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.5D+00, 10.0D+00, &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
    10.0D+00,  2.0D+00 /), (/ dim_num, 3, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2115'
  write ( *, '(a)' ) '  For a triangle in 2D:'
  write ( *, '(a)' ) '  TRIANGLE_QUALITY_2D computes the quality.'

  do test = 1, test_num

    t(1:dim_num,1:3) = t_test(1:dim_num,1:3,test)

    call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

    call triangle_quality_2d ( t, quality )

    write ( *, '(a,g14.6)' ) '  Quality = ', quality

  end do

  return
end
subroutine test212 ( )

!*****************************************************************************80
!
!! TEST212 tests TRIANGLE_SAMPLE, TRIANGLE_XY_TO_XSI_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 10

  real ( kind = 8 ) p(dim_num)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ), dimension(dim_num,3) :: t = reshape ( (/ &
     4.0D+00, 2.0D+00, &
     1.0D+00, 5.0D+00, &
    -2.0D+00, 2.0D+00 /), (/ dim_num, 3 /) )
  integer ( kind = 4 ) test
  real ( kind = 8 ) xsi(dim_num+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST212'
  write ( *, '(a)' ) '  TRIANGLE_SAMPLE samples a triangle.'
  write ( *, '(a)' ) '  TRIANGLE_XY_TO_XSI_2D converts XY to XSI coordinates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We are computing the XSI coordinates just to verify'
  write ( *, '(a)' ) '  that the points are inside the triangle.'

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sample points (X,Y) and (XSI1,XSI2,XSI3) coordinates:'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    call triangle_sample ( t, 1, seed, p )
    call triangle_xy_to_xsi_2d ( t, p, xsi )
    write ( *, '(2x,2f8.4,4x,3f8.4)' ) p(1:dim_num), xsi(1:dim_num+1)
  end do

  return
end
subroutine test213 ( )

!*****************************************************************************80
!
!! TEST213 tests TRIANGLE_SAMPLE, TRIANGLE_XY_TO_XSI_2D, TRIANGLE_XSI_TO_XY_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 10

  integer ( kind = 4 ) i
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) p2(dim_num)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ), dimension(dim_num,3) :: t = reshape ( (/ &
     4.0D+00, 2.0D+00, &
     1.0D+00, 5.0D+00, &
    -2.0D+00, 2.0D+00 /), (/ dim_num, 3 /) )
  integer ( kind = 4 ) test
  real ( kind = 8 ) xsi(dim_num+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST213'
  write ( *, '(a)' ) '  TRIANGLE_SAMPLE samples a triangle.'
  write ( *, '(a)' ) '  TRIANGLE_XY_TO_XSI_2D converts XY to XSI coordinates.'
  write ( *, '(a)' ) '  TRIANGLE_XSI_TO_XY_2D converts XSI to XY coordinates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We verify that (X,Y) -> (XSI1,XSI2,XSI3) -> (X,Y)'
  write ( *, '(a)' ) '  works properly.'

  call r8mat_transpose_print ( dim_num, 3, t, '  Triangle vertices:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sample points:'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    if ( test == 1 ) then
      do i = 1, dim_num
        p(i) = sum ( t(i,1:3) ) / 3.0D+00
      end do
    else if ( test == 2 ) then
      p(1) = 3.0D+00
      p(2) = 0.0D+00
    else
      call triangle_sample ( t, 1, seed, p )
    end if

    call triangle_xy_to_xsi_2d ( t, p, xsi )

    call triangle_xsi_to_xy_2d ( t, xsi, p2 )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,2f8.4,4x,3f8.4)' ) p(1:dim_num), xsi(1:dim_num+1)
    write ( *, '(2x,2f8.4)' ) p2(1:dim_num)

  end do

  return
end
subroutine test219 ( )

!*****************************************************************************80
!
!! TEST219 tests TUBE_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 July 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) dist
  real ( kind = 8 ), dimension ( test_num ) :: dist_test = (/ &
    0.5D+00, 0.5D+00, 1.0D+00, 1.0D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nhi
  integer ( kind = 4 ) nlo
  integer ( kind = 4 ), dimension ( test_num ) :: n_test = (/ 4, 5, 5, 5 /)
  real ( kind = 8 ), allocatable, dimension(:,:) :: p
  real ( kind = 8 ), dimension ( dim_num, 19 ) :: p_test = reshape ( (/ &
     0.0D+00,  0.0D+00, &
     4.0D+00,  3.0D+00, &
     4.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00, &
     2.0D+00,  0.0D+00, &
     2.0D+00,  1.0D+00, &
     0.0D+00,  1.0D+00, &
     0.0D+00,  0.0D+00, &
    10.0D+00, 20.0D+00, &
    20.0D+00, 20.0D+00, &
    10.0D+00, 10.0D+00, &
    20.0D+00, 10.0D+00, &
    10.0D+00, 20.0D+00, &
     0.0D+00,  0.0D+00, &
    10.0D+00,  0.0D+00, &
    10.0D+00, 10.0D+00, &
    10.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00 /), (/ dim_num, 19 /) )
  real ( kind = 8 ), allocatable, dimension(:,:) :: p1
  real ( kind = 8 ), allocatable, dimension(:,:) :: p2
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST219'
  write ( *, '(a)' ) '  TUBE_2D computes corners of a tube of radius'
  write ( *, '(a)' ) '  DIST surrounding a sequence of points.'
 
  do test = 1, test_num
 
    n = n_test ( test )
    dist = dist_test(test)

    allocate ( p(1:dim_num,1:n) )
    nlo = sum ( n_test(1:test-1) ) + 1
    nhi = nlo + n - 1
  
    p(1:dim_num,1:n) = p_test(1:dim_num,nlo:nhi)

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Test ', test
    write ( *, '(a,i8)' ) '  Number of points N = ', n
    write ( *, '(a,g14.6)' ) '  Tube radius DIST = ', dist

    call r8mat_transpose_print ( dim_num, n, p, '  Points to surround:' )
 
    allocate ( p1(1:dim_num,1:n) )
    allocate ( p2(1:dim_num,1:n) )

    call tube_2d ( dist, n, p, p1, p2 )
 
    call r8mat_transpose_print ( dim_num, n, p1, '  P1:' )

    call r8mat_transpose_print ( dim_num, n, p2, '  P2:' )

    deallocate ( p )
    deallocate ( p1 )
    deallocate ( p2 )
 
  end do
 
  return
end
subroutine test220 ( )

!*****************************************************************************80
!
!! TEST220 tests VECTOR_DIRECTIONS_ND;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 5

  real ( kind = 8 ) angle(dim_num)
  real ( kind = 8 ) angle_degrees(dim_num)
  integer ( kind = 4 ) j
  real ( kind = 8 ) radians_to_degrees
  integer ( kind = 4 ) test
  real ( kind = 8 ) v(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: v_test = reshape ( (/ &
     1.0D+00,        0.0D+00, &
     1.7320508D+00,  1.0D+00, &
    -1.7320508D+00,  1.0D+00, &
    -1.7320508D+00, -1.0D+00, &
     1.7320508D+00, -1.0D+00 /), (/ dim_num, test_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST220'
  write ( *, '(a)' ) '  VECTOR_DIRECTIONS_ND computes the angles'
  write ( *, '(a)' ) '  that a vector makes with the axes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X        Y       AX       AY   ' // &
    '    AX       AY'
  write ( *, '(a)' ) '                     (__Radians___)' // &
    '  (___Degrees___)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    v(1:dim_num) = v_test(1:dim_num,test)

    call vector_directions_nd ( dim_num, v, angle )

    do j = 1, dim_num
      angle_degrees(j) = radians_to_degrees ( angle(j) )
    end do

    write ( *, '(2x,6f9.3)' ) &
      v(1:dim_num), angle(1:dim_num), angle_degrees(1:dim_num)
 
  end do
 
  return
end
subroutine test221 ( )

!*****************************************************************************80
!
!! TEST221 tests VECTOR_DIRECTIONS_ND;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) angle(dim_num)
  real ( kind = 8 ) angle_degrees(dim_num)
  integer ( kind = 4 ) j
  real ( kind = 8 ) radians_to_degrees
  integer ( kind = 4 ) test
  real ( kind = 8 ) v(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: v_test = reshape ( (/ &
    1.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 2.0D+00, 3.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00 /), (/ dim_num, test_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST221'
  write ( *, '(a)' ) '  VECTOR_DIRECTIONS_ND computes the angles'
  write ( *, '(a)' ) '  that a vector makes with the axes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X       Y       Z      AX      AY      AZ   ' // &
     '   AX      AY      AZ   '
  write ( *, '(a)' ) '                         (_____Radians_______)' // &
     ' (_______Degrees_______)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    v(1:dim_num) = v_test(1:dim_num,test)

    call vector_directions_nd ( dim_num, v, angle )

    do j = 1, dim_num
      angle_degrees(j) = radians_to_degrees ( angle(j) )
    end do

    write ( *, '(2x,9f8.3)' ) &
      v(1:dim_num), angle(1:dim_num), angle_degrees(1:dim_num)
 
  end do
 
  return
end
subroutine test222 ( )

!*****************************************************************************80
!
!! TEST222 tests VECTOR_ROTATE_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension(test_num) :: a_test = (/ &
    30.0D+00, -45.0D+00, 270.0D+00 /)
  real ( kind = 8 ) degrees_to_radians
  integer ( kind = 4 ) test
  real ( kind = 8 ) v(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: v_test = reshape ( (/ &
    1.0D+00, 0.0D+00, &
    0.0D+00, 2.0D+00, &
    1.0D+00, 1.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) w(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST222'
  write ( *, '(a)' ) '  VECTOR_ROTATE_2D rotates a vector through'
  write ( *, '(a)' ) '  a given angle around the origin.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X1      Y1   Angle      X2      Y2'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    v(1:dim_num) = v_test(1:dim_num,test)

    angle = degrees_to_radians ( a_test(test) )

    call vector_rotate_2d ( v, angle, w )

    write ( *, '(2x,5f8.3)') v(1:dim_num), a_test(test), w(1:dim_num)
 
  end do
 
  return
end
subroutine test2225 ( )

!*****************************************************************************80
!
!! TEST2225 tests VECTOR_ROTATE_3D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 5

  real ( kind = 8 ), dimension(dim_num) :: axis = (/ &
    1.0D+00, 1.0D+00, 1.0D+00 /)
  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( test_num ) :: a_test = (/ &
    30.0D+00, -45.0D+00, 90.0D+00, 270.0D+00, 30.0D+00 /)
  real ( kind = 8 ) degrees_to_radians
  integer ( kind = 4 ) test
  real ( kind = 8 ) v1(dim_num)
  real ( kind = 8 ), dimension(dim_num,test_num) :: v1_test = reshape ( (/ &
    1.0D+00,  0.0D+00,  0.0D+00, &
    0.0D+00,  2.0D+00,  0.0D+00, &
    0.0D+00,  0.0D+00,  3.0D+00, &
    1.0D+00,  1.0D+00,  1.0D+00, &
    1.0D+00,  1.0D+00, -2.0D+00 /), (/ dim_num, test_num /) )
  real ( kind = 8 ) v2(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2225'
  write ( *, '(a)' ) '  VECTOR_ROTATE_3D rotates a vector through'
  write ( *, '(a)' ) '  a given angle around the origin.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rotations will be about the following axis:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,3f8.3)' ) axis(1:dim_num)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '               V1            Angle             V2'
  write ( *, '(a)' ) &
    '    ----------------------  ------  ----------------------'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    v1(1:dim_num) = v1_test(1:dim_num,test)

    angle = degrees_to_radians ( a_test(test) )

    call vector_rotate_3d ( v1, axis, angle, v2 )

    write ( *, '(2x,7f8.3)') v1(1:dim_num), a_test(test), v2(1:dim_num)
 
  end do
!
!  Test using an axis that is not of unit length!
!
  axis(1:3) = (/ 0.0D+00, 0.0D+00, 2.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rotations will be about the following axis:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,3f8.3)' ) axis(1:dim_num)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '               V1            Angle             V2'
  write ( *, '(a)' ) &
    '    ----------------------  ------  ----------------------'
  write ( *, '(a)' ) ' '

  v1(1:3) = (/ 1.0D+00, 1.0D+00, 1.0D+00 /)

  angle = 90.0D+00
  angle = degrees_to_radians ( angle )

  call vector_rotate_3d ( v1, axis, angle, v2 )

  write ( *, '(2x,7f8.3)') v1(1:dim_num), 90.0, v2(1:dim_num)
 
  return
end
subroutine test223 ( )

!*****************************************************************************80
!
!! TEST223 tests VECTOR_ROTATE_BASE_2D;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) angle
  real ( kind = 8 ), dimension ( test_num ) :: a_test = (/ &
    30.0D+00, -45.0D+00, 270.0D+00, 20.0D+00 /) 
  real ( kind = 8 ) degrees_to_radians
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ), dimension(dim_num) :: pb = (/ 10.0D+00, 5.0D+00 /)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: p_test = reshape ( (/ &
    11.0D+00, 5.0D+00, &
    10.0D+00, 7.0D+00, &
    11.0D+00, 6.0D+00, &
    10.0D+00, 5.0D+00 /), (/ dim_num, test_num /) )
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST223'
  write ( *, '(a)' ) '  VECTOR_ROTATE_BASE_2D rotates a vector (X1,Y1)'
  write ( *, '(a)' ) '  through an angle around a base point (XB,YB).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        P1              PB       Angle          P2'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    p1 = p_test(1:dim_num,test)

    angle = degrees_to_radians ( a_test(test) )

    call vector_rotate_base_2d ( p1, pb, angle, p2 )

    write ( *, '(2x,7f8.3)' ) &
      p1(1:dim_num), pb(1:dim_num), a_test(test), p2(1:dim_num)
 
  end do
 
  return
end
subroutine test224 ( )

!*****************************************************************************80
!
!! TEST224 tests VECTOR_SEPARATION_ND;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: test_num = 5

  real ( kind = 8 ) radians_to_degrees
  integer ( kind = 4 ) test1
  integer ( kind = 4 ) test2
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta_deg
  real ( kind = 8 ) v1(dim_num)
  real ( kind = 8 ) v2(dim_num)
  real ( kind = 8 ), dimension ( dim_num, test_num ) :: v_test = reshape ( (/ &
    1.0D+00,  0.0D+00,  0.0D+00, &
    1.0D+00,  2.0D+00,  3.0D+00, &
    0.0D+00,  0.0D+00,  1.0D+00, &
   -3.0D+00,  2.0D+00, -1.0D+00, &
   -2.0D+00, -4.0D+00, -6.0D+00 /), (/ dim_num, test_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST224'
  write ( *, '(a)' ) '  VECTOR_SEPARATION_3D computes the separation angle'
  write ( *, '(a)' ) '  between two vectors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    -----Vector 1-----      -----Vector 2-----  ' // &
     '   Radians    Degrees'
  write ( *, '(a)' ) ' '

  do test1 = 1, test_num

    v1(1:dim_num) = v_test(1:dim_num,test1)

    do test2 = test1 + 1, test_num

      v2(1:dim_num) = v_test(1:dim_num,test2)

      call vector_separation_nd ( dim_num, v1, v2, theta )

      theta_deg = radians_to_degrees ( theta )

      write ( *, '(2x,6f8.3,f8.3,5x,f8.3)') &
        v1(1:dim_num), v2(1:dim_num), theta, theta_deg

    end do

  end do

  return
end
subroutine test2245 ( )

!*****************************************************************************80
!
!! TEST2245 tests VOXELS_DIST_L1_ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) dist
  integer ( kind = 4 ), dimension(dim_num) :: p1 = (/ 1, 1, 5 /)
  integer ( kind = 4 ), dimension(dim_num) :: p2 = (/ 9, 4, 4 /)
  integer ( kind = 4 ) voxels_dist_l1_nd

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST2245'
  write ( *, '(a)' ) '  VOXELS_DIST_L1_ND prints the voxels on a line in 3D.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  P1:'
  write ( *, '(4x,3i8)' ) p1(1:dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  P2:'
  write ( *, '(4x,3i8)' ) p2(1:dim_num)

  dist = voxels_dist_l1_nd ( dim_num, p1, p2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  L1 distance = ', dist

  return
end
subroutine test225 ( )

!*****************************************************************************80
!
!! TEST225 tests VOXELS_LINE_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 May 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) n
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: v
  integer ( kind = 4 ) voxels_dist_l1_nd
  integer ( kind = 4 ), dimension ( dim_num ) :: p1 = (/ 1, 1, 5 /)
  integer ( kind = 4 ), dimension ( dim_num ) :: p2 = (/ 9, 4, 4 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST225'
  write ( *, '(a)' ) '  VOXELS_LINE_3D computes the voxels on a line in 3D'
  write ( *, '(a)' ) '  starting at the first voxel, and heading towards'
  write ( *, '(a)' ) '  the second one.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Starting voxel:'
  write ( *, '(4x,3i8)' ) p1(1:dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  "Heading" voxel:'
  write ( *, '(4x,3i8)' ) p2(1:dim_num)

  n = voxels_dist_l1_nd ( dim_num, p1, p2 ) + 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of voxels we will compute is ', n

  allocate ( v(1:3,n) )

  call voxels_line_3d ( p1, p2, n, v )

  call i4mat_transpose_print ( 3, n, v, '  The voxels:' )

  deallocate ( v )

  return
end
subroutine test226 ( )

!*****************************************************************************80
!
!! TEST226 tests VOXELS_REGION_3D.
!
!  Discussion:
!
!    The test region is 8 by 9 by 1 voxels:
!
!    123456789
!  1 .........
!  2 ...11.1..
!  3 ..11111..
!  4 ...11.1..
!  5 ......1..
!  6 .11..11..
!  7 ..1......
!  8 .......1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: list_max = 100
  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: nx = 8 
  integer ( kind = 4 ), parameter :: ny = 9 
  integer ( kind = 4 ), parameter :: nz = 1 

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ishow(nx,ny,nz)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) list(list_max)
  integer ( kind = 4 ) list_num
  integer ( kind = 4 ) nelements
  integer ( kind = 4 ) region
  integer ( kind = 4 ) region_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST226'
  write ( *, '(a)' ) '  VOXELS_REGION_3D groups voxels into regions.'

  ishow(1:nx,1:ny,1:nz) = 0

  ishow(2,4,1) = 1
  ishow(2,5,1) = 1
  ishow(2,7,1) = 1

  ishow(3,3,1) = 1
  ishow(3,4,1) = 1
  ishow(3,5,1) = 1
  ishow(3,6,1) = 1
  ishow(3,7,1) = 1

  ishow(4,4,1) = 1
  ishow(4,5,1) = 1
  ishow(4,7,1) = 1

  ishow(5,7,1) = 1

  ishow(6,2,1) = 1
  ishow(6,3,1) = 1
  ishow(6,6,1) = 1
  ishow(6,7,1) = 1

  ishow(7,3,1) = 1

  ishow(8,8,1) = 1

  call voxels_region_3d ( list_max, nx, ny, nz, ishow, list_num, list, &
    region_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of regions found = ', region_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The nonzero ISHOW array elements are:'
  write ( *, '(a)' ) ' '

  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        l = ishow(i,j,k)
        if ( l /= 0 ) then
          write ( *, '(2x,4i8)' ) i, j, k, l
        end if
      end do
    end do
  end do

  if ( list_max < list_num ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The stack-based list of regions is unusable.'

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The stack-based list of regions is:'
    write ( *, '(a)' ) ' '

    region = region_num

    do while ( 0 < list_num )

      nelements = list(list_num)
      list_num = list_num - 1

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a,i8,a)' ) &
        '  Region ', region, ' includes ', nelements, ' voxels:'
      write ( *, '(a)' ) ' '

      do l = 1, nelements
        k = list(list_num)
        list_num = list_num - 1
        j = list(list_num)
        list_num = list_num - 1
        i = list(list_num)
        list_num = list_num - 1
        write ( *, '(2x,3i8)' ) i, j, k
      end do

      region = region - 1

    end do

  end if

  return
end
subroutine test227 ( )

!*****************************************************************************80
!
!! TEST227 tests VOXELS_STEP_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) i
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) jnc
  integer ( kind = 4 ) knc
  integer ( kind = 4 ) v1(dim_num)
  integer ( kind = 4 ) v2(dim_num)
  integer ( kind = 4 ) v3(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST227'
  write ( *, '(a)' ) '  VOXELS_STEP_3D steps along a line from'
  write ( *, '(a)' ) '  one voxel to another.'

  v1(1:dim_num) = (/ 1, 1, 5 /)
  v2(1:dim_num) = v1(1:dim_num)

  inc = 7
  jnc = 3
  knc = -1

  write ( *, '(a)' ) ' '
  write ( *, '(2x,i4,2x,3i8)' ) 0, v2(1:dim_num)

  do i = 1, 10
    call voxels_step_3d ( v1, v2, inc, jnc, knc, v3 )
    write ( *, '(2x,i4,2x,3i8)' ) i, v3(1:dim_num)
    v2(1:dim_num) = v3(1:dim_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now, as a check, reverse direction and return.'
  write ( *, '(a)' ) ' '

  v1(1:dim_num) = v2(1:dim_num)

  inc = -inc
  jnc = -jnc
  knc = -knc

  v2(1:dim_num) = v1(1:dim_num)

  write ( *, '(2x,i4,2x,3i8)' ) 0, v2(1:dim_num)
  do i = 1, 10
    call voxels_step_3d ( v1, v2, inc, jnc, knc, v3 )
    write ( *, '(2x,i4,2x,3i8)' ) i, v3(1:dim_num)
    v2(1:dim_num) = v3(1:dim_num)
  end do

  return
end
