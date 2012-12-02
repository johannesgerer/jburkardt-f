program main

!*****************************************************************************80
!
!! MAIN is the main program for COLORS_PRB.
!
!  Discussion:
!
!    COLORS_PRB calls the test programs for COLORS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COLORS_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Simple tests for the COLORS routines.'
 
  call test01
  call test02
  call test025
  call test027
  call test03
  call test035
  call test04
  call test05
  call test06
  call test07
  call test08
  call test09
  call test10

  call test11
  call test12
  call test13
  call test131
  call test133
  call test135
  call test14
  call test15
  call test16
  call test17
  call test18
  call test185
  call test19
  call test20

  call test21
  call test22
  call test23
  call test24
  call test25
  call test26
  call test27
  call test28

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COLORS_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
 
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01

!*****************************************************************************80
!
!! TEST01 tests CMY_TO_RGB, RGB_TO_CMY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) c
  real ( kind = 8 ) g
  real ( kind = 8 ) g2
  integer itest
  real ( kind = 8 ) m
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  CMY_TO_RGB converts CMY to RGB colors.'
  write ( *, '(a)' ) '  RGB_TO_CMY converts RGB to CMY colors;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Rin     Gin     Bin     C       M       Y    ' &
    // '   Rout    Gout    Bout'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call rgb_test ( itest, r, g, b )

    if ( r < 0.0D+00 ) then
      exit
    end if

    call rgb_to_cmy ( r, g, b, c, m, y )

    call cmy_to_rgb ( c, m, y, r2, g2, b2 )

    write ( *, '(9f8.3)' ) r, g, b, c, m, y, r2, g2, b2
 
  end do

  return
end
subroutine test02

!*****************************************************************************80
!
!! TEST02 tests CMYK_TO_RGB, RGB_TO_CMYK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) c
  real ( kind = 8 ) g
  real ( kind = 8 ) g2
  integer itest
  real ( kind = 8 ) k
  real ( kind = 8 ) m
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  CMYK_TO_RGB converts CMYK to RGB colors.'
  write ( *, '(a)' ) '  RGB_TO_CMYK converts RGB to CMYK colors;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Rin     Gin     Bin     C       M       Y    ' &
    // '   K        Rout    Gout    Bout'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call rgb_test ( itest, r, g, b )

    if ( r < 0.0D+00 ) then
      exit
    end if

    call rgb_to_cmyk ( r, g, b, c, m, y, k )

    call cmyk_to_rgb ( c, m, y, k, r2, g2, b2 )

    write ( *, '(10f8.3)' ) r, g, b, c, m, y, k, r2, g2, b2
 
  end do

  return
end
subroutine test025

!*****************************************************************************80
!
!! TEST025 tests GRAYSCALE_LUV, LUV_TO_XYZ_CAP, XYZ_CAP_TO_RGB709.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 11

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  integer i
  real ( kind = 8 ) l(n)
  real ( kind = 8 ) r
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) xcap
  real ( kind = 8 ) xcapn
  real ( kind = 8 ) xn
  real ( kind = 8 ) ycap
  real ( kind = 8 ) ycapn
  real ( kind = 8 ) yn
  real ( kind = 8 ) zcap
  real ( kind = 8 ) zcapn
  real ( kind = 8 ) zn

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST025'
  write ( *, '(a)' ) '  GRAYSCALE_LUV determines the appropriate CIE LUV'
  write ( *, '(a)' ) '    coordinates for a grayscale.'
  write ( *, '(a)' ) '  LUV_TO_XYZ_CAP converts LUV to XYZ colors;'
  write ( *, '(a)' ) '  XYZ_CAP_TO_RGB709 converts XYZ to RGB709 colors;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i3)' ) '  The number of grays is ', n

  call grayscale_luv ( n, l, u, v )
!
!  Determine the (X,Y,Z) coordinates of the illuminant.
!
  ycapn = 1.0D+00

  call name_to_xyz ( 'D65', xn, yn, zn )

  call xyy_to_xyz_cap ( xn, yn, xcapn, ycapn, zcapn )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       L       u       v       ' // &
    'X       Y       Z       R       G       B'
  write ( *, '(a)' ) ' '

  do i = 1, n

    call luv_to_xyz_cap ( l(i), u(i), v(i), xcap, ycap, zcap, xcapn, &
      ycapn, zcapn )

    call xyz_cap_to_rgb709 ( xcap, ycap, zcap, r, g, b )

    write ( *, '(9f8.3)' ) l(i), u(i), v(i), xcap, ycap, zcap, r, g, b
 
  end do

  return
end
subroutine test027

!*****************************************************************************80
!
!! TEST027 tests GRAYSCALE_RGB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 11

  real ( kind = 8 ) b(n)
  real ( kind = 8 ) g(n)
  integer i
  real ( kind = 8 ) r(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST027'
  write ( *, '(a)' ) '  GRAYSCALE_RGB determines the appropriate RGB'
  write ( *, '(a)' ) '    coordinates for a grayscale.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i3)' ) '  The number of grays is ', n

  call grayscale_rgb ( n, r, g, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       R       G       B'
  write ( *, '(a)' ) ' '

  do i = 1, n

    write ( *, '(3f8.3)' ) r(i), g(i), b(i)
 
  end do

  return
end
subroutine test03

!*****************************************************************************80
!
!! TEST03 tests HLS_TO_RGB, RGB_TO_HLS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) g
  real ( kind = 8 ) g2
  real ( kind = 8 ) h
  integer itest
  real ( kind = 8 ) l
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  HLS_TO_RGB converts HLS to RGB colors.'
  write ( *, '(a)' ) '  RGB_TO_HLS converts RGB to HLS colors;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Rin     Gin     Bin     H       L       S    ' &
    // '   Rout    Gout    Bout'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call rgb_test ( itest, r, g, b )

    if ( r < 0.0D+00 ) then
      exit
    end if

    call rgb_to_hls ( r, g, b, h, l, s )

    call hls_to_rgb ( h, l, s, r2, g2, b2 )

    write ( *, '(9f8.3)' ) r, g, b, h, l, s, r2, g2, b2
 
  end do

  return
end
subroutine test04

!*****************************************************************************80
!
!! TEST04 tests HSV_TO_RGB, RGB_TO_HSV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) g
  real ( kind = 8 ) g2
  real ( kind = 8 ) h
  integer itest
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) s
  real ( kind = 8 ) v

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  HSV_TO_RGB converts HSV to RGB colors.'
  write ( *, '(a)' ) '  RGB_TO_HSV converts RGB to HSV colors;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Rin     Gin     Bin     H       S       V    ' &
    // '   Rout    Gout    Bout'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call rgb_test ( itest, r, g, b )

    if ( r < 0.0D+00 ) then
      exit
    end if

    call rgb_to_hsv ( r, g, b, h, s, v )

    call hsv_to_rgb ( h, s, v, r2, g2, b2 )

    write ( *, '(9f8.3)' ) r, g, b, h, s, v, r2, g2, b2
 
  end do

  return
end
subroutine test035

!*****************************************************************************80
!
!! TEST035 tests HSI_TO_RGB, RGB_TO_HSI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) g
  real ( kind = 8 ) g2
  real ( kind = 8 ) h
  real ( kind = 8 ) i
  integer itest
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST035'
  write ( *, '(a)' ) '  RGB_TO_HSI converts RGB to HSI colors;'
  write ( *, '(a)' ) '  HSI_TO_RGB converts HSI to RGB colors;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       Rin     Gin     Bin     ' // &
    'H       S       I       Rout    Gout    Bout'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call rgb_test ( itest, r, g, b )

    if ( r < 0.0D+00 ) then
      exit
    end if

    call rgb_to_hsi ( r, g, b, h, s, i )

    call hsi_to_rgb ( h, s, i, r2, g2, b2 )

    write ( *, '(9f8.3)' ) r, g, b, h, s, i, r2, g2, b2
 
  end do

  return
end
subroutine test05

!*****************************************************************************80
!
!! TEST05 tests LAB_TO_XYZ_CAP, XYZ_CAP_TO_LAB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ntest = 17

  real ( kind = 8 ) astar
  real ( kind = 8 ) bstar
  integer itest
  real ( kind = 8 ) lstar
  real ( kind = 8 ) nm
  real ( kind = 8 ), parameter :: nmhi = 700.0D+00
  real ( kind = 8 ), parameter :: nmlo = 380.0D+00
  real ( kind = 8 ) temp
  real ( kind = 8 ) x
  real ( kind = 8 ) xcap
  real ( kind = 8 ) xcap2
  real ( kind = 8 ) xcapn
  real ( kind = 8 ) xn
  real ( kind = 8 ) y
  real ( kind = 8 ) ycap
  real ( kind = 8 ) ycap2
  real ( kind = 8 ) ycapn
  real ( kind = 8 ) yn
  real ( kind = 8 ) z
  real ( kind = 8 ) zcap
  real ( kind = 8 ) zcap2
  real ( kind = 8 ) zcapn
  real ( kind = 8 ) zn

  ycapn = 1.0D+00

  call name_to_xyz ( 'D65', xn, yn, zn )

  call xyy_to_xyz_cap ( xn, yn, xcapn, ycapn, zcapn )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  LAB_TO_XYZ_CAP converts L*a*b* to XYZ colors.'
  write ( *, '(a)' ) '  XYZ_CAP_TO_LAB converts XYZ to L*a*b* colors;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Illuminant XYZ color coordinates:'
  write ( *, '(3g14.6)' ) xcapn, ycapn, zcapn 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NM  Xin     Yin     Zin       L*       a*' &
    // '       b*      Xout   Yout   Zout'
  write ( *, '(a)' ) ' '

  do itest = 1, ntest

    nm = ( real ( ntest - itest,     kind = 8 ) * nmlo &
         + real (         itest - 1, kind = 8 ) * nmhi ) &
         / real ( ntest         - 1, kind = 8 )

    call nm_to_xyz ( nm, x, y, z )
    ycap = 0.95D+00 * ycapn
    call xyy_to_xyz_cap ( x, y, xcap, ycap, zcap )
!
!  Make sure every component of (XCAP,YCAP,ZCAP) is smaller than
!  (XCAPN,YCAPN,ZCAPN).
!
    temp = 1.0D+00
    if ( 0.0D+00 < xcap .and. 0.0D+00 < xcapn ) then
      temp = min ( temp, 0.95D+00 * xcapn / xcap )
    end if
    if ( 0.0D+00 < ycap .and. 0.0D+00 < ycapn ) then
      temp = min ( temp, 0.95D+00 * ycapn / ycap )
    end if
    if ( 0.0D+00 < zcap .and. 0.0D+00 < zcapn ) then
      temp = min ( temp, 0.95D+00 * zcapn / zcap )
    end if

    xcap = xcap * temp
    ycap = ycap * temp
    zcap = zcap * temp

    call xyz_cap_to_lab ( xcap, ycap, zcap, xcapn, ycapn, zcapn, &
      lstar, astar, bstar )

    call lab_to_xyz_cap ( lstar, astar, bstar, xcap2, ycap2, zcap2, &
      xcapn, ycapn, zcapn )

    write ( *, '(f5.1,3f8.3,3f9.3,3f8.3)' ) &
      nm, xcap, ycap, zcap, lstar, astar, bstar, xcap2, ycap2, zcap2
 
  end do

  return
end
subroutine test06

!*****************************************************************************80
!
!! TEST06 tests LCC_TO_RGBPRIME, RGBPRIME_TO_LCC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) bprime
  real ( kind = 8 ) bprime2
  real ( kind = 8 ) chroma1
  real ( kind = 8 ) chroma2
  real ( kind = 8 ) g
  real ( kind = 8 ) gprime
  real ( kind = 8 ) gprime2
  integer itest
  real ( kind = 8 ) luma
  real ( kind = 8 ) r
  real ( kind = 8 ) rprime
  real ( kind = 8 ) rprime2
  real ( kind = 8 ) yb
  real ( kind = 8 ) yg
  real ( kind = 8 ) yr

  yr = 0.299D+00
  yg = 0.587D+00
  yb = 0.114D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  LCC_TO_RGBPRIME converts LCC to R''G''B'' colors;'
  write ( *, '(a)' ) '  RGBPRIME_TO_LCC converts R''G''B'' to LCC colors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    R''in    G''in    B''in    Luma    Chroma1 ' &
    // 'Chroma2 R''out   G''out   B''out'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call rgb_test ( itest, r, g, b )

    if ( r < 0.0D+00 ) then
      exit
    end if

    call rgb_to_rgbprime ( r, g, b, rprime, gprime, bprime )

    call rgbprime_to_lcc ( rprime, gprime, bprime, yr, yg, yb, luma, chroma1, &
      chroma2 )

    call lcc_to_rgbprime ( luma, chroma1, chroma2, yr, yg, yb, &
      rprime2, gprime2, bprime2 )

    write ( *, '(9f8.3)' ) rprime, gprime, bprime, luma, chroma1, &
      chroma2, rprime2, gprime2, bprime2

  end do

  return
end
subroutine test07

!*****************************************************************************80
!
!! TEST07 tests LCC_TO_YCBCR, YCBCR_TO_LCC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) chroma1in
  real ( kind = 8 ) chroma1out
  real ( kind = 8 ) chroma2in
  real ( kind = 8 ) chroma2out
  real ( kind = 8 ) cb
  real ( kind = 8 ) cr
  integer itest
  real ( kind = 8 ) lumain
  real ( kind = 8 ) lumaout
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  LCC_TO_YCBCR converts LCC to YCBCR colors;'
  write ( *, '(a)' ) '  YCBCR_TO_LCC converts YCBCR to LCC colors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Luma   Chroma1    Chroma2     Y''     Cb' &
    // '      Cr      Luma  Chroma1   Chroma2'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call ycc_test ( itest, lumain, chroma1in, chroma2in )

    if ( lumain < 0.0D+00 ) then
      exit
    end if

    call lcc_to_ycbcr ( lumain, chroma1in, chroma2in, y, cb, cr )

    call ycbcr_to_lcc ( y, cb, cr, lumaout, chroma1out, chroma2out )

    write ( *, '(9f8.3)' ) lumain, chroma1in, chroma2in, y, cb, cr, lumaout, &
      chroma1out, chroma2out
 
  end do

  return
end
subroutine test08

!*****************************************************************************80
!
!! TEST08 tests LCC_TO_YCC, YCC_TO_LCC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) c1in
  real ( kind = 8 ) c1out
  real ( kind = 8 ) c2in
  real ( kind = 8 ) c2out
  real ( kind = 8 ) chroma1
  real ( kind = 8 ) chroma2
  integer itest
  real ( kind = 8 ) luma
  real ( kind = 8 ) yin
  real ( kind = 8 ) yout

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  LCC_TO_YCC converts LCC to PhotoYCC colors.'
  write ( *, '(a)' ) '  YCC_TO_LCC converts PhotoYCC to LCC colors;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Yin    C1in     C2in     Luma     Chroma1' &
    // '      Chroma2      Yout    C1out    C2out'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call ycc_test ( itest, yin, c1in, c2in )

    if ( yin < 0.0D+00 ) then
      exit
    end if

    call ycc_to_lcc ( yin, c1in, c2in, luma, chroma1, chroma2 )

    call lcc_to_ycc ( luma, chroma1, chroma2, yout, c1out, c2out )

    write ( *, '(9f8.3)' ) yin, c1in, c2in, luma, chroma1, chroma2, yout, &
      c1out, c2out
 
  end do

  return
end
subroutine test09

!*****************************************************************************80
!
!! TEST09 tests LIN_TO_NONLIN, NONLIN_TO_LIN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) rprime

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  LIN_TO_NONLIN converts linear to nonlinear RGB;'
  write ( *, '(a)' ) '  NONLIN_TO_LIN converts nonlinear to linear RGB.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Rin     R''      Rout'
  write ( *, '(a)' ) ' '

  do i = -1, 26

    r = real ( i, kind = 8 ) / 25.0D+00

    call lin_to_nonlin ( r, rprime )
    call nonlin_to_lin ( rprime, r2 )

    write ( *, '(9f8.3)' ) r, rprime, r2 
 
  end do

  return
end
subroutine test10

!*****************************************************************************80
!
!! TEST10 tests LUV_TO_XYZ_CAP, XYZ_CAP_TO_LUV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ntest = 17

  integer itest
  real ( kind = 8 ) lstar
  real ( kind = 8 ) nm
  real ( kind = 8 ), parameter :: nmhi = 700.0D+00
  real ( kind = 8 ), parameter :: nmlo = 380.0D+00
  real ( kind = 8 ) temp
  real ( kind = 8 ) unprime
  real ( kind = 8 ) ustar
  real ( kind = 8 ) vnprime
  real ( kind = 8 ) vstar
  real ( kind = 8 ) wnprime
  real ( kind = 8 ) x
  real ( kind = 8 ) xcap
  real ( kind = 8 ) xcap2
  real ( kind = 8 ) xcapn
  real ( kind = 8 ) xn
  real ( kind = 8 ) y
  real ( kind = 8 ) ycap
  real ( kind = 8 ) ycap2
  real ( kind = 8 ) ycapn
  real ( kind = 8 ) yn
  real ( kind = 8 ) z
  real ( kind = 8 ) zcap
  real ( kind = 8 ) zcap2
  real ( kind = 8 ) zcapn
  real ( kind = 8 ) zn

  ycapn = 1.0D+00
  call name_to_xyz ( 'D65', xn, yn, zn )
  call xyy_to_xyz_cap ( xn, yn, xcapn, ycapn, zcapn )
  call xyz_cap_to_uvwprime ( xcapn, ycapn, zcapn, unprime, vnprime, wnprime )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  LUV_TO_XYZ_CAP converts L*u*v* to XYZ colors.'
  write ( *, '(a)' ) '  XYZ_CAP_TO_LUV converts XYZ to L*u*v* colors;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Illuminant XYZ color coordinates:'
  write ( *, '(3g14.6)') xcapn, ycapn, zcapn 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Illuminant u''v''w'' color coordinates:'
  write ( *, '(3g14.6)' ) unprime, vnprime, wnprime
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   NM Xin     Yin     Zin       L*       u*' &
    // '       v*      Xout   Yout   Zout'
  write ( *, '(a)' ) ' '

  do itest = 1, ntest

    nm = ( real ( ntest - itest,     kind = 8 ) * nmlo &
         + real (         itest - 1, kind = 8 ) * nmhi ) &
         / real ( ntest         - 1, kind = 8 )
    call nm_to_xyz ( nm, x, y, z )
    ycap = 0.95D+00 * ycapn
    call xyy_to_xyz_cap ( x, y, xcap, ycap, zcap )
!
!  Make sure every component of (XCAP,YCAP,ZCAP) is smaller than
!  (XCAPN,YCAPN,ZCAPN).
!
    temp = 1.0D+00
    if ( 0.0D+00 < xcap .and. 0.0D+00 < xcapn ) then
      temp = min ( temp, 0.95D+00 * xcapn / xcap )
    end if
    if ( 0.0D+00 < ycap .and. 0.0D+00 < ycapn ) then
      temp = min ( temp, 0.95D+00 * ycapn / ycap )
    end if
    if ( 0.0D+00 < zcap .and. 0.0D+00 < zcapn ) then
      temp = min ( temp, 0.95D+00 * zcapn / zcap )
    end if

    xcap = xcap * temp
    ycap = ycap * temp
    zcap = zcap * temp

    call xyz_cap_to_luv ( xcap, ycap, zcap, xcapn, ycapn, zcapn, lstar, ustar, &
      vstar )

    call luv_to_xyz_cap ( lstar, ustar, vstar, xcap2, ycap2, zcap2, &
      xcapn, ycapn, zcapn )

    write ( *, '(f5.1,3f8.3,3f9.3,3f8.3)' ) &
      nm, xcap, ycap, zcap, lstar, ustar, vstar, xcap2, ycap2, zcap2
 
  end do

  return
end
subroutine test11

!*****************************************************************************80
!
!! TEST11 tests NAME_TO_PRIMARIES, PRIMARIES_TO_Y.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bx
  real ( kind = 8 ) by
  real ( kind = 8 ) gx
  real ( kind = 8 ) gy
  real ( kind = 8 ) rx
  real ( kind = 8 ) ry
  real ( kind = 8 ) wx
  real ( kind = 8 ) wy
  real ( kind = 8 ) yb
  real ( kind = 8 ) yg
  real ( kind = 8 ) yr

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  NAME_TO_PRIMARIES returns the CIE xy'
  write ( *, '(a)' ) '    chromaticities of the primaries and'
  write ( *, '(a)' ) '    reference white used for various'
  write ( *, '(a)' ) '    television standards;'
  write ( *, '(a)' ) '  PRIMARIES_TO_Y computes the coefficients in'
  write ( *, '(a)' ) '    the luminance function, given the'
  write ( *, '(a)' ) '    chromaticities of the three primaries,'
  write ( *, '(a)' ) '    and the reference white.'
!
!  Get primary values for NTSC..
!
  call name_to_primaries ( 'NTSC', rx, ry, gx, gy, bx, by, wx, wy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Primary definition:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  R primary        ', rx, ry
  write ( *, '(a,2g14.6)' ) '  G primary        ', gx, gy
  write ( *, '(a,2g14.6)' ) '  B primary	 ', bx, by
  write ( *, '(a,2g14.6)' ) '  Reference white: ', wx, wy

  call primaries_to_y ( rx, ry, gx, gy, bx, by, wx, wy, yr, yg, yb )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R luminance coefficient ', yr
  write ( *, '(a,g14.6)' ) '  G luminance coefficient ', yg
  write ( *, '(a,g14.6)' ) '  B luminance coefficient ', yb

  return
end
subroutine test12

!*****************************************************************************80
!
!! TEST12 tests NAME_TEST, NAME_TO_RGB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b2
  real ( kind = 8 ) g2
  integer itest
  character ( len = 30 ) name
  real ( kind = 8 ) r2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  NAME_TO_RGB converts a name to RGB colors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Name    Rout    Gout    Bout'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call name_test ( itest, name )

    if ( len_trim ( name ) == 0 ) then
      exit
    end if

    call name_to_rgb ( name, r2, g2, b2 )

    write ( *, '(a10,2x,3f8.3)' ) name, r2, g2, b2
 
  end do

  return
end
subroutine test13

!*****************************************************************************80
!
!! TEST13 tests NAME_TO_RGB, RGB_TO_NAME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ntest = 10

  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) g
  real ( kind = 8 ) g2
  integer itest
  character ( len = 30 ) name
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  NAME_TO_RGB converts a name to RGB colors.'
  write ( *, '(a)' ) '  RGB_TO_NAME converts RGB colors to a name;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Pick RGB at random.'
  write ( *, '(a)' ) '   R       G       B      Nearest Name   Rout  Gout  Bout'
  write ( *, '(a)' ) ' '

  do itest = 1, ntest

    call rgb_uniform ( seed, r, g, b )

    call rgb_to_name ( r, g, b, name )

    call name_to_rgb ( name, r2, g2, b2 )
 
    write ( *, '(3f8.3,2x,a10,2x,3f8.3)' ) r, g, b, name, r2, g2, b2

  end do

  return
end
subroutine test131

!*****************************************************************************80
!
!! TEST131 tests NM_TO_XYZ, XYZ_CAP_TO_LUV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  real ( kind = 8 ) l
  real ( kind = 8 ) nm
  real ( kind = 8 ), parameter :: nm_hi = 780.0D+00
  real ( kind = 8 ), parameter :: nm_lo = 380.0D+00
  integer, parameter :: nm_step = 81
  integer, parameter :: purple_step = 21
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) xcap
  real ( kind = 8 ) xcapn
  real ( kind = 8 ) xcap_hi
  real ( kind = 8 ) xcap_lo
  real ( kind = 8 ) xn
  real ( kind = 8 ) ycap
  real ( kind = 8 ) ycap_hi
  real ( kind = 8 ) ycap_lo
  real ( kind = 8 ) ycapn
  real ( kind = 8 ) yn
  real ( kind = 8 ) zcap
  real ( kind = 8 ) zcap_hi
  real ( kind = 8 ) zcap_lo
  real ( kind = 8 ) zcapn
  real ( kind = 8 ) zn

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST131'
  write ( *, '(a)' ) '  NM_TO_XYZ_CAP converts pure light wavelengths to XYZ.'
  write ( *, '(a)' ) '  XYZ_CAP_TO_LUV converts XYZ to LUV.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (Then slide from red to blue along "nonpure" purple...)'
  write ( *, '(a)' ) ' '

  ycapn = 1.0D+00
  call name_to_xyz ( 'D65', xn, yn, zn )

  call xyy_to_xyz_cap ( xn, yn, xcapn, ycapn, zcapn )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   NM         X       Y       Z' // &
    '           L         u         v'
  write ( *, '(a)' ) ' '

  do i = 1, nm_step
    nm = ( real ( nm_step - i,     kind = 8 ) * nm_lo &
         + real (           i - 1, kind = 8 ) * nm_hi ) &
         / real ( nm_step     - 1, kind = 8 )
    call nm_to_xyz_cap ( nm, xcap, ycap, zcap )
    call xyz_cap_to_luv ( xcap, ycap, zcap, xcapn, ycapn, zcapn, l, u, v )
    write ( *, '(f5.1,2x,3f8.4,2x,3f10.4)' ) nm, xcap, ycap, zcap, l, u, v
  end do
!
!  Now fill in values along the "purple trail" from red back to blue.
!

  call nm_to_xyz_cap ( nm_hi, xcap_lo, ycap_lo, zcap_lo )
  call nm_to_xyz_cap ( nm_lo, xcap_hi, ycap_hi, zcap_hi )

  do i = 2, purple_step

    xcap = ( real ( purple_step - i,     kind = 8 ) * xcap_lo &
           + real (               i - 1, kind = 8 ) * xcap_hi ) &
           / real ( purple_step - 1, kind = 8 )

    ycap = ( real ( purple_step - i,     kind = 8 ) * ycap_lo &
           + real (               i - 1, kind = 8 ) * ycap_hi ) &
           / real ( purple_step     - 1, kind = 8 )

    zcap = ( real ( purple_step - i,     kind = 8 ) * zcap_lo &
           + real (               i - 1, kind = 8 ) * zcap_hi ) &
           / real ( purple_step     - 1, kind = 8 )

    call xyz_cap_to_luv ( xcap, ycap, zcap, xcapn, ycapn, zcapn, l, u, v )
    write ( *, '(5x,2x,3f8.4,2x,3f10.4)' )     xcap, ycap, zcap, l, u, v

  end do

  return
end
subroutine test133

!*****************************************************************************80
!
!! TEST133 tests RGB_NAMED_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ntest = 10

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  integer itest
  character ( len = 30 ) name
  real ( kind = 8 ) r
  integer seed

  seed = 123456789
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  RGB_NAMED_UNIFORM picks a random named color;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   R       G       B      Name'
  write ( *, '(a)' ) ' '

  do itest = 1, ntest

    call rgb_named_uniform ( seed, r, g, b, name )
 
    write ( *, '(3f8.3,2x,a)' ) r, g, b, name

  end do

  return
end
subroutine test135

!*****************************************************************************80
!
!! TEST135 tests NCS_TO_RGB, RGB_TO_NCS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  integer c
  character c1
  character c2
  real ( kind = 8 ) g
  real ( kind = 8 ) g2
  integer itest
  integer n
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  integer s
  integer w

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST135'
  write ( *, '(a)' ) '  NCS_TO_RGB converts NCS to RGB colors.'
  write ( *, '(a)' ) '  RGB_TO_NCS converts RGB to NCS colors;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  WARNING:'
  write ( *, '(a)' ) '    These routines are NOT worked out yet!'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Rin     Gin     Bin   C1 C2   N   C   S   W   ' &
    // 'Rout    Gout    Bout'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1

    call rgb_test ( itest, r, g, b )

    if ( r < 0.0D+00 ) then
      exit
    end if

    call rgb_to_ncs ( r, g, b, c1, c2, n, c, s )

    w = 100 - c - s

    call ncs_to_rgb ( c1, c2, n, c, s, r2, g2, b2 )

    write ( *, '(3f8.3,1x,a1,2x,a1,1x,i4,i4,i4,i4,3f8.3)' ) &
      r, g, b, c1, c2, n, c, s, w, r2, g2, b2
 
  end do

  return
end
subroutine test14

!*****************************************************************************80
!
!! TEST14 tests NM_TO_XYZ, XYY_TO_XYZ_CAP, XYZ_CAP_TO_XYY.
!
!  Discussion:
!
!    Thanks to Harald Anlauf, of the Technical University of Darmstadt,
!    for pointing out an error in an output format, 30 April 2002.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ntest = 17

  integer itest
  real ( kind = 8 ) nm
  real ( kind = 8 ), parameter :: nmhi = 700.0D+00
  real ( kind = 8 ), parameter :: nmlo = 380.0D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) x2
  real ( kind = 8 ) xcap
  real ( kind = 8 ) y
  real ( kind = 8 ) y2
  real ( kind = 8 ) ycap
  real ( kind = 8 ) zcap
  real ( kind = 8 ) z

  ycap = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  NM_TO_XYZ converts wavelengths to xyz colors;'
  write ( *, '(a)' ) '  XYY_TO_XYZ_CAP converts xyY to XYZ colors;'
  write ( *, '(a)' ) '  XYZ_CAP_TO_XYY converts XYZ to xyY colors.'
  write ( *, '(a,g14.6)' ) '  (Assume a luminosity of YCAP = ', ycap
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    NM      xin     yin     Yin     X       Y     ' &
    // '   Z      xout    yout    Yout'
  write ( *, '(a)' ) ' '

  do itest = 1, ntest

    nm = ( real ( ntest - itest,     kind = 8 ) * nmlo &
         + real (         itest - 1, kind = 8 ) * nmhi ) &
         / real ( ntest         - 1, kind = 8 )

    call nm_to_xyz ( nm, x, y, z )

    call xyy_to_xyz_cap ( x, y, xcap, ycap, zcap )

    call xyz_cap_to_xyy ( xcap, ycap, zcap, x2, y2 )

    write ( *, '(10f8.3)' ) nm, x, y, ycap, xcap, ycap, zcap, x2, y2, ycap
 
  end do

  return
end
subroutine test15

!*****************************************************************************80
!
!! TEST15 tests RGB_TO_HUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  real ( kind = 8 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  RGB_TO_HUE computes a hue between 0 and 1'
  write ( *, '(a)' ) '  corresponding to a given (R,G,B) color.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   R       G       B       H'
  write ( *, '(a)' ) ' '

  r = 1.00D+00
  g = 0.00D+00
  b = 0.00D+00

  call rgb_to_hue ( r, g, b, h )
  write ( *, '(4f8.4)' ) r, g, b, h

  r = 1.00D+00
  g = 1.00D+00
  b = 0.00D+00

  call rgb_to_hue ( r, g, b, h )
  write ( *, '(4f8.4)' ) r, g, b, h

  r = 0.00D+00
  g = 1.00D+00
  b = 0.00D+00

  call rgb_to_hue ( r, g, b, h )
  write ( *, '(4f8.4)' ) r, g, b, h

  r = 0.00D+00
  g = 1.00D+00
  b = 1.00D+00

  call rgb_to_hue ( r, g, b, h )
  write ( *, '(4f8.4)' ) r, g, b, h

  r = 0.00D+00
  g = 0.00D+00
  b = 1.00D+00

  call rgb_to_hue ( r, g, b, h )
  write ( *, '(4f8.4)' ) r, g, b, h

  r = 1.00D+00
  g = 0.00D+00
  b = 1.00D+00

  call rgb_to_hue ( r, g, b, h )
  write ( *, '(4f8.4)' ) r, g, b, h

  r = 0.00D+00
  g = 0.00D+00
  b = 0.00D+00

  call rgb_to_hue ( r, g, b, h )
  write ( *, '(4f8.4)' ) r, g, b, h

  r = 0.50D+00
  g = 0.50D+00
  b = 0.50D+00

  call rgb_to_hue ( r, g, b, h )
  write ( *, '(4f8.4)' ) r, g, b, h

  r = 1.00D+00
  g = 1.00D+00
  b = 1.00D+00

  call rgb_to_hue ( r, g, b, h )
  write ( *, '(4f8.4)' ) r, g, b, h

  r = 0.94D+00
  g = 0.70D+00
  b = 0.15D+00

  call rgb_to_hue ( r, g, b, h )
  write ( *, '(4f8.4)' ) r, g, b, h

  r = 0.24D+00
  g = 0.70D+00
  b = 0.85D+00

  call rgb_to_hue ( r, g, b, h )
  write ( *, '(4f8.4)' ) r, g, b, h

  r = 0.24D+00
  g = 0.24D+00
  b = 0.85D+00

  call rgb_to_hue ( r, g, b, h )
  write ( *, '(4f8.4)' ) r, g, b, h

  return
end
subroutine test16

!*****************************************************************************80
!
!! TEST16 tests RGB_TO_RGBPRIME, RGBPRIME_TO_RGB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) bprime
  real ( kind = 8 ) g
  real ( kind = 8 ) g2
  real ( kind = 8 ) gprime
  integer itest
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) rprime

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  RGB_TO_RGBPRIME: RGB => R''G''B'' colors;'
  write ( *, '(a)' ) '  RGBPRIME_TO_RGB: R''G''B'' => RGB colors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Rin     Gin     Bin     R''    G''      B''   ' &
    // '   Rout    Gout    Bout'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call rgb_test ( itest, r, g, b )

    if ( r < 0.0D+00 ) then
      exit
    end if

    call rgb_to_rgbprime ( r, g, b, rprime, gprime, bprime )

    call rgbprime_to_rgb ( rprime, gprime, bprime, r2, g2, b2 )

    write ( *, '(9f8.3)' ) r, g, b, rprime, gprime, bprime, r2, g2, b2
 
  end do

  return
end
subroutine test17

!*****************************************************************************80
!
!! TEST17 tests RGB_TO_YCBCR, YCBCR_TO_RGB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) cb
  real ( kind = 8 ) cr
  real ( kind = 8 ) g
  real ( kind = 8 ) g2
  integer itest
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) yb
  real ( kind = 8 ) yg
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yr

  yr = 0.299D+00
  yg = 0.587D+00
  yb = 0.114D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  RGB_TO_YCBCR converts RGB to Y''CbCr colors;'
  write ( *, '(a)' ) '  YCBCR_TO_RGB converts Y''CbCr to RGB colors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Rin     Gin     Bin     Yprime  Cb      Cr   ' &
    // '   Rout    Gout    Bout'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call rgb_test ( itest, r, g, b )

    if ( r < 0.0D+00 ) then
      exit
    end if

    call rgb_to_ycbcr ( r, g, b, yr, yg, yb, yprime, cb, cr )

    call ycbcr_to_rgb ( yprime, cb, cr, yr, yg, yb, r2, g2, b2 )

    write ( *, '(9f8.3)' ) r, g, b, yprime, cb, cr, r2, g2, b2
 
  end do

  return
end
subroutine test18

!*****************************************************************************80
!
!! TEST18 tests RGB_TO_YIQ, YIQ_TO_RGB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) g
  real ( kind = 8 ) g2
  real ( kind = 8 ) i
  integer itest
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) yb
  real ( kind = 8 ) yg
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yr

  yr = 0.299D+00
  yg = 0.587D+00
  yb = 0.114D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  RGB_TO_YIQ converts RGB to Y''IQ colors;'
  write ( *, '(a)' ) '  YIQ_TO_RGB converts Y''IQ to RGB colors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Rin     Gin     Bin     Y''      I       Q    ' &
    // '   Rout    Gout    Bout'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call rgb_test ( itest, r, g, b )

    if ( r < 0.0D+00 ) then
      exit
    end if

    call rgb_to_yiq ( r, g, b, yr, yg, yb, yprime, i, q )

    call yiq_to_rgb ( yprime, i, q, yr, yg, yb, r2, g2, b2 )

    write ( *, '(9f8.3)' ) r, g, b, yprime, i, q, r2, g2, b2
 
  end do

  return
end
subroutine test185

!*****************************************************************************80
!
!! TEST185 tests RGB_TO_YCBCR, YCBCR_TO_RGB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) g
  real ( kind = 8 ) g2
  integer itest
  real ( kind = 8 ) pb
  real ( kind = 8 ) pr
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST185'
  write ( *, '(a)' ) '  RGB_TO_YPBPR converts RGB to YPbPr colors;'
  write ( *, '(a)' ) '  YPBPR_TO_RGB converts YPbPr to RGB colors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Rin     Gin     Bin     Y       Pb      Pr   ' &
    // '   Rout    Gout    Bout'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call rgb_test ( itest, r, g, b )

    if ( r < 0.0D+00 ) then
      exit
    end if

    call rgb_to_ypbpr ( r, g, b, y, pb, pr )

    call ypbpr_to_rgb ( y, pb, pr, r2, g2, b2 )

    write ( *, '(9f8.3)' ) r, g, b, y, pb, pr, r2, g2, b2
 
  end do

  return
end
subroutine test19

!*****************************************************************************80
!
!! TEST19 tests RGB_TO_YUV, YUV_TO_RGB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) g
  real ( kind = 8 ) g2
  integer itest
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) yb
  real ( kind = 8 ) yg
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yr

  yr = 0.299D+00
  yg = 0.587D+00
  yb = 0.114D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  RGB_TO_YUV converts RGB to Y''UV colors;'
  write ( *, '(a)' ) '  YUV_TO_RGB converts Y''UV to RGB colors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Rin     Gin     Bin     Y''      U       V    ' &
    // '   Rout    Gout    Bout'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call rgb_test ( itest, r, g, b )

    if ( r < 0.0D+00 ) then
      exit
    end if

    call rgb_to_yuv ( r, g, b, yr, yg, yb, yprime, u, v )

    call yuv_to_rgb ( yprime, u, v, yr, yg, yb, r2, g2, b2 )

    write ( *, '(9f8.3)' ) r, g, b, yprime, u, v, r2, g2, b2
 
  end do

  return
end
subroutine test20

!*****************************************************************************80
!
!! TEST20 tests RGB709_TO_XYZ_CAP, XYZ_CAP_TO_RGB709.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) g
  real ( kind = 8 ) g2
  integer itest
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) xcap
  real ( kind = 8 ) ycap
  real ( kind = 8 ) zcap

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  RGB709_TO_XYZ_CAP: RGB709 => CIE XYZ colors;'
  write ( *, '(a)' ) '  XYZ_CAP_TO_RGB709: CIE XYZ => RGB709 colors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Rin     Gin     Bin     X       Y       Z    ' &
    // '   Rout    Gout    Bout'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1

    call rgb_test ( itest, r, g, b )

    if ( r < 0.0D+00 ) then
      exit
    end if

    call rgb709_to_xyz_cap ( r, g, b, xcap, ycap, zcap )

    call xyz_cap_to_rgb709 ( xcap, ycap, zcap, r2, g2, b2 )

    write ( *, '(9f8.3)' ) r, g, b, xcap, ycap, zcap, r2, g2, b2
 
  end do

  return
end
subroutine test21

!*****************************************************************************80
!
!! TEST21 tests RGBCIE_TO_XYZ_CAP, XYZ_CAP_TO_RGBCIE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  real ( kind = 8 ) g
  real ( kind = 8 ) g2
  integer itest
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  real ( kind = 8 ) xcap
  real ( kind = 8 ) ycap
  real ( kind = 8 ) zcap

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  RGBCIE_TO_XYZ_CAP converts CIE RGB to XYZ colors;'
  write ( *, '(a)' ) '  XYZ_CAP_TO_RGBCIE converts XYZ to CIE RGB colors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Rin     Gin     Bin     X       Y       Z    ' &
    // '   Rout    Gout    Bout'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call rgb_test ( itest, r, g, b )

    if ( r < 0.0D+00 ) then
      exit
    end if

    call rgbcie_to_xyz_cap ( r, g, b, xcap, ycap, zcap )

    call xyz_cap_to_rgbcie ( xcap, ycap, zcap, r2, g2, b2 )

    write ( *, '(9f8.3)' ) r, g, b, xcap, ycap, zcap, r2, g2, b2
 
  end do

  return
end
subroutine test22

!*****************************************************************************80
!
!! TEST22 tests SRGB_TO_XYZ_CAP, XYZ_CAP_TO_SRGB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  integer itest
  real ( kind = 8 ) r
  integer sb
  integer sb2
  integer sg
  integer sg2
  integer sr
  integer sr2
  real ( kind = 8 ) xcap
  real ( kind = 8 ) xcap2
  real ( kind = 8 ) ycap
  real ( kind = 8 ) ycap2
  real ( kind = 8 ) zcap
  real ( kind = 8 ) zcap2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  SRGB_TO_XYZ_CAP converts sRGB to XYZ colors;'
  write ( *, '(a)' ) '  XYZ_CAP_TO_SRGBCIE converts XYZ to sRGB colors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '   Xin     Yin     Zin    sR  sG  sB   Xout    Yout    Zout'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call rgb_test ( itest, r, g, b )

    if ( r < 0.0D+00 ) then
      exit
    end if

    call rgbcie_to_xyz_cap ( r, g, b, xcap, ycap, zcap )

    call xyz_cap_to_srgb ( xcap, ycap, zcap, sr, sg, sb )

    call srgb_to_xyz_cap ( sr, sg, sb, xcap2, ycap2, zcap2 )

    write ( *, '(3f8.3,3i4,3f8.3)' ) &
      xcap, ycap, zcap, sr, sg, sb, xcap2, ycap2, zcap2
 
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Repeat test, but study sRGB->XYZ->sRGB'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  sR  sG  sB   X       Y       Z     sR2 sG2 sB2'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call rgb_test ( itest, r, g, b )

    if ( r < 0.0D+00 ) then
      exit
    end if

    call rgbcie_to_xyz_cap ( r, g, b, xcap, ycap, zcap )

    call xyz_cap_to_srgb ( xcap, ycap, zcap, sr, sg, sb )

    call srgb_to_xyz_cap ( sr, sg, sb, xcap, ycap, zcap )

    call xyz_cap_to_srgb ( xcap, ycap, zcap, sr2, sg2, sb2 )

    write ( *, '(3i4,3f8.3,3i4)' ) &
      sr, sg, sb, xcap, ycap, zcap, sr2, sg2, sb2
 
  end do
  return
end
subroutine test23

!*****************************************************************************80
!
!! TEST23 tests T_TO_SPD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  integer ihi
  integer j
  real ( kind = 8 ) lambda
  real ( kind = 8 ) power
  real ( kind = 8 ) t
  real ( kind = 8 ) thi
  real ( kind = 8 ) tlo

  ihi = 10D+00
  tlo = 1000.0D+00
  thi = 10000.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  T_TO_SPD evaluates the black body spectral'
  write ( *, '(a)' ) '  power distribution function SPD(T,LAMBDA).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   T        Lambda      SPD(T,LAMBDA)'
  write ( *, '(a)' ) ' '

  do i = 1, ihi

    t = ( real ( ihi - i, kind = 8 ) * tlo &
        + real (       i, kind = 8 ) * thi ) &
        / real ( ihi,     kind = 8 )

    write ( *, * ) ' '

    do j = 380, 780, 40

      lambda = real ( j, kind = 8 )
    
      call t_to_spd ( t, lambda, power )

      write ( *, '(7g12.5)' ) t, lambda, power

    end do

  end do

  return
end
subroutine test24

!*****************************************************************************80
!
!! TEST24 tests T_TO_XY.
!
!  Discussion:
!
!    Thanks to Harald Anlauf, of the Technical University of Darmstadt,
!    for pointing out an error in an output format, 30 April 2002.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  integer ihi
  real ( kind = 8 ) t
  real ( kind = 8 ) thi
  real ( kind = 8 ) tlo
  real ( kind = 8 ) x
  real ( kind = 8 ) xcap
  real ( kind = 8 ) y
  real ( kind = 8 ) ycap
  real ( kind = 8 ) z
  real ( kind = 8 ) zcap

  ihi = 20
  tlo = 1000.0D+00
  thi = 1400.0D+00
  ycap = 100.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  T_TO_XYZ returns the CIE xyz chromaticities of'
  write ( *, '(a)' ) '  a black body radiator at temperature T.'
  write ( *, '(a,g14.6)' ) '  Assume constant luminosity YCAP = ', ycap
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     T         x         y         z        X' &
    // '         Y          Z'
  write ( *, '(a)' ) ' '

  do i = -2, ihi

    t = ( real ( ihi - i, kind = 8 ) * tlo &
        + real (       i, kind = 8 ) * thi ) &
        / real ( ihi,     kind = 8 )

    call t_to_xyz ( t, x, y, z )

    call xyy_to_xyz_cap ( x, y, xcap, ycap, zcap )

    write ( *, '(7f10.4)' ) t, x, y, z, xcap, ycap, zcap

  end do

  return
end
subroutine test25

!*****************************************************************************80
!
!! TEST25 tests UVPRIMEY_TO_XYZ_CAP, XYZ_CAP_TO_UVWPRIME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ntest = 17

  integer itest
  real ( kind = 8 ) nm
  real ( kind = 8 ), parameter :: nmhi = 700.0D+00
  real ( kind = 8 ), parameter :: nmlo = 380.0D+00
  real ( kind = 8 ) temp
  real ( kind = 8 ) uprime
  real ( kind = 8 ) vprime
  real ( kind = 8 ) wprime
  real ( kind = 8 ) x
  real ( kind = 8 ) xcap
  real ( kind = 8 ) xcap2
  real ( kind = 8 ) xcapn
  real ( kind = 8 ) xn
  real ( kind = 8 ) y
  real ( kind = 8 ) ycap
  real ( kind = 8 ) ycapn
  real ( kind = 8 ) yn
  real ( kind = 8 ) z
  real ( kind = 8 ) zcap
  real ( kind = 8 ) zcap2
  real ( kind = 8 ) zcapn
  real ( kind = 8 ) zn

  ycapn = 1.0D+00
  call name_to_xyz ( 'D65', xn, yn, zn )
  call xyy_to_xyz_cap ( xn, yn, xcapn, ycapn, zcapn )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) '  UVPRIMEY_TO_XYZ_CAP converts u''v''Y to XYZ colors.'
  write ( *, '(a)' ) '  XYZ_CAP_TO_UVWPRIME converts XYZ to u''v''w'' colors;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Illuminant XYZ color coordinates:'
  write ( *, '(3g14.6)' ) xcapn, ycapn, zcapn 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Xin     Yin     Zin       u''       v''' &
    // '       w''      Xout   Yout   Zout'
  write ( *, '(a)' ) ' '

  do itest = 1, ntest

    nm = ( real ( ntest - itest,     kind = 8 ) * nmlo &
         + real (         itest - 1, kind = 8 ) * nmhi ) &
         / real ( ntest         - 1, kind = 8 )

    call nm_to_xyz ( nm, x, y, z )
    ycap = 0.95D+00 * ycapn
    call xyy_to_xyz_cap ( x, y, xcap, ycap, zcap )
!
!  Make sure every component of (XCAP,YCAP,ZCAP) is smaller than
!  (XCAPN,YCAPN,ZCAPN).
!
    temp = 1.0D+00
    if ( 0.0D+00 < xcap .and. 0.0D+00 < xcapn ) then
      temp = min ( temp, 0.95D+00 * xcapn / xcap )
    end if
    if ( 0.0D+00 < ycap .and. 0.0D+00 < ycapn ) then
      temp = min ( temp, 0.95D+00 * ycapn / ycap )
    end if
    if ( 0.0D+00 < zcap .and. 0.0D+00 < zcapn ) then
      temp = min ( temp, 0.95D+00 * zcapn / zcap )
    end if

    xcap = xcap * temp
    ycap = ycap * temp
    zcap = zcap * temp

    call xyz_cap_to_uvwprime ( xcap, ycap, zcap, uprime, vprime, wprime )

    call uvprimey_to_xyz_cap ( uprime, vprime, xcap2, ycap, zcap2 )

    write ( *, '(3f8.3,3f9.3,3f8.3)' ) &
      xcap, ycap, zcap, uprime, vprime, wprime, xcap2, ycap, zcap2
 
  end do

  return
end
subroutine test26

!*****************************************************************************80
!
!! TEST26 tests XYZ_CAP_TO_YCC, YCC_TO_XYZ_CAP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ntest = 17

  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  integer itest
  real ( kind = 8 ) nm
  real ( kind = 8 ), parameter :: nmhi = 700.0D+00
  real ( kind = 8 ), parameter :: nmlo = 380.0D+00
  real ( kind = 8 ) temp
  real ( kind = 8 ) x
  real ( kind = 8 ) xcap
  real ( kind = 8 ) xcap2
  real ( kind = 8 ) xcapn
  real ( kind = 8 ) xn
  real ( kind = 8 ) y
  real ( kind = 8 ) yb
  real ( kind = 8 ) ycap
  real ( kind = 8 ) ycap2
  real ( kind = 8 ) ycapn
  real ( kind = 8 ) yg
  real ( kind = 8 ) yn
  real ( kind = 8 ) yr
  real ( kind = 8 ) yval
  real ( kind = 8 ) z
  real ( kind = 8 ) zcap
  real ( kind = 8 ) zcap2
  real ( kind = 8 ) zcapn
  real ( kind = 8 ) zn

  yr = 0.299D+00
  yg = 0.587D+00
  yb = 0.114D+00

  ycapn = 1.0D+00
  call name_to_xyz ( 'D65', xn, yn, zn )
  call xyy_to_xyz_cap ( xn, yn, xcapn, ycapn, zcapn )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) '  XYZ_CAP_TO_YCC converts XYZ to PhotoYCC colors;'
  write ( *, '(a)' ) '  YCC_TO_XYZ_CAP converts PhotoYCC to XYZ colors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Illuminant XYZ color coordinates:'
  write ( *, '(3g14.6)' ) xcapn, ycapn, zcapn
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' NM     Xin     Yin     Zin     Y       C1' &
    // '      C2      Xout    Yout    Zout'
  write ( *, '(a)' ) ' '

  do itest = 1, ntest

    nm = ( real ( ntest - itest,     kind = 8 ) * nmlo &
         + real (         itest - 1, kind = 8 ) * nmhi ) &
         / real ( ntest         - 1, kind = 8 )

    call nm_to_xyz ( nm, x, y, z )
    ycap = 0.95D+00 * ycapn
    call xyy_to_xyz_cap ( x, y, xcap, ycap, zcap )
!
!  Make sure every component of (XCAP,YCAP,ZCAP) is smaller than
!  (XCAPN,YCAPN,ZCAPN).
!
    temp = 1.0D+00
    if ( 0.0D+00 < xcap .and. 0.0D+00 < xcapn ) then
      temp = min ( temp, 0.95D+00 * xcapn / xcap )
    end if
    if ( 0.0D+00 < ycap .and. 0.0D+00 < ycapn ) then
      temp = min ( temp, 0.95D+00 * ycapn / ycap )
    end if
    if ( 0.0D+00 < zcap .and. 0.0D+00 < zcapn ) then
      temp = min ( temp, 0.95D+00 * zcapn / zcap )
    end if

    xcap = xcap * temp
    ycap = ycap * temp
    zcap = zcap * temp

    call xyz_cap_to_ycc ( xcap, ycap, zcap, yr, yg, yb, yval, c1, c2 )

    call ycc_to_xyz_cap ( yval, c1, c2, yr, yg, yb, xcap2, ycap2, zcap2 )

    write ( *, '(f5.1,9f8.3)' ) nm, xcap, ycap, zcap, yval, c1, c2, xcap2, &
      ycap2, zcap2
 
  end do

  return
end
subroutine test27

!*****************************************************************************80
!
!! TEST27 tests YCBCR_TO_YCC, YCC_TO_YCBCR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) c1in
  real ( kind = 8 ) c1out
  real ( kind = 8 ) c2in
  real ( kind = 8 ) c2out
  real ( kind = 8 ) cb
  real ( kind = 8 ) cr
  integer itest
  real ( kind = 8 ) y
  real ( kind = 8 ) yin
  real ( kind = 8 ) yout

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)' ) '  YCBCR_TO_YCC converts YCBCR to PhotoYCC colors.'
  write ( *, '(a)' ) '  YCC_TO_YCBCR converts PhotoYCC to YCBCR colors;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Yin    C1in     C2in     Y''     Cb' &
    // '      Cr      Yout    C1out    C2out'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1
    call ycc_test ( itest, yin, c1in, c2in )

    if ( yin < 0.0D+00 ) then
      exit
    end if

    call ycc_to_ycbcr ( yin, c1in, c2in, y, cb, cr )

    call ycbcr_to_ycc ( y, cb, cr, yout, c1out, c2out )

    write ( *, '(9f8.3)' ) yin, c1in, c2in, y, cb, cr, yout, c1out, c2out
 
  end do

  return
end
subroutine test28

!*****************************************************************************80
!
!! TEST28 tests YIQ_TO_YUV, YUV_TO_YIQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ) i
  real ( kind = 8 ) i2
  integer itest
  real ( kind = 8 ) q
  real ( kind = 8 ) q2
  real ( kind = 8 ) r
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) yb
  real ( kind = 8 ) yg
  real ( kind = 8 ) yprime
  real ( kind = 8 ) yprime2
  real ( kind = 8 ) yprime3
  real ( kind = 8 ) yr

  yr = 0.299D+00
  yg = 0.587D+00
  yb = 0.114D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST28'
  write ( *, '(a)' ) '  YIQ_TO_YUV converts Y''IQ to Y''UV colors;'
  write ( *, '(a)' ) '  YUV_TO_YIQ converts Y''UV to Y''IQ colors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Y''in    Iin     Qin     Y''      U' &
    // '       V       Y''out   Iout    Qout'
  write ( *, '(a)' ) ' '

  itest = 0

  do

    itest = itest + 1

    call rgb_test ( itest, r, g, b )

    if ( r < 0.0D+00 ) then
      exit
    end if

    call rgb_to_yiq ( r, g, b, yr, yg, yb, yprime, i, q )

    call yiq_to_yuv ( yprime, i, q, yprime2, u, v )

    call yuv_to_yiq ( yprime2, u, v, yprime3, i2, q2 )

    write ( *, '(9f8.3)' ) yprime, i, q, yprime2, u, v, yprime3, i2, q2
 
  end do

  return
end
