program main

!*****************************************************************************80
!
!! MAIN is the main program for BLEND_PRB.
!
!  Discussion:
!
!    BLEND_PRB tests routines from BLEND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BLEND_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the BLEND library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )
  call test10 ( )
  call test11 ( )
  call test12 ( )
  call test13 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BLEND_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 checks for a gross error in the blend coefficients.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 3

  external identity_r
  external identity_rs
  external identity_rst
  integer ( kind = 4 ) n
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x(nmax)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Simple identity test to detect gross errors.'
!
!  Set N to 1.
!
  n = 1
!
!  Test BLEND_R_0DN on identity.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Identity test for BLEND_R_0DN:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  call blend_r_0dn ( r, x, n, identity_r )
  write ( *, '(2f8.4)' ) r, x(1:n)

  r = 1.0D+00
  call blend_r_0dn ( r, x, n, identity_r )
  write ( *, '(2f8.4)' ) r, x(1:n)

  r = 0.5D+00
  call blend_r_0dn ( r, x, n, identity_r )
  write ( *, '(2f8.4)' ) r, x(1:n)
!
!  Set N to 2.
!
  n = 2
!
!  Test BLEND_RS_0DN on identity.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Identity test for BLEND_RS_0DN:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  call blend_rs_0dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  call blend_rs_0dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  call blend_rs_0dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  call blend_rs_0dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  call blend_rs_0dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)
!
!  Test BLEND_RS_1DN on identity.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Identity test for BLEND_RS_1DN:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  call blend_rs_1dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  call blend_rs_1dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  call blend_rs_1dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  call blend_rs_1dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  call blend_rs_1dn ( r, s, x, n, identity_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)
!
!  Set N to 3.
!
  n = 3
!
!  Test BLEND_RST_0DN on identity.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Identity test for BLEND_RST_0DN:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  t = 0.0D+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 0.0D+00
  t = 1.0D+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  t = 1.0D+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  t = 0.5D+00
  call blend_rst_0dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)
!
!  Test BLEND_RST_1DN on identity.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Identity test for BLEND_RST_1DN:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  t = 0.0D+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 0.0D+00
  t = 1.0D+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  t = 1.0D+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  t = 0.5D+00
  call blend_rst_1dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)
!
!  Test BLEND_RST_2DN on identity.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Identity test for BLEND_RST_2DN:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  t = 0.0D+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 0.0D+00
  t = 1.0D+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  t = 1.0D+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  t = 0.5D+00
  call blend_rst_2dn ( r, s, t, x, n, identity_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 checks for simple errors in the blend coefficients.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 3

  integer ( kind = 4 ) n
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  external stretch_r
  external stretch_rs
  external stretch_rst
  real ( kind = 8 ) t
  real ( kind = 8 ) x(nmax)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Shift and stretch test to detect simple errors.'
!
!  Set N to 1.
!
  n = 1
!
!  Test BLEND_R_0DN on shift by 1, stretch by 2.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Shift and stretch test for BLEND_R_0DN:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  call blend_r_0dn ( r, x, n, stretch_r )
  write ( *, '(2f8.4)' ) r, x(1:n)

  r = 1.0D+00
  call blend_r_0dn ( r, x, n, stretch_r )
  write ( *, '(2f8.4)' ) r, x(1:n)

  r = 0.5D+00
  call blend_r_0dn ( r, x, n, stretch_r )
  write ( *, '(2f8.4)' ) r, x(1:n)
!
!  Set N to 2.
!
  n = 2
!
!  Test BLEND_RS_0DN on shift by (1,2), stretch by (3,4).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Shift and stretch test for BLEND_RS_0DN:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  call blend_rs_0dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  call blend_rs_0dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  call blend_rs_0dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  call blend_rs_0dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  call blend_rs_0dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)
!
!  Test BLEND_RS_1D on shift by (1,2), stretch by (3,4).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Shift and stretch test for BLEND_RS_1DN:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  call blend_rs_1dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  call blend_rs_1dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  call blend_rs_1dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  call blend_rs_1dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  call blend_rs_1dn ( r, s, x, n, stretch_rs )
  write ( *, '(4f8.4)' ) r, s, x(1:n)
!
!  Set N to 3.
!
  n = 3
!
!  Test BLEND_RST_0DN on shift by (1,2,3), stretch by (4,5,6).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Shift and stretch test for BLEND_RST_0DN:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  t = 0.0D+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 0.0D+00
  t = 1.0D+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  t = 1.0D+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  t = 0.5D+00
  call blend_rst_0dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)
!
!  Test BLEND_RST_1DN on shift by (1,2,3), stretch by (4,5,6).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Shift and stretch test for BLEND_RST_1DN:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  t = 0.0D+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 0.0D+00
  t = 1.0D+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  t = 1.0D+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)
  r = 0.5D+00
  s = 0.5D+00
  t = 0.5D+00
  call blend_rst_1dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)
!
!  Test BLEND_RST_2DN on shift by (1,2,3), stretch by (4,5,6).
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Shift and stretch test for BLEND_RST_2DN:'
  write ( *, '(a)' ) ' '
  r = 0.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 0.0D+00
  t = 0.0D+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 1.0D+00
  t = 0.0D+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.0D+00
  s = 0.0D+00
  t = 1.0D+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 1.0D+00
  s = 1.0D+00
  t = 1.0D+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  r = 0.5D+00
  s = 0.5D+00
  t = 0.5D+00
  call blend_rst_2dn ( r, s, t, x, n, stretch_rst )
  write ( *, '(6f8.4)' ) r, s, t, x(1:n)

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 checks out BLEND_I_0D1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(m)

  x(1) = 100.0D+00
  x(m) = 100.0 + real ( ( m - 1 ) * 5, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  BLEND_I_0D1 interpolates data in a vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  X(1) = ', x(1)
  write ( *, '(a,g14.6)' ) '  X(', m, ')= ', x(m)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Interpolated values:'
  write ( *, '(a)' ) ' '

  call blend_i_0d1 ( x, m )

  do i = 1, m
    write ( *, '(2x,i8,2x,g14.6)' ) i, x(i)
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 checks out BLEND_IJ_0D1 and BLEND_IJ_1D1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m1 = 5
  integer ( kind = 4 ), parameter :: m2 = 4

  external cubic_rs
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) x(m1,m2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  BLEND_IJ_0D1 interpolates data in a table,'
  write ( *, '(a)' ) '  from corner data.'
  write ( *, '(a)' ) '  BLEND_IJ_1D1 interpolates data in a table,'
  write ( *, '(a)' ) '  from edge data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8,a)' ) '  The table is ', m1, ' rows by ', &
    m2, ' columns.'
!
!  Load data in the corners only.
!
  i = 1
  j = 1
  r = real (  i - 1, kind = 8 ) &
    / real ( m1 - 1, kind = 8 )
  s = real (  j - 1, kind = 8 ) &
    / real ( m2 - 1, kind = 8 )
  call cubic_rs ( r, s, 1, x(i,j) )

  i = m1
  j = 1
  r = real (  i - 1, kind = 8 ) &
    / real ( m1 - 1, kind = 8 )
  s = real (  j - 1, kind = 8 ) &
    / real ( m2 - 1, kind = 8 )
  call cubic_rs ( r, s, 1, x(i,j) )

  i = 1
  j = m2
  r = real (  i - 1, kind = 8 ) &
    / real ( m1 - 1, kind = 8 )
  s = real (  j - 1, kind = 8 ) &
    / real ( m2 - 1, kind = 8 )
  call cubic_rs ( r, s, 1, x(i,j) )

  i = m1
  j = m2
  r = real (  i - 1, kind = 8 ) &
    / real ( m1 - 1, kind = 8 )
  s = real (  j - 1, kind = 8 ) &
    / real ( m2 - 1, kind = 8 )
  call cubic_rs ( r, s, 1, x(i,j) )

  call blend_ij_0d1 ( x, m1, m2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Values interpolated by BLEND_IJ_0D1:'
  write ( *, '(a)' ) ' '

  do i = 1, m1
    write ( *, '(5g14.6)' ) x(i,1:m2)
  end do
!
!  Load data in the edges.
!
  j = 1
  do i = 1, m1
    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )
    s = real (  j - 1, kind = 8 ) &
      / real ( m2 - 1, kind = 8 )
    call cubic_rs ( r, s, 1, x(i,j) )
  end do

  j = m2
  do i = 1, m1
    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )
    s = real (  j - 1, kind = 8 ) &
      / real ( m2 - 1, kind = 8 )
    call cubic_rs ( r, s, 1, x(i,j) )
  end do

  i = 1
  do j = 2, m2 - 1
    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )
    s = real (  j - 1, kind = 8 ) &
      / real ( m2 - 1, kind = 8 )
    call cubic_rs ( r, s, 1, x(i,j) )
  end do

  i = m1
  do j = 2, m2 - 1
    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )
    s = real (  j - 1, kind = 8 ) &
      / real ( m2 - 1, kind = 8 )
    call cubic_rs ( r, s, 1, x(i,j) )
  end do

  call blend_ij_1d1 ( x, m1, m2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Values interpolated by BLEND_IJ_1D1:'
  write ( *, '(a)' ) ' '

  do i = 1, m1
    write ( *, '(5g14.6)' ) x(i,1:m2)
  end do
!
!  Compare with BLEND_RS_1D1
!
  do i = 1, m1
    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )
    do j = 1, m2
      s = real (  j - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )
      call blend_rs_1dn ( r, s, x(i,j), 1, cubic_rs )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data blended by BLEND_RS_1DN:'
  write ( *, '(a)' ) ' '

  do i = 1, m1
    write ( *, '(5g14.6)' ) x(i,1:m2)
  end do
!
!  Load all data.
!
  do i = 1, m1
    do j = 1, m2
      r = real (  i - 1, kind = 8 ) &
        / real ( m1 - 1, kind = 8 )
      s = real (  j - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )
      call cubic_rs ( r, s, 1, x(i,j) )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact data:'
  write ( *, '(a)' ) ' '

  do i = 1, m1
    write ( *, '(5g14.6)' ) x(i,1:m2)
  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 checks out BLEND_IJK_0D1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m1 = 4
  integer ( kind = 4 ), parameter :: m2 = 3
  integer ( kind = 4 ), parameter :: m3 = 3

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) num_extreme
  external quad_rst
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x(m1,m2,m3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  BLEND_IJK_0D1 interpolates data in a table,'
  write ( *, '(a)' ) '  from corner data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8,a,i8,a)' ) '  The table is ', m1, ' rows by ', &
    m2, ' columns by ', m3, ' layers.'
!
!  Load data on the faces.
!
  do i = 1, m1
    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )
    do j = 1, m2
      s = real (  j - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )
      do k = 1, m3
        t = real (  k - 1, kind = 8 ) &
          / real ( m3 - 1, kind = 8 )

        num_extreme = 0
        if ( i == 1 .or. i == m1 ) then
          num_extreme = num_extreme + 1
        end if
        if ( j == 1 .or. j == m2 ) then
          num_extreme = num_extreme + 1
        end if
        if ( k == 1 .or. k == m3 ) then
          num_extreme = num_extreme + 1
        end if

        if ( num_extreme == 3 ) then
          call quad_rst ( r, s, t, 1, x(i,j,k) )
        else
          x(i,j,k) = 0.0D+00
        end if

      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data given to BLEND_IJK_0D1:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  call blend_ijk_0d1 ( x, m1, m2, m3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Values interpolated by BLEND_IJK_0D1:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do
!
!  Load all data.
!
  do i = 1, m1
    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )
    do j = 1, m2
      s = real (  j - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )
      do k = 1, m3
        t = real (  k - 1, kind = 8 ) &
          / real ( m3 - 1, kind = 8 )
        call quad_rst ( r, s, t, 1, x(i,j,k) )
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact data:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 checks out BLEND_IJK_1D1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m1 = 4
  integer ( kind = 4 ), parameter :: m2 = 3
  integer ( kind = 4 ), parameter :: m3 = 3

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) num_extreme
  external quad_rst
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x(m1,m2,m3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  BLEND_IJK_1D1 interpolates data in a table,'
  write ( *, '(a)' ) '  from edge data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8,a,i8,a)' ) '  The table is ', m1, ' rows by ', &
    m2, ' columns by ', m3, ' layers.'
!
!  Load data on the faces.
!
  do i = 1, m1
    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )
    do j = 1, m2
      s = real (  j - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )
      do k = 1, m3
        t = real (  k - 1, kind = 8 ) &
          / real ( m3 - 1, kind = 8 )

        num_extreme = 0
        if ( i == 1 .or. i == m1 ) then
          num_extreme = num_extreme + 1
        end if
        if ( j == 1 .or. j == m2 ) then
          num_extreme = num_extreme + 1
        end if
        if ( k == 1 .or. k == m3 ) then
          num_extreme = num_extreme + 1
        end if

        if ( 2 <= num_extreme ) then
          call quad_rst ( r, s, t, 1, x(i,j,k) )
        else
          x(i,j,k) = 0.0D+00
        end if

      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data given to BLEND_IJK_1D1:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  call blend_ijk_1d1 ( x, m1, m2, m3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Values interpolated by BLEND_IJK_1D1:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do
!
!  Load all data.
!
  do i = 1, m1
    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )
    do j = 1, m2
      s = real (  j - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )
      do k = 1, m3
        t = real (  k - 1, kind = 8 ) &
          / real ( m3 - 1, kind = 8 )
        call quad_rst ( r, s, t, 1, x(i,j,k) )
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact data:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 checks out BLEND_IJK_2D1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m1 = 4
  integer ( kind = 4 ), parameter :: m2 = 3
  integer ( kind = 4 ), parameter :: m3 = 3

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) num_extreme
  external quad_rst
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x(m1,m2,m3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  BLEND_IJK_2D1 interpolates data in a table,'
  write ( *, '(a)' ) '  from face data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8,a,i8,a)' ) '  The table is ', m1, ' rows by ', &
    m2, ' columns by ', &
    m3, ' layers.'
!
!  Load data on the faces.
!
  do i = 1, m1
    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )
    do j = 1, m2
      s = real (  j - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )
      do k = 1, m3
        t = real (  k - 1, kind = 8 ) &
          / real ( m3 - 1, kind = 8 )

        num_extreme = 0
        if ( i == 1 .or. i == m1 ) then
          num_extreme = num_extreme + 1
        end if
        if ( j == 1 .or. j == m2 ) then
          num_extreme = num_extreme + 1
        end if
        if ( k == 1 .or. k == m3 ) then
          num_extreme = num_extreme + 1
        end if

        if ( 1 <= num_extreme ) then
          call quad_rst ( r, s, t, 1, x(i,j,k) )
        else
          x(i,j,k) = 0.0D+00
        end if

      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data given to BLEND_IJK_2D1:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  call blend_ijk_2d1 ( x, m1, m2, m3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Values interpolated by BLEND_IJK_2D1:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do
!
!  Load all data.
!
  do i = 1, m1
    r = real (  i - 1, kind = 8 ) &
      / real ( m1 - 1, kind = 8 )
    do j = 1, m2
      s = real (  j - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )
      do k = 1, m3
        t = real (  k - 1, kind = 8 ) &
          / real ( m3 - 1, kind = 8 )
        call quad_rst ( r, s, t, 1, x(i,j,k) )
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact data:'
  write ( *, '(a)' ) ' '

  do k = 1, m3
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Layer K = ', k
    write ( *, '(a)' ) ' '
    do i = 1, m1
      write ( *, '(5g14.6)' ) x(i,1:m2,k)
    end do
  end do

  return
end
subroutine cubic_rs ( r, s, i, x )

!*****************************************************************************80
!
!! CUBIC_RS evaluates a function of R and S used for some tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, the (local) coordinates of a point.
!
!    Input, integer I, the component of X to be returned.
!
!    Output, real X, the value of the I-th component of X at the point whose
!    local coordinates are (R,S).
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) x

  x = 20.0D+00 * ( r**2 * s**3 )

  return
end
subroutine quad_rst ( r, s, t, i, x )

!*****************************************************************************80
!
!! QUAD_RST evaluates a function of (R,S,T) used for some tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, T, the (local) coordinates of a point.
!
!    Input, integer I, the component of X to be returned.
!
!    Output, real X, the value of the I-th component of X at the point whose
!    local coordinates are (R,S,T).
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x

  x = 18.0D+00 * ( r**2 + s + t )

  return
end
subroutine identity_r ( r, i, x )

!*****************************************************************************80
!
!! IDENTITY_R returns a data component given (R).
!
!  Discussion:
!
!    This is a dummy routine, which simply returns (R).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the coordinate of a point that lies on the
!    boundary of the cube.
!
!    Input, integer I, the component of the data which is to be returned.
!
!    Output, real X, the I-th component of the data vector X, evaluated
!    at the point (R), which is on an endpoint of the unit line segment.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = r
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IDENTITY_R - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal component index I = ', i
    stop
  end if

  return
end
subroutine identity_rs ( r, s, i, x )

!*****************************************************************************80
!
!! IDENTITY_RS returns a data component given (R,S).
!
!  Discussion:
!
!    This is a dummy routine, which simply returns (R,S).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, the coordinates of a point that lies on the
!    boundary of the square.
!
!    Input, integer I, the component of the data which is to be returned.
!
!    Output, real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, or edge, of the unit square.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = r
  else if ( i == 2 ) then
    x = s
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IDENTITY_RS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal component index I = ', i
    stop
  end if

  return
end
subroutine identity_rst ( r, s, t, i, x )

!*****************************************************************************80
!
!! IDENTITY_RST returns a data component given (R,S,T).
!
!  Discussion:
!
!    This is a dummy routine, which simply returns (R,S,T).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, T, the coordinates of a point that lies on the
!    boundary of the cube.
!
!    Input, integer I, the component of the data which is to be returned.
!
!    Output, real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, edge or face of the unit cube.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = r
  else if ( i == 2 ) then
    x = s
  else if ( i == 3 ) then
    x = t
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IDENTITY_RST - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal component index I = ', i
    stop
  end if

  return
end
subroutine stretch_r ( r, i, x )

!*****************************************************************************80
!
!! STRETCH_R returns a data component given (R).
!
!  Discussion:
!
!    This routine shifts by 1 and stretches by 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the coordinate of a point that lies on the
!    boundary of the cube.
!
!    Input, integer ( kind = 4 ) I, the component of the data which is to be returned.
!
!    Output, real X, the I-th component of the data vector X, evaluated
!    at the point (R), which is on an endpoint of the unit line segment.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 2.0D+00 * r + 1.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRETCH_R - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal component index I = ', i
    stop
  end if

  return
end
subroutine stretch_rs ( r, s, i, x )

!*****************************************************************************80
!
!! STRETCH_RS returns a data component given (R,S).
!
!  Discussion:
!
!    This routine shifts by (1,2) and stretches by (3,4).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, the coordinates of a point that lies on the
!    boundary of the square.
!
!    Input, integer ( kind = 4 ) I, the component of the data which is to be returned.
!
!    Output, real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, or edge, of the unit square.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 3.0D+00 * r + 1.0D+00
  else if ( i == 2 ) then
    x = 4.0D+00 * s + 2.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRETCH_RS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal component index I = ', i
    stop
  end if

  return
end
subroutine stretch_rst ( r, s, t, i, x )

!*****************************************************************************80
!
!! STRETCH_RST returns a data component given (R,S,T).
!
!  Discussion:
!
!    This routine shifts by (1,2,3) and stretches by (4,5,6)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, T, the coordinates of a point that lies on the
!    boundary of the cube.
!
!    Input, integer ( kind = 4 ) I, the component of the data which is to be returned.
!
!    Output, real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, edge or face of the unit cube.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x

  if ( i == 1 ) then
    x = 4.0D+00 * r + 1.0D+00
  else if ( i == 2 ) then
    x = 5.0D+00 * s + 2.0D+00
  else if ( i == 3 ) then
    x = 6.0D+00 * t + 3.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STRETCH_RST - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal component index I = ', i
    stop
  end if

  return
end
subroutine ellipse_rs ( r, s, i, x )

!*****************************************************************************80
!
!! ELLIPSE_RS maps the boundary of the unit square to an ellipse.
!
!  Discussion:
!
!    The ellipse is ( 3 * cos ( THETA ), 2 * sin ( THETA ) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, the coordinates of a point that lies on the
!    boundary of the square.
!
!    Input, integer ( kind = 4 ) I, the component of the data which is to be returned.
!
!    Output, real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, or edge, of the unit square.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.14159265D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) theta
  real ( kind = 8 ) x

  if ( r == 0.0D+00 ) then
    theta = 0.25D+00 * pi * ( 5.0D+00 * ( 1.0D+00 - s ) + 3.0D+00 * s  )
  else if ( r == 1.0D+00 ) then
    theta = 0.25D+00 * pi * ( - 1.0D+00 * ( 1.0D+00 - s ) + 1.0D+00 * s )
  else if ( s == 0.0D+00 ) then
    theta = 0.25D+00 * pi * ( 5.0D+00 * ( 1.0D+00 - r ) + 7.0D+00 * r )
  else if ( s == 1.0D+00 ) then
    theta = 0.25D+00 * pi * ( 3.0D+00 * ( 1.0D+00 - r ) + 1.0D+00 * r )
  end if

  if ( i == 1 ) then

    x = 3.0D+00 * cos ( theta )

  else if ( i == 2 ) then

    x = 2.0D+00 * sin ( theta )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ELLIPSE_RS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal component index I = ', i
    stop

  end if

  return
end
subroutine sphere_rst ( r, s, t, i, x )

!*****************************************************************************80
!
!! SPHERE_RST maps the boundary of the unit cube to a sphere.
!
!  Discussion:
!
!    The sphere is
!      x = cos ( theta ) * cos ( phi )
!      y = sin ( theta ) * cos ( phi )
!      z = sin ( phi )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, the coordinates of a point that lies on the
!    boundary of the square.
!
!    Input, integer ( kind = 4 ) I, the component of the data which is to be returned.
!
!    Output, real X, the I-th component of the data vector X, evaluated
!    at the point (R,S), which is on a corner, or edge, of the unit square.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) norm
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x
!
!  Compute length of vector from ( 0.5, 0.5, 0.5 ) to ( r, s, t )
!
  norm = sqrt ( ( r - 0.5D+00 )**2 + ( s - 0.5D+00 )**2 + ( t - 0.5D+00 )**2 )
!
!  Compute ( x, y, z ) coordinates of a point on the sphere
!  ( x - 0.5 )**2 + ( y - 0.5 )**2 + ( z - 0.5 )**2 = 0.25 that is
!  the projection of the point ( r, s, t ) on the unit cube.
!
  if ( i == 1 ) then

    x = 0.5D+00 + 0.5D+00 * ( r - 0.5D+00 ) / norm

  else if ( i == 2 ) then

    x = 0.5D+00 + 0.5D+00 * ( s - 0.5D+00 ) / norm

  else if ( i == 3 ) then

    x = 0.5D+00 + 0.5D+00 * ( t - 0.5D+00 ) / norm

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE_RST - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal component index I = ', i
    stop

  end if

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests BLEND_IJ_W_1D1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 5
  integer ( kind = 4 ), parameter :: n2 = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: pi = 3.14159265D+00
  real ( kind = 8 ) r(n1)
  real ( kind = 8 ) rad
  real ( kind = 8 ) rr
  real ( kind = 8 ) s(n2)
  real ( kind = 8 ) ss
  real ( kind = 8 ) x(n1,n2)
  real ( kind = 8 ) y(n1,n2)

  rad = 3.0D+00

  x(1:n1,1:n2) = 0.0D+00
  y(1:n1,1:n2) = 0.0D+00
!
!  Set the boundary values.
!
!  It turns out that our values correspond to the X and Y
!  coordinates of a quarter circle of radius 3, although
!  it is by no means necessary that a formula for the data
!  be known.
!
  do i = 1, n1
    rr = ( real ( i - 1, kind = 8 ) / real ( n1 - 1, kind = 8 ) )**2
    r(i) = rr
    x(i,1) = 0.0D+00
    y(i,1) = 0.0D+00
    x(i,n2) = rad * cos ( 0.5D+00 * pi * ( 1.0D+00 - rr ) )
    y(i,n2) = rad * sin ( 0.5D+00 * pi * ( 1.0D+00 - rr ) )
  end do

  do j = 1, n2
    ss = ( real ( j - 1, kind = 8 ) / real ( n2 - 1, kind = 8 ) )**2
    s(j) = ss
    x(1,j) = 0.0D+00
    y(1,j) = rad * ss
    x(n1,j) = rad * ss
    y(n1,j) = 0.0D+00
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  BLEND_IJ_W_1D1 uses blending to fill in the'
  write ( *, '(a)' ) '    interior of a table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     R           S           X           Y'
  write ( *, '(a)' ) ' '
  call blend_ij_w_1d1 ( x, r, s, n1, n2 )
  call blend_ij_w_1d1 ( y, r, s, n1, n2 )

  do i = 1, n1
    write ( *, '(a)' ) ' '
    do j = 1, n2
      write ( *, '(4g12.4)' ) r(i), s(j), x(i,j), y(i,j)
    end do
  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests BLEND_102.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 5
  integer ( kind = 4 ), parameter :: n2 = 5

  real ( kind = 8 ) d(n1,n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r
  real ( kind = 8 ) s

  d(1:n1,1:n2) = 0.0D+00

  do i = 1, n1
    do j = 1, n2
      if ( &
        ( i == 1 .or. i == n1 ) .and. &
        ( j == 1 .or. j == n2 ) ) then
        d(i,j) = real ( i + j, kind = 8 )
      end if
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  BLEND_102 blends corner values into a table.'
  write ( *, '(a)' ) ' '

  call r8mat_print ( n1, n2, d, '  Initial data array' )

  do i = 1, n1

    r = real (  i - 1, kind = 8 ) &
      / real ( n1 - 1, kind = 8 )

    do j = 1, n2

      s = real (  j - 1, kind = 8 ) &
        / real ( n2 - 1, kind = 8 )

      if ( &
        ( i == 1 .or. i == n1 ) .and. &
        ( j == 1 .or. j == n2 ) ) then
        cycle
      end if

      call blend_102 ( r, s, d(1,1), d(1,n2), d(n1,1), d(n1,n2), d(i,j) )

    end do

  end do

  call r8mat_print ( n1, n2, d, '  Interpolated data array' )

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests BLEND_112.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 5
  integer ( kind = 4 ), parameter :: n2 = 5

  real ( kind = 8 ) d(n1,n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r
  real ( kind = 8 ) s

  d(1:n1,1:n2) = 0.0D+00

  do i = 1, n1
    do j = 1, n2
      if ( i == 1 .or. i == n1 .or. j == 1 .or. j == n2 ) then
        d(i,j) = real ( i + j, kind = 8 )
      end if
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  BLEND_112 blends side values into a table.'
  write ( *, '(a)' ) ' '

  call r8mat_print ( n1, n2, d, '  Initial data array' )

  do i = 1, n1

    r = real (  i - 1, kind = 8 ) &
      / real ( n1 - 1, kind = 8 )

    do j = 1, n2

      s = real (  j - 1, kind = 8 ) &
        / real ( n2 - 1, kind = 8 )

      if ( i == 1 .or. i == n1 .or. j == 1 .or. j == n2 ) then
        cycle
      end if

      call blend_112 ( r, s, d(1,1), d(1,n2), d(n1,1), d(n1,n2), &
        d(i,1), d(i,n2), d(1,j), d(n1,j), d(i,j) )

    end do

  end do

  call r8mat_print ( n1, n2, d, '  Interpolated data array' )

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests BLEND_103.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 3
  integer ( kind = 4 ), parameter :: n2 = 5
  integer ( kind = 4 ), parameter :: n3 = 4

  real ( kind = 8 ) d(n1,n2,n3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t

  d(1:n1,1:n2,1:n3) = 0.0D+00

  do i = 1, n1
    do j = 1, n2
      do k = 1, n3
        if ( &
          ( i == 1 .or. i == n1 ) .and. &
          ( j == 1 .or. j == n2 ) .and. &
          ( k == 1 .or. k == n3  ) ) then
          d(i,j,k) = real ( i + j + k, kind = 8 )
        end if
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  BLEND_103 blends corner values into a table.'
  write ( *, '(a)' ) ' '

  call r8block_print ( n1, n2, n3, d, '  Initial data array' )

  do i = 1, n1

    r = real (  i - 1, kind = 8 ) &
      / real ( n1 - 1, kind = 8 )

    do j = 1, n2

      s = real (  j - 1, kind = 8 ) &
        / real ( n2 - 1, kind = 8 )

      do k = 1, n3

        t = real (  k - 1, kind = 8 ) &
          / real ( n3 - 1, kind = 8 )

        if ( &
          ( i == 1 .or. i == n1 ) .and. &
          ( j == 1 .or. j == n2 ) .and. &
          ( k == 1 .or. k == n3  ) ) then
          cycle
        end if

        call blend_103 ( r, s, t, &
          d(1,1,1), d(1,1,n3), d(1,n2,1), d(1,n2,n3), &
          d(n1,1,1), d(n1,1,n3), d(n1,n2,1), d(n1,n2,n3), d(i,j,k) )

      end do

    end do

  end do

  call r8block_print ( n1, n2, n3, d, '  Interpolated data array' )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests BLEND_113.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 3
  integer ( kind = 4 ), parameter :: n2 = 5
  integer ( kind = 4 ), parameter :: n3 = 4

  real ( kind = 8 ) d(n1,n2,n3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t

  d(1:n1,1:n2,1:n3) = 0.0D+00

  do i = 1, n1
    do j = 1, n2
      do k = 1, n3
        if ( &
          ( ( i == 1 .or. i == n1 ) .and. ( j == 1 .or. j == n2 ) ) .or. &
          ( ( i == 1 .or. i == n1 ) .and. ( k == 1 .or. k == n3 ) ) .or. &
          ( ( j == 1 .or. j == n2 ) .and. ( k == 1 .or. k == n3 ) ) ) then
          d(i,j,k) = real ( i + j + k, kind = 8 )
        end if
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  BLEND_113 blends edge values into a table.'
  write ( *, '(a)' ) ' '

  call r8block_print ( n1, n2, n3, d, '  Initial data array' )

  do i = 1, n1

    r = real (  i - 1, kind = 8 ) &
      / real ( n1 - 1, kind = 8 )

    do j = 1, n2

      s = real (  j - 1, kind = 8 ) &
        / real ( n2 - 1, kind = 8 )

      do k = 1, n3

        t = real (  k - 1, kind = 8 ) &
          / real ( n3 - 1, kind = 8 )

        if ( &
          ( ( i == 1 .or. i == n1 ) .and. ( j == 1 .or. j == n2 ) ) .or. &
          ( ( i == 1 .or. i == n1 ) .and. ( k == 1 .or. k == n3 ) ) .or. &
          ( ( j == 1 .or. j == n2 ) .and. ( k == 1 .or. k == n3 ) ) ) then
          cycle
        end if

        call blend_113 ( r, s, t, &
          d(1,1,1), d(1,1,n3), d(1,n2,1), d(1,n2,n3), &
          d(n1,1,1), d(n1,1,n3), d(n1,n2,1), d(n1,n2,n3), &
          d(i,1,1), d(i,1,n3), d(i,n2,1), d(i,n2,n3), &
          d(1,j,1), d(1,j,n3), d(n1,j,1), d(n1,j,n3), &
          d(1,1,k), d(1,n2,k), d(n1,1,k), d(n1,n2,k), &
          d(i,j,k) )

      end do

    end do

  end do

  call r8block_print ( n1, n2, n3, d, '  Interpolated data array' )

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests BLEND_123.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n1 = 3
  integer ( kind = 4 ), parameter :: n2 = 5
  integer ( kind = 4 ), parameter :: n3 = 4

  real ( kind = 8 ) d(n1,n2,n3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t

  d(1:n1,1:n2,1:n3) = 0.0D+00

  do i = 1, n1
    do j = 1, n2
      do k = 1, n3
        if ( &
          ( i == 1 .or. i == n1 ) .or. &
          ( j == 1 .or. j == n2 ) .or. &
          ( k == 1 .or. k == n3  ) ) then
          d(i,j,k) = real ( i + j + k, kind = 8 )
        end if
      end do
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  BLEND_123 blends face values into a table.'
  write ( *, '(a)' ) ' '

  call r8block_print ( n1, n2, n3, d, '  Initial data array' )

  do i = 1, n1

    r = real (  i - 1, kind = 8 ) &
      / real ( n1 - 1, kind = 8 )

    do j = 1, n2

      s = real (  j - 1, kind = 8 ) &
        / real ( n2 - 1, kind = 8 )

      do k = 1, n3

        t = real (  k - 1, kind = 8 ) &
          / real ( n3 - 1, kind = 8 )

        if ( &
          ( i == 1 .or. i == n1 ) .or. &
          ( j == 1 .or. j == n2 ) .or. &
          ( k == 1 .or. k == n3  ) ) then
          cycle
        end if

        call blend_123 ( r, s, t, &
          d(1,1,1), d(1,1,n3), d(1,n2,1), d(1,n2,n3), &
          d(n1,1,1), d(n1,1,n3), d(n1,n2,1), d(n1,n2,n3), &
          d(i,1,1), d(i,1,n3), d(i,n2,1), d(i,n2,n3), &
          d(1,j,1), d(1,j,n3), d(n1,j,1), d(n1,j,n3), &
          d(1,1,k), d(1,n2,k), d(n1,1,k), d(n1,n2,k), &
          d(1,j,k), d(n1,j,k), d(i,1,k), d(i,n2,k), d(i,j,1), d(i,j,n3), &
          d(i,j,k) )

      end do

    end do

  end do

  call r8block_print ( n1, n2, n3, d, '  Interpolated data array' )

  return
end
