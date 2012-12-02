program main

!*****************************************************************************80
!
!! TOMS446_PRB tests algorithm 446.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS446_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS446 library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS446_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests CHEBY, which computes Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nf = 5
  integer ( kind = 4 ), parameter :: npl = 10

  external functn
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(npl,nf)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Test CHEBY, which computes the '
  write ( *, '(a)' ) '  Chebyshev series for several functions.'

  call cheby ( nf, npl, functn, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        Sin(x)      Cos(x)    Sin(2x)     Cos(2x)       X^5'
  write ( *, '(a)' ) ' '

  do i = 1, npl
    write ( *, '(5(2x,f10.4))' ) x(i,1:nf)
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests MULTPLY, which multiplies two Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nf = 5
  integer ( kind = 4 ), parameter :: npl = 10

  external functn
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(npl,nf)
  real ( kind = 8 ) x1(npl)
  real ( kind = 8 ) x2(npl)
  real ( kind = 8 ) x3(npl)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Test MLTPLY, which computes the '
  write ( *, '(a)' ) '  product of two Chebyshev series.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Multiply series for SIN(X) and COS(X)'
  write ( *, '(a)' ) '  and compare with series for 1/2*SIN(2X).'

  call cheby ( nf, npl, functn, x )

  x1(1:npl) = x(1:npl,1)
  x2(1:npl) = x(1:npl,2)
  x(1:npl,3) = 0.5D+00 * x(1:npl,3)

  call mltply ( x1, x2, npl, x3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        Sin(x)      Cos(x)   1/2*Sin(2x)     RESULT'
  write ( *, '(a)' ) ' '

  do i = 1, npl
    write ( *, '(5(2x,f10.4))' ) x(i,1:3), x3(i)
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests ECHEB, which evaluates a Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nf = 5
  integer ( kind = 4 ), parameter :: npl = 10

  external functn
  real ( kind = 8 ) fval
  real ( kind = 8 ) fxj(nf)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nx
  real ( kind = 8 ) x(npl,nf)
  real ( kind = 8 ) x2(npl)
  real ( kind = 8 ) xval

  nx = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Test ECHEB, which evaluates a '
  write ( *, '(a)' ) '  Chebyshev series.'

  call cheby ( nf, npl, functn, x )

  do j = 1, nf

    x2(1:npl) = x(1:npl,j)

    write ( *, '(a)' ) ' '
    if ( j == 1 ) then
      write ( *, '(a)' ) '  Sin(x)'
    else if ( j == 2 ) then
      write ( *, '(a)' ) '  Cos(x)'
    else if ( j == 3 ) then
      write ( *, '(a)' ) '  Sin(2x)'
    else if ( j == 4 ) then
      write ( *, '(a)' ) '  Cos(2x)'
    else if ( j == 5 ) then
      write ( *, '(a)' ) '  x^5'
    end if

    write ( *, '(a)' ) ' '

    do k = 1, nx

      xval = 2.0D+00 * real ( k - 1, kind = 8 ) / real ( nx - 1, kind = 8 ) &
        - 1.0D+00

      call functn ( xval, fxj )

      call echeb ( xval, x2, npl, fval )

      write ( *, '(5(2x,f10.4))' ) xval, fxj(j), fval

    end do

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests EDCHEB, which evaluates the derivative of a Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nf = 5
  integer ( kind = 4 ), parameter :: npl = 10

  external functn
  real ( kind = 8 ) fval
  real ( kind = 8 ) fxj(nf)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nx
  real ( kind = 8 ) x(npl,nf)
  real ( kind = 8 ) x2(npl)
  real ( kind = 8 ) xval

  nx = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Test EDCHEB, which evaluates the  '
  write ( *, '(a)' ) '  derivative of a Chebyshev series.'

  call cheby ( nf, npl, functn, x )

  do j = 1, nf

    x2(1:npl) = x(1:npl,j)

    write ( *, '(a)' ) ' '
    if ( j == 1 ) then
      write ( *, '(a)' ) '  d/dx Sin(x)'
    else if ( j == 2 ) then
      write ( *, '(a)' ) '  d/dx Cos(x)'
    else if ( j == 3 ) then
      write ( *, '(a)' ) '  d/dx Sin(2x)'
    else if ( j == 4 ) then
      write ( *, '(a)' ) '  d/dx Cos(2x)'
    else if ( j == 5 ) then
      write ( *, '(a)' ) '  d/dx x^5'
    end if

    write ( *, '(a)' ) ' '

    do k = 1, nx

      xval = 2.0D+00 * real ( k - 1, kind = 8 ) / real ( nx - 1, kind = 8 ) &
        - 1.0D+00

      call functn_d ( xval, fxj )

      call edcheb ( xval, x2, npl, fval )

      write ( *, '(5(2x,f10.4))' ) xval, fxj(j), fval

    end do

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests DFRNT, which computes the Chebyshev series of a derivative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nf = 5
  integer ( kind = 4 ), parameter :: npl = 10

  external functn
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(npl,nf)
  real ( kind = 8 ) x2(npl)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Test DFRNT, which computes the '
  write ( *, '(a)' ) '  Chebyshev series for the derivative'
  write ( *, '(a)' ) '  of several functions.'

  call cheby ( nf, npl, functn, x )

  do j = 1, nf
    x2(1:npl) = x(1:npl,j)
    call dfrnt ( x2, npl, x2 )
    x(1:npl,j) = x2(1:npl)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Chebyshev series for d/dx of:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        Sin(x)      Cos(x)    Sin(2x)     Cos(2x)       X^5'
  write ( *, '(a)' ) ' '

  do i = 1, npl
    write ( *, '(5(2x,f10.4))' ) x(i,1:nf)
  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests NTGRT, which computes the Chebyshev series of an indefinite integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nf = 5
  integer ( kind = 4 ), parameter :: npl = 10

  external functn
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(npl,nf)
  real ( kind = 8 ) x2(npl)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Test NTGRT, which computes the '
  write ( *, '(a)' ) '  Chebyshev series for the indefinite'
  write ( *, '(a)' ) '  integral of several functions.'

  call cheby ( nf, npl, functn, x )

  do j = 1, nf
    x2(1:npl) = x(1:npl,j)
    call ntgrt ( x2, npl, x2 )
    x(1:npl,j) = x2(1:npl)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Series for indefinite integral of:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        Sin(x)      Cos(x)    Sin(2x)     Cos(2x)       X^5'
  write ( *, '(a)' ) ' '

  do i = 1, npl
    write ( *, '(5(2x,f10.4))' ) x(i,1:nf)
  end do

  return
end
subroutine functn ( x, fxj )

!*****************************************************************************80
!
!! FUNCTN evaluates several functions at X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) FXJ(5), the derivative values.
!
  implicit none

  integer ( kind = 4 ), parameter :: nf = 5

  real ( kind = 8 ) fxj(nf)
  real ( kind = 8 ) x

  fxj(1) = sin ( x )
  fxj(2) = cos ( x )
  fxj(3) = sin ( 2.0D+00 * x )
  fxj(4) = cos ( 2.0D+00 * x )
  fxj(5) = x**5

  return
end
subroutine functn_d ( x, fxj )

!*****************************************************************************80
!
!! FUNCTN_D evaluates the derivatives of several functions at X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) FXJ(5), the derivative values.
!
  implicit none

  integer ( kind = 4 ), parameter :: nf = 5

  real ( kind = 8 ) fxj(nf)
  real ( kind = 8 ) x

  fxj(1) =  cos ( x )
  fxj(2) = -sin ( x )
  fxj(3) =  2.0D+00 * cos ( 2.0D+00 * x )
  fxj(4) = -2.0D+00 * sin ( 2.0D+00 * x )
  fxj(5) =  5.0D+00 * x**4

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
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
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8  ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
