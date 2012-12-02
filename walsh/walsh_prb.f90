program main

!*****************************************************************************80
!
!! MAIN is the main program for WALSH_PRB.
!
!  Discussion:
!
!    WALSH_PRB tests the WALSH library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WALSH_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the WALSH library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WALSH_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests FWT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 16

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) work(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  FWT computes a fast Walsh transform.'

  do j = 1, 2

    if ( j == 1 ) then
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, w )
    else
      do i = 1, n
        w(i) = real ( i, kind = 8 )
      end do
    end if

    x(1:n) = w(1:n)
    call fwt ( n, w, work )
    y(1:n) = w(1:n) / real ( n, kind = 8 )
    call fwt ( n, w, work )
    z(1:n) = w(1:n) / real ( n, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     I        X(I)    Y=FWT(X)/N   Z=FWT(Y)/N'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) i, x(i), y(i), z(i)
    end do

  end do
    
  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests WALSH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 16

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) work(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  WALSH computes a fast Walsh transform.'

  do j = 1, 2

    if ( j == 1 ) then
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, w )
    else
      do i = 1, n
        w(i) = real ( i, kind = 8 )
      end do
    end if

    x(1:n) = w(1:n)
    call walsh ( n, w, work )
    y(1:n) = w(1:n) / real ( n, kind = 8 )
    call walsh ( n, w, work )
    z(1:n) = w(1:n) / real ( n, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     I        X(I)    Y=FWT(X)/N   Z=FWT(Y)/N'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) i, x(i), y(i), z(i)
    end do

  end do
    
  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests HAAR, HAARIN and HNORM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 16

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) work(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  HAAR computes a Haar transform.'
  write ( *, '(a)' ) '  HNORM normalizes the transformed data.'
  write ( *, '(a)' ) '  HAARIN computes an inverse Haar transform.'

  do j = 1, 2

    if ( j == 1 ) then
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, w )
    else
      do i = 1, n
        w(i) = real ( i, kind = 8 )
      end do
    end if

    x(1:n) = w(1:n)

    call haar ( n, w, work )

    y(1:n) = w(1:n)

    call hnorm ( n, w )

    z(1:n) = w(1:n)

    call haarin ( n, w, work )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     I        X(I)    Y=HAAR(X)  Z=HNORM(Y)  W=HAARIN(Z)'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
        i, x(i), y(i), z(i), w(i)
    end do

  end do
    
  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests FFWT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 16

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  FFWT computes a fast Walsh transform.'

  do j = 1, 2

    if ( j == 1 ) then
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, w )
    else
      do i = 1, n
        w(i) = real ( i, kind = 8 )
      end do
    end if

    x(1:n) = w(1:n)
    call ffwt ( n, w )
    y(1:n) = w(1:n) / real ( n, kind = 8 )
    call ffwt ( n, w )
    z(1:n) = w(1:n) / real ( n, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     I        X(I)   Y=FFWT(X)/N  Z=FFWT(Y)/N'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
        i, x(i), y(i), z(i)
    end do

  end do
    
  return
end
