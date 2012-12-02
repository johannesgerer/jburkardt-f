program main

!*****************************************************************************80
!
!! MAIN is the main program for LAMBERT_PRB.
!
!  Discussion:
!
!    LAMBERT_PRB tests the LAMBERT sequence routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 December 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAMBERT_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LAMBERT library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAMBERT_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests LAMBERT1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 December 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 50

  real ( kind = 8 ) eta(1,n)
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  LAMBERT1 computes entries of the 1D Lambert sequence.'

  call lambert1 ( n, eta )

  do j = 1, n
    write ( *, '(i4,2f10.4)' ) j, eta(1,j)
  end do

  call lambert_write ( 1, n, eta, 'lambert_01_00050.txt' )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests LAMBERT2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 December 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 50

  real ( kind = 8 ) eta(2,n)
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  LAMBERT2 computes entries of the 2D Lambert sequence.'

  call lambert2 ( n, eta )

  do j = 1, n
    write ( *, '(i4,2f10.4)' ) j, eta(1:2,j)
  end do

  call lambert_write ( 2, n, eta, 'lambert_02_00050.txt' )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests LAMBERT3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 December 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 50

  real ( kind = 8 ) eta(3,n)
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  LAMBERT3 computes entries of the 3D Lambert sequence.'

  call lambert3 ( n, eta )

  do j = 1, n
    write ( *, '(i4,3f10.4)' ) j, eta(1:3,j)
  end do

  call lambert_write ( 3, n, eta, 'lambert_03_00050.txt' )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests LAMBERT4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 December 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 50

  real ( kind = 8 ) eta(4,n)
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  LAMBERT4 computes entries of the 4D Lambert sequence.'

  call lambert4 ( n, eta )

  do j = 1, n
    write ( *, '(i4,4f10.4)' ) j, eta(1:4,j)
  end do

  call lambert_write ( 3, n, eta, 'lambert_04_00050.txt' )

  return
end

