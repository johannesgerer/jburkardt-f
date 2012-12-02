program main

!*****************************************************************************80
!
!! MAIN is the main program for IHS_PRB.
!
!  Discussion:
!
!    IHS_PRB tests the improved hypercube sampling algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'IHS_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the IHS library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'IHS_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests the improved distributed hypercube sampling algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 8 ) average
  real    ( kind = 8 ) covc
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), parameter :: duplication = 5
  integer ( kind = 4 ) j
  real    ( kind = 8 ) opt
  integer ( kind = 4 ), parameter :: point_num = 10
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) std
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  IHS implements the IHS Algorithm'
  write ( *, '(a)' ) '  (Improved Distributed Hypercube Sampling)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate the code for a fixed number of points'
  write ( *, '(a)' ) '  and an increasing dimension.'

  do dim_num = 1, 4

    allocate ( x(1:dim_num,1:point_num) )

    seed = 17
    opt = real ( point_num, kind = 8 ) / &
      ( real ( point_num, kind = 8 ) )**( 1.0D+00 / real ( dim_num, kind = 8 ) )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)'   ) '  Random number seed =       ', seed
    write ( *, '(a,i8 )'   ) '  Spatial dimension =        ', dim_num
    write ( *, '(a,i8 )'   ) '  Number of points =         ', point_num
    write ( *, '(a,i8 )'   ) '  Duplication factor =       ', duplication
    write ( *, '(a,g14.6)' ) '  Desired minimum distance = ', opt
!
!  Get the points.
!
    call ihs ( dim_num, point_num, duplication, seed, x )
!
!  Compute the covariance.
!
    call covariance ( dim_num, point_num, x, average, std, covc )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Average minimum distance ', average
    write ( *, '(a,g14.6)' ) '  Standard deviation:      ', std
    write ( *, '(a,g14.6)' ) '  Covariance:              ', covc

    write ( *, '(a)' ) ' '

    do j = 1, point_num
      write ( *, '(i4,4x,10i4)' ) j, x(1:dim_num,j)
    end do

    deallocate ( x )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests the improved distributed hypercube sampling algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: point_num = 10

  real    ( kind = 8 ) average
  real    ( kind = 8 ) covc
  integer ( kind = 4 ) duplication
  integer ( kind = 4 ) j
  real    ( kind = 8 ) opt
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) std
  integer ( kind = 4 ) x(dim_num,point_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  IHS implements the IHS Algorithm'
  write ( *, '(a)' ) '  (Improved Distributed Hypercube Sampling)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate the code for a fixed number of points'
  write ( *, '(a)' ) '  and dimension, but vary the duplication value.'

  opt = real ( point_num, kind = 8 ) / &
    ( real ( point_num, kind = 8 ) )**( 1.0D+00 / real ( dim_num, kind = 8 ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8 )'   ) '  Spatial dimension =  ', dim_num
  write ( *, '(a,i8 )'   ) '  Number of points =   ', point_num
  write ( *, '(a,g14.6)' ) '  Desired minimum distance = ', opt

  do duplication = 1, 5

    seed = 17

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)'   ) '  Random number seed = ', seed
    write ( *, '(a,i8 )'   ) '  Duplication factor = ', duplication
!
!  Get the points.
!
    call ihs ( dim_num, point_num, duplication, seed, x )
!
!  Compute the covariance.
!
    call covariance ( dim_num, point_num, x, average, std, covc )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Average minimum distance ', average
    write ( *, '(a,g14.6)' ) '  Standard deviation:      ', std
    write ( *, '(a,g14.6)' ) '  Covariance:              ', covc

    write ( *, '(a)' ) ' '

    do j = 1, point_num
      write ( *, '(i4,4x,10i4)' ) j, x(1:dim_num,j)
    end do

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests the improved distributed hypercube sampling algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 8 ) average
  real    ( kind = 8 ) covc
  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: duplication = 5
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) opt
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) std
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  IHS implements the IHS Algorithm'
  write ( *, '(a)' ) '  (Improved Distributed Hypercube Sampling)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate the code for a fixed dimension'
  write ( *, '(a)' ) '  and duplication value, and increasing number of points.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8 )'   ) '  Spatial dimension =  ', dim_num
  write ( *, '(a,i8 )'   ) '  Duplication factor = ', duplication

  do i = 1, 5

    point_num = 10 * 2**(i-1)

    allocate ( x(1:dim_num,1:point_num) )

    opt = real ( point_num, kind = 8 ) / &
      ( real ( point_num, kind = 8 ) )**( 1.0D+00 / real ( dim_num, kind = 8 ) )

    seed = 17

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)'   ) '  Random number seed = ', seed
    write ( *, '(a,i8 )'   ) '  Number of points =   ', point_num
    write ( *, '(a,g14.6)' ) '  Desired minimum distance = ', opt
!
!  Get the points.
!
    call ihs ( dim_num, point_num, duplication, seed, x )
!
!  Compute the covariance.
!
    call covariance ( dim_num, point_num, x, average, std, covc )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Average minimum distance ', average
    write ( *, '(a,g14.6)' ) '  Standard deviation:      ', std
    write ( *, '(a,g14.6)' ) '  Covariance:              ', covc

    write ( *, '(a)' ) ' '

    do j = 1, point_num
      if ( j <= 10 .or. point_num-10 <= j ) then
        write ( *, '(i4,4x,10i4)' ) j, x(1:dim_num,j)
      else if ( j == 11 ) then
        write ( *, '(a)' ) '....    ........'
      end if

    end do

    deallocate ( x )

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests the improved distributed hypercube sampling algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: point_num = 10

  real    ( kind = 8 ) average
  real    ( kind = 8 ) covc
  integer ( kind = 4 ), parameter :: duplication = 5
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) opt
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) std
  integer ( kind = 4 ) x(dim_num,point_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  IHS implements the IHS Algorithm'
  write ( *, '(a)' ) '  (Improved Distributed Hypercube Sampling)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate the code for a fixed number of points,'
  write ( *, '(a)' ) '  dimension, and duplication factor, but with a'
  write ( *, '(a)' ) '  varying random number seed.'

  opt = real ( point_num, kind = 8 ) / &
    ( real ( point_num, kind = 8 ) )**( 1.0D+00 / real ( dim_num, kind = 8 ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8 )'   ) '  Spatial dimension =        ', dim_num
  write ( *, '(a,i8 )'   ) '  Number of points =         ', point_num
  write ( *, '(a,i8 )'   ) '  Duplication factor =       ', duplication
  write ( *, '(a,g14.6)' ) '  Desired minimum distance = ', opt

  seed = 17

  do i = 1, 4

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)'   ) '  Random number seed =       ', seed
!
!  Get the points.
!
    call ihs ( dim_num, point_num, duplication, seed, x )
!
!  Compute the covariance.
!
    call covariance ( dim_num, point_num, x, average, std, covc )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Average minimum distance ', average
    write ( *, '(a,g14.6)' ) '  Standard deviation:      ', std
    write ( *, '(a,g14.6)' ) '  Covariance:              ', covc

    write ( *, '(a)' ) ' '

    do j = 1, point_num
      write ( *, '(i4,4x,10i4)' ) j, x(1:dim_num,j)
    end do

  end do

  return
end
