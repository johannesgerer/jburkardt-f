program main

!*****************************************************************************80
!
!! MAIN tests the PINK_NOISE routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PINK_NOISE_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PINK_NOISE library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) 'PINK_NOISE_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests WRAP2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  integer ( kind = 4 ) q_in

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  WRAP2 performs a circular wrap.'
  write ( *, '(a)' ) '  Q is expected to range between 0 and M.'
  write ( *, '(a)' ) '  WRAP2 takes an input value of Q, and either'
  write ( *, '(a)' ) '  increments it by M+1 until in the range, or'
  write ( *, '(a)' ) '  decrements it by M+1 until in the range,'
  write ( *, '(a)' ) '  and returns the result as the function value.'

  do m = 2, 4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '   M  Qin  Qout'
    write ( *, '(a)' ) ' '
    do i = -5, 3 * m - 1
      q = i
      q_in = q
      call wrap2 ( m, q )
      write ( *, '(2x,i2,2x,i2,2x,i2)' ) m, q_in, q
    end do
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests CDELAY2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) q
  integer ( kind = 4 ) q_in

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  CDELAY2 is a circular buffer implementation'
  write ( *, '(a)' ) '  of an M-fold delay.  Q is a counter'
  write ( *, '(a)' ) '  which is decremented by CDELAY2, but reset to M'
  write ( *, '(a)' ) '  after it reaches 0.'

  do m = 2, 4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '   I   M  Qin  Qout'
    write ( *, '(a)' ) ' '
    q = m
    do i = 1, 3 * ( m + 1 )
      q_in = q
      call cdelay2 ( m, q )
      write ( *, '(2x,i2,2x,i2,2x,i2,2x,i2)' ) i, m, q_in, q
    end do
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests RANH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) q
  real    ( kind = 8 ) ranh
  real    ( kind = 8 ) u
  real    ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  RANH is a random hold function.'
  write ( *, '(a)' ) '  Given a value U and a delay D, it returns the value'
  write ( *, '(a)' ) '  U for D calls, then resets U.'

  do d = 5, 1, -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '   I   D   Q      U           Y'
    write ( *, '(a)' ) ' '
    u = 0.5D+00
    q = 3
    do i = 1, 20
      y = ranh ( d, u, q )
      write ( *, '(2x,i2,2x,i2,2x,i2,2x,f10.6,2x,f10.6)' ) i, d, q, u, y
    end do
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests RAN1F.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ), allocatable :: q(:)
  real    ( kind = 8 ) ran1f
  integer ( kind = 4 ) rep
  real    ( kind = 8 ), allocatable :: u(:)
  real    ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  RAN1F generates random values with an approximate'
  write ( *, '(a)' ) '  1/F distribution.'

  b = 1

  do while ( b < 32 )

    allocate ( u(1:b) )
    allocate ( q(1:b) )

    do rep = 1, 4

      call random_number ( harvest = u(1:b) )
      q(1:b) = 0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   B   I      Y'
      write ( *, '(a)' ) ' '

      do i = 1, 20
        y = ran1f ( b, u, q )
        write ( *, '(2x,i2,2x,i2,2x,f10.6)' ) b, i, y 
      end do

    end do

    deallocate ( q )
    deallocate ( u )

    b = b * 2

  end do

  return
end
