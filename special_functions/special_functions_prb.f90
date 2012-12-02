program main

!*****************************************************************************80
!
!! MAIN tests the SPECIAL_FUNCTIONS library.
!
!  Modified:
!
!    03 July 2012
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPECIAL_FUNCTIONS_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SPECIAL_FUNCTIONS library.'

  call mbeta ( )
  call mcisia ( )
  call mcisib ( )
  call mhygfx ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPECIAL_FUNCTIONS_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine mbeta ( )

!*****************************************************************************80
!
!! MBETA tests BETA.
!
!  Modified:
!
!    12 March 2012
!
!  Example:
!
!                 p       q           B(p,q)
!               ---------------------------------
!                1.5     2.0     .2666666667D+00
!                2.5     2.0     .1142857143D+00
!                1.5     3.0     .1523809524D+00
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) bt
  real ( kind = 8 ) p
  real ( kind = 8 ), dimension ( test_num ) :: p_test = (/ &
    1.5D+00, 2.5D+00, 1.5D+00 /)
  real ( kind = 8 ) q
  real ( kind = 8 ), dimension ( test_num ) :: q_test = (/ &
    2.0D+00, 2.0D+00, 3.0D+00 /)
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MBETA:'
  write ( *, '(a)' ) '  Test BETA.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    p       q           B(p,q)'
  write ( *, '(a)' ) '  ---------------------------------'

  do test = 1, test_num

    p = p_test(test)
    q = q_test(test)

    call beta ( p, q, bt )
    write ( *, '(2x,f5.1,3x,f5.1,d20.10)' ) p, q, bt

  end do

  return
end
subroutine mcisia ( )

!*****************************************************************************80
!
!! MCISIA tests CISIA.
!
!  Modified:
!
!    03 July 2012
!
!  Example:
!
!      x        Ci(x)           Si(x)
!    ------------------------------------
!     0.0    - oo                 0
!     5.0    -.190030D+00      1.549931
!    10.0    -.454563D-01      1.658348
!    20.0     .444201D-01      1.548241
!    30.0    -.330326D-01      1.566757
!    40.0     .190201D-01      1.586985
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 6

  double precision ci
  double precision si
  integer test
  double precision x
  double precision, save, dimension ( test_num ) :: x_test = (/ &
    0.0D+00, 5.0D+00, 10.0D+00, 20.0D+00, 30.0D+00, 40.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MCISIA'
  write ( *, '(a)' ) '  Test CISIA, which computes the'
  write ( *, '(a)' ) '  cosine and sine integrals.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   x        ci(x)           si(x)'
  write ( *, '(a)' ) '------------------------------------'

  do test = 1, test_num

    x = x_test(test)

    call cisia ( x, ci, si )

    write ( *, '(1x,f5.1,g16.8,g16.8)' ) x, ci, si

  end do

  return
end
subroutine mcisib ( )

!*****************************************************************************80
!
!! MCISIB tests CISIB.
!
!  Modified:
!
!    01 July 2012
!
!  Example:
!
!      x        Ci(x)           Si(x)
!    ------------------------------------
!     0.0    - oo                 0
!     5.0    -.190030D+00      1.549931
!    10.0    -.454563D-01      1.658348
!    20.0     .444201D-01      1.548241
!    30.0    -.330326D-01      1.566757
!    40.0     .190201D-01      1.586985
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 6

  double precision ci
  double precision si
  integer test
  double precision x
  double precision, save, dimension ( test_num ) :: x_test = (/ &
    0.0D+00, 5.0D+00, 10.0D+00, 20.0D+00, 30.0D+00, 40.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MCISIB'
  write ( *, '(a)' ) '  Test CISIB, which computes the'
  write ( *, '(a)' ) '  cosine and sine integrals.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   x        ci(x)           si(x)'
  write ( *, '(a)' ) '------------------------------------'

  do test = 1, test_num

    x = x_test(test)

    call cisib ( x, ci, si )

    write ( *, '(1x,f5.1,g16.8,g16.8)' ) x, ci, si

  end do

  return
end
subroutine mhygfx ( )

!*****************************************************************************80
!
!! MHYGFX tests HYGFX.
!
!  Modified:
!
!    12 March 2012
!
!  Example:
!
!     A    B     C     X     F(A,B,C,X)
!
!   -2.5  3.3   6.7  0.25  0.72356129D+00
!   -0.5  3.3   6.7  0.25  0.93610145D+00    
!    0.5  3.3   6.7  0.25  0.10689695D+01    
!    2.5  3.3   6.7  0.25  0.14051563D+01
!
!   -2.5  3.3   6.7  0.55  0.46961432D+00
!   -0.5  3.3   6.7  0.55  0.85187390D+00
!    0.5  3.3   6.7  0.55  0.11795358D+01
!    2.5  3.3   6.7  0.55  0.23999063D+01
!
!   -2.5  3.3   6.7  0.85  0.29106096D+00
!   -0.5  3.3   6.7  0.85  0.75543187D+00
!    0.5  3.3   6.7  0.85  0.13510497D+00
!    2.5  3.3   6.7  0.85  0.57381566D+01
!
!    3.3  6.7  -5.5  0.25  0.15090670D+05
!    3.3  6.7  -0.5  0.25 -0.21631479D+04
!    3.3  6.7   0.5  0.25  0.26451677D+03
!    3.3  6.7   4.5  0.25  0.41946916D+01
!
!    3.3  6.7  -5.5  0.55  0.10170778D+11
!    3.3  6.7  -0.5  0.55 -0.30854772D+07
!    3.3  6.7   0.5  0.55  0.11967860D+06
!    3.3  6.7   4.5  0.55  0.58092729D+02
!
!    3.3  6.7  -5.5  0.85  0.58682088D+19
!    3.3  6.7  -0.5  0.85 -0.10217370D+13
!    3.3  6.7   0.5  0.85  0.92370648D+10
!    3.3  6.7   4.5  0.85  0.20396914D+05 
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a_test(4)
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) c_test(4)
  real ( kind = 8 ) hf
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) x
  real ( kind = 8 ) x_test(3)

  save a_test
  save c_test
  save x_test

  data a_test / -2.5D+00, -0.5D+00, 0.5D+00, 2.5D+00 /
  data c_test / -5.5D+00, -0.5D+00, 0.5D+00, 4.5D+00/
  data x_test /  0.25D+00, 0.55D+00, 0.85D+00 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MHYGFX:'
  write ( *, '(a)' ) '  HYGFX evaluates the hypergeometric function 2F1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     A              B            C            X             F(A,B,C,X)'

  do l = 1, 3
    x = x_test(l)
    c = 6.7D+00
    b = 3.3D+00
    write ( *, '(a)' ) ' '
    do i = 1, 4
      a = a_test(i)
      call hygfx ( a, b, c, x, hf )
      write ( *, '(4g14.6,g24.16)' ) a, b, c, x, hf
    end do
  end do

  do l = 1, 3
    x = x_test(l)
    write ( *, '(a)' ) ' '
    do k = 1, 4
      c = c_test(k)
      b = 6.7D+00
      a = 3.3D+00
      call hygfx ( a, b, c, x, hf )
      write ( *, '(4g14.6,g24.16)' ) a, b, c, x, hf
    end do
  end do

  return
end
