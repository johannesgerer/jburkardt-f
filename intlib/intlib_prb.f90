program main

!*****************************************************************************80
!
!! MAIN is the main program for INTLIB_PRB.
!
!  Discussion:
!
!    INTLIB_PRB runs the INTLIB tests.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INTLIB_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the INTLIB library.'

  call test1d ( )

  call test27 ( )
  call test28 ( )
  call test29 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INTLIB_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test1d ( )

!*****************************************************************************80
!
!! TEST1D tests all the 1D integration codes.
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

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: f1d1
  real ( kind = 8 ), external :: fbd1
  real ( kind = 8 ), external :: fed1
  real ( kind = 8 ), external :: fqd1
  real ( kind = 8 ), external :: fxd1
  real ( kind = 8 ), external :: fx2d1
  real ( kind = 8 ), external :: fx3d1

  a = 0.0D+00
  b = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1D'
  write ( *, '(a)' ) '  Test 1D quadrature codes'
  write ( *, '(a)' ) '  for integral of F(X) on [A,B].'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g24.16)' ) '  A = ', a
  write ( *, '(a,g24.16)' ) '  B = ', b

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  F(X)=1'
  write ( *, '(a)' ) ' '

  call tst1d ( f1d1, a, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  F(X)=X'
  write ( *, '(a)' ) ' '

  call tst1d ( fxd1, a, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  F(X)=X*2'
  write ( *, '(a)' ) ' '

  call tst1d ( fx2d1, a, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  F(X)=X*3'
  write ( *, '(a)' ) ' '

  call tst1d ( fx3d1, a, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  F(X)=EXP(X)'
  write ( *, '(a)' ) ' '

  call tst1d ( fed1, a, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  F(X)=SQRT(X)'
  write ( *, '(a)' ) ' '

  call tst1d ( fqd1, a, b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  F(X)=1/(1+X*X)'
  write ( *, '(a)') ' '

  call tst1d ( fbd1, a, b )

  return
end
subroutine tst1d ( func, a, b )

!*****************************************************************************80
!
!! TST1D calls all the 1D integrators.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: func

  call test01 ( func, a, b )
  call test02 ( func, a, b )
  call test03 ( func, a, b )
  call test04 ( func, a, b )
  call test05 ( func, a, b )
  call test06 ( func, a, b )
  call test07 ( func, a, b )
  call test08 ( func, a, b )
  call test09 ( func, a, b )
  call test10 ( func, a, b )
  call test11 ( func, a, b )
  call test12 ( func, a, b )
  call test13 ( func, a, b )
  call test14 ( func, a, b )
  call test15 ( func, a, b )
  call test16 ( func, a, b )
  call test17 ( func, a, b )
  call test88 ( func, a, b )

  return
end
subroutine test01 ( func, a, b )

!*****************************************************************************80
!
!! TEST01 tests PLINT.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ntab = 11

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) ftab(ntab)
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) result
  real ( kind = 8 ) xtab(ntab)

  call r8vec_even ( ntab, a, b, xtab )

  do i = 1, ntab
    ftab(i) = func ( xtab(i) )
  end do

  call plint ( ntab, xtab, ftab, a, b, result )

  write ( *, '(a,g24.16)' ) '  PLINT ', result

  return
end
subroutine test02 ( func, a, b )

!*****************************************************************************80
!
!! TEST02 tests AVINT.
!
!  Modified:
!
!    10 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ntab = 11

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) ftab(ntab)
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) result
  real ( kind = 8 ) xtab(ntab)

  call r8vec_even ( ntab, a, b, xtab )

  do i = 1, ntab
    ftab(i) = func ( xtab(i) )
  end do

  call avint ( ntab, xtab, ftab, a, b, result )

  write ( *, '(a,g24.16)' ) '  AVINT ', result

  return
end
subroutine test03 ( func, a, b )

!*****************************************************************************80
!
!! TEST03 tests CUBINT.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ntab = 11

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) error
  real ( kind = 8 ) ftab(ntab)
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) ib
  real ( kind = 8 ) result
  real ( kind = 8 ) xtab(ntab)

  ia = 1
  ib = ntab

  call r8vec_even ( ntab, a, b, xtab )

  do i = 1, ntab
    ftab(i) = func ( xtab(i) )
  end do

  call cubint ( ntab, xtab, ftab, ia, ib, result, error )

  write ( *, '(a,g24.16)' ) '  CUBINT', result

  return
end
subroutine test04 ( func, a, b )

!*****************************************************************************80
!
!! TEST04 tests WEDINT.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2
  integer ( kind = 4 ), parameter :: ntab = 6 * n + 1

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) ftab(ntab)
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  real ( kind = 8 ) result
  real ( kind = 8 ) xtab(ntab)

  h = ( b - a ) / real ( ntab - 1, kind = 8 )

  call r8vec_even ( ntab, a, b, xtab )

  do i = 1, ntab
    ftab(i) = func ( xtab(i) )
  end do

  call wedint ( ntab, h, ftab, result )

  write ( *, '(a,g24.16)' ) '  WEDINT', result

  return
end
subroutine test05 ( func, a, b )

!*****************************************************************************80
!
!! TEST05 tests CSPINT.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ntab = 13

  real ( kind = 8 ) a
  real ( kind = 8 ) aleft
  real ( kind = 8 ) b
  real ( kind = 8 ) brite
  real ( kind = 8 ) e(ntab)
  real ( kind = 8 ) ftab(ntab)
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) result
  real ( kind = 8 ) xtab(ntab)
  real ( kind = 8 ) work(ntab)
  real ( kind = 8 ) y(3,ntab)
!
!  Note that, for accuracy, it is useful to have two data points
!  outside the interval (A,B).
!
  aleft = a - ( b - a ) / real ( ntab - 3, kind = 8 )
  brite = b + ( b - a ) / real ( ntab - 3, kind = 8 )

  call r8vec_even ( ntab, aleft, brite, xtab )

  do i = 1, ntab
    ftab(i) = func ( xtab(i) )
  end do

  call cspint ( ntab, xtab, ftab, a, b, y, e, work, result )

  write ( *, '(a,g24.16)' ) '  CSPINT', result

  return
end
subroutine test06 ( func, a, b )

!*****************************************************************************80
!
!! TEST06 tests GAUS8.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) err
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) ier
  real ( kind = 8 ) result

  err = 0.001D+00

  call gaus8 ( func, a, b, err, result, ier )

  write ( *, '(a,g24.16)' ) '  GAUS8 ', result

  return
end
subroutine test07 ( func, a, b )

!*****************************************************************************80
!
!! TEST07 tests QNC79.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) err
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) k
  real ( kind = 8 ) result

  err = 0.001D+00

  call qnc79 ( func, a, b, err, result, ier, k )

  write ( *, '(a,g24.16)' ) '  QNC79 ', result

  return
end
subroutine test08 ( func, a, b )

!*****************************************************************************80
!
!! TEST08 tests QUAD.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nleast = 3
  integer ( kind = 4 ), parameter :: nmost = nleast + 2

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) relerr
  real ( kind = 8 ) result
  real ( kind = 8 ) work(nmost+1)

  abserr = 0.001D+00
  relerr = 0.001D+00

  call quad ( func, a, b, abserr, relerr, nleast, nmost, work, result )

  write ( *, '(a,g24.16)' ) '  QUAD  ', result

  return
end
subroutine test09 ( func, a, b )

!*****************************************************************************80
!
!! TEST09 tests RMINSP.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) epsin
  real ( kind = 8 ) epsout
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) iop
  real ( kind = 8 ) result

  epsin = 0.001D+00
  iop = 1

  call rminsp ( func, a, b, epsin, epsout, iop, result )

  write ( *, '(a,g24.16)' ) '  RMINSP', result

  return
end
subroutine test10 ( func, a, b )

!*****************************************************************************80
!
!! TEST10 tests RMINSP with cosine transformation.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) epsin
  real ( kind = 8 ) epsout
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) iop
  real ( kind = 8 ) result

  epsin = 0.001D+00
  iop = 2

  call rminsp ( func, a, b, epsin, epsout, iop, result )

  write ( *, '(a,g24.16)' ) '  RMINSP', result

  return
end
subroutine test11 ( func, a, b )

!*****************************************************************************80
!
!! TEST11 tests IRATEX.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) epsin
  real ( kind = 8 ) epsout
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) ind
  real ( kind = 8 ) result

  epsin = 0.001D+00

  call iratex ( func, a, b, epsin, epsout, result, ind )

  write ( *, '(a,g24.16)' ) '  IRATEX', result

  return
end
subroutine test12 ( func, a, b )

!*****************************************************************************80
!
!! TEST12 tests CADRE.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) error
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) ind
  real ( kind = 8 ) relerr
  real ( kind = 8 ) result

  abserr = 0.001D+00
  relerr = 0.001D+00

  call cadre ( func, a, b, abserr, relerr, error, result, ind )

  write ( *, '(a,g24.16)' ) '  CADRE ', result

  return
end
subroutine test13 ( func, a, b )

!*****************************************************************************80
!
!! TEST13 tests CHINSP.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) epsin
  real ( kind = 8 ) epsout
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) result

  epsin = 0.001D+00

  call chinsp ( func, a, b, epsin, epsout, result )

  write ( *, '(a,g24.16)' ) '  CHINSP', result

  return
end
subroutine test14 ( func, a, b )

!*****************************************************************************80
!
!! TEST14 tests SIMP.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) eps
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) result

  eps = 0.001D+00

  call simp ( func, a, b, eps, result )

  write ( *, '(a,g24.16)' ) '  SIMP  ', result

  return
end
subroutine test15 ( func, a, b )

!*****************************************************************************80
!
!! TEST15 tests HIORDQ.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ntab = 11
  integer ( kind = 4 ), parameter :: nwork = 2 * ( ntab - 1 )

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  real ( kind = 8 ) result
  real ( kind = 8 ) work(nwork)
  real ( kind = 8 ) x(ntab)
  real ( kind = 8 ) y(ntab)

  call r8vec_even ( ntab, a, b, x )

  do i = 1, ntab
    y(i) = func ( x(i) )
  end do

  h = ( b - a ) / real ( ntab - 1, kind = 8 )

  call hiordq ( ntab, h, y, work, result )

  write ( *, '(a,g24.16)' ) '  HIORDQ', result

  return
end
subroutine test16 ( func, a, b )

!*****************************************************************************80
!
!! TEST16 tests SIMPSN.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ntab = 11

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  real ( kind = 8 ) result
  real ( kind = 8 ) x(ntab)
  real ( kind = 8 ) y(ntab)

  call r8vec_even ( ntab, a, b, x )

  do i = 1, ntab
    y(i) = func ( x(i) )
  end do

  h = ( b - a ) / real ( ntab - 1, kind = 8 )

  call simpsn ( ntab, h, y, result )

  write ( *, '(a,g24.16)' ) '  SIMPSN', result

  return
end
subroutine test17 ( func, a, b )

!*****************************************************************************80
!
!! TEST17 tests SIMPNE.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ntab = 11

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) result
  real ( kind = 8 ) x(ntab)
  real ( kind = 8 ) y(ntab)

  call r8vec_even ( ntab, a, b, x )

  do i = 1, ntab
    y(i) = func ( x(i) )
  end do

  call simpne ( ntab, x, y, result )

  write ( *, '(a,g24.16)' ) '  SIMPNE', result

  return
end
subroutine test88 ( func, a, b )

!*****************************************************************************80
!
!! TEST88 tests MONTE_CARLO.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) n
  real ( kind = 8 ) result

  n = 1000

  call monte_carlo ( func, a, b, n, result )

  write ( *, '(a,g24.16)' ) '  MONTE ', result

  return
end
subroutine test27 ( )

!*****************************************************************************80
!
!! TEST27 tests FILON_COS.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable, dimension ( : ) :: ftab
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ntab
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) result
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable, dimension ( : ) :: xtab

  a = 0.0D+00
  b = 2.0D+00 * pi

  ntab = 11
  allocate ( ftab(1:ntab) )
  allocate ( xtab(1:ntab) )

  call r8vec_even ( ntab, a, b, xtab )

  h = ( b - a ) / real ( ntab - 1, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)') '  FILON_COS estimates the integral of.'
  write ( *, '(a)' ) '  F(X) * COS ( T * X )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrate F(X)*COS(T*X):'
  write ( *, '(a)' ) '  with F(X)=1, X, X**2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g24.16)' ) '  A = ', a
  write ( *, '(a,g24.16)' ) '  B = ', b
  write ( *, '(a,i6)' ) '  NTAB = ', ntab
  write ( *, '(a,g24.16)' ) '  H = ', h
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '       T                      Approximate             Exact'
  write ( *, '(a)' ) ' '

  do k = 1, 3

    if ( k == 1 ) then
      t = 1.0D+00
    else if ( k == 2 ) then
      t = 2.0D+00
    else if ( k == 3 ) then
      t = 10.0D+00
    end if

    do i = 1, 3

      do j = 1, ntab

        if ( i == 1 )then
          ftab(j) = 1.0D+00
        else
          ftab(j) = xtab(j)**( i - 1 )
        end if

      end do

      call filon_cos ( ntab, ftab, a, b, t, result )

      if ( i == 1 ) then
        exact = ( sin ( t * b ) - sin ( t * a ) ) / t
      else if ( i == 2 ) then
        exact = ( ( cos ( t * b ) + t * b * sin ( t * b ) ) &
                - ( cos ( t * a ) + t * a * sin ( t * a ) ) ) / t**2
      else if ( i == 3 ) then
        exact = ( ( 2.0D+00 * t * b * cos ( t * b ) &
              + ( t * t * b**2 - 2.0D+00 ) * sin ( t * b ) ) &
                - ( 2.0D+00 * t * a * cos ( t * a ) &
              + ( t * t * a**2 - 2.0D+00 ) * sin ( t * a ) ) ) / t**3
      end if

      write ( *, '(2x,3g24.16)' ) t, result, exact

    end do

    write ( *, '(a)' ) ' '

  end do

  deallocate ( ftab )
  deallocate ( xtab )
!
!  Example suggested by James Roedder.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrate log(1+X)*cos(T*X):'
  write ( *, '(a)' ) '  T = 10, and NTAB increases'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '    NTAB    H       Approximate             Exact                   Error'
  write ( *, '(a)' ) ' '

  do i = 1, 6

    ntab = 2**i * 10 + 1
    allocate ( ftab(1:ntab) )
    allocate ( xtab(1:ntab) )

    call r8vec_even ( ntab, a, b, xtab )

    h = ( b - a ) / real ( ntab - 1, kind = 8 )

    do j = 1, ntab
      ftab(j) = log ( 1.0D+00 + xtab(j) )
    end do

    t = 10.0D+00

    call filon_cos ( ntab, ftab, a, b, t, result )

    exact = -0.008446594405D+00
    error = result - exact

    write ( *, '(2x,i6,f8.4,g24.16,g24.16,g16.8)' ) &
      ntab, h, result, exact, error

    deallocate ( ftab )
    deallocate ( xtab )

  end do

  return
end
subroutine test28 ( )

!*****************************************************************************80
!
!! TEST28 tests FILON_SIN.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable, dimension ( : ) :: ftab
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ntab
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) result
  real ( kind = 8 ) t
  real ( kind = 8 ), allocatable, dimension ( : ) :: xtab

  a = 0.0D+00
  b = 2.0D+00 * pi

  ntab = 11
  allocate ( ftab(1:ntab) )
  allocate ( xtab(1:ntab) )

  call r8vec_even ( ntab, a, b, xtab )

  h = ( b - a ) / real ( ntab - 1, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST28'
  write ( *, '(a)' ) '  FILON_SIN estimates the integral of.'
  write ( *, '(a)' ) '  F(X) * SIN ( T * X )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g24.16)' ) '  A = ', a
  write ( *, '(a,g24.16)' ) '  B = ', b
  write ( *, '(a,i6)' ) '  NTAB = ', ntab
  write ( *, '(a,g24.16)' ) '  H = ', h
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrate F(X)*SIN(T*X)'
  write ( *, '(a)' ) '  with F(X)=1, X, X**2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '       T                      Approximate             Exact'
  write ( *, '(a)' ) ' '

  do k = 1, 3

    if ( k == 1 ) then
      t = 1.0D+00
    else if ( k == 2 ) then
      t = 2.0D+00
    else if ( k == 3 ) then
      t = 10.0D+00
    end if

    do i = 1, 3

      do j = 1, ntab

        if ( i == 1 )then
          ftab(j) = 1.0D+00
        else
          ftab(j) = xtab(j)**( i - 1 )
        end if

      end do

      call filon_sin ( ntab, ftab, a, b, t, result )

      if ( i == 1 ) then
        exact = ( - cos ( t * b ) + cos ( t * a ) ) / t
      else if ( i == 2 ) then
        exact = ( ( sin ( t * b ) - t * b * cos ( t * b ) ) &
                - ( sin ( t * a ) - t * a * cos ( t * a ) ) ) / t**2
      else if ( i == 3 ) then
        exact = ( ( 2.0D+00 * t * b * sin ( t * b ) &
                + ( 2.0D+00 - t**2 * b**2 ) * cos ( t * b ) ) &
                - ( 2.0D+00 * t * a * sin ( t * a ) &
                + ( 2.0D+00 - t**2 * a**2 ) * cos ( t * a ) ) ) / t**3
      end if

      write ( *, '(2x,3g24.16)' ) t, result, exact

    end do

    write ( *, '(a)' ) ' '

  end do

  deallocate ( ftab )
  deallocate ( xtab )
!
!  Example suggested by James Roedder.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrate log(1+X)*sin(T*X):'
  write ( *, '(a)' ) '  T = 10, and NTAB increases'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '    NTAB    H       Approximate             Exact                   Error'
  write ( *, '(a)' ) ' '

  do i = 1, 6

    ntab = 2**i * 10 + 1
    allocate ( ftab(1:ntab) )
    allocate ( xtab(1:ntab) )

    call r8vec_even ( ntab, a, b, xtab )

    h = ( b - a ) / real ( ntab - 1, kind = 8 )

    do j = 1, ntab
      ftab(j) = log ( 1.0D+00 + xtab(j) )
    end do

    t = 10.0D+00

    call filon_sin ( ntab, ftab, a, b, t, result )

    exact = -0.19762680771872D+00
    error = result - exact

    write ( *, '(2x,i6,f8.4,g24.16,g24.16,g16.8)' ) &
      ntab, h, result, exact, error

    deallocate ( ftab )
    deallocate ( xtab )

  end do

  return
end
subroutine test29 ( )

!*****************************************************************************80
!
!! TEST29 tests AVINT and PLINT.
!
!  Modified:
!
!    10 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ntab = 25

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), dimension (ntab) :: ftab = (/ &
    1.01D-03, &
    7.88D-04, &
    5.63D-04, &
    3.38D-04, &
    2.25D-04, &
    1.13D-04, &
    5.63D-05, &
    2.25D-05, &
    1.13D-05, &
    5.63D-06, &
    2.25D-06, &
    1.13D-06, &
    5.63D-07, &
    2.25D-07, &
    1.13D-07, &
    5.63D-08, &
    1.81D-08, &
    1.13D-08, &
    5.63D-09, &
    2.25D-09, &
    1.13D-09, &
    5.63D-10, &
    2.25D-10, &
    1.13D-10, &
    5.63D-11 /)
  integer ( kind = 4 ) i
  real ( kind = 8 ) result1
  real ( kind = 8 ) result2
  real ( kind = 8 ), dimension (ntab) :: xtab = (/ &
    9.821221D+00, &
    64.28147D+00, &
    216.5835D+00, &
    730.615D+00, &
    1638.551D+00, &
    5870.479D+00, &
    16380.7D+00, &
    52411.82D+00, &
    117519.1D+00, &
    255870.0D+00, &
    705120.8D+00, &
    1525173.0D+00, &
    3339822.0D+00, &
    9766004.0D+00, &
    2.30D+07, &
    5.84D+07, &
    5.21D+08, &
    7.78D+08, &
    1.02D+09, &
    1.28D+09, &
    1.44D+09, &
    1.58D+09, &
    1.75D+09, &
    1.87D+09, &
    1.98D+09 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST29'
  write ( *, '(a)' ) '  Compare PLINT and AVINT on badly scaled data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Y values are all POSITIVE, so we expect'
  write ( *, '(a)' ) '  every partial integral to be positive as well.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      B                PLINT           AVINT'
  write ( *, '(a)' ) ' '

  do i = 3, ntab

    a = xtab(1)
    b = xtab(i)

    call plint ( ntab, xtab, ftab, a, b, result1 )

    call avint ( ntab, xtab, ftab, a, b, result2 )

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) b, result1, result2

  end do


  return
end
function f1d1 ( x )

!*****************************************************************************80
!
!! F1D1(X) = 1.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f1d1
  real ( kind = 8 ) x

  f1d1 = 1.0D+00

  return
end
function f1d2 ( x, y )

!*****************************************************************************80
!
!! F1D2(X,Y) = 1.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f1d2
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  f1d2 = 1.0D+00

  return
end
function fbd1 ( x )

!*****************************************************************************80
!
!! FBD1(X) = 1 / ( 1 + X**2 )
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fbd1
  real ( kind = 8 ) x

  fbd1 = 1.0D+00 / ( 1.0D+00 + x**2 )

  return
end
function fchby ( x )

!*****************************************************************************80
!
!! FCHBY ...
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fchby
  real ( kind = 8 ) g
  integer ( kind = 4 ) iprob
  real ( kind = 8 ) x

  common /problm/ iprob

  if ( iprob == 1 )then
    g = 1.0D+00
  else
    g = x**(iprob-1)
  end if

  fchby = g / sqrt ( 1.0D+00 - x * x )

  return
end
function fed1 ( x )

!*****************************************************************************80
!
!! FED1(X) = EXP(X).
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fed1
  real ( kind = 8 ) x

  fed1 = exp ( x )

  return
end
function fqd1 ( x )

!*****************************************************************************80
!
!! FQD1(X) = SQRT ( X ).
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fqd1
  real ( kind = 8 ) x

  fqd1 = sqrt ( abs ( x ) )

  return
end
function fxd1 ( x )

!*****************************************************************************80
!
!! FXD1(X) = X.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fxd1
  real ( kind = 8 ) x

  fxd1 = x

  return
end
function fx2d1 ( x )

!*****************************************************************************80
!
!! FX2D1(X) = X**2.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx2d1
  real ( kind = 8 ) x

  fx2d1 = x**2

  return
end
function fx3d1 ( x )

!*****************************************************************************80
!
!! FX3D1(X) = X**3.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx3d1
  real ( kind = 8 ) x

  fx3d1 = x**3

  return
end
function f1sd1 ( x, y )

!*****************************************************************************80
!
!! F1SD1(X,Y) = 1 / SQRT ( 1 - X**2 )
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f1sd1
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  f1sd1 = 1.0D+00 / sqrt ( 1.0D+00 - x**2)

  return
end
function fxd2 ( x, y )

!*****************************************************************************80
!
!! FXD2(X,Y) = X + Y.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fxd2
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  fxd2 = x + y

  return
end
function fx2d2 ( x, y )

!*****************************************************************************80
!
!! FX2D2(X,Y) = X**2 + Y**2
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx2d2
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  fx2d2 = x**2 + y**2

  return
end
function fx3d2 ( x, y )

!*****************************************************************************80
!
!! FX3D2(X,Y) = X**3 + Y**3
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx3d2
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  fx3d2 = x**3 + y**3

  return
end
function fxsd1 ( x )

!*****************************************************************************80
!
!! FXSD1(X) = X / SQRT ( 1 - X**2 )
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fxsd1
  real ( kind = 8 ) x

  fxsd1 = x / sqrt ( 1.0D+00 - x**2 )

  return
end
function fx2sd1 ( x, y )

!*****************************************************************************80
!
!! FX2SD1(X,Y) = X**2 / SQRT ( 1 - X**2 )
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx2sd1
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  fx2sd1 = x**2 / sqrt ( 1.0D+00 - x**2 )

  return
end
function fed2 ( x, y )

!*****************************************************************************80
!
!! FED2(X,Y) = EXP(X+Y).
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fed2
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  fed2 = exp ( x + y )

  return
end
function fbd2 ( x, y )

!*****************************************************************************80
!
!! FBD2(X,Y) = 1 / ( 1 + X**2 + Y**2 )
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fbd2
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  fbd2 = 1.0D+00 / ( 1.0D+00 + x**2 + y**2 )

  return
end
