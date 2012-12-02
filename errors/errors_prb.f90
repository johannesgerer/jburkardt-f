program main

!*****************************************************************************80
!
!! MAIN is the main program for ERRORS_PRB.
!
!  Discussion:
!
!    ERRORS_PRB runs the ERRORS examples.
!
!    This program demonstrates how unreliable certain
!    calculations can be when performed with a finite accuracy
!    on a computer.  Standard software is used, where possible,
!    for the solution of these problems.
!
!  Quoting from the reference:
!
!    Despite the increasing spread of computers and calculators,
!    many users still suffer from a certain computational naivete.
!    Whatever the computer produces as a result is regarded as
!    almost mathematically proved.  It is often quite astonishing
!    for a user to discover that even in simple calculations with
!    just a few operations it is possible to produce completely
!    incorrect results, and that this is a necessary consequence
!    of current computer techniques (not just the user's method,
!    but the arithmetic unit of the computer as well).
!
!    The following examples are offered to show how few operations
!    it takes to cause a computer to return meaningless results.
!    This fact applies not just to hand calculators
!    but also to the biggest supercomputers.  The particular
!    problem involved in these effects is that one never knows
!    whether, when or where these sorts of errors
!    will enter into a long, involved calculation and destroy the
!    results.
!
!    Only in the last few years have methods and solution
!    processes been developed which can mathematically describe
!    these effects, and numerically control them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    U Kulisch, C Ullrich, Editors,
!    Wissenschaftliches Rechnen und Programmiersprachen,
!    (Scientific Computing and Programming Languages),
!    Berichte des German Chapter of the ACM,
!    (Reports of the German Chapter of the ACM)
!    Volume 10, Teubner Verlag, 1982.
!
!    Yves Nievergelt,
!    Numerical Linear Algebra on the HP-28, or How to Lie with Supercalculators,
!    The American Mathematical Monthly,
!    Volume 98, Number 6, June-July 1991, pages 539-544.
!
!    S M Rump,
!    Wie Zuverlaessig Sind die Ergebnisse Unserer Rechenanlagen?
!    (How Reliable are the Results of our Computations?)
!    Jahrbuch Ueberblicke Mathematik 1983, pages 163-168.
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ERRORS_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ERRORS library.'

  call header ( )

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test065 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )
  call test10 ( )

  call test11 ( )
  call test12 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
  call test16 ( )
  call test17 ( )
  call test18 ( )
  call test19 ( )
  call test20 ( )

  call test21 ( )
  call test215 ( )
  call test22 ( )
  call test23 ( )
  call test24 ( )
  call test25 ( )
  call test26 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ERRORS_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine header ( )

!*****************************************************************************80
!
!! HEADER prints out a header.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HOW RELIABLE ARE OUR CALCULATIONS?'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'This program demonstrates how unreliable certain'
  write ( *, '(a)' ) 'calculations are, when performed with a finite'
  write ( *, '(a)' ) 'accuracy on a computer. '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Standard software is used where possible for the '
  write ( *, '(a)' ) 'solution of these problems.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Paraphrasing the reference:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' "Despite the spread of computers and calculators,'
  write ( *, '(a)' ) '  many users suffer a certain computational '
  write ( *, '(a)' ) '  gullibility.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Whatever the computer produces is regarded as'
  write ( *, '(a)' ) '  almost mathematically proven.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  It is often quite astonishing to discover that'
  write ( *, '(a)' ) '  even in simple calculations with just a few '
  write ( *, '(a)' ) '  operations it is possible to produce completely '
  write ( *, '(a)' ) '  incorrect results, and that this is a necessary '
  write ( *, '(a)' ) '  consequence of current computer techniques - not '
  write ( *, '(a)' ) '  just the user''s method, but the arithmetic '
  write ( *, '(a)' ) '  processes of the computer as well.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The following examples show how few operations '
  write ( *, '(a)' ) '  are needed for a computer to return meaningless '
  write ( *, '(a)' ) '  results.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This fact applies not just to hand calculators'
  write ( *, '(a)' ) '  but also to the biggest supercomputers.  The '
  write ( *, '(a)' ) '  worst aspect of these effects is that one never'
  write ( *, '(a)' ) '  knows whether, when or where these sorts of '
  write ( *, '(a)' ) '  errors will enter into a long, involved '
  write ( *, '(a)' ) '  calculation and destroy the results.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Only in the last few years have methods and '
  write ( *, '(a)' ) '  solution processes been developed to '
  write ( *, '(a)' ) '  mathematically describe these effects, and '
  write ( *, '(a)' ) '  numerically control them."'

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), parameter :: dp = 665857.0D+00
  real ( kind = 8 ), parameter :: dq = 470832.0D+00
  real ( kind = 8 ) dr
  real ( kind = 4 ), parameter :: p = 665857.0E+00
  real ( kind = 4 ), parameter :: q = 470832.0E+00
  real ( kind = 4 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 1                                 |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    P = 665857                               |'
  write ( *, '(a)' ) '|    Q = 470832                               |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Compute:                                   |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    R = P**2 - 2 * Q**2                      |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct value:                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    R = 1                                    |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '

  r = p**2 - 2.0E+00 * q**2

  write ( *, '(a,g14.6)' ) '  Single precision R = ', r

  dr = dp**2 - 2.0D+00 * dq**2

  write ( *, '(a,g14.6)' ) '  Double precision R = ', dr

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), parameter :: dp = 10864.0D+00
  real ( kind = 8 ), parameter :: dq = 18817.0D+00
  real ( kind = 8 ) dr
  real ( kind = 4 ), parameter :: p = 10864.0E+00
  real ( kind = 4 ), parameter :: q = 18817.0E+00
  real ( kind = 4 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 2                                 |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    P = 10864                                |'
  write ( *, '(a)' ) '|    Q = 18817                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Compute:                                   |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    R = 9 * P**4 - Q**4 + 2 * Q**2           |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct value:                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    R = 1                                    |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '

  r = 9.0E+00 * p**4 - q**4 + 2.0E+00 * q**2

  write ( *, '(a,g14.6)' ) '  Single precision R = ', r

  dr = 9.0D+00 * dp**4 - dq**4 + 2.0D+00 * dq**2

  write ( *, '(a,g14.6)' ) '  Double precision R = ', dr

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) dp
  real ( kind = 8 ) dq
  real ( kind = 8 ) dr
  real ( kind = 4 ) p
  real ( kind = 4 ) q
  real ( kind = 4 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 3                                 |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    P = 10**34                               |'
  write ( *, '(a)' ) '|    Q = - 2                                  |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Compute:                                   |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    R = ( P + Q ) - P                        |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct value:                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    R = - 2                                  |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '

  p = 1.0E+34
  q = - 2.0E+00

  r = ( p + q ) - p

  write ( *, '(a,g14.6)' ) '  Single precision R = ', r

  dp = 1.0D+34
  dq = - 2.0D+00

  dr = ( dp + dq ) - dp

  write ( *, '(a,g14.6)' ) '  Double precision R = ', dr

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 4 ) c(0:n)
  real ( kind = 8 ) dc(0:n)
  real ( kind = 8 ) dq
  real ( kind = 8 ) dr
  real ( kind = 4 ) q
  real ( kind = 4 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 4                                 |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    Q = 1.091608                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Compute:                                   |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    R =                                      |'
  write ( *, '(a)' ) '|        170.4  * Q**3                        |'
  write ( *, '(a)' ) '|      - 356.41 * Q**2                        |'
  write ( *, '(a)' ) '|      + 168.97 * Q                           |'
  write ( *, '(a)' ) '|      +  18.601                              |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct value:                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    R = 0.821248E-13                         |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '

  q = 1.091608E+00

  c(0) =    18.601E+00
  c(1) =   168.97E+00
  c(2) = - 356.41E+00
  c(3) =   170.4E+00

  call rpoly_val ( n, c, q, r )

  write ( *, '(a,g14.6)' ) '  Single precision R (direct evaluation) = ', r

  call rpoly_val_horner ( n, c, q, r )

  write ( *, '(a,g14.6)' ) '  Single precision R (Horner''s method) =  ', r

  dq = 1.091608D+00

  dc(0) =    18.601D+00
  dc(1) =   168.97D+00
  dc(2) = - 356.41D+00
  dc(3) =   170.4D+00

  call dpoly_val ( n, dc, dq, dr )

  write ( *, '(a,g14.6)' ) '  Double precision R (direct evaluation) = ', dr

  call dpoly_val_horner ( n, dc, dq, dr )

  write ( *, '(a,g14.6)' ) '  Double precision R (Horner''s method) =  ', dr

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 4 ) c(0:n)
  real ( kind = 4 ) q
  real ( kind = 4 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 5                                 |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    Q = 0.707107                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Compute:                                   |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    R =                                      |'
  write ( *, '(a)' ) '|         8118 * Q**4                         |'
  write ( *, '(a)' ) '|      - 11482 * Q**3                         |'
  write ( *, '(a)' ) '|      +         Q**2                         |'
  write ( *, '(a)' ) '|      +  5741 * Q                            |'
  write ( *, '(a)' ) '|      -  2030                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct value:                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    R = - 0.191527325270E-10                 |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '

  q = 0.707107E+00

  c(0) =  - 2030.0E+00
  c(1) =    5741.0E+00
  c(2) =       1.0E+00
  c(3) = - 11482.0E+00
  c(4) =    8118.0E+00

  call rpoly_val ( n, c, q, r )

  write ( *, '(a,g14.6)' ) '  R (direct evaluation) = ', r

  call rpoly_val_horner ( n, c, q, r )

  write ( *, '(a,g14.6)' ) '  R (Horner''s method) =  ', r

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) ab(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot(n)
  integer ( kind = 4 ) job
  real ( kind = 4 ) qraux(n)
  real ( kind = 4 ) qty(n)
  real ( kind = 4 ) qy(n)
  real ( kind = 4 ) r(n)
  real ( kind = 4 ) rsd(n)
  real ( kind = 4 ) work(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 6                                 |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    A = [ 64919121    - 159018721 ]          |'
  write ( *, '(a)' ) '|        [ 41869520.5  - 102558961 ]          |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    B = [ 1 ]                                |'
  write ( *, '(a)' ) '|        [ 0 ]                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Solve:                                     |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    A * X = B                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct X:                                 |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    205,117,922.0                            |'
  write ( *, '(a)' ) '|     83,739,041.0                            |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R1: Gauss elimination using SGEFA and SGESL'
  write ( *, '(a)' ) '  R2: QR factorization using SQRDC and SQRSL.'

  a(1,1) =   64919121.0E+00
  a(1,2) = -159018721.0E+00
  a(2,1) =   41869520.5E+00
  a(2,2) = -102558961.0E+00

  r(1) = 1.0E+00
  r(2) = 0.0E+00

  call sgefa ( a, lda, n, ipivot, info )

  if ( info == 0 ) then

    job = 0

    call sgesl ( a, lda, n, ipivot, r, job )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Gauss elimination solution:'
    write ( *, '(a)' ) ' '

    do i = 1, n
      write ( *, '(i4,g14.6)' ) i, r(i)
    end do

  else

    write ( *, '(a,i4)' ) '  SGEFA failed, returning INFO = ', info

  end if

  a(1,1) =   64919121.0E+00
  a(1,2) = -159018721.0E+00
  a(2,1) =   41869520.5E+00
  a(2,2) = -102558961.0E+00

  r(1) = 1.0E+00
  r(2) = 0.0E+00

  job = 0

  call sqrdc ( a, lda, n, n, qraux, ipivot, work, job )

  job = 110

  call sqrsl ( a, lda, n, n, qraux, r, qy, qty, r, rsd, ab, job, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Least Squares solution:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(g14.6)' ) r(i)
  end do

  return
end
subroutine test065 ( )

!*****************************************************************************80
!
!! TEST065
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) a_save(n,n)
  real ( kind = 4 ) ab(n)
  real ( kind = 4 ) det(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot(n)
  integer ( kind = 4 ) job
  real ( kind = 4 ) qraux(n)
  real ( kind = 4 ) qty(n)
  real ( kind = 4 ) qy(n)
  real ( kind = 4 ) resid(n)
  real ( kind = 4 ) rhs_save(n)
  real ( kind = 4 ) rsd(n)
  real ( kind = 4 ) work(n)
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) x_save(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 6.5                               |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    A = [ 888445  887112 ]                   |'
  write ( *, '(a)' ) '|        [ 887112  885781 ]                   |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    det ( A ) = 1                            |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    B = [ 1 ]                                |'
  write ( *, '(a)' ) '|        [ 0 ]                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Solve:                                     |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    A * X = B                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct X:                                 |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    885781                                   |'
  write ( *, '(a)' ) '|   -887112                                   |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R1: Gauss elimination using SGEFA and SGESL'
  write ( *, '(a)' ) '  R2: QR factorization using SQRDC and SQRSL.'
  write ( *, '(a)' ) ' '

  a_save(1,1) = 888445.0E+00
  a_save(1,2) = 887112.0E+00
  a_save(2,1) = 887112.0E+00
  a_save(2,2) = 885781.0E+00

  rhs_save(1:2) = (/ 1.0E+00, 0.0E+00 /)

  x_save(1:2) = (/ 885781.0E+00, -887112.0E+00 /)

  a(1:2,1:2) = a_save(1:2,1:2)
  x(1:2) = rhs_save(1:2)

  call sgefa ( a, lda, n, ipivot, info )

  if ( info == 0 ) then

    job = 0

    call sgesl ( a, lda, n, ipivot, x, job )

    resid(1:n) = matmul ( a_save(1:n,1:n), x(1:n) ) - rhs_save(1:n)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Gauss elimination solution and residual:'
    write ( *, '(a)' ) ' '

    do i = 1, n
      write ( *, '(i4,2g14.6)' ) i, x(i), resid(i)
    end do

    job = 10

    call sgedi ( a, lda, n, ipivot, det, work, job )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6,a,g14.6)' ) 'SGEDI estimates determinant = ', &
      det(1), ' * 10** ', det(2)

  else

    write ( *, '(a,i4)' ) '  SGEFA failed, returning INFO = ', info

  end if

  a(1:2,1:2) = a_save(1:2,1:2)
  x(1:2) = rhs_save(1:2)

  job = 0

  call sqrdc ( a, lda, n, n, qraux, ipivot, work, job )

  job = 110

  call sqrsl ( a, lda, n, n, qraux, x, qy, qty, x, rsd, ab, job, info )

  resid(1:n) = matmul ( a_save(1:n,1:n), x(1:n) ) - rhs_save(1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Least Squares solution and residual:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(i4,2g14.6)' ) i, x(i), resid(i)
  end do

  resid(1:n) = matmul ( a_save(1:n,1:n), x_save(1:n) ) - rhs_save(1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact solution and residual:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(i4,2g14.6)' ) i, x_save(i), resid(i)
  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) ab(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot(n)
  integer ( kind = 4 ) job
  real ( kind = 4 ) qraux(n)
  real ( kind = 4 ) qty(n)
  real ( kind = 4 ) qy(n)
  real ( kind = 4 ) r(n)
  real ( kind = 4 ) rsd(n)
  real ( kind = 4 ) work(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 7                                 |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    A = [ -367296  -43199   519436  -954302 ]|'
  write ( *, '(a)' ) '|        [  259718 -477151  -367295 -1043199 ]|'
  write ( *, '(a)' ) '|        [  886731   88897 -1254026 -1132096 ]|'
  write ( *, '(a)' ) '|        [  627013  566048  -886732   911103 ]|'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    B = [ 1 ]                                |'
  write ( *, '(a)' ) '|        [ 1 ]                                |'
  write ( *, '(a)' ) '|        [ 1 ]                                |'
  write ( *, '(a)' ) '|        [ 0 ]                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Solve:                                     |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    A * X = B                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct X:                                 |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    8.86731088897E+17                        |'
  write ( *, '(a)' ) '|    8.86731088897E+11                        |'
  write ( *, '(a)' ) '|    6.27013566048E+17                        |'
  write ( *, '(a)' ) '|    6.27013566048E+11                        |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'

  a(1,1) =   -367296.0E+00
  a(1,2) =    -43199.0E+00
  a(1,3) =    519436.0E+00
  a(1,4) =   -954302.0E+00
  a(2,1) =    259718.0E+00
  a(2,2) =   -477151.0E+00
  a(2,3) =   -367295.0E+00
  a(2,4) =  -1043199.0E+00
  a(3,1) =    886731.0E+00
  a(3,2) =     88897.0E+00
  a(3,3) =  -1254026.0E+00
  a(3,4) =  -1132096.0E+00
  a(4,1) =    627013.0E+00
  a(4,2) =    566048.0E+00
  a(4,3) =   -886732.0E+00
  a(4,4) =    911103.0E+00

  r(1:4) = (/ 1.0E+00, 1.0E+00, 1.0E+00, 0.0E+00 /)

  call sgefa ( a, lda, n, ipivot, info )

  if ( info /= 0 ) then
    write ( *, '(a,i4)' ) '  SGEFA failed, returning INFO = ', info
  else
    job = 0

    call sgesl ( a, lda, n, ipivot, r, job )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  SGEFA/SGESL solution (Gauss elimination):'
    write ( *, '(a)' ) ' '

    do i = 1, n
      write ( *, '(g14.6)' ) r(i)
    end do

  end if

  a(1,1) =  -367296.0E+00
  a(1,2) =   -43199.0E+00
  a(1,3) =   519436.0E+00
  a(1,4) =  -954302.0E+00
  a(2,1) =   259718.0E+00
  a(2,2) =  -477151.0E+00
  a(2,3) =  -367295.0E+00
  a(2,4) = -1043199.0E+00
  a(3,1) =   886731.0E+00
  a(3,2) =    88897.0E+00
  a(3,3) = -1254026.0E+00
  a(3,4) = -1132096.0E+00
  a(4,1) =   627013.0E+00
  a(4,2) =   566048.0E+00
  a(4,3) =  -886732.0E+00
  a(4,4) =   911103.0E+00

  r(1:4) = (/ 1.0E+00, 1.0E+00, 1.0E+00, 0.0E+00 /)

  job = 0

  call sqrdc ( a, lda, n, n, qraux, ipivot, work, job )

  job = 110

  call sqrsl ( a, lda, n, n, qraux, r, qy, qty, r, rsd, ab, job, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SQRDC/SQRSL solution (QR factorization):'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(g14.6)' ) r(i)
  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  real ( kind = 4 ) f1
  real ( kind = 4 ) f2
  real ( kind = 4 ) f3
  real ( kind = 4 ) f4
  real ( kind = 4 ) fmin
  real ( kind = 4 ), external :: p8horn
  real ( kind = 4 ), external :: p8reg
  real ( kind = 4 ) r1
  real ( kind = 4 ) r2
  real ( kind = 4 ) r3
  real ( kind = 4 ) r4
  real ( kind = 4 ), parameter :: tol = 1.0E-14

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 8                                 |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    P(X) =                                   |'
  write ( *, '(a)' ) '|        170.4  * x**3                        |'
  write ( *, '(a)' ) '|      - 356.41 * x**2                        |'
  write ( *, '(a)' ) '|      + 168.97 * x                           |'
  write ( *, '(a)' ) '|      +  18.601                              |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Seek the minimizer X* of P(X)              |'
  write ( *, '(a)' ) '|  over the interval 0 <= X <= 2.             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  The true minimizer lies in the interval    |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  [1.091607978, 1.091607981 ]                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'

  a = 0.0E+00
  b = 2.0E+00

  r1 = fmin ( a, b, p8reg, tol )
  f1 = p8reg ( r1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Evaluating P(X) in the usual way'
  write ( *, '(a,g14.6)' ) '  FMIN computes the minimizer at ', r1
  write ( *, '(a,g14.6)' ) '  with function value ', f1

  r2 = fmin ( a, b, p8horn, tol )
  f2 = p8horn ( r2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Evaluating P(X) using Horner''s method,'
  write ( *, '(a,g14.6)' ) '  FMIN computes the minimizer at ',r2
  write ( *, '(a,g14.6)' ) '  with function value ', f2

  r3 = 1.091607978E+00
  r4 = 1.091607981E+00
  f3 = p8horn ( r1 )
  f4 = p8horn ( r2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g24.16,a,g24.16)' ) '  The true value lies between ', r3, &
    ' and ', r4
  write ( *, '(a,2g24.16)' ) '  with function value below ', f3, f4

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 4 ), dimension (0:n) :: c = (/ -945804881.0E+00, 1753426039.0E+00, &
    -1083557822.0E+00, 223200658.0E+00 /)
  real ( kind = 8 ), dimension (0:n) :: dc = (/ -945804881.0D+00, &
    1753426039.0D+00, -1083557822.0D+00, 223200658.0D+00 /)
  real ( kind = 8 ) dphorn
  real ( kind = 8 ) dpreg
  real ( kind = 8 ) dx
  integer ( kind = 4 ) i
  real ( kind = 4 ) phorn
  real ( kind = 4 ) preg
  real ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 9                                 |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    P(X) =                                   |'
  write ( *, '(a)' ) '|         223200658 * X**3                    |'
  write ( *, '(a)' ) '|      - 1083557822 * X**2                    |'
  write ( *, '(a)' ) '|      + 1753426039 * X                       |'
  write ( *, '(a)' ) '|       - 945804881                           |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Evaluate at 11 equally spaced values       |'
  write ( *, '(a)' ) '|  from 1.61801916 to 1.61801917              |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct table:                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  ===== X =====    ========= P(X) ========   |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    1.618019160  - 0.17081105112320e-11      |'
  write ( *, '(a)' ) '|    1.618019161  - 0.89804011575510e-12      |'
  write ( *, '(a)' ) '|    1.618019162  - 0.34596536943057e-12      |'
  write ( *, '(a)' ) '|    1.618019163  - 0.51884933054474e-13      |'
  write ( *, '(a)' ) '|    1.618019164  - 0.15797467422848e-13      |'
  write ( *, '(a)' ) '|    1.618019165  - 0.23770163333175e-12      |'
  write ( *, '(a)' ) '|    1.618019166  - 0.71759609157723e-12      |'
  write ( *, '(a)' ) '|    1.618019167  - 0.14554795029553e-11      |'
  write ( *, '(a)' ) '|    1.618019168  - 0.24513505282621e-11      |'
  write ( *, '(a)' ) '|    1.618019169  - 0.37052078282936e-11      |'
  write ( *, '(a)' ) '|    1.618019170  - 0.52170500638460e-11      |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R1: compute P(X) the regular way.'
  write ( *, '(a)' ) '  R2: compute P(X) using Horner''s method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Single precision:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, R1, R2'
  write ( *, '(a)' ) ' '

  do i = 0, 10

    x = ( real ( 10 - i ) * 1.61801916E+00 &
        + real (      i ) * 1.61801917E+00 ) / 10.0E+00

    call rpoly_val ( n, c, x, preg )

    call rpoly_val_horner ( n, c, x, phorn )

    write ( *, '(f14.10,3g14.6)' ) x, preg, phorn

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Double precision:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X, DR1, DR2'
  write ( *, '(a)' ) ' '
  do i = 0, 10

    dx = ( real ( 10 - i, kind = 8 ) * 1.61801916D+00 &
      +  real ( i, kind = 8 ) * 1.61801917D+00 ) / 10.0D+00

    dpreg = dc(3) * dx**3 + dc(2) * dx**2 + dc(1) * dx + dc(0)

    call dpoly_val_horner ( n, dc, dx, dphorn )

    write ( *, '(f14.10,3g14.6)' ) dx, dpreg, dphorn

  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 seeks the nearest straight line to given data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) drb
  real ( kind = 8 ) drm
  real ( kind = 8 ) drm_bot
  real ( kind = 8 ) drm_top
  real ( kind = 8 ) drv
  real ( kind = 8 ) dx
  real ( kind = 8 ) dx1
  real ( kind = 8 ) dx2
  real ( kind = 8 ) dx3
  real ( kind = 8 ) dy1
  real ( kind = 8 ) dy2
  real ( kind = 8 ) dy3
  real ( kind = 4 ) rb
  real ( kind = 4 ) rm
  real ( kind = 4 ) rm_bot
  real ( kind = 4 ) rm_top
  real ( kind = 4 ) rv
  real ( kind = 4 ) x
  real ( kind = 4 ) x1
  real ( kind = 4 ) x2
  real ( kind = 4 ) x3
  real ( kind = 4 ) y1
  real ( kind = 4 ) y2
  real ( kind = 4 ) y3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 10                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  For the data:                              |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|       X       Y                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    5201477   99999                          |'
  write ( *, '(a)' ) '|    5201478  100000                          |'
  write ( *, '(a)' ) '|    5201479  100001                          |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Find the "nearest" straight line:          |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    Y(X) = A * X + B                         |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  in the least squares sense.                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct values:                            |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    A = 1.0                                  |'
  write ( *, '(a)' ) '|    B = - 5101478                            |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Evaluate:                                  |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    Y(5201480)                               |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct value:                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    Y(5201480) = 100002.0                    |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                    Slope A    Intercept B    Y(5201480) '
  write ( *, '(a)' ) ' '

  x1 = 5201477.0E+00
  x2 = 5201478.0E+00
  x3 = 5201479.0E+00

  y1 =  99999.0E+00
  y2 = 100000.0E+00
  y3 = 100001.0E+00

  rm_top = ( x1 * y1 + x2 * y2 + x3 * y3 &
         - ( x1 + x2 + x3 ) * ( y1 + y2 + y3 ) / 3.0E+00 )
  rm_bot = ( x1**2 + x2**2 + x3**2 - ( x1 + x2 + x3 )**2 / 3.0E+00 )

  if ( rm_bot == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Single precision slope calculation fails.'
    write ( *, '(a)' ) '  The divisor was computed as 0.'
  else

    rm = rm_top / rm_bot

    rb = ( y1 + y2 + y3 ) / 3.0E+00 - rm * ( x1 + x2 + x3 ) / 3.0E+00

    x = 5201480.0E+00
    rv = rm * x + rb

    write ( *, '(a,3g14.6)' ) '  Single precision ', rm, rb, rv

  end if

  dx1 = 5201477.0D+00
  dx2 = 5201478.0D+00
  dx3 = 5201479.0D+00

  dy1 = 99999.0D+00
  dy2 = 100000.0D+00
  dy3 = 100001.0D+00

  drm_top = ( dx1 * dy1 + dx2 * dy2 + dx3 * dy3 &
          - ( dx1 + dx2 + dx3 ) * ( dy1 + dy2 + dy3 ) / 3.0D+00 )
  drm_bot = ( dx1**2 + dx2**2 + dx3**2 - ( dx1 + dx2 + dx3 )**2 / 3.0D+00 )

  if ( drm_bot == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Double precision slope calculation fails.'
    write ( *, '(a)' ) '  The divisor was computed as 0.'
  else
    drm = drm_top / drm_bot
    drb = ( dy1 + dy2 + dy3 ) / 3.0D+00 - drm * ( dx1 + dx2 + dx3 ) / 3.0D+00

    dx = 5201480.0D+00
    drv = drm * dx + drb

    write ( *, '(a,3g14.6)' ) '  Double precision ', drm, drb, drv
  end if

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  real ( kind = 4 ) f1
  real ( kind = 4 ) f2
  real ( kind = 4 ), external :: fhorn11
  real ( kind = 4 ) fmin
  real ( kind = 4 ), external :: freg11
  real ( kind = 4 ) r1
  real ( kind = 4 ) r2
  real ( kind = 4 ) :: tol = 1.0E-14

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 11                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  P(X) =                                     |'
  write ( *, '(a)' ) '|      2124476931 * X**4                      |'
  write ( *, '(a)' ) '|    - 1226567328 * X**3                      |'
  write ( *, '(a)' ) '|    -  708158977 * X**2                      |'
  write ( *, '(a)' ) '|    +  408855776 * X                         |'
  write ( *, '(a)' ) '|    +    1.0E-27                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  P(X) is positive for X >= 0.               |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Use a minimizer to seek the minimum        |'
  write ( *, '(a)' ) '|  value of P(X) in (0.56, 0.59).             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '

  a = 0.56E+00
  b = 0.59E+00

  r1 = fmin ( a, b, freg11, tol )
  f1 = freg11 ( r1 )

  r2 = fmin ( a, b, fhorn11, tol )
  f2 = fhorn11 ( r2 )

  write ( *, '(a)' ) '  Using standard evaluation,'
  write ( *, '(a,g18.10)' ) '  subroutine FMIN finds a minimizer at X = ', r1
  write ( *, '(a,g18.10)' ) '  with function value P(X) = ', f1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Using Horner''s rule,'
  write ( *, '(a,g18.10)' ) '  subroutine FMIN finds a minimizer at X = ', r2
  write ( *, '(a,g18.10)' ) '  with function value P(X) = ', f2

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 9
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) ab(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lcm_12n
  integer ( kind = 4 ) mult
  real ( kind = 4 ) qraux(n)
  real ( kind = 4 ) qty(n)
  real ( kind = 4 ) qy(n)
  real ( kind = 4 ) r1(n)
  real ( kind = 4 ) r2(n)
  real ( kind = 4 ) rsd(n)
  real ( kind = 4 ) work(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 12                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Solve a linear system A * x = b            |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  A is the Hilbert matrix of order 9,        |'
  write ( *, '(a)' ) '|  multiplied by the least common multiple of |'
  write ( *, '(a)' ) '|  1 through 9, so all entries are integers.  |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  B = (1,0,0,...,0)                          |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  We are interested in the first component   |'
  write ( *, '(a)' ) '|  X(1) of the solution vector.               |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct X(1):                              |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    6.611036022800E-06                       |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R1: Gauss elimination using SGEFA and SGESL'
  write ( *, '(a)' ) '  R2: QR factorization using SQRDC and SQRSL.'
  write ( *, '(a)' ) ' '

  mult = lcm_12n ( n )

  do i = 1, n
    do j = 1, n
      a(i,j) = real ( mult ) / real ( i + j - 1 )
    end do
  end do

  r1(1) = 1.0E+00
  r1(2:n) = 0.0E+00

  call sgefa ( a, lda, n, ipivot, info )

  if ( info == 0 ) then

    job = 0

    call sgesl ( a, lda, n, ipivot, r1, job )

  else

    write ( *, '(a,i4)' ) '  SGEFA failed, returning INFO = ', info

  end if

  do i = 1, n
    do j = 1, n
      a(i,j) = real ( mult ) / real ( i + j - 1 )
    end do
  end do

  r2(1) = 1.0E+00
  r2(2:n) = 0.0E+00

  job = 0

  call sqrdc ( a, lda, n, n, qraux, ipivot, work, job )

  job = 110

  call sqrsl ( a, lda, n, n, qraux, r2, qy, qty, r2, rsd, ab, job, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First solution component only:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, R1(I), R2(I)'
  write ( *, '(a)' ) ' '

  write ( *, '(i6,2g14.6)' ) i, r1(1), r2(1)

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 21
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) ab(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) lcm_12n
  integer ( kind = 4 ) mult
  real ( kind = 4 ) qraux(n)
  real ( kind = 4 ) qty(n)
  real ( kind = 4 ) qy(n)
  real ( kind = 4 ) r1(n)
  real ( kind = 4 ) r2(n)
  real ( kind = 4 ) rsd(n)
  real ( kind = 4 ) work(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 13                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Solve a linear system A * x = b            |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  A is the Hilbert matrix of order 21        |'
  write ( *, '(a)' ) '|  multiplied by the least common multiple of |'
  write ( *, '(a)' ) '|  1 through 25 so all entries are integral.  |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  B is (1,0,0,...,0)                         |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  We are interested in the first component   |'
  write ( *, '(a)' ) '|  X(1) of the solution vector.               |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct X(1):                              |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    2.013145339298E-15                       |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R1: Gauss elimination using SGEFA and SGESL'
  write ( *, '(a)' ) '  R2: QR factorization using SQRDC and SQRSL.'
  write ( *, '(a)' ) ' '

  mult = lcm_12n ( n )

  do i = 1, n
    do j = 1, n
      a(i,j) = real ( mult ) / real ( i + j - 1 )
    end do
  end do

  r1(1) = 1.0E+00
  r1(2:n) = 0.0E+00

  call sgefa ( a, lda, n, ipivot, info )

  if ( info == 0 )then

    job = 0

    call sgesl ( a, lda, n, ipivot, r1, job )

  else

    write ( *, '(a,i4)' ) '  SGEFA failed, returning INFO = ', info

  end if

  do i = 1, n
    do j = 1, n
      a(i,j) = real ( mult ) / real ( i + j - 1 )
    end do
  end do

  r2(1) = 1.0E+00
  r2(2:n) = 0.0E+00

  job = 0

  call sqrdc ( a, lda, n, n, qraux, ipivot, work, job )

  job = 110

  call sqrsl ( a, lda, n, n, qraux, r2, qy, qty, r2, rsd, ab, job, info )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First solution component only:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, R1(I), R2(I)'
  write ( *, '(a)' ) ' '
  write ( *, '(i6,2g14.6)' )   i, r1(1), r2(1)

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2
  integer ( kind = 4 ), parameter :: lda = n

  real ( kind = 4 ) a(lda,n)
  real ( kind = 4 ) ab(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot(n)
  integer ( kind = 4 ) job
  real ( kind = 4 ) qraux(n)
  real ( kind = 4 ) qty(n)
  real ( kind = 4 ) qy(n)
  real ( kind = 4 ) r1(n)
  real ( kind = 4 ) r2(n)
  real ( kind = 4 ) rsd(n)
  real ( kind = 4 ) work(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 14                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    A = [ 64079  57314 ]                     |'
  write ( *, '(a)' ) '|        [ 51860  46385 ]                     |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    B = [   2 ]                              |'
  write ( *, '(a)' ) '|        [ 305 ]                              |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Solve:                                     |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    A * X = B                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct X:                                 |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    - 46368.0                                |'
  write ( *, '(a)' ) '|      51841.0                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '            X(1)      X(2)'
  write ( *, '(a)' ) ' '

  a(1,1) = 64079.0E+00
  a(1,2) = 57314.0E+00
  a(2,1) = 51860.0E+00
  a(2,2) = 46385.0E+00

  r1(1) = 2.0E+00
  r1(2) = 305.0E+00

  call sgefa ( a, lda, n, ipivot, info )

  if ( info == 0 ) then

    job = 0

    call sgesl ( a, lda, n, ipivot, r1, job )

    write ( *, '(a,2g14.6)' ) '  SGEFA/SGESL  ', r1(1), r1(2)

  else

    write ( *, '(a)' ) '  SGEFA failed.'

  end if

  a(1,1) = 64079.0E+00
  a(1,2) = 57314.0E+00
  a(2,1) = 51860.0E+00
  a(2,2) = 46385.0E+00

  r2(1) = 2.0E+00
  r2(2) = 305.0E+00

  job = 0

  call sqrdc ( a, lda, n, n, qraux, ipivot, work, job )

  job = 110

  call sqrsl ( a, lda, n, n, qraux, r2, qy, qty, r2, rsd, ab, job, info )

  write ( *, '(a,2g14.6)' ) '  SQRDC/SQRSL  ', r2(1), r2(2)

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) fx15
  real ( kind = 4 ) h
  integer ( kind = 4 ) i
  real ( kind = 4 ) r
  real ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 15                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    F(X) =       ( 4970 * X - 4923 )         |'
  write ( *, '(a)' ) '|           --------------------------------  |'
  write ( *, '(a)' ) '|           ( 4970 * X**2 - 9799 * X + 4830 ) |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Approximate F"(X) at X = 1 by:             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    Del2(X)(F,H) =                           |'
  write ( *, '(a)' ) '|      ( F(X-H) - 2*F(X) + F(X+H) ) / H**2    |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct values:                            |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    F"(1.0)             = 94.0               |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    Del2(1.0)(F,1.0E-4) = 70.78819           |'
  write ( *, '(a)' ) '|    Del2(1.0)(F,1.0E-5) = 93.76790           |'
  write ( *, '(a)' ) '|    Del2(1.0)(F,1.0E-8) = 94.0               |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  H, Del2(X,H)'

  x = 1.0E+00

  do i = 1, 3

    if ( i == 1 ) then
      h = 1.0E-04
    else if ( i == 2 ) then
      h = 1.0E-05
    else if ( i == 3 ) then
      h = 1.0E-08
    end if

    r = ( fx15 ( x + h ) - 2.0E+00 * fx15 ( x ) + fx15 ( x - h ) ) / h**2

    write ( *, '(2g14.6)' ) h, r

  end do

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) scale
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  real ( kind = 4 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 16                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  P(X,Y) =                                   |'
  write ( *, '(a)' ) '|     83521 *        y**8                     |'
  write ( *, '(a)' ) '|    +  578 * x**2 * y**4                     |'
  write ( *, '(a)' ) '|    -    2 * x**4                            |'
  write ( *, '(a)' ) '|    +    2 * x**6                            |'
  write ( *, '(a)' ) '|    -        x**8                            |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Evaluate P(X,Y) for                        |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    X = 9478657                              |'
  write ( *, '(a)' ) '|    Y = 2298912                              |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct value:                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    P(X,Y) = - 179,689,877,047,297           |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'

  scale = 1000000.0E+00
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The unscaled calculation will overflow.'
  write ( *, '(a,g14.6)' ) '  We use a scale factor on X and Y of ', scale

  x = 9478657.0E+00 / scale
  y = 2298912.0E+00 / scale

  z =     83521.0E+00 *        y**8 &
        +   578.0E+00 * x**2 * y**4 / scale**2 &
        -     2.0E+00 * x**4 / scale**4 &
        +     2.0E+00 * x**6 * scale**2 &
        -           x**8
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Scaled value of P(X,Y) = ', z
  write ( *, '(a,g14.6,a,g14.6,a)' ) '  Computed value of P(X,Y) = ', z , &
    ' * ', scale, ' ** 8'

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 4 ), dimension ( n ) :: a = (/ 2.718281828E+00, -3.141592654E+00, &
    1.414213562E+00, 0.5772156649E+00, 0.3010299957E+00 /)
  real ( kind = 4 ) b(n)
  real ( kind = 8 ), dimension ( n ) :: da = (/ 2.718281828D+00, &
    -3.141592654D+00, 1.414213562D+00, 0.5772156649D+00, &
    0.3010299957D+00 /)
  real ( kind = 8 ) db(n)
  real ( kind = 8 ) dr
  integer ( kind = 4 ) i
  real ( kind = 4 ) r
  real ( kind = 4 ) sdot
  real ( kind = 4 ) sdsdot

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 17                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  A =                                        |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|     2.718281828                             |'
  write ( *, '(a)' ) '|    -3.141592654                             |'
  write ( *, '(a)' ) '|     1.414213562                             |'
  write ( *, '(a)' ) '|     0.5772156649                            |'
  write ( *, '(a)' ) '|     0.3010299957                            |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  B =                                        |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|        1486.2497                            |'
  write ( *, '(a)' ) '|      878366.9879                            |'
  write ( *, '(a)' ) '|        - 22.37492                           |'
  write ( *, '(a)' ) '|     4773714.647                             |'
  write ( *, '(a)' ) '|           0.000185049                       |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Compute the scalar product:                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    R = A dot B                              |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct value:                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    R = - 0.100657107E-10                    |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '

  b(1) =    1486.2497E+00
  b(2) =  878366.9879E+00
  b(3) =    - 22.37492E+00
  b(4) = 4773714.647E+00
  b(5) =       0.000185049E+00

  r = 0.0E+00
  do i = 1, n
    r = r + a(i) * b(i)
  end do

  write ( *, '(a,g14.6)' ) '  Standard loop     R = ', r

  r = sdot ( n, a, 1, b, 1 )
  write ( *, '(a,g14.6)' ) '  SDOT:             R = ', r

  r = sdsdot ( n, a, 1, b, 1 )
  write ( *, '(a,g14.6)' ) '  SDSDOT:           R = ', r

  r = dot_product ( a, b )
  write ( *, '(a,g14.6)' ) '  DOT_PRODUCT:      R = ', r

  db(1) =    1486.2497D+00
  db(2) =  878366.9879D+00
  db(3) =    - 22.37492D+00
  db(4) = 4773714.647D+00
  db(5) =       0.000185049D+00

  dr = 0.0E+00
  do i = 1, n
    dr = dr + da(i) * db(i)
  end do

  write ( *, '(a,g14.6)' ) '  Standard loop    DR = ', dr
  dr = dot_product ( da, db )

  write ( *, '(a,g14.6)' ) '  DOT_PRODUCT:     DR = ', dr

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) scale
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  real ( kind = 4 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 18                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  X = 192119201                              |'
  write ( *, '(a)' ) '|  Y = 35675640                               |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Z(X,Y) =                                   |'
  write ( *, '(a)' ) '|    (                                        |'
  write ( *, '(a)' ) '|         1682 * x    * y**4                  |'
  write ( *, '(a)' ) '|       +    3 * x**3                         |'
  write ( *, '(a)' ) '|       +   29 * x    * y**2                  |'
  write ( *, '(a)' ) '|       -    2 * x**5                         |'
  write ( *, '(a)' ) '|       +  832                                |'
  write ( *, '(a)' ) '|      ) / 107751                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct Z:                                 |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    1783.0                                   |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'

  scale = 1000000.0E+00
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The unscaled calculation will overflow.'
  write ( *, '(a,g14.6)' ) '  We use a scale factor on X and Y of ', scale

  x = 192119201.0E+00 / scale
  y = 35675640.0E+00 / scale

  z = (      1682.0E+00 * x    * y**4 &
           +    3.0E+00 * x**3 / scale**2 &
           +   29.0E+00 * x    * y**2 / scale**2 &
           -    2.0E+00 * x**5 &
           +  832.0E+00 / scale**5 &
         ) / 107751.0E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Scaled value of Z = ', z
  write ( *, '(a,g14.6,a,g14.6,a,i5)' ) '  Computed value of Z = ', z, &
    ' * ', scale, ' ** ', 5

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 4 ) fx3
  real ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 19                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    F1(X) = 1 - Cos(X)                       |'
  write ( *, '(a)' ) '|    F2(X) = Sin**2(X) / ( 1 + Cos(X) )       |'
  write ( *, '(a)' ) '|    F3(X) = Taylor series for 1 - Cos(X)     |'
  write ( *, '(a)' ) '|      around X = 0, up to X**6 terms.        |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  For all X                                  |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    F1(X) = F2(X)                            |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  For X near 0,                              |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    F3(X) should approximate F1(X)           |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X          1-Cos        Sin**2/(1+Cos)  Taylor'
  write ( *, '(a)' ) ' '

  x = 0.5E+00

  do

    fx1 = 1.0E+00 - cos ( x )
    fx2 = ( sin ( x ) )**2 / ( 1.0E+00 + cos ( x ) )
    fx3 = 0.5E+00 * x**2 * ( 1.0E+00 - ( x**2 / 12.0E+00 ) + ( x**4 / 360.0E+00 ) )

    write ( *, '(4g14.6)' ) x, fx1, fx2, fx3

    x = 0.25E+00 * x

    if ( x <= 1.0E-08 ) then
      exit
    end if

  end do

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2

  real ( kind = 8 ), dimension (0:n) :: dp = &
    (/ 0.005D+00, 100.0D+00, 0.00005D+00 /)
  real ( kind = 8 ) dpval(2)
  complex ( kind = 8 ) dr(2)
  integer ( kind = 4 ) ierror
  real ( kind = 4 ), dimension (0:n) :: p = (/ 0.005E+00, 100.0E+00, 0.00005E+00 /)
  real ( kind = 4 ) pval(2)
  complex ( kind = 4 ) r(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 20                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  P(X) =                                     |'
  write ( *, '(a)' ) '|        0.00005 * x**2                       |'
  write ( *, '(a)' ) '|    + 100        * x                         |'
  write ( *, '(a)' ) '|    +   0.005                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Evaluate the standard quadratic formula:   |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    X1 = 0.5 * ( -B + sqrt(B**2-4*A*C) ) / A |'
  write ( *, '(a)' ) '|    X2 = 0.5 * ( -B - sqrt(B**2-4*A*C) ) / A |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Evaluate an alternate quadratic formula:   |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    X3 = 2 * C / ( - B + sqrt(B**2-4*A*C) )  |'
  write ( *, '(a)' ) '|    X4 = 2 * C / ( - B - sqrt(B**2-4*A*C) )  |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct results:                           |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    Root1 =       -0.00005...                |'
  write ( *, '(a)' ) '|    Root2 = -1999999.99995...                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'

  call rpoly2_roots ( p, r )

  call rpoly_val_horner ( n, p, real ( r(1) ), pval(1) )
  call rpoly_val_horner ( n, p, real ( r(2) ), pval(2) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  X1 = ', real ( r(1) ), &
    ' with P(X1) = ', pval(1)
  write ( *, '(a,g14.6,a,g14.6)' ) '  X2 = ', real ( r(2) ), &
    ' with P(X2) = ', pval(2)

  call rpoly2_roots2 ( p, r, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  RPOLY2_ROOTS2 returned an error flag.'
  else
    call rpoly_val_horner ( n, p, real ( r(1) ), pval(1) )
    call rpoly_val_horner ( n, p, real ( r(2) ), pval(2) )
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6,a,g14.6)' ) '  X3 = ', real ( r(1) ), &
      ' with P(X3) = ', pval(1)
    write ( *, '(a,g14.6,a,g14.6)' ) '  X4 = ', real ( r(2) ), &
      ' with P(X4) = ', pval(2)
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat calculations in real ( kind = 8 ):'
  write ( *, '(a)' ) ' '

  call dpoly2_roots ( dp, dr )

  call dpoly_val_horner ( n, dp, real ( dr(1), kind = 8 ), dpval(1) )
  call dpoly_val_horner ( n, dp, real ( dr(2), kind = 8 ), dpval(2) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g20.12,a,g14.6)' ) '  X1 = ', real ( dr(1), kind = 8 ), &
    ' with P(X1) = ', dpval(1)
  write ( *, '(a,g20.12,a,g14.6)' ) '  X2 = ', real ( dr(2), kind = 8 ), &
    ' with P(X2) = ', dpval(2)

  call dpoly2_roots2 ( dp, dr, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  DPOLY2_ROOTS2 returned an error flag.'
  else
    call dpoly_val_horner ( n, dp, real ( dr(1), kind = 8 ), dpval(1) )
    call dpoly_val_horner ( n, dp, real ( dr(2), kind = 8 ), dpval(2) )
    write ( *, '(a)' ) ' '
    write ( *, '(a,g20.12,a,g14.6)' ) '  X3 = ', real ( dr(1), kind = 8 ), &
      ' with P(X3) = ', dpval(1)
    write ( *, '(a,g20.12,a,g14.6)' ) '  X4 = ', real ( dr(2), kind = 8 ), &
      ' with P(X4) = ', dpval(2)
  end if

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2

  integer ( kind = 4 ) ierror
  real ( kind = 4 ), dimension(0:n) :: p = (/ &
    3.9999999E+00, -4.0E+00, 1.0E+00 /)
  real ( kind = 4 ) pval(2)
  complex ( kind = 4 ) r(2)
  real ( kind = 4 ) temp
  real ( kind = 4 ) :: x1val = 1.999683772E+00
  real ( kind = 4 ) :: x2val = 2.000316228E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 21                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  P(X) =                                     |'
  write ( *, '(a)' ) '|       1         * x**2                      |'
  write ( *, '(a)' ) '|    -  4         * x                         |'
  write ( *, '(a)' ) '|    +  3.9999999                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Evaluate the standard quadratic formula:   |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    X1 = 0.5 * ( -B + sqrt(B**2-4*A*C) ) / A |'
  write ( *, '(a)' ) '|    X2 = 0.5 * ( -B - sqrt(B**2-4*A*C) ) / A |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Evaluate an alternate quadratic formula:   |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    X3 = 2 * C / ( - B + sqrt(B**2-4*A*C) )  |'
  write ( *, '(a)' ) '|    X4 = 2 * C / ( - B - sqrt(B**2-4*A*C) )  |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct results:                           |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    1.999683772                              |'
  write ( *, '(a)' ) '|    2.000316228                              |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'

  call rpoly2_roots ( p, r )

  call rpoly_val_horner ( n, p, real ( r(1) ), pval(1) )
  call rpoly_val_horner ( n, p, real ( r(2) ), pval(2) )

  write ( *, '(a,g18.10,a,g14.6)' ) '  X1 = ', real ( r(1) ), &
    ' with P(X1) = ', pval(1)
  write ( *, '(a,g18.10,a,g14.6)' ) '  X2 = ', real ( r(2) ), &
    ' with P(X2) = ', pval(2)

  call rpoly2_roots2 ( p, r, ierror )

  if ( ierror /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  RPOLY2_ROOTS2 returned an error flag.'

  else

    call rpoly_val_horner ( n, p, real ( r(1) ), pval(1) )
    call rpoly_val_horner ( n, p, real ( r(2) ), pval(2) )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g18.10,a,g14.6)' ) '  X3 = ', real ( r(1) ), &
      ' with P(X3) = ', pval(1)
    write ( *, '(a,g18.10,a,g14.6)' ) '  X4 = ', real ( r(2) ), &
      ' with P(X4) = ', pval(2)

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,g18.10)' ) '  True value = ', x1val
  write ( *, '(a,g18.10)' ) '  True value = ', x2val

  return
end
subroutine test215 ( )

!*****************************************************************************80
!
!! TEST215
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2
  integer ( kind = 4 ), parameter :: ntest = 3

  real ( kind = 4 ) aval(ntest)
  real ( kind = 4 ) bval(ntest)
  real ( kind = 4 ) cval(ntest)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  real ( kind = 4 ) p(0:n)
  real ( kind = 4 ) pval(2)
  complex ( kind = 4 ) r(2)
  real ( kind = 4 ) temp
  real ( kind = 4 ) x1val(ntest)
  real ( kind = 4 ) x2val(ntest)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 215                               |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Evaluate the standard quadratic formula:   |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    X1 = 0.5 * ( -B + sqrt(B**2-4*A*C) ) / A |'
  write ( *, '(a)' ) '|    X2 = 0.5 * ( -B - sqrt(B**2-4*A*C) ) / A |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Evaluate an alternate quadratic formula:   |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    X3 = 2 * C / ( - B + sqrt(B**2-4*A*C) )  |'
  write ( *, '(a)' ) '|    X4 = 2 * C / ( - B - sqrt(B**2-4*A*C) )  |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct results:                           |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    yabba                                    |'
  write ( *, '(a)' ) '|    dabba                                    |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'

  aval(1) = 6.0E+00
  bval(1) = 5.0E+00
  cval(1) = -4.0E+00
  x1val(1) = 0.5E+00
  x2val(1) = -4.0E+00 / 3.0E+00
!
!  Causes overflow.  We give up.
!
  temp = 1.0E+30

  aval(2) = 6.0E+00 * temp
  bval(2) = 5.0E+00 * temp
  cval(2) = -4.0E+00 * temp
  x1val(2) = 0.5E+00
  x2val(2) = -4.0E+00 / 3.0E+00
!
!  Roots were not given for this example.
!  Forsythe requests that only the root near 1 be solved for.
!
  aval(3) = 1.0E-30
  bval(3) = -1.0E+30
  cval(3) = 1.0E+30
  x1val(3) = 1.0E+00
!
!  Evaluation of the constant 1.0E+60 causes a compile time error
!  on the Alpha!  So we'll fudge and just set it to the maximal value
!  for now.
!
! x2val(3) = 1.0E+60
  x2val(3) = huge ( x2val(3) )

  do i = 1, ntest

    p(2) = aval(i)
    p(1) = bval(i)
    p(0) = cval(i)

    write ( *, '(a)' ) ' '
    write ( *, '(a,3g14.6)' ) '  A, B, C  = ', aval(i), bval(i), cval(i)

    if ( i == 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SKIP THIS TEST!.'
      write ( *, '(a)' ) '  It''s guaranteed to cause overflow for us!'
      cycle
    else if ( i == 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SKIP THIS TEST!.'
      write ( *, '(a)' ) '  It''s guaranteed to cause overflow for us!'
      cycle
    end if

    call rpoly2_roots ( p, r )

    call rpoly_val_horner ( n, p, real ( r(1) ), pval(1) )
    call rpoly_val_horner ( n, p, real ( r(2) ), pval(2) )

    write ( *, '(a,g18.10,a,g14.6)' ) '  X1 = ', real ( r(1) ), &
      ' with P(X1) = ', pval(1)
    write ( *, '(a,g18.10,a,g14.6)' ) '  X2 = ', real ( r(2) ), &
      ' with P(X2) = ', pval(2)

    call rpoly2_roots2 ( p, r, ierror )

    if ( ierror /= 0 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  RPOLY2_ROOTS2 returned an error flag.'

    else

      call rpoly_val_horner ( n, p, real ( r(1) ), pval(1) )
      call rpoly_val_horner ( n, p, real ( r(2) ), pval(2) )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g18.10,a,g14.6)' ) '  X3 = ', real ( r(1) ), &
        ' with P(X3) = ', pval(1)
      write ( *, '(a,g18.10,a,g14.6)' ) '  X4 = ', real ( r(2) ), &
        ' with P(X4) = ', pval(2)

    end if

    write ( *, '(a,g18.10)' ) '  True value = ', x1val(i)
    write ( *, '(a,g18.10)' ) '  True value = ', x2val(i)

  end do

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) facn
  real ( kind = 4 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 4 ) fx2old
  real ( kind = 4 ) fx3
  integer ( kind = 4 ) n
  real ( kind = 4 ) sum2
  real ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 22                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Define F(X,N) to be the difference between |'
  write ( *, '(a)' ) '|  EXP(X) and its Taylor series, multiplied   |'
  write ( *, '(a)' ) '|  by N!.                                     |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Formula 1:                                 |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    F(X,N) = N! *                            |'
  write ( *, '(a)' ) '|      ( EXP(X) - (1+X+X**2/2+...+X**N/N!) )  |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  But we can also determine a second         |'
  write ( *, '(a)' ) '|  definition, namely:                        |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Formula 2:                                 |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    F(X,N) = N*F(X,N-1) - X**N               |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Compare formula 1 with formula 2 at X = 1  |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    N  Formula 1  Formula 2  x**n/(n+1)'
  write ( *, '(a)' ) ' '

  x = 1.0E+00
  facn = 1.0E+00
  sum2 = 1.0E+00
  fx2old = exp ( x ) - 1.0E+00

  do n = 0, 20

    fx1 = facn * ( exp ( x ) - sum2 )

    if ( n == 0 ) then
      fx2 = fx1
      fx3 = 0.0E+00
    else
      fx2 = n * fx2old - x**n
      fx3 =  x**n / real ( n )
    end if

    fx2old = fx2
    facn = facn * ( n + 1 )
    sum2 = sum2 + x**( n + 1 ) / facn

    write ( *, '(i4,3g14.6)' ) n, fx1, fx2, fx3

  end do

  return
end
subroutine test23 ( )

!*****************************************************************************80
!
!! TEST23
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) fac
  integer ( kind = 4 ) i
  real ( kind = 4 ) sum2
  real ( kind = 4 ) sumold
  real ( kind = 4 ) temp
  real ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 23                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Estimate EXP(-5.5) in two ways:            |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    R1(N) = sum ( I = 0 to N ) X**I/I        |'
  write ( *, '(a)' ) '|    R2(N) = 1 / sum ( I = 0 to N ) (-X)**I/I |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct value:                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    EXP(-5.5) = 0.00408677143846406699       |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '

  sum2 = 0.0E+00
  i = 0
  fac = 1.0E+00
  x = - 5.5E+00

  do

    temp = x**i / fac
    sumold = sum2
    sum2 = sum2 + temp
    i = i + 1
    fac = fac * i
    temp = sum2 - sumold

    if ( temp == 0.0E+00 ) then
      exit
    end if

  end do

  write ( *, '(a,g14.6)' ) '  Computed R1 = ', sum2

  sum2 = 0.0E+00
  i = 0
  fac = 1.0E+00
  x = 5.5E+00

  do

    temp = x**i / fac
    sumold = sum2
    sum2 = sum2 + temp
    i = i + 1
    fac = fac * i
    temp = sum2 - sumold

    if ( temp == 0.0E+00 ) then
      exit
    end if

  end do

  sum2 = 1.0E+00 / sum2

  write ( *, '(a,g14.6)' ) '  Computed R2 = ', sum2

  return
end
subroutine test24 ( )

!*****************************************************************************80
!
!! TEST24
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2

  real ( kind = 4 ) a(n,n)
  real ( kind = 4 ) a_exp(n,n)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 24                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Estimate the matrix exponential e^A        |'
  write ( *, '(a)' ) '|  using a naive Taylor series approach.      |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  The matrix A is                            |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    -49 24                                   |'
  write ( *, '(a)' ) '|    -64 31                                   |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct value of e^A:                      |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|     -0.735759  0.551819                     |'
  write ( *, '(a)' ) '|     -1.471518  1.103638                     |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Reported computed value of e^A:            |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    -22.25880  -1.432766                     |'
  write ( *, '(a)' ) '|    -61.49931  -3.474280                     |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'

  a(1:2,1:2) = reshape ( (/ -49, -64, 24, 31 /), (/ 2, 2 /) )

  call matrix_exponential_taylor ( n, a, a_exp )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Computed value of e^A:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2f14.6)' ) a_exp(i,1:n)
  end do

  return
end
subroutine test25 ( )

!*****************************************************************************80
!
!! TEST25
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Fred Gruenberger,
!    Computer Recreations,
!    Scientific American,
!    April, 1984, pages 19-26.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: n = 27
  real ( kind = 4 ), parameter :: x = 1.0000001E+00
  real ( kind = 4 ) x_1
  real ( kind = 4 ) x_2
  real ( kind = 4 ) x_3
  real ( kind = 4 ) x_4
  real ( kind = 4 ), parameter :: x_n = 674530.4707E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 25                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Estimate an integer power of a             |'
  write ( *, '(a)' ) '|  real number.                               |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  The computation is                         |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    (1.0000001)^(2^27)                       |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct value is:                          |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    674,530.4707                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '

  x_1 = x**(2**n)

  x_2 = x
  do i = 1, n
    x_2 = x_2 * x_2
  end do

  x_3 = exp ( log ( x ) * real ( 2**n ) )

  x_4 = exp ( log ( x ) * exp ( log ( 2.0E+00 ) * real ( n ) ) )

  write ( *, '(a,g16.8)' ) '  X**(2**N) =            ', x_1
  write ( *, '(a,g16.8)' ) '  X*X*...*X =            ', x_2
  write ( *, '(a,g16.8)' ) '  E^(logX*27**N)) =      ', x_3
  write ( *, '(a,g16.8)' ) '  E^(logX*E^(27*logN)) = ', x_4

  return
end
subroutine test26 ( )

!*****************************************************************************80
!
!! TEST26
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: n = 27
  real ( kind = 8 ), parameter :: x = 1.0000001D+00
  real ( kind = 8 ) x_1
  real ( kind = 8 ) x_2
  real ( kind = 8 ) x_3
  real ( kind = 8 ) x_4
  real ( kind = 8 ), parameter :: x_n = 674530.4707D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Exercise 26                                |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Estimate an integer power of a             |'
  write ( *, '(a)' ) '|  real number.                               |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  In this version, use DOUBLE PRECISION.     |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  The computation is                         |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    (1.0000001)^(2^27)                       |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|  Correct value is:                          |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '|    674,530.4707                             |'
  write ( *, '(a)' ) '|                                             |'
  write ( *, '(a)' ) '+=============================================+'
  write ( *, '(a)' ) ' '

  x_1 = x**(2**n)

  x_2 = x
  do i = 1, n
    x_2 = x_2 * x_2
  end do

  x_3 = exp ( log ( x ) * real ( 2**n, kind = 8 ) )

  x_4 = exp ( log ( x ) * exp ( log ( 2.0D+00 ) * real ( n, kind = 8 ) ) )

  write ( *, '(a,g16.8)' ) '  X**(2**N) =            ', x_1
  write ( *, '(a,g16.8)' ) '  X*X*...*X =            ', x_2
  write ( *, '(a,g16.8)' ) '  E^(logX*27**N)) =      ', x_3
  write ( *, '(a,g16.8)' ) '  E^(logX*E^(27*logN)) = ', x_4

  return
end
function p8horn ( x )

!*****************************************************************************80
!
!! P8HORN evaluates a certain polynomial using Horner's method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the polynomial.
!
!    Output, real P8HORN, the value of the polynomial.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 4 ), dimension ( 0:n ) :: c = (/ 18.601E+00, 168.97E+00, &
    -356.41E+00, 170.4E+00 /)
  real ( kind = 4 ) p8horn
  real ( kind = 4 ) x

  call rpoly_val_horner ( n, c, x, p8horn )

  return
end
function p8reg ( x )

!*****************************************************************************80
!
!! P8REG evaluates a certain polynomial using the simple method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 4 ) c(0:n)
  real ( kind = 4 ) p8reg
  real ( kind = 4 ) x

  c(0) =    18.601E+00
  c(1) =   168.97E+00
  c(2) = - 356.41E+00
  c(3) =   170.4E+00

  call rpoly_val ( n, c, x, p8reg )

  return
end
function fhorn11 ( x )

!*****************************************************************************80
!
!! FHORN11 evaluates a certain polynomial using Horner's method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 4 ) fhorn11
  real ( kind = 4 ) p(0:n)
  real ( kind = 4 ) x

  p(0) = 1.0E-27
  p(1) = 408855776.0E+00
  p(2) = -708158977.0E+00
  p(3) = -1226567328.0E+00
  p(4) = 2124476931.0E+00

  call rpoly_val_horner ( n, p, x, fhorn11 )

  return
end
function freg11 ( x )

!*****************************************************************************80
!
!! FREG11 evaluates a certain polynomial using the simple method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 4 ) freg11
  real ( kind = 4 ) p(0:n)
  real ( kind = 4 ) x

  p(0) = 1.0E-27
  p(1) = 408855776.0E+00
  p(2) = -708158977.0E+00
  p(3) = -1226567328.0E+00
  p(4) = 2124476931.0E+00

  call rpoly_val ( n, p, x, freg11 )

  return
end
function fx15 ( x )

!*****************************************************************************80
!
!! FX15 evaluates a certain function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) fx15
  real ( kind = 4 ) x

  fx15 = ( 4970.0E+00 * x - 4923.0E+00 ) / &
    ( 4970.0E+00 * x**2 - 9799.0E+00 * x + 4830.0E+00 )

  return
end

