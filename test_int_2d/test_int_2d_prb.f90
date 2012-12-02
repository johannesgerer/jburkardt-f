program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_INT_2D_PRB.
!
!  Discussion:
!
!    TEST_INT_2D_PRB demonstrates the TEST_INT_2D integration test functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INT_2D_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_INT_2D library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INT_2D_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 applies a Monte Carlo rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a(2)
  real ( kind = 8 ) b(2)
  integer ( kind = 4 ) dim
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable :: fx(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num
  real ( kind = 8 ) quad
  integer ( kind = 4 ) seed
  real ( kind = 8 ) volume
  real ( kind = 8 ), allocatable :: x(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use a Monte Carlo rule.'
  write ( *, '(a)' ) '  Repeatedly multiply the number of points by 4.'

  call p00_problem_num ( problem_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '   Problem      Points       Approx         Error'

  do problem = 1, problem_num

    write ( *, '(a)' ) ' '

    n = 1

    do i = 1, 12

      seed = 123456789

      allocate ( x(2,n) )
      allocate ( fx(n) )

      call r8mat_uniform_01 ( 2, n, seed, x )

      call p00_lim ( problem, a, b )

      do dim = 1, 2
        x(dim,1:n) = ( 1.0D+00 - x(dim,1:n) ) * a(dim) &
                   +             x(dim,1:n)   * b(dim)
      end do

      volume = product ( b(1:2) - a(1:2) )

      call p00_fun ( problem, n, x, fx )

      quad = volume * sum ( fx(1:n) ) / real ( n, kind = 8 )

      call p00_exact ( problem, exact )

      error = abs ( quad - exact )

      write ( *, '(2x,i8,2x,i10,2x,g14.6,2x,g14.6)' ) &
        problem, n, quad, error

      deallocate ( fx )
      deallocate ( x )

      n = n * 4

    end do

    write ( *, '(2x,i8,2x,a10,2x,g14.6)' ) &
      problem, '     Exact', exact

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 applies a product of composite midpoint rules.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a(2)
  real ( kind = 8 ) b(2)
  integer ( kind = 4 ) dim
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable :: fx(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num
  real ( kind = 8 ) quad
  real ( kind = 8 ) volume
  real ( kind = 8 ), allocatable :: x(:,:)
  real ( kind = 8 ) xval
  real ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Use a product of composite midpoint rules.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeatedly multiply the number of points by 4.'

  call p00_problem_num ( problem_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '   Problem      Points       Approx         Error'

  do problem = 1, problem_num

    write ( *, '(a)' ) ' '

    nx = 1
    ny = 1

    do i = 1, 12

      n = nx * ny

      allocate ( x(2,n) )
      allocate ( fx(n) )

      call p00_lim ( problem, a, b )

      k = 0

      do ix = 1, nx

        xval = ( real ( 2 * nx - 2 * ix + 1, kind = 8 ) * a(1)   &
               + real (          2 * ix - 1, kind = 8 ) * b(1) ) &
               / real ( 2 * nx,              kind = 8 )

        do iy = 1, ny

          yval = ( real ( 2 * ny - 2 * iy + 1, kind = 8 ) * a(2)   &
                 + real (          2 * iy - 1, kind = 8 ) * b(2) ) &
                 / real ( 2 * ny,              kind = 8 )

          k = k + 1
          x(1,k) = xval
          x(2,k) = yval

        end do

      end do

      volume = product ( b(1:2) - a(1:2) )

      call p00_fun ( problem, n, x, fx )

      quad = volume * sum ( fx(1:n) ) / real ( n, kind = 8 )

      call p00_exact ( problem, exact )

      error = abs ( quad - exact )

      write ( *, '(2x,i8,2x,i10,2x,g14.6,2x,g14.6)' ) &
        problem, n, quad, error

      deallocate ( fx )
      deallocate ( x )

      nx = nx * 2
      ny = ny * 2

    end do

    write ( *, '(2x,i8,2x,a10,2x,g14.6)' ) &
      problem, '     Exact', exact

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 applies a product of Gauss-Legendre rules.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a(2)
  real ( kind = 8 ) b(2)
  integer ( kind = 4 ) dim
  real ( kind = 8 ) error
  real ( kind = 8 ) exact
  real ( kind = 8 ), allocatable :: fxy(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) nxy
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) problem
  integer ( kind = 4 ) problem_num
  real ( kind = 8 ) quad
  real ( kind = 8 ) volume
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: wxy(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: xy(:,:)
  real ( kind = 8 ) xval
  real ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Use a product of Gauss-Legendre rules.'
  write ( *, '(a)' ) '  The 1D rules essentially double in order.'

  call p00_problem_num ( problem_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '   Problem      Points       Approx         Error'

  do problem = 1, problem_num

    write ( *, '(a)' ) ' '

    nx = 1
    ny = 1

    do i = 1, 8

      allocate ( x(1:nx) )
      allocate ( w(1:nx) )

      call legendre_set ( nx, x, w )

      nxy = nx * ny

      allocate ( wxy(nxy) )
      allocate ( xy(2,nxy) )
      allocate ( fxy(nxy) )

      call p00_lim ( problem, a, b )

      k = 0

      do ix = 1, nx

        xval = ( ( 1.0D+00 + x(ix) ) * a(1)   &
               + ( 1.0D+00 - x(ix) ) * b(1) ) &
               /   2.0D+00

        do iy = 1, ny

          yval = ( ( 1.0D+00 + x(iy) ) * a(2)   &
                 + ( 1.0D+00 - x(iy) ) * b(2) ) &
                 /   2.0D+00

          k = k + 1
          xy(1,k) = xval
          xy(2,k) = yval
          wxy(k) = w(ix) * w(iy)

        end do

      end do

      volume = product ( b(1:2) - a(1:2) )

      call p00_fun ( problem, nxy, xy, fxy )

      quad = volume * dot_product ( wxy(1:nxy), fxy(1:nxy) ) / 4.0D+00

      call p00_exact ( problem, exact )

      error = abs ( quad - exact )

      write ( *, '(2x,i8,2x,i10,2x,g14.6,2x,g14.6)' ) &
        problem, nxy, quad, error

      deallocate ( fxy )
      deallocate ( w )
      deallocate ( wxy )
      deallocate ( x )
      deallocate ( xy )

      nx = 2 * nx + 1
      ny = nx

    end do

    write ( *, '(2x,i8,2x,a10,2x,g14.6)' ) &
      problem, '     Exact', exact

  end do

  return
end
