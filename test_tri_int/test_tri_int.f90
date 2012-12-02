subroutine get_prob_num ( prob_num )

!*****************************************************************************80
!
!! GET_PROB_NUM returns the number of test integration problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) PROB_NUM, the number of problems.
!
  implicit none

  integer ( kind = 4 ) prob_num

  prob_num = 22

  return
end
subroutine p00_fun ( problem, n, p, f )

!*****************************************************************************80
!
!! P00_FUN evaluates the integrand for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)
  integer ( kind = 4 ) problem

  if ( problem == 1 ) then
    call p01_fun ( n, p, f )
  else if ( problem == 2 ) then
    call p02_fun ( n, p, f )
  else if ( problem == 3 ) then
    call p03_fun ( n, p, f )
  else if ( problem == 4 ) then
    call p04_fun ( n, p, f )
  else if ( problem == 5 ) then
    call p05_fun ( n, p, f )
  else if ( problem == 6 ) then
    call p06_fun ( n, p, f )
  else if ( problem == 7 ) then
    call p07_fun ( n, p, f )
  else if ( problem == 8 ) then
    call p08_fun ( n, p, f )
  else if ( problem == 9 ) then
    call p09_fun ( n, p, f )
  else if ( problem == 10 ) then
    call p10_fun ( n, p, f )
  else if ( problem == 11 ) then
    call p11_fun ( n, p, f )
  else if ( problem == 12 ) then
    call p12_fun ( n, p, f )
  else if ( problem == 13 ) then
    call p13_fun ( n, p, f )
  else if ( problem == 14 ) then
    call p14_fun ( n, p, f )
  else if ( problem == 15 ) then
    call p15_fun ( n, p, f )
  else if ( problem == 16 ) then
    call p16_fun ( n, p, f )
  else if ( problem == 17 ) then
    call p17_fun ( n, p, f )
  else if ( problem == 18 ) then
    call p18_fun ( n, p, f )
  else if ( problem == 19 ) then
    call p19_fun ( n, p, f )
  else if ( problem == 20 ) then
    call p20_fun ( n, p, f )
  else if ( problem == 21 ) then
    call p21_fun ( n, p, f )
  else if ( problem == 22 ) then
    call p22_fun ( n, p, f )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_FUN - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_monte_carlo ( problem, n, seed, result )

!*****************************************************************************80
!
!! P00_MONTE_CARLO applies the Monte Carlo rule to integrate a function.
!
!  Discussion:
!
!    The function f(x,y) is to be integrated over a triangle T.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) N, the number of sample points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) area
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)
  integer ( kind = 4 ) problem
  real ( kind = 8 ) result
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t(2,3)

  call p00_vertices ( problem, t )

  call triangle_sample ( t, n, seed, p )

  call p00_fun ( problem, n, p, f )

  call triangle_area ( t, area )

  result = area * sum ( f(1:n) ) / real ( n, kind = 8 )

  return
end
subroutine p00_singularity ( problem, singularity )

!*****************************************************************************80
!
!! P00_SINGULARITY warns of common singularities for any problem.
!
!  Discussion:
!
!    This routine can be used to check whether the integrand function
!    for a given problem has singularities at the vertices or along
!    the edges of the triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Output, integer ( kind = 4 ) SINGULARITY.
!    0, there are no vertex or edge singularities.
!    1, there are singularities at one or more vertices, but not on edges.
!    2, there are singularities on one or more edges, possibly
!       including vertices.
!    3, there are singularities somewhere inside or on the triangle.
!
  implicit none

  integer ( kind = 4 ) problem
  integer ( kind = 4 ) singularity

  if ( problem == 1 ) then
    singularity = 0
  else if ( problem == 2 ) then
    singularity = 0
  else if ( problem == 3 ) then
    singularity = 0
  else if ( problem == 4 ) then
    singularity = 0
  else if ( problem == 5 ) then
    singularity = 0
  else if ( problem == 6 ) then
    singularity = 0
  else if ( problem == 7 ) then
    singularity = 0
  else if ( problem == 8 ) then
    singularity = 0
  else if ( problem == 9 ) then
    singularity = 0
  else if ( problem == 10 ) then
    singularity = 2
  else if ( problem == 11 ) then
    singularity = 1
  else if ( problem == 12 ) then
    singularity = 2
  else if ( problem == 13 ) then
    singularity = 2
  else if ( problem == 14 ) then
    singularity = 2
  else if ( problem == 15 ) then
    singularity = 2
  else if ( problem == 16 ) then
    singularity = 2
  else if ( problem == 17 ) then
    singularity = 3
  else if ( problem == 18 ) then
    singularity = 1
  else if ( problem == 19 ) then
    singularity = 0
  else if ( problem == 20 ) then
    singularity = 0
  else if ( problem == 21 ) then
    singularity = 1
  else if ( problem == 22 ) then
    singularity = 1
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_SINGULARITY - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_title ( problem, title )

!*****************************************************************************80
!
!! P00_TITLE returns the title for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Output, character ( len = * ) TITLE, the title.
!
  implicit none

  integer ( kind = 4 ) problem
  character ( len = * ) title

  if ( problem == 1 ) then
    call p01_title ( title )
  else if ( problem == 2 ) then
    call p02_title ( title )
  else if ( problem == 3 ) then
    call p03_title ( title )
  else if ( problem == 4 ) then
    call p04_title ( title )
  else if ( problem == 5 ) then
    call p05_title ( title )
  else if ( problem == 6 ) then
    call p06_title ( title )
  else if ( problem == 7 ) then
    call p07_title ( title )
  else if ( problem == 8 ) then
    call p08_title ( title )
  else if ( problem == 9 ) then
    call p09_title ( title )
  else if ( problem == 10 ) then
    call p10_title ( title )
  else if ( problem == 11 ) then
    call p11_title ( title )
  else if ( problem == 12 ) then
    call p12_title ( title )
  else if ( problem == 13 ) then
    call p13_title ( title )
  else if ( problem == 14 ) then
    call p14_title ( title )
  else if ( problem == 15 ) then
    call p15_title ( title )
  else if ( problem == 16 ) then
    call p16_title ( title )
  else if ( problem == 17 ) then
    call p17_title ( title )
  else if ( problem == 18 ) then
    call p18_title ( title )
  else if ( problem == 19 ) then
    call p19_title ( title )
  else if ( problem == 20 ) then
    call p20_title ( title )
  else if ( problem == 21 ) then
    call p21_title ( title )
  else if ( problem == 22 ) then
    call p22_title ( title )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_vertex_sub ( problem, level, n, result )

!*****************************************************************************80
!
!! P00_VERTEX_SUB approximates an integral in a triangle by subdivision.
!
!  Discussion:
!
!    The function f(x,y) is to be integrated over a triangle T.
!
!    The first approximation averages the values at the vertices.
!
!    If a second approximation is requested, the routine subdivides each
!    existing triangle into 4, evaluates the function at the new vertices,
!    and returns an improved estimate.
!
!    The routine may be called repeatedly in this way, to get an improved
!    estimate of the integral.
!
!    Note that this routine will fail in the case that there
!    are singularities at the vertices or along the sides of the triangle.
!
!    Moreover, since the number of new vertices grows as a power of 4,
!    the use of an automatic array to store all the new vertices at one
!    time may fail when a memory limit is reached.
!
!    Finally, note that the rule has a very low order of convergence.
!    We're just exhibiting this routine as an EXAMPLE of how one
!    would use the test integrands.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) LEVEL, the level of subdivision.  The first
!    call should be with LEVEL = 0.  For successive refinement, the routine
!    may be called repeatedly.  Each time, the user should increase the
!    value of LEVEL by 1, and also input the value of RESULT that was
!    output on the previous call.
!
!    Input/output, integer ( kind = 4 ) N, the number of function evaluations
!    used.  If LEVEL = 0, the input value is ignored.  Otherwise, the input
!    value is assumed to be the output value from the previous call.
!
!    Input/output, real ( kind = 8 ) RESULT, the approximate integral.
!    If LEVEL = 0, then the input value is ignored.  Otherwise, the
!    input value is assumed to be the result from the previous call,
!    at the previous level.  The output value is based on the input
!    value, adjusted by information determined at the new level.
!
  implicit none

  real ( kind = 8 ), save :: area
  real ( kind = 8 ), allocatable, dimension ( : ) :: f
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) level
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_new
  integer ( kind = 4 ) order_max_1d
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: p
  integer ( kind = 4 ) problem
  real ( kind = 8 ) result
  real ( kind = 8 ) result_new
  real ( kind = 8 ), save, dimension ( 2, 3 ) :: t
  real ( kind = 8 ) x
  real ( kind = 8 ) xsi(3)
  real ( kind = 8 ) y
!
!  Compute first level.
!
  if ( level == 0 ) then

    call p00_vertices ( problem, t )
    call triangle_area ( t, area )

    n = 3
    allocate ( f(1:n) )
    allocate ( p(1:2,1:n) )

    call p00_fun ( problem, n, t, f )

    result = sum ( f(1:n) ) * area / 3.0D+00
!
!  Compute next level.
!
  else

    if ( level == 1 ) then
      n_new = 3
    else
      n_new = ( 2**(level-1) + 1 ) * 2**(level-1) * 3
    end if

    allocate ( f(1:n_new) )
    allocate ( p(1:2,1:n_new) )

    order_max_1d = 2**level

    n_new = 0

    do i = 0, order_max_1d - 1, 2
      do j = 0, order_max_1d - 1 - i, 2

        n_new = n_new + 1

        xsi(1) = real ( order_max_1d - i - 1 - j, kind = 8 ) &
               / real ( order_max_1d,             kind = 8 )
        xsi(2) = real (                i + 1,     kind = 8 ) &
               / real ( order_max_1d,             kind = 8 )
        xsi(3) = real (                        j, kind = 8 ) &
               / real ( order_max_1d,             kind = 8 )

        p(1:2,n_new) = matmul ( t(1:2,1:3), xsi(1:3) )

        n_new = n_new + 1

        xsi(1) = real ( order_max_1d - i - j - 1, kind = 8 ) &
               / real ( order_max_1d,             kind = 8 )
        xsi(2) = real (                i,         kind = 8 ) &
               / real ( order_max_1d,             kind = 8 )
        xsi(3) = real (                    j + 1, kind = 8 ) &
               / real ( order_max_1d,             kind = 8 )

        p(1:2,n_new) = matmul ( t(1:2,1:3), xsi(1:3) )

        n_new = n_new + 1

        xsi(1) = real ( order_max_1d - i - 1 - j - 1, kind = 8 ) &
               / real ( order_max_1d,                 kind = 8 )
        xsi(2) = real (                i + 1,         kind = 8 ) &
               / real ( order_max_1d,                 kind = 8 )
        xsi(3) = real (                        j + 1, kind = 8 ) &
               / real ( order_max_1d,                 kind = 8 )

        p(1:2,n_new) = matmul ( t(1:2,1:3), xsi(1:3) )

      end do
    end do

    call p00_fun ( problem, n_new, p, f )

    result_new = sum ( f(1:n_new) ) * area / real ( n_new, kind = 8 )

    result = ( 3.0D+00 * result_new + result ) / 4.0D+00

    n = n + n_new

  end if

  deallocate ( f )
  deallocate ( p )

  return
end
subroutine p00_vertices ( problem, t )

!*****************************************************************************80
!
!! P00_VERTICES returns the vertices for any problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem number.
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  integer ( kind = 4 ) problem
  real ( kind = 8 ) t(2,3)

  if ( problem == 1 ) then
    call p01_vertices ( t )
  else if ( problem == 2 ) then
    call p02_vertices ( t )
  else if ( problem == 3 ) then
    call p03_vertices ( t )
  else if ( problem == 4 ) then
    call p04_vertices ( t )
  else if ( problem == 5 ) then
    call p05_vertices ( t )
  else if ( problem == 6 ) then
    call p06_vertices ( t )
  else if ( problem == 7 ) then
    call p07_vertices ( t )
  else if ( problem == 8 ) then
    call p08_vertices ( t )
  else if ( problem == 9 ) then
    call p09_vertices ( t )
  else if ( problem == 10 ) then
    call p10_vertices ( t )
  else if ( problem == 11 ) then
    call p11_vertices ( t )
  else if ( problem == 12 ) then
    call p12_vertices ( t )
  else if ( problem == 13 ) then
    call p13_vertices ( t )
  else if ( problem == 14 ) then
    call p14_vertices ( t )
  else if ( problem == 15 ) then
    call p15_vertices ( t )
  else if ( problem == 16 ) then
    call p16_vertices ( t )
  else if ( problem == 17 ) then
    call p17_vertices ( t )
  else if ( problem == 18 ) then
    call p18_vertices ( t )
  else if ( problem == 19 ) then
    call p19_vertices ( t )
  else if ( problem == 20 ) then
    call p20_vertices ( t )
  else if ( problem == 21 ) then
    call p21_vertices ( t )
  else if ( problem == 22 ) then
    call p22_vertices ( t )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'P00_VERTICES - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal problem number = ', problem
    stop
  end if

  return
end
subroutine p00_wandzura05_sub ( problem, level, n, result )

!*****************************************************************************80
!
!! P00_WANDZURA05_SUB uses subdivision and a Wandzura rule.
!
!  Discussion:
!
!    The Wandzura rule is a seven point rule of polynomial exactness 5.
!
!    The function f(x,y) is to be integrated over a triangle T.
!
!    The triangle is subdivided by subdividing each side into LEVEL sections,
!    which produces LEVEL*LEVEL subtriangles.  The Wandzura rule is then
!    applied to each subtriangle, and the result is summed.
!
!    The abscissas of this Wandzura rule do not lie on the vertices
!    or sides of the reference triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Stephen Wandzura, Hong Xiao,
!    Symmetric Quadrature Rules on a Triangle,
!    Computers and Mathematics with Applications,
!    Volume 45, pages 1829-1840, 2003.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PROBLEM, the problem index.
!
!    Input, integer ( kind = 4 ) LEVEL, the level of subdivision.  This
!    indicates the number of equally spaced subedges into which each edge of
!    the triangle is to be divided.  This will result in a total of
!    LEVEL*LEVEL subtriangles being used.
!
!    Output, integer ( kind = 4 ) N, the number of function evaluations used.
!
!    Output, real ( kind = 8 ) RESULT, the approximate integral.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_3 = 3
  integer ( kind = 4 ), parameter :: order = 7

  real ( kind = 8 ) area
  real ( kind = 8 ) f(order)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) level
  logical more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) problem
  real ( kind = 8 ) result
  real ( kind = 8 ) sub_tri_phys(2,3)
  real ( kind = 8 ) tri_phys(2,3)
  real ( kind = 8 ), dimension ( order ) :: w = (/ &
    0.22500000000000D+00, &
    0.13239415278851D+00, &
    0.13239415278851D+00, &
    0.13239415278851D+00, &
    0.12593918054483D+00, &
    0.12593918054483D+00, &
    0.12593918054483D+00  &
    /)
  real ( kind = 8 ) xsi(3,3)
  real ( kind = 8 ) xy_phys(2,order)
  real ( kind = 8 ), dimension ( 3, order ) :: xy_ref = reshape ( (/ &
      0.33333333333333D+00, 0.33333333333333D+00, 0.33333333333333D+00, &
      0.05971587178977D+00, 0.47014206410512D+00, 0.47014206410512D+00, &
      0.47014206410512D+00, 0.05971587178977D+00, 0.47014206410512D+00, &
      0.47014206410512D+00, 0.47014206410512D+00, 0.05971587178977D+00, &
      0.79742698535309D+00, 0.10128650732346D+00, 0.10128650732346D+00, &
      0.10128650732346D+00, 0.79742698535309D+00, 0.10128650732346D+00, &
      0.10128650732346D+00, 0.10128650732346D+00, 0.79742698535309D+00  &
    /), (/ i4_3, order /) )

  result = 0.0D+00

  call p00_vertices ( problem, tri_phys )

  call triangle_area ( tri_phys, area )

  more = .false.

  do
!
!  Get the integer indices of the next reference subtriangle.
!
    call subtriangle_next ( level, more, i1, j1, i2, j2, i3, j3 )
!
!  Get the barycentric coordinates of the vertices of the reference subtriangle.
!
    xsi(1,1) = real (         i1,      kind = 8 ) / real ( level, kind = 8 )
    xsi(2,1) = real (              j1, kind = 8 ) / real ( level, kind = 8 )
    xsi(3,1) = real ( level - i1 - j1, kind = 8 ) / real ( level, kind = 8 )

    xsi(1,2) = real (         i2,      kind = 8 ) / real ( level, kind = 8 )
    xsi(2,2) = real (              j2, kind = 8 ) / real ( level, kind = 8 )
    xsi(3,2) = real ( level - i2 - j2, kind = 8 ) / real ( level, kind = 8 )

    xsi(1,3) = real (         i3,      kind = 8 ) / real ( level, kind = 8 )
    xsi(2,3) = real (              j3, kind = 8 ) / real ( level, kind = 8 )
    xsi(3,3) = real ( level - i3 - j3, kind = 8 ) / real ( level, kind = 8 )
!
!  Map the reference subtriangle to the physical subtriangle.
!
    sub_tri_phys(1:2,1:3) = matmul ( tri_phys(1:2,1:3), xsi(1:3,1:3) )
!
!  Now map the integration abscissas to the physical subtriangle.
!
    xy_phys(1:2,1:order) = matmul ( sub_tri_phys(1:2,1:3), xy_ref(1:3,1:order) )
!
!  Evaluate the function.
!
    call p00_fun ( problem, order, xy_phys, f )
!
!  Update the quadrature estimate.
!
    result = result + dot_product ( w(1:order), f(1:order) )

    if ( .not. more ) then
      exit
    end if

  end do
!
!  Scale by area and number of subtriangles.
!
  n = level * level * order

  result = result * area / real ( level * level, kind = 8 )

  return
end
subroutine p01_fun ( n, p, f )

!*****************************************************************************80
!
!! P01_FUN evaluates the integrand for problem 1.
!
!  Integrand:
!
!    f(x,y) = 2
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = 2.0D+00

  return
end
subroutine p01_title ( title )

!*****************************************************************************80
!
!! P01_TITLE returns the title of problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = 2'

  return
end
subroutine p01_vertices ( t )

!*****************************************************************************80
!
!! P01_VERTICES returns the vertices for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p02_fun ( n, p, f )

!*****************************************************************************80
!
!! P02_FUN evaluates the integrand for problem 2.
!
!  Integrand:
!
!    f(x,y) = 6 * x
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = 6.0D+00 * p(1,1:n)

  return
end
subroutine p02_title ( title )

!*****************************************************************************80
!
!! P02_TITLE returns the title of problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = 6 * x'

  return
end
subroutine p02_vertices ( t )

!*****************************************************************************80
!
!! P02_VERTICES returns the vertices for problem 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p03_fun ( n, p, f )

!*****************************************************************************80
!
!! P03_FUN evaluates the integrand for problem 3.
!
!  Integrand:
!
!    f(x,y) = 6 * y
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = 6.0D+00 * p(2,1:n)

  return
end
subroutine p03_title ( title )

!*****************************************************************************80
!
!! P03_TITLE returns the title of problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = 6 * y'

  return
end
subroutine p03_vertices ( t )

!*****************************************************************************80
!
!! P03_VERTICES returns the vertices for problem 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p04_fun ( n, p, f )

!*****************************************************************************80
!
!! P04_FUN evaluates the integrand for problem 4.
!
!  Integrand:
!
!    f(x,y) = 12 * x^2
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = 12.0D+00 * p(1,1:n)**2

  return
end
subroutine p04_title ( title )

!*****************************************************************************80
!
!! P04_TITLE returns the title of problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = 12 * x^2'

  return
end
subroutine p04_vertices ( t )

!*****************************************************************************80
!
!! P04_VERTICES returns the vertices for problem 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p05_fun ( n, p, f )

!*****************************************************************************80
!
!! P05_FUN evaluates the integrand for problem 5.
!
!  Integrand:
!
!    f(x,y) = 24 * x*y
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = 24.0D+00 * p(1,1:n) * p(2,1:n)

  return
end
subroutine p05_title ( title )

!*****************************************************************************80
!
!! P05_TITLE returns the title of problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = 24 * x*y'

  return
end
subroutine p05_vertices ( t )

!*****************************************************************************80
!
!! P05_VERTICES returns the vertices for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p06_fun ( n, p, f )

!*****************************************************************************80
!
!! P06_FUN evaluates the integrand for problem 6.
!
!  Integrand:
!
!    f(x,y) = 12 * y^2
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = 12.0D+00 * p(2,1:n)**2

  return
end
subroutine p06_title ( title )

!*****************************************************************************80
!
!! P06_TITLE returns the title of problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = 12 * y^2'

  return
end
subroutine p06_vertices ( t )

!*****************************************************************************80
!
!! P06_VERTICES returns the vertices for problem 6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p07_fun ( n, p, f )

!*****************************************************************************80
!
!! P07_FUN evaluates the integrand for problem 7.
!
!  Integrand:
!
!    f(x,y) = 20 * x^3
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = 20.0D+00 * p(1,1:n)**3

  return
end
subroutine p07_title ( title )

!*****************************************************************************80
!
!! P07_TITLE returns the title of problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = 20 * x^3'

  return
end
subroutine p07_vertices ( t )

!*****************************************************************************80
!
!! P07_VERTICES returns the vertices for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p08_fun ( n, p, f )

!*****************************************************************************80
!
!! P08_FUN evaluates the integrand for problem 8.
!
!  Integrand:
!
!    f(x,y) = 30 * x^4
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = 30.0D+00 * p(1,1:n)**4

  return
end
subroutine p08_title ( title )

!*****************************************************************************80
!
!! P08_TITLE returns the title of problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = 30 * x^4'

  return
end
subroutine p08_vertices ( t )

!*****************************************************************************80
!
!! P08_VERTICES returns the vertices for problem 8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p09_fun ( n, p, f )

!*****************************************************************************80
!
!! P09_FUN evaluates the integrand for problem 9.
!
!  Integrand:
!
!    f(x,y) = 42 * x^5
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = 42.0D+00 * p(1,1:n)**5

  return
end
subroutine p09_title ( title )

!*****************************************************************************80
!
!! P09_TITLE returns the title of problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = 42 * x^5'

  return
end
subroutine p09_vertices ( t )

!*****************************************************************************80
!
!! P09_VERTICES returns the vertices for problem 9.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p10_fun ( n, p, f )

!*****************************************************************************80
!
!! P10_FUN evaluates the integrand for problem 10.
!
!  Discussion:
!
!    The integral has been transformed from the integral of G(X,Y)
!    over the unit reference triangle.
!
!  Integrand:
!
!    PA = 1
!    PB = 5
!    PC = 0
!    PD = 0
!    PG = 0.25
!    PH = -0.25
!    D = PB - PA
!    U(X) = ( X - PA ) / D
!    V1(X) = ( 1 - ( X - PA ) / D ) / ( ( PG - PC ) * X + PH - PD )
!    V(X,Y) = V1(X) * ( Y - PC * X - PD )
!
!    G(X,Y) = X^(-0.2)
!
!    c = 36/25
!
!    f(x,y) = c * g ( u(x), v(x,y) ) * v1(x) / d
!
!  Vertices:
!
!    (1,0), (5,0), (5,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Elise deDoncker, Ian Robinson,
!    Algorithm 612:
!    Integration over a Triangle Using Nonlinear Extrapolation,
!    ACM Transactions on Mathematical Software,
!    Volume 10, Number 1, March 1984, pages 17-22.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)
  real ( kind = 8 ), parameter :: pa = 1.0D+00
  real ( kind = 8 ), parameter :: pb = 5.0D+00
  real ( kind = 8 ), parameter :: pc = 0.0D+00
  real ( kind = 8 ), parameter :: pd = 0.0D+00
  real ( kind = 8 ), parameter :: pg = 0.25D+00
  real ( kind = 8 ), parameter :: ph = -0.25D+00
  real ( kind = 8 ), parameter :: power = -0.2D+00
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) v1(n)

  d = pb - pa
  u(1:n) = ( p(1,1:n) - pa ) / d
  v1(1:n) = ( 1.0D+00 - u(1:n) ) / ( ( pg - pc ) * p(1,1:n) + ph - pd )
  v(1:n) = v1(1:n) * ( p(2,1:n) - pc * p(1,1:n) - pd )
  f(1:n) = ( 36.0D+00 / 25.0D+00 ) * u(1:n)**power * v1(1:n) / d

  return
end
subroutine p10_title ( title )

!*****************************************************************************80
!
!! P10_TITLE returns the title of problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = x^(-0.2) on ((1,0),(5,0),(5,1))'

  return
end
subroutine p10_vertices ( t )

!*****************************************************************************80
!
!! P10_VERTICES returns the vertices for problem 10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
    1.0D+00, 0.0D+00, &
    5.0D+00, 0.0D+00, &
    5.0D+00, 1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p11_fun ( n, p, f )

!*****************************************************************************80
!
!! P11_FUN evaluates the integrand for problem 11.
!
!  Discussion:
!
!    The integral has been transformed from the integral of G(X,Y)
!    over the unit reference triangle.
!
!  Integrand:
!
!    PA = 0
!    PB = 1
!    PC = 0
!    PD = 0
!    PG = -1.0
!    PH = 1.0
!    D = PB - PA
!    U(X) = ( X - PA ) / D
!    V1(X) = ( 1 - ( X - PA ) / D ) / ( ( PG - PC ) * X + PH - PD )
!    V(X,Y) = V1(X) * ( Y - PC * X - PD )
!
!    G(X,Y) = (X+Y)**-0.2
!    c = 9 / 5
!
!    f(x,y) = c * g ( u(x), v(x,y) ) * v1(x) / d
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Elise deDoncker, Ian Robinson,
!    Algorithm 612:
!    Integration over a Triangle Using Nonlinear Extrapolation,
!    ACM Transactions on Mathematical Software,
!    Volume 10, Number 1, March 1984, pages 17-22.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)
  real ( kind = 8 ), parameter :: pa = 0.0D+00
  real ( kind = 8 ), parameter :: pb = 1.0D+00
  real ( kind = 8 ), parameter :: pc = 0.0D+00
  real ( kind = 8 ), parameter :: pd = 0.0D+00
  real ( kind = 8 ), parameter :: pg = -1.0D+00
  real ( kind = 8 ), parameter :: ph =  1.0D+00
  real ( kind = 8 ), parameter :: power = -0.2D+00
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) v1(n)

  d = pb - pa
  u(1:n) = ( p(1,1:n) - pa ) / d
  v1(1:n) = ( 1.0D+00 - u(1:n) ) / ( ( pg - pc ) * p(1,1:n) + ph - pd )
  v(1:n) = v1(1:n) * ( p(2,1:n) - pc * p(1,1:n) - pd )

  f(1:n) = ( 9.0D+00 / 5.0D+00 ) * ( u(1:n) + v(1:n) )**power * v1(1:n) / d

  return
end
subroutine p11_title ( title )

!*****************************************************************************80
!
!! P11_TITLE returns the title of problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = (x+y)^(-0.2)'

  return
end
subroutine p11_vertices ( t )

!*****************************************************************************80
!
!! P11_VERTICES returns the vertices for problem 11.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p12_fun ( n, p, f )

!*****************************************************************************80
!
!! P12_FUN evaluates the integrand for problem 12.
!
!  Discussion:
!
!    The integral has been transformed from the integral of G(X,Y)
!    over the unit reference triangle.
!
!  Integrand:
!
!    PA = -1
!    PB = 3
!    PC = 0.25
!    PD = -2.75
!    PG = -1.0
!    PH = 1.0
!    D = PB - PA
!    U(X) = ( X - PA ) / D
!    V1(X) = ( 1 - ( X - PA ) / D ) / ( ( PG - PC ) * X + PH - PD )
!    V(X,Y) = V1(X) * ( Y - PC * X - PD )
!
!    G(X,Y) = (1-X-Y)**-0.2
!    c = 36 / 25
!
!    f(x,y) = c * g ( u(x), v(x,y) ) * v1(x) / d
!
!  Vertices:
!
!    (-1,-3), (3,-2), (-1,2)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Elise deDoncker, Ian Robinson,
!    Algorithm 612:
!    Integration over a Triangle Using Nonlinear Extrapolation,
!    ACM Transactions on Mathematical Software,
!    Volume 10, Number 1, March 1984, pages 17-22.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)
  real ( kind = 8 ), parameter :: pa = -1.0D+00
  real ( kind = 8 ), parameter :: pb = 3.0D+00
  real ( kind = 8 ), parameter :: pc = 0.25D+00
  real ( kind = 8 ), parameter :: pd = -2.75D+00
  real ( kind = 8 ), parameter :: pg = -1.0D+00
  real ( kind = 8 ), parameter :: ph =  1.0D+00
  real ( kind = 8 ), parameter :: power = -0.2D+00
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) v1(n)

  d = pb - pa
  u(1:n) = ( p(1,1:n) - pa ) / d
  v1(1:n) = ( 1.0D+00 - u(1:n) ) / ( ( pg - pc ) * p(1,1:n) + ph - pd )
  v(1:n) = v1(1:n) * ( p(2,1:n) - pc * p(1,1:n) - pd )
  f(1:n) = ( 36.0D+00 / 25.0D+00 ) * &
    ( 1.0D+00 - u(1:n) - v(1:n) )**power * v1(1:n) / d

  return
end
subroutine p12_title ( title )

!*****************************************************************************80
!
!! P12_TITLE returns the title of problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = (1-x-y)^(-0.2) on ((-1,-3),(3,-2),(-1,2))'

  return
end
subroutine p12_vertices ( t )

!*****************************************************************************80
!
!! P12_VERTICES returns the vertices for problem 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
    -1.0D+00, -3.0D+00, &
     3.0D+00, -2.0D+00, &
    -1.0D+00,  2.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p13_fun ( n, p, f )

!*****************************************************************************80
!
!! P13_FUN evaluates the integrand for problem 13.
!
!  Discussion:
!
!    The integral has been transformed from the integral of G(X,Y)
!    over the unit reference triangle.
!
!  Integrand:
!
!    PA = 0
!    PB = -7
!    PC = 0
!    PD = 0
!    PG = -3/7
!    PH = -3
!    D = PB - PA
!    U(X) = ( X - PA ) / D
!    V1(X) = ( 1 - ( X - PA ) / D ) / ( ( PG - PC ) * X + PH - PD )
!    V(X,Y) = V1(X) * ( Y - PC * X - PD )
!
!    G(X,Y) = (X*Y)**-0.2
!    c = 0.94810264549557699446
!
!    f(x,y) = g ( u(x), v(x,y) ) * v1(x) / d / c
!
!  Vertices:
!
!    (0,0), (-7,0), (0,-3)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Elise deDoncker, Ian Robinson,
!    Algorithm 612:
!    Integration over a Triangle Using Nonlinear Extrapolation,
!    ACM Transactions on Mathematical Software,
!    Volume 10, Number 1, March 1984, pages 17-22.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: c = 0.94810264549557699446D+00
  real ( kind = 8 ) d
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)
  real ( kind = 8 ), parameter :: pa = 0.0D+00
  real ( kind = 8 ), parameter :: pb = -7.0D+00
  real ( kind = 8 ), parameter :: pc = 0.0D+00
  real ( kind = 8 ), parameter :: pd = 0.0D+00
  real ( kind = 8 ), parameter :: pg = ( -3.0D+00 / 7.0D+00 )
  real ( kind = 8 ), parameter :: ph =  -3.0D+00
  real ( kind = 8 ), parameter :: power = -0.2D+00
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) v1(n)

  d = pb - pa
  u(1:n) = ( p(1,1:n) - pa ) / d
  v1(1:n) = ( 1.0D+00 - u(1:n) ) / ( ( pg - pc ) * p(1,1:n) + ph - pd )
  v(1:n) = v1(1:n) * ( p(2,1:n) - pc * p(1,1:n) - pd )
  f(1:n) = ( u(1:n) * v(1:n) )**power * v1(1:n) / d / c

  return
end
subroutine p13_title ( title )

!*****************************************************************************80
!
!! P13_TITLE returns the title of problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = (x*y)^(-0.2) on ((0,0),(-7,0),(0,-3))'

  return
end
subroutine p13_vertices ( t )

!*****************************************************************************80
!
!! P13_VERTICES returns the vertices for problem 13.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
     0.0D+00,  0.0D+00, &
    -7.0D+00,  0.0D+00, &
     0.0D+00, -3.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p14_fun ( n, p, f )

!*****************************************************************************80
!
!! P14_FUN evaluates the integrand for problem 14.
!
!  Integrand:
!
!    f(x,y) = 3/10 * ( 1 / sqrt ( X ) + 1 / sqrt ( Y ) + 1 / sqrt ( X + Y ) )
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Elise deDoncker, Ian Robinson,
!    Algorithm 612:
!    Integration over a Triangle Using Nonlinear Extrapolation,
!    ACM Transactions on Mathematical Software,
!    Volume 10, Number 1, March 1984, pages 17-22.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = 3.0D+00 / 10.0D+00 * (      &
           1.0D+00 / sqrt ( p(1,1:n) ) &
         + 1.0D+00 / sqrt ( p(2,1:n) ) &
         + 1.0D+00 / sqrt ( p(1,1:n) + p(2,1:n) ) )

  return
end
subroutine p14_title ( title )

!*****************************************************************************80
!
!! P14_TITLE returns the title of problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = 1/sqrt(x) + 1/sqrt(y) + 1/sqrt(x+y)'

  return
end
subroutine p14_vertices ( t )

!*****************************************************************************80
!
!! P14_VERTICES returns the vertices for problem 14.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.0D+00,  1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p15_fun ( n, p, f )

!*****************************************************************************80
!
!! P15_FUN evaluates the integrand for problem 15.
!
!  Integrand:
!
!    f(x,y) = (3/4) / sqrt ( 1 - X - Y )
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Elise deDoncker, Ian Robinson,
!    Algorithm 612:
!    Integration over a Triangle Using Nonlinear Extrapolation,
!    ACM Transactions on Mathematical Software,
!    Volume 10, Number 1, March 1984, pages 17-22.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = 3.0D+00 / 4.0D+00 &
         / sqrt ( 1.0D+00 - p(1,1:n) - p(2,1:n) )

  return
end
subroutine p15_title ( title )

!*****************************************************************************80
!
!! P15_TITLE returns the title of problem 15.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = 1/sqrt(1-x-y)'

  return
end
subroutine p15_vertices ( t )

!*****************************************************************************80
!
!! P15_VERTICES returns the vertices for problem 15.
!
!  Modified:
!
!    28 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.0D+00,  1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p16_fun ( n, p, f )

!*****************************************************************************80
!
!! P16_FUN evaluates the integrand for problem 16.
!
!  Integrand:
!
!    f(x,y) = (-2/3) * log ( x * y )
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Elise deDoncker, Ian Robinson,
!    Algorithm 612:
!    Integration over a Triangle Using Nonlinear Extrapolation,
!    ACM Transactions on Mathematical Software,
!    Volume 10, Number 1, March 1984, pages 17-22.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = ( - 2.0D+00 / 3.0D+00 ) * log ( p(1,1:n) * p(2,1:n) )

  return
end
subroutine p16_title ( title )

!*****************************************************************************80
!
!! P16_TITLE returns the title of problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = log(x*y)'

  return
end
subroutine p16_vertices ( t )

!*****************************************************************************80
!
!! P16_VERTICES returns the vertices for problem 16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.0D+00,  1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p17_fun ( n, p, f )

!*****************************************************************************80
!
!! P17_FUN evaluates the integrand for problem 17.
!
!  Integrand:
!
!    f(x,y) = ( 1/sqrt(|x-1/4|) + 1/sqrt(|y-1/2|) ) / 3.11357229949
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Elise deDoncker, Ian Robinson,
!    Algorithm 612:
!    Integration over a Triangle Using Nonlinear Extrapolation,
!    ACM Transactions on Mathematical Software,
!    Volume 10, Number 1, March 1984, pages 17-22.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = ( 1.0D+00 / 3.11357229949D+00 ) * &
    ( 1.0D+00 / sqrt ( abs ( p(1,1:n) - 0.25D+00 ) ) &
    + 1.0D+00 / sqrt ( abs ( p(2,1:n) - 0.50D+00 ) ) )

  return
end
subroutine p17_title ( title )

!*****************************************************************************80
!
!! P17_TITLE returns the title of problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = 1/sqrt(|x-1/4|) + 1/sqrt(|y-1/2|)'

  return
end
subroutine p17_vertices ( t )

!*****************************************************************************80
!
!! P17_VERTICES returns the vertices for problem 17.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.0D+00,  1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p18_fun ( n, p, f )

!*****************************************************************************80
!
!! P18_FUN evaluates the integrand for problem 18.
!
!  Integrand:
!
!    f(x,y) = -4 * log ( x + y )
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Elise deDoncker, Ian Robinson,
!    Algorithm 612:
!    Integration over a Triangle Using Nonlinear Extrapolation,
!    ACM Transactions on Mathematical Software,
!    Volume 10, Number 1, March 1984, pages 17-22.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = -4.0D+00 * log ( p(1,1:n) + p(2,1:n) )

  return
end
subroutine p18_title ( title )

!*****************************************************************************80
!
!! P18_TITLE returns the title of problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = log ( x + y )'

  return
end
subroutine p18_vertices ( t )

!*****************************************************************************80
!
!! P18_VERTICES returns the vertices for problem 18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.0D+00,  1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p19_fun ( n, p, f )

!*****************************************************************************80
!
!! P19_FUN evaluates the integrand for problem 19.
!
!  Integrand:
!
!    f(x,y) = ( sin ( x ) * cos ( 5 * y ) ) / 0.043052326655855175018
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Elise deDoncker, Ian Robinson,
!    Algorithm 612:
!    Integration over a Triangle Using Nonlinear Extrapolation,
!    ACM Transactions on Mathematical Software,
!    Volume 10, Number 1, March 1984, pages 17-22.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) :: c = 0.043052326655855175018D+00
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = sin ( p(1,1:n) ) * cos ( 5.0D+00 * p(2,1:n) ) / c

  return
end
subroutine p19_title ( title )

!*****************************************************************************80
!
!! P19_TITLE returns the title of problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = sin ( x ) cos ( 5 y )'

  return
end
subroutine p19_vertices ( t )

!*****************************************************************************80
!
!! P19_VERTICES returns the vertices for problem 19.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.0D+00,  1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p20_fun ( n, p, f )

!*****************************************************************************80
!
!! P20_FUN evaluates the integrand for problem 20.
!
!  Integrand:
!
!    f(x,y) = ( sin ( 11 x ) * cos ( y ) ) / 0.085468091995313041919D+00
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Elise deDoncker, Ian Robinson,
!    Algorithm 612:
!    Integration over a Triangle Using Nonlinear Extrapolation,
!    ACM Transactions on Mathematical Software,
!    Volume 10, Number 1, March 1984, pages 17-22.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) :: c = 0.085468091995313041919D+00
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = sin ( 11.0D+00 * p(1,1:n) ) * cos ( p(2,1:n) ) / c

  return
end
subroutine p20_title ( title )

!*****************************************************************************80
!
!! P20_TITLE returns the title of problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = sin ( 11 x ) cos ( y )'

  return
end
subroutine p20_vertices ( t )

!*****************************************************************************80
!
!! P20_VERTICES returns the vertices for problem 20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.0D+00,  1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p21_fun ( n, p, f )

!*****************************************************************************80
!
!! P21_FUN evaluates the integrand for problem 21.
!
!  Discussion:
!
!    To do this integral by hand, convert to polar coordinates:
!
!    Integral ( 0 <= t <= Pi/2 )
!      Integral ( 0 <= r <= 1/(cos(t)+sin(t)) ) 1/r * r dr dt
!
!  Integrand:
!
!    f(x,y) = 1 / sqrt ( x * x + y * y ) / 1.246450480...
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Elise deDoncker, Ian Robinson,
!    Algorithm 612:
!    Integration over a Triangle Using Nonlinear Extrapolation,
!    ACM Transactions on Mathematical Software,
!    Volume 10, Number 1, March 1984, pages 17-22.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: c = 1.2464504802804610268D+00
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = 1.0D+00 / sqrt ( p(1,1:n)**2 + p(2,1:n)**2 ) / c

  return
end
subroutine p21_title ( title )

!*****************************************************************************80
!
!! P21_TITLE returns the title of problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = 1 / r = 1 / sqrt ( x^2 + y^2 )'

  return
end
subroutine p21_vertices ( t )

!*****************************************************************************80
!
!! P21_VERTICES returns the vertices for problem 21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.0D+00,  1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine p22_fun ( n, p, f )

!*****************************************************************************80
!
!! P22_FUN evaluates the integrand for problem 22.
!
!  Discussion:
!
!    To do this integral by hand, convert to polar coordinates:
!
!    Integral ( 0 <= t <= Pi/2 )
!      Integral ( 0 <= r <= 1/(cos(t)+sin(t))) Log(r)/r * r dr dt
!
!  Integrand:
!
!    f(x,y) = log ( r ) / r / (-1.5280234546641884580)
!
!  Vertices:
!
!    (0,0), (1,0), (0,1)
!
!  Integral:
!
!    1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Elise deDoncker, Ian Robinson,
!    Algorithm 612:
!    Integration over a Triangle Using Nonlinear Extrapolation,
!    ACM Transactions on Mathematical Software,
!    Volume 10, Number 1, March 1984, pages 17-22.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) F(N), the function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ), parameter :: c = -1.5280234546641884580D+00
  real ( kind = 8 ) f(n)
  real ( kind = 8 ) p(2,n)

  f(1:n) = log ( sqrt ( p(1,1:n)**2 + p(2,1:n)**2 ) ) &
               / sqrt ( p(1,1:n)**2 + p(2,1:n)**2 ) / c

  return
end
subroutine p22_title ( title )

!*****************************************************************************80
!
!! P22_TITLE returns the title of problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) TITLE, the title of the problem.
!
  implicit none

  character ( len = * ) title

  title = 'f(x,y) = log ( r ) / r'

  return
end
subroutine p22_vertices ( t )

!*****************************************************************************80
!
!! P22_VERTICES returns the vertices for problem 22.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) T(2,3), the vertices.
!
  implicit none

  real ( kind = 8 ) t(2,3)

  t(1:2,1:3) = reshape ( (/ &
     0.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00, &
     0.0D+00,  1.0D+00 /), (/ 2, 3 /) )

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of real ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + huge ( seed )
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine random_initialize ( seed_input )

!*****************************************************************************80
!
!! RANDOM_INITIALIZE initializes the FORTRAN90 random number seed.
!
!  Discussion:
!
!    If you don't initialize the FORTRAN90 random number generator
!    routine RANDOM_NUMBER, which is used by calls like
!
!      call random_number ( harvest = x )
!
!    then its behavior is not specified.  That may be OK for you.  But
!    you may want to be able to force the same sequence to be generated
!    each time, or to force a different sequence.
!
!    To control the sequence of random numbers, you need to set the seed.
!    In FORTRAN90, this is done by calling the RANDOM+SEED routine.
!    You can call it with no arguments, in fact.  But if you call
!    it with no arguments:
!
!      call random_seed ( )
!
!    then its behavior (or more particularly, the behavior of RANDOM_NUMBER)
!    is still not specified.  You might hope that the system will go out
!    and pick a nice random seed for you, but there's no guarantee.
!
!
!    For example, on the DEC ALPHA, if you compile a program that calls
!    RANDOM_NUMBER, then every time you run it, you get the same sequence
!    of "random" values.  If you compile a program that calls RANDOM_SEED
!    with no arguments, and then calls RANDOM_NUMBER, you still get the
!    same sequence each time you run the program.
!
!    In order to actually try to scramble up the random number generator
!    a bit, this routine goes through the tedious process of getting the
!    size of the random number seed, making up values based on the current
!    time, and setting the random number seed.
!
!    Unfortunately, the RANDOM_SEED routine has a very elastic definition.
!    It does not use a single scalar integer SEED.  Instead, it communicates
!    with the user through an integer vector whose size is not specified.
!    You actually have to "ask" the routine to tell you the size of this
!    vector.  Then you can fill up the vector with values to be used to
!    influence the seeding of the random number routine.  The details of
!    how the seed affects the sequence are also unspecified, but the only
!    thing we are reasonably confident about is that specifying the same
!    seed should result in the same sequence, and specifying different
!    seeds should result in different sequences!
!
!    I assume this is far more than you wanted to know.  (It's certainly
!    more than I wanted to know!)
!
!    The argument SEED is an input quantity only, so it is legal
!    to type
!
!      call random_initialize ( 0 )
!
!    or
!
!      call random_initialize ( 18867 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED_INPUT, the user's "suggestion" for a seed.
!    However, if the input value is 0, the routine will come up with
!    its own "suggestion", based on the system clock.
!
  implicit none

  integer ( kind = 4 ) count
  integer ( kind = 4 ) count_max
  integer ( kind = 4 ) count_rate
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_input
  integer ( kind = 4 ), allocatable :: seed_vector(:)
  integer ( kind = 4 ) seed_size
  real ( kind = 8 ) t
  integer ( kind = 4 ), parameter :: warm_up = 100

  seed = seed_input
!
!  Initialize the random seed routine.
!
  call random_seed ( )
!
!  Determine the size of the random number seed vector.
!
  call random_seed ( size = seed_size )
!
!  Allocate a vector of the right size to be used as a random seed.
!
  allocate ( seed_vector(seed_size) )
!
!  If the user supplied a SEED value, use that.
!
!  Otherwise, use the system clock value to make up a value that is
!  likely to change based on when this routine is called.
!
  if ( seed /= 0 ) then

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, user SEED = ', seed
    end if

  else

    call system_clock ( count, count_rate, count_max )

    seed = count

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, arbitrary SEED = ', &
        seed
    end if

  end if
!
!  Set the seed vector.  We don't know the significance of the
!  individual entries of the internal seed vector, so we'll just set
!  all entries to SEED.
!
  seed_vector(1:seed_size) = seed
!
!  Now call RANDOM_SEED, and tell it to use this seed vector.
!
  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Call the random number routine a bunch of times just to "warm it up".
!
  do i = 1, warm_up
    call random_number ( harvest = t )
  end do

  return
end
subroutine subtriangle_next ( n, more, i1, j1, i2, j2, i3, j3 )

!*****************************************************************************80
!
!! SUBTRIANGLE_NEXT computes the next subtriangle of a triangle.
!
!  Discussion:
!
!    The three sides of a triangle have been subdivided into N segments,
!    inducing a natural subdivision of the triangle into N*N subtriangles.
!    It is desired to consider each subtriangle, one at a time, in some
!    definite order.  This routine can produce information defining each
!    of the subtriangles, one after another.
!
!    The subtriangles are described in terms of the integer coordinates
!    (I,J) of their vertices.  These coordinates both range from 0 to N,
!    with the additional restriction that I + J <= N.
!
!    The vertices of each triangle are listed in counterclockwise order.
!
!  Example:
!
!    N = 4
!
!    4  *
!       |\
!       16\
!    3  *--*
!       |14|\
!       13\15\
!    2  *--*--*
!       |\9|11|\
!       |8\10\12\
!    1  *--*--*--*
!       |\2|\4|\6|\
!       |1\|3\|5\|7\
!   0   *--*--*--*--*
!
!       0  1  2  3  4
!
!    Rank  I1 J1  I2 J2  I3 J3
!    ----  -----  -----  -----
!       1   0  0   1  0   0  1
!       2   1  1   0  1   1  0
!       3   1  0   2  0   1  1
!       4   2  1   1  1   2  0
!       5   2  0   3  0   2  1
!       6   3  1   1  1   3  0
!       7   3  0   4  0   3  1
!       8   0  1   1  1   0  2
!       9   1  2   0  2   1  1
!      10   1  1   2  1   1  2
!      11   2  2   1  2   2  1
!      12   2  1   3  1   2  2
!      13   0  2   1  2   0  3
!      14   1  3   0  3   1  2
!      15   1  2   2  2   1  3
!      16   0  3   1  3   0  4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, indicates the number of subdivisions of each side
!    of the original triangle.
!
!    Input/output, logical MORE.
!    On first call, set MORE to FALSE.  Thereafter, the output value of MORE
!    will be TRUE if there are more subtriangles that can be generated by
!    further calls.  However, if MORE is returned as FALSE, the accompanying
!    subtriangle information refers to the last subtriangle that can be
!    generated.
!
!    Input/output, integer ( kind = 4 ) I1, J1, I2, J2, I3, J3, the indices of the
!    vertices of the subtriangle.
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  logical more
  integer ( kind = 4 ) n

  if ( n <= 0 ) then
    more = .false.
    return
  end if

  if ( .not. more ) then

    i1 = 0
    j1 = 0
    i2 = 1
    j2 = 0
    i3 = 0
    j3 = 1

    if ( n == 1 ) then
      more = .false.
    else
      more = .true.
    end if
!
!  We last generated a triangle like:
!
!    2---1
!     \  |
!      \ |
!       \|
!        3
!
  else if ( i2 < i3 ) then

    i1 = i3
    j1 = j3
    i2 = i1 + 1
    j2 = j1
    i3 = i1
    j3 = j1 + 1
!
!  We last generated a triangle like
!
!    3
!    |\
!    | \
!    |  \
!    1---2
!
  else if ( i1 + 1 + j1 + 1 <= n ) then

    i1 = i1 + 1
    j1 = j1 + 1
    i2 = i1 - 1
    j2 = j1
    i3 = i1
    j3 = j1 - 1
!
!  We must be at the end of a row.
!
  else

    i1 = 0
    j1 = j1 + 1
    i2 = i1 + 1
    j2 = j1
    i3 = i1
    j3 = j1 + 1

    if ( n <= j1 + 1 ) then
      more = .false.
    end if

  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2001
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

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine triangle_area ( t, area )

!*****************************************************************************80
!
!! TRIANGLE_AREA computes the area of a triangle in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = 8 ) AREA, the absolute area of the triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) area
  real ( kind = 8 ) t(dim_num,3)

  area = 0.5D+00 * abs ( &
      t(1,1) * ( t(2,2) - t(2,3) ) &
    + t(1,2) * ( t(2,3) - t(2,1) ) &
    + t(1,3) * ( t(2,1) - t(2,2) ) )

  return
end
subroutine triangle_sample ( t, n, seed, p )

!*****************************************************************************80
!
!! TRIANGLE_SAMPLE returns random points in a triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Input, integer ( kind = 4 ) N, the number of points to generate.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) P(2,N), random points in the triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha(n)
  integer ( kind = 4 ) dim
  real ( kind = 8 ) p(dim_num,n)
  real ( kind = 8 ) p12(dim_num,n)
  real ( kind = 8 ) p13(dim_num,n)
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t(dim_num,3)
!
!  For comparison between F90, C++ and MATLAB codes, call R8VEC_UNIFORM_01.
!
! call r8vec_uniform_01 ( n, seed, alpha )
!
!  For faster execution, call RANDOM_NUMBER.
!
  call random_number ( harvest = alpha(1:n) )
!
!  Interpret R as a percentage of the triangle's area.
!
!  Imagine a line L, parallel to side 1, so that the area between
!  vertex 1 and line L is R percent of the full triangle's area.
!
!  The line L will intersect sides 2 and 3 at a fraction
!  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
!
  alpha(1:n) = sqrt ( alpha(1:n) )
!
!  Determine the coordinates of the points on sides 2 and 3 intersected
!  by line L.
!
  do dim = 1, dim_num

    p12(dim,1:n) = ( 1.0D+00 - alpha(1:n) ) * t(dim,1) &
                             + alpha(1:n)   * t(dim,2)

    p13(dim,1:n) = ( 1.0D+00 - alpha(1:n) ) * t(dim,1) &
                             + alpha(1:n)   * t(dim,3)

  end do
!
!  Now choose, uniformly at random, a point on the line L.
!
!  For comparison between F90, C++ and MATLAB codes, call R8VEC_UNIFORM_01.
!
! call r8vec_uniform_01 ( n, seed, alpha )
!
!  For faster execution, call RANDOM_NUMBER.
!
  call random_number ( harvest = alpha(1:n) )

  do dim = 1, dim_num

    p(dim,1:n) = ( 1.0D+00 - alpha(1:n) ) * p12(dim,1:n) &
                           + alpha(1:n)   * p13(dim,1:n)

  end do

  return
end
subroutine triangle_sample_old ( t, n, seed, p )

!*****************************************************************************80
!
!! TRIANGLE_SAMPLE_OLD returns random points in a triangle.
!
!  Discussion:
!
!    This version is not optimized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Input, integer ( kind = 4 ) N, the number of points to generate.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) P(2,N), random points in the triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(dim_num,n)
  real ( kind = 8 ) p12(dim_num)
  real ( kind = 8 ) p13(dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t(dim_num,3)

  do j = 1, n

    r = r8_uniform_01 ( seed )
!
!  Interpret R as a percentage of the triangle's area.
!
!  Imagine a line L, parallel to side 1, so that the area between
!  vertex 1 and line L is R percent of the full triangle's area.
!
!  The line L will intersect sides 2 and 3 at a fraction
!  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
!
    alpha = sqrt ( r )
!
!  Determine the coordinates of the points on sides 2 and 3 intersected
!  by line L.
!
    p12(1:dim_num) = ( 1.0D+00 - alpha ) * t(1:dim_num,1) &
                               + alpha   * t(1:dim_num,2)

    p13(1:dim_num) = ( 1.0D+00 - alpha ) * t(1:dim_num,1) &
                               + alpha   * t(1:dim_num,3)
!
!  Now choose, uniformly at random, a point on the line L.
!
    beta = r8_uniform_01 ( seed )

    p(1:dim_num,j) = ( 1.0D+00 - beta ) * p12(1:dim_num) &
                               + beta   * p13(1:dim_num)

  end do

  return
end
