program main

!*****************************************************************************80
!
!! MAIN is the main program for SGMGA_SIZE_TABLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 July 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) ctime1
  real ( kind = 8 ) ctime2
  integer ( kind = 4 ) dim_max
  integer ( kind = 4 ) dim_min
  integer ( kind = 4 ) growth_1d
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) np_1d
  real ( kind = 8 ), allocatable :: p_1d(:)
  integer ( kind = 4 ) rule_1d

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_SIZE_TABLE:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Make tables of point counts.'
  write ( *, '(a)' ) '  Measure the CPU time for each table.'
!
!  Clenshaw-Curtis Grid (1), slow exponential growth (4).
!
  rule_1d = 1
  growth_1d = 4
  np_1d = 0
  allocate ( p_1d(1:np_1d) )
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 7
  call cpu_time ( ctime1 )
  call sgmga_size_tabulate ( rule_1d, growth_1d, np_1d, p_1d, dim_min, &
    dim_max, level_max_min, level_max_max )
  call cpu_time ( ctime2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  CPU Time = ', ctime2 - ctime1
  deallocate ( p_1d )
!
!  Clenshaw-Curtis Grid (1), exponential growth (6).
!
  rule_1d = 1
  growth_1d = 6
  np_1d = 0
  allocate ( p_1d(1:np_1d) )
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 7
  call cpu_time ( ctime1 )
  call sgmga_size_tabulate ( rule_1d, growth_1d, np_1d, p_1d, dim_min, &
    dim_max, level_max_min, level_max_max )
  call cpu_time ( ctime2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  CPU Time = ', ctime2 - ctime1
  deallocate ( p_1d )
!
!  Clenshaw-Curtis Grid (1), exponential growth (6).
!
  rule_1d = 1
  growth_1d = 6
  np_1d = 0
  allocate ( p_1d(1:np_1d) )
  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 7
  call cpu_time ( ctime1 )
  call sgmga_size_tabulate ( rule_1d, growth_1d, np_1d, p_1d, dim_min, &
    dim_max, level_max_min, level_max_max )
  call cpu_time ( ctime2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  CPU Time = ', ctime2 - ctime1
  deallocate ( p_1d )
!
!  Clenshaw-Curtis Grid (1), exponential growth (6).
!
  rule_1d = 1
  growth_1d = 6
  np_1d = 0
  allocate ( p_1d(1:np_1d) )
  dim_min = 100
  dim_max = 100
  level_max_min = 0
  level_max_max = 2
  call cpu_time ( ctime1 )
  call sgmga_size_tabulate ( rule_1d, growth_1d, np_1d, p_1d, dim_min, &
   dim_max, level_max_min, level_max_max )
  call cpu_time ( ctime2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  CPU Time = ', ctime2 - ctime1
  deallocate ( p_1d )
!
!  Gauss Patterson Grid (3), slow exponential growth (4).
!
  rule_1d = 3
  growth_1d = 4
  np_1d = 0
  allocate ( p_1d(1:np_1d) )
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 7
  call cpu_time ( ctime1 )
  call sgmga_size_tabulate ( rule_1d, growth_1d, np_1d, p_1d, dim_min, &
    dim_max, level_max_min, level_max_max )
  call cpu_time ( ctime2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  CPU Time = ', ctime2 - ctime1
  deallocate ( p_1d )
!
!  Gauss Patterson Grid (3), moderate exponential growth (5).
!
  rule_1d = 3
  growth_1d = 5
  np_1d = 0
  allocate ( p_1d(1:np_1d) )
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 7
  call cpu_time ( ctime1 )
  call sgmga_size_tabulate ( rule_1d, growth_1d, np_1d, p_1d, dim_min, &
    dim_max, level_max_min, level_max_max )
  call cpu_time ( ctime2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  CPU Time = ', ctime2 - ctime1
  deallocate ( p_1d )
!
!  Gauss Patterson Grid (3), exponential growth (6).
!
  rule_1d = 3
  growth_1d = 6
  np_1d = 0
  allocate ( p_1d(1:np_1d) )
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 6
  call cpu_time ( ctime1 )
  call sgmga_size_tabulate ( rule_1d, growth_1d, np_1d, p_1d, dim_min, &
    dim_max, level_max_min, level_max_max )
  call cpu_time ( ctime2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  CPU Time = ', ctime2 - ctime1
  deallocate ( p_1d )
!
!  Gauss Legendre grid (4), slow linear odd growth (2)
!
  rule_1d = 4
  growth_1d = 2
  np_1d = 0
  allocate ( p_1d(1:np_1d) )
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 7
  call cpu_time ( ctime1 )
  call sgmga_size_tabulate ( rule_1d, growth_1d, np_1d, p_1d, dim_min, &
    dim_max, level_max_min, level_max_max )
  call cpu_time ( ctime2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  CPU Time = ', ctime2 - ctime1
  deallocate ( p_1d )
!
!  Gauss Legendre grid (4), moderate linear growth (3).
!
  rule_1d = 4
  growth_1d = 3
  np_1d = 0
  allocate ( p_1d(1:np_1d) )
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 7
  call cpu_time ( ctime1 )
  call sgmga_size_tabulate ( rule_1d, growth_1d, np_1d, p_1d, dim_min, &
    dim_max, level_max_min, level_max_max )
  call cpu_time ( ctime2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  CPU Time = ', ctime2 - ctime1
  deallocate ( p_1d )
!
!  Gauss Laguerre grid (7), moderate linear growth (3).
!
  rule_1d = 7
  growth_1d = 3
  np_1d = 0
  allocate ( p_1d(1:np_1d) )
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 7
  call cpu_time ( ctime1 )
  call sgmga_size_tabulate ( rule_1d, growth_1d, np_1d, p_1d, dim_min, &
    dim_max, level_max_min, level_max_max )
  call cpu_time ( ctime2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  CPU Time = ', ctime2 - ctime1
  deallocate ( p_1d )
!
!  Hermite Genz Keister (10), slow exponential growth (4).
!
  rule_1d = 10
  growth_1d = 4
  np_1d = 0
  allocate ( p_1d(1:np_1d) )
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 7
  call cpu_time ( ctime1 )
  call sgmga_size_tabulate ( rule_1d, growth_1d, np_1d, p_1d, dim_min, &
    dim_max, level_max_min, level_max_max )
  call cpu_time ( ctime2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  CPU Time = ', ctime2 - ctime1
  deallocate ( p_1d )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_SIZE_TABLE:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine sgmga_size_tabulate ( rule_1d, growth_1d, np_1d, p_1d, &
   dim_min, dim_max, level_max_min, level_max_max )

!****************************************************************************80
!
!! SGMGA_SIZE_TABULATE tests SGMGA_SIZE.
!
!  Discussion:
!
!    We do NOT consider mixed rules.  Instead, we are looking at sparse grid
!    rules for which all dimensions use the same 1D rule family.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE_1D, the 1D rule.
!     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
!     2, "F2",  Fejer Type 2, Open Fully Nested.
!     3, "GP",  Gauss Patterson, Open Fully Nested.
!     4, "GL",  Gauss Legendre, Open Weakly Nested.
!     5, "GH",  Gauss Hermite, Open Weakly Nested.
!     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
!     7, "LG",  Gauss Laguerre, Open Non Nested.
!     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
!     9, "GJ",  Gauss Jacobi, Open Non Nested.
!    10, "HGK", Hermite Genz-Keister, Open Fully Nested.
!    11, "UO",  User supplied Open, presumably Non Nested.
!    12, "UC",  User supplied Closed, presumably Non Nested.
!
!    Input, integer ( kind = 4 ) GROWTH_1D, the growth in each dimension.
!    0, "DF", default growth associated with this quadrature rule;
!    1, "SL", slow linear, L+1;
!    2  "SO", slow linear odd, O=1+2((L+1)/2)
!    3, "ML", moderate linear, 2L+1;
!    4, "SE", slow exponential;
!    5, "ME", moderate exponential;
!    6, "FE", full exponential.
!
!    Input, integer ( kind -= 4 ) NP_1D, the number of parameters in the
!    1D rule.
!
!    Input, real ( kind = 8 ) P_1D(NP_1D), the parameters.
!
!    Input, integer ( kind = 4 ) DIM_MIN, the minimum spatial dimension.
!
!    Input, integer ( kind = 4 ) DIM_MAX, the maximum spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX_MIN, the minimum value of LEVEL_MAX.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX_MAX, the maximum value of LEVEL_MAX.
!
  implicit none

  integer ( kind = 4 ) dim_max
  integer ( kind = 4 ) dim_min
  integer ( kind = 4 ) np_1d

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), allocatable :: growth(:)
  integer ( kind = 4 ) growth_1d
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  real ( kind = 8 ), allocatable :: level_weight(:)
  integer ( kind = 4 ) np_sum
  integer ( kind = 4 ), allocatable :: np(:)
  real ( kind = 8 ), allocatable :: p(:)
  real ( kind = 8 ) p_1d(np_1d)
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ), allocatable :: point_vec(:)
  real ( kind = 8 ) r8_epsilon
  integer ( kind = 4 ) rule_1d
  integer ( kind = 4 ), allocatable :: rule(:)
  real ( kind = 8 ) tol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_SIZE_TABULATE'
  write ( *, '(a)' ) '  SGMGA_SIZE returns the number of distinct'
  write ( *, '(a)' ) '  points in a sparse grid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We use the same rule in all dimensions, and count the points'
  write ( *, '(a)' ) '  for a range of dimensions and levels.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  1D rule index is  ', rule_1d
  write ( *, '(a,i8)' ) '  1D growth rule is ', growth_1d
  write ( *, '(a)' ) ' '

  tol = sqrt ( r8_epsilon ( ) )

  allocate ( point_vec(dim_min:dim_max) )

  do dim_num = dim_min, dim_max
    point_vec(dim_num) = dim_num
  end do

  write ( *, '(a8)', ADVANCE = 'NO' ) '   DIM: '
  write ( *, '(5(2x,i8))' ) point_vec(dim_min:dim_max)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max

    do dim_num = dim_min, dim_max

      allocate ( growth(1:dim_num) )
      allocate ( level_weight(1:dim_num) )
      allocate ( rule(1:dim_num) )
      allocate ( np(1:dim_num) )
      np_sum = dim_num * np_1d
      allocate ( p(1:np_sum) )

      level_weight(1:dim_num) = 1.0D+00
      rule(1:dim_num) = rule_1d
      growth(1:dim_num) = growth_1d
      np(1:dim_num) = np_1d
      do dim = 1, dim_num
        p(1+(dim-1)*np_1d:np_1d+(dim-1)*np_1d) = p_1d(1:np_1d)
      end do

      call sgmga_size ( dim_num, level_weight, level_max, rule, growth, &
        np, p, tol, point_num )

      point_vec(dim_num) = point_num

      deallocate ( growth )
      deallocate ( level_weight )
      deallocate ( np )
      deallocate ( p )
      deallocate ( rule )

    end do

    write ( *, '(4x,i4)', ADVANCE = 'NO' ) level_max
    write ( *, '(5(2x,i8))' ) point_vec(dim_min:dim_max)

  end do

  deallocate ( point_vec )

  return
end
