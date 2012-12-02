program main

!*****************************************************************************80
!
!! MAIN is the main program for SPARSE_GRID_OPEN_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
  implicit none

  integer ( kind = 4 ) dim_max
  integer ( kind = 4 ) dim_min
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_PURE_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SPARSE_GRID_PURE library.'
!
!  CC_ME
!
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10
  call test005 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 10
  call test005 ( dim_min, dim_max, level_max_min, level_max_max )
!
!  CC_SE
!
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10
  call test01 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 10
  call test01 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 11
  dim_max = 15
  level_max_min = 0
  level_max_max = 5
  call test01 ( dim_min, dim_max, level_max_min, level_max_max )
!
!  CFN_E
!
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10
  call test02 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 10
  call test02 ( dim_min, dim_max, level_max_min, level_max_max )
!
!  F2_SE
!
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10
  call test03 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 10
  call test03 ( dim_min, dim_max, level_max_min, level_max_max )
!
!  GP_ME
!
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10
  call test035 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 10
  call test035 ( dim_min, dim_max, level_max_min, level_max_max )
!
!  GP_SE
!
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10
  call test04 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 10
  call test04 ( dim_min, dim_max, level_max_min, level_max_max )
!
!  OFN_E
!
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10
  call test05 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 10
  call test05 ( dim_min, dim_max, level_max_min, level_max_max )
!
!  ONN_E
!
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10
  call test06 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 10
  call test06 ( dim_min, dim_max, level_max_min, level_max_max )
!
!  ONN_L
!
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10
  call test07 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 10
  call test07 ( dim_min, dim_max, level_max_min, level_max_max )
!
! OWN_E
!
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10
  call test08 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 10
  call test08 ( dim_min, dim_max, level_max_min, level_max_max )
!
!  OWN_L
!
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10
  call test09 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 10
  call test09 ( dim_min, dim_max, level_max_min, level_max_max )
!
!  OWN_LS
!
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 10
  call test10 ( dim_min, dim_max, level_max_min, level_max_max )

  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 10
  call test10 ( dim_min, dim_max, level_max_min, level_max_max )
!
!  Terminate
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_PURE_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test005 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST005 tests SPARSE_GRID_CC_ME_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005'
  write ( *, '(a)' ) '  SPARSE_GRID_CC_ME_SIZE returns the number of '
  write ( *, '(a)' ) '  distinct points in a sparse grid, made from'
  write ( *, '(a)' ) '  * CC_ME, Clenshaw Curtis Moderate Exponential Growth family.'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i10))' ) '   DIM: ', point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max
    do dim_num = dim_min, dim_max
      call sparse_grid_cc_me_size ( dim_num, level_max, point_num(dim_num) )
    end do
    write ( *, '(a4,i4,6(2x,i10))' ) '    ', level_max, point_num(dim_min:dim_max)
  end do

  return
end
subroutine test01 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST01 tests SPARSE_GRID_CC_SE_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  SPARSE_GRID_CC_SE_SIZE returns the number of '
  write ( *, '(a)' ) '  distinct points in a sparse grid, made from'
  write ( *, '(a)' ) '  * CC_SE, Clenshaw Curtis Slow Exponential Growth family.'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i10))' ) '   DIM: ', point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max
    do dim_num = dim_min, dim_max
      call sparse_grid_cc_se_size ( dim_num, level_max, point_num(dim_num) )
    end do
    write ( *, '(a4,i4,6(2x,i10))' ) '    ', level_max, point_num(dim_min:dim_max)
  end do

  return
end
subroutine test02 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST02 tests SPARSE_GRID_CFN_E_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  SPARSE_GRID_CFN_E_SIZE returns the number of '
  write ( *, '(a)' ) '  distinct points in a CFN_E sparse grid, made from'
  write ( *, '(a)' ) '  any closed fully nested family of 1D quadrature'
  write ( *, '(a)' ) '  rules with exponential growth, including:'
  write ( *, '(a)' ) '  * CC_E, the Clenshaw Curtis Exponential Growth family;'
  write ( *, '(a)' ) '  * NCC_E, the Newton Cotes Closed Exponential Growth family.'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i10))' ) '   DIM: ', point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max
    do dim_num = dim_min, dim_max
      call sparse_grid_cfn_e_size ( dim_num, level_max, point_num(dim_num) )
    end do
    write ( *, '(a4,i4,6(2x,i10))' ) '    ', level_max, point_num(dim_min:dim_max)
  end do

  return
end
subroutine test03 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST03 tests SPARSE_GRID_F2_SE_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  SPARSE_GRID_F2_SE_SIZE returns the number of '
  write ( *, '(a)' ) '  distinct points in a sparse grid, made from'
  write ( *, '(a)' ) '  * F2_SE, the Fejer Type 2 Slow Exponential Growth family.'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i10))' ) '   DIM: ', point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max
    do dim_num = dim_min, dim_max
      call sparse_grid_f2_se_size ( dim_num, level_max, point_num(dim_num) )
    end do
    write ( *, '(a4,i4,6(2x,i10))' ) '    ', level_max, point_num(dim_min:dim_max)
  end do

  return
end
subroutine test035 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST035 tests SPARSE_GRID_GP_ME_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST035'
  write ( *, '(a)' ) '  SPARSE_GRID_GP_ME_SIZE returns the number of '
  write ( *, '(a)' ) '  distinct points in a sparse grid made from'
  write ( *, '(a)' ) '  * GP_ME, Gauss Patterson Moderate Exponential Growth family.'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i10))' ) '   DIM: ', point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max
    do dim_num = dim_min, dim_max
      call sparse_grid_gp_me_size ( dim_num, level_max, point_num(dim_num) )
    end do
    write ( *, '(a4,i4,6(2x,i10))' ) '    ', level_max, point_num(dim_min:dim_max)
  end do

  return
end
subroutine test04 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST04 tests SPARSE_GRID_GP_SE_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  SPARSE_GRID_GP_SE_SIZE returns the number of '
  write ( *, '(a)' ) '  distinct points in a sparse grid made from'
  write ( *, '(a)' ) '  * GP_SE, Gauss Patterson Slow Exponential Growth family.'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i10))' ) '   DIM: ', point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max
    do dim_num = dim_min, dim_max
      call sparse_grid_gp_se_size ( dim_num, level_max, point_num(dim_num) )
    end do
    write ( *, '(a4,i4,6(2x,i10))' ) '    ', level_max, point_num(dim_min:dim_max)
  end do

  return
end
subroutine test05 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST05 tests SPARSE_GRID_OFN_E_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  SPARSE_GRID_OFN_E_SIZE returns the number of '
  write ( *, '(a)' ) '  distinct points in an OFN_E sparse grid made up of all '
  write ( *, '(a)' ) '  product grids formed from open fully nested '
  write ( *, '(a)' ) '  quadrature rules with Exponential Growth, including'
  write ( *, '(a)' ) '  * F2_E, the Fejer Type 2 Exponential Growth Family;'
  write ( *, '(a)' ) '  * GP_E, the Gauss Patterson Exponential Growth Family;'
  write ( *, '(a)' ) '  * NCO_E, the Newton Cotes Open Exponential Growth Family.'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i10))' ) '   DIM: ', point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max
    do dim_num = dim_min, dim_max
      call sparse_grid_ofn_e_size ( dim_num, level_max, point_num(dim_num) )
    end do
    write ( *, '(a4,i4,6(2x,i10))' ) '    ', level_max, point_num(dim_min:dim_max)
  end do

  return
end
subroutine test06 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST06 tests SPARSE_GRID_ONN_E_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  SPARSE_GRID_ONN_E_SIZE returns the number of '
  write ( *, '(a)' ) '  distinct points in an ONN_E sparse grid made up of all '
  write ( *, '(a)' ) '  product grids formed from open non-nested '
  write ( *, '(a)' ) '  quadrature rules with exponential growth, including:'
  write ( *, '(a)' ) '  * LG_E, the Gauss Laguerre Exponential Growth Family;'
  write ( *, '(a)' ) '  * GJ_E, the Gauss Jacobi Exponential Growth Family;'
  write ( *, '(a)' ) '  * GLG_E, the Generalized Gauss Laguerre Exponential Growth Family;'
  write ( *, '(a)' ) '  * GW_E, any Golub Welsch Exponential Growth Family;'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i10))' ) '   DIM: ', point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max
    do dim_num = dim_min, dim_max
      call sparse_grid_onn_e_size ( dim_num, level_max, point_num(dim_num) )
    end do
    write ( *, '(a4,i4,6(2x,i10))' ) '    ', level_max, point_num(dim_min:dim_max)
  end do

  return
end
subroutine test07 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST07 tests SPARSE_GRID_ONN_L_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  SPARSE_GRID_ONN_L_SIZE returns the number of '
  write ( *, '(a)' ) '  distinct points in an ONN_L sparse grid made up of all '
  write ( *, '(a)' ) '  product grids formed from open non-nested '
  write ( *, '(a)' ) '  quadrature rules with linear growth, including:'
  write ( *, '(a)' ) '  * LG_L, the Gauss Laguerre Linear Growth Family;'
  write ( *, '(a)' ) '  * GJ_L, the Gauss Jacobi Linear Growth Family;'
  write ( *, '(a)' ) '  * GLG_L, the Generalized Gauss Laguerre Linear Growth Family;'
  write ( *, '(a)' ) '  * GW_L, any Golub Welsch Linear Growth Family;'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i10))' ) '   DIM: ', point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max
    do dim_num = dim_min, dim_max
      call sparse_grid_onn_l_size ( dim_num, level_max, point_num(dim_num) )
    end do
    write ( *, '(a4,i4,6(2x,i10))' ) '    ', level_max, point_num(dim_min:dim_max)
  end do

  return
end
subroutine test08 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST08 tests SPARSE_GRID_OWN_E_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  SPARSE_GRID_OWN_E_SIZE returns the number of '
  write ( *, '(a)' ) '  distinct points in an OWN_E sparse grid made up of all '
  write ( *, '(a)' ) '  product grids formed from open weakly nested '
  write ( *, '(a)' ) '  quadrature rules with exponential growth, including:'
  write ( *, '(a)' ) '  * GGH_E, the Generalized Gauss-Hermite Exponential Growth Family;'
  write ( *, '(a)' ) '  * GH_E, the Gauss-Hermite Exponential Growth Family;'
  write ( *, '(a)' ) '  * LG_E, the Gauss-Legendre Exponential Growth Family;'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i10))' ) '   DIM: ', point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max
    do dim_num = dim_min, dim_max
      call sparse_grid_own_e_size ( dim_num, level_max, point_num(dim_num) )
    end do
    write ( *, '(a4,i4,6(2x,i10))' ) '    ', level_max, point_num(dim_min:dim_max)
  end do

  return
end
subroutine test09 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST09 tests SPARSE_GRID_OWN_L_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  SPARSE_GRID_OWN_L_SIZE returns the number of '
  write ( *, '(a)' ) '  distinct points in an OWN_L sparse grid made up of all '
  write ( *, '(a)' ) '  product grids formed from open weakly nested '
  write ( *, '(a)' ) '  quadrature rules with linear growth, including:'
  write ( *, '(a)' ) '  * GGH_L, the Generalized Gauss-Hermite Linear Growth Family;'
  write ( *, '(a)' ) '  * GH_L, the Gauss-Hermite Linear Growth Family;'
  write ( *, '(a)' ) '  * LG_L, the Gauss-Legendre Linear Growth Family;'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i10))' ) '   DIM: ', point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max
    do dim_num = dim_min, dim_max
      call sparse_grid_own_l_size ( dim_num, level_max, point_num(dim_num) )
    end do
    write ( *, '(a4,i4,6(2x,i10))' ) '    ', level_max, point_num(dim_min:dim_max)
  end do

  return
end
subroutine test10 ( dim_min, dim_max, level_max_min, level_max_max )

!*****************************************************************************80
!
!! TEST10 tests SPARSE_GRID_OWN_LS_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
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
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  SPARSE_GRID_OWN_LS_SIZE returns the number of '
  write ( *, '(a)' ) '  distinct points in an OWN_LS sparse grid made up of all '
  write ( *, '(a)' ) '  product grids formed from open weakly nested '
  write ( *, '(a)' ) '  quadrature rules with linear (slow) growth, including:'
  write ( *, '(a)' ) '  * GGH_LS, the Generalized Gauss-Hermite Linear (Slow) Growth Family;'
  write ( *, '(a)' ) '  * GH_LS, the Gauss-Hermite Linear (Slow) Growth Family;'
  write ( *, '(a)' ) '  * LG_LS, the Gauss-Legendre Linear (Slow) Growth Family;'
  write ( *, '(a)' ) ' '

  do dim_num = dim_min, dim_max
    point_num(dim_num) = dim_num
  end do

  write ( *, '(a8,6(2x,i10))' ) '   DIM: ', point_num(dim_min:dim_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   LEVEL_MAX'
  write ( *, '(a)' ) ' '

  do level_max = level_max_min, level_max_max
    do dim_num = dim_min, dim_max
      call sparse_grid_own_ls_size ( dim_num, level_max, point_num(dim_num) )
    end do
    write ( *, '(a4,i4,6(2x,i10))' ) '    ', level_max, point_num(dim_min:dim_max)
  end do

  return
end
