program main

!*****************************************************************************80
!
!! MAIN is the main program for SGMGA_VCN_COEF_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_VCN_COEF_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SGMA_VCN_COEF function.'

  call sgmga_vcn_coef_tests ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_VCN_COEF_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine sgmga_vcn_coef_tests ( )

!*****************************************************************************80
!
!! SGMGA_VCN_COEF_TESTS tests SGMGA_VCN_COEF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  real    ( kind = 8 ), allocatable :: importance(:)
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  real    ( kind = 8 ), allocatable :: level_weight(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_VCN_COEF_TESTS'
  write ( *, '(a)' ) '  calls SGMGA_VCN_COEF_TEST.'

  dim_num = 2
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = 1.0D+00
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max_min = 0
  level_max_max = 4
  call sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max )
  deallocate ( importance )
  deallocate ( level_weight )

  dim_num = 2
  allocate ( importance(1:dim_num) )
  importance(1:dim) = (/ 2.0D+00, 1.0D+00 /)
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max_min = 0
  level_max_max = 8
  call sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max )
  deallocate ( importance )
  deallocate ( level_weight )

  dim_num = 3
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = 1.0D+00
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max_min = 0
  level_max_max = 4
  call sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max )
  deallocate ( importance )
  deallocate ( level_weight )

  dim_num = 3
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = real ( dim, kind = 8 )
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max_min = 0
  level_max_max = 4
  call sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max )
  deallocate ( importance )
  deallocate ( level_weight )

  dim_num = 4
  allocate ( importance(1:dim_num) )
  do dim = 1, dim_num
    importance(dim) = real ( dim, kind = 8 )
  end do
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max_min = 0
  level_max_max = 3
  call sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max )
  deallocate ( importance )
  deallocate ( level_weight )
!
!  Try a case which includes a dimension of "0 importance".
!
  dim_num = 3
  allocate ( importance(1:dim_num) )
  importance = (/ 1.0D+00, 0.0D+00, 1.0D+00 /)
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max_min = 0
  level_max_max = 3
  call sgmga_vcn_coef_test ( dim_num, importance, level_weight, level_max_min, &
    level_max_max )
  deallocate ( importance )
  deallocate ( level_weight )

  return
end
subroutine sgmga_vcn_coef_test ( dim_num, importance, level_weight, &
  level_max_min, level_max_max )

!*****************************************************************************80
!
!! SGMGA_VCN_COEF_TEST tests SGMGA_VCN_COEF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 April 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num

  real    ( kind = 8 ) coef1
  real    ( kind = 8 ) coef1_sum
  real    ( kind = 8 ) coef2
  real    ( kind = 8 ) coef2_sum
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) i
  real    ( kind = 8 ) importance(dim_num)
  integer ( kind = 4 ) level_1d(dim_num)
  integer ( kind = 4 ) level_1d_max(dim_num)
  integer ( kind = 4 ) level_1d_min(dim_num)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  real    ( kind = 8 ) level_weight(dim_num)
  real    ( kind = 8 ) level_weight_min_pos
  logical more_grids
  real    ( kind = 8 ) q
  real    ( kind = 8 ) q_max
  real    ( kind = 8 ) q_min
  integer ( kind = 4 ) r8_floor
  real    ( kind = 8 ) r8vec_min_pos

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_VCN_COEF_TEST'
  write ( *, '(a)' ) '  For anisotropic problems, each product grid in a sparse'
  write ( *, '(a)' ) '  grid has an associated "combinatorial coefficient".'
  write ( *, '(a)' ) '  SGMGA_VCN_COEF_NAIVE uses a naive algorithm.'
  write ( *, '(a)' ) '  SGMGA_VCN_COEF attempts a more efficient method.'
  write ( *, '(a)' ) '  Here, we simply compare COEF1 and COEF2, the same'
  write ( *, '(a)' ) '  coefficient computed by the naive and efficient ways.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  IMPORTANCE:  '
  write ( *, '(5g14.6)' ) importance(1:dim_num)
  write ( *, '(a)' ) '  LEVEL_WEIGHT:'
  write ( *, '(5g14.6)' ) level_weight(1:dim_num)

  do level_max = level_max_min, level_max_max

    i = 0
    coef1_sum = 0.0D+00
    coef2_sum = 0.0D+00
!
!  Initialization.
!
    level_1d_min(1:dim_num) = 0

    level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight )
    q_min = real ( level_max, kind = 8 ) * level_weight_min_pos &
      - sum ( level_weight(1:dim_num) )
    q_max = real ( level_max, kind = 8 ) * level_weight_min_pos
    do dim = 1, dim_num
      if ( 0.0D+00 < level_weight(dim) ) then
        level_1d_max(dim) = r8_floor ( q_max / level_weight(dim) ) + 1
        if ( q_max <= ( level_1d_max(dim) - 1 ) * level_weight(dim) ) then
          level_1d_max(dim) = level_1d_max(dim) - 1
        end if
      else
        level_1d_max(dim) = 0
      end if
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     I           Q       Coef1       Coef2       X'
    write ( *, '(2x,a4,2x,g14.6,24x,10(2x,i2))' ) &
      ' MIN', q_min, level_1d_min(1:dim_num)
!
!  Seek all vectors LEVEL_1D which satisfy the constraint:
!
!    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT )
!      < sum ( 0 <= I < DIM_NUM ) LEVEL_WEIGHT(I) * LEVEL_1D(I)
!      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
!
    more_grids = .false.

    do

      call sgmga_vcn_ordered_naive ( dim_num, level_weight, level_1d_max, &
        level_1d, q_min, q_max, more_grids )

      if ( .not. more_grids ) then
        exit
      end if
!
!  Compute the combinatorial coefficient.
!
      call sgmga_vcn_coef_naive ( dim_num, level_weight, level_1d_max, level_1d, &
        q_min, q_max, coef1 )

      call sgmga_vcn_coef ( dim_num, level_weight,level_1d, q_max, coef2 )

      i = i + 1

      q = dot_product ( level_weight(1:dim_num), &
        real ( level_1d(1:dim_num), kind = 8 ) )

      coef1_sum = coef1_sum + coef1
      coef2_sum = coef2_sum + coef2

      write ( *, '(2x,i4,2x,g14.6,2x,g10.2,2x,g10.2,10(2x,i2))' ) &
        i, q, coef1, coef2, level_1d(1:dim_num)

    end do

    write ( *, '(2x,a4,2x,g14.6,24x,10(2x,i2))' ) &
      ' MAX', q_max, level_1d_max(1:dim_num)
    write ( *, '(a,g10.2,2x,g10.2)' ) &
      '   SUM                  ', coef1_sum, coef2_sum

  end do

  return
end

