program main

!*****************************************************************************80
!
!! MAIN is the main program for SGMGA_VCN_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_VCN_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SGMGA_VCN and SGMGA_VCN_ORDERED functions.'

  call sgmga_vcn_tests ( )

  call sgmga_vcn_timing_tests ( )

  call sgmga_vcn_ordered_tests ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_VCN_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine sgmga_vcn_tests ( )

!*****************************************************************************80
!
!! SGMGA_VCN_TESTS tests SGMGA_VCN.
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

  integer ( kind = 4 ), parameter :: test_num = 12

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) :: dim_num_array(test_num) = (/ &
    2, 2, 2, 2, 2, &
    3, 3, 3, 3, 3, &
    4, 4 /)
  real ( kind = 8 ), allocatable :: importance(:)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) :: level_max_array(test_num) = (/ &
    0, 1, 2, 3, 4, &
    0, 1, 2, 3, 4, &
    2, 3 /)
  real ( kind = 8 ), allocatable :: level_weight(:)
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_VCN_TESTS'
  write ( *, '(a)' ) '  calls SGMGA_VCN_TEST.'
!
!  Isotropic examples.
!
  do test = 1, test_num

    dim_num = dim_num_array(test)
    allocate ( importance(1:dim_num) )
    importance(1:dim_num) = 1.0D+00
    allocate ( level_weight(1:dim_num) )
    call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
    level_max = level_max_array(test)
    q_min = real ( level_max, kind = 8 ) - sum ( level_weight(1:dim_num) )
    q_max = real ( level_max, kind = 8 )

    call sgmga_vcn_test ( dim_num, importance, level_weight, q_min, q_max )

    deallocate ( importance )
    deallocate ( level_weight )

  end do
!
!  Zero weight example.
!
  dim_num = 3
  allocate ( importance(1:dim_num) )
  importance(1:3) = (/ 1.0D+00, 0.0D+00, 1.0D+00 /)
  allocate ( level_weight(1:dim_num) )
  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
  level_max = 2
  q_min = real ( level_max, kind = 8 ) - sum ( level_weight(1:dim_num) )
  q_max = real ( level_max, kind = 8 )

  call sgmga_vcn_test ( dim_num, importance, level_weight, q_min, q_max )

  deallocate ( importance )
  deallocate ( level_weight )
!
!  Anisotropic examples.
!
  do test = 1, test_num

    dim_num = dim_num_array(test)
    allocate ( importance(1:dim_num) )
    do dim = 1, dim_num
      importance(dim) = real ( dim, kind = 8 )
    end do
    allocate ( level_weight(1:dim_num) )
    call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
    level_max = level_max_array(test)
    q_min = real ( level_max, kind = 8 ) - sum ( level_weight(1:dim_num) )
    q_max = real ( level_max, kind = 8 )

    call sgmga_vcn_test ( dim_num, importance, level_weight, q_min, q_max )

    deallocate ( importance )
    deallocate ( level_weight )

  end do

  return
end
subroutine sgmga_vcn_test ( dim_num, importance, level_weight, q_min, q_max )

!*****************************************************************************80
!
!! SGMGA_VCN_TEST tests SGMGA_VCN_NAIVE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) i
  real ( kind = 8 ) importance(dim_num)
  integer ( kind = 4 ) level_1d(dim_num)
  integer ( kind = 4 ) level_1d_max(dim_num)
  integer ( kind = 4 ) level_1d_min(dim_num)
  real ( kind = 8 ) level_weight(dim_num)
  logical              more_grids
  real ( kind = 8 ) q
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  integer ( kind = 4 ) r8_floor
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_VCN_TEST'
  write ( *, '(a)' ) '  Consider vectors 0 <= LEVEL_1D(1:N) <= LEVEL_1D_MAX(1:N),'
  write ( *, '(a)' ) '  Set Q = sum ( LEVEL_1D_WEIGHT(1:N) * LEVEL_1D(1:N) )'
  write ( *, '(a)' ) '  Accept vectors for which Q_MIN < Q <= Q_MAX'
  write ( *, '(a)' ) '  No particular order is imposed on the LEVEL_1D values.'
  write ( *, '(a)' ) '  SGMGA_VCN_NAIVE uses a naive approach;'
  write ( *, '(a)' ) '  SGMGA_VCN tries to be more efficient.'
  write ( *, '(a)' ) '  Here, we just compare the results.'

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
  level_1d_min(1:dim_num) = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  IMPORTANCE:'
  write ( *, '(5g14.6)' ) importance(1:dim_num)
  write ( *, '(a)' ) '  LEVEL_WEIGHT:'
  write ( *, '(5g14.6)' ) level_weight(1:dim_num)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SGMGA_VCN_NAIVE:'
  write ( *, '(a)' ) '     I       Q           X'
  write ( *, '(2x,a4,2x,g14.6,10(2x,i2))' ) &
    ' MIN', q_min, level_1d_min(1:dim_num)

  i = 0
  more_grids = .false.

  do

    call sgmga_vcn_naive ( dim_num, level_weight, level_1d_max, level_1d, &
      q_min, q_max, more_grids )

    if ( .not. more_grids ) then
      exit
    end if

    q = dot_product ( level_weight(1:dim_num), &
      real ( level_1d(1:dim_num), kind = 8 ) )
    i = i + 1
    write ( *, '(2x,i4,2x,g14.6,10(2x,i2))' ) i, q, level_1d(1:dim_num)

  end do

  write ( *, '(2x,a4,2x,g14.6,10(2x,i2))' ) &
    ' MAX', q_max, level_1d_max(1:dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SGMGA_VCN:'
  write ( *, '(a)' ) '     I       Q           X'
  write ( *, '(2x,a4,2x,g14.6,10(2x,i2))' ) &
    ' MIN', q_min, level_1d_min(1:dim_num)

  i = 0
  more_grids = .false.

  do

    call sgmga_vcn ( dim_num, level_weight, level_1d, q_min, q_max, more_grids )

    if ( .not. more_grids ) then
      exit
    end if

    q = dot_product ( level_weight(1:dim_num), &
      real ( level_1d(1:dim_num), kind = 8 ) )
    i = i + 1
    write ( *, '(2x,i4,2x,g14.6,10(2x,i2))' ) i, q, level_1d(1:dim_num)

  end do

  write ( *, '(2x,a4,2x,g14.6,10(2x,i2))' ) &
    ' MAX', q_max, level_1d_max(1:dim_num)

  return
end
subroutine sgmga_vcn_timing_tests ( )

!*****************************************************************************80
!
!! SGMGA_VCN_TIMING_TESTS times SGMGA_VCN and SGMGA_VCN_NAIVE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ), allocatable :: importance(:)
  integer ( kind = 4 ) level_max
  real ( kind = 8 ), allocatable :: level_weight(:)
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_VCN_TIMING_TESTS'
  write ( *, '(a)' ) '  calls SGMGA_VCN_TIMING_TEST.'
!
!  Isotropic examples.
!
  dim_num = 2

  do test = 1, 2

    dim_num = dim_num * 2
    allocate ( importance(1:dim_num) )
    importance(1:dim_num) = 1.0D+00
    allocate ( level_weight(1:dim_num) )
    call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
    level_max = 2
    q_min = real ( level_max, kind = 8 ) - sum ( level_weight(1:dim_num) )
    q_max = real ( level_max, kind = 8 )

    call sgmga_vcn_timing_test ( dim_num, importance, level_weight, q_min, &
      q_max )

    deallocate ( importance )
    deallocate ( level_weight )

  end do
!
!  Anisotropic examples.
!
  dim_num = 2

  do test = 1, 2

    dim_num = dim_num * 2
    allocate ( importance(1:dim_num) )
    do dim = 1, dim_num
      importance(dim) = real ( dim, kind = 8 )
    end do
    allocate ( level_weight(1:dim_num) )
    call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
    level_max = 2
    q_min = real ( level_max, kind = 8 ) - sum ( level_weight(1:dim_num) )
    q_max = real ( level_max, kind = 8 )

    call sgmga_vcn_timing_test ( dim_num, importance, level_weight, q_min, &
      q_max )

    deallocate ( importance )
    deallocate ( level_weight )

  end do

  return
end
subroutine sgmga_vcn_timing_test ( dim_num, importance, level_weight, q_min, &
  q_max )

!*****************************************************************************80
!
!! SGMGA_VCN_TIMING_TEST times SGMGA_VCN_NAIVE and SGMGA_VCN.
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

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) i
  real ( kind = 8 ) importance(dim_num)
  integer ( kind = 4 ) level_1d(dim_num)
  integer ( kind = 4 ) level_1d_max(dim_num)
  integer ( kind = 4 ) level_1d_min(dim_num)
  real ( kind = 8 ) level_weight(dim_num)
  logical              more_grids
  real ( kind = 8 ) q
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  integer ( kind = 4 ) r8_floor
  integer ( kind = 4 ) test
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_VCN_TEST'
  write ( *, '(a)' ) '  Consider vectors 0 <= LEVEL_1D(1:N) <= LEVEL_1D_MAX(1:N),'
  write ( *, '(a)' ) '  Set Q = sum ( LEVEL_1D_WEIGHT(1:N) * LEVEL_1D(1:N) )'
  write ( *, '(a)' ) '  Accept vectors for which Q_MIN < Q <= Q_MAX'
  write ( *, '(a)' ) '  No particular order is imposed on the LEVEL_1D values.'
  write ( *, '(a)' ) '  SGMGA_VCN_NAIVE uses a naive approach;'
  write ( *, '(a)' ) '  SGMGA_VCN tries to be more efficient.'
  write ( *, '(a)' ) '  Here, we compare the timings.'

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
  level_1d_min(1:dim_num) = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  IMPORTANCE:'
  write ( *, '(5g14.6)' ) importance(1:dim_num)
  write ( *, '(a)' ) '  LEVEL_WEIGHT:'
  write ( *, '(5g14.6)' ) level_weight(1:dim_num)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SGMGA_VCN_NAIVE:'
  write ( *, '(a)' ) '     I       Q           X'
  write ( *, '(2x,a4,2x,g14.6,10(2x,i2))' ) &
    ' MIN', q_min, level_1d_min(1:dim_num)

  i = 0
  more_grids = .false.

  call cpu_time ( t1 )
  do

    call sgmga_vcn_naive ( dim_num, level_weight, level_1d_max, level_1d, &
      q_min, q_max, more_grids )

    if ( .not. more_grids ) then
      exit
    end if

!   q = dot_product ( level_weight(1:dim_num), &
!     real ( level_1d(1:dim_num), kind = 8 ) )
!   i = i + 1
!   write ( *, '(2x,i4,2x,g14.6,10(2x,i2))' ) i, q, level_1d(1:dim_num)

  end do

  call cpu_time ( t2 )

  write ( *, '(2x,a4,2x,g14.6,10(2x,i2))' ) &
    ' MAX', q_max, level_1d_max(1:dim_num)
  write ( *, '(2x,a4,2x,g14.6)' ) 'TIME', t2 - t1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SGMGA_VCN:'
  write ( *, '(a)' ) '     I       Q           X'
  write ( *, '(2x,a4,2x,g14.6,10(2x,i2))' ) &
    ' MIN', q_min, level_1d_min(1:dim_num)

  i = 0
  more_grids = .false.

  call cpu_time ( t1 )

  do

    call sgmga_vcn ( dim_num, level_weight, level_1d, q_min, q_max, more_grids )

    if ( .not. more_grids ) then
      exit
    end if

!   q = dot_product ( level_weight(1:dim_num), &
!     real ( level_1d(1:dim_num), kind = 8 ) )
!   i = i + 1
!   write ( *, '(2x,i4,2x,g14.6,10(2x,i2))' ) i, q, level_1d(1:dim_num)

  end do

  call cpu_time ( t2 )

  write ( *, '(2x,a4,2x,g14.6,10(2x,i2))' ) &
    ' MAX', q_max, level_1d_max(1:dim_num)
  write ( *, '(2x,a4,2x,g14.6)' ) 'TIME', t2 - t1


  return
end
subroutine sgmga_vcn_ordered_tests ( )

!*****************************************************************************80
!
!! SGMGA_VCN_TESTS tests SGMGA_VCN_ORDERED_NAIVE.
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

  integer ( kind = 4 ), parameter :: test_num = 12

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) :: dim_num_array(test_num) = (/ &
    2, 2, 2, 2, 2, &
    3, 3, 3, 3, 3, &
    4, 4 /)
  real ( kind = 8 ), allocatable :: importance(:)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) :: level_max_array(test_num) = (/ &
    0, 1, 2, 3, 4, &
    0, 1, 2, 3, 4, &
    2, 3 /)
  real ( kind = 8 ), allocatable :: level_weight(:)
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_VCN_ORDERED_TESTS'
  write ( *, '(a)' ) '  calls SGMGA_VCN_ORDERED_TEST.'
!
!  Isotropic examples.
!
  do test = 1, test_num

    dim_num = dim_num_array(test)
    allocate ( importance(1:dim_num) )
    importance(1:dim_num) = 1.0D+00
    allocate ( level_weight(1:dim_num) )
    call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
    level_max = level_max_array(test)
    q_min = real ( level_max, kind = 8 ) - sum ( level_weight(1:dim_num) )
    q_max = real ( level_max, kind = 8 )

    call sgmga_vcn_ordered_test ( dim_num, importance, level_weight, q_min, &
      q_max )

    deallocate ( importance )
    deallocate ( level_weight )

  end do
!
!  Anisotropic examples.
!
  do test = 1, test_num

    dim_num = dim_num_array(test)
    allocate ( importance(1:dim_num) )
    do dim = 1, dim_num
      importance(dim) = real ( dim, kind = 8 )
    end do
    allocate ( level_weight(1:dim_num) )
    call sgmga_importance_to_aniso ( dim_num, importance, level_weight )
    level_max = level_max_array(test)
    q_min = real ( level_max, kind = 8 ) - sum ( level_weight(1:dim_num) )
    q_max = real ( level_max, kind = 8 )

    call sgmga_vcn_ordered_test ( dim_num, importance, level_weight, q_min, &
      q_max )

    deallocate ( importance )
    deallocate ( level_weight )

  end do

  return
end
subroutine sgmga_vcn_ordered_test ( dim_num, importance, level_weight, q_min, &
  q_max )

!*****************************************************************************80
!
!! SGMGA_VCN_ORDERED_TEST tests SGMGA_VCN_ORDERED_NAIVE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) i
  real ( kind = 8 ) importance(dim_num)
  integer ( kind = 4 ) level_1d(dim_num)
  integer ( kind = 4 ) level_1d_max(dim_num)
  integer ( kind = 4 ) level_1d_min(dim_num)
  real ( kind = 8 ) level_weight(dim_num)
  logical more_grids
  real ( kind = 8 ) q
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  integer ( kind = 4 ) r8_floor
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_VCN_ORDERED_TEST'
  write ( *, '(a)' ) '  Consider vectors 0 <= LEVEL_1D(1:N) <= LEVEL_1D_MAX(1:N),'
  write ( *, '(a)' ) '  Set Q = sum ( LEVEL_WEIGHT(1:N) * LEVEL_1D(1:N) )'
  write ( *, '(a)' ) '  Accept only vectors for which Q_MIN < Q <= Q_MAX'
  write ( *, '(a)' ) '  The solutions are weakly ordered by the value of Q.'
  write ( *, '(a)' ) '  SGMGA_VCN_ORDERED_NAIVE calls SGMGA_VCN_NAIVE;'
  write ( *, '(a)' ) '  SGMGA_VCN_ORDERED calls SGMGA_VCN.'

  do dim = 1, dim_num
    if ( 0.0D+00 < level_weight(dim) ) then
      level_1d_max(dim) = r8_floor ( q_max / level_weight(dim) ) + 1
    else
      level_1d_max(dim) = 0
    end if
  end do

  level_1d_min(1:dim_num) = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  IMPORTANCE:'
  write ( *, '(5g14.6)' ) importance(1:dim_num)
  write ( *, '(a)' ) '  LEVEL_WEIGHT:'
  write ( *, '(5g14.6)' ) level_weight(1:dim_num)
!-------------------------------------------------
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SGMGA_VCN_ORDERED_NAIVE:'
  write ( *, '(a)' ) '     I       Q           X'
  write ( *, '(2x,a4,2x,g14.6,10(2x,i2))' ) &
    ' MIN', q_min, level_1d_min(1:dim_num)

  i = 0
  more_grids = .false.

  do

    call sgmga_vcn_ordered_naive ( dim_num, level_weight, level_1d_max, &
      level_1d, q_min, q_max, more_grids )

    if ( .not. more_grids ) then
      exit
    end if

    q = dot_product ( level_weight(1:dim_num), &
      real ( level_1d(1:dim_num), kind = 8 ) )
    i = i + 1
    write ( *, '(2x,i4,2x,g14.6,10(2x,i2))' ) i, q, level_1d(1:dim_num)

  end do

  write ( *, '(2x,a4,2x,g14.6,10(2x,i2))' ) &
    ' MAX', q_max, level_1d_max(1:dim_num)
!-------------------------------------------------
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SGMGA_VCN_ORDERED:'
  write ( *, '(a)' ) '     I       Q           X'
  write ( *, '(2x,a4,2x,g14.6,10(2x,i2))' ) &
    ' MIN', q_min, level_1d_min(1:dim_num)

  i = 0
  more_grids = .false.

  do

    call sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, &
      level_1d, q_min, q_max, more_grids )

    if ( .not. more_grids ) then
      exit
    end if

    q = dot_product ( level_weight(1:dim_num), &
      real ( level_1d(1:dim_num), kind = 8 ) )
    i = i + 1
    write ( *, '(2x,i4,2x,g14.6,10(2x,i2))' ) i, q, level_1d(1:dim_num)

  end do

  write ( *, '(2x,a4,2x,g14.6,10(2x,i2))' ) &
    ' MAX', q_max, level_1d_max(1:dim_num)

  return
end
