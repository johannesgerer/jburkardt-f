program main

!*****************************************************************************80
!
!! MAIN is the main program for SGMGA_ANISO_NORMALIZE_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_ANISO_NORMALIZE_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SGMGA_ANISO_NORMALIZE and'
  write ( *, '(a)' ) '  SGMGA_IMPORTANCE_TO_ANISO functions.'

  call sgmga_aniso_balance_tests ( )

  call sgmga_aniso_normalize_tests ( )

  call sgmga_importance_to_aniso_tests ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_ANISO_NORMALIZE_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine sgmga_aniso_balance_tests ( )

!*****************************************************************************80
!
!! SGMGA_ANISO_BALANCE_TESTS call SGMGA_ANISO_BALANCE_TEST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) alpha_max
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ), allocatable :: level_weight(:)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_ANISO_BALANCE_TESTS'
  write ( *, '(a)' ) '  Call SGMGA_ANISO_BALANCE_TEST with various arguments.'

  alpha_max = 10.0D+00
  seed = 123456789
  dim_num = 5
  allocate ( level_weight(1:dim_num) )

  do test = 1, test_num
    call r8vec_uniform_01 ( dim_num, seed, level_weight )
    level_weight(1:dim_num) = 10.0D+00 * level_weight(1:dim_num)
    call sgmga_aniso_balance_test ( alpha_max, dim_num, level_weight )
  end do

  alpha_max = 5.0D+00
  seed = 123456789
  dim_num = 5

  do test = 1, test_num
    call r8vec_uniform_01 ( dim_num, seed, level_weight )
    level_weight(1:dim_num) = 10.0D+00 * level_weight(1:dim_num)
    call sgmga_aniso_balance_test ( alpha_max, dim_num, level_weight )
  end do

  deallocate ( level_weight )

  return
end
subroutine sgmga_aniso_balance_test ( alpha_max, dim_num, level_weight )

!*****************************************************************************80
!
!! SGMGA_ANISO_BALANCE_TEST calls SGMGA_ANISO_BALANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) alpha_max
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) option
  real ( kind = 8 ) level_weight(dim_num)
  real ( kind = 8 ) level_weight2(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_ANISO_BALANCE_TEST'
  write ( *, '(a,g14.6)' ) '  ALPHA_MAX = ', alpha_max
  write ( *, '(a,g14.6)' ) &
    '  Input weight sum: ', sum ( level_weight(1:dim_num) )
  write ( *, '(5g14.6)' ) level_weight(1:dim_num)

  call sgmga_aniso_balance ( alpha_max, dim_num, level_weight, level_weight2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) &
    '  Output weight sum:', sum ( level_weight2(1:dim_num) )
  write ( *, '(5g14.6)' ) level_weight2(1:dim_num)

  return
end
subroutine sgmga_aniso_normalize_tests ( )

!*****************************************************************************80
!
!! SGMGA_ANISO_NORMALIZE_TESTS call SGMGA_ANISO_NORMALIZE_TEST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num
  real ( kind = 8 ), allocatable :: level_weight(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_ANISO_NORMALIZE_TESTS'
  write ( *, '(a)' ) '  Call SGMGA_ANISO_NORMALIZE_TEST with various arguments.'

  dim_num = 2
  allocate ( level_weight(1:dim_num) )
  level_weight(1:dim_num) = (/ 1.0D+00, 1.0D+00 /)
  call sgmga_aniso_normalize_test ( dim_num, level_weight )
  deallocate ( level_weight )

  dim_num = 2
  allocate ( level_weight(1:dim_num) )
  level_weight(1:dim_num) = (/ 10.0D+00, 10.0D+00 /)
  call sgmga_aniso_normalize_test ( dim_num, level_weight )
  deallocate ( level_weight )

  dim_num = 2
  allocate ( level_weight(1:dim_num) )
  level_weight(1:dim_num) = (/ 10.0D+00, 2.0D+00 /)
  call sgmga_aniso_normalize_test ( dim_num, level_weight )
  deallocate ( level_weight )

  dim_num = 2
  allocate ( level_weight(1:dim_num) )
  level_weight(1:dim_num) = (/ 1.0D+00, 2.0D+00 /)
  call sgmga_aniso_normalize_test ( dim_num, level_weight )
  deallocate ( level_weight )

  dim_num = 3
  allocate ( level_weight(1:dim_num) )
  level_weight(1:dim_num) = (/ 1.0D+00, 2.0D+00, 3.0D+00 /)
  call sgmga_aniso_normalize_test ( dim_num, level_weight )
  deallocate ( level_weight )
!
!  Try a case in which one variable has 0 weight.
!
  dim_num = 3
  allocate ( level_weight(1:dim_num) )
  level_weight(1:dim_num) = (/ 2.0D+00, 0.0D+00, 1.5D+00 /)
  call sgmga_aniso_normalize_test ( dim_num, level_weight )
  deallocate ( level_weight )

  dim_num = 4
  allocate ( level_weight(1:dim_num) )
  level_weight(1:dim_num) = (/ 1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00 /)
  call sgmga_aniso_normalize_test ( dim_num, level_weight )
  deallocate ( level_weight )

  return
end
subroutine sgmga_aniso_normalize_test ( dim_num, level_weight )

!*****************************************************************************80
!
!! SGMGA_ANISO_NORMALIZE_TEST calls SGMGA_ANISO_NORMALIZE.
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

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) option
  real ( kind = 8 ) level_weight(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_ANISO_NORMALIZE_TEST'
  write ( *, '(a,g14.6)' ) &
    '  Input weight sum: ', sum ( level_weight(1:dim_num) )
  write ( *, '(5g14.6)' ) level_weight(1:dim_num)

  do option = 0, 2

    call sgmga_aniso_normalize ( option, dim_num, level_weight )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i4,a,g14.6)' ) &
    '  For OPTION = ', option, &
    '  Normalized weight sum:', sum ( level_weight(1:dim_num) )
    write ( *, '(5g14.6)' ) level_weight(1:dim_num)

  end do

  return
end
subroutine sgmga_importance_to_aniso_tests ( )

!*****************************************************************************80
!
!! SGMGA_IMPORTANCE_TO_ANISO_TESTS call SGMGA_IMPORTANCE_TO_ANISO_TEST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num
  real ( kind = 8 ), allocatable :: importance(:)
  real ( kind = 8 ), allocatable :: level_weight(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_IMPORTANCE_TO_ANISO_TESTS'
  write ( *, '(a)' ) '  Call SGMGA_IMPORTANCE_TO_ANISO_TEST with various arguments.'

  dim_num = 2
  allocate ( importance(1:dim_num) )
  allocate ( level_weight(1:dim_num) )
  importance(1:dim_num) = (/ 1.0D+00, 1.0D+00 /)
  call sgmga_importance_to_aniso_test ( dim_num, importance, level_weight )
  deallocate ( importance )
  deallocate ( level_weight )

  dim_num = 2
  allocate ( importance(1:dim_num) )
  allocate ( level_weight(1:dim_num) )
  importance(1:dim_num) = (/ 10.0D+00, 10.0D+00 /)
  call sgmga_importance_to_aniso_test ( dim_num, importance, level_weight )
  deallocate ( importance )
  deallocate ( level_weight )

  dim_num = 2
  allocate ( importance(1:dim_num) )
  allocate ( level_weight(1:dim_num) )
  importance(1:dim_num) = (/ 10.0D+00, 2.0D+00 /)
  call sgmga_importance_to_aniso_test ( dim_num, importance, level_weight )
  deallocate ( importance )
  deallocate ( level_weight )

  dim_num = 2
  allocate ( importance(1:dim_num) )
  allocate ( level_weight(1:dim_num) )
  importance(1:dim_num) = (/ 1.0D+00, 2.0D+00 /)
  call sgmga_importance_to_aniso_test ( dim_num, importance, level_weight )
  deallocate ( importance )
  deallocate ( level_weight )

  dim_num = 3
  allocate ( importance(1:dim_num) )
  allocate ( level_weight(1:dim_num) )
  importance(1:dim_num) = (/ 1.0D+00, 2.0D+00, 3.0D+00 /)
  call sgmga_importance_to_aniso_test ( dim_num, importance, level_weight )
  deallocate ( importance )
  deallocate ( level_weight )
!
!  Try a case in which one variable has 0 importance.
!
  dim_num = 3
  allocate ( importance(1:dim_num) )
  allocate ( level_weight(1:dim_num) )
  importance(1:dim_num) = (/ 2.0D+00, 0.0D+00, 1.5D+00 /)
  call sgmga_importance_to_aniso_test ( dim_num, importance, level_weight )
  deallocate ( importance )
  deallocate ( level_weight )

  dim_num = 4
  allocate ( importance(1:dim_num) )
  allocate ( level_weight(1:dim_num) )
  importance(1:dim_num) = (/ 1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00 /)
  call sgmga_importance_to_aniso_test ( dim_num, importance, level_weight )
  deallocate ( importance )
  deallocate ( level_weight )

  return
end
subroutine sgmga_importance_to_aniso_test ( dim_num, importance, level_weight )

!*****************************************************************************80
!
!! SGMGA_IMPORTANCE_TO_ANISO_TEST calls SGMGA_IMPORTANCE_TO_ANISO.
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

  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) importance(dim_num)
  real ( kind = 8 ) level_weight(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_IMPORTANCE_TO_ANISO_TEST'
  write ( *, '(a)' ) '  Importances:'
  write ( *, '(5g14.6)' ) importance(1:dim_num)

  call sgmga_importance_to_aniso ( dim_num, importance, level_weight )

  write ( *, '(a)' ) '  Anisotropic coefficients:'
  write ( *, '(5g14.6)' ) level_weight(1:dim_num)

  return
end
