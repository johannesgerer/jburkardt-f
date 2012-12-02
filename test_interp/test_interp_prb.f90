program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_INTERP_PRB.
!
!  Discussion:
!
!    TEST_INTERP_PRB calls the TEST_INTERP tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp (  )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INTERP_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_INTERP library.'
  write ( *, '(a)' ) '  This test also requires the R8LIB library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_INTERP_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 shows how P00_STORY can be called.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  P00_STORY prints the problem "story".'

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', prob

    call p00_story ( prob )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 prints the data for each problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 May 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ), allocatable :: p(:,:)
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  P00_DATA_NUM returns N, the number of data points.'
  write ( *, '(a)' ) '  P00_DIM_NUM returns M, the dimension of data.'
  write ( *, '(a)' ) '  P00_DATA returns the actual (MxN) data.'

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem  ', prob

    call p00_data_num ( prob, data_num )
    write ( *, '(a,i8)' ) '  DATA_NUM ', data_num
    call p00_dim_num ( prob, dim_num )
    write ( *, '(a,i8)' ) '  DIM_NUM  ', dim_num
    allocate ( p(1:dim_num,1:data_num) )

    call p00_data ( prob, dim_num, data_num, p )

    call r8mat_transpose_print ( dim_num, data_num, p, '  Data array:' )

    deallocate ( p )

  end do

  return
end
