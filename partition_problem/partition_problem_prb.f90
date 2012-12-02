program main

!*****************************************************************************80
!
!! MAIN is the main program for PARTITION_PROBLEM_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5
  integer ( kind = 4 ), allocatable :: w(:)

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PARTITION_PROBLEM_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PARTITION_PROBLEM library.'
!
!  Find individual solutions.
!
  do test = 1, test_num

    if ( test == 1 ) then
      n = 5
      allocate ( w(1:n) )
      w = (/ 19, 17, 13, 9, 6 /)
    else if ( test == 2 ) then
      n = 9
      allocate ( w(1:n) )
      w = (/ 484, 114, 205, 288, 506, 503, 201, 127, 410 /)
    else if ( test == 3 ) then
      n = 10
      allocate ( w(1:n) )
      w = (/ 771, 121, 281, 854, 885, 734, 486, 1003, 83, 62 /)
    else if ( test == 4 ) then
      n = 10
      allocate ( w(1:n) )
      w = (/ 2, 10, 3, 8, 5, 7, 9, 5, 3, 2 /)
    else if ( test == 5 ) then
      n = 9
      allocate ( w(1:n) )
      w = (/ 3, 4, 3, 1, 3, 2, 3, 2, 1 /)
    end if

    call test01 ( n, w )

    deallocate ( w )

  end do
!
!  Count solutions.
!
  do test = 1, test_num

    if ( test == 1 ) then
      n = 5
      allocate ( w(1:n) )
      w = (/ 19, 17, 13, 9, 6 /)
    else if ( test == 2 ) then
      n = 9
      allocate ( w(1:n) )
      w = (/ 484, 114, 205, 288, 506, 503, 201, 127, 410 /)
    else if ( test == 3 ) then
      n = 10
      allocate ( w(1:n) )
      w = (/ 771, 121, 281, 854, 885, 734, 486, 1003, 83, 62 /)
    else if ( test == 4 ) then
      n = 10
      allocate ( w(1:n) )
      w = (/ 2, 10, 3, 8, 5, 7, 9, 5, 3, 2 /)
    else if ( test == 5 ) then
      n = 9
      allocate ( w(1:n) )
      w = (/ 3, 4, 3, 1, 3, 2, 3, 2, 1 /)
    end if

    call test02 ( n, w )

    deallocate ( w )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PARTITION_PROBLEM_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine test01 ( n, w )

!*****************************************************************************80
!
!! TEST01 tests PARTITION_BRUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of weights.
!
!    Input, integer ( kind = 4 ) W(N), a set of weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) discrepancy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) w(n)
  integer ( kind = 4 ) w0_sum
  integer ( kind = 4 ) w1_sum

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Partition a set of N integers W so that the subsets'
  write ( *, '(a)' ) '  have equal sums.'

  call partition_brute ( n, w, c, discrepancy )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I        W0        W1'
  write ( *, '(a)' ) ' '
  w0_sum = 0
  w1_sum = 0
  do i = 1, n
    if ( c(i) == 0 ) then
      w0_sum = w0_sum + w(i)
      write ( *, '(2x,i4,2x,i8,2x,8x)' ) i, w(i)
    else
      w1_sum = w1_sum + w(i)
      write ( *, '(2x,i4,2x,8x,2x,i8)' ) i, w(i)
    end if
  end do
  write ( *, '(a)' ) '        --------  --------'
  write ( *, '(2x,4x,2x,i8,2x,i8)' ) w0_sum, w1_sum
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Discrepancy = ', discrepancy

  return
end
subroutine test02 ( n, w )

!*****************************************************************************80
!
!! TEST02 tests PARTITION_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of weights.
!
!    Input, integer ( kind = 4 ) W(N), a set of weights.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) count
  integer ( kind = 4 ) i
  integer ( kind = 4 ) w(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  PARTITION_COUNT counts the number of exact solutions'
  write ( *, '(a)' ) '  of the partition problem.'

  call partition_count ( n, w, count )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I        W'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i4,2x,i8)' ) i, w(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(2x,a,i4)' ) 'Number of solutions = ', count

  return
end
