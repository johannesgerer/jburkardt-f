program main

!*****************************************************************************80
!
!! TASK_DICISION_TEST tests TASK_DIVISION.
!
!  Discussion:
!
!    This program simply demonstrates how one might automate the
!    assignment of T tasks to P processors, assuming that the assignment
!    is to be beforehand.
!
!    In that case, we just want to make sure that we assign each task
!    to a processor, that we assign about the same number of tasks
!    to each processor, and that we assign each processor a contiguous
!    range of tasks, say tasks I_LO to I_HI.
!
!    The routine that is called simulates this process.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 October 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) proc_first
  integer ( kind = 4 ) proc_last
  integer ( kind = 4 ) task_number

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TASK_DIVISION_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate how to automate the division of'
  write ( *, '(a)' ) '  T tasks among a range of P processors'
  write ( *, '(a)' ) '  indexed from PROC_FIRST to PROC_LAST.'

  task_number = 23
  proc_first = 0
  proc_last = 3
  call task_division ( task_number, proc_first, proc_last )

  task_number = 17
  proc_first = 1
  proc_last = 6
  call task_division ( task_number, proc_first, proc_last )

  task_number = 17
  proc_first = 4
  proc_last = 6
  call task_division ( task_number, proc_first, proc_last )

  task_number = 5
  proc_first = -2
  proc_last = 6
  call task_division ( task_number, proc_first, proc_last )

  task_number = 5
  proc_first = 0
  proc_last = 4
  call task_division ( task_number, proc_first, proc_last )

  task_number = 5
  proc_first = 0
  proc_last = 0
  call task_division ( task_number, proc_first, proc_last )

  task_number = 1000
  proc_first = 1
  proc_last = 17
  call task_division ( task_number, proc_first, proc_last )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TASK_DIVISION_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
