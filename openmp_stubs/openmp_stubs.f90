subroutine omp_destroy_lock ( lock )

!*****************************************************************************80
!
!! OMP_DESTROY_LOCK destroys a simple lock.
!
!  Discussion:
!
!    The routine is intended to return the state of the lock to the 
!    uninitialized state.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, integer ( kind = omp_lock_kind ) LOCK, the simple lock.
!
  implicit none

  include 'omp_lib_kinds.h'

  integer ( kind = omp_lock_kind ) lock

  lock = 0

  return 
end subroutine
subroutine omp_destroy_nest_lock ( nlock )

!*****************************************************************************80
!
!! OMP_DESTROY_NEST_LOCK destroys a nestable lock.
!
!  Discussion:
!
!    The routine is intended to return the state of the lock to the 
!    uninitialized state.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July, 2011.
!
!  Parameters:
!
!    Output, integer ( kind = omp_nest_lock_kind ) NLOCK, the nestable lock.
!
  implicit none

  include 'omp_lib_kinds.h'

  integer ( kind = omp_nest_lock_kind ) nlock

  nlock = 0

  return 
end subroutine
function omp_get_active_level ( )

!*****************************************************************************80
!
!! OMP_GET_ACTIVE_LEVEL returns the number of nested active parallel regions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, integer OMP_GET_ACTIVE_LEVEL, the number of nested active parallel 
!    regions enclosing the task that contains this call.
!
  implicit none

  integer omp_get_active_level

  omp_get_active_level = 0

  return
end function
function omp_get_ancestor_thread_num ( level )

!*****************************************************************************80
!
!! OMP_GET_ANCESTOR_THREAD_NUM returns the thread number of the ancestor of this thread.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Input, integer LEVEL, the nested level.
!
!    Output, integer OMP_GET_ANCESTOR_THREAD_NUM, the thread number of the 
!    ancestor of this thread.
!
  implicit none

  integer level
  integer omp_get_ancestor_thread_num

  if ( level == 0 ) then
    omp_get_ancestor_thread_num = 0
  else
    omp_get_ancestor_thread_num = -1
  end if

  return 
end function
function omp_get_dynamic ( )

!*****************************************************************************80
!
!! OMP_GET_DYNAMIC reports if dynamic adjustment of thread number is allowed.
!
!  Discussion:
!
!    The user can request dynamic thread adjustment by calling OMP_SET_DYNAMIC.
!
!    For this stub library, the value FALSE is always returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, logical OMP_GET_DYNAMIC, is TRUE if dynamic adjustment of thread
!    number has been enable, by default or by a user call.
!
  implicit none

  logical omp_get_dynamic

  omp_get_dynamic = .false.

  return
end function
function omp_get_level ( )

!*****************************************************************************80
!
!! OMP_GET_LEVEL returns the number of nested parallel regions enclosing this task.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, integer OMP_GET_LEVEL, the number of nested parallel regions 
!    enclosing this task.
!
  implicit none

  integer omp_get_level

  omp_get_level = 0

  return
end function
function omp_get_max_active_levels ( )

!*****************************************************************************80
!
!! OMP_GET_MAX_ACTIVE_LEVELS gets the maximum number of nested active parallel regions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, integer OMP_GET_MAX_ACTIVE_LEVELS gets the maximum number of 
!    nested active parallel regions.
!
  implicit none

  integer omp_get_max_active_levels

  omp_get_max_active_levels = 0

  return
end function
function omp_get_max_threads ( )

!*****************************************************************************80
!
!! OMP_GET_MAX_THREADS returns the default number of threads.
!
!  Discussion:
!
!    If a parallel region is reached, and no number of threads has been
!    specified explicitly, there is a default number of threads that will
!    be used to form the new team.  That value is returned by this function.
!
!    For this stub library, the value 1 is always returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, integer OMP_GET_MAX_THREADS, the default number of threads.
!
  implicit none

  integer omp_get_max_threads

  omp_get_max_threads = 1

  return
end
function omp_get_nested ( )

!*****************************************************************************80
!
!! OMP_GET_NESTED reports if nested parallelism has been enabled.
!
!  Discussion:
!
!    The user can request nested parallelism by calling OMP_SET_NESTED.
!
!    For this stub library, the value FALSE is always returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, logical OMP_GET_NESTED, is TRUE if nested parallelism has been
!    enable by default or by a user call.
!
  implicit none

  logical omp_get_nested

  omp_get_nested = .false.

  return
end function
function omp_get_num_procs ( )

!*****************************************************************************80
!
!! OMP_GET_NUM_PROCS returns the number of processors available to the program.
!
!  Discussion:
!
!    For this stub library, the value 1 is always returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, integer GET_NUM_PROCS, the number of processors available.
!
  implicit none

  integer omp_get_num_procs 

  omp_get_num_procs = 1

  return
end function
function omp_get_num_threads ( )

!*****************************************************************************80
!
!! OMP_GET_NUM_THREADS returns the number of threads in the current team.
!
!  Discussion:
!
!    For this stub library, the value 1 is always returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, integer OMP_GET_NUM_THREADS, the number of threads in the 
!    current team.
!
  implicit none

  integer omp_get_num_threads

  omp_get_num_threads = 1

  return
end function
subroutine omp_get_schedule ( kind, modifier )

!*****************************************************************************80
!
!! OMP_GET_SCHEDULE returns information about the "runtime" schedule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!    For the stub library, the static schedule is returned.
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, integer ( kind = omp_sched_kind ) KIND, may be
!    1, omp_sched_static,
!    2, omp_sched_dynamic,
!    3, omp_sched_guided,
!    4, omp_sched_auto.
!
!    Output, integer MODIFIER; this contains the "chunk_size" information for
!    static, dynamic, or guided schedules, and is ignored for the auto schedule.
!    If the chunk_size is less than 1, then the default value is used instead.
!
  implicit none

  include 'omp_lib_kinds.h'
  integer ( kind = omp_sched_kind ) kind
  integer modifier
 
  kind = omp_sched_static
  modifier = 0

  return
end subroutine
function omp_get_team_size ( level )

!*****************************************************************************80
!
!! OMP_GET_TEAM_SIZE returns the size of the thread team for a given level.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Input, integer LEVEL, the nested level.
!
!    Output, integer OMP_GET_TEAM_SIZE, the size of the thread team for 
!    this level.
!
  implicit none

  integer level
  integer omp_get_team_size

  if ( level == 0 ) then
    omp_get_team_size = 1
  else
    omp_get_team_size = -1
  end if

  return
end function
function omp_get_thread_limit ( )

!*****************************************************************************80
!
!! OMP_GET_THREAD_LIMIT returns the maximum number of OpenMP threads available.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, integer OMP_GET_THREAD_LIMIT, the maximum number of OpenMP
!    threads available.
!
  implicit none

  integer omp_get_thread_limit

  omp_get_thread_limit = 1

  return
end function
function omp_get_thread_num ( )

!*****************************************************************************80
!
!! OMP_GET_THREAD_NUM returns the thread number of a thread in a team.
!
!  Discussion:
!
!    Thread numbers start at 0.
!
!    If this function is not called from a parallel region, then only one
!    thread is executing, so the value returned is 0.
!
!    For this stub library, the value 0 is always returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, integer OMP_GET_THREAD_NUM, the thread number.
!
  implicit none

  integer omp_get_thread_num

  omp_get_thread_num = 0 

  return
end function
function omp_get_wtick ( )

!*****************************************************************************80
!
!! OMP_GET_WTICK returns the precision of the timer used by OMP_GET_WTIME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, double precision OMP_GET_WTICK, the number of seconds between
!    successive "ticks" of the wall clock timer.
!
  implicit none

  integer count
  integer count_max
  integer count_rate
  double precision omp_get_wtick

  call system_clock ( count, count_rate, count_max )

  omp_get_wtick = 1.0D+00 / dble ( count_rate )

  return
end
function omp_get_wtime ( )

!*****************************************************************************80
!
!! OMP_GET_WTIME returns elapsed wall clock time in seconds.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, double precision OMP_GET_WTIME, the current reading of the
!    wall clock timer.
!
  implicit none

  integer count
  integer count_max
  integer count_rate
  double precision omp_get_wtime

  call system_clock ( count, count_rate, count_max )

  omp_get_wtime = dble ( count ) / dble ( count_rate )

  return
end
function omp_in_final ( )

!*****************************************************************************80
!
!! OMP_IN_FINAL is true if the routine is executed in a final task region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, logical OMP_IN_FINAL, is true if the routine is executed in a
!    final task region.
!
  implicit none

  logical omp_in_final

  omp_in_final = .true.

  return
end function
function omp_in_parallel ( )

!*****************************************************************************80
!
!! OMP_IN_PARALLEL returns TRUE if the call is made from a parallel region.
!
!  Discussion:
!
!    For this stub library, the value FALSE is always returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, logical OMP_IN_PARALLEL, is TRUE if the routine was called
!    from a parallel region.
!
  implicit none

  logical omp_in_parallel

  omp_in_parallel = .false.

  return
end function
subroutine omp_init_lock ( lock )

!*****************************************************************************80
!
!! OMP_INIT_LOCK initializes a simple lock.
!
!  Discussion:
!
!    This routine is intended to initialize the lock to the unlocked state.
!
!    For this stub library, the lock is set to -1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, integer ( kind = omp_lock_kind ) LOCK, the lock.
!    0, if the simple lock is not initialized 
!	-1, if the simple lock is initialized but not set 
!	 1, if the simple lock is set 
!
  implicit none

  include 'omp_lib_kinds.h' 

  integer ( kind = omp_lock_kind ) lock

  lock = -1

  return
end subroutine
subroutine omp_init_nest_lock ( nlock ) 

!*****************************************************************************80
!
!! OMP_INIT_NEST_LOCK initializes a nestable lock.
!
!  Discussion:
!
!    This routine is intended to initialize the lock to the unlocked state.
!
!    For this stub library, the lock is set to -1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Output, integer ( kind = omp_nest_lock_kind ) NLOCK, the lock.
!    0, if the nestable lock is not initialized;
!   -1, if the nestable lock is initialized but not set;
!    1, if the nestable lock is set 
!
  implicit none

  include 'omp_lib_kinds.h'
  integer ( kind = omp_nest_lock_kind ) nlock

  nlock = -1

  return
end subroutine
subroutine omp_set_dynamic ( dynamic_threads )

!*****************************************************************************80
!
!! OMP_SET_DYNAMIC enables dynamic adjustment of the number of threads.
!
!  Discussion:
!
!    If DYNAMIC_THREADS is TRUE, then the number of threads available for
!    execution in a parallel region may be dynamically adjusted.
!
!    For this stub library, the input value is ignored.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Input, logical DYNAMIC_THREADS, is TRUE if the user wishes to allow
!    dynamic adjustment of the number of threads available for execution
!    in any parallel region.
!
  implicit none

  logical dynamic_threads

  return
end subroutine
subroutine omp_set_lock ( lock )

!*****************************************************************************80
!
!! OMP_SET_LOCK sets a simple lock.
!
!  Discussion:
!
!    The lock must already have been initialized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Input/output, integer ( kind = omp_lock_kind ) LOCK, the simple lock.
! 
  implicit none

  include 'omp_lib_kinds.h'
  integer ( kind = omp_lock_kind ) lock

  if ( lock == -1 ) then
    lock = 1
  else if ( lock == 1 ) then 
    print *, 'error in omp_set_lock: deadlock in using lock variable.'
    stop
  else
    print *, 'error in omp_set_lock: lock not initialized.' 
    stop
  end if

  return
end subroutine
subroutine omp_set_max_active_levels ( max_levels )

!*****************************************************************************80
!
!! OMP_SET_MAX_ACTIVE_LEVELS limits the number of nested active parallel regions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Input, integer MAX_LEVELS, the maximum number of nested active parallel
!    regions.
!
  implicit none

  integer max_levels

  return
end subroutine
subroutine omp_set_nest_lock ( nlock )

!*****************************************************************************80
!
!! OMP_SET_NEST_LOCK sets a nestable lock.
!
!  Discussion:
!
!    The lock must already have been initialized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Input/output, integer ( kind = omp_nest_lock_kind ) NLOCK, the nestable lock.
!
  implicit none

  include 'omp_lib_kinds.h'
  integer ( kind = omp_nest_lock_kind ) nlock

  if ( nlock == -1 ) then
    nlock = 1
  else if ( nlock == 0 ) then
    print *, 'error in omp_set_nest_lock: nested lock not initialized'
    stop
  else
    print *, 'error in omp_set_nest_lock: deadlock using nested lock variable.'
    stop
  end if

  return 
end subroutine
subroutine omp_set_nested ( nested )

!*****************************************************************************80
!
!! OMP_SET_NESTED controls the use of nested parallelism.
!
!  Discussion:
!
!    For this stub library, the input value is ignored.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Input, logical NESTED, is TRUE if nested parallelism is to be enabled.
!
  implicit none
 
  logical nested

  return
end subroutine
subroutine omp_set_num_threads ( num_threads )

!*****************************************************************************80
!
!! OMP_SET_NUM_THREADS sets the number of threads.
!
!  Discussion:
!
!    This routine sets the number of threads to be used in all subsequent
!    parallel regions for which an explicit number of threads is not
!    specified.
!
!    For this stub library, the input value is ignored.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Input, integer NUM_THREADS, the number of threads to be used in all
!    subsequent parallel regions for which an explicit number of threads
!    is not specified.  0 < NUM_THREADS.
!
  implicit none

  integer num_threads

  return
end subroutine
subroutine omp_set_schedule ( kind, modifier )

!*****************************************************************************80
!
!! OMP_SET_SCHEDULE chooses the schedule when "runtime" is the schedule kind.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Input, integer ( kind = omp_sched_kind ) KIND, may be
!    1, omp_sched_static,
!    2, omp_sched_dynamic,
!    3, omp_sched_guided,
!    4, omp_sched_auto.
!
!    Input, integer MODIFIER; this contains the "chunk_size" information for
!    static, dynamic, or guided schedules, and is ignored for the auto schedule.
!    If the chunk_size is less than 1, then the default value is used instead.
!
  implicit none

  include 'omp_lib_kinds.h'
  integer ( kind = omp_sched_kind ) kind
  integer modifier

  return
end subroutine
function omp_test_lock ( lock )

!*****************************************************************************80
!
!! OMP_TEST_LOCK tests a simple lock.
!
!  Discussion:
!
!    Calling this routine with an uninitialized lock causes an error.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Input/output, integer ( kind = omp_lock_kind ) LOCK, the simple lock.
!    If the lock was initialized but not set, it is set by this call.
!
!    Output, logical OMP_TEST_LOCK:
!    TRUE, on input, the lock was initialized and not set;
!    FALSE, on input the lock was initialized, and set.
!
  implicit none

  include 'omp_lib_kinds.h'
  integer ( kind = omp_lock_kind ) lock
  logical omp_test_lock

  if ( lock == -1 ) then
    lock = 1
    omp_test_lock = .true. 
  else if ( lock == 1 ) then
    omp_test_lock = .false.
  else 
    print *, 'error in omp_test_lock: lock not initialized' 
    stop
  end if

  return 
end function
function omp_test_nest_lock ( nlock ) 

!*****************************************************************************80
!
!! OMP_TEST_NEST_LOCK tests a nestable lock.
!
!  Discussion:
!
!    Calling this routine with an uninitialized lock causes an error.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Input/output, integer ( kind = omp_nest_lock_kind ) NLOCK, the nestable lock.
!    If the lock was initialized but not set, it is set by this call.
!
!    Output, integer OMP_TEST_NEST_LOCK, returns the new nesting count,
!    if the call was successful.  Otherwise, the value 0 is returned.
!
  implicit none

  include 'omp_lib_kinds.h'
  integer ( kind = omp_nest_lock_kind ) nlock
  integer omp_test_nest_lock

  if ( nlock == -1 ) then
    nlock = 1
    omp_test_nest_lock = 1
  else if ( nlock == 1 ) then
    omp_test_nest_lock = 0
  else
    print *, 'error in omp_test_nest_lock: nested lock not initialized'
    stop
  end if

  return 
end function
subroutine omp_unset_lock ( lock )

!*****************************************************************************80
!
!! OMP_UNSET_LOCK unsets a simple lock.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Input/output, integer ( kind = omp_lock_kind ) LOCK, the simple lock.
!
  implicit none

  include 'omp_lib_kinds.h'
  integer ( kind = omp_lock_kind ) lock

  if ( lock == 1 ) then
    lock = -1
  else if ( lock == -1 ) then 
    print *, 'error in omp_unset_lock: lock not set'
    stop
  else
    print *, 'error in omp_unset_lock: lock not initialized'
    stop
  end if

  return 
end subroutine
subroutine omp_unset_nest_lock ( nlock )

!*****************************************************************************80
!
!! OMP_UNSET_NEST_LOCK unsets a nestable lock.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    OpenMP Application Program Interface,
!    Version 3.1,
!    July 2011.
!
!  Parameters:
!
!    Input/output, integer ( kind = omp_nest_lock_kind ) NLOCK, the nestable lock.
!
  implicit none

  include 'omp_lib_kinds.h'
  integer ( kind = omp_nest_lock_kind ) nlock

  if ( nlock == 1 ) then
    nlock = -1
  elseif ( nlock == 0 ) then 
    print *, 'error in omp_unset_nest_lock: nested lock not initialized'
    stop
  else
    print *, 'error in omp_unset_nest_lock: nested lock not set'
    stop
  end if

  return 
end subroutine
