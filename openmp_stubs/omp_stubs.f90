module omp_lib_kinds
  integer, parameter :: omp_lock_kind = selected_int_kind( 10 )
  integer, parameter :: omp_nest_lock_kind = selected_int_kind( 10 )
  integer, parameter :: omp_sched_kind = selected_int_kind( 8 )
  integer(kind=omp_sched_kind), parameter :: omp_sched_static = 1
  integer(kind=omp_sched_kind), parameter :: omp_sched_dynamic = 2
  integer(kind=omp_sched_kind), parameter :: omp_sched_guided = 3 
  integer(kind=omp_sched_kind), parameter :: omp_sched_auto = 4 
end module omp_lib_kinds

module omp_lib

  use omp_lib_kinds

  integer, parameter :: openmp_version = 201107

  interface
    subroutine omp_set_num_threads (number_of_threads_expr)
      integer, intent(in) :: number_of_threads_expr
    end subroutine omp_set_num_threads
    function omp_get_num_threads ()
      integer :: omp_get_num_threads
    end function omp_get_num_threads
    function omp_get_max_threads ()
      integer :: omp_get_max_threads
    end function omp_get_max_threads
    function omp_get_thread_num ()
      integer :: omp_get_thread_num
    end function omp_get_thread_num
    function omp_get_num_procs ()
      integer :: omp_get_num_procs
    end function omp_get_num_procs
    function omp_in_parallel ()
      logical :: omp_in_parallel 
    end function omp_in_parallel
    subroutine omp_set_dynamic (enable_expr)
      logical, intent(in) :: enable_expr
    end subroutine omp_set_dynamic
    function omp_get_dynamic ()
      logical :: omp_get_dynamic
    end function omp_get_dynamic
    subroutine omp_set_nested (enable_expr)
      logical, intent(in) :: enable_expr
    end subroutine omp_set_nested
    function omp_get_nested ()
      logical :: omp_get_nested
    end function omp_get_nested
    subroutine omp_set_schedule (kind, modifier)
      use omp_lib_kinds
      integer(kind=omp_sched_kind), intent(in) :: kind
      integer, intent(in) :: modifier
    end subroutine omp_set_schedule
    subroutine omp_get_schedule (kind, modifier)
      use omp_lib_kinds
      integer(kind=omp_sched_kind), intent(out) :: kind
      integer, intent(out)::modifier
    end subroutine omp_get_schedule
    function omp_get_thread_limit()
      integer :: omp_get_thread_limit
    end function omp_get_thread_limit
    subroutine omp_set_max_active_levels(var)
      integer, intent(in) :: var
    end subroutine omp_set_max_active_levels
    function omp_get_max_active_levels()
      integer :: omp_get_max_active_levels
    end function omp_get_max_active_levels
    function omp_get_level()
      integer :: omp_get_level
    end function omp_get_level
    function omp_get_ancestor_thread_num(level)
      integer, intent(in) :: level
      integer :: omp_get_ancestor_thread_num
    end function omp_get_ancestor_thread_num
    function omp_get_team_size(level)
      integer, intent(in) :: level
      integer :: omp_get_team_size
    end function omp_get_team_size
    function omp_get_active_level()
      integer :: omp_get_active_level
    end function omp_get_active_level
    function omp_in_final()
      logical omp_in_final
    end function omp_in_final
    subroutine omp_init_lock (var)
      use omp_lib_kinds
      integer (kind=omp_lock_kind), intent(out) :: var
    end subroutine omp_init_lock
    subroutine omp_destroy_lock (var)
      use omp_lib_kinds
      integer (kind=omp_lock_kind), intent(inout) :: var
    end subroutine omp_destroy_lock
    subroutine omp_set_lock (var)
      use omp_lib_kinds
      integer (kind=omp_lock_kind), intent(inout) :: var
    end subroutine omp_set_lock
    subroutine omp_unset_lock (var)
      use omp_lib_kinds
      integer (kind=omp_lock_kind), intent(inout) :: var
    end subroutine omp_unset_lock
    function omp_test_lock (var)
      use omp_lib_kinds
      logical :: omp_test_lock
      integer (kind=omp_lock_kind), intent(inout) :: var
    end function omp_test_lock
    subroutine omp_init_nest_lock (var)
      use omp_lib_kinds
      integer (kind=omp_nest_lock_kind), intent(out) :: var
    end subroutine omp_init_nest_lock
    subroutine omp_destroy_nest_lock (var)
      use omp_lib_kinds
      integer (kind=omp_nest_lock_kind), intent(inout) :: var
    end subroutine omp_destroy_nest_lock
    subroutine omp_set_nest_lock (var)
      use omp_lib_kinds
      integer (kind=omp_nest_lock_kind), intent(inout) :: var
    end subroutine omp_set_nest_lock
    subroutine omp_unset_nest_lock (var)
      use omp_lib_kinds
      integer (kind=omp_nest_lock_kind), intent(inout) :: var
    end subroutine omp_unset_nest_lock
    function omp_test_nest_lock (var)
      use omp_lib_kinds
      integer :: omp_test_nest_lock
      integer (kind=omp_nest_lock_kind), intent(inout) :: var
    end function omp_test_nest_lock
    function omp_get_wtick ()
      double precision :: omp_get_wtick
    end function omp_get_wtick
    function omp_get_wtime ()
      double precision :: omp_get_wtime
    end function omp_get_wtime

  end interface 

end module omp_lib
