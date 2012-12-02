!  omp_lib_kinds.h
!
  integer, parameter :: omp_lock_kind = selected_int_kind ( 10 )
  integer, parameter :: omp_nest_lock_kind = selected_int_kind ( 10 )
  integer, parameter :: omp_sched_kind = selected_int_kind ( 8 )

  integer, parameter :: omp_sched_static = 1
  integer, parameter :: omp_sched_dynamic = 2
  integer, parameter :: omp_sched_guided = 3
  integer, parameter :: omp_sched_auto = 4
