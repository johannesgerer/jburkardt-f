  include 'omp_lib_kinds.h'

  parameter, integer         :: openmp_version = 201107

  external                   :: omp_set_num_threads
  external, integer          :: omp_get_num_threads
  external, integer          :: omp_get_max_threads
  external, integer          :: omp_get_thread_num
  external, integer          :: omp_get_num_procs
  external, logical          :: omp_in_parallel
  external                   :: omp_set_dynamic
  external, logical          :: omp_get_dynamic
  external                   :: omp_set_nested
  external, logical          :: omp_get_nested
  external                   :: omp_set_schedule
  external                   :: omp_get_schedule
  external, integer          :: omp_get_thread_limit
  external                   :: omp_set_max_active_levels
  external, integer          :: omp_get_max_active_levels
  external, integer          :: omp_get_level
  external, integer          :: omp_get_ancestor_thread_num
  external, integer          :: omp_get_team_size
  external, integer          :: omp_get_active_level
  external, logical          :: omp_in_final
  external                   :: omp_init_lock
  external                   :: omp_destroy_lock
  external                   :: omp_set_lock
  external                   :: omp_unset_lock
  external, logical          :: omp_test_lock
  external                   :: omp_init_nest_lock
  external                   :: omp_destroy_nest_lock
  external                   :: omp_set_nest_lock
  external                   :: omp_unset_nest_lock
  external, integer          :: omp_test_nest_lock
  external, double precision :: omp_get_wtick
  external, double precision :: omp_get_wtime