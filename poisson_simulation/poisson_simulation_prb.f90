program main

!*****************************************************************************80
!
!! POISSON_SIMULATION_TEST tests POISSON_SIMULATION.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'POISSON_SIMULATION_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the POISSON_SIMULATION library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'POISSON_SIMULATION_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 simulates waiting for a given number of events.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: bin_num = 30
  integer ( kind = 4 ), parameter :: event_num = 1000

  character ( len = 80 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 80 ) data_filename
  integer ( kind = 4 ) data_unit
  integer ( kind = 4 ) f_bin(bin_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t(0:event_num)
  real ( kind = 8 ) w(0:event_num)
  real ( kind = 8 ) w_bin(bin_num)
  real ( kind = 8 ) w_max
  real ( kind = 8 ) w_min
  real ( kind = 8 ) width

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  POISSON_FIXED_EVENTS simulates a Poisson process'
  write ( *, '(a)' ) '  until a given number of events have occurred.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Simulate a Poisson process, for which, on average,'
  write ( *, '(a)' ) '  LAMBDA events occur per unit time.'
  write ( *, '(a)' ) '  Run until you have observed EVENT_NUM events.'
 
  lambda = 0.5D+00
  seed = 123456789

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  LAMBDA = ', lambda
  write ( *, '(a,i6)' ) '  EVENT_NUM = ', event_num

  call poisson_fixed_events ( lambda, event_num, seed, t, w )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Minimum wait = ', minval ( w(1:event_num) )
  write ( *, '(a,g14.6)' ) '  Average wait = ', &
    sum ( w(1:event_num) ) / real ( event_num, kind = 8 )
  write ( *, '(a,g14.6)' ) '  Maximum wait = ', maxval ( w(1:event_num) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) ' Count            Time            Wait'
  write ( *, '(a)' ) ''
  do i = 0, 5
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, t(i), w(i)
  end do
  write ( *, '(a)' ) '  ....  ..............  ..............'
  do i = event_num - 5, event_num 
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, t(i), w(i)
  end do
!
!  Create the data file.
!
  call get_unit ( data_unit )

  data_filename = 'poisson_timeline_data.txt'

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do i = 0, event_num
    write ( data_unit, '(2x,g14.6,2x,i8)' ) t(i), i
  end do
  close ( unit = data_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data stored in "' // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  call get_unit ( command_unit )

  command_filename = 'poisson_timeline_commands.txt'

  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# poisson_timeline_commands.txt'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < poisson_timeline_commands.txt'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "poisson_timeline.png"'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'set xlabel "Time"'
  write ( command_unit, '(a)' ) 'set ylabel "Number of Events"'
  write ( command_unit, '(a)' ) 'set title "Observation of Fixed Number of Poisson Events"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a,f8.2,a)' ) &
    'plot "poisson_timeline_data.txt" using 1:2 lw 2'
  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  Plot commands stored in "' &
    // trim ( command_filename ) // '".'
!
!  Determine bin information.
!
  w_min = minval ( w )
  w_max = maxval ( w )

  call r8vec_midspace ( bin_num, w_min, w_max, w_bin )

  f_bin(1:bin_num) = 0
  do i = 0, event_num
    j = 1 + int ( real ( bin_num, kind = 8 ) * ( w(i) - w_min ) &
      / ( w_max - w_min ) )
    j = min ( j, bin_num )
    f_bin(j) = f_bin(j) + 1
  end do
!
!  Create the data file.
!
  call get_unit ( data_unit )

  data_filename = 'poisson_times_data.txt'

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do i = 1, bin_num
    write ( data_unit, '(2x,g14.6,2x,i6)' ) w_bin(i), f_bin(i)
  end do
  close ( unit = data_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data stored in "' // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  call get_unit ( command_unit )

  command_filename = 'poisson_times_commands.txt'

  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# poisson_times_commands.txt'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < poisson_times_commands.txt'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "poisson_times.png"'
  write ( command_unit, '(a)' ) 'set xlabel "Waiting Time"'
  write ( command_unit, '(a)' ) 'set ylabel "Frequency"'
  write ( command_unit, '(a)' ) 'set title "Waiting Times Observed Over Fixed Time"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style fill solid'
  width = 0.85D+00 * ( w_max - w_min ) / real ( bin_num, kind = 8 )
  write ( command_unit, '(a,f8.2,a)' ) &
    'plot "poisson_times_data.txt" using 1:2:(', width, ') with boxes'
  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  Plot commands stored in "' &
    // trim ( command_filename ) // '".'

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 simulates waiting for a given length of time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: bin_num = 30
  integer ( kind = 4 ), parameter :: test_num = 20000

  character ( len = 80 ) command_filename
  integer ( kind = 4 ) command_unit
  character ( len = 80 ) data_filename
  integer ( kind = 4 ) data_unit
  real ( kind = 8 ) f_bin(bin_num)
  integer ( kind = 4 ) i
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) n(test_num)
  real ( kind = 8 ) n_bin(bin_num)
  real ( kind = 8 ) n_max
  real ( kind = 8 ) n_mean
  real ( kind = 8 ) n_min
  real ( kind = 8 ) n_var
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  integer ( kind = 4 ) test
  real ( kind = 8 ) w

  lambda = 0.5D+00
  t = 1000.0D+00
  seed = 123456789

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  POISSON_FIXED_EVENTS simulates a Poisson process'
  write ( *, '(a)' ) '  counting the number of events that occur during'
  write ( *, '(a)' ) '  a given time.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Simulate a Poisson process, for which, on average,'
  write ( *, '(a)' ) '  LAMBDA events occur per unit time.'
  write ( *, '(a,g14.6,a)' ) '  Run for a total of ', t, ' time units.'
  write ( *, '(a,g14.6)' ) '  LAMBDA = ', lambda

  do test = 1, test_num
    call poisson_fixed_time ( lambda, t, seed, n(test) )
  end do

  call i4vec_mean ( test_num, n, n_mean )
  call i4vec_variance ( test_num, n, n_var )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Mean number of events = ', n_mean
  write ( *, '(a,g14.6)' ) '  Variance = ', n_var  
  write ( *, '(a,g14.6)' ) '  STD = ', sqrt ( n_var )

  n_min = real ( minval ( n ), kind = 8 )
  n_max = real ( maxval ( n ), kind = 8 )

  call r8vec_midspace ( bin_num, n_min, n_max, n_bin )

  f_bin(1:bin_num) = 0
  do test = 1, test_num
    i = 1 + int ( real ( bin_num * ( n(test) - n_min ), kind = 8 ) &
      / real ( n_max - n_min, kind = 8 ) )
    i = min ( i, bin_num )
    f_bin(i) = f_bin(i) + 1
  end do
!
!  Create the data file.
!
  call get_unit ( data_unit )

  data_filename = 'poisson_events_data.txt'

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do i = 1, bin_num
    write ( data_unit, '(2(2x,g14.6))' ) n_bin(i), f_bin(i)
  end do
  close ( unit = data_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data stored in "' // trim ( data_filename ) // '".'
!
!  Create the command file.
!
  call get_unit ( command_unit )

  command_filename = 'poisson_events_commands.txt'

  open ( unit = command_unit, file = command_filename, status = 'replace' )

  write ( command_unit, '(a)' ) '# poisson_events_commands.txt'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < poisson_events_commands.txt'
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "poisson_events.png"'
  write ( command_unit, '(a)' ) 'set xlabel "Number of Events"'
  write ( command_unit, '(a)' ) 'set ylabel "Frequency"'
  write ( command_unit, '(a)' ) 'set title "Number of Poisson Events Over Fixed Time"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style fill solid'
  w = 0.85D+00 * ( n_max - n_min ) / real ( bin_num, kind = 8 )
  write ( command_unit, '(a,f8.2,a)' ) &
    'plot "poisson_events_data.txt" using 1:2:(', w, ') with boxes'
  write ( command_unit, '(a)' ) 'quit'

  close ( unit = command_unit )

  write ( *, '(a)' ) '  Plot commands stored in "' &
    // trim ( command_filename ) // '".'
  
  return
end
