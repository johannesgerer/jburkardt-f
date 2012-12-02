program main

!*****************************************************************************80
!
!! MAIN is the main program for MPI_MULTITASK.
!
!  Discussion:
!
!    Message tag 1: P0 sends input to P1
!    Message tag 2: P0 sends input to P2
!    Message tag 3: P1 sends output to P0.
!    Message tag 4: P2 sends output to P0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2011
!
!  Author:
!
!    John Burkardt
!
  use mpi

  implicit none

  integer ( kind = 4 ) id
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) input1
  integer ( kind = 4 ) input2
  integer ( kind = 4 ) my_id
  integer ( kind = 4 ) output1
  integer ( kind = 4 ) output2
  integer ( kind = 4 ) p
  real ( kind = 8 ) wtime
!
!  Process 0 is the "monitor".
!  It chooses the inputs, and sends them to the workers.
!  It waits for the outputs.
!  It plots the outputs.
!
  call MPI_Init ( ierr )

  call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )

  call MPI_Comm_size ( MPI_COMM_WORLD, p, ierr )
!
!  Make sure we have enough processes.
!
  if ( p < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MPI_MULTITASK - Fatal error!'
    write ( *, '(a)' ) '  Number of available processes must be at least 3!'
    call MPI_Finalize ( ierr )
    stop
  end if
!
!  Run program P0 on process 0, and so on.
!
  if ( id == 0 ) then

    call timestamp ( )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MPI_MULTITASK:'
    write ( *, '(a)' ) '  FORTRAN90 / MPI version'

    wtime = MPI_Wtime ( )

    call p0_set_input ( input1, input2 )
    call p0_send_input ( input1, input2 )
    call p0_receive_output ( output1, output2 )
!
!  Get timing.
!
    wtime = MPI_Wtime ( ) - wtime
    write ( *, '(a,g14.6)' ) '  Process 0 time = ', wtime

    call MPI_Finalize ( ierr )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MPI_MULTITASK:'
    write ( *, '(a)' ) '  Normal end of execution.'

    call timestamp ( )

    stop
!
!  Process 1 works on task 1.
!  It receives input from process 0.
!  It computes the output.
!  It sends the output to process 0.
!
  else if ( id == 1 ) then

    wtime = MPI_Wtime ( )
    call p1_receive_input ( input1 )
    call p1_compute_output ( input1, output1 )
    call p1_send_output ( output1 )
    wtime = MPI_Wtime ( ) - wtime
    write ( *, '(a,g14.6)' ) '  Process 1 time = ', wtime
    call MPI_Finalize ( ierr )
    stop
!
!  Process 2 works on task 2.
!  It receives input from process 0.
!  It computes the output.
!  It sends the output to process 0.
!
  else if ( id == 2 ) then

    wtime = MPI_Wtime ( )
    call p2_receive_input ( input2 )
    call p2_compute_output ( input2, output2 )
    call p2_send_output ( output2 )
    wtime = MPI_Wtime ( ) - wtime
    write ( *, '(a,g14.6)' ) '  Process 2 time = ', wtime
    call MPI_Finalize ( ierr )
    stop

  end if

end
subroutine p0_set_input ( input1, input2 )

!*****************************************************************************80
!
!! P0_SET_INPUT sets input.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) INPUT1, INPUT2, the values of two
!    inputs used by tasks 1 and 2.
!
  implicit none

  integer ( kind = 4 ) input1
  integer ( kind = 4 ) input2

  input1 = 10000000
  input2 = 100000

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) 'P0_SET_PARAMETERS:'
  write ( *, '(a,i12)' ) '  Set INPUT1 = ', input1
  write ( *, '(a,i12)' ) '      INPUT2 = ', input2

  return
end
subroutine p0_send_input ( input1, input2 )

!*****************************************************************************80
!
!! P0_SEND_INPUT sends input to processes 1 and 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT1, INPUT2, the values of two
!    inputs used by tasks 1 and 2.
!
  use mpi

  implicit none

  integer ( kind = 4 ) id
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) input1
  integer ( kind = 4 ) input2
  integer ( kind = 4 ) tag

  id = 1
  tag = 1
  call MPI_Send ( input1, 1, MPI_INTEGER, id, tag, MPI_COMM_WORLD, ierr )

  id = 2
  tag = 2
  call MPI_Send ( input2, 1, MPI_INTEGER, id, tag, MPI_COMM_WORLD, ierr )

  return
end
subroutine p0_receive_output ( output1, output2 )

!*****************************************************************************80
!
!! P0_RECEIVE_OUTPUT receives output from processes 1 and 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) OUTPUT1, OUTPUT2, the values of the
!    outputs of tasks 1 and 2.
!
  use mpi

  implicit none

  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) output
  integer ( kind = 4 ) output_received
  integer ( kind = 4 ) output1
  integer ( kind = 4 ) output2
  integer ( kind = 4 ) source
  integer status(MPI_STATUS_SIZE)

  output_received = 0
!
!  Loop until every worker has checked in.
!
  do while ( output_received < 2 )
!
!  Receive the next message that arrives.
!
    call MPI_Recv ( output, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, &
      MPI_COMM_WORLD, status, ierr )
!
!  The actual source of the message is saved in STATUS.
!
    source = status(MPI_SOURCE)
!
!  Save the value in OUTPUT1 or OUTPUT2.
!
    if ( source == 1 ) then
      output1 = output
    else
      output2 = output
    end if

    output_received = output_received + 1

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Process 1 returned OUTPUT1 = ', output1
  write ( *, '(a,i8)' ) '  Process 2 returned OUTPUT2 = ', output2

  return
end
subroutine p1_receive_input ( input1 )

!*****************************************************************************80
!
!! P1_RECEIVE_INPUT receives input from process 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) INPUT1, the value of the parameter.
!
  use mpi

  implicit none

  integer ( kind = 4 ) id
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) input1
  integer ( kind = 4 ) status(MPI_STATUS_SIZE)
  integer ( kind = 4 ) tag

  id = 0
  tag = 1
  call MPI_Recv ( input1, 1, MPI_INTEGER, id, tag, MPI_COMM_WORLD, status, ierr )

  return
end
subroutine p1_compute_output ( input1, output1 )

!*****************************************************************************80
!
!! P1_COMPUTE_OUTPUT carries out computation number 1.
!
!  Discussion:
!
!    No MPI calls occur in this function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT1, the problem input.
!
!    Output, integer ( kind = 4 ) OUTPUT1, the problem output.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) input1
  integer ( kind = 4 ) output1

  output1 = 0

  do i = 2, input1

    j = i
    k = 0

    do while ( 1 < j )

      if ( mod ( j, 2 ) == 0 ) then
        j = j / 2
      else
        j = 3 * j + 1
      end if

      k = k + 1

    end do

    output1 = max ( output1, k )

  end do

  return
end
subroutine p1_send_output ( output1 )

!*****************************************************************************80
!
!! P1_SEND_OUTPUT sends output to process 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OUTPUT1, the problem output.
!
  use mpi

  implicit none

  integer ( kind = 4 ) id
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) output1
  integer ( kind = 4 ) tag

  id = 0
  tag = 3
  call MPI_Send ( output1, 1, MPI_INTEGER, id, tag, MPI_COMM_WORLD, ierr )

  return
end
subroutine p2_receive_input ( input2 )

!*****************************************************************************80
!
!! P2_RECEIVE_INPUT receives input from process 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) INPUT2, the value of the parameter.
!
  use mpi

  implicit none

  integer ( kind = 4 ) id
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) input2
  integer ( kind = 4 ) status(MPI_STATUS_SIZE)
  integer ( kind = 4 ) tag

  id = 0
  tag = 2
  call MPI_Recv ( input2, 1, MPI_INTEGER, id, tag, MPI_COMM_WORLD, status, ierr )

  return
end
subroutine p2_compute_output ( input2, output2 )

!*****************************************************************************80
!
!! P2_COMPUTE_OUTPUT carries out computation number 2.
!
!  Discussion:
!
!    No MPI calls occur in this function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INPUT2, the problem input.
!
!    Output, integer ( kind = 4 ) OUTPUT2, the problem output.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) input2
  integer ( kind = 4 ) output2
  logical prime

  output2 = 0

  do i = 2, input2

    prime = .true.
    do j = 2, i - 1
      if ( mod ( i, j ) == 0 ) then
        prime = .false.
        exit
      end if
    end do

    if ( prime ) then
      output2 = output2 + 1
    end if

  end do

  return
end
subroutine p2_send_output ( output2 )

!*****************************************************************************80
!
!! P2_SEND_OUTPUT sends output to process 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OUTPUT2, the problem output.
!
  use mpi

  implicit none

  integer ( kind = 4 ) id
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) output2
  integer ( kind = 4 ) tag

  id = 0
  tag = 4
  call MPI_Send ( output2, 1, MPI_INTEGER, id, tag, MPI_COMM_WORLD, ierr )

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 )  ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
