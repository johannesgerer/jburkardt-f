program main

!*****************************************************************************80
!
!! MAIN is the main program for BLACS_PRB.
!
!  Discussion:
!
!    BLACS_PRB is a test program for the BLACS.
!
!    The BLACS is a library of Basic Linear Algebra Communications
!    Subroutines which can facilitate the solution of linear algebra
!    computations that use message passing.
!
!  Modified:
!
!    20 November 2003
!
!  Reference:
!
!    Jack Dongarra, Clint Whaley,
!    LAPACK Working Note 94:
!    A User's Guide to the BLACS v1.1,
!    pages 60-62.
!
  implicit none

  integer blacs_pnum
  integer caller
  integer contxt
  integer i
  integer ios
  integer j
  integer pcol_me
  integer pcol_num
  integer proc_me
  integer proc_num
  integer prow_me
  integer prow_num
  integer you
  integer yourcol
  integer yourrow
!
!  The call to BLACS_PINFO will return the process number assigned to
!  this process, and the total number of processes.
!
  call blacs_pinfo ( proc_me, proc_num )

  if ( proc_me == 0 ) then
    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Process: ', proc_me
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BLACS_PRB'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) '  Test the BLACS library.'
    write ( *, '(a)') ' '
    write ( *, '(a)' ) '  A sample program for the BLACS.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  In this simple example, we begin with a'
    write ( *, '(a)' ) '  certain number of processes.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  These processes are arranged into a 2D'
    write ( *, '(a)' ) '  computational grid (and any left over processes'
    write ( *, '(a)' ) '  will exit.)'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Then BLACS is initialized, and the (0,0) process'
    write ( *, '(a)' ) '  expects to receive a "check-in" message from all'
    write ( *, '(a)' ) '  other active processes.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  The number of processes is ', proc_num
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Process: ', proc_me
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BLACS_PRB'
    write ( *, '(a)' ) '  Process beginning.'
  end if
!
!  If in PVM, create the virtual machine if it doesn't exist.
!
  if ( proc_num < 1 ) then
    if ( proc_me == 0 ) then
      write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'BLACS_PRB - Process ', proc_me
      write ( *, '(a)' ) '  Please enter the number of processes.'
      read ( *, *, iostat = ios ) proc_num
      if ( ios /= 0 ) then 
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) 'BLACS_PRB - Process ', proc_me
        write ( *, '(a)' ) '  Abnormal end of execution.'
        write ( *, '(a)' ) '  Could not input number of processes.'
        stop
      end if
      call blacs_setup ( proc_me, proc_num )
    end if
  end if
!
!  Set up the process grid.
!
  prow_num = int ( sqrt ( real ( proc_num ) ) )
  pcol_num = proc_num / prow_num

  if ( proc_me == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'BLACS_PRB - Process ', proc_me
    write ( *, '(a)' ) '  Setting up the 2D process grid.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i3)' ) '  prow_num = ', prow_num
    write ( *, '(a,i3)' ) '  pcol_num = ', pcol_num
  end if
!
!  Get the default system context.
!
  call blacs_get ( 0, 0, contxt )
!
!  Define the grid.
!
  if ( proc_me == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'BLACS_PRB - Process ', proc_me
    write ( *, '(a)' ) '  Calling BLACS_GRIDINIT to define the grid.'
  end if

  call blacs_gridinit ( contxt, 'ROW', prow_num, pcol_num )
!
!  Get this process's row and column grid coordinates.
!
  if ( proc_me == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'BLACS_PRB - Process ', proc_me
    write ( *, '(a)' ) '  Calling BLACS_GRIDINFO for process grid coordinates.'
  end if

  call blacs_gridinfo ( contxt, prow_num, pcol_num, prow_me, pcol_me )
!
!  If this process is not in the grid, exit.
!
  if ( prow_num <= prow_me .or. pcol_num <= pcol_me ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BLACS_PRB - Process ', proc_me
    write ( *, '(a)' ) '  Not part of the grid, exiting.'
    call blacs_exit ( 0 )
    stop
  end if
!
!  Get the process id from the grid coordinates.
!
  if ( proc_me == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'BLACS_PRB - Process ', proc_me
    write ( *, '(a)' ) '  Call BLACS_PNUM for process id from grid coordinates.'
  end if

  caller = blacs_pnum ( contxt, prow_me, pcol_me )
!
!  Process (0,0) RECEIVES check-in messages from all other processes.
!
  if ( prow_me == 0 .and. pcol_me == 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Process(0,0):'
    write ( *, '(a)' ) '  All other processes must send me a check in message.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    Process     BLACS'
    write ( *, '(a)' ) '   Row   Col    ID'
    write ( *, '(a)' ) ' '

    do i = 0, prow_num - 1
      do j = 0, pcol_num - 1

        if ( i == 0 .and. j == 0 ) then
          you = caller
        else
          call igerv2d ( contxt, 1, 1, you, 1, i, j )
        end if
!
!  From the remote process's rank, determine its grid coordinates.
!
        call blacs_pcoord ( contxt, you, yourrow, yourcol )
!
!  If the grid coordinates are not what we expect, fail.
!
        if ( yourrow /= i .or. yourcol /= j ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'BLACS_PRB - Master (0,0) process:'
          write ( *, '(a,i6)' ) '  A grid error has occurred with process ', you
          stop
        end if

        write ( *, '(3i6)' ) i, j, you

      end do
    end do
!
!  Non-master processes SEND their process number to BLACS process (0,0).
!
  else

    call igesd2d ( contxt, 1, 1, caller, 1, 0, 0 )

  end if

  call blacs_exit ( 0 )

  if ( prow_me == 0 .and. pcol_me == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BLACS_PRB - Master (0,0) process:'
    write ( *, '(a)' ) '  Normal end of execution.'
    write ( *, '(a)' ) ' '
    call timestamp ( )
  end if

  stop
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
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

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
