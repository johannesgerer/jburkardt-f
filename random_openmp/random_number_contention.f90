program main

!*****************************************************************************80
!
!! RANDOM_NUMBER_CONTENTION demonstrates a memory contention problem.
!
!  Discussion:
!
!    This example illustrates a problem that can occur when using OpenMP.
!
!    The problem can crop up in a variety of places, but here it shows up
!    in the use of the FORTRAN90 random number generator.
!
!    The FORTRAN90 random number generator includes internal status 
!    variables.  When called in a shared memory environment, these variables
!    cause contention between the various threads of execution, as they
!    try to access and modify these internal status variables.
!
!    The surprising symptom of this contention is that the program can run
!    much more slowly when more than one thread is used.
!
!    This memory contention should be expected, because of the likely
!    structure of the random number generator, no matter what compiler
!    is used.  The contention is not guaranteed, though, and the severity
!    of the effect would vary with compiler and computer.
!
!    This example was inspired by a similar example posted on the 
!    CMISS wiki page.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 April 2008
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  integer i
  integer j
  integer :: nt = 10000
  integer proc_num
  double precision r
  integer thread_max
  integer thread_num
  double precision wtime

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDOM_NUMBER_CONTENTION'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This program demonstrates a case where the use'
  write ( *, '(a)' ) '  of multiple threads in OpenMP may DEGRADE the'
  write ( *, '(a)' ) '  performance of a program.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this program, we simply call the FORTRAN90'
  write ( *, '(a)' ) '  random number generator many times.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Because that function has internal saved variables,'
  write ( *, '(a)' ) '  when it is accessed by multile threads, there can'
  write ( *, '(a)' ) '  be huge delays because the threads fight to gain'
  write ( *, '(a)' ) '  access to the saved variables.'
!
!  How many processors are available?
!
  proc_num = omp_get_num_procs ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The number of processors available:'
  write ( *, '(a,i8)' ) '  OMP_GET_NUM_PROCS () = ', proc_num
  
  thread_max = proc_num
  if ( 8 < thread_max ) then
    thread_max = 8
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Threads    Wallclock time'
  write ( *, '(a)' ) '              (seconds)'
  write ( *, '(a)' ) ' '

  do thread_num = 1, thread_max

    call omp_set_num_threads ( thread_num )

    wtime = omp_get_wtime ( )

!$omp parallel shared ( nt ) private ( i, j, r )

!$omp do
    do i = 1, 100
      do j = 1, nt
        call random_number ( r )
      end do
    end do
!$omp end do

!$omp end parallel

    wtime = omp_get_wtime ( ) - wtime

    write ( *, '(2x,i8,2x,g14.6)' ) thread_num, wtime

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDOM_NUMBER_CONTENTION'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

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

