program main

!*****************************************************************************80
!
!! MAIN is the main program for MXV_OPENMP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 April 2008
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) proc_num
  integer ( kind = 4 ) thread_num

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MXV_OPENMP:'
  write ( *, '(a)' ) '  FORTRAN90/OpenMP version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compute matrix vector products y = A*x.'

  proc_num = omp_get_num_procs ( )
  thread_num = omp_get_max_threads ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of processors available = ', proc_num
  write ( *, '(a,i8)' ) '  The number of threads available    = ', thread_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare various algorithms:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  MXV_PLAIN          - "plain vanilla" FORTRAN.'
  write ( *, '(a)' ) '  MXV_PLAIN_OPENMP  - PLAIN + OpenMP.'
  write ( *, '(a)' ) '  MXV_MATMUL         - the FORTRAN90 MATMUL function.'
  write ( *, '(a)' ) '  MXV_MATMUL_OPENMP - MATMUL + OpenMP.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Algorithm                  M         N      Seconds'
!
!  N = M
!
  m = 10

  do i = 1, 3

    write ( *, '(a)' ) ' '

    n = m
    call test01 ( m, n )

    m = m * 10

  end do
!
!  N = 10 * M
!
  m = 1

  do i = 1, 4

    write ( *, '(a)' ) ' '

    n = 10 * m
    call test01 ( m, n )

    m = m * 10

  end do
!
!  M = 10 * N
!
  n = 1

  do i = 1, 4

    write ( *, '(a)' ) ' '

    m = 10 * n
    call test01 ( m, n )

    n = n * 10

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MXV_OPENMP:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( m, n )

!*****************************************************************************80
!
!! TEST01 compares various algorithms for a given matrix size MxN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    of the matrix.
!
  implicit none

  real    ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real    ( kind = 8 ) omp_get_wtime
  real    ( kind = 8 ) seconds
  real    ( kind = 8 ), allocatable, dimension ( : ) :: x
  real    ( kind = 8 ), allocatable, dimension ( : ) :: y

  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:m) )

  call matgen ( m, n, a, x )

  seconds = omp_get_wtime ( )
  call mxv_plain ( m, n, a, x, y )
  seconds = omp_get_wtime ( ) - seconds
  write ( *, '(2x,a18,2x,i8,2x,i8,2x,g14.6)' ) &
    'MXV_PLAIN         ', m, n, seconds

  seconds = omp_get_wtime ( )
  call mxv_plain_openmp ( m, n, a, x, y )
  seconds = omp_get_wtime ( ) - seconds
  write ( *, '(2x,a18,2x,i8,2x,i8,2x,g14.6)' ) &
    'MXV_PLAIN_OPENMP ', m, n, seconds

  seconds = omp_get_wtime ( )
  call mxv_matmul ( m, n, a, x, y )
  seconds = omp_get_wtime ( ) - seconds
  write ( *, '(2x,a18,2x,i8,2x,i8,2x,g14.6)' ) &
    'MXV_MATMUL        ', m, n, seconds

  seconds = omp_get_wtime ( )
  call mxv_matmul ( m, n, a, x, y )
  seconds = omp_get_wtime ( ) - seconds
  write ( *, '(2x,a18,2x,i8,2x,i8,2x,g14.6)' ) &
    'MXV_MATMUL_OPENMP', m, n, seconds

  deallocate ( a )
  deallocate ( x )
  deallocate ( y )

  return
end
subroutine matgen ( m, n, a, x )

!*****************************************************************************80
!
!! MATGEN generates a random matrix A and vector X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    of the matrix.
!
!    Output, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, real ( kind = 8 ) X(N), the vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) x(n)

  seed = 1325
!
! Set the matrix A.
!
  do j = 1, n
    do i = 1, m
      seed = mod ( ( 3125 * seed ), 65536 )
      a(i,j) = ( seed - 32768.0 ) / 16384.0
    end do
  end do
!
!  Set X.
!
  do i = 1, n
    x(i) = i
  end do

  return
end
subroutine mxv_matmul ( m, n, a, x, y )

!*****************************************************************************80
!
!! MXV_MATMUL computes y = A * x, using the FORTRAN90 MATMUL function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) Y(M), the product vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) y(m)

  y(1:m) = matmul ( a(1:m,1:n), x(1:n) )

  return
end
subroutine mxv_matmul_openmp ( m, n, a, x, y )

!*****************************************************************************80
!
!! MXV_MATMUL_OPENMP computes y = A * x, using MATMUL + OpenMP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 July 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) Y(M), the product vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) y(m)

!$omp parallel 

!$omp workshare

  y(1:m) = matmul ( a(1:m,1:n), x(1:n) )

!$omp end workshare

!$omp end parallel

  return
end
subroutine mxv_plain ( m, n, a, x, y )

!*****************************************************************************80
!
!! MXV_PLAIN computes y = A * x, using "plain" code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) Y(M), the product vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) y(m)

  do i = 1, m
    y(i) = 0.0
    do j = 1, n
      y(i) = y(i) + a(i,j) * x(j)
    end do
  end do

  return
end
subroutine mxv_plain_openmp ( m, n, a, x, y )

!*****************************************************************************80
!
!! MXV_PLAIN_OPENMP computes y = A * x, using OpenMP parallel directives.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Barbara Chapman, Gabriele Jost, Ruud vanderPas, David Kuck,
!    Using OpenMP: Portable Shared Memory Parallel Processing,
!    MIT Press, 2007,
!    ISBN13: 978-0262533027,
!    LC: QA76.642.C49.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) Y(M), the product vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) y(m)

!$omp parallel &
!$omp shared ( m, n, a, x, y ) &
!$omp private ( i, j )

!$omp do

  do i = 1, m
    y(i) = 0.0
    do j = 1, n
      y(i) = y(i) + a(i,j) * x(j)
    end do
  end do
!$omp end do

!$omp end parallel

  return
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
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
