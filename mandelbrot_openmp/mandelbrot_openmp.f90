program main

!*****************************************************************************80
!
!! MAIN is the main program for MANDELBROT_OPENMP.
!
!  Discussion:
!
!    MANDELBROT_OPENMP computes an image of the Mandelbrot set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Local Parameters:
!
!    Local, integer COUNT_MAX, the maximum number of iterations taken
!    for a particular pixel.
!
  use omp_lib

  implicit none

  integer ( kind = 4 ), parameter :: m = 500
  integer ( kind = 4 ), parameter :: n = 500

  integer ( kind = 4 ) b(m,n)
  integer ( kind = 4 ) c
  integer ( kind = 4 ) c_max
  integer ( kind = 4 ) count(m,n)
  integer ( kind = 4 ) :: count_max = 2000
  integer ( kind = 4 ) g(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  character ( len = 80 ) output_filename
  integer ( kind = 4 ) output_unit
  integer ( kind = 4 ) r(m,n)
  real ( kind = 8 ) wtime
  real ( kind = 8 ) wtime_total
  real ( kind = 8 ) :: x_max =   1.25D+00
  real ( kind = 8 ) :: x_min = - 2.25D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) :: y_max =   1.75D+00
  real ( kind = 8 ) :: y_min = - 1.75D+00
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MANDELBROT_OPENMP'
  write ( *, '(a)' ) '  FORTRAN90/OpenMP version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Create an ASCII PPM image of the Mandelbrot set.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For each point C = X + i*Y'
  write ( *, '(a,g14.6,a,g14.6,a)' ) '  with X range [', x_min, ',', x_max, ']'
  write ( *, '(a,g14.6,a,g14.6,a)' ) '  and  Y range [', y_min, ',', y_max, ']'
  write ( *, '(a,i8,a)' ) '  carry out ', count_max, ' iterations of the map'
  write ( *, '(a)' ) '  Z(n+1) = Z(n)^2 + C.'
  write ( *, '(a)' ) '  If the iterates stay bounded (norm less than 2)'
  write ( *, '(a)' ) '  then C is taken to be a member of the set.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  An ASCII PPM image of the set is created using'
  write ( *, '(a,i8,a)' ) '    M = ', m, ' pixels in the X direction and'
  write ( *, '(a,i8,a)' ) '    N = ', n, ' pixels in the Y direction.'

  wtime = omp_get_wtime ( )
!
!  Carry out the iteration for each pixel, determining COUNT.
!
!$omp parallel &
!$omp shared ( b, count, count_max, g, r, x_max, x_min, y_max, y_min ) &
!$omp private ( i, j, k, x, x1, x2, y, y1, y2 )

!$omp do

  do i = 1, m

    do j = 1, n

      x = ( real (     j - 1, kind = 8 ) * x_max   &
          + real ( m - j,     kind = 8 ) * x_min ) &
          / real ( m     - 1, kind = 8 )

      y = ( real (     i - 1, kind = 8 ) * y_max   &
          + real ( n - i,     kind = 8 ) * y_min ) &
          / real ( n     - 1, kind = 8 )

      count(i,j) = 0

      x1 = x
      y1 = y

      do k = 1, count_max

        x2 = x1 * x1 - y1 * y1 + x
        y2 = 2 * x1 * y1 + y

        if ( x2 < -2.0D+00 .or. 2.0D+00 < x2 .or. &
             y2 < -2.0D+00 .or. 2.0D+00 < y2 ) then
          count(i,j) = k
          exit
        end if

        x1 = x2
        y1 = y2

      end do

      if ( mod ( count(i,j), 2 ) == 1 ) then
        r(i,j) = 255
        g(i,j) = 255
        b(i,j) = 255
      else
        c = int ( 255.0D+00 * sqrt ( sqrt ( sqrt ( &
          ( real ( count(i,j), kind = 8 ) / real ( count_max, kind = 8 ) ) ) ) ) )
        r(i,j) = 3 * c / 5
        g(i,j) = 3 * c / 5
        b(i,j) = c
      end if

    end do

  end do
!$omp end do

!$omp end parallel

  wtime = omp_get_wtime ( ) - wtime
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Time = ', wtime
!
!  Write data to an ASCII PPM file.
!
  output_unit = 10
  output_filename = 'mandelbrot.ppm'

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  write ( output_unit, '(a2)' ) 'P3'
  write ( output_unit, '(i5,2x,i5)' ) n, m
  write ( output_unit, '(i3)' ) 255
  do i = 1, m
    do jlo = 1, n, 4
      jhi = min ( jlo + 3, n )
      write ( output_unit, '(12i5)' ) ( r(i,j), g(i,j), b(i,j), j = jlo, jhi )
    end do
  end do

  close ( unit = output_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Graphics data written to "' // trim ( output_filename ) // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MANDELBROT_OPENMP'
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
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
