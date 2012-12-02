program main

!*****************************************************************************80
!
!! MAIN is the main program for IMAGE_DENOISE_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'IMAGE_DENOISE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the IMAGE_DENOISE library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'IMAGE_DENOISE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests GRAY_MEDIAN_NEWS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: g
  integer ( kind = 4 ) g_max
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: g2
  character ( len = 80 ) :: input_filename = "glassware_noisy.ascii.pgm"
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  character ( len = 80 ) :: output_filename = "glassware_median_news.ascii.pgm"

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  GRAY_MEDIAN_NEWS uses a NEWS median filter'
  write ( *, '(a)' ) '  on a noisy grayscale image.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The input file is "' // trim ( input_filename ) // '".'
!
!  Open the input file and read the data.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call pgma_read_header ( input_unit, m, n, g_max )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Number of rows =          ', m
  write ( *, '(a,i4)' ) '  Number of columns =       ', n
  write ( *, '(a,i4)' ) '  Maximum pixel intensity = ', g_max

  allocate ( g(1:m,1:n) )

  call pgma_read_data ( input_unit, m, n, g );

  close ( unit = input_unit )

  allocate ( g2(1:m,1:n) )

  call gray_median_news ( m, n, g, g2 )
!
!  Write the denoised images.
!
  call pgma_write ( output_filename, m, n, g2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Wrote denoised image to "' // trim ( output_filename ) // '".'

  deallocate ( g )
  deallocate ( g2 )

  return
end
