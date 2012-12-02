program main

!*****************************************************************************80
!
!! MAIN is the main program for IMAGE_EDGE_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 December 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'IMAGE_EDGE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the IMAGE_EDGE library.'

  call test01 ( 'coins.ascii.pgm', 'coin_edges.ascii.pbm' )

  call test01 ( 'pepper.ascii.pgm', 'pepper_edges.ascii.pbm' )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'IMAGE_EDGE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( input_filename, output_filename )

!*****************************************************************************80
!
!! MAIN is the main program for NEWS_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), allocatable :: e(:,:)
  integer   ( kind = 4 ), allocatable :: g(:,:)
  integer   ( kind = 4 )  g_histo(0:255)
  integer   ( kind = 4 )  g_max
  integer   ( kind = 4 )  i
  character ( len = * )   input_filename
  integer   ( kind = 4 )  input_unit
  integer   ( kind = 4 )  ios
  integer   ( kind = 4 )  m
  integer   ( kind = 4 )  n
  character ( len = * )   output_filename
  integer   ( kind = 4 )  output_unit

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The input file is "' // trim ( input_filename ) // '"'
!
!  Open the input file and read the data.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  call pgma_read_header ( input_unit, m, n, g_max )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of rows =          ', m
  write ( *, '(a,i8)' ) '  Number of columns =       ', n
  write ( *, '(a,i8)' ) '  Maximum pixel intensity = ', g_max

  allocate ( g(1:m,1:n) )

  call pgma_read_data ( input_unit, m, n, g )

  close ( unit = input_unit )

  call i4mat_histogram ( m, n, g, 255, g_histo )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Gray     Count'
  write ( *, '(a)' ) ' '
  do i = 0, 255
    write ( *, '(2x,i3,2x,i8)' ) i, g_histo(i)
  end do

  allocate ( e(1:m,1:n) )

  call news ( m, n, g, e )
!
!  Write the Edge information as a portable BIT map (0/1).
!
  call pbma_write ( output_filename, m, n, e )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Wrote edge information to "' &
    // trim ( output_filename ) // '".'

  deallocate ( e )
  deallocate ( g )

  return
end
