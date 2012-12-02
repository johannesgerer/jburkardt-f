program main

!*****************************************************************************80
!
!! MAIN is the main program for HEX_GRID_DATASET.
!
!  Discussion:
!
!    HEX_GRID_DATASET generates a hexagonal grid dataset and writes it to a file.
!
!    This program is meant to be used interactively.  It's also
!    possible to prepare a simple input file beforehand and use it
!    in batch mode.
!
!    The program requests input values from the user:
!
!    * X1, Y1, the lower left corner of a bounding box.
!    * X2, Y2, the upper right corner of a bounding box.
!    * NODES_PER_LAYER, the number of nodes per layer.
!
!    The program will now define the dataset, and write it to a file.
!
!    The program will now request that you type "Y" if you want to 
!    set up another dataset.  Otherwise the program terminates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: dim_num = 2

  real ( kind = 8 ) box(dim_num,2)
  character c
  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer ios
  integer layers
  integer n
  integer nodes_per_layer
  character ( len = 80 ) :: output_file_name
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: P

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HEX_GRID_DATASET'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Generate a hexagonal grid dataset.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This program is meant to be used interactively.'
  write ( *, '(a)' ) '  It is also possible to prepare a simple input '
  write ( *, '(a)' ) '  file beforehand and use it in batch mode.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program requests input values from the user:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  * X1, Y1, the lower left corner of the region,'
  write ( *, '(a)' ) '  * X2, Y2, the upper right corner of the region,'
  write ( *, '(a)' ) '  * NODES_PER_LAYER, the number of nodes per layer,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After the dataset of nodes is computed, it is'
  write ( *, '(a)' ) '  written to a file, and another dataset may be made.'
  write ( *, '(a)' ) ' '

  do

    write ( *, '(a)' ) '  *'
    write ( *, '(a)' ) ' *'
    write ( *, '(a)' ) '*  Ready to generate a new dataset:'
    write ( *, '(a)' ) ' *'
    write ( *, '(a)' ) '  *'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter X1, Y1, the lower left corner of the region:'
    write ( *, '(a)' ) '  (Try ''0, 0'' if you do not have a preference.)'
    write ( *, '(a)' ) '  (Enter "*" or "QUIT" to terminate execution.)'

    read ( *, *, iostat = ios ) box(1:dim_num,1)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HEX_GRID_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for X1, Y1.'
      exit
    end if

    write ( *, '(a,2g14.6)' ) '  User input X1, Y1 = ', box(1:dim_num,1)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter X2, Y2, the upper right corner of the region:'
    write ( *, '(a)' ) '  (Try ''10, 10'' if you do not have a preference.)'
    write ( *, '(a)' ) '  (Enter "*" or "QUIT" to terminate execution.)'

    read ( *, *, iostat = ios ) box(1:dim_num,2)

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HEX_GRID_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for X2, Y2.'
      exit
    end if

    write ( *, '(a,2g14.6)' ) '  User input X2, Y2 = ', box(1:dim_num,2)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '  Enter NODES_PER_LAYER, the number of nodes in a layer.'
    write ( *, '(a)' ) '  (Try ''10'' if you do not have a preference.)'
    write ( *, '(a)' ) '  (1 or any smaller value terminates execution).'

    read ( *, *, iostat = ios ) nodes_per_layer

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HEX_GRID_DATASET - Warning!'
      write ( *, '(a)' ) '  Terminating abnormally because of an I/O error'
      write ( *, '(a)' ) '  while expecting input for NODES_PER_LAYER.'
      exit
    end if

    write ( *, '(a,i12)' ) '  User input NODES_PER_LAYER = ', nodes_per_layer

    if ( nodes_per_layer <= 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HEX_GRID_DATASET'
      write ( *, '(a,i12)' ) &
        '  The input value of NODES_PER_LAYER = ', nodes_per_layer
      write ( *, '(a)' ) '  is interpreted as a request for termination.'
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
    end if

    call hex_grid_layers ( nodes_per_layer, box, layers )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  The number of layers will be ', layers

    call hex_grid_h ( nodes_per_layer, box, hx, hy )

    write ( *, '(a,g14.6)' ) '  The X spacing will be ', hx
    write ( *, '(a,g14.6)' ) '  The Y spacing will be ', hy

    call hex_grid_n ( nodes_per_layer, box, n )

    write ( *, '(a,i12)' ) '  The number of nodes ', n

    allocate ( p(1:2,1:n) )

    call hex_grid_points ( nodes_per_layer, layers, n, box, p )

    write ( output_file_name, '(a,i3,a,i3,a,i6,a)' ) &
      'hex_grid_', nodes_per_layer, '_', layers, '_', n, '.txt'

    call s_blank_delete ( output_file_name )

    call hex_grid_write ( n, nodes_per_layer, layers, hx, hy, box, &
      p, output_file_name )

    deallocate ( p )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The data was written to the file "' &
       // trim ( output_file_name ) // '".'

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter "Y" if you want to define another dataset.'

    read ( *, '(a)' ) c

    if ( c /= 'y' .and. c /= 'Y' ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HEX_GRID_DATASET:'
      write ( *, '(a)' ) '  Normal end of execution.'
      exit
    end if

  end do

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
