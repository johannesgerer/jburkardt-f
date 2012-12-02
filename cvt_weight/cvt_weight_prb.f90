program main

!*****************************************************************************80
!
!! MAIN is the main program for CVT_WEIGHT_PRB.
!
!  Discussion:
!
!    CVT_WEIGHT_PRB tests the CVT_WEIGHT routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 November 2002
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: cell_num = 3
  integer, parameter :: dim_num = 1
  integer, parameter :: sample_num = 10000

  integer, parameter :: density_num = cell_num + 2

  real box_max(dim_num)
  real box_min(dim_num)
  integer cell
  real cell_gen_coord(dim_num,cell_num)
  real cell_volume_desired(cell_num)
  real cell_volume(cell_num)
  integer count(cell_num)
  integer cvt_it
  integer, parameter :: cvt_it_max = 5
  real density_coord(density_num-1)
  real density_max
  real density_value(density_num)
  real density_value_max
  integer i
! integer, parameter :: density_it_max = 10
  real, parameter :: region_volume = 10.0
  real sample_coord(dim_num,sample_num)
  integer seed
  real x

  call timestamp()

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_WEIGHT_PRB'
  write ( *, '(a)' ) '  Test the CVT_WEIGHT library.'

  seed = 1952
  call random_initialize ( seed )

  box_min(1) = 0.0
  box_max(1) = 10.0

  cell_volume_desired(1:cell_num) = &
    (/ 2.0, 3.0, 5.0 /)

  cell_volume_desired(1:cell_num) = cell_volume_desired(1:cell_num) &
    * region_volume / sum ( cell_volume_desired(1:cell_num) )
!
!  Set the initial generators.
!  SORT THEM.
!
  cell_gen_coord(1,1:cell_num) = &
    (/ 1.0, 4.0, 8.0 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cell, Generator'
  write ( *, '(a)' ) ' '
  do cell = 1, cell_num
    write ( *, '(2x,i3,3f10.4)' ) cell, cell_gen_coord(1:dim_num,cell)
  end do
!
!  Determine the cell volumes.
!
  call cell_volume_computation ( dim_num, box_min, box_max, cell_num, &
    cell_gen_coord, sample_num, cell_volume, region_volume )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cell, CELL VOLUME, CELL_VOLUME_DESIRED'
  write ( *, '(a)' ) ' '
  do cell = 1, cell_num
    write ( *, '(2x,i4,3f10.4)' ) &
      cell, cell_volume(cell), cell_volume_desired(cell)
  end do
!
!  Set the density function.
!
  call rvec_even ( box_min(1), box_max(1), density_num-1, density_coord )

  density_value(1) = 0.0
  do cell = 1, cell_num
    density_value(cell+1) = 1.0 / cell_volume_desired(cell)**(1.0/3.0)
  end do
  density_value(cell_num+2) = 0.0

  density_value(2:4) = (/ 10.0, 1.0, 1.0 /)

  density_coord(1) = box_min(1)
  do cell = 1, cell_num-1
    density_coord(cell+1) = 0.5 * &
      ( cell_gen_coord(1,cell) + cell_gen_coord(1,cell+1) )
  end do
  density_coord(cell_num+1) = box_max(1)

  density_max = maxval ( density_value(1:density_num) )
  density_value(1:density_num) = density_value(1:density_num) / density_max
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Piecewise constant density function:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Breaks, Values'
  write ( *, '(a)' ) ' '
  do i = 1, density_num
    write ( *, '(14x,g14.6)' ) density_value(i)
    if ( i < density_num ) then
      write ( *, '(g14.6)' ) density_coord(i)
    end if
  end do
!
!  Find the CVT, given the current density function.
!
  do cvt_it = 1, cvt_it_max

    call cvt_dense_iteration ( density_num, dim_num, box_min, box_max, &
      cell_num, cell_gen_coord, density_value, density_coord, sample_num )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Cell, Cell_Generator'
    write ( *, '(a)' ) ' '
    do cell = 1, cell_num
      write ( *, '(2x,i3,3f10.4)' ) cell, cell_gen_coord(1:dim_num,cell)
    end do
!
!  Determine the cell volumes.
!
    call cell_volume_computation ( dim_num, box_min, box_max, cell_num, &
      cell_gen_coord, sample_num, cell_volume, region_volume )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Cell, CELL VOLUME, CELL_VOLUME_DESIRED'
    write ( *, '(a)' ) ' '
    do cell = 1, cell_num
      write ( *, '(2x,i4,3f10.4)' ) &
        cell, cell_volume(cell), cell_volume_desired(cell)
    end do

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CVT_WEIGHT_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
