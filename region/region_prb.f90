program main

!*****************************************************************************80
!
!! MAIN is the main program for REGION_PRB.
!
!  Discussion:
!
!    REGION_PRB reads voxel data, processes it, and writes it back out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxlist = 8000
  integer ( kind = 4 ), parameter :: region_max = 20
  integer ( kind = 4 ), parameter :: nx = 64
  integer ( kind = 4 ), parameter :: ny = 64
  integer ( kind = 4 ), parameter :: nz = 26

  real ( kind = 8 ) ave_pos
  integer ( kind = 4 ) c(3)
  integer ( kind = 4 ) center(4,region_max)
  character ( len = 80 ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iregion
  integer ( kind = 4 ) i4voxel(nx,ny,nz)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) list(maxlist)
  integer ( kind = 4 ) nelements
  integer ( kind = 4 ) nlist
  integer ( kind = 4 ) num_pos
  real ( kind = 8 ) r8voxel(nx,ny,nz)
  integer ( kind = 4 ) region_num
  integer ( kind = 4 ) total
  integer ( kind = 4 ) thresh

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'REGION_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the REGION library.'

  call test01 ( )
!
!  Read the MRI data.
!
  filename = 'roi.ascii'

  call i4voxel_read ( nx, ny, nz, i4voxel, filename )
!
!  Find the average nonzero value.
!
  call i4voxel_sum ( nx, ny, nz, i4voxel, total )

  call i4voxel_count_positive ( nx, ny, nz, i4voxel, num_pos )

  ave_pos = real ( total ) / real ( num_pos )

! call i4voxel_plot ( nx, ny, nz, i4voxel )
!
!  Zero out all voxels below a given threshold value.
!
  thresh = nint ( 0.75D+00 * ave_pos )

  call i4voxel_thresh ( nx, ny, nz, i4voxel, thresh )
!
!  Thicken the voxels.
!
  call i4voxel_thicken ( nx, ny, nz, i4voxel )
!
!  Now find the regions.
!
  call i4voxel_to_region ( nx, ny, nz, i4voxel, list, maxlist, nlist, &
    region_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) &
    '  The voxels are grouped into ', region_num, ' regions.'
  write ( *, '(a)' ) ' '

  if ( .false. ) then
    write ( *, '(a)' ) 'The nonzero array elements are:'
    write ( *, '(a)' ) ' '

    do i = 1, nx
      do j = 1, ny
        do k = 1, nz
          l = i4voxel(i,j,k)
          if ( l /= 0 ) then
            write ( *, '(4i6)' ) i, j, k, l
          end if
        end do
      end do
    end do
  end if

  if ( maxlist < nlist ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The stack-based list of regions is unusable.'

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The stack-based list of regions:'
    write ( *, '(a)' ) ' '

    iregion = region_num + 1

    do

      iregion = iregion - 1

      if ( nlist <= 0 ) then
        exit
      end if

      nelements = list(nlist)
      nlist = nlist - 1

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a,i8,a)' ) &
        '  Region ', iregion, ' includes ', nelements, ' voxels:'
      write ( *, '(a)' ) ' '

      do l = 1, nelements

        k = list(nlist)
        nlist = nlist - 1
        j = list(nlist)
        nlist = nlist - 1
        i = list(nlist)
        nlist = nlist - 1
        write ( *, '(3i6)' ) i, j, k

      end do

    end do

  end if
!
!  Make a typewriter plot of the regions in Z slices.
!
  if ( .true. ) then
    call i4voxel_plot2 ( nx, ny, nz, i4voxel )
  end if
!
!  Compute center of mass of each region.
!
  call region_center ( nx, ny, nz, i4voxel, region_max, region_num, center )
!
!  Save the center of mass of the central region.
!
  iregion = 4
  do i = 1, 3
    c(i) = center(i,iregion)
  end do
!
!  Zero out region #4, slide other region numbers down 1.
!
  iregion = 4
  call region_blank ( nx, ny, nz, i4voxel, center, iregion, region_max, &
    region_num )
!
!  Write out an ASCII MRI file containing the voxels, marked by region only.
!
  filename = 'roi_region.ascii'

  call i4voxel_write ( nx, ny, nz, i4voxel, filename )
!
!  Build a new array, R8VOXEL.
!
!  For each region IREGION, construct the line from the center of
!  the central region through its own center.
!
!  For each voxel in region IREGION, move outward along this line
!  until you reach the boundary.  Add PERCENT to each voxel in ARRAY
!  through which the voxel passes, where PERCENT is the percentage of
!  the total region volume represented by one voxel.  And if this
!  is the first time any voxel has passed through this voxel,
!  add 100 * IREGION to it.
!
  call transport ( nx, ny, nz, r8voxel, c, center, i4voxel, region_max )
!
!  Copy the real ARRAY into the integer IVOXEL.
!
  call r8voxel_to_i4voxel ( nx, ny, nz, r8voxel, i4voxel )
!
!  Do a 0/nonzero plot of the data.
!
  if ( .false. ) then
    call i4voxel_plot3 ( nx, ny, nz, i4voxel )
  end if
!
!  Write data to a file.
!
  filename = 'roi_colors.ascii'

  call i4voxel_write ( nx, ny, nz, i4voxel, filename )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'REGION_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates I4VOXEL_TO_OBJ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 December 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nx = 64
  integer ( kind = 4 ), parameter :: ny = 64
  integer ( kind = 4 ), parameter :: nz = 26

  character ( len = 80 ) ascii_filename
  integer ( kind = 4 ), allocatable :: i4voxel(:,:,:)
  character ( len = 80 ) obj_filename

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Read an ASCII MRI data file.'
  write ( *, '(a)' ) '  Write a corresponding OBJ 3D graphics file.'
!
!  Read the MRI data.
!
  allocate ( i4voxel(1:nx,1:ny,1:nz) )

  ascii_filename = 'roi.ascii'
  call i4voxel_read ( nx, ny, nz, i4voxel, ascii_filename )
!
!  Make an OBJ file of the data.
!
  obj_filename = 'roi.obj'
  call i4voxel_to_obj ( nx, ny, nz, i4voxel, obj_filename )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created the OBJ graphics file "' &
    // trim ( obj_filename ) // '".'

  deallocate ( i4voxel )

  return
end
