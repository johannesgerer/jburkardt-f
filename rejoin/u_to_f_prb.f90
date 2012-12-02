program u_to_f_prb

!*****************************************************************************80
!
!! MAIN is the main program for U_TO_F_PRB.
!
!  Discussion:
!
!    U_TO_F_PRB converts unformatted DNS files to formatted versions.
!
!    DNS was run on the Cray.  However, the development of the REJOIN
!    routines was carried out on another machine.  The unformatted
!    Cray files were converted to formatted versions, so that the
!    REJOIN routines could be tested on a more convenient platform.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 May 1999
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: npx = 4
  integer ( kind = 4 ), parameter :: npy = 2
  integer ( kind = 4 ), parameter :: nsc = 8
  integer ( kind = 4 ), parameter :: nx_global = 100
  integer ( kind = 4 ), parameter :: ny_global = 100

  integer ( kind = 4 ), parameter :: nx_local = nx_global / npx
  integer ( kind = 4 ), parameter :: ny_local = ny_global / npy

  integer ( kind = 4 ) i_status
  integer ( kind = 4 ) john
  real time

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'U_TO_F_PRB'
  write ( *, '(a)' ) '  Make formatted copies of DNS unformatted files.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data files are associated with'
  write ( *, '(a)' ) '  a decomposition of a rectangular array of data'
  write ( *, '(a)' ) '  into rectangular subarrays associated with'
  write ( *, '(a)' ) '  NPX by NPY processors.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The original region had a grid point shape of '
  write ( *, '(i6,a,i6)' ) nx_global, ' by ', ny_global
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The processor array has the shape:'
  write ( *, '(i6,a,i6)' ) npx, ' by ', npy
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Therefore, the subarrays have shape:'
  write ( *, '(i6,a,i6)' ) nx_local, ' by ', ny_local
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The original files are unformatted.'
  write ( *, '(a)' ) '  The new copies will be formatted.'
!
!  Set JOHN to
!    0, do nothing
!    1, convert L files only.
!    2, convert       T files only.
!    3, convert L and T files.
!
  john = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Requested conversion:'

  if ( john == 0 ) then
    write ( *, '(a)' ) '    NOTHING.'
  else if ( john == 1 ) then
    write ( *, '(a)' ) '    Convert FIELD_L files.'
  else if ( john == 2 ) then
    write ( *, '(a)' ) '    Convert FIELD_T files.'
  else if ( john == 3 ) then
    write ( *, '(a)' ) '    Convert FIELD_L and FIELD_T files.'
  end if
!
!  Convert the FIELD_L_****_ files.
!
  if ( john == 1 .or. john == 3 ) then

    i_status = 1

    time = 0.0E+00

    call u_to_f_all ( i_status, npx, npy, nsc, nx_local, ny_local, time )

  end if
!
!  Convert the FIELD_T_****_tttttt files.
!
  if ( john == 2 .or. john == 3 ) then

    i_status = 0

    time = 0.0E+00

    call u_to_f_all ( i_status, npx, npy, nsc, nx_local, ny_local, time )

  end if
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'U_TO_F_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
