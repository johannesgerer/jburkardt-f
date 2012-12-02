program main

!*****************************************************************************80
!
!! MAIN is the main program for HB_IO_PRB.
!
!  Discussion:
!
!    HB_IO_PRB runs the HB_IO tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HB_IO_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the HB_IO library.'

  call test01 ( 'rua_32.txt' )
  call test01 ( 'rse_5.txt' )
  call test02 ( )
  call test03 ( 'rua_32.txt' )
  call test03 ( 'rse_5.txt' )
  call test04 ( )
  call test05 ( 'rua_32.txt' )
  call test05 ( 'rse_5.txt' )
  call test06 ( )
  call test07 ( 'rua_32.txt' )
  call test08 ( )
  call test09 ( 'rua_32.txt' )
  call test10 ( )
  call test11 ( )
  call test12 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HB_IO_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( input_file )

!*****************************************************************************80
!
!! TEST01 tests HB_HEADER_READ;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer indcrd
  character ( len = 16 ) indfmt
  character ( len = * ) :: input_file
  integer input_unit
  integer ios
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ncol
  integer neltvl
  integer nnzero
  integer nrhs
  integer nrhsix
  integer nrow
  integer ptrcrd
  character ( len = 16 ) ptrfmt
  integer rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  character ( len = 72 ) title
  integer totcrd
  integer valcrd
  character ( len = 20 ) valfmt

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  HB_HEADER_READ reads the header of an HB file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reading the file "' // trim ( input_file ) // '".'

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a)' ) '  Error opening the file.'
    return
  end if

  call hb_header_read ( input_unit, title, key, totcrd, ptrcrd, indcrd, &
    valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, &
    valfmt, rhsfmt, rhstyp, nrhs, nrhsix )

  close ( unit = input_unit )
!
!  Print out the  header information.
!
  call hb_header_print ( title, key, totcrd, ptrcrd, indcrd, valcrd, &
    rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, valfmt, &
    rhsfmt, rhstyp, nrhs, nrhsix )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests HB_HEADER_WRITE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer                :: indcrd = 8
  character ( len = 16 ) :: indfmt = '(16I5)'
  integer ios
  character ( len = 8 )  :: key = 'RUA_32'
  character ( len = 3 )  :: mxtype = 'PUA'
  integer                :: ncol = 32
  integer                :: neltvl = 0
  integer                :: nnzero = 126
  integer                :: nrhs = 0
  integer                :: nrhsix = 0
  integer                :: nrow = 32
  character ( len = 80 ) :: output_file = 'rua_32_header.txt'
  integer output_unit
  integer                :: ptrcrd = 3
  character ( len = 16 ) :: ptrfmt = '(16I5)'
  integer                :: rhscrd = 0
  character ( len = 20 ) :: rhsfmt = ' '
  character ( len = 3 )  :: rhstyp = '   '
  character ( len = 72 ) :: title = &
    '1Real unsymmetric assembled matrix based on IBM32'
  integer                :: totcrd = 11
  integer                :: valcrd = 0
  character ( len = 20 ) :: valfmt = ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  HB_HEADER_WRITE writes the header of an HB file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Writing the file "' // trim ( output_file ) // '".'

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Fatal error!'
    write ( *, '(a)' ) '  Error opening the file.'
    return
  end if

  call hb_header_write ( output_unit, title, key, totcrd, ptrcrd, indcrd, &
    valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, &
    valfmt, rhsfmt, rhstyp, nrhs, nrhsix )

  close ( unit = output_unit )

  return
end
subroutine test03 ( input_file )

!*****************************************************************************80
!
!! TEST03 tests HB_STRUCTURE_READ;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, allocatable, dimension ( : ) :: colptr
  integer indcrd
  character ( len = 16 ) indfmt
  character ( len = * ) :: input_file
  integer input_unit
  integer ios
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ncol
  integer neltvl
  integer nnzero
  integer nrhs
  integer nrhsix
  integer nrow
  integer ptrcrd
  character ( len = 16 ) ptrfmt
  integer rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  integer, allocatable, dimension ( : ) :: rowind
  character ( len = 72 ) title
  integer totcrd
  integer valcrd
  character ( len = 20 ) valfmt

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  HB_STRUCTURE_READ reads the structure of an HB file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reading the file "' // trim ( input_file ) // '".'

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST03 - Fatal error!'
    write ( *, '(a)' ) '  Error opening the file.'
    return
  end if

  write ( *, '(a)' ) '  Reading the header.'

  call hb_header_read ( input_unit, title, key, totcrd, ptrcrd, indcrd, &
    valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, &
    valfmt, rhsfmt, rhstyp, nrhs, nrhsix )

  allocate ( colptr(ncol+1) )

  if ( mxtype(3:3) == 'A' ) then
    allocate ( rowind(nnzero) )
  else if ( mxtype(3:3) == 'E' ) then
    allocate ( rowind(neltvl) )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST03 - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of MXTYPE(3:3) = ' // mxtype(3:3)
    stop
  end if

  write ( *, '(a)' ) '  Reading the structure.'

  call hb_structure_read ( input_unit, ncol, mxtype, nnzero, neltvl, &
    ptrcrd, ptrfmt, indcrd, indfmt, colptr, rowind )

  close ( unit = input_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) '  KEY =    "' // key // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  NROW =   ', nrow
  write ( *, '(a,i6)' ) '  NCOL =   ', ncol
  write ( *, '(a,i6)' ) '  NNZERO = ', nnzero
  write ( *, '(a,i6)' ) '  NELTVL = ', neltvl

  call hb_structure_print ( ncol, mxtype, nnzero, neltvl, colptr, rowind )

  deallocate ( colptr )
  deallocate ( rowind )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests HB_STRUCTURE_WRITE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ncol = 32
  integer, parameter :: neltvl = 0
  integer, parameter :: nnzero = 126

  integer, dimension ( ncol + 1 ) :: colptr = (/ &
      1,   7,  12,  18,  22,  26,  29,  34,  39,  46, &
     53,  58,  61,  63,  65,  68,  71,  74,  79,  82, &
     85,  88,  90,  94,  97, 102, 106, 110, 112, 117, &
    121, 124, 127 /)
  character ( len = 16 ) :: indfmt = '(16I5)'
  integer ios
  character ( len = 3 ) :: mxtype = 'RUA'
  character ( len = 80 ) :: output_file = 'rua_32_structure.txt'
  integer output_unit
  character ( len = 16 ) :: ptrfmt = '(16I5)'
  integer, dimension ( nnzero ) :: rowind = (/ &
    1,    2,    3,    4,    7,   26,    1,    2,    9,   21, &
   28,    2,    3,    6,    8,    9,   29,    3,    4,    5, &
   12,    3,    5,   23,   27,    1,    6,   16,    3,    7, &
   14,   21,   31,    1,    8,   12,   17,   27,    7,    9, &
   10,   13,   19,   23,   27,    1,   10,   11,   21,   23, &
   25,   27,    2,   11,   15,   18,   29,    6,   12,   24, &
   11,   13,    3,   14,    2,   15,   20,    4,   16,   22, &
    4,   16,   17,    6,   10,   18,   20,   30,    1,   19, &
   26,    8,   16,   20,    3,   21,   32,   11,   22,    2, &
   17,   21,   23,   12,   24,   26,    6,   15,   18,   24, &
   25,   13,   18,   22,   26,    5,   24,   26,   27,    9, &
   28,    3,    5,   27,   29,   32,   12,   17,   23,   30, &
   13,   14,   31,   24,   28,   32 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) &
    '  HB_STRUCTURE_WRITE writes the structure of an HB file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Writing the file "' // trim ( output_file ) // '".'

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST04 - Fatal error!'
    write ( *, '(a)' ) '  Error opening the file.'
    return
  end if

  call hb_structure_write ( output_unit, ncol, mxtype, nnzero, neltvl, &
    ptrfmt, indfmt, colptr, rowind )

  close ( unit = output_unit )

  return
end
subroutine test05 ( input_file )

!*****************************************************************************80
!
!! TEST05 tests HB_VALUES_READ;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, allocatable, dimension ( : ) :: colptr
  integer indcrd
  character ( len = 16 ) indfmt
  character ( len = * ) :: input_file
  integer input_unit
  integer ios
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ncol
  integer neltvl
  integer nnzero
  integer nrhs
  integer nrhsix
  integer nrow
  integer ptrcrd
  character ( len = 16 ) ptrfmt
  integer rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  integer, allocatable, dimension ( : ) :: rowind
  character ( len = 72 ) title
  integer totcrd
  integer valcrd
  character ( len = 20 ) valfmt
  real ( kind = 8 ), allocatable, dimension ( : ) :: values

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  HB_VALUES_READ reads the values of an HB file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reading the file "' // trim ( input_file ) // '".'

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05 - Fatal error!'
    write ( *, '(a)' ) '  Error opening the file.'
    return
  end if

  write ( *, '(a)' ) '  Reading the header.'

  call hb_header_read ( input_unit, title, key, totcrd, ptrcrd, indcrd, &
    valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, &
    valfmt, rhsfmt, rhstyp, nrhs, nrhsix )

  allocate ( colptr(ncol+1) )

  if ( mxtype(3:3) == 'A' ) then
    allocate ( rowind(nnzero) )
  else if ( mxtype(3:3) == 'E' ) then
    allocate ( rowind(neltvl) )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05 - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of MXTYPE(3:3).'
    stop
  end if

  write ( *, '(a)' ) '  Reading the structure.'

  call hb_structure_read ( input_unit, ncol, mxtype, nnzero, neltvl, &
    ptrcrd, ptrfmt, indcrd, indfmt, colptr, rowind )

  if ( mxtype(3:3) == 'A' ) then
    allocate ( values(nnzero) )
  else if ( mxtype(3:3) == 'E' ) then
    allocate ( values(neltvl) )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05 - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of MXTYPE(3:3) = ' // mxtype(3:3)
    stop
  end if

  write ( *, '(a)' ) '  Reading the values.'

  call hb_values_read ( input_unit, valcrd, mxtype, nnzero, neltvl, &
    valfmt, values )

  close ( unit = input_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) '  KEY =    "' // key // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  NROW =   ', nrow
  write ( *, '(a,i6)' ) '  NCOL =   ', ncol
  write ( *, '(a,i6)' ) '  NNZERO = ', nnzero
  write ( *, '(a,i6)' ) '  NELTVL = ', neltvl

  call hb_values_print ( ncol, colptr, mxtype, nnzero, neltvl, values )

  deallocate ( colptr )
  deallocate ( rowind )
  deallocate ( values )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests HB_VALUES_WRITE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: neltvl = 0
  integer, parameter :: nnzero = 126

  integer ios
  character ( len = 3 ) :: mxtype = 'RUA'
  character ( len = 80 ) :: output_file = 'rua_32_values.txt'
  integer output_unit
  integer :: valcrd = 13
  character ( len = 20 ) :: valfmt = '(10F7.1)'
  real ( kind = 8 ), dimension ( nnzero ) :: values = (/ &
  101.0,  102.0,  103.0,  104.0,  107.0, &
  126.0,  201.0,  202.0,  209.0,  221.0, &
  228.0,  302.0,  303.0,  306.0,  308.0, &
  309.0,  329.0,  403.0,  404.0,  405.0, &
  412.0,  503.0,  505.0,  523.0,  527.0, &
  601.0,  606.0,  616.0,  703.0,  707.0, &
  714.0,  721.0,  731.0,  801.0,  808.0, &
  812.0,  817.0,  827.0,  907.0,  909.0, &
  910.0,  913.0,  919.0,  923.0,  927.0, &
 1001.0, 1010.0, 1011.0, 1021.0, 1023.0, &
 1025.0, 1027.0, 1102.0, 1111.0, 1115.0, &
 1118.0, 1129.0, 1206.0, 1212.0, 1224.0, &
 1311.0, 1313.0, 1403.0, 1414.0, 1502.0, &
 1515.0, 1520.0, 1604.0, 1616.0, 1622.0, &
 1704.0, 1716.0, 1717.0, 1806.0, 1810.0, &
 1818.0, 1820.0, 1830.0, 1901.0, 1919.0, &
 1926.0, 2008.0, 2016.0, 2020.0, 2103.0, &
 2121.0, 2132.0, 2211.0, 2222.0, 2302.0, &
 2317.0, 2321.0, 2323.0, 2412.0, 2424.0, &
 2426.0, 2506.0, 2515.0, 2518.0, 2524.0, &
 2525.0, 2613.0, 2618.0, 2622.0, 2626.0, &
 2705.0, 2724.0, 2726.0, 2727.0, 2809.0, &
 2828.0, 2903.0, 2905.0, 2927.0, 2929.0, &
 2932.0, 3012.0, 3017.0, 3023.0, 3030.0, &
 3113.0, 3114.0, 3131.0, 3224.0, 3228.0, &
 3232.0 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) &
    '  HB_VALUES_WRITE writes the values of an HB file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Writing the file "' // trim ( output_file ) // '".'

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST06 - Fatal error!'
    write ( *, '(a)' ) '  Error opening the file.'
    return
  end if

  call hb_values_write ( output_unit, valcrd, mxtype, nnzero, neltvl, &
    valfmt, values )

  close ( unit = output_unit )

  return
end
subroutine test07 ( input_file )

!*****************************************************************************80
!
!! TEST07 tests HB_RHS_READ, HB_GUESS_READ, TEST07 tests HB_EXACT_READ;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, allocatable, dimension ( : ) :: colptr
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: exact
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: guess
  integer indcrd
  character ( len = 16 ) indfmt
  character ( len = * ) :: input_file
  integer input_unit
  integer ios
  character ( len = 8 ) key
  character ( len = 3 ) mxtype
  integer ncol
  integer neltvl
  integer nnzero
  integer nrhs
  integer nrhsix
  integer nrow
  integer ptrcrd
  character ( len = 16 ) ptrfmt
  integer rhscrd
  character ( len = 20 ) rhsfmt
  integer, allocatable, dimension ( : ) :: rhsind
  integer, allocatable, dimension ( : ) :: rhsptr
  character ( len = 3 ) rhstyp
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: rhsval
  real ( kind = 8 ), allocatable, dimension ( : )    :: rhsvec
  integer, allocatable, dimension ( : ) :: rowind
  character ( len = 72 ) title
  integer totcrd
  integer valcrd
  character ( len = 20 ) valfmt
  real ( kind = 8 ), allocatable, dimension ( : ) :: values

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  HB_RHS_READ reads right hand sides from an HB file.'
  write ( *, '(a)' ) '  HB_GUESS_READ reads starting guesses from an HB file.'
  write ( *, '(a)' ) '  HB_EXACT_READ reads exact solutions from an HB file.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reading the file "' // trim ( input_file ) // '".'

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST07 - Fatal error!'
    write ( *, '(a)' ) '  Error opening the file.'
    return
  end if

  write ( *, '(a)' ) '  Reading the header.'

  call hb_header_read ( input_unit, title, key, totcrd, ptrcrd, indcrd, &
    valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, &
    valfmt, rhsfmt, rhstyp, nrhs, nrhsix )

  allocate ( colptr(ncol+1) )

  if ( mxtype(3:3) == 'A' ) then
    allocate ( rowind(nnzero) )
  else if ( mxtype(3:3) == 'E' ) then
    allocate ( rowind(neltvl) )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST07 - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of MXTYPE(3:3).'
    stop
  end if

  write ( *, '(a)' ) '  Reading the structure.'

  call hb_structure_read ( input_unit, ncol, mxtype, nnzero, neltvl, &
    ptrcrd, ptrfmt, indcrd, indfmt, colptr, rowind )

  if ( mxtype(3:3) == 'A' ) then
    allocate ( values(nnzero) )
  else if ( mxtype(3:3) == 'E' ) then
    allocate ( values(neltvl) )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST07 - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of MXTYPE(3:3) = ' // mxtype(3:3)
    stop
  end if

  write ( *, '(a)' ) '  Reading the values.'

  call hb_values_read ( input_unit, valcrd, mxtype, nnzero, neltvl, &
    valfmt, values )

  if ( 0 < rhscrd ) then
!
!  Read the right hand sides.
!
    if ( rhstyp(1:1) == 'F' ) then
      allocate ( rhsval(nrow,nrhs) )
    else if ( rhstyp(1:1) == 'M' ) then
      if ( mxtype(3:3) == 'A' ) then
        allocate ( rhsptr(1:nrhs+1) )
        allocate ( rhsind(1:nrhsix) )
        allocate ( rhsvec(1:nrhsix) )
      else if ( mxtype(3:3) == 'E' ) then
        allocate ( rhsval(nnzero,nrhs) )
      end if
    end if

    write ( *, '(a)' ) '  Reading the right hand side.'

    call hb_rhs_read ( input_unit, nrow, nnzero, nrhs, nrhsix, &
      rhscrd, ptrfmt, indfmt, rhsfmt, mxtype, rhstyp, rhsval, &
      rhsind, rhsptr, rhsvec )

    if ( rhstyp(1:1) == 'F' ) then

      call r8mat_print_some ( nrow, nrhs, rhsval, 1, 1, 5, 5, &
        '  Part of RHS array' )

    else if ( rhstyp(1:1) == 'M' .and. mxtype(3:3) == 'A' ) then

      call i4vec_print_some ( nrhs+1, rhsptr, 1, 10, '  Part of RHSPTR' )
      call i4vec_print_some ( nrhsix, rhsind, 1, 10, '  Part of RHSIND' )
      call r8vec_print_some ( nrhsix, rhsvec, 1, 10, '  Part of RHSVEC' )

    else if ( rhstyp(1:1) == 'M' .and. mxtype(3:3) == 'E' ) then

      call r8mat_print_some ( nnzero, nrhs, rhsval, 1, 1, 5, 5, &
        '  Part of RHS array' )

    end if
!
!  Read the starting guesses.
!
    if ( rhstyp(2:2) == 'G' ) then

      allocate ( guess(1:nrow,1:nrhs) )

      write ( *, '(a)' ) '  Reading the starting guesses.'

      call hb_guess_read ( input_unit, nrow, nrhs, rhscrd, rhsfmt, rhstyp, &
        guess )

      call r8mat_print_some ( nrow, nrhs, guess, 1, 1, 5, 5, &
        '  Part of GUESS array' )

    end if
!
!  Read the exact solutions.
!
    if ( rhstyp(3:3) == 'X' ) then

      allocate ( exact(1:nrow,1:nrhs) )

      write ( *, '(a)' ) '  Reading the exact solutions.'

      call hb_exact_read ( input_unit, nrow, nrhs, rhscrd, rhsfmt, rhstyp, &
        exact )

      call r8mat_print_some ( nrow, nrhs, exact, 1, 1, 5, 5, &
        '  Part of EXACT array' )

    end if

  end if

  close ( unit = input_unit )

  if ( allocated ( colptr ) ) then
    deallocate ( colptr )
  end if
  if ( allocated ( exact ) ) then
    deallocate ( exact )
  end if
  if ( allocated ( guess ) ) then
    deallocate ( guess )
  end if
  if ( allocated ( rhsind ) ) then
    deallocate ( rhsind )
  end if
  if ( allocated ( rhsptr ) ) then
    deallocate ( rhsptr )
  end if
  if ( allocated ( rhsval ) ) then
    deallocate ( rhsval )
  end if
  if ( allocated ( rhsvec ) ) then
    deallocate ( rhsvec )
  end if
  if ( allocated ( rowind ) ) then
    deallocate ( rowind )
  end if
  if ( allocated ( values ) ) then
    deallocate ( values )
  end if

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests HB_RHS_WRITE, HB_GUESS_WRITE, HB_EXACT_WRITE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: neltvl = 0
  integer, parameter :: nnzero = 126
  integer, parameter :: nrhs = 1
  integer, parameter :: nrhsix = 0
  integer, parameter :: nrow = 32

  real ( kind = 8 ), dimension ( nrow, nrhs ) :: exact = reshape ( (/ &
    1.0,   2.0,   3.0,   4.0,   5.0,   6.0,   7.0,   8.0,   9.0,  10.0, &
   11.0,  12.0,  13.0,  14.0,  15.0,  16.0,  17.0,  18.0,  19.0,  20.0, &
   21.0,  22.0,  23.0,  24.0,  25.0,  26.0,  27.0,  28.0,  29.0,  30.0, &
   31.0,  32.0 /), (/ nrow, nrhs /) )
  real ( kind = 8 ), dimension ( nrow, nrhs ) :: guess = reshape ( (/ &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0 /), (/ nrow, nrhs /) )
  character ( len = 16 ) :: indfmt = '(16I5)'
  integer ios
  character ( len = 3 ) :: mxtype = 'RUA'
  character ( len = 80 ) :: output_file = 'rua_32_rhs.txt'
  integer output_unit
  character ( len = 16 ) :: ptrfmt = '(16I5)'
  integer :: rhscrd = 12
  character ( len = 20 ) :: rhsfmt = '(10F7.1)'
  integer, dimension ( 0 ) :: rhsind
  integer, dimension ( 0 ) :: rhsptr
  real ( kind = 8 ), dimension ( nrow, nrhs ) :: rhsval = reshape ( (/ &
    101.0, 102.0, 103.0, 104.0, 107.0, 126.0, 201.0, 202.0, 209.0, 221.0, &
    228.0, 302.0, 303.0, 306.0, 308.0, 309.0, 329.0, 403.0, 404.0, 405.0, &
    412.0, 503.0, 505.0, 523.0, 527.0, 601.0, 606.0, 616.0, 703.0, 707.0, &
    714.0, 721.0 /), (/ nrow, nrhs /) )
  real ( kind = 8 ), dimension (0) :: rhsvec
  character ( len = 3 ) :: rhstyp = 'FGX'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) &
    '  HB_RHS_WRITE writes the right hand sides to an HB file.'
  write ( *, '(a)' ) &
    '  HB_GUESS_WRITE writes starting guesses to an HB file.'
  write ( *, '(a)' ) &
    '  HB_EXACT_WRITE writes exact solutions to an HB file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Writing the file "' // trim ( output_file ) // '".'

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST08 - Fatal error!'
    write ( *, '(a)' ) '  Error opening the file.'
    return
  end if
!
!  Write the right hand sides.
!
  call hb_rhs_write ( output_unit, nrow, nnzero, nrhs, nrhsix, &
    rhscrd, ptrfmt, indfmt, rhsfmt, mxtype, rhstyp, rhsval, &
    rhsind, rhsptr, rhsvec )
!
!  Write the right hand sides.
!
  call hb_guess_write ( output_unit, nrow, nrhs, rhscrd, rhsfmt, rhstyp, &
    guess )
!
!  Write the right hand sides.
!
  call hb_exact_write ( output_unit, nrow, nrhs, rhscrd, rhsfmt, rhstyp, &
    exact )

  close ( unit = output_unit )

  return
end
subroutine test09 ( input_file )

!*****************************************************************************80
!
!! TEST09 tests HB_FILE_READ;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2008
!
!  Author:
!
!    John Burkardt
!
  use hb_file_module

  implicit none

  character ( len = * ) :: input_file
  integer input_unit
  integer ios

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  HB_FILE_READ reads all the data in an HB file.'
  write ( *, '(a)' ) '  HB_FILE_MODULE is the module that stores the data.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reading the file "' // trim ( input_file ) // '".'

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST09 - Fatal error!'
    write ( *, '(a)' ) '  Error opening the file.'
    return
  end if

  call hb_file_read ( input_unit )
!
!  Print out the header information.
!
  call hb_header_print ( title, key, totcrd, ptrcrd, indcrd, valcrd, &
    rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, valfmt, &
    rhsfmt, rhstyp, nrhs, nrhsix )
!
!  Print the structure information.
!
  call hb_structure_print ( ncol, mxtype, nnzero, neltvl, colptr, rowind )
!
!  Print the values.
!
  call hb_values_print ( ncol, colptr, mxtype, nnzero, neltvl, values )

  if ( 0 < rhscrd ) then
!
!  Print a bit of the right hand sides.
!
    if ( rhstyp(1:1) == 'F' ) then

      call r8mat_print_some ( nrow, nrhs, rhsval, 1, 1, 5, 5, &
        '  Part of RHS array' )

    else if ( rhstyp(1:1) == 'M' .and. mxtype(3:3) == 'A' ) then

      call i4vec_print_some ( nrhs+1, rhsptr, 1, 10, '  Part of RHSPTR' )
      call i4vec_print_some ( nrhsix, rhsind, 1, 10, '  Part of RHSIND' )
      call r8vec_print_some ( nrhsix, rhsvec, 1, 10, '  Part of RHSVEC' )

    else if ( rhstyp(1:1) == 'M' .and. mxtype(3:3) == 'E' ) then

      call r8mat_print_some ( nnzero, nrhs, rhsval, 1, 1, 5, 5, &
        '  Part of RHS array' )

    end if
!
!  Print a bit of the starting guesses.
!
    if ( rhstyp(2:2) == 'G' ) then

      call r8mat_print_some ( nrow, nrhs, guess, 1, 1, 5, 5, &
        '  Part of GUESS array' )

    end if
!
!  Print a bit of the exact solutions.
!
    if ( rhstyp(3:3) == 'X' ) then

      call r8mat_print_some ( nrow, nrhs, exact, 1, 1, 5, 5, &
        '  Part of EXACT array' )

    end if

  end if

  close ( unit = input_unit )

  if ( allocated ( colptr ) ) then
    deallocate ( colptr )
  end if
  if ( allocated ( exact ) ) then
    deallocate ( exact )
  end if
  if ( allocated ( guess ) ) then
    deallocate ( guess )
  end if
  if ( allocated ( rhsind ) ) then
    deallocate ( rhsind )
  end if
  if ( allocated ( rhsptr ) ) then
    deallocate ( rhsptr )
  end if
  if ( allocated ( rhsval ) ) then
    deallocate ( rhsval )
  end if
  if ( allocated ( rhsvec ) ) then
    deallocate ( rhsvec )
  end if
  if ( allocated ( rowind ) ) then
    deallocate ( rowind )
  end if
  if ( allocated ( values ) ) then
    deallocate ( values )
  end if

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests HB_FILE_WRITE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter       :: ncol = 32
  integer, parameter       :: neltvl = 0
  integer, parameter       :: nnzero = 126
  integer, parameter       :: nrhs = 1
  integer, parameter       :: nrhsix = 0
  integer, parameter       :: nrow = 32

  integer, dimension ( ncol + 1 ) :: colptr = (/ &
      1,   7,  12,  18,  22,  26,  29,  34,  39,  46, &
     53,  58,  61,  63,  65,  68,  71,  74,  79,  82, &
     85,  88,  90,  94,  97, 102, 106, 110, 112, 117, &
    121, 124, 127 /)
  real ( kind = 8 ), dimension ( nrow, nrhs ) :: exact = reshape ( (/ &
    1.0,   2.0,   3.0,   4.0,   5.0,   6.0,   7.0,   8.0,   9.0,  10.0, &
   11.0,  12.0,  13.0,  14.0,  15.0,  16.0,  17.0,  18.0,  19.0,  20.0, &
   21.0,  22.0,  23.0,  24.0,  25.0,  26.0,  27.0,  28.0,  29.0,  30.0, &
   31.0,  32.0 /), (/ nrow, nrhs /) )
  real ( kind = 8 ), dimension ( nrow, nrhs ) :: guess = reshape ( (/ &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0 /), (/ nrow, nrhs /) )
  integer                  :: indcrd = 8
  character ( len = 16 )   :: indfmt = '(16I5)'
  integer ios
  character ( len = 8 )    :: key = 'RUA_32'
  character ( len = 3 )    :: mxtype = 'RUA'
  character ( len = 80 )   :: output_file = 'rua_32_file.txt'
  integer output_unit
  integer                  :: ptrcrd = 3
  character ( len = 16 )   :: ptrfmt = '(16I5)'
  integer                  :: rhscrd = 12
  character ( len = 20 )   :: rhsfmt = '(10F7.1)'
  integer, dimension ( 0 ) :: rhsind
  integer, dimension ( 0 ) :: rhsptr
  real ( kind = 8 ), dimension ( nrow, nrhs ) :: rhsval = reshape ( (/ &
    101.0, 102.0, 103.0, 104.0, 107.0, 126.0, 201.0, 202.0, 209.0, 221.0, &
    228.0, 302.0, 303.0, 306.0, 308.0, 309.0, 329.0, 403.0, 404.0, 405.0, &
    412.0, 503.0, 505.0, 523.0, 527.0, 601.0, 606.0, 616.0, 703.0, 707.0, &
    714.0, 721.0 /), (/ nrow, nrhs /) )
  character ( len = 3 )  :: rhstyp = 'FGX'
  real ( kind = 8 ), allocatable, dimension ( : ) :: rhsvec
  integer, dimension ( nnzero ) :: rowind = (/ &
    1,    2,    3,    4,    7,   26,    1,    2,    9,   21, &
   28,    2,    3,    6,    8,    9,   29,    3,    4,    5, &
   12,    3,    5,   23,   27,    1,    6,   16,    3,    7, &
   14,   21,   31,    1,    8,   12,   17,   27,    7,    9, &
   10,   13,   19,   23,   27,    1,   10,   11,   21,   23, &
   25,   27,    2,   11,   15,   18,   29,    6,   12,   24, &
   11,   13,    3,   14,    2,   15,   20,    4,   16,   22, &
    4,   16,   17,    6,   10,   18,   20,   30,    1,   19, &
   26,    8,   16,   20,    3,   21,   32,   11,   22,    2, &
   17,   21,   23,   12,   24,   26,    6,   15,   18,   24, &
   25,   13,   18,   22,   26,    5,   24,   26,   27,    9, &
   28,    3,    5,   27,   29,   32,   12,   17,   23,   30, &
   13,   14,   31,   24,   28,   32 /)
  character ( len = 72 ) :: title = &
    '1Real unsymmetric assembled matrix based on IBM32'
  integer                :: totcrd = 36
  integer                :: valcrd = 13
  character ( len = 20 ) :: valfmt = '(10F7.1)'
  real ( kind = 8 ), dimension ( nnzero ) :: values = (/ &
  101.0,  102.0,  103.0,  104.0,  107.0, &
  126.0,  201.0,  202.0,  209.0,  221.0, &
  228.0,  302.0,  303.0,  306.0,  308.0, &
  309.0,  329.0,  403.0,  404.0,  405.0, &
  412.0,  503.0,  505.0,  523.0,  527.0, &
  601.0,  606.0,  616.0,  703.0,  707.0, &
  714.0,  721.0,  731.0,  801.0,  808.0, &
  812.0,  817.0,  827.0,  907.0,  909.0, &
  910.0,  913.0,  919.0,  923.0,  927.0, &
 1001.0, 1010.0, 1011.0, 1021.0, 1023.0, &
 1025.0, 1027.0, 1102.0, 1111.0, 1115.0, &
 1118.0, 1129.0, 1206.0, 1212.0, 1224.0, &
 1311.0, 1313.0, 1403.0, 1414.0, 1502.0, &
 1515.0, 1520.0, 1604.0, 1616.0, 1622.0, &
 1704.0, 1716.0, 1717.0, 1806.0, 1810.0, &
 1818.0, 1820.0, 1830.0, 1901.0, 1919.0, &
 1926.0, 2008.0, 2016.0, 2020.0, 2103.0, &
 2121.0, 2132.0, 2211.0, 2222.0, 2302.0, &
 2317.0, 2321.0, 2323.0, 2412.0, 2424.0, &
 2426.0, 2506.0, 2515.0, 2518.0, 2524.0, &
 2525.0, 2613.0, 2618.0, 2622.0, 2626.0, &
 2705.0, 2724.0, 2726.0, 2727.0, 2809.0, &
 2828.0, 2903.0, 2905.0, 2927.0, 2929.0, &
 2932.0, 3012.0, 3017.0, 3023.0, 3030.0, &
 3113.0, 3114.0, 3131.0, 3224.0, 3228.0, &
 3232.0 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  HB_FILE_WRITE writes an HB file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Writing the file "' // trim ( output_file ) // '".'

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST10 - Fatal error!'
    write ( *, '(a)' ) '  Error opening the file.'
    return
  end if

  call hb_file_write ( output_unit, title, key, totcrd, ptrcrd, indcrd, &
    valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, &
    valfmt, rhsfmt, rhstyp, nrhs, nrhsix, colptr, rowind, values, &
    rhsval, rhsptr, rhsind, rhsvec, guess, exact )

  close ( unit = output_unit )

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests HB_MATVEC_A_MEM;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 September 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter       :: ncol = 32
  integer, parameter       :: neltvl = 0
  integer, parameter       :: nnzero = 126
  integer, parameter       :: nrhs = 2
  integer, parameter       :: nrhsix = 0
  integer, parameter       :: nrow = 32

  integer, dimension ( ncol + 1 ) :: colptr = (/ &
      1,   7,  12,  18,  22,  26,  29,  34,  39,  46, &
     53,  58,  61,  63,  65,  68,  71,  74,  79,  82, &
     85,  88,  90,  94,  97, 102, 106, 110, 112, 117, &
    121, 124, 127 /)
  real ( kind = 8 ), dimension ( ncol, nrhs ) :: exact = reshape ( (/ &
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   1.0, &
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, &
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, &
    0.0,   0.0, &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0 /), (/ nrow, nrhs /) )
  real ( kind = 8 ), dimension ( nrow, nrhs ) :: guess = reshape ( (/ &
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, &
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, &
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, &
    0.0,   0.0, &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0 /), (/ nrow, nrhs /) )
  integer                  :: indcrd = 8
  character ( len = 16 )   :: indfmt = '(16I5)'
  integer ios
  character ( len = 8 )    :: key = 'RUA_32'
  character ( len = 3 )    :: mxtype = 'RUA'
  character ( len = 80 )   :: output_file = 'rua_32_ax.txt'
  integer output_unit
  integer                  :: ptrcrd = 3
  character ( len = 16 )   :: ptrfmt = '(16I5)'
  integer                  :: rhscrd = 12
  character ( len = 20 )   :: rhsfmt = '(10F7.1)'
  integer, dimension ( 0 ) :: rhsind
  integer, dimension ( 0 ) :: rhsptr
  real ( kind = 8 ), dimension ( nrow, nrhs ) :: rhsval
  character ( len = 3 )  :: rhstyp = 'FGX'
  real ( kind = 8 ), allocatable, dimension ( : ) :: rhsvec
  integer, dimension ( nnzero ) :: rowind = (/ &
    1,    2,    3,    4,    7,   26,    1,    2,    9,   21, &
   28,    2,    3,    6,    8,    9,   29,    3,    4,    5, &
   12,    3,    5,   23,   27,    1,    6,   16,    3,    7, &
   14,   21,   31,    1,    8,   12,   17,   27,    7,    9, &
   10,   13,   19,   23,   27,    1,   10,   11,   21,   23, &
   25,   27,    2,   11,   15,   18,   29,    6,   12,   24, &
   11,   13,    3,   14,    2,   15,   20,    4,   16,   22, &
    4,   16,   17,    6,   10,   18,   20,   30,    1,   19, &
   26,    8,   16,   20,    3,   21,   32,   11,   22,    2, &
   17,   21,   23,   12,   24,   26,    6,   15,   18,   24, &
   25,   13,   18,   22,   26,    5,   24,   26,   27,    9, &
   28,    3,    5,   27,   29,   32,   12,   17,   23,   30, &
   13,   14,   31,   24,   28,   32 /)
  character ( len = 72 ) :: title = &
    '1Real unsymmetric assembled matrix based on IBM32'
  integer                :: totcrd = 36
  integer                :: valcrd = 13
  character ( len = 20 ) :: valfmt = '(10F7.1)'
  real ( kind = 8 ), dimension ( nnzero ) :: values = (/ &
  101.0,  102.0,  103.0,  104.0,  107.0, &
  126.0,  201.0,  202.0,  209.0,  221.0, &
  228.0,  302.0,  303.0,  306.0,  308.0, &
  309.0,  329.0,  403.0,  404.0,  405.0, &
  412.0,  503.0,  505.0,  523.0,  527.0, &
  601.0,  606.0,  616.0,  703.0,  707.0, &
  714.0,  721.0,  731.0,  801.0,  808.0, &
  812.0,  817.0,  827.0,  907.0,  909.0, &
  910.0,  913.0,  919.0,  923.0,  927.0, &
 1001.0, 1010.0, 1011.0, 1021.0, 1023.0, &
 1025.0, 1027.0, 1102.0, 1111.0, 1115.0, &
 1118.0, 1129.0, 1206.0, 1212.0, 1224.0, &
 1311.0, 1313.0, 1403.0, 1414.0, 1502.0, &
 1515.0, 1520.0, 1604.0, 1616.0, 1622.0, &
 1704.0, 1716.0, 1717.0, 1806.0, 1810.0, &
 1818.0, 1820.0, 1830.0, 1901.0, 1919.0, &
 1926.0, 2008.0, 2016.0, 2020.0, 2103.0, &
 2121.0, 2132.0, 2211.0, 2222.0, 2302.0, &
 2317.0, 2321.0, 2323.0, 2412.0, 2424.0, &
 2426.0, 2506.0, 2515.0, 2518.0, 2524.0, &
 2525.0, 2613.0, 2618.0, 2622.0, 2626.0, &
 2705.0, 2724.0, 2726.0, 2727.0, 2809.0, &
 2828.0, 2903.0, 2905.0, 2927.0, 2929.0, &
 2932.0, 3012.0, 3017.0, 3023.0, 3030.0, &
 3113.0, 3114.0, 3131.0, 3224.0, 3228.0, &
 3232.0 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  HB_MATVEC_A_MEM multiplies a matrix times a vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This particular version assumes:'
  write ( *, '(a)' ) '  * the matrix is in "A" format (assembled),'
  write ( *, '(a)' ) '  * the matrix and vectors can fit in memory,'
  write ( *, '(a)' ) '  * the matrix and multiplicand have been read into'
  write ( *, '(a)' ) '    memory before the routine is called.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this example, the first vector X is zero except'
  write ( *, '(a)' ) '  for a 1 in row 10.  This means A*X should return'
  write ( *, '(a)' ) '  column 10 of A.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The second vector X is all 1''s.  A*X should be'
  write ( *, '(a)' ) '  the sum of the entries of each row.'

  call hb_matvec_a_mem ( nrow, ncol, nnzero, nrhs, colptr, rowind, values, &
    exact, rhsval )

  call r8mat_print ( nrow, nrhs, rhsval,  '  The product vectors A*X' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Writing the file "' // trim ( output_file ) // '".'

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST11 - Fatal error!'
    write ( *, '(a)' ) '  Error opening the file.'
    return
  end if

  call hb_file_write ( output_unit, title, key, totcrd, ptrcrd, indcrd, &
    valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, &
    valfmt, rhsfmt, rhstyp, nrhs, nrhsix, colptr, rowind, values, &
    rhsval, rhsptr, rhsind, rhsvec, guess, exact )

  close ( unit = output_unit )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests HB_VECMAT_A_MEM;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2004
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter       :: ncol = 32
  integer, parameter       :: neltvl = 0
  integer, parameter       :: nnzero = 126
  integer, parameter       :: nrhs = 2
  integer, parameter       :: nrhsix = 0
  integer, parameter       :: nrow = 32

  integer, dimension ( ncol + 1 ) :: colptr = (/ &
      1,   7,  12,  18,  22,  26,  29,  34,  39,  46, &
     53,  58,  61,  63,  65,  68,  71,  74,  79,  82, &
     85,  88,  90,  94,  97, 102, 106, 110, 112, 117, &
    121, 124, 127 /)
  real ( kind = 8 ), dimension ( nrow, nrhs ) :: exact = reshape ( (/ &
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   1.0, &
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, &
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, &
    0.0,   0.0, &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0 /), (/ nrow, nrhs /) )
  real ( kind = 8 ), dimension ( nrow, nrhs ) :: guess = reshape ( (/ &
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, &
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, &
    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, &
    0.0,   0.0, &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0, &
    1.0,   1.0 /), (/ nrow, nrhs /) )
  integer                  :: indcrd = 8
  character ( len = 16 )   :: indfmt = '(16I5)'
  integer ios
  character ( len = 8 )    :: key = 'RUA_32'
  character ( len = 3 )    :: mxtype = 'RUA'
  character ( len = 80 )   :: output_file = 'rua_32_xa.txt'
  integer output_unit
  integer                  :: ptrcrd = 3
  character ( len = 16 )   :: ptrfmt = '(16I5)'
  integer                  :: rhscrd = 12
  character ( len = 20 )   :: rhsfmt = '(10F7.1)'
  integer, dimension ( 0 ) :: rhsind
  integer, dimension ( 0 ) :: rhsptr
  real ( kind = 8 ), dimension ( ncol, nrhs ) :: rhsval
  character ( len = 3 )  :: rhstyp = 'FGX'
  real ( kind = 8 ), allocatable, dimension ( : ) :: rhsvec
  integer, dimension ( nnzero ) :: rowind = (/ &
    1,    2,    3,    4,    7,   26,    1,    2,    9,   21, &
   28,    2,    3,    6,    8,    9,   29,    3,    4,    5, &
   12,    3,    5,   23,   27,    1,    6,   16,    3,    7, &
   14,   21,   31,    1,    8,   12,   17,   27,    7,    9, &
   10,   13,   19,   23,   27,    1,   10,   11,   21,   23, &
   25,   27,    2,   11,   15,   18,   29,    6,   12,   24, &
   11,   13,    3,   14,    2,   15,   20,    4,   16,   22, &
    4,   16,   17,    6,   10,   18,   20,   30,    1,   19, &
   26,    8,   16,   20,    3,   21,   32,   11,   22,    2, &
   17,   21,   23,   12,   24,   26,    6,   15,   18,   24, &
   25,   13,   18,   22,   26,    5,   24,   26,   27,    9, &
   28,    3,    5,   27,   29,   32,   12,   17,   23,   30, &
   13,   14,   31,   24,   28,   32 /)
  character ( len = 72 ) :: title = &
    '1Real unsymmetric assembled matrix based on IBM32'
  integer                :: totcrd = 36
  integer                :: valcrd = 13
  character ( len = 20 ) :: valfmt = '(10F7.1)'
  real ( kind = 8 ), dimension ( nnzero ) :: values = (/ &
  101.0,  102.0,  103.0,  104.0,  107.0, &
  126.0,  201.0,  202.0,  209.0,  221.0, &
  228.0,  302.0,  303.0,  306.0,  308.0, &
  309.0,  329.0,  403.0,  404.0,  405.0, &
  412.0,  503.0,  505.0,  523.0,  527.0, &
  601.0,  606.0,  616.0,  703.0,  707.0, &
  714.0,  721.0,  731.0,  801.0,  808.0, &
  812.0,  817.0,  827.0,  907.0,  909.0, &
  910.0,  913.0,  919.0,  923.0,  927.0, &
 1001.0, 1010.0, 1011.0, 1021.0, 1023.0, &
 1025.0, 1027.0, 1102.0, 1111.0, 1115.0, &
 1118.0, 1129.0, 1206.0, 1212.0, 1224.0, &
 1311.0, 1313.0, 1403.0, 1414.0, 1502.0, &
 1515.0, 1520.0, 1604.0, 1616.0, 1622.0, &
 1704.0, 1716.0, 1717.0, 1806.0, 1810.0, &
 1818.0, 1820.0, 1830.0, 1901.0, 1919.0, &
 1926.0, 2008.0, 2016.0, 2020.0, 2103.0, &
 2121.0, 2132.0, 2211.0, 2222.0, 2302.0, &
 2317.0, 2321.0, 2323.0, 2412.0, 2424.0, &
 2426.0, 2506.0, 2515.0, 2518.0, 2524.0, &
 2525.0, 2613.0, 2618.0, 2622.0, 2626.0, &
 2705.0, 2724.0, 2726.0, 2727.0, 2809.0, &
 2828.0, 2903.0, 2905.0, 2927.0, 2929.0, &
 2932.0, 3012.0, 3017.0, 3023.0, 3030.0, &
 3113.0, 3114.0, 3131.0, 3224.0, 3228.0, &
 3232.0 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  HB_VECMAT_A_MEM multiplies a vector times a matrix'
  write ( *, '(a)' ) '  b'' = x'' * A,'
  write ( *, '(a)' ) '  or, equivalently, computes b = A''*x.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This particular version assumes:'
  write ( *, '(a)' ) '  * the matrix is in "A" format (assembled),'
  write ( *, '(a)' ) '  * the matrix and vectors can fit in memory,'
  write ( *, '(a)' ) '  * the matrix and multiplicand have been read into'
  write ( *, '(a)' ) '    memory before the routine is called.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this example, the first vector X is zero except'
  write ( *, '(a)' ) '  for a 1 in row 10.  This means A''*X should return'
  write ( *, '(a)' ) '  row 10 of A.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The second vector X is all 1''s.  A''*X should be'
  write ( *, '(a)' ) '  the sum of the entries of each column.'

  call hb_vecmat_a_mem ( nrow, ncol, nnzero, nrhs, colptr, rowind, values, &
    exact, rhsval )

  call r8mat_print ( ncol, nrhs, rhsval,  '  The product vectors A''*X' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Writing the file "' // trim ( output_file ) // '".'

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST12 - Fatal error!'
    write ( *, '(a)' ) '  Error opening the file.'
    return
  end if

  call hb_file_write ( output_unit, title, key, totcrd, ptrcrd, indcrd, &
    valcrd, rhscrd, mxtype, nrow, ncol, nnzero, neltvl, ptrfmt, indfmt, &
    valfmt, rhsfmt, rhstyp, nrhs, nrhsix, colptr, rowind, values, &
    rhsval, rhsptr, rhsind, rhsvec, guess, exact )

  close ( unit = output_unit )

  return
end
