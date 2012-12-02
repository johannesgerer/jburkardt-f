program main

!*****************************************************************************80
!
!! MAIN is the main program for MM_IO_PRB.
!
!  Discussion:
!
!    MM_IO_PRB calls the MM_IO test problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MM_IO_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the MM_IO library.'

  call test01 ( 'matrix_05_05_crg.txt' )
  call test02 ( 'matrix_05_05_crg.txt' )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MM_IO_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( input_file )

!*****************************************************************************80
!
!! TEST01 tests MM_HEADER_READ
!
!  Discussion:
!
!    The size information can be handled either by
!    MM_SIZE_READ_FILE (by backspacing once the comments have been read)
!    or by MM_SIZE_READ_STRING (by passing the most recently read line,
!    which is NOT a comment).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) comment
  integer ( kind = 4 ) nnz
  character ( len = 7 ) field
  character ( len = 14 ) id
  character ( len = * ) :: input_file
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  character ( len = 10 ) rep
  character ( len = 19 ) symm
  character ( len = 6 ) type

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  MM_HEADER_READ reads the header line of'
  write ( *, '(a)' ) '    a Matrix Market file.'
  write ( *, '(a)' ) '  MM_SIZE_READ_FILE or MM_SIZE_READ_STRING'
  write ( *, '(a)' ) '    reads the size line of a Matrix Market file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reading "' // trim ( input_file ) // '".'

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file.'
    stop
  end if

  call mm_header_read ( input_unit, id, type, rep, field, symm )

  call mm_header_check ( id, type, rep, field, symm )

  call mm_header_print ( input_file, id, type, rep, field, symm )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Comment lines:'
  write ( *, '(a)' ) ' '

  do

    call mm_comment_read ( input_unit, comment )

    if ( comment(1:1) /= '%' ) then

      if ( .true. ) then

      else
        backspace ( input_unit )
      end if

      exit

    end if

    call mm_comment_print ( comment )

  end do

  if ( .true. ) then
    call mm_size_read_string ( comment, rep, symm, nrow, ncol, nnz )
  else
    call mm_size_read_file ( input_unit, rep, symm, nrow, ncol, nnz )
  end if

  call mm_size_print ( input_file, rep, symm, nrow, ncol, nnz )

  close ( unit = input_unit )

  return
end
subroutine test02 ( input_file )

!*****************************************************************************80
!
!! TEST02 tests MM_FILE_READ
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnzmax = 100

  complex ( kind = 4 ) cval(nnzmax)
  real ( kind = 8 ) dval(nnzmax)
  integer ( kind = 4 ) nnz
  character ( len = 7 ) field
  character ( len = 14 ) id
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) indx(nnzmax)
  character ( len = * ) :: input_file
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ival(nnzmax)
  integer ( kind = 4 ) jndx(nnzmax)
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  character ( len = 10 ) rep
  real ( kind = 4 ) rval(nnzmax)
  character ( len = 19 ) symm
  character ( len = 6 ) type

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  MM_FILE_READ reads a Matrix Market file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reading "' // trim ( input_file ) // '".'

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file.'
    stop
  end if

  call mm_file_read ( input_unit, id, type, rep, field, symm, nrow, ncol, &
    nnz, nnzmax, indx, jndx, ival, rval, dval, cval )

  close ( unit = input_unit )

  call mm_header_print ( input_file, id, type, rep, field, symm )

  call mm_size_print ( input_file, rep, symm, nrow, ncol, nnz )

  ilo = 1
  ihi = 5

  call mm_values_print_some ( rep, field, nnz, indx, jndx, ival, rval, &
    dval, cval, ilo, ihi )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests MM_FILE_WRITE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nnz = 10

  complex ( kind = 4 ), dimension ( nnz ) :: cval
  real ( kind = 8 ), dimension ( nnz ) :: dval
  character ( len = 7 ) :: field = 'real'
  character ( len = 14 ) :: id = '%%MatrixMarket'
  integer ( kind = 4 ), dimension ( nnz ) :: indx = (/ &
    1, 1, 2, 2, 3, 3, 4, 5, 5, 5 /)
  integer ( kind = 4 ) ios
  integer ( kind = 4 ), dimension ( nnz ) :: ival
  integer ( kind = 4 ), dimension ( nnz ) :: jndx = (/ &
    1, 5, 3, 4, 2, 5, 1, 2, 4, 5 /)
  integer ( kind = 4 ) :: ncol = 5
  integer ( kind = 4 ) :: nrow = 5
  character ( len = 80 ) :: output_file = 'matrix_05_05_crg.txt'
  integer ( kind = 4 ) output_unit
  character ( len = 10 ) :: rep = 'coordinate'
  real ( kind = 4 ), dimension ( nnz ) :: rval = (/ &
    11.0E+00, 15.0E+00, 23.0E+00, 24.0E+00, 32.0E+00, &
    35.0E+00, 41.0E+00, 52.0E+00, 54.0E+00, 55.0E+00 /)
  character ( len = 19 ) :: symm = 'general'
  character ( len = 6 ) :: type = 'matrix'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  MM_FILE_WRITE writes the header and data of'
  write ( *, '(a)' ) '  a Matrix Market file.'

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST03 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if

  call mm_file_write ( output_unit, id, type, rep, field, symm, nrow, &
    ncol, nnz, indx, jndx, ival, rval, dval, cval )

  close ( unit = output_unit )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests MM_FILE_WRITE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 4 ), dimension ( 0 ) :: cval
  real (  kind = 8 ), allocatable, dimension ( :, : ) :: dval
  character ( len = 7 ) :: field = 'double'
  character ( len = 14 ) :: id = '%%MatrixMarket'
  integer ( kind = 4 ), dimension ( 0 ) :: indx
  integer ( kind = 4 ) ios
  integer ( kind = 4 ), dimension ( 0 ) :: ival
  integer ( kind = 4 ), dimension ( 0 ) :: jndx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nnz
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ), parameter :: nx = 3
  integer ( kind = 4 ), parameter :: ny = 2
  character ( len = 80 ) :: output_file = 'wathen_29_29_adg.txt'
  integer ( kind = 4 ) output_unit
  character ( len = 10 ) :: rep = 'array'
  real ( kind = 4 ), dimension ( 0 ) :: rval
  character ( len = 19 ) :: symm = 'general'
  character ( len = 6 ) :: type = 'matrix'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04:'
  write ( *, '(a)' ) '  MM_FILE_WRITE writes the header and data of'
  write ( *, '(a)' ) '  a Matrix Market file.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we generate the Wathen matrix'
  write ( *, '(a)' ) '  of order 29 x 29, and store it (inefficiently)'
  write ( *, '(a)' ) '  as an array.'
!
!  Generate the WATHEN matrix.
!
  call wathen_size ( nx, ny, n )

  allocate ( dval(1:n,1:n) )

  call wathen ( nx, ny, n, dval )
!
!  Write the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST04 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if

  nrow = n
  ncol = n
  nnz = nrow * ncol

  call mm_file_write ( output_unit, id, type, rep, field, symm, nrow, &
    ncol, nnz, indx, jndx, ival, rval, dval, cval )

  deallocate ( dval )

  close ( unit = output_unit )

  return
end
subroutine wathen ( nx, ny, n, a )

!*****************************************************************************80
!
!! WATHEN returns a finite element matrix.
!
!  Discussion:
!
!    A is the consistent mass matrix for a regular NX by NY grid
!    of 8 node serendipity elements.  Here is an illustration
!    for NX = 3, NX = 2:
!
!     23-24-25-26-27-28-29
!      |     |     |     |
!     19    20    21    22
!      |     |     |     |
!     12-13-14-15-16-17-18
!      |     |     |     |
!      8     9    10    11
!      |     |     |     |
!      1--2--3--4--5--6--7
!
!    For this example, the total number of nodes is, as expected,
!
!      N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29
!
!  Properties:
!
!    A is symmetric positive definite for any positive values of the
!    density RHO(NX,NY), which is here given the value 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    A J Wathen,
!    Realistic eigenvalue bounds for the Galerkin mass matrix,
!    IMA Journal of Numerical Analysis,
!    Volume 7, 1987, pages 449-457.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, values which determine the size of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix, as determined
!    by NX and NY.
!
!    Output, real ( kind = 8 ) A(N,N), the matrix.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ), save, dimension ( 8, 8 ) :: em =  reshape ( (/ &
     6.0, -6.0,  2.0, -8.0,  3.0, -8.0,  2.0, -6.0, &
    -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, &
     2.0, -6.0,  6.0, -6.0,  2.0, -8.0,  3.0, -8.0, &
    -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, &
     3.0, -8.0,  2.0, -6.0,  6.0, -6.0,  2.0, -8.0, &
    -8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, &
     2.0, -8.0,  3.0, -8.0,  2.0, -6.0,  6.0, -6.0, &
    -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0 /), &
    (/ 8, 8 /) )
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kcol
  integer ( kind = 4 ) krow
  integer ( kind = 4 ) node(8)
  real ( kind = 8 ) rho

  a(1:n,1:n) = 0.0D+00

  do j = 1, ny

    do i = 1, nx
!
!  For the element (I,J), determine the indices of the 8 nodes.
!
      node(1) = 3 * j * nx + 2 * j + 2 * i + 1
      node(2) = node(1) - 1
      node(3) = node(1) - 2

      node(4) = ( 3 * j - 1 ) * nx + 2 * j + i - 1
      node(8) = node(4) + 1

      node(5) = ( 3 * j - 3 ) * nx + 2 * j + 2 * i - 3
      node(6) = node(5) + 1
      node(7) = node(5) + 2
!
!  The density RHO can also be set to a random positive value.
!
      rho = 1.0D+00

      do krow = 1, 8
        do kcol = 1, 8

          if ( node(krow) < 1 .or. n < node(krow) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'WATHEN - Fatal error!'
            write ( *, '(a)' ) '  Index NODE(KROW) out of bounds.'
            write ( *, '(a)' ) '  I = ', i
            write ( *, '(a)' ) '  J = ', j
            write ( *, '(a)' ) '  KROW = ', krow
            write ( *, '(a)' ) '  NODE(KROW) = ', node(krow)
            return
          else if ( node(kcol) < 1 .or. n < node(kcol) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'WATHEN - Fatal error!'
            write ( *, '(a)' ) '  Index NODE(KCOL) out of bounds.'
            write ( *, '(a)' ) '  I = ', i
            write ( *, '(a)' ) '  J = ', j
            write ( *, '(a)' ) '  KCOL = ', kcol
            write ( *, '(a)' ) '  NODE(KCOL) = ', node(kcol)
            return
          end if

          a(node(krow),node(kcol)) = a(node(krow),node(kcol)) &
            + 20.0D+00 * rho * em(krow,kcol) / 9.0D+00

        end do
      end do

    end do
  end do

  return
end
subroutine wathen_size ( nx, ny, n )

!*****************************************************************************80
!
!! WATHEN_SIZE returns the size of a finite element matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    A J Wathen,
!    Realistic eigenvalue bounds for the Galerkin mass matrix,
!    IMA Journal of Numerical Analysis,
!    Volume 7, 1987, pages 449-457.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NX, NY, values which determine the size of A.
!
!    Output, integer ( kind = 4 ) N, the order of the matrix, as determined
!    by NX and NY.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny

  n = 3 * nx * ny + 2 * nx + 2 * ny + 1

  return
end
