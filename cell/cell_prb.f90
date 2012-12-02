program main

!*****************************************************************************80
!
!! MAIN tests CELL.
!
!  Discussion:
!
!    An R8CVV is a "cell vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CELL_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the CELL library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CELL_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 stores some of Pascal's triangle in an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "cell array vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 November 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5

  real ( kind = 8 ), allocatable :: a(:)
  real ( kind = 8 ), allocatable :: ai(:)
  real ( kind = 8 ) aij
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ), allocatable :: in(:)
  integer ( kind = 4 ) j
  integer ( kind = 4 ), allocatable :: jn(:)
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) nn
  integer ( kind = 4 ), dimension ( 5 ) :: nr = (/ &
    4, 5, 6, 7, 8 /)
  integer ( kind = 4 ) nv
  integer ( kind = 4 ), allocatable :: roff(:)
  integer ( kind = 4 ) row
  real ( kind = 8 ), allocatable :: vn(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use a cell array (vector of vectors) to store rows 3:7'
  write ( *, '(a)' ) '  of Pascal''s triangle.'

  call i4vec_print ( m, nr, '  The row sizes:' )
!
!  From the NR information:
!  * determine the total size, MN
!
  call r8cvv_size ( m, nr, mn )
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  The storage for the cell array is ', mn
!
!  Allocate the cell array.
!
  allocate ( a(1:mn) )
!
!  Zero out the cell array.
!
  a(1:mn) = 0.0D+00
!
!  Allocate a vector big enough to hold any single row.
!
  nr_max = maxval ( nr(1:m) )
  allocate ( ai(1:nr_max) )
!
!  From the NR information:
!  * determine the offsets.
!
  allocate ( roff(1:m+1) )
  call r8cvv_offset ( m, nr, roff )
  call i4vec_print ( m + 1, roff, '  The row offsets:' )
!
!  Rows 1 through 5 of A will contain rows 3 through 7 of Pascal's triangle.
!  Set these values one row at a time.
!
  ai(1) = 1.0D+00

  do row = 1, 7

    col = row + 1
    ai(col) = ai(col-1)
    do col = row, 2, -1
      ai(col) = ai(col) + ai(col-1)
    end do

    if ( 3 <= row ) then
      i = row - 2
      call r8cvv_rset ( mn, a, m, roff, i, ai )
    end if

  end do
!
!  Print the cell array.
!
  call r8cvv_print ( mn, a, m, roff, '  Rows 3:7 of Pascal''s Triangle:' )
!
!  Retrieve the entry from cell array row 3, column 4:
!
  i = 3
  j = 4
  call r8cvv_iget ( mn, a, m, roff, i, j, aij )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i1,a,i1,a,g14.6)' ) '  A(', i, ',', j, ') = ', aij
!
!  Retrieve row 4:
!
  i = 4
  call r8cvv_rget ( mn, a, m, roff, i, ai )
  nv = roff(i+1) - roff(i)
  call r8vec_transpose_print ( nv, ai, '  A(4,*):' )
!
!  Retrieve a list of entries.
!
  nn = 4
  allocate ( in(1:nn) )
  allocate ( jn(1:nn) )
  allocate ( vn(1:nn) )
  in = (/ 1, 2, 5, 5 /)
  jn = (/ 2, 3, 4, 8 /)
  call r8cvv_nget ( mn, a, m, roff, nn, in, jn, vn )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Retrieve an arbitrary list of items:'
  write ( *, '(a)' ) ' '
  do i = 1, nn
    write ( *, '(a,i1,a,i1,a,g14.6)' ) '  A(', in(i), ',', jn(i), ') = ', vn(i)
  end do
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( ai )
  deallocate ( in )
  deallocate ( jn )
  deallocate ( roff )
  deallocate ( vn )

  return
end

