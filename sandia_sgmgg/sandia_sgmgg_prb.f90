program main

!*****************************************************************************80
!
!! MAIN is the main program for SANDIA_SGMGG_PRB.
!
!  Discussion:
!
!    SANDIA_SGMGG_PRB tests generalized sparse grid routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2011
!
!  Author:
!
!    John Burkardt
!
  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SANDIA_SGMGG_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SANDIA_SGMGG library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SANDIA_SGMGG_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates the naive coefficient calculations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 August 2010
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), allocatable :: sparse_index(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Demonstrate naive coefficient calculations.'
!
!  Isotropic grid in 2D.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  1) Isotropic grid in 2D'

  dim_num = 2
  point_num = 7
  allocate ( sparse_index(dim_num,point_num) )

  sparse_index = reshape ( (/ &
    0, 2, &
    0, 3, &
    1, 1, &
    1, 2, &
    2, 0, &
    2, 1, &
    3, 0 /), (/ dim_num, point_num /) )
    
  call sandia_sgmgg_coef_naive_test ( dim_num, point_num, sparse_index  )

  deallocate ( sparse_index )
!
!  Isotropic grid in 3D.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  2) Isotropic grid in 3D'

  dim_num = 3
  point_num = 19
  allocate ( sparse_index(dim_num,point_num) )

  sparse_index = reshape ( (/ &
    0, 1, 0, &
    0, 2, 0, &
    0, 3, 0, &
    1, 0, 0, &
    1, 1, 0, &
    1, 2, 0, &
    2, 0, 0, &
    2, 1, 0, &
    3, 0, 0, &
    0, 0, 1, &
    0, 1, 1, &
    0, 2, 1, &
    1, 0, 1, &
    1, 1, 1, &
    2, 0, 1, &
    0, 0, 2, &
    0, 1, 2, &
    1, 0, 2, &
    0, 0, 3 /), (/ dim_num, point_num /) )
    
  call sandia_sgmgg_coef_naive_test ( dim_num, point_num, sparse_index  )

  deallocate ( sparse_index )
!
!  Anisotropic grid in 2D.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  3) Anisotropic grid in 2D'

  dim_num = 2
  point_num = 8
  allocate ( sparse_index(dim_num,point_num) )

  sparse_index = reshape ( (/ &
    0, 2, &
    1, 1, &
    1, 2, &
    2, 1, &
    3, 0, &
    3, 1, &
    4, 0, &
    5, 0 /), (/ dim_num, point_num /) )
    
  call sandia_sgmgg_coef_naive_test ( dim_num, point_num, sparse_index  )

  deallocate ( sparse_index )
!
!  Generalized grid in 2D.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  4) Generalized grid in 2D'

  dim_num = 2
  point_num = 8
  allocate ( sparse_index(dim_num,point_num) )

  sparse_index = reshape ( (/ &
    0, 0, &
    0, 1, &
    0, 2, &
    0, 3, &
    1, 0, &
    1, 1, &
    2, 0, &
    3, 0 /), (/ dim_num, point_num /) )
    
  call sandia_sgmgg_coef_naive_test ( dim_num, point_num, sparse_index  )

  deallocate ( sparse_index )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 demonstrates the naive neighbor calculations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2010
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), allocatable :: sparse_index(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  Demonstrate naive neighbor calculations.'
!
!  Isotropic grid in 2D.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  1) Isotropic grid in 2D'

  dim_num = 2
  point_num = 7
  allocate ( sparse_index(dim_num,point_num) )

  sparse_index = reshape ( (/ &
    0, 2, &
    0, 3, &
    1, 1, &
    1, 2, &
    2, 0, &
    2, 1, &
    3, 0 /), (/ dim_num, point_num /) )
    
  call sandia_sgmgg_neighbor_naive_test ( dim_num, point_num, sparse_index  )

  deallocate ( sparse_index )
!
!  Isotropic grid in 3D.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  2) Isotropic grid in 3D'

  dim_num = 3
  point_num = 19
  allocate ( sparse_index(dim_num,point_num) )

  sparse_index = reshape ( (/ &
    0, 1, 0, &
    0, 2, 0, &
    0, 3, 0, &
    1, 0, 0, &
    1, 1, 0, &
    1, 2, 0, &
    2, 0, 0, &
    2, 1, 0, &
    3, 0, 0, &
    0, 0, 1, &
    0, 1, 1, &
    0, 2, 1, &
    1, 0, 1, &
    1, 1, 1, &
    2, 0, 1, &
    0, 0, 2, &
    0, 1, 2, &
    1, 0, 2, &
    0, 0, 3 /), (/ dim_num, point_num /) )
    
  call sandia_sgmgg_neighbor_naive_test ( dim_num, point_num, sparse_index  )

  deallocate ( sparse_index )
!
!  Anisotropic grid in 2D.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  3) Anisotropic grid in 2D'

  dim_num = 2
  point_num = 8
  allocate ( sparse_index(dim_num,point_num) )

  sparse_index = reshape ( (/ &
    0, 2, &
    1, 1, &
    1, 2, &
    2, 1, &
    3, 0, &
    3, 1, &
    4, 0, &
    5, 0 /), (/ dim_num, point_num /) )
    
  call sandia_sgmgg_neighbor_naive_test ( dim_num, point_num, sparse_index  )

  deallocate ( sparse_index )
!
!  Generalized grid in 2D.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  4) Generalized grid in 2D'

  dim_num = 2
  point_num = 8
  allocate ( sparse_index(dim_num,point_num) )

  sparse_index = reshape ( (/ &
    0, 0, &
    0, 1, &
    0, 2, &
    0, 3, &
    1, 0, &
    1, 1, &
    2, 0, &
    3, 0 /), (/ dim_num, point_num /) )
    
  call sandia_sgmgg_neighbor_naive_test ( dim_num, point_num, sparse_index  )

  deallocate ( sparse_index )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 sets up the GG data structure.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ) d1
  integer   ( kind = 4 ) d2
  integer   ( kind = 4 ) d3
  integer   ( kind = 4 ), allocatable :: gg_a(:)
  integer   ( kind = 4 ), allocatable :: gg_b(:,:)
  integer   ( kind = 4 ), allocatable :: gg_f(:,:)
  real      ( kind = 8 ), allocatable :: gg_g(:)
  integer   ( kind = 4 ), allocatable :: gg_i(:,:)
  integer   ( kind = 4 ) gg_ma
  integer   ( kind = 4 ) gg_mi
  integer   ( kind = 4 ) gg_mo
  integer   ( kind = 4 ) gg_na
  integer   ( kind = 4 ) gg_nd
  integer   ( kind = 4 ) gg_ni
  integer   ( kind = 4 ) gg_no
  integer   ( kind = 4 ), allocatable :: gg_o(:)
  integer   ( kind = 4 ), allocatable :: gg_s(:)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i4_wrap
  integer   ( kind = 4 ) indx_max
  integer   ( kind = 4 ) indx_nab
  integer   ( kind = 4 ) indx_new
  integer   ( kind = 4 ) nb
  integer   ( kind = 4 ) nbf
  real      ( kind = 8 ) r8_uniform_01
  integer   ( kind = 4 ) seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  Set up examples of the GG (Gerstner-Griebel)'
  write ( *, '(a)' ) '  data structure for generalized sparse grids.'
!
!  Isotropic grid in 2D.
!
!  3 | 4
!  2 | 3  7
!  1 | 2  6  9 
!  0 | 1  5  8 10
!    +-----------
!      0  1  2  3
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  1) Isotropic grid in 2D'

  gg_na = 4
  gg_nd = 2
  gg_ni = 10
  gg_no = 6

  gg_ma = 20
  gg_mi = 20
  gg_mo = 20

  allocate ( gg_a(gg_ma) )
  allocate ( gg_b(gg_nd,gg_mi) )
  allocate ( gg_f(gg_nd,gg_mi) )
  allocate ( gg_g(gg_mi) )
  allocate ( gg_i(gg_nd,gg_mi) )
  allocate ( gg_o(gg_mo) )
  allocate ( gg_s(gg_mi) )

  gg_a(1:gg_na) = (/ 4, 7, 9, 10 /)

  gg_b(1:gg_nd,1:gg_ni) = reshape ( (/ &
    -1, -1, &
    -1,  1, &
    -1,  2, &
    -1,  3, &
     1, -1, &
     2,  5, &
     3,  6, &
     5, -1, &
     6,  8, &
     8, -1 /), (/ gg_nd, gg_ni /) )

  gg_f(1:gg_nd,1:gg_ni) = reshape ( (/ &
     5,  2, &
     6,  3, &
     7,  4, &
    -1, -1, &
     8,  6, &
     9,  7, &
    -1, -1, &
    10,  9, &
    -1, -1, &
    -1, -1 /), (/ gg_nd, gg_ni /)  )

  gg_g(1:gg_ni) = (/ &
    0.1D+00, 1.1D+00, 2.2D+00, 3.0D+00, 1.0D+00, &
    2.1D+00, 3.2D+00, 2.0D+00, 3.3D+00, 3.1D+00 /)

  gg_i(1:gg_nd,1:gg_ni) = reshape ( (/ &
    0, 0, &
    0, 1, & 
    0, 2, &
    0, 3, &
    1, 0, &
    1, 1, &
    1, 2, &
    2, 0, &
    2, 1, &
    3, 0 /), (/ gg_nd, gg_ni /) )

  gg_o(1:gg_no) = (/ 1, 2, 3, 5, 6, 8 /)

  gg_s(1:gg_ni) = (/ 0, 0, 0, 1, 0, 0, 1, 0, 1, 1 /)
!
!  Implicit decreasing heap sort on GG_G and GG_A
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Before Heap:'
  write ( *, '(a)' ) '     I     A      G'
  write ( *, '(a)' ) ' '

  do i = 1, gg_na
    write ( *, '(2x,i4,2x,i4,2x,g14.6)' ) i, gg_a(i), gg_g(gg_a(i))
  end do

  call r8vec_indexed_heap_d ( gg_na, gg_g, gg_a )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After Heap:'
  write ( *, '(a)' ) '     I     A      G'
  write ( *, '(a)' ) ' '

  do i = 1, gg_na
    write ( *, '(2x,i4,2x,i4,2x,g14.6)' ) i, gg_a(i), gg_g(gg_a(i))
  end do
!
!  Print out the current data structure:
!
  call sgmgg_print ( gg_ma, gg_mi, gg_mo, gg_na, gg_nd, gg_ni, gg_no, &
    gg_a, gg_b, gg_f, gg_g, gg_i, gg_o, gg_s )
!
!  Now identify the active index with maximum error indicator.
!
  call r8vec_indexed_heap_d_max ( gg_na, gg_g, gg_a, indx_max )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4,a,g14.6)' ) &
    '  Maximum error indicator G(',indx_max, ') = ', gg_g(indx_max)
!
!  We are going to move INDX_MAX from the active to the old list.
!
!  1) Extract INDX_MAX from the active heap.
!     This automatically decrements GG_NA by 1.
!
  call r8vec_indexed_heap_d_extract ( gg_na, gg_g, gg_a, indx_max )
!
!  2) Add INDX_MAX to the old list.
!
  gg_no = gg_no + 1
  gg_o(gg_no) = indx_max
  gg_s(indx_max) = 0
!
!  3) Determine new index sets which are forward neighbors of INDX_MAX.
!     For each index set J found:
!     a) update forward neighbor array of INDX_MAX;
!     b) update backward neighbor array of J;
!     c) add J to A, increment NA
!     d) compute G(J) and add to G heap
!
  do d1 = 1, gg_nd

    d2 = i4_wrap ( d1 + 1, 1, gg_nd )

    nb = gg_b(d2,indx_max)
    indx_nab = gg_f(d1,nb)

    if ( indx_nab == -1 ) then
      cycle
    end if

    write ( *, * ) ' '
    write ( *, * ) '  Extension in direction D1 = ', d1, ' is legal'
    write ( *, * ) '  NB = ', nb
    write ( *, * ) '  INDX_NAB = ', indx_nab

    gg_ni = gg_ni + 1
    indx_new = gg_ni

    gg_i(1:gg_nd,indx_new) = gg_i(1:gg_nd,indx_max)
    gg_i(d1,indx_new) = gg_i(d1,indx_new) + 1

    gg_b(1:gg_nd,indx_new) = -1
    gg_f(1:gg_nd,indx_new) = -1

    gg_b(d1,indx_new) = indx_max
    gg_f(d1,indx_max) = indx_new

    gg_b(d2,indx_new) = indx_nab
    gg_f(d2,indx_nab) = indx_new
!
!  Check whether there are any new neighbor relations in the remaining
!  directions D + 2, D + 3, ..., D + ND - 1.
!
    do d2 = d1 + 2, d1 + gg_nd - 1

      d3 = i4_wrap ( d2, 1, gg_nd )

      nb = gg_b(d3,indx_max)
      indx_nab = gg_f(d1,nb)

      if ( indx_nab == -1 ) then
        cycle
      end if

      gg_b(d3,indx_new) = indx_nab
      gg_f(d3,indx_nab) = indx_new

    end do
    gg_s(indx_new) = 1
    gg_g(indx_new) = r8_uniform_01 ( seed )
    call r8vec_indexed_heap_d_insert ( gg_na, gg_g, gg_a, indx_new )
!
!  Print out the current data structure:
!
  call sgmgg_print ( gg_ma, gg_mi, gg_mo, gg_na, gg_nd, gg_ni, gg_no, &
    gg_a, gg_b, gg_f, gg_g, gg_i, gg_o, gg_s )

  end do

  deallocate ( gg_a )
  deallocate ( gg_b )
  deallocate ( gg_f )
  deallocate ( gg_g )
  deallocate ( gg_i )
  deallocate ( gg_o )
  deallocate ( gg_s )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 simulates a set of incremental coefficient calculations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 August 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: point_num = 8

  integer ( kind = 4 ) point_num2
  integer ( kind = 4 ), dimension ( 2, 8 ) :: sparse_index = reshape ( (/ &
    0, 0, &
    0, 1, &
    0, 2, &
    0, 3, &
    1, 0, &
    1, 1, &
    2, 0, &
    3, 0 /), (/ 2, 8 /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04:';
  write ( *, '(a)' ) '  Simulate incremental coefficient calculations.'
!
!  Generalized grid in 2D.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Generalized grid in 2D:'

  do point_num2 = 1, point_num
    call sandia_sgmgg_coef_naive_test ( dim_num, point_num2, sparse_index )
  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 demonstrates SANDIA_SGMGG_COEF_INC2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 August 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ), parameter :: n1 = 5

  integer ( kind = 4 ) :: c1(n1) = (/ &
    +1, &
    +1, &
    +1, &
    -1, &
    -1 /)
  integer ( kind = 4 ) c3(n1+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4vec_sum
  integer ( kind = 4 ) j
  integer ( kind = 4 ), dimension ( m, n1 ) :: s1 = reshape ( (/ &
    0, 2, &
    1, 1, &
    2, 0, &
    0, 1, &
    1, 0 /), (/ m, n1 /) )
  integer ( kind = 4 ) :: s2(m) = (/ 1, 2 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05:'
  write ( *, '(a)' ) '  Predict new coefficients given candidate index.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Spatial dimension M = ', m
  write ( *, '(a,i4)' ) '  Number of items in active set N1 = ', n1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Index  Coef   Indices'
  write ( *, '(a)' ) ' '

  do j = 1, n1
    write ( *, '(4x,i2,2x,i4,2(2x,i2))' ) j, c1(j), s1(1:m,j)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,2(2x,i2))' ) '  Candidate:', s2(1:m)
!
!  Generalized grid in 2D.
!
  call sandia_sgmgg_coef_inc2 ( m, n1, s1, c1, s2, c3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Index  Coef  Coef'
  write ( *, '(a)' ) '          Old   New'
  write ( *, '(a)' ) ' '
  do i = 1, n1
    write ( *, '(4x,i2,2x,i4,2x,i4)' ) i, c1(i), c3(i)
  end do

  write ( *, '(4x,i2,2x,4x,2x,i4)' ) n1 + 1, c3(n1+1)
  write ( *, '(a)' ) '    --   ----  ----'
  write ( *, '(a,2x,i4,2x,i4)' ) &
    '  Sum:', sum ( c1(1:n1) ), sum ( c3(1:n1+1) )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 demonstrates SANDIA_SGMGG_COEF_INC2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 August 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2
  integer ( kind = 4 ), parameter :: n1 = 5

  integer ( kind = 4 ) :: c1(n1) = (/ &
    +1, &
    +1, &
    +1, &
    -1, &
    -1 /)
  integer ( kind = 4 ) c3(n1+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4vec_sum
  integer ( kind = 4 ) j
  integer ( kind = 4 ), dimension ( m, n1 ) :: s1 = reshape ( (/ &
    2, 0, &
    1, 1, &
    0, 2, &
    1, 0, &
    0, 1 /), (/ m, n1 /) )
  integer ( kind = 4 ) :: s2(m) = (/ 3, 0 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06:'
  write ( *, '(a)' ) '  Predict new coefficients given candidate index.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Spatial dimension M = ', m
  write ( *, '(a,i4)' ) '  Number of items in active set N1 = ', n1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Index  Coef   Indices'
  write ( *, '(a)' ) ' '

  do j = 1, n1
    write ( *, '(4x,i2,2x,i4,2(2x,i2))' ) j, c1(j), s1(1:m,j)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,2(2x,i2))' ) '  Candidate:', s2(1:m)
!
!  Generalized grid in 2D.
!
  call sandia_sgmgg_coef_inc2 ( m, n1, s1, c1, s2, c3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Index  Coef  Coef'
  write ( *, '(a)' ) '          Old   New'
  write ( *, '(a)' ) ' '
  do i = 1, n1
    write ( *, '(4x,i2,2x,i4,2x,i4)' ) i, c1(i), c3(i)
  end do

  write ( *, '(4x,i2,2x,4x,2x,i4)' ) n1 + 1, c3(n1+1)
  write ( *, '(a)' ) '    --   ----  ----'
  write ( *, '(a,2x,i4,2x,i4)' ) &
    '  Sum:', sum ( c1(1:n1) ), sum ( c3(1:n1+1) )

  return
end
subroutine sandia_sgmgg_coef_naive_test ( dim_num, point_num, sparse_index  )

!*****************************************************************************80
!
!! SANDIA_SGMGG_COEF_NAIVE_TEST demonstrates "naive" coefficient computations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) coef(point_num)
  integer ( kind = 4 ) sparse_index(dim_num,point_num)

  call i4mat_transpose_print ( dim_num, point_num, sparse_index, &
    '  SPARSE_INDEX:' )

  call sandia_sgmgg_coef_naive ( dim_num, point_num, sparse_index, coef )

  call i4vec_print ( point_num, coef, '  COEF:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  COEF Sum = ', sum ( coef(1:point_num) )

  return
end
subroutine sandia_sgmgg_neighbor_naive_test ( dim_num, point_num, sparse_index  )

!*****************************************************************************80
!
!! SANDIA_SGMGG_NEIGHBOR_NAIVE_TEST demonstrates "naive" neighbor computations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) neighbor(dim_num,point_num)
  integer ( kind = 4 ) sparse_index(dim_num,point_num)

  call i4mat_transpose_print ( dim_num, point_num, sparse_index, &
    '  SPARSE_INDEX:' )

  call sandia_sgmgg_neighbor_naive ( dim_num, point_num, sparse_index, neighbor )

  call i4mat_transpose_print ( dim_num, point_num, neighbor, &
    '  NEIGHBOR:' )

  return
end
