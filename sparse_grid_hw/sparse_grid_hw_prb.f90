program main

!*****************************************************************************80
!
!! MAIN is the main program for SPARSE_GRID_HW_PRB.
!
!  Discussion:
!
!    SPARSE_GRID_HW_PRB tests the SPARSE_GRID_HW library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 December 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_HW_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SPARSE_GRID_HW library.'

  call ccu_test ( )
! call ccu_sparse_test ( )
  call get_seq_test ( )
! call gqn_sparse_test ( )
  call gqn_test ( )
! call gqu_sparse_test ( )
  call gqu_test ( )
! call kpn_sparse_test ( )
  call kpn_test ( )
! call kpu_sparse_test ( )
  call kpu_test ( )
  call nwspgr_size_test ( )
! call nwspgr_test ( )
  call order_report ( )
  call pack_rules_test ( )
  call symmetric_sparse_size_test ( )
  call tensor_product_test ( )
  call tensor_product_cell_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_HW_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine ccu_test ( )

!*****************************************************************************80
!
!! CCU_TEST uses the CCU function for 1D quadrature over [0,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) exact
  real ( kind = 8 ) fu_integral
  real ( kind = 8 ), allocatable :: fx(:)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n
  real ( kind = 8 ) q
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CCU_TEST:'
  write ( *, '(a)' ) '  Clenshaw Curtis quadrature over [0,1]:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Level   Nodes    Estimate  Error'
  write ( *, '(a)' ) ' '

  d = 1
  exact = fu_integral ( d );

  do l = 1, 5

    if ( l == 1 ) then
      n = 1
    else
      n = 2 ** ( l - 1 ) + 1
    end if

    allocate ( x(1:n) )
    allocate ( w(1:n) )

    call ccu ( n, x, w )

    allocate ( fx(1:n) )
    call fu_value ( d, n, x, fx )

    q = dot_product ( w(1:n), fx(1:n) )

    e = sqrt ( ( q - exact )**2 ) / exact

    write ( *, '(2x,i2,4x,i6,2x,g14.6,2x,g14.6)' ) l, n, q, e

    deallocate ( fx )
    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine get_seq_test ( )

!*****************************************************************************80
!
!! GET_SEQ_TEST tests GET_SEQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d

  integer ( kind = 4 ), allocatable :: fs(:,:)
  integer ( kind = 4 ) norm
  integer ( kind = 4 ) seq_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GET_SEQ_TEST'
  write ( *, '(a)' ) '  GET_SEQ returns all D-dimensional vectors that sum to NORM.'

  d = 3
  norm = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  D = ', d
  write ( *, '(a,i4)' ) '  NORM = ', norm

  call num_seq ( norm - d, d, seq_num )

  allocate ( fs(1:seq_num,1:d) )

  call get_seq ( d, norm, seq_num, fs )

  call i4mat_print ( seq_num, d, fs, '  The compositions' )

  deallocate ( fs )

  return
end
subroutine gqn_test ( )

!*****************************************************************************80
!
!! GQN_TEST uses the GQN function for 1D quadrature over (-oo,+oo).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) exact
  real ( kind = 8 ) fn_integral
  real ( kind = 8 ), allocatable :: fx(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nh
  real ( kind = 8 ) q
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: wh(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: xh(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GQN_TEST:'
  write ( *, '(a)' ) '  Gauss-Hermite quadrature over (-oo,+oo):'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Level   Nodes    Estimate  Error'
  write ( *, '(a)' ) ' '

  d = 1
  exact = fn_integral ( d );

  do l = 1, 5

    nh = ( l + 1 ) / 2
    allocate ( xh(1:nh) )
    allocate ( wh(1:nh) )

    call gqn ( l, xh, wh )

    if ( mod ( l, 2 ) == 1 ) then
      n = 2 * nh - 1
    else
      n = 2 * nh
    end if

    allocate ( x(1:n) )
    allocate ( w(1:n) )

    if ( mod ( l, 2 ) == 1 ) then
      x(1:nh-1) = - xh(nh:2:-1)
      x(nh:n) = xh(1:nh)
      w(1:nh-1) = wh(nh:2:-1)
      w(nh:n) = wh(1:nh)
    else
      x(1:nh) = - xh(nh:1:-1)
      x(nh+1:n) = xh(1:nh)
      w(1:nh) = wh(nh:1:-1)
      w(nh+1:n) = wh(1:nh)
    end if

    allocate ( fx(1:n) )
    call fn_value ( d, n, x, fx )

    q = dot_product ( w(1:n), fx(1:n) )

    e = sqrt ( ( q - exact )**2 ) / exact

    write ( *, '(2x,i2,4x,i6,2x,g14.6,2x,g14.6)' ) l, n, q, e

    deallocate ( fx )
    deallocate ( w )
    deallocate ( wh )
    deallocate ( x )
    deallocate ( xh )

  end do

  return
end
subroutine gqu_test ( )

!*****************************************************************************80
!
!! GQU_TEST uses the GQU function for 1D quadrature over [0,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) exact
  real ( kind = 8 ) fu_integral
  real ( kind = 8 ), allocatable :: fx(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nh
  real ( kind = 8 ) q
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: wh(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: xh(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GQU_TEST:'
  write ( *, '(a)' ) '  Gauss-Legendre quadrature over [0,1]:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Level   Nodes    Estimate  Error'
  write ( *, '(a)' ) ' '

  d = 1
  exact = fu_integral ( d );

  do l = 1, 5

    nh = ( l + 1 ) / 2
    allocate ( xh(1:nh) )
    allocate ( wh(1:nh) )

    call gqu ( l, xh, wh )

    if ( mod ( l, 2 ) == 1 ) then
      n = 2 * nh - 1
    else
      n = 2 * nh
    end if

    allocate ( x(1:n) )
    allocate ( w(1:n) )

    if ( mod ( l, 2 ) == 1 ) then
      x(1:nh-1) = 1.0D+00 - xh(nh:2:-1)
      x(nh:n) = xh(1:nh)
      w(1:nh-1) = wh(nh:2:-1)
      w(nh:n) = wh(1:nh)
    else
      x(1:nh) = 1.0D+00 - xh(nh:1:-1)
      x(nh+1:n) = xh(1:nh)
      w(1:nh) = wh(nh:1:-1)
      w(nh+1:n) = wh(1:nh)
    end if

    allocate ( fx(1:n) )
    call fu_value ( d, n, x, fx )

    q = dot_product ( w(1:n), fx(1:n) )

    e = sqrt ( ( q - exact )**2 ) / exact

    write ( *, '(2x,i2,4x,i6,2x,g14.6,2x,g14.6)' ) l, n, q, e

    deallocate ( fx )
    deallocate ( w )
    deallocate ( wh )
    deallocate ( x )
    deallocate ( xh )

  end do

  return
end
subroutine kpn_test ( )

!*****************************************************************************80
!
!! KPN_TEST uses the KPN function for 1D quadrature over (-oo,+oo).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) exact
  real ( kind = 8 ) fn_integral
  real ( kind = 8 ), allocatable :: fx(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nh
  real ( kind = 8 ) q
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: wh(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: xh(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'KPN_TEST:'
  write ( *, '(a)' ) '  Kronrod-Patterson-Hermite quadrature over (-oo,+oo):'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Level   Nodes    Estimate  Error'
  write ( *, '(a)' ) ' '

  d = 1
  exact = fn_integral ( d );

  do l = 1, 5

    call kpn_order ( l, n )

    nh = ( n + 1 ) / 2
    allocate ( xh(1:nh) )
    allocate ( wh(1:nh) )

    call kpn ( n, xh, wh )

    allocate ( x(1:n) )
    allocate ( w(1:n) )

    if ( mod ( n, 2 ) == 1 ) then
      x(1:nh-1) = - xh(nh:2:-1)
      x(nh:n) = xh(1:nh)
      w(1:nh-1) = wh(nh:2:-1)
      w(nh:n) = wh(1:nh)
    else
      x(1:nh) = - xh(nh:1:-1)
      x(nh+1:n) = xh(1:nh)
      w(1:nh) = wh(nh:1:-1)
      w(nh+1:n) = wh(1:nh)
    end if

    allocate ( fx(1:n) )
    call fn_value ( d, n, x, fx )

    q = dot_product ( w(1:n), fx(1:n) )

    e = sqrt ( ( q - exact )**2 ) / exact

    write ( *, '(2x,i2,4x,i6,2x,g14.6,2x,g14.6)' ) l, n, q, e

    deallocate ( fx )
    deallocate ( w )
    deallocate ( wh )
    deallocate ( x )
    deallocate ( xh )

  end do

  return
end
subroutine kpu_test ( )

!*****************************************************************************80
!
!! KPU_TEST uses the KPU function for 1D quadrature over [0,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) exact
  real ( kind = 8 ) fu_integral
  real ( kind = 8 ), allocatable :: fx(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nh
  real ( kind = 8 ) q
  real ( kind = 8 ), allocatable :: w(:)
  real ( kind = 8 ), allocatable :: wh(:)
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: xh(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'KPU_TEST:'
  write ( *, '(a)' ) '  Kronrod-Patterson quadrature over [0,1]:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Level   Nodes    Estimate  Error'
  write ( *, '(a)' ) ' '

  d = 1
  exact = fu_integral ( d );

  do l = 1, 5

    call kpu_order ( l, n )

    nh = ( n + 1 ) / 2
    allocate ( xh(1:nh) )
    allocate ( wh(1:nh) )

    call kpu ( n, xh, wh )

    allocate ( x(1:n) )
    allocate ( w(1:n) )

    if ( mod ( n, 2 ) == 1 ) then
      x(1:nh-1) = 1.0D+00 - xh(nh:2:-1)
      x(nh:n) = xh(1:nh)
      w(1:nh-1) = wh(nh:2:-1)
      w(nh:n) = wh(1:nh)
    else
      x(1:nh) = 1.0D+00 - xh(nh:1:-1)
      x(nh+1:n) = xh(1:nh)
      w(1:nh) = wh(nh:1:-1)
      w(nh+1:n) = wh(1:nh)
    end if

    allocate ( fx(1:n) )
    call fu_value ( d, n, x, fx )

    q = dot_product ( w(1:n), fx(1:n) )

    e = sqrt ( ( q - exact )**2 ) / exact

    write ( *, '(2x,i2,4x,i6,2x,g14.6,2x,g14.6)' ) l, n, q, e

    deallocate ( fx )
    deallocate ( w )
    deallocate ( wh )
    deallocate ( x )
    deallocate ( xh )

  end do

  return
end
subroutine nwspgr_size_test ( )

!*****************************************************************************80
!
!! NWSPGR_SIZE_TEST tests NWSPGR_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 December 2012
!
!  Author:
!
!    John Burkardt.
!
  implicit none

  external ccu_order
  integer ( kind = 4 ) dim
  external gqn_order
  external gqu_order
  integer ( kind = 4 ) k
  external kpn_order
  external kpu_order
  integer ( kind = 4 ) r_size

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'NWSPGR_SIZE_TEST:'
  write ( *, '(a)' ) '  NWSPGR_SIZE returns the size of a sparse grid, based on either:'
  write ( *, '(a)' ) '  one of the built-in 1D rules, or a family of 1D rules'
  write ( *, '(a)' ) '  supplied by the user.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Kronrod-Patterson, [0,1], Dim 2, Level 3, Symmetric'
  write ( *, '(a)' ) ''
  call nwspgr_size ( kpu_order, 2, 3, r_size )
  write ( *, '(a,i6)' ) '  Full          ', r_size

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Kronrod-Patterson, (-oo,+oo), Dim 2, Level 3, Symmetric'
  write ( *, '(a)' ) ''
  call nwspgr_size ( kpn_order, 2, 3, r_size )
  write ( *, '(a,i6)' ) '  Full          ', r_size

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Gauss-Legendre, [0,1], Dim 2, Level 3, Symmetric'
  write ( *, '(a)' ) ''
  call nwspgr_size ( gqu_order, 2, 3, r_size )
  write ( *, '(a,i6)' ) '  Full          ', r_size

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Gauss Hermite, (-oo,+oo), [0,1], Dim 2, Level 3, Symmetric'
  write ( *, '(a)' ) ''
  call nwspgr_size ( gqn_order, 2, 3, r_size )
  write ( *, '(a,i6)' ) '  Full          ', r_size

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Clenshaw Curtis, [-1,+1], [0,1], Dim 2, Level 3, Unsymmetric'
  write ( *, '(a)' ) ''
  call nwspgr_size ( ccu_order, 2, 3, r_size )
  write ( *, '(a,i6)' ) '  Full          ', r_size
!
!  Do a table.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Dimension / Level table for Clenshaw Curtis Uncompressed'
  write ( *, '(a)' ) ''
  write ( *, '(a)', advance = 'no' ) ' Dim: '
  do dim = 1, 10
    write ( *, '(2x,i6)', advance = 'no' ) dim
  end do
  write ( *, '(a)', advance = 'yes' ) ''
  write ( *, '(a)' ) 'Level'
  do k = 1, 5
    write ( *, '(2x,i2,2x)', advance = 'no' ) k
    do dim = 1, 10
      call nwspgr_size ( ccu_order, dim, k, r_size )
      write ( *, '(2x,i6)', advance = 'no' ) r_size
    end do
    write ( *, '(a)', advance = 'yes' ) ''
  end do

  return
end
subroutine order_report ( )

!*****************************************************************************80
!
!! ORDER_REPORT reports on the order of each family of rules.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) ap
  integer ( kind = 4 ) k
  integer ( kind = 4 ), dimension ( 5 ) :: kpn_order = (/ &
    1, 3, 9, 19, 35 /)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) o
  integer ( kind = 4 ) rp

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ORDER_REPORT'
  write ( *, '(a)' ) '  For each family of rules, report:'
  write ( *, '(a)' ) '  L,  the level index,'
  write ( *, '(a)' ) '  RP, the required polynomial precision,'
  write ( *, '(a)' ) '  AP, the actual polynomial precision,'
  write ( *, '(a)' ) '  O,  the rule order (number of points).'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  GQN family'
  write ( *, '(a)' ) '  Gauss quadrature, exponential weight, (-oo,+oo)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   L  RP  AP   O'
  write ( *, '(a)' ) ' '

  do l = 1, 25
    rp = 2 * l - 1
    o = l
    ap = 2 * o - 1
    write ( *, '(2x,i2,2x,i2,2x,i2,2x,i2)' ) l, rp, ap, o
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  GQU family'
  write ( *, '(a)' ) '  Gauss quadrature, unit weight, [0,1]'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   L  RP  AP   O'
  write ( *, '(a)' ) ' '

  do l = 1, 25
    rp = 2 * l - 1
    o = l
    ap = 2 * o - 1
    write ( *, '(2x,i2,2x,i2,2x,i2,2x,i2)' ) l, rp, ap, o
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  KPN family'
  write ( *, '(a)' ) '  Gauss-Kronrod-Patterson quadrature, exponential weight, (-oo,+oo)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   L  RP  AP   O'
  write ( *, '(a)' ) ' '

  k = 1
  o = 1
  ap = 1

  do l = 1, 25

    rp = 2 * l - 1

    do while ( ap < rp )

      if ( k == 5 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  No higher order rule is available!'
        exit
      end if
!
!  Can we use a simple rule?
!
      if ( rp < kpn_order(k+1) ) then
        o = rp
        ap = rp
!
!  Otherwise, move to next higher rule.
!
      else
        k = k + 1
        ap = 2 * kpn_order(k) - kpn_order(k-1)
        o = kpn_order(k)
      end if

    end do

    write ( *, '(2x,i2,2x,i2,2x,i2,2x,i2)' ) l, rp, ap, o

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  KPU family'
  write ( *, '(a)' ) '  Gauss-Kronrod-Patterson quadrature, unit weight, [0,1]'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   L  RP  AP   O'
  write ( *, '(a)' ) ' '

  do l = 1, 25
    rp = 2 * l - 1
    o = 1
    ap = 1
    do while ( ap < rp )
      o = 2 * ( o + 1 ) - 1
      ap = ( 3 * o + 1 ) / 2
    end do
    write ( *, '(2x,i2,2x,i2,2x,i2,2x,i2)' ) l, rp, ap, o
  end do

  return
end
subroutine pack_rules_test ( )

!*****************************************************************************80
!
!! PACK_RULES_TEST tests RULES_1D_SIZE and RULES_1D_SET.
!
!  Discussion:
!
!    These two procedures create packed arrays of information defining
!    the 1D rules needed for the sparse grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 June 2012
!
!  Author:
!
!    John Burkardt.
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  external ccu
  external ccu_order
  integer ( kind = 4 ) r1d(k+1)
  integer ( kind = 4 ) r1d_size
  logical sym
  real ( kind = 8 ), allocatable :: w1d(:)
  real ( kind = 8 ), allocatable :: x1d(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PACK_RULES_TEST'
  write ( *, '(a)' ) '  Given a sparse grid level K, the code must collect'
  write ( *, '(a)' ) '  the nodes and weights for the 1D rule of levels 1 through K.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RULES_1D_SIZE determines the size of the packed vectors.'
  write ( *, '(a)' ) '  RULES_1D_SET creates the packed vectors.'

  sym = .true.

  call rules_1d_size ( k, sym, ccu_order, r1d_size, r1d )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  R1D_SIZE = ', r1d_size

  call i4vec_print ( k + 1, r1d, '  R1D pointer vector:' )
!
!  Now call CCU to set values...
!
  allocate ( x1d(1:r1d_size) )
  allocate ( w1d(1:r1d_size) )

  call rules_1d_set ( k, sym, ccu_order, ccu, r1d_size, r1d, x1d, w1d )

  call r8vecs_print ( k, r1d, r1d_size, x1d, '  X vectors' )
  call r8vecs_print ( k, r1d, r1d_size, w1d, '  W vectors' )

  return
end
subroutine symmetric_sparse_size_test ( )

!*****************************************************************************80
!
!! SYMMETRIC_SPARSE_SIZE_TEST tests SYMMETRIC_SPARSE_SIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2012
!
!  Author:
!
!    Original MATLAB version by Florian Heiss, Viktor Winschel.
!    FORTRAN90 version by John Burkardt.
!
!  Local parameters:
!
!    Local, integer D, the spatial dimension.
!
!    Local, integer MAXK, the maximum level to check.
!

  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  integer ( kind = 4 ) dim
  integer ( kind = 4 ), dimension ( test_num ) :: dim_test = (/ 5, 5, 3 /)
  real ( kind = 8 ), dimension ( 6, 5 ) :: nodes1 = reshape ( (/ &
   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, &
   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, &
   0.0, 0.0, 0.0, 1.0, 0.0, 0.0, &
   0.0, 0.0, 1.0, 0.0, 0.0, 0.0, &
   0.0, 1.0, 0.0, 0.0, 0.0, 0.0 /), (/ 6, 5 /) )
  real ( kind = 8 ), dimension ( 21, 5 ) :: nodes2 = reshape ( (/ &
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.73205, &
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.73205, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.73205, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 1.0, 1.0, 1.73205, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, &
    0.0, 1.0, 1.73205, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 /), (/ 21, 5 /) )
  real ( kind = 8 ), dimension ( 23, 3 ) :: nodes3 = reshape ( (/ &
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.741964, 1.0, &
    1.0, 1.0, 1.0, 1.0, 1.0, 1.73205, 1.73205, 1.73205, 2.33441, &
    0.0, 0.0, 0.0, 0.0, 0.0, 0.741964, 1.0, 1.0, 1.0, 1.73205, 1.73205, 2.33441, &
    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.73205, 0.0, 0.0, 1.0, 0.0, &
    0.0, 0.741964, 1.0, 1.73205, 2.33441, 0.0, 0.0, 1.0, 1.73205, 0.0, 1.0, 0.0, &
    0.0, 0.0, 1.0, 1.73205, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0 /), (/ 23, 3 /) )
  integer ( kind = 4 ) r
  integer ( kind = 4 ), dimension ( test_num ) :: r_test = (/ 6, 21, 23 /)
  integer ( kind = 4 ) r2
  integer ( kind = 4 ) test
  real ( kind = 8 ) x0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SYMMETRIC_SPARSE_SIZE_TEST'
  write ( *, '(a)' ) '  Given a symmetric sparse grid rule represented only by'
  write ( *, '(a)' ) '  the points with positive values, determine the total number'
  write ( *, '(a)' ) '  of points in the grid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For dimension DIM, we report '
  write ( *, '(a)' ) '  R, the number of points in the positive orthant, and '
  write ( *, '(a)' ) '  R2, the total number of points.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       DIM         R        R2'
  write ( *, '(a)' ) ' '

  x0 = 0.0

  do test = 1, test_num

    r = r_test(test)
    dim = dim_test(test)
    if ( test == 1 ) then
      call symmetric_sparse_size ( r, dim, nodes1, x0, r2 )
    else if ( test == 2 ) then
      call symmetric_sparse_size ( r, dim, nodes2, x0, r2 )
    else if ( test == 3 ) then
      call symmetric_sparse_size ( r, dim, nodes3, x0, r2 )
    end if

    write ( *, '(2x,i8,2x,i8,2x,i8)' ) dim, r, r2

  end do

  return
end
subroutine tensor_product_test ( )

!*****************************************************************************80
!
!! TENSOR_PRODUCT_TEST tests TENSOR_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order1 = 2
  integer ( kind = 4 ), parameter :: order2 = 3
  integer ( kind = 4 ), parameter :: order3 = 2

  integer ( kind = 4 ) d
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i4vec_product
  integer ( kind = 4 ) i4vec_sum
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n1d
  integer ( kind = 4 ), allocatable :: order1d(:)
  real ( kind = 8 ), allocatable :: w1d(:)
  real ( kind = 8 ), allocatable :: wnd(:)
  real ( kind = 8 ) :: w1_1d(order1) = (/ 1.0D+00, 1.0D+00 /)
  real ( kind = 8 ) :: w2_1d(order2) = (/ 0.25D+00, 0.50D+00, 0.25D+00 /)
  real ( kind = 8 ) :: w3_1d(order3) = (/ 2.50D+00, 2.50D+00 /)
  real ( kind = 8 ) :: x1_1d(order1) = (/ -1.0D+00, +1.0D+00 /)
  real ( kind = 8 ) :: x2_1d(order2) = (/ 2.0D+00, 2.5D+00, 3.0D+00 /)
  real ( kind = 8 ) :: x3_1d(order3) = (/ 10.0D+00, 15.0D+00 /)
  real ( kind = 8 ), allocatable :: x1d(:)
  real ( kind = 8 ), allocatable :: xnd(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TENSOR_PRODUCT_TEST:'
  write ( *, '(a)' ) '  Given a sequence of 1D quadrature rules, construct the'
  write ( *, '(a)' ) '  tensor product rule.'
!
!  1D rule.
!
  d = 1
  allocate ( order1d(1:d) )

  order1d(1) = order1

  n1d = i4vec_sum ( d, order1d )
  allocate ( x1d(1:n1d) )
  allocate ( w1d(1:n1d) )

  n = i4vec_product ( d, order1d )
  allocate ( xnd(1:d,1:n) )
  allocate ( wnd(1:n) )

  i1 = 1
  i2 = order1
  x1d(i1:i2) = x1_1d(1:order1)
  w1d(i1:i2) = w1_1d(1:order1)
 
  call tensor_product ( d, order1d,n1d, x1d, w1d, n, xnd, wnd )
  
  call quad_rule_print ( d, n, xnd, wnd, '  A 1D rule over [-1,+1]:' )

  deallocate ( order1d )
  deallocate ( w1d )
  deallocate ( wnd )
  deallocate ( x1d )
  deallocate ( xnd )
!
!  2D rule.
!
  d = 2
  allocate ( order1d(1:d) )

  order1d(1:2) = (/ order1, order2 /)

  n1d = i4vec_sum ( d, order1d )
  allocate ( x1d(1:n1d) )
  allocate ( w1d(1:n1d) )

  n = i4vec_product ( d, order1d )
  allocate ( xnd(1:d,1:n) )
  allocate ( wnd(1:n) )

  i1 = 1
  i2 = order1
  x1d(i1:i2) = x1_1d(1:order1)
  w1d(i1:i2) = w1_1d(1:order1)
  i1 = i2 + 1
  i2 = i2 + order2
  x1d(i1:i2) = x2_1d(1:order2)
  w1d(i1:i2) = w2_1d(1:order2)

  call tensor_product ( d, order1d, n1d, x1d, w1d, n, xnd, wnd )
  
  call quad_rule_print ( d, n, xnd, wnd, '  A 2D rule over [-1,+1] x [2.0,3.0]:' )

  deallocate ( order1d )
  deallocate ( w1d )
  deallocate ( wnd )
  deallocate ( x1d )
  deallocate ( xnd )
!
!  3D rule.
!
  d = 3
  allocate ( order1d(1:d) )

  order1d(1:3) = (/ order1, order2, order3 /)

  n1d = i4vec_sum ( d, order1d )
  allocate ( x1d(1:n1d) )
  allocate ( w1d(1:n1d) )

  n = i4vec_product ( d, order1d )
  allocate ( xnd(1:d,1:n) )
  allocate ( wnd(1:n) )

  i1 = 1
  i2 = order1
  x1d(i1:i2) = x1_1d(1:order1)
  w1d(i1:i2) = w1_1d(1:order1)
  i1 = i2 + 1
  i2 = i2 + order2
  x1d(i1:i2) = x2_1d(1:order2)
  w1d(i1:i2) = w2_1d(1:order2)
  i1 = i2 + 1
  i2 = i2 + order3
  x1d(i1:i2) = x3_1d(1:order3)
  w1d(i1:i2) = w3_1d(1:order3)

  call tensor_product ( d, order1d, n1d, x1d, w1d, n, xnd, wnd )

  call quad_rule_print ( d, n, xnd, wnd, &
    '  A 3D rule over [-1,+1] x [2.0,3.0] x [10.0,15.0]:' )

  deallocate ( order1d )
  deallocate ( w1d )
  deallocate ( wnd )
  deallocate ( x1d )
  deallocate ( xnd )

  return
end
subroutine tensor_product_cell_test ( )

!*****************************************************************************80
!
!! TENSOR_PRODUCT_CELL_TEST tests TENSOR_PRODUCT_CELL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 December 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: order1 = 2
  integer ( kind = 4 ), parameter :: order2 = 3
  integer ( kind = 4 ), parameter :: order3 = 2

  integer ( kind = 4 ) d
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i4vec_product
  integer ( kind = 4 ) i4vec_sum
  integer ( kind = 4 ) n1d
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) np
  integer ( kind = 4 ) :: nr(3) = (/ 2, 3, 2 /)
  integer ( kind = 4 ), allocatable :: order1d(:)
  integer ( kind = 4 ) roff(4)
  real ( kind = 8 ), allocatable :: w1d(:)
  real ( kind = 8 ), allocatable :: wc(:)
  real ( kind = 8 ), allocatable :: wp(:)
  real ( kind = 8 ) :: w1_1d(order1) = (/ 1.0D+00, 1.0D+00 /)
  real ( kind = 8 ) :: w2_1d(order2) = (/ 0.25D+00, 0.50D+00, 0.25D+00 /)
  real ( kind = 8 ) :: w3_1d(order3) = (/ 2.50D+00, 2.50D+00 /)
  real ( kind = 8 ) :: x1_1d(order1) = (/ -1.0D+00, +1.0D+00 /)
  real ( kind = 8 ) :: x2_1d(order2) = (/ 2.0D+00, 2.5D+00, 3.0D+00 /)
  real ( kind = 8 ) :: x3_1d(order3) = (/ 10.0D+00, 15.0D+00 /)
  real ( kind = 8 ), allocatable :: x1d(:)
  real ( kind = 8 ), allocatable :: xc(:)
  real ( kind = 8 ), allocatable :: xp(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TENSOR_PRODUCT_TEST_CELL:'
  write ( *, '(a)' ) '  Given a set of 1D quadrature rules stored in a cell array,'
  write ( *, '(a)' ) '  construct the tensor product rule.'
!
!  We can construct ROFF once and for all.
!
  call r8cvv_offset ( 3, nr, roff )
!
!  1D rule.
!
  d = 1
  nc = sum ( nr(1:d) )
  allocate ( xc(1:nc) )
  call r8cvv_rset ( nc, xc, d, roff, 1, x1_1d )
  allocate ( wc(1:nc) )
  call r8cvv_rset ( nc, wc, d, roff, 1, w1_1d )
  np = product ( nr(1:d) )
  allocate( xp(1:d,1:np) )
  allocate( wp(1:np) )

  call tensor_product_cell ( nc, xc, wc, d, nr, roff, np, xp, wp )

  call quad_rule_print ( d, np, xp, wp, '  A 1D rule over [-1,+1]:' )

  deallocate ( wc )
  deallocate ( wp )
  deallocate ( xc )
  deallocate ( xp )
!
!  2D rule.
!
  d = 2
  nc = sum ( nr(1:d) )
  allocate ( xc(1:nc) )
  call r8cvv_rset ( nc, xc, d, roff, 1, x1_1d )
  call r8cvv_rset ( nc, xc, d, roff, 2, x2_1d )
  allocate ( wc(1:nc) )
  call r8cvv_rset ( nc, wc, d, roff, 1, w1_1d )
  call r8cvv_rset ( nc, wc, d, roff, 2, w2_1d )
  np = product ( nr(1:d) )
  allocate( xp(1:d,1:np) )
  allocate( wp(1:np) )

  call tensor_product_cell ( nc, xc, wc, d, nr, roff, np, xp, wp )

  call quad_rule_print ( d, np, xp, wp, '  A 1D rule over [-1,+1]:' )

  deallocate ( wc )
  deallocate ( wp )
  deallocate ( xc )
  deallocate ( xp )
!
!  3D rule.
!
  d = 3
  nc = sum ( nr(1:d) )
  allocate ( xc(1:nc) )
  call r8cvv_rset ( nc, xc, d, roff, 1, x1_1d )
  call r8cvv_rset ( nc, xc, d, roff, 2, x2_1d )
  call r8cvv_rset ( nc, xc, d, roff, 3, x3_1d )
  allocate ( wc(1:nc) )
  call r8cvv_rset ( nc, wc, d, roff, 1, w1_1d )
  call r8cvv_rset ( nc, wc, d, roff, 2, w2_1d )
  call r8cvv_rset ( nc, wc, d, roff, 3, w3_1d )
  np = product ( nr(1:d) )
  allocate( xp(1:d,1:np) )
  allocate( wp(1:np) )

  call tensor_product_cell ( nc, xc, wc, d, nr, roff, np, xp, wp )

  call quad_rule_print ( d, np, xp, wp, '  A 1D rule over [-1,+1]:' )

  deallocate ( wc )
  deallocate ( wp )
  deallocate ( xc )
  deallocate ( xp )

  return
end

