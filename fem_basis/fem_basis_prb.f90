program main

!*****************************************************************************80
!
!! FEM_BASIS_PRB tests FEM_BASIS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM_BASIS_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the FEM_BASIS library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Repeat 1D, 2D, 3D tests but now call FEM_BASIS_MD in each case.
!
  call test04 ( )
  call test05 ( )
  call test06 ( )
!
!  Test the function for triangular prisms.
!
  call test07 ( )
  call test08 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM_BASIS_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests FEM_BASIS_1D
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) lij
  real ( kind = 8 ) r8_fraction
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  FEM_BASIS_1D evaluates an arbitrary'
  write ( *, '(a)' ) '  basis function over an interval.'

  i1 = 2
  j1 = 1
  d = i1 + j1
  x1 = r8_fraction ( i1, d )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I   J        X      L(I,J)(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i2,2x,i2,2x,f10.4,2x,g14.6)' ) i1, j1, x1, 1.0D+00
  write ( *, '(a)' ) ' '
  do i2 = 0, d
    j2 = d - i2
    x2 = r8_fraction ( i2, d )
    call fem_basis_1d ( i1, j1, x2, lij )
    write ( *, '(2x,i2,2x,i2,2x,f10.4,2x,g14.6)' ) &
      i2, j2, x2, lij
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests FEM_BASIS_2D
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) lijk
  real ( kind = 8 ) r8_fraction
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  FEM_BASIS_2D evaluates an arbitrary triangular'
  write ( *, '(a)' ) '  basis function.'

  i1 = 1
  j1 = 0
  k1 = 2
  d = i1 + j1 + k1
  x1 = r8_fraction ( i1, d )
  y1 = r8_fraction ( j1, d )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I   J   K        X           Y      L(I,J,K)(X,Y)'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i2,2x,i2,2x,i2,2x,f10.4,2x,f10.4,2x,g14.6)' ) &
    i1, j1, k1, x1, y1, 1.0D+00
  write ( *, '(a)' ) ' '
  do j2 = 0, d
    do i2 = 0, d - j2
      k2 = d - i2 - j2
      x2 = r8_fraction ( i2, d )
      y2 = r8_fraction ( j2, d )
      call fem_basis_2d ( i1, j1, k1, x2, y2, lijk )
      write ( *, '(2x,i2,2x,i2,2x,i2,2x,f10.4,2x,f10.4,2x,g14.6)' ) &
        i2, j2, k2, x2, y2, lijk
    end do
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests FEM_BASIS_3D
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  real ( kind = 8 ) lijkl
  real ( kind = 8 ) r8_fraction
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) z1
  real ( kind = 8 ) z2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  FEM_BASIS_3D evaluates an arbitrary tetrahedral'
  write ( *, '(a)' ) '  basis function.'

  i1 = 1
  j1 = 0
  k1 = 2
  l1 = 1
  d = i1 + j1 + k1 + l1
  x1 = r8_fraction ( i1, d )
  y1 = r8_fraction ( j1, d )
  z1 = r8_fraction ( k1, d )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I   J   K   L        X           Y           Z      L(I,J,K,L)(X,Y,Z)'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i2,2x,i2,2x,i2,2x,i2,2x,f10.4,2x,f10.4,2x,f10.4,2x,g14.6)' ) &
    i1, j1, k1, l1, x1, y1, z1, 1.0D+00
  write ( *, '(a)' ) ' '
  do k2 = 0, d
    do j2 = 0, d - k2
      do i2 = 0, d - j2 - k2
        l2 = d - i2 - j2 - k2
        x2 = r8_fraction ( i2, d )
        y2 = r8_fraction ( j2, d )
        z2 = r8_fraction ( k2, d )
        call fem_basis_3d ( i1, j1, k1, l1, x2, y2, z2, lijkl )
        write ( *, '(2x,i2,2x,i2,2x,i2,2x,i2,2x,f10.4,2x,f10.4,2x,f10.4,2x,g14.6)' ) &
          i2, j2, k2, l2, x2, y2, z2, lijkl
      end do
    end do
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests FEM_BASIS_MD, repeating TEST01.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 1

  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1(m+1)
  integer ( kind = 4 ) i2(m+1)
  real ( kind = 8 ) l
  integer ( kind = 4 ) p1
  real ( kind = 8 ) r8_fraction
  real ( kind = 8 ) x1(m)
  real ( kind = 8 ) x2(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  FEM_BASIS_MD evaluates an arbitrary'
  write ( *, '(a)' ) '  basis function over an M-dimensional simplex.'

  i1(1) = 2
  i1(2) = 1
  d = sum ( i1(1:m+1) )
  do i = 1, m
    x1(i) = r8_fraction ( i1(i), d )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I   J        X      L(I,J)(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i2,2x,i2,2x,f10.4,2x,g14.6)' ) i1(1:m+1), x1(1:m), 1.0D+00
  write ( *, '(a)' ) ' '
  do p1 =  0, d
    i2(1) = p1
    i2(2) = d - i2(1)
    do i = 1, m
      x2(i) = r8_fraction ( i2(i), d )
    end do
    call fem_basis_md ( m, i1, x2, l )
    write ( *, '(2x,i2,2x,i2,2x,f10.4,2x,g14.6)' ) i2(1:m+1), x2(1:m), l
  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests FEM_BASIS_MD, repeating TEST02.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 2

  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1(m+1)
  integer ( kind = 4 ) i2(m+1)
  real ( kind = 8 ) l
  integer ( kind = 4 ) p1
  integer ( kind = 4 ) p2
  real ( kind = 8 ) r8_fraction
  real ( kind = 8 ) x1(m)
  real ( kind = 8 ) x2(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  FEM_BASIS_MD evaluates an arbitrary'
  write ( *, '(a)' ) '  basis function over an M-dimensional simplex.'

  i1(1) = 1
  i1(2) = 0
  i1(3) = 2
  d = sum ( i1(1:m+1) )
  do i = 1, m
    x1(i) = r8_fraction ( i1(i), d )
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I   J   K        X           Y      L(I,J,K)(X,Y)'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i2,2x,i2,2x,i2,2x,f10.4,2x,f10.4,2x,g14.6)' ) &
    i1(1:m+1), x1(1:m), 1.0D+00
  write ( *, '(a)' ) ' '
  do p2 = 0, d
    i2(2) = p2
    do p1 = 0, d - p2
      i2(1) = p1
      i2(3) = d - sum ( i2(1:m) )
      do i = 1, m
        x2(i) = r8_fraction ( i2(i), d )
      end do
      call fem_basis_md ( m, i1, x2, l )
      write ( *, '(2x,i2,2x,i2,2x,i2,2x,f10.4,2x,f10.4,2x,g14.6)' ) &
        i2(1:m+1), x2(1:m), l
    end do
  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests FEM_BASIS_MD, repeating TEST03.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 3

  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1(m+1)
  integer ( kind = 4 ) i2(m+1)
  real ( kind = 8 ) l
  integer ( kind = 4 ) p1
  integer ( kind = 4 ) p2
  integer ( kind = 4 ) p3
  real ( kind = 8 ) r8_fraction
  real ( kind = 8 ) x1(m)
  real ( kind = 8 ) x2(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  FEM_BASIS_MD evaluates an arbitrary'
  write ( *, '(a)' ) '  basis function over an M-dimensional simplex.'

  i1(1) = 1
  i1(2) = 0
  i1(3) = 2
  i1(4) = 1
  d = sum ( i1(1:m+1) )
  do i = 1, m
    x1(i) = r8_fraction ( i1(i), d )
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I   J   K   L        X           Y           Z      L(I,J,K,L)(X,Y,Z)'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i2,2x,i2,2x,i2,2x,i2,2x,f10.4,2x,f10.4,2x,f10.4,2x,g14.6)' ) &
    i1(1:m+1), x1(1:m), 1.0D+00
  write ( *, '(a)' ) ' '
  do p3 = 0, d
    i2(3) = p3
    do p2 = 0, d - p3
      i2(2) = p2
      do p1 = 0, d - p2 - p3
        i2(1) = p1
        i2(4) = d - sum ( i2(1:m) )
        do i = 1, m
          x2(i) = r8_fraction ( i2(i), d )
        end do
        call fem_basis_md ( m, i1, x2, l )
        write ( *, '(2x,i2,2x,i2,2x,i2,2x,i2,2x,f10.4,2x,f10.4,2x,f10.4,2x,g14.6)' ) &
          i2(1:m+1), x2(1:m), l
      end do
    end do
  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests FEM_BASIS_PRISM_TRIANGLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  integer ( kind = 4 ) di
  integer ( kind = 4 ) dj
  integer ( kind = 4 ) i1(3)
  integer ( kind = 4 ) i2(3)
  integer ( kind = 4 ) i_1
  integer ( kind = 4 ) i_2
  integer ( kind = 4 ) j1(2)
  integer ( kind = 4 ) j2(2)
  integer ( kind = 4 ) j_1
  real ( kind = 8 ) r8_fraction
  real ( kind = 8 ) xyz1(3)
  real ( kind = 8 ) xyz2(3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  FEM_BASIS_PRISM_TRIANGLE evaluates an arbitrary'
  write ( *, '(a)' ) '  basis function over a right triangular prism.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, we generate basis functions which can be'
  write ( *, '(a)' ) '  up to degree 2 in X and Y, and up to degree 2 in Z.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Choose a node N1, define the basis function associated'
  write ( *, '(a)' ) '  with that node, and then evaluate it at all other nodes.'

  i1(1:3) = (/ 2, 0, 0 /)
  di = sum ( i1(1:3) )
  xyz1(1) = r8_fraction ( i1(1), di )
  xyz1(2) = r8_fraction ( i1(2), di )

  j1(1:2) = (/ 1, 1 /)
  dj = sum ( j1(1:2) )
  xyz1(3) = r8_fraction ( j1(1), dj )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  I1  I2  I3  J1  J2        X           Y           Z          B(X,Y,Z)'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i2,2x,i2,2x,i2, &
               2x,i2,2x,i2, &
               2x,f10.4,2x,f10.4,2x,f10.4,2x,g14.6)' ) &
               i1(1:3), j1(1:2), xyz1(1:3), 1.0D+00
  write ( *, '(a)' ) ' '

  do i_1 = 0, di
    i2(1) = i_1
    xyz2(1) = r8_fraction ( i2(1), di )
    do i_2 = 0, di - i2(1)
      i2(2) = i_2
      xyz2(2) = r8_fraction ( i2(2), di )
      i2(3) = di - i2(1) - i2(2)
      do j_1 = 0, dj
        j2(1) = j_1
        j2(2) = dj - j2(1)
        xyz2(3) = r8_fraction ( j2(1), dj )

        call fem_basis_prism_triangle ( i1, j1, xyz2, b )

        write ( *, '(2x,i2,2x,i2,2x,i2, &
                     2x,i2,2x,i2, &
                     2x,f10.4,2x,f10.4,2x,f10.4,2x,g14.6)' ) &
                     i2(1:3), j2(1:2), xyz2(1:3), b
      end do
    end do
  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests FEM_BASIS_PRISM_TRIANGLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  integer ( kind = 4 ) di
  integer ( kind = 4 ) dj
  integer ( kind = 4 ) i1(3)
  integer ( kind = 4 ) i2(3)
  integer ( kind = 4 ) i_1
  integer ( kind = 4 ) i_2
  integer ( kind = 4 ) j1(2)
  integer ( kind = 4 ) j2(2)
  integer ( kind = 4 ) j_1
  real ( kind = 8 ) r8_fraction
  real ( kind = 8 ) xyz1(3)
  real ( kind = 8 ) xyz2(3)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  FEM_BASIS_PRISM_TRIANGLE evaluates an arbitrary'
  write ( *, '(a)' ) '  basis function over a right triangular prism.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, we generate basis functions which can be'
  write ( *, '(a)' ) '  up to degree 3 in X and Y, and up to degree 1 in Z.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Choose a node N1, define the basis function associated'
  write ( *, '(a)' ) '  with that node, and then evaluate it at all other nodes.'

  i1(1:3) = (/ 2, 0, 1 /)
  di = sum ( i1(1:3) )
  xyz1(1) = r8_fraction ( i1(1), di )
  xyz1(2) = r8_fraction ( i1(2), di )

  j1(1:2) = (/ 1, 0 /)
  dj = sum ( j1(1:2) )
  xyz1(3) = r8_fraction ( j1(1), dj )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  I1  I2  I3  J1  J2        X           Y           Z          B(X,Y,Z)'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,i2,2x,i2,2x,i2, &
               2x,i2,2x,i2, &
               2x,f10.4,2x,f10.4,2x,f10.4,2x,g14.6)' ) &
               i1(1:3), j1(1:2), xyz1(1:3), 1.0D+00
  write ( *, '(a)' ) ' '

  do i_1 = 0, di
    i2(1) = i_1
    xyz2(1) = r8_fraction ( i2(1), di )
    do i_2 = 0, di - i2(1)
      i2(2) = i_2
      xyz2(2) = r8_fraction ( i2(2), di )
      i2(3) = di - i2(1) - i2(2)
      do j_1 = 0, dj
        j2(1) = j_1
        j2(2) = dj - j2(1)
        xyz2(3) = r8_fraction ( j2(1), dj )

        call fem_basis_prism_triangle ( i1, j1, xyz2, b )

        write ( *, '(2x,i2,2x,i2,2x,i2, &
                     2x,i2,2x,i2, &
                     2x,f10.4,2x,f10.4,2x,f10.4,2x,g14.6)' ) &
                     i2(1:3), j2(1:2), xyz2(1:3), b
      end do
    end do
  end do

  return
end
