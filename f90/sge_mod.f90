module sge_module

!*****************************************************************************80
!
!! SGE_MODULE is a module of information about a single SGE matrix.
!
!  Discussion:
!
!    Using a module, it is possible to set up an allocatable array
!    that is shared by various routines.  In particular, a user's 
!    program, by including the statement "use sge_module", will be
!    able to examine or alter any of the variables in this module.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2007
!
!   Author:
!
!    John Burkardt
!
  implicit none

  real, save :: det_sge = 0.0E+00
  logical, save :: a_created = .false.
  logical, save :: lu_computed = .false.
  logical, save :: inv_computed = .false.
  logical, save :: det_computed = .false.
  integer, save :: n_sge = -2
  integer, save :: info_sge = -2
  integer, save, allocatable, dimension ( : ) :: pivot_sge
  real, save, allocatable, dimension ( :, : ) :: lu_sge

end
subroutine sge_create ( n, a )

!*****************************************************************************80
!
!! SGE_CREATE stores information about a single SGE matrix.
!
!  Discussion:
!
!    The user calls this routine with input values of N and A.
!    A copy of N is stored.  A copy of A is stored in the LU factor
!    array, but the LU factors are not computed by this routine.
!    Any old information in the module is cancelled.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the array.
!
!    Input, real A(N,N), the array.
!
  use sge_module

  implicit none

  integer n

  real a(n,n)
!
!  Make space for the LU factors and pivot vector.
!
  allocate ( lu_sge(1:n,1:n) )
  allocate ( pivot_sge(1:n) )
!
!  "Remember" the value of N.
!
  n_sge = n
!
!  Note that we have not computed the LU factor yet.
!
  info_sge = -1
!
!  Initialize the LU factor.
!
  lu_sge(1:n,1:n) = a(1:n,1:n)
!
!  Note what we have and have not done.
!
  a_created = .true.
  lu_computed = .false.
  inv_computed = .false.
  det_computed = .false.

  return
end
subroutine sge_delete ( )

!*****************************************************************************80
!
!! SGE_DELETE deletes information about a single SGE matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2007
!
!  Author:
!
!    John Burkardt
!
  use sge_module

  implicit none

  deallocate ( lu_sge )
  deallocate ( pivot_sge )

  n_sge = -1
  info_sge = -1

  a_created = .false.
  lu_computed = .false.
  inv_computed = .false.
  det_computed = .false.

  return
end
function sge_det ( )

!*****************************************************************************80
!
!! SGE_DET computes the determinant of the current matrix.
!
!  Discussion:
!
!    In order to use this routine, the user must first have
!    specified the values of N and A, and passed them to SGE_CREATE,
!    and then called SGE_FA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real SGE_DET, the determinant of the matrix.
!
  use sge_module

  implicit none

  real sge_det
  integer i

  if ( .not. a_created ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_DET - Fatal error!'
    write ( *, '(a)' ) '  The matrix has not been defined yet.'
    write ( *, '(a)' ) '  Call SGE_CREATE first!'
    stop
  end if

  if ( .not. lu_computed ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_DET - Fatal error!'
    write ( *, '(a)' ) '  The matrix has not been factored yet.'
    write ( *, '(a)' ) '  Call SGE_FA first!'
    stop
  end if

  det_sge = 1.0E+00

  do i = 1, n_sge
    det_sge = det_sge * lu_sge(i,i)
  end do

  do i = 1, n_sge
    if ( pivot_sge(i) /= i ) then
      det_sge = - det_sge
    end if
  end do

  det_computed = .true.

  sge_det = det_sge

  return
end
subroutine sge_fa ( )

!*****************************************************************************80
!
!! SGE_FA factors a general matrix.
!
!  Discussion:
!
!    In order to use this routine, the user must first have
!    specified the values of N and A, and passed them to SGE_CREATE.
!
!    SGE_FA is a simplified version of the LINPACK routine SGEFA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
  use sge_module

  implicit none

  integer i
  integer j
  integer k
  integer l

  if ( .not. a_created ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_FA - Fatal error!'
    write ( *, '(a)' ) '  The matrix has not been defined yet.'
    write ( *, '(a)' ) '  Call SGE_CREATE first!'
    stop
  end if

  info_sge = 0

  do k = 1, n_sge-1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k+1, n_sge
      if ( abs ( lu_sge(l,k) ) < abs ( lu_sge(i,k) ) ) then
        l = i
      end if
    end do

    pivot_sge(k) = l
!
!  If the pivot index is zero, the algorithm has failed.
!
    if ( lu_sge(l,k) == 0.0E+00 ) then
      info_sge = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SGE_FA - Fatal error!'
      write ( *, '(a,i6)' ) '  Zero pivot on step ', info_sge
      return
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      call r4_swap ( lu_sge(l,k), lu_sge(k,k) )
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    lu_sge(k+1:n_sge,k) = -lu_sge(k+1:n_sge,k) / lu_sge(k,k)
!
!  Row elimination with column indexing.
!
    do j = k+1, n_sge

      if ( l /= k ) then
        call r4_swap ( lu_sge(l,j), lu_sge(k,j) )
      end if

      lu_sge(k+1:n_sge,j) = lu_sge(k+1:n_sge,j) &
        + lu_sge(k+1:n_sge,k) * lu_sge(k,j)

    end do

  end do

  pivot_sge(n_sge) = n_sge

  if ( lu_sge(n_sge,n_sge) == 0.0E+00 ) then
    info_sge = n_sge
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_FA - Fatal error!'
    write ( *, '(a,i6)' ) '  Zero pivot on step ', info_sge
  end if

  lu_computed = .true.

  return
end
subroutine sge_sl ( b, x )

!*****************************************************************************80
!
!! SGE_SL solves a system factored by SGE_FA.
!
!  Discussion:
!
!    In order to use this routine, the user must first have
!    specified the values of N and A, and passed them to SGE_CREATE,
!    and then called SGE_FA.
!
!    SGE_SL is a simplified version of the LINPACK routine SGESL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real B(N), the right hand side vector.
!
!    Output, real X(N), the solution vector.
!
  use sge_module

  implicit none

  real b(n_sge)
  integer k
  integer l
  real x(n_sge)

  if ( .not. a_created ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_SL - Fatal error!'
    write ( *, '(a)' ) '  The matrix has not been defined yet.'
    write ( *, '(a)' ) '  Call SGE_CREATE first!'
    stop
  end if

  if ( .not. lu_computed ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_SL - Fatal error!'
    write ( *, '(a)' ) '  The matrix has not been factored yet.'
    write ( *, '(a)' ) '  Call SGE_FA first!'
    stop
  end if

  if ( info_sge /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_SL - Fatal error!'
    write ( *, '(a)' ) '  The matrix is singular.'
    stop
  end if

  x(1:n_sge) = b(1:n_sge)
!
!  Solve PL * Y = B.
!
  do k = 1, n_sge-1

    l = pivot_sge(k)

     if ( l /= k ) then
      call r4_swap ( x(l), x(k) )
    end if

    x(k+1:n_sge) = x(k+1:n_sge) + lu_sge(k+1:n_sge,k) * x(k)

  end do
!
!  Solve U * X = Y.
!
  do k = n_sge, 1, -1
    x(k) = x(k) / lu_sge(k,k)
    x(1:k-1) = x(1:k-1) - lu_sge(1:k-1,k) * x(k)
  end do

  return
end
subroutine sge_slt ( b, x )

!*****************************************************************************80
!
!! SGE_SLT solves a transposed system factored by SGE_FA.
!
!  Discussion:
!
!    In order to use this routine, the user must first have
!    specified the values of N and A, and passed them to SGE_CREATE,
!    and then called SGE_FA.
!
!    SGE_SL is a simplified version of the LINPACK routine SGESL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real B(N), the right hand side vector.
!
!    Output, real X(N), the solution vector.
!
  use sge_module

  implicit none

  real b(n_sge)
  integer k
  integer l
  real x(n_sge)

  if ( .not. a_created ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_SL - Fatal error!'
    write ( *, '(a)' ) '  The matrix has not been defined yet.'
    write ( *, '(a)' ) '  Call SGE_CREATE first!'
    stop
  end if

  if ( .not. lu_computed ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_SL - Fatal error!'
    write ( *, '(a)' ) '  The matrix has not been factored yet.'
    write ( *, '(a)' ) '  Call SGE_FA first!'
    stop
  end if

  if ( info_sge /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGE_SL - Fatal error!'
    write ( *, '(a)' ) '  The matrix is singular.'
    stop
  end if

  x(1:n_sge) = b(1:n_sge)
!
!  Solve U' * Y = B.
!
  do k = 1, n_sge
    x(k) = ( x(k) - dot_product ( x(1:k-1), lu_sge(1:k-1,k) ) ) / lu_sge(k,k)
  end do
!
!  Solve ( PL )' * X = Y.
!
  do k = n_sge-1, 1, -1

    x(k) = x(k) + dot_product ( x(k+1:n_sge), lu_sge(k+1:n_sge,k) )

    l = pivot_sge(k)

    if ( l /= k ) then
      call r4_swap ( x(l), x(k) )
    end if

  end do

  return
end
subroutine r4_swap ( x, y )

!*****************************************************************************80
!
!! R4_SWAP swaps two R4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real x
  real y
  real z

  z = x
  x = y
  y = z

  return
end
subroutine r4vec_print ( n, a, title )

!*****************************************************************************80
!
!! R4VEC_PRINT prints an R4VEC.
!
!  Discussion:
!
!    If all the entries are integers, the data is printed
!    in integer format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer n

  real a(n)
  integer i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  if ( all ( a(1:n) == aint ( a(1:n) ) ) ) then
    do i = 1, n
      write ( *, '(i6,i6)' ) i, int ( a(i) )
    end do
  else if ( all ( abs ( a(1:n) ) < 1000000.0E+00 ) ) then
    do i = 1, n
      write ( *, '(i6,f14.6)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(i6,g14.6)' ) i, a(i)
    end do
  end if

  return
end
