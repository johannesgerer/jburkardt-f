program main

!*****************************************************************************80
!
!! MAIN is the main program for SVD_TRUNCATED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SVD_TRUNCATED:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the use of the truncated or economy-size'
  write ( *, '(a)' ) '  Singular Value Decomposition (SVD) for cases where'
  write ( *, '(a)' ) '  the sizes of M and N are very different.'

  m = 4
  n = 3
  call svd_truncated_u_test ( m, n )

  m = 3
  n = 4
  call svd_truncated_v_test ( m, n )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SVD_TRUNCATED:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine svd_truncated_u_test ( m, n )

!*****************************************************************************80
!
!! SVD_TRUNCATED_U_TEST tests SVD_TRUNCATED_U.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ), allocatable :: a_save(:,:)
  real ( kind = 8 ) err
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: sn(:)
  real ( kind = 8 ), allocatable :: un(:,:)
  real ( kind = 8 ), allocatable :: v(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SVD_TRUNCATED_U_TEST'
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  allocate ( a(m,n) )
  allocate ( a_save(m,n) )
  allocate ( un(m,n) )
  allocate ( sn(n) )
  allocate ( v(n,n) )

  call random_number ( harvest = a_save(1:m,1:n) )
 
  call r8mat_print ( m, n, a_save, '  A:' )

  a(1:m,1:n) = a_save(1:m,1:n)

  call svd_truncated_u ( m, n, a, un, sn, v )
!
!  Check the factorization by computing A = U * S * V'
!
  do j = 1, n
    do i = 1, m
      a(i,j) = 0.0D+00
      do k = 1, n
        a(i,j) = a(i,j) + un(i,k) * sn(k) * v(k,j)
      end do
    end do
  end do

  err = maxval ( abs ( a(1:m,1:n) - a_save(1:m,1:n) ) )

  write ( *, '(a)' ) ' '

  write ( *, '(a,g14.6)' ) '  Maximum error |A - U*S*V''| = ', err

  call r8mat_print ( m, n, a, '  Recomputed A = U * S * V'':' )

  deallocate ( a )
  deallocate ( a_save )
  deallocate ( sn )
  deallocate ( un )
  deallocate ( v )

  return
end
subroutine svd_truncated_u ( m, n, a, un, sn, v )

!*****************************************************************************80
!
!! SVD_TRUNCATED_U computes the SVD when N <= M.
!
!  Discussion:
!
!    A(mxn) = U(mxm)  * S(mxn)  * V(nxn)'
!           = Un(mxn) * Sn(nxn) * V(nxn)'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix to be decomposed.
!
!    Ouptut, real ( kind = 8 ) UN(M,N), the first N left singular vectors.
!
!    Output, real ( kind = 8 ) SN(N), the first N singular values.
!
!    Output, real ( kind = 8 ) V(N,N), the right singular vectors.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) info
  character jobu
  character jobv
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) lwork
  real ( kind = 8 ) sn(n)
  real ( kind = 8 ) un(m,n)
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ) work(5*n+m)

  if ( m < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVD_TRUNCATED_U - Fatal error!'
    write ( *, '(a)' ) '  Only call this function when N <= M.'
    write ( *, '(a,i8)' ) '  But we have M = ', m
    write ( *, '(a,i8)' ) '  and N = ', n
    stop
  end if

  jobu = 's'
  jobv = 'a'
  lda = m
  ldu = m
  ldv = n
  lwork = 5 * n + m

  call dgesvd ( jobu, jobv, m, n, a, lda, sn, un, ldu, v, ldv, work, lwork, info )

  if ( info == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVD_TRUNCATED_U:'
    write ( *, '(a)' ) '  DGESVD computation was successful.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVD_TRUNCATED_U - Warning!'
    write ( *, '(a,i8)' ) '  DGESVD returned INFO = ', info
  end if

  return
end
subroutine svd_truncated_v_test ( m, n )

!*****************************************************************************80
!
!! SVD_TRUNCATED_V_TEST tests SVD_TRUNCATED_V.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ), allocatable :: a_save(:,:)
  real ( kind = 8 ) err
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: sm(:)
  real ( kind = 8 ), allocatable :: u(:,:)
  real ( kind = 8 ), allocatable :: vm(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SVD_TRUNCATED_V_TEST'
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n

  allocate ( a(m,n) )
  allocate ( a_save(m,n) )
  allocate ( u(m,m) )
  allocate ( sm(m) )
  allocate ( vm(m,n) )

  call random_number ( harvest = a_save(1:m,1:n) )
 
  call r8mat_print ( m, n, a_save, '  A:' )

  a(1:m,1:n) = a_save(1:m,1:n)

  call svd_truncated_v ( m, n, a, u, sm, vm )
!
!  Check the factorization by computing A = U * S * V'
!
  do j = 1, n
    do i = 1, m
      a(i,j) = 0.0D+00
      do k = 1, m
        a(i,j) = a(i,j) + u(i,k) * sm(k) * vm(k,j)
      end do
    end do
  end do

  err = maxval ( abs ( a(1:m,1:n) - a_save(1:m,1:n) ) )

  write ( *, '(a)' ) ' '

  write ( *, '(a,g14.6)' ) '  Maximum error |A - U*S*V''| = ', err

  call r8mat_print ( m, n, a, '  Recomputed A = U * S * V'':' )

  deallocate ( a )
  deallocate ( a_save )
  deallocate ( sm )
  deallocate ( u )
  deallocate ( vm )

  return
end
subroutine svd_truncated_v ( m, n, a, u, sm, vm )

!*****************************************************************************80
!
!! SVD_TRUNCATED_V computes the SVD when M <= N.
!
!  Discussion:
!
!    A(mxn) = U(mxm) * S(mxn)  * V(nxn)'
!           = U(mxm) * Sm(mxm) * Vm(mxn)'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix to be decomposed.
!
!    Ouptut, real ( kind = 8 ) U(M,M), the left singular vectors.
!
!    Output, real ( kind = 8 ) SM(M), the first M singular values.
!
!    Output, real ( kind = 8 ) VM(M,N), the first M right singular vectors.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) info
  character jobu
  character jobv
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) lwork
  real ( kind = 8 ) sm(m)
  real ( kind = 8 ) u(m,m)
  real ( kind = 8 ) vm(m,n)
  real ( kind = 8 ) work(5*m+n)

  if ( n < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVD_TRUNCATED_V - Fatal error!'
    write ( *, '(a)' ) '  Only call this function when M <= N.'
    write ( *, '(a,i8)' ) '  But we have M = ', m
    write ( *, '(a,i8)' ) '  and N = ', n
    stop
  end if

  jobu = 'a'
  jobv = 's'
  lda = m
  ldu = m
  ldv = m
  lwork = 5 * m + n

  call dgesvd ( jobu, jobv, m, n, a, lda, sm, u, ldu, vm, ldv, work, lwork, info )

  if ( info == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVD_TRUNCATED_V:'
    write ( *, '(a)' ) '  DGESVD computation was successful.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVD_TRUNCATED_V - Warning!'
    write ( *, '(a,i8)' ) '  DGESVD returned INFO = ', info
  end if

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
