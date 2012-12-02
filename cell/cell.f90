subroutine i4cvv_iget ( mn, a, m, roff, i, j, aij )

!*****************************************************************************80
!
!! I4CVV_IGET gets item J from row I in an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input, integer ( kind = 4 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) I, the row of the item.
!    1 <= I <= M.
!
!    Input, integer ( kind = 4 ) J, the column of the item.
!    1 <= J <= NR(I).
!
!    Output, integer ( kind = 4 ) AIJ, the value of item A(I,J).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn

  integer ( kind = 4 ) a(mn)
  integer ( kind = 4 ) aij
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) roff(m+1)

  k = roff(i) + j
  aij = a(k)

  return
end
subroutine i4cvv_iinc ( mn, a, m, roff, i, j, daij )

!*****************************************************************************80
!
!! I4CVV_IINC increments item J from row I in an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input/output, integer ( kind = 4 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) I, the row of the item.
!    1 <= I <= M.
!
!    Input, integer ( kind = 4 ) J, the column of the item.
!    1 <= J <= NR(I).
!
!    Input, integer ( kind = 4 ) DAIJ, the increment to item A(I,J).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn

  integer ( kind = 4 ) a(mn)
  integer ( kind = 4 ) daij
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) roff(m+1)

  k = roff(i) + j
  a(k) = a(k) + daij

  return
end
subroutine i4cvv_iset ( mn, a, m, roff, i, j, aij )

!*****************************************************************************80
!
!! I4CVV_ISET sets item J from row I in an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input/output, integer ( kind = 4 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) I, the row of the item.
!    1 <= I <= M.
!
!    Input, integer ( kind = 4 ) J, the column of the item.
!    1 <= J <= NR(I).
!
!    Input, integer ( kind = 4 ) AIJ, the new value of item A(I,J).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn

  integer ( kind = 4 ) a(mn)
  integer ( kind = 4 ) aij
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) roff(m+1)

  k = roff(i) + j
  a(k) = aij

  return
end
subroutine i4cvv_nget ( mn, a, m, roff, nn, in, jn, vn )

!*****************************************************************************80
!
!! I4CVV_NGET gets N items JN(*) from row IN(*) in an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input, integer ( kind = 4 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) NN, the number of items.
!
!    Input, integer ( kind = 4 ) IN(NN), the rows of the items.
!    1 <= IN(*) <= M.
!
!    Input, integer ( kind = 4 ) JN(NN), the columns of the items.
!    1 <= JN(*) <= NR(IN(*)).
!
!    Output, integer ( kind = 4 ) VN(NN), the value of items A(IN(*),JN(*)).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) nn

  integer ( kind = 4 ) a(mn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) in(nn)
  integer ( kind = 4 ) jn(nn)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) roff(m+1)
  integer ( kind = 4 ) vn(nn)

  do i = 1, nn
    k = roff(in(i)) + jn(i)
    vn(i) = a(k)
  end do

  return
end
subroutine i4cvv_ninc ( mn, a, m, roff, nn, in, jn, dvn )

!*****************************************************************************80
!
!! I4CVV_NINC increments items JN(*) from row IN(*) in an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input/output, integer ( kind = 4 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) NN, the number of items.
!
!    Input, integer ( kind = 4 ) IN(NN), the rows of the items.
!    1 <= IN(*) <= M.
!
!    Input, integer ( kind = 4 ) JN(NN), the columns of the items.
!    1 <= JN(*) <= NR(IN(*)).
!
!    Input, integer ( kind = 4 ) DVN(NN), the increments of items A(IN(*),JN(*)).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) nn

  integer ( kind = 4 ) a(mn)
  integer ( kind = 4 ) dvn(nn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) in(nn)
  integer ( kind = 4 ) jn(nn)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) roff(m+1)

  do i = 1, nn
    k = roff(in(i)) + jn(i)
    a(k) = a(k) + dvn(i)
  end do

  return
end
subroutine i4cvv_nset ( mn, a, m, roff, nn, in, jn, vn )

!*****************************************************************************80
!
!! I4CVV_NSET sets items JN(*) from row IN(*) in an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input/output, integer ( kind = 4 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) NN, the number of items.
!
!    Input, integer ( kind = 4 ) IN(NN), the rows of the items.
!    1 <= IN(*) <= M.
!
!    Input, integer ( kind = 4 ) JN(NN), the columns of the items.
!    1 <= JN(*) <= NR(IN(*)).
!
!    Input, integer ( kind = 4 ) VN(NN), the new value of items A(IN(*),JN(*)).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) nn

  integer ( kind = 4 ) a(mn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) in(nn)
  integer ( kind = 4 ) jn(nn)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) roff(m+1)
  integer ( kind = 4 ) vn(nn)

  do i = 1, nn
    k = roff(in(i)) + jn(i)
    a(k) = vn(i)
  end do

  return
end
subroutine i4cvv_offset ( m, nr, roff )

!*****************************************************************************80
!
!! I4CVV_OFFSET determines the row offsets of an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) NR(M), the row sizes.
!
!    Output, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  integer ( kind = 4 ) roff(m+1)
  integer ( kind = 4 ) nr(m)

  roff(1) = 0
  do i = 1, m
    roff(i+1) = roff(i) + nr(i)
  end do

  return
end
subroutine i4cvv_print ( mn, a, m, roff, title )

!*****************************************************************************80
!
!! I4CVV_PRINT prints an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input, integer ( kind = 4 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn

  integer ( kind = 4 ) a(mn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  integer ( kind = 4 ) roff(m+1)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, m

    k1 = roff(i) + 1
    k2 = roff(i+1)

    do klo = k1, k2, 5
      khi = min ( klo + 5 - 1, k2 )
      if ( klo == k1 ) then
        write ( *, '(i5,2x, 10i7)' ) i, a(klo:khi)
      else
        write ( *, '(5x,2x, 10i7)' )    a(klo:khi)
      end if
    end do

  end do

  return
end
subroutine i4cvv_rget ( mn, a, m, roff, i, ai )

!*****************************************************************************80
!
!! I4CVV_RGET gets row I from an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input, integer ( kind = 4 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) I, the row.
!    1 <= I <= M.
!
!    Output, integer ( kind = 4 ) AI(NR(I)), the value of A(I,*).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn

  integer ( kind = 4 ) a(mn)
  integer ( kind = 4 ) ai(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) roff(m+1)

  k1 = roff(i) + 1
  k2 = roff(i+1)
  nv = k2 + 1 - k1
  ai(1:nv) = a(k1:k2)

  return
end
subroutine i4cvv_rinc ( mn, a, m, roff, i, dai )

!*****************************************************************************80
!
!! I4CVV_RINC increments row I in an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input/output, integer ( kind = 4 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) I, the row.
!    1 <= I <= M.
!
!    Input, integer ( kind = 4 ) DAI(NR(I)), the increment for A(I,*).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn

  integer ( kind = 4 ) a(mn)
  integer ( kind = 4 ) dai(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) roff(m+1)

  k1 = roff(i) + 1
  k2 = roff(i+1)
  nv = k2 + 1 - k1
  a(k1:k2) = a(k1:k2) + dai(1:nv)

  return
end
subroutine i4cvv_rset ( mn, a, m, roff, i, ai )

!*****************************************************************************80
!
!! I4CVV_RSET sets row I from an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input/output, integer ( kind = 4 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) I, the row.
!    1 <= I <= M.
!
!    Input, integer ( kind = 4 ) AI(NR(I)), the new value of A(I,*).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn

  integer ( kind = 4 ) a(mn)
  integer ( kind = 4 ) ai(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) roff(m+1)

  k1 = roff(i) + 1
  k2 = roff(i+1)
  nv = k2 + 1 - k1
  a(k1:k2) = ai(1:nv)

  return
end
subroutine i4cvv_size ( m, nr, mn )

!*****************************************************************************80
!
!! I4CVV_SIZE determines the size of an I4CVV.
!
!  Discussion:
!
!    An I4CVV is a "vector of vectors" of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) NR(M), the size of each row.
!
!    Output, integer ( kind = 4 ) MN, the size of the cell array.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) mn
  integer ( kind = 4 ) nr(m)

  mn = sum ( nr(1:m) )

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end
subroutine i4vec_transpose_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Example:
!
!    A = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 /)
!    TITLE = 'My vector:  '
!
!    My vector:
!
!        1    2    3    4    5
!        6    7    8    9   10
!       11
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do ilo = 1, n, 5
    ihi = min ( ilo + 5 - 1, n )
    write ( *, '(5i12)' ) a(ilo:ihi)
  end do

  return
end
subroutine r8cvv_iget ( mn, a, m, roff, i, j, aij )

!*****************************************************************************80
!
!! R8CVV_IGET gets item J from row I in an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input, real ( kind = 8 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) I, the row of the item.
!    1 <= I <= M.
!
!    Input, integer ( kind = 4 ) J, the column of the item.
!    1 <= J <= NR(I).
!
!    Output, real ( kind = 8 ) AIJ, the value of item A(I,J).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn

  real ( kind = 8 ) a(mn)
  real ( kind = 8 ) aij
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) roff(m+1)

  k = roff(i) + j
  aij = a(k)

  return
end
subroutine r8cvv_iinc ( mn, a, m, roff, i, j, daij )

!*****************************************************************************80
!
!! R8CVV_IINC increments item J from row I in an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input/output, real ( kind = 8 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) I, the row of the item.
!    1 <= I <= M.
!
!    Input, integer ( kind = 4 ) J, the column of the item.
!    1 <= J <= NR(I).
!
!    Input, real ( kind = 8 ) DAIJ, the increment to the value of item A(I,J).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn

  real ( kind = 8 ) a(mn)
  real ( kind = 8 ) daij
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) roff(m+1)

  k = roff(i) + j
  a(k) = a(k) + daij

  return
end
subroutine r8cvv_iset ( mn, a, m, roff, i, j, aij )

!*****************************************************************************80
!
!! R8CVV_ISET sets item J from row I in an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input/output, real ( kind = 8 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) I, the row of the item.
!    1 <= I <= M.
!
!    Input, integer ( kind = 4 ) J, the column of the item.
!    1 <= J <= NR(I).
!
!    Input, real ( kind = 8 ) AIJ, the new value of item A(I,J).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn

  real ( kind = 8 ) a(mn)
  real ( kind = 8 ) aij
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) roff(m+1)

  k = roff(i) + j
  a(k) = aij

  return
end
subroutine r8cvv_nget ( mn, a, m, roff, nn, in, jn, vn )

!*****************************************************************************80
!
!! R8CVV_NGET gets N items JN(*) from row IN(*) in an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input, real ( kind = 8 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) NN, the number of items.
!
!    Input, integer ( kind = 4 ) IN(NN), the rows of the items.
!    1 <= IN(*) <= M.
!
!    Input, integer ( kind = 4 ) JN(NN), the columns of the items.
!    1 <= JN(*) <= NR(IN(*)).
!
!    Output, real ( kind = 8 ) VN(NN), the value of items A(IN(*),JN(*)).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) nn

  real ( kind = 8 ) a(mn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) in(nn)
  integer ( kind = 4 ) jn(nn)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) roff(m+1)
  real ( kind = 8 ) vn(nn)

  do i = 1, nn
    k = roff(in(i)) + jn(i)
    vn(i) = a(k)
  end do

  return
end
subroutine r8cvv_ninc ( mn, a, m, roff, nn, in, jn, dvn )

!*****************************************************************************80
!
!! R8CVV_NINC increments items JN(*) from row IN(*) in an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input/output, real ( kind = 8 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) NN, the number of items.
!
!    Input, integer ( kind = 4 ) IN(NN), the rows of the items.
!    1 <= IN(*) <= M.
!
!    Input, integer ( kind = 4 ) JN(NN), the columns of the items.
!    1 <= JN(*) <= NR(IN(*)).
!
!    Input, real ( kind = 8 ) DVN(NN), the increments of items A(IN(*),JN(*)).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) nn

  real ( kind = 8 ) a(mn)
  real ( kind = 8 ) dvn(nn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) in(nn)
  integer ( kind = 4 ) jn(nn)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) roff(m+1)

  do i = 1, nn
    k = roff(in(i)) + jn(i)
    a(k) = a(k) + dvn(i)
  end do

  return
end
subroutine r8cvv_nset ( mn, a, m, roff, nn, in, jn, vn )

!*****************************************************************************80
!
!! R8CVV_NSET sets items JN(*) from row IN(*) in an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input/output, real ( kind = 8 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) NN, the number of items.
!
!    Input, integer ( kind = 4 ) IN(NN), the rows of the items.
!    1 <= IN(*) <= M.
!
!    Input, integer ( kind = 4 ) JN(NN), the columns of the items.
!    1 <= JN(*) <= NR(IN(*)).
!
!    Input, real ( kind = 8 ) VN(NN), the new value of items A(IN(*),JN(*)).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) nn

  real ( kind = 8 ) a(mn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) in(nn)
  integer ( kind = 4 ) jn(nn)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) roff(m+1)
  real ( kind = 8 ) vn(nn)

  do i = 1, nn
    k = roff(in(i)) + jn(i)
    a(k) = vn(i)
  end do

  return
end
subroutine r8cvv_offset ( m, nr, roff )

!*****************************************************************************80
!
!! R8CVV_OFFSET determines the row offsets of an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) NR(M), the row sizes.
!
!    Output, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  integer ( kind = 4 ) roff(m+1)
  integer ( kind = 4 ) nr(m)

  roff(1) = 0
  do i = 1, m
    roff(i+1) = roff(i) + nr(i)
  end do

  return
end
subroutine r8cvv_print ( mn, a, m, roff, title )

!*****************************************************************************80
!
!! R8CVV_PRINT prints an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input, real ( kind = 8 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn

  real ( kind = 8 ) a(mn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  integer ( kind = 4 ) roff(m+1)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, m

    k1 = roff(i) + 1
    k2 = roff(i+1)

    do klo = k1, k2, 5
      khi = min ( klo + 5 - 1, k2 )
      if ( klo == k1 ) then
        write ( *, '(i5,2x, 5g14.6)' ) i, a(klo:khi)
      else
        write ( *, '(5x,2x, 5g14.6)' )    a(klo:khi)
      end if
    end do

  end do

  return
end
subroutine r8cvv_rget ( mn, a, m, roff, i, ai )

!*****************************************************************************80
!
!! R8CVV_RGET gets row I from an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input, real ( kind = 8 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) I, the row.
!    1 <= I <= M.
!
!    Output, real ( kind = 8 ) AI(NR(I)), the value of A(I,*).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn

  real ( kind = 8 ) a(mn)
  real ( kind = 8 ) ai(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) roff(m+1)

  k1 = roff(i) + 1
  k2 = roff(i+1)
  nv = k2 + 1 - k1
  ai(1:nv) = a(k1:k2)

  return
end
subroutine r8cvv_rinc ( mn, a, m, roff, i, dai )

!*****************************************************************************80
!
!! R8CVV_RINC increments row I in an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input/output, real ( kind = 8 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) I, the row.
!    1 <= I <= M.
!
!    Input, real ( kind = 8 ) DAI(NR(I)), the increment for A(I,*).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn

  real ( kind = 8 ) a(mn)
  real ( kind = 8 ) dai(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) roff(m+1)

  k1 = roff(i) + 1
  k2 = roff(i+1)
  nv = k2 + 1 - k1
  a(k1:k2) = a(k1:k2) + dai(1:nv)

  return
end
subroutine r8cvv_rset ( mn, a, m, roff, i, ai )

!*****************************************************************************80
!
!! R8CVV_RSET sets row I from an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MN, the size of the cell array.
!
!    Input/output, real ( kind = 8 ) A(MN), the cell array.
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) ROFF(M+1), the row offsets.
!
!    Input, integer ( kind = 4 ) I, the row.
!    1 <= I <= M.
!
!    Input, real ( kind = 8 ) AI(NR(I)), the new value of A(I,*).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) mn

  real ( kind = 8 ) a(mn)
  real ( kind = 8 ) ai(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) roff(m+1)

  k1 = roff(i) + 1
  k2 = roff(i+1)
  nv = k2 + 1 - k1
  a(k1:k2) = ai(1:nv)

  return
end
subroutine r8cvv_size ( m, nr, mn )

!*****************************************************************************80
!
!! R8CVV_SIZE determines the size of an R8CVV.
!
!  Discussion:
!
!    An R8CVV is a "vector of vectors" of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in the array.
!
!    Input, integer ( kind = 4 ) NR(M), the size of each row.
!
!    Output, integer ( kind = 4 ) MN, the size of the cell array.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) mn
  integer ( kind = 4 ) nr(m)

  mn = sum ( nr(1:m) )

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vec_transpose_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_TRANSPOSE_PRINT prints an R8VEC "transposed".
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Example:
!
!    A = (/ 1.0, 2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0 /)
!    TITLE = 'My vector:  '
!
!    My vector:
!
!        1.0    2.1    3.2    4.3    5.4
!        6.5    7.6    8.7    9.8   10.9
!       11.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do ilo = 1, n, 5
    ihi = min ( ilo + 5 - 1, n )
    write ( *, '(5g14.6)' ) a(ilo:ihi)
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
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
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
