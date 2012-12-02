subroutine amask ( nrow, ncol, a, ja, ia, jmask, imask, c, jc, ic, iw, &
  nzmax, ierr )

!*****************************************************************************80
!
!! AMASK extracts a sparse matrix from a masked input matrix.
!
!  Discussion:
!
!    The routine looks at the positions defined by MASK, JMASK and IMASK.
!
!    The algorithm is "in place": C, JC, IC can be the same as
!    A, JA, IA.  
!
!  Modified:
!
!    08 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Reference:
!
!    Youcef Saad,
!    Sparsekit: a basic tool kit for sparse matrix computations,
!    Technical Report, Computer Science Department,
!    University of Minnesota, June 1994
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NCOL, the column dimension of the matrix.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, integer JMASK(*), IMASK((NROW+1), defining mask (pattern only) 
!    stored in compressed sparse row format.
!
!    Input, integer NZMAX, the length of arrays C and JC.
!
!    Output, C, JC, IC, the output matrix in Compressed Sparse Row format.
!
!    Workspace, logical IW(NCOL).
!
!    Input, integer NZMAX, the dimension of C.
!
!    Output, integer IERR, serving as error message.
!    ierr = 1  means normal return
!    ierr > 1 means that amask stopped when processing
!    row number ierr, because there was not enough space in
!    c, jc according to the value of nzmax.
!
  implicit none

  integer ncol
  integer nrow
  integer nzmax

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) c(nzmax)
  integer ia(nrow+1)
  integer ic(nrow+1)
  integer ierr
  integer ii
  integer imask(nrow+1)
  logical iw(ncol)
  integer j
  integer ja(*)
  integer jc(nzmax)
  integer jmask(*)
  integer k
  integer k1
  integer k2
  integer len

  ierr = 0
  len = 0

  iw(1:ncol) = .false.
!
!  Unpack the mask for row II in IW.
!
  do ii = 1, nrow
!
!  Save pointer in order to be able to do things in place.
!
    do k = imask(ii), imask(ii+1)-1
      iw(jmask(k)) = .true.
    end do
!
!  Add unmasked elemnts of row II.
!
    k1 = ia(ii)
    k2 = ia(ii+1)-1
    ic(ii) = len+1

    do k = k1, k2
      j = ja(k)
      if ( iw(j) ) then
        len = len + 1
        if ( nzmax < len ) then
          ierr = ii
          return
        end if
        jc(len) = j
        c(len) = a(k)
      end if
    end do

    do k = imask(ii), imask(ii+1)-1
      iw(jmask(k)) = .false.
    end do

  end do

  ic(nrow+1) = len + 1

  return
end
subroutine amub ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, &
  iw, ierr )

!*****************************************************************************80
!
!! AMUB performs the matrix product C = A * B.
!
!  Discussion:
!
!    The column dimension of B is not needed.
!
!  Modified:
!
!    08 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NCOL, the column dimension of the matrix.
!
!    Input, integer JOB, job indicator.  When JOB = 0, only the structure
!    is computed, that is, the arrays JC and IC, but the real values
!    are ignored.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, b, jb, ib, matrix B in compressed sparse row format.
!
!    Input, integer NZMAX, the length of the arrays c and jc.
!    The routine will stop if the result matrix C  has a number
!    of elements that exceeds exceeds NZMAX.
!
! on return:
!
! c,
! jc,
! ic    = resulting matrix C in compressed sparse row sparse format.
!
! ierr      = integer. serving as error message.
!         ierr = 0 means normal return,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!
!  iw      = integer work array of length equal to the number of
!         columns in A.
!
  implicit none

  integer ncol
  integer nrow
  integer nzmax

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) c(nzmax)
  integer ia(nrow+1)
  integer ib(ncol+1)
  integer ic(ncol+1)
  integer ierr
  integer ii
  integer iw(ncol)
  integer ja(*)
  integer jb(*)
  integer jc(nzmax)
  integer jcol
  integer jj
  integer job
  integer jpos
  integer k
  integer ka
  integer kb
  integer len
  real ( kind = 8 ) scal
  logical values

  values = ( job /= 0 )
  len = 0
  ic(1) = 1
  ierr = 0
!
!  Initialize IW.
!
  iw(1:ncol) = 0

  do ii = 1, nrow
!
!  Row I.
!
    do ka = ia(ii), ia(ii+1)-1

      if ( values ) then
        scal = a(ka)
      end if

      jj = ja(ka)

      do kb = ib(jj), ib(jj+1)-1

           jcol = jb(kb)
           jpos = iw(jcol)

           if ( jpos == 0 ) then
              len = len + 1
              if ( nzmax < len ) then
                 ierr = ii
                 return
              end if
              jc(len) = jcol
              iw(jcol)= len
              if ( values ) then
                c(len) = scal * b(kb)
              end if
           else
              if ( values ) then
                c(jpos) = c(jpos) + scal * b(kb)
              end if
           end if

         end do

    end do

    do k = ic(ii), len
      iw(jc(k)) = 0
    end do

    ic(ii+1) = len + 1

  end do

  return
end
subroutine amubdg ( nrow, ncol, ncolb, ja, ia, jb, ib, ndegr, nnz, iw )

!*****************************************************************************80
!
!! AMUBDG gets the number of nonzero elements in each row of A * B.
!
!  Discussion:
!
!    The routine also computes the total number of nonzero elements in A * B.
!
!    Method: A' * A = sum [over i = 1, nrow]  a(i)^T a(i)
!    where a(i) = i-th row of  A.  We must be careful not to add  the
!    elements already accounted for.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix A.
!
!    Input, integer NCOL, the column dimension of the matrix A,
!    (and the row dimension of B).
!
!    Input, integer NCOLB, the column dimension of the matrix B.
!
!    Input, ja, ia= row structure of input matrix A: ja = column indices of
!    the nonzero elements of A stored by rows.
!    ia = pointer to beginning of each row in ja.
!
!    Input, jb, ib, the row structure of input matrix B: jb = column indices of
!    the nonzero elements of A stored by rows.
!    ib is a pointer to beginning of each row in jb.
!
!    Output, integer NDEGR(NROW), contains the degrees (the number of 
!    nonzeros in each row of the matrix A * B.
!
!    Output, integer NNZ, the number of nonzero elements found in A * B.
!
!   Workspace, integer IW(NCOLB).
!
  implicit none

  integer ncol
  integer ncolb
  integer nrow

  integer ia(nrow+1)
  integer ib(ncol+1)
  integer ii
  integer iw(ncolb)
  integer j
  integer ja(*)
  integer jb(*)
  integer jc
  integer jr
  integer k
  integer last
  integer ldg
  integer ndegr(nrow)
  integer nnz

  iw(1:ncolb) = 0
  ndegr(1:nrow) = 0

  do ii = 1, nrow
!
!  For each row of A.
!
    ldg = 0
!
!  End-of-linked list.
!
    last = -1

    do j = ia(ii), ia(ii+1)-1
!
!  Row number to be added.
!
        jr = ja(j)

        do k = ib(jr), ib(jr+1)-1
           jc = jb(k)
!
!  Add one element to the linked list.
!
           if ( iw(jc) == 0 ) then
              ldg = ldg + 1
              iw(jc) = last
              last = jc
           end if

         end do

    end do

    ndegr(ii) = ldg
!
!  Reset IW to zero.
!
    do k = 1, ldg
      j = iw(last)
      iw(last) = 0
      last = j
     end do

  end do

  nnz = sum ( ndegr(1:nrow) )

  return
end
subroutine amudia ( nrow, job, a, ja, ia, diag, b, jb, ib )

!*****************************************************************************80
!
!! AMUDIA performs the matrix by matrix product B = A * Diag  (in place)
!
!  Discussion:
!
!    The column dimension of A is not needed.
!    The algorithm is "in place", so B can take the place of A.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer JOB, job indicator. Job=0 means get array b only
!    job = 1 means get b, and the integer arrays ib, jb.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real DIAG(NROW), the diagonal matrix stored as a vector.
!
!    Output, B(*), JB(*), IB(NROW+1), the resulting matrix B in 
!    compressed sparse row sparse format.
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) diag(nrow)
  integer ia(nrow+1)
  integer ib(nrow+1)
  integer ii
  integer ja(*)
  integer jb(*)
  integer job
  integer k
  integer k1
  integer k2

  do ii = 1, nrow
!
!  Scale each element.
!
    k1 = ia(ii)
    k2 = ia(ii+1) - 1
    do k = k1, k2
      b(k) = a(k) * diag(ja(k))
    end do

  end do

  if ( job == 0 ) then
    return
  end if

  ib(1) = ia(1)
  do ii = 1, nrow
    ib(ii) = ia(ii)
    do k = ia(ii), ia(ii+1)-1
      jb(k) = ja(k)
    end do
  end do

  return
end
subroutine amux ( n, x, y, a, ja, ia )

!*****************************************************************************80
!
!! AMUX multiplies a CSR matrix A times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix by a vector using the dot product form.
!    Matrix A is stored in compressed sparse row storage.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, real X(*), and array of length equal to the column dimension 
!    of A.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real Y(N), the product A * X.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  integer i
  integer ia(*)
  integer ja(*)
  integer k
  real ( kind = 8 ) t
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(n)

  do i = 1, n
!
!  Compute the inner product of row I with vector X.
!
    t = 0.0D+00
    do k = ia(i), ia(i+1)-1
      t = t + a(k) * x(ja(k))
    end do

    y(i) = t

  end do

  return
end
subroutine amuxd ( n, x, y, diag, ndiag, idiag, ioff )

!*****************************************************************************80
!
!! AMUXD multiplies a DIA matrix times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix by a vector when the original matrix 
!    is stored in the DIA diagonal storage format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, real X(*), array of length equal to the column dimension of
!    the A matrix.
!
!    Output, real Y(N), the product A * X.
!
!    Input, real DIAG(NDIAG,IDIAG), the diagonals.
!
!    Input, integer NDIAG, the first dimension of array adiag as declared in
!    the calling program.
!
!    Input, integer IDIAG, the number of diagonals in the matrix.
!
!    Input, integer IOFF(IDIAG), the offsets of the diagonals of the matrix:
!    diag(i,k) contains the element a(i,i+ioff(k)) of the matrix.
!
  implicit none

  integer idiag
  integer n
  integer ndiag

  real ( kind = 8 ) diag(ndiag,idiag)
  integer i1
  integer i2
  integer io
  integer ioff(idiag)
  integer j
  integer k
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  y(1:n) = 0.0D+00

  do j = 1, idiag
    io = ioff(j)
    i1 = max ( 1, 1 - io )
    i2 = min ( n, n - io )
    do k = i1, i2
      y(k) = y(k) + diag(k,j) * x(k+io)
    end do
  end do

  return
end
subroutine amuxe ( n, x, y, na, ncol, a, ja )

!*****************************************************************************80
!
!! AMUXE multiplies an ELL matrix times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix by a vector, where the matrix is stored
!    in the ELL Ellpack/Itpack sparse format.
!
!  Modified:
!
!    09 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, real X(*), array of length equal to the column dimension of
!    the A matrix.
!
!    Input, integer NA, the first dimension of arrays A and JA
!    as declared by the calling program.
!
!    Input, integer NCOL, the number of active columns in array a.
!    (i.e., the number of generalized diagonals in matrix.)
!
!    a, ja = the real and integer arrays of the Ellpack/Itpack format
!    (a(i,k),k = 1,ncol contains the elements of row i in matrix
!    ja(i,k),k = 1,ncol contains their column numbers)
!
!    Output, real Y(N), the product A * X.
!
  implicit none

  integer n
  integer na
  integer ncol

  real ( kind = 8 ) a(na,ncol)
  integer i
  integer j
  integer ja(na,ncol)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  y(1:n) = 0.0D+00

  do j = 1, ncol
    do i = 1, n
      y(i) = y(i) + a(i,j) * x(ja(i,j))
    end do
  end do

  return
end
subroutine amuxj ( n, x, y, jdiag, a, ja, ia )

!*****************************************************************************80
!
!! AMUXJ multiplies a JAD matrix times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix A times a vector, where A is
!    stored in JAD Jagged-Diagonal storage format.
!
!    Permutation related to the JAD format is not performed.
!    This can be done by:
!
!      call permvec ( n, y, y, iperm )
!
!    after the call to AMUXJ.  Here IPERM is the permutation produced
!    by CSRJAD.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, real X(*), an array of length equal to the column dimension of
!    the A matrix.
!
!    Input, integer JDIAG, the number of jagged-diagonals in the data
!    structure.
!
! a      = real array containing the jagged diagonals of A stored
!          in succession (in decreasing lengths)
!
! j      = integer array containing the column indices of the
!          corresponding elements in a.
!
! ia     = integer array containing the lengths of the  jagged diagonals
!
!    Output, real Y(N), the product A*X.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  integer ia(*)
  integer ii
  integer j
  integer ja(*)
  integer jdiag
  integer k1
  integer len
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  y(1:n) = 0.0D+00

  do ii = 1, jdiag
    k1 = ia(ii) - 1
    len = ia(ii+1) - k1 - 1
    do j = 1, len
      y(j) = y(j) + a(k1+j) * x(ja(k1+j))
    end do
  end do

  return
end
subroutine aplb ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, &
  iw, ierr )

!*****************************************************************************80
!
!! APLB performs the CSR matrix sum C = A + B.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of A and B.
!
!    Input, integer NCOL, the column dimension of A and B.
!
!    Input, integer JOB.  When JOB = 0, only the structure
!    (i.e. the arrays jc, ic) is computed and the
!    real values are ignored.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! b,
! jb,
! ib      =  Matrix B in compressed sparse row format.
!
! nzmax      = integer. The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!
! c,
! jc,
! ic      = resulting matrix C in compressed sparse row sparse format.
!
! ierr      = integer. serving as error message.
!         ierr = 0 means normal return,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!
! iw      = integer work array of length equal to the number of
!         columns in A.
!
  implicit none

  integer ncol
  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) c(*)
  integer ia(nrow+1)
  integer ib(nrow+1)
  integer ic(nrow+1)
  integer ierr
  integer ii
  integer iw(ncol)
  integer ja(*)
  integer jb(*)
  integer jc(*)
  integer jcol
  integer job
  integer jpos
  integer k
  integer ka
  integer kb
  integer len
  integer nzmax
  logical values

  values = ( job /= 0 )
  ierr = 0
  len = 0
  ic(1) = 1
  iw(1:ncol) = 0

  do ii = 1, nrow
!
!  Row I.
!
     do ka = ia(ii), ia(ii+1)-1

        len = len + 1
        jcol = ja(ka)

        if ( nzmax < len ) then
          ierr = ii
          return
        end if

        jc(len) = jcol
        if ( values ) then
          c(len) = a(ka)
        end if
        iw(jcol) = len
     end do

     do kb = ib(ii), ib(ii+1)-1

        jcol = jb(kb)
        jpos = iw(jcol)

        if ( jpos == 0 ) then

           len = len + 1

           if ( nzmax < len ) then
             ierr = ii
             return
           end if

           jc(len) = jcol
           if ( values ) then
             c(len) = b(kb)
           end if
           iw(jcol)= len
        else
           if ( values ) then
             c(jpos) = c(jpos) + b(kb)
           end if
        end if

     end do

     do k = ic(ii), len
       iw(jc(k)) = 0
     end do

     ic(ii+1) = len+1
  end do

  return
end
subroutine aplb1 ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, &
  nzmax, ierr )

!*****************************************************************************80
!
!! APLB1 performs the sum C = A + B for sorted CSR matrices.
!
!  Discussion:
!
!    The difference between this routine and APLB is that here the 
!    resulting matrix is such that the elements of each row are sorted,
!    with increasing column indices in each row, provided the original
!    matrices are sorted in the same way.
!
!    This routine will not work if either of the two input matrices is 
!    not sorted.
!
!  Modified:
!
!    11 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of A and B.
!
!    Input, integer NCOL, the column dimension of A and B.
!
!    Input, integer JOB.  When JOB = 0, only the structure
!    (i.e. the arrays jc, ic) is computed and the
!    real values are ignored.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format with entries sorted.
!
! b,
! jb,
! ib      =  Matrix B in compressed sparse row format with entries sorted
!        ascendly in each row
!
! nzmax      = integer. The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!
! c,
! jc,
! ic      = resulting matrix C in compressed sparse row sparse format
!         with entries sorted ascendly in each row.
!
! ierr      = integer. serving as error message.
!         ierr = 0 means normal return,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) c(*)
  integer i
  integer ia(nrow+1)
  integer ib(nrow+1)
  integer ic(nrow+1)
  integer ierr
  integer j1
  integer j2
  integer ja(*)
  integer jb(*)
  integer jc(*)
  integer job
  integer ka
  integer kamax
  integer kb
  integer kbmax
  integer kc
  integer ncol
  integer nzmax
  logical values

  values = ( job /= 0 )
  ierr = 0
  kc = 1
  ic(1) = kc

  do i = 1, nrow

    ka = ia(i)
    kb = ib(i)
    kamax = ia(i+1) - 1
    kbmax = ib(i+1) - 1

    do

      if ( ka <= kamax ) then
        j1 = ja(ka)
      else
        j1 = ncol + 1
      end if

      if ( kb <= kbmax ) then
        j2 = jb(kb)
      else
        j2 = ncol + 1
      end if
!
!  Three cases
!
      if ( j1 == j2 ) then
        if ( values ) then
          c(kc) = a(ka) + b(kb)
        end if
        jc(kc) = j1
        ka = ka + 1
        kb = kb + 1
        kc = kc + 1
      else if ( j1 < j2 ) then
        jc(kc) = j1
        if ( values ) then
          c(kc) = a(ka)
        end if
        ka = ka + 1
        kc = kc + 1
      else if ( j2 < j1 ) then
        jc(kc) = j2
        if ( values ) then
          c(kc) = b(kb)
        end if
        kb = kb + 1
        kc = kc + 1
      end if

      if ( nzmax < kc ) then
        ierr = i
        return
      end if

      if ( kamax < ka .and. kbmax < kb ) then
        exit
      end if

     end do

     ic(i+1) = kc

  end do

  return
end
subroutine aplbdg ( nrow, ncol, ja, ia, jb, ib, ndegr, nnz, iw )

!*****************************************************************************80
!
!! APLBDG gets the number of nonzero elements in each row of A + B.
!
!  Discussion:
!
!    It also reports the total number of nonzero elements in A + B.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of A and B.
!
!    Input, integer NCOL, the column dimension of A and B.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, b, jb, ib, matrix B in compressed sparse row format.
!
!    Output, integer NDEGR(NROW), the number of nonzeros in each row 
!    of the matrix A + B.
!
!    Output, integer NNZ, the total number of nonzero elements found 
!    in A * B.
!
!    Workspace, integer IW(NCOL).
!
  implicit none

  integer ncol
  integer nrow

  integer ia(nrow+1)
  integer ib(nrow+1)
  integer ii
  integer iw(ncol)
  integer j
  integer ja(*)
  integer jb(*)
  integer jc
  integer jr
  integer k
  integer last
  integer ldg
  integer ndegr(nrow)
  integer nnz

  iw(1:ncol) = 0

  ndegr(1:nrow) = 0

  do ii = 1, nrow

     ldg = 0
!
!  End-of-linked list.
!
     last = -1
!
!  Row of A.
!
     do j = ia(ii), ia(ii+1)-1
        jr = ja(j)
!
!  Add element to the linked list.
!
        ldg = ldg + 1
        iw(jr) = last
        last = jr
     end do
!
!  Row of B.
!
     do j = ib(ii), ib(ii+1)-1

        jc = jb(j)
!
!  Add one element to the linked list.
!
        if ( iw(jc) == 0 ) then
           ldg = ldg + 1
           iw(jc) = last
           last = jc
        end if

     end do
!
!  Done with row II.
!
     ndegr(ii) = ldg
!
!  Reset IW to zero.
!
     do k = 1, ldg
        j = iw(last)
        iw(last) = 0
        last = j
     end do

  end do

  nnz = sum ( ndegr(1:nrow) )

  return
end
subroutine apldia ( nrow, job, a, ja, ia, diag, b, jb, ib, iw )

!*****************************************************************************80
!
!! APLDIA adds a diagonal matrix to a general sparse matrix:  B = A + Diag.
!
!  Discussion:
!
!    The column dimension of A is not needed.
!
!    The algorithm is in place (b, jb, ib, can be the same as
!    a, ja, ia, on entry). See comments for parameter job.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer JOB, job indicator. Job=0 means get array b only
!    (i.e. assume that a has already been copied into array b,
!    or that algorithm is used in place. ) For all practical
!    puposes enter job=0 for an in-place call and job=1 otherwise.
!    In case there are missing diagonal elements in A,
!    then the option job =0 will be ignored, since the algorithm
!    must modify the data structure (i.e. jb, ib) in this
!    situation.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real DIAG(NROW), a diagonal matrix.
!
! on return:
!
! b,
! jb,
! ib      = resulting matrix B in compressed sparse row sparse format.
!
!
! iw    = integer work array of length n. On return iw will
!         contain  the positions of the diagonal entries in the
!         output matrix. (i.e., a(iw(k)), ja(iw(k)), k = 1,...n,
!         are the values/column indices of the diagonal elements
!         of the output matrix. ).
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) diag(nrow)
  integer ia(nrow+1)
  integer ib(nrow+1)
  integer icount
  integer ii
  integer iw(*)
  integer j
  integer ja(*)
  integer jb(*)
  integer job
  integer k
  integer k1
  integer k2
  integer ko
  integer nnz
  logical test
!
!  Copy integer arrays into B's data structure if required.
!
  if ( job /= 0 ) then
    nnz = ia(nrow+1)-1
    jb(1:nnz) = ja(1:nnz)
    ib(1:nrow+1) = ia(1:nrow+1)
  end if
!
!  Get positions of diagonal elements in data structure.
!
  call diapos ( nrow, ja, ia, iw )
!
!  Count number of holes in diagonal and add DIAG elements to
!  valid diagonal entries.
!
  icount = 0

  do j = 1, nrow

     if ( iw(j) == 0 ) then
        icount = icount + 1
     else
        b(iw(j)) = a(iw(j)) + diag(j)
     end if

  end do
!
!  If no diagonal elements to insert, return.
!
  if ( icount == 0 ) then
    return
  end if
!
!  Shift the nonzero elements if needed, to allow for created
!  diagonal elements.
!
  ko = ib(nrow+1) + icount
!
!  Copy rows backward.
!
  do ii = nrow, 1, -1
!
!  Go through row II.
!
     k1 = ib(ii)
     k2 = ib(ii+1) - 1
     ib(ii+1) = ko
     test = ( iw(ii) == 0 )

     do k = k2, k1, -1

        j = jb(k)

         if ( test .and. j < ii ) then
           test = .false.
           ko = ko - 1
           b(ko) = diag(ii)
           jb(ko) = ii
           iw(ii) = ko
        end if

        ko = ko - 1
        b(ko) = a(k)
        jb(ko) = j

      end do
!
!  The diagonal element has not been added yet.
!
     if ( test ) then
        ko = ko - 1
        b(ko) = diag(ii)
        jb(ko) = ii
        iw(ii) = ko
     end if

  end do

  ib(1) = ko

  return
end
subroutine aplsb ( nrow, ncol, a, ja, ia, s, b, jb, ib, c, jc, ic, nzmax, &
  iw, ierr )

!*****************************************************************************80
!
!! APLSB performs the matrix linear combination C = A + s * B.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NCOL, the column dimension of the matrix B.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real S, scalar factor for B.
!
! b,
! jb,
! ib      =  Matrix B in compressed sparse row format.
!
! nzmax      = integer. The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!
! c,
! jc,
! ic      = resulting matrix C in compressed sparse row sparse format.
!
! ierr      = integer. serving as error message.
!         ierr = 0 means normal return,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!
! iw      = integer work array of length equal to the number of
!         columns in A.
!
  implicit none

  integer ncol
  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) c(*)
  integer ia(nrow+1)
  integer ib(nrow+1)
  integer ic(nrow+1)
  integer ierr
  integer ii
  integer iw(ncol)
  integer ja(*)
  integer jb(*)
  integer jc(*)
  integer jcol
  integer jpos
  integer k
  integer ka
  integer kb
  integer len
  integer nzmax
  real ( kind = 8 ) s

  ierr = 0
  len = 0
  ic(1) = 1

  iw(1:ncol) = 0

  do ii = 1, nrow
!
!  Row I.
!
     do ka = ia(ii), ia(ii+1)-1

        len = len + 1
        jcol = ja(ka)

        if ( nzmax < len ) then
          ierr = ii
          return
        end if

        jc(len) = jcol
        c(len) = a(ka)
        iw(jcol)= len

     end do

     do kb = ib(ii), ib(ii+1)-1

        jcol = jb(kb)
        jpos = iw(jcol)
        if ( jpos == 0 ) then
           len = len + 1
           if ( nzmax < len ) then
             ierr = ii
             return
           end if
           jc(len) = jcol
           c(len) = s * b(kb)
           iw(jcol)= len
        else
           c(jpos) = c(jpos) + s * b(kb)
        end if

     end do

     do k = ic(ii), len
        iw(jc(k)) = 0
     end do

     ic(ii+1) = len + 1
 
  end do

  return
end
subroutine aplsb1 ( nrow, ncol, a, ja, ia, s, b, jb, ib, c, jc, ic, &
  nzmax, ierr )

!*****************************************************************************80
!
!! APLSB1 performs the operation C = A + s * B for sorted CSR matrices.
!
!  Discussion:
!
!    The difference with aplsb is that the resulting matrix is such that
!    the elements of each row are sorted with increasing column indices in
!    each row, provided the original matrices are sorted in the same way.
!
!    This will not work if any of the two input matrices is not sorted
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of A and B.
!
!    Input, integer NCOL, the column dimension of A and B.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format with entries sorted.
!
!    Input, real S, a scale factor for B.
!
! b,
! jb,
! ib      =  Matrix B in compressed sparse row format with entries sorted
!        ascendly in each row
!
! nzmax      = integer. The  length of the arrays c and jc.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!
! c,
! jc,
! ic      = resulting matrix C in compressed sparse row sparse format
!         with entries sorted ascendly in each row.
!
! ierr      = integer. serving as error message.
!         ierr = 0 means normal return,
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) c(*)
  integer i
  integer ia(nrow+1)
  integer ib(nrow+1)
  integer ic(nrow+1)
  integer ierr
  integer j1
  integer j2
  integer ja(*)
  integer jb(*)
  integer jc(*)
  integer ka
  integer kamax
  integer kb
  integer kbmax
  integer kc
  integer ncol
  integer nzmax
  real ( kind = 8 ) s

  ierr = 0
  kc = 1
  ic(1) = kc

  do i = 1, nrow

    ka = ia(i)
    kb = ib(i)
    kamax = ia(i+1) - 1
    kbmax = ib(i+1) - 1

    do

      if ( ka <= kamax ) then
        j1 = ja(ka)
      else
        j1 = ncol + 1
      end if

      if ( kb <= kbmax ) then
        j2 = jb(kb)
      else
        j2 = ncol + 1
      end if
!
!  Three cases.
!
      if ( j1 == j2 ) then
        c(kc) = a(ka) + s * b(kb)
        jc(kc) = j1
        ka = ka + 1
        kb = kb + 1
        kc = kc + 1
      else if ( j1 < j2 ) then
        jc(kc) = j1
        c(kc) = a(ka)
        ka = ka + 1
        kc = kc + 1
      else if ( j2 < j1 ) then
        jc(kc) = j2
        c(kc) = s * b(kb)
        kb = kb + 1
        kc = kc + 1
      end if

      if ( nzmax < kc ) then
        ierr = i
        return
      end if

      if ( kamax < ka .and. kbmax < kb ) then
        exit
      end if

    end do

    ic(i+1) = kc

  end do

  return
end
subroutine aplsbt ( nrow, ncol, a, ja, ia, s, b, jb, ib, c, jc, ic, nzmax, &
  iw, ierr )

!*****************************************************************************80
!
!! APLSBT performs the matrix sum C = A + B'.
!
!  Discussion:
!
!    It is important to note that here all of three arrays c, ic,
!    and jc are assumed to be of length nnz(c). This is because
!    the matrix is internally converted to coordinate "COO" format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of A.  This must also be
!    the column dimension of B.
!
!    Input, integer NCOL, the column dimension of the matrix A.
!    This must also be the row dimension of B.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real S, the scalar factor for B.
!
! b,
! jb,
! ib      =  Matrix B in compressed sparse row format.
!
! nzmax      = integer. The  length of the arrays c, jc, and ic.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!
! c,
! jc,
! ic      = resulting matrix C in compressed sparse row format.
!
! ierr      = integer. serving as error message.
!         ierr = 0 means normal return.
!         ierr = -1 means that nzmax was < either the number of
!         nonzero elements of A or the number of nonzero elements in B.
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!
! iw      = integer work array of length equal to the number of
!         columns in A.
!
  implicit none

  integer ncol
  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) c(*)
  integer ia(nrow+1)
  integer ib(ncol+1)
  integer ic(*)
  integer ierr
  integer ii
  integer ipos
  integer iw(ncol)
  integer ja(*)
  integer jb(*)
  integer jc(*)
  integer jcol
  integer jpos
  integer k
  integer ka
  integer len
  integer ljob
  integer nnza
  integer nnzb
  integer nzmax
  real ( kind = 8 ) s

  ierr = 0
  iw(1:ncol) = 0

  nnza = ia(nrow+1) - 1
  nnzb = ib(ncol+1) - 1
  len = nnzb

  if ( nzmax < nnzb .or. nzmax < nnza ) then
    ierr = -1
    return
  end if
!
!  Transpose matrix B into C.
!
  ljob = 1
  ipos = 1
  call csrcsc ( ncol, ljob, ipos, b, jb, ib, c, jc, ic )
  c(1:len) = c(1:len) * s
!
!  The main loop. 
!  Add rows from 1 through NROW.
!
  do ii = 1, nrow
!
!  IW is used as a system to recognize whether there
!  was a nonzero element in C.
!
    do k = ic(ii), ic(ii+1)-1
      iw(jc(k)) = k
    end do

    do ka = ia(ii), ia(ii+1)-1

      jcol = ja(ka)
      jpos = iw(jcol)
!
!  If fill-in, append in coordinate format to matrix.
!
      if ( jpos == 0 ) then

        len = len + 1

        if ( nzmax < len ) then
          ierr = ii
          return
        end if

        jc(len) = jcol
        ic(len) = ii
        c(len) = a(ka)

      else
!
!  else, do addition.
!
        c(jpos) = c(jpos) + a(ka)
      end if

    end do

    do k = ic(ii), ic(ii+1)-1
      iw(jc(k)) = 0
    end do

  end do
!
!  Convert matrix without fill-ins into coordinate format.
!
  ljob = 3

  call csrcoo ( nrow, ljob, nnzb, c, jc, ic, nnzb, c, ic, jc, ierr )

  if ( ierr /= 0 ) then
    ierr = -ierr
  end if
!
!  Convert the whole thing back to CSR format.
!
  ljob = 1

  call coocsr_inplace ( nrow, len, 1, c, jc, ic, iw )

  return
end
subroutine aplsca ( nrow, a, ja, ia, scal, iw )

!*****************************************************************************80
!
!! APLSCA adds a scalar to the diagonal entries of a sparse matrix A :=A + s I
!
!  Discussion:
!
!    The column dimension of A is not needed.
!
!    important: the matrix A may be expanded slightly to allow for
!    additions of nonzero elements to previously nonexisting diagonals.
!    The is no checking as to whether there is enough space appended
!    to the arrays a and ja. if not sure allow for n additional
!    elements.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real SCAL, a scalar to be added to the diagonal entries.
!
! on return:
!
!
! a,
! ja,
! ia      = matrix A with diagonal elements shifted (or created).
!
! iw    = integer work array of length n. On return iw will
!         contain  the positions of the diagonal entries in the
!         output matrix. (i.e., a(iw(k)), ja(iw(k)), k = 1,...n,
!         are the values/column indices of the diagonal elements
!         of the output matrix. ).
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  integer ia(nrow+1)
  integer icount
  integer ii
  integer iw(*)
  integer j
  integer ja(*)
  integer k
  integer k1
  integer k2
  integer ko
  real ( kind = 8 ) scal
  logical test

  call diapos ( nrow, ja, ia, iw )
  icount = 0

  do j = 1, nrow

     if ( iw(j) == 0 ) then
        icount = icount + 1
     else
        a(iw(j)) = a(iw(j)) + scal
     end if

  end do
!
!  If no diagonal elements to insert in data structure, return.
!
  if ( icount == 0 ) then
    return
  end if
!
!  Shift the nonzero elements if needed, to allow for created
!  diagonal elements.
!
  ko = ia(nrow+1) + icount
!
!  Copy rows backward.
!
  do ii = nrow, 1, -1
!
!  Go through row II.
!
     k1 = ia(ii)
     k2 = ia(ii+1) - 1
     ia(ii+1) = ko
     test = ( iw(ii) == 0 )

     do k = k2, k1, -1

        j = ja(k)

        if ( test .and. j < ii ) then
           test = .false.
           ko = ko - 1
           a(ko) = scal
           ja(ko) = ii
           iw(ii) = ko
        end if

        ko = ko - 1
        a(ko) = a(k)
        ja(ko) = j

    end do
!
!  The diagonal element has not been added yet.
!
     if ( test ) then
        ko = ko - 1
        a(ko) = scal
        ja(ko) = ii
        iw(ii) = ko
     end if

  end do

  ia(1) = ko

  return
end
subroutine apmbt ( nrow, ncol, job, a, ja, ia, b, jb, ib, c, jc, ic, nzmax, &
  iw, ierr )

!*****************************************************************************80
!
!! APMBT performs the matrix sum  C = A + B' or C = A - B'.
!
!  Discussion:
!
!    It is important to note that here all of three arrays c, ic,
!    and jc are assumed to be of length nnz(c).  This is because
!    the matrix is internally converted to coordinate "COO" format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of A, which must also be
!    the column dimension of B.
!
!    Input, integer NCOL, the column dimension of the matrix A, which
!    must also be the row dimension of B.
!
! job      = integer. if job = -1, apmbt will compute C= A - transp(B)
!         (structure + values)
!         if job == 1, it will compute C=A+transp(A)
!         (structure+ values)
!         if job == 0, it will compute the structure of
!         C= A+/-transp(B) only (ignoring all real values).
!         any other value of job will be treated as  job=1
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! b,
! jb,
! ib      =  Matrix B in compressed sparse row format.
!
! nzmax      = integer. The  length of the arrays c, jc, and ic.
!         amub will stop if the result matrix C  has a number
!         of elements that exceeds exceeds nzmax. See ierr.
!
! on return:
!
! c,
! jc,
! ic      = resulting matrix C in compressed sparse row format.
!
! ierr      = integer. serving as error message.
!         ierr = 0 means normal return.
!         ierr = -1 means that nzmax was < either the number of
!         nonzero elements of A or the number of nonzero elements in B.
!         ierr > 0 means that amub stopped while computing the
!         i-th row  of C with i = ierr, because the number
!         of elements in C exceeds nzmax.
!
! work arrays:
!
! iw      = integer work array of length equal to the number of
!         columns in A.
!
  implicit none

  integer ncol
  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) c(*)
  integer ia(nrow+1)
  integer ib(ncol+1)
  integer ic(*)
  integer ierr
  integer ii
  integer ipos
  integer iw(ncol)
  integer ja(*)
  integer jb(*)
  integer jc(*)
  integer jcol
  integer job
  integer jpos
  integer k
  integer ka
  integer len
  integer ljob
  integer nnza
  integer nnzb
  integer nzmax
  logical values

  values = ( job /= 0 )

  ierr = 0
  iw(1:ncol) = 0

  nnza = ia(nrow+1) - 1
  nnzb = ib(ncol+1) - 1
  len = nnzb

  if ( nzmax < nnzb .or. nzmax < nnza ) then
    ierr = -1
    return
  end if
!
!  Transpose matrix B into C.
!
  ljob = 0
  if ( values ) then
    ljob = 1
  end if

  ipos = 1
  call csrcsc ( ncol, ljob, ipos, b, jb, ib, c, jc, ic )

  if ( job == -1 ) then
    c(1:len) = -c(1:len)
  end if
!
!  The main loop.
!
  do ii = 1, nrow

     do k = ic(ii), ic(ii+1)-1
        iw(jc(k)) = k
     end do

     do ka = ia(ii), ia(ii+1)-1

        jcol = ja(ka)
        jpos = iw(jcol)
!
!  If fill-in, append in coordinate format to matrix.
!
        if ( jpos == 0 ) then

           len = len + 1

           if ( nzmax < len ) then
             ierr = ii
             return
           end if

           jc(len) = jcol

           ic(len) = ii
           if ( values ) then
             c(len) = a(ka)
           end if
        else
!
!  else do addition.
!
           if ( values ) then
             c(jpos) = c(jpos) + a(ka)
           end if

        end if

     end do

     do k = ic(ii), ic(ii+1)-1
       iw(jc(k)) = 0
     end do

  end do
!
!  Convert first part of matrix (without fill-ins) into COO format.
!
  ljob = 2
  if ( values ) then
    ljob = 3
  end if

  call csrcoo ( nrow, ljob, nnzb, c, jc, ic, nnzb ,c, ic, jc, ierr )

  if ( ierr /= 0 ) then
    ierr = -ierr
  end if
!
!  Convert the whole thing back to CSR format.
!
  ljob = 0
  if ( values ) then
    ljob = 1
  end if

  call coocsr_inplace ( nrow, len, ljob, c, jc, ic, iw )

  return
end
subroutine assmb1 ( u, a, ja, ia, fu, f, node_num, element_num, element_node, &
  node_code, npe )

!*****************************************************************************80
!
!! ASSMB1 assembles a finite element matrix in the CSR format.
!
!  Discussion:
!
!    The routine receives all the unassembled local finite element
!    matrices and right hand sides, and assembles them into a global
!    matrix and right hand side.
!
!    The JA and IA arrays are constructed based on the information
!    in the element connectivity array ELEMENT_NODE.
!
!  Modified:
!
!    01 July 2005
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, real U(ELEMENT_NUM,NPE,NPE), the unassembled local matrices.
!
!    Output, real A(*), integer JA(*), IA(NODE_NUM+1), the assembled global
!    matrix in CSR (Compressed Sparse Row) format.
!
!    Input, real FU(ELEMENT_NUM,NPE), the unassembled right hand sides.
!
!    Output, real F(NODE_NUM), the assembled global right hand side.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(NPE,ELEMENT_NUM), the connectivity matrix.
!    ELEMENT_NODE(I,J) is the global index of the I-th local node in element J.
!
!    Input, integer NODE_CODE(NODE_NUM), boundary information for each node 
!    with the following meaning:
!    * 0, node I is internal;
!    * 1, node I is a boundary but not a corner point;
!    * 2, node I is a corner point.
!
!    Input, integer NPE, the number of nodes per element.
!
!  Local parameters:
!
!    Workspace, integer IWK(NODE_NUM).
!
!    Workspace, integer JWK(NODE_NUM+1).
!
!    Integer LOCAL, LOCAL1, LOCAL2, local node numbers.
!
!    Integer NODE, NODE1, NODE2, global node numbers.
!
  implicit none

  integer element_num
  integer node_num
  integer npe

  real ( kind = 8 ) a(*)
  integer element
  integer element_node(npe,element_num)
  real ( kind = 8 ) f(node_num)
  real ( kind = 8 ) fu(element_num,npe)
  integer ia(node_num+1)
  integer iwk(node_num)
  integer j
  integer ja(*)
  integer jwk(node_num+1)
  integer k
  integer k1
  integer k2
  integer local
  integer local1
  integer local2
  integer node
  integer node_code(node_num)
  integer node1
  integer node2
  integer row_last
  integer row_start
  real ( kind = 8 ) u(element_num,npe,npe)
!
!  Initialize.
!
  f(1:node_num) = 0.0D+00
!
!  Initialize the pointer arrays.
!
  ia(1:node_num+1) = 1
  jwk(1:node_num+1) = 0
!
!  Count the number of elements (or boundary conditions) where a given
!  node occurs.  Put this into the IA vector.  Then replace this count
!  by an incremental count of all the entries preceding it.
!
  do element = 1, element_num
    do local = 1, npe
      node = element_node(local,element)
      ia(node) = ia(node) + 1
    end do
  end do

  do node = 1, node_num
    if ( 1 <= node_code(node) ) then
      ia(node) = ia(node) + 1
    end if
  end do

  k1 = ia(1)
  ia(1) = 1
  do node = 2, node_num + 1
    k2 = ia(node)
    ia(node) = ia(node-1) + k1
    iwk(node-1) = ia(node-1) - 1
    k1 = k2
  end do
!
!  The assembly loop.
!
  do element = 1, element_num
!
!  The local row number is LOCAL1.
!  The global row number is NODE1.
!
    do local1 = 1, npe

      node1 = element_node(local1,element)
      f(node1) = f(node1) + fu(element,local1)
!
!  Unpack the row into JWK1.
!
      row_start = ia(node1)
      row_last = iwk(node1)

      do k = row_start, row_last
        jwk(ja(k)) = k
      end do
!
!  The local column is LOCAL2.
!  The global column number is JJ.
!
      do local2 = 1, npe

        node2 = element_node(local2,element)
        k = jwk(node2)

        if ( k == 0 ) then
          row_last = row_last + 1
          jwk(node2) = row_last
          ja(row_last) = node2
          a(row_last) = u(element,local1,local2)
        else
          a(k) = a(k) + u(element,local1,local2)
        end if

      end do
!
!  Refresh JWK.
!
      jwk(ja(row_start:row_last)) = 0
      iwk(node1) = row_last

    end do

  end do

  return
end
subroutine assmbo ( nx, nelx, node, ijk, nodcode, x, y, a, ja, ia, f, iwk, &
  jwk, ierr, xyk )

!*****************************************************************************80
!
!! ASSMBO assembles a finite element matrix.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NX, the number of nodes.
!
!    Input, integer NELX, the number of elements.
!
!    Input, integer NODE, the number of nodes per element.
!
!    Input, integer IJK(NODE,NELX), lists the nodes that make up 
!    each element.
!
!    Output, real A(*), integer JA(*), IA(NX+1), the assembled matrix in CSR
!    Compressed Sparse Row format.
!
! f      = right hand side (global load vector)
!
! nodcode= boundary information list for each node with the
!         following meaning:
!      nodcode(i) = 0 -->  node i is internal
!      nodcode(i) = 1 -->  node i is a boundary but not a corner point
!      nodcode(i) = 2 -->  node i is a corner point (corner points
!
! mp      = group material number for each element.
!
! x,y   = real arrays containing the x and y coordinates
!        resp. of the nodes.
!
! iwk,jwk = two integer work arrays.
! ierr      = error message integer .
!        ierr = 0 --> normal return
!        ierr = 1 --> negative area encountered (due to bad
!                 numbering of nodes of an element- see
!               message printed in unit iout ).
! iout      = output unit (not used here).
!
! xyk      = routine defining the material properties at each
!         element. Form:
!       call xyk(nel,xyke,x,y,ijk,node) with on return
!         xyke =  material constant matrices.
!         for each element nel, xyke(1,nel),xyke(2,nel)
!         and xyke(3,nel) represent the constants
!         K11, K22, and K12 at that element.
!
  implicit none

  integer node
  integer nx

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) det
  real ( kind = 8 ) f(*)
  real ( kind = 8 ) fe(3)
  integer i
  integer ia(nx+1)
  integer ierr
  integer ii
  integer ijk(node,*)
  integer ilast
  integer irowst
  integer iwk(*)
  integer j
  integer ja(*)
  integer jj
  integer jwk(*)
  integer k
  integer ka
  integer kb
  integer knod
  integer ksav
  integer ksavn
  integer nel
  integer nelx
  integer nodcode(*)
  real ( kind = 8 ) ske(3,3)
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xe(3)
  external xyk
  real ( kind = 8 ) xyke(2,2)
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) ye(3)
!
!  Maximum number of nonzeros allowed  = 200
!
!  Initialize.
!
  f(1:nx) = 0.0D+00
!
!  Initialize the pointer arrays.
!
  ia(1:nx+1) = 1
  jwk(1:nx+1) = 0

  do k = 1, nelx
    do j = 1, node
      knod = ijk(j,k)
      ia(knod) = ia(knod) + 1
    end do
  end do

  do k = 1, nx
    if ( 1 <= nodcode(k) ) then
      ia(k) = ia(k) + 1
    end if
  end do

  ksav = ia(1)
  ia(1) = 1
  do j = 2, nx+1
    ksavn = ia(j)
    ia(j) = ia(j-1) +  ksav
    iwk(j-1) = ia(j-1)-1
    ksav = ksavn
  end do
!
!  The main loop.
!
  do nel = 1, nelx
!
!  Get coordinates of nodal points.
!
    do i = 1, node
      j = ijk(i,nel)
      xe(i) = x(j)
      ye(i) = y(j)
    end do
!
!  Compute determinant.
!
    det = xe(2) * ( ye(3) - ye(1) ) &
        + xe(3) * ( ye(1) - ye(2) ) &
        + xe(1) * ( ye(2) - ye(3) )
!
!  Set material properties.
!
    call xyk ( nel, xyke, x, y, ijk, node )
!
!  Construct element stiffness matrix.
!
    ierr = 0

    call estif3 ( nel, ske, fe, det, xe, ye, xyke, ierr )

    if ( ierr /= 0 ) then
      return
    end if
!
!  Assemble: add element stiffness matrix to global matrix.
!
    do ka = 1, node

      ii = ijk(ka,nel)
      f(ii) = f(ii) + fe(ka)
!
!  Unpack row into JWK1.
!
      irowst = ia(ii)
      ilast  = iwk(ii)
      do k = irowst, ilast
       jwk(ja(k)) = k
      end do

      do kb = 1, node
!
!  Column number = JJ.
!
        jj = ijk(kb,nel)
        k = jwk(jj)

        if ( k == 0 ) then
          ilast = ilast + 1
          jwk(jj) = ilast
          ja(ilast) = jj
          a(ilast) = ske(ka,kb)
        else
          a(k) = a(k) + ske(ka,kb)
        end if

      end do
!
!  Refresh JWK.
!
      do k = irowst, ilast
        jwk(ja(k)) = 0
      end do

      iwk(ii) = ilast

    end do

  end do

  return
end
subroutine atmux ( n, x, y, a, ja, ia )

!*****************************************************************************80
!
!! ATMUX computes A' * x for a CSR matrix A.
!
!  Discussion:
!
!    This routine multiplies the transpose of a matrix by a vector when the
!    original matrix is stored in compressed sparse row storage. Can also be
!    viewed as the product of a matrix by a vector when the original
!    matrix is stored in the compressed sparse column format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, real X(*), an array whose length is equal to the
!    column dimension of A.
!
!    Output, real Y(N), the product A' * X.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  integer i
  integer ia(*)
  integer ja(*)
  integer k
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(n)

  y(1:n) = 0.0D+00

  do i = 1, n
    do k = ia(i), ia(i+1)-1
      y(ja(k)) = y(ja(k)) + x(i) * a(k)
    end do
  end do

  return
end
subroutine blkchk ( nrow, ja, ia, nblk, imsg )

!*****************************************************************************80
!
!! BLKCHK checks whether the input matrix is a block matrix.
!
!  Discussion:
!
!    This routine checks whether the input matrix is a block matrix with block 
!    size of NBLK.  A block matrix is one which is comprised of small square 
!    dense blocks.  If there are zero elements within the square blocks and the 
!    data structure takes them into account then blkchk may fail to find the
!    correct block size.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, JA(*), IA(NROW+1), the matrix information (but no values) in CSR
!    Compressed Sparse Row format.
!
!    Input, integer NBLK, the block value to be checked.
!
!    Output, integer IMSG, a message  with the following meaning:
!     0 : the output value of NBLK is a correct block size. 
!    -1 : NBLK does not divide NROW;
!    -2 : a starting element in a row is at wrong position
!         (j /= mult*nblk +1 );
!    -3 : NBLK does not divide a row length;
!    -4 : an element is isolated outside a block or two rows in same 
!         group have different lengths
!
  implicit none

  integer nrow

  integer i
  integer i1
  integer ia(nrow+1)
  integer ii
  integer imsg
  integer irow
  integer j
  integer j2
  integer ja(*)
  integer jstart
  integer k
  integer len
  integer lena
  integer nblk
  integer nr
!
!  First part of code will find candidate block sizes.
!  This is not guaranteed to work, so a check is done at the end.
!  The criterion used here is a simple one:
!  scan rows and determine groups of rows that have the same length
!  and such that the first column number and the last column number
!  are identical.
!
  imsg = 0

  if ( nblk <= 1 ) then
    return
  end if

  nr = nrow / nblk

  if ( nr * nblk /= nrow ) then
    imsg = -1
    return
  end if
!
!  The main loop.
!
  irow = 1

  do ii = 1, nr
!
!  I1 = starting position for group of NBLK rows in original matrix.
!
    i1 = ia(irow)
    j2 = i1
!
!  LENA = length of each row in that group in the original matrix.
!
    lena = ia(irow+1) - i1
!
!  LEN = length of each block-row in that group in the output matrix.
!
    len = lena / nblk
    if ( len * nblk /= lena ) then
      imsg = -3
      return
    end if
!
!  For each row.
!
    do i = 1, nblk

      irow = irow + 1

      if ( ia(irow) - ia(irow-1) /= lena ) then
        imsg = -4
        return
      end if
!
!  For each block.
!
      do k = 0, len-1

        jstart = ja(i1+nblk*k) - 1

        if ( ( jstart / nblk ) * nblk /= jstart ) then
          imsg = -2
          return
        end if
!
!  For each column.
!
        do j = 1, nblk
          if ( jstart + j /= ja(j2) ) then
            imsg = -4
          end if
          j2 = j2 + 1
        end do

      end do

    end do

  end do

  return
end
subroutine blkfnd ( nrow, ja, ia, nblk )

!*****************************************************************************80
!
!! BLKFND determines the block structure of a matrix.
!
!  Discussion:
!
!    If the matrix has a block structure, this routine finds the block 
!    size.  A block matrix is one which is comprised of small square 
!    dense blocks.  If there are zero elements within the square blocks 
!    and the original data structure takes these zeros into account 
!    then this routine may fail to find the correct block size.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, real JA(*), IA(NROW+1), the matrix information (but not the
!    values) in CSR Compressed Sparse Row format.
!
!    Output, integer NBLK, the block value that was found.
!
  implicit none

  integer nrow

  integer i
  integer i1
  integer i2
  integer ia(nrow+1)
  integer iblk
  integer imsg
  integer irow
  integer ja(*)
  integer jf
  integer jfirst
  integer jl
  integer jlast
  integer jrow
  integer len
  integer len0
  integer minlen
  integer nblk
!
!  The first part of code will find candidate block sizes.
!  The criterion used here is a simple one: scan rows and  determine groups
!  of rows that have the same length and such that the first column
!  number and the last column number are identical.
!
  minlen = ia(2) - ia(1)
  irow = 1

  do i = 2, nrow
    len = ia(i+1) - ia(i)
    if ( len < minlen ) then
      minlen = len
      irow = i
    end if
  end do
!
!  Candidates are all dividers of MINLEN.
!
  nblk = 1
  if ( minlen <= 1 ) then
    return
  end if

  do iblk = minlen, 1, -1

    if ( mod ( minlen, iblk ) /= 0 ) then
      cycle
    end if

     len = ia(2) - ia(1)
     len0 = len
     jfirst = ja(1)
     jlast = ja(ia(2)-1)

     do jrow = irow+1, irow+nblk-1

        i1 = ia(jrow)
        i2 = ia(jrow+1) - 1
        len = i2 + 1 - i1
        jf = ja(i1)
        jl = ja(i2)

        if ( len /= len0 .or. jf /= jfirst .or. jl /= jlast ) then
          go to 99
        end if

     end do
!
!  Check for this candidate.
!
     call blkchk ( nrow, ja, ia, iblk, imsg )
!
!  Block size found.
!
     if ( imsg == 0 ) then
        nblk = iblk
        return
     end if

 99   continue

  end do

  return
end
subroutine bndcsr ( n, abd, nabd, lowd, ml, mu, a, ja, ia, len, ierr )

!*****************************************************************************80
!
!! BNDCSR converts Banded Linpack format to Compressed Sparse Row format.
!
!  Discussion:
!
!    The matrix values found to be equal to zero
!    (actual test: if (abd(...) == 0.0) are removed.
!
!    The resulting may not be identical to a CSR matrix
!    originally transformed to a BND format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, integer NABD, the first dimension of ABD.
!
! abd   = real array containing the values of the matrix stored in
!         banded form. The j-th column of abd contains the elements
!         of the j-th column of  the original matrix,comprised in the
!         band ( i in (j-ml,j+mu) ) with the lowest diagonal located
!         in row lowd (see below).
!
! lowd  = integer. this should be set to the row number in abd where
!         the lowest diagonal (leftmost) of A is located.
!         lowd should be s.t.  ( 1  <=  lowd  <= nabd).
!         The routines dgbco, ... of linpack use lowd=2*ml+mu+1.
!
! ml      = integer. equal to the bandwidth of the strict lower part of A
! mu      = integer. equal to the bandwidth of the strict upper part of A
!         thus the total bandwidth of A is ml+mu+1.
!         if ml+mu+1 is found to be larger than nabd then an error
!         message is set. see ierr.
!
! len   = integer. length of arrays a and ja. bndcsr will stop if the
!         length of the arrays a and ja is insufficient to store the
!         matrix. see ierr.
!
! on return:
!
!    Output, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!
! lowd  = if on entry lowd was zero then lowd is reset to the default
!         value ml+mu+l.
!
! ierr  = integer. used for error message output.
!         ierr == 0 :means normal return
!         ierr == -1 : means invalid value for lowd.
!        ierr > 0 : means that there was not enough storage in a and ja
!         for storing the ourput matrix. The process ran out of space
!         (as indicated by len) while trying to fill row number ierr.
!         This should give an idea of much more storage might be required.
!         Moreover, the first irow-1 rows are correctly filled.
!
  implicit none

  integer n
  integer nabd

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) abd(nabd,*)
  integer i
  integer ia(n+1)
  integer ierr
  integer irow
  integer j
  integer ja(*)
  integer ko
  integer len
  integer lowd
  integer ml
  integer mu
  real ( kind = 8 ) t

  ierr = 0

  if ( nabd < lowd .or. lowd <= 0 ) then
    ierr = -1
    return
  end if

  ko = 1
  ia(1) = 1

  do irow = 1, n

    i = lowd

    do j = irow-ml, irow+mu

      if ( j <= 0 ) then
        go to 19
      end if

      if ( n < j ) then
        go to 21
      end if

      t = abd(i,j)

      if ( t /= 0.0D+00 ) then

        if ( len < ko ) then
          ierr = irow
          return
        end if

        a(ko) = t
        ja(ko) = j
        ko = ko + 1

      end if

19    continue

      i = i - 1

    end do

 21 continue

    ia(irow+1) = ko

  end do

  return
end
subroutine bound ( nx, nelx, ijk, nodcode, node, n_int, iperm, &
  x, y, wk, iwk )

!*****************************************************************************80
!
!! BOUND counts the number of boundary points.
!
!  Discussion:
!
!    It also reorders the points in such a way that the boundary nodes
!    are last.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NX, the number of nodes.
!
!    Input, integer NELX, the number of elements.
!
!    Input/output, integer IJK(NODE,NELX), lists the nodes that make up 
!    each element.  On output, IJK has been updated.
!
! nodcode, node: see other routines
!
!    Output, integer N_INT,  the number of points on the boundary.
!
! iperm = permutation array from old ordering to new ordering,
!
! iwk   = reverse permutation array or return.
! wk      = real work array
! On return
! x, y, nodecode, are permuted
! ijk  is updated according to new oerdering.
! n_int = number of interior points.
!
  implicit none

  integer node

  integer ijk(node,*)
  integer iperm(*)
  integer iwk(*)
  integer j
  integer k
  integer knod
  integer n_int
  integer nbound
  integer nel
  integer nelx
  integer nodcode(*)
  integer nx
  real ( kind = 8 ) wk(*)
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)
!
!  Maximum number of nonzeros allowed  = 200
!
!  Put all boundary points at the end, backwards.
!
  n_int = 1
  nbound = nx

  do j = 1, nx
    if ( nodcode(j) == 0 ) then
      iperm(n_int) = j
      n_int = n_int + 1
    else
      iperm(nbound) = j
      nbound = nbound - 1
    end if
  end do

  n_int = n_int - 1
!
!  Permute X's.
!
  wk(1:nx) = x(1:nx)
  x(1:nx) = wk(iperm(1:nx))
!
!  Permute the Y's.
!
  wk(1:nx) = y(1:nx)
  y(1:nx) = wk(iperm(1:nx))
!
!  Permute the boundary information.
!
  iwk(1:nx) = nodcode(1:nx)

  do k = 1, nx
    nodcode(k) = iwk(iperm(k))
  end do
!
!  Get the reverse permutation.
!
  do k = 1, nx
    iwk(iperm(k)) = k
  end do
!
!  Update the element connectivity matrix.
!
  do nel = 1, nelx
    do j = 1, node
      knod = ijk(j,nel)
      ijk(j,nel) = iwk(knod)
    end do
  end do

  return
end
subroutine bsort2 ( w, ind, n, ncut )

!*****************************************************************************80
!
!! BSORT2 returns the NCUT largest elements of an array, using bubble sort.
!
!  Discussion:
!
!    This routine carries out a simple bubble sort for getting the NCUT largest
!    elements in modulus, in array W.  IND is sorted accordingly.
!    (Ought to be replaced by a more efficient sort especially
!    if NCUT is not that small).
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
  implicit none

  integer n

  integer i
  integer ind(*)
  integer iswp
  integer j
  integer ncut
  logical test
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) wswp

  i = 1

  do

    test = .false.

    do j = n-1, i, -1

      if ( abs ( w(j) ) < abs ( w(j+1) ) ) then
!
!  Swap.
!
        wswp = w(j)
        w(j) = w(j+1)
        w(j+1) = wswp
!
!  Reorder the original ind array accordingly.
!
        iswp = ind(j)
        ind(j) = ind(j+1)
        ind(j+1) = iswp
!
!  Set indicator that sequence is still unsorted.
!
        test = .true.

      end if

    end do

    i = i + 1

    if ( .not. test .or. ncut < i ) then
      exit
    end if

  end do

  return
end
subroutine bsrcsr ( n, nblk, na, a, ja, ia, ao, jao, iao )

!*****************************************************************************80
!
!! BSRCSR converts Block Sparse Row to Compressed Sparse Row (CSR) format.
!
!  Discussion:
!
!    This routine converts a matrix stored in block-reduced
!    a, ja, ia format to the general sparse row a, ja, ia format.
!    A matrix that has a block structure is a matrix whose entries
!    are blocks of the same size nblk (e.g. 3 x 3).  Then it is often
!    preferred to work with the reduced graph of the matrix, i.e.,
!    Instead of storing one element at a time one can store the whole
!    block.  In this storage scheme a row of the array a will
!    hold the nblk**2 entries of a block.
!
!    This code is not "in place".
!
! general picture: (nblk = 2)
!     --- A ---                                --- JA --  -- IA --
! A=  x x x x   1st block in block row 1           x         x
!     x x x x  2-nd block in block row 1           x
!     . . . .                                      .
!     x x x x  last block in block row 1           x
!     -------                                     ---
!     x x x x   1st block in block row 2           x          x
!     x x x x  2-nd block in block row 2           x
!     . . . .                                      x
!     x x x x   last block in block row 2          x
!     -------                                     ---
!     .......                                     ...         .
!     -------                                     ---
!     x x x x   1st block in block row n/nblk      x          x
!     x x x x  2-nd block in block row n/nblk      x
!     . . . .                                      x
!     x x x x  last block in block row n/nblk      x
!     -------                                     ---
!                                               end + 1       x
!
!
! example  with nblk = 2:
!
!
!             1   2   0   0   3   4
!             5   6   0   0   7   8
!             0   0   9  10  11  12
!             0   0  13  14  15  16
!             17 18   0   0   0   0
!             22 23   0   0   0   0
! THEN:
!
!  ---- A ----                                     -- JA --   -- IA --
!-
!  1   5   2  6  Block row 1 (2 block matrices)      | 1  <--- | 1
!  3   7   4  8                                      | 5       |
!  ------------                                      | --      |
!  9  13  10 14  block row 2 (2 block matrices)      | 3  <--- | 3
! 11  15  12 16                                      | 5       |
!  ------------                                      | --      |
! 17  22  18 23  Block row 3 (1 block matrix)        | 1  <--- | 5
!  ------------                                      | --      |
!                                                   end+1 <--- | 6
!
! JA  =  1  5 | 3  5 | 1       column numbers of (1,1) entries of blocks
! IA  =  1      3      5  6    pointers to beginnings of BLOCK-rows
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
! nblk  = integer equal to the dimension of each block.
!         nblk must divide n.
!
! na      = first dimension of array a as declared in calling program
!
! a      = real array containing the values of the matrix. For details
!         on the format see below. Each row of a contains the nblk x nblk
!         block matrix unpacked column-wise (this allows the user to
!         declare the array a as a(na,nblk,nblk) on entry if desired).
!         the block rows are stored in sequence just as for the compressed
!         sparse row format.
!
! ja      = integer array of length n/nblk. ja(k) contains the column index
!         of the leading element, i.e., the element (1,1) of the block
!         that is held in the row a(k,*) of the value array.
!
! ia    = integer array of length n/nblk+1. ia(i) points to the beginning
!        of block row number i in the arrays a and ja.
!
!    Output, real AO(*), JAO(*), IAO(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer n
  integer na

  real ( kind = 8 ) a(na,*)
  real ( kind = 8 ) ao(*)
  integer i
  integer i1
  integer i2
  integer ia(*)
  integer iao(n+1)
  integer ii
  integer ij
  integer irow
  integer j
  integer ja(*)
  integer jao(*)
  integer jj
  integer jst
  integer k
  integer krow
  integer nblk
  integer nr
!
!  Get the IA, JA data structure for output matrix
!
  nr = n / nblk

  iao(1:n+1) = 0

  irow = 0
  krow = 1

  do ii = 1, nr
!
!  NR is the dimension of the reduced matrix.
!
     i1 = ia(ii)
     i2 = ia(ii+1) - 1
!
!  Create NBLK rows for each K.
!
     do i = 1, nblk
        do k = i1, i2
           jst = ja(k) - 1
           do j = 1, nblk
              ij = ( j - 1 ) * nblk + i
              ao(krow) = a(k,ij)
              jao(krow) = jst + j
              krow = krow + 1
           end do
         end do
      iao(irow+i) = krow
     end do

     irow = irow + nblk

  end do

  do jj = 1, n
     j = n - jj + 1
     iao(j+1) = iao(j)
  end do

  iao(1) = 1

  return
end
subroutine bsten ( nx, ny, nz, kx, ky, kz, nfree, stencil, h )

!*****************************************************************************80
!
!! BSTEN calculates block stencil values.
!
!  Discussion:
!
!    This routine calculates the correct block-stencil values for
!    a centered difference discretization of the elliptic operator
!    (block version of stencil)
!
! L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) +
!       d delx ( u ) + e dely (u) + f delz( u ) + g u
!
!   For 2-D problems the discretization formula that is used is:
!
! h**2 * Lu == a(i+1/2,j) * {u(i+1,j) - u(i,j)} +
!             a(i-1/2,j) * {u(i-1,j) - u(i,j)} +
!              b(i,j+1/2) * {u(i,j+1) - u(i,j)} +
!              b(i,j-1/2) * {u(i,j-1) - u(i,j)} +
!              (h/2) * d(i,j) * {u(i+1,j) - u(i-1,j)} +
!              (h/2) * e(i,j) * {u(i,j+1) - u(i,j-1)} +
!              (h/2) * e(i,j) * {u(i,j+1) - u(i,j-1)} +
!              (h**2) * g(i,j) * u(i,j)
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
  implicit none

  real ( kind = 8 ) cntr(225)
  real ( kind = 8 ) coeff(225)
  real ( kind = 8 ) h
  real ( kind = 8 ) h2
  real ( kind = 8 ) hhalf
  integer k
  integer kx
  integer ky
  integer kz
  integer nfree
  integer nfree2
  integer nx
  integer ny
  integer nz
  real ( kind = 8 ) stencil(7,*)
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  if ( 15 < nfree ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BSTEN  - FATAL ERROR'
    write ( *, '(a)' ) '  Input value of NFREE is greater than 15.'
    stop
  end if

  nfree2 = nfree * nfree

  cntr(1:nfree2) = 0.0D+00
  stencil(1:7,1:nfree2) = 0.0D+00

  hhalf = h * 0.5D+00
  h2 = h * h
  x = h * real ( kx, kind = 8 )
  y = h * real ( ky, kind = 8 )
  z = h * real ( kz, kind = 8 )
!
!  Differentiation with respect to X:
!
  call afunbl ( nfree, x + hhalf, y, z, coeff )

  do k = 1, nfree2
    stencil(3,k) = stencil(3,k) + coeff(k)
    cntr(k) = cntr(k) + coeff(k)
  end do

  call afunbl ( nfree, x - hhalf, y, z, coeff )

  do k = 1, nfree2
    stencil(2,k) = stencil(2,k) + coeff(k)
    cntr(k) = cntr(k) + coeff(k)
  end do

  call dfunbl ( nfree, x, y, z, coeff )

  do k = 1, nfree2
    stencil(3,k) = stencil(3,k) + coeff(k) * hhalf
    stencil(2,k) = stencil(2,k) - coeff(k) * hhalf
  end do

  if ( ny <= 1 ) then
    go to 99
  end if
!
!  Differentiation with respect to Y:
!
  call bfunbl ( nfree, x, y + hhalf, z, coeff )

  do k = 1, nfree2
     stencil(5,k) = stencil(5,k) + coeff(k)
     cntr(k) = cntr(k) + coeff(k)
  end do

  call bfunbl ( nfree, x, y - hhalf, z, coeff )

  do k = 1, nfree2
     stencil(4,k) = stencil(4,k) + coeff(k)
     cntr(k) = cntr(k) + coeff(k)
  end do

  call efunbl ( nfree, x, y, z, coeff )

  do k = 1, nfree2
     stencil(5,k) = stencil(5,k) + coeff(k) * hhalf
     stencil(4,k) = stencil(4,k) - coeff(k) * hhalf
  end do
!
!  Differentiation with respect to Z:
!
  if ( 1 < nz ) then

    call cfunbl ( nfree, x, y, z + hhalf, coeff )

    do k = 1, nfree2
      stencil(7,k) = stencil(7,k) + coeff(k)
      cntr(k) = cntr(k) + coeff(k)
    end do

    call cfunbl ( nfree, x, y, z - hhalf, coeff )

    do k = 1, nfree2
      stencil(6,k) = stencil(6,k) + coeff(k)
      cntr(k) = cntr(k) + coeff(k)
    end do

    call ffunbl ( nfree, x, y, z, coeff )

    do k = 1, nfree2
      stencil(7,k) = stencil(7,k) + coeff(k) * hhalf
      stencil(6,k) = stencil(6,k) - coeff(k) * hhalf
    end do

  end if
!
!  Discretization of product by G:
!
 99   continue

  call gfunbl ( nfree, x, y, z, coeff )

  do k = 1, nfree2
     stencil(1,k) = h2*coeff(k) - cntr(k)
  end do

  return
end
subroutine checkref ( nx, nelx, ijk, node, nodcode, nbound, nxnew, nelxnew )

!*****************************************************************************80
!
!! CHECKREF returns the expected number of new nodes and elements.
!
!  Discussion:
!
!    These numbers indicate the number of new nodes and elements that may
!    be expected if routine REFALL is applied once to the current grid.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NX, the number of nodes.
!
!    Input, integer NELX, the number of elements.
!
!    Input, integer IJK(NODE,NELX), lists the nodes that make up 
!    each element.
!
! nbound  = number of boundary points on entry - enter zero if
!           unknown
!
! nodcode= boundary information list for each node with the
!         following meaning:
!      nodcode(i) = 0 -->  node i is internal
!      nodcode(i) = 1 -->  node i is a boundary but not a corner point
!      nodcode(i) = 2 -->  node i is a corner point.
!
! nxnew  = new number of nodes if refall were to be applied
! nelxnew = same for nelx.
!
  implicit none

  integer node
  integer nx

  integer ijk(node,*)
  integer j
  integer nbound
  integer nelx
  integer nelxnew
  integer nodcode(nx)
  integer nxnew

  nelxnew = nelx * 4
!
!  Count the boundary nodes.
!
  if ( nbound == 0 ) then

    do j = 1, nx
      if ( 1 <= nodcode(j) ) then
        nbound = nbound + 1
      end if
    end do

  end if
!
!  Number of edges = ( 3 * ( number of elements ) + number of bound nodes ) / 2
!
  nxnew = nx + ( 3 * nelx + nbound ) / 2
  nbound = 2 * nbound

  return
end
subroutine chkelmt ( nx, x, y, nelx, ijk, node )

!*****************************************************************************80
!
!! CHKELMT checks the labeling within each element and reorders if necessary.
!
!  Discussion:
!
!    If the nodes are not correctly ordered, this routine reorders them.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NX, the number of nodes.
!
!    Input, real X(*), Y(*), the coordinates of the nodes.
!
!    Input, integer NELX, the number of elements.
!
!    Input, integer IJK(NODE,NELX), lists the nodes that make up 
!    each element.
!
!    Input, integer NODE, the number of nodes per element.
!
  implicit none

  integer node

  real ( kind = 8 ) det
  integer ijk(node,*)
  integer j
  integer nel
  integer nelx
  integer nx
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  do nel = 1, nelx

    det = x(ijk(2,nel)) * ( y(ijk(3,nel)) - y(ijk(1,nel)) ) + &
          x(ijk(3,nel)) * ( y(ijk(1,nel)) - y(ijk(2,nel)) ) + &
          x(ijk(1,nel)) * ( y(ijk(2,nel)) - y(ijk(3,nel)) )
!
!  If the determinant is negative, switch the last two nodes of the element.
!
    if ( det < 0.0D+00 ) then
      j = ijk(2,nel)
      ijk(2,nel) = ijk(3,nel)
      ijk(3,nel) = j
    end if

  end do

  return
end
subroutine cnrms ( nrow, nrm, a, ja, ia, diag )

!*****************************************************************************80
!
!! CNRMS gets the norms of each column of A.
!
!  Discussion:
!
!    There is a choice of three norms.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NRM, choosed the norm:
!    1, means 1-norm,
!    2, means the 2-nrm,
!    0, means max norm
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!   Output, real ( kind = 8 ) DIAG(NROW), the row norms.
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) diag(nrow)
  integer ia(nrow+1)
  integer ii
  integer j
  integer ja(*)
  integer k
  integer k1
  integer k2
  integer nrm

  diag(1:nrow) = 0.0D+00

  do ii = 1, nrow

     k1 = ia(ii)
     k2 = ia(ii+1) - 1

     do k = k1, k2

        j = ja(k)
!
!  Update the norm of each column.
!
        if ( nrm == 0 ) then
           diag(j) = max ( diag(j), abs ( a(k) ) )
        else if ( nrm == 1 ) then
           diag(j) = diag(j) + abs ( a(k) )
        else
           diag(j) = diag(j) + a(k)**2
        end if

    end do

  end do

  if ( nrm /= 2 ) then
    return
  end if

  do k = 1, nrow
     diag(k) = sqrt ( diag(k) )
  end do

  return
end
subroutine coocsr_inplace ( n, nnz, job, a, ja, ia, iwk )

!*****************************************************************************80
!
!! COOCSR_INPLACE converts COO to CSR in place.
!
!  Discussion:
!
!    This routine converts a matrix stored in coordinate format into
!    the CSR format.  The conversion is done in place in that the arrays
!    a,ja,ia of the result are overwritten onto the original arrays.
!
!    The entries of the output matrix are not sorted (the column
!    indices in each are not in increasing order) use COOCSR
!    if you want them sorted.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, integer NNZ, the number of nonzero elements in A.
!
!    Input, integer JOB.  When JOB = 1, the real values in A are
!    filled.  Otherwise A is not touched and the structure of the
!    array only (i.e. JA, IA)  is obtained.
!
!    Input/output, real A(NNZ).  On input, the matrix numeric values,
!    stored in the COO format.  On output, the numeric values, stored
!    in CSR format.
!
! ja      = integer array of length nnz containing the column positions
!         of the corresponding elements in a.
!
! ia      = integer array of length nnz containing the row positions
!         of the corresponding elements in a.
!
! iwk      = integer work array of length n.
!
! on return:
!
!    Output, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer n
  integer nnz

  real ( kind = 8 ) a(*)
  integer i
  integer ia(nnz)
  integer inext
  integer init
  integer ipos
  integer iwk(n)
  integer j
  integer ja(nnz)
  integer jnext
  integer job
  integer k
  real ( kind = 8 ) t
  real ( kind = 8 ) tnext
  logical values

  values = (job == 1)
!
!  Find pointer array for resulting matrix.
!
  iwk(1:n+1) = 0

  do k = 1, nnz
    i = ia(k)
    iwk(i+1) = iwk(i+1) + 1
  end do

  iwk(1) = 1
  do i = 2, n
    iwk(i) = iwk(i-1) + iwk(i)
  end do
!
!  Loop for a cycle in chasing process.
!
  init = 1
  k = 0

 5    continue

  if ( values ) then
    t = a(init)
  end if

  i = ia(init)
  j = ja(init)
  ia(init) = -1

 6 continue
   k = k + 1
!
!  Current row number is I.  Determine where to go.
!
  ipos = iwk(i)
!
!  Save the chased element.
!
  if ( values ) then
    tnext = a(ipos)
  end if

  inext = ia(ipos)
  jnext = ja(ipos)
!
!  Then occupy its location.
!
  if ( values ) then
    a(ipos) = t
  end if

  ja(ipos) = j
!
!  Update pointer information for next element to come in row I.
!
  iwk(i) = ipos + 1
!
!  Determine the next element to be chased.
!
  if ( ia(ipos) < 0 ) then
    go to 65
  end if

  t = tnext
  i = inext
  j = jnext
  ia(ipos) = -1

  if ( k < nnz ) then
    go to 6
  end if

  go to 70

 65 continue

  init = init + 1

  if ( nnz < init ) then
    go to 70
  end if

  if ( ia(init) < 0 ) then
    go to 65
  end if
!
!  Restart chasing.
!
  go to 5

 70   continue

  ia(1) = 1
  ia(2:n+1) = iwk(1:n)

  return
end
subroutine coocsr ( nrow, nnz, a, ir, jc, ao, jao, iao )

!*****************************************************************************80
!
!! COOCSR converts COO to CSR.
!
!  Discussion:
!
!    This routine converts a matrix that is stored in COO coordinate format
!    a, ir, jc into a CSR row general sparse ao, jao, iao format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NNZ, the number of nonzero elements in the matrix.
!
! a,
! ir,
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
!         the elements, ir(k) = its row number and jc(k) = its column
!        number. The order of the elements is arbitrary.
!
! on return:
!
! ir       is destroyed
!
!    Output, real AO(*), JAO(*), IAO(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer i
  integer iad
  integer iao(nrow+1)
  integer ir(*)
  integer j
  integer jao(*)
  integer jc(*)
  integer k
  integer k0
  integer nnz
  real ( kind = 8 ) x

  iao(1:nrow+1) = 0
!
!  Determine the row lengths.
!
  do k = 1, nnz
    iao(ir(k)) = iao(ir(k)) + 1
  end do
!
!  The starting position of each row.
!
  k = 1
  do j = 1, nrow+1
     k0 = iao(j)
     iao(j) = k
     k = k + k0
  end do
!
!  Go through the structure once more.  Fill in output matrix.
!
  do k = 1, nnz
     i = ir(k)
     j = jc(k)
     x = a(k)
     iad = iao(i)
     ao(iad) = x
     jao(iad) = j
     iao(i) = iad + 1
  end do
!
!  Shift back IAO.
!
  do j = nrow, 1, -1
    iao(j+1) = iao(j)
  end do
  iao(1) = 1

  return
end
subroutine cooell ( n, nnz, a, ja, ia, ac, jac, nac, ner, ncmax, ierr )

!*****************************************************************************80
!
!! COOELL converts coordinate format to Ellpack/Itpack format.
!
!  Discussion:
!
!    This routine takes a sparse matrix in coordinate format and
!    converts it into the Ellpack/Itpack storage.
!
!  Example:
!
!       (   11   0   13    0     0     0  )
!       |   21  22    0   24     0     0  |
!       |    0  32   33    0    35     0  |
!   A = |    0   0   43   44     0    46  |
!       |   51   0    0   54    55     0  |
!       (   61  62    0    0    65    66  )
!
!   Coordinate storage scheme:
!
!    A  = (11,22,33,44,55,66,13,21,24,32,35,43,46,51,54,61,62,65)
!    IA = (1, 2, 3, 4, 5, 6, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6 )
!    JA = ( 1, 2, 3, 4, 5, 6, 3, 1, 4, 2, 5, 3, 6, 1, 4, 1, 2, 5)
!
!   Ellpack/Itpack storage scheme:
!
!       (   11  13    0    0   )          (   1   3   *    *  )
!       |   22  21   24    0   |          |   2   1   4    *  |
!  AC = |   33  32   35    0   |    JAC = |   3   2   5    *  |
!       |   44  43   46    0   |          |   4   3   6    *  |
!       |   55  51   54    0   |          |   5   1   4    *  |
!       (   66  61   62   65   )          (   6   1   2    5  )
!
!   Note: * means that you can store values from 1 to 6 (1 to n, where
!         n is the order of the matrix) in that position in the array.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Ernest E. Rothman, Cornell Theory Center
!
!  Reference:
!
!    David Kincaid, T C Oppe, J R Respess, D M Young,
!    ITPACKV 2C User's Guide, 
!    Technical Report CNA-191. 
!    Center for Numerical Analysis,
!    University of Texas at Austin, 1984.
!
!    Engineering and Scientific Subroutine Library; 
!    Guide and Reference; 
!    Release 3 (SC23-0184-3), pages 79-86.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer NNZ, the number of nonzero elements in the sparse matrix.
!
!    Input, integer NCA, the first dimension of output arrays CA and JAC.
!
!  A(NNZ)  - Real array.
!            Stored entries of the sparse matrix A.
!            NNZ is the number of nonzeros.
!
!  IA(NNZ) - Integer array.
!            Pointers to specify rows for the stored nonzero entries
!            in A.
!
!  JA(NNZ) - Integer array.
!            Pointers to specify columns for the stored nonzero
!            entries in A.
!
!  NER     - Integer. Must be set greater than or equal to the maximum
!            number of nonzeros in any row of the sparse matrix.
!
!  OUTPUT PARAMETERS
!
!  AC(NAC,*)  - Real array.
!               Stored entries of the sparse matrix A in compressed
!               storage mode.
!
!  JAC(NAC,*) - Integer array.
!               Contains the column numbers of the sparse matrix
!               elements stored in the corresponding positions in
!               array AC.
!
!  NCMAX   -  Integer. Equals the maximum number of nonzeros in any
!             row of the sparse matrix.
!
!  IERR    - Error parameter is returned as zero on successful
!             execution of the subroutin<e.
!             Error diagnostics are given by means of positive values
!             of this parameter as follows:
!
!             IERR = -1   -  NER is too small and should be set equal
!                            to NCMAX. The array AC may not be large
!                            enough to accomodate all the non-zeros of
!                            of the sparse matrix.
!             IERR =  1   -  The array AC has a zero column. (Warning)
!             IERR =  2   -  The array AC has a zero row.    (Warning)
!
  implicit none

  integer nac
  integer ner
  integer nnz

  real ( kind = 8 ) a(nnz)
  real ( kind = 8 ) ac(nac,ner)
  integer ia(nnz)
  integer icount
  integer ierr
  integer ii
  integer in
  integer inn
  integer is
  integer ja(nnz)
  integer jac(nac,ner)
  integer k
  integer n
  integer ncmax
!
!  Initialize the error parameter to zero.
!
  ierr = 0
!
!  Initialize the output arrays.
!
  jac(1:n,1:ner) = n
  ac(1:n,1:ner) = 0.0D+00
!
!  Assign nonzero elements of the sparse matrix (stored in the one
!  dimensional array A to the two dimensional array AC.
!  Also, assign the correct values with information about their
!  column indices to the two dimensional array KA. And at the same
!  time count the number of nonzeros in each row so that the
!  parameter NCMAX equals the maximum number of nonzeros in any row
!  of the sparse matrix.
!
  ncmax = 1

  do is = 1, n

     k = 0

     do ii = 1, nnz
        if ( ia(ii) == is ) then
           k = k + 1
           if ( k <= ner ) then
              ac(is,k) = a(ii)
              jac(is,k) = ja(ii)
           end if
        end if
     end do

     if ( ncmax <= k ) then
       ncmax = k
     end if

  end do
!
!  Perform some simple error checks:
!
!  Check maximum number of nonzeros in each row:
!
  if ( ncmax == ner ) then
    ierr = 0
  end if

  if ( ner < ncmax ) then
     ierr = -1
     return
  end if
!
!  Check if there are any zero columns in AC:
!
  do in = 1, ncmax

     icount = 0

     do inn = 1, n
        if ( ac(inn,in) /= 0.0D+00 ) then
          icount = 1
        end if
     end do

     if ( icount == 0 ) then
        ierr = 1
        return
     end if

  end do
!
!  Check if there are any zero rows in AC:
!
  do inn = 1, n

     icount = 0

     do in = 1, ncmax
        if ( ac(inn,in) /= 0.0D+00 ) then
          icount = 1
        end if
     end do

     if ( icount == 0 ) then
        ierr = 2
        return
     end if

  end do

  return
end
subroutine copmat ( nrow, a, ja, ia, ao, jao, iao, ipos )

!*****************************************************************************80
!
!! COPMAT copies the matrix A, JA, IA, into the matrix AO, JAO, IAO.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, integer IPOS, indicates the position in the array AO, JAO
!    where the first element should be copied.  Thus IAO(1) = IPOS on return.
!
!    Output, real AO(*), JAO(*), IAO(NROW+1), the copied matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer i
  integer ia(*)
  integer iao(nrow+1)
  integer ipos
  integer ja(*)
  integer jao(*)
  integer k
  integer kst

  kst = ipos - ia(1)

  do i = 1, nrow+1
    iao(i) = ia(i) + kst
  end do

  do k = ia(1), ia(nrow+1)-1
    ao(kst+k) = a(k)
    jao(kst+k)= ja(k)
  end do

  return
end
subroutine cperm ( nrow, a, ja, ia, ao, jao, iao, perm, job )

!*****************************************************************************80
!
!! CPERM permutes the columns of a matrix.
!
!  Discussion:
!
!    This routine permutes the columns of a matrix a, ja, ia.
!    The result is written in the output matrix  ao, jao, iao.
!    cperm computes B = A P, where  P is a permutation matrix
!    that maps column j into column perm(j), i.e., on return
!    a(i,j) becomes a(i,perm(j)) in new matrix
!
!    if job=1 then ao, iao are not used.
!    This routine is in place: ja, jao can be the same.
!    If the matrix is initially sorted (by increasing column number)
!    then ao,jao,iao  may not be on return.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! perm      = integer array of length ncol (number of columns of A
!         containing the permutation array  the columns:
!         a(i,j) in the original matrix becomes a(i,perm(j))
!         in the output matrix.
!
! job      = integer indicating the work to be done:
!             job = 1      permute a, ja, ia into ao, jao, iao
!                       (including the copying of real values ao and
!                       the array iao).
!             job /= 1 :  ignore real values ao and ignore iao.
!
!
! on return:
!
! ao, jao, iao = input matrix in a, ja, ia format (array ao not needed)
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer ia(nrow+1)
  integer iao(nrow+1)
  integer ja(*)
  integer jao(*)
  integer job
  integer k
  integer nnz
  integer perm(*)

  nnz = ia(nrow+1)-1

  do k = 1, nnz
    jao(k) = perm(ja(k))
  end do
!
!  Done with the JA array.  Return if no need to touch values.
!
  if ( job /= 1 ) then
    return
  end if
!
!  Else get new pointers, and copy values too.
!
  iao(1:nrow+1) = ia(1:nrow+1)

  ao(1:nnz) = a(1:nnz)

  return
end
subroutine cscal ( nrow, job, nrm, a, ja, ia, diag, b, jb, ib )

!*****************************************************************************80
!
!! CSCAL scales the columns of A such that their norms are one.
!
!  Discussion:
!
!    The result matrix is written on B, or overwritten on A.
!
!    3 choices of norms: 1-norm, 2-norm, max-norm. in place.
!
!    The column dimension of A is not needed.
!
!    The algorithm in place (B can take the place of A).
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
! job   = integer. job indicator. Job=0 means get array b only
!         job = 1 means get b, and the integer arrays ib, jb.
!
! nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
!                  means the 2-nrm, nrm = 0 means max norm
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! on return:
!
!
! diag = diagonal matrix stored as a vector containing the matrix
!        by which the columns have been scaled, i.e., on return
!        we have B = A * Diag
!
!    Output, real B(*), integer JB(*), IB(NROW+1), the scaled matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) diag(nrow)
  integer ia(nrow+1)
  integer ib(nrow+1)
  integer ja(*)
  integer jb(*)
  integer job
  integer nrm

  call cnrms ( nrow, nrm, a, ja, ia, diag )

  diag(1:nrow) = 1.0D+00 / diag(1:nrow)

  call amudia ( nrow, job, a, ja, ia, diag, b, jb, ib )

  return
end
subroutine csort ( n, a, ja, ia, iwork, values )

!*****************************************************************************80
!
!! CSORT sorts the elements of a CSR matrix.
!
!  Discussion:
!
!    This routine sorts the elements of a CSR matrix (stored in Compressed
!    Sparse Row Format) in increasing order of their column indices within
!    each row. It uses a form of bucket sort with a cost of O(nnz) where
!    nnz = number of nonzero elements.
!
!    Requires an integer work array of size length 2*nnz.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! iwork = integer work array of length max ( n+1, 2*nnz )
!         where nnz = 2* (ia(n+1)-ia(1))  ) .
!
! values= logical indicating whether or not the real values a(*) must
!         also be permuted. if (.not. values) then the array a is not
!         touched by csort and can be a dummy array.
!
! on return:
!
! the matrix stored in the structure a, ja, ia is permuted in such a
! way that the column indices are in increasing order within each row.
! iwork(1:nnz) contains the permutation used  to rearrange the elements.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  integer i
  integer ia(n+1)
  integer ifirst
  integer irow
  integer iwork(*)
  integer j
  integer ja(*)
  integer k
  integer ko
  integer next
  integer nnz
  logical values
!
!  Count the number of elements in each column.
!
  iwork(1:n+1) = 0

  do i = 1, n
     do k = ia(i), ia(i+1)-1
        j = ja(k) + 1
        iwork(j) = iwork(j) + 1
     end do
  end do
!
!  Compute pointers from lengths.
!
  iwork(1) = 1

  do i = 1, n
     iwork(i+1) = iwork(i) + iwork(i+1)
  end do
!
!  Get the positions of the nonzero elements in order of columns.
!
  ifirst = ia(1)
  nnz = ia(n+1)-ifirst

  do i = 1, n
    do k = ia(i), ia(i+1)-1
      j = ja(k)
      next = iwork(j)
      iwork(nnz+next) = k
      iwork(j) = next + 1
    end do
  end do
!
!  Convert to coordinate format.
!
  do i = 1, n
    do k = ia(i), ia(i+1)-1
      iwork(k) = i
    end do
  end do
!
!  Loop to find permutation: for each element find the correct
!  position in (sorted) arrays A, JA.  Record this in IWORK.
!
  do k = 1, nnz
     ko = iwork(nnz+k)
     irow = iwork(ko)
     next = ia(irow)
!
!  The current element should go in next position in row. IWORK
!  records this position.
!
     iwork(ko) = next
     ia(irow) = next + 1
  end do
!
!  Perform an in-place permutation of the arrays.
!
     call ivperm ( nnz, ja(ifirst), iwork )

     if ( values ) then
       call dvperm ( nnz, a(ifirst), iwork )
     end if
!
!  Reshift the pointers of the original matrix back.
!
  do i = n, 1, -1
    ia(i+1) = ia(i)
  end do

  ia(1) = ifirst

  return
end
subroutine csrbnd ( n, a, ja, ia, job, abd, nabd, lowd, ml, mu, ierr )

!*****************************************************************************80
!
!! CSRBND converts Compressed Sparse Row to Banded Linpack format.
!
!  Discussion:
!
!    This routine converts a general sparse matrix stored in
!    compressed sparse row format into the banded format.  For the
!    banded format,the Linpack conventions are assumed (see below).
!
!    Additional details on banded format.  (this closely follows the      
!    format used in linpack. may be useful for converting a matrix into   
!    this storage format in order to use the linpack  banded solvers).    
!
!    Band storage format  for matrix abd             
!    uses ml+mu+1 rows of abd(nabd,*) to store the diagonals of           
!    a in rows of abd starting from the lowest (sub)-diagonal  which  is  
!    stored in row number lowd of abd. the minimum number of rows needed  
!    in abd is ml+mu+1, i.e., the minimum value for lowd is ml+mu+1. the  
!    j-th  column  of  abd contains the elements of the j-th column of a, 
!    from bottom to top: the element a(j+ml,j) is stored in  position     
!    abd(lowd,j), then a(j+ml-1,j) in position abd(lowd-1,j) and so on.   
!    Generally, the element a(j+k,j) of original matrix a is stored in    
!    position abd(lowd+k-ml,j), for k = ml,ml-1,..,0,-1, -mu.               
!    The first dimension nabd of abd must be >= lowd                    
!                                                                      
!    example [from linpack ]:   if the original matrix is             
!                                                                      
!              11 12 13  0  0  0                                       
!              21 22 23 24  0  0                                       
!               0 32 33 34 35  0     original banded matrix            
!               0  0 43 44 45 46                                       
!               0  0  0 54 55 56                                       
!               0  0  0  0 65 66                                       
!                                                                      
!    then  n = 6, ml = 1, mu = 2. lowd should be >= 4 (=ml+mu+1)  and   
!    if lowd = 5 for example, abd  should be:                             
!                                                                      
!    untouched --> x  x  x  x  x  x                                       
!                  *  * 13 24 35 46                                       
!                  * 12 23 34 45 56    resulting abd matrix in banded     
!                 11 22 33 44 55 66    format                             
!     row lowd--> 21 32 43 54 65  *                                       
!                                                                      
!    * = not used                                                         
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! job      = integer. if job=1 then the values of the lower bandwith ml
!         and the upper bandwidth mu are determined internally.
!         otherwise it is assumed that the values of ml and mu
!         are the correct bandwidths on input. See ml and mu below.
!
! nabd  = integer. first dimension of array abd.
!
! lowd  = integer. this should be set to the row number in abd where
!         the lowest diagonal (leftmost) of A is located.
!         lowd should be  ( 1  <=  lowd  <= nabd).
!         if it is not known in advance what lowd should be
!         enter lowd = 0 and the default value lowd = ml+mu+1
!         will be chosen. Alternative: call routine getbwd from unary
!         first to detrermione ml and mu then define lowd accordingly.
!         (Note: the banded solvers in linpack use lowd=2*ml+mu+1. )
!
! ml      = integer. equal to the bandwidth of the strict lower part of A
! mu      = integer. equal to the bandwidth of the strict upper part of A
!         thus the total bandwidth of A is ml+mu+1.
!         if ml+mu+1 is found to be larger than lowd then an error
!         flag is raised (unless lowd = 0). see ierr.
!
! note:   ml and mu are assumed to have       the correct bandwidth values
!         as defined above if job is set to zero on entry.
!
! on return:
!
! abd   = real array of dimension abd(nabd,n).
!         on return contains the values of the matrix stored in
!         banded form. The j-th column of abd contains the elements
!         of the j-th column of  the original matrix comprised in the
!         band ( i in (j-ml,j+mu) ) with the lowest diagonal at
!         the bottom row (row lowd). See details below for this format.
!
! ml      = integer. equal to the bandwidth of the strict lower part of A
! mu      = integer. equal to the bandwidth of the strict upper part of A
!         if job=1 on entry then these two values are internally computed.
!
! lowd  = integer. row number in abd where the lowest diagonal
!         (leftmost) of A is located on return. In case lowd = 0
!         on return, then it is defined to ml+mu+1 on return and the
!         lowd will contain this value on return. `
!
! ierr  = integer. used for error messages. On return:
!         ierr == 0  :means normal return
!         ierr == -1 : means invalid value for lowd. (either < 0
!         or larger than nabd).
!         ierr == -2 : means that lowd is not large enough and as
!         result the matrix cannot be stored in array abd.
!         lowd should be at least ml+mu+1, where ml and mu are as
!         provided on output.
!
  implicit none

  integer n
  integer nabd

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) abd(nabd,n)
  integer i
  integer ia(n+1)
  integer ierr
  integer ii
  integer j
  integer ja(*)
  integer job
  integer k
  integer lowd
  integer m
  integer mdiag
  integer ml
  integer mu
!
!  Determine ML and MU.
!
  ierr = 0

  if ( job == 1 ) then
    call getbwd ( n, a, ja, ia, ml, mu )
  end if

  m = ml + mu + 1

  if ( lowd == 0 ) then
    lowd = m
  end if

  if ( lowd < m ) then
    ierr = -2
  end if

  if ( nabd < lowd .or. lowd < 0 ) then
    ierr = -1
  end if

  if ( ierr < 0 ) then
    return
  end if

  do i = 1, m
    ii = lowd - i + 1
    abd(ii,1:n) = 0.0D+00
  end do

  mdiag = lowd - ml

  do i = 1, n
    do k = ia(i), ia(i+1)-1
      j = ja(k)
      abd(i-j+mdiag,j) = a(k)
    end do
  end do

  return
end
subroutine csrbsr ( n, nblk, na, a, ja, ia, ao, jao, iao )

!*****************************************************************************80
!
!! CSRBSR converts Compressed Sparse Row to Block Sparse Row.
!
!  Discussion:
!
!    This routine does the reverse of BSRCSR.  It converts
!    a matrix stored in a general compressed a, ja, ia format into a
!    a block reduced matrix a(*,*),ja(*),ia(*) format. The code
!    assumes that the original matrix is indeed a block format
!    and that the elements are ordered in such a way that their
!    column numbers are increasing. (This can be achieved
!    by transposing a, ja, ia twice, putting the resulting matrix
!    into a, ja, ia).
!
!    See routine bsrcsr for more details on data structure for blocked
!    matrices. The input matrix is a, ja, ia (in compressed format) and
!    the output matrix is the matrix ao, jao, iao in block-reduced
!    format.
!
!    This code is not in place.
!
!    See routine bsrcsr for details on data sctructure for
!    block sparse row format.
!
!    The routine assumes that  the input matrix has been
!    sorted in such a way that the column indices are always
!    in increasing order for the same row.
!    for all k "in the SAME ROW."
!
!    THERE IS NO CHECKING AS TO WHETHER the input is correct.
!    it is recommended to use the routine blchk to check
!    if the matrix is a block-matrix before calling csrbsr.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, integer NBLK, the dimension of each block.
!    NBLK must divide N.
!
! na      = first dimension of array ao as declared in calling program
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! on return:
!
! ao    = real array containing the values of the matrix. For details
!         on the format see below. Each row of a contains the nblk x nblk
!         block matrix unpacked column-wise (this allows the user to
!         declare the array a as a(na,nblk,nblk) on entry if desired).
!         the block rows are stored in sequence just as for the compressed
!         sparse row format.
! jao   = integer array of length n/nblk. ja(k) contains the column index
!         of the leading element, i.e., the element (1,1) of the block
!         that is held in the row a(k,*) of the value array.
! iao   = integer array of length n/nblk+1. ia(i) points to the beginning
!        of block row number i in the arrays a and ja.
!
  implicit none

  integer n
  integer na

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(na,*)
  integer i
  integer i1
  integer ia(n+1)
  integer iao(*)
  integer ibrow
  integer ii
  integer irow
  integer j
  integer j1
  integer j2
  integer ja(*)
  integer jao(*)
  integer k
  integer k1
  integer len
  integer lena
  integer nblk
  integer nr
!
!  NR is the dimension of the reduced matrix.
!
  nr = n / nblk
  iao(1) = 1
  ibrow = 1
  irow  = 1
!
!  The main loop.
!
  do ii = 1, nr
!
!  I1 = starting position for group of nblk rows in original matrix.
!
     i1 = ia(irow)
!
!  LENA = length of each row in that group in the original matrix.
!
     lena = ia(irow+1) - i1
!
!  LEN = length of each block-row in that group in the output matrix.
!
     len = lena / nblk
     k1 = iao(ibrow)
!
!  Copy the real values of A.
!
!  For each block.
!
     do k = 0, len-1
!
!  Store column positions of the (1,1) elements of each block.
!
        jao(k1+k) = ja(i1+nblk*k)
!
!  For each column.
!
        do j = 1, nblk

           j1 = ( j - 1 ) * nblk
           j2 = i1 + k * nblk + j - 1
!
!  For each row.
!
           do i = 1, nblk
              ao(k1+k,j1+i) = a(j2+(i-1)*lena)
           end do

         end do

    end do
!
!  Done with a whole block row. 
!  Update IAO, IBROW and IROW.
!
     iao(ibrow+1) = iao(ibrow) + len
     ibrow = ibrow + 1
     irow  = irow + nblk

  end do

  return
end
subroutine csrcoo ( nrow, job, nzmax, a, ja, ia, nnz, ao, ir, jc, ierr )

!*****************************************************************************80
!
!! CSRCOO converts Compressed Sparse Row to Coordinate format.
!
!  Discussion:
!
!   This routine converts a matrix that is stored in row general sparse 
!   A, JA, IA format into coordinate format AO, IR, JC. 
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
! job   = integer serving as a job indicator.
!         if job = 1 fill in only the array ir, ignore jc, and ao.
!         if job = 2 fill in ir, and jc but not ao
!         if job = 3 fill in everything.
!         The reason why these options are provided is that on return
!         ao and jc are the same as a, ja. So when job = 3, a and ja are
!         simply copied into ao, jc.  When job=2, only jc and ir are
!         returned. With job=1 only the array ir is returned. Moreover,
!         the algorithm is in place:
!           call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr)
!         will write the output matrix in coordinate format on a, ja,ia.
!         (Important: note the order in the output arrays a, ja, ia. )
!         i.e., ao can be the same as a, ir can be the same as ia
!         and jc can be the same as ja.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! nzmax = length of space available in ao, ir, jc.
!         the code will stop immediatly if the number of
!         nonzero elements found in input matrix exceeds nzmax.
!
! on return:
!-
! ao, ir, jc = matrix in coordinate format.
!
! nnz        = number of nonzero elements in matrix.
!
! ierr       = integer error indicator.
!         ierr == 0 means normal retur
!         ierr == 1 means that the the code stopped
!         because there was no space in ao, ir, jc
!         (according to the value of  nzmax).
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer i
  integer ia(nrow+1)
  integer ierr
  integer ir(*)
  integer ja(*)
  integer jc(*)
  integer job
  integer k
  integer k1
  integer k2
  integer nnz
  integer nzmax

  ierr = 0
  nnz = ia(nrow+1)-1

  if ( nzmax < nnz ) then
    ierr = 1
    return
  end if

  if ( 3 <= job ) then
    ao(1:nnz) = a(1:nnz)
  end if

  if ( 2 <= job ) then
    jc(1:nnz) = ja(1:nnz)
  end if
!
!  Copy backward.
!
  do i = nrow, 1, -1
    k1 = ia(i+1) - 1
    k2 = ia(i)
    do k = k1, k2, -1
      ir(k) = i
    end do
  end do

  return
end
subroutine csrcsc ( n, job, ipos, a, ja, ia, ao, jao, iao )
 
!*****************************************************************************80
!
!! CSRCSC converts Compressed Sparse Row to Compressed Sparse Column.
!
!  Discussion:
!
!    This is essentially a transposition operation.  
!
!    It is NOT an in-place algorithm.
!
!    This routine transposes a matrix stored in a, ja, ia format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer JOB, indicates whether or not to fill the values of the
!    matrix AO or only the pattern (IA, and JA).  Enter 1 for yes.
!
! ipos  = starting position in ao, jao of the transposed matrix.
!         the iao array takes this into account (thus iao(1) is set to ipos.)
!         Note: this may be useful if one needs to append the data structure
!         of the transpose to that of A. In this case use
!                call csrcsc (n,1,n+2,a,ja,ia,a,ja,ia(n+2))
!        for any other normal usage, enter ipos=1.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real AO(*), JAO(*), IAO(N+1), the matrix in CSC
!    Compressed Sparse Column format.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer i
  integer ia(n+1)
  integer iao(n+1)
  integer ipos
  integer j
  integer ja(*)
  integer jao(*)
  integer job
  integer k
  integer next
!
!  Compute lengths of rows of A'.
!
  iao(1:n+1) = 0

  do i = 1, n
    do k = ia(i), ia(i+1)-1
      j = ja(k) + 1
      iao(j) = iao(j) + 1
    end do
  end do
!
!  Compute pointers from lengths.
!
  iao(1) = ipos
  do i = 1, n
    iao(i+1) = iao(i) + iao(i+1)
  end do
!
!  Do the actual copying.
!
  do i = 1, n
    do k = ia(i), ia(i+1)-1
      j = ja(k)
      next = iao(j)
      if ( job == 1 ) then
        ao(next) = a(k)
      end if
      jao(next) = i
      iao(j) = next + 1
    end do
  end do
!
!  Reshift IAO and leave.
!
  do i = n, 1, -1
    iao(i+1) = iao(i)
  end do
  iao(1) = ipos

  return
end
subroutine csrdia ( n, idiag, job, a, ja, ia, ndiag, diag, ioff, ao, &
  jao, iao, ind )

!*****************************************************************************80
!
!! CSRDIA converts Compressed Sparse Row to diagonal format.
!
!  Discussion:
!
!    This routine extracts IDIAG diagonals from the input matrix A,
!    JA, IA, and puts the rest of the matrix in the output matrix AO,
!    JAO, IAO.  The diagonals to be extracted depend on the value of JOB.
!
!    In  the first case, the diagonals to be
!    extracted are simply identified by their offsets provided in ioff
!    by the caller.  In the second case, the code internally determines
!    the idiag most significant diagonals, i.e., those diagonals of the
!    matrix which have the largest number of nonzero elements, and
!    extracts them.
!
!    The algorithm is in place: ao, jao, iao can be overwritten on
!    a, ja, ia if desired.
!
!    When the code is required to select the diagonals (job >= 10)
!    the selection of the diagonals is done from left to right
!    as a result if several diagonals have the same weight (number
!    of nonzero elemnts) the leftmost one is selected first.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, integer IDIAG.  On intput, the number of diagonals 
!    to be extracted.  On output, IDIAG may be modified to the
!    actual number of diagonals found.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! job      = integer. serves as a job indicator.  Job is better thought
!         of as a two-digit number job=xy. If the first (x) digit
!         is one on entry then the diagonals to be extracted are
!         internally determined. In this case csrdia exctracts the
!         idiag most important diagonals, i.e. those having the largest
!         number on nonzero elements. If the first digit is zero
!         then csrdia assumes that ioff(*) contains the offsets
!         of the diagonals to be extracted. there is no verification
!         that ioff(*) contains valid entries.
!         The second (y) digit of job determines whether or not
!         the remainder of the matrix is to be written on ao,jao,iao.
!         If it is zero  then ao, jao, iao is not filled, i.e.,
!         the diagonals are found  and put in array diag and the rest is
!         is discarded. if it is one, ao, jao, iao contains matrix
!         of the remaining elements.
!         Thus:
!         job= 0 means do not select diagonals internally (pick those
!                defined by ioff) and do not fill ao,jao,iao
!         job= 1 means do not select diagonals internally
!                      and fill ao,jao,iao
!         job=10 means  select diagonals internally
!                      and do not fill ao,jao,iao
!         job=11 means select diagonals internally
!                      and fill ao,jao,iao
!
!  Input, integer NDIAG, the first dimension of array DIAG.
!
! on return:
!
! diag  = real array of size (ndiag x idiag) containing the diagonals
!         of A on return
!
! ioff  = integer array of length idiag, containing the offsets of the
!           diagonals to be extracted.
! ao, jao
!  iao  = remainder of the matrix in a, ja, ia format.
!
! work arrays:
!
! ind   = integer array of length 2*n-1 used as integer work space.
!         needed only when job>=10 i.e., in case the diagonals are to
!         be selected internally.
!
  implicit none

  integer idiag
  integer ndiag

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  real ( kind = 8 ) diag(ndiag,idiag)
  integer i
  integer ia(*)
  integer iao(*)
  integer idum
  integer ii
  integer ind(*)
  integer ioff(*)
  integer j
  integer ja(*)
  integer jao(*)
  integer jmax
  integer job
  integer job1
  integer job2
  integer k
  integer ko
  integer l
  integer n
  integer n2

  job1 = job / 10
  job2 = job - job1 * 10

  if ( job1 /= 0 ) then

    n2 = n + n - 1
    call infdia ( n, ja, ia, ind, idum )
!
!  Determine the diagonals to accept.
!
    ii = 0

    do

      ii = ii + 1
      jmax = 0

      do k = 1, n2

        j = ind(k)

        if ( jmax < j ) then
          i = k
          jmax = j
        end if

      end do

      if ( jmax <= 0 ) then
        ii = ii - 1
        exit
      end if

      ioff(ii) = i - n
      ind(i) = - jmax

      if ( idiag <= ii ) then
        exit
      end if

    end do

    idiag = ii

  end if
!
!  Initialize DIAG to zero.
!
  diag(1:n,1:idiag) = 0.0D+00

  ko = 1
!
!  Extract diagonals and accumulate remaining matrix.
!
  do i = 1, n

     do k = ia(i), ia(i+1)-1

        j = ja(k)

        do l = 1, idiag
           if ( j - i == ioff(l) ) then
             diag(i,l) = a(k)
             go to 51
           end if
        end do
!
!  Append element not in any diagonal to AO, JAO, IAO.
!
        if ( job2 /= 0 ) then
          ao(ko) = a(k)
          jao(ko) = j
          ko = ko + 1
        end if

51      continue

     end do

     if ( job2 /= 0 ) then
       ind(i+1) = ko
     end if

  end do
!
!  Finish with IAO.
!
  if ( job2 /= 0 ) then
    iao(1) = 1
    iao(2:n+1) = ind(2:n+1)
  end if

  return
end
subroutine csrdns ( nrow, ncol, a, ja, ia, dns, ndns, ierr )

!*****************************************************************************80
!
!! CSRDNS converts Compressed Sparse Row to Dense format.
!
!  Discussion:
!
!    This routine converts a row-stored sparse matrix into a densely stored one.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NCOL, the column dimension of the matrix.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real DNS(NDNS,NDNS), the dense array containing a
!    copy of the matrix.
!
!    Input, integer NDNS, the dimension of the DNS array.
!
!    Output, integer IERR, error indicator.
!    0, means normal return
!    i, means that the code has stopped when processing
!       row number i, because it found a column number > ncol.
!
  implicit none

  integer ncol
  integer ndns

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) dns(ndns,ncol)
  integer i
  integer ia(*)
  integer ierr
  integer j
  integer ja(*)
  integer k
  integer nrow
  
  ierr = 0
  dns(1:nrow,1:ncol) = 0.0D+00

  do i = 1, nrow
    do k = ia(i), ia(i+1)-1
      j = ja(k)
      if ( ncol < j ) then
        ierr = i
        return
      end if
      dns(i,j) = a(k)
    end do
  end do

  return
end
subroutine csrell ( nrow, a, ja, ia, maxcol, coef, jcoef, ncoef, &
  ndiag, ierr )

!*****************************************************************************80
!
!! CSRELL converts Compressed Sparse Row to Ellpack/Itpack format
!
!  Discussion:
!
!    This routine converts a matrix stored in the general A, JA, IA
!    format into the COEF, JCOEF Ellpack/Itpack format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix A.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix, stored in
!    compressed sparse row format.
!
!    Input, integer NCOEF, the first dimension of arrays COEF, and JCOEF.
!
!    Input, integer MAXCOL, the number of columns available in COEF.
!
!    Output, real COEF(NCOEF,MAXCOL), the values of the matrix A in
!    Ellpack/Itpack format.
!
!    Output, integer JCOEF(NCOEF,MAXCOL), the column indices of each entry
!    in COEF.
!
!    Output, integer NDIAG, the number of active 'diagonals' found.
!
!    Output, integer IERR, an error flag. 
!    0 = correct return. 
!    nonzero means that NDIAG, the number of diagonals found, exceeds
!    the limit of MAXCOL.
!
  implicit none

  integer maxcol
  integer ncoef
  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) coef(ncoef,maxcol)
  integer i
  integer ia(nrow+1)
  integer ierr
  integer ja(*)
  integer jcoef(ncoef,maxcol)
  integer k
  integer k1
  integer k2
  integer ndiag
!
!  Determine the length of each row of the lower part of A.
!
  ierr = 0
  ndiag = 0
  do i = 1, nrow
    k = ia(i+1) - ia(i)
    ndiag = max ( ndiag, k )
  end do
!
!  Check that sufficient columns are available.
!
  if ( maxcol < ndiag ) then
    ierr = 1
    return
  end if
!
!  Initialize COEF and JCOEF.
!
  coef(1:nrow,1:ndiag) = 0.0D+00
  jcoef(1:nrow,1:ndiag) = 1
!
!  Copy elements by row.
!
  do i = 1, nrow

    k1 = ia(i)
    k2 = ia(i+1)-1

    do k = k1, k2
      coef(i,k-k1+1) = a(k)
      jcoef(i,k-k1+1) = ja(k)
    end do

  end do

  return
end
subroutine csrjad ( nrow, a, ja, ia, idiag, iperm, ao, jao, iao )

!*****************************************************************************80
!
!! CSRJAD converts Compressed Sparse Row to Jagged Diagonal storage.
!
!  Discussion:
!
!    This routine converts a matrix stored in the compressed sparse
!    row format to the jagged diagonal format.  The data structure
!    for the JAD (Jagged Diagonal storage) is as follows.  The rows of
!    the matrix are implicitly permuted so that their lengths are in
!    decreasing order.  The real entries AO(*) and their column indices
!    JAO(*) are stored in succession.  The number of such diagonals is IDIAG.
!    The lengths of each of these diagonals is stored in IAO(*).
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Reference:
!
!    E. Anderson, Youcef Saad,
!    Solving sparse triangular systems on parallel computers,
!    International Journal of High Speed Computing, 
!    Volume 1, pages 73-96, 1989.
!
!    Youcef Saad, 
!    Krylov Subspace Methods on Supercomputers,
!    SIAM Journal on Statistical and Scientific Computing, 
!    Volume 10, pages 1200-1232, 1989.
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! on return:
!
!    Output, integer IDIAG, the number of jagged diagonals in the data
!    structure A, JA, IA.
!
!    Output, integer IPERM(NROW), the permutation of the rows that leads 
!    to a decreasing order of the number of nonzero elements.
!
! ao    = real array containing the values of the matrix A in
!         jagged diagonal storage. The j-diagonals are stored
!         in ao in sequence.
!
!    Output, integer JAO(*), the column indices of the entries in ao.
!
! iao   = integer array containing pointers to the beginning
!         of each j-diagonal in ao, jao. iao is also used as
!         a work array and it should be of length n at least.
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer i
  integer ia(nrow+1)
  integer iao(nrow)
  integer idiag
  integer ilo
  integer iperm(nrow)
  integer j
  integer ja(*)
  integer jao(*)
  integer jj
  integer k
  integer k0
  integer k1
  integer len
!
!  Define initial IPERM and get lengths of each row.
!  JAO is used a work vector to store tehse lengths.
!
  idiag = 0
  ilo = nrow
  do j = 1, nrow
    iperm(j) = j
    len = ia(j+1) - ia(j)
    ilo = min ( ilo, len )
    idiag = max ( idiag, len )
    jao(j) = len
  end do
!
!  Call the sorter to get permutation.  Use IAO as a work array.
!
  call dcsort ( jao, nrow, iao, iperm, ilo, idiag )
!
!  Define the output data structure.  First lengths of the J-diagonals.
!
  iao(1:nrow) = 0

  do k = 1, nrow
    len = jao(iperm(k))
    do i = 1, len
      iao(i) = iao(i) + 1
    end do
  end do
!
!  Get the output matrix itself.
!
  k1 = 1
  k0 = k1
  do jj = 1, idiag
    len = iao(jj)
    do k = 1, len
      i = ia(iperm(k)) + jj -1
      ao(k1) = a(i)
      jao(k1) = ja(i)
      k1 = k1 + 1
    end do
    iao(jj) = k0
    k0 = k1
  end do

  iao(idiag+1) = k1

  return
end
subroutine csrlnk ( n, a, ja, ia, link )

!*****************************************************************************80
!
!! CSRLNK converts Compressed Sparse Row to Linked storage format.
!
!  Discussion:
!
!    This routine translates a matrix stored in compressed sparse
!    row into one with a linked list storage format.  Only the link
!    array needs to be obtained since the arrays A, JA, and IA may
!    be unchanged and have carry the same meaning for the output matrix.
!
!    In other words a, ja, ia, link   ia the output linked list data
!    structure with a, ja, ia being the same.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! on return:
!
! a     = nonzero elements.
!
! ja    = column positions.
!
! ia    = points to the first row of matrix in structure.
!
! link      = integer array of size containing the linked list information.
!         link(k) points to the next element of the row after element
!         ao(k), jcol(k). if link(k) = 0, then there is no next element,
!         i.e., ao(k), jcol(k) is the last element of the current row.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  integer i
  integer ia(n+1)
  integer ja(*)
  integer k
  integer link(*)
!
!  Loop through all rows.
!
  do i = 1, n
    do k = ia(i), ia(i+1)-2
      link(k) = k + 1
    end do
    link(ia(i+1)-1) = 0
  end do

  return
end
subroutine csrmsr ( n, a, ja, ia, ao, jao, wk, iwk )

!*****************************************************************************80
!
!! CSRMSR converts Compressed Sparse Row to Modified Sparse Row.
!
!  Discussion:
!
!    This routine converts a general sparse matrix a, ja, ia into
!    a compressed matrix using a separated diagonal (referred to as
!    the bell-labs format as it is used by bell labs semi conductor
!    group. We refer to it here as the modified sparse row format.
!
!    This has been coded in such a way that one can overwrite
!    the output matrix onto the input matrix if desired by a call of
!    the form
!
!     call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
!
!    In case ao, jao, are different from a, ja, then one can
!    use ao, jao as the work arrays in the calling sequence:
!
!     call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
!
!    Algorithm is in place.  i.e. both:
!
!          call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
!          (in which  ao, jao, are different from a, ja)
!           and
!          call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
!          (in which  wk, jwk, are different from a, ja)
!        are OK.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! on return :
!
! ao, jao  = sparse matrix in modified sparse row storage format:
!         +  ao(1:n) contains the diagonal of the matrix.
!         +  ao(n+2:nnz) contains the nondiagonal elements of the
!             matrix, stored rowwise.
!         +  jao(n+2:nnz) : their column indices
!         +  jao(1:n+1) contains the pointer array for the nondiagonal
!             elements in ao(n+1:nnz) and jao(n+2:nnz).
!             i.e., for i <= n+1 jao(i) points to beginning of row i
!            in arrays ao, jao.
!             here nnz = number of nonzero elements+1
!
!    Work array, real WK(N).
!
!    Work array, integer IWK(N+1).
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer i
  integer ia(n+1)
  integer icount
  integer ii
  integer iptr
  integer iwk(n+1)
  integer j
  integer ja(*)
  integer jao(*)
  integer k
  real ( kind = 8 ) wk(n)

  icount = 0
!
!  Store away diagonal elements and count nonzero diagonal elements.
!
  do i = 1, n
    wk(i) = 0.0D+00
    iwk(i+1) = ia(i+1) - ia(i)
    do k = ia(i), ia(i+1)-1
      if ( ja(k) == i ) then
        wk(i) = a(k)
        icount = icount + 1
        iwk(i+1) = iwk(i+1) - 1
      end if
    end do
  end do
!
!  Compute total length.
!
  iptr = n + ia(n+1) - icount
!
!  Copy backwards, to avoid collisions.
!
  do ii = n, 1, -1
    do k = ia(ii+1)-1, ia(ii), -1
      j = ja(k)
      if ( j /= ii ) then
        ao(iptr) = a(k)
        jao(iptr) = j
        iptr = iptr - 1
      end if
    end do
  end do
!
!  Compute the pointer values and copy WK.
!
  jao(1) = n + 2
  do i = 1, n
    ao(i) = wk(i)
    jao(i+1) = jao(i) + iwk(i+1)
  end do

  return
end
subroutine csrncf ( nrow, a, ja, ia, maxnz, nonz, coef, jcoef, ierr )

!*****************************************************************************80
!
!! CSRNCF converts CSR to NSPCG NCF format.
!
!  Discussion:
!
!    This routine converts a matrix stored in the general A, JA, IA
!    compressed sparse row format into the Nonsymmetric Coordinate Format
!    used as storage format 5 by NSPCG.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix A.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix, stored in
!    compressed sparse row format.
!
!    Input, integer MAXNZ, the maximum number of nonzeros allowed for
!    in the storage of COEF and JCOEF.
!
!    Output, integer NONZ, the actual number of nonzeros encountered.
!
!    Output, real COEF(MAXNZ), the values of the matrix A in NCF format.
!
!    Output, integer JCOEF(MAXNZ,2), the row and column indices of each 
!    entry in COEF.
!
!    Output, integer IERR, an error flag. 
!    0 = correct return. 
!    nonzero means that MAXNZ < NONZ.
!
  implicit none

  integer maxnz
  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) coef(maxnz)
  integer i
  integer ia(nrow+1)
  integer ierr
  integer ja(*)
  integer jcoef(maxnz,2)
  integer k
  integer k1
  integer k2
  integer nonz

  ierr = 0
!
!  Initialize COEF and JCOEF.
!
  coef(1:maxnz) = 0.0D+00
  jcoef(1:maxnz,1:2) = 0
!
!  The first N entries are reserved for the diagonals.
!
  do i = 1, nrow
    jcoef(i,1:2) = i
  end do

  nonz = nrow

  if ( maxnz < nonz ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSRNCF - Fatal error!'
    write ( *, '(a)' ) '  MAXNZ < NONZ.'
    ierr = 1
    return
  end if

  do i = 1, nrow

    k1 = ia(i)
    k2 = ia(i+1) - 1

    do k = k1, k2

      if ( ja(k) == i ) then

        coef(i) = coef(i) + a(k)

      else if ( 0.0D+00 /= a(k) ) then

        nonz = nonz + 1

        if ( maxnz < nonz ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CSRNCF - Fatal error!'
          write ( *, '(a)' ) '  MAXNZ < NONZ.'
          ierr = 1
          return
        end if

        coef(nonz) = a(k)
        jcoef(nonz,1) = i
        jcoef(nonz,2) = ja(k)

      end if

    end do

  end do

  return
end
subroutine csrssk ( n, imod, a, ja, ia, asky, isky, nzmax, ierr )

!*****************************************************************************80
!
!! CSRSSK converts Compressed Sparse Row to Symmetric Skyline Format.
!
!  Discussion:
!
!    This routine translates a compressed sparse row or a symmetric
!    sparse row format into a symmetric skyline format.
!    the input matrix can be in either compressed sparse row or the
!    symmetric sparse row format.  The output matrix is in a symmetric
!    skyline format: a real array containing the (active portions) of the
!    rows in  sequence and a pointer to the beginning of each row.
!
!    This module is NOT in place.
!
!    Even when imod = 2, length of  isky is  n+1, not n.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
! imod  = integer indicating the variant of skyline format wanted:
!         imod = 0 means the pointer isky points to the `zeroth'
!         element of the row, i.e., to the position of the diagonal
!         element of previous row (for i = 1, isky(1)= 0)
!         imod = 1 means that itpr points to the beginning of the row.
!         imod = 2 means that isky points to the end of the row (diagonal
!                  element)
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, integer NZMAX, the amount of storage available in ASKY.
!
! on return:
!
! asky    = real array containing the values of the matrix stored in skyline
!         format. asky contains the sequence of active rows from
!         i = 1, to n, an active row being the row of elemnts of
!         the matrix contained between the leftmost nonzero element
!         and the diagonal element.
!
! isky      = integer array of size n+1 containing the pointer array to
!         each row. The meaning of isky depends on the input value of
!         imod (see above).
!
! ierr  =  integer.  Error message. If the length of the
!         output array asky exceeds nzmax. ierr returns the minimum value
!         needed for nzmax. otherwise ierr=0 (normal return).
!
  implicit none

  integer n
  integer nzmax

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) asky(nzmax)
  integer i
  integer ia(n+1)
  integer ierr
  integer imod
  integer isky(n+1)
  integer j
  integer ja(*)
  integer k
  integer kend
  integer ml
  integer nnz
!
!  Determine the individual bandwidths and pointers.
!
  ierr = 0
  isky(1) = 0
  do i = 1, n
    ml = 0
    do k = ia(i), ia(i+1)-1
      ml = max ( ml, i-ja(k)+1 )
    end do
    isky(i+1) = isky(i) + ml
  end do
!
!  Test if there is enough space in ASKY to do the copying.
!
  nnz = isky(n+1)

  if ( nzmax < nnz ) then
    ierr = nnz
    return
  end if
!
!  Fill ASKY with zeros.
!
  asky(1:nnz) = 0.0D+00
!
!  Copy the nonzero elements.
!
  do i = 1, n
    kend = isky(i+1)
    do k = ia(i), ia(i+1)-1
      j = ja(k)
      if ( j <= i ) then
        asky(kend+j-i) = a(k)
      end if
    end do
  end do
!
!  Modify the pointer according to IMOD if necessary.
!
  if ( imod == 1 ) then
    do k = 1, n+1
      isky(k) = isky(k) + 1
    end do
  else if ( imod == 2 ) then
    do k = 1, n
      isky(k) = isky(k+1)
    end do
  end if

  return
end
subroutine csrssr ( nrow, a, ja, ia, nzmax, ao, jao, iao, ierr )
 
!*****************************************************************************80
!
!! CSRSSR converts Compressed Sparse Row to Symmetric Sparse Row.
!
!  Discussion:
!
!    This routine extracts the lower triangular part of a matrix.
!
!    It can be used as a means for converting a symmetric matrix for
!    which all the entries are stored in sparse format into one
!    in which only the lower part is stored. The routine uses an 
!    in place algorithm, in that the output matrix ao, jao, iao can 
!    be overwritten on the input matrix  a, ja, ia if desired.  
!
!    This routine has been coded to
!    put the diagonal elements of the matrix in the last position in
!    each row (i.e. in position  ao(ia(i+1)-1   of ao and jao)
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, integer NZMAX, the length of AO and JAO.
!
! On return:
!
! ao, jao,
!     iao = lower part of input matrix (a,ja,ia) stored in compressed sparse
!          row format format.
!
! ierr   = integer error indicator.
!          ierr == 0  means normal return
!          ierr == i  means that the code has stopped when processing
!          row number i, because there is not enough space in ao, jao
!          (according to the value of nzmax)
!
  implicit none

  integer nrow
  integer nzmax

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(nzmax)
  integer i
  integer ia(nrow+1)
  integer iao(*)
  integer ierr
  integer ja(*)
  integer jao(nzmax)
  integer k
  integer kdiag
  integer ko
  integer kold
  real ( kind = 8 ) t

  ierr = 0
  ko = 0

  do i = 1, nrow

    kold = ko
    kdiag = 0

    do k = ia(i), ia(i+1) -1

      if ( ja(k) <= i ) then
        ko = ko + 1
        if ( nzmax < ko ) then
          ierr = i
          return
        end if
        ao(ko) = a(k)
        jao(ko) = ja(k)
        if ( ja(k) == i ) then 
          kdiag = ko
        end if
      end if

    end do
!
!  Exchange.
!
    if ( kdiag /= 0 .and. kdiag /= ko ) then

       t = ao(kdiag)
       ao(kdiag) = ao(ko)
       ao(ko) = t

       k = jao(kdiag)
       jao(kdiag) = jao(ko)
       jao(ko) = k

     end if

     iao(i) = kold + 1

  end do
!
!  Redefine IAO(N+1).
!
  iao(nrow+1) = ko + 1

  return
end
subroutine daxpy ( n, da, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DAXPY computes constant times a vector plus a vector.
!
!  Discussion:
!
!    Uses unrolled loops for increments equal to one.
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Dongarra, Moler, Bunch, Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979.
!
!    Lawson, Hanson, Kincaid, Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of elements in DX and DY.
!
!    Input, real ( kind = 8 ) DA, the multiplier of DX.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer INCX, the increment between successive entries of DX.
!
!    Input/output, real ( kind = 8 ) DY(*), the second vector.
!    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
!
!    Input, integer INCY, the increment between successive entries of DY.
!
  implicit none

  real ( kind = 8 ) da
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
  integer i
  integer incx
  integer incy
  integer ix
  integer iy
  integer m
  integer n

  if ( n <= 0 ) then
    return
  end if

  if ( da  == 0.0D+00 ) then
    return
  end if
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
  if ( incx /= 1 .or. incy /= 1 ) then

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      dy(iy) = dy(iy) + da * dx(ix)
      ix = ix + incx
      iy = iy + incy
    end do
!
!  Code for both increments equal to 1.
!
  else

    m = mod ( n, 4 )

    do i = 1, m
      dy(i) = dy(i) + da * dx(i)
    end do

    do i = m+1, n, 4
      dy(i  ) = dy(i  ) + da * dx(i  )
      dy(i+1) = dy(i+1) + da * dx(i+1)
      dy(i+2) = dy(i+2) + da * dx(i+2)
      dy(i+3) = dy(i+3) + da * dx(i+3)
    end do

  end if

  return
end
subroutine dcn ( ar, ia, ja, n, ne, ic, nn, ierr )

!*****************************************************************************80
!
!! DCN generates sparse square matrices of type D(N,C).
!
!  Discussion:
!
!    The routine generates sparse square matrices of the type D(N,C).
!
!    This type of matrix has the following characteristics:
!
!    * 1's on the diagonal, 
!
!    * three bands at the distance C above the diagonal and reappearing 
!      cyclicly under it, 
!
!    * a 10 x 10 triangle of elements in the upper right hand corner.
!
!    This routine generates the matrix in the storage by indices mode.
!
!    If A is a sparse matrix of type D(N,C), then
!
!      min|A(i,j)| = 1,
!      max|A(i,j)| = max ( 1000, N + 1 )
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Ernest E Rothman,
!    Cornell Theory Center
!
!  Reference:
!
!    Zahari Zlatev, Kjeld Schaumburg, Jerzy Wasniewski,
!    A testing Scheme for Subroutines Solving Large Linear Problems,
!    Computers and Chemistry, 
!    Volume 5, Number 2-3, pages 91-100, 1981.
!
!    Ole Osterby, Zahari Zlatev,
!    Direct Methods for Sparse Matrices,
!    Springer-Verlag 1983.
!
!  Parameters:
!
!    Output, real AR(NN), the numerical values of the sparse matrix.
!
!    Output, integer IA(NN), the corresponding rows of the entries of AR.
!
!    Output, integer JA(NN), the corresponding columns of the entries of AR.
!
!    Input, integer N, the order of the matrix.  N must be at least 14.
!
!    Input, integer NE, the number of nonzero elements in the matrix.
!    NE = 4*N + 55.
!
!    Input, integer IC, sets the sparsity pattern.  
!    0 < IC < N-12 is required.
!
!    Input, integer NN, the dimension of AR, IA, and JA.  NN must be at
!    least NE.
!
!    Output, integer IERR, an error flag.
!    0, no error.
!    1, N is out of range.
!    2, IC is out of range.
!    3, NN is out of range.
!
  implicit none

  integer nn

  real ( kind = 8 ) ar(nn)
  integer i
  integer ia(nn)
  integer ic
  integer icount
  integer ierr
  integer ilast
  integer it
  integer j
  integer ja(nn)
  integer n
  integer ne

  ierr = 0
!
!  Check the input parameters.
!
  if ( n <= 13 ) then
    ierr = 1
    return
  end if

  if ( ic <= 0 .or. n - 12 <= ic ) then
    ierr = 2
    return
  end if

  ne = 4 * n + 55

  if ( nn < ne ) then
    ierr = 3
    return
  end if
!
!  Begin to generate the nonzero elements as well as the row and column
!  pointers:
!
  ar(1:n) = 1.0D+00

  do i = 1, n
    ia(i) = i
    ja(i) = i
  end do

  ilast = n

  do i = 1, n-ic
    it = ilast + i
    ar(it) = 1.0D+00 + real ( i, kind = 8 )
    ia(it) = i
    ja(it) = i + ic
  end do

  ilast = ilast + n - ic

  do i = 1, n-ic-1
    it = ilast + i
    ar(it) = - real ( i, kind = 8 )
    ia(it) = i
    ja(it) = i + ic + 1
  end do

  ilast = ilast + n - ic - 1

  do i = 1, n-ic-2
    it = ilast + i
    ar(it) = 16.0D+00
    ia(it) = i
    ja(it) = i + ic + 2
  end do

  ilast = ilast + n - ic - 2
  icount = 0

  do j = 1, 10
    do i = 1, 11-j
      icount = icount + 1
      it = ilast + icount
      ar(it) = 100.0D+00 * real ( j, kind = 8 )
      ia(it) = i
      ja(it) = n - 11 + i + j
    end do
  end do

  icount = 0
  ilast = 55 + ilast

  do i = n-ic+1, n
    icount = icount + 1
    it = ilast + icount
    ar(it) = 1.0D+00 + real ( i, kind = 8 )
    ia(it) = i
    ja(it) = i - n + ic
  end do

  ilast = ilast + ic
  icount = 0

  do i = n-ic, n
    icount = icount + 1
    it = ilast + icount
    ar(it) = - real ( i, kind = 8 )
    ia(it) = i
    ja(it) = i - n + ic + 1
  end do

  ilast = ilast + ic + 1
  icount = 0

  do i = n-ic-1, n
    icount = icount + 1
    it = ilast + icount
    ar(it) = 16.0D+00
    ia(it) = i
    ja(it) = i - n + ic + 2
  end do

  return
end
subroutine dcsort ( ival, n, icnt, index, ilo, ihi )

!*****************************************************************************80
!
!! DCSORT computes a sorting permutation for a vector.
!
!  Discussion:
!
!    This routine computes a permutation which, when applied to the
!    input vector IVAL, sorts the integers in ival in descending
!    order.  The permutation is represented by the vector INDEX.  The
!    permuted IVAL can be interpreted as follows:
!
!      ival(index(i-1)) >= ival(index(i)) >= ival(index(i+1))
!
!    A specialized sort, the distribution counting sort, is used
!    which takes advantage of the knowledge that
!    1)  The values are in the (small) range [ ilo, ihi ]
!    2)  Values are likely to be repeated often
!
!    The permutation is NOT applied to the vector IVAL.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Michael Heroux, Sandra Carney
!    Mathematical Software Research Group
!    Cray Research, Inc.
!
!  Reference:
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 3: Sorting and Searching,
!    Addison-Wesley, 1973, pages 78-79.
!
!  Parameters:
!
!    Input, integer IVAL(N), the values to be sorted.
!
!    Input, integer N, the number of values to be sorted.
!
!    Workspace, integer ICNT(IHI-ILO+1).
!
!    Output, integer INDEX(N), the permutation which sorts the IVAL.
!
!    Input, integer ILO, IHI, the minimum and maximum values in IVAL
!    to be sorted.
!
  implicit none

  integer ihi
  integer ilo
  integer n

  integer i
  integer icnt(ilo:ihi)
  integer index(n)
  integer ival(n)
  integer ivalj
  integer j

  icnt(ilo:ihi) = 0

  do i = 1, n
    icnt(ival(i)) = icnt(ival(i)) + 1
  end do

  do i = ihi-1, ilo, -1
    icnt(i) = icnt(i) + icnt(i+1)
  end do

  do j = n, 1, -1
    ivalj = ival(j)
    index(icnt(ivalj)) = j
    icnt(ivalj) = icnt(ivalj) - 1
  end do

  return
end
function ddot ( n, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DDOT forms the dot product of two vectors.
!
!  Discussion:
!
!    This routine uses unrolled loops for increments equal to one.
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Dongarra, Moler, Bunch, Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) DX(*), the first vector.
!
!    Input, integer INCX, the increment between successive entries in X.
!
!    Input, real ( kind = 8 ) DY(*), the second vector.
!
!    Input, integer INCY, the increment between successive entries in Y.
!
!    Output, real DDOT, the sum of the product of the corresponding
!    entries of X and Y.
!
  implicit none

  real ( kind = 8 ) ddot
  real ( kind = 8 ) dtemp
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
  integer i
  integer incx
  integer incy
  integer ix
  integer iy
  integer m
  integer n

  ddot = 0.0D+00
  dtemp = 0.0D+00

  if ( n <= 0 ) then
    return
  end if
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
  if ( incx /= 1 .or. incy /= 1 ) then

    if ( 0 <= incx ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( 0 <= incy ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      dtemp = dtemp + dx(ix) * dy(iy)
      ix = ix + incx
      iy = iy + incy
    end do
!
!  Code for both increments equal to 1.
!
  else

    m = mod ( n, 5 )

    do i = 1, m
      dtemp = dtemp + dx(i) * dy(i)
    end do

    do i = m+1, n, 5

      dtemp = dtemp + dx(i  ) * dy(i  ) &
                    + dx(i+1) * dy(i+1) &
                    + dx(i+2) * dy(i+2) &
                    + dx(i+3) * dy(i+3) &
                    + dx(i+4) * dy(i+4)
    end do

  end if

  ddot = dtemp

  return
end
subroutine diacsr ( n, job, idiag, diag, ndiag, ioff, a, ja, ia )

!*****************************************************************************80
!
!! DIACSR converts diagonal format to compressed sparse row
!
!  Discussion:
!
!    This routine extract the IDIAG most important diagonals from the
!    input matrix a, ja, ia, that is, those diagonals of the matrix which have
!    the largest number of nonzero elements.  If requested (see job),
!    the rest of the matrix is put in a the output matrix ao, jao, iao
!
!    The arrays A and JA should be of length n*idiag.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer JOB, if 0, indicates that entries in DIAG that
!    are exactly zero are not to be included in the output matrix.
!    0, then check for each entry in DIAG
!
!    Input, integer IDIAG, the number of diagonals to be extracted.
!
!    Output, real DIAG(NDIAG,IDIAG), the diagonals of A.
!
!    Input, integer NDIAG, the first dimension of DIAG.
!
!    Input, integer IOFF(IDIAG), the offsets of the diagonals to be 
!    extracted.
!
!    Output, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer idiag
  integer n
  integer ndiag

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) diag(ndiag,idiag)
  integer i
  integer ia(n+1)
  integer ioff(*)
  integer j
  integer ja(*)
  integer jj
  integer job
  integer ko
  real ( kind = 8 ) t

  ia(1) = 1
  ko = 1

  do i = 1, n

    do jj = 1, idiag

      j = i + ioff(jj)
      if ( j < 1 .or. n < j ) then
        cycle
      end if

      t = diag(i,jj)

      if ( job == 0 .and. t == 0.0D+00 ) then
        cycle
      end if

      a(ko) = t
      ja(ko) = j
      ko = ko + 1

    end do

    ia(i+1) = ko

  end do

  return
end
subroutine diamua ( nrow, job, a, ja, ia, diag, b, jb, ib )

!*****************************************************************************80
!
!! DIAMUA performs the matrix by matrix product B = Diag * A.
!
!  Discussion:
!
!    The column dimension of A is not needed.
!
!    The algorithm is in-place; that is, B can take the place of A.
!    in this case use job=0.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer JOB, indicates the job to be done.
!    0, means get array B only;
!    1, means get B, and the integer arrays IB and JB.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real DIAG(N), a diagonal matrix stored as a vector.
!
!    Output, real B(*), integer JB(*), integer IB(NROW+1), the resulting 
!    matrix B in compressed sparse row sparse format.
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) diag(nrow)
  integer ia(nrow+1)
  integer ib(nrow+1)
  integer ii
  integer ja(*)
  integer jb(*)
  integer job
  integer k
  integer k1
  integer k2
  real ( kind = 8 ) scal

  do ii = 1, nrow
!
!  Normalize each row.
!
     k1 = ia(ii)
     k2 = ia(ii+1) - 1
     scal = diag(ii)
     b(k1:k2) = a(k1:k2) * scal

  end do

  if ( job == 0 ) then
    return
  end if

  ib(1) = ia(1)

  do ii = 1, nrow
    ib(ii) = ia(ii)
    do k = ia(ii), ia(ii+1)-1
      jb(k) = ja(k)
    end do
  end do

  return
end
subroutine diapos ( n, ja, ia, idiag )

!*****************************************************************************80
!
!! DIAPOS returns the positions of the diagonal elements of a sparse matrix.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, JA(*), IA(N+1), the matrix information, (but no values) 
!    in CSR Compressed Sparse Row format.
!
!    Output, integer IDIAG(N); the I-th entry of IDIAG points to the 
!    diagonal element A(I,I) in the arrays A and JA.  That is,
!    A(IDIAG(I)) = element A(I,I) of matrix A.  If no diagonal element 
!    is found, the entry is set to 0.
!
  implicit none

  integer n

  integer i
  integer ia(n+1)
  integer idiag(n)
  integer ja(*)
  integer k

  idiag(1:n) = 0
!
!  Sweep through the data structure.
!
  do i = 1, n
    do k = ia(i), ia(i+1) -1
      if ( ja(k) == i ) then
        idiag(i) = k
      end if
    end do
  end do

  return
end
subroutine dinfo1 ( n, iout, a, ja, ia, valued, title, key, type, ao, jao, iao )

!*****************************************************************************80
!
!! DINFO1 computes and prints matrix statistics.
!
!  Discussion:
!
!    This routine obtains a number of statistics on a sparse matrix and writes 
!    it into the output unit iout.  The matrix is assumed
!    to be stored in the compressed sparse COLUMN format sparse a, ja, ia
!
!    On return, elementary statistics on the matrix are written on output unit
!    iout, and the entries of a, ja, ia are sorted.
!
!    title, key, type are the same paramaters as those
!    used for Harwell-Bowing matrices.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the column dimension of the matrix.
!
!    Input, integer IOUT, the FORTRAN unit number where the information
!    is to be output.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.  If values are not provided,
!    then A may be just a dummy array.
!
!    Input, logical VALUED, is TRUE if values are provided.
!
!    Input, character ( len = 72 ) TITLE, a title describing the matrix
!    The first character in title is ignored (it is often a one).
!
!    Input, character ( len = 8 ) KEY, an 8-character key for the matrix
!
! type  = a 3-character string to describe the type of the matrix.
!         see harwell/Boeing documentation for more details on the
!         above three parameters.
!
! on return
!
! ao      = real array of length nnz used as work array.
! if values are not provided, then AO may be a dummy array.
!
! jao      = integer work array of length max ( 2*n+1, nnz )
!
! iao   = integer work array of length n+1
!
! Output description:
!
! The following info needs to be updated.
!
! + A header containing the Title, key, type of the matrix and, if values
!   are not provided a message to that effect.
!
!    SYMMETRIC STRUCTURE MEDIEVAL RUSSIAN TOWNS
!                       Key = RUSSIANT , Type = SSA
!    No values provided - Information of pattern only
!
!
! +  dimension n, number of nonzero elements nnz, average number of
!    nonzero elements per column, standard deviation for this average.
! +  if the matrix is upper or lower triangular a message to that effect
!    is printed. Also the number of nonzeros in the strict upper
!    (lower) parts and the main diagonal are printed.
! +  weight of longest column. This is the largest number of nonzero
!    elements in a column encountered. Similarly for weight of
!    largest/smallest row.
! +  lower dandwidth as defined by
!          ml = max ( i-j, / all  a(i,j)/= 0 )
! +  upper bandwidth as defined by
!          mu = max ( j-i, / all  a(i,j)/= 0 )
!    NOTE that ml or mu can be negative. ml < 0 would mean
!    that A is confined to the strict upper part above the diagonal
!    number -ml. Similarly for mu.
!
! +  maximun bandwidth as defined by
!    Max (  Max [ j ; a(i,j) /= 0 ] - Min [ j ; a(i,j) /= 0 ] )
!     i
! +  average bandwidth = average over all columns of the widths each column.
!
! +  If there are zero columns /or rows a message is printed
!    giving the number of such columns/rows.
!
! +  matching elements in A and transp(A) :this counts the number of
!    positions (i,j) such that if a(i,j) /= 0 then a(j,i) /= 0.
!    if this number is equal to nnz then the matrix is symmetric.
! +  Relative symmetry match : this is the ratio of the previous integer
!    over nnz. If this ratio is equal to one then the matrix has a
!    symmetric structure.
!
! +  average distance of a given element from the diagonal, standard dev.
!    the distance of a(i,j) is defined as abs ( j-i ).
!
! +  Frobenious norm of A
!    Frobenious norm of 0.5*(A + transp(A))
!    Frobenious norm of 0.5*(A - transp(A))
!    these numbers provide information on the degree of symmetry
!    of the matrix. If the norm of the nonsymmetric part is
!    zero then the matrix is symmetric.
!
! + 90% of matrix is in the band of width k, means that
!   by moving away and in a symmetric manner from the main
!   diagonal you would have to include exactly k diagonals
!   (k is always odd), in order to include 90% of the nonzero
!   elements of A.  The same thing is then for 80%.
!
! + The total number of nonvoid diagonals, i.e., among the
!   2n-1 diagonals of the matrix which have at least one nonxero
!   element.
!
! +  Most important diagonals. The code selects a number of k
!    (k <= 10) diagonals that are the most important ones, i.e.
!    that have the largest number of nonzero elements. Any diagonal
!    that has fewer than 1% of the nonzero elements of A is dropped.
!    the numbers printed are the offsets with respect to the
!    main diagonal, going from left tp right.
!    Thus 0 means the main diagonal -1 means the subdiagonal, and
!    +10 means the 10th upper diagonal.
! +  The accumulated percentages in the next line represent the
!    percentage of the nonzero elements represented by the diagonals
!    up the current one put together.
!    Thus:
!    *  The 10 most important diagonals are (offsets)    :             *
!    *     0     1     2    24    21     4    23    22    20    19     *
!    *  The accumulated percentages they represent are   :             *
!    *  40.4  68.1  77.7  80.9  84.0  86.2  87.2  88.3  89.4  90.4     *
!    *-----------------------------------------------------------------*
!    shows the offsets of the most important  diagonals and
!    40.4 represent ratio of the number of nonzero elements in the
!    diagonal zero (main diagonal) over the total number of nonzero
!    elements. the second number indicates that the diagonal 0 and the
!    diagonal 1 together hold 68.1% of the matrix, etc..
!
! +  Block structure:
!    if the matrix has a block structure then the block size is found
!    and printed. Otherwise the info1 will say that the matrix
!    does not have a block structure. Note that block structure has
!    a very specific meaning here. the matrix has a block structure
!    if it consists of square blocks that are dense. even if there
!    are zero elements in the blocks  they should be represented
!    otherwise it would be possible to determine the block size.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) aii
  real ( kind = 8 ) amx
  real ( kind = 8 ) ao(*)
  real ( kind = 8 ) av
  real ( kind = 8 ) bndav
  real ( kind = 8 ) dcount(20)
  real ( kind = 8 ) ddomc
  real ( kind = 8 ) ddomr
  real ( kind = 8 ) dianrm
  real ( kind = 8 ) dist
  real ( kind = 8 ) dsumc
  real ( kind = 8 ) dsumr
  real ( kind = 8 ) eps
  integer i
  integer ia(n+1)
  integer iacc
  integer iao(n+1)
  integer iband
  integer idiag
  integer ii
  integer imatch
  integer indiag
  integer ioff(20)
  integer iout
  integer, parameter :: ipar1 = 1
  integer ipos
  integer itot
  integer j
  integer j0
  integer j0r
  integer j1
  integer j1r
  integer j2
  integer ja(*)
  integer jao(*)
  integer jb1
  integer jb2
  integer jmax
  integer jmaxc
  integer jminc
  integer jminr
  integer job
  integer k
  integer k1
  integer k1max
  integer k2
  integer k2max
  character ( len = 8 ) key
  integer lenc
  integer lenr
  integer ml
  integer mu
  integer n2
  integer nblk
  integer nddomc
  integer nddomr
  integer ndiag
  integer nlower
  integer nnz
  integer nsky
  integer nskyl
  integer nskyu
  integer nupper
  integer nzcol
  integer nzdiag
  integer nzmaxc
  integer nzmaxr
  integer nzminc
  integer nzminr
  integer nzrow
  real ( kind = 8 ) st
  real ( kind = 8 ) std
  logical sym
  real ( kind = 8 ) ta
  real ( kind = 8 ) tan
  real ( kind = 8 ) tas
  character ( len = 72 ) title
  character( len = 61 ) tmpst
  character ( len = 3 ) type
  logical valued

  write (iout,99)
  write (iout,97) title(2:72), key, type
 97     format(2x,' * ',a71,' *'/, &
           2x,' *',20x,'Key = ',a8,' , Type = ',a3,25x,' *')

  if ( .not. valued ) then
    write ( iout, '(a)' ) '  No values provided - Information on pattern only'
  end if

  nnz = ia(n+1) - ia(1)
  sym = ( type(2:2) == 'S' )

  write ( iout, 99)
  write ( iout, 100) n, nnz
!
!  Average and standard deviation.
!
  av = real ( nnz, kind = 8 ) / real ( n, kind = 8 )
!
!  AV will be corrected later.
!
  if ( sym ) then
    av = 2.0D+00 * av - 1.0D+00
  end if

  job = 0
  if ( valued ) then
    job = 1
  end if

  ipos = 1
  call csrcsc ( n, job, ipos, a, ja, ia, ao, jao, iao )
  call csrcsc(n, job, ipos, ao, jao, iao, a, ja, ia)
!
!  Bandwidth.
!
  iband = 0
!
!  Number of nonzero elements in lower part
!
  nupper = 0
!
!  Number of nonzero elements in diagonal.
!
  ndiag  = 0
!
!  Distance of an element from diagonal.
!
  dist = 0.0D+00
!
!  Number of diagonally dominant columns
!
  nddomc = 0
!
!  Number of diagonally dominant rows
!
  nddomr = 0
!
!  Maximum length of columns.
!
  nzmaxc = 0
!
!  Minimum length of column.
!
  nzminc = n
!
!  Maximum length of rows.
!
  nzmaxr = 0
!
!  Minimum length of rows.
!
  nzminr = n
!
!  Number of zero columns.
!
  nzcol  = 0
!
!  Number of zero columns.
!
  nzrow  = 0
!
!  Lower and upper bandwidths.
!
  ml = -n
  mu = -n
!
!  Average bandwidth.
!
  bndav = 0.0D+00
!
!  Standard deviation for average nonzero elements.
!
  st = 0.0D+00
!
!  DIANRM = Frobenius norm of the diagonal (used only in symmetric case).
!
  dianrm = 0.0D+00
!
!  NSKYU = skyline storage for upper part.
!
  nskyu = 0
!
!  NSKYL = skyline storage for lower part.
!
  nskyl = 0
!
!  Computing maximum bandwith, maximum number of nonzero elements per column,
!  minimum nonzero elements per column, column and column diagonal dominance
!  occurrences, average distance of an element from diagonal, number of
!  elemnts in lower and upper parts.
!
  do i = 1, n

    j0 = ia(i)
    j1 = ia(i+1) - 1
    j0r = iao(i)
    j1r = iao(i+1) - 1
!
!  Bandwidth info:
!
    jminc = ja(j0)
    jmaxc = ja(j1)
    jminr = jao(j0r)
    if ( sym ) then
      jminc = jminr
    end if

    nskyl = nskyl + i - jminr + 1
    nskyu = nskyu + i - jminc + 1

    ml = max ( ml, i-jminc )
    mu = max ( mu, jmaxc-i )
    iband = max ( iband, jmaxc-jminc+1 )
    bndav = bndav + real ( jmaxc-jminc+1, kind = 8 )
!
!  Maximum and minimum number of nonzero elements per column.
!
    lenr = j1r + 1 - j0r
    if ( lenr <= 0 ) then
      nzrow = nzrow + 1
    end if
    nzmaxr = max ( nzmaxr, lenr )
    nzminr = min ( nzminr, lenr )
!
!  INDIAG = nonzero diagonal element indicator.
!
    indiag = 0

    do k = j0, j1
      j = ja(k)
      if ( j < i ) then
        nupper = nupper + 1
      end if
      if ( j == i ) then
        indiag = 1
      end if
      dist = dist + real ( abs ( j - i ), kind = 8 )
    end do

    ndiag = ndiag + indiag
!
!  Maximum and minimum number of nonzero elements per column.
!
    lenc = j1 + 1 - j0
    if ( sym ) then
      lenc = lenc + lenr - indiag
    end if
    if ( lenc <= 0 ) then
      nzcol = nzcol + 1
    end if
    nzmaxc = max ( nzmaxc, lenc )
    nzminc = min ( nzminc, lenc )
    st = st + ( real ( lenc, kind = 8 ) - av )**2
!
!  Diagonal dominance.
!
    if ( valued ) then

      dsumc = 0.0D+00
      aii = 0.0D+00

      do k = j0, j1
        j = ja(k)
        if ( j == i ) then
          aii = abs ( a(k) )
        end if
        dsumc = dsumc + abs ( a(k) )
      end do

      dianrm = dianrm + aii * aii

      dsumr = 0.0D+00
      do k = iao(i), iao(i+1)-1
        dsumr = dsumr + abs ( ao(k) )
      end do

      if ( sym ) then
        if ( dsumr + dsumc <= 3.0D+00 * aii ) then
          nddomc = nddomc + 1
        end if
      else
        if ( dsumc <= 2.0D+00 * aii ) then
          nddomc = nddomc + 1
        end if
        if ( dsumr <= 2.0D+00 * aii ) then
          nddomr = nddomr + 1
        end if
      end if

    end if

  end do

    if ( sym ) then
      nddomr = nddomc
    end if

    nlower = nnz - nupper - ndiag
!
!  Write bandwidth info.
!
    dist = dist / real ( nnz, kind = 8 )
!
!  If NDIAG /= N, then we should correct AV and STD in symmetric case.
!
    if ( sym .and. ndiag /= n ) then
      eps = real ( ndiag - n, kind = 8 ) / real ( n, kind = 8 )
      av = av + eps
      st = st - eps * eps
    end if

    st = sqrt ( st / real ( n, kind = 8 ) )
    bndav = bndav / real ( n, kind = 8 )
!
!  Write out info.
!
    if ( sym ) then
      nupper = nlower
    end if

    write(iout, 101) av, st

    if ( nlower == 0 ) then
      write(iout, 105)
    end if

 1  continue

    if ( nupper == 0 ) then
      write(iout, 106)
    end if

    write(iout, 107) nlower
    write(iout, 108) nupper
    write(iout, 109) ndiag

    write(iout, 1020) nzmaxc, nzminc
    if ( .not. sym ) then
      write(iout, 1021) nzmaxc, nzminc
    end if

    if ( nzcol /= 0 ) then
      write(iout,116) nzcol
    end if

    if ( nzrow /= 0 ) then
      write(iout,115) nzrow
    end if
!
!  Normalize various results of above loop.
!
    ddomr = real ( nddomc, kind = 8 ) / real ( n, kind = 8 )
    ddomc = real ( nddomr, kind = 8 ) / real ( n, kind = 8 )
!
!  Symmetry and near symmetry, Frobenius norms.
!
      st = 0.0D+00
      tan = 0.0D+00
      tas = 0.0D+00
      std = 0.0D+00
      imatch = 0
!
!  Main loop for symmetry detection and Frobenius norms.
!
      do i = 1, n

         k1 = ia(i)
         k2 = iao(i)
         k1max = ia(i+1) - 1
         k2max = iao(i+1) - 1

         do k = k1, k1max
            std = std + ( dist - real ( abs ( ja(k) - i ), kind = 8 ) )**2
         end do

         if ( sym ) then
           go to 6
         end if

 5       continue

         if ( k1max < k1 .or. k2max < k2 ) then
           go to 6
         end if

         j1 = ja(k1)
         j2 = jao(k2)

         if ( j1 == j2 ) then

           imatch = imatch + 1

           if ( valued ) then
             tas = tas + ( a(k1) + ao(k2) )**2
             tan = tan + ( a(k1) - ao(k2) )**2
             st = st + a(k1)**2
           end if

         end if

         k1 = k1 + 1
         k2 = k2 + 1

         if ( j1 < j2 ) then
           k2 = k2 - 1
         end if

         if ( j2 < j1 ) then
           k1 = k1 - 1
         end if

         go to 5

 6        continue

      end do

      if ( sym ) then
        imatch = nnz
      end if

      av = real ( imatch, kind = 8 ) / real ( nnz, kind = 8 )
      std = sqrt ( std / real ( nnz, kind = 8 ) )
!
!  Maximum absolute value in A.
!
      if ( valued ) then
         amx = 0.0D+00
         ta = 0.0D+00

         do k = 1, nnz
            ta = ta + a(k)**2
            amx = max ( amx, abs ( a(k) ) )
         end do

         if ( sym ) then
            ta = sqrt ( 2.0D+00 * ta - dianrm )
            tas = ta
            tan = 0.0D+00
         else
            st = ta - st
            tas = 0.5D+00 * sqrt ( tas + st )
            tan = 0.5D+00 * sqrt ( tan + st )
            ta = sqrt ( ta )
         end if
      end if

      write (iout,103) imatch, av, dist, std
      write(iout,96)
      if ( valued ) then
         write(iout,104) ta, tas, tan, amx, ddomr, ddomc
         write (iout,96)
      end if
!
!  Bandedness- main diagonals.
!
  n2 = n + n - 1
  jao(1:n2) = 0

      do i = 1, n
         k1 = ia(i)
         k2 = ia(i+1) -1
         do k = k1, k2
            j = ja(k)
            jao(n+i-j) = jao(n+i-j) + 1
         end do
      end do

      iacc = jao(n)
      jb1 = 0
      jb2 = 0
      j = 0

 92   continue

      j = j + 1
      iacc = iacc + jao(n+j) + jao(n-j)

      if ( iacc * 100 <= nnz * 80 ) then
        jb1 = jb1 + 1
      end if

 93   continue

      if ( iacc * 100 <= nnz * 90 ) then
         jb2 = jb2 + 1
         go to 92
      end if
!
!  Write bandwidth information.
!
      write(iout,117)  ml, mu, iband, bndav

      nsky = nskyl + nskyu - n

      if ( sym ) then
        nsky = nskyl
      end if

      write(iout,1175) nsky

      write (iout,112) 2*jb2+1, 2*jb1+1
!
!  Count the number of nonzero diagonals.
!
      nzdiag = 0
      do i = 1, n2
         if ( jao(i) /= 0 ) then
           nzdiag = nzdiag + 1
         end if
      end do

      ndiag = 10
      ndiag = min ( n2, ndiag )
      itot  = 0
      ii = 0
      idiag = 0
!
!  Sort diagonals by decreasing order of weights.
!
 40       jmax = 0
      i = 1
      do k = 1, n2
         j = jao(k)
         if ( jmax <= j ) then
           i = k
           jmax = j
         end if
      end do
!
!  Permute.
!  Save offsets and accumulated count if diagonal is acceptable.
!  (if it has at least IPAR1*NNZ/100 nonzero elements)
!  Quit if no more acceptable diagonals.
!
      if ( jmax * 100 < ipar1 * nnz ) then
        go to 4
      end if

      ii = ii + 1
      ioff(ii) = i - n
      jao(i) = -jmax
      itot = itot + jmax
      dcount(ii) = real ( 100 * itot, kind = 8 ) / real ( nnz, kind = 8 )

      if ( ii < ndiag ) then
        go to 40
      end if

 4        continue
      ndiag = ii
!
!     t = real ( icount, kind = 8 ) / real ( nnz, kind = 8 )
!
      write (iout,118) nzdiag
      write (tmpst,'(10i6)') ioff(1:ndiag)
      write (iout,110) ndiag,tmpst
      write (tmpst,'(10f6.1)') dcount(1:ndiag)
      write (iout,111) tmpst
      write (iout, 96)
!
!  Determine block size if matrix is a block matrix.
!
       call blkfnd ( n, ja, ia, nblk )

      if (nblk <= 1) then
         write(iout,113)
      else
         write(iout,114) nblk
      end if
      write (iout,96)
!
!  Done.  Define all the formats.
!
 99    format (2x,38(2h *))
 96    format (6x,' *',65(1h-),'*')

 100   format( &
   6x,' *  Dimension N                                      = ', &
   i10,'  *'/ &
   6x,' *  Number of nonzero elements                       = ', &
   i10,'  *')
 101   format( &
   6x,' *  Average number of nonzero elements/Column        = ', &
   f10.4,'  *'/ &
   6x,' *  Standard deviation for above average             = ', &
   f10.4,'  *')

 1020       format( &
   6x,' *  Weight of longest column                         = ', &
   i10,'  *'/ &
   6x,' *  Weight of shortest column                        = ', &
   i10,'  *')
 1021       format( &
   6x,' *  Weight of longest row                            = ', &
   i10,'  *'/ &
   6x,' *  Weight of shortest row                           = ', &
   i10,'  *')
 117        format( &
   6x,' *  Lower bandwidth  (max: i-j, a(i,j) /= 0)       = ', &
   i10,'  *'/ &
   6x,' *  Upper bandwidth  (max: j-i, a(i,j) /= 0)       = ', &
   i10,'  *'/ &
   6x,' *  Maximum Bandwidth                                = ', &
   i10,'  *'/ &
   6x,' *  Average Bandwidth                                = ', &
   e10.3,'  *')
 1175       format( &
   6x,' *  Number of nonzeros in skyline storage            = ', &
   i10,'  *')
 103   format( &
   6x,' *  Matching elements in symmetry                    = ', &
   i10,'  *'/ &
   6x,' *  Relative Symmetry Match (symmetry=1)             = ', &
   f10.4,'  *'/ &
   6x,' *  Average distance of a(i,j)  from diag.           = ', &
   e10.3,'  *'/ &
   6x,' *  Standard deviation for above average             = ', &
   e10.3,'  *')
 104   format( &
   6x,' *  Frobenius norm of A                              = ', &
   e10.3,'  *'/ &
   6x,' *  Frobenius norm of symmetric part                 = ', &
   e10.3,'  *'/ &
   6x,' *  Frobenius norm of nonsymmetric part              = ', &
   e10.3,'  *'/ &
   6x,' *  Maximum element in A                             = ', &
   e10.3,'  *'/ &
   6x,' *  Percentage of weakly diagonally dominant rows    = ', &
   e10.3,'  *'/ &
   6x,' *  Percentage of weakly diagonally dominant columns = ', &
   e10.3,'  *')
 105        format( &
   6x,' *  The matrix is lower triangular ...       ',21x,' *')
 106        format( &
   6x,' *  The matrix is upper triangular ...       ',21x,' *')
 107        format( &
   6x,' *  Nonzero elements in strict lower part            = ', &
   i10,'  *')
 108       format( &
   6x,' *  Nonzero elements in strict upper part            = ', &
   i10,'  *')
 109       format( &
   6x,' *  Nonzero elements in main diagonal                = ', &
   i10,'  *')
 110   format(6x,' *  The ', i2, ' most important', &
           ' diagonals are (offsets)    : ',10x,'  *',/, &
   6x,' *',a61,3x,' *')

 111   format(6x,' *  The accumulated percentages they represent are ', &
   '  : ', 10x,'  *',/, &
   6x,' *',a61,3x,' *')
 112      format( &
   6x,' *  90% of matrix is in the band of width            = ', &
   i10,'  *',/, &
   6x,' *  80% of matrix is in the band of width            = ', &
   i10,'  *')
 113       format( &
   6x,' *  The matrix does not have a block structure ',19x, &
      ' *')
 114       format( &
   6x,' *  Block structure found with block size            = ', &
   i10,'  *')
 115       format( &
   6x ' *  There are zero rows. Number of such rows         = ', &
   i10,'  *')
 116       format( &
   6x ' *  There are zero columns. Number of such columns   = ', &
   i10,'  *')
 118       format( &
   6x ' *  The total number of nonvoid diagonals is         = ', &
   i10,'  *')

  return
end
subroutine diric ( nx, nint, a, ja, ia, f )

!*****************************************************************************80
!
!! DIRIC accounts for Dirichlet boundary conditions.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, real A(*), integer JA(*), IA(?+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) f(*)
  integer ia(*)
  integer ja(*)
  integer nint
  integer nc
  integer nr
  integer nx
!
!  Call extract from unary.
!
  call submat ( nx, 1, 1, nint, 1, nint, a, ja, ia, nr, nc, a, ja, ia )

  return
end
subroutine dlauny ( x, y, nodes, elmnts, nemax, nelmnt )

!*****************************************************************************80
!
!! DLAUNY is a simple, nonoptimal Delaunay triangulation code.
!
!  Discussion:
!
!    This is a simple, noptimal routine for the Delaunay triangulation
!    of a set of points in 2D.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    P K Sweby
!
!  Parameters:
!
!    Input, real X(NODES+3), Y(NODES+3), hold the coordinates of the
!    nodes, with at least 3 extra entries for workspace.
!
!    Input, integer NODES, the number of nodes.
!              
!    Output, integer ELMNTS(NEMAX,3), the nodes that form each triangular
!    element.                                    
!
!    Input, integer NEMAX, the maximum number of elements.
!
!    Output, integer NELMNT, the number of elements.                             !
  implicit none
                                                                  
  integer nemax
  integer nodes

  real ( kind = 8 ) cx
  real ( kind = 8 ) cy
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  integer elmnts(nemax,3)
  integer i1
  integer i2
  integer i3
  integer ie
  integer in
  integer j
  integer je
  integer k
  integer l
  integer match
  integer nart
  integer ndel
  integer nelmnt
  integer newel
  integer nn
  real ( kind = 8 ) pi
  real ( kind = 8 ) r2
  real ( kind = 8 ) rn2
  real ( kind = 8 ) x(nodes)
  real ( kind = 8 ) xl
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xr
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y(nodes)
  real ( kind = 8 ) yl
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yr
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) z

  pi = 4.0D+00 * atan ( 1.0D+00 )
!
!  Calculate artificial nodes NODES+i i = 1,2,3,4 and construct first
!  two (artificial) elements.
!
  xmin = minval ( x(1:nodes) )
  xmax = maxval ( x(1:nodes) )
  ymin = minval ( y(1:nodes) )
  ymax = maxval ( y(1:nodes) )

  dx = xmax - xmin
  dy = ymax - ymin

  xl = xmin - 4.0D+00 * dx
  xr = xmax + 4.0D+00 * dx
  yl = ymin - 4.0D+00 * dy
  yr = ymax + 4.0D+00 * dy

  x(nodes+1) = xl
  y(nodes+1) = yl
  x(nodes+2) = xl
  y(nodes+2) = yr
  x(nodes+3) = xr
  y(nodes+3) = yr
  x(nodes+4) = xr
  y(nodes+4) = yl

  elmnts(1,1) = nodes + 1
  elmnts(1,2) = nodes + 2
  elmnts(1,3) = nodes + 3
  elmnts(2,1) = nodes + 3
  elmnts(2,2) = nodes + 4
  elmnts(2,3) = nodes + 1

  nelmnt = 2

  do in = 1, nodes
!
!  Add one mesh point at a time and remesh locally if necessary.
!
    ndel = 0
    newel = 0

    do ie = 1, nelmnt
!
!  Is point IN inside the circumcircle of element IE?
!
      i1 = elmnts(ie,1)
      i2 = elmnts(ie,2)
      i3 = elmnts(ie,3)

      x2 = x(i2) - x(i1)
      x3 = x(i3) - x(i1)
      y2 = y(i2) - y(i1)
      y3 = y(i3) - y(i1)

      z = ( x2 * ( x2 - x3 ) + y2 * ( y2 - y3 ) ) / ( y2 * x3 - y3 * x2 )
      cx = 0.5D+00 * ( x3 - z * y3 )
      cy = 0.5D+00 * ( y3 + z * x3 )
      r2 = cx**2 + cy**2
      rn2 = ( ( x(in) - x(i1) - cx )**2 + ( y(in) - y(i1) - cy )**2 )
!
!  It is inside.  Create new elements and mark old for deletion.
!
      if ( rn2 <= r2 ) then

        do j = 1, 3
          do k = 1, 3
            elmnts(nelmnt+newel+j,k) = elmnts(ie,k)
          end do
          elmnts(nelmnt+newel+j,j) = in
        end do

        newel = newel + 3
        elmnts(ie,1)=0
        ndel = ndel + 1

      end if

    end do
!
!  If IN was inside circumcircle of more than 1 element, then we will
!  have created 2 identical new elements: delete them both.
!
    if ( 1 < ndel ) then

      do ie = nelmnt+1, nelmnt+newel-1
        do je = ie+1, nelmnt+newel

          match = 0

          do k = 1, 3
            do l = 1, 3
              if ( elmnts(ie,k) == elmnts(je,l) ) then
                match = match + 1
              end if
            end do
          end do

          if ( match == 3 ) then
            elmnts(ie,1) = 0
            elmnts(je,1) = 0
            ndel = ndel + 2
          end if

        end do
      end do

    end if
!
!  Delete any elements.
!
    nn = nelmnt + newel
    ie = 1

    do

      if ( elmnts(ie,1) == 0 ) then
        do j = ie, nn-1
          do k = 1, 3
            elmnts(j,k) = elmnts(j+1,k)
          end do
        end do
        nn = nn - 1
        ie = ie - 1
      end if

      ie = ie + 1

      if ( nn < ie ) then
        exit
      end if

    end do

    nelmnt = nn

  end do
!
!  Remove elements containing artificial nodes.
!
  ie = 1

  do

    nart = 0
    do l = 1, 3
      if ( nodes < elmnts(ie,l) ) then
        nart = nart + 1
      end if
    end do

    if ( 0 < nart ) then
      do j = ie, nn - 1
        do k = 1, 3
          elmnts(j,k) = elmnts(j+1,k)
        end do
      end do
      nelmnt = nelmnt - 1
      ie = ie - 1
    end if

    ie = ie + 1

    if ( nelmnt < ie ) then
      exit
    end if

  end do

  return
end
subroutine dnscsr ( nrow, ncol, nzmax, dns, ndns, a, ja, ia, ierr )

!*****************************************************************************80
!
!! DNSCSR converts Dense to Compressed Row Sparse format.
!
!  Discussion:
!
!    This routine converts a densely stored matrix into a row orientied
!    compactly sparse matrix.  It is the reverse of CSRDNS.
!
!    This routine does not check whether an element is small.  It considers 
!    that A(I,J) is zero only if it is exactly equal to zero.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NCOL, the column dimension of the matrix.
!
!    Input, integer NZMAX, the maximum number of nonzero elements 
!    allowed.  This should be set to be the lengths of the arrays A and JA.
!
!    Input, real DNS(NDNS,NCOL), an NROW by NCOL dense matrix.
!
!    Input, integer NDNS, the first dimension of DNS, which must be
!    at least NROW.
!
!    Output, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, integer IERR, error indicator.
!    0 means normal return;
!    I, means that the the code stopped while processing row I, because
!       there was no space left in A and JA, as defined by NZMAX.
!
  implicit none

  integer ncol
  integer ndns
  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) dns(ndns,ncol)
  integer i
  integer ia(nrow+1)
  integer ierr
  integer j
  integer ja(*)
  integer next
  integer nzmax

  ierr = 0
  next = 1
  ia(1) = 1

  do i = 1, nrow

    do j = 1, ncol

      if ( dns(i,j) /= 0.0D+00 ) then

        if ( nzmax < next ) then
          ierr = i
          return
        end if

        ja(next) = j
        a(next) = dns(i,j)
        next = next + 1

      end if

    end do

    ia(i+1) = next

  end do

  return
end
subroutine dperm ( nrow, a, ja, ia, ao, jao, iao, perm, qperm, job )

!*****************************************************************************80
!
!! DPERM permutes the rows and columns of a matrix stored in CSR format. 
!
!  Discussion:
!
!    This routine computes P*A*Q, where P and Q are permutation matrices.
!    P maps row i into row perm(i) and Q maps column j into column qperm(j).
!    A(I,J) becomes A(perm(i),qperm(j)) in the new matrix.
!
!    In the particular case where Q is the transpose of P (symmetric
!    permutation of A) then qperm is not needed.
!    note that qperm should be of length ncol (number of columns) but this
!    is not checked.
!
!    The algorithm is "in place".
!
!    The column indices may not be sorted on return even if they are
!    sorted on entry.
!
!    In case job == 2 or job == 4, a and ao are never referred to
!    and can be dummy arguments.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the order of the matrix.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, integer PERM(NROW), the permutation array for the rows: PERM(I)
!    is the destination of row I in the permuted matrix; also the destination
!    of column I in case the permutation is symmetric (JOB <= 2).
!
!    Input, integer QPERM(NROW), the permutation array for the columns.
!    This should be provided only if JOB=3 or JOB=4, that is, only in 
!    the case of a nonsymmetric permutation of rows and columns. 
!    Otherwise QPERM is a dummy argument.
!
! job      = integer indicating the work to be done:
! * job = 1,2 permutation is symmetric  Ao :== P * A * transp(P)
!             job = 1      permute a, ja, ia into ao, jao, iao
!             job = 2 permute matrix ignoring real values.
! * job = 3,4 permutation is non-symmetric  Ao :== P * A * Q
!             job = 3      permute a, ja, ia into ao, jao, iao
!             job = 4 permute matrix ignoring real values.
!
!    Output, real AO(*), JAO(*), IAO(NROW+1), the permuted matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer ia(nrow+1)
  integer iao(nrow+1)
  integer ja(*)
  integer jao(*)
  integer job
  integer locjob
  integer perm(nrow)
  integer qperm(nrow)
!
!  LOCJOB indicates whether or not real values must be copied.
!
  locjob = mod ( job, 2 )
!
!  Permute the rows first.
!
  call rperm ( nrow, a, ja, ia, ao, jao, iao, perm, locjob )
!
!  Permute the columns.
!
  locjob = 0

  if ( job <= 2 ) then
    call cperm ( nrow, ao, jao, iao, ao, jao, iao, perm, locjob )
  else
    call cperm ( nrow, ao, jao, iao, ao, jao, iao, qperm, locjob )
  end if

  return
end
subroutine dscaldg ( n, a, ja, ia, diag, job )

!*****************************************************************************80
!
!! DSCALDG scales rows by a diagonal factor.
!
!  Discussion:
!
!    This routine scales rows of a matrix by a diagonal factor DIAG.
!    DIAG is either given or to be computed.
!
!    If job = 1, we scale row I by by  +/- max |a(i,j) | and put the
!    inverse of the scaling factor in DIAG(i), where +/- is the sign of a(i,i).
!
!    If job = 2, we scale by the 2-norm of each row.
!
!    If DIAG(I) = 0, then DIAG(I) is replaced by 1.0.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, integer JOB, describes the task to be performed.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) diag(*)
  integer i
  integer ia(n+1)
  integer j
  integer ja(*)
  integer job
  integer k
  integer k1
  integer k2
  real ( kind = 8 ) t

  if ( job == 2 ) then

    do j = 1, n
      k1 = ia(j)
      k2 = ia(j+1) - 1
      t = 0.0D+00
      do k = k1, k2
        t = t + a(k) * a(k)
      end do
      diag(j) = sqrt ( t )
    end do

  else if ( job == 1 ) then

    call retmx ( n, a, ja, ia, diag )

  end if

   do j = 1, n

     if ( diag(j) /= 0.0D+00 ) then
        diag(j) = 1.0D+00 / diag(j)
     else
        diag(j) = 1.0D+00
     end if

  end do

  do i = 1, n
    t = diag(i)
    do k = ia(i), ia(i+1) -1
      a(k) = a(k) * t
    end do
  end do

  return
end
subroutine dump ( n, a, ja, ia, iout )

!*****************************************************************************80
!
!! DUMP writes the matrix to a file.
!
!  Discussion:
!
!    This routine writes the matrix to a file, one row at a time in a nice 
!    readable format.  This is a simple routine which is useful for debugging.
!
!    The output unit iout will have written in it the matrix in
!    one of two possible formats (depending on the max number of
!    elements per row. the values are output with only two digits
!    of accuracy (D9.2).
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, integer IOUT, the FORTRAN output unit number.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  integer i
  integer ia(n+1)
  integer iout
  integer ja(*)
  integer k
  integer k1
  integer k2
  integer maxr
!
!  Select mode horizontal or vertical.
!
  maxr = 0
  do i = 1, n
    maxr = max ( maxr, ia(i+1) - ia(i) )
  end do

  if ( maxr <= 8 ) then
!
!  Able to print one row across line.
!
    do i = 1, n
      write(iout,100) i
      k1 = ia(i)
      k2 = ia(i+1) - 1
      write (iout,101) ja(k1:k2)
      write (iout,102) a(k1:k2)
    end do

  else
!
!  Unable to print one row acros line.  Do three items at a time acros line.
!
     do i = 1, n
        write(iout,200) i
        k1 = ia(i)
        k2 = ia(i+1) - 1
        write (iout,201) (ja(k),a(k), k = k1, k2)
    end do

  end if

 100  format(1X,35(1h-),' row',i3,1x,35(1h-) )
 101  format(' col:',8(i5,6h     :))
 102  format(' val:',8(E9.2,2h :) )
 200  format(1h ,31(1h-),' row',i3,1x,31(1h-),/ &
         3('  columns :   values   *') )
 201  format(3(1h ,i5,6h    : ,D9.2,3h  *) )
  return
end
subroutine dvperm ( n, x, perm )

!*****************************************************************************80
!
!! DVPERM performs an in-place permutation of a real vector.
!
!  Discussion:
!
!    This routine permutes a real vector X using a permutation PERM.
!
!    On return, the vector X satisfies,
!
!      x(perm(j)) :== x(j), j = 1,2,.., n
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the length of X.
!
!    Input/output, real X(N), the vector to be permuted.
!
!    Input, integer PERM(N), the permutation.
!
  implicit none

  integer n

  integer ii
  integer init
  integer k
  integer next
  integer perm(n)
  real ( kind = 8 ) tmp
  real ( kind = 8 ) tmp1
  real ( kind = 8 ) x(n)

  init = 1
  tmp = x(init)
  ii = perm(init)
  perm(init)= -perm(init)
  k = 0
!
!  The main loop.
!
 6  continue

   k = k + 1
!
!  Save the chased element.
!
  tmp1 = x(ii)
  x(ii) = tmp
  next = perm(ii)

  if ( next < 0 ) then
    go to 65
  end if
!
!  Test for end.
!
  if ( n < k ) then
    perm(1:n) = -perm(1:n)
    return
  end if

  tmp = tmp1
  perm(ii) = -perm(ii)
  ii = next
!
!  End of the loop.
!
  go to 6
!
!  Reinitialize cycle.
!
 65   continue

  init = init + 1

  if ( n < init ) then 
    perm(1:n) = -perm(1:n)
    return
  end if

  if ( perm(init) < 0 ) then
    go to 65
  end if

  tmp = x(init)
  ii = perm(init)
  perm(init) = -perm(init)
  go to 6

end
subroutine ecn ( n, ic, ne, ia, ja, ar, nn, ierr )

!*****************************************************************************80
!
!! ECN generates sparse (square) matrices of the type E(N,C).  
!
!  Discussion:
!
!    This type of matrix has the following characteristics:
!    Symmetric, positive definite, N x N matrices with 4 in the diagonal
!    and -1 in the two sidediagonal and in the two bands at the distance
!    C from the diagonal. These matrices are similar to matrices obtained
!    from using the five point formula in the discretization of the
!    elliptic PDE.
!
!    If A is the sparse matrix of type E(N,C), then
!
!       min|A(i,j)| = 1,     max|A(i,j)| = 4
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Ernest Rothman, Cornell Theory Center
!
!  Reference:
!
!    Zahari Zlatev, Kjeld Schaumburg, Jerzy Wasniewski,
!    A testing Scheme for Subroutines Solving Large Linear Problems,
!    Computers and Chemistry, 
!    Volume 5, Number 2-3, pages 91-100, 1981.
!
!    Ole Osterby, Zahari Zlatev,
!    Direct Methods for Sparse Matrices;
!    Springer-Verlag, 1983.
!
!  Parameters:
!
!    Input, integer N, the size of the matrix.
!
!    Input, integer IC, controls the sparsity pattern.  1 < IC < N
!    is required.
!
!    Input, integer NN, the dimension of IA, JA and AR.  NN must
!    be at least NE.
!
!    Output, integer NE, the number of nonzero elements in the sparse matrix
!    of the type E(N,C). NE = 5*N - 2*IC - 2 .
!
!    Output, real AR(NN), the stored entries of the sparse matrix A.
!    NE is the number of nonzeros including a mandatory
!    diagonal entry for each row.
!
!    Output, integer IA(NN), pointers to specify rows for the stored nonzero
!    entries in AR.
!
!    Output, integer JA(NN), pointers to specify columns for the stored 
!    nonzero entries in AR.
!
!    Output, integer IERR, an error parameter, returned as zero on successful
!    execution of the subroutine.  Error diagnostics are given by means of
!    positive values of this parameter as follows:
!    1: N is out of range.
!    2: IC is out of range.
!    3: NN is out of range.
!
  implicit none

  integer nn

  real ( kind = 8 ) ar(nn)
  integer i
  integer ia(nn)
  integer ic
  integer ierr
  integer ilast
  integer it
  integer ja(nn)
  integer n
  integer ne

  ierr = 0
!
!  Check the input parameters.
!
  if ( n <= 2 ) then
    ierr = 1
    return
  end if

  if ( ic <= 1 .or. n <= ic ) then
    ierr = 2
    return
  end if

  ne = 5 * n - 2 * ic - 2

  if ( nn < ne ) then
    ierr = 3
    return
  end if
!
!  Generate the nonzero elements as well as the row and column pointers.
!
  ar(1:n) = 4.0D+00

  do i = 1, n
    ia(i) = i
    ja(i) = i
  end do

  ilast = n

  do i = 1, n-1
    it = ilast + i
    ar(it) = -1.0D+00
    ia(it) = i + 1
    ja(it) = i
  end do

  ilast = ilast + n - 1

  do i = 1, n-1
    it = ilast + i
    ar(it) = -1.0D+00
    ia(it) = i
    ja(it) = i + 1
  end do

  ilast = ilast + n - 1

  do i = 1, n-ic
    it = ilast + i
    ar(it) = -1.0D+00
    ia(it) = i + ic
    ja(it) = i
  end do

  ilast = ilast + n - ic

  do i = 1, n-ic
    it = ilast + i
    ar(it) = -1.0D+00
    ia(it) = i
    ja(it) = i + ic
  end do

  return
end
subroutine ellcsr ( nrow, coef, jcoef, ncoef, ndiag, a, ja, ia, nzmax, ierr )

!*****************************************************************************80
!
!! ELLCSR converts Ellpack/Itpack to Compressed Sparse Row.
!
!  Discussion:
!
!    This routine converts a matrix stored in Ellpack/Itpack format
!    coef-jcoef into the compressed sparse row format.  It actually checks
!    whether an entry in the input matrix is a nonzero element before
!    putting it in the output matrix.  The test does not account for small
!    values but only for exact zeros.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix A.
!
!    Input, real COEF(NCOEF,NDIAG), the values of the matrix A in 
!    Ellpack/Itpack format.
!
!    Input, integer JCOEF(NCOEF,NDIAG), the column indices of the 
!    corresponding elements in COEF.
!
!    Input, integer NCOEF, the maximum number of coefficients per diagonal.
!
!    Input, integer NDIAG, the number of active columns in COEF and JCOEF.
!    and the number of columns made available in coef.
!
!    Output, real A(NZMAX), JA(NZMAX), IA(NROW+1), the matrix, stored
!    in compressed sparse row format.
!
!    Input, integer NZMAX, the size of the arrays A and JA.
!
!    Output, integer IERR, an error flag.
!    0, means normal return.
!    nonzero, means that NZMAX is too small, and there is not enough 
!    space in A and JA to store output matrix.
!
  implicit none

  integer ncoef
  integer ndiag
  integer nrow
  integer nzmax

  real ( kind = 8 ) a(nzmax)
  real ( kind = 8 ) coef(ncoef,ndiag)
  integer i
  integer ia(nrow+1)
  integer ierr
  integer ja(nzmax)
  integer jcoef(ncoef,ndiag)
  integer k
  integer kpos

  ierr = 0
!
!  Copy elements by row.
!
  kpos = 1

  do i = 1, nrow

    do k = 1, ndiag

      if ( coef(i,k) /= 0.0D+00 ) then

        if ( nzmax < kpos ) then
          ierr = kpos
          return
        end if

        a(kpos) = coef(i,k)
        ja(kpos) = jcoef(i,k)
        kpos = kpos + 1

      end if

    end do

    ia(i+1) = kpos

  end do

  return
end
subroutine estif3 ( nel, ske, fe, det, xe, ye, xyke, ierr )

!*****************************************************************************80
!
!! ESTIF3 constructs an element stiffness matrix using 3 node triangles.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NEL, the index of the element.
!
!    Output, real SKE(3,3), the element stiffness matrix.
!
!    Output, real FE(3), the element load vector.
!
!    Input, real DET, twice the area of the triangle.
!
!    Input, real XE(3), YE(3), the coordinates of the vertices of the
!    triangle.
!
!    Input, real XYKE(2,2), the material constants KXX, KXY, KYX and KYY.
!
!    Output, integer IERR, an error flag, which is nonzero if
!    an error was detected.
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ) det
  real ( kind = 8 ) dn(3,2)
  real ( kind = 8 ) fe(3)
  integer i
  integer ierr
  integer j
  integer k
  integer l
  integer nel
  real ( kind = 8 ) ske(3,3)
  real ( kind = 8 ) t
  real ( kind = 8 ) xe(3)
  real ( kind = 8 ) xyke(2,2)
  real ( kind = 8 ) ye(3)
!
!  Initialize.
!
  area = 0.5D+00 * det
  fe(1:3) = 0.0D+00
  ske(1:3,1:3) = 0.0D+00
!
!  Get the first gradient of the shape function.
!
  call gradi3 ( nel, xe, ye, dn, det, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  do i = 1, 3
    do j = 1, 3
      t = 0.0D+00
      do k = 1, 2
        do l = 1, 2
          t = t + xyke(k,l) * dn(i,k) * dn(j,l)
        end do
      end do
      ske(i,j) = t * area
    end do
  end do

  return
end
subroutine exphes ( n, m, dt, eps, u, w, job, z, wkc, beta, errst, hh, ih, &
  x, y, indic, ierr )

!*****************************************************************************80
!
!! EXPHES computes the Arnoldi basis.
!
!  Discussion:
!
!    This routine computes the Arnoldi basis and the corresponding
!    coefficient vector in the approximation
!
!      w  ::= beta  Vm  ym
!
!    where ym = exp(- Hm *dt) * e1
!
!    To the vector exp(-A dt) w where A is an arbitary matrix and
!    w is a given input vector. In case job = 0 the arnoldi basis
!    is recomputed. Otherwise the
!    code assumes assumes that  u(*) contains an already computed
!    Arnoldi basis and computes only the y-vector (which is stored in v(*))
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer M, the dimension of the Krylov subspace.  This is also
!    the degree of the polynomial approximation to the exponential.
!
!    Input, real ( kind = 8 ) DT, a scalar by which to multiply the matrix.
!    DT can be viewed as a time step.  DT must be positive.
!
! eps   = scalar indicating the relative error tolerated for the result.
!         the code will try to compute an answer such that
!         norm2(exactanswer-approximation) / norm2(w) <= eps
!
! u      = work array of size n*(m+1) to contain the Arnoldi basis
!
! w      = real array of length n = input vector to  which exp(-A) is
!         to be applied.
!
! y     = real work array of  size (m+1)
! wkc   = COMPLEX work array of size (m+1)
!
! job      = integer. job indicator. If job <  0 then the Arnoldi
!         basis is recomputed. If 0 < job then it is assumed
!         that the user wants to use a previously computed Krylov
!         subspace but a different dt. Thus the Arnoldi basis and
!         the Hessenberg matrix Hm are not recomputed.
!        In that case the user should not modify the values of beta
!         and the matrices hh and u(n,*) when recalling phipro.
!         job = -1 : recompute basis and get an initial estimate for
!                    time step dt to be used.
!         job = 0  : recompute basis and do not alter dt.
!         job = 1  : do not recompute arnoldi basis.
!
! hh    = work array of size size at least (m+1) * m
!
! ih      = first dimension of hh as declared in the calling program.
!         m <= ih is required.
!
!  Entries specific to the matrix 
!
! diagonal storage is used :
!         a(n,ndiag) is a rectangular array with a(*,k) containing the
!         the diagonal offset by ioff(k) (negative or positive or zero)
!         i.e.,
!        a(i,jdiag) contains the element A(i,i+ioff(jdiag)) in the
!         usual dense storage scheme.
!
! a      = matrix in diagonal storage form
! ioff      = offsets  of diagonals
! ndiag = number of diagonals.
!
! on return:
!
! w2      = resulting vector w2 = exp(-A *dt) * w
! beta  = real equal to the 2-norm of w. Needed if exppro will
!         be recalled with the same Krylov subspace and a different
!         dt.
! errst = rough estimates of the 2-norm of the error.
! hh      = work array of dimension at least (m+1) x m
!
  implicit none

  integer ih
  integer m
  integer n
  integer, parameter :: ndmx = 20

  complex alp(ndmx+1)
  real ( kind = 8 ) alp0
  real ( kind = 8 ) beta
  real ( kind = 8 ) ddot
  real ( kind = 8 ) dt
  real ( kind = 8 ) eps
  real ( kind = 8 ) errst
  real ( kind = 8 ) fnorm
  integer i
  integer i0
  integer i1
  integer ierr
  integer indic
  integer job
  integer k
  integer ldg
  real ( kind = 8 ) hh(ih,ih)
  integer m1
  complex rd(ndmx+1)
  real ( kind = 8 ) rm
  real ( kind = 8 ) t
  real ( kind = 8 ) u(n,*)
  real ( kind = 8 ) w(*)
  complex wkc(ih)
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) z(m+1)

  save
!
!  Use degree 14 chebyshev all the time.
!
  if ( 3 <= indic ) then
    go to 60
  end if
!
!  Input fraction expansion of rational function.
!
  ldg = 7

  alp0 = 0.183216998528140087D-11

  alp(1)=( 0.557503973136501826D+02,-0.204295038779771857D+03)
  alp(2)=(-0.938666838877006739D+02, 0.912874896775456363D+02)
  alp(3)=( 0.469965415550370835D+02,-0.116167609985818103D+02)
  alp(4)=(-0.961424200626061065D+01,-0.264195613880262669D+01)
  alp(5)=( 0.752722063978321642D+00, 0.670367365566377770D+00)
  alp(6)=(-0.188781253158648576D-01,-0.343696176445802414D-01)
  alp(7)=( 0.143086431411801849D-03, 0.287221133228814096D-03)

  rd(1)=(-0.562314417475317895D+01, 0.119406921611247440D+01)
  rd(2)=(-0.508934679728216110D+01, 0.358882439228376881D+01)
  rd(3)=(-0.399337136365302569D+01, 0.600483209099604664D+01)
  rd(4)=(-0.226978543095856366D+01, 0.846173881758693369D+01)
  rd(5)=( 0.208756929753827868D+00, 0.109912615662209418D+02)
  rd(6)=( 0.370327340957595652D+01, 0.136563731924991884D+02)
  rd(7)=( 0.889777151877331107D+01, 0.166309842834712071D+02)
!
!  if 0 < job, skip Arnoldi process:
!
  if ( 0 < job ) then
    go to 2
  end if
!
!  Normalize vector W and put in first column of U.
!
  beta = sqrt ( ddot ( n, w, 1, w, 1 ) )

  if ( beta == 0.0D+00 ) then
    ierr = -1
    indic = 1
    return
  end if

  t = 1.0D+00 / beta
  u(1:n,1) = w(1:n) * t
!
!  The Arnoldi loop.
!
  i1 = 1

 58 continue

  i = i1
  i1 = i + 1

  x(1:n) = u(1:n,i)

  indic = 3
  return

 60   continue

  u(1:n,i1) = y(1:n)

  i0 = 1
!
!  Switch for Lanczos version.
!
! i0 = max ( 1, i-1 )

  call mgsr ( n, i0, i1, u, hh(1,i) )

  fnorm = fnorm + ddot ( i1, hh(1,i), 1, hh(1,i), 1 )

  if ( hh(i1,i) == 0.0D+00 ) then
    m = i
  end if

  if  ( i < m ) go to 58
!
!  Done with the Arnoldi loop.
!
  rm = real ( m, kind = 8 )
  fnorm = sqrt ( fnorm / rm )
!
!  Get BETA * E1 into Z.
!
  m1 = m + 1
  hh(1:m1,m1) = 0.0D+00
!
!  Compute initial DT when  0 <= JOB.
!
!
!  T = eps / beta
!
  if ( job < 0 ) then

    t = eps
    do k = 1, m-1
      t = t * ( 1.0D+00 - real ( m - k, kind = 8 ) / rm )
    end do

    t = 2.0D+00 * rm * ( t**( 1.0D+00 / rm ) )  / fnorm
    t = min ( abs ( dt ), t )
    dt = sign ( t, dt )

  end if

 2    continue

  z(1) = beta
  z(2:m1) = 0.0D+00
!
!  Get exp ( H ) * BETA * E1
!
  call hes ( ldg, m1, hh, ih, dt, z, rd, alp, alp0, wkc )
!
!  Error estimate.
!
  errst = abs ( z(m1) )

  indic = 2

  return
end
subroutine exppro ( n, m, eps, tn, u, w, x, y, indic, ierr )

!*****************************************************************************80
!
!! EXPPRO computes an approximation to the vector
!
!              w :=  exp ( - A * tn ) * w
!
! where A is an arbitary matrix and w is a given input vector
! uses a dynamic estimation of internal time advancement (dt)
!
! THIS IS A REVERSE COMMUNICATION IMPLEMENTATION.
!
!      indic = 0
!
!      do
!
!        call exppro ( n, m, eps, tn, u, w, x, y, indic )
!
!        if ( indic == 1 ) then
!          exit
!        end if
!
!        call matvec(n, x, y)     <--- user's matrix-vec. product
!                                    with x = input vector, and
!                                     y = result = A * x.
!      end do
!      .....
!
!    IM should not exceed 60 in this version  (see ih0 below)
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Reference:
!
!    E. Gallopoulos, Youcef Saad,
!    Efficient solution of parabolic equations by Krylov approximation methods,
!    RIACS Technical Report, 90-14.
!
!    Youcef Saad,
!    Analysis of some Krylov subspace approximations to the
!    matrix exponential operator,
!    RIACS Technical Report, 90-14.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer M, the dimension of the Krylov subspace.  This is also
!    the degree of the polynomial approximation to the exponential.
!
! eps   = scalar indicating the relative error tolerated for the result.
!         the code will try to compute an answer such that
!         norm2(exactanswer-approximation) / norm2(w) <= eps
!
! tn      = scalar by which to multiply matrix. (may be < 0)
!         the code will compute an approximation to exp(- tn * A) w
!         and overwrite the result onto w.
!
! u      = work array of size n*(m+1) (used to hold the Arnoldi basis )
!
! w      = real array of length n = input vector to  which exp(-A) is
!         to be applied. this is also an output argument
!
! x, y  = two real work vectors of length at least  n each.
!         see indic for usage.
!
! indic = integer used as indicator for the reverse communication.
!         in the first call enter indic = 0. See below for more.
!
! on return:
!
! w     = contains the resulting vector exp(-A * tn ) * w when
!         exppro has finished (see indic)
!
! indic = indicator for the reverse communication protocole.
!       * INDIC == 1  means that exppro has finished and w contains the
!         result.
!       * 1 < INDIC,  means that exppro has not finished and that
!         it is requesting another matrix vector product before
!         continuing. The user must compute Ax where A is the matrix
!         and x is the vector provided by exppro, and return the
!         result in y. Then exppro must be called again without
!         changing any other argument. typically this must be
!         implemented in a loop with exppro being called as long
!         indic is returned with a value /= 1.
!
! ierr  = error indicator.
!         ierr = 1 means phipro was called with indic=1 (not allowed)
!         ierr = -1 means that the input is zero the solution has been
!         unchanged.
!
  implicit none

  integer, parameter :: ih0 = 60
  integer n

  real ( kind = 8 ) beta
  real ( kind = 8 ) dtl
  real ( kind = 8 ) eps
  real ( kind = 8 ) errst
  real ( kind = 8 ) hh
  integer ierr
  integer ih
  integer indic
  integer job
  integer m
  real ( kind = 8 ) red
  real ( kind = 8 ) tcur
  real ( kind = 8 ) tn
  real ( kind = 8 ) told
  real ( kind = 8 ) u(*)
  real ( kind = 8 ) w(n)
  complex wkc(ih0)
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) z(ih0)

  save
!
!  INDIC = 3  means passing through only with result of y = A*x to EXPHES.
!  INDIC = 2  means EXPHES has finished its job.
!  INDIC = 1  means EXPPRO has finished its job (real end).
!
  ierr = 0

  if ( indic == 3 ) then
    go to 101
  end if

  if ( indic == 1 ) then
    ierr = 1
    return
  end if

  ih = ih0
  m  = min ( m, ih0 )
  tcur = 0.0D+00
  dtl = tn - tcur
  job = -1
!
!  Outer loop.
!
 100  continue
!
!  Call the exponential propagator.
!
  told = tcur

 101  continue

  call exphes ( n, m, dtl, eps, u, w, job, z, wkc, beta, errst, hh, ih, &
    x, y, indic, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  if ( 3 <= indic ) then
    return
  end if

  tcur = told + dtl
!
!  Relative error.
!
  errst = errst / beta

  if ( errst <= eps .and. ( eps / 100.0D+00 < errst .or. tcur == tn ) ) then
    go to 102
  end if
!
!  Use approximation :  new error = fact**m  * current error.
!
  red = ( 0.5D+00 * eps / errst )**( 1.0D+00 / real ( m, kind = 8 ) )
  dtl = dtl * red

  if ( abs ( tn ) < abs ( told + dtl ) ) then
    dtl = tn - told
  end if

  job = 1
  go to 101

 102  continue

  call project ( n, m, u, z, w )
  job = 0
  dtl = min ( dtl, tn-tcur )

  if ( abs ( tn ) < abs ( tcur + dtl ) ) then
    dtl = tn-tcur
  end if

  if ( abs ( tcur ) < abs ( tn ) ) then
    go to 100
  end if

  indic = 1

  return
end
subroutine expprod ( n, m, eps, tn, u, w, x, y, a, ioff, ndiag )

!*****************************************************************************80
!
!! EXPPROD computes an approximation to the vector
!
!              w :=  exp( - A * tn ) * w
!
! for matrices stored in diagonal (DIA) format.
!
! This routine constitutes an interface for the routine exppro for
! matrices stored in diagonal (DIA) format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer M, the dimension of the Krylov subspace.  This is also
!    the degree of the polynomial approximation to the exponential.
!
! see exppro for meaning of parameters eps, tn, u, w, x, y.
!
! a, ioff, and ndiag are the arguments of the matrix:
!
! a(n,ndiag) = a rectangular array with a(*,k) containing the diagonal
!              offset by ioff(k) (negative or positive or zero), i.e.,
!              a(i,jdiag) contains the element A(i,i+ioff(jdiag)) in
!              the usual dense storage scheme.
!
! ioff           = integer array containing the offsets  of the ndiag diagonals
! ndiag      = integer. the number of diagonals.
!
  implicit none

  integer n
  integer ndiag

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) eps
  integer ierr
  integer indic
  integer ioff(ndiag)
  integer m
  real ( kind = 8 ) tn
  real ( kind = 8 ) u(*)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  indic = 0

  do

    call exppro ( n, m, eps, tn, u, w, x, y, indic, ierr )

    if ( indic == 1 ) then
      exit
    end if
!
!  Matrix vector-product for diagonal storage.
!
    call oped ( n, x, y, a, ioff, ndiag )

  end do

  return
end
subroutine extbdg ( n, a, ja, ia, bdiag, nblk, ao, jao, iao )

!*****************************************************************************80
!
!! EXTBDG extracts the main diagonal blocks of a matrix.
!
!  Discussion:
!
!    The matrix is stored in compressed sparse row format.  This routine
!    puts the result into the array bdiag and the remainder in ao,jao,iao.
!
!    This version is sequential.  There is a more parallel version
!    that goes through the structure twice.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, integer NBLK, the dimension of each diagonal block.  The diagonal
!    blocks are stored in compressed format rowwise.  We store in
!    succession the I nonzeros of the I-th row after those of
!    row number I-1.
!
!    Output, real BDIAG(N,NBLK), the diagonal blocks of A.
!
!    Output, real AO(*), JAO(*), IAO(N+1), the remainder of the
!    matrix, in CSR Compressed Sparse Row format.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  real ( kind = 8 ) bdiag(*)
  integer i
  integer ia(n+1)
  integer iao(n+1)
  integer j
  integer j1
  integer j2
  integer ja(*)
  integer jao(*)
  integer jj
  integer k
  integer kb
  integer ko
  integer l
  integer ltr
  integer m
  integer nblk

  m = 1 + ( n - 1 ) / nblk

  ltr = ( ( nblk - 1 ) * nblk ) / 2
  l = m * ltr
  bdiag(1:l) = 0.0D+00
  ko = 0
  kb = 1
  iao(1) = 1

  do jj = 1, m

     j1 = ( jj - 1 ) * nblk + 1
     j2 = min ( n, j1 + nblk - 1 )

     do j = j1, j2

        do i = ia(j), ia(j+1) -1

          k = ja(i)

          if ( k < j1 ) then
            ko = ko + 1
            ao(ko) = a(i)
            jao(ko) = k
          else if ( k < j ) then
            bdiag(kb+k-j1) = a(i)
          end if

        end do

        kb = kb + j - j1
        iao(j+1) = ko + 1
      end do
  end do

  return
end
subroutine filter ( n, job, drptol, a, ja, ia, b, jb, ib, len, ierr )

!*****************************************************************************80
!
!! FILTER copies a matrix, dropping small elements.
!
!  Discussion:
!
!    The input parameter job selects a definition of small.
!
!    This module is in place. (b,jb,ib can ne the same as
!       a, ja, ia in which case the result will be overwritten).
!
!    Contributed by David Day,  Sep 19, 1989.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, integer JOB, determines strategy chosen by caller to
!    drop elements from matrix A.
!    * 1, Elements whose absolute value is less than the drop tolerance 
!    are removed.
!    * 2, Elements whose absolute value is less than the product of the 
!    drop tolerance and the Euclidean norm of the row are removed.
!    * 3, Elements whose absolute value is less that the product of the 
!    drop tolerance and the largest element in the row are removed.
!
!    Input, real DRPTOL, the drop tolerance.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, integer LEN, the amount of space in A and JA.
!
!    Output, real B(*), integer JB(*), IB(N+1), the filtered matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, integer IERR, error flag.
!    0 indicates normal return
!    0 < IERR indicates that there is'nt enough
!    space is a and ja to store the resulting matrix.
!    IERR then contains the row number where filter stopped.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) drptol
  integer ia(n+1)
  integer ib(n+1)
  integer ierr
  integer index
  integer ja(*)
  integer jb(*)
  integer job
  integer k
  integer k1
  integer k2
  integer len
  real ( kind = 8 ) loctol
  real ( kind = 8 ) norm
  integer row

  index = 1

  do row = 1, n

    k1 = ia(row)
    k2 = ia(row+1) - 1
    ib(row) = index

    if ( job == 1 ) then

      norm = 1.0D+00

    else if ( job == 2 ) then

      norm = sqrt ( sum ( a(k1:k2)**2 ) )

    else

      norm = 0.0D+00

      do k = k1, k2
        if ( norm < abs ( a(k) ) ) then
          norm = abs ( a(k) )
        end if
      end do

    end if

    loctol = drptol * norm

    do k = k1, k2

      if ( loctol < abs ( a(k) ) ) then
        if ( len < index ) then
          ierr = row
          return
        end if
        b(index) = a(k)
        jb(index) = ja(k)
        index = index + 1
      end if
    end do

  end do

  ib(n+1) = index

  return
end
subroutine gen57bl ( nx, ny, nz, nfree, na, n, a, ja, ia, iau, stencil )

!*****************************************************************************80
!
!! GEN57BL computes the sparse matrix for an elliptic operator.
!
!  Discussion:
!
!    This routine computes the sparse matrix, in compressed
!    format, associated with the discretization of the elliptic operator:
!
!    L u = delx( a . delx u ) + dely ( b . dely u) + delz ( c . delz u ) 
!        + delx ( d . u ) + dely (e . u) + delz( f . u ) + g . u
!
!    Here u is a vector of nfree componebts and each of the functions
!    a, b, c, d, e, f, g   is an (nfree x nfree) matrix depending of
!    the coordinate (x,y,z).
!    with Dirichlet Boundary conditions, on a rectangular 1-D,
!    2-D or 3-D grid using centered difference schemes.
!
!    The functions a, b, ..., g are known through the
!    routines  afunbl, bfunbl, ..., gfunbl. (user supplied) .
!
!    uses natural ordering, first x direction, then y, then z
!    mesh size h is uniform and determined by grid points
!    in the x-direction.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NX, NY, NZ, the number of nodes in the X, Y and Z
!    directions.
!
!    Input, integer NFREE, the number of degrees of freedom per node.
!
!    Output, integer N, the dimension of the matrix.
!
!    Input, integer NA, the first dimension of A as declared in the calling
!    program.  We require NFREE**2 <= NA.
!
! a, ja, ia = resulting matrix in  row-sparse block-reduced format
!           a(1:nfree**2, j ) contains a nonzero block.
!           ja(j) contains the column number of (1,1) entry of the block.
!
! iau     = integer*n containing the position of the diagonal element
!           in the a, ja, ia structure
!
! stencil =  work array of size (7,nfree**2), used to store
!            local stencils.
!
!     stencil (1:7,*) has the following meaning:
!
!     center point = stencil(1)
!     west point   = stencil(2)
!     east point   = stencil(3)
!     south point  = stencil(4)
!     north point  = stencil(5)
!     front point  = stencil(6)
!     back point   = stencil(7)
!
!
!                           st(5)
!                            |
!                            |
!                            |
!                            |          .st(7)
!                            |     .
!                            | .
!         st(2) ----------- st(1) ---------- st(3)
!                       .    |
!                   .        |
!               .            |
!            st(6)           |
!                            |
!                            |
!                           st(4)
!
!
  implicit none

  integer na

  real ( kind = 8 ) a(na,*)
  real ( kind = 8 ) h
  integer ia(*)
  integer iau(*)
  integer iedge
  integer ix
  integer iy
  integer iz
  integer ja(*)
  integer k
  integer kx
  integer ky
  integer kz
  integer n
  integer nfree
  integer nfree2
  integer node
  integer nx
  integer ny
  integer nz
  real ( kind = 8 ) stencil(7,*)

  h = 1.0D+00 / real ( nx + 1, kind = 8 )
  kx = 1
  ky = nx
  kz = nx * ny
  nfree2 = nfree * nfree
  iedge = 1
  node = 1

  do iz = 1, nz
    do iy = 1, ny
      do ix = 1, nx

        ia(node) = iedge
        call bsten ( nx, ny, nz, ix, iy, iz, nfree, stencil, h )
!
!  West
!
        if ( 1 < ix ) then
          ja(iedge) = node - kx
          do k = 1, nfree2
            a(iedge,k) = stencil(2,k)
          end do
          iedge = iedge + 1
        end if
!
!  South
!
        if ( 1 < iy ) then
          ja(iedge) = node - ky
          do k = 1, nfree2
            a(iedge,k) = stencil(4,k)
          end do
          iedge = iedge + 1
        end if
!
!  Front plane
!
        if ( 1 < iz ) then
          ja(iedge) = node - kz
          do k = 1, nfree2
            a(iedge,k) = stencil(6,k)
          end do
          iedge = iedge + 1
        end if
!
!  Center node
!
        ja(iedge) = node
        iau(node) = iedge
        a(iedge,1:nfree2) = stencil(1,1:nfree2)
        iedge = iedge + 1
!
!  Upper part
!  East
!
        if ( ix < nx ) then
          ja(iedge) = node + kx
          do k = 1, nfree2
            a(iedge,k) = stencil(3,k)
          end do
          iedge = iedge + 1
        end if
!
!  North
!
        if ( iy < ny ) then
          ja(iedge) = node + ky
          do k = 1, nfree2
            a(iedge,k) = stencil(5,k)
          end do
          iedge = iedge + 1
        end if
!
!  Back plane
!
        if ( iz < nz ) then
          ja(iedge) = node + kz
          do k = 1, nfree2
            a(iedge,k) = stencil(7,k)
          end do
          iedge = iedge + 1
        end if
!
!  Next node.
!
        node = node + 1

      end do
    end do
  end do
!
!  Change numbering of nodes so that each JA(K) will contain the
!  actual column number in the original matrix of entry (1,1) of each
!  block (K).
!
  do k = 1, iedge - 1
    ja(k) = (ja(k)-1) * nfree + 1
  end do

  n = ( node - 1 ) * nfree
  ia(node) = iedge

  return
end
subroutine gen57pt ( nx, ny, nz, a, ja, ia, iau, stencil )

!*****************************************************************************80
!
!! GEN57PT computes the compressed sparse matrix for an elliptic operator.
!
!  Discussion:
!
!    This routine computes the compressed sparse matrix discretization
!    for the elliptic operator:
!
!    L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) 
!        + d delx ( u ) + e dely (u) + f delz( u ) + g u
!
!    with Dirichlet Boundary conditions, on a rectangular 1-D,
!    2-D or 3-D grid using centered difference schemes.
!
!    The functions a, b, ..., g are known through the
!    routines  afun, bfun, ..., gfun.
!    note that to obtain the correct matrix, any function that is not
!    needed should be set to zero. For example for two-dimensional
!    problems, nz should be set to 1 and the functions cfun and ffun
!    should be zero functions.
!
!    uses natural ordering, first x direction, then y, then z
!    mesh size h is uniform and determined by grid points
!    in the x-direction.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NX, NY, NZ, the number of points in the X, Y and Z
!    directions.
!
!    Output, real A(*), integer JA(*), IA(?+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! iau     = integer IAU(N) containing the position of the diagonal element
!           in the a, ja, ia structure
!
!    Output, real STENCIL(7), used to store local stencils.
!     center point = stencil(1)
!     west point = stencil(2)
!     east point = stencil(3)
!     south point = stencil(4)
!     north point = stencil(5)
!     front point = stencil(6)
!     back point = stencil(7)
!
!
!                           st(5)
!                            |
!                            |
!                            |
!                            |          .st(7)
!                            |     .
!                            | .
!         st(2) ----------- st(1) ---------- st(3)
!                       .    |
!                   .        |
!               .            |
!            st(6)           |
!                            |
!                            |
!                           st(4)
!
!
  implicit none

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) h
  integer ia(*)
  integer iau(*)
  integer iedge
  integer ix
  integer iy
  integer iz
  integer ja(*)
  integer kx
  integer ky
  integer kz
  integer node
  integer nx
  integer ny
  integer nz
  real ( kind = 8 ) stencil(7)

  h = 1.0D+00 / real ( nx + 1, kind = 8 )
  kx = 1
  ky = nx
  kz = nx * ny
  iedge = 1
  node = 1

  do iz = 1, nz
     do iy = 1, ny
        do ix = 1, nx

           ia(node) = iedge
           call getsten ( nx, ny, nz, ix, iy, iz, stencil, h )
!
!  West
!
           if ( 1 < ix ) then
              ja(iedge)=node - kx
          a(iedge) = stencil(2)
              iedge = iedge + 1
           end if
!
!  South
!
           if ( 1 < iy ) then
              ja(iedge)=node - ky
          a(iedge) = stencil(4)
              iedge = iedge + 1
           end if
!
!  Front plane
!
           if ( 1 < iz ) then
              ja(iedge)=node - kz
          a(iedge) = stencil(6)
              iedge=iedge + 1
           end if
!
!  Center node
!
           ja(iedge) = node
           iau(node) = iedge
           a(iedge) = stencil(1)
           iedge = iedge + 1
!
!  Upper part
!  East
!
           if ( ix < nx ) then
              ja(iedge)=node + kx
          a(iedge) = stencil(3)
              iedge=iedge + 1
           end if
!
!  North
!
           if ( iy < ny ) then
              ja(iedge)=node + ky
          a(iedge) = stencil(5)
              iedge=iedge + 1
           end if
!
!  Back plane
!
           if ( iz < nz ) then
              ja(iedge)=node + kz
              a(iedge) = stencil(7)
              iedge=iedge + 1
           end if
!
!  Next node.
!
           node=node + 1

      end do
    end do
  end do

  ia(node)=iedge

  return
end
subroutine genfea ( nx, nelx, node, job, x, y, ijk, nodcode, fs, nint, &
  a, ja, ia, f, iwk, jwk, ierr, xyk )

!*****************************************************************************80
!
!! GENFEA generates finite element matrices for heat conduction problems.
!
!  Discussion:
!
!    This routine generates finite element matrices for the
!    heat conduction problem:
!
!      -Div ( K(x,y) Grad u ) = f
!      u = 0 on boundary
!
!    (with Dirichlet boundary conditions). The matrix is returned
!    assembled in compressed sparse row format. See genfeu for
!    matrices in unassembled form. The user must provide the grid,
!    (coordinates x, y and connectivity matrix ijk) as well as some
!    information on the nodes (nodcode) and the material properties
!    (the function K(x,y) above) in the form of a routine xyk.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NX, the number of nodes in the grid.
!
!    Input, integer NELX, the number of elements.
!
!    Input, integer NODE, the number of nodes per element, which
!    should be 3 for this routine.
!
! job          = integer. If job=0, it is assumed that there is no heat
!             source (i.e. fs = 0) and the right hand side
!             produced will therefore be a zero vector.
!             If job = 1 on entry then the contributions from the
!             heat source in each element are taken into account.
!
! x, y      = two real arrays containing the coordinates of the nodes.
!
!    Input, integer IJK(NODE,NELX), lists the nodes that make up 
!    each element.
!
! nodcode   = an integer array containing the boundary information for
!             each node with the following meaning.
!      nodcode(i) = 0 -->  node i is internal
!      nodcode(i) = 1 -->  node i is a boundary but not a corner point
!      nodcode(i) = 2 -->  node i is a corner node. [This node and the
!             corresponmding element are discarded.]
!
! fs          = real array of length nelx on entry containing the heat
!             source for each element (job = 1 only)
!
! xyk          = routine defining the material properties at each
!            element. Form:
!             call xyk(nel,xyke,x,y,ijk,node) with on return
!             xyke =  material constant matrices.
!            for each element nel, xyke(1,nel),xyke(2,nel)
!             and xyke(3,nel) represent the constants
!             K11, K22, and K12 at that element.
!
! on return
!
! nint          = integer. The number of active (nonboundary) nodes. Also
!             equal to the dimension of the assembled matrix.
!
!    Input, real A(*), integer JA(*), IA(?+1), the assembled matrix in CSR
!    Compressed Sparse Row format.
!
! f          = real array containing the right hand for the linears
!             system to solve.
!
!    Workspace, integer IWK(NX), JWK(NX).
!
! ierr          = integer. Error message. If (ierr /= 0) on return
!             it means that one of the elements has a negative or zero
!             area probably because of a bad ordering of the nodes
!             (see ijk above). Use the routine chkelmt to reorder
!             the nodes properly if necessary.
!
  implicit none

  integer node

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) f(*)
  real ( kind = 8 ) fs(*)
  integer ia(*)
  integer ierr
  integer ijk(node,*)
  integer indic
  integer iwk(*)
  integer ja(*)
  integer job
  integer jwk(*)
  integer nelx
  integer nint
  integer nodcode(*)
  integer nx
  real ( kind = 8 ) x(*)
  external xyk
  real ( kind = 8 ) y(*)

  ierr = 0
!
!  Take into boundary conditions to remove boundary nodes.
!
  call bound ( nx, nelx, ijk, nodcode, node, nint, jwk, x, y, f, iwk )
!
!  Assemble the matrix.
!
  call assmbo ( nx, nelx, node, ijk, nodcode, x, y, &
    a, ja, ia, f, iwk, jwk, ierr, xyk )
!
!  If applicable (JOB == 1), get heat source function.
!
  indic = 1

  if ( job == 1 ) then
    call hsourc ( indic, nx, nelx, node, x, y, ijk, fs, f )
  end if
!
!  Get the Dirichlet conditions.
!
  call diric ( nx, nint, a, ja, ia, f )

  return
end
subroutine genfeu ( nx, nelx, node, job, x, y, ijk, nodcode, fs, nint, &
  a, na, f, iwk, jwk, ierr, xyk )

!*****************************************************************************80
!
!! GENFEU generates finite element matrices for heat conduction problems.
!
!  Discussion:
!
!    This routine generates finite element matrices for the 
!    heat conduction problem:
!
!      - Div ( K(x,y) Grad u ) = f
!      u = 0 on boundary
!
!    (with Dirichlet boundary conditions). The matrix is returned
!    in unassembled form. The user must provide the grid,
!    (coordinates x, y and connectivity matrix ijk) as well as some
!    information on the nodes (nodcode) and the material properties
!    (the function K(x,y) above) in the form of a routine xyk.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NX, the number of nodes in the grid.
!
!    Input, integer NELX, the number of elements.
!
!    Input, integer NODE, the number of nodes per element, which
!    should be 3 for this routine.
!
! job          = integer. If job=0, it is assumed that there is no heat
!             source (i.e. fs = 0) and the right hand side
!             produced will therefore be a zero vector.
!             If job = 1 on entry then the contributions from the
!             heat source in each element are taken into account.
!
! na          = integer. The first dimension of the array a.
!             a is declared as an array of dimension a(na,node,node).
!
! x, y      = two real arrays containing the coordinates of the nodes.
!
!    Input, integer IJK(NODE,NELX), lists the nodes that make up 
!    each element.
!
! xyk          = routine defining the material properties at each
!            element. Form:
!             call xyk(nel,xyke,x,y,ijk,node) with on return
!             xyke = material constant matrices.
!            for each element nel, xyke(1,nel),xyke(2,nel)
!             and xyke(3,nel) represent the constants
!             K11, K22, and K12 at that element.
!
! nodcode   = an integer array containing the boundary information for
!             each node with the following meaning.
!      nodcode(i) = 0 -->  node i is internal
!      nodcode(i) = 1 -->  node i is a boundary but not a corner point
!      nodcode(i) = 2 -->  node i is a corner node. [This node and the
!             corresponmding element are discarded.]
!
! fs          = real array of length nelx on entry containing the heat
!             source for each element (job = 1 only)
!
! on return
!
! nint          = integer. The number of active (nonboundary) nodes. Also
!             equal to the dimension of the assembled matrix.
!
! a         = matrix in unassembled form. a(nel,*,*) contains the
!             element matrix for element nel.
!
! f          = real array containing the right hand for the linears
!             system to solve, in assembled form.
!
!    Workspace, integer IWK(NX), JWK(NX).
!
! ierr          = integer. Error message. If (ierr /= 0) on return
!             it means that one of the elements has a negative or zero
!             area probably because of a bad ordering of the nodes
!             (see ijk above). Use the routine chkelmt to reorder
!             the nodes properly if necessary.
!
  implicit none

  integer na
  integer node

  real ( kind = 8 ) a(na,node,node)
  real ( kind = 8 ) f(*)
  real ( kind = 8 ) fs(*)
  integer ierr
  integer ijk(node,*)
  integer indic
  integer iwk(*)
  integer job
  integer jwk(*)
  integer nelx
  integer nint
  integer nodcode(*)
  integer nx
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)
  external xyk

  ierr = 0
!
!  Take boundary conditions into account to move boundary nodes to the end.
!
  call bound ( nx, nelx, ijk, nodcode, node, nint, jwk, &
    x, y, f, iwk )
!
!  Assemble the matrix.
!
  call unassbl ( a, na, f, nx, nelx, ijk, nodcode, &
    node, x, y, ierr, xyk )
!
!  If applicable (JOB == 1), get heat source function.
!
  indic = 0

  if ( job == 1 ) then
    call hsourc ( indic, nx, nelx, node, x, y, ijk, fs, f )
  end if

  return
end
subroutine getbwd ( n, a, ja, ia, ml, mu )

!*****************************************************************************80
!
!! GETBWD gets the bandwidth of lower part and upper part of A.
!
!  Discussion:
!
!    This routine does not assume that the matrix is sorted.
!
!    ml and mu are allowed to be negative on return. This may be
!    useful since it will tell us whether a band is confined
!    in the strict  upper/lower triangular part.
!    indeed the definitions of ml and mu are
!
!      ml = max ( (i-j)  s.t. a(i,j) /= 0  )
!      mu = max ( (j-i)  s.t. a(i,j) /= 0  )
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, integer ML, MU, the lower and upper bandwidths of
!    the matrix.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  integer i
  integer ia(n+1)
  integer ja(*)
  integer k
  integer ldist
  integer ml
  integer mu

  ml = -n
  mu = -n

  do i = 1, n
     do k = ia(i), ia(i+1)-1
        ldist = i - ja(k)
        ml = max ( ml,  ldist )
        mu = max ( mu, -ldist )
     end do
  end do

  return
end
subroutine getdia ( nrow, ncol, job, a, ja, ia, len, diag, idiag, ioff )

!*****************************************************************************80
!
!! GETDIA extracts a given diagonal from a matrix stored in CSR format. 
!
!  Discussion:
!
!    The output matrix may be transformed with the diagonal removed
!    from it if desired (as indicated by job.)
!
!    Our definition of a diagonal of matrix is a vector of length nrow
!    (always) which contains the elements in rows 1 to nrow of
!    the matrix that are contained in the diagonal offset by ioff
!    with respect to the main diagonal. If the diagonal element
!    falls outside the matrix then it is defined as a zero entry.
!    Thus the proper definition of diag(*) with offset ioff is
!
!    diag(k) = a(k,ioff+k) k = 1,2,...,nrow
!    with elements falling outside the matrix being defined as zero.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NCOL, the column dimension of the matrix.
!
! job   = integer. Job indicator.  If job = 0 then
!         the matrix a, ja, ia, is not altered on return.
!         if job/=1  then getdia will remove the entries
!         collected in diag from the original matrix.
!         This is done in place.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! ioff  = integer,containing the offset of the wanted diagonal
!        the diagonal extracted is the one corresponding to the
!        entries a(i,j) with j-i = ioff.
!        thus ioff = 0 means the main diagonal
!
! on return:
!
! len   = number of nonzero elements found in diag.
!         (len <= min ( nrow, ncol-ioff ) - max ( 1, 1-ioff) + 1 )
!
! diag  = real array of length nrow containing the wanted diagonal.
!        diag contains the diagonal (a(i,j),j-i = ioff ) as defined
!         above.
!
! idiag = integer array of  length len, containing the poisitions
!         in the original arrays a and ja of the diagonal elements
!         collected in diag. A zero entry in idiag(i) means that
!         there was no entry found in row i belonging to the diagonal.
!
! a, ja,
!    ia = if job /= 0 the matrix is unchanged. otherwise the nonzero
!         diagonal entries collected in diag are removed from the
!         matrix. the structure is modified since the diagonal elements
!        are removed from a,ja,ia. Thus, the  returned matrix will
!         have len fewer elements if the diagonal is full.
!
  implicit none

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) diag(*)
  integer i
  integer ia(*)
  integer idiag(*)
  integer iend
  integer ioff
  integer istart
  integer ja(*)
  integer job
  integer k
  integer kdiag
  integer ko
  integer kold
  integer len
  integer ncol
  integer nrow

  istart = max ( 0, -ioff )
  iend = min ( nrow, ncol-ioff )
  len = 0
  idiag(1:nrow) = 0
  diag(1:nrow) = 0.0D+00
!
!  Extract the diagonal elements.
!
  do i = istart+1, iend

     do k = ia(i), ia(i+1) -1
        if ( ja(k) - i == ioff ) then
           diag(i) = a(k)
           idiag(i) = k
           len = len + 1
           exit
        end if
     end do

  end do

  if ( job == 0 .or. len == 0 ) then
    return
  end if
!
!  Rewind the structure.
!
  ko = 0

  do i = istart+1, iend

    kold = ko
    kdiag = idiag(i)

    if ( kdiag /= 0 ) then

      do k = ia(i), ia(i+1)-1
        if ( ja(k) /= kdiag ) then
          ko = ko + 1
          a(ko) = a(k)
          ja(ko) = ja(k)
        end if
      end do
      ia(i) = kold + 1
    end if

  end do
!
!  Redefine IA(NROW+1).
!
  ia(nrow+1) = ko + 1

  return
end
function getelm ( i, j, a, ja, ia, iadd, sorted )

!*****************************************************************************80
!
!! GETELM returns the element A(I,J) of a CSR matrix A. 
!
!  Discussion:
!
!    The matrix is assumed to be stored in Compressed Sparse Row (CSR) format.
!    This routine performs a binary search in the case where it is known 
!    that the elements are sorted, so that the column indices are in 
!    increasing order.  It also returns, in IADD, the address of the 
!    element A(I,J) in arrays A and JA when the search is successsful.
!    IADD is 0 if the element could not be found.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Noel Nachtigal, MIT
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer I, J, the row and column indices of the element.
!
!    Input, real A(*), integer JA(*), IA(?+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, integer IADD, the address of element A(I,J) in arrays A, JA 
!    if found, zero if not found.
!
!    Input, logical SORTED, is true if the matrix is known to have its 
!    column indices sorted in increasing order.
!
!    Output, real GETELM, the value of A(I,J).
!
  implicit none

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) getelm
  integer i
  integer ia(*)
  integer iadd
  integer ibeg
  integer iend
  integer imid
  integer j
  integer ja(*)
  integer k
  logical sorted
!
!  Initialization.
!
  iadd = 0
  getelm = 0.0D+00
  ibeg = ia(i)
  iend = ia(i+1) - 1
!
!  Case where the matrix is not necessarily sorted.
!
  if ( .not. sorted ) then
!
!  Scan the row, and exit as soon as A(I,J) is found.
!
    do k = ibeg, iend
      if ( ja(k) == j ) then
        iadd = k
        go to 20
      end if
    end do

  else
!
!  Begin binary search.  Compute the middle index.
!
 10 continue

   imid = ( ibeg + iend ) / 2
!
!  Test if found.
!
     if ( ja(imid) == j ) then
        iadd = imid
        go to 20
     end if

     if ( iend <= ibeg ) then
       go to 20
     end if
!
!  else update the interval bounds.
!
     if ( j < ja(imid) ) then
        iend = imid - 1
     else
        ibeg = imid + 1
     end if

     go to 10

  end if

 20 continue

  if ( iadd /= 0 ) then
    getelm = a(iadd)
  end if

  return
end
subroutine getl ( n, a, ja, ia, ao, jao, iao )

!*****************************************************************************80
!
!! GETL extracts the lower triangular part of a matrix.
!
!  Discussion:
!
!    This routine extracts the lower triangle of a matrix and writes the result 
!    as ao, jao, iao.  The routine is "in place" in that ao, jao, iao can be 
!    the same as a, ja, ia if desired.
!
!    The diagonal element is the last element in each row.
!    That is, the diagonal element of row I is in A(IA(I+1)-1).
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real AO(*), JAO(*), IAO(N+1), the lower triangular
!    part of the input matrix, in CSR Compressed Sparse Row format.
!
  implicit none

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer i
  integer ia(*)
  integer iao(*)
  integer ja(*)
  integer jao(*)
  integer k
  integer kdiag
  integer ko
  integer kold
  integer n
  real ( kind = 8 ) t
!
!  Inititialize KO, the pointer for the output matrix.
!
  ko = 0

  do i = 1, n

    kold = ko
    kdiag = 0

    do k = ia(i), ia(i+1) -1

      if ( ja(k) <= i ) then
        ko = ko+1
        ao(ko) = a(k)
        jao(ko) = ja(k)
        if ( ja(k) == i ) then
          kdiag = ko
        end if
      end if

    end do
!
!  Exchange.
!
    if ( kdiag /= 0 .and. kdiag /= ko ) then

      t = ao(kdiag)
      ao(kdiag) = ao(ko)
      ao(ko) = t

      k = jao(kdiag)
      jao(kdiag) = jao(ko)
      jao(ko) = k

    end if

    iao(i) = kold + 1

  end do
!
!  Redefine IAO(N+1).
!
  iao(n+1) = ko + 1

  return
end
subroutine getsten ( nx, ny, nz, kx, ky, kz, stencil, h )

!*****************************************************************************80
!
!! GETSTEN calculates the stencil for centered elliptic discretization.
!
!  Discussion:
!
!    This routine calculates the stencil for a centered difference 
!    discretization of the elliptic operator:
!
!    L u = delx ( a delx u ) 
!        + dely ( b dely u ) 
!        + delz ( c delz u ) 
!        + delx ( d u ) 
!        + dely ( e u ) 
!        + delz ( f u ) 
!        + g u
!
!    For 2-D problems, the discretization formula that is used is:
!
! h**2 * Lu == a(i+1/2,j) * {u(i+1,j) - u(i,j)} +
!             a(i-1/2,j) * {u(i-1,j) - u(i,j)} +
!              b(i,j+1/2) * {u(i,j+1) - u(i,j)} +
!              b(i,j-1/2) * {u(i,j-1) - u(i,j)} +
!              (h/2)*d(i,j) * {u(i+1,j) - u(i-1,j)} +
!              (h/2)*e(i,j) * {u(i,j+1) - u(i,j-1)} +
!              (h/2)*e(i,j) * {u(i,j+1) - u(i,j-1)} +
!               (h**2) * g(i,j)*u(i,j)
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NX, NY, NZ, the number of nodes in the X, Y and Z
!    directions.
!
!    ?, integer KX, KY, KZ, ?
!
!    Output, real STENCIL(7), ?
!
!    ?, real H, ?
!
  implicit none

  real ( kind = 8 ) afun
  real ( kind = 8 ) bfun
  real ( kind = 8 ) cfun
  real ( kind = 8 ) cntr
  real ( kind = 8 ) coeff
  real ( kind = 8 ) dfun
  real ( kind = 8 ) efun
  real ( kind = 8 ) ffun
  real ( kind = 8 ) gfun
  real ( kind = 8 ) h
  real ( kind = 8 ) hhalf
  integer kx
  integer ky
  integer kz
  integer nx
  integer ny
  integer nz
  real ( kind = 8 ) stencil(7)
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  stencil(1:7) = 0.0D+00

  hhalf = h * 0.5D+00
  x = h * real ( kx, kind = 8 )
  y = h * real ( ky, kind = 8 )
  z = h * real ( kz, kind = 8 )
  cntr = 0.0D+00
!
!  Differentiation with respect to X.
!
  coeff = afun(x+hhalf,y,z)
  stencil(3) = stencil(3) + coeff
  cntr = cntr + coeff

  coeff = afun(x-hhalf,y,z)
  stencil(2) = stencil(2) + coeff
  cntr = cntr + coeff

  coeff = dfun(x,y,z) * hhalf
  stencil(3) = stencil(3) + coeff
  stencil(2) = stencil(2) - coeff

  if ( 1 < ny ) then
!
!  Differentiation with respect to Y.
!
    coeff = bfun(x,y+hhalf,z)
    stencil(5) = stencil(5) + coeff
    cntr = cntr + coeff

    coeff = bfun(x,y-hhalf,z)
    stencil(4) = stencil(4) + coeff
    cntr = cntr + coeff

    coeff = efun(x,y,z) * hhalf
    stencil(5) = stencil(5) + coeff
    stencil(4) = stencil(4) - coeff
!
!  Differentiation with respect to Z.
!
    if ( 1 < nz ) then

      coeff = cfun(x,y,z+hhalf)
      stencil(7) = stencil(7) + coeff
      cntr = cntr + coeff

      coeff = cfun(x,y,z-hhalf)
      stencil(6) = stencil(6) + coeff
      cntr = cntr + coeff

      coeff = ffun(x,y,z) * hhalf
      stencil(7) = stencil(7) + coeff
      stencil(6) = stencil(6) - coeff

    end if

  end if
!
!  Discretization of the product by G.
!
  coeff = gfun(x,y,z)
  stencil(1) = h * h * coeff - cntr

  return
end
subroutine getu ( n, a, ja, ia, ao, jao, iao )

!*****************************************************************************80
!
!! GETU extracts the upper triangular part of a matrix.
!
!  Discussion:
!
!    The routine writes the result ao, jao, iao. 
!
!    The routine is in place in that ao, jao, iao can be the same 
!    as a, ja, ia if desired.
!
!    The diagonal element is the last element in each row.
!    i.e. in  a(ia(i+1)-1 )
!    ao, jao, iao may be the same as a, ja, ia on entry -- in which case
!    getu will overwrite the result on a, ja, ia.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Input, real AO(*), JAO(*), IAO(N+1), the upper triangular
!    part of the input matrix, in CSR Compressed Sparse Row format.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer i
  integer ia(n+1)
  integer iao(n+1)
  integer ja(*)
  integer jao(*)
  integer k
  integer kdiag
  integer kfirst
  integer ko
  real ( kind = 8 ) t

  ko = 0

  do i = 1, n

     kfirst = ko + 1
     kdiag = 0

     do k = ia(i), ia(i+1)-1

       if ( i <= ja(k) ) then
         ko = ko + 1
         ao(ko) = a(k)
         jao(ko) = ja(k)
         if ( ja(k) == i ) then
           kdiag = ko
         end if
       end if

     end do
!
!  Exchange.
!
     if ( kdiag /= 0 .and. kdiag /= kfirst ) then

       t = ao(kdiag)
       ao(kdiag) = ao(kfirst)
       ao(kfirst) = t

       k = jao(kdiag)
       jao(kdiag) = jao(kfirst)
       jao(kfirst) = k

    end if

    iao(i) = kfirst

  end do
!
!  Redefine IAO(N+1).
!
  iao(n+1) = ko + 1

  return
end
subroutine gradi3 ( nel, xe, ye, dn, det, ierr )

!*****************************************************************************80
!
!! GRADI3 constructs the first derivative of the shape functions.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NEL, the element number.
!
!    Input, real XE(3), YE(3), the coordinates of the three nodal points
!    in an element.
!
!    Output, real DN(3,2), the gradients of the shape functions.
!
!    Input, real DET, the determinant of the triangle.
!
!    Output, integer IERR, error flag, which is nonzero if an error occurred.
!
  implicit none

  real ( kind = 8 ) det
  real ( kind = 8 ) dn(3,2)
  integer ierr
  integer nel
  real ( kind = 8 ), parameter :: tol = 1.0D-17
  real ( kind = 8 ) xe(3)
  real ( kind = 8 ) ye(3)

  if ( det <= tol ) then

    ierr = 3

  else

    ierr = 0

    dn(1,1) = ( ye(2) - ye(3) ) / det
    dn(2,1) = ( ye(3) - ye(1) ) / det
    dn(3,1) = ( ye(1) - ye(2) ) / det
    dn(1,2) = ( xe(3) - xe(2) ) / det
    dn(2,2) = ( xe(1) - xe(3) ) / det
    dn(3,2) = ( xe(2) - xe(1) ) / det

  end if

  return
end
subroutine hes ( ndg, m, hh, ih, dt, y, root, coef, coef0, w2 )

!*****************************************************************************80
!
!! HES computes exp ( H dt) * y  where H = Hessenberg matrix (hh)
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, real Y(M), an arbitrary vector.
! 
!    Input, integer NDG, the number of poles as determined by GETRAT.
!
!    Input, integer M, the dimension of the Hessenberg matrix.
!
!    hh      = hessenberg matrix (real)
!
!    ih      = first dimenbsion of hh
!
!    dt      = scaling factor used for hh (see (1))
!
!    y      = real vector. on return exp(H dt ) y is computed
!         and overwritten on y.
!
!    ROOT(NDG)  = poles of the rational approximation to exp as
!         computed by getrat
!
!    coef,
!    coef0 = coefficients of partial fraction expansion
!
!    exp(t) ~ coef0 +  sum     Real [   coef(i) / (t - root(i)  ]
!                  i = 1,ndg
!
! valid for real t.
! coef0 is real, coef(*) is a complex array.
!
  implicit none

  integer ih
  integer, parameter :: mmax = 70
  integer ndg

  complex coef(*)
  real ( kind = 8 ) coef0
  real ( kind = 8 ) dt
  real ( kind = 8 ) hh(ih,*)
  complex hloc(mmax+1,mmax)
  integer i
  integer ii
  integer j
  integer m
  complex root(ndg)
  complex t
  complex w2(*)
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) yloc(mmax)
  complex zpiv
!
!  Loop associated with the poles.
!
  yloc(1:m) = y(1:m)
  y(1:m) = y(1:m) * coef0

  do ii = 1, ndg
!
!  Copy the Hessenberg matrix into a temporary array.
!
     do j = 1, m
        do i = 1, j+1
           hloc(i,j) = cmplx ( dt * hh(i,j) )
        end do
        hloc(j,j) = hloc(j,j) - root(ii)
        w2(j) = cmplx ( yloc(j) )
     end do
!
!  Forward solve.
!
     do i = 2, m
        zpiv  = hloc(i,i-1) / hloc(i-1,i-1)
        do j = i, m
           hloc(i,j) = hloc(i,j) - zpiv * hloc(i-1,j)
        end do
        w2(i) = w2(i) - zpiv * w2(i-1)
     end do
!
!  Backward solve.
!
     do i = m, 1, -1
        t = w2(i)
        do j = i+1, m
           t = t - hloc(i,j) * w2(j)
        end do
        w2(i) = t / hloc(i,i)
     end do
!
!  Accumulate result in Y.
!
     do i = 1, m
        y(i) = y(i) + coef(ii) * w2(i)
     end do

  end do

  return
end
subroutine hsourc ( indic, nx, nelx, node, x, y, ijk, fs, f )

!*****************************************************************************80
!
!! HSOURC assembles the load vector F from element contributions in FS.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    indic = indicates if f is to be assembled (1) or not (zero)
!    note: f(*) not initilazed. because might use values from boundary
!    conditions.
!
!    Input, integer NX, the number of nodes in the grid.
!
!    Input, integer NELX, the number of elements.
!
!    Input, integer NODE, the number of nodes per element, which
!    should be 3 for this routine.
!
!    Input, real X(NX), Y(NX), the coordinates of the nodes.
!
!    Input, integer IJK(NODE,NELX), lists the nodes that make up 
!    each element.
!
!    Output, real F(*), ?
!
!    Input, real FS(*), ?
!
  implicit none

  integer node

  real ( kind = 8 ) areao3
  real ( kind = 8 ) det
  real ( kind = 8 ) f(*)
  real ( kind = 8 ) fs(*)
  integer i
  integer ii
  integer ijk(node,*)
  integer indic
  integer j
  integer jnod
  integer ka
  integer nel
  integer nelx
  integer nx
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xe(3)
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) ye(3)

  jnod = 0

  do nel = 1, nelx
!
!  Get coordinates of nodal points.
!
    do i = 1, node
      j = ijk(i,nel)
      xe(i) = x(j)
      ye(i) = y(j)
    end do
!
!  Compute the determinant.
!
    det = xe(2) * ( ye(3) - ye(1) ) &
        + xe(3) * ( ye(1) - ye(2) ) &
        + xe(1) * ( ye(2) - ye(3) )

    areao3 = det / 6.0D+00
!
!  Contributions to nodes in the element.
!
    if ( indic == 0 ) then

      do ka = 1, node
        jnod = jnod + 1
        f(jnod) = fs(nel) * areao3
      end do

    else

      do ka = 1, node
        ii = ijk(ka,nel)
        f(ii) = f(ii) + fs(nel) * areao3
      end do

    end if

  end do

  return
end
subroutine ilu0 ( n, a, ja, ia, alu, jlu, ju, iw, ierr )

!*****************************************************************************80
!
!! ILU0 is an ILU(0) preconditioner.
!
!  Discussion:
!
!    Note that this has been coded in such a way that it can be used
!    with PGMRES.  Normally, since the data structure of a, ja, ia is
!    the same as that of a, ja, ia, savings can be made. In fact with
!    some definitions (not correct for general sparse matrices) all we
!    need in addition to a, ja, ia is an additional diagonal.
!    Ilu0 is not recommended for serious problems. It is only provided
!    here for comparison purposes.
!
!    It is assumed that the the elements in the input matrix are stored
!    in such a way that in each row the lower part comes first and
!    then the upper part. To get the correct ILU factorization, it is
!    also necessary to have the elements of L sorted by increasing
!    column number. It may therefore be necessary to sort the
!    elements of a, ja, ia prior to calling ilu0. This can be
!    achieved by transposing the matrix twice using csrcsc.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! on return:
!
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju        = pointer to the diagonal elements in alu, jlu.
!
! ierr        = integer indicating error code on return
!           ierr = 0 --> normal return
!           ierr = k --> code encountered a zero pivot at step k.
! work arrays:
!
! iw          = integer work array of length n.
!
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) alu(*)
  integer i
  integer ia(n+1)
  integer ierr
  integer ii
  integer iw(n)
  integer j
  integer ja(*)
  integer jcol
  integer jf
  integer jj
  integer jlu(*)
  integer jm
  integer jrow
  integer js
  integer ju(*)
  integer ju0
  integer jw
  real ( kind = 8 ) tl

  ju0 = n + 2
  jlu(1) = ju0
!
!  Initialize the work vector.
!
  iw(1:n) = 0
!
!  The main loop.
!
  do ii = 1, n

    js = ju0
!
!  Generating row II of L and U.
!
    do j = ia(ii), ia(ii+1)-1
!
!  Copy row II of A, JA, IA into row II of ALU, JLU (L/U) matrix.
!
      jcol = ja(j)

      if ( jcol == ii ) then
        alu(ii) = a(j)
        iw(jcol) = ii
        ju(ii) = ju0
      else
        alu(ju0) = a(j)
        jlu(ju0) = ja(j)
        iw(jcol) = ju0
        ju0 = ju0 + 1
      end if

    end do

    jlu(ii+1) = ju0
    jf = ju0 - 1
    jm = ju(ii) - 1
!
!  Exit if the diagonal element is reached.
!
    do j = js, jm

      jrow = jlu(j)
      tl = alu(j) * alu(jrow)
      alu(j) = tl
!
!  Perform linear combination.
!
      do jj = ju(jrow), jlu(jrow+1)-1
        jw = iw(jlu(jj))
        if ( jw /= 0 ) then
          alu(jw) = alu(jw) - tl * alu(jj)
        end if
      end do

    end do
!
!  Invert and store the diagonal element.
!
    if ( alu(ii) == 0.0D+00 ) then
      ierr = ii
      return
    end if

    alu(ii) = 1.0D+00 / alu(ii)
!
!  Reset pointer IW to zero.
!
    iw(ii) = 0
    do i = js, jf
      iw(jlu(i)) = 0
    end do

  end do

  ierr = 0
  return
end
subroutine ilut ( n, a, ja, ia, lfil, tol, alu, jlu, ju, iwk, wu, wl, jr, &
  jwl, jwu, ierr )

!*****************************************************************************80
!
!! ILUT is an ILUT preconditioner.
!
!  Discussion:
!
!    This routine carries ouot incomplete LU factorization with dual
!    truncation mechanism.  Sorting is done for both L and U. 
!
!    The dual drop-off strategy works as follows:
!
!    1) Theresholding in L and U as set by TOL.  Any element whose size
!       is less than some tolerance (relative to the norm of current
!       row in u) is dropped.
!         
!    2) Keeping only the largest lenl0+lfil elements in L and the
!       largest lenu0+lfil elements in U, where lenl0=initial number 
!       of nonzero elements in a given row of lower part of A 
!       and lenlu0 is similarly defined.
!  
!    Flexibility: one can use tol=0 to get a strategy based on keeping the
!    largest elements in each row of L and U. Taking tol /= 0 but lfil=n
!    will give the usual threshold strategy (however, fill-in is then     
!    unpredictible).                                                      
!
!    A must have all nonzero diagonal elements.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!                                                                      
!  Parameters:
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! lfil    = integer. The fill-in parameter. Each row of L and
!           each row of U will have a maximum of lfil elements
!           in addition to the original number of nonzero elements.
!           Thus storage can be determined beforehand.
!           lfil must be >= 0.
!
! iwk     = integer. The minimum length of arrays alu and jlu
!
! On return:
!
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju      = integer array of length n containing the pointers to
!           the beginning of each row of U in the matrix alu,jlu.
!
! ierr    = integer. Error message with the following meaning.
!           ierr  = 0    --> successful return.
!           ierr > 0  --> zero pivot encountered at step number ierr.
!           ierr  = -1   --> Error. input matrix may be wrong.
!                            (The elimination process has generated a
!                            row in L or U whose length is >  n.)
!           ierr  = -2   --> The matrix L overflows the array al.
!           ierr  = -3   --> The matrix U overflows the array alu.
!           ierr  = -4   --> Illegal value for lfil.
!           ierr  = -5   --> zero pivot encountered.
!
! work arrays:
!
! jr,jwu,jwl, integer work arrays of length n.
! wu, wl, real work arrays of length n+1, and n resp.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) alu(*)
  real ( kind = 8 ) fact
  integer ia(n+1)
  integer idiag
  integer ierr
  integer ii
  integer iwk
  integer j 
  integer j1
  integer j2
  integer ja(*)
  integer jj
  integer jlu(*)
  integer jpos
  integer jr(*)
  integer jrow
  integer ju(*)
  integer ju0
  integer jwl(n)
  integer jwu(n)
  integer k
  integer len
  integer lenl
  integer lenl0
  integer lenu
  integer lenu0
  integer lfil
  integer nl
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) tnorm
  real ( kind = 8 ) tol
  real ( kind = 8 ) wl(n)
  real ( kind = 8 ) wu(n)

  if ( lfil < 0 ) then
    ierr = -4
    return
  end if
!
!  Initialize JU0 (points to next element to be added to ALU, JLU)
!  and pointer.
!
  ju0 = n + 2
  jlu(1) = ju0
!
!  Integer double pointer array.
!
  jr(1:n) = 0
!
!  The main loop.
!
  do ii = 1, n

    j1 = ia(ii)
    j2 = ia(ii+1) - 1
    lenu = 0
    lenl = 0

    tnorm = 0.0D+00
    do k = j1, j2
      tnorm = tnorm + abs ( a(k) )
    end do
    tnorm = tnorm / real ( j2-j1+1, kind = 8 )
!
!  Unpack L-part and U-part of row of A in arrays WL, WU.
!
    do j = j1, j2

      k = ja(j)
      t = a(j)

      if ( tol * tnorm <= abs ( t ) ) then

        if ( k < ii ) then
          lenl = lenl + 1
          jwl(lenl) = k
          wl(lenl) = t
          jr(k) = lenl
        else
          lenu = lenu+1
          jwu(lenu) = k
          wu(lenu) = t
          jr(k) = lenu
        end if

      end if

    end do

    lenl0 = lenl
    lenu0 = lenu
    jj = 0
    nl = 0
!
!  Eliminate previous rows.
!
150 continue

    jj = jj + 1

    if ( lenl < jj ) then
      go to 160
    end if
!
!  In order to do the elimination in the correct order we need to
!  exchange the current row number with the one that has
!  smallest column number, among JJ, JJ+1, ..., LENL.
!
    jrow = jwl(jj)
    k = jj
!
!  Determine the smallest column index.
!
    do j = jj+1, lenl
       if ( jwl(j) < jrow ) then
          jrow = jwl(j)
          k = j
       end if
    end do
!
!  Exchange in JWL.
!
    j = jwl(jj)
    jwl(jj) = jrow
    jwl(k) = j
!
!  Exchange in JR.
!
    jr(jrow) = jj
    jr(j) = k
!
!  Exchange in WL.
!
    s = wl(k)
    wl(k) = wl(jj)
    wl(jj) = s

    if ( ii <= jrow ) then
      go to 160
    end if
!
!  Get the multiplier for row to be eliminated: JROW.
!
    fact = wl(jj) * alu(jrow)
    jr(jrow) = 0

    if ( abs ( fact ) * wu(n+2-jrow) <= tol * tnorm ) then
      go to 150
    end if
!
!  Combine current row and row JROW.
!
    do k = ju(jrow), jlu(jrow+1)-1
       s = fact * alu(k)
       j = jlu(k)
       jpos = jr(j)
!
!  If fill-in element and small disregard.
!
       if ( abs ( s ) < tol * tnorm .and. jpos == 0 ) then
         cycle
       end if

       if ( ii <= j ) then
!
!  Dealing with upper part.
!
          if ( jpos == 0 ) then
!
!  This is a fill-in element.
!
             lenu = lenu + 1

             if ( n < lenu ) then
               go to 995
             end if

             jwu(lenu) = j
             jr(j) = lenu
             wu(lenu) = - s
          else
!
!  No fill-in element.
!
             wu(jpos) = wu(jpos) - s
          end if
       else
!
!  Dealing with lower part.
!
          if ( jpos == 0 ) then
!
!  This is a fill-in element.
!
             lenl = lenl + 1

             if ( n < lenl ) then
               go to 995
             end if

             jwl(lenl) = j
             jr(j) = lenl
             wl(lenl) = -s
          else
!
!  No fill-in element.
!
             wl(jpos) = wl(jpos) - s
          end if
       end if

  end do

    nl = nl + 1
    wl(nl) = fact
    jwl(nl) = jrow
  go to 150
!
!  Update the L matrix.
!
 160 continue

    len = min ( nl, lenl0 + lfil )

    call bsort2 ( wl, jwl, nl, len )

    do k = 1, len

       if ( iwk < ju0 ) then
         ierr = -2
         return
       end if

       alu(ju0) =  wl(k)
       jlu(ju0) =  jwl(k)
       ju0 = ju0 + 1

    end do
!
!  Save pointer to beginning of row II of U.
!
    ju(ii) = ju0
!
!  Reset double pointer JR to zero (L-part - except first
!  JJ-1 elements which have already been reset).
!
  do k = jj, lenl
    jr(jwl(k)) = 0
  end do
!
!  Be sure that the diagonal element is first in W and JW.
!
    idiag = jr(ii)

    if ( idiag == 0 ) then
      go to 900
    end if

    if ( idiag /= 1 ) then

       s = wu(1)
       wu(j) = wu(idiag)
       wu(idiag) = s

       j = jwu(1)
       jwu(1) = jwu(idiag)
       jwu(idiag) = j

    end if

    len = min ( lenu, lenu0 + lfil )

    call bsort2 ( wu(2), jwu(2), lenu-1, len )
!
! Update the U-matrix.
!
    t = 0.0D+00

    do k = 2, len

       if ( iwk < ju0 ) then
         ierr = -3
         return
       end if

       jlu(ju0) = jwu(k)
       alu(ju0) = wu(k)
       t = t + abs ( wu(k) )
       ju0 = ju0 + 1

    end do
!
!  Save norm in WU (backwards). Norm is in fact average absolute value.
!
    wu(n+2-ii) = t / real ( len + 1, kind = 8 )
!
!  Store inverse of diagonal element of U.
!
    if ( wu(1) == 0.0D+00 ) then
      ierr = -5
      return
    end if

    alu(ii) = 1.0D+00 / wu(1)
!
!  Update pointer to beginning of next row of U.
!
  jlu(ii+1) = ju0
!
!  Reset double pointer JR to zero (U-part).
!
  do k = 1, lenu
    jr(jwu(k)) = 0
  end do

  end do

  ierr = 0

  return
!
!  Zero pivot :
!
 900    ierr = ii
    return
!
!  Incomprehensible error. Matrix must be wrong.
!
 995    ierr = -1
    return
end
subroutine infdia ( n, ja, ia, ind, idiag )

!*****************************************************************************80
!
!! INFDIA obtains information on the diagonals of A.
!
!  Discussion:
!
!    This routine finds the lengths of each of the 2*N-1 diagonals of A
!
!    It also outputs the number of nonzero diagonals found.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer JA(*), IA(N+1), the matrix information (but
!    no values) in CSR Compressed Sparse Row format.
!
!    Output, integer IND(2*N-1); The K-th entry in IND contains the number 
!    of nonzero elements in diagonal K, the numbering being from the 
!    lowermost diagonal (bottom-left).  In other words IND(K) = length 
!    of diagonal whose offset with respect to the main diagonal is = - N + K.
!
!    Output, integer IDIAG, the number of nonzero diagonals found.
!
  implicit none

  integer n

  integer i
  integer ia(n+1)
  integer idiag
  integer ind(*)
  integer j
  integer ja(*)
  integer k
  integer n2

  n2 = n+n-1
  ind(1:n2) = 0

  do i = 1, n
     do k = ia(i), ia(i+1)-1
        j = ja(k)
        ind(n+j-i) = ind(n+j-i) + 1
     end do
  end do
!
!  Count the nonzero ones.
!
  idiag = 0

  do k = 1, n2
    if ( ind(k) /= 0 ) then
      idiag = idiag + 1
    end if
  end do

  return
end
subroutine ivperm ( n, ix, perm )

!*****************************************************************************80
!
!! IVPERM performs an in-place permutation of an integer vector.
!
!  Discussion:
!
!    The integer vector ix is permuted according to the permutation 
!    array perm(*), i.e., on return, the vector x satisfies,
!
!      ix(perm(j)) :== ix(j), j = 1,2,.., n
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the length of the vector.
!
!    Input/output, integer IX(N), the vector to be permuted.
!
!    Input, integer PERM(N), the permutation.
!
  implicit none

  integer n

  integer ii
  integer init
  integer ix(n)
  integer k
  integer next
  integer perm(n)
  integer tmp
  integer tmp1

  init = 1
  tmp = ix(init)
  ii = perm(init)
  perm(init)= -perm(init)
  k = 0
!
!  Loop.
!
 6    continue

  k = k + 1
!
!  Save the chased element.
!
  tmp1 = ix(ii)
  ix(ii) = tmp
  next = perm(ii)

  if ( next < 0 ) then
    go to 65
  end if
!
!  Test for end.
!
  if ( n < k ) then
    perm(1:n) = -perm(1:n)
    return
  end if

  tmp = tmp1
  perm(ii) = -perm(ii)
  ii = next
!
!  End of loop.
!
  go to 6
!
!  Reinitilaize cycle.
!
 65   continue

  init = init + 1

  if ( n < init ) then
    perm(1:n) = -perm(1:n)
    return
  end if

  if ( perm(init) < 0 ) then
    go to 65
  end if

  tmp = ix(init)
  ii = perm(init)
  perm(init)=-perm(init)
  go to 6

end
subroutine jadcsr ( nrow, idiag, a, ja, ia, iperm, ao, jao, iao )

!*****************************************************************************80
!
!! JADSCR converts Jagged Diagonal Storage to Compressed Sparse Row.
!
!  Discussion:
!
!    This routine converts a matrix stored in the jagged diagonal format
!    to the compressed sparse row format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer IDIAG, the number of jagged diagonals in the data
!    structure A, JA, IA.
!
! a,
! ja,
! ia, input matrix in jagged diagonal format.
!
!    Input, integer IPERM(NROW), the row permutation used to obtain
!    the JAD ordering.
!
!    Output, real AO(*), integer JAO(*), IAO(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer idiag
  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer i
  integer ia(idiag+1)
  integer iao(nrow+1)
  integer iperm(nrow)
  integer j
  integer ja(*)
  integer jao(*)
  integer jj
  integer k
  integer k1
  integer kpos
  integer len
!
!  Determine first the pointers for output matrix.  Go through the
!  structure once:
!
  jao(1:nrow) = 0
!
!  Compute the lengths of each row of the output matrix.
!
  do i = 1, idiag
    len = ia(i+1) - ia(i)
    do k = 1, len
      jao(iperm(k)) = jao(iperm(k)) + 1
    end do
  end do
!
!  Permute.
!
  kpos = 1
  iao(1) = 1
  do i = 1, nrow
    kpos = kpos + jao(i)
    iao(i+1) = kpos
  end do
!
!  Copy elemnts one at a time.
!
  do jj = 1, idiag

     k1 = ia(jj) - 1
     len = ia(jj+1) - k1 - 1

     do k = 1, len
        kpos = iao(iperm(k))
        ao(kpos) = a(k1+k)
        jao(kpos) = ja(k1+k)
        iao(iperm(k)) = kpos + 1
     end do

  end do
!
!  Rewind the pointers.
!
  do j = nrow, 1, -1
    iao(j+1) = iao(j)
  end do
  iao(1) = 1

  return
end
subroutine ldsol ( n, x, y, al, jal )

!*****************************************************************************80
!
!! LDSOL solves L * x = y, for L a triangular matrix in MSR format.
!
!  Discussion:
!
!    This routine solves a non-unit lower triangular system by standard
!    (sequential) forward elimination, with the matrix stored in the
!    MSR format, with diagonal elements already inverted. 
!
!    (Otherwise do inversion, al(1:n) = 1.0/al(1:n),  before calling ldsol).
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real Y(N), the right hand side of the linear system.
!
! al,
! jal,   = Lower triangular matrix stored in Modified Sparse Row
!          format.
!
!    Output, real X(N), the solution.
!
  implicit none

  integer n

  real ( kind = 8 ) al(*)
  integer j
  integer jal(*)
  integer k
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  x(1) = y(1) * al(1)

  do k = 2, n

     t = y(k)
     do j = jal(k), jal(k+1)-1
       t = t - al(j) * x(jal(j))
     end do

     x(k) = al(k) * t

  end do

  return
end
subroutine ldsolc ( n, x, y, al, jal )

!*****************************************************************************80
!
!! LDSOLC solves L*x = y;    L = nonunit Low. Triang. MSC format
!
!  Discussion:
!
!    This routine solves a non-unit lower triangular system by standard
!    sequential forward elimination, with the matrix stored in Modified
!    Sparse Column format with diagonal elements already inverted. 
!    (otherwise do inversion, al(1:n) = 1.0/al(1:n),  before calling ldsol).
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real Y(N), the right hand side of the linear system.
!
! al,
! jal,
! ial,    = Lower triangular matrix stored in Modified Sparse Column
!           format.
!
!    Output, real X(N), the solution.
!
  implicit none

  integer n

  real ( kind = 8 ) al(*)
  integer j
  integer jal(*)
  integer k
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) t

  x(1:n) = y(1:n)

  do k = 1, n

     x(k) = x(k) * al(k)
     t = x(k)
     do j = jal(k), jal(k+1)-1
       x(jal(j)) = x(jal(j)) - t * al(j)
     end do

  end do

  return
end
subroutine ldsoll ( n, x, y, al, jal, nlev, lev, ilev )

!*****************************************************************************80
!
!! LDSOLL solves L*x = y; L = triangular. 
!
!  Discussion:
!
!    This routine uses LEVEL SCHEDULING/MSR format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real Y(N), the right hand side of the linear system.
!
! al,
! jal,   = Lower triangular matrix stored in Modified Sparse Row
!          format.
! nlev   = number of levels in matrix
! lev    = integer array of length n, containing the permutation
!          that defines the levels in the level scheduling ordering.
! ilev   = pointer to beginning of levels in lev.
!          the numbers lev(i) to lev(i+1)-1 contain the row numbers
!          that belong to level number i, in the level shcheduling
!          ordering.
!
!    Output, real X(N), the solution.
!
  implicit none

  integer n
  integer nlev

  real ( kind = 8 ) al(*)
  integer i
  integer ii
  integer ilev(nlev+1)
  integer jal(*)
  integer jrow
  integer k
  integer lev(n)
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
!
!  Outer loop goes through the levels. (SEQUENTIAL loop)
!
  do ii = 1, nlev
!
!  Next loop executes within the same level. PARALLEL loop
!
     do i = ilev(ii), ilev(ii+1)-1

        jrow = lev(i)
!
!  Compute inner product of row JROW with X.
!
        t = y(jrow)
        do k = jal(jrow), jal(jrow+1)-1
          t = t - al(k) * x(jal(k))
        end do
        x(jrow) = t * al(jrow)

    end do

  end do

  return
end
subroutine levels ( n, jal, ial, nlev, lev, ilev, levnum )

!*****************************************************************************80
!
!! LEVELS gets the level structure of a lower triangular matrix.
!
!  Discussion:
!
!    The level structure is used for level scheduling in the parallel
!    solution of triangular systems.  Strict lower matrices (e.g. unit) 
!    as well matrices with their main diagonal are accepted.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
! jal, ial =
!
! on return:
!
!    Output, integer NLEV, the number of levels found.
!
! lev      = integer array of length n containing the level
!            scheduling permutation.
! ilev     = integer array. pointer to beginning of levels in lev.
!            the numbers lev(i) to lev(i+1)-1 contain the row numbers
!            that belong to level number i, in the level scheduling
!            ordering. The equations of the same level can be solved
!            in parallel, once those of all the previous levels have
!            been solved.
! work arrays:
!
! levnum   = integer array of length n (containing the level numbers
!            of each unknown on return)
!
  implicit none

  integer n

  integer i
  integer ial(*)
  integer ilev(*)
  integer j
  integer jal(*)
  integer lev(*)
  integer levi
  integer levnum(n)
  integer nlev

  levnum(1:n) = 0
!
!  Compute level of each node.
!
  nlev = 0
  do i = 1, n
     levi = 0
     do j = ial(i), ial(i+1) - 1
        levi = max ( levi, levnum(jal(j)) )
     end do
     levi = levi + 1
     levnum(i) = levi
     nlev = max ( nlev, levi )
  end do
!
!  Set data structure.
!
  ilev(1:nlev+1) = 0
!
!  Count number of elements in each level.
!
  do j = 1, n
     i = levnum(j) + 1
     ilev(i) = ilev(i) + 1
  end do
!
!  Set up pointer for each level.
!
  ilev(1) = 1
  do j = 1, nlev
    ilev(j+1) = ilev(j) + ilev(j+1)
  end do
!
!  Determine elements of each level.
!
  do j = 1, n
     i = levnum(j)
     lev(ilev(i)) = j
     ilev(i) = ilev(i)+1
  end do
!
!  Reset pointers backwards.
!
  do j = nlev, 1, -1
    ilev(j+1) = ilev(j)
  end do

  return
end
subroutine lnkcsr ( n, a, jcol, istart, link, ao, jao, iao )

!*****************************************************************************80
!
!! LNKCSR converts linked list storage to Compressed Sparse Row format.
!
!  Discussion:
!
!    This routine translates a matrix stored in linked list storage
!    format into the compressed sparse row format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
! a      = real array of size nna containing the nonzero elements
!
! jcol      = integer array of size      nnz containing the column positions
!         of the corresponding elements in a.
!
! istart= integer array of size n poiting to the beginning of the rows.
!         istart(i) contains the position of the first element of
!         row i in data structure. (a, jcol, link).
!         if a row is empty istart(i) must be zero.
!
! link      = integer array of size nnz containing the links in the linked
!         list data structure. link(k) points to the next element
!         of the row after element ao(k), jcol(k). if link(k) = 0,
!         then there is no next element, i.e., ao(k), jcol(k) is
!         the last element of the current row.
!
!    Output, real AO(*), integer JAO(*), IAO(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer iao(n+1)
  integer ipos
  integer irow
  integer istart(n)
  integer jao(*)
  integer jcol(*)
  integer link(*)
  integer next
!
!  Determine individual bandwidths and pointers.
!
  ipos = 1
  iao(1) = ipos
!
!  Loop through all rows.
!
  do irow = 1, n
!
!  Unroll I-th row.
!
     next = istart(irow)

     do

       if ( next == 0 ) then
         exit
       end if

       jao(ipos) = jcol(next)
       ao(ipos) = a(next)
       ipos = ipos + 1
       next = link(next)

     end do

     iao(irow+1) = ipos

  end do

  return
end
subroutine lsol ( n, x, y, al, jal, ial )

!*****************************************************************************80
!
!! LSOL solves L*x = y ; L = lower unit triang. /  CSR format
!
!  Discussion:
!
!    This routine solves a unit lower triangular system by standard 
!    (sequential ) forward elimination - matrix stored in CSR format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
! y      = real array containg the right side.
!
! al,
! jal,
! ial,    = Lower triangular matrix stored in compressed sparse row
!          format.
!
!    Output, real X(N), the solution.
!
  implicit none

  integer n

  real ( kind = 8 ) al(*)
  integer ial(n+1)
  integer j
  integer jal(*)
  integer k
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  x(1) = y(1)
  do k = 2, n
    t = y(k)
    do j = ial(k), ial(k+1)-1
      t = t-al(j) * x(jal(j))
    end do
    x(k) = t
  end do

  return
end
subroutine lsolc ( n, x, y, al, jal, ial )

!*****************************************************************************80
!
!! LSOLC solves L*x = y where L = unit lower triang. CSC format
!
!  Discussion:
!
!    This routine solves a unit lower triangular system by standard 
!    (sequential ) forward elimination - matrix stored in CSC format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real Y(N), the right hand side of the linear system.
!
! al,
! jal,
! ial,    = Lower triangular matrix stored in compressed sparse column
!          format.
!
!    Output, real X(N), the solution.
!
  implicit none

  integer n

  real ( kind = 8 ) al(*)
  integer ial(*)
  integer j
  integer jal(*)
  integer k
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  x(1:n) = y(1:n)

  do k = 1, n-1
     t = x(k)
     do j = ial(k), ial(k+1)-1
       x(jal(j)) = x(jal(j)) - t * al(j)
     end do
  end do

  return
end
subroutine lusol0 ( n, y, x, alu, jlu, ju )

!*****************************************************************************80
!
!! LUSOL0 performs a forward followed by a backward solve
! for LU matrix as produced by  ILUT
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real Y(N), the right hand side of the linear system.
!
!    Output, real X(N), the solution.
!
!    ALU, JLU, JU, ...
!
  implicit none

  integer n

  real ( kind = 8 ) alu(*)
  integer i
  integer jlu(*)
  integer ju(*)
  integer k
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
!
!  Forward solve
!
  do i = 1, n
    x(i) = y(i)
    do k = jlu(i), ju(i)-1
      x(i) = x(i) - alu(k) * x(jlu(k))
    end do
  end do
!
!  Backward solve.
!
  do i = n, 1, -1
    do k = ju(i), jlu(i+1)-1
      x(i) = x(i) - alu(k) * x(jlu(k))
    end do
    x(i) = alu(i) * x(i)
  end do

  return
end
subroutine markgen ( m, n, a, ja, ia )

!*****************************************************************************80
!
!! MARKGEN is a matrix generator for a Markov random walk on a triang. grid
!
!  Discussion:
!
!    This routine generates a test matrix that models a random
!    walk on a triangular grid. This test example was used by
!    G. W. Stewart ["{SRRIT} - a FORTRAN subroutine to calculate the
!    dominant invariant subspaces of a real matrix",
!    Tech. report. TR-514, University of Maryland (1978).] and in a few
!    papers on eigenvalue problems by Y. Saad [see e.g. LAA, vol. 34,
!    pp. 269-295 (1980) ]. These matrices provide reasonably easy
!    test problems for eigenvalue algorithms. The transpose of the
!    matrix  is stochastic and so it is known that one is an exact
!    eigenvalue. One seeks the eigenvector of the transpose associated
!    with the eigenvalue unity. The problem is to calculate the
!    steady state probability distribution of the system, which is
!    the eigevector associated with the eigenvalue one and scaled in
!    such a way that the sum all the components is equal to one.
!
!    1) the code will actually compute the transpose of the
!    stochastic matrix that contains the transition probibilities.
!
!    2) It should also be possible to have a matrix generator
!    with an additional parameter (basically redefining `half' below
!    to be another parameter and changing the rest accordingly, but
!    this is not as simple as it sounds). This is not likely to provide
!    any more interesting matrices.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer M, the number of points in each direction.
!
!    Output, integer N, the dimension of the matrix (which is
!    ( M * ( M + 1 ) ) / 2.
!
!    Output, real AO(*), integer JAO(*), IAO(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) cst
  integer i
  integer ia(*)
  integer ix
  integer j
  integer ja(*)
  integer jax
  integer jmax
  integer m
  integer n
  real ( kind = 8 ) pd
  real ( kind = 8 ) pu

  cst = 0.5D+00 / real ( m - 1, kind = 8 )
!
!  IX counts the grid point (natural ordering used), i.e.,
!  the row number of the matrix.
!
  ix = 0
  jax = 1
  ia(1) = jax
!
!  Sweep the Y coordinates.
!
  do i = 1, m

     jmax = m - i + 1
!
!  Sweep X coordinates.
!
     do j = 1, jmax

        ix = ix + 1
        if ( j == jmax ) then
          go to 2
        end if
        pd = cst * real ( i+j-1, kind = 8 )
!
!  north
!
        a(jax) = pd
        if ( i == 1 ) then
          a(jax) = a(jax) + pd
        end if
        ja(jax) = ix + 1
        jax = jax+1
!
!  east
!
        a(jax) = pd
        if ( j == 1 ) then
          a(jax) = a(jax) + pd
        end if
        ja(jax) = ix + jmax
        jax = jax + 1
!
!  south
!
 2      continue

        pu = 0.5D+00 - cst * real ( i+j-3, kind = 8 )

        if ( 1 < j ) then
           a(jax) = pu
           ja(jax) = ix - 1
           jax = jax + 1
        end if
!
!  west
!
        if ( 1 < i ) then
           a(jax) = pu
           ja(jax) = ix - jmax - 1
           jax = jax + 1
        end if

        ia(ix+1) = jax

    end do

  end do

  n = ix

  return
end
subroutine matrf2 ( m, n, c, index, alpha, nn, nz, a, snr, rnr, fejlm )

!*****************************************************************************80
!
!! MATRF2 generates sparse (rectangular or square) matrices.
!
!  Discussion:
!
!    The dimensions of the matrix and the average number of nonzero
!    elements per row can be specified by the user. Moreover, the user
!    can also change the sparsity pattern and the condition number of the
!    matrix. The non-zero elements of the desired matrix will be
!    accumulated (in an arbitrary order) in the first NZ positions of
!    array A. The column and the row numbers of the non-zero element
!    stored in A(I), I = 1,...,NZ, will be found in SNR(I) and RNR(I),
!    respectively. The matrix generated by this routine is of the
!    class F(M,N,C,R,ALPHA) (see reference).
!
!    If A is the sparse matrix of type F(M,N,C,R,ALPHA), then
!
!      min |A(i,j)| = 1/ALPHA,
!
!      max |A(i,j)| = max ( INDEX*N - N, 10*ALPHA ).
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Zahari Zlatev, Kjeld Schaumburg, Jerzy Wasniewski
!
!  Reference:
!
!    Zahari Zlatev, Kjeld Schaumburg, Jerzy Wasniewski,
!    A testing Scheme for Subroutines Solving Large Linear Problems,
!    Computers and Chemistry, 
!    Volume 5, Number 2-3, pages 91-100, 1981.
!
!  Parameters:
!
!   INPUT PARAMETERS
!
!   M    - Integer. The number of rows in the desired matrix.
!          N < M+1 < 9000001 must be specified.
!
!   N    - Integer. The number of columns in the desired matrix.
!          21 < N < 9000001 must be specified.
!
!   C    - Integer. The sparsity pattern can be changed by means of this
!          parameter.  10 < C < N-10  must be specified.
!
!   INDEX - Integer.  The average number of non-zero elements per row in
!           the matrix will be equal to INDEX.
!           1 < INDEX < N-C-8 must be specified.
!
!   ALPHA - Real. The condition number of the matrix can be changed
!           BY THIS PARAMETER. ALPHA > 0.0 MUST BE SPECIFIED.
!           If ALPHA is approximately equal to 1.0 then the generated
!           matrix is well-conditioned. Large values of ALPHA will
!           usually produce ill-conditioned matrices. Note that no
!           round-off errors during the computations in this routine
!           are made if ALPHA = 2**I (where I is an arbitrary integer
!           which produces numbers in the machine range).
!
!   Input, integer NN, the length of arrays A, RNR, and SNR.
!   INDEX*M+109 < NN < 9000001 must be specified.
!
!   Output, integer NZ, the number of nonzero elements in the matrix.
!
!   Output, real A(NN), the nonzero elements of the matrix,
!   accumulated in the first NZ locations of array A.
!
!   Output, integer SNR(NN), the column number of the non-zero element
!   kept in A(I), I = 1,...NZ.
!
!   Output, integer RNR(NN), the row number of the non-zero element
!   kept in A(I).
!
!   Output, integer FEJLM, error indicator.
!   0, indicates that the call is successful.
!   1, N is out of range.
!   2, M is out of range.
!   3, C is out of range.
!   4, INDEX is out of range.
!   5, NN is out of range.
!   7, ALPHA is out of range.
!
  implicit none

  integer nn

  real ( kind = 8 ) a(nn)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) alpha1
  integer c
  integer fejlm
  integer i
  integer index
  integer index1
  integer j
  integer j1
  integer k
  integer m
  integer m1
  integer m2
  integer n
  integer n2
  integer nz
  integer nz1
  integer rnr(nn)
  integer rr1
  integer rr2
  integer rr3
  integer snr(nn)

  m1 = m
  fejlm = 0
  nz1 = index * m + 110
  k = 1
  alpha1 = alpha
  index1 = index - 1
!
!  Check the parameters.
!
  if ( n < 22 ) then
    fejlm = 1
    return
  end if

  if ( 9000000 < n ) then
    fejlm = 1
    return
  end if

  if ( n < m ) then
    fejlm = 2
    return
  end if

  if ( 9000000 < m ) then
    fejlm = 2
    return
  end if

  if ( c < 11 ) then
    fejlm = 3
    return
  end if

  if ( n-c < 11 ) then
    fejlm = 3
    return
  end if

  if ( index < 1 ) then
    fejlm = 4
    return
  end if

  if ( n - c - index < 9 ) then
    fejlm = 4
    return
  end if

  if ( nn < nz1 ) then
    fejlm = 5
    return
  end if

  if ( 9000000 < nn ) then
    fejlm = 5
    return
  end if

  if ( alpha <= 0.0D+00 ) then
    fejlm = 6
    return
  end if
!
!  End of the error check.  Begin to generate the non-zero elements of
!  the required matrix.
!
  a(1:n) = 1.0D+00

  do i = 1, n
    snr(i) = i
  end do

  do i = 1, n
    rnr(i) = i
  end do

  nz = n
  j1 = 1

  do j = 1, index1

    j1 = -j1

    do i = 1, n

      a(nz+i) = real ( j1 * j * i, kind = 8 )

      if ( i + c + j - 1 <= n ) then
        snr(nz+i) = i + c + j - 1
      end if

      if ( i + c + j - 1 > n ) then
        snr(nz+i) = c + i + j - 1 - n
      end if

      rnr(nz+i) = i

    end do

    nz = nz + n

  end do

  rr1 = 10
  rr2 = nz
  rr3 = 1

  do

    do i = 1, rr1
      a(rr2+i) = alpha * real ( i, kind = 8 )
      snr(rr2+i) = n - rr1 + i
      rnr(rr2+i) = rr3
    end do

    if ( rr1 == 1 ) then
      exit
    end if

    rr2 = rr2 + rr1
    rr1 = rr1 - 1
    rr3 = rr3 + 1

  end do

  nz = nz + 55

  do

    m1 = m1 - n
    alpha = 1.0D+00 / alpha

    if ( m1 <= 0 ) then
      exit
    end if

    n2 = k * n

    if ( n <= m1 ) then
      m2 = n
    end if

    if ( m1 < n ) then
      m2 = m1
    end if

    do i = 1, m2
      a(nz+i) = alpha * real ( k + 1, kind = 8 )
      snr(nz+i) = i
      rnr(nz+i) = n2 + i
    end do

    nz = nz + m2

    j1 = 1

    do j = 1, index1
 
      j1 = -j1

      do i = 1, m2

        a(nz+i) = alpha * real ( j * j1, kind = 8 ) &
          * ( real ( ( k + 1 ) * i, kind = 8 ) + 1.0D+00 )

        if ( i + c + j - 1 <= n ) then
          snr(nz+i) = i + c + j - 1
        end if

        if ( n < i + c + j - 1 ) then
          snr(nz+i) = c + i + j - 1 - n
        end if

        rnr(nz+i) = n2 + i

      end do
      nz = nz + m2
    end do

    k = k + 1
  end do


  alpha = 1.0D+00 / alpha1
  rr1 = 1
  rr2 = nz

  do

    do i = 1, rr1
      a(rr2+i) = alpha * real ( rr1 + 1 - i, kind = 8 )
      snr(rr2+i) = i
      rnr(rr2+i) = n - 10 + rr1
    end do

    if ( rr1 == 10 ) then
      exit
    end if

    rr2 = rr2 + rr1
    rr1 = rr1 + 1

  end do

  nz = nz + 55
  alpha = alpha1

  return
end
subroutine mgsr ( n, i0, i1, ss, r )

!*****************************************************************************80
!
!! MGSR is a modified Gram - Schmidt with partial reorthogonalization.
!
!  Discussion:
!
!    the vector ss(*,i1) is
!    orthogonalized against the first i vectors  of ss  (which  are  already
!    orthogonal).  the coefficients of the orthogonalization are returned in
!    the array r
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
  implicit none

  integer n

  real ( kind = 8 ) ddot
  real ( kind = 8 ) hinorm
  integer i
  integer i0
  integer i1
  integer it
  integer j
  real ( kind = 8 ) r(*)
  real ( kind = 8 ) ss(n,*)
  real ( kind = 8 ) t

  r(1:i1) = 0.0D+00
  i = i1 - 1

  do it = 1, 2

    hinorm = 0.0D+00

    if ( 0 < i ) then

      do j = i0, i
        t = ddot ( n, ss(1,j), 1, ss(1,i1), 1 )
        hinorm = hinorm + t**2
        r(j) = r(j) + t
        call daxpy ( n, -t, ss(1,j), 1, ss(1,i1), 1 )
      end do

      t = ddot ( n, ss(1,i1), 1, ss(1,i1), 1 )

    end if
!
!  Test for reorthogonalization.  See Daniel et. al.
!  Two reorthogonalizations allowed.
!
    if ( hinorm < t * 10.0D+00 ) then
      exit
    end if

  end do

  t = sqrt ( t )
  r(i1)= t

  if ( t == 0.0D+00 ) then
    return
  end if

  t = 1.0D+00 / t

  ss(1:n,i1) = ss(1:n,i1) * t

  return
end
subroutine milu0 ( n, a, ja, ia, alu, jlu, ju, iw, ierr )

!*****************************************************************************80
!
!! MILU0 is a simple milu(0) preconditioner.
!
!  Discussion:
!
!    Note that this has been coded in such a way that it can be used
!    with pgmres. Normally, since the data structure of a, ja, ia is
!    the same as that of a, ja, ia, savings can be made. In fact with
!    some definitions (not correct for general sparse matrices) all we
!    need in addition to a, ja, ia is an additional diagonal.
!    Ilu0 is not recommended for serious problems. It is only provided
!    here for comparison purposes.
!
!    It is assumed that the the elements in the input matrix are ordered
!    in such a way that in each row the lower part comes first and
!    then the upper part. To get the correct ILU factorization, it is
!    also necessary to have the elements of L ordered by increasing
!    column number. It may therefore be necessary to sort the
!    elements of a, ja, ia prior to calling milu0. This can be
!    achieved by transposing the matrix twice using csrcsc.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! on return:
!
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju        = pointer to the diagonal elements in alu, jlu.
!
!    Workspace, integer IW(N).
!
! ierr        = integer indicating error code on return
!           ierr = 0 --> normal return
!           ierr = k --> code encountered a zero pivot at step k.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) alu(*)
  integer i
  integer ia(n+1)
  integer ierr
  integer ii
  integer iw(n)
  integer j
  integer ja(*)
  integer jcol
  integer jf
  integer jj
  integer jlu(*)
  integer jm
  integer jrow
  integer js
  integer ju(*)
  integer ju0
  integer jw
  real ( kind = 8 ) s
  real ( kind = 8 ) tl

  ju0 = n + 2
  jlu(1) = ju0
!
!  Initialize work vector.
!
  iw(1:n) = 0
!
!  The main loop.
!
  do ii = 1, n

    js = ju0
!
!  Generating row II or L and U.
!
    do j = ia(ii), ia(ii+1)-1
!
!  Copy row II of A, JA, IA into row II of ALU, JLU (L/U) matrix.
!
      jcol = ja(j)

      if ( jcol == ii ) then
        alu(ii) = a(j)
        iw(jcol) = ii
        ju(ii) = ju0
      else
        alu(ju0) = a(j)
        jlu(ju0) = ja(j)
        iw(jcol) = ju0
        ju0 = ju0 + 1
      end if

    end do

    jlu(ii+1) = ju0
    jf = ju0 - 1
    jm = ju(ii) - 1
!
!  S accumulates fill-in values.
!
    s = 0.0D+00

    do j = js, jm

      jrow = jlu(j)
      tl = alu(j) * alu(jrow)
      alu(j) = tl
!
!  Perform linear combination.
!
      do jj = ju(jrow), jlu(jrow+1)-1
        jw = iw(jlu(jj))
        if ( jw /= 0 ) then
          alu(jw) = alu(jw) - tl * alu(jj)
        else
          s = s + tl * alu(jj)
        end if
      end do

    end do
!
!  Invert and store diagonal element.
!
    alu(ii) = alu(ii) - s

    if ( alu(ii) == 0.0D+00 ) then
      ierr = ii
      return
    end if

    alu(ii) = 1.0D+00 / alu(ii)
!
!  Reset pointer IW to zero.
!
    iw(ii) = 0
    do i = js, jf
      iw(jlu(i)) = 0
    end do

  end do

  ierr = 0

  return
end
subroutine msrcsr ( n, a, ja, ao, jao, iao, wk )

!*****************************************************************************80
!
!! MSRCSR converts Modified Sparse Row to Compressed Sparse Row.
!
!  Discussion:
!
!    This routine converts a compressed matrix using a separated diagonal
!    (modified sparse row format) in the Compressed Sparse Row format.
!
!    does not check for zero elements in the diagonal.
!
!    This is an "in place" algorithm (see a, ja, ia).
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
! ao, jao  = sparse matrix in msr sparse storage format
!           see routine csrmsr for details
!
! on return :
!
! a, ja, ia = matrix in csr format. note that the
!           algorithm is in place: ao, jao can be the same
!            as a, ja, in which case it will be overwritten on it
!            upon return.
!
!             here nnz = number of nonzero elements+1
!
!    Workspace, real WK(N).
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  logical added
  real ( kind = 8 ) ao(*)
  integer iao(n+1)
  integer idiag
  integer ii
  integer iptr
  integer j
  integer ja(*)
  integer jao(*)
  integer k
  real ( kind = 8 ) wk(n)

  wk(1:n) = a(1:n)

  iao(1) = 1
  iptr = 1

  do ii = 1, n

    added = .false.
    idiag = iptr + ( ja(ii+1) - ja(ii) )

    do k = ja(ii), ja(ii+1)-1

      j = ja(k)

      if ( j < ii ) then
        ao(iptr) = a(k)
        jao(iptr) = j
        iptr = iptr + 1
      else if ( added ) then
        ao(iptr) = a(k)
        jao(iptr) = j
        iptr = iptr + 1
      else
!
!  Add diagonal element.  Only reserve a position for it.
!
        idiag = iptr
        iptr = iptr + 1
        added = .true.
!
!  Then other elements.
!
        ao(iptr) = a(k)
        jao(iptr) = j
        iptr = iptr + 1
      end if

    end do

    ao(idiag) = wk(ii)
    jao(idiag) = ii
    if ( .not. added ) then
      iptr = iptr + 1
    end if
    iao(ii+1) = iptr

  end do

  return
end
subroutine ope ( n, x, y, a, ja, ia )

!*****************************************************************************80
!
!! OPE sparse matrix * vector multiplication
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real X(N), the vector to be multiplied.
!
!    Output, real Y(N), the product A * X.
!
!    Input, real A(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  integer i
  integer ia(n+1)
  integer ja(*)
  integer k
  integer k1
  integer k2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) -1
    y(i) = 0.0D+00
    do k = k1, k2
      y(i) = y(i) + a(k) * x(ja(k))
    end do
  end do

  return
end
subroutine pgmres ( n, im, rhs, sol, vv, eps, maxits, iout, &
  aa, ja, ia, alu, jlu, ju, ierr )

!*****************************************************************************80
!
!! PGMRES is an ILUT - Preconditioned GMRES solver.
!                                                                      
!  Discussion:
!
!    This is a simple version of the ILUT preconditioned GMRES algorithm. 
!    The ILUT preconditioner uses a dual strategy for dropping elements   
!    instead  of the usual level of-fill-in approach. See details in ILUT 
!    subroutine documentation. PGMRES uses the L and U matrices generated 
!    from the subroutine ILUT to precondition the GMRES algorithm.        
!    The preconditioning is applied to the right. The stopping criterion  
!    utilized is based simply on reducing the residual norm by epsilon.   
!    This preconditioning is more reliable than ilu0 but requires more    
!    storage. It seems to be much less prone to difficulties related to   
!    strong nonsymmetries in the matrix. We recommend using a nonzero tol 
!    (tol=.005 or .001 usually give good results) in ILUT. Use a large    
!    lfil whenever possible (e.g. lfil = 5 to 10). The higher lfil the    
!    more reliable the code is. Efficiency may also be much improved.     
!    Note that lfil=n and tol=0.0 in ILUT  will yield the same factors as 
!    Gaussian elimination without pivoting.                               
!                                                                      
!    ILU(0) and MILU(0) are also provided for comparison purposes         
!    USAGE: first call ILUT or ILU0 or MILU0 to set up preconditioner and 
!    then call pgmres.                                                    
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad              
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer IM, the size of the Krylov subspace.  IM should not 
!    exceed 50 in this version.  This restriction can be reset by changing
!    the parameter command for KMAX below.
!                                                
!    Input/output, real RHS(N), on input, the right hand side vector.
!    On output, the information in this vector has been destroyed.                                      
!
! sol   == real vector of length n containing an initial guess to the  
!          solution on input. approximate solution on output           
!
! eps   == tolerance for stopping criterion. process is stopped        
!          as soon as ( ||.|| is the euclidean norm):                  
!          || current residual||/||initial residual|| <= eps           
!
! maxits== maximum number of iterations allowed                        
!
! iout  == output unit number number for printing intermediate results 
!          if (iout <= 0) nothing is printed out.                    
!                                                                      
!    Input, real AA(*), integer JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!                                                                      
! alu,jlu== A matrix stored in Modified Sparse Row format containing   
!           the L and U factors, as computed by routine ilut.       
!                                                                      
! ju     == integer array of length n containing the pointers to       
!           the beginning of each row of U in alu, jlu as computed     
!           by routine ILUT.                                        
!                                                                      
! on return:                                                           
!                                                          
! sol   == contains an approximate solution (upon successful return).  
! ierr  == integer. Error message with the following meaning.          
!          ierr = 0 --> successful return.                            
!          ierr = 1 --> convergence not achieved in itmax iterations. 
!          ierr =-1 --> the initial guess seems to be the exact        
!                       solution (initial residual computed was zero) 
!                                                                      
! work arrays:                                                        
!                                                       
! vv    == work array of length  n x (im+1) (used to store the Arnoli  
!          basis)                                                      
!
  implicit none

  integer, parameter :: kmax = 50
  integer n

  real ( kind = 8 ) aa(*)
  real ( kind = 8 ) alu(*)
  real ( kind = 8 ) c(kmax)
  real ( kind = 8 ) ddot
  real ( kind = 8 ) eps
  real ( kind = 8 ) eps1
  real ( kind = 8 ), parameter :: epsmac = 1.0D-16
  real ( kind = 8 ) gam
  real ( kind = 8 ) hh(kmax+1,kmax)
  integer i
  integer i1
  integer ia(n+1)
  integer ierr
  integer ii
  integer im
  integer iout
  integer its
  integer j
  integer ja(*)
  integer jj
  integer jlu(*)
  integer ju(*)
  integer k
  integer k1
  integer maxits
  integer n1
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) ro
  real ( kind = 8 ) rs(kmax+1)
  real ( kind = 8 ) s(kmax)
  real ( kind = 8 ) sol(n)
  real ( kind = 8 ) t
  real ( kind = 8 ) vv(n,*)
!
!  Arnoldi size should not exceed KMAX=50 in this version.
!  To reset modify parameter KMAX accordingly.
!
  n1 = n + 1
  its = 0
!
!  Outer loop starts here.
!  Compute initial residual vector.
!
  call ope ( n, sol, vv, aa, ja, ia )

  vv(1:n,1) = rhs(1:n) - vv(1:n,1)

  do

    ro = sqrt ( ddot ( n, vv, 1, vv, 1 ) )

    if ( 0 < iout .and. its == 0 ) then
      write(iout, 199) its, ro
    end if

    if ( ro == 0.0D+00 ) then
      ierr = -1
      exit
    end if

    t = 1.0D+00 / ro
    vv(1:n,1) = vv(1:n,1) * t

    if ( its == 0 ) then
      eps1 = eps * ro
    end if
!
!  Initialize first term of RHS of Hessenberg system.
!
     rs(1) = ro
     i = 0

 4   continue

     i = i + 1
     its = its + 1
     i1 = i + 1
     call lusol0 ( n, vv(1,i), rhs, alu, jlu, ju )
     call ope ( n, rhs, vv(1,i1), aa, ja, ia )
!
!  Modified Gram - Schmidt.
!
     do j = 1, i
       t = ddot ( n, vv(1,j), 1, vv(1,i1), 1 )
       hh(j,i) = t
       call daxpy ( n, -t, vv(1,j), 1, vv(1,i1), 1 )
     end do

     t = sqrt ( ddot ( n, vv(1,i1), 1, vv(1,i1), 1 ) )
     hh(i1,i) = t

     if ( t /= 0.0D+00 ) then
       t = 1.0D+00 / t
       vv(1:n,i1) = vv(1:n,i1) * t
     end if
!
!  Update factorization of HH.
!
    if ( i == 1 ) then
      go to 121
    end if
!
!  Perform previous transformations on I-th column of H.
!
    do k = 2, i
       k1 = k-1
       t = hh(k1,i)
       hh(k1,i) = c(k1) * t + s(k1) * hh(k,i)
       hh(k,i) = -s(k1) * t + c(k1) * hh(k,i)
    end do

121 continue

    gam = sqrt ( hh(i,i)**2 + hh(i1,i)**2 )
!
!  If GAMMA is zero then any small value will do.
!  It will affect only residual estimate.
!
    if ( gam == 0.0D+00 ) then
      gam = epsmac
    end if
!
!  Get the next plane rotation.
!
    c(i) = hh(i,i) / gam
    s(i) = hh(i1,i) / gam
    rs(i1) = -s(i) * rs(i)
    rs(i) = c(i) * rs(i)
!
!  Determine residual norm and test for convergence.
!
    hh(i,i) = c(i) * hh(i,i) + s(i) * hh(i1,i)
    ro = abs ( rs(i1) )
131 format(1h ,2e14.4)

    if ( 0 < iout ) then
      write(iout, 199) its, ro
    end if

    if ( i < im .and. eps1 < ro ) then
      go to 4
    end if
!
!  Now compute solution.  First solve upper triangular system.
!
    rs(i) = rs(i) / hh(i,i)

    do ii = 2, i
      k = i - ii + 1
      k1 = k + 1
      t = rs(k)
      do j = k1, i
        t = t - hh(k,j) * rs(j)
      end do
      rs(k) = t / hh(k,k)
    end do
!
!  Form linear combination of V(*,i)'s to get solution.
!
    t = rs(1)
    rhs(1:n) = vv(1:n,1) * t

    do j = 2, i
      t = rs(j)
      rhs(1:n) = rhs(1:n) + t * vv(1:n,j)
    end do
!
!  Call preconditioner.
!
    call lusol0 ( n, rhs, rhs, alu, jlu, ju )

    sol(1:n) = sol(1:n) + rhs(1:n)
!
!  Restart outer loop when necessary.
!
    if ( ro <= eps1 ) then
      ierr = 0
      exit
    end if

    if ( maxits < its ) then
      ierr = 1
      exit
    end if
!
!  Else compute residual vector and continue.
!
    do j = 1, i
      jj = i1 - j + 1
      rs(jj-1) = -s(jj-1) * rs(jj)
      rs(jj) = c(jj-1) * rs(jj)
    end do

    do j = 1, i1
      t = rs(j)
      if ( j == 1 ) then
        t = t - 1.0D+00
      end if
      call daxpy ( n, t, vv(1,j), 1,  vv, 1 )
    end do

199 format(' its =', i4, ' res. norm =', G14.6)

  end do

  return
end
subroutine pltmt ( nrow, ncol, mode, ja, ia, title, key, type, job, iounit )

!*****************************************************************************80
!
!! PLTMT creates a 'pic' plot of a matrix.
!
!  Discussion:
!
!    This routine creates a PIC file for plotting the pattern of
!    a sparse matrix stored in general sparse format. It is not intended
!    to be a means of plotting large matrices (It is very inefficient).
!
!    It is however useful for small matrices and can be used for example
!    for inserting matrix plots in a text. The size of the plot can be
!    7in x 7in or 5 in x 5in .. There is also an option for writing a
!    3-line header in troff (see description of parameter job).
!    See SPARSKIT/UNSUPP/ for a version of this to produce a post-script
!    file.
!
! example of usage .
!
! In the fortran code:
!  a) read a Harwell/Boeing matrix
!          call readmt (.....)
!         iout = 13
!  b) generate pic file:
!          call  pltmt (nrow,ncol,mode,ja,ia,title,key,type,iout)
!         stop
!
! Then in a unix environment plot the matrix by the command
!
!      pic FOR013.DAT | troff -me | lpr -Ppsx
!
!      1) Plots square as well as rectangular matrices.
!            (however not as much tested with rectangular matrices.)
!        2) the dot-size is adapted according to the size of the
!            matrix.
!        3) This is not meant at all as a way of plotting large
!            matrices. The pic file generaled will have one line for
!            each nonzero element. It is  only meant for use in
!           such things as document poreparations etc..
!         4) The caption written will print the 71 character long
!            title. This may not be centered correctly if the
!            title has trailing blanks (a problem with Troff).
!            if you want the title centered then you can center
!            the string in title before calling pltmt.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NCOL, the column dimension of the matrix.
!
!    Input, integer MODE, indicates matrix storage mode:
!    0, by rows,
!    1, by columns.
!
! ja     = column indices of nonzero elements when matrix is
!         stored rowise. Row indices if stores column-wise.
!
! ia     = integer array of containing the pointers to the
!         beginning of the columns in arrays a, ja.
!
! title  = character*71 = title of matrix test ( character a*71 ).
! key    = character*8  = key of matrix
! type   = character*3  = type of matrix.
!
! job    = this integer parameter allows to set a few minor
!          options. First it tells pltmt whether or not to
!          reduce the plot. The standard size of 7in is then
!          replaced by a 5in plot. It also tells pltmt whether or
!          not to append to the pic file a few 'troff' lines that
!          produce a centered caption includingg the title, key and
!          types as well as the size and number of nonzero elements.
!          job = 0 : do not reduce and do not make caption.
!          job = 1 : reduce and do not make caption.
!          job = 10 : do not reduce and make caption
!          job = 11 : reduce and make caption.
!          (i.e. trailing digit for reduction, leading digit for caption)
!
! iounit = logical unit number where to write the matrix into.
!
  implicit none

  integer ncol

  real ( kind = 8 ) hscale
  integer ia(ncol+1)
  integer ii
  integer ilast
  integer iounit
  integer ips
  integer istart
  integer ja(*)
  integer job
  integer k
  character ( len = 8 ) key
  integer maxdim
  integer mode
  integer n
  integer nnz
  integer nrow
  real ( kind = 8 ) ptsize
  real ( kind = 8 ) tiny
  character ( len = 72 ) title
  character ( len = 3 ) type
  real ( kind = 8 ) vscale
  real ( kind = 8 ) x
  real ( kind = 8 ) xht
  real ( kind = 8 ) xncol
  real ( kind = 8 ) xnrow
  real ( kind = 8 ) xshift
  real ( kind = 8 ) xwid
  real ( kind = 8 ) y
  real ( kind = 8 ) yshift

  n = ncol
  if ( mode == 0 ) then
    n = nrow
  end if

  nnz = ia(n+1) - ia(1)
  maxdim = max ( nrow, ncol )
  xnrow = real ( nrow, kind = 8 )
  xncol = real ( ncol, kind = 8 )
  ptsize = 0.08D+00
  hscale = ( 7.0D+00 - 2.0D+00 * ptsize ) / real ( maxdim - 1, kind = 8 )
  vscale = hscale
  xwid = ptsize + real ( ncol - 1, kind = 8 ) * hscale + ptsize
  xht = ptsize + real ( nrow - 1, kind = 8 ) * vscale + ptsize
  xshift = ( 7.0D+00 - xwid ) / 2.0D+00
  yshift = ( 7.0D+00 - xht ) / 2.0D+00

  if ( mod ( job, 10 ) == 1 ) then
    write (iounit,88)
  else
    write (iounit,89)
  end if

 88   format('.PS 5in',/,'.po 1.8i')
 89   format('.PS',/,'.po 0.7i')
  write(iounit,90)
 90   format('box invisible wid 7.0 ht 7.0 with .sw at (0.0,0.0) ')
  write(iounit,91) xwid, xht, xshift, yshift
 91   format('box wid ',f5.2,' ht ',f5.2, &
       ' with .sw at (',f5.2,',',f5.2,')' )
!
!  Shift points slightly to account for size of dot.
!
  tiny = 0.03D+00
  if ( mod ( job, 10 ) == 1 ) then
    tiny = 0.05D+00
  end if

  xshift = xshift + ptsize - tiny
  yshift = yshift + ptsize + tiny

  ips = 8
  if ( maxdim <= 500 ) then
    ips = 10
  end if
  if ( maxdim <= 300 ) then
    ips = 12
  end if
  if ( maxdim <= 100 ) then
    ips = 16
  end if
  if ( maxdim < 50 ) then
    ips = 24
  end if

  write(iounit,92) ips
 92   format('.ps ',i2)
!
!  Plotting loop
!
  do ii = 1, n

     istart = ia(ii)
     ilast = ia(ii+1)-1

     if ( mode /= 0 ) then
        x = real ( ii - 1, kind = 8 )
        do k = istart, ilast
           y = xnrow - real ( ja(k), kind = 8 )
           write(iounit,128) xshift+x*hscale, yshift+y*vscale
        end do
     else
        y = xnrow - real ( ii, kind = 8 )
        do k = istart, ilast
           x = real ( ja(k) - 1, kind = 8 )
           write(iounit,128) xshift+x*hscale, yshift+y*vscale
        end do
     end if

  end do

 128  format(7h"." at ,f6.3,',',f6.3,8h ljust  )
  write (iounit, 129)
 129  format('.PE')
!
!  Quit if caption not desired.
!
  if ( job / 10 /= 1 ) then
    return
  end if

  write(iounit,127) key, type, title
  write(iounit,130) nrow,ncol,nnz
 127  format('.sp 4'/'.ll 7i'/'.ps 12'/'.po 0.7i'/'.ce 3'/, &
       'Matrix:  ',a8,',  Type:  ',a3,/,a72)
 130  format('Dimension: ',i4,' x ',i4,',  Nonzero elements: ',i5)
  return
end
subroutine pltmtps ( nrow, ncol, mode, ja, ia, title, key, type, job, iounit )

!*****************************************************************************80
!
!! PLTMTPS creates a PostScript plot of a sparse matrix.
!
!  Discussion:
!
!    This routine creates a 'PS' file for plotting the pattern of
!    a sparse matrix stored in general sparse format. It can be used
!    for inserting matrix plots in a text. The size of the plot can be
!    7in x 7in or 5 in x 5in ..
!
!    1) Plots square as well as rectangular matrices.
!    2) Does not writer a caption yet.
!    3) No bounding box put in yet
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Paul Frederickson
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NCOL, the column dimension of the matrix.
!
!    Input, integer MODE, indicates the matrix storage mode:
!    0, by rows;
!    1, by columns.
!
! ja     = column indices of nonzero elements when matrix is
!         stored rowise. Row indices if stores column-wise.
! ia     = integer array of containing the pointers to the
!         beginning of the columns in arrays a, ja.
!
! title  = character*72 = title of matrix test ( character a*72 ).
! key    = character*8  = key of matrix
! type   = character*3  = type of matrix.
!
! job, integer. tells pltmt whether or not to reduce the plot.
!           if enabled then the standard size of 7in will be
!           replaced by a 5in plot.
!          job = 0 : do not reduce
!          job = 1 : reduce plot to 5 inches.
!
! iounit = logical unit number where to write the matrix into.
!
  implicit none

  real ( kind = 8 ) delta
  integer ia(*)
  integer ii
  integer ilast
  integer iounit
  integer istart
  integer ja(*)
  integer job
  integer k
  character ( len = 8 ) key
  integer m
  integer maxdim
  integer mode
  integer n
  integer ncol
  integer nrow
  integer nnz
  character ( len = 72 ) title
  character ( len = 3 ) type

  if ( mode == 0 ) then
    n = nrow
  else
    n = ncol
  end if

  nnz = ia(n+1) - ia(1)
  maxdim = max ( nrow, ncol )
  m = 1 + maxdim
!
!  Keep this test as in old pltmt (for future changes).
!
   if ( mod ( job, 10 ) == 1 ) then
     delta = 72.0D+00 * 5.0D+00 / ( 2.0D+00 + maxdim )
   else
     delta = 72.0D+00 * 7.0D+00 / (2.0D+00 + maxdim )
   end if

  write(iounit,*)'%!PS'
  write(iounit,*)' gsave 50 50 translate'
  write(iounit,*) delta, delta, ' scale'
  write(iounit,*) ' 0.25 setlinewidth'

  if ( mod ( job, 10 ) == 1 ) then
    write (iounit,*) ' 23 55 translate'
  else
    write (iounit,*) ' 2 35 translate'
  end if

  write(iounit,*) ' newpath'
  write(iounit,*) 0,0,' moveto'
  write(iounit,*) m,0,' lineto'
  write(iounit,*) m,m,' lineto'
  write(iounit,*) 0,m,' lineto'
  write(iounit,*) ' closepath stroke'
  write(iounit,*) ' 1 1 translate'
  write(iounit,*) ' 0.5 setlinewidth'
  write(iounit,*) ' /p {moveto 0 -.25 rmoveto '
  write(iounit,*) '            0  .50 rlineto stroke} def'
!
!  Plotting loop
!
  do ii = 1, n

    istart = ia(ii)
    ilast  = ia(ii+1)-1

    if ( mode /= 0 ) then

      do k = istart, ilast
        write(iounit,*) ii-1, nrow-ja(k), ' p'
      end do

    else
!             y = xnrow - real ( ii, kind = 8 )
      do k = istart, ilast
!               x = real ( ja(k) - 1, kind = 8 )
        write(iounit,*) ja(k)-1, nrow-ii, ' p'
      end do

    end if

  end do

  write(iounit,*)' showpage grestore'
 130    format('Dimension: ',i4,' x ',i4',  Nonzero elements: ',i5)
  return
end
subroutine project ( n, m, u, v, w )

!*****************************************************************************80
!
!! PROJECT computes the matrix-vector product w = U * v.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
!    Input, integer M, the column dimension of the matrix.
!
!    Input, real ( kind = 8 ) U(N,M), the matrix.
!
!    Input, real ( kind = 8 ) V(M), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) W(N), the product U*V.
!
  implicit none

  integer m
  integer n

  real ( kind = 8 ) u(n,m)
  real ( kind = 8 ) v(m)
  real ( kind = 8 ) w(n)

  w(1:n) = matmul ( u(1:n,1:m), v(1:m) )

  return
end
subroutine prtmt ( nrow, ncol, a, ja, ia, rhs, guesol, title, key, type, &
  ifmt, job, iounit )

!*****************************************************************************80
!
!! PRTMT writes a matrix in Harwell-Boeing format into a file.
!
!  Discussion:
!
!    This routine assumes that the matrix is stored in CSC format
!    (Compressed Sparse Column format).
!    There is some limited functionality for right hand sides.
!
!    This code attempts to pack as many elements as possible per
!    80-character line.
!
!    This code attempts to avoid as much as possible to put
!    blanks in the formats that are written in the 4-line header
!    This is done for purely esthetical reasons since blanks
!    are ignored in format descriptors.
!
!    sparse formats for right hand sides and guesses not supported.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NCOL, the column dimension of the matrix.
!
!    Input, real A(*), integer JA(*), IA(NCOL+1), the matrix in CSC
!    Compressed Sparse Column format.
!
!    Input, real RHS(*), contains the right hand sides and optionally
!    the associated initial guesses and/or exact solutions
!    in this order.  See also GUESOL for details.   RHS will
!    be used only if 2 < JOB.  Only full storage for the right hand 
!    sides is supported.
!
! guesol = a 2-character string indicating whether an initial guess
!          (1-st character) and / or the exact solution (2-nd)
!          character) is provided with the right hand side.
!         if the first character of guesol is 'G' it means that an
!          an intial guess is provided for each right hand sides.
!          These are assumed to be appended to the right hand sides in
!          the array rhs.
!         if the second character of guesol is 'X' it means that an
!          exact solution is provided for each right hand side.
!          These are assumed to be appended to the right hand sides
!          and the initial guesses (if any) in the array rhs.
!
! title  = character*71 = title of matrix test ( character a*71 ).
! key    = character*8  = key of matrix
! type   = charatcer*3  = type of matrix.
!
! ifmt       = integer specifying the format chosen for the real values
!         to be output (i.e., for a, and for rhs-guess-sol if
!          applicable). the meaning of ifmt is as follows.
!        * if (ifmt < 100) then the E descriptor is used,
!           format Ed.m, in which the length (m) of the mantissa is
!           precisely the integer ifmt (and d = ifmt+6)
!        * if (ifmt > 100) then prtmt will use the
!           F- descriptor (format Fd.m) in which the length of the
!           mantissa (m) is the integer mod(ifmt,100) and the length
!           of the integer part is k = ifmt/100 (and d = k+m+2)
!          Thus  ifmt= 4   means  E10.4  +.xxxxD+ee    while
!                ifmt=104  means  F7.4   +x.xxxx
!                ifmt=205  means  F9.5   +xx.xxxxx
!          Note: formats for ja, and ia are internally computed.
!
! job       = integer to indicate whether matrix values and
!         a right hand side is available to be written
!          job = 1   write srtucture only, i.e., the arrays ja and ia.
!          job = 2   write matrix including values, i.e., a, ja, ia
!          job = 3   write matrix and one right hand side: a,ja,ia,rhs.
!         job = nrhs+2 write matrix and nrhs successive right hand sides
!         Note that there cannot be any right hand side if the matrix
!         has no values. Also the initial guess and exact solutions when
!          provided are for each right hand side. For example if nrhs=2
!          and guesol='GX' there are 6 vectors to write.
!
!
! iounit = logical unit number where to write the matrix into.
!
! on return:
!
! the matrix a, ja, ia will be written in output unit iounit
! in the Harwell-Boeing format. Noe of the inputs is modofied.
!
  implicit none

  integer ncol

  real ( kind = 8 ) a(*)
  character ( len = 2 ) guesol
  integer i
  integer ia(ncol+1)
  integer iend
  integer ifmt
  integer ihead
  integer indcrd
  character ( len = 16 ) indfmt
  integer iounit
  integer ix
  integer ja(*)
  integer job
  character ( len = 8 ) key
  integer len
  integer next
  integer nnz
  integer nperli
  integer nrhs
  integer nrow
  integer ptrcrd
  character ( len = 16 ) ptrfmt
  real ( kind = 8 ) rhs(*)
  integer rhscrd
  character ( len = 3 ) rhstyp
  character ( len = 72 ) title
  integer totcrd
  character ( len = 3 ) type
  integer valcrd
  character ( len = 20 ) valfmt
!
!  Compute pointer format.
!
  nnz = ia(ncol+1) - 1
  len = int ( log10 ( 0.1D+00 + real ( nnz + 1, kind = 8 ) ) ) + 1
  nperli = 80 / len
  ptrcrd = ncol / nperli + 1

  if ( 9 < len ) then
     assign 101 to ix
  else
     assign 100 to ix
  end if

  write (ptrfmt,ix) nperli,len
 100  format(1h(,i2,1HI,i1,1h) )
 101  format(1h(,i2,1HI,i2,1h) )
!
!  Compute the ROW index format.
!
  len = int ( log10 ( 0.1D+00 + real ( nrow, kind = 8 ) ) ) + 1
  nperli = min ( 80 / len, nnz )
  indcrd = ( nnz - 1 ) / nperli + 1
  write (indfmt,100) nperli,len
!
!  Compute values and RHS format (using the same for both).
!
  valcrd = 0
  rhscrd = 0
!
!  Skip this part if no values provided.
!
  if ( job <= 1 ) then
    go to 20
  end if

  if ( 100 <= ifmt ) then
     ihead = ifmt / 100
     ifmt = ifmt - 100 * ihead
     len = ihead + ifmt + 2
     nperli = 80 / len

     if ( len <= 9 ) then
        assign 102 to ix
     elseif ( ifmt <= 9 ) then
        assign 103 to ix
     else
        assign 104 to ix
     end if

     write(valfmt,ix) nperli,len,ifmt
 102     format(1h(,i2,1hF,i1,1h.,i1,1h) )
 103     format(1h(,i2,1hF,i2,1h.,i1,1h) )
 104     format(1h(,i2,1hF,i2,1h.,i2,1h) )

  else
     len = ifmt + 6
     nperli = 80 / len
!
!  Try to minimize the blanks in the format strings.
!
     if ( nperli <= 9 ) then
      if ( len <= 9 ) then
         assign 105 to ix
      else if ( ifmt <= 9 ) then
         assign 106 to ix
      else
         assign 107 to ix
      end if
   else
      if ( len <= 9 ) then
         assign 108 to ix
      else if ( ifmt <= 9 ) then
         assign 109 to ix
      else
           assign 110 to ix
        end if
     end if

     write(valfmt,ix) nperli,len,ifmt
 105     format(1h(,i1,1hE,i1,1h.,i1,1h) )
 106     format(1h(,i1,1hE,i2,1h.,i1,1h) )
 107     format(1h(,i1,1hE,i2,1h.,i2,1h) )
 108     format(1h(,i2,1hE,i1,1h.,i1,1h) )
 109     format(1h(,i2,1hE,i2,1h.,i1,1h) )
 110     format(1h(,i2,1hE,i2,1h.,i2,1h) )

  end if
  valcrd = ( nnz - 1 ) / nperli + 1
  nrhs = job - 2

  if ( 1 <= nrhs ) then
     i = ( nrhs * nrow - 1 ) / nperli + 1
     rhscrd = i
     if ( guesol(1:1) == 'G' ) then
       rhscrd = rhscrd + i
     end if
     if ( guesol(2:2) == 'X' ) then
       rhscrd = rhscrd + i
     end if
     rhstyp = 'F' // guesol
  end if

 20   continue

  totcrd = ptrcrd + indcrd + valcrd + rhscrd
!
!  Write four line or five line header.
!
  write(iounit,10) title,key,totcrd,ptrcrd,indcrd,valcrd, &
       rhscrd,type,nrow,ncol,nnz,nrhs,ptrfmt,indfmt,valfmt,valfmt

  if ( 1 <= nrhs ) then
    write (iounit,11) rhstyp, nrhs
  end if

 10   format (a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
 11   format(A3,11x,i4)

  write(iounit,ptrfmt) ia(1:ncol+1)
  write(iounit,indfmt) ja(1:nnz)

  if ( job <= 1 ) then
    return
  end if

  write(iounit,valfmt) (a(i), i = 1, nnz)
  if ( job <= 2 ) then
    return
  end if
  len = nrow * nrhs
  next = 1
  iend = len
  write(iounit,valfmt) (rhs(i), i = next, iend)
!
!  Write initial guesses if available
!
  if ( guesol(1:1) == 'G' ) then
     next = next + len
     iend = iend + len
     write(iounit,valfmt) (rhs(i), i = next, iend)
  end if
!
!  Write exact solutions if available.
!
  if ( guesol(2:2) == 'X' ) then
     next = next + len
     iend = iend + len
     write(iounit,valfmt) (rhs(i), i = next, iend)
  end if

  return
end
subroutine readmt ( nmax, nzmax, job, iounit, a, ja, ia, rhs, nrhs, &
  guesol, nrow, ncol, nnz, title, key, type, ierr )

!*****************************************************************************80
!
!! READMT reads a Harwell/Boeing sparse matrix file. 
!
!  Discussion:
!
!    The routine handles right hand sides in full format only (no 
!    sparse right hand sides).
!
!    The file inout must be open (and possibly rewound if necessary)
!    prior to calling readmt.
!
!    Refer to the documentation on the Harwell-Boeing formats
!    for details on the format assumed by readmt.
!    We summarize the format here for convenience.
!
!    a) all lines in inout are assumed to be 80 character long.
!    b) the file consists of a header followed by the block of the
!       column start pointers followed by the block of the
!       row indices, followed by the block of the real values and
!       finally the numerical values of the right hand side if a
!       right hand side is supplied.
!    c) the file starts by a header which contains four lines if no
!       right hand side is supplied and five lines otherwise.
!       * first line contains the title (72 characters long) followed by
!         the 8-character identifier (name of the matrix, called key)
!        [ A72,A8 ]
!       * second line contains the number of lines for each
!         of the following data blocks (4 of them) and the total number
!         of lines excluding the header.
!        [5i4]
!       * the third line contains a three character string identifying
!         the type of matrices as they are referenced in the Harwell
!         Boeing documentation [e.g., rua, rsa,..] and the number of
!         rows, columns, nonzero entries.
!         [A3,11X,4I14]
!       * The fourth line contains the variable fortran format
!         for the following data blocks.
!         [2A16,2A20]
!       * The fifth line is only present if right hand sides are
!         supplied. It consists of three one character-strings containing
!         the storage format for the right hand sides
!         ('F'= full,'M'=sparse=same as matrix), an initial guess
!         indicator ('G' for yes), an exact solution indicator
!         ('X' for yes), followed by the number of right hand sides
!         and then the number of row indices.
!         [A3,11X,2I14]
!     d) The three following blocks follow the header as described
!        above.
!     e) In case the right hand side are in sparse formats then
!        the fourth block uses the same storage format as for the matrix
!        to describe the NRHS right hand sides provided, with a column
!        being replaced by a right hand side.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
! nmax, max column dimension  allowed for matrix. The array ia should
!          be of length at least ncol+1 (see below) if job>0
! nzmax       = max number of nonzeros elements allowed. the arrays a,
!          and ja should be of length equal to nnz (see below) if these
!          arrays are to be read (see job).
!
! job       = integer to indicate what is to be read. (note: job is an
!          input and output parameter, it can be modified on return)
!          job = 0    read the values of ncol, nrow, nnz title, key,
!                     type and return. matrix is not read and arrays
!                     a, ja, ia, rhs are not touched.
!          job = 1    read srtucture only, i.e., the arrays ja and ia.
!          job = 2    read matrix including values, i.e., a, ja, ia
!          job = 3    read matrix and right hand sides: a,ja,ia,rhs.
!                  rhs may contain initial guesses and exact
!                     solutions appended to the actual right hand sides.
!                  this will be indicated by the output parameter
!                     guesol [see below].
!
! nrhs   = integer. nrhs is an input as well as ouput parameter.
!          at input nrhs contains the total length of the array rhs.
!          See also ierr and nrhs in output parameters.
!
! iounit = logical unit number where to read the matrix from.
!
! on return:
!
! job    = on return job may be modified to the highest job it could
!          do: if job=2 on entry but no matrix values are available it
!          is reset to job=1 on return. Similarly of job=3 but no rhs
!          is provided then it is rest to job=2 or job=1 depending on
!          whether or not matrix values are provided.
!          Note that no error message is triggered (i.e. ierr = 0
!          on return in these cases. It is therefore important to
!          compare the values of job on entry and return ).
!
!    Output, real A(*), JA(*), IA(NCOL+1), the matrix in CSC
!    Compressed Sparse Column format.
!
! rhs    = real array of size nrow + 1 if available (see job)
!
! nrhs   = integer containing the number of right hand sides found
!          each right hand side may be accompanied with an intial guess
!          and also the exact solution.
!
! guesol = a 2-character string indicating whether an initial guess
!          (1-st character) and / or the exact solution (2-nd
!          character) is provided with the right hand side.
!         if the first character of guesol is 'G' it means that an
!          an intial guess is provided for each right hand side.
!          These are appended to the right hand sides in the array rhs.
!         if the second character of guesol is 'X' it means that an
!          exact solution is provided for each right hand side.
!          These are  appended to the right hand sides
!          and the initial guesses (if any) in the array rhs.
!
!    Output, integer NROW, the row dimension of the matrix.
!
!    Output, integer NCOL, the column dimension of the matrix.
!
! nnz       = number of nonzero elements in A. This info is returned
!          even if there is not enough space in a, ja, ia, in order
!          to determine the minimum storage needed.
!
! title  = character*72 = title of matrix test ( character a*72).
! key    = character*8  = key of matrix
! type   = charatcer*3  = type of matrix.
!          for meaning of title, key and type refer to documentation
!          Harwell/Boeing matrices.
!
! ierr   = integer used for error messages
!         * ierr  = 0 means that  the matrix has been read normally.
!         * ierr  = 1 means that  the array matrix could not be read
!         because ncol+1 > nmax
!         * ierr  = 2 means that  the array matrix could not be read
!         because nnz > nzmax
!         * ierr  = 3 means that  the array matrix could not be read
!         because both (ncol+1 > nmax) and  (nnz > nzmax )
!         * ierr  = 4 means that  the right hand side (s) initial
!         guesse (s) and exact solution (s)   could  not be
!         read because they are stored in sparse format (not handled
!         by this routine ...)
!         * ierr  = 5 means that the right hand sides, initial guesses
!         and exact solutions could not be read because the length of
!         rhs as specified by the input value of nrhs is not
!         insufficient to store them. The rest of the matrix may have
!         been read normally.
!
  implicit none

  integer nmax
  integer nzmax

  real ( kind = 8 ) a(nzmax)
  character ( len = 2 ) guesol
  integer ia(nmax+1)
  integer iend
  integer ierr
  integer indcrd
  character ( len = 16 ) indfmt
  integer iounit
  integer ja(nzmax)
  integer job
  character ( len = 8 ) key
  integer len
  integer lenrhs
  integer n
  integer ncol
  integer neltvl
  integer next
  integer nnz
  integer nrhs
  integer nrow
  integer nvec
  integer ptrcrd
  character ( len = 16 ) ptrfmt
  real ( kind = 8 ) rhs(*)
  integer rhscrd
  character ( len = 20 ) rhsfmt
  character ( len = 3 ) rhstyp
  character ( len = 72 ) title
  integer totcrd
  character ( len = 3 ) type
  integer valcrd
  character ( len = 20 ) valfmt

  lenrhs = nrhs

  read(iounit,2010,end=10)title,key
  read(iounit,2011,end=10)totcrd,ptrcrd,indcrd,valcrd,rhscrd
  read(iounit,2012,end=10)type,nrow,ncol,nnz,neltvl
  read(iounit,2013,end=10)ptrfmt,indfmt,valfmt,rhsfmt
2010  format(a72,a8)
2011  format(5i14)
2012  format(a3,11x,4i14)
2013  format(2a16, 2a20)

  if ( 0 < rhscrd ) then
    read (iounit,2014,end=10) rhstyp, nrhs
  end if

2014  format (a3,11x,i4)
!
!  Anything else to read?
!
  if ( job <= 0 ) then
    return
  end if

  ierr = 0
!
!  Check whether matrix is readable.
!
  n = ncol
  if ( nmax < ncol ) then
    ierr = 1
  end if

  if ( nzmax < nnz ) then
    ierr = ierr + 2
  end if

  if ( ierr /= 0 ) then
    return
  end if
!
!  Read pointer and row numbers.
!
  read (iounit,ptrfmt,end=10) ia(1:n+1)
  read (iounit,indfmt,end=10) ja(1:nnz)
!
!  Reading values of matrix if required...
!
  if ( job <= 1 ) then
    return
  end if
!
!  ...and if available.
!
  if ( valcrd <= 0 ) then
    job = 1
    return
  end if

  read (iounit,valfmt,end=10) a(1:nnz)
!
!  Reading RHS if required...
!
  if ( job <= 2 ) then
    return
  end if
!
!  ...and if available.
!
  if ( rhscrd <= 0 ) then
    job = 2
    return
  end if
!
!  Read right hand side.
!
  if ( rhstyp(1:1) == 'M' ) then
     ierr = 4
     return
  end if

  guesol = rhstyp(2:3)

  nvec = 1
  if ( guesol(1:1) == 'G' ) then
    nvec=nvec+1
  end if
  if ( guesol(2:2) == 'X' ) then
    nvec=nvec+1
  end if

  len = nrhs * nrow

  if ( lenrhs < len * nvec ) then
     ierr = 5
     return
  end if
!
!  Read right hand sides.
!
  next = 1
  iend = len
  read(iounit,rhsfmt,end=10) rhs(next:iend)
!
!  Read initial guesses if available.
!
  if ( guesol(1:1) == 'G' ) then
     next = next + len
     iend = iend + len
     read(iounit,valfmt,end=10) rhs(next:iend)
  end if
!
!  Read exact solutions if available.
!
  if ( guesol(2:2) == 'X' ) then
     next = next + len
     iend = iend + len
     read(iounit,valfmt,end=10) rhs(next:iend)
  end if

  return
10    continue
  WRITE(*,*)' '
  WRITE(*,*)'READMT - Fatal error.'
  WRITE(*,*)'         End of file while reading information!'
  WRITE(*,*)'         Results are unreliable!'

  return
end
subroutine refall ( nx, nelx, ijk, node, ndeg, x, y, ichild, iparnts, &
  nodcode, nxmax, nelmax, ierr )

!*****************************************************************************80
!
!! REFALL refines a finite element grid using triangular elements.
!
!  Discussion:
!
!    The routine uses midpoints to refine all the elements of the grid.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NX, the number of nodes at input.
!
!    Input, integer NELX, the number of elements.
!
!    Input, integer IJK(NODE,NELX), lists the nodes that make up 
!    each element.
!
!    Input, integer NODE, the number of nodes per element.
!
! ndeg      = first dimension of array ichild which is at least as large
!         as the max degree of each node
!
! x,y   = real arrays containing the x(*) and y(*) coordinates
!        resp. of the nodes.
! ichild= list of the children of a node: ichild(1,k) stores
!         the position in ichild(*,k)  of the last child so far.
!         (local use)
! iparnts= list of the 2 parents of each node.
!         (local use)
! nodcode= boundary information list for each node with the
!         following meaning:
!      nodcode(i) = 0 -->  node i is internal
!      nodcode(i) = 1 -->  node i is a boundary but not a corner point
!      nodcode(i) = 2 -->  node i is a corner point.
! corner elements are used only to generate the grid by refinement
! since they do not  correspond to real elements.
! nxmax  = maximum number of nodes allowed. If during the algorithm
!          the number of nodes being created exceeds nxmax then
!         refall  quits without modifying the (x,y) xoordinates
!         and nx, nelx. ijk is modified. Also ierr is set to 1.
! nelmax = same as above for number of elements allowed. See ierr..
! ierr       = error message:
!         0 --> normal return
!         1 --> refall quit because nxmax  was exceeded.
!         2 --> refall quit because nelmax was exceeded.
!
  implicit none

  integer ndeg
  integer node
  integer nx

  integer i
  integer ichild(ndeg,*)
  integer ierr
  integer ii
  integer ipar1
  integer ipar2
  integer iparnts(2,nx)
  integer ijk(node,*)
  integer jchild
  integer jj
  integer jnod
  integer k
  integer k1
  integer k2
  integer last
  integer nodcode(nx)
  integer midnode(10)
  integer inod(10)
  integer nel
  integer nelmax
  integer nelx
  integer nelxnew
  integer nxmax
  integer nxnew
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)
!
!  Initialize the lists of children and parents.
!
!  The data structure is as follows:
!   ICHILD(1,K) stores the position of last child of node K so far in list
!   ICHILD(J,K), J >= 2 = list of children of node K.
!   IPARNTS(1,K) and IPARNTS(2,K) are the two parents of node K.
!
!  First check:
!
  if ( nxmax <= nx ) then
    ierr = 1
    return
  end if

  if ( nelmax <= nelx ) then
    ierr = 2
    return
  end if
!
!  Initialize.
!
  iparnts(1:2,1:nx)= 0

  ichild(1,1:nx) = 1
  ichild(2:ndeg,1:nx) = 0

  nelxnew = nelx
  nxnew = nx
  ierr = 0
!
!  The main loop: scan all elements.
!
  do nel = 1, nelx
!
!  Interesting question which order is best for parallelism?
!  alternative order: do nel = nelx, 1, -1
!
!  Unpack nodes of element.
!
    do i = 1, node

      inod(i) = ijk(i,nel)
!
!  Convention: node after last node is first node.
!
      inod(node+i) = inod(i)
      midnode(i) = 0

    end do
!
!  For each new potential node determine if it has already been
!  numbered.  A potential node is the middle of any two nodes.
!
    do 80 ii = 1, node

          k1 = inod(ii)
          k2 = inod(ii+1)
!
!  Test for current pair.
!
    last = ichild(1,k1)

    do k = 2, last

       jchild = ichild(k,k1)
        ipar1 = iparnts(1,jchild)
        ipar2 = iparnts(2,jchild)
          if ( (ipar1 == k1 .and. ipar2 == k2) .or. &
                  (ipar2 == k1 .and. ipar1 == k2)) then
!
!  The node has already been created and numbered.
!
      midnode(ii) = jchild
!
!  Therefore it must be an internal node...
!
      nodcode(jchild) = 0
!
!  ...and there is no new node to create.
!
      go to 80
          end if

    end do
!
!  Else create a new node.
!
      nxnew = nxnew + 1

      if ( nxmax < nxnew ) then
        ierr = 1
        return
      end if

      x(nxnew) = (x(k1) + x(k2)) * 0.5D+00
      y(nxnew) = (y(k1) + y(k2)) * 0.5D+00
      midnode(ii) = nxnew
!
!  Update NODCODE information, normally min ( NODCODE(K1), NODCODE(K2) ).
!
   nodcode(nxnew) = min ( 1, nodcode(k1), nodcode(k2) )
!
!  Update parents and children's lists.
!
      iparnts(1,nxnew) = k1
      iparnts(2,nxnew) = k2

      last = last+1
      ichild(last,k1) = nxnew
      ichild(1,k1) = last

      last = ichild(1,k2)+1
      ichild(last,k2) = nxnew
      ichild(1,k2) = last

80  continue
!
!  Replace current element by new one.
!
    do i = 1, node
      jnod = midnode(i)
      ijk(i,nel) = jnod
    end do
!
!  Create new elements.
!
    do ii = 1, node

      nelxnew = nelxnew + 1

      if ( nelmax < nelxnew ) then
        ierr = 2
        return
      end if

      ijk(1,nelxnew) = inod(ii)
      k = ii

      do jj = 2, node
        ijk(jj,nelxnew) = midnode(k)
        k = k + 2
        if ( node < k ) then
          k = k - node
        end if
      end do

    end do
!
!  Done!
!
  end do

  nx = nxnew
  nelx = nelxnew

  return
end
subroutine retmx ( n, a, ja, ia, dd )

!*****************************************************************************80
!
!! RETMX returns in dd(*) the max absolute value of elements in row *.
!
!  Discussion:
!
!    This routine is used for scaling.  It has been superseded by RNRMS.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real A(*), JA(*), IA(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real DD(N), the element of each row that has the largest absolute
!    value.  The sign of DD is modified such that it is the same as that 
!    of the diagonal element in its row.
!
  implicit none

  integer n

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) dd(*)
  integer i
  integer ia(n+1)
  integer ja(*)
  integer k
  integer k1
  integer k2
  real ( kind = 8 ) t
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
!
!  Initialize.
!
  k2 = 1

  do i = 1, n

     k1 = k2
     k2 = ia(i+1) - 1
     t = 0.0D+00

     do k = k1, k2

        t1 = abs ( a(k) )
        if ( t < t1 ) then
          t = t1
        end if

        if ( ja(k) == i ) then

          if ( a(k) < 0.0D+00 ) then
            t2 = -1.0D+00
          else if ( a(k) == 0.0D+00 ) then
            t2 = 0.0D+00
          else
            t2 = 1.0D+00
            end if

          end if

     end do

     dd(i) = t2 * t
!
!  We do not invert the diagonal entries here.
!
  end do

  return
end
subroutine rnrms ( nrow, nrm, a, ja, ia, diag )

!*****************************************************************************80
!
!! RNRMS gets the norms of each row of A.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
! nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
!                  means the 2-nrm, nrm = 0 means max norm
!
!    Input, real A(*), integer, JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real DIAG(NROW), the row norms.
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) diag(nrow)
  integer ia(nrow+1)
  integer ii
  integer ja(*)
  integer k
  integer k1
  integer k2
  integer nrm
  real ( kind = 8 ) scal
!
!  Compute the norm of each element.
!
  do ii = 1, nrow

    scal = 0.0D+00
    k1 = ia(ii)
    k2 = ia(ii+1) - 1

    if ( nrm == 0 ) then
      do k = k1, k2
        scal = max ( scal, abs ( a(k) ) )
      end do
    else if ( nrm == 1 ) then
      do k = k1, k2
        scal = scal + abs ( a(k) )
      end do
    else
      do k = k1, k2
        scal = scal + a(k)**2
      end do
    end if

    if ( nrm == 2 ) then
      scal = sqrt ( scal )
    end if

    diag(ii) = scal

  end do

  return
end
subroutine rperm ( nrow, a, ja, ia, ao, jao, iao, perm, job )

!*****************************************************************************80
!
!! RPERM permutes the rows of a matrix in CSR format.
!
!  Discussion:
!
!    This routine computes B = P*A  where P is a permutation matrix.
!    the permutation P is defined through the array perm: for each j,
!    perm(j) represents the destination row number of row number j.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! perm       = integer array of length nrow containing the permutation arrays
!        for the rows: perm(i) is the destination of row i in the
!         permuted matrix.
!         ---> a(i,j) in the original matrix becomes a(perm(i),j)
!         in the output  matrix.
!
! job      = integer indicating the work to be done:
!             job = 1      permute a, ja, ia into ao, jao, iao
!                       (including the copying of real values ao and
!                       the array iao).
!             job /= 1 :  ignore real values.
!                     (in which case arrays a and ao are not needed nor
!                      used).
!
!    Output, real AO(*), integer JAO(*), IAO(NROW+1), the permuted matrix
!    in CSR Compressed Sparse Row format.
!
! note :
!        if (job/=1)  then the arrays a and ao are not used.
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer i
  integer ia(nrow+1)
  integer iao(nrow+1)
  integer ii
  integer j
  integer ja(*)
  integer jao(*)
  integer job
  integer k
  integer ko
  integer perm(nrow)
  logical values

  values = ( job == 1 )
!
!  Determine pointers for output matrix.
!
  do j = 1, nrow
    i = perm(j)
    iao(i+1) = ia(j+1) - ia(j)
  end do
!
!  Get pointers from lengths.
!
  iao(1) = 1
  do j = 1, nrow
    iao(j+1) = iao(j+1) + iao(j)
  end do
!
!  Copying.
!
  do ii = 1, nrow
!
!  Old row = II is new row IPERM(II), and KO is the new pointer.
!
     ko = iao(perm(ii))

     do k = ia(ii), ia(ii+1)-1
       jao(ko) = ja(k)
       if ( values ) then
         ao(ko) = a(k)
       end if
       ko = ko + 1
     end do

  end do

  return
end
subroutine rscal ( nrow, job, nrm, a, ja, ia, diag, b, jb, ib )

!*****************************************************************************80
!
!! RSCAL normalizes the rows of A.
!
!  Discussion:
!
!    There are three choices for the norm to use, the 1-norm, 2-norm 
!    or max-norm.
!
!    The column dimension of A is not needed.
!
!    The algorithm is in-place, so the A and B information can share
!    the same memory.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
! job   = integer. job indicator. Job=0 means get array b only
!         job = 1 means get b, and the integer arrays ib, jb.
!
! nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
!                  means the 2-nrm, nrm = 0 means max norm
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! on return:
!
! diag = diagonal matrix stored as a vector containing the matrix
!        by which the rows have been scaled, i.e., on return
!        we have B = Diag*A.
!
!    Output, real B(*), integer JB(*), IB(NROW+1), the scaled matrix in CSR
!    Compressed Sparse Row format.
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) b(*)
  real ( kind = 8 ) diag(nrow)
  integer ia(nrow+1)
  integer ib(nrow+1)
  integer ja(*)
  integer jb(*)
  integer job
  integer nrm

  call rnrms ( nrow, nrm, a, ja, ia, diag )

  diag(1:nrow) = 1.0D+00 / diag(1:nrow)

  call diamua ( nrow, job, a, ja, ia, diag, b, jb, ib )

  return
end
subroutine sskssr ( n, imod, asky, isky, ao, jao, iao, nzmax, ierr )

!*****************************************************************************80
!
!! SSKSSR converts Symmetric Skyline Format to Symmetric Sparse Row format.
!
!  Discussion:
!
!    This routine translates a symmetric skyline format into a
!    symmetric sparse row format.  Each element is tested to see if it is
!    a zero element.  Only the actual nonzero elements are retained. Note
!    that the test used is simple and does take into account the smallness
!    of a value.  The routine FILTER can be used
!    for this purpose.
!
!    This routine is an in-place algorithm, so ASKY and ISKY can be
!    the same as AO and IAO.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
! imod  = integer indicating the variant of skyline format used:
!         imod = 0 means the pointer iao points to the `zeroth'
!         element of the row, i.e., to the position of the diagonal
!         element of previous row (for i = 1, iao(1)= 0)
!         imod = 1 means that itpr points to the beginning of the row.
!         imod = 2 means that iao points to the end of the row
!                  (diagonal element)
! asky  = real array containing the values of the matrix. asky contains
!         the sequence of active rows from i = 1, to n, an active row
!         being the row of elemnts of the matrix contained between the
!         leftmost nonzero element and the diagonal element.
!
! isky       = integer array of size n+1 containing the pointer array to
!         each row. isky (k) contains the address of the beginning of the
!         k-th active row in the array asky.
!
! nzmax = integer. equal to the number of available locations in the
!         output array ao.
!
! on return:
!
!    Output, real AO(*), integer JAO(*), IAO(N+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! ierr  = integer. Serving as error message. If the length of the
!         output arrays ao, jao exceeds nzmax then ierr returns
!         the row number where the algorithm stopped: rows
!         i, to ierr-1 have been processed succesfully.
!         ierr = 0 means normal return.
!         ierr = -1  : illegal value for imod
!
  implicit none

  integer n
  integer nzmax

  real ( kind = 8 ) ao(nzmax)
  real ( kind = 8 ) asky(*)
  integer i
  integer iao(n+1)
  integer ierr
  integer imod
  integer isky(n+1)
  integer j
  integer jao(nzmax)
  integer k
  integer kend
  integer kstart
  integer next

  ierr = 0
!
!  Check for validity of IMOD.
!
  if ( imod /= 0 .and. imod /= 1 .and. imod /= 2 ) then
     ierr =-1
     return
  end if
!
!  NEXT = pointer to next available position in output matrix
!  KEND = pointer to end of current row in skyline matrix.
!
  next = 1
!
!  Set KEND = start position -1 in skyline matrix.
!
  kend = 0
  if ( imod == 1 ) then
    kend = isky(1)-1
  end if
  if ( imod == 0 ) then
    kend = isky(1)
  end if
!
!  Loop through all rows
!
  do i = 1, n
!
!  Save value of pointer to I-th row in output matrix.
!
     iao(i) = next
!
!  Get beginnning and end of skyline row.
!
     kstart = kend + 1

     if ( imod == 0 ) then
       kend = isky(i+1)
     else if ( imod == 1 ) then
       kend = isky(i+1)-1
     else if ( imod == 2 ) then
       kend = isky(i)
     end if
!
!  Copy element into output matrix unless it is a zero element.
!
     do k = kstart, kend
       if ( asky(k) /= 0.0D+00 ) then
         j = i - ( kend - k )
         jao(next) = j
         ao(next) = asky(k)
         next = next + 1
         if ( nzmax+1 < next ) then
           ierr = i
           return
         end if
       end if
    end do

  end do

  iao(n+1) = next

  return
end
subroutine ssrcsr ( nrow, a, ja, ia, nzmax, ao, jao, iao, indu, ierr )

!*****************************************************************************80
!
!! SSRCSR converts Symmetric Sparse Row to (regular) Compressed Sparse Row.
!
!  Discussion:
!
!    This routine converts a symmetric  matrix in which only the lower
!    part is  stored in compressed sparse row format, i.e.,
!    a matrix stored in symmetric sparse format, into a fully stored matrix
!    i.e., a matrix where both the lower and upper parts are stored in
!    compressed sparse row format. the algorithm is in place (i.e. result
!    may be overwritten onto the input matrix a, ja, ia ----- ).
!
!    the output matrix delivered by ssrcsr is such that each row starts with
!    the elements of the lower part followed by those of the upper part.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
! a,
! ia,
! ja    = matrix in compressed sparse row format. This is assumed to be
!         a lower triangular matrix.
!
! nzmax      = size of arrays ao and jao. ssrcsr will abort if the storage
!         provided in a, ja is not sufficient to store A. See ierr.
!
! on return:
!
! ao, iao,
!   jao = output matrix in compressed sparse row format. The resulting
!         matrix is symmetric and is equal to A+A**T - D, if
!         A is the original lower triangular matrix. ao, jao, iao,
!         can be the same as a, ja, ia in the calling sequence.
!
! indu  = integer array of length nrow+1. If the input matrix is such
!         that the last element in each row is its diagonal element then
!         on return, indu will contain the pointers to the diagonal
!         element in each row of the output matrix. Otherwise used as
!         work array.
! ierr  = integer. Serving as error message. If the length of the arrays
!         ao, jao exceeds nzmax, ierr returns the minimum value
!         needed for nzmax. otherwise ierr=0 (normal return).
!
  implicit none

  integer nrow
  integer nzmax

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(nzmax)
  integer i
  integer ia(nrow+1)
  integer iao(nrow+1)
  integer ierr
  integer indu(nrow+1)
  integer ipos
  integer j
  integer ja(*)
  integer jao(nzmax)
  integer k
  integer ko
  integer kfirst
  integer klast
  integer kosav
  integer lenrow
  integer nnz

  ierr = 0
  indu(1:nrow+1) = 0
!
!  Compute number of elements in each row of strict upper part.
!  Put result in INDU(I+1) for row I.
!
  do i = 1, nrow
    do k = ia(i), ia(i+1)-1
       j = ja(k)
       if ( j < i ) then
         indu(j+1) = indu(j+1) + 1
       end if
    end do
  end do
!
!  Find addresses of first elements of ouput matrix.  Result in INDU.
!
  indu(1) = 1
  do i = 1, nrow
    lenrow = ia(i+1)-ia(i)
    indu(i+1) = indu(i) + indu(i+1) + lenrow
  end do
!
!  Enough storage in A, JA?
!
  nnz = indu(nrow+1) - 1

  if ( nzmax < nnz ) then
    ierr = nnz
    return
  end if
!
!  Now copy lower part (backwards).
!
  kosav = indu(nrow+1)

  do i = nrow, 1, -1

     klast = ia(i+1) - 1
     kfirst = ia(i)
     iao(i+1) = kosav
     ko = indu(i)
     kosav = ko

     do k = kfirst, klast
       ao(ko) = a(k)
       jao(ko) = ja(k)
       ko = ko+1
     end do

     indu(i) = ko
  end do

  iao(1) = 1
!
!  Copy upper part.  Go through the structure of AO, JAO, IAO
!  that has already been copied (lower part).  INDU(I) is the address
!  of the next free location in row I for AO, JAO.
!
!  I-th row is now in AO, JAO, IAO structure, lower half part.
!
  do i = 1, nrow

    do k = iao(i), iao(i+1)-1

      j = jao(k)

      if ( i <= j ) then
        exit
      end if

      ipos = indu(j)
      ao(ipos) = ao(k)
      jao(ipos) = i
      indu(j) = indu(j) + 1

    end do

  end do

  return
end
subroutine submat ( n, job, i1, i2, j1, j2, a, ja, ia, nr, nc, ao, jao, iao )

!*****************************************************************************80
!
!! SUBMAT extracts the submatrix A(i1:i2,j1:j2).
!
!  Discussion:
!
!    This routine extracts a submatrix and puts the result in
!    matrix ao,iao,jao.  It is an "in place" routine, so ao,jao,iao may be 
!    the same as a,ja,ia.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the row dimension of the matrix.
!
! i1,i2 = two integers with i2 >= i1 indicating the range of rows to be
!          extracted.
!
! j1,j2 = two integers with j2 >= j1 indicating the range of columns
!         to be extracted.
!         * There is no checking whether the input values for i1, i2, j1,
!           j2 are between 1 and n.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
! job      = job indicator: if job /= 1 then the real values in a are NOT
!         extracted, only the column indices (i.e. data structure) are.
!         otherwise values as well as column indices are extracted...
!
! on output
!
! nr      = number of rows of submatrix
! nc      = number of columns of submatrix
!        * if either of nr or nc is nonpositive the code will quit.
!
! ao,
! jao,iao = extracted matrix in general sparse format with jao containing
!      the column indices,and iao being the pointer to the beginning
!      of the row,in arrays a,ja.
!
  implicit none

  real ( kind = 8 ) a(*)
  real ( kind = 8 ) ao(*)
  integer i
  integer i1
  integer i2
  integer ia(*)
  integer iao(*)
  integer ii
  integer j
  integer j1
  integer j2
  integer ja(*)
  integer jao(*)
  integer job
  integer k
  integer k1
  integer k2
  integer klen
  integer n
  integer nc
  integer nr

  nr = i2 - i1 + 1
  nc = j2 - j1 + 1

  if ( nr <= 0 .or. nc <= 0 ) then
    return
  end if

  klen = 0
!
!  Simple procedure that proceeds row-wise.
!
  do i = 1, nr

    ii = i1 + i - 1
    k1 = ia(ii)
    k2 = ia(ii+1) - 1
    iao(i) = klen + 1

    do k = k1, k2
      j = ja(k)
      if ( j1 <= j .and. j <= j2 ) then
        klen = klen + 1
        if ( job == 1 ) then
          ao(klen) = a(k)
        end if
        jao(klen) =j
      end if
    end do

  end do

  iao(nr+1) = klen + 1

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    15 March 2003
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

  character ( len = 40 ) string

  call timestring ( string )

  write ( *, '(a)' ) trim ( string )

  return
end
subroutine timestring ( string )

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = 'May 31 2001   9:45:54.872 AM'
!
!  Modified:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = * ) string
  character ( len = 10 ) time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( string, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine transp ( nrow, ncol, a, ja, ia, iwk, ierr )

!*****************************************************************************80
!
!! TRANSP carries out in-place transposition routine.
!
!  Discussion:
!
!    This routine transposes a matrix stored in compressed sparse row
!    format.  The transposition is done in place in that the arrays 
!    A, JA, and IA of the transpose are overwritten onto the original arrays.
!
!    If you do not need the transposition to be done in place
!    it is preferrable to use the conversion routine csrcsc
!    (see conversion routines in formats).
!
!    The entries of the output matrix are not sorted (the column
!    indices in each are not in increasing order).  Use CSRCSC
!    if you want them sorted.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NCOL, the column dimension of the matrix.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Workspace, integer IWK(*), of the same length as JA.
!
! on return:
!
!
! ncol      = actual row dimension of the transpose of the input matrix.
!         Note that this may be <= the input value for ncol, in
!         case some of the last columns of the input matrix are zero
!         columns. In the case where the actual number of rows found
!         in transp(A) exceeds the input value of ncol, transp will
!         return without completing the transposition. see ierr.
!
!    Input, real A(*), integer JA(*), IA(NCOL+1), the transposed matrix
!    in CSR Compressed Sparse Row format.
!
! ierr      = integer. error message. If the number of rows for the
!         transposed matrix exceeds the input value of ncol,
!         then ierr is  set to that number and transp quits.
!         Otherwise ierr is set to 0 (normal return).
!
  implicit none

  integer nrow

  real ( kind = 8 ) a(*)
  integer i
  integer ia(nrow+1)
  integer ierr
  integer inext
  integer init
  integer iwk(*)
  integer j
  integer ja(*)
  integer jcol
  integer k
  integer l
  integer ncol
  integer nnz
  real ( kind = 8 ) t
  real ( kind = 8 ) t1

  ierr = 0
  nnz = ia(nrow+1) - 1
!
!  Determine the column dimension.
!
  jcol = 0
  do k = 1, nnz
    jcol = max ( jcol, ja(k) )
  end do

  if ( ncol < jcol ) then
     ierr = jcol
     return
  end if
!
!  Convert to coordinate format.  Use IWK for row indices.
!
  ncol = jcol

  do i = 1, nrow
    do k = ia(i), ia(i+1)-1
      iwk(k) = i
    end do
  end do
!
!  Find pointer array for transpose.
!
  ia(1:ncol+1) = 0

  do k = 1, nnz
    i = ja(k)
    ia(i+1) = ia(i+1) + 1
  end do
  ia(1) = 1

  do i = 1, ncol
    ia(i+1) = ia(i) + ia(i+1)
  end do
!
!  Loop for a cycle in chasing process.
!
  init = 1
  k = 0

 5    continue

  t = a(init)
  i = ja(init)
  j = iwk(init)
  iwk(init) = -1

 6 continue

   k = k + 1
!
!  Current row number is I.  Determine where to go.
!
  l = ia(i)
!
!  Save the chased element.
!
  t1 = a(l)
  inext = ja(l)
!
!  Then occupy its location.
!
  a(l) = t
  ja(l) = j
!
!  Update pointer information for next element to be put in row i.
!
  ia(i) = l + 1
!
!  Determine next element to be chased.
!
  if ( iwk(l) < 0 ) then
    go to 65
  end if

  t = t1
  i = inext
  j = iwk(l)
  iwk(l) = -1

  if ( k < nnz ) then
    go to 6
  end if

  do i = ncol, 1, -1
    ia(i+1) = ia(i)
  end do

  ia(1) = 1

  return

 65   continue

  init = init + 1

  if ( nnz < init ) then

    do i = ncol, 1, -1
      ia(i+1) = ia(i)
    end do

    ia(1) = 1

    return
  end if

  if ( iwk(init) < 0 ) then
    go to 65
  end if
!
!  Restart chasing.
!
  go to 5
end
subroutine udsol ( n, x, y, au, jau )

!*****************************************************************************80
!
!! UDSOL solves U*x = y;   U = upper triangular in MSR format
!
!  Discussion:
!
!    This routine solves a non-unit upper triangular matrix by standard 
!    sequential backward elimination.  The matrix is stored in MSR format,
!    with diagonal elements already inverted (otherwise do inversion,
!    au(1:n) = 1.0/au(1:n),  before calling).
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real Y(N), the right hand side.
!
! au,
! jau,    = Lower triangular matrix stored in modified sparse row
!          format.
!
!    Output, real X(N), the solution.
!
  implicit none

  integer n

  real ( kind = 8 ) au(*)
  integer j
  integer jau(*)
  integer k
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  x(n) = y(n) * au(n)

  do k = n-1, 1, -1
     t = y(k)
     do j = jau(k), jau(k+1)-1
       t = t - au(j) * x(jau(j))
     end do
     x(k) = au(k) * t
  end do

  return
end
subroutine udsolc ( n, x, y, au, jau )

!*****************************************************************************80
!
!! UDSOLC solves U * x = y, for upper triangular U in MSC format.
!
!  Discussion:
!
!    This routine solves a non-unit upper triangular system by standard
!    sequential forward elimination.  The matrix is stored in Modified
!    Sparse Column format with diagonal elements already inverted 
!    (otherwise do inversion,
!
!      au(1:n) = 1.0 / au(1:n),  
!
!    before calling this routine.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, real Y(N), contains the right hand side of the linear system.
!
! au,
! jau,   = Upper triangular matrix stored in Modified Sparse Column
!          format.
!
!    Output, real X(N), the solution of  U x = y .
!
  implicit none

  integer n

  real ( kind = 8 ) au(*)
  integer j
  integer jau(*)
  integer k
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  x(1:n) = y(1:n)

  do k = n, 1, -1
    x(k) = x(k) * au(k)
    t = x(k)
    do j = jau(k), jau(k+1)-1
      x(jau(j)) = x(jau(j)) - t * au(j)
    end do
  end do

  return
end
subroutine unassbl ( a, na, f, nx, nelx, ijk, nodcode, node, x, y, ierr, xyk )

!*****************************************************************************80
!
!! UNASSBL ???
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    a = un-assembled matrix on output
!
!    na = 1-st dimension of a.  a(na,node,node)
!
!    f = right hand side (global load vector) in un-assembled form
!
!    nx = number of nodes at input
!
!    nelx = number of elements at input
!
!    ijk = connectivity matrix: for node k, ijk(*,k) point to the
!    nodes of element k.
!
!    node = total number of nodal points in each element
!    also second dimension of a.
!
!    nodcode= boundary information list for each node with the
!    following meaning:
!    nodcode(i) = 0 -->  node i is internal
!    nodcode(i) = 1 -->  node i is a boundary but not a corner point
!    nodcode(i) = 2 -->  node i is a corner point (corner points
!
!    x,y = real arrays containing the $x$ and $y$ coordinates
!    resp. of the nodes.
!    K11, K22, and K12 at that element.
!
!    ierr = error message integer .
!    ierr = 0 --> normal return
!    ierr = 1 --> negative area encountered (due to bad
!    numbering of nodes of an element-
!    message printed in unit iout).
!
!    iout = output unit (not used here).
!
!    xyk  = routine defining the material properties at each
!    element. Form:
!    call xyk(nel,xyke,x,y,ijk,node)
!
  implicit none

  integer na
  integer node

  real ( kind = 8 ) a(na,node,node)
  real ( kind = 8 ) det
  real ( kind = 8 ) f(node,*)
  real ( kind = 8 ) fe(3)
  integer i
  integer ierr
  integer ijk(node,*)
  integer j
  integer ka
  integer kb
  integer nel
  integer nelx
  integer nodcode(*)
  integer nx
  real ( kind = 8 ) ske(3,3)
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xe(3)
  external xyk
  real ( kind = 8 ) xyke(2,2)
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) ye(3)
!
!  The maximum number of nonzeros allowed  = 200
!
!  Initialize.
!
  f(1:node,1:nx) = 0.0D+00
!
!  The main loop.
!
  do nel = 1, nelx
!
!  Get coordinates of nodal points.
!
    do i = 1, node
      j = ijk(i,nel)
      xe(i) = x(j)
      ye(i) = y(j)
    end do
!
!  Compute determinant.
!
     det = xe(2) * ( ye(3) - ye(1) ) &
         + xe(3) * ( ye(1) - ye(2) ) &
         + xe(1) * ( ye(2) - ye(3) )
!
!  Set material properties
!
    call xyk ( nel, xyke, x, y, ijk, node )
!
!  Construct element stiffness matrix
!
    ierr = 0
    call estif3 ( nel, ske, fe, det, xe, ye, xyke, ierr )

    if ( ierr /= 0 ) then
      return
    end if
!
!  Assemble: add element stiffness matrix to global matrix.
!
    do ka = 1, node
      f(ka,nel) = fe(ka)
      do kb = 1, node
        a(nel,ka,kb) = ske(ka,kb)
      end do
    end do

  end do

  return
end
subroutine usol ( n, x, y, au, jau, iau )

!*****************************************************************************80
!
!! USOL solves   U x = y    U = unit upper triangular.
!
!  Discussion:
!
!    This routine solves a unit upper triangular system by standard 
!    sequential backward elimination.  The matrix is stored in CSR format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
! y      = real array containg the right side.
!
! au,
! jau,
! iau,    = Lower triangular matrix stored in compressed sparse row
!          format.
!
! On return:
!
!      x = The solution of  U x = y .
!
  implicit none

  integer n

  real ( kind = 8 ) au(*)
  integer iau(n+1)
  integer j
  integer jau(*)
  integer k
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  x(n) = y(n)

  do k = n-1, 1, -1
    t = y(k)
    do j = iau(k), iau(k+1)-1
      t = t - au(j) * x(jau(j))
    end do
    x(k) = t
  end do

  return
end
subroutine usolc ( n, x, y, au, jau, iau )

!*****************************************************************************80
!
!! USOLC solves U * x = y for unit upper triangular U in CSC format.
!
!  Discussion:
!
!    This routine solves a unit upper triangular system by standard 
!    sequential forward elimination.  The matrix is stored in the CSC format.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
! y      = real array containg the right side.
!
! au,
! jau,
! iau,    = Uower triangular matrix stored in compressed sparse column
!          format.
!
! On return:
!
!      x  = The solution of  U x  = y.
!
  implicit none

  integer n

  real ( kind = 8 ) au(*)
  integer iau(*)
  integer j
  integer jau(*)
  integer k
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  x(1:n) = y(1:n)

  do k = n, 1, -1
    t = x(k)
    do j = iau(k), iau(k+1)-1
      x(jau(j)) = x(jau(j)) - t * au(j)
    end do
  end do

  return
end
