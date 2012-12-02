program main

!*****************************************************************************80
!
!! MAIN is a test suite for the Formats routines.
!
! tests all of the routines in the module formats.
!
! Note: the comments may not have been updated.
!
! Here is the sequence of what is done by this program.
! 1) generates a block matrix associated with a simple 5-point
!    matrix on a 4 x 2 grid (2-D) with 2 degrees of freedom per
!    grid point. Thus N = 16. This is produced in block format
!    using the generation routine genpbl.
! 2) the block format is translated into a compressed sparse row
!    format by bsrcsr. The result is dumped in file csr.mat
! 3) the matrix is translated in dense format by csrdns.
!    result a 16 x 16 matrix is written in unit dns.mat.
!    This is a good file to look at to see what the matrix is
!    and to compare results of other formats with.
! 4) the dense matrix obtained in 3) is reconverted back to
!    csr format using dnscsr. Result appended to file csr.mat
! 5) The matrix obtained in 4) is converted in coordinate format
!    and the resulting matrix is written in file coo.mat
! 6) the result is converted back to csr format. matrix
!    appended to csr.mat.
! 7) result of 6) is converted to symmetric sparse row storage
!    (ssr) and the result is appended to csr.mat
! 8) result of 7) converted back to csr format and result is
!    appended to csr.mat
! 9) matrix resulting from 8) is converted to modified sparse
!    row format using csrmsr and result is written in msr.mat.
!10) the resulting matrix is converted back to csrformat and
!    result is appended to csr.mat
!11) result is converted to ellpack-itpack format with
!    csrell and result is printed in ell.mat
!12) result is converted back to csr format and appended to csr.mat
!12) result converted to csc format (transposition) using csrcsc
!    which should produce the same matrix here. result appended
!    to csr.mat. A second call to csrcsc is made on resulting
!    matrix.
!13) the subroutine csrdia is used to extract two diagonals
!    (offsets -1 and 0) and then all the diagonals of matrix.
!    results in dia.mat
!14) diacsr is then called to convert the diagonally stored matrix
!    back to csr format. result appended to csr.mat
!15) result is converted to band format (bnd) by calling
!    csrbnd. result dumped to bnd.mat
!16) result is converted back to csr format and appended to csr.mat
!17) result sorted by a call to csrcsc and then converted to
!    block format (csrbsr) and then back to csr format again.
!    result appedned to csr.mat.
!18) matrix converted to symmetric skyline format. result appended
!    to file band.mat
!19) matrix converted back to csr format and result appended to
!    csr.mat.
!20) result converted to jad format. result output in jad.mat
!21) result concverted back to csr fromat. appended to csr.mat
!
  implicit none

  integer, parameter :: nxmax = 10
  integer, parameter :: nmx = nxmax * nxmax
  integer, parameter :: nnzmax = 10 * nmx
  integer, parameter :: ndns = 20

  real ( kind = 8 ) a(nnzmax)
  real ( kind = 8 ) a1(nnzmax)
  real ( kind = 8 ) dns(ndns,ndns)
  integer i
  integer ia(nmx+1)
  integer ia1(nnzmax)
  integer idiag
  integer idiag0
  integer ierr
  integer imod
  integer ioff(20)
  integer iout
  integer iwk(nmx*2+1)
  integer j
  integer ja(nnzmax)
  integer ja1(nnzmax)
  integer job
  integer k
  integer k1
  integer k2
  integer kdiag
  integer kend
  integer kstart
  integer len
  integer lowd
  integer maxcol
  integer ml
  integer mu
  integer n
  integer na
  integer ndiag
  integer nel
  integer nfree
  integer nnz
  integer nx
  integer ny
  integer nz
  real ( kind = 8 ) stencil(7,100)
  real ( kind = 8 ) wk(nmx)

  call timestamp ( )

  WRITE(*,*)' '
  WRITE(*,*)'SPARSEKIT_PRB03'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  WRITE(*,*)'A set of tests for SPARSKIT'
  WRITE(*,*)' '

  open (unit=7,file='csr.mat',STATUS='replace')
  open (unit=8,file='dns.mat',STATUS='replace')
  open (unit=9,file='coo.mat',STATUS='replace')
  open (unit=10,file='msr.mat',STATUS='replace')
  open (unit=11,file='ell.mat',STATUS='replace')
  open (unit=12,file='dia.mat',STATUS='replace')
  open (unit=13,file='bnd.mat',STATUS='replace')
  open (unit=14,file='jad.mat',STATUS='replace')
!
!  dimension of grid
!
  nx = 4
  ny = 2
  nz = 1
  nfree = 2
!
!  generate grid problem.
!
  na = nx*ny*nz*5
  call gen57bl (nx,ny,nz,nfree,na,n,a1,ja1,ia1,iwk,stencil)
!
!  write out the matrix
!
      call bsrcsr (n,nfree,na,a1,ja1,ia1,a,ja,ia)
      iout = 7
      nnz = ia(n+1)-1
      write (iout,*) '-----------------------------------------'
      write (iout,*) '  +++  initial matrix in CSR format +++ '
      write (iout,*) '-----------------------------------------'
        call dump (n,a,ja,ia,6)
!
!  call csrdns
!
       call csrdns(n,n,a,ja,ia,dns,ndns,ierr)
       iout = iout+1
       write (iout,*) '-----------------------------------------'
         write (iout,*) '  +++  initial matrix in DENSE format+++ '
       write (iout,*) '-----------------------------------------'
       write (iout,'(4x,16i4)') (j,j=1,n)
       write (iout,'(3x,65(1h-))')
       do i=1, n
          write (8,102) i,(dns(i,j), j=1,n)
 102          format(1h ,i2,1h|,16f4.1)
       end do
!
!  convert it back to sparse format.
!
        call dnscsr(n,n,nnzmax,dns,ndns,a1,ja1,ia1,ierr)
          write(*,*) '-----------------------------------------'
        write(*,*) '  +++ matrix after conversion from dnscsr +++ '
         write(*,*) '-----------------------------------------'
        if (ierr .ne. 0) write(*,*) ' ***** ERROR FROM DNSCSR'
        if (ierr .ne. 0) write(*,*)  '     IERR = ', ierr
        call dump (n,a1,ja1,ia1,6)
!
!  convert it to coordinate format.
!
          call csrcoo(n,3,nnzmax,a,ja,ia,nnz,a1,ia1,ja1,ierr)
        iout = iout+ 1
        if (ierr .ne. 0) write (iout,*) ' ***** ERROR IN CSRCOO'
        if (ierr .ne. 0) write (iout,*)   '     IERR = ', ierr
          write (iout,*) '-----------------------------------------'
        write(iout,*) ' +++ Matrix in coordinate format +++ '
         write (iout,*) '-----------------------------------------'
          write(iout,103) (ia1(j),ja1(j),a1(j),j=1,nnz)
 103      format (' i =', i3,'    j = ',i3,'     a(i,j) = ',f4.1)
!
!  convert it back again to csr format
!
          call coocsr(n,nnz,a1,ia1,ja1,a,ja,ia)

         write(*,*) '-----------------------------------------'
        write(*,*) '  +++ matrix after conversion from coocsr +++ '
         write(*,*) '-----------------------------------------'
          call dump(n,a,ja,ia,6)
!
!  going to srs format
!
          call csrssr(n,a,ja,ia,nnzmax,a1,ja1,ia1,ierr)
         write(*,*) '-----------------------------------------'
        write(*,*) '  +++ matrix after conversion to ssr format +++ '
        write(*,*) '      (lower part only stored in csr format)    '
         write(*,*) '-----------------------------------------'
         call dump(n,a1,ja1,ia1,6)
!
!  back to csr
!
          call ssrcsr (n,a1,ja1,ia1,nnzmax,a,ja,ia,iwk,ierr)
          if (ierr .ne. 0) write(7,*) ' error in ssrcsr-IERR=',ierr
         write(*,*) '-----------------------------------------'
        write(*,*) '  +++ matrix after conversion from ssrcsr +++ '
         write(*,*) '-----------------------------------------'
         call dump(n,a,ja,ia,6)
!
!  msr format
!
         iout = iout+1
         call csrmsr (n,a,ja,ia,a1,ja1,a1,ja1)
        write (iout,*) '-----------------------------------------'
         write (iout,*) '  +++ matrix in modified sparse row format +++'
        write (iout,*) '-----------------------------------------'
         write (iout,*) ' ** MAIN DIAGONAL '
         write (iout,'(16f4.1)') (a1(k),k=1,n)
            write (iout,*) ' ** POINTERS: '
         write (iout,'(17i4)') (ja1(k),k=1,n+1)
             write (iout,*) ' ** REMAINDER :'
         call dump(n,a1,ja1,ja1,iout)

         call msrcsr (n,a1,ja1,a,ja,ia,wk)
        write(*,*) '-----------------------------------------'
        write(*,*) ' +++ matrix after conversion from msrcsr +++'
        write(*,*) '-----------------------------------------'

           call dump(n,a,ja,ia,6)

        maxcol = 13
!
        call csrell (n,a,ja,ia,maxcol,a1,ja1,n,ndiag,ierr)
        iout = iout+1
        if (ierr .ne. 0) write (iout,*) ' ***** ERROR IN CSRELL'
        if (ierr .ne. 0) write (iout,*)   '     IERR = ', ierr
        write (iout,*) '-----------------------------------------'
        write (iout,*) '  +++ matrix in ELLPACK-ITPACK format +++ '
        write (iout,*) '-----------------------------------------'
         do i=1,ndiag
         write (iout,*) ' Column number: ', i
         write (iout,104) (a1(n*(i-1)+k),k=1,n)
 104         format(9h COEF  = ,16f4.0)
         write (iout,105) (ja1(n*(i-1)+k),k=1,n)
 105         format (9h JCOEF = ,16i4)
        end do
        call ellcsr (n,a1,ja1,n,ndiag,a,ja,ia,nnzmax,ierr)
        if (ierr .ne. 0) write(*,*) ' ***** ERROR IN ELLCSR'
        if (ierr .ne. 0) write(*,*)   '     IERR = ', ierr
        write(*,*) '-----------------------------------------'
         write(*,*) '  +++ matrix after conversion from ellcsr +++'
        write(*,*) '-----------------------------------------'
         call dump(n,a,ja,ia,6)

       call csrcsc(n,1,1,a,ja,ia, a1,ja1,ia1)
        write(*,*) '-----------------------------------------'
         write(*,*) '  +++ matrix after conversion from csrcsc  +++ '
        write(*,*) '-----------------------------------------'
       call dump(n,a1,ja1,ia1,6)
       call csrcsc(n,1,1,a1,ja1,ia1, a,ja,ia)
!
!  test 1:
!  get main diagonal and subdiagonal
!  get some info on diagonals
!
       call infdia(n,ja,ia,iwk,idiag0)
       job = 0
       ioff(1) = 0
       ioff(2) = -1
       idiag = 2
       call csrdia (n,idiag,job,a,ja,ia,ndns,dns,ioff,a1,ja1,ia1,iwk)
       iout = iout+1
              write (iout,*) '-----------------------------------------'
         write (iout,*) '  +++  diagonal format +++ '
        write (iout,*) '-----------------------------------------'
         write (iout,*) '  diagonals ioff = 0 and ioff = -1 '
       write (iout,*) ' number of diag.s returned from csrdia=',idiag
       do kdiag = 1, idiag
         write (iout,*) ' diagonal offset = ', ioff(kdiag)
         write (iout,'(16f4.1)') (dns(k,kdiag),k=1,n)
       end do
!
! reverse conversion
!
       ndiag = ndns
       idiag = idiag0
       job = 10
       call csrdia (n,idiag,job,a,ja,ia,ndns,dns,ioff,a1,ja1,ia1,iwk)
       write (iout,*) '-----------------------------------------'
         write (iout,*) '  +++  second test diagonal format +++ '
         write (iout,*) '         ** all diagonals of A  ** '
       write (iout,*) '-----------------------------------------'
       write (iout,*) ' number of diagonals on return from csrdia=', &
                        idiag
       do kdiag = 1, idiag
         write (iout,*) ' diagonal offset = ', ioff(kdiag)
         write (iout,'(16f4.1)') (dns(k,kdiag),k=1,n)
       end do
!
!  reverse conversion
!
       job = 0
       call diacsr (n,job,idiag,dns,ndns,ioff,a,ja,ia)

        write(*,*) '-----------------------------------------'
         write(*,*) '  +++ matrix after conversion from diacsr  +++ '
        write(*,*) '-----------------------------------------'
                call dump(n,a,ja,ia,6)
!
!  checking the banded format
!
      lowd = 0
      job = 1
      call csrbnd(n,a,ja,ia,job,dns,ndns,lowd,ml,mu,ierr)
        iout = iout+1
        if (ierr .ne. 0) write (iout,*) ' ***** ERROR IN CSRBND'
        if (ierr .ne. 0) write (iout,*)   '     IERR = ', ierr
        write (iout,*) '-----------------------------------------'
         write (iout,*) '       +++  banded  format +++ '
         write (iout,*) ' bandwidth values found ml=',ml,'  mu=',mu
       write (iout,*) '-----------------------------------------'
       write (iout,'(4x,16i4)') (j,j=1,n)
       write (iout,'(3x,65(1h-))')
       do i=1, lowd
         write (iout,102) i, (dns(i,j), j=1,n)
       end do
!
!  convert back to a, ja, ia format.
!
      len = nnzmax
      call bndcsr(n,dns,ndns,lowd,ml,mu,a,ja,ia,len,ierr)
        write(*,*) ' IERR IN BNDCSR = ', ierr
        write(*,*) '-----------------------------------------'
        write(*,*) '  +++ matrix after conversion from bndcsr +++'
        write(*,*) '-----------------------------------------'
      call dump(n,a,ja,ia,6)
!
!  make sure it is sorted
!
        call csrcsc(n,1,1,a,ja,ia, a1,ja1,ia1)
!
!  checking skyline format.
!
         imod = 1
       call csrssk (n,imod,a1,ja1,ia1,a,ia,nnzmax,ierr)
       if (ierr .ne. 0) write (iout,*)   '     IERR = ', ierr
        write (iout,*) '-----------------------------------------'
         write (iout,*) '    +++ Sym. Skyline format +++ '
       write (iout,*) '-----------------------------------------'
       write (iout,'(3x,65(1h-))')
!
!  create column values.
!
       do i=1, n
         kend = ia(i+1)-1
         kstart = ia(i)
         do k=kstart,kend
           ja(k) =  i-(kend-k)
         end do
      end do

         call dump(n,a,ja,ia,iout)
!
!  back to ssr format.
!
      call sskssr (n,imod,a,ia,a1,ja1,ia1,nnzmax,ierr)
        write(*,*) '-----------------------------------------'
        write(*,*) '  +++ matrix after conversion from sskcsr +++'
       write(*,*) '-----------------------------------------'
      call dump(n,a1,ja1,ia1,6)
!
!  checking jad format.
!
!  first go back to the csr format ----
!
       call ssrcsr (n,a1,ja1,ia1,nnzmax,a,ja,ia,iwk,ierr)

       call csrjad (n, a, ja, ia, ndiag, iwk, a1, ja1, ia1)

       iout = iout+1
        write (iout,*) '-----------------------------------------'
        write (iout,*) '   +++   matrix in JAD format +++ '
        write (iout,*) '-----------------------------------------'
!
!  permutation array
!
                write (iout,*) ' ** PERMUTATION ARRAY '
         write (iout,'(17i4)') (iwk(k),k=1,n)
!
!  diagonals.
!
         do i=1,ndiag
           write (iout,*) ' J-diagonal number: ', i
           k1 = ia1(i)
           k2 = ia1(i+1)-1
           write (iout,104) (a1(k),k=k1,k2)
           write (iout,105) (ja1(k),k=k1,k2)
         end do
!
!  back to csr format.
!
        call jadcsr (n, ndiag, a1, ja1, ia1, iwk, a, ja, ia)

        write(*,*) '-----------------------------------------'
        write(*,*) '  +++ matrix after conversion from jadcsr +++'
       write(*,*) '-----------------------------------------'
      call dump(n,a,ja,ia,6)
!
! checking the linked list format
!
       nnz = ia(n+1) - ia(1)
       call csrlnk (n, a, ja, ia, iwk)
!
!  print links in file 7 (no need for another file)
!
       iout = 7
       write (iout,*) '-----------------------------------------'
       write (iout,*) '   +++   matrix in LNK format +++ '
       write (iout,*) '-----------------------------------------'
!
!  permutation array
!
                write (iout,*) ' LINK ARRAY '
           write (iout,*) ' ---------- '
         write (iout,'(17i4)') (iwk(k),k=1,nnz)
!
!  back to csr format..
!
        call lnkcsr (n, a, ja, ia, iwk, a1, ja1, ia1)

        write(*,*) '-----------------------------------------'
        write(*,*) '  +++ matrix after conversion from lnkcsr +++'
       write(*,*) '-----------------------------------------'
      call dump(n,a,ja,ia,6)

      CLOSE(UNIT=7)
      CLOSE(UNIT=8)
      CLOSE(UNIT=9)
      CLOSE(UNIT=10)
      CLOSE(UNIT=11)
      CLOSE(UNIT=12)
      CLOSE(UNIT=13)
      CLOSE(UNIT=14)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB03'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function afun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) afun, x,y, z

  afun = -1.0

  return
end
function bfun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) bfun, x,y, z

  bfun = -1.0

  return
end
function cfun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) cfun, x,y, z

  cfun = -1.0

  return
end
function dfun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) dfun
  real ( kind = 8 ) gamma
  real ( kind = 8 ) x,y, z
  data gamma /100.0/
!        dfun = gamma*dexp(x*y)

  dfun = 10.0

  return
end
function efun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) efun
  real ( kind = 8 ) gamma
  real ( kind = 8 ) x,y, z

  data gamma /100.0/
!        efun = gamma*dexp(-x*y)

  efun = 0.0

  return
end
function ffun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) ffun, x,y, z

  ffun = 0.0

  return
end
function gfun (x,y,z)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) gfun, x,y, z

  gfun = 0.0

  return
end
subroutine afunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0
    end do
    coeff((j-1)*nfree+j) = -1.0
  end do

  return
end
subroutine bfunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0
    end do
    coeff((j-1)*nfree+j) = -1.0
  end do

  return
end
subroutine cfunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0
    end do
    coeff((j-1)*nfree+j) = -1.0
  end do

  return
end
subroutine dfunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0
    end do
  end do

  return
end
subroutine efunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0
    end do
  end do

  return
end
subroutine ffunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0
    end do
  end do

  return
end
subroutine gfunbl (nfree,x,y,z,coeff)

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) coeff(100)
  integer i
  integer j
  integer nfree
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  do j=1, nfree
    do i=1, nfree
      coeff((j-1)*nfree+i) = 0.0
    end do
  end do

  return
end
subroutine xyk(nel,xyke,x,y,ijk,node)

!*****************************************************************************80
!
! The material property function xyk for the
! finite element problem
!
  implicit none

  integer node

  integer ijk(node,*)
  integer nel
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xyke(2,2)
  real ( kind = 8 ) y(*)
!
!  this is the identity matrix.
!
  xyke(1,1) = 1.0
  xyke(2,2) = 1.0
  xyke(1,2) = 0.0
  xyke(2,1) = 0.0

  return
end
