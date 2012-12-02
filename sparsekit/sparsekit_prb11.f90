program main

!*****************************************************************************80
!
!! MAIN is a finite element matrix generator.
!
!  Discussion:
!
!    This driver will generate a finite element matrix for the
!    conduction problem
!
!      -Div ( K(x,y) Grad u ) = f
!      u = 0 on boundary
!
!    (Dirichlet boundary conditions). The matrix is returned
!    assembled in compressed sparse row format. Unassembled matrices
!    can be generated (using genfeu) but this is not supported yet.
!
!    This driver will provide a few grids if wanted, with an
!    arbitrary number of levels of refinement (as permitted by the
!    sizes of the arrays as declared below).
!
!  Modified:
!
!    02 July 2005
!
!  Reference:
!
!    Noborou Kikuchi
!    Finite element methods in mechanics,
!    Cambridge University Press, 1986.
!
  implicit none

  integer, parameter :: maxnx = 2000
  integer, parameter :: maxnel = 4000

  real ( kind = 8 ) a(7*maxnx)
  real ( kind = 8 ) f(3*maxnx)
  real ( kind = 8 ) fs(3*maxnx)
  character ( len = 50 ) gridfile 
  character ( len = 2 ) guesol
  integer ia(maxnx)
  integer ichild(8,maxnel)
  integer ierr
  integer ifmt
  integer ii
  integer :: iin = 7
  integer ijk(3,maxnel)
  integer :: iout = 8
  integer iparnts(2,maxnx)
  integer iu
  integer iwk(maxnx)
  integer ja(7*maxnx)
  integer job
  integer jwk(maxnx)
  character ( len = 8 ) key
  character ( len = 50 ) matfile
  integer n
  integer n2
  integer :: na = 3000
  integer nb
  integer :: ndeg = 8
  integer nelmax
  integer nelx
  integer nelxnew
  integer ngrid
  integer nodcode(maxnx)
  integer :: node = 3
  integer nref
  integer nx
  integer nxmax
  integer nxnew
  character ( len = 72 ) title
  character ( len = 3 ) type
  real ( kind = 8 ) x(maxnx)
  external xyk
  real ( kind = 8 ) y(maxnx)      
!
!  Choose starting grid.
!
  iu = 10

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB11'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Demonstrate the SPARSKIT routines that can generate'
  write ( *, '(a)' ) '  test matrices based on finite element grids.'
  write ( *, '(a)' ) ' '
!
!  Force NGRID to be 1.  Normally, this would be chosen by the user.
!
  ngrid = 1
!
!  Generate the grid.
!
  if ( ngrid == 7 ) then
    write(*,*)'Grid type 7 : user provided initial grid '
    write(*,*)'Filename for initial grid :'
    read(*,'(A)') gridfile
    open ( unit = iin, file = gridfile, STATUS = 'old' )
  end if

  nx = 0
  nelx = 0
  call ingrid(ngrid,iin,nx,nelx,node,x,y,nodcode,ijk)

  if ( ngrid == 7 ) then
    close ( unit = iin )
  end if
!
!  Refine the grid.  Here we choose NREF=2, although this would
!  normally be interactively chosen.
!
  nref = 2
  nxmax = maxnx
  nelmax= maxnel
  nb = 0

  do ii = 1, nref
!
! estimate the number nx and nelx at next refinement level.
!
    call checkref(nx,nelx,ijk,node,nodcode,nb,nxnew,nelxnew)

    if ( nxmax < nxnew .or. nelmax < nelxnew ) then
      WRITE ( *, * ) 'Was able to do only ', ii-1 ,'  refinements'
      exit
    end if

    call refall(nx,nelx,ijk,node,ndeg,x,y,ichild,iparnts,nodcode, &
      nxmax, nelmax, ierr)

    if (ierr /= 0) then 
      WRITE ( *, * ) '** ERROR IN REFALL : ierr =',ierr
    end if

  end do

  job = 0

  call genfea ( nx, nelx, node, job, x, y, ijk, nodcode, fs, n2, &
    a, ja, ia, f, iwk, jwk, ierr, xyk )
!
!  Store matrix as a Harwell Boeing Matrix, by a call to prtmt.
!
  guesol = 'NN'
  title = ' FINITE ELEMENT TEST MATRIX FROM SPARSKIT            '
  type = 'RSA'
  key = 'SPARSKIT'
  ifmt = 104
  job = 2
  n = n2
!
!  Set the filename for the matrix data.
!
  matfile = 'test.mat'

  open ( unit = iout, file = matfile, STATUS = 'replace' )

  call prtmt(n,n,a,ja,ia,f,guesol,title,key,type, &
    ifmt,job,iout)

  close ( UNIT = IOUT )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The information about the matrix generated has been'
  write ( *, '(a)' ) '  stored in a file, using Harwell-Boeing format. '
  write ( *, '(a)' ) '  The file is "' // trim ( matfile ) // '".'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSEKIT_PRB11'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine ingrid ( ngrid, iin, nx, nelx, node, x, y, nodcode, ijk )

!*****************************************************************************80
!
!! INGRID initializes the grid according to the choice NGRID.
!
!  Discussion:
!
!    There are 6 initial grids provided and the user can
!   also enter his own grid as a seventh option.
!
! on entry:
!
! ngrid          = integer indicating the grid chosen. ngrid=1...6
!             corresponds to one of the 6 examples supplied by
!             SPARSKIT. ngrid = 7 is a user supplied initial grid.
!             see below for additional information for the format.
! iin       = integer containing the I/O unit number where to read
!             the data from in case ngrid = 7. A dummy integer
!             otherwise.
! node      = integer = the number of nodes per element (should be
!             set to three in this version). also the first dimension
!             of array ijk.
!
! on return
! 
! nx          = integer . the number of nodes
! nelx          = integer . the number of elements
! x, y      = two real arrays containing the coordinates of the nodes.
! nodcode   = an integer array containing the boundary information for
!             each node with the following meaning.
!      nodcode(i) = 0 -->  node i is internal
!      nodcode(i) = 1 -->  node i is a boundary but not a corner point
!      nodcode(i) = 2 -->  node i is a corner node.
!
! ijk(node,*)= an integer array containing the connectivity matrix.
!
!
! format for user supplied grid (when ngrid = 7)
!
! option 7 is a user defined initial grid.
!
! format is as follows:
! line 1: two integers, the first containing the number of nodes
!         the second the number of elements.
! line 2 to line nx+1:  node information
!        enter the following one line per node:
!        * the number of the node in the numbering chosen (integer from
!        taking the values 1 to nx), followed by
!        * the coordinates of the nodes (2 reals)  followed by
!       the boundary information, an integer taking one of the
!        values 0, 1, or 2,  with the meaning explained above.
!
! line nx+2 to nx+nelx+1: connectivity matrix
!       enter the following one line per element:
!       * number of the element in the numbering chosen, followed by
!       * The three numbers of the nodes (according to the above numbering
!       of the nodes) that constitute the element, in a counter clock-wise
!       order (this is in fact not important since it is checked by the
!       subroutine chkelemt).
!
! AN EXAMPLE: consisting on one single element (a triangle)
!------------
!    3    1
!    1    0.0000    0.0000    2
!    2    4.0000    0.0000    2
!    3    0.0000    4.0000    2
!    1    1    2    3
!
  implicit none

  integer node

  integer i
  integer ii
  integer iin
  integer ijk(node,*)
  integer j
  integer nelx
  integer ngrid
  integer nodcode(*)
  integer nx
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  if ( ngrid == 1 ) then

    call fgrid1 ( nx, nelx, node, x, y, nodcode, ijk )

  else if ( ngrid == 2 ) then

    call fgrid2 ( nx, nelx, node, x, y, nodcode, ijk )

  else if ( ngrid == 3 ) then

    call fgrid3 ( nx, nelx, node, x, y, nodcode, ijk )

  else if ( ngrid == 4 ) then

    call fgrid4 ( nx, nelx, node, x, y, nodcode, ijk )

  else if ( ngrid == 5 ) then

    call fgrid5 ( nx, nelx, node, x, y, nodcode, ijk )

  else if ( ngrid == 6 ) then

    call fgrid6 ( nx, nelx, node, x, y, nodcode, ijk )

  else if ( ngrid == 7 ) then

    read (iin,*) nx, nelx

    do i = 1, nx
      read(iin,*) ii, x(ii), y(ii), nodcode(ii)
    end do

    do i = 1, nelx

      read(iin,*) ii, (ijk(j,ii),j=1,node)
      nelx = max ( nelx, ii )

    end do

  end if

  call chkelmt ( nx, x, y, nelx, ijk, node )

  return
end
subroutine fgrid1 (nx,nelx,node,x,y,nodcode,ijk)

!*****************************************************************************80
!
!! FGRID1: initial grid for a simple square with two elements.
!
!      3             4
!       --------------
!       |          . |
!       |   2    .   |
!       |      .     |
!       |   .    1   |
!       | .          |
!       --------------
!      1              2
!
! input parameters: node = first dimensoin of ijk (must be .ge. 3)
! output parameters:
!    nx    = number of nodes
!    nelx = number of elemnts
!    (x(1:nx), y(1:nx)) = coordinates of nodes
!    nodcode(1:nx) = integer code for each node with the
!          following meening:
!      nodcode(i) = 0 -->  node i is internal
!      nodcode(i) = 1 -->  node i is a boundary but not a corner point
!      nodcode(i) = 2 -->  node i is a corner point.
!   ijk(1:3,1:nelx) = connectivity matrix. for a given element
!          number nel, ijk(k,nel), k=1,2,3 represent the nodes
!          composing the element nel.
!
  implicit none

  integer node

  integer ijk(node,*)
  integer, dimension ( 2 ) :: ijk1 = (/ 1, 1 /)
  integer, dimension ( 2 ) :: ijk2 = (/ 2, 4 /)
  integer, dimension ( 2 ) :: ijk3 = (/ 4, 3 /)
  integer k
  integer nelx
  integer nodcode(*)
  integer nx
  real ( kind = 8 ) x(*)
  real ( kind = 8 ), dimension ( 4 ) :: x1 = (/ 0.0, 1.0, 0.0, 1.0 /)
  real ( kind = 8 ) y(*)
  real ( kind = 8 ), dimension ( 4 ) :: y1 = (/ 0.0, 0.0, 1.0, 1.0 /)

  nx = 4

  do k = 1, nx
    x(k) = x1(k)
    y(k) = y1(k)
    nodcode(k) = 1
  end do

  nodcode(2) = 2
  nodcode(3) = 2

  nelx = 2

  do k = 1, nelx
    ijk(1,k) = ijk1(k)
    ijk(2,k) = ijk2(k)
    ijk(3,k) = ijk3(k)
  end do

  return
end
subroutine fgrid2 (nx,nelx,node,x,y,nodcode,ijk)

!*****************************************************************************80
!
!! FGRID2: initial grid for a simple D-shaped region with 4 elemnts
!       6
!       | .
!       |    .
!       |      .
!       |   4     .
!       |           .
!     4 -------------- 5
!       |          . |
!       |   3    .   |
!       |      .     |
!       |   .    2   |
!       | .          |
!       --------------
!       | 2         . 3
!       |         .
!       |   1   .
!       |     .
!       |  .
!       |.
!       1
!--------------------------------------------------------------
! input parameters: node = first dimensoin of ijk (must be .ge. 3)
! output parameters:
!    nx    = number of nodes
!    nelx = number of elemnts
!    (x(1:nx), y(1:nx)) = coordinates of nodes
!    nodcode(1:nx) = integer code for each node with the
!          following meening:
!      nodcode(i) = 0 -->  node i is internal
!      nodcode(i) = 1 -->  node i is a boundary but not a corner point
!      nodcode(i) = 2 -->  node i is a corner point.
!   ijk(1:3,1:nelx) = connectivity matrix. for a given element
!          number nel, ijk(k,nel), k=1,2,3 represent the nodes
!          composing the element nel.
!
  implicit none

  integer node

  integer ijk(node,*)
  integer ijk1(4)
  integer ijk2(4)
  integer ijk3(4)
  integer k
  integer nelx
  integer nodcode(*)
  integer nx
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) x1(6)
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) y1(6)
!
! coordinates of nodal points
!
        data x1/0.0, 0.0, 1.0, 0.0, 1.0, 0.0/
        data y1/0.0, 1.0, 1.0, 2.0, 2.0, 3.0/
!
!------------------|--|--|--|
! elements         1  2  3  4
!------------------|--|--|--|
      data ijk1   /1, 2, 2, 4/
      data ijk2   /3, 3, 5, 5/
      data ijk3   /2, 5, 4, 6/

  nx = 6

  do k = 1, nx
    x(k) = x1(k)
    y(k) = y1(k)
    nodcode(k) = 1
  end do

  nelx = 4

  do k = 1, nelx
    ijk(1,k) = ijk1(k)
    ijk(2,k) = ijk2(k)
    ijk(3,k) = ijk3(k)
  end do

  return
end
subroutine fgrid3 (nx,nelx,node,x,y,nodcode,ijk)

!*****************************************************************************80
!
!! FGRID3: initial grid for a C-shaped region composed of 10 elements --
!
!
!      10           11            12
!       ---------------------------
!       |          . |          . |
!       |  7     .   |   9    .   |
!       |      .     |      .     |
!       |   .    8   |   .   10   |
!       | .          | .          |
!     7 ---------------------------
!       |          . |8           9
!       |   5    .   |
!       |      .     |
!       |   .    6   |
!     4 | .          |5           6
!       ---------------------------
!       |          . |          . |
!       |   1    .   |  3     .   |
!       |      .     |      .     |
!       |   .    2   |   .   4    |
!       | .          | .          |
!       ---------------------------
!      1             2            3
!
!
! input parameters: node = first dimensoin of ijk (must be .ge. 3)
!    nx    = number of nodes
!    nelx = number of elemnts
!    (x(1:nx), y(1:nx)) = coordinates of nodes
!    nodcode(1:nx) = integer code for each node with the
!          following meening:
!      nodcode(i) = 0 -->  node i is internal
!      nodcode(i) = 1 -->  node i is a boundary but not a corner point
!      nodcode(i) = 2 -->  node i is a corner point.
!   ijk(1:3,1:nelx) = connectivity matrix. for a given element
!          number nel, ijk(k,nel), k=1,2,3 represent the nodes
!          composing the element nel.
!
  implicit none

  integer node

  integer ijk(node,*)
  integer ijk1(10)
  integer ijk2(10)
  integer ijk3(10)
  integer k
  integer nelx
  integer nodcode(*)
  integer nx
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) x1(12)
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) y1(12)
!
! coordinates of nodal points
!
        data x1/0.0,1.0,2.0,0.0,1.0,2.0,0.0,1.0,2.0,0.0,1.0,2.0/
        data y1/0.0,0.0,0.0,1.0,1.0,1.0,2.0,2.0,2.0,3.0,3.0,3.0/
!
!------------------|--|--|--|--|--|--|---|---|---|
! elements         1  2  3  4  5  6  7   8   9  10
!------------------|--|--|--|--|--|--|---|---|---|
      data ijk1   /1, 1, 2, 2, 4, 4, 7,  7,  8, 8/
      data ijk2   /5, 2, 6, 3, 8, 5, 11, 8, 12, 9/
      data ijk3   /4, 5, 5, 6, 7, 8, 10, 11,11, 12/

  nx = 12

  do k = 1, nx
    x(k) = x1(k)
    y(k) = y1(k)
    nodcode(k) = 1
  end do

  nodcode(3) = 2
  nodcode(10) = 2
  nodcode(9) = 2

  nelx = 10

  do k = 1, nelx
    ijk(1,k) = ijk1(k)
    ijk(2,k) = ijk2(k)
    ijk(3,k) = ijk3(k)
  end do

  return
end
subroutine fgrid4 (nx,nelx,node,x,y,nodcode,ijk)

!*****************************************************************************80
!
!! FGRID4: initial grid for a C-shaped region composed of 10 elements --
!      10                   11
!       +------------------+ .
!       | .                |    .
!       |    .       8     |       . 12
!       |        .         |  9   . |
!       |     7      .     |   .    |
!     7 |                . | .   10 |
!       -------------------+--------+ 9
!       |                 .| 8
!       |     5       .    |
!       |         .        |
!       |    .       6     |
!       |.                 | 5      6
!    4  +------------------+--------+
!       |               .  | .   4  |
!       |    1       .     |    .   |
!       |        .         |  3    .| 3
!       |    .        2    |    .
!       | .                | .
!       --------------------
!       1                  2
!
! input parameters: node = first dimensoin of ijk (must be .ge. 3)
!    nx    = number of nodes
!    nelx = number of elemnts
!    (x(1:nx), y(1:nx)) = coordinates of nodes
!    nodcode(1:nx) = integer code for each node with the
!          following meening:
!      nodcode(i) = 0 -->  node i is internal
!      nodcode(i) = 1 -->  node i is a boundary but not a corner point
!      nodcode(i) = 2 -->  node i is a corner point.
!   ijk(1:3,1:nelx) = connectivity matrix. for a given element
!          number nel, ijk(k,nel), k=1,2,3 represent the nodes
!          composing the element nel.
!
  implicit none

  integer node

  integer ijk(node,*)
  integer ijk1(10)
  integer ijk2(10)
  integer ijk3(10)
  integer k
  integer nelx
  integer nodcode(*)
  integer nx
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) x1(12)
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) y1(12)
!
! coordinates of nodal points
!
        data x1/0.0,1.0,1.5,0.0,1.0,1.5,0.0,1.0,1.5,0.0,1.0,1.5/
        data y1/0.0,0.0,0.5,1.0,1.0,1.0,2.0,2.0,2.0,3.0,3.0,2.5/
!
!------------------|--|--|--|--|--|--|---|---|---|
! elements         1  2  3  4  5  6  7   8   9  10
!------------------|--|--|--|--|--|--|---|---|---|
      data ijk1   /1, 1, 2, 5, 4, 4, 7, 10,  8, 8/
      data ijk2   /5, 2, 3, 3, 8, 5, 8,  8, 12, 9/
      data ijk3   /4, 5, 5, 6, 7, 8, 10, 11,11, 12/

  nx = 12

  do k = 1, nx
    x(k) = x1(k)
    y(k) = y1(k)
    nodcode(k) = 1
  end do

  nodcode(6) = 2
  nodcode(9) = 2

  nelx = 10

  do k = 1, nelx
    ijk(1,k) = ijk1(k)
    ijk(2,k) = ijk2(k)
    ijk(3,k) = ijk3(k)
  end do

  return
end
subroutine fgrid5 (nx,nelx,node,x,y,nodcode,ijk)

!*****************************************************************************80
!
!! FGRID5: initial grid for a whrench shaped region composed of 14 elements --
!
!                                      13            15
!                                        . ----------.           |-3
!                                      .   .   13  .   .         |
!                                   .   12   .   .  14    .      |
! 9        10        11       12  .            . 14        . 16  |
! ----------------------------------------------------------     |-2
! |       . |       . |       . |            . |                 |
! | 1   .   |  3  .   |  5  .   |    7   .     |                 |
! |   .  2  |   .  4  |   .  6  |     .    8   |                 |
! |.        |.        |.        | .            |                 |
! -----------------------------------------------------------    |-1
! 1         2         3       4  .           6 .           . 8   |
!                                   .   9    .   .   11   .      |
!                                      .   .  10    .   .        |
!                                        .___________.           |-0
!                                       5             7
!
! 0---------1--------2----------3--------------4-------------5
!
! input parameters: node = first dimensoin of ijk (must be .ge. 3)
!    nx    = number of nodes
!    nelx = number of elemnts
!    (x(1:nx), y(1:nx)) = coordinates of nodes
!    nodcode(1:nx) = integer code for each node with the
!          following meening:
!      nodcode(i) = 0 -->  node i is internal
!      nodcode(i) = 1 -->  node i is a boundary but not a corner point
!      nodcode(i) = 2 -->  node i is a corner point.
!   ijk(1:3,1:nelx) = connectivity matrix. for a given element
!          number nel, ijk(k,nel), k=1,2,3 represent the nodes
!          composing the element nel.
!
  implicit none

  integer node

  integer ijk(node,*)
  integer ijk1(14)
  integer ijk2(14)
  integer ijk3(14)
  integer k
  integer nelx
  integer nodcode(*)
  integer nx
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) x1(16)
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) y1(16)
!
! coordinates of nodal points
!
        data x1/0.,1.,2.,3.,3.5,4.,4.5,5.,0.,1.,2.,3.,3.5,4.,4.5,5./
        data y1/1.,1.,1.,1.,0.,1.,0.,1.,2.,2.,2.,2.,3.,2.,3.,2./
!
!------------------|--|--|--|--|--|--|---|---|---|--|---|---|---|
! elements         1  2  3  4  5  6  7   8   9  10  11  12  13  14
!------------------|--|--|--|--|--|--|---|---|---|--|---|---|---|
      data ijk1   /1, 1, 2, 2, 3, 3, 4,  4,  4,  5, 6, 12, 14, 14/
      data ijk2   /10,2,11, 3,12, 4,14,  6,  5,  7, 7, 14, 15, 16/
      data ijk3   /9,10,10,11,11,12,12, 14,  6,  6, 8, 13, 13, 15/

  nx = 16

  do k=1, nx
    x(k) = x1(k)
    y(k) = y1(k)
    nodcode(k) = 1
  end do

  nodcode(9) = 2
  nodcode(8) = 2
  nodcode(16) = 2

  nelx = 14

  do k=1,nelx
    ijk(1,k) = ijk1(k)
    ijk(2,k) = ijk2(k)
    ijk(3,k) = ijk3(k)
  end do

  return
end
subroutine fgrid6 (nx,nelx,node,x,y,nodcode,ijk)

!*****************************************************************************80
!
!! FGRID6 generates a random finite element grid. 
!
!  Discussion:
!
!    Random coordinates are generated by using the library random number
!    generator and then a Delauney triangulation is used to generate the grid.
!
!    The algorithm used for the triangulation is coded by Sweby.
!
  implicit none

  integer node

  integer adjlist(200,12)
  integer i
  integer i1
  integer i2
  integer ijk(node,*)
  integer ijktr(200,3)
  integer il(6)
  integer j
  integer jj
  integer k
  integer nadj(200)
  integer nbr
  integer nel(200)
  integer nelx
  integer nemax
  integer nod
  integer nodcode(*)
  integer nx
  real ( kind = 8 ) random
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  nx = 20

  do j = 1, nx
    x(j) = random()
  end do

  do j = 1, nx
    y(j) = random()
  end do

  nemax = 200
  call dlauny ( x, y, nx, ijktr, nemax, nelx )

  print *, ' delauny -- nx, nelx ', nx, nelx

  do j = 1, nx
    nel(j) = 0
    nadj(j) = 0
  end do
!
!  transpose ijktr into ijk and count the number of
!  elemnts to which each node belongs.
!
  do j = 1, nelx
    do k = 1, node
      i = ijktr(j,k)
      ijk(k,j) = i
      nel(i) = nel(i)+1
    end do
  end do
!
!  Take care of ordering within each element.
!
  call chkelmt (nx, x, y, nelx, ijk, node)
!
!  The next blocks are to determine  the nature of each point.
!  (interior point, boundary point, corner point.
!
!  List and count the neighbors of each node.
!
  do j = 1, nelx

    do k=1, node
      il(k) = ijk(k,j)
      il(k+node) = il(k)
    end do
!
!  neighbors of node il(k) are il(k+1), il(k+2)
!
    do k = 1, node

      nod = il(k)
      i1 = il(k+1)
      i2 = il(k+2)
!
!  see if already there
!
      nbr = nadj(nod)

      do jj = 1, nbr

        if ( adjlist(nod,nbr) == i1 ) then 
          i1 = 0
        end if

        if ( adjlist(nod,nbr) == i2 ) then
          i2 = 0
        end if

      end do

      if ( i1 /= 0 ) then
        nbr = nbr + 1
        adjlist(nod,nbr) = i1
      end if

      if ( i2 /= 0 ) then
        nbr = nbr + 1
        adjlist(nod,nbr) = i2
      end if

      nadj(nod) = nbr
 
    end do

  end do
!
!  Boundary info:
!  if number of neighbors = number of elemnts to which it belongs then it is
!  an internal point
!  if not but number of neighnors >= 2 then boundary point
!  if nadj(k) = 2 then corner point.
!
  do j=1, nx

    nodcode(j) = 0
    nbr = nadj(j)

    if (nel(j) < nbr) then
      nodcode(j) = 1
    end if

    if ( nbr == 2 ) then
      nodcode(j) = 2
    end if

  end do

  return
end
function random ()

!*****************************************************************************80
!
!! RANDOM returns a pseudorandom value.
!
!  This routine was extracted from ELEFUNT.
!
  implicit none

  integer, save :: iy = 100001
  real ( kind = 8 ) random

  iy = iy * 125
  iy = iy - ( iy / 2796203 ) * 2796203
  random = real ( iy, kind = 8 ) / 2796203.0

  return
end
subroutine xyk ( nel, xyke, x, y, ijk, node )

!*****************************************************************************80
!
!! XYK evaluates the material property matrix function.
!
!  Discussion:
!
!    In this example, the function is just the identity matrix.
!
  implicit none

  integer node

  integer ijk(node,*)
  integer nel
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xyke(2,2)
  real ( kind = 8 ) y(*)

  xyke(1,1) = 1.0
  xyke(2,2) = 1.0
  xyke(1,2) = 0.0
  xyke(2,1) = 0.0

  return
end
