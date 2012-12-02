program main

!*****************************************************************************80
!
!! MAIN is the main program for TABLE_TET_MESH.
!
!  Discussion:
!
!    TABLE_TET_MESH computes a tet mesh of a 3D TABLE dataset.
!
!  Usage:
!
!    table_tet_mesh input_file_name
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Local, integer ( kind = 4 ), parameter BF_MAX, the maximum number of boundary 
!    faces.  I don't know a reasonable formula for this quantity.
!    If it's not large enough, the program will print a warning message.
!
!    Local, integer ( kind = 4 ), parameter FC_MAX, the maximum number of faces.  
!    I don't know a reasonable formula for this quantity.
!    If it's not large enough, the program will print a warning message.
!
  implicit none

  integer   ( kind = 4 ), parameter :: bf_max = 8000
  integer   ( kind = 4 ), parameter :: fc_max = 40000

  integer   ( kind = 4 ) arg_num
  integer   ( kind = 4 ), dimension (1:3,1:bf_max) :: bf
  integer   ( kind = 4 ) bf_num
  integer   ( kind = 4 ) dim_num
  integer   ( kind = 4 ) face_num
  integer   ( kind = 4 ), dimension (1:7,1:fc_max) :: fc
  integer   ( kind = 4 ) fc_num
  integer   ( kind = 4 ), allocatable, dimension ( : ) :: ht
  integer   ( kind = 4 ) ht_num
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) iarg
  integer   ( kind = 4 ) iargc
  integer   ( kind = 4 ) :: ierror = 0
  character ( len = 255 ) :: node_filename = ' '
  integer   ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xyz
  character ( len = 255 ) :: element_filename = ' '
  integer   ( kind = 4 ) tetra_num
  integer   ( kind = 4 ) tetra_num2
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: tetra_node
  integer   ( kind = 4 ), allocatable, dimension ( : ) :: vm

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TABLE_TET_MESH'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Read a real TABLE dataset of N points in 3 dimensions,'
  write ( *, '(a)' ) '  Compute the Delaunay tet mesh.'
  write ( *, '(a)' ) '  Write an integer TABLE dataset of the tet mesh.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  If at least one command line argument, it's the input file name.
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, node_filename )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TABLE_TET_MESH:'
    write ( *, '(a)' ) '  Please enter the name of the input file.'

    read ( *, '(a)' ) node_filename

  end if
!
!  Create the output file name from the input file name.
!
  element_filename = node_filename
  call file_name_ext_swap ( element_filename, 'tetra.txt' )
!
!  Read the point coordinates.
!
  call r8mat_header_read (  node_filename, dim_num, node_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the header of "' // trim ( node_filename ) //'".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points NODE_NUM = ', node_num

  allocate ( node_xyz(1:dim_num,1:node_num) )

  call r8mat_data_read ( node_filename, dim_num, node_num, node_xyz )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the data in "' // trim ( node_filename ) //'".'

  call r8mat_transpose_print_some ( dim_num, node_num, node_xyz, 1, 1, 5, 5, &
    '  5 by 5 portion of node data read from file:' )
!
!  Determine the tet mesh.
!
  ht_num = ( 3 * node_num ) / 2

  allocate ( ht(ht_num) )
  allocate ( vm(1:node_num) )

  do i = 1, node_num
    vm(i) = i
  end do

  call dtris3 ( node_num, ht_num, bf_max, fc_max, node_xyz, vm, bf_num, &
    fc_num, face_num, tetra_num, bf, fc, ht, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TABLE_TET_MESH - Fatal error!'
    write ( *, '(a,i8)' ) '  DTRIS3 returned IERROR = ', ierror
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TABLE_TET_MESH'
    write ( *, '(a)' ) '  ABNORMAL end of execution!'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  BF_MAX = ', bf_max
  write ( *, '(a,i8)' ) '  BF_NUM = ', bf_num
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  FC_MAX = ', fc_max
  write ( *, '(a,i8)' ) '  FC_NUM = ', fc_num
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  HT_NUM = ', ht_num
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  TETRA_NUM = ', tetra_num

  allocate ( tetra_node(1:4,1:tetra_num) )

  call tetlst ( fc_max, fc_num, vm, fc, tetra_num, tetra_num2, tetra_node )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  TETRA_NUM2 = ', tetra_num2
!
!  Print a portion.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed the tet mesh.'
 
  call i4mat_transpose_print_some ( 4, tetra_num, tetra_node, 1, 1, &
    4, 5, '  4 by 5 portion of tetra data:' )
!
!  Write the tet mesh to a file.
!
  call i4mat_write ( element_filename, 4, tetra_num, tetra_node )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Wrote the tetra data to "' &
    // trim ( element_filename ) //'".'
!
!  Free memory.
!
  deallocate ( ht )
  deallocate ( node_xyz )
  deallocate ( tetra_node )
  deallocate ( vm )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TABLE_TET_MESH'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine availf ( hdavfc, fc_num, fc_max, fc, ind, ierr )

!*****************************************************************************80
!
!! AVAILF returns the index of the next available record in the FC array.
!
!  Discussion: 
!
!    This routine returns the index of the next available record in the
!    FC array, either HDAVFC or FC_NUM+1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2005
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) HDAVFC, head pointer of available records in FC.
!
!    Input/output, integer ( kind = 4 ) FC_NUM, current number of records used in FC.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum number of records available in FC.
!
!    Input, integer ( kind = 4 ) FC(1:7,1:*), array of face records; see routine DTRIS3.
!
!    Output, integer ( kind = 4 ) IND, the index of available record (if FC not full).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) fc_num
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) fc_max

  ierr = 0

  if ( hdavfc /= 0 ) then
    ind = hdavfc
    hdavfc = -fc(1,hdavfc)
  else if ( fc_max <= fc_num ) then
    ierr = 11
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AVAILF - Fatal error!'
    write ( *, '(a)' ) '  Memory requirements for array FC exceed the'
    write ( *, '(a,i12)' ) '  current limit of FC_MAX = ', fc_max
  else
    fc_num = fc_num + 1
    ind = fc_num
  end if

  return
end
subroutine baryth ( a, b, c, d, e, alpha, degen )

!*****************************************************************************80
!
!! BARYTH computes barycentric coordinates of a point in 3D.
!
!  Discussion: 
!
!    This routine computes the barycentric coordinates of a 3D point with
!    respect to the four vertices of a tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2009
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(1:3), B(1:3), C(1:3), D(1:3), 4 vertices
!    of tetrahedron.
!
!    Input, real ( kind = 8 ) E(1:3), fifth point for which 
!    barycentric coordinates found
!
!    Output, real ( kind = 8 ) ALPHA(1:4), the scaled barycentric coordinates
!    (if DEGEN = .FALSE.) such that 
!      E = (ALPHA(1)*A + ALPHA(2)*B + ALPHA(3)*C +ALPHA(4)*D)/DET 
!    where DET = 6 * (volume of tetra ABCD);  an ALPHA(I) may be set to 0 
!    after tolerance test to indicate that E is coplanar with a face, so 
!    sum of ALPHA(I)/DET may not be 1; if the actual barycentric
!    coordinates rather than just their signs are needed,
!    modify this routine to divide ALPHA(I) by DET.
!
!    Output, logical DEGEN, TRUE iff A, B, C, D are coplanar.
!
  implicit none

  real ( kind = 8 ) a(3)
  real ( kind = 8 ) alpha(4)
  real ( kind = 8 ) amax
  real ( kind = 8 ) b(3)
  real ( kind = 8 ) bmax
  real ( kind = 8 ) c(3)
  real ( kind = 8 ) cmax
  real ( kind = 8 ) cp1
  real ( kind = 8 ) cp2
  real ( kind = 8 ) cp3
  real ( kind = 8 ) d(3)
  real ( kind = 8 ) da(3)
  real ( kind = 8 ) db(3)
  real ( kind = 8 ) dc(3)
  real ( kind = 8 ) de(3)
  logical degen
  real ( kind = 8 ) det
  real ( kind = 8 ) dmax
  real ( kind = 8 ) e(3)
  real ( kind = 8 ) ea(3)
  real ( kind = 8 ) eb(3)
  real ( kind = 8 ) ec(3)
  real ( kind = 8 ) emax
  real ( kind = 8 ) tol

  tol = 100.0D+00 * epsilon ( tol )
  degen = .false.

  da(1:3) = a(1:3) - d(1:3)
  db(1:3) = b(1:3) - d(1:3)
  dc(1:3) = c(1:3) - d(1:3)

  amax = max ( abs ( a(1) ), abs ( a(2) ), abs ( a(3) ) )
  bmax = max ( abs ( b(1) ), abs ( b(2) ), abs ( b(3) ) )
  cmax = max ( abs ( c(1) ), abs ( c(2) ), abs ( c(3) ) )
  dmax = max ( abs ( d(1) ), abs ( d(2) ), abs ( d(3) ) )

  cp1 = db(2) * dc(3) - db(3) * dc(2)
  cp2 = db(3) * dc(1) - db(1) * dc(3)
  cp3 = db(1) * dc(2) - db(2) * dc(1)
  det = da(1) * cp1 + da(2) * cp2 + da(3) * cp3

  if ( abs ( det ) <= 0.01D+00 * tol * max ( amax, bmax, cmax, dmax ) ) then
    degen = .true. 
    return
  end if

  de(1:3) = e(1:3) - d(1:3)
  ea(1:3) = a(1:3) - e(1:3)
  eb(1:3) = b(1:3) - e(1:3)
  ec(1:3) = c(1:3) - e(1:3)

  alpha(1) = de(1) * cp1 + de(2) * cp2 + de(3) * cp3

  cp1 = da(2) * de(3) - da(3) * de(2)
  cp2 = da(3) * de(1) - da(1) * de(3)
  cp3 = da(1) * de(2) - da(2) * de(1)

  alpha(2) = dc(1) * cp1 + dc(2) * cp2 + dc(3) * cp3
  alpha(3) = db(1) * cp1 + db(2) * cp2 + db(3) * cp3

  alpha(4) = ea(1) * ( eb(2) * ec(3) - eb(3) * ec(2) ) &
           + ea(2) * ( eb(3) * ec(1) - eb(1) * ec(3) ) &
           + ea(3) * ( eb(1) * ec(2) - eb(2) * ec(1) )

  if ( det < 0.0D+00 ) then
    alpha(1) = -alpha(1)
    alpha(2) = -alpha(2)
    alpha(4) = -alpha(4)
  else
    alpha(3) = -alpha(3)
  end if

  emax = max ( abs ( e(1) ), abs ( e(2) ), abs ( e(3) ) )

  if ( abs ( alpha(1) ) <= tol * max ( bmax, cmax, dmax, emax ) ) then
    alpha(1) = 0.0D+00
  end if

  if ( abs ( alpha(2) ) <= tol * max ( amax, cmax, dmax, emax ) ) then
    alpha(2) = 0.0D+00
  end if

  if ( abs ( alpha(3) ) <= tol * max ( amax, bmax, dmax, emax ) ) then
    alpha(3) = 0.0D+00
  end if

  if ( abs ( alpha(4) ) <= tol * max ( amax, bmax, cmax, emax ) ) then
    alpha(4) = 0.0D+00
  end if

  return
end
subroutine ccsph ( intest, a, b, c, d, e, center, radsq, in )

!*****************************************************************************80
!
!! CCSPH finds the circumsphere through the vertices of a tetrahedron.
!
!  Discussion: 
!
!    This routine finds the center and the square of the radius of 
!    the circumsphere through four vertices of a tetrahedron, and 
!    possibly determines whether a fifth 3D point is inside the sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2009
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, logical INTEST, is TRUE, if and only if the test for fifth point 
!    in sphere is to be made.
!
!    Input, real ( kind = 8 ) A(1:3), B(1:3), C(1:3), D(1:3), vertices 
!    of tetrahedron.
!
!    Input, real ( kind = 8 ) E(1:3), a fifth point; referenced if 
!    and only if INTEST is TRUE.
!
!    Output, real ( kind = 8 ) CENTER(1:3), center of sphere; undefined 
!    if A,B,C,D coplanar.
!
!    Output, real ( kind = 8 ) RADSQ, the square of radius of sphere; 
!    -1 if A,B,C,D coplanar.
!
!    Output, integer ( kind = 4 ) IN, contains following value if INTEST is .TRUE.:
!     2 if A,B,C,D coplanar; 
!     1 if E inside sphere;
!     0 if E on sphere; 
!    -1 if E outside sphere
!
  implicit none

  real ( kind = 8 ) a(3)
  real ( kind = 8 ) b(3)
  real ( kind = 8 ) c(3)
  real ( kind = 8 ) center(3)
  real ( kind = 8 ) cmax
  real ( kind = 8 ) cp1
  real ( kind = 8 ) cp2
  real ( kind = 8 ) cp3
  real ( kind = 8 ) d(3)
  real ( kind = 8 ) da(3)
  real ( kind = 8 ) db(3)
  real ( kind = 8 ) dc(3)
  real ( kind = 8 ) det
  real ( kind = 8 ) dsq
  real ( kind = 8 ) e(3)
  integer ( kind = 4 ) in
  logical intest
  real ( kind = 8 ) radsq
  real ( kind = 8 ) rhs(3)
  real ( kind = 8 ) tol

  tol = 100.0D+00 * epsilon ( tol )

  da(1:3) = a(1:3) - d(1:3)
  db(1:3) = b(1:3) - d(1:3)
  dc(1:3) = c(1:3) - d(1:3)

  rhs(1) = 0.5D+00 * sum ( da(1:3)**2 )
  rhs(2) = 0.5D+00 * sum ( db(1:3)**2 )
  rhs(3) = 0.5D+00 * sum ( dc(1:3)**2 )

  cmax = max ( &
    abs ( a(1) ), abs ( a(2) ), abs ( a(3) ), &
    abs ( b(1) ), abs ( b(2) ), abs ( b(3) ), &
    abs ( c(1) ), abs ( c(2) ), abs ( c(3) ), &
    abs ( d(1) ), abs ( d(2) ), abs ( d(3) ) )

  cp1 = db(2) * dc(3) - dc(2) * db(3)
  cp2 = dc(2) * da(3) - da(2) * dc(3)
  cp3 = da(2) * db(3) - db(2) * da(3)

  det = da(1) * cp1 + db(1) * cp2 + dc(1) * cp3

  if ( abs ( det ) <= 0.01D+00 * tol * cmax ) then
    radsq = -1.0D+00
    in = 2
    return
  end if

  center(1) = ( rhs(1) * cp1 + rhs(2) * cp2 + rhs(3) * cp3 ) / det

  cp1 = db(1) * rhs(3) - dc(1) * rhs(2)
  cp2 = dc(1) * rhs(1) - da(1) * rhs(3)
  cp3 = da(1) * rhs(2) - db(1) * rhs(1)

  center(2) =  ( da(3) * cp1 + db(3) * cp2 + dc(3) * cp3 ) / det
  center(3) = -( da(2) * cp1 + db(2) * cp2 + dc(2) * cp3 ) / det

  radsq = sum ( center(1:3)**2 )

  center(1:3) = center(1:3) + d(1:3)

  if ( intest ) then

    dsq = sum ( ( e(1:3) - center(1:3) )**2 )

    if ( ( 1.0D+00 + tol ) * radsq < dsq ) then
      in = -1
    else if ( dsq < ( 1.0D+00 - tol ) * radsq ) then
      in = 1
    else
      in = 0
    end if

  end if

  return
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  logical ch_eqi
  character c1
  character c1_cap
  character c2
  character c2_cap

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine dhpsrt ( k, n, lda, a, map )

!*****************************************************************************80
!
!! DHPSRT sorts a list of double precision points in KD.
!
!  Discussion: 
!
!    This routine uses heapsort to obtain the permutation of N K-dimensional
!    double precision points so that the points are in lexicographic
!    increasing order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 August 2005
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the dimension of points.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling 
!    routine; K <= LDA.
!
!    Input, real ( kind = 8 ) A(1:K,1:*), array of points.
!
!    Input/output, integer ( kind = 4 ) MAP(1:N).  On input, he points of A with indices 
!    MAP(1), MAP(2), ..., MAP(N) are to be sorted.  On output, the elements 
!    are permuted so that A(*,MAP(1)) <= A(*,MAP(2)) <= ... <= A(*,MAP(N)).
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) map(n)
  integer ( kind = 4 ) t

  do i = n/2, 1, -1
    call dsftdw ( i, n, k, lda, a, map )
  end do

  do i = n, 2, -1
    t = map(1)
    map(1) = map(i)
    map(i) = t
    call dsftdw ( 1, i-1, k, lda, a, map )
  end do

  return
end
function dless ( k, p, q )

!*****************************************************************************80
!
!! DLESS determines the lexicographically lesser of two double precision values.
!
!  Discussion: 
!
!    This routine determines whether P is lexicographically less than Q in
!    floating point arithmetic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2009
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, dimension of points.
!
!    Input, real ( kind = 8 ) P(1:K), Q(1:K), two points.
!
!    Output, logical DLESS, TRUE if P < Q, FALSE otherwise.
!
  implicit none

  real ( kind = 8 ) cmax
  logical dless
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind = 8 ) p(k)
  real ( kind = 8 ) q(k)
  real ( kind = 8 ) tol

  tol = 100.0D+00 * epsilon ( tol )

  do i = 1, k

    cmax = max ( abs ( p(i) ), abs ( q(i) ) )

    if ( abs ( p(i) - q(i) ) <= tol * cmax .or. cmax <= tol ) then
      cycle
    end if

     if ( p(i) < q(i) ) then
       dless = .true.
     else
       dless = .false.
     end if

     return

  end do

  dless = .false.

  return
end
subroutine dsftdw ( l, u, k, lda, a, map )

!*****************************************************************************80
!
!! DSFTDW does one step of the heap sort algorithm for double precision data.
!
!  Discussion: 
!
!    This routine sifts A(*,MAP(L)) down a heap of size U.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2009
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, U, the lower and upper index of part of heap.
!
!    Input, integer ( kind = 4 ) K, the dimension of points.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling routine.
!
!    Input, real ( kind = 8 ) A(1:K,1:*), see routine DHPSRT.
!
!    Input/output, integer ( kind = 4 ) MAP(1:*), see routine DHPSRT.
!
  implicit none

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) map(*)
  integer ( kind = 4 ) u

  real ( kind = 8 ) a(lda,*)
  logical dless
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) t

  i = l
  j = 2 * i
  t = map(i)

  do

    if ( u < j ) then
      exit
    end if

    if ( j < u ) then
      if ( dless ( k, a(1,map(j)), a(1,map(j+1)) ) ) then
        j = j + 1
      end if
    end if

    if ( dless ( k, a(1,map(j)), a(1,t)) ) then
      exit
    end if

    map(i) = map(j)
    i = j
    j = 2 * i

  end do

  map(i) = t

  return
end
subroutine dtris3 ( npt, sizht, bf_max, fc_max, vcl, vm, bf_num, fc_num, nface, &
  ntetra, bf, fc, ht, ierr )

!*****************************************************************************80
!
!! DTRIS3 constructs a Delaunay triangulation of vertices in 3D.
!
!  Discussion:
!
!    This routine constructs a Delaunay triangulation of 3D vertices using
!    an incremental approach and local transformations.  Vertices are
!    first sorted in lexicographically increasing (x,y,z) order, and
!    then are inserted one at a time from outside the convex hull.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 September 2005
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPT, the number of 3D vertices.
!
!    Input, integer ( kind = 4 ) SIZHT, the size of the hash table HT; a good choice is
!    a prime number which is about 1/8 * NFACE (or 3/2 * NPT for random
!    points from the uniform distribution).
!
!    Input, integer ( kind = 4 ) BF_MAX, the maximum size available for BF array.
!    This needs to be at least as big as the number of boundary faces.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!    This needs to be at least as big as the number of faces.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:NPT), the vertex coordinates.
!    In the general case, VCL may contain the coordinates for more
!    than NPT vertices, and the VM array is used to select them.
!
!    Input/output, integer ( kind = 4 ) VM(1:NPT), the vertices of VCL to be triangulated.
!    On output, these indices are permuted, so that VCL(*,VM(1)), ... ,
!    VCL(*,VM(NPT)) are in lexicographic increasing order,
!    with possible slight reordering so first 4 vertices are
!    non-coplanar.  Typically, the input value of VM might be 1 through
!    NPT.
!
!    Output, integer ( kind = 4 ) BF_NUM, the number of positions used in BF array;
!    BF_NUM <= BF_MAX.
!
!    Output, integer ( kind = 4 ) FC_NUM, the number of positions used in FC array;
!    FC_NUM <= FC_MAX.
!
!    Output, integer ( kind = 4 ) NFACE, the number of faces in triangulation; 
!    NFACE <= FC_NUM.
!
!    Output, integer ( kind = 4 ) NTETRA, the number of tetrahedra in the triangulation.
!
!    Output, integer ( kind = 4 ) BF(1:3,1:BF_NUM), boundary face records containing pointers
!    (indices) to FC; if FC(5,I) = -J < 0 and FC(1:3,I) = ABC,
!    then BF(1,J) points to other boundary face with edge BC,
!    BF(2,J) points to other boundary face with edge AC, and
!    BF(3,J) points to other boundary face with edge AB;
!    if BF(1,J) <= 0, record is not used and is in avail list.
!
!    Output, integer ( kind = 4 ) FC(1:7,1:FC_NUM), face records which are in linked lists
!    in hash table with direct chaining. Fields are:
!    FC(1:3,*) - A,B,C with 1<=A<B<C<=NPT; indices in VM of 3
!    vertices of face; if A <= 0, record is not used (it is
!    in linked list of available records with indices <= FC_NUM);
!    internal use: if B <= 0, face in queue, not in triangulation.
!    FC(4:5,*) - D,E; indices in VM of 4th vertex of 1 or 2
!    tetrahedra with face ABC; if ABC is boundary face
!    then E < 0 and |E| is an index of BF array
!    FC(6,*) - HTLINK; pointer (index in FC) of next element
!    in linked list (or NULL = 0)
!    FC(7,*) - used internally for QLINK (link for queues or
!    stacks); pointer (index in FC) of next face in queue/
!    stack (or NULL = 0); QLINK = -1 indicates face is not
!    in any queue/stack, and is output value (for records
!    not in avail list), except:
!    FC(7,1:2) - HDAVBF,HDAVFC : head pointers of avail list in BF, FC.
!
!    Output, integer ( kind = 4 ) HT(0:SIZHT-1), a hash table using direct chaining;
!    entries are head pointers of linked lists (indices of FC array)
!    containing the faces and tetrahedra of the triangulation.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) bf_max
  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(3,bf_max)
  integer ( kind = 4 ) bfi
  real ( kind = 8 ) ctr(3)
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) fc_num
  integer ( kind = 4 ) front
  integer ( kind = 4 ) hdavbf
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) ntetra
  integer ( kind = 4 ) op
  integer ( kind = 4 ) opside
  integer ( kind = 4 ) ptr
  integer ( kind = 4 ) topnv
  integer ( kind = 4 ) top
  integer ( kind = 4 ) va
  integer ( kind = 4 ) vb
  integer ( kind = 4 ) vc
  real ( kind = 8 ) vcl(3,*)
  integer ( kind = 4 ) vi
  integer ( kind = 4 ) vm(npt)

  ierr = 0
!
!  Permute elements of VM so that vertices are in lexicographic order.
!
  call dhpsrt ( 3, npt, 3, vcl, vm )
!
!  Reorder points so that first four points are in general position.
!
  call frstet ( .true., npt, vcl, vm, i3, i4, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTRIS3 - Error!'
    write ( *, '(a,i8)' ) '  FRSTET returned IERR = ', ierr
    return
  end if
!
!  Initialize data structures.
!
  do i = 1, 3
    ctr(i) = sum ( vcl(i,vm(1:4)) ) / 4.0D+00
  end do

  ht(0:sizht-1) = 0
  hdavbf = 0
  hdavfc = 0
  bf_num = 4
  fc_num = 4
  ntetra = 1

  call htins ( 1, 1, 2, 3, 4, -1, npt, sizht, fc, ht )
  call htins ( 2, 1, 2, 4, 3, -2, npt, sizht, fc, ht )
  call htins ( 3, 1, 3, 4, 2, -3, npt, sizht, fc, ht )
  call htins ( 4, 2, 3, 4, 1, -4, npt, sizht, fc, ht )

  bf(1:3,1) = (/ 4, 3, 2 /)
  bf(1:3,2) = (/ 4, 3, 1 /)
  bf(1:3,3) = (/ 4, 2, 1 /)
  bf(1:3,4) = (/ 3, 2, 1 /)

  if ( msglvl == 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTRIS3:'
    write ( *, '(a)' ) '  First tetrahedron:'
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) vm(1:4)
    write ( *, '(a,i8,a,i8)' ) '  I3 = ', i3, '  I4 = ', i4
  end if
!
!  Insert the I-th vertex into Delaunay triangle of first I-1 vertices.
!
  do i = 5, npt

    vi = vm(i)

    if ( msglvl == 4 ) then
      write ( *, '(a,i8,a,i8)' ) '  Step: ', i, '  Vertex: ', vi
    end if

    if ( i == 5 ) then
      ip = 2
    else
      ip = i - 1
    end if

    if ( i == i3 + 2 ) then
      ip = 3
    end if

    if ( i == i4 + 1 ) then
      ip = 4
    end if
!
!  Form stacks of boundary faces involving vertex IP.
!  TOP is for stack of boundary faces to be tested for visibility.
!  FRONT is for stack of boundary faces visible from vertex I.
!  TOPNV is for stack of boundary faces not visible from I.
!
    front = 0
    topnv = 0

    if ( i == 5 ) then

      top = 4

      if ( ip == 2 ) then
        a = 2
      else
        a = 3
      end if

      if ( ip <= 3 ) then
        b = 1
      else
        b = 2
      end if

      fc(7,top) = a
      fc(7,a) = b
      fc(7,b) = 0

    else if ( ip == i - 1 ) then

      top = bfi
      fc(7,bfi) = 0
      b = fc(2,bfi)
      ptr = bf(1,-fc(5,bfi))

      do

        if ( fc(1,ptr) == b ) then
          b = fc(2,ptr)
          j = 1
        else
          b = fc(1,ptr)
          j = 2
        end if

        fc(7,ptr) = top
        top = ptr
        ptr = bf(j,-fc(5,ptr))

        if ( ptr == bfi ) then
          exit
        end if

      end do

    else

      j = 0

      do k = 1, bf_num

        if ( bf(1,k) <= 0 ) then
          cycle
        end if

        do e = 1, 3

          ptr = bf(e,k)

          if ( fc(1,ptr) == ip ) then
            b = fc(2,ptr)
            j = 3
            exit
          else if ( fc(2,ptr) == ip ) then
            b = fc(1,ptr)
            j = 3
            exit
          else if ( fc(3,ptr) == ip ) then
            b = fc(1,ptr)
            j = 2
            exit
          end if

        end do

        if ( j /= 0 ) then
          exit
        end if

      end do

      bfi = ptr
      top = bfi
      fc(7,bfi) = 0
      ptr = bf(j,-fc(5,bfi))

      do

        if ( fc(1,ptr) == b ) then
          j = 1
          if ( fc(2,ptr) == ip ) then
            b = fc(3,ptr)
          else
            b = fc(2,ptr)
          end if
        else if ( fc(2,ptr) == b ) then
          j = 2
          if ( fc(1,ptr) == ip ) then
            b = fc(3,ptr)
          else
            b = fc(1,ptr)
          end if
        else
          j = 3
          if ( fc(1,ptr) == ip ) then
            b = fc(2,ptr)
          else
            b = fc(1,ptr)
          end if
        end if

        fc(7,ptr) = top
        top = ptr
        ptr = bf(j,-fc(5,ptr))

        if ( ptr == bfi ) then
          exit
        end if

      end do

    end if
!
!  Find a boundary face visible from vertex I.
!
    do while ( top /= 0 )

      ptr = top
      top = fc(7,ptr)
      va = vm(fc(1,ptr))
      vb = vm(fc(2,ptr))
      vc = vm(fc(3,ptr))
      op = opside ( vcl(1,va), vcl(1,vb), vcl(1,vc), ctr, vcl(1,vi) )

      if ( op == 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DTRIS3 - Error!'
        write ( *, '(a)' ) '  Unexpected return value from OPSIDE.'
        ierr = 301
        return
      end if

      if ( op == 1 ) then

        front = ptr

        do while ( top /= 0 )

          ptr = top
          top = fc(7,ptr)
          fc(7,ptr) = -1

        end do

      else

        fc(7,ptr) = topnv
        topnv = ptr

      end if

    end do

    if ( front == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS3 - Error!'
      write ( *, '(a)' ) '  FRONT = 0.'
      ierr = 306
      return
    end if
!
!  Find remaining visible boundary faces, add new tetrahedra with
!  vertex I, apply local transformation based on empty sphere criterion.
!
    call vbfac ( vcl(1,vi), ctr, vcl, vm, bf, fc, front, topnv )

    call nwthou ( i, npt, sizht, bf_num, fc_num, bf_max, fc_max, bf, fc, ht, &
      ntetra, hdavbf, hdavfc, front, back, bfi, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS3 - Error!'
      write ( *, '(a,i8)' ) '  NWTHOU returned IERR = ', ierr
      return
    end if

    call swapes ( .false., i, npt, sizht, fc_num, fc_max, vcl, vm, bf, fc, ht, &
      ntetra, hdavfc, front, back, j, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS3 - Error!'
      write ( *, '(a,i8)' ) '  SWAPES returned IERR = ', ierr
      return
    end if

  end do

  nface = fc_num
  ptr = hdavfc

  do while ( ptr /= 0 )
    nface = nface - 1
    ptr = -fc(1,ptr)
  end do

  fc(7,1) = hdavbf
  fc(7,2) = hdavfc

  return
end
subroutine file_column_count ( input_file_name, column_num )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of COLUMN_NUM words,
!    separated by spaces.  There may also be some blank lines, and some
!    comment lines,
!    which have a "#" in column 1.
!
!    The routine tries to find the first non-comment non-blank line and
!    counts the number of words in that line.
!
!    If all lines are blanks or comments, it goes back and tries to analyze
!    a comment line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer ( kind = 4 ) column_num
  logical got_one
  character ( len = * ) input_file_name
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    ' // trim ( input_file_name )
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    got_one = .true.
    exit

  end do

  if ( .not. got_one ) then

    rewind ( input_unit )

    do

      read ( input_unit, '(a)', iostat = ios ) line

      if ( ios /= 0 ) then
        exit
      end if

      if ( len_trim ( line ) == 0 ) then
        cycle
      end if

      got_one = .true.
      exit

    end do

  end if

  close ( unit = input_unit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    column_num = -1
    return
  end if

  call s_word_count ( line, column_num )

  return
end
subroutine file_name_ext_get ( file_name, i, j )

!*****************************************************************************80
!
!! FILE_NAME_EXT_GET determines the "extension" of a file name.
!
!  Discussion:
!
!    The "extension" of a file name is the string of characters
!    that appears after the LAST period in the name.  A file
!    with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!    Blanks are unusual in file names.  This routine ignores all
!    trailing blanks, but will treat initial or internal blanks
!    as regular characters acceptable in a file name.
!
!  Example:
!
!    FILE_NAME   I  J
!
!    bob.for     4  7
!    N.B.C.D     6  7
!    Naomi.      6  6
!    Arthur      0  0
!    .com        1  1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, a file name to be examined.
!
!    Output, integer ( kind = 4 ) I, J, the indices of the first and last characters
!    in the file extension.
!
!    If no period occurs in FILE_NAME, then
!      I = J = 0;
!    Otherwise,
!      I is the position of the LAST period in FILE_NAME, and J is the
!      position of the last nonblank character following the period.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s_index_last

  i = s_index_last ( file_name, '.' )

  if ( i /= 0 ) then

    j = len_trim ( file_name )

  else

    j = 0

  end if

  return
end
subroutine file_name_ext_swap ( file_name, ext )

!*****************************************************************************80
!
!! FILE_NAME_EXT_SWAP replaces the current "extension" of a file name.
!
!  Discussion:
!
!    The "extension" of a file name is the string of characters
!    that appears after the LAST period in the name.  A file
!    with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!  Example:
!
!          Input           Output
!    ================     =========
!    FILE_NAME    EXT     FILE_NAME
!
!    bob.for      obj     bob.obj
!    bob.bob.bob  txt     bob.bob.txt
!    bob          yak     bob.yak
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME, a file name.
!    On output, the extension of the file has been changed.
!
!    Input, character ( len = * ) EXT, the extension to be used on the output
!    copy of FILE_NAME, replacing the current extension if any.
!
  implicit none

  character ( len = * ) ext
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) len_max
  integer ( kind = 4 ) len_name

  len_max = len ( file_name )
  len_name = len_trim ( file_name )

  call file_name_ext_get ( file_name, i, j )

  if ( i == 0 ) then

    if ( len_max < len_name + 1 ) then
      return
    end if

    len_name = len_name + 1
    file_name(len_name:len_name) = '.'
    i = len_name + 1

  else

    i = i + 1
    file_name(i:j) = ' '

  end if

  file_name(i:) = ext

  return
end
subroutine file_row_count ( input_file_name, row_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of row records in a file.
!
!  Discussion:
!
!    It does not count lines that are blank, or that begin with a
!    comment symbol '#'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
!
  implicit none

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) comment_num
  integer ( kind = 4 ) ierror
  character ( len = * ) input_file_name
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  integer ( kind = 4 ) record_num
  integer ( kind = 4 ) row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( input_file_name )
    stop
  end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = record_num
      exit
    end if

    record_num = record_num + 1

    if ( line(1:1) == '#' ) then
      comment_num = comment_num + 1
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      comment_num = comment_num + 1
      cycle
    end if

    row_num = row_num + 1

  end do

  close ( unit = input_unit )

  return
end
subroutine frstet ( shift, nv, vcl, map, i3, i4, ierr )

!*****************************************************************************80
!
!! FRSTET shifts vertices so the first 4 vertices are in general position in 3D.
!
!  Discussion: 
!
!    This routine shifts or swaps vertices if necessary so first 3 vertices
!    (according to MAP) are not collinear and first 4 vertices are
!    not coplanar (so that first tetrahedron is valid).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2009
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms, 
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, logical SHIFT, if TRUE, MAP(3), MAP(4) may be updated due to shift,
!    else they may be updated due to swaps; in former case,
!    it is assumed MAP gives vertices in lexicographic order.
!
!    Input, integer ( kind = 4 ) NV, the number of vertices.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) MAP(1:NV), on input, contains vertex indices of VCL.
!    On output, shifted or 2 swaps applied if necessary so that vertices
!    indexed by MAP(1), MAP(2), MAP(3), MAP(4) not coplanar.
!
!    Output, integer ( kind = 4 ) I3, I4, the indices such that MAP_in(I3) = MAP_out(3) and
!    MAP_in(I4) = MAP_out(4).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) nv

  real ( kind = 8 ) cmax
  real ( kind = 8 ) cp1
  real ( kind = 8 ) cp2
  real ( kind = 8 ) cp3
  real ( kind = 8 ) dmax
  real ( kind = 8 ) dotp
  real ( kind = 8 ) dv2(3)
  real ( kind = 8 ) dvk(3)
  real ( kind = 8 ) dvl(3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) map(nv)
  logical shift
  real ( kind = 8 ) tol
  real ( kind = 8 ) vcl(3,*)
!
!  First check that consecutive vertices are not identical.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )

  if ( shift ) then
    l = nv - 1
  else
    l = 1
  end if

  m1 = map(1)

  do i = 1, l

    m = m1
    m1 = map(i+1)

    do k = 1, 3
      cmax = max ( abs ( vcl(k,m) ), abs ( vcl(k,m1) ) )
      if ( tol * cmax < abs ( vcl(k,m) - vcl(k,m1) ) .and. tol < cmax ) then
        go to 20
      end if
    end do

    ierr = 302
    return

20  continue

  end do
!
!  Find index K = I3 and L = I4.
!
  m1 = map(1)
  m2 = map(2)

  dv2(1:3) = vcl(1:3,m2) - vcl(1:3,m1)

  cmax = max ( abs ( vcl(1,m1) ), abs ( vcl(2,m1) ), abs ( vcl(3,m1) ), &
    abs ( vcl(1,m2) ), abs ( vcl(2,m2) ), abs ( vcl(3,m2) ) )
  k = 2

  do

    k = k + 1

    if ( nv < k ) then
      ierr = 303
      return
    end if

    m = map(k)

    dvk(1:3) = vcl(1:3,m) - vcl(1:3,m1)

    dmax = max ( cmax, abs ( vcl(1,m) ), abs ( vcl(2,m) ), abs ( vcl(3,m) ) )

    cp1 = dv2(2) * dvk(3) - dv2(3) * dvk(2)
    cp2 = dv2(3) * dvk(1) - dv2(1) * dvk(3)
    cp3 = dv2(1) * dvk(2) - dv2(2) * dvk(1)

    if ( tol * dmax < max ( abs ( cp1 ), abs ( cp2 ), abs ( cp3 ) ) ) then
      exit
    end if

  end do

  cmax = dmax
  l = k

  do

    l = l + 1

    if ( nv < l ) then
      ierr = 304
      return
    end if

    m = map(l)

    dvl(1:3) = vcl(1:3,m) - vcl(1:3,m1)

    dmax = max ( cmax, abs ( vcl(1,m) ), abs ( vcl(2,m) ), abs ( vcl(3,m) ) )

    dotp = dvl(1) * cp1 + dvl(2) * cp2 + dvl(3) * cp3

    if ( tol * dmax < abs ( dotp ) ) then
      exit
    end if

  end do
!
!  Shift or swap elements of MAP if necessary.
!
  if ( shift ) then

    if ( 3 < k ) then
      m1 = map(k)
    end if

    if ( 4 < l ) then
      m2 = map(l)
      do i = l, k+2, -1
         map(i) = map(i-1)
      end do
      do i = k+1, 5, -1
        map(i) = map(i-2)
      end do
      map(4) = m2
    end if

    if ( 3 < k ) then
      map(3) = m1
    end if

  else

    if ( 3 < k ) then
      m = map(3)
      map(3) = map(k)
      map(k) = m
    end if

    if ( 4 < l ) then
      m = map(4)
      map(4) = map(l)
      map(l) = m
    end if

  end if

  i3 = k
  i4 = l

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine htdel ( ind, n, p, fc, ht )

!*****************************************************************************80
!
!! HTDEL deletes a record from the hash table.
!
!  Discussion: 
!
!    This routine deletes record FC(1:7,IND) from the hash table HT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2009
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IND, the index of FC array.
!
!    Input, integer ( kind = 4 ) N, the upper bound on vertex indices.
!
!    Input, integer ( kind = 4 ) P, the size of hash table.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:*), the array of face records; 
!    see routine DTRIS3.  On output, one link in FC is updated.
!
!    Input/output, integer ( kind = 4 ) HT(0:P-1), the hash table using direct chaining.
!    On output, one link in HT is updated.
!
  implicit none

  integer ( kind = 4 ) p

  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ptr

  k = mod ( fc(1,ind) * n + fc(2,ind), p )
  k = mod ( k * n + fc(3,ind), p )
  ptr = ht(k)

  if ( ptr == ind ) then

    ht(k) = fc(6,ind)

  else

    do while ( fc(6,ptr) /= ind )
      ptr = fc(6,ptr)
    end do
    fc(6,ptr) = fc(6,ind)

  end if

  return
end
subroutine htins ( ind, a, b, c, d, e, n, p, fc, ht )

!*****************************************************************************80
!
!! HTINS inserts a record into the hash table.
!
!  Discussion: 
!
!    This routine inserts record FC(1:7,IND) containing A,B,C,D,E,HTLINK,-1
!    into hash table HT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2009
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IND, the index of FC array.
!
!    Input, integer ( kind = 4 ) A, B, C, D, E, the first 5 fields of FC record (or column).
!
!    Input, integer ( kind = 4 ) N, the upper bound on vertex indices.
!
!    Input, integer ( kind = 4 ) P, the size of hash table.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:*), the array of face records; 
!    see routine DTRIS3.
!
!    Input/output, integer ( kind = 4 ) HT(0:P-1), the hash table using direct chaining
!
  implicit none

  integer ( kind = 4 ) p

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bb
  integer ( kind = 4 ) c
  integer ( kind = 4 ) cc
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  aa = a
  bb = b
  cc = c
  call order3 ( aa, bb, cc )
  k = mod ( aa * n + bb, p )
  k = mod ( k * n + cc, p )
  fc(1,ind) = aa
  fc(2,ind) = bb
  fc(3,ind) = cc
  fc(4,ind) = d
  fc(5,ind) = e
  fc(6,ind) = ht(k)
  fc(7,ind) = -1
  ht(k) = ind

  return
end
function htsrc ( a, b, c, n, p, fc, ht )

!*****************************************************************************80
!
!! HTSRC searches for a record in the hash table.
!
!  Discussion: 
!
!    This routine searches for record FC(1:7,IND) containing key A,B,C
!    in hash table HT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2009
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A,B,C, first 3 fields of FC record (in any order).
!
!    Input, integer ( kind = 4 ) N, upper bound on vertex indices.
!
!    Input, integer ( kind = 4 ) P, size of hash table.
!
!    Input, integer ( kind = 4 ) FC(1:7,1:*), array of face records; see routine DTRIS3.
!
!    Input, integer ( kind = 4 ) HT(0:P-1), hash table using direct chaining.
!
!    Output, integer ( kind = 4 ) HTSRC, index of FC record with key A,B,C if found,
!    or 0 if not found.
!
  implicit none

  integer ( kind = 4 ) p

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bb
  integer ( kind = 4 ) c
  integer ( kind = 4 ) cc
  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  aa = a
  bb = b
  cc = c
  call order3 ( aa, bb, cc )
  k = mod ( aa * n + bb, p )
  k = mod ( k * n + cc, p )
  ind = ht(k)

  do

    if ( ind == 0 ) then
      exit
    end if

    if ( fc(1,ind) == aa .and. fc(2,ind) == bb .and. fc(3,ind) == cc ) then
      exit
    end if

    ind = fc(6,ind)

  end do

  htsrc = ind

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of integer division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!  Example:
!
!        I     J     MOD  I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an integer to lie between given limits by wrapping.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I4_WRAP
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, an integer value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds for the integer value.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    i4_wrap = jlo
  else
    i4_wrap = jlo + i4_modp ( ival - jlo, wide )
  end if

  return
end
subroutine i4mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = * ) title

  call i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 10
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = 7 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i7)') i
    end do

    write ( *, '(''  Row '',10a7)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i7)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a7)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine i4mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! I4MAT_WRITE writes an I4MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) j
  character ( len = * )  output_filename
  integer   ( kind = 4 ) output_status
  integer   ( kind = 4 ) output_unit
  character ( len = 30 ) string
  integer   ( kind = 4 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
  write ( string, '(a1,i8,a4)' ) '(', m, 'i10)'
!
!  Write the data.
!
  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine nwthou ( i, npt, sizht, bf_num, fc_num, bf_max, fc_max, bf, fc, ht, &
  ntetra, hdavbf, hdavfc, front, back, bfi, ierr )

!*****************************************************************************80
!
!! NWTHOU creates new tetrahedra outside the current convex hull.
!
!  Discussion: 
!
!    This routine creates new tetrahedra in a 3D triangulation outside the
!    convex hull by joining vertex I to visible boundary faces.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2009
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the (local) index of next vertex inserted in
!    triangulation.
!
!    Input, integer ( kind = 4 ) NPT, the number of 3D points to be triangulated.
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT.
!
!    Input/output, integer ( kind = 4 ) BF_NUM, the number of positions used in BF array.
!
!    Input/output, integer ( kind = 4 ) FC_NUM, the number of positions used in FC array.
!
!    Input, integer ( kind = 4 ) BF_MAX, the maximum size available for BF array.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input/output, integer ( kind = 4 ) BF(1:3,1:BF_MAX), the array of boundary face 
!    records; see DTRIS3.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:FC_MAX), the array of face records;
!    see routine DTRIS3.
!
!    Input/output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining.
!
!    Input/output, integer ( kind = 4 ) NTETRA, the number of tetrahedra in triangulation.
!
!    Input/output, integer ( kind = 4 ) HDAVBF, the head pointer to available BF records.
!
!    Input/output, integer ( kind = 4 ) HDAVFC, the head pointer to available FC records.
!
!    Input, integer ( kind = 4 ) FRONT, the index of front of queue (or top of stack) 
!    of visible boundary faces.
!
!    Output, integer ( kind = 4 ) BACK, the index of back of queue (or bottom of stack)
!    of visible boundary faces (which become interior faces).
!
!    Output, integer ( kind = 4 ) BFI, the index of FC of a boundary face containing 
!    vertex I.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) bf_max
  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(3,bf_max)
  integer ( kind = 4 ) bfi
  integer ( kind = 4 ) bfnew
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) fc_num
  integer ( kind = 4 ) front
  integer ( kind = 4 ) hdavbf
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) ntetra
  integer ( kind = 4 ) ptr
!
!  For ABC in queue, form tetrahedron ABCI + add faces ABI, ACI, BCI.
!  PTR, NBR, IND are indices of FC; K, L, BFNEW indices of BF.
!
  ierr = 0
  bfi = 0
  ptr = front

  do

    back = ptr
    a = fc(1,ptr)
    b = fc(2,ptr)
    c = fc(3,ptr)
    k = -fc(5,ptr)
    fc(5,ptr) = i
    ntetra = ntetra + 1

    if ( msglvl == 4 ) then
      write ( *,600) a,b,c,i
    end if

    do e = 1, 3

      if ( e == 2 ) then
        call i4_swap ( a, b )
      else if ( e == 3 ) then
        call i4_swap ( a, c )
      end if

      nbr = bf(e,k)

      if ( fc(7,nbr) /= -1 ) then
        if ( fc(5,nbr) == i ) then
          cycle
        end if
      end if

      call availf ( hdavfc, fc_num, fc_max, fc, ind, ierr )

      if ( ierr /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'NWTHOU - Error!'
        write ( *, '(a,i8)' ) '  AVAILF returned IERR = ', ierr
        return
      end if

      l = -fc(5,nbr)

      if ( bf(1,l) == ptr ) then
        j = 1
      else if ( bf(2,l) == ptr ) then
        j = 2
      else
        j = 3
      end if

      if ( fc(7,nbr) /= -1 ) then

        call htins ( ind, b, c, i, a, fc(j,nbr), npt, sizht, fc, ht )

      else

        if ( hdavbf /= 0 ) then

          bfnew = hdavbf
          hdavbf = -bf(1,hdavbf)

        else

          if ( bf_max <= bf_num ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'NWTHOU - Error!'
            write ( *, '(a)' ) '  BF_MAX <= BF_NUM.'
            write ( *, '(a)' ) '  BF_MAX must be increased to proceed.'
            ierr = 12
            return
          end if

          bf_num = bf_num + 1
          bfnew = bf_num

        end if

        if ( bfi == 0 ) then
          bfi = ind
        end if

        call htins ( ind, b, c, i, a, -bfnew, npt, sizht, fc, ht )
        bf(j,l) = ind
        bf(3,bfnew) = nbr

      end if

    end do

    if ( k == bf_num ) then
      bf_num = bf_num - 1
    else
      bf(1,k) = -hdavbf
      hdavbf = k
    end if

    ptr = fc(7,ptr)

    if ( ptr == 0 ) then
      exit
    end if

  end do
!
!  Set BF(1:2,BFNEW) fields for new boundary faces.
!
  ptr = bfi
  a = fc(1,ptr)
  j = 2

  do

    b = fc(j,ptr)
    c = fc(4,ptr)

    do

      nbr = htsrc ( a, c, i, npt, sizht, fc, ht )
 
      if ( nbr <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'NWTHOU - Error!'
        write ( *, '(a,i8)' ) '  HTSRC returned NBR <= 0.'
        ierr = 300
        return
      end if

      if ( fc(5,nbr) <= 0 ) then
        exit
      end if

      if ( fc(4,nbr) == b ) then
        d = fc(5,nbr)
      else
        d = fc(4,nbr)
      end if

      b = c
      c = d

    end do

    k = -fc(5,ptr)
    l = -fc(5,nbr)

    if ( fc(1,ptr) == a ) then
      bf(2,k) = nbr
    else
      bf(1,k) = nbr
    end if

    if ( fc(1,nbr) == a ) then
      j = 1
    else
      j = 2
    end if

    bf(3-j,l) = ptr
    a = fc(3-j,nbr)
    ptr = nbr

    if ( ptr == bfi ) then
      exit
    end if

  end do

  600 format ( '  New tetra: ',4i7)

  return
end
function opside ( a, b, c, d, e )

!*****************************************************************************80
!
!! OPSIDE tests if points are on opposite sides of a triangular face.
!
!  Discussion: 
!
!    This routine tests if points D, E are on opposite sides of triangular
!    face with vertices A, B, C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 August 2005
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(1:3), B(1:3), C(1:3), D(1:3), E(1:3),
!    five 3D points.
!
!    Output, integer ( kind = 4 ) OPSIDE, the result of the test:
!    +1 if D, E on opposite sides; 
!    -1 if on same side;
!     2 if D is coplanar with face ABC (ABCD is a degenerate tetrahedron); 
!     0 if E is coplanar with face ABC
!
  implicit none

  real ( kind = 8 ) a(3)
  real ( kind = 8 ) ab(3)
  real ( kind = 8 ) ac(3)
  real ( kind = 8 ) b(3)
  real ( kind = 8 ) c(3)
  real ( kind = 8 ) d(3)
  real ( kind = 8 ) ddp
  real ( kind = 8 ) dmax
  real ( kind = 8 ) e(3)
  real ( kind = 8 ) edp
  real ( kind = 8 ) emax
  real ( kind = 8 ) nrml1
  real ( kind = 8 ) nrml2
  real ( kind = 8 ) nrml3
  integer ( kind = 4 ) opside
  real ( kind = 8 ) tol

  tol = 100.0D+00 * epsilon ( tol )

  ab(1:3) = b(1:3) - a(1:3)
  ac(1:3) = c(1:3) - a(1:3)

  emax = max ( &
    abs ( a(1) ), abs ( a(2) ), abs ( a(3) ), &
    abs ( b(1) ), abs ( b(2) ), abs ( b(3) ), &
    abs ( c(1) ), abs ( c(2) ), abs ( c(3) ) )

  dmax = max ( emax, abs ( d(1) ), abs ( d(2) ), abs ( d(3) ) )

  nrml1 = ab(2) * ac(3) - ab(3) * ac(2)
  nrml2 = ab(3) * ac(1) - ab(1) * ac(3)
  nrml3 = ab(1) * ac(2) - ab(2) * ac(1)

  ddp = ( d(1) - a(1) ) * nrml1 &
      + ( d(2) - a(2) ) * nrml2 &
      + ( d(3) - a(3) ) * nrml3

  if ( abs ( ddp ) <= tol * dmax ) then
    opside = 2
    return
  end if

  emax = max ( emax, abs ( e(1) ), abs ( e(2) ), abs ( e(3) ) )

  edp = ( e(1) - a(1) ) * nrml1 &
      + ( e(2) - a(2) ) * nrml2 &
      + ( e(3) - a(3) ) * nrml3

  if ( abs ( edp ) <= tol * emax ) then
    opside = 0
  else if ( ddp * edp < 0.0D+00 ) then
    opside = 1
  else
    opside = -1
  end if

  return
end
subroutine order3 ( i, j, k )

!*****************************************************************************80
!
!! ORDER3 reorders 3 integers into ascending order.
!
!  Discussion: 
!
!    This routine reorders I, J, K so that I <= J <= K.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2009
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J, K, on output are sorted into
!    nondecreasing order.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) t

  if ( j < i ) then
    if ( k < j ) then
      call i4_swap ( i, k )
    else if ( k < i ) then
      t = i
      i = j
      j = k
      k = t
    else
      call i4_swap ( i, j )
    end if
  else
    if ( k < i ) then
      t = i
      i = k
      k = j
      j = t
    else if ( k < j ) then
      call i4_swap ( j, k )
    end if
  end if

  return
end
subroutine r8mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_DATA_READ reads data from an R8MAT file.
!
!  Discussion:
!
!    The file may contain more than N points, but this routine will
!    return after reading N of them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 )   m
  integer   ( kind = 4 )   n

  integer   ( kind = 4 )   ierror
  character ( len = * )    input_filename
  integer   ( kind = 4 )   input_status
  integer   ( kind = 4 )   input_unit
  integer   ( kind = 4 )   j
  character ( len = 255 )  line
  real ( kind = 8 )   table(m,n)
  real ( kind = 8 )   x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i8)' ) '  Number of values expected per line M = ', m
      write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

  end do

  close ( unit = input_unit )

  return
end
subroutine r8mat_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! R8MAT_HEADER_READ reads the header from an R8MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) M, spatial dimension.
!
!    Output, integer ( kind = 4 ) N, the number of points.
!
  implicit none

  character ( len = * )  input_filename
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  call file_column_count ( input_filename, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 June 2004
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
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer   ( kind = 4 ), parameter :: incx = 5
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2
  integer   ( kind = 4 ) i2hi
  integer   ( kind = 4 ) i2lo
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j2hi
  integer   ( kind = 4 ) j2lo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character ( len = * )  title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,1x,5a14)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) get
  integer ( kind = 4 ) put
  integer ( kind = 4 ) nchar
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  put = 0
  nchar = len_trim ( s )

  do get = 1, nchar

    c = s(get:get)

    if ( c /= ' ' .and. c /= TAB ) then
      put = put + 1
      s(put:put) = c
    end if

  end do

  s(put+1:nchar) = ' '

  return
end
function s_index_last ( s, sub )

!*****************************************************************************80
!
!! S_INDEX_LAST finds the LAST occurrence of a given substring.
!
!  Discussion:
!
!    It returns the location in the string at which the substring SUB is
!    first found, or 0 if the substring does not occur at all.
!
!    The routine is also trailing blank insensitive.  This is very
!    important for those cases where you have stored information in
!    larger variables.  If S is of length 80, and SUB is of
!    length 80, then if S = 'FRED' and SUB = 'RED', a match would
!    not be reported by the standard FORTRAN INDEX, because it treats
!    both variables as being 80 characters long!  This routine assumes that
!    trailing blanks represent garbage!
!
!    This means that this routine cannot be used to find, say, the last
!    occurrence of a substring 'A ', since it assumes the blank space
!    was not specified by the user, but is, rather, padding by the
!    system.  However, as a special case, this routine can properly handle
!    the case where either S or SUB is all blanks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be searched.
!
!    Input, character ( len = * ) SUB, the substring to search for.
!
!    Output, integer ( kind = 4 ) S_INDEX_LAST.  0 if SUB does not occur in
!    the string.  Otherwise S_INDEX_LAST = I, where S(I:I+LENS-1) = SUB,
!    where LENS is the length of SUB, and is the last place
!    this happens.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) llen1
  integer ( kind = 4 ) llen2
  character ( len = * ) s
  integer ( kind = 4 ) s_index_last
  character ( len = * ) sub

  s_index_last = 0

  llen1 = len_trim ( s )
  llen2 = len_trim ( sub )
!
!  In case S or SUB is blanks, use LEN
!
  if ( llen1 == 0 ) then
    llen1 = len ( s )
  end if

  if ( llen2 == 0 ) then
    llen2 = len ( sub )
  end if

  if ( llen1 < llen2 ) then
    return
  end if

  do j = 1, llen1+1-llen2

    i = llen1 + 2 - llen2 - j

    if ( s(i:i+llen2-1) == sub ) then
      s_index_last = i
      return
    end if

  end do

  return
end
subroutine s_to_r8 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Example:
!
!    S                 R
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) R, the real value that was read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical ch_eqi
  character c
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) r
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0D+00
  lchar = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    lchar = lchar + 1
    c = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        lchar = lchar + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = - 1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = - 1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if ( ihave < 11 .and. lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 .or. nchar <= lchar+1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave

    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0D+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0D+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

  return
end
subroutine s_to_r8vec ( s, n, rvec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer ( kind = 4 ) N, the number of values expected.
!
!    Output, real ( kind = 8 ) RVEC(N), the values read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) lchar
  real ( kind = 8 ) rvec(n)
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + lchar

  end do

  return
end
subroutine s_word_count ( s, nword )

!*****************************************************************************80
!
!! S_WORD_COUNT counts the number of "words" in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, integer ( kind = 4 ) NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens
  integer ( kind = 4 ) nword
  character ( len = * ) s

  nword = 0
  lens = len ( s )

  if ( lens <= 0 ) then
    return
  end if

  blank = .true.

  do i = 1, lens

    if ( s(i:i) == ' ' ) then
      blank = .true.
    else if ( blank ) then
      nword = nword + 1
      blank = .false.
    end if

  end do

  return
end
subroutine swapes ( bndcon, i, npt, sizht, fc_num, fc_max, vcl, vm, bf, fc, ht, &
  ntetra, hdavfc, front, back, ifac, ierr )

!*****************************************************************************80
!
!! SWAPES swaps faces in a 3D triangulation.
!
!  Discussion: 
!
!    This routine swaps faces, applying local transformations, in a 3D 
!    triangulation based on the empty circumsphere criterion until (nearly)
!    all faces are locally optimal.  I is the index of the new vertex
!    added to the triangulation, or 0 if an initial triangulation is given.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2005
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
! 
!    Input, logical BNDCON, TRUE iff boundary faces are constrained (i.e. not
!    swapped by local transformations).
!
!    Input, integer ( kind = 4 ) I, the local index of next vertex inserted in
!    triangulation, or 0; if positive, it is assumed I is largest index so far.
!
!    Input, integer ( kind = 4 ) NPT, the number of 3D points to be triangulated.
!
!    Input, integer ( kind = 4 ) SIZHT, the size of hash table HT.
!
!    Input/output, integer ( kind = 4 ) FC_NUM, the number of positions used in FC array.
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum size available for FC array.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:NPT), the indices of vertices of VCL being
!    triangulated.
!
!    Input/output, integer ( kind = 4 ) BF(1:3,1:*), the  array of boundary face records;
!    see DTRIS3.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:FC_MAX), the array of face records;
!    see routine DTRIS3.
!
!    Input/output, integer ( kind = 4 ) HT(0:SIZHT-1), the hash table using direct chaining.
!
!    Input/output, integer ( kind = 4 ) NTETRA, the number of tetrahedra in triangulation.
!
!    Input/output, integer ( kind = 4 ) HDAVFC, the head pointer to available FC records.
!
!    Input/output, integer ( kind = 4 ) FRONT, BACK, the indices of front and back of
!    queue of interior faces for which sphere test is applied.
!
!    Output, integer ( kind = 4 ) IFAC, the index of last face for which sphere test 
!    applied, or 0.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) sizht

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  real ( kind = 8 ) alpha(4)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) bf(3,*)
  integer ( kind = 4 ) bfx(2)
  logical bndcon
  integer ( kind = 4 ) c
  real ( kind = 8 ) center(3)
  integer ( kind = 4 ) d
  integer ( kind = 4 ) dd
  logical degen
  integer ( kind = 4 ) e
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) fc_num
  integer ( kind = 4 ) front
  integer ( kind = 4 ) g
  integer ( kind = 4 ) hdavfc
  integer ( kind = 4 ) ht(0:sizht-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) in
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) ind1
  integer ( kind = 4 ) ind2
  integer ( kind = 4 ) indx(2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kneg
  integer ( kind = 4 ) kzero
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) nbr(2,2)
  integer ( kind = 4 ) ntetra
  real ( kind = 8 ) radsq
  real ( kind = 8 ) vcl(3,*)
  integer ( kind = 4 ) va
  integer ( kind = 4 ) vb
  integer ( kind = 4 ) vc
  integer ( kind = 4 ) vd
  integer ( kind = 4 ) ve
  integer ( kind = 4 ) vm(npt)

  ierr = 0
  ifac = 0

  do

    do

      if ( front == 0 ) then
        return
      end if

      ind = front
      front = fc(7,ind)

      if ( fc(2,ind) /= 0 ) then
        exit
      end if

      if ( ind == fc_num ) then
        fc_num = fc_num - 1
      else
        fc(1,ind) = -hdavfc
        hdavfc = ind
      end if

    end do

    ifac = ind
    fc(7,ind) = -1
    a = fc(1,ind)
    b = fc(2,ind)
    c = fc(3,ind)
    d = fc(4,ind)
    e = fc(5,ind)
    va = vm(a)
    vb = vm(b)
    vc = vm(c)
    vd = vm(d)
    ve = vm(e)

    if ( msglvl == 4 ) then
      write ( *,600) ind,a,b,c,d,e
    end if

    call ccsph ( .true., vcl(1,va), vcl(1,vb), vcl(1,vc), vcl(1,vd), &
      vcl(1,ve), center, radsq, in )

    if ( in == 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SWAPES - Fatal error!'
      write ( *, '(a)' ) '  CCSPH returned IN = 2.'
      ierr = 301
      return
    end if

    if ( 1 <= in ) then

      call baryth ( vcl(1,va), vcl(1,vb), vcl(1,vc), vcl(1,vd), &
        vcl(1,ve), alpha, degen )

      if ( degen ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SWAPES - Fatal error!'
        write ( *, '(a)' ) '  BARYTH detected a degenerate tetrahedron.'
        ierr = 301
        return
      else if ( 0.0D+00 < alpha(4) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SWAPES - Fatal error!'
        write ( *, '(a)' ) '  BARYTH detected 0 < ALPHA(4).'
        ierr = 309
        return
      end if

      kneg = 1
      kzero = 0

      do j = 1, 3
        if ( alpha(j) < 0.0D+00 ) then
          kneg = kneg + 1
        else if ( alpha(j) == 0.0D+00 ) then
          kzero = kzero + 1
        end if
      end do
!
!  Swap 2 tetrahedra for 3.
!
      if ( kneg == 1 .and. kzero == 0 ) then

        call updatf ( a, b, d, c, e, i, npt, sizht, front, back, fc, ht, ierr )
        call updatf ( a, c, d, b, e, i, npt, sizht, front, back, fc, ht, ierr )
        call updatf ( b, c, d, a, e, i, npt, sizht, front, back, fc, ht, ierr )
        call updatf ( a, b, e, c, d, i, npt, sizht, front, back, fc, ht, ierr )
        call updatf ( a, c, e, b, d, i, npt, sizht, front, back, fc, ht, ierr )
        call updatf ( b, c, e, a, d, i, npt, sizht, front, back, fc, ht, ierr )

        if ( ierr /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a, i8)' ) '  UPDATF returned IERR = ', ierr
          return
        end if

        call htdel ( ind, npt, sizht, fc, ht )
        call htins ( ind, a, d, e, b, c, npt, sizht, fc, ht )
        call availf ( hdavfc, fc_num, fc_max, fc, ind, ierr )

        if ( ierr /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a,i8)' ) '  AVAILF returned IERR = ', ierr
          return
        end if

        call htins ( ind, b, d, e, a, c, npt, sizht, fc, ht )
        call availf ( hdavfc, fc_num, fc_max, fc, ind, ierr )

        if ( ierr /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a,i8)' ) '  AVAILF returned IERR = ', ierr
          return
        end if

        call htins ( ind, c, d, e, a, b, npt, sizht, fc, ht )
        ntetra = ntetra + 1

        if ( msglvl == 4 ) then
          write ( *,610)
        end if
!
!  Swap 3 tetrahedra for 2 if possible. Relabel so edge
!  AB would be deleted. Swap if ABDE is in current triangulation.
!
      else if ( kneg == 2 .and. kzero == 0 ) then

        if ( alpha(1) < 0.0D+00 ) then
          call i4_swap ( a, c )
        else if ( alpha(2) < 0.0D+00 ) then
          call i4_swap ( b, c )
        end if

        ind1 = htsrc ( a, b, d, npt, sizht, fc, ht )

        if ( ind1 <= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a)' ) '  HTSRC returned IND1 <= 0.'
          ierr = 300
          return
        end if

        if ( fc(4,ind1) == e .or. fc(5,ind1) == e ) then

          call updatf ( a, c, d, b, e, i, npt, sizht, front, back, fc, ht, ierr )
          call updatf ( a, c, e, b, d, i, npt, sizht, front, back, fc, ht, ierr )
          call updatf ( a, d, e, b, c, i, npt, sizht, front, back, fc, ht, ierr )
          call updatf ( b, c, d, a, e, i, npt, sizht, front, back, fc, ht, ierr )
          call updatf ( b, c, e, a, d, i, npt, sizht, front, back, fc, ht, ierr )
          call updatf ( b, d, e, a, c, i, npt, sizht, front, back, fc, ht, ierr )

          if ( ierr /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SWAPES - Fatal error!'
            write ( *, '(a,i8)' ) '  UPDATF returned IERR = ', ierr
            return
          end if

          call htdel ( ind, npt, sizht, fc, ht )
          call htins ( ind, c, d, e, a, b, npt, sizht, fc, ht )
          call htdel ( ind1, npt, sizht, fc, ht )

          if ( 0 <= fc(7,ind1) ) then
            fc(2,ind1) = 0
          else
            if ( ind1 == fc_num ) then
              fc_num = fc_num - 1
            else
              fc(1,ind1) = -hdavfc
              hdavfc = ind1
            end if
          end if

          ind1 = htsrc ( a, b, e, npt, sizht, fc, ht )

          if ( ind1 <= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SWAPES - Fatal error!'
            write ( *, '(a)' ) '  HTSRC returned IND1 <= 0.'
            ierr = 300
            return
          end if

          call htdel ( ind1, npt, sizht, fc, ht )

          if ( 0 <= fc(7,ind1) ) then

            fc(2,ind1) = 0

          else

            if ( ind1 == fc_num ) then
              fc_num = fc_num - 1
            else
              fc(1,ind1) = -hdavfc
              hdavfc = ind1
            end if

          end if

          ntetra = ntetra - 1

          if ( msglvl == 4 ) then
            write ( *,620) c,d,e
          end if

        else

          if ( msglvl == 4 ) then
            write ( *,630) a,b,d,e
          end if

        end if
!
!  Coplanar faces: swap 2 tetrahedra for 2 if boundary faces
!  (and BNDCON is .FALSE.), else do pair of 2 for 2 swaps if
!  possible.  Relabel vertices so that DE intersects AB.
!  Also swap if necessary to make A < B and D < E.
!
      else if ( kneg == 1 .and. kzero == 1 ) then

        if ( alpha(1) == 0.0D+00 ) then

          call i4_swap ( a, c )

        else if ( alpha(2) == 0.0D+00 ) then

          call i4_swap ( b, c )

        end if

        if ( b < a ) then
          call i4_swap ( a, b )
        end if

        if ( e < d ) then
          call i4_swap ( d, e )
        end if

        ind1 = htsrc ( a, b, d, npt, sizht, fc, ht )

        if ( ind1 <= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a)' ) '  HTSRC returned IND1 <= 0.'
          ierr = 300
          return
        end if

        ind2 = htsrc ( a, b, e, npt, sizht, fc, ht )

        if ( ind2 <= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SWAPES - Fatal error!'
          write ( *, '(a)' ) '  HTSRC returned IND2 <= 0.'
          ierr = 300
          return
        end if

        if ( fc(4,ind1) == c ) then
          f = fc(5,ind1)
        else
          f = fc(4,ind1)
        end if

        if ( fc(4,ind2) == c ) then
          g = fc(5,ind2)
        else
          g = fc(4,ind2)
        end if

        if ( f <= 0 .and. g <= 0 ) then

          if ( .not. bndcon ) then

            call updatf(a,c,d,b,e,i,npt,sizht,front,back,fc,ht, ierr )
            call updatf(a,c,e,b,d,i,npt,sizht,front,back,fc,ht, ierr )
            call updatf(b,c,d,a,e,i,npt,sizht,front,back,fc,ht, ierr )
            call updatf(b,c,e,a,d,i,npt,sizht,front,back,fc,ht, ierr )

            if ( ierr /= 0 ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SWAPES - Fatal error!'
              write ( *, '(a,i8)' ) '  UPDATF returned IERR = ', ierr
              return
            end if

            call htdel(ind,npt,sizht,fc,ht)
            call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)
            call htdel(ind1,npt,sizht,fc,ht)
            call htins(ind1,a,d,e,c,fc(5,ind1),npt,sizht,fc,ht)
            call htdel(ind2,npt,sizht,fc,ht)
            call htins(ind2,b,d,e,c,fc(5,ind2),npt,sizht,fc,ht)
            indx(1) = ind1
            indx(2) = ind2
            bfx(1) = -fc(5,ind1)
            bfx(2) = -fc(5,ind2)
            dd = d

            do j = 1, 2

              if ( j == 2 ) then
                dd = e
              end if

              if ( dd < a ) then
                nbr(j,1) = bf(3,bfx(j))
                nbr(j,2) = bf(2,bfx(j))
              else if ( dd < b ) then
                nbr(j,1) = bf(3,bfx(j))
                nbr(j,2) = bf(1,bfx(j))
              else
                nbr(j,1) = bf(2,bfx(j))
                nbr(j,2) = bf(1,bfx(j))
              end if

            end do

            aa = a
            k = -fc(5,nbr(1,2))
  
            do j = 1, 2

              if ( j == 2 ) then
                aa = b
                k = -fc(5,nbr(2,1))
              end if

              if ( aa < d ) then
                bf(1,bfx(j)) = indx(3-j)
                bf(2,bfx(j)) = nbr(2,j)
                bf(3,bfx(j)) = nbr(1,j)
              else if ( aa < e ) then
                bf(1,bfx(j)) = nbr(2,j)
                bf(2,bfx(j)) = indx(3-j)
                bf(3,bfx(j)) = nbr(1,j)
              else
                bf(1,bfx(j)) = nbr(2,j)
                bf(2,bfx(j)) = nbr(1,j)
                bf(3,bfx(j)) = indx(3-j)
              end if

              if ( bf(1,k) == indx(j) ) then
                bf(1,k) = indx(3-j)
              else if ( bf(2,k) == indx(j) ) then
                bf(2,k) = indx(3-j)
              else
                bf(3,k) = indx(3-j)
              end if

            end do

            if ( msglvl == 4 ) then
              write ( *,640) a,b,d,e
            end if
  
          end if

        else if ( f == g ) then

          call updatf(a,c,d,b,e,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(a,c,e,b,d,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(b,c,d,a,e,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(b,c,e,a,d,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(a,d,f,b,e,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(a,e,f,b,d,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(b,d,f,a,e,i,npt,sizht,front,back,fc,ht, ierr )
          call updatf(b,e,f,a,d,i,npt,sizht,front,back,fc,ht, ierr )

          if ( ierr /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SWAPES - Fatal error!'
            write ( *, '(a,i8)' ) '  UPDATF returned IERR = ', ierr
            return
          end if

          call htdel(ind,npt,sizht,fc,ht)
          call htins(ind,c,d,e,a,b,npt,sizht,fc,ht)

          ind = htsrc ( a, b, f, npt, sizht, fc, ht )

          if ( ind <= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SWAPES - Fatal error!'
            write ( *, '(a)' ) '  HTSRC returned IND <= 0.'
            ierr = 300
            return
          end if

          call htdel(ind,npt,sizht,fc,ht)

          if ( 0 <= fc(7,ind) ) then
            fc(2,ind) = 0
            call availf ( hdavfc, fc_num, fc_max, fc, ind, ierr )
            if ( ierr /= 0 ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SWAPES - Fatal error!'
              write ( *, '(a,i8)' ) '  AVAILF returned IERR = ', ierr
              return
            end if
          end if

          call htins(ind,d,e,f,a,b,npt,sizht,fc,ht)
          call htdel(ind1,npt,sizht,fc,ht)
          j = fc(7,ind1)
          call htins(ind1,a,d,e,c,f,npt,sizht,fc,ht)
          fc(7,ind1) = j
          call htdel(ind2,npt,sizht,fc,ht)
          j = fc(7,ind2)
          call htins(ind2,b,d,e,c,f,npt,sizht,fc,ht)
          fc(7,ind2) = j

          if ( i <= 0 .and. fc(7,ind1) == -1 ) then
            fc(7,ind1) = 0
            if ( front == 0 ) then
              front = ind1
            else
              fc(7,back) = ind1
            end if
            back = ind1
          end if

          if ( i <= 0 .and. fc(7,ind2) == -1 ) then
            fc(7,ind2) = 0
            if ( front == 0 ) then
              front = ind2
            else
              fc(7,back) = ind2
            end if
            back = ind2
          end if

          if ( msglvl == 4 ) then
            write ( *,650) a,b,d,e,f
          end if

        else

          if ( msglvl == 4 ) then
            write ( *,660) a,b,d,e,f,g
          end if

        end if

      end if

    end if

  end do

  600 format (1x,'index =',i7,' : ',5i7)
  610 format (4x,'swap 2-3')
  620 format (4x,'swap 3-2 with new common face:',3i7)
  630 format (4x,'swap 3-2 not poss, tetra missing:',4i7)
  640 format (4x,'swap 2-2: edge ',2i7,' repl by ',2i7)
  650 format (4x,'swap 4-4: edge ',2i7,' repl by ',2i7,'   f =',i7)
  660 format (4x,'swap 4-4 not poss: a,b,d,e,f,g =',6i7)

  return
end
subroutine tetlst ( fc_max, fc_num, vm, fc, tetra_num, tetra_num2, tetra )

!*****************************************************************************80
!
!! TETLST constructs a list of tetrahedra from the FC array.
!
!  Discussion: 
!
!    This routine constructs a list of tetrahedra from the FC array. 
!
!    Global vertex indices from VM are produced.  The vertex indices for each
!    tetrahedron are sorted in increasing order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 September 2005
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FC_MAX, the maximum number of positions in the FC array.
!
!    Input, integer ( kind = 4 ) FC_NUM, the number of positions used in the FC array.
!
!    Input, integer ( kind = 4 ) VM(1:*), the indices of vertices of VCL that are
!    triangulated.
!
!    Input, integer ( kind = 4 ) FC(7,FC_MAX), array of face records; see routine DTRIS3.
!
!    Input, integer ( kind = 4 ) TETRA_NUM, the number of tetrahedrons expected.
!
!    Output, integer ( kind = 4 ) TETRA_NUM2, the number of tetrahedrons found.
!
!    Output, integer ( kind = 4 ) TETRA(4,TETRA_NUM), contains global tetrahedron indices;
!    it is assumed there is enough space.
!
  implicit none

  integer ( kind = 4 ) fc_max
  integer ( kind = 4 ) tetra_num

  integer ( kind = 4 ) a
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) fc_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) t(4)
  integer ( kind = 4 ) tetra(4,tetra_num)
  integer ( kind = 4 ) tetra_num2
  integer ( kind = 4 ) vm(*)

  tetra_num2 = 0

  do i = 1, fc_num

    if ( fc(1,i) <= 0 )  then
      cycle
    end if

    do k = 4, 5

      if ( fc(3,i) < fc(k,i) ) then

        tetra_num2 = tetra_num2 + 1

        if ( tetra_num2 <= tetra_num ) then
          tetra(1,tetra_num2) = fc(1,i)
          tetra(2,tetra_num2) = fc(2,i)
          tetra(3,tetra_num2) = fc(3,i)
          tetra(4,tetra_num2) = fc(k,i)
        end if

      end if

    end do

  end do

  if ( tetra_num2 /= tetra_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETLST - Warning!'
    write ( *, '(a)' ) '  Inconsistent tetrahedron information.'
    write ( *, '(a,i8)' ) '  Was expecting TETRA_NUM = ', tetra_num
    write ( *, '(a,i8)' ) '  Found TETRA_NUM2 = ', tetra_num2
    return
  end if

  do k = 1, tetra_num2

    t(1:4) = vm(tetra(1:4,k))

    do i = 1, 3
      l = i
      do j = i+1, 4
        if ( t(j) < t(l) ) then
          l = j
        end if
      end do
      a = t(i)
      t(i) = t(l)
      t(l) = a
    end do

    tetra(1:4,k) = t(1:4)

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
subroutine updatf ( a, b, c, d, e, i, n, p, front, back, fc, ht, ierr )

!*****************************************************************************80
!
!! UPDATF updates a record in FC after a local transformation.
!
!  Discussion: 
!
!    This routine updates a record in FC due to a local transformation.
!
!    Tetrahedron ABCD becomes ABCE. Add face ABC to queue if it is
!    interior face, not yet in queue, and its largest index isn't I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2009
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, C, the first 3 fields of FC record 
!    (in any order).
!
!    Input, integer ( kind = 4 ) D, E, the fourth vertex indices of old and 
!    new tetrahedrons.
!
!    Input, integer ( kind = 4 ) I, the vertex index determining whether 
!    face put on queue.
!
!    Input, integer ( kind = 4 ) N, the upper bound on vertex indices.
!
!    Input, integer ( kind = 4 ) P, the size of hash table.
!
!    Input/output, integer ( kind = 4 ) FRONT, BACK, the front and back 
!    pointers of queue.
!
!    Input, integer ( kind = 4 ) FC(1:7,1:*), the array of face records; 
!    see routine DTRIS3.
!
!    Input, integer ( kind = 4 ) HT(0:P-1), the hash table using direct 
!    chaining.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless 
!    an error occurred.
!
  implicit none

  integer ( kind = 4 ) p

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) n

  ierr = 0
  ind = htsrc ( a, b, c, n, p, fc, ht )

  if ( ind <= 0 ) then
    ierr = 300
    return
  end if

  if ( fc(4,ind) == d ) then
    fc(4,ind) = e
  else
    fc(5,ind) = e
  end if

  if ( fc(7,ind) == -1 .and. fc(3,ind) /= i .and. 0 < fc(5,ind) ) then

    fc(7,ind) = 0

    if ( front == 0 ) then
      front = ind
    else
      fc(7,back) = ind
    end if

    back = ind
  end if

  return
end
subroutine vbfac ( pt, ctr, vcl, vm, bf, fc, topv, topnv )

!*****************************************************************************80
!
!! VBFAC determines the boundary faces of a 3D triangulation.
!
!  Discussion: 
!
!    This routine determines boundary faces of a 3D triangulation visible
!    from point PT, given a starting visible boundary face.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2009
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) PT(1:3), the 3D point.
!
!    Input, real ( kind = 8 ) CTR(1:3), the 3D point in interior of
!    triangulation.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:*), the indices of vertices of VCL being triangulated.
!
!    Input, integer ( kind = 4 ) BF(1:3,1:*), the array of boundary face records; see DTRIS3.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:*), the array of face records; see routine
!    DTRIS3; row 7 is used for links of 3 stacks in this routine.  On output, 
!    FC(7,*) has been updated, so that only stack of visible boundary 
!    faces remains.
!
!    Input/output, integer ( kind = 4 ) TOPV.  On input, index of FC of visible boundary
!    face.  On output, index of top of stack of visible boundary faces.
!
!    Input, integer ( kind = 4 ) TOPNV, the index of top of stack of boundary faces 
!    already found to be not visible from PT, or 0 for empty stack.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) bf(3,*)
  real ( kind = 8 ) ctr(3)
  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) op
  integer ( kind = 4 ) opside
  real ( kind = 8 ) pt(3)
  integer ( kind = 4 ) ptr
  integer ( kind = 4 ) topn
  integer ( kind = 4 ) topnv
  integer ( kind = 4 ) topt
  integer ( kind = 4 ) topv
  integer ( kind = 4 ) va
  integer ( kind = 4 ) vb
  integer ( kind = 4 ) vc
  integer ( kind = 4 ) vm(*)
  real ( kind = 8 ) vcl(3,*)
!
!  TOPN is index of top of stack of non-visible boundary faces.
!  TOPT is index of top of stack of boundary faces to be tested.
!
  topn = topnv
  topt = 0
  fc(7,topv) = 0
  k = -fc(5,topv)

  do j = 1, 3
    nbr = bf(j,k)
    if ( fc(7,nbr) == -1 ) then
      fc(7,nbr) = topt
      topt = nbr
    end if
  end do

  do while ( topt /= 0 )

    ptr = topt
    topt = fc(7,ptr)
    va = vm(fc(1,ptr))
    vb = vm(fc(2,ptr))
    vc = vm(fc(3,ptr))
    op = opside ( vcl(1,va), vcl(1,vb), vcl(1,vc), ctr, pt )

    if ( op == 2 ) then
      ierr = 301
      return
    end if

    if ( op == 1 ) then

      fc(7,ptr) = topv
      topv = ptr
      k = -fc(5,ptr)

      do j = 1, 3
        nbr = bf(j,k)
        if ( fc(7,nbr) == -1 ) then
          fc(7,nbr) = topt
          topt = nbr
        end if
      end do

    else

      fc(7,ptr) = topn
      topn = ptr

    end if

  end do
!
!  For boundary faces not visible from PT, set FC(7,*) = -1.
!
  do while ( topn /= 0 ) 
    ptr = topn
    topn = fc(7,ptr)
    fc(7,ptr) = -1
  end do

  return
end
