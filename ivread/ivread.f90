program main

!*****************************************************************************80
!
!! MAIN is the main program for IVREAD.
!
!  Discussion:
!
!    IVREAD is the main program for converting a graphics file.
!
!    This program was written for the Studio for Creative Inquiry.
!
!    This program was originally intended to read an Inventor 3D graphics
!    file, and write a corresponding VLA graphics file.  
!
!    However, it can also try to read files of several different types,
!    convert the data internally to line data, and write out files
!    of various types.
!
!    DXF and Inventor file formats are sophisticated and extensive.
!    This program was only designed to work on a few simple cases, and
!    may easily be confused by unexpected (though legal) data.
!
!    A sketch of the file formats and data structures is included
!    in most of the READ and WRITE routines.
!
!    The sizes of COR3_MAX, FACE_MAX and LINE_MAX control how much 
!    information this program can handle. 
!
!    The "head.25.iv" file has 240146 faces and 119,736 points.
!    The "fish.iv" file has about 76000 faces.
!    The "brain.iv" file has between 40000 and 75000 faces.
!
!  Development:
!
!    27 August 2003: Added OFF_READ and OFF_WRITE.
!
!    18 October 2001: TS_READ and TS_WRITE handle ambientlight by assigning
!    it to color of default material #1.  IV_READ and IV_WRITE go along
!    with this.  See if other formats have a single ambient light setting.
!
!    Also, trying to get line colors to work.  I see that IV_READ doesn't
!    even read line colors in...trying to fix that.  Then I can read in
!    MATERIALS.IV and output MATERIALS.TS.
!
!    16 October 2001: Added draft TS_READ.
!
!    14 October 2001: Added POINT data.
!
!    13 October 2001: Added TS_WRITE.
!
!    05 June 2001: Restored XYZ_WRITE.
!
!    26 May 1999: Added LINE_PRUNE switch, which will try to cut down
!    (by about half) the number of superfluous lines created when
!    faces are turned into lines by FACE_TO_LINE for VLA_WRITE output.
!
!    22 May 1999: For VLA output files, the program now will automatically 
!    try to temporarily convert all face information to line information.  No
!    sophisticated attempt is made to delete superfluous lines (the
!    way FACE_TO_EDGE tries.)
!
!    The "<<" merge command:
!
!      On 20 April 1999, the "<<" command was added, to allow data
!      from two or more files to be merged.  It works, on a simple example,
!      for OBJ files.  However, some tuning of OBJ_READ was necessary.
!      Similar testing and tuning must be done to the other READ routines
!      before they will work with this option.
!      On 21 April 1999, the "<<" command worked on a simple example
!      using two IV files as input.  SMF_READ was also updated for the "<<"
!      command, but not tested.  ASE_READ, DXF_READ, HRC_READ,
!      STLA_READ and VLA_READ may already be OK.
!      On 22 April 1999, "<<" command works with ASE_READ, HRC_READ,
!      SMF_READ and STLA_READ.
!
!    The "MatrixTransform" field in IV_READ/IV_WRITE.
!      I'm having problems because I am reading in an IV file that has
!      a matrix transform that is not the identity.  I just ignore it,
!      and so my data is not rotated and scaled, when it should be.
!      As a start toward addressing this issue, I have IV_WRITE 
!      writing out the current transform matrix.  One problem with
!      Inventor is that the transform matrix can be specified on
!      every level, and the actual transform matrix that applies
!      has to be deduced from where you are in the tree.
!      Right now, all I've done is have IV_READ read the matrix,
!      and IV_WRITE write it out.  No concatenation is possible right now,
!      but my kludgy code will at least apply ONE transformation matrix
!      to the data, in IV_READ, anyway...
!
!    Adding material/normal/texture binding stubs, because 
!    A) new SMF format allows it;
!    B) Inventor uses it;
!    C) SCI wants to do textures eventually.
!
!    SMF_READ and SMF_WRITE can now read and write face and vertex colors 
!    of SMF2.0 files.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Parameter, integer COR3_MAX, the maximum number of points.
!
!    Parameter, integer EDGE_MAX, the maximum number of edges.
!
!    Parameter, integer FACE_MAX, the maximum number of faces.
!
!    Parameter, integer LINE_MAX, the maximum number of line definition items.
!
!    Parameter, integer MATERIAL_MAX, the maximum number of materials.
!
!    Parameter, integer ORDER_MAX, the maximum number of vertices per face.
!
!    Parameter, integer POINT_MAX, the maximum number of points to display.
!
!    Parameter, integer TEXTURE_MAX, the maximum number of textures.
!
  implicit none

  integer ( kind = 4 ), parameter :: cor3_max = 50000
  integer ( kind = 4 ), parameter :: edge_max = 100
  integer ( kind = 4 ), parameter :: face_max = 50000
  integer ( kind = 4 ), parameter :: line_max = 50000
  integer ( kind = 4 ), parameter :: material_max = 200
  integer ( kind = 4 ), parameter :: order_max = 6
  integer ( kind = 4 ), parameter :: point_max = 1000
  integer ( kind = 4 ), parameter :: texture_max = 10

  integer ( kind = 4 ) arg_num
  logical byte_swap
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_material(cor3_max)
  real ( kind = 4 ) cor3_new(3,cor3_max)
  real ( kind = 4 ) cor3_normal(3,cor3_max)
  real ( kind = 4 ) cor3_tex_uv(2,cor3_max)
  logical debug
  integer ( kind = 4 ) edge(4,edge_max)
  integer ( kind = 4 ) face(order_max,face_max)
  real ( kind = 4 ) face_area(face_max)
  integer ( kind = 4 ) face_material(face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_object(face_max)
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) face_rank(face_max)
  real ( kind = 4 ) face_tex_uv(2,face_max)
  integer ( kind = 4 ) face_tier(face_max)
  character ( len = 255 ) filein_name
  character ( len = 255 ) fileout_name
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ipxfargc
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_prune
  integer ( kind = 4 ) list(cor3_max)
  character ( len = 255 ) material_name(material_max)
  real ( kind = 4 ) material_rgba(4,material_max)
  character ( len = 255 ) object_name
  integer ( kind = 4 ) point(point_max)
  integer ( kind = 4 ) point_num
  character ( len = 255 ) texture_name(texture_max)
  real ( kind = 4 ) texture_temp(2,order_max*face_max)
  real ( kind = 4 ) transform_matrix(4,4)
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  real ( kind = 4 ) vertex_tex_uv(2,order_max,face_max)

  call timestamp ( )
!
!  Initialize a few program variables.
!
  byte_swap = .true.
  debug = .false.
  filein_name = ' '
  fileout_name = ' '
  ierror = 0
  line_prune = 1
!
!  Get the number of command line arguments.
!
!  Old style:
!
  arg_num = iargc ( )
!
!  New style:
!
! arg_num = ipxfargc ( )

  if ( 2 <= arg_num ) then

    call command_line ( arg_num, cor3, cor3_material, cor3_max, cor3_normal, &
      cor3_tex_uv, debug, face, face_area, face_material, face_max, &
      face_normal, face_order, face_tex_uv, filein_name, fileout_name, &
      ierror, line_dex, line_material, line_prune, material_name, &
      material_rgba, line_max, material_max, order_max, texture_max, &
      object_name, point, point_max, point_num, &
      texture_name, texture_temp, &
      transform_matrix, vertex_material, vertex_normal, vertex_tex_uv )

  else

    call interact ( byte_swap, cor3, cor3_material, cor3_max, cor3_new, &   
      cor3_normal, cor3_tex_uv, debug, edge, edge_max, face, face_area, &
      face_material, face_max, face_normal, face_object, face_order, &
      face_rank, face_tex_uv, face_tier, filein_name, fileout_name, ierror, &
      line_dex, line_material, line_max, line_prune, list, material_max, &
      material_name, material_rgba, order_max, object_name, &
      point, point_max, point_num, texture_max, texture_name, texture_temp, &
      transform_matrix, vertex_material, vertex_normal, vertex_tex_uv )

  end if

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IVREAD - Error exit.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IVREAD - Normal exit.'
  end if

  stop
end
function angle_rad_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )

!*****************************************************************************80
!
!! ANGLE_RAD_3D returns the angle in radians between two vectors in 3D.
!
!  Discussion:
!
!    The routine always computes the SMALLER of the two angles between
!    two vectors.  Thus, if the vectors make an (exterior) angle of 
!    1.5 radians, the (interior) angle of 0.5 radians will be reported.
!
!    X dot Y = Norm(X) * Norm(Y) * Cos ( Angle(X,Y) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are three points
!    which define the vectors.  The vectors are:
!    ( X1-X2, Y1-Y2, Z1-Z2 ) and ( X3-X2, Y3-Y2, Z3-Z2 ).
!
!    Output, real ( kind = 4 ) ANGLE_RAD_3D, the angle between the two vectors, in radians.
!    This value will always be between 0 and PI.  If either vector has 
!    zero length, then the angle is returned as zero.
!
  implicit none

  real ( kind = 4 ) angle_rad_3d
  real ( kind = 4 ) dot
  real ( kind = 4 ) dot0_3d
  real ( kind = 4 ) enorm0_3d
  real ( kind = 4 ) v1norm
  real ( kind = 4 ) v2norm
  real ( kind = 4 ) x1
  real ( kind = 4 ) x2
  real ( kind = 4 ) x3
  real ( kind = 4 ) y1
  real ( kind = 4 ) y2
  real ( kind = 4 ) y3
  real ( kind = 4 ) z1
  real ( kind = 4 ) z2
  real ( kind = 4 ) z3

  dot = dot0_3d ( x2, y2, y2, x1, y1, z1, x3, y3, z3 )
  v1norm = enorm0_3d ( x1, y1, z1, x2, y2, z2 )
  v2norm = enorm0_3d ( x3, y3, z3, x2, y2, z2 )
 
  if ( v1norm == 0.0E+00 .or. v2norm == 0.0E+00 ) then
    angle_rad_3d = 0.0E+00
  else
    angle_rad_3d = acos ( dot / ( v1norm * v2norm ) )
  end if
 
  return
end
subroutine ase_read ( bad_num, cor3, cor3_material, cor3_max, cor3_num, &
  face, face_material, face_max, face_normal, face_num, face_order, &
  filein_name, ierror, iunit, material_max, material_name, &
  material_num, material_rgba, order_max, text_num, &
  vertex_material, vertex_normal )

!*****************************************************************************80
!
!! ASE_READ reads graphics information from an ASE file.
!
!  Problems:
!
!    Processing of the MESH_FACELIST assumes faces are always triangles
!    or quadrilaterals.
!
!  Discussion:
!
!    It is intended that the information read from the file can
!    either start a whole new graphics object, or simply be added
!    to a current graphics object via the '<<' command.  
!
!    This is controlled by whether the input values have been zeroed 
!    out or not.  This routine simply tacks on the information it 
!    finds to the current graphics object.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input/output, integer ( kind = 4 ) COR3_MATERIAL(COR3_MAX), the material index 
!    of each node.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input/output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input/output, real ( kind = 4 ) FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices 
!    per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which data is read.
!
!    Input/output, character ( len = * ) MATERIAL_NAME(MATERIAL_MAX), 
!    material names.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input/output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input/output, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Output, real ( kind = 4 ) TRANSFORM_MATRIX(4,4), the transformation matrix.
!
!    Input/output, integer ( kind = 4 ) VERTEX_MATERIAL(ORDER_MAX,FACE_MAX), 
!    vertex materials.
!
!    Input/output, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals 
!    at vertices.
!
  implicit none

  logical, parameter :: debug = .false.
  integer ( kind = 4 ), parameter :: level_max = 10
  integer ( kind = 4 ), parameter :: OFFSET = 1

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max

  integer ( kind = 4 ) bad_num
  real ( kind = 4 ) bval
  character ( len = 4 ) char4
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_material(cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) cor3_num_old
  logical done
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_num_old
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  real ( kind = 4 ) gval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) imat
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) iword
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) level
  character ( len = 256 ) level_name(0:level_max)
  character ( len = 256 ) line
  character ( len = * ) material_name(material_max)
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ) nlbrack
  integer ( kind = 4 ) node
  integer ( kind = 4 ) nrbrack
  real ( kind = 4 ) rgba(4)
  real ( kind = 4 ) rval
  logical s_eqi
  real ( kind = 4 ) temp
  integer ( kind = 4 ) text_num
  real ( kind = 4 ) transform_matrix(4,4)
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  character ( len = 256 ) word
  character ( len = 256 ) word1
  character ( len = 256 ) wordm1
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  real ( kind = 4 ) z
!
  ierror = 0
  level = 0
  level_name(0) = 'Top'
  cor3_num_old = cor3_num
  face_num_old = face_num
  nlbrack = 0
  nrbrack = 0
  call tmat_init ( transform_matrix )
 
  word = ' '
  wordm1 = ' '
!
!  Read a line of text from the file.
!
  do
 
    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    text_num = text_num + 1
!
!  Replace any control characters (in particular, TAB's) by blanks.
!
    call s_control_blank ( line )
 
    done = .true.
    iword = 0
!
!  Read a word from the line.
!
20  continue
 
    if ( word /= ' ' ) then
      wordm1 = word
    end if
 
    call word_next_read ( line, word, done )
!
!  If no more words in this line, read a new line.
!
    if ( done ) then
      cycle
    end if
 
    iword = iword + 1
    if ( iword == 1 ) then
      word1 = word
    end if
!
!  In cases where the word is a left bracket, record the level name,
!  and for right brackets, do a parity check.
!
    if ( word == '{' ) then
 
      nlbrack = nlbrack + 1

      level = nlbrack - nrbrack
      if ( level < 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ASE_READ - Error!'
        write ( *, '(a)' ) '  Too many right brackets!'
        level = 0
      else if ( level_max < level ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ASE_READ - Error!'
        write ( *, '(a)' ) '  Too many left brackets!'
        level = level_max
      end if

      level_name(level) = wordm1
 
      if ( debug ) then
        do i = 0, level
          write ( *, '(i6,2x,a)' ) i, trim ( level_name(i) )
        end do
      end if

    else if ( word == '}' ) then
 
      nrbrack = nrbrack + 1
 
      if ( nlbrack < nrbrack ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ASE_READ - Fatal error!'
        write ( *, '(a,i6)' ) '  Extraneous right bracket on line ', text_num
        write ( *, '(a)' ) trim ( line )
        write ( *, '(a)' ) '  Currently processing field:'
        write ( *, '(a)' ) trim ( level_name(level) )
        ierror = 1
        return
      end if
 
    end if
!
!  *3DSMAX_ASCIIEXPORT  200
!
    if ( word1 == '*3DSMAX_ASCIIEXPORT' ) then
 
      cycle
!
!  *COMMENT
!
    else if ( word1 == '*COMMENT' ) then
 
      cycle
!
!  *GEOMOBJECT
!
    else if ( level_name(level) == '*GEOMOBJECT' ) then
 
      if ( word == '{' ) then
        go to 20
      else if ( word == '}' ) then

        level = nlbrack - nrbrack

      else if ( word == '*NODE_NAME' ) then
        cycle
      else if ( word == '*NODE_TM' ) then
        go to 20
      else if ( word == '*MESH' ) then
        go to 20
      else if ( word == '*PROP_CASTSHADOW' ) then
        cycle
      else if ( word == '*PROP_MOTIONBLUR' ) then
        cycle
      else if ( word == '*PROP_RECVSHADOW' ) then
        cycle
      else
        go to 99
      end if
!
!  *MESH
!
    else if ( level_name(level) == '*MESH' ) then
 
      if ( word == '{' ) then
        go to 20
      else if ( word == '}' ) then
        level = nlbrack - nrbrack
      else if ( word == '*MESH_CFACELIST' ) then
        go to 20
      else if ( word == '*MESH_CVERTLIST' ) then
        go to 20
      else if ( word == '*MESH_FACE_LIST' ) then
        go to 20
      else if ( word == '*MESH_NORMALS' ) then
        go to 20
      else if ( word == '*MESH_NUMCVERTEX' ) then
        cycle
      else if ( word == '*MESH_NUMCVFACES' ) then
        cycle
      else if ( word == '*MESH_NUMFACES' ) then
        cycle
      else if ( word == '*MESH_NUMTVERTEX' ) then
        cycle
      else if ( word == '*MESH_NUMTVFACES' ) then
        cycle
      else if ( word == '*MESH_NUMVERTEX' ) then
        cycle
      else if ( word == '*MESH_TFACELIST' ) then
        go to 20
      else if ( word == '*MESH_TVERTLIST' ) then
        go to 20
      else if ( word == '*MESH_VERTEX_LIST' ) then
        go to 20
      else if ( word == '*TIMEVALUE' ) then
        cycle
      else
        bad_num = bad_num + 1
        if ( bad_num < 10 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'ASE_READ - Error!'
          write ( *, '(a)' ) '  Bad data while reading *MESH.'
          write ( *, '(a)' ) trim ( line )
        end if
        cycle
      end if
!
!  *MESH_CFACELIST
!
    else if ( level_name(level) == '*MESH_CFACELIST' ) then
 
      if ( word == '{' ) then
        go to 20
      else if ( word == '}' ) then
        level = nlbrack - nrbrack
      else if ( word == '*MESH_CFACE' ) then
        cycle
      else
        go to 99
      end if
!
!  *MESH_CVERTLIST
!
!  Mesh vertex indices must be incremented by COR3_NUM_OLD before 
!  being stored in the internal array.
!
    else if ( level_name(level) == '*MESH_CVERTLIST' ) then
 
      if ( word == '{' ) then
 
        go to 20
 
      else if ( word == '}' ) then
 
        level = nlbrack - nrbrack

      else if ( word == '*MESH_VERTCOL' ) then
 
        call word_next_read ( line, word, done )
        call s_to_i4 ( word, i, ierror, lchar )
        i = i + cor3_num_old + OFFSET
 
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, rval, ierror, lchar )
 
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, gval, ierror, lchar )
 
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, bval, ierror, lchar )

        rgba(1) = rval
        rgba(2) = gval
        rgba(3) = bval
        rgba(4) = 1.0E+00

        if ( material_num <= 1000 ) then
          call r4col_find ( 4, 4, material_num, material_rgba, rgba, imat )
        else
          imat = 0
        end if

        if ( imat == 0 ) then

          material_num = material_num + 1

          if ( material_num <= material_max ) then

            call i4_to_s_zero ( material_num, char4 )
            material_name(material_num) = 'Material_' // char4
            material_rgba(1:4,material_num) = rgba(1:4)
            imat = material_num

          else

            imat = 0

          end if

        end if

        cor3_material(i) = imat
   
      else
        go to 99
      end if
!
!  *MESH_FACE_LIST
!
!  WARNING:
!  The following coding assumes that the faces are always triangles
!  or quadrilaterals, but not higher order.
!
    else if ( level_name(level) == '*MESH_FACE_LIST' ) then
 
      if ( word == '{' ) then
 
        go to 20
 
      else if ( word == '}' ) then
 
        level = nlbrack - nrbrack

      else if ( word == '*MESH_FACE' ) then
 
        face_num = face_num + 1

        if ( face_num <= face_max ) then

          call word_next_read ( line, word, done )
          call s_to_i4 ( word, i, ierror, lchar )
          face_material(face_num) = material_num
          face_order(face_num) = 0

          call word_next_read ( line, word, done )
          call word_next_read ( line, word, done )
          call s_to_i4 ( word, i, ierror, lchar )
          face(1,face_num) = i + cor3_num_old + OFFSET
          vertex_material(1,face_num) = material_num
          face_order(face_num) = face_order(face_num) + 1

          call word_next_read ( line, word, done )
          call word_next_read ( line, word, done )
          call s_to_i4 ( word, i, ierror, lchar )
          face(2,face_num) = i + cor3_num_old + OFFSET
          vertex_material(2,face_num) = material_num
          face_order(face_num) = face_order(face_num) + 1

          call word_next_read ( line, word, done )
          call word_next_read ( line, word, done )
          call s_to_i4 ( word, i, ierror, lchar )
          face(3,face_num) = i + cor3_num_old + OFFSET
          vertex_material(3,face_num) = material_num
          face_order(face_num) = face_order(face_num) + 1

          call word_next_read ( line, word, done )
          if ( s_eqi ( word, 'D:' ) ) then
            call word_next_read ( line, word, done )
            call s_to_i4 ( word, i, ierror, lchar )
            face(4,face_num) = i + cor3_num_old + OFFSET
            vertex_material(4,face_num) = material_num
            face_order(face_num) = face_order(face_num) + 1
          end if

        end if
 
        cycle
 
      else
        go to 99
      end if
!
!  *MESH_NORMALS
!
    else if ( level_name(level) == '*MESH_NORMALS' ) then
 
      if ( word == '{' ) then
        go to 20
      else if ( word == '}' ) then
        level = nlbrack - nrbrack
      else if ( word == '*MESH_FACENORMAL' ) then

        call word_next_read ( line, word, done )
        call s_to_i4 ( word, iface, ierror, lchar )
        iface = iface + face_num_old + OFFSET
        ivert = 0

        call word_next_read ( line, word, done )
        call s_to_r4 ( word, x, ierror, lchar )

        call word_next_read ( line, word, done )
        call s_to_r4 ( word, y, ierror, lchar )

        call word_next_read ( line, word, done )
        call s_to_r4 ( word, z, ierror, lchar )

        face_normal(1,iface) = x
        face_normal(2,iface) = y
        face_normal(3,iface) = z

        cycle

      else if ( word == '*MESH_VERTEXNORMAL' ) then

        call word_next_read ( line, word, done )
        call s_to_i4 ( word, node, ierror, lchar )
        ivert = ivert + 1

        call word_next_read ( line, word, done )
        call s_to_r4 ( word, x, ierror, lchar )

        call word_next_read ( line, word, done )
        call s_to_r4 ( word, y, ierror, lchar )

        call word_next_read ( line, word, done )
        call s_to_r4 ( word, z, ierror, lchar )

        vertex_normal(1,ivert,iface) = x
        vertex_normal(2,ivert,iface) = y
        vertex_normal(3,ivert,iface) = z

        cycle

      else
        go to 99
      end if
!
!  *MESH_TFACELIST
!
    else if ( level_name(level) == '*MESH_TFACELIST' ) then
 
      if ( word == '{' ) then
        go to 20
      else if ( word == '}' ) then

        level = nlbrack - nrbrack

      else if ( word1 == '*MESH_TFACE' ) then
        cycle
      else
        go to 99
      end if
!
!  *MESH_TVERTLIST
!
    else if ( level_name(level) == '*MESH_TVERTLIST' ) then
 
      if ( word == '{' ) then
        go to 20
      else if ( word == '}' ) then
        level = nlbrack - nrbrack
      else if ( word1 == '*MESH_TVERT' ) then
        cycle
      else
        go to 99
      end if
!
!  *MESH_VERTEX_LIST
!
    else if ( level_name(level) == '*MESH_VERTEX_LIST' ) then
 
      if ( word == '{' ) then 
        cor3_num_old = cor3_num
        go to 20
       else if ( word == '}' ) then
        level = nlbrack - nrbrack
      else if ( word1 == '*MESH_VERTEX' ) then
 
        call word_next_read ( line, word, done ) 
        call s_to_i4 ( word, i, ierror, lchar )
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, x, ierror, lchar )
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, y, ierror, lchar )
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, z, ierror, lchar )
 
        i = i + cor3_num_old + OFFSET
        cor3_num = max ( cor3_num, i )
        if ( i <= cor3_max ) then
          cor3(1,i) =  transform_matrix(1,1) * x + transform_matrix(1,2) * y &
                     + transform_matrix(3,1) * z + transform_matrix(4,1)

          cor3(2,i) =  transform_matrix(2,1) * x + transform_matrix(2,2) * y &
                     + transform_matrix(2,3) * z + transform_matrix(2,4)

          cor3(3,i) =  transform_matrix(3,1) * x + transform_matrix(3,2) * y &
                     + transform_matrix(3,3) * z + transform_matrix(3,4)
        end if

        cycle

      else
        go to 99
      end if
!
!  *NODE_TM
!
!  Each node should start out with a default transformation matrix.
!
    else if ( level_name(level) == '*NODE_TM' ) then
 
      if ( word == '{' ) then

        call tmat_init ( transform_matrix )
        go to 20

      else if ( word == '}' ) then
        level = nlbrack - nrbrack
      else if ( word == '*INHERIT_POS' ) then
        cycle
      else if ( word == '*INHERIT_ROT' ) then
        cycle
      else if ( word == '*INHERIT_SCL' ) then
        cycle
      else if ( word == '*NODE_NAME' ) then
        cycle
      else if ( word == '*TM_POS' ) then
        cycle
      else if ( word == '*TM_ROTANGLE' ) then
        cycle
      else if ( word == '*TM_ROTAXIS' ) then
        cycle
      else if ( word == '*TM_ROW0' ) then
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, temp, ierror, lchar )
        transform_matrix(1,1) = temp
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, temp, ierror, lchar )
        transform_matrix(2,1) = temp
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, temp, ierror, lchar )
        transform_matrix(3,1) = temp
        cycle
      else if ( word == '*TM_ROW1' ) then
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, temp, ierror, lchar )
        transform_matrix(1,2) = temp
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, temp, ierror, lchar )
        transform_matrix(2,2) = temp
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, temp, ierror, lchar )
        transform_matrix(3,2) = temp
        cycle
      else if ( word == '*TM_ROW2' ) then
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, temp, ierror, lchar )
        transform_matrix(1,3) = temp
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, temp, ierror, lchar )
        transform_matrix(2,3) = temp
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, temp, ierror, lchar )
        transform_matrix(3,3) = temp
        cycle
      else if ( word == '*TM_ROW3' ) then
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, temp, ierror, lchar )
        transform_matrix(1,4) = temp
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, temp, ierror, lchar )
        transform_matrix(2,4) = temp
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, temp, ierror, lchar )
        transform_matrix(3,4) = temp
        cycle
      else if ( word == '*TM_SCALE' ) then
        cycle
      else if ( word == '*TM_SCALEAXIS' ) then
        cycle
      else if ( word == '*TM_SCALEAXISANG' ) then
        cycle
      else
        go to 99
      end if
!
!  *SCENE
!
    else if ( level_name(level) == '*SCENE' ) then

      if ( word == '{' ) then
        go to 20
      else if ( word == '}' ) then
        level = nlbrack - nrbrack
      else if ( word == '*SCENE_AMBIENT_STATIC' ) then
        cycle
      else if ( word == '*SCENE_BACKGROUND_STATIC' ) then
        cycle
      else if ( word == '*SCENE_FILENAME' ) then
        cycle
      else if ( word == '*SCENE_FIRSTFRAME' ) then
        cycle
      else if ( word == '*SCENE_FRAMESPEED' ) then
        cycle
      else if ( word == '*SCENE_LASTFRAME' ) then
        cycle
      else if ( word == '*SCENE_TICKSPERFRAME' ) then
        cycle
      else
        go to 99
      end if
 
    end if
 
    go to 20
!
!  Bad data
!
  99    continue

    bad_num = bad_num + 1

    if ( bad_num <= 10 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASE_READ - Warning!'
      write ( *, '(a)' ) '  Bad data on level ' // trim ( level_name(level) )
      write ( *, '(a)' ) '  Bad word: ' // trim ( word )
      write ( *, '(a,i6)' ) '  Line number: ', text_num
      write ( *, '(a)' ) trim ( line )
    end if

  end do
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'ASE_READ - Read ', text_num, ' text lines from ' &
    // trim ( filein_name )

  return
end
subroutine ase_write ( cor3, cor3_max, cor3_num, face, face_max, &
  face_normal, face_num, face_order, filein_name, fileout_name, &
  iunit, order_max, vertex_normal )

!*****************************************************************************80
!
!! ASE_WRITE writes graphics information to an ASE file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, real ( kind = 4 ) FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  character ( len = 10 ) chrtmp
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) face(order_max,face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  character ( len = * ) fileout_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: OFFSET = 1
  integer ( kind = 4 ) text_num
  character ( len = 200 ) text
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  real ( kind = 4 ) z

  text_num = 0
!
!  Write the header.
!
  write ( iunit, '(a)' ) '*3DSMAX_ASCIIEXPORT 200'
  write ( iunit, '(a)' ) '*COMMENT "' // trim ( fileout_name ) // &
    ', created by IVREAD."'
  write ( iunit, '(a)' ) '*COMMENT "Original data in ' // trim ( filein_name ) &
    // '."'

  text_num = text_num + 3
!
!  Write the scene block.
!
  write ( iunit, '(a)' ) '*SCENE {'
  write ( iunit, '(a)' ) '  *SCENE_FILENAME ""'
  write ( iunit, '(a)' ) '  *SCENE_FIRSTFRAME 0'
  write ( iunit, '(a)' ) '  *SCENE_LASTFRAME 100'
  write ( iunit, '(a)' ) '  *SCENE_FRAMESPEED 30'
  write ( iunit, '(a)' ) '  *SCENE_TICKSPERFRAME 160'
  write ( iunit, '(a)' ) '  *SCENE_BACKGROUND_STATIC 0.0000 0.0000 0.0000'
  write ( iunit, '(a)' ) '  *SCENE_AMBIENT_STATIC 0.0431 0.0431 0.0431'
  write ( iunit, '(a)' ) '}'

  text_num = text_num + 9
!
!  Begin the big geometry block.
!
  write ( iunit, '(a)' ) '*GEOMOBJECT {'
  write ( iunit, '(a)' ) '  *NODE_NAME "Object01"'

  text_num = text_num + 2
!
!  Sub block NODE_TM:
!
  write ( iunit, '(a)' ) '  *NODE_TM {'
  write ( iunit, '(a)' ) '    *NODE_NAME "Object01"'
  write ( iunit, '(a)' ) '    *INHERIT_POS 0 0 0'
  write ( iunit, '(a)' ) '    *INHERIT_ROT 0 0 0'
  write ( iunit, '(a)' ) '    *INHERIT_SCL 0 0 0'
  write ( iunit, '(a)' ) '    *TM_ROW0 1.0000 0.0000 0.0000'
  write ( iunit, '(a)' ) '    *TM_ROW1 0.0000 1.0000 0.0000'
  write ( iunit, '(a)' ) '    *TM_ROW2 0.0000 0.0000 1.0000'
  write ( iunit, '(a)' ) '    *TM_ROW3 0.0000 0.0000 0.0000'
  write ( iunit, '(a)' ) '    *TM_POS 0.0000 0.0000 0.0000'
  write ( iunit, '(a)' ) '    *TM_ROTAXIS 0.0000 0.0000 0.0000'
  write ( iunit, '(a)' ) '    *TM_ROTANGLE 0.0000'
  write ( iunit, '(a)' ) '    *TM_SCALE 1.0000 1.0000 1.0000'
  write ( iunit, '(a)' ) '    *TM_SCALEAXIS 0.0000 0.0000 0.0000'
  write ( iunit, '(a)' ) '    *TM_SCALEAXISANG 0.0000'
  write ( iunit, '(a)' ) '  }'

  text_num = text_num + 16
!
!  Sub block MESH:
!    Items
!
  write ( iunit, '(a)' ) '  *MESH {'
  write ( iunit, '(a)' ) '    *TIMEVALUE 0'
  write ( chrtmp, '(i8)' ) cor3_num
  write ( iunit, '(a)' ) '    *MESH_NUMVERTEX ' // trim ( chrtmp )
  write ( chrtmp, '(i8)' ) face_num
  write ( iunit, '(a)' ) '    *MESH_NUMFACES ' // trim ( chrtmp )

  text_num = text_num + 4
!
!  Sub sub block MESH_VERTEX_LIST
!
  write ( iunit, '(a)' ) '    *MESH_VERTEX_LIST {'

  do j = 1, cor3_num
    write ( text, '(a,i8,3g12.4)' ) '*MESH_VERTEX ', j - OFFSET, cor3(1:3,j)
    call s_blanks_delete ( text )
    write ( iunit, '(6x,a)' ) trim ( text )
  end do

  write ( iunit, '(a)' ) '    }'

  text_num = text_num + cor3_num + 2
!
!  Sub sub block MESH_FACE_LIST
!    Items MESH_FACE
!
!  ???  What do you do when the face has more than 4 vertices?
!
  write ( iunit, '(a)' ) '    *MESH_FACE_LIST {'

  do iface = 1, face_num

    i1 = face(1,iface) - OFFSET
    i2 = face(2,iface) - OFFSET
    i3 = face(3,iface) - OFFSET

    if ( face_order(iface) == 3 ) then

      write ( text, '(a,i8,a,i8,a,i8,a,i8,a)' )  '*MESH_FACE ', &
        iface - OFFSET, ': A: ', i1, ' B: ', i2, ' C: ', i3, &
        ' AB: 1 BC: 1 CA: 1 *MESH_SMOOTHING *MESH_MTLID 1'

      call s_blanks_delete ( text )
      write ( iunit, '(6x,a)' ) trim ( text )

    else if ( face_order(iface) == 4 ) then

      i4 = face(4,iface) - OFFSET
      write ( text, '(a,i8,a,i8,a,i8,a,i8,a,i8,a)' ) '*MESH_FACE ', &
        iface - OFFSET, ': A: ', i1, ' B: ', i2, ' C: ', i3, ' D: ', i4, &
        ' AB: 1 BC: 1 CD: 1 DA: 1 *MESH_SMOOTHING *MESH_MTLID 1'
      call s_blanks_delete ( text )
      write ( iunit, '(6x,a)' ) trim ( text )

    end if

  end do

  write ( iunit, '(a)' ) '    }'

  text_num = text_num + face_num + 2
!
!  Item MESH_NUMTVERTEX
!
  write ( iunit, '(a)' ) '    *MESH_NUMTVERTEX 0'
  text_num = text_num + 1
!
!  Item NUMCVERTEX
!
  write ( iunit, '(a)' ) '    *MESH_NUMCVERTEX 0'
  text_num = text_num + 1
!
!  Sub block MESH_NORMALS
!    Items MESH_FACENORMAL, MESH_VERTEXNORMAL (repeated)
!
  write ( iunit, '(a)' ) '    *MESH_NORMALS {'
  text_num = text_num + 1

  do iface = 1, face_num

    x = face_normal(1,iface)
    y = face_normal(2,iface)
    z = face_normal(3,iface)

    write ( text, '(a,i8,3g12.4)' ) '*MESH_FACENORMAL ', iface-OFFSET, x, y, z
    call s_blanks_delete ( text )

    write ( iunit, '(6x,a)' ) trim ( text )
    text_num = text_num + 1

    do ivert = 1, face_order(iface)

      x = vertex_normal(1,ivert,iface)
      y = vertex_normal(2,ivert,iface)
      z = vertex_normal(3,ivert,iface)

      write ( text, '(a,i8,3g12.4)' ) '*MESH_VERTEXNORMAL ', ivert-OFFSET, &
        x, y, z 
      call s_blanks_delete ( text )
      write ( iunit, '(6x,a)' ) trim ( text )
      text_num = text_num + 1

    end do

  end do

  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 1
!
!  Close the MESH object.
!
  write ( iunit, '(a)' ) '  }'
!
!  A few closing parameters.
!
  write ( iunit, '(a)' ) '  *PROP_MOTIONBLUR 0'
  write ( iunit, '(a)' ) '  *PROP_CASTSHADOW 1'
  write ( iunit, '(a)' ) '  *PROP_RECVSHADOW 1'
!
!  Close the GEOM object.
!
  write ( iunit, '(a)' ) '}'

  text_num = text_num + 5
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'ASE_WRITE - Wrote ', text_num, ' text lines to ' &
    // trim ( fileout_name )

  return
end
subroutine byu_read ( cor3, cor3_max, cor3_num, face, face_max, face_num, &
  face_order, filein_name, ierror, iunit, order_max )

!*****************************************************************************80
!
!! BYU_READ reads graphics data from a Movie.BYU surface geometry file.
!
!  Discussion:
!
!    A Movie.BYU surface geometry file contains 4 groups of data.
!
!    The first group of data is a single line, containing 4 integers,
!    each one left justified in 8 columns.  The integers are:
!
!      PART_NUM, VERTEX_NUM, POLY_NUM, EDGE_NUM,
! 
!    that is, the number of parts or objects, the number of vertices or nodes,
!    the number of polygons or faces, and the number of edges.
!
!    The second group of data is a single line, containing 2 integers,
!    each one left justified in 8 columnes.  The integers are:
!
!      POLY1, POLY2,
!
!    the starting and ending polygon numbers.  Presumably, this means
!    that the polygons are labeled POLY1, POLY1+1, ..., POLY2, comprising
!    a total of POLY_NUM polygons.
!
!    The third group is the X, Y and Z coordinates of all the vertices.
!    These may be written using a FORTRAN format of 6E12.5, which
!    crams two sets of (X,Y,Z) data onto each line, with each real value
!    written in an exponential format with 5 places after the decimal.
!    However, it is generally possible to write the XYZ coordinate data
!    for each vertex on a separate line.
!
!    The fourth group defines the polygons in terms of the vertex indices.
!    For each polygon, the vertices that make up the polygon are listed in
!    counterclockwise order.  The last vertex listed is given with a negative
!    sign to indicate the end of the list.  All the vertices for all the
!    polygons are listed one after the other, using a format that puts
!    up to 10 left-justified integers on a line, with each integer occupying
!    8 spaces.
!
!    This code will certainly read a BYU file created by BYU_WRITE, but
!    it will not handle more general files.  In particular, an object
!    can have several parts, the coordinate data can be grouped so
!    that there are 2 sets of (x,y,z) data per line, and so on.
!
!  Example:
!
!          1       8       6      24
!          1       6
!    0.00000E+00 0.00000E+00 0.00000E+00
!    1.00000E+00 0.00000E+00 0.00000E+00
!    1.00000E+00 2.00000E+00 0.00000E+00
!    0.00000E+00 2.00000E+00 0.00000E+00
!    0.00000E+00 0.00000E+00 1.00000E+00
!    1.00000E+00 0.00000E+00 1.00000E+00
!    1.00000E+00 2.00000E+00 1.00000E+00
!    0.00000E+00 2.00000E+00 1.00000E+00
!          4       3       2      -1
!          5       6       7      -8
!          1       5       8      -4
!          4       8       7      -3
!          3       7       6      -2
!          2       6       5      -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input/output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices 
!    per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which data is read.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input/output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) cor3_num_new
  logical done
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_num_new
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) np2
  integer ( kind = 4 ) part_num
  character ( len = 256 ) text
  integer ( kind = 4 ) text_num

  ierror = 0
  text_num = 0

  read ( iunit, *, iostat = ios ) part_num, cor3_num_new, face_num_new, edge_num

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BYU_READ - Fatal error!'
    write ( *, '(a)' ) '  Unexpected end of file.'
    return
  end if

  text_num = text_num + 1

  read ( iunit, *, iostat = ios ) np1, np2

  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BYU_READ - Fatal error!'
    write ( *, '(a)' ) '  Unexpected end of file.'
    return
  end if

  text_num = text_num + 1

  do j = cor3_num + 1, cor3_num + cor3_num_new

    read ( iunit, *, iostat = ios ) cor3(1:3,j)

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BYU_READ - Fatal error!'
      write ( *, '(a)' ) '  Unexpected end of file.'
      return
    end if

    text_num = text_num + 1

  end do

  do iface = face_num + 1, face_num + face_num_new

    read ( iunit, '(a)', iostat = ios ) text
    text_num = text_num + 1

    ivert = 0
    done = .true.

    do

      call intnex ( text, ival, done )

      if ( done ) then
        exit
      end if

      ivert = ivert + 1
      face(ivert,iface) = abs ( ival ) + cor3_num

      if ( ival <= 0 ) then
        exit
      end if

    end do

    face_order(iface) = ivert

  end do

  cor3_num = cor3_num + cor3_num_new
  face_num = face_num + face_num_new
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'BYU_READ - Read ', text_num, ' text lines from ' &
    // trim ( filein_name )

  return
end
subroutine byu_write ( cor3, cor3_max, cor3_num, face, face_max, face_num, &
  face_order, fileout_name, iunit, order_max )

!*****************************************************************************80
!
!! BYU_WRITE writes out the graphics data as a Movie.BYU surface geometry file.
!
!  Discussion:
!
!    A Movie.BYU surface geometry file contains 4 groups of data.
!
!    The first group of data is a single line, containing 4 integers,
!    each one left justified in 8 columns.  The integers are:
!
!      PART_NUM, VERTEX_NUM, POLY_NUM, EDGE_NUM,
! 
!    that is, the number of parts or objects, the number of vertices or nodes,
!    the number of polygons or faces, and the number of edges.
!
!    The second group of data is a single line, containing 2 integers,
!    each one left justified in 8 columnes.  The integers are:
!
!      POLY1, POLY2,
!
!    the starting and ending polygon numbers.  Presumably, this means
!    that the polygons are labeled POLY1, POLY1+1, ..., POLY2, comprising
!    a total of POLY_NUM polygons.
!
!    The third group is the X, Y and Z coordinates of all the vertices.
!    These may be written using a FORTRAN format of 6E12.5, which
!    crams two sets of (X,Y,Z) data onto each line, with each real value
!    written in an exponential format with 5 places after the decimal.
!    However, it is generally possible to write the XYZ coordinate data
!    for each vertex on a separate line.
!
!    The fourth group defines the polygons in terms of the vertex indices.
!    For each polygon, the vertices that make up the polygon are listed in
!    counterclockwise order.  The last vertex listed is given with a negative
!    sign to indicate the end of the list.  All the vertices for all the
!    polygons are listed one after the other, using a format that puts
!    up to 10 left-justified integers on a line, with each integer occupying
!    8 spaces.
!
!    This code will certainly read a BYU file created by BYU_WRITE, but
!    it will not handle more general files.  In particular, an object
!    can have several parts, the coordinate data can be grouped so
!    that there are 2 sets of (x,y,z) data per line, and so on.
!
!  Example:
!
!          1       8       6      24
!          1       6
!    0.00000E+00 0.00000E+00 0.00000E+00
!    1.00000E+00 0.00000E+00 0.00000E+00
!    1.00000E+00 2.00000E+00 0.00000E+00
!    0.00000E+00 2.00000E+00 0.00000E+00
!    0.00000E+00 0.00000E+00 1.00000E+00
!    1.00000E+00 0.00000E+00 1.00000E+00
!    1.00000E+00 2.00000E+00 1.00000E+00
!    0.00000E+00 2.00000E+00 1.00000E+00
!          4       3       2      -1
!          5       6       7      -8
!          1       5       8      -4
!          4       8       7      -3
!          3       7       6      -2
!          2       6       5      -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which data is read.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) fileout_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp(8)
  integer ( kind = 4 ) part_num
  integer ( kind = 4 ) text_num

  text_num = 0
!
!  EDGE_NUM is the total number of edges.
!
  edge_num = sum ( face_order(1:face_num) )

  part_num = 1
  write ( iunit, '(10i8)' ) part_num, cor3_num, face_num, edge_num
  text_num = text_num + 1

  write ( iunit, '(10i8)' ) 1, face_num
  text_num = text_num + 1

  do j = 1, cor3_num
    write ( iunit, '(1p6e12.5)' ) cor3(1:3,j)
    text_num = text_num + 1
  end do
!
!  It takes a little mangling in order to print out all the edges in a
!  single list, one face at a time, with the last node in each face
!  negative, and written in groups of 8.
!
  ihi = 0

  do iface = 1, face_num

    do ivert = 1, face_order(iface)

      ihi = ihi + 1
      jp(ihi) = face(ivert,iface)
      if ( ivert == face_order(iface) ) then
        jp(ihi) = - jp(ihi)
      end if

      if ( ihi == 8 .or. ivert == face_order(iface) ) then
        write ( iunit, '(10i8)' ) jp(1:ihi)
        text_num = text_num + 1
        ihi = 0
      end if

    end do
  end do
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'BYU_WRITE - Wrote ', text_num, ' text lines to ' &
    // trim ( fileout_name )

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
!    14 August 1999
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
  character c2
  character cc1
  character cc2

  cc1 = c1
  cc2 = c2

  call ch_cap ( cc1 )
  call ch_cap ( cc2 )

  if ( cc1 == cc2 ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
function ch_is_control ( c )

!*****************************************************************************80
!
!! CH_IS_CONTROL reports whether a character is a control character or not.
!
!  Discussion:
!
!    A "control character" has ASCII code <= 31 or 127 <= ASCII code.
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
!    Input, character C, the character to be tested.
!
!    Output, logical CH_IS_CONTROL, TRUE if C is a control character, and
!    FALSE otherwise.
!
  implicit none

  character c
  logical ch_is_control

  if ( ichar ( c ) <= 31 .or. 127 <= ichar ( c ) ) then
    ch_is_control = .true.
  else
    ch_is_control = .false.
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
subroutine command_line ( arg_num, cor3, cor3_material, cor3_max, cor3_normal, &
  cor3_tex_uv, debug, face, face_area, face_material, face_max, face_normal, &
  face_order, face_tex_uv, filein_name, fileout_name, ierror, &
  line_dex, line_material, line_prune, material_name, material_rgba, &
  line_max, material_max, order_max, texture_max, &
  object_name, point, point_max, point_num, texture_name, &
  texture_temp, transform_matrix, vertex_material, vertex_normal, &
  vertex_tex_uv )

!*****************************************************************************80
!
!! COMMAND_LINE works with command line parameters.
!
!  Discussion:
!
!    This routine is invoked when the user command is something like
!
!      ivread filein_name fileout_name
!
!    or
!
!      ivread -rn filein_name fileout_name
!
!    where "-rn" signals the "reverse normals" option, or
!
!      ivread -rf filein_name fileout_name
!
!    where "-rf" signals the "reverse faces" option.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ARG_NUM, the number of command-line arguments.
!
!    Workspace, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Workspace, integer COR3_MATERIAL(COR3_MAX), the material of each node.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Workspace, real COR3_NORMAL(3,COR3_MAX), normals at nodes.
!
!    Workspace, real COR3_TEX_UV(2,COR3_MAX), UV texture coordinates for nodes.
!
!    Input, logical DEBUG, debugging switch.
!
!    Workspace, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Workspace, real FACE_AREA(FACE_MAX), area of each face.
!
!    Workspace, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Workspace, real FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Workspace, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Workspace, real FACE_TEX_UV(2,FACE_MAX), face texture coordinates.
!
!    Workspace, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Workspace, character ( len = * ) FILEOUT_NAME, the name of the 
!    output file.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!
!    Workspace, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated 
!    by -1.
!
!    Workspace, integer LINE_MATERIAL(LINE_MAX), material index for line.
!
!    Input, integer ( kind = 4 ) LINE_PRUNE, pruning option.
!    0, no pruning, draw every line.
!    nonzero, prune.  Only draw the line from node I to node J if I < J.
!    This should cut down on repeated drawing of lines in the common
!    case of a face mesh where each line is drawn twice, once with positive
!    and once with negative orientation.  In other cases, pruning
!    may omit some lines that only occur with negative orientation.
!
!    Workspace, character ( len = * ) MATERIAL_NAME(MATERIAL_MAX), 
!    material names.
!
!    Workspace, real MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and
!    A values.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer ( kind = 4 ) TEXTURE_MAX, the maximum number of textures.
!
!    Workspace, character ( len = * ) OBJECT_NAME, object name.
!
!    Workspace, character ( len = * ) TEXTURE_NAME(TEXTURE_MAX), 
!    texture names.
!
!    Workspace, real TRANSFORM_MATRIX(4,4), the transformation matrix.
!
!    Workspace, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Workspace, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
!    Workspace, real VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture
!    coordinates.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) point_max
  integer ( kind = 4 ) texture_max

  integer ( kind = 4 ) arg_num
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_material(cor3_max)
  real ( kind = 4 ) cor3_normal(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  real ( kind = 4 ) cor3_tex_uv(2,cor3_max)
  logical debug
  integer ( kind = 4 ) face(order_max,face_max)
  real ( kind = 4 ) face_area(face_max)
  integer ( kind = 4 ) face_material(face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  real ( kind = 4 ) face_tex_uv(2,face_max)
  character ( len = * ) filein_name
  character ( len = 10 ) filein_type
  character ( len = * ) fileout_name
  character ( len = 10 ) fileout_type
  integer ( kind = 4 ) group_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) ilen
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) line_prune
  character ( len = * ) material_name(material_max)
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  character ( len = * ) object_name
  integer ( kind = 4 ) object_num
  integer ( kind = 4 ) point(point_max)
  integer ( kind = 4 ) point_num
  logical reverse_faces
  logical reverse_normals
  logical s_eqi
  character ( len = * ) texture_name(texture_max)
  integer ( kind = 4 ) texture_num
  real ( kind = 4 ) texture_temp(2,order_max*face_max)
  real ( kind = 4 ) transform_matrix(4,4)
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  real ( kind = 4 ) vertex_tex_uv(2,order_max,face_max)

  reverse_faces = .false.
  reverse_normals = .false.
!
!  Initialize the graphics data.
!
  call data_init ( cor3, cor3_material, cor3_max, cor3_normal, &
    cor3_num,  cor3_tex_uv, face, face_area, face_material, face_max, &
    face_normal, face_num, face_order, face_tex_uv, group_num, line_dex, &
    line_material, line_max, line_num, material_max, material_name, &
    material_num, material_rgba, object_name, object_num, order_max, &
    point, point_max, point_num, texture_max, texture_name, texture_num, &
    texture_temp, transform_matrix, vertex_material, vertex_normal, &
    vertex_tex_uv )
!
!  Sort out the command line arguments.
!
  iarg = 1
!
!  Old style:
!
  call getarg ( iarg, filein_name )
!
!  New style:
!
! call pxfgetarg ( iarg, filein_name, ilen, ierror )
!
! if ( ierror /= 0 ) then
!   write ( *, '(a)' ) ' '
!   write ( *, '(a)' ) 'COMMAND_LINE - Fatal error!'
!   write ( *, '(a)' ) '  Could not read command line argument.'
!   stop
! end if

  if ( s_eqi ( filein_name, '-RN' ) ) then

    reverse_normals = .true.
    iarg = iarg + 1
!
!  Old style:
!
    call getarg ( iarg, filein_name )
!
!  New style:
!
!   call pxfgetarg ( iarg, filein_name, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'COMMAND_LINE - Fatal error!'
!     write ( *, '(a)' ) '  Could not read command line argument.'
!    stop
! end  if

  else if ( s_eqi ( filein_name, '-RF' ) ) then

    reverse_faces = .true.
    iarg = iarg + 1
!
!  Old style:
!
    call getarg ( iarg, filein_name )
!
!  New style:
!
!   call pxfgetarg ( iarg, filein_name, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'COMMAND_LINE - Fatal error!'
!     write ( *, '(a)' ) '  Could not read command line argument.'
!     stop
!   end if

  end if
 
  iarg = iarg + 1

  if ( iarg <= arg_num ) then
!
!  Old style:
!
    call getarg ( iarg, fileout_name )
!
!  New style:
!
!   call pxfgetarg ( iarg, fileout_name, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'COMMAND_LINE - Fatal error!'
!     write ( *, '(a)' ) '  Could not read command line argument.'
!     stop
!   end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMMAND_LINE - Fatal error!'
    write ( *, '(a)' ) '  Not enough arguments on the command line.'
    write ( *, '(a)' ) '  Could not find the output file name.'
    return

  end if
!
!  Check the input file name.
!
  call infile ( filein_name, ierror, filein_type )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMMAND_LINE - Fatal error!'
    write ( *, '(a)' ) '  Improper input file name.'
    return
  end if
!
!  Read the input file.
!
  call data_read ( cor3, cor3_material, cor3_max, &
    cor3_normal, cor3_num, cor3_tex_uv, debug, face, face_area, &
    face_material, face_max, face_normal, face_num, face_order, &
    face_tex_uv, filein_name, filein_type, group_num, ierror, line_dex, &
    line_material, line_max, line_num, material_max, material_name, &
    material_num, material_rgba, object_num, order_max, &
    point, point_max, point_num, texture_max, texture_name, texture_num, &
    texture_temp, vertex_material, vertex_normal, vertex_tex_uv )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMMAND_LINE - Fatal error!'
    write ( *, '(a)' ) '  The input file could not be opened or read properly.'
    return
  end if
!
!  Check the output file name.
!
  call outfile ( filein_name, fileout_name, ierror, fileout_type )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMMAND_LINE - Fatal error!'
    write ( *, '(a)' ) '  Improper output file name.'
    return
  end if
!
!  -RN: Reverse the normal vectors at points, vertices, and faces.
!
  if ( reverse_normals ) then

    cor3_normal(1:3,1:cor3_num) = - cor3_normal(1:3,1:cor3_num)
    vertex_normal(1:3,1:order_max,1:face_num) = &
      - vertex_normal(1:3,1:order_max,1:face_num)

    face_normal(1:3,1:face_num) = - face_normal(1:3,1:face_num)
!
!  -RF: Reverse the order of the nodes defining each face.
!
  else if ( reverse_faces ) then
 
    call face_reverse_order ( cor3_max, cor3_normal, cor3_num, face, &
      face_max, face_normal, face_num, face_order, order_max, &
      vertex_material, vertex_normal, vertex_tex_uv )

  end if
!
!  Write the output file.
!
  call data_write ( cor3, cor3_material, cor3_max, cor3_normal, cor3_num, &
    cor3_tex_uv, debug, face, face_material, face_max, face_normal, &
    face_num, face_order, face_tex_uv, filein_name, fileout_name, &
    fileout_type, ierror, line_dex, line_material, line_max, line_num, &
    line_prune, material_name, material_max, material_num, material_rgba, &
    object_name, order_max, point, point_max, point_num, texture_max, &
    texture_name, texture_num, vertex_material, vertex_normal, &
    vertex_tex_uv )

  return
end
subroutine cor3_normal_set ( cor3_max, cor3_normal, face, face_area, &
  face_max, face_num, face_order, order_max, vertex_normal )

!*****************************************************************************80
!
!! COR3_NORMAL_SET recomputes zero node normal vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input/output, real ( kind = 4 ) COR3_NORMAL(3,COR3_MAX), normals at nodes.  
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3_normal(3,cor3_max)
  real ( kind = 4 ) cor3_temp(cor3_max)
  integer ( kind = 4 ) face(order_max,face_max)
  real ( kind = 4 ) face_area(face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  real ( kind = 4 ) norm
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
!
!  Determine which nodes need to have their normals recomputed.
!
  do icor3 = 1, cor3_max
    cor3_temp(icor3) = sum ( cor3_normal(1:3,icor3)**2 )
  end do
!
!  Add up the vertex normals from all the faces to which the node belongs,
!  weighted by area.
!
  do iface = 1, face_num

    do ivert = 1, face_order(iface)

      icor3 = face(ivert,iface)

      if ( cor3_temp(icor3) == 0.0E+00 ) then

        do j = 1, 3
          cor3_normal(j,icor3) = cor3_normal(j,icor3) &
               + face_area(iface) * vertex_normal(j,ivert,iface)
        end do

      end if

    end do

  end do
!
!  Renormalize.
!
  do icor3 = 1, cor3_max

    norm = sum ( cor3_normal(1:3,icor3)**2 )

    if ( norm == 0.0E+00 ) then

      norm = 3.0E+00
      cor3_normal(1:3,icor3) = 1.0E+00

    end if

    norm = sqrt ( norm )

    cor3_normal(1:3,icor3) = cor3_normal(1:3,icor3) / norm

  end do

  return
end
subroutine cor3_range ( cor3, cor3_max, cor3_num )

!*****************************************************************************80
!
!! COR3_RANGE computes and prints the coordinate minima and maxima.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
  implicit none

  integer ( kind = 4 ) cor3_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) i
  real ( kind = 4 ) xave
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin
  real ( kind = 4 ) yave
  real ( kind = 4 ) ymax
  real ( kind = 4 ) ymin
  real ( kind = 4 ) zave
  real ( kind = 4 ) zmax
  real ( kind = 4 ) zmin

  if ( cor3_num <= 0 ) then

    xave = 0.0E+00
    xmax = 0.0E+00
    xmin = 0.0E+00

    yave = 0.0E+00
    ymax = 0.0E+00
    ymin = 0.0E+00

    zave = 0.0E+00
    zmax = 0.0E+00
    zmin = 0.0E+00

  else

    xave = cor3(1,1)
    xmax = cor3(1,1)
    xmin = cor3(1,1)

    yave = cor3(2,1)
    ymax = cor3(2,1)
    ymin = cor3(2,1)

    zave = cor3(3,1)
    zmax = cor3(3,1)
    zmin = cor3(3,1)

    do i = 2, cor3_num

      xave = xave + cor3(1,i)
      xmin = min ( xmin, cor3(1,i) )
      xmax = max ( xmax, cor3(1,i) )
  
      yave = yave + cor3(2,i)
      ymin = min ( ymin, cor3(2,i) )
      ymax = max ( ymax, cor3(2,i) )

      zave = zave + cor3(3,i)
      zmin = min ( zmin, cor3(3,i) )
      zmax = max ( zmax, cor3(3,i) )

    end do

    xave = xave / real ( cor3_num )
    yave = yave / real ( cor3_num )
    zave = zave / real ( cor3_num )

  end if

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COR3_RANGE - Nodal coordinate range:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     Minimum     Average    Maximum       Range'
    write ( *, '(a)' ) ' '
    write ( *, '(a1,4g12.4)' ) 'X', xmin, xave, xmax, xmax-xmin
    write ( *, '(a1,4g12.4)' ) 'Y', ymin, yave, ymax, ymax-ymin
    write ( *, '(a1,4g12.4)' ) 'Z', zmin, zave, zmax, zmax-zmin
  end if
 
  return
end
subroutine cross0_3d ( x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3 )

!*****************************************************************************80
!
!! CROSS0_3D computes the cross product of (P1-P0) and (P2-P0) in 3D.
!
!  Discussion:
!
!    The vectors are specified with respect to a basis point P0.
!    We are computing the normal to the triangle (P0,P1,P2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2, the coordinates of
!    three points.  The basis point is (X0,Y0,Z0).
!
!    Output, real ( kind = 4 ) X3, Y3, Z3, the cross product (P1-P0) x (P2-P0), or
!    (X1-X0,Y1-Y0,Z1-Z0) x (X2-X0,Y2-Y0,Z2-Z0).
!
  implicit none

  real ( kind = 4 ) x0
  real ( kind = 4 ) x1
  real ( kind = 4 ) x2
  real ( kind = 4 ) x3
  real ( kind = 4 ) y0
  real ( kind = 4 ) y1
  real ( kind = 4 ) y2
  real ( kind = 4 ) y3
  real ( kind = 4 ) z0
  real ( kind = 4 ) z1
  real ( kind = 4 ) z2
  real ( kind = 4 ) z3

  x3 = ( y1 - y0 ) * ( z2 - z0 ) - ( z1 - z0 ) * ( y2 - y0 )
  y3 = ( z1 - z0 ) * ( x2 - x0 ) - ( x1 - x0 ) * ( z2 - z0 )
  z3 = ( x1 - x0 ) * ( y2 - y0 ) - ( y1 - y0 ) * ( x2 - x0 )

  return
end
subroutine data_check ( cor3_max, cor3_num, face_max, face_num, face_order, &
  line_max, line_num, material_max, material_name, material_num, order_max, &
  texture_max, texture_num, texture_name )

!*****************************************************************************80
!
!! DATA_CHECK checks the input data, and enforces limits.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input/output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per
!    face.
!
!    Input/output, character ( len = * ) MATERIAL_NAME(MATERIAL_MAX), 
!    material names.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer ( kind = 4 ) TEXTURE_MAX, the maximum number of textures.
!
!    Input/output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input/output, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Input/output, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Input/output, integer ( kind = 4 ) TEXTURE_NUM, the number of textures.
!
!    Input/output, character ( len = * ) TEXTURE_NAME(TEXTURE_MAX), 
!    texture names.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) texture_max

  character ( len = 4 ) char4
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) material_num
  character ( len = * ) material_name(material_max)
  integer ( kind = 4 ) nfix
  character ( len = * ) texture_name(texture_max)
  integer ( kind = 4 ) texture_num

  if ( cor3_max < cor3_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_CHECK - Warning!'
    write ( *, '(a,i6,a)' ) '  The input data requires ', cor3_num, ' points.'
    write ( *, '(a,i6)' ) '  There was only room for ', cor3_max
    cor3_num = cor3_max
  else if ( cor3_num < 0 ) then
    cor3_num = 0
  end if
 
  if ( face_max < face_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_CHECK - Warning!'
    write ( *, '(a,i6,a)' ) '  The input data requires ', face_num, ' faces.'
    write ( *, '(a,i6)' ) '  There was only room for ', face_max
    face_num = face_max
  else if ( face_num <0 ) then
    face_num = 0
  end if

  if ( line_max < line_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_CHECK - Warning!'
    write ( *, '(a,i6,a)' ) '  The input data requires ', line_num, ' line items.'
    write ( *, '(a,i6)' ) '  There was only room for ', line_max
    line_num = line_max
  else if ( line_num < 0 ) then 
    line_num = 0
  end if

  if ( material_max < material_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_CHECK - Warning!'
    write ( *, '(a,i6,a)' ) '  The input data requires ', material_num, ' materials.'
    write ( *, '(a,i6)' ) '  There was only room for ', material_max
    material_num = material_max
  else if ( material_num < 0 ) then 
    material_num = 0
  end if

  nfix = 0

  do iface = 1, face_num

    if ( face_order(iface) < 3 ) then
      nfix = nfix + 1
      face_order(iface) = 0
    end if

  end do

  if ( 0 < nfix ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_CHECK - Warning!'
    write ( *, '(a,i6,a)' ) '  Corrected ', nfix, &
      ' faces using less than 3 vertices per face.'
  end if
  nfix = 0

  do iface = 1, face_num

    if ( order_max < face_order(iface) ) then
      nfix = nfix + 1
      face_order(iface) = order_max
    end if

  end do

  if ( 0 < nfix ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_CHECK - Warning!'
    write ( *, '(a,i6,a,i6,a)' ) '  Corrected ', nfix, &
      ' faces using more than ', order_max, ' vertices per face.'
  end if

  do i = 1, material_num
    if ( material_name(i) == ' ' ) then
      call i4_to_s_zero ( i, char4 )
      texture_name(i) = 'Material_' // char4
    end if
  end do

  do i = 1, texture_num
    if ( texture_name(i) == ' ' ) then
      call i4_to_s_zero ( i, char4 )
      texture_name(i) = 'Texture_' // char4
    end if
  end do

  return
end
subroutine data_init ( cor3, cor3_material, cor3_max, cor3_normal, &
  cor3_num,  cor3_tex_uv, face, face_area, face_material, face_max, &
  face_normal, face_num, face_order, face_tex_uv, group_num, line_dex, &
  line_material, line_max, line_num, material_max, material_name, &
  material_num, material_rgba, object_name, object_num, order_max, &
  point, point_max, point_num, texture_max, texture_name, texture_num, &
  texture_temp, transform_matrix, vertex_material, vertex_normal, &
  vertex_tex_uv )

!*****************************************************************************80
!
!! DATA_INIT initializes internal graphics data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Output, integer ( kind = 4 ) COR3_MATERIAL(COR3_MAX), the material of each node.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Output, real ( kind = 4 ) COR3_TEX_UV(2,COR3_MAX), UV texture coordinates for nodes.
!
!    Output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Output, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Output, real ( kind = 4 ) FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Output, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Output, integer ( kind = 4 ) LINE_MATERIAL(LINE_MAX), material index for line.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of line data items.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Output, character ( len = * ) MATERIAL_NAME(MATERIAL_MAX), 
!    material names.
!
!    Output, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Output, character ( len = * ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
!    Output, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Output, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
!    Output, real ( kind = 4 ) VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture
!    coordinates.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) point_max
  integer ( kind = 4 ) texture_max

  character ( len = 4 ) char4
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_material(cor3_max)
  real ( kind = 4 ) cor3_normal(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  real ( kind = 4 ) cor3_tex_uv(2,cor3_max)
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) face(order_max,face_max)
  real ( kind = 4 ) face_area(face_max)
  integer ( kind = 4 ) face_material(face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  real ( kind = 4 ) face_tex_uv(2,face_max)
  integer ( kind = 4 ) group_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_num
  character ( len = * ) material_name(material_max)
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  character ( len = * ) object_name
  integer ( kind = 4 ) object_num
  integer ( kind = 4 ), parameter :: OFFSET = 1
  integer ( kind = 4 ) point(point_max)
  integer ( kind = 4 ) point_num
  character ( len = * ) texture_name(texture_max)
  integer ( kind = 4 ) texture_num
  real ( kind = 4 ) texture_temp(2,order_max*face_max)
  real ( kind = 4 ) transform_matrix(4,4)
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  real ( kind = 4 ) vertex_tex_uv(2,order_max,face_max)

  cor3(1:3,1:cor3_max) = 0.0E+00
  cor3_material(1:cor3_max) = -1
  cor3_normal(1:3,1:cor3_max) = 0.0E+00
  cor3_num = 0
  cor3_tex_uv(1:2,1:cor3_max) = 0.0E+00
  face(1:order_max,1:face_max) = -1 + OFFSET
  face_area(1:face_max) = 0.0E+00
  face_material(1:face_max) = 0 + OFFSET
  face_normal(1:3,1:face_max) = 0.0E+00
  face_num = 0
  face_order(1:face_max) = 0
  face_tex_uv(1:2,1:face_max) = 0.0E+00
  group_num = 0
  line_dex(1:line_max) = -1 + OFFSET
  line_material(1:line_max) = 0 + OFFSET
  line_num = 0
!
!  There is one special default material, which is gray.
!
  material_name(1) = 'Material_Default'
  do i = 2, material_max
    call i4_to_s_zero ( i, char4 )
    material_name(i) = 'Material_' // char4
  end do
  material_num = 1
  material_rgba(1:4,1:material_max) = 0.5E+00

  object_name = 'IVCON'
  object_num = 0
  point(1:point_max) = 0
  point_num = 0

  do i = 1, texture_max
    call i4_to_s_zero ( i, char4 )
    texture_name(i) = 'Texture_' // char4
  end do

  texture_num = 0
  texture_temp(1:2,1:order_max*face_max) = 0.0E+00

  call tmat_init ( transform_matrix )

  vertex_material(1:order_max,1:face_max) = 0 + OFFSET
  vertex_normal(1:3,1:order_max,1:face_max) = 0.0E+00
  vertex_tex_uv(1:2,1:order_max,1:face_max) = 0.0E+00

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_INIT: Initialized graphics data.'
  end if

  return
end
subroutine data_read ( cor3, cor3_material, cor3_max, &
  cor3_normal, cor3_num, cor3_tex_uv, debug, face, face_area, &
  face_material, face_max, face_normal, face_num, face_order, &
  face_tex_uv, filein_name, filein_type, group_num, ierror, line_dex, &
  line_material, line_max, line_num, material_max, material_name, &
  material_num, material_rgba, object_num, order_max, &
  point, point_max, point_num, texture_max, texture_name, texture_num, &
  texture_temp, vertex_material, vertex_normal, vertex_tex_uv )

!*****************************************************************************80
!
!! DATA_READ reads a file into internal graphics data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Output, integer ( kind = 4 ) COR3_MATERIAL(COR3_MAX), the material of each node.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Output, real ( kind = 4 ) COR3_TEX_UV(2,COR3_MAX), UV texture coordinates for nodes.
!
!    Output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Output, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Output, real ( kind = 4 ) FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input/output, character ( len = * ) FILEIN_NAME, the name of the 
!    input file.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Output, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Output, integer ( kind = 4 ) LINE_MATERIAL(LINE_MAX), material index for line.
!
!    Output, character ( len = * ) MATERIAL_NAME(MATERIAL_MAX), material
!    names.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Output, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Output, character ( len = * ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
!    Output, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Output, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
!    Output, real ( kind = 4 ) VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture
!    coordinates.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) point_max
  integer ( kind = 4 ) texture_max

  integer ( kind = 4 ) bad_num
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_material(cor3_max)
  real ( kind = 4 ) cor3_normal(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  real ( kind = 4 ) cor3_tex_uv(2,cor3_max)
  logical debug
  integer ( kind = 4 ) dup_num
  integer ( kind = 4 ) face(order_max,face_max)
  real ( kind = 4 ) face_area(face_max)
  integer ( kind = 4 ) face_material(face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  real ( kind = 4 ) face_tex_uv(2,face_max)
  character ( len = * ) filein_name
  character ( len = 10 ) filein_type
  integer ( kind = 4 ) group_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iline
  integer ( kind = 4 ) ios
  integer ( kind = 4 ), parameter :: iunit = 1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_num
  character ( len = * ) material_name(material_max)
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ) ncol_oogl
  integer ( kind = 4 ) nrow_oogl
  integer ( kind = 4 ) ntemp
  integer ( kind = 4 ) object_num
  integer ( kind = 4 ), parameter :: OFFSET = 1
  integer ( kind = 4 ) order_num
  integer ( kind = 4 ) point(point_max)
  integer ( kind = 4 ) point_num
  logical s_eqi
  integer ( kind = 4 ) text_num
  character ( len = * ) texture_name(texture_max)
  integer ( kind = 4 ) texture_num
  real ( kind = 4 ) texture_temp(2,order_max*face_max)
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  real ( kind = 4 ) vertex_tex_uv(2,order_max,face_max)

  ierror = 0
  bad_num = 0
  dup_num = 0
  text_num = 0
!
!  Open the file.
!
  open ( unit = iunit, file = filein_name, status = 'old', iostat = ios )
 
  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_READ - Error!'
    write ( *, '(a)' ) '  The input file "' // trim ( filein_name ) &
      // '"  could not be opened!'
    ierror = ios
    return
  end if
!
!  Read the information in the file.
!
  if ( s_eqi ( filein_type, 'ASE' ) ) then
 
    call ase_read ( bad_num, cor3, cor3_material, cor3_max, cor3_num, &
      face, face_material, face_max, face_normal, face_num, face_order, &
      filein_name, ierror, iunit, material_max, material_name, &
      material_num, material_rgba, order_max, text_num, &
      vertex_material, vertex_normal )

    call node_to_vertex_material ( cor3_material, cor3_max, face, &
      face_max, face_num, face_order, order_max, vertex_material )

    face_material(1:face_num) = vertex_material(1,1:face_num)

  else if ( s_eqi ( filein_type, 'BYU' ) ) then
 
    call byu_read ( cor3, cor3_max, cor3_num, face, face_max, face_num, &
      face_order, filein_name, ierror, iunit, order_max )

  else if ( s_eqi ( filein_type, 'DXF' ) ) then
 
    call dxf_read ( bad_num, cor3, cor3_max, cor3_num, dup_num, face, &
      face_material, face_max, face_num, face_order, filein_name, &
      ierror, iunit, line_dex, line_material, line_max, line_num, &
      material_num, order_max, text_num )

  else if ( s_eqi ( filein_type, 'HRC' ) ) then
 
    call hrc_read ( bad_num, cor3, cor3_max, cor3_num, dup_num, face, &
      face_material, face_max, face_num, face_order, filein_name, &
      ierror, iunit, line_dex, line_material, line_max, line_num, &
      material_max, material_name, material_num, material_rgba, &
      order_max, text_num, texture_max, texture_name, texture_num, &
      vertex_material, vertex_normal, vertex_tex_uv )

  else if ( s_eqi ( filein_type, 'IV' ) ) then
 
    call iv_read ( bad_num, cor3, cor3_max, cor3_num, &
      debug, face, face_max, face_num, face_order, filein_name, ierror, &
      iunit, line_dex, line_material, line_max, line_num, material_max, &
      material_num, material_rgba, order_max, text_num, &
      texture_max, texture_name, texture_num, texture_temp, &
      vertex_material, vertex_normal, vertex_tex_uv )

  else if ( s_eqi ( filein_type, 'OBJ' ) ) then

    call obj_read ( bad_num, cor3, cor3_max, cor3_num, face, &
      face_material, face_max, face_num, face_order, filein_name, &
      group_num, ierror, iunit, line_dex, line_material, line_max, &
      line_num, material_max, material_name, material_num, material_rgba, &
      object_num, order_max, text_num, vertex_material, vertex_normal )

  else if ( s_eqi ( filein_type, 'OFF' ) ) then

    call off_read ( cor3, cor3_max, cor3_num, face, face_max, &
      face_num, face_order, filein_name, ierror, iunit, order_max, text_num )

  else if ( s_eqi ( filein_type, 'OOGL' ) ) then

    call oogl_read ( cor3, cor3_material, cor3_max, cor3_normal, cor3_num, &
      face, face_area, face_material, face_max, face_normal, face_num, &
      face_order, filein_name, ierror, iunit, material_max, material_name, &
      material_num, material_rgba, ncol_oogl, nrow_oogl, order_max, text_num, &
      vertex_material, vertex_normal )

  else if ( s_eqi ( filein_type, 'SMF' ) ) then

    call smf_read ( bad_num, cor3, cor3_material, cor3_max, cor3_normal, &
      cor3_num, cor3_tex_uv, debug, face, face_material, face_max, &
      face_normal, face_num, face_order, face_tex_uv, filein_name, &
      group_num, ierror, iunit, material_max, material_name, material_num, &
      material_rgba, order_max, text_num, texture_max, texture_num, &
      texture_name, vertex_material )

  else if ( s_eqi ( filein_type, 'STL' ) .or. &
            s_eqi ( filein_type, 'STLA' ) ) then

    call stla_read ( bad_num, cor3, cor3_max, cor3_num, dup_num, face, &
      face_material, face_max, face_normal, face_num, face_order, &
      filein_name, ierror, iunit, material_num, object_num, order_max, &
      text_num, vertex_material, vertex_normal )

  else if ( s_eqi ( filein_type, 'TRI' ) .or. &
            s_eqi ( filein_type, 'TRIA' ) ) then

    call tria_read ( cor3, cor3_max, cor3_num, dup_num, face, &
      face_material, face_max, face_num, face_order, filein_name, ierror, &
      iunit, order_max, text_num, vertex_material, vertex_normal )

  else if ( s_eqi ( filein_type, 'TS' ) .or. &
            s_eqi ( filein_type, '3S' ) ) then

    call ts_read ( bad_num, cor3, cor3_max, cor3_num, face, face_material, &
      face_max, face_num, face_order, filein_name, ierror, iunit, line_dex, &
      line_material, line_max, line_num, material_max, material_num, &
      material_rgba, order_max, point, point_max, point_num )

  else if ( s_eqi ( filein_type, 'VLA' ) ) then
 
    call vla_read ( bad_num, cor3, cor3_max, cor3_num, dup_num, &
      filein_name, ierror, iunit, line_dex, line_material, line_max, &
      line_num, material_num, text_num )

  else if ( s_eqi ( filein_type, 'WRL' ) ) then

    call vrml_read ( bad_num, cor3, cor3_material, cor3_max, cor3_num, &
      face, face_material, face_max, face_num, face_order, filein_name, &
      ierror, iunit, line_dex, line_material, line_max, line_num, &
      material_max, material_name, material_num, material_rgba, &
      order_max, text_num, texture_max, texture_name, texture_num, &
      vertex_material, vertex_normal, vertex_tex_uv )

  else if ( s_eqi ( filein_type, 'XYZ' ) ) then

    call xyz_read ( cor3, cor3_max, cor3_num, filein_name, ierror, iunit, &
      line_dex, line_max, line_num, text_num )

  else 

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized input file type:'
    write ( *, '(a)' ) trim ( filein_type )
    ierror = 1
    return

  end if
 
  close ( unit = iunit )
!
!  Make a report on what we saw in the file.
!
  ntemp = min ( face_num, face_max )
  call i4vec_max ( ntemp, face_order, order_num )

  call data_report ( bad_num, cor3_num, dup_num, face_num, &
    group_num, line_num, material_num, object_num, order_num, text_num )
!
!  Warn about any errors that occurred during reading.
!
  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_READ - Warning!'
    write ( *, '(a)' ) '  An error occurred reading the input file.'
    return
  end if
!
!  Check the data.
!
!  You MUST wait until after this check before doing other computations,
!  since COR3_NUM, and the other variables could be much larger than
!  the legal maximums, until corrected by this routine.
!
  call data_check ( cor3_max, cor3_num, face_max, face_num, face_order, &
    line_max, line_num, material_max, material_name, material_num, order_max, &
    texture_max, texture_num, texture_name )
!
!  MATERIALS FIXUPS:
!
!  If there are no materials, define one.
!
  if ( material_num <= 0 ) then

    material_num = 1

    material_name(material_num) = 'Default_Material'

    material_rgba(1,material_num) = 0.7E+00
    material_rgba(2,material_num) = 0.7E+00
    material_rgba(3,material_num) = 0.7E+00
    material_rgba(4,material_num) = 1.0E+00

  end if
!
!  If a point hasn't been assigned a material, set it to material 1.
!
  do i = 1, cor3_num
    if ( cor3_material(i) < 1 .or. material_num < cor3_material(i) ) then
      cor3_material(i) = 1
    end if
  end do
!
!  If a vertex hasn't been assigned a material, set it to material 1.
!
  do i = 1, face_num
    do j = 1, face_order(i)
      if ( vertex_material(j,i) < 1 .or. &
           material_num < vertex_material(j,i) ) then
        vertex_material(j,i) = 1
      end if
    end do
  end do
!
!  If a face hasn't been assigned a material, set it to material 1.
!
  do i = 1, face_num
    if ( face_material(i) < 1 .or. material_num < face_material(i) ) then
      face_material(i) = 1
    end if
  end do
!
!  If a line item hasn't been assigned a material, set it to material 1.
!
  do iline = 1, line_num
    if ( line_dex(iline) == -1 + OFFSET ) then
      line_material(iline) = -1 + OFFSET
    else if ( line_material(iline) < 1 .or. &
      material_num < line_material(iline) ) then
      line_material(iline) = material_num
    end if
  end do
!
!  NULL EDGE DELETION.
!
  call edge_null_delete ( cor3, cor3_max, face, face_max, face_num, &
    face_order, order_max, vertex_normal )
!
!  COMPUTE FACE AREAS.
!
  call face_area_set ( cor3, cor3_max, face, face_area, face_max, &
    face_num, face_order, order_max )
!
!  NULL FACE DELETION.
!
  call face_null_delete ( face, face_area, face_material, face_max, &
    face_num, face_order, order_max, vertex_material, vertex_normal )
!
!  NORMAL VECTOR FIXUPS:
!
!  Recompute zero vertex normals from vertex positions.
!
  call vertex_normal_set ( cor3, cor3_max, face, face_max, &
    face_num, face_order, order_max, vertex_normal )
!
!  Recompute zero node normals by averaging vertex normals.
!
  call cor3_normal_set ( cor3_max, cor3_normal, face, face_area, &
    face_max, face_num, face_order, order_max, vertex_normal )
!
!  Recompute zero face normals by averaging vertex normals.
!
  call face_normal_ave ( face_max, face_normal, face_num, face_order, &
    order_max, vertex_normal )
!
!  Report the range of the nodal coordinates.
!
  call cor3_range ( cor3, cor3_max, cor3_num )

  return
end
subroutine data_report ( bad_num, cor3_num, dup_num, &
  face_num, group_num, line_num, material_num, object_num, order_num, text_num )

!*****************************************************************************80
!
!! DATA_REPORT gives a summary of the contents of the data file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Input, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
  implicit none

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) dup_num
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) group_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) material_num
  integer ( kind = 4 ) object_num
  integer ( kind = 4 ) order_num
  integer ( kind = 4 ) text_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DATA_REPORT - The input file contains:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Bad data items             ', bad_num
  write ( *, '(a,i6)' ) '  Text lines                 ', text_num
  write ( *, '(a,i6)' ) '  Duplicate points           ', dup_num
  write ( *, '(a,i6)' ) '  Faces                      ', face_num
  write ( *, '(a,i6)' ) '  Groups                     ', group_num
  write ( *, '(a,i6)' ) '  Vertices per face, maximum ', order_num
  write ( *, '(a,i6)' ) '  Line items                 ', line_num
  write ( *, '(a,i6)' ) '  Materials                  ', material_num
  write ( *, '(a,i6)' ) '  Points                     ', cor3_num
  write ( *, '(a,i6)' ) '  Objects                    ', object_num

  return
end
subroutine data_write ( cor3, cor3_material, cor3_max, cor3_normal, cor3_num, &
  cor3_tex_uv, debug, face, face_material, face_max, face_normal, face_num, &
  face_order, face_tex_uv, filein_name, fileout_name, fileout_type, ierror, &
  line_dex, line_material, line_max, line_num, line_prune, material_name, &
  material_max, material_num, material_rgba, object_name, order_max, &  
  point, point_max, point_num, texture_max, texture_name, texture_num, &
  vertex_material, vertex_normal, vertex_tex_uv )

!*****************************************************************************80
!
!! DATA_WRITE writes the internal graphics data to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MATERIAL(COR3_MAX), the material of each node.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of 3D points.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, real ( kind = 4 ) COR3_TEX_UV(2,COR3_MAX), UV texture coordinates for nodes.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, real ( kind = 4 ) FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Input/output, character ( len = * ) FILEOUT_NAME, the name of the 
!    output file.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Input, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer ( kind = 4 ) LINE_MATERIAL(LINE_MAX), material index for line.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, character ( len = * ) MATERIAL_NAME(MATERIAL_MAX), material names.
!
!    Input, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer ( kind = 4 ) TEXTURE_MAX, the maximum number of textures.
!
!    Input, character ( len = * ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
!    Input, real ( kind = 4 ) TRANSFORM_MATRIX(4,4), the transformation matrix.
!
!    Input, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
!    Input, real ( kind = 4 ) VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture
!    coordinates.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) point_max
  integer ( kind = 4 ) texture_max

  real ( kind = 4 ) cor2(2,cor3_max)
  integer ( kind = 4 ) cor2_max
  integer ( kind = 4 ) cor2_num
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_material(cor3_max)
  real ( kind = 4 ) cor3_normal(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  real ( kind = 4 ) cor3_tex_uv(2,cor3_max)
  logical debug
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  real ( kind = 4 ) face_tex_uv(2,face_max)
  character ( len = * ) filein_name
  character ( len = * ) fileout_name
  character ( len = 10 ) fileout_type
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ), parameter :: iunit = 1
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) line_num_save
  integer ( kind = 4 ) line_prune
  character ( len = * ) material_name(material_max)
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  character ( len = * ) object_name
  integer ( kind = 4 ) point(point_max)
  integer ( kind = 4 ) point_num
  logical s_eqi
  character ( len = * ) texture_name(texture_max)
  integer ( kind = 4 ) texture_num
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  real ( kind = 4 ) vertex_tex_uv(2,order_max,face_max)

  ierror = 0
!
!  Open the file.
!
  open ( unit = iunit, file = fileout_name, status = 'replace', iostat = ios )
  
  if ( ios /= 0 ) then
    ierror = ios
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_WRITE - Warning!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if
!
!  Write an Autodesk file...
!
  if ( s_eqi ( fileout_type, 'ASE' ) ) then
 
    call ase_write ( cor3, cor3_max, cor3_num, face, face_max, &
      face_normal, face_num, face_order, filein_name, fileout_name, &
      iunit, order_max, vertex_normal )
!
!  ...or a BYU file...
!
  else if ( s_eqi ( fileout_type, 'BYU' ) ) then

    call byu_write ( cor3, cor3_max, cor3_num, face, face_max, face_num, &
      face_order, fileout_name, iunit, order_max )
!
!  ...or a DXF file...
!
  else if ( s_eqi ( fileout_type, 'DXF' ) ) then
 
    call dxf_write ( cor3, cor3_max, face, face_max, face_num, &
      face_order, filein_name, fileout_name, iunit, line_dex, &
      line_max, line_num, order_max )
!
!  ...or an HRC SoftImage file...
!
  else if ( s_eqi ( fileout_type, 'HRC' ) ) then
 
    call hrc_write ( cor3, cor3_max, cor3_num, face, face_material, &
      face_max, face_num, face_order, fileout_name, iunit, line_dex, &
      line_max, line_num, material_max, material_name, material_num, &
      material_rgba, order_max, texture_max, texture_name, texture_num, &
      vertex_normal, vertex_tex_uv )
!
!  ...or an IV Inventor file...
!
  else if ( s_eqi ( fileout_type, 'IV' ) ) then
 
!   if ( face_num == 0 .and. line_num == 0 ) then

    if ( face_num == 0 ) then

      call iv_point_write ( cor3, cor3_max, cor3_num, filein_name, &
        fileout_name, line_dex, line_max, line_num, iunit )

    else

      call iv_write ( cor3, cor3_max, cor3_normal, cor3_num, face, &
        face_max, face_num, face_order, filein_name, fileout_name, &
        iunit, line_dex, line_material, line_max, line_num, material_max, &
        material_num, material_rgba, order_max, texture_max, &
        texture_num, texture_name, vertex_material, vertex_tex_uv )

    end if
!
!  ...or a WaveFront OBJ file...
!
  else if ( s_eqi ( fileout_type, 'OBJ' ) ) then
 
    call obj_write ( cor3, cor3_max, cor3_num, face, face_max, face_num, &
      face_order, filein_name, fileout_name, iunit, line_dex, line_max, &
      line_num, order_max, vertex_normal )
!
!  ...or a GEOMVIEW OFF file...
!
  else if ( s_eqi ( fileout_type, 'OFF' ) ) then

    call off_write ( cor3, cor3_max, cor3_num, face, face_max, face_num, &
      face_order, fileout_name, iunit, order_max )
!
!  ...or a POV file...
!
  else if ( s_eqi ( fileout_type, 'POV' ) ) then
 
    call pov_write ( cor3, cor3_max, face, face_material, face_max, &
      face_num, face_order, filein_name, fileout_name, iunit, &
      material_max, material_num, material_rgba, order_max, vertex_normal )
!
!  ...or a PS file...
!
  else if ( s_eqi ( fileout_type, 'PS' ) ) then
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WATCH OUT!'
    write ( *, '(a)' ) 'PS_WRITE not ready yet!'

    cor2_num = cor3_num
    cor2_max = cor3_max

    call project_2d ( cor2, cor3, ierror, cor2_max, cor3_max, cor2_num, &
      cor3_num )

    if ( ierror == 0 ) then

      call ps_write ( cor2, cor2_max, cor2_num, face, face_material, &
        face_max, face_num, face_order, fileout_name, iunit, line_dex, &
        line_material, line_max, line_num, material_max, &
        material_rgba, order_max )

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DATA_WRITE - Error!'
      write ( *, '(a)' ) '  2D projection canceled.'
      return

    end if
!
!  ...or an SMF file...
!
  else if ( s_eqi ( fileout_type, 'SMF' ) ) then
 
    call smf_write ( cor3, cor3_material, cor3_max, cor3_normal, cor3_num, &
      cor3_tex_uv, face, face_max, face_num, face_order, filein_name, &
      fileout_name, iunit, material_max, material_rgba, order_max, &
      texture_max, texture_name, texture_num )
!
!  ...or an ASCII STL file...
!
  else if ( s_eqi ( fileout_type, 'STL' ) .or. &
            s_eqi ( fileout_type, 'STLA' ) ) then
 
    call stla_write ( cor3, cor3_max, face, face_max, face_normal, &
      face_num, face_order, filein_name, fileout_name, iunit, order_max )
!
!  ...or a TEC file...
!
  else if ( s_eqi ( fileout_type, 'TEC' ) ) then
 
    call tec_write ( cor3, cor3_material, cor3_max, cor3_num, face, &
      face_max, face_num, face_order, fileout_name, iunit, material_max, &
      material_rgba, order_max )
!
!  ...or a TRI/TRIA file...
!
  else if ( s_eqi ( fileout_type, 'TRI' ) .or. &
            s_eqi ( fileout_type, 'TRIA' ) ) then

    call tria_write ( cor3, cor3_max, cor3_normal, face, face_max, face_num, &
      face_order, fileout_name, iunit, order_max )
!
!  ...or a Mathematics TS file...
!
  else if ( s_eqi ( fileout_type, 'TS' ) .or. &
            s_eqi ( fileout_type, '3S' ) ) then
 
    call ts_write ( cor3, cor3_max, cor3_num, face, face_material, &
      face_max, face_num, face_order, filein_name, fileout_name, iunit, &
      line_dex, line_material, line_max, line_num, material_max, &
      material_num, material_rgba, order_max, point, point_max, point_num )
!
!  ...or a TXT text file...
!
  else if ( s_eqi ( fileout_type, 'TXT' ) ) then

    call txt_write ( cor3, cor3_material, cor3_max, cor3_normal, cor3_num, &
      cor3_tex_uv, face, face_material, face_max, face_normal, face_num, &
      face_order, face_tex_uv, filein_name, fileout_name, iunit, line_dex, &
      line_material, line_max, line_num, material_max, material_name, &
      material_num, material_rgba, object_name, order_max, point, point_max, &
      point_num, texture_max, texture_name, texture_num, vertex_material, &
      vertex_normal, vertex_tex_uv )
!
!  ...or a UCD file...
!
  else if ( s_eqi ( fileout_type, 'UCD' ) ) then
 
    call ucd_write ( cor3, cor3_material, cor3_max, cor3_num, face, &
      face_material, face_max, face_num, face_order, fileout_name, iunit, &
      material_max, material_num, material_rgba, order_max )
!
!  ...or a VLA file...
!
  else if ( s_eqi ( fileout_type, 'VLA' ) ) then
 
    line_num_save = line_num

    if ( 0 < face_num ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DATA_WRITE - Note:'
      write ( *, '(a)' ) '  Face information will be temporarily'
      write ( *, '(a)' ) '  converted to line information for '
      write ( *, '(a)' ) '  the VLA output.'

      call face_to_line ( debug, face, face_max, face_num, face_order, &
        line_dex, line_material, line_max, line_num, line_prune, &
        order_max, vertex_material )

      if ( line_max < line_num ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DATA_WRITE - Note:'
        write ( *, '(a)' ) '  Some face information was lost.'
        write ( *, '(a,i6)' ) '  The maximum number of lines is ', line_max
        write ( *, '(a,i6)' ) '  but we would need at least ', line_num
        line_num = line_max

      end if

    end if

    call vla_write ( cor3, cor3_max, filein_name, fileout_name, iunit, &
      line_dex, line_max, line_num )

    line_num = line_num_save
!
!  ...or a VRML file...
!
  else if ( s_eqi ( fileout_type, 'WRL' ) ) then
 
    call vrml_write ( cor3, cor3_max, cor3_num, face, face_max, face_num, &
      face_order, filein_name, fileout_name, iunit, line_dex, line_material, &
      line_max, line_num, material_max, material_num, material_rgba, &
      order_max, vertex_material )
!
!  ...or an XGL file...
!
  else if ( s_eqi ( fileout_type, 'XGL' ) ) then
 
    call xgl_write ( cor3, cor3_max, cor3_normal, cor3_num, face, &
      face_material, face_max, face_num, face_order, fileout_name, iunit, &
      material_max, material_num, material_rgba, order_max )
!
!  ...or an XYZ file...
!
  else if ( s_eqi ( fileout_type, 'XYZ' ) ) then
 
    line_num_save = line_num

    if ( 0 < face_num ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DATA_WRITE - Note:'
      write ( *, '(a)' ) '  Face information will be temporarily'
      write ( *, '(a)' ) '  converted to line information for '
      write ( *, '(a)' ) '  the XYZ output.'

      call face_to_line ( debug, face, face_max, face_num, face_order, &
        line_dex, line_material, line_max, line_num, line_prune, &
        order_max, vertex_material )

      if ( line_max < line_num ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DATA_WRITE - Note:'
        write ( *, '(a)' ) '  Some face information was lost.'
        write ( *, '(a,i6)' ) '  The maximum number of lines is ', line_max
        write ( *, '(a,i6)' ) '  but we would need at least ', line_num
        line_num = line_max

      end if

    end if

    call xyz_write ( cor3, cor3_max, cor3_num, filein_name, fileout_name, &
      iunit, line_dex, line_max, line_num )

    line_num = line_num_save
!
!  ...or what the hell happened?
!
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_WRITE - Warning!'
    write ( *, '(a)' ) '  Unrecognized output file type.'
    ierror = 1
 
  end if
!
!  Close the file.
!
  close ( unit = iunit )
 
  return
end
function degrees_to_radians ( angle )

!*****************************************************************************80
!
!! DEGREES_TO_RADIANS converts an angle from degrees to radians.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) ANGLE, an angle in degrees.
!
!    Output, real ( kind = 4 ) DEGREES_TO_RADIANS, the equivalent angle
!    in radians.
!
  implicit none

  real ( kind = 4 ) angle
  real ( kind = 4 ) degrees_to_radians
  real ( kind = 4 ), parameter :: pi = &
    3.14159265358979323846264338327950288419716939937510E+00

  degrees_to_radians = ( angle / 180.0E+00 ) * pi

  return
end
subroutine digit_to_ch ( digit, c )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
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
!    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
function dot0_3d ( x0, y0, z0, x1, y1, z1, x2, y2, z2 )

!*****************************************************************************80
!
!! DOT0_3D computes the dot product of (P1-P0) and (P2-P0) in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X0, Y0, Z0, the coordinates of the point P0.
!
!    Input, real ( kind = 4 ) X1, Y1, Z1, the coordinates of the point P1.
!
!    Input, real ( kind = 4 ) X2, Y2, Z2, the coordinates of the point P2.
!
!    Output, real ( kind = 4 ) DOT0_3D, the dot product of (P1-P0) and (P2-P0).
!
  implicit none

  real ( kind = 4 ) dot0_3d
  real ( kind = 4 ) x0
  real ( kind = 4 ) x1
  real ( kind = 4 ) x2
  real ( kind = 4 ) y0
  real ( kind = 4 ) y1
  real ( kind = 4 ) y2
  real ( kind = 4 ) z0
  real ( kind = 4 ) z1
  real ( kind = 4 ) z2

  dot0_3d = ( x1 - x0 ) * ( x2 - x0 ) + ( y1 - y0 ) * ( y2 - y0 ) + &
    ( z1 - z0 ) * ( z2 - z0 )
 
  return
end
subroutine dxf_read ( bad_num, cor3, cor3_max, cor3_num, dup_num, face, &
  face_material, face_max, face_num, face_order, filein_name, &
  ierror, iunit, line_dex, line_material, line_max, line_num, &
  material_num, order_max, text_num )

!*****************************************************************************80
!
!! DXF_READ reads graphics information from an AutoCAD DXF file.
!
!  Discussion:
!
!    It is intended that the information read from the file can
!    either start a whole new graphics object, or simply be added
!    to a current graphics object via the '<<' command.  
!
!    This is controlled by whether the input values have been zeroed 
!    out or not.  This routine simply tacks on the information it 
!    finds to the current graphics object.
!
!  Example:
!
!      0
!    SECTION
!      2
!    HEADER
!    999
!    diamond.dxf created by IVREAD.
!    999
!    Original data in diamond.obj.
!      0
!    ENDSEC
!      0
!    SECTION
!      2
!    TABLES
!      0
!    ENDSEC
!      0
!    SECTION
!      2
!    BLOCKS
!      0
!    ENDSEC
!      0
!    SECTION
!      2
!    ENTITIES
!      0
!    LINE
!      8
!    0
!     10
!      0.00  (X coordinate of beginning of line.)
!     20
!      0.00  (Y coordinate of beginning of line.)
!     30
!      0.00  (Z coordinate of beginning of line.)
!     11
!      1.32  (X coordinate of end of line.)
!     21
!      1.73  (Y coordinate of end of line.)
!     31
!      2.25  (Z coordinate of end of line.)
!      0
!    3DFACE
!      8
!     Cube
!    10
!    -0.50  (X coordinate of vertex 1)
!    20
!     0.50  (Y coordinate of vertex 1)   
!    30
!      1.0  (Z coordinate of vertex 1)  
!    11
!     0.50  (X coordinate of vertex 2)  
!    21
!     0.50  (Y coordinate of vertex 2)
!    31
!      1.0  (Z coordinate of vertex 2)
!    12
!     0.50  (X coordinate of vertex 3) 
!    22
!     0.50  (Y coordinate of vertex 3)
!    32
!     0.00  (Z coordinate of vertex 3)
!     0
!    ENDSEC
!      0
!    EOF
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BAD_NUM, the number of "bad" lines of input text.
!
!    Input/output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input/output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Output, integer ( kind = 4 ) DUP_NUM, the number of duplicate nodes discovered.
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input/output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices 
!    per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which data is read.
!
!    Input/output, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, 
!    terminated by -1.
!
!    Input/output, integer ( kind = 4 ) LINE_MATERIAL(LINE_MAX), material index for line.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line items.
!
!    Input/output, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Output, integer ( kind = 4 ) TEXT_NUM, the number of lines of input text.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) order_max

  integer ( kind = 4 ) bad_num
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  real ( kind = 4 ) cvec(3)
  integer ( kind = 4 ) dup_num
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icode
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) ixyz
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) material_num
  character ( len = 256 ) mode
  integer ( kind = 4 ), parameter :: OFFSET = 1
  real ( kind = 4 ) rval
  logical s_eqi
  integer ( kind = 4 ) text_num
  character ( len = 256 ) text1
  character ( len = 256 ) text2

  ierror = 0
  mode = ' '
!
!  Read the next pair of lines.  
!  TEXT1 is a numeric tag, TEXT2 contains data.
!
  do
 
    read ( iunit, '(a)', iostat = ios ) text1

    if ( ios /= 0 ) then
      exit
    end if

    text_num = text_num + 1
 
    call s_to_i4 ( text1, icode, ierror, lchar )
 
    if ( ierror /= 0 ) then
      ierror = 1
      bad_num = bad_num + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DXF_READ - Fatal error!'
      write ( *, '(a)' ) '  Could not interpret line:'
      write ( *, '(a)' ) trim ( text1 )
      return
    end if
!
!  Read the second item, which might be a label or numeric value.
!
    read ( iunit, '(a)', iostat = ios ) text2

    if ( ios /= 0 ) then

      if ( text1 /= ' ' ) then
        ierror = 3
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DXF_READ - Warning!'
        write ( *, '(a)' ) '  The last code was not followed by data:'
        write ( *, '(a)' ) trim ( text1 )
      end if

      exit

    end if

    text_num = text_num + 1
!
!  Codes 0 through 9 are followed by a string.
!
!  All we want to do is know when we can expect to be reading
!  LINE data and when we can expect to be reading FACE data.
!
    if ( 0 <= icode .and. icode <= 9 ) then

      if ( s_eqi ( text2(1:6), '3DFACE' ) ) then
        mode = '3DFACE'
        ivert = 0
      else if ( s_eqi ( text2(1:4), 'LINE' ) ) then
        mode = 'LINE'
      end if
!
!  Codes 10 through 59 are followed by a real value.
!
!    10, 11, 12, ... are followed by a line of X data;
!    20, 21, 22, ... are followed by a line of Y data;
!    30, 31, 33, ... are followed by a line of Z data.
!
    else if ( 10 <= icode .and. icode <= 59 ) then
 
      call s_to_r4 ( text2, rval, ierror, lchar )
 
      if ( ierror /= 0 ) then

        ierror = 2
        rval = 0.0E+00
        bad_num = bad_num + 1

        if ( bad_num == 1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DXF_READ - Fatal error!'
          write ( *, '(a)' ) '  Could not interpret line:'
          write ( *, '(a)' ) trim ( text2 )
        end if

      end if
 
      if ( mode == 'LINE' ) then

        if ( icode == 10 .and. 0 < line_num ) then

          line_num = line_num + 1
          if ( line_num <= line_max ) then
            line_dex(line_num) = -1 + OFFSET
            line_material(line_num) = -1 + OFFSET
          end if

        end if

      else if ( mode == '3DFACE' ) then

        if ( icode == 10 ) then

          face_num = face_num + 1
          face_order(face_num) = 0
          face_material(face_num) = material_num

        end if

      end if

      if ( 10 <= icode .and. icode <= 19 ) then
        ixyz = 1
      else if ( icode >= 20 .and. icode <= 29 ) then
        ixyz = 2
      else if ( icode >= 30 .and. icode <= 39 ) then
        ixyz = 3
      end if

      cvec(ixyz) = rval
!
!  Once the entire (X,Y,Z) triple has been read, check to see if the
!  values in CVEC already exist in COR3.  If so, save space by using 
!  the index of a previous copy.
!
!  Otherwise, add CVEC to COR3, and increment COR3_NUM.
!
      if ( ixyz == 3 ) then

        if ( cor3_num <= 1000 ) then
          call r4col_find ( 3, 3, cor3_num, cor3, cvec, icor3 )
        else
          icor3 = 0
        end if
 
        if ( icor3 == 0 ) then

          cor3_num = cor3_num + 1
          icor3 = cor3_num

          if ( cor3_num <= cor3_max ) then
            cor3(1:3,cor3_num) = cvec(1:3)
          end if

        else

          dup_num = dup_num + 1

        end if

        if ( mode == 'LINE' ) then

          line_num = line_num + 1

          if ( line_num <= line_max ) then
            line_dex(line_num) = icor3 - 1 + OFFSET
            line_material(line_num) = material_num
          end if

        else if ( mode == '3DFACE' ) then

          ivert = ivert + 1
          face(ivert,face_num) = icor3
          face_order(face_num) = face_order(face_num) + 1

        end if

      end if
!
!  Codes 60 through 79 are followed by an integer.
!
    else if ( 60 <= icode .and. icode <= 79 ) then

      call s_to_i4 ( text2, ival, ierror, lchar )
!
!  Codes 140 through 147 are followed by a real.
!
    else if ( icode >= 140 .and. icode <= 147 ) then

      call s_to_r4 ( text2, rval, ierror, lchar )
!
!  Codes 170 through 175 are followed by an integer.
!
    else if ( icode >= 170 .and. icode <= 175 ) then

      call s_to_i4 ( text2, ival, ierror, lchar )
!
!  Codes 210 through 239 are followed by a real.
!
    else if ( icode >= 210 .and. icode <= 239 ) then

      call s_to_r4 ( text2, rval, ierror, lchar )
!
!  Code 999 is followed by a (comment) string.
!
    else if ( icode == 999 ) then
!
!  Codes 1000 through 1009 are followed by a string.
!
    else if ( icode >= 1000 .and. icode <= 1009 ) then
!
!  Codes 1010 through 1059 are followed by a real.
!
    else if ( icode >= 1010 .and. icode <= 1059 ) then

      call s_to_r4 ( text2, rval, ierror, lchar )
!
!  Codes 1060 through 1079 are followed by an integer.
!
    else if ( icode >= 1060 .and. icode <= 1079 ) then

      call s_to_i4 ( text2, ival, ierror, lchar )
!
!  Unrecognized code.
!
    else

      bad_num = bad_num + 1

    end if
 
  end do
!
!  END OF INPUT.
!
!  Slap a trailing "-1" on the end of the line indices.
!
  if ( 0 < line_num ) then
    line_num = line_num + 1
    if ( line_num <= line_max ) then
      line_dex(line_num) = -1 + OFFSET
      line_material(line_num) = -1 + OFFSET
    end if
  end if
 
  return
end
subroutine dxf_write ( cor3, cor3_max, face, face_max, face_num, &
  face_order, filein_name, fileout_name, iunit, line_dex, line_max, &
  line_num, order_max )

!*****************************************************************************80
!
!! DXF_WRITE writes graphics data to an AutoCAD DXF file.
!
!  Example:
!
!      0
!    SECTION
!      2
!    HEADER
!    999
!    diamond.dxf created by IVREAD.
!    999
!    Original data in diamond.obj.
!      0
!    ENDSEC
!      0
!    SECTION
!      2
!    TABLES
!      0
!    ENDSEC
!      0
!    SECTION
!      2
!    BLOCKS
!      0
!    ENDSEC
!      0
!    SECTION
!      2
!    ENTITIES
!      0
!    LINE
!      8
!    0
!     10
!      0.00  (X coordinate of beginning of line.)
!     20
!      0.00  (Y coordinate of beginning of line.)
!     30
!      0.00  (Z coordinate of beginning of line.)
!     11
!      1.32  (X coordinate of end of line.)
!     21
!      1.73  (Y coordinate of end of line.)
!     31
!      2.25  (Z coordinate of end of line.)
!      0
!    3DFACE
!      8
!     Cube
!    10
!    -0.50  (X coordinate of vertex 1)
!    20
!     0.50  (Y coordinate of vertex 1)   
!    30
!      1.0  (Z coordinate of vertex 1)  
!    11
!     0.50  (X coordinate of vertex 2)  
!    21
!     0.50  (Y coordinate of vertex 2)
!    31
!      1.0  (Z coordinate of vertex 2)
!    12
!     0.50  (X coordinate of vertex 3) 
!    22
!     0.50  (Y coordinate of vertex 3)
!    32
!     0.00  (Z coordinate of vertex 3)
!     0
!    ENDSEC
!      0
!    EOF
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  character ( len = * ) fileout_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) jcor3
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_num
  logical newline
  integer ( kind = 4 ), parameter :: OFFSET = 1
  integer ( kind = 4 ) text_num

  text_num = 0

  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'SECTION'
  write ( iunit, '(i3)' ) 2
  write ( iunit, '(a)' ) 'HEADER'
  write ( iunit, '(i3)' ) 999
  write ( iunit, '(a,a)' ) trim ( fileout_name ), ' created by IVREAD.'
  write ( iunit, '(i3)' ) 999
  write ( iunit, '(a,a)' ) 'Original data in ', trim ( filein_name ) // '.'
  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'ENDSEC'
  text_num = text_num + 10

  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'SECTION'
  write ( iunit, '(i3)' ) 2
  write ( iunit, '(a)' ) 'TABLES'
  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'ENDSEC'
  text_num = text_num + 6

  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'SECTION'
  write ( iunit, '(i3)' ) 2
  write ( iunit, '(a)' ) 'BLOCKS'
  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'ENDSEC'
  text_num = text_num + 6

  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'SECTION'
  write ( iunit, '(i3)' ) 2
  write ( iunit, '(a)' ) 'ENTITIES'
  text_num = text_num + 4

  jcor3 = 0
  newline = .true.

  do i = 1, line_num

    icor3 = line_dex(i)

    if ( line_dex(i) - OFFSET == -1 ) then

      newline = .true.
!
!  LINE_DEX(I) is the index of a new point that begins or continues a line.
!
    else
!
!  LINE_DEX(I) is the index of a new point that continues a line. 
!    Output the pair of points that define this segment of the line.
!
      if ( .not. newline ) then

        write ( iunit, '(i3)' ) 0
        write ( iunit, '(a)' ) 'LINE'
        write ( iunit, '(i3)' ) 8
        write ( iunit, '(i3)' ) 0

        write ( iunit, '(i3)' ) 10
        write ( iunit, '(g12.4)' ) cor3(1,jcor3)
        write ( iunit, '(i3)' ) 20
        write ( iunit, '(g12.4)' ) cor3(2,jcor3)
        write ( iunit, '(i3)' ) 30
        write ( iunit, '(g12.4)' ) cor3(3,jcor3)

        write ( iunit, '(i3)' ) 11
        write ( iunit, '(g12.4)' ) cor3(1,icor3)
        write ( iunit, '(i3)' ) 21
        write ( iunit, '(g12.4)' ) cor3(2,icor3)
        write ( iunit, '(i3)' ) 31
        write ( iunit, '(g12.4)' ) cor3(3,icor3)

        text_num = text_num + 16
 
      end if
!
!  Save the index of this new point, and note that a line is in progress.
!
      jcor3 = icor3
      newline = .false.

    end if
 
  end do
!
!  Handle faces.
!  This is going to fail bigtime if FACE_ORDER is larger than 9
!
  do iface = 1, face_num
    
    write ( iunit, '(a)' ) '  0'
    write ( iunit, '(a)' ) '3DFACE'
    write ( iunit, '(a)' ) '  8'
    write ( iunit, '(a)' ) '  Cube'
    text_num = text_num + 4
    
    do ivert = 1, face_order(iface)

      icor3 = face(ivert,iface)
  
      write ( iunit, '(i3)' ) 10 + ivert - 1
      write ( iunit, '(g12.4)' ) cor3(1,icor3)
      write ( iunit, '(i3)' ) 20 + ivert - 1
      write ( iunit, '(g12.4)' ) cor3(2,icor3)
      write ( iunit, '(i3)' ) 30 + ivert - 1
      write ( iunit, '(g12.4)' ) cor3(3,icor3)

      text_num = text_num + 6

    end do
  end do

  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'ENDSEC'
  write ( iunit, '(i3)' ) 0
  write ( iunit, '(a)' ) 'EOF'
 
  text_num = text_num + 4
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DXF_WRITE - Wrote ', text_num, ' text lines to ' &
    // trim ( fileout_name )
 
  return
end
subroutine edge_add_nodes ( edge, edge_max, edge_num, iface, n1, n2, ierror )

!*****************************************************************************80
!
!! EDGE_ADD_NODES adds the edge defined by two nodes to the edge list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EDGE(4,EDGE_MAX), contains edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer ( kind = 4 ) EDGE_MAX, the maximum number of edges.
!
!    Input/output, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
!    Input, integer ( kind = 4 ) IFACE, the face to which the nodes belong.
!
!    Input, integer ( kind = 4 ) N1, N2, two nodes which form an edge.
!
!    Output, integer ( kind = 4 ) IERROR, error flag, 0 = no error, nonzero = error.
!
  implicit none

  integer ( kind = 4 ) edge_max

  integer ( kind = 4 ) edge(4,edge_max)
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  if ( edge_num < edge_max ) then
    edge_num = edge_num + 1
    edge(1,edge_num) = n1
    edge(2,edge_num) = n2
    edge(3,edge_num) = iface
    edge(4,edge_num) = 0
    ierror = 0
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EDGE_ADD_NODES - Fatal error!'
    write ( *, '(a,i6)' ) '  Exceeding EDGE_MAX = ', edge_max
    ierror = 1
  end if

  return
end
subroutine edge_bound ( edge, edge_max, edge_num )

!*****************************************************************************80
!
!! EDGE_BOUND reports the edges which are part of the boundary.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EDGE(4,EDGE_MAX), contains edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer ( kind = 4 ) EDGE_MAX, the maximum number of edges.
!
!    Input, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
  implicit none

  integer ( kind = 4 ) edge_max

  integer ( kind = 4 ) bound_num
  logical, parameter :: DEBUG = .false.
  integer ( kind = 4 ) edge(4,edge_max)
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) iedge

  if ( DEBUG ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Boundary edges:'
    write ( *, '(a)' ) ' '
  end if

  bound_num = 0

  do iedge = 1, edge_num
    if ( edge(4,iedge) == 0 ) then
      bound_num = bound_num + 1
      if ( DEBUG ) then
        write ( *, '(2i6)' ) edge(2,iedge), edge(1,iedge)
      end if
    end if
  end do

  write ( *, '(a,i6,a)' ) 'EDGE_BOUND found ', bound_num, ' boundary edges.'

  return
end
subroutine edge_count ( face_max, order_max, face_num, face, face_order, &
  edge_num )

!*****************************************************************************80
!
!! EDGE_COUNT determines the number of edges in a graph.
!
!  Discussion:
!
!    The routine extracts the successive pairs of vertices that
!    define each edge of a face.  It reorders each pair so that 
!    the lesser element is listed first.  It sorts the entire list.
!    Then it counts the unique entries.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Output, integer ( kind = 4 ) EDGE_NUM, the number of unique edges.
!
  implicit none

  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: edge
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) edge_num_old
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) vert
  integer ( kind = 4 ) vert2
!
!  First count the number of edges with duplication.
!
  edge_num = 0
  do iface = 1, face_num
    edge_num = edge_num + face_order(iface)
  end do
!
!  Allocate space, and store the edges.
!
  allocate ( edge(edge_num,2) )

  edge_num = 0
  do iface = 1, face_num
    do vert = 1, face_order(iface)
      edge_num = edge_num + 1
      i = face(vert,iface)
      vert2 = i4_wrap ( vert+1, 1, face_order(iface) )
      j = face(vert2,iface)
      edge(edge_num,1) = min ( i, j )
      edge(edge_num,2) = max ( i, j )
    end do
  end do
!
!  Sort the edges.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( edge_num, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call r4_swap ( edge(i,1), edge(j,1) )
      call r4_swap ( edge(i,2), edge(j,2) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      if ( edge(i,1) < edge(j,1) ) then
        isgn = -1
      else if ( edge(i,1) == edge(j,1) ) then
        if ( edge(i,2) < edge(j,2) ) then
          isgn = -1
        else if ( edge(i,2) == edge(j,2) ) then
          isgn = 0
        else
          isgn = 1
        end if
      else
        isgn = 1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do
!
!  Count the unique entries.
!
  edge_num_old = edge_num

  edge_num = 0

  do i = 1, edge_num_old
 
    if ( i == 1 ) then
      edge_num = 1
    else
      if ( edge(i-1,1) /= edge(i,1) .or. &
           edge(i-1,2) /= edge(i,2) ) then
        edge_num = edge_num + 1
      end if
    end if
 
  end do

  deallocate ( edge )

  return
end
subroutine edge_match_face ( edge, edge_max, edge_num, face_list, n, index )

!*****************************************************************************80
!
!! EDGE_MATCH_FACE seeks an edge common to a face and the edge list.
!
!  Discussion:
!
!    If a common edge is found, then the information in the face node
!    list is adjusted so that the first two entries correspond to the
!    matching edge in EDGE, but in reverse order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EDGE(4,EDGE_MAX), contains edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer ( kind = 4 ) EDGE_MAX, the maximum number of edges.
!
!    Input, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
!    Input/output, integer ( kind = 4 ) FACE_LIST(N), the list of nodes making a face.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the face.
!
!    Output, integer ( kind = 4 ) INDEX, the results of the search.
!    0, there is no edge common to the face and the EDGE array.
!    nonzero, edge INDEX is common to the face and the EDGE array.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) edge_max

  integer ( kind = 4 ) edge(4,edge_max)
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) face_list(n)
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) index
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  index = 0

  if ( n <= 0 ) then
    return
  end if

  if ( edge_num <= 0 ) then
    return
  end if

  do j = 1, n

    if ( j == n ) then
      jp1 = 1
    else
      jp1 = j + 1
    end if

    n1 = face_list(j)
    n2 = face_list(jp1)

    do iedge = 1, edge_num

      if ( edge(1,iedge) == n2 .and. edge(2,iedge) == n1 ) then

        call i4vec_rotate ( n, 1 - j, face_list )

        index = iedge
        return

      else if ( edge(1,iedge) == n1 .and. edge(2,iedge) == n2 ) then

        call i4vec_rotate ( n, n - jp1, face_list )

        call i4vec_reverse ( n, face_list )

        index = iedge
        return

      end if

    end do
   
  end do

  return
end
subroutine edge_match_nodes ( edge, edge_max, edge_num, n1, n2, iedge )

!*****************************************************************************80
!
!! EDGE_MATCH_NODES seeks an edge of the form (N1,N2) or (N2,N1) in EDGE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EDGE(4,EDGE_MAX), contains edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer ( kind = 4 ) EDGE_MAX, the maximum number of edges.
!
!    Input, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
!    Input, integer ( kind = 4 ) N1, N2, two nodes that form an edge.
!
!    Output, integer ( kind = 4 ) IEDGE, the results of the search.
!    0, no matching edge was found.
!    nonzero, edge IEDGE of the EDGE array matches (N1,N2) or (N2,N1).
!
  implicit none

  integer ( kind = 4 ) edge_max

  integer ( kind = 4 ) edge(4,edge_max)
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  iedge = 0
  do i = 1, edge_num

    if ( ( n1 == edge(1,i) .and. n2 == edge(2,i) ) .or. &
         ( n2 == edge(1,i) .and. n1 == edge(2,i) ) ) then
      iedge = i
      return
    end if

  end do

  return
end
subroutine edge_null_delete ( cor3, cor3_max, face, face_max, face_num, &
  face_order, order_max, vertex_normal )

!*****************************************************************************80
!
!! EDGE_NULL_DELETE deletes face edges with zero length.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3(3,cor3_max)
  logical, parameter :: debug = .false.
  real ( kind = 4 ) distsq
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) edge_num_del
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) face_order2
  integer ( kind = 4 ) face2(100)
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) inode
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode
  integer ( kind = 4 ) jvert
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  real ( kind = 4 ) vertex_normal2(3,100)

  edge_num = 0
  edge_num_del = 0
!
!  Consider each face.
!
  do iface = 1, face_num
!
!  Consider each pair of consecutive vertices.
!
    face_order2 = 0

    do ivert = 1, face_order(iface)

      edge_num = edge_num + 1

      jvert = ivert + 1
      if ( face_order(iface) < jvert ) then
        jvert = 1
      end if

      inode = face(ivert,iface)
      jnode = face(jvert,iface)

      distsq = sum (  ( cor3(1:3,inode) - cor3(1:3,jnode) )**2 )

      if ( distsq /= 0.0E+00 ) then
        face_order2 = face_order2 + 1
        face2(face_order2) = face(ivert,iface)
        vertex_normal2(1:3,face_order2) = vertex_normal(1:3,ivert,iface)
      else
        edge_num_del = edge_num_del + 1
      end if

    end do

    face_order(iface) = face_order2
    do ivert = 1, face_order(iface)
      face(ivert,iface) = face2(ivert)
      vertex_normal(1:3,ivert,iface) = vertex_normal2(1:3,ivert)
    end do

  end do

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EDGE_NULL_DELETE:'
    write ( *, '(a,i6,a)' ) '  There are a total of ', edge_num, ' edges.'
    write ( *, '(a,i6,a)' ) '  Of these, ', edge_num_del, &
      ' were of zero length, and deleted.'
  end if

  return
end
function enorm_nd ( n, x )

!*****************************************************************************80
!
!! ENORM_ND computes the Euclidean norm of a vector in ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input, real ( kind = 4 ) X(N), the coordinates of the vector.
!
!    Output, real ( kind = 4 ) ENORM_ND, the Euclidean norm of the vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 4 ) enorm_nd
  integer ( kind = 4 ) i
  real ( kind = 4 ) x(n)

  enorm_nd = 0.0E+00

  do i = 1, n
    enorm_nd = enorm_nd + x(i) * x(i)
  end do

  enorm_nd = sqrt ( enorm_nd )
 
  return
end
function enorm0_3d ( x0, y0, z0, x1, y1, z1 )

!*****************************************************************************80
!
!! ENORM0_3D computes the Euclidean norm of (P1-P0) in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X0, Y0, Z0, X1, Y1, Z1, the coordinates of the points 
!    P0 and P1.
!
!    Output, real ( kind = 4 ) ENORM0_3D, the Euclidean norm of (P1-P0).
!
  implicit none

  real ( kind = 4 ) enorm0_3d
  real ( kind = 4 ) x0
  real ( kind = 4 ) x1
  real ( kind = 4 ) y0
  real ( kind = 4 ) y1
  real ( kind = 4 ) z0
  real ( kind = 4 ) z1

  enorm0_3d = sqrt ( ( x1 - x0 )**2 + ( y1 - y0 )**2 + ( z1 - z0 )**2 )
 
  return
end
subroutine face_area_set ( cor3, cor3_max, face, face_area, face_max, &
  face_num, face_order, order_max )

!*****************************************************************************80
!
!! FACE_AREA_SET computes the area of the faces.
!
!  Discussion:
!
!    The area is the sum of the areas of the triangles formed by
!    node N with consecutive pairs of nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Output, real ( kind = 4 ) FACE_AREA(FACE_MAX), the area of each face.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices allowed per face.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) alpha
  real ( kind = 4 ) area_max
  real ( kind = 4 ) area_min
  real ( kind = 4 ) area_tri
  real ( kind = 4 ) base
  real ( kind = 4 ) cor3(3,cor3_max)
  logical, parameter :: debug = .false.
  real ( kind = 4 ) dot
  integer ( kind = 4 ) face(order_max,face_max)
  real ( kind = 4 ) face_area(face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_num_del
  integer ( kind = 4 ) face_order(face_max)
  real ( kind = 4 ) height
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) iface
  real ( kind = 4 ) tol
  real ( kind = 4 ) x1
  real ( kind = 4 ) x2
  real ( kind = 4 ) x3
  real ( kind = 4 ) y1
  real ( kind = 4 ) y2
  real ( kind = 4 ) y3
  real ( kind = 4 ) z1
  real ( kind = 4 ) z2
  real ( kind = 4 ) z3

  do iface = 1, face_num

    face_area(iface) = 0.0E+00

    do i = 1, face_order(iface) - 2

      i1 = face(i,iface)
      i2 = face(i+1,iface)
      i3 = face(i+2,iface)

      x1 = cor3(1,i1)
      y1 = cor3(2,i1)
      z1 = cor3(3,i1)

      x2 = cor3(1,i2)
      y2 = cor3(2,i2)
      z2 = cor3(3,i2)

      x3 = cor3(1,i3)
      y3 = cor3(2,i3)
      z3 = cor3(3,i3)
!
!  Find the projection of (P3-P1) onto (P2-P1).
!
      dot = ( x2 - x1 ) * ( x3 - x1 ) + &
            ( y2 - y1 ) * ( y3 - y1 ) + &
            ( z2 - z1 ) * ( z3 - z1 )

      base = sqrt ( ( x2 - x1 )**2 + ( y2 - y1 )**2 + ( z2 - z1 )**2 )
!
!  The height of the triangle is the length of (P3-P1) after its
!  projection onto (P2-P1) has been subtracted.
!
      if ( base == 0.0E+00 ) then

        height = 0.0E+00

      else

        alpha = dot / base**2
  
        height = sqrt ( ( x3 - x1 - alpha * ( x2 - x1 ) )**2 + &
                        ( y3 - y1 - alpha * ( y2 - y1 ) )**2 + &
                        ( z3 - z1 - alpha * ( z2 - z1 ) )**2 )

      end if

      area_tri = 0.5E+00 * base * height

      face_area(iface) = face_area(iface) + area_tri

    end do

  end do

  area_min = minval ( face_area(1:face_num) )
  area_max = maxval ( face_area(1:face_num) )

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FACE_AREA_SET:'
    write ( *, '(a,g14.6)' ) '  Minimum face area is ', area_min
    write ( *, '(a,g14.6)' ) '  Maximum face area is ', area_max
  end if

  tol = area_max / 10000.0E+00

  if ( area_min < tol ) then

    face_num_del = 0

    do iface = 1, face_num
      if ( face_area(iface) < tol ) then
        face_order(iface) = 0
        face_num_del = face_num_del + 1
      end if
    end do

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FACE_AREA_SET:'
      write ( *, '(a,i6,a)' ) '  Marked ', face_num_del, &
        ' tiny faces for deletion.'
    end if

  end if

  return
end
subroutine face_check ( edge, edge_max, edge_num, face, face_material, &
  face_max, face_normal, face_num, face_object, face_order, face_rank, &
  face_tier, object_num, order_max, vertex_material, vertex_normal )

!*****************************************************************************80
!
!! FACE_CHECK checks and analyzes a set of faces.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) EDGE(4,EDGE_MAX), contains edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, or 0 if the edge is used once.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_NUM), the nodes making faces
!
!    Input, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Output, integer ( kind = 4 ) FACE_OBJECT(FACE_NUM), describes the objects.
!    FACE_OBJECT(I) is the index of the edge-connected "object" to 
!    which face I belongs.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_NUM), the number of nodes per face.
!
!    Output, integer ( kind = 4 ) FACE_RANK(FACE_NUM), is an ordered list of faces.
!    FACE_RANK(1) is the index of the face in the first tier of the 
!    first object, followed by second tier faces, and so on until
!    object one is complete.  Object two follows, and so on.
!
!    Output, integer ( kind = 4 ) FACE_TIER(FACE_NUM).  FACE_TIER(I) is the "tier"
!    of face I in its object.  The seed of the object has tier 1,
!    the neighbors of the seed have tier 2, and so on.
!
!    Input, integer ( kind = 4 ) EDGE_MAX, the maximum number of edges.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, is the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Output, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Output, integer ( kind = 4 ) OBJECT_NUM, the number of objects.
!
!    Input, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  implicit none

  integer ( kind = 4 ) edge_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) face_max

  integer ( kind = 4 ) edge(4,edge_max)
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_fixed
  integer ( kind = 4 ) face_material(face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_object(face_max)
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) face_rank(face_max)
  integer ( kind = 4 ) face_tier(face_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) object_num
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
!
!  Organize the faces into layered objects.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FACE_CHECK:'
  write ( *, '(a)' ) '  Determine edge-connected objects.'

  call object_build ( face, face_num, face_object, face_order, face_rank, &
    face_tier, object_num, order_max )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of objects = ', object_num

  if ( face_num <= 20 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Face, Object, Tier'
    write ( *, '(a)' ) ' '

    do i = 1, face_num
      write ( *, '(3i6)' ) i, face_object(i), face_tier(i)
    end do

  end if

  if ( face_num <= 20 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Preferred order:'
    write ( *, '(a)' ) '  Order, Face'
    write ( *, '(a)' ) ' '
    do i = 1, face_num
      write ( *, '(i6,i6)' ) i, face_rank(i)
    end do
  end if
!
!  Reorder the faces by object and tier.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Reorder the faces.'

  call face_sort ( face, face_material, face_max, face_normal, face_num, &
    face_object, face_order, face_tier, order_max, &
    vertex_material, vertex_normal )

  if ( face_num <= 20 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Face, Label, Object, Tier'
    write ( *, '(a)' ) ' '
    do i = 1, face_num
      write ( *, '(4i6)' ) i, face_rank(i), face_object(i), face_tier(i)
    end do
  end if
!
!  Construct the edge list.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Construct the edge list.'
  write ( *, '(a)' ) '(While doing so, check for edges used more'
  write ( *, '(a)' ) 'than twice.)'

  call face_to_edge ( edge, edge_max, edge_num, face, face_num, &
    face_order, ierror, order_max )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FACE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  FACE_TO_EDGE failed.'
    return
  end if

  if ( face_num <= 20 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Edge, Node1, Node2, Face1, Face2, Tier, Object'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' I, node1(i), node2(i), face1(i), face2(i)'
    write ( *, '(a)' ) ' '

    do i = 1, edge_num
      write ( *, '(10i3)' ) i, edge(1:4,i)
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Face, Order, Nodes'
    write ( *, '(a)' ) ' '
    do i = 1, face_num
      write ( *, '(10i3)' ) i, face_order(i), &
        ( face(j,i), j = 1, face_order(i) )
    end do
  end if
!
!  Now force faces to have a consistent orientation.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Force faces to consistent orientation.'
  
  call face_flip ( edge, edge_max, edge_num, face, face_fixed, face_num, &
    face_order, order_max )

  if ( face_num <= 20 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Face, Order, Nodes'
    write ( *, '(a)' ) ' '
    do i = 1, face_num
      write ( *, '(10i3)' ) i, face_order(i), &
        ( face(j,i), j = 1, face_order(i) )
    end do
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'List boundary edges.'

  call edge_bound ( edge, edge_max, edge_num )

  return
end
subroutine face_flip ( edge, edge_max, edge_num, face, face_fixed, face_num, &
  face_order, order_max )

!*****************************************************************************80
!
!! FACE_FLIP flips faces to achieve a consistent orientation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EDGE(4,EDGE_MAX), contains edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_NUM), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_NUM), the number of nodes per face.
!
!    Input, integer ( kind = 4 ) EDGE_MAX, the maximum number of edges.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Output, integer ( kind = 4 ) FACE_FIXED, the number of bad faces that were found.
!
!    Input, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
  implicit none

  integer ( kind = 4 ) edge_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) face_num

  integer ( kind = 4 ) edge(4,edge_max)
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) f1
  integer ( kind = 4 ) f2
  integer ( kind = 4 ) face(order_max,face_num)
  integer ( kind = 4 ) face_fixed
  integer ( kind = 4 ) face_order(face_num)
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
!
  face_fixed = 0

  do iedge = 1, edge_num

    n1 = edge(1,iedge)
    n2 = edge(2,iedge)
    f1 = edge(3,iedge)
    f2 = edge(4,iedge)
!
!  For now, just whine unless (N1,N2) is positive in F1 and negative in F2.
!
    if ( f1 /= 0 ) then

      do j = 1, face_order(f1)

        if ( j < face_order(f1) ) then
          jp1 = j + 1
        else
          jp1 = j
        end if

        m1 = face(j,f1)
        m2 = face(jp1,f1)

        if ( m1 == n1 .and. m2 == n2 ) then
          exit
        else if ( m1 == n2 .and. m2 == n1 ) then
          face_fixed = face_fixed + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FACE_FLIP - Warning!'
          write ( *, '(a)' ) '  Bad orientation:'
          write ( *, '(a,i6)' ) '  Face ', f1
          write ( *, '(a,i6)' ) '  Side ', j
          exit
        end if

      end do

    end if

    if ( f2 /= 0 ) then

      do j = 1, face_order(f2)

        if ( j < face_order(f2) ) then
          jp1 = j + 1
        else
          jp1 = j
        end if

        m1 = face(j,f2)
        m2 = face(jp1,f2)

        if ( m1 == n2 .and. m2 == n1 ) then
          exit
        else if ( m1 == n1 .and. m2 == n2 ) then
          face_fixed = face_fixed + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FACE_FLIP - Warning!'
          write ( *, '(a)' ) 'Bad orientation:'
          write ( *, '(a,i6)' ) '  Face ', f2
          write ( *, '(a,i6)' ) '  Side ', j
          exit
        end if

      end do

    end if

  end do

  write ( *, '(a,i6,a)' ) 'FACE_FLIP found ', face_fixed, &
    ' badly oriented faces.'

  return
end
subroutine face_normal_ave ( face_max, face_normal, face_num, face_order, &
  order_max, vertex_normal )

!*****************************************************************************80
!
!! FACE_NORMAL_AVE sets face normals as average of face vertex normals.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Output, real ( kind = 4 ) FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  implicit none

  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  logical, parameter :: debug = .false.
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) nfix
  real ( kind = 4 ) norm
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)

  if ( face_num <= 0 ) then
    return
  end if

  nfix = 0

  do iface = 1, face_num

    norm = sqrt ( sum ( face_normal(1:3,iface)**2 ) )

    if ( norm == 0.0E+00 ) then

      nfix = nfix + 1

      face_normal(1:3,iface) = 0.0E+00

      do ivert = 1, face_order(iface)
        face_normal(1:3,iface) = face_normal(1:3,iface) + &
          vertex_normal(1:3,ivert,iface)
      end do

      norm = sqrt ( sum ( face_normal(1:3,iface)**2 ) )

      if ( norm == 0.0E+00 ) then
        face_normal(1:3,iface) = 1.0E+00 / sqrt ( 3.0E+00 )
      else
        face_normal(1:3,iface) = face_normal(1:3,iface) / norm
      end if

    end if

  end do

  if ( 0 < nfix ) then
    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6,a)' ) 'FACE_NORMAL_AVE: Recomputed ', nfix, &
        ' face normals by averaging face vertex normals.'
    end if
  end if

  return
end
subroutine face_null_delete ( face, face_area, face_material, face_max, &
  face_num, face_order, order_max, vertex_material, vertex_normal )

!*****************************************************************************80
!
!! FACE_NULL_DELETE deletes faces of order less than 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input/output, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normal vectors
!    at vertices.  
!
  implicit none

  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  logical, parameter :: debug = .false.
  integer ( kind = 4 ) face(order_max,face_max)
  real ( kind = 4 ) face_area(face_max)
  integer ( kind = 4 ) face_material(face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_num2
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
!
!  Drop faces of order 0, 1 or 2.
!
  face_num2 = 0

  do iface = 1, face_num

    if ( 3 <= face_order(iface) ) then

      face_num2 = face_num2 + 1

      if ( face_num2 /= iface ) then

        face_area(face_num2) = face_area(iface)
        face_material(face_num2) = face_material(iface)
        face_order(face_num2) = face_order(iface)
        face(1:order_max,face_num2) = face(1:order_max,iface)
        vertex_material(1:order_max,face_num2) = &
          vertex_material(1:order_max,iface)
        vertex_normal(1:3,1:order_max,face_num2) = &
          vertex_normal(1:3,1:order_max,iface)

      end if

    end if

  end do

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FACE_NULL_DELETE'
    write ( *, '(a,i6,a)' ) '  There are a total of ', face_num, ' faces.'
    write ( *, '(a,i6,a)' ) '  Of these, ', face_num2, ' passed the order test.'
  end if

  face_num = face_num2 

  return
end
subroutine face_print ( cor3, cor3_max, face, face_index, face_material, &
  face_max, face_normal, face_num, face_order, order_max, vertex_material, &
  vertex_normal )

!*****************************************************************************80
!
!! FACE_PRINT prints out information about a face.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, real ( kind = 4 ) FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer ( kind = 4 ) IFACE, the face about which information is desired.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_index
  integer ( kind = 4 ) face_material(face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(order_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)

  if ( face_index < 1 .or. face_num < face_index ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FACE_PRINT - Error!'
    write ( *, '(a,i6)' ) '  Face indices must be between 1 and ', face_num
    write ( *, '(a,i6)' ) '  But your requested value was ', face_index
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FACE_PRINT'
  write ( *, '(a,i6)' ) '  Information about face ', face_index
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of vertices is ', face_order(face_index)
  write ( *, '(a,i6)' ) '  Face material is ', face_material(face_index)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vertex list:'
  write ( *, '(a)' ) '    Vertex #, Node #, Material #, X, Y, Z'
  write ( *, '(a)' ) ' '
  do ivert = 1, face_order(face_index)
    j = face(ivert,face_index)
    k = vertex_material(ivert,face_index)
    write ( *, '(3i8,3f10.4)' ) ivert, j, k, cor3(1,j), cor3(2,j), cor3(3,j)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Face normal vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(3f10.4)' ) face_normal(1:3,face_index)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vertex face normals:'
  write ( *, '(a)' ) ' '
  do ivert = 1, face_order(face_index)
    write ( *, '(i8,3f10.4)' ) ivert, vertex_normal(1:3,ivert,face_index)
  end do

  return
end
subroutine face_reverse_order ( cor3_max, cor3_normal, cor3_num, face, &
  face_max, face_normal, face_num, face_order, order_max, &
  vertex_material, vertex_normal, vertex_tex_uv )

!*****************************************************************************80
!
!! FACE_REVERSE_ORDER reverses the order of the nodes in each face.
!
!  Discussion:
!
!    Reversing the order of the nodes requires that the normal vectors
!    be reversed as well, so this routine will automatically reverse
!    the normals associated with nodes, vertices and faces.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of nodes.
!
!    Input/output, real ( kind = 4 ) COR3_NORMAL(3,COR3_MAX), normals at nodes.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of nodes.
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input/output, real ( kind = 4 ) FACE_NORMAL(3,FACE_MAX), normals at faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input/output, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input/output, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals 
!    at vertices.  
!
!    Input/output, real ( kind = 4 ) VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex textures.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3_normal(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) face(order_max,face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 4 ) temp
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  real ( kind = 4 ) vertex_tex_uv(2,order_max,face_max)
!
  do iface = 1, face_num

    m = face_order(iface)

    do ivert = 1, m/2

      call i4_swap ( face(ivert,iface), face(m+1-ivert,iface) )
      call i4_swap ( vertex_material(ivert,iface), &
                    vertex_material(m+1-ivert,iface) )

      do j = 1, 3
        call r4_swap ( vertex_normal(j,ivert,iface), &
                      vertex_normal(j,m+1-ivert,iface) )
      end do

      do j = 1, 2
        call r4_swap ( vertex_tex_uv(j,ivert,iface), &
                      vertex_tex_uv(j,m+1-ivert,iface) )
      end do

    end do

  end do

  cor3_normal(1:3,1:cor3_num) = - cor3_normal(1:3,1:cor3_num)
  face_normal(1:3,1:face_num) = - face_normal(1:3,1:face_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FACE_REVERSE_ORDER'
  write ( *, '(a)' ) '  Each list of nodes defining a face'
  write ( *, '(a)' ) '  has been reversed; related information,'
  write ( *, '(a)' ) '  including normal vectors, was also updated.'

  return
end
subroutine face_sort ( face, face_material, face_max, face_normal, face_num, &
  face_object, face_order, face_tier, order_max, vertex_material, &
  vertex_normal )

!*****************************************************************************80
!
!! FACE_SORT renumbers the faces in order of object and tier.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_NUM), the nodes making faces.
!
!    Input/output, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input/output, integer ( kind = 4 ) FACE_OBJECT(FACE_NUM), describes the objects.
!    FACE_OBJECT(I) is the index of the edge-connected "object" to 
!    which face I belongs.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(FACE_NUM), the number of nodes per face.
!
!    Input/output, integer ( kind = 4 ) FACE_TIER(FACE_NUM).  FACE_TIER(I) is the "tier"
!    of face I in its object.  The seed of the object has tier 1,
!    the neighbors of the seed have tier 2, and so on.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, is the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input/output, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  implicit none

  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) face_max

  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_object(face_max)
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) face_tier(face_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) jface
  real ( kind = 4 ) temp
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)

  iface = 0
  jface = 0
  indx = 0
  isgn = 0

  do

    call sort_heap_external ( face_num, indx, iface, jface, isgn )
!
!  Interchange faces IFACE and JFACE.
!
    if ( 0 < indx ) then

      do i = 1, order_max
        call i4_swap ( face(i,iface), face(i,jface) )
      end do

      call i4_swap ( face_material(iface),    face_material(jface) )
      call i4_swap ( face_object(iface), face_object(jface) )
      call i4_swap ( face_order(iface),  face_order(jface) )
      call i4_swap ( face_tier(iface),   face_tier(jface) )

      do i = 1, 3
        call r4_swap ( face_normal(i,iface), face_normal(i,jface) )
      end do

      do ivert = 1, order_max
        call i4_swap ( vertex_material(ivert,iface), vertex_material(ivert,jface) )
      end do

      do i = 1, 3
        do ivert = 1, order_max
          call r4_swap ( vertex_normal(i,ivert,iface), vertex_normal(i,ivert,jface) )
        end do
      end do
!
!  Compare faces IFACE and JFACE.
!
    else if ( indx < 0 ) then

      if ( ( face_object(iface) < face_object(jface) ) .or. &
           ( face_object(iface) == face_object(jface) .and. &
             face_tier(iface) < face_tier(jface) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine face_subset ( cor3, cor3_max, cor3_num, face, face_material, &
  face_max, face_normal, face_num, face_order, ierror, line_num, list, &
  order_max, vertex_material, vertex_normal )

!*****************************************************************************80
!
!! FACE_SUBSET selects a subset of the current faces as the new object.
!
!  Warning:
!
!    The original graphic object is overwritten by the new one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input/output, real ( kind = 4 ) FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Workspace, integer LIST(COR3_MAX), contains the indices of the points
!    to be copied from the old graphics object to the new one.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input/output, integer ( kind = 4 ) LINE_NUM, the number of lines.  
!    This routine resets LINE_NUM to zero, since we will be dropping
!    as many points as possible.
!
!    Input/output, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input/output, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals 
!    at vertices.  
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) cor3_num2
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) iface1
  integer ( kind = 4 ) iface2
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) list(cor3_max)
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)

  ierror = 0

  line_num = 0
!
!  Get the first and last faces to save, IFACE1 and IFACE2.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Enter lowest face number to save,'
  write ( *, '(a,i6)' ) 'between 1 and ', face_num
  read ( *, * ) iface1

  if ( iface1 < 1 .or. face_num < iface1 ) then
    write ( *, '(a)' ) 'Illegal choice!'
    ierror = 1
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Enter highest face number to save'
  write ( *, '(a,i6,a,i6)' ) 'between ', iface1, ' and ', face_num
  read ( *, * ) iface2

  if ( iface2 < iface1 .or. face_num < iface2 ) then
    write ( *, '(a)' ) 'Illegal choice!'
    ierror = 1
    return
  end if

  inc = iface1 - 1
!
!  "Slide" the data for the saved faces down the face arrays.
!
  do iface = 1, iface2 + 1 - iface1

    face_material(iface) = face_material(iface+inc)
    face_order(iface) = face_order(iface+inc)

    do ivert = 1, order_max
      face(ivert,iface) = face(ivert,iface+inc)
      vertex_material(ivert,iface) = vertex_material(ivert,iface+inc)
      vertex_normal(1:3,ivert,iface) = vertex_normal(1:3,ivert,iface+inc)
    end do

    face_normal(1:3,iface) = face_normal(1:3,iface+inc)

  end do
!
!  Now reset the number of faces.
!
  face_num = iface2 + 1 - iface1
!
!  Now, for each point I, set LIST(I) = J if point I is the J-th
!  point we are going to save, and 0 otherwise.  Then J will be
!  the new label of point I.
!
  list(1:cor3_num) = 0
  cor3_num2 = 0

  do iface = 1, face_num

    do ivert = 1, face_order(iface)

      j = face(ivert,iface)

      if ( list(j) == 0 ) then
        cor3_num2 = cor3_num2 + 1
        list(j) = cor3_num2
      end if

    end do

  end do
!
!  Now make the nonzero list entries rise in order, so that
!  we can compress the COR3 data in a minute.
!
  cor3_num2 = 0
  do i = 1, cor3_num
    if ( list(i) /= 0 ) then
      cor3_num2 = cor3_num2 + 1
      list(i) = cor3_num2
    end if
  end do
!
!  Relabel the FACE array with the new node indices.
!
  do iface = 1, face_num
    do ivert = 1, face_order(iface)
      j = face(ivert,iface)
      face(ivert,iface) = list(j)
    end do
  end do
!
!  Rebuild the COR3 array by sliding data down.
!
  do i = 1, cor3_num
    k = list(i)
    if ( k /= 0 ) then
      cor3(1:3,k) = cor3(1:3,i)
    end if
  end do

  cor3_num = cor3_num2

  return
end
subroutine face_to_edge ( edge, edge_max, edge_num, face, face_num, &
  face_order, ierror, order_max )

!*****************************************************************************80
!
!! FACE_TO_EDGE converts face data to edge data.
!
!  Discussion:
!
!    The computation will fail if:
!
!      More than two faces claim to share an edge (Node1,Node2).
!      Not enough storage is set aside by EDGE_MAX.
!
!    If is expected that the edge (Node1,Node2) in Face1 is traversed in
!    the opposite sense, as (Node2,Node1), in Face2.  If this is not the
!    case, then some faces may need to be reoriented, but that will not
!    affect the computation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) EDGE(4,EDGE_MAX), contains edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, or 0 if the edge is used once.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_NUM), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_NUM), the number of nodes per face.
!
!    Output, integer ( kind = 4 ) IERROR, error flag: 0 = no error, nonzero = error.
!
!    Input, integer ( kind = 4 ) EDGE_MAX, the maximum number of edges.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Output, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
  implicit none

  integer ( kind = 4 ) edge_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) face_num

  integer ( kind = 4 ) edge(4,edge_max)
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) face(order_max,face_num)
  integer ( kind = 4 ) face_order(face_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) index
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
!
!  Initialize.
!
  ierror = 0

  edge(1:4,1:edge_max) = 0

  edge_num = 0
!
!  Consider face #I.
!
  do iface = 1, face_num
!
!  Seek an edge of face IFACE that already occurs in the edge list.
!  If there is one, then slide and reverse the entries in FACE(*,IFACE)
!  so that that edge occurs first, and in the opposite sense to its
!  occurrence in the edge list.
!
    call edge_match_face ( edge, edge_max, edge_num, face(1,iface), &
      face_order(iface), index )
!
!  Now, in any case, we know that the first two nodes in FACE(*,IFACE)
!  are the negative of an existing edge, or no nodes in FACE(*,IFACE)
!  occur in any existing edge.
!
    do j = 1, face_order(iface)

      n1 = face(j,iface)
      if ( j == face_order(iface) ) then
        jp1 = 1
      else
        jp1 = j + 1
      end if

      n2 = face(jp1,iface)

      call edge_match_nodes ( edge, edge_max, edge_num, n1, n2, iedge )

      if ( iedge == 0 ) then

        call edge_add_nodes ( edge, edge_max, edge_num, iface, n1, n2, ierror )

        if ( ierror /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FACE_TO_EDGE - Fatal error!'
          write ( *, '(a)' ) '  EDGE_ADD_NODES failed.'
          ierror = 1
          return
        end if

      else if ( edge(4,iedge) == 0 ) then

        edge(4,iedge) = iface

      else 

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FACE_TO_EDGE - Fatal error!'
        write ( *, '(a,2i6)' ) '  Edge between nodes ', edge(1,iedge), edge(2,iedge)
        write ( *, '(a)' ) '  is used at least 3 times, by faces:'
        write ( *, '(3i6)' ) edge(3,iedge), edge(4,iedge), iface
        ierror = 1
        return

      end if

    end do
  end do

  return
end
subroutine face_to_line ( debug, face, face_max, face_num, face_order, &
  line_dex, line_material, line_max, line_num, line_prune, &
  order_max, vertex_material )

!*****************************************************************************80
!
!! FACE_TO_LINE converts face information to line information.
!
!  Discussion:
!
!    In some cases, the graphic information represented by polygonal faces
!    must be converted to a representation based solely on line segments.
!    This is particularly true if a VLA file is being written.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input/output, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, 
!    terminated by -1.
!
!    Output, integer ( kind = 4 ) LINE_MATERIAL(LINE_MAX), material index for line.
!
!    Input, integer ( kind = 4 ) LINE_PRUNE, pruning option.
!    0, no pruning, draw every line.
!    nonzero, prune.  Only draw the line from node I to node J if I < J.
!    This should cut down on repeated drawing of lines in the common
!    case of a face mesh where each line is drawn twice, once with positive
!    and once with negative orientation.  In other cases, pruning
!    may omit some lines that only occur with negative orientation.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input/output, integer ( kind = 4 ) LINE_NUM, the number of line data items.
!
!    Input, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  implicit none

  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) order_max

  logical debug
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) jcor3
  integer ( kind = 4 ) jvert
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) line_prune
  integer ( kind = 4 ), parameter :: OFFSET = 1
  integer ( kind = 4 ) vertex_material(order_max,face_max)
!
!  Case 1: 
!  No line pruning.
!
  if ( line_prune == 0 ) then

    do iface = 1, face_num

      do ivert = 1, face_order(iface)

        icor3 = face(ivert,iface)
 
        line_num = line_num + 1
        if ( line_num <= line_max ) then
          line_dex(line_num) = icor3
          line_material(line_num) = vertex_material(ivert,iface)
        end if

      end do

      ivert = 1
      icor3 = face(ivert,iface)

      line_num = line_num + 1
      if ( line_num <= line_max ) then
        line_dex(line_num) = icor3
        line_material(line_num) = vertex_material(ivert,iface)
      end if

      line_num = line_num + 1
      if ( line_num <= line_max ) then
        line_dex(line_num) = -1 + OFFSET
        line_material(line_num) = -1 + OFFSET
      end if

    end do
!
!  Case 2: 
!    Simple-minded line pruning.
!    Only draw line (I,J) if I < J.
!
  else

    do iface = 1, face_num

      do ivert = 1, face_order(iface)

        icor3 = face(ivert,iface)

        if ( ivert < face_order(iface) ) then
          jvert = ivert + 1
        else
          jvert = 1
        end if

        jcor3 = face(jvert,iface)

        if ( icor3 < jcor3 ) then

          if ( line_num + 3 <= line_max ) then

            line_num = line_num + 1
            line_dex(line_num) = icor3
            line_material(line_num) = vertex_material(ivert,iface)
 
            line_num = line_num + 1
            line_dex(line_num) = jcor3
            line_material(line_num) = vertex_material(jvert,iface)

            line_num = line_num + 1
            line_dex(line_num) = -1 + OFFSET
            line_material(line_num) = -1 + OFFSET

          end if

        end if

      end do

    end do

  end if

  if ( debug ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FACE_TO_LINE:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I, LINE_DEX(I), LINE_MATERIAL(I)'
    write ( *, '(a)' ) ' '

    do i = 1, line_num
      write ( *, '(i6,2x,i6,2x,i6)' ) i, line_dex(i), line_material(i)
    end do

  end if

  return
end
subroutine face_touch ( face, face_order, order_max, face_num, iface, jface, &
  touch )

!*****************************************************************************80
!
!! FACE_TOUCH reports whether two polygonal faces touch.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_NUM), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_NUM), the number of nodes per face.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, is the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) IFACE, JFACE, the faces to be checked.
!
!    Output, integer ( kind = 4 ) TOUCH:
!     0, the faces do not touch;
!    +1, the faces touch, both using an arc in the same direction;
!    -1, the faces touch, using an arc in opposite directions.
!
  implicit none

  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) face_num

  integer ( kind = 4 ) face(order_max,face_num)
  integer ( kind = 4 ) face_order(face_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jface
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mp1
  integer ( kind = 4 ) mm1
  integer ( kind = 4 ) n
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) touch

  touch = 0
!
!  Arc N1-N2 on IFACE must be matched by arc N1-N2 or N2-N1 on JFACE.
!
  do i = 1, face_order(iface)

    n = face(i,iface)
    if ( i < face_order(iface) ) then
      np1 = face(i+1,iface)
    else
      np1 = face(1,iface)
    end if

    do j = 1, face_order(jface)

      m = face(j,jface)
      if ( j < face_order(jface) ) then
        mp1 = face(j+1,jface)
      else
        mp1 = face(1,jface)
      end if
      if ( 1 < j ) then
        mm1 = face(j-1,jface)
      else
        mm1 = face(face_order(jface),jface)
      end if

      if ( n == m ) then
        if ( np1 == mp1 ) then
          touch = + 1
          return
        else if ( np1 == mm1 ) then
          touch = - 1
          return
        end if
      end if

    end do
  end do

  return
end
subroutine file_get_next_word ( iunit, word, text, num_text, ierror )

!*****************************************************************************80
!
!! FILE_GET_NEXT_WORD returns the next word and trailing context from a file.
!
!  Discussion:
!
!    The file should have been opened before calling this routine.
!    The file should contain ASCII text, which can be thought of as
!    words separated by one or more blanks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) IUNIT, the unit number associated with the file.
!
!    Output, character ( len = * ) WORD, the next word in the file.  If the
!    current line of the file is blank, or if the file has been exhausted,
!    WORD will be set to ' '.
!
!    Input/output, character ( len = * ) TEXT, the remaining text of the line
!    that contains the information in WORD.  On each call, the next word
!    in TEXT is extracted until TEXT is empty, when it is refilled by
!    reading another line from the file.  Because TEXT contains information
!    needed by this routine, it should not be altered by the user
!    between calls.
!
!    Input/output, integer ( kind = 4 ) NUM_TEXT, the number of lines read from the file.
!    Before the first call to this routine, the user should set NUM_TEXT
!    to 0.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error, another word was read, and returned in WORD.
!    1, end of file.  WORD and TEXT were set to ' '.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) num_text
  character ( len = * ) text
  character ( len = * ) word

  ierror = 0
!
!  If NUM_TEXT is zero, then initialize TEXT.
!
  if ( num_text <= 0 ) then
    num_text = 0
    text = ' '
  end if
!
!  If TEXT is blank, try to read a new line from the file.
!
  if ( text == ' ' ) then

    read ( iunit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      ierror = 1
      word = ' '
      text = ' '
      return
    end if

    num_text = num_text + 1

    if ( text == ' ' ) then
      word = ' '
      return
    end if

  end if
!
!  Extract the next word from TEXT into WORD and return.
!
  lenc = len_trim ( text )
!
!  Find ILO, the index of the first nonblank in TEXT.
!
  ilo = 1

  do while ( text(ilo:ilo) == ' ' )
    ilo = ilo + 1
  end do
!
!  Find IHI, the index of the last consecutive nonblank after the one at ILO.
!
  ihi = ilo

  do while ( ihi+1 <= lenc )
    if ( text(ihi+1:ihi+1) == ' ' ) then
      exit
    end if
    ihi = ihi + 1
  end do
!
!  Set WORD.
!
  word = text(ilo:ihi)
!
!  Slide TEXT to the left.
!
  if ( ihi+1 <= lenc ) then
    text = text(ihi+1:)
  else
    text = ' '
  end if

  return
end
subroutine file_name_ext_get ( filnam, i, j )

!*****************************************************************************80
!
!! FILE_NAME_EXT_GET determines the "extension" of a file name.
!
!  Discussion:
!
!    The "extension" of a filename is the string of characters
!    that appears after the LAST period in the name.  A file
!    with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!    Blanks are unusual in filenames.  This routine ignores all
!    trailing blanks, but will treat initial or internal blanks
!    as regular characters acceptable in a file name.
!
!  Example:
!
!    FILNAM   I  J
!
!    bob.for    5  7
!    N.B.C.D    7  7
!    Naomi.     0  0
!    Arthur     0  0
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
!    Input, character ( len = * ) FILNAM, a file name to be examined.
!
!    Output, integer ( kind = 4 ) I, J, the indices of the first and last characters
!    in the file extension.
!
!    If at least one period occurs in the filename, and at least one
!    nonblank character follows that period, then I will be the index
!    of the first character after the period, and J the index of the
!    last nonblank character after the period.  The extension is
!    therefore equal to FILNAM(I:J).
!
!    Otherwise, I and J will be returned as 0, indicating that the file
!    has no extension.
!
  implicit none

  character ( len = * ) filnam
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s_index_last

  i = s_index_last ( filnam, '.' )

  if ( i /= 0 ) then

    j = len_trim ( filnam )

    if ( i == j ) then
      i = 0
      j = 0
    else
      i = i + 1
    end if

  else

    j = 0

  end if

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
subroutine hello ( cor3_max, face_max, line_max, material_max, order_max, &
  point_max, texture_max )

!*****************************************************************************80
!
!! HELLO prints out a message about the program.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of lines.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum vertices per face.
!
!    Input, integer ( kind = 4 ) POINT_MAX, the maximum number of points to display.
!
!    Input, integer ( kind = 4 ) TEXTURE_MAX, the maximum number of textures.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) point_max
  integer ( kind = 4 ) texture_max

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Hello:  This is IVRead,'
  write ( *, '(a)' ) '  a program which can convert some files from'
  write ( *, '(a)' ) '  some 3D graphics format to some others:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    ".ase"  3D Studio Max ASCII export;'
  write ( *, '(a)' ) '    ".byu"  Movie.BYU surface geometry;'
  write ( *, '(a)' ) '    ".dxf"  AutoCAD DXF;'
  write ( *, '(a)' ) '    ".hrc"  SoftImage hierarchy;'
  write ( *, '(a)' ) '    ".iv"   SGI Open Inventor;'
  write ( *, '(a)' ) '    ".obj"  WaveFront Advanced Visualizer;'
  write ( *, '(a)' ) '    ".off"  Geomview OFF file;'  
  write ( *, '(a)' ) '    ".oogl" OOGL file (input only);'  
  write ( *, '(a)' ) '    ".pov"  Persistence of Vision (output only);'
  write ( *, '(a)' ) '    ".ps"   PostScript (output only)(NOT READY);'
  write ( *, '(a)' ) '    ".smf"  Michael Garland''s format;'
  write ( *, '(a)' ) '    ".stl"  ASCII StereoLithography;'
  write ( *, '(a)' ) '    ".stla" ASCII StereoLithography;'
  write ( *, '(a)' ) '    ".tec"  TECPLOT (output only);'
  write ( *, '(a)' ) '    ".tri"  [Greg Hood triangles, ASCII];'
  write ( *, '(a)' ) '    ".tria" [Greg Hood triangles, ASCII];'
  write ( *, '(a)' ) '    ".ts"   Mathematica ThreeScript (output only);'
  write ( *, '(a)' ) '    ".3s"   Mathematica ThreeScript (output only);'
  write ( *, '(a)' ) '    ".txt"  Text (output only);'
  write ( *, '(a)' ) '    ".ucd"  AVS unstructured cell data (output only);'
  write ( *, '(a)' ) '    ".vla"  VLA; (points and lines);' 
  write ( *, '(a)' ) '    ".wrl"  VRML;'
  write ( *, '(a)' ) '    ".xgl"  XGL (output only) (DEVELOPMENT)'
  write ( *, '(a)' ) '    ".xyz"  XYZ (points and lines);'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Current limits:'
  write ( *, '(a)') ' '
  write ( *, '(i8,a)' ) cor3_max,  ' points;'
  write ( *, '(i8,a)' ) line_max,  ' line items;'
  write ( *, '(i8,a)' ) face_max,  ' faces.'
  write ( *, '(a)' ) ' '
  write ( *, '(i8,a)' ) order_max, ' vertices per face;'
  write ( *, '(i8,a)' ) point_max, ' points to display;'
  write ( *, '(i8,a)' ) material_max,   ' materials;'
  write ( *, '(i8,a)' ) texture_max,   ' textures.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Last modified: 04 June 2002.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Send problem reports to burkardt@iastate.edu.'
 
  return
end
subroutine help

!*****************************************************************************80
!
!! HELP prints out a help message about the interactive commands.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 June 1999
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HELP:'
  write ( *, '(a)' ) '  Batch commands to convert IN_FILE to OUT_FILE:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ivread  in_file  out_file'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ivread  -rn  in_file  out_file'
  write ( *, '(a)' ) '    Reverse normals before output.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ivread  -rf  in_file  out_file'
  write ( *, '(a)' ) '    Reverse faces and normals before output.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HELP:'
  write ( *, '(a)' ) '  These are legal interactive commands:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  < in_file    Read data from a file.'
  write ( *, '(a)' ) '  << out_file  Append data from a file'
  write ( *, '(a)' ) '                 to current information.'
  write ( *, '(a)' ) '  > out_file   Write data out to a file.'
  write ( *, '(a)' ) '  B            Switch byte swapping option.'
  write ( *, '(a)' ) '  D            Switch debug option.'
  write ( *, '(a)' ) '  F            Print info about one face.'
  write ( *, '(a)' ) '  H	     Print this help message.'
  write ( *, '(a)' ) '  I	     Info, print out recent changes.'
  write ( *, '(a)' ) '  LINE_PRUNE   Set FACE_TO_LINE pruning option.'
  write ( *, '(a)' ) '  LINES        Convert faces to lines.'
  write ( *, '(a)' ) '  N            Recompute normal vectors.'
  write ( *, '(a)' ) '  O            Use an average node normal.'
  write ( *, '(a)' ) '  Q            Quit.'
  write ( *, '(a)' ) '  REVERSE      Reverse the normal vectors.'
  write ( *, '(a)' ) '  RELAX        Smooth surface via relaxation.'
  write ( *, '(a)' ) '  S            Select face subset.'
  write ( *, '(a)' ) '  T            Transform data.'
  write ( *, '(a)' ) '  U            Renumber faces and analyze.'
  write ( *, '(a)' ) '  V            Convert polygons to triangles.'
  write ( *, '(a)' ) '  W            Reverse faces and normals.'

  return
end
subroutine hrc_read ( bad_num, cor3, cor3_max, cor3_num, dup_num, face, &
  face_material, face_max, face_num, face_order, filein_name, &
  ierror, iunit, line_dex, line_material, line_max, line_num, &
  material_max, material_name, material_num, material_rgba, &
  order_max, text_num, texture_max, texture_name, texture_num, &
  vertex_material, vertex_normal, vertex_tex_uv )

!*****************************************************************************80
!
!! HRC_READ reads graphics information from a SoftImage HRC file.
!
!  Discussion:
!
!    It is intended that the information read from the file can
!    either start a whole new graphics object, or simply be added
!    to a current graphics object via the '<<' command.  
!
!    This is controlled by whether the input values have been zeroed 
!    out or not.  This routine simply tacks on the information it 
!    finds to the current graphics object.
!
!  Example:
!
!    HRCH: Softimage 4D Creative Environment v3.00
!
!
!    model
!    {
!      name         "cube_10x10"
!      scaling      1.000 1.000 1.000
!      rotation     0.000 0.000 0.000
!      translation  0.000 0.000 0.000
!
!      mesh
!      {
!        flag    ( PROCESS )
!        discontinuity  60.000
!
!        vertices   8
!        {
!          [0] position  -5.000  -5.000  -5.000
!          [1] position  -5.000  -5.000  5.000
!          [2] position  -5.000  5.000  -5.000
!          [3] position  -5.000  5.000  5.000
!          [4] position  5.000  -5.000  -5.000
!          [5] position  5.000  -5.000  5.000
!          [6] position  5.000  5.000  -5.000
!          [7] position  5.000  5.000  5.000
!        }
!
!        polygons   6
!        {
!          [0] nodes  4
!              {
!                [0] vertex  0
!                    normal  -1.000  0.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [1] vertex  1
!                    normal  -1.000  0.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [2] vertex  3
!                    normal  -1.000  0.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [3] vertex  2
!                    normal  -1.000  0.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!              }
!              material  0
!          [1] nodes  4
!             {
!                [0] vertex  1
!                    normal  0.000  0.000  1.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [1] vertex  5
!
!    ...etc.....
!
!          [5] nodes  4
!              {
!                [0] vertex  2
!                    normal  0.000  1.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [1] vertex  3
!                    normal  0.000  1.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [2] vertex  7
!                    normal  0.000  1.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [3] vertex  6
!                    normal  0.000  1.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!              }
!              material  0
!        }
!
!        edges   12
!        {
!          [1] vertices  3  2
!          [2] vertices  2  0
!          [3] vertices  0  1
!          [4] vertices  1  3
!          [5] vertices  7  3
!          [6] vertices  1  5
!          [7] vertices  5  7
!          [8] vertices  6  7
!          [9] vertices  5  4
!          [10] vertices  4  6
!          [11] vertices  2  6
!          [12] vertices  4  0
!        }
!      }
!
!      material [0]
!      {
!      name           "kazoo"
!      type           PHONG
!      ambient        0.0  1.0  0.0E+00
!      diffuse        1.0  0.0  0.0E+00
!      specular       0.0  0.0  1.0E+00
!      exponent      50.0E+00
!      reflectivity   0.0E+00
!      transparency   0.0E+00
!      refracIndex    1.0E+00
!      glow           0
!      coc            0.0E+00
!      }
!
!      texture [0]
!      {
!      name          "/usr/users/foss/HOUSE/PICTURES/mellon"
!      glbname       "t2d1"
!      anim          STATIC
!      method        XY
!      repeat        1  1
!      scaling       1.000  1.000
!      offset        0.000  0.000
!      pixelInterp
!      effect        INTENSITY
!      blending      1.000
!      ambient       0.977
!      diffuse       1.000
!      specular      0.966
!      reflect       0.000
!      transp        0.000
!      roughness     0.000
!      reflMap       1.000
!      rotation      0.000
!      txtsup_rot    0.000  0.000  0.000
!      txtsup_trans  0.000  0.000  0.000
!      txtsup_scal   1.000  1.000  1.000
!      }
!
!    }
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BAD_NUM, the number of bad text lines.
!
!    Input/output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input/output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Output, integer ( kind = 4 ) DUP_NUM, the number of duplicate points that were dropped.
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input/output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices 
!    per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which data is read.
!
!    Input/output, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated 
!    by -1.
!
!    Input/output, integer ( kind = 4 ) LINE_MATERIAL(LINE_MAX), material index for line.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input/output, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input/output, character ( len = * ) MATERIAL_NAME(MATERIAL_MAX), 
!    material names.
!
!    Input/output, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Input/output, real ( kind = 4 ) MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B
!    and A values.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Output, integer ( kind = 4 ) TEXT_NUM, the number of text lines.
!
!    Input, integer ( kind = 4 ) TEXTURE_MAX, the maximum number of textures.
!
!    Input/output, character ( len = * ) TEXTURE_NAME(TEXTURE_MAX), 
!    texture names.
!
!    Input/output, integer ( kind = 4 ) TEXTURE_NUM, the number of textures.
!
!    Input/output, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input/output, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals 
!    at vertices.
!
!    Input/output, real ( kind = 4 ) VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex 
!    texture coordinates.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ), parameter :: level_max = 10
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) texture_max

  integer ( kind = 4 ) bad_num
  character ( len = 4 ) char4
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) cor3_num_old
  logical, parameter :: debug = .FALSE.
  logical done
  integer ( kind = 4 ) dup_num
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ival1
  integer ( kind = 4 ) ival2
  integer ( kind = 4 ) ival3
  integer ( kind = 4 ) ival4
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) iword
  integer ( kind = 4 ) jval
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) level
  character ( len = 256 ) level_name(0:level_max)
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_num
  logical lval
  character ( len = * ) material_name(material_max)
  integer ( kind = 4 ) material_num
  integer ( kind = 4 ) material_num_old
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ) nlbrack
  integer ( kind = 4 ) nrbrack
  integer ( kind = 4 ), parameter :: OFFSET = 1
  real ( kind = 4 ) rval
  logical s_eqi
  logical s_is_i4
  real ( kind = 4 ) temp(3)
  character ( len = 256 ) text
  integer ( kind = 4 ) text_num
  character ( len = * ) texture_name(texture_max)
  integer ( kind = 4 ) texture_num
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  real ( kind = 4 ) vertex_tex_uv(2,order_max,face_max)
  character ( len = 256 ) word
  character ( len = 256 ) word2
  character ( len = 256 ) wordm1
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  real ( kind = 4 ) z

  ierror = 0
  ival1 = 0
  ival2 = 0
  ival3 = 0
  ival4 = 0
  jval = 0
  level = 0
  level_name(0) = 'Top'
  nlbrack = 0
  nrbrack = 0
  cor3_num_old = cor3_num
  material_num_old = material_num
  word = ' '
  wordm1 = ' '
!
!  Read a line of text from the file.
!
10    continue

  read ( iunit, '(a)', iostat = ios ) text

  if ( ios /= 0 ) then
    go to 50
  end if

  text_num = text_num + 1

  if ( text == ' ' ) then
    go to 10
  end if
!
!  The first line of the file must be the header.
!
  if ( text_num == 1 ) then

    if ( .not. s_eqi ( text(1:5), 'HRCH:' ) ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HRC_READ - Fatal error!'
      write ( *, '(a)' ) '  The input file has a bad header.'
      write ( *, '(a)' ) trim ( text )
      return
    else
      go to 10
    end if

  end if

  done = .true.
  iword = 0
!
!  Save the previous word read.  
!  It helps when a word depends on its cotext_num.
!
20    continue

  if ( word /= ' ' .and. word .ne. ',' ) then
    wordm1 = word
  end if
!
!  Read a word from the line.
!
  call word_next_read ( text, word, done )
!
!  If no more words in this line, read in a whole new line.
!
  if ( done ) then
    go to 10
  end if
!
!  Ignore blanks and commas.
!
  if ( word == ' ' .or. word == ',' ) then
    go to 20
  end if
!
!  Count the words in the current line, and the total.
!
  iword = iword + 1
!
!  If the word is a curly bracket, count it.
!
  if ( word == '{' ) then

    nlbrack = nlbrack + 1

  else if ( word .eq. '}' ) then

    nrbrack = nrbrack + 1

    if ( nlbrack < nrbrack ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HRC_READ - Fatal error!'
      write ( *, '(a,i6)' ) '  Extraneous right bracket, line ', text_num
      write ( *, '(a)' ) trim ( text )
      write ( *, '(a)' ) 'Currently processing field:'
      write ( *, '(a)' ) trim ( level_name(level) )
      ierror = 1
      return
    end if

  end if
!
!  If the word is a left bracket, the previous word is the name of a node.
!
  if ( word == '{' ) then

    level = nlbrack - nrbrack
    if ( level < 0 ) then
      write ( *, '(a)' ) 'Too many right brackets!'
      level = 0
    else if ( level_max < level ) then
      write ( *, '(a)' ) 'Too many left brackets!'
      level = level_max
    end if
    level_name(level) = wordm1

    if ( debug ) then
      write ( *, '(a)' ) ' '
      do i = 0, level
        write ( *, '(i3,2x,a)' ) i, trim ( level_name(i) )
      end do
    end if

  end if
!
!  CONTROLPOINTS
!
  if ( s_eqi ( level_name(level), 'CONTROLPOINTS' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      line_num = line_num + 1
      if ( line_num <= line_max ) then
        line_dex(line_num) = -1 + OFFSET
        line_material(line_num) = material_num
      end if

      level = nlbrack - nrbrack

    else if ( word(1:1) == '[' ) then

      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'POSITION' ) ) then

      do i = 1, 3
        call word_next_read ( text, word2, done )
        call s_to_r4 ( word2, rval, ierror, lchar )
        temp(i) = rval
      end do
!
!  If the coordinate values already exist in COR3, then  
!  save space by using the index of a previous copy.
!
      if ( cor3_num <= 1000 ) then
        call r4col_find ( 3, 3, cor3_num, cor3, temp, icor3 )
      else
        icor3 = 0
      end if

      if ( icor3 == 0 ) then

        cor3_num = cor3_num + 1
        if ( cor3_num <= cor3_max ) then
          cor3(1,cor3_num) = temp(1)
          cor3(2,cor3_num) = temp(2)
          cor3(3,cor3_num) = temp(3)
        end if

        icor3 = cor3_num 

      else

        dup_num = dup_num + 1

      end if

      line_num = line_num + 1
      if ( line_num <= line_max ) then
        line_dex(line_num) = icor3
        line_material(line_num) = material_num
      end if

    else

      go to 99

    end if
!
!  EDGES
!
  else if ( s_eqi ( level_name(level), 'EDGES' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word(1:1) == '[' ) then

      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'VERTICES' ) ) then

      call word_next_read ( text, word2, done )
      lval = s_is_i4 ( word2, jval )
      line_num = line_num + 1

      if ( line_num <= line_max ) then
        line_dex(line_num) = jval + cor3_num_old + OFFSET
        line_material(line_num) = material_num
      end if

      call word_next_read ( text, word2, done )
      lval = s_is_i4 ( word2, jval )
      line_num = line_num + 1

      if ( line_num <= line_max ) then
        line_dex(line_num) = jval + cor3_num_old + OFFSET
        line_material(line_num) = material_num
      end if

      line_num = line_num + 1

      if ( line_num <= line_max ) then
        line_dex(line_num) = -1 + OFFSET
        line_material(line_num) = -1 + OFFSET
      end if

    else

      go to 99

    end if
!
!  MATERIAL
!
  else if ( s_eqi ( level_name(level), 'MATERIAL' ) ) then

    if ( word == '{' ) then

      material_num = material_num + 1
      if ( material_num <= material_max ) then
        call i4_to_s_zero ( material_num, char4 )
        material_name(material_num) = 'Material_' // char4
      end if

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word(1:1) == '[' ) then

    else if ( s_eqi ( word, 'AMBIENT' ) ) then

      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'COC' ) ) then

      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'DIFFUSE' ) ) then

      call word_next_read ( text, word2, done )
      call s_to_r4 ( word2, rval, ierror, lchar )
      material_rgba(1,material_num) = rval

      call word_next_read ( text, word2, done )
      call s_to_r4 ( word2, rval, ierror, lchar )
      material_rgba(2,material_num) = rval

      call word_next_read ( text, word2, done )
      call s_to_r4 ( word2, rval, ierror, lchar )
      material_rgba(3,material_num) = rval

    else if ( s_eqi ( word, 'EXPONENT' ) ) then

      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'GLOW' ) ) then

      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'NAME' ) ) then

      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )
      material_name(material_num) = word2
      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'REFLECTIVITY' ) ) then

      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'REFRACINDEX' ) ) then

      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'SPECULAR' ) ) then

      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'TRANSPARENCY' ) ) then

      call word_next_read ( text, word2, done )
      call s_to_r4 ( word2, rval, ierror, lchar )
      material_rgba(4,material_num) = 1.0E+00 - rval

    else if ( s_eqi ( word, 'TYPE' ) ) then

      call word_next_read ( text, word2, done )

    else

      go to 99

    end if
!
!  MESH
!
  else if ( s_eqi ( level_name(level), 'MESH' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then
      level = nlbrack - nrbrack
    else if ( s_eqi ( word, 'DISCONTINUITY' ) ) then
      go to 10
    else if ( s_eqi ( word, 'EDGES' ) ) then
      go to 10
    else if ( s_eqi ( word, 'FLAG' ) ) then
      go to 10
    else if ( s_eqi ( word, 'POLYGONS' ) ) then
      go to 10
    else if ( s_eqi ( word, 'VERTICES' ) ) then
      go to 10
    else
      go to 99
    end if
!
!  MODEL
!
  else if ( s_eqi ( level_name(level), 'MODEL' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'MATERIAL' ) ) then
      go to 10
    else if ( s_eqi ( word, 'MESH' ) ) then
      go to 10
    else if ( s_eqi ( word, 'NAME' ) ) then
      go to 10
    else if ( s_eqi ( word, 'PATCH' ) ) then
      go to 10
    else if ( s_eqi ( word, 'ROTATION' ) ) then
      go to 10
    else if ( s_eqi ( word, 'SCALING' ) ) then
      go to 10
    else if ( s_eqi ( word, 'SPLINE' ) ) then
      go to 10
    else if ( s_eqi ( word, 'TRANSLATION' ) ) then
      go to 10
    else
      go to 99
    end if
!
!  NODES
!
  else if ( s_eqi ( level_name(level), 'NODES' ) ) then

    if ( word == '{' ) then

      face_num = face_num + 1
      ivert = 0
      face_order(face_num) = 0

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word(1:1) == '[' ) then

      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'NORMAL' ) ) then

      call word_next_read ( text, word2, done )
      call s_is_r4 ( word2, x, lval )

      call word_next_read ( text, word2, done )
      call s_is_r4 ( word2, y, lval )

      call word_next_read ( text, word2, done )
      call s_is_r4 ( word2, z, lval )

      if ( face_num <= face_max .and. ivert <= order_max ) then
        vertex_normal(1,ivert,face_num) = x
        vertex_normal(2,ivert,face_num) = y
        vertex_normal(3,ivert,face_num) = z
      end if

    else if ( s_eqi ( word, 'UVTEXTURE' ) ) then

      call word_next_read ( text, word2, done )
      call s_is_r4 ( word2, x, lval )

      call word_next_read ( text, word2, done )
      call s_is_r4 ( word2, y, lval )

      if ( face_num <= face_max .and. ivert <= order_max ) then
        vertex_tex_uv(1,ivert,face_num) = x
        vertex_tex_uv(2,ivert,face_num) = y
      end if

    else if ( s_eqi ( word, 'VERTEX' ) ) then

      call word_next_read ( text, word2, done )
      lval = s_is_i4 ( word2, jval )
      ivert = ivert + 1

      if ( ivert <= order_max .and. face_num <= face_max ) then
        face_order(face_num) = face_order(face_num) + 1
        face(ivert,face_num) = jval + cor3_num_old + OFFSET
      end if
!
!  What do I do with this?  Define a vertex material?
!
    else if ( s_eqi ( word, 'VERTEXCOLOR' ) ) then

      call word_next_read ( text, word2, done )
      lval = s_is_i4 ( word2, ival1 )

      call word_next_read ( text, word2, done )
      lval = s_is_i4 ( word2, ival2 )

      call word_next_read ( text, word2, done )
      lval = s_is_i4 ( word2, ival3 )

      call word_next_read ( text, word2, done )
      lval = s_is_i4 ( word2, ival4 )

    else
      go to 99
    end if
!
!  PATCH
!
!  JVB: I don't know what to do with this yet.
!
  else if ( s_eqi ( level_name(level), 'PATCH' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then
      level = nlbrack - nrbrack
    else if ( s_eqi ( word, 'APPROX_TYPE' ) ) then
      call word_next_read ( text, word2, done )
    else if ( s_eqi ( word, 'CONTROLPOINTS' ) ) then
      go to 10
    else if ( s_eqi ( word, 'CURV_U' ) ) then
      call word_next_read ( text, word2, done )
    else if ( s_eqi ( word, 'CURV_V' ) ) then
      call word_next_read ( text, word2, done )
    else if ( s_eqi ( word, 'RECMIN' ) ) then
      call word_next_read ( text, word2, done )
    else if ( s_eqi ( word, 'RECMAX' ) ) then
      call word_next_read ( text, word2, done )
    else if ( s_eqi ( word, 'RECURSION' ) ) then
      call word_next_read ( text, word2, done )
    else if ( s_eqi ( word, 'SPACIAL' ) ) then
      call word_next_read ( text, word2, done )
    else if ( s_eqi ( word, 'TAGGEDPOINTS' ) ) then
      go to 10
    else if ( s_eqi ( word, 'UCURVE' ) ) then
      call word_next_read ( text, word2, done )
    else if ( s_eqi ( word, 'UPOINT' ) ) then
      call word_next_read ( text, word2, done )
    else if ( s_eqi ( word, 'USTEP' ) ) then
      call word_next_read ( text, word2, done )
    else if ( s_eqi ( word, 'UTENSION' ) ) then
      call word_next_read ( text, word2, done )
    else if ( s_eqi ( word, 'UTYPE' ) ) then
      call word_next_read ( text, word2, done )
    else if ( s_eqi ( word, 'VCLOSE' ) ) then

    else if ( s_eqi ( word, 'VCURVE' ) ) then
      call word_next_read ( text, word2, done )
    else if ( s_eqi ( word, 'VIEWDEP' ) ) then
      call word_next_read ( text, word2, done )
    else if ( s_eqi ( word, 'VPOINT' ) ) then
      call word_next_read ( text, word2, done )
    else if ( s_eqi ( word, 'VSTEP' ) ) then
      call word_next_read ( text, word2, done )
    else if ( s_eqi ( word, 'VTENSION' ) ) then
      call word_next_read ( text, word2, done )
    else if ( s_eqi ( word, 'VTYPE' ) ) then
      call word_next_read ( text, word2, done )
    else
      go to 99
    end if
!
!  POLYGONS
!
  else if ( s_eqi ( level_name(level), 'POLYGONS' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word == '[' ) then

      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'MATERIAL' ) ) then

      call word_next_read ( text, word2, done )
      lval = s_is_i4 ( word2, jval )

      face_material(face_num) = jval + material_num_old + OFFSET

      do i = 1, order_max
        vertex_material(i,face_num) = jval + material_num_old + OFFSET
      end do

    else if ( s_eqi ( word, 'NODES' ) ) then
      call word_next_read ( text, word2, done )
    else
      go to 99
    end if
!
!  SPLINE
!
  else if ( s_eqi ( level_name(level), 'SPLINE' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'CONTROLPOINTS' ) ) then
      go to 10
    else if ( s_eqi ( word, 'NAME' ) ) then
      go to 10
    else if ( s_eqi ( word, 'NBKEYS' ) ) then
      go to 10
    else if ( s_eqi ( word, 'STEP' ) ) then
      go to 10
    else if ( s_eqi ( word, 'TENSION' ) ) then
      go to 10
    else if ( s_eqi ( word, 'TYPE' ) ) then
      go to 10
    else
      go to 99
    end if
!
!  TAGGEDPOINTS
!
  else if ( s_eqi ( level_name(level), 'TAGGEDPOINTS' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word(1:1) == '[' ) then

      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'TAGGED' ) ) then
      call word_next_read ( text, word2, done )
    else
      go to 99
    end if
!
!  TEXTURE
!
  else if ( s_eqi ( level_name(level), 'TEXTURE' ) ) then

    if ( word == '{' ) then

      texture_num = texture_num + 1

      if ( texture_num <= texture_max ) then
        call i4_to_s_zero ( texture_num, char4 )
        texture_name(texture_num) = 'Texture_' // char4
      end if

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word(1:1) == '[' ) then

    else if ( s_eqi ( word, 'AMBIENT' ) ) then

      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'ANIM' ) ) then

      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'BLENDING' ) ) then

      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'DIFFUSE' ) ) then

      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'EFFECT' ) ) then

      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'GLBNAME' ) ) then

      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'METHOD' ) ) then

      call word_next_read ( text, word2, done )
!
!  (I assume there are initial and trailing quotes in the NAME field,
!  which are treated as separate words.)
!
    else if ( s_eqi ( word, 'NAME' ) ) then

      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )
      texture_name(texture_num) = word2
      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'OFFSET' ) ) then

      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'PIXELINTERP' ) ) then

    else if ( s_eqi ( word, 'REFLECT' ) ) then

      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'REFLMAP' ) ) then

      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'REPEAT' ) ) then

      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'ROTATION' ) ) then

      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'ROUGHNESS' ) ) then

      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'SCALING' ) ) then

      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'SPECULAR' ) ) then

      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'TRANSP' ) ) then

      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'TXTSUP_ROT' ) ) then

      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'TXTSUP_SCAL' ) ) then

      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'TXTSUP_TRANS' ) ) then

      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )

    else

      go to 99

    end if
!
!  VERTICES
!
  else if ( s_eqi ( level_name(level), 'VERTICES' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word(1:1) == '[' ) then

      call word_next_read ( text, word2, done )
      call word_next_read ( text, word2, done )

    else if ( s_eqi ( word, 'POSITION' ) ) then

      cor3_num = cor3_num + 1

      call word_next_read ( text, word2, done )
      call s_is_r4 ( word2, x, lval )

      call word_next_read ( text, word2, done )
      call s_is_r4 ( word2, y, lval )

      call word_next_read ( text, word2, done )
      call s_is_r4 ( word2, z, lval )

      if ( cor3_num <= cor3_max ) then
        cor3(1,cor3_num) = x
        cor3(2,cor3_num) = y
        cor3(3,cor3_num) = z
      end if

    else
      go to 99
    end if
!
!  Any other word:
!
  else

  end if

  go to 20
!
!  Bad data
!
99    continue

  bad_num = bad_num + 1

  if ( bad_num <= 10 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HRC_READ - Warning!'
    write ( *, '(a)' ) '  Bad data on level ' // trim ( level_name(level) )
    write ( *, '(a)' ) '  Bad word: ' // trim ( word )
    write ( *, '(a,i6)' ) '  Line number: ', text_num
    write ( *, '(a)' ) trim ( text )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HRC_READ - Fatal error!'
    write ( *, '(a)' ) '  Too many warnings!'
    return
  end if

  go to 10
!
!  Normal end of information in file.
!
50    continue
!
!  Check the "materials" defining a line.
!
!  If COORDINDEX is -1, so should be the MATERIALINDEX.
!  If COORDINDEX is not -1, then the MATERIALINDEX shouldn't be either.
!
  do i = 1, line_num

    if ( line_dex(i) == -1 + OFFSET ) then
      line_material(i) = -1 + OFFSET
    else if ( line_material(i) == -1 + OFFSET ) then
      line_material(i) = material_num
    end if

  end do
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'HRC_READ - Read ', text_num, ' text lines from ' &
    // trim ( filein_name )

  return
end
subroutine hrc_write ( cor3, cor3_max, cor3_num, face, face_material, &
  face_max, face_num, face_order, fileout_name, iunit, line_dex, &
  line_max, line_num, material_max, material_name, material_num, &
  material_rgba, order_max, texture_max, texture_name, texture_num, &
  vertex_normal, vertex_tex_uv )

!*****************************************************************************80
!
!! HRC_WRITE writes graphics data to an HRC SoftImage file.
!
!  Example:
!
!    HRCH: Softimage 4D Creative Environment v3.00
!
!
!    model
!    {
!      name         "cube_10x10"
!      scaling      1.000 1.000 1.000
!      rotation     0.000 0.000 0.000
!      translation  0.000 0.000 0.000
!
!      mesh
!      {
!        flag    ( PROCESS )
!        discontinuity  60.000
!
!        vertices   8
!        {
!          [0] position  -5.000  -5.000  -5.000
!          [1] position  -5.000  -5.000  5.000
!          [2] position  -5.000  5.000  -5.000
!          [3] position  -5.000  5.000  5.000
!          [4] position  5.000  -5.000  -5.000
!          [5] position  5.000  -5.000  5.000
!          [6] position  5.000  5.000  -5.000
!          [7] position  5.000  5.000  5.000
!        }
!
!        polygons   6
!        {
!          [0] nodes  4
!              {
!                [0] vertex  0
!                    normal  -1.000  0.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [1] vertex  1
!                    normal  -1.000  0.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [2] vertex  3
!                    normal  -1.000  0.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [3] vertex  2
!                    normal  -1.000  0.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!              }
!              material  0
!          [1] nodes  4
!             {
!                [0] vertex  1
!                    normal  0.000  0.000  1.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [1] vertex  5
!
!    ...etc.....
!
!          [5] nodes  4
!              {
!                [0] vertex  2
!                    normal  0.000  1.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [1] vertex  3
!                    normal  0.000  1.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [2] vertex  7
!                    normal  0.000  1.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!                [3] vertex  6
!                    normal  0.000  1.000  0.000
!                    uvTexture  0.000  0.000
!                    vertexColor 255 178 178 178
!              }
!              material  0
!        }
!
!        edges   12
!        {
!          [1] vertices  3  2
!          [2] vertices  2  0
!          [3] vertices  0  1
!          [4] vertices  1  3
!          [5] vertices  7  3
!          [6] vertices  1  5
!          [7] vertices  5  7
!          [8] vertices  6  7
!          [9] vertices  5  4
!          [10] vertices  4  6
!          [11] vertices  2  6
!          [12] vertices  4  0
!        }
!      }
!
!      material [0]
!      {
!      name           "kazoo"
!      type           PHONG
!      ambient        0.0  1.0  0.0E+00
!      diffuse        1.0  0.0  0.0E+00
!      specular       0.0  0.0  1.0E+00
!      exponent      50.0E+00
!      reflectivity   0.0E+00
!      transparency   0.0E+00
!      refracIndex    1.0E+00
!      glow           0
!      coc            0.0E+00
!      }
!
!      texture [0]
!      {
!      name          "/usr/users/foss/HOUSE/PICTURES/mellon"
!      glbname       "t2d1"
!      anim          STATIC
!      method        XY
!      repeat        1  1
!      scaling       1.000  1.000
!      offset        0.000  0.000
!      pixelInterp
!      effect        INTENSITY
!      blending      1.000
!      ambient       0.977
!      diffuse       1.000
!      specular      0.966
!      reflect       0.000
!      transp        0.000
!      roughness     0.000
!      reflMap       1.000
!      rotation      0.000
!      txtsup_rot    0.000  0.000  0.000
!      txtsup_trans  0.000  0.000  0.000
!      txtsup_scal   1.000  1.000  1.000
!      }
!
!    }
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), face materials.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, character ( len = * ) MATERIAL_NAME(MATERIAL_MAX), material names.
!
!    Input, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, character ( len = * ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
!    Input, integer ( kind = 4 ) TEXTURE_NUM, the number of textures.
!
!    Input, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
!    Input, real ( kind = 4 ) VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture
!    coordinates.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) texture_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) fileout_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) jrel
  integer ( kind = 4 ) k
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_num
  character ( len = * ) material_name(material_max)
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ) npts
  integer ( kind = 4 ) nseg
  integer ( kind = 4 ), parameter :: OFFSET = 1
  character ( len = 100 ) text
  integer ( kind = 4 ) text_num
  character ( len = * ) texture_name(texture_max)
  integer ( kind = 4 ) texture_num
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  real ( kind = 4 ) vertex_tex_uv(2,order_max,face_max)
  character ( len = 10 ) word
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  real ( kind = 4 ) z

  nseg = 0
  text_num = 0

  write ( iunit, '(a)' ) 'HRCH: Softimage 4D Creative Environment v3.00'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) ' '
  text_num = text_num + 3

  write ( iunit, '(a)' ) 'model'
  write ( iunit, '(a)' ) '{'
  write ( iunit, '(a)' ) '  name         "' // trim ( fileout_name ) // '"'
  write ( iunit, '(a)' ) '  scaling      1.000 1.000 1.000'
  write ( iunit, '(a)' ) '  rotation     0.000 0.000 0.000'
  write ( iunit, '(a)' ) '  translation  0.000 0.000 0.000'
  text_num = text_num + 6

  if ( 0 < face_num ) then

    write ( iunit, '(a)' ) ' '
    write ( iunit, '(a)' ) '  mesh'
    write ( iunit, '(a)' ) '  {'
    write ( iunit, '(a)' ) '    flag    ( PROCESS )'
    write ( iunit, '(a)' ) '    discontinuity  60.000'
    text_num = text_num + 5
!
!  Point coordinates.
!
    if ( 0 < cor3_num ) then

      write ( iunit, '(a)' ) ' '
      write ( text, '(a, i8)' ) 'vertices ', cor3_num
      call s_blanks_delete ( text )
      write ( iunit, '(4x,a)' ) trim ( text )
      write ( iunit, '(a)' )     '    {'
      text_num = text_num + 3
 
      do j = 1, cor3_num

        write ( word, '( ''['', i8, '']'' )' ) j-OFFSET
        call s_blank_delete ( word )

        write ( text, '(a,'' position '',3f12.3)' ) trim ( word ), cor3(1:3,j)
        call s_blanks_delete ( text )
        write ( iunit, '(6x,a)' ) trim ( text )
        text_num = text_num + 1
      end do

      write ( iunit, '(a)' )     '    }'
      text_num = text_num + 1

    end if
!
!  Faces.
!
    write ( iunit, '(a)' ) ' '
    write ( text, '(a,i8)' ) 'polygons ', face_num
    call s_blanks_delete ( text )
    write ( iunit, '(4x,a)' ) trim ( text )
    write ( iunit, '(a)' ) '    {'
    text_num = text_num + 3

    do iface = 1, face_num

      write ( word, '( ''['', i8, '']'' )' ) iface-OFFSET
      call s_blank_delete ( word )
      write ( text, '(a,'' nodes '',i8 )' ) trim ( word ), face_order(iface)
      call s_blanks_delete ( text )
      write ( iunit, '(6x,a)' ) trim ( text )
      write ( iunit, '(a)' ) '      {'
      text_num = text_num + 2

      do ivert = 1, face_order(iface)

        write ( word, '( ''['', i8, '']'' )' ) ivert-OFFSET
        call s_blank_delete ( word )
        write ( text, '( a,'' vertex '',i8 )' ) trim ( word ), &
          face(ivert,iface) - OFFSET
        call s_blanks_delete ( text )
        write ( iunit, '(8x,a)' ) trim ( text )

        x = vertex_normal(1,ivert,iface)
        y = vertex_normal(2,ivert,iface)
        z = vertex_normal(3,ivert,iface)
        write ( text, '(a,3f12.3)' ) 'normal ', x, y, z
        call s_blanks_delete ( text )
        write ( iunit, '(12x,a)' ) trim ( text )

        x = vertex_tex_uv(1,ivert,iface)
        y = vertex_tex_uv(2,ivert,iface)
        write ( text, '(a,2f12.3)' ) 'uvTexture ', x, y
        call s_blanks_delete ( text )
        write ( iunit, '(12x,a)' ) trim ( text )

        text = 'vertexColor 255 178 178 178'
        call s_blanks_delete ( text )
        write ( iunit, '(12x,a)' ) trim ( text )

        text_num = text_num + 4

      end do

      write ( iunit, '(a)' ) '      }'
      write ( text, '(''material '',i8)' ) face_material(iface) - OFFSET
      call s_blanks_delete ( text )
      write ( iunit, '(6x,a)' ) trim ( text )
      text_num = text_num + 2

    end do
 
    write ( iunit, '(a)' ) '    }'
    write ( iunit, '(a)' ) '  }'
 
    text_num = text_num + 2

  end if
!
!  IndexedLineSet.
!
  if ( 0 < line_num ) then

    nseg = 0

    jhi = 0

10  continue

    jlo = jhi
!
!  Look for the next index JLO that is not -1.
!
    do

      jlo = jlo + 1

      if ( line_num < jlo ) then
        go to 20
      end if

      if ( line_dex(jlo) /= -1+OFFSET ) then
        exit
      end if

    end do
!
!  Look for the highest following index JHI that is not -1.
!
    jhi = jlo + 1

    if ( line_num < jhi ) then
      go to 20
    end if

    if ( line_dex(jhi) == -1+OFFSET ) then
      go to 10
    end if

    do while ( jhi < line_num )
   
      if ( line_dex(jhi+1) == -1+OFFSET ) then
        exit
      end if

      jhi = jhi + 1

    end do
!
!  Our next line segment involves LINE_DEX indices JLO through JHI.
!
    nseg = nseg + 1
    write ( text, '(''spl'', i8 )' ) nseg
    call s_blank_delete ( text )

    npts = jhi + 1 - jlo

    write ( iunit, '(a)' ) ' '
    write ( iunit, '(a)' ) '  spline'
    write ( iunit, '(a)' ) '  {'
    write ( iunit, '(a)' ) '    name     "' // trim ( text ) // '"'
    write ( iunit, '(a)' ) '    type     LINEAR'
    write ( iunit, '(a,i8)' ) '    nbKeys   ', npts
    write ( iunit, '(a)' ) '    tension  0.000'
    write ( iunit, '(a)' ) '    step     1'
    write ( iunit, '(a)' ) ' '
    text_num = text_num + 9

    write ( iunit, '(a)' ) '    controlpoints'
    write ( iunit, '(a)' ) '    {'
    text_num = text_num + 2

    do j = jlo, jhi
 
      jrel = j - jlo
      k = line_dex(j)
      write ( word, '( ''['', i8, '']'')' ) jrel
      call s_blank_delete ( word )
      write ( text, '( a, '' position '', 3f12.4)' ) &
        trim ( word ), cor3(1,k), cor3(2,k), cor3(3,k)
      call s_blanks_delete ( text )
      write ( iunit, '(6x,a)' ) trim ( text )
      text_num = text_num + 1

    end do

    write ( iunit, '(a)' ) '    }'
    write ( iunit, '(a)' ) '  }'
    text_num = text_num + 2

    go to 10

20  continue

  end if
!
!  MATERIALS
!
  do i = 1, material_num

    write ( text, '(''['', i8, '']'' )' ) i-OFFSET
    call s_blank_delete ( text )
    write ( iunit, '(a)' ) '  material ' // trim ( text )

    write ( iunit, '(a)' ) '  {'
    write ( iunit, '(a)' ) 'name          "' // trim ( material_name(i) ) // '"'
    write ( iunit, '(a)' ) '    type           PHONG'
    write ( iunit, '(a,3f10.4)' ) '    ambient        ', &
      material_rgba(1,i), material_rgba(2,i), material_rgba(3,i)
    write ( iunit, '(a,3f10.4)' ) '    diffuse        ', &
      material_rgba(1,i), material_rgba(2,i), material_rgba(3,i)
    write ( iunit, '(a,3f10.4)' ) '    specular       ', &
      material_rgba(1,i), material_rgba(2,i), material_rgba(3,i)
    write ( iunit, '(a)' ) '    exponent      50.0'
    write ( iunit, '(a)' ) '    reflectivity   0.0'
    write ( iunit, '(a,f10.4)' ) '    transparency   ', 1.0E+00 - material_rgba(4,i)
    write ( iunit, '(a)' ) '    refracIndex    1.0'
    write ( iunit, '(a)' ) '    glow           0'
    write ( iunit, '(a)' ) '    coc            0.0'
    write ( iunit, '(a)' ) '  }'

    text_num = text_num + 14

  end do
!
!  TEXTURES
!
  do i = 1, texture_num

    write ( text, '(''['', i8, '']'' )' ) i-OFFSET
    call s_blank_delete ( text )
    write ( iunit, '(a)' ) '  texture [' // trim ( text ) //']'

    write ( iunit, '(a)' ) '{'
    write ( iunit, '(a)' ) 'name          "' // trim ( texture_name(i) ) // '"'
    write ( iunit, '(a)' ) 'glbname       "t2d1"'
    write ( iunit, '(a)' ) 'anim          STATIC'
    write ( iunit, '(a)' ) 'method        XY'
    write ( iunit, '(a)' ) 'repeat        1  1'
    write ( iunit, '(a)' ) 'scaling       1.000  1.000'
    write ( iunit, '(a)' ) 'offset        0.000  0.000'
    write ( iunit, '(a)' ) 'pixelInterp'
    write ( iunit, '(a)' ) 'effect        INTENSITY'
    write ( iunit, '(a)' ) 'blending      1.000'
    write ( iunit, '(a)' ) 'ambient       0.977'
    write ( iunit, '(a)' ) 'diffuse       1.000'
    write ( iunit, '(a)' ) 'specular      0.966'
    write ( iunit, '(a)' ) 'reflect       0.000'
    write ( iunit, '(a)' ) 'transp        0.000'
    write ( iunit, '(a)' ) 'roughness     0.000'
    write ( iunit, '(a)' ) 'reflMap       1.000'
    write ( iunit, '(a)' ) 'rotation      0.000'
    write ( iunit, '(a)' ) 'txtsup_rot    0.000  0.000  0.000'
    write ( iunit, '(a)' ) 'txtsup_trans  0.000  0.000  0.000'
    write ( iunit, '(a)' ) 'txtsup_scal   1.000  1.000  1.000'
    write ( iunit, '(a)' ) '}'
    write ( iunit, '(1x)' ) 

    text_num = text_num + 25

  end do

  write ( iunit, '(a)' ) '}'
  text_num = text_num + 1
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'HRC_WRITE - Wrote ', text_num, ' text lines to ' // trim ( fileout_name )


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
!        I         J     MOD  I4_MODP   I4_MODP Factorization
!
!      107        50       7       7    107 =  2 *  50 + 7
!      107       -50       7       7    107 = -2 * -50 + 7
!     -107        50      -7      43   -107 = -3 *  50 + 43
!     -107       -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 199
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
    write ( *, '(a,i6)' ) '  I4_MODP ( I, J ) called with J = ', j
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
!! I4_SWAP switches two integer values.
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
subroutine i4_to_s_zero ( intval, s )

!*****************************************************************************80
!
!! I4_TO_S_ZERO converts an integer to a string, with zero padding.
!
!  Example:
!
!    Assume that S is 6 characters long:
!
!    INTVAL  S
!
!         1  000001
!        -1  -00001
!         0  000000
!      1952  001952
!    123456  123456
!   1234567  ******  <-- Not enough room!
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
!    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be right justified, and zero padded.
!    If there is not enough space, the string will be filled with stars.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) ival
  character ( len = * ) s

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  ival = intval
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = - ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  Working from right to left, strip off the digits of the integer
!  and place them into S(ILO:IHI).
!
  ipos = ihi

  do while ( ival /= 0 .or. ipos == ihi )

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do i = 1, ihi
        s(i:i) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )

    s(ipos:ipos) = c
    ipos = ipos - 1

  end do
!
!  Fill the empties with zeroes.
!
  do i = ilo, ipos
    s(i:i) = '0'
  end do

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
subroutine infile ( filein_name, ierror, filein_type )

!*****************************************************************************80
!
!! INFILE determines the input filename and type.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Output, character ( len = 10 ) FILEIN_TYPE, the type of the file, which is
!    set to the filename extension. 
!
  implicit none

  character ( len = * ) filein_name
  character ( len = 10 ) filein_type
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  logical s_eqi

  ierror = 0

  if ( filein_name == ' ' ) then
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INFILE:'
    write ( *, '(a)' ) '  Enter the name of a graphics file to be read:'
 
    read ( *, '(a)', iostat = ios ) filein_name
 
    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INFILE - Error!'
      write ( *, '(a)' ) '  The input file was not specified correctly.'
      ierror = ios
      return
    end if
 
  end if
!
!  Determine the input file type.
!
  call file_name_ext_get ( filein_name, i1, i2 )

  if ( i1 /= 0 ) then
    filein_type = filein_name(i1:i2)
  else
    filein_type = ' '
  end if
 
  if ( filein_type == ' ' ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INFILE - Warning!'
    write ( *, '(a)' ) '  Could not determine the input file type.'
    write ( *, '(a)' ) '  The input file name is:'
    write ( *, '(a)' ) trim ( filein_name )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The file type should occur after the period.'
    write ( *, '(a)' ) '  Please specify the file type you are using:'
    write ( *, '(a)' ) '    3ds,  ase,  byu, dxf,  hrc, ' // &
                           'iv,   obj,  off, oogl, smf, ' // &
                           'stl,  stla, tri, tria, ts,  ' // &
                           'vla,  wrl,  xyz:'
    read ( *, '(a)', iostat = ios ) filein_type

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INFILE - Fatal error!'
      write ( *, '(a)' ) '  Input failure while reading input file type.'
      stop
    end if

  end if
 
  if ( .not. ( &
    s_eqi ( filein_type, '3S' ) .or. &
    s_eqi ( filein_type, 'ASE' ) .or. &
    s_eqi ( filein_type, 'BYU' ) .or. &
    s_eqi ( filein_type, 'DXF' ) .or. &
    s_eqi ( filein_type, 'HRC' ) .or. &
    s_eqi ( filein_type, 'IV' )  .or. &
    s_eqi ( filein_type, 'OBJ' ) .or. &
    s_eqi ( filein_type, 'OFF' ) .or. &
    s_eqi ( filein_type, 'OOGL' ) .or. &
    s_eqi ( filein_type, 'SMF' ) .or. &
    s_eqi ( filein_type, 'STL' ) .or. &
    s_eqi ( filein_type, 'STLA' )  .or. &
    s_eqi ( filein_type, 'TRI' )  .or. &
    s_eqi ( filein_type, 'TRIA' )  .or. &
    s_eqi ( filein_type, 'TS' ) .or. &
    s_eqi ( filein_type, 'VLA' ) .or. &
    s_eqi ( filein_type, 'WRL' ) .or. &
    s_eqi ( filein_type, 'XYZ' ) ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INFILE - Error!'
    write ( *, '(a)' ) '  The input file type was not acceptable!'
    return
  end if

  return
end
subroutine interact ( byte_swap, cor3, cor3_material, cor3_max, cor3_new, &   
  cor3_normal, cor3_tex_uv, debug, edge, edge_max, face, face_area, &
  face_material, face_max, face_normal, face_object, face_order, face_rank,  &
  face_tex_uv, face_tier, filein_name, fileout_name, ierror, line_dex, &  
  line_material, line_max, line_prune, list, material_max, material_name, &
  material_rgba, order_max, object_name, point, point_max, &
  point_num, texture_max, texture_name, texture_temp, transform_matrix, &
  vertex_material, vertex_normal, vertex_tex_uv )

!*****************************************************************************80
!
!! INTERACT interacts with the user to specify input and output files.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Workspace, real COR3(3,COR3_MAX), the coordinates of points.
!
!    Workspace, integer COR3_MATERIAL(COR3_MAX), the material of each node.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Workspace, real COR3_TEX_UV(2,COR3_MAX), UV texture coordinates for nodes.
!
!    Workspace, integer FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Workspace, integer FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Workspace, real FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Workspace, integer FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  0 = no error, nonzero = error.
!
!    Workspace, integer LINE_DEX(LINE_MAX), nodes forming a line, terminated 
!    by -1.
!
!    Workspace, integer LINE_MATERIAL(LINE_MAX), material index for line.
!
!    Workspace, character ( len = * ) MATERIAL_NAME(MATERIAL_MAX), 
!    material names.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Workspace, character ( len = * ) TEXTURE_NAME(TEXTURE_MAX), 
!    texture names.
!
!    Workspace, integer VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Workspace, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
!    Workspace, real VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture
!    coordinates.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) edge_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) point_max
  integer ( kind = 4 ) texture_max

  logical byte_swap
  character ( len = 100 ) command
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_material(cor3_max)
  real ( kind = 4 ) cor3_new(3,cor3_max)
  real ( kind = 4 ) cor3_normal(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  real ( kind = 4 ) cor3_tex_uv(2,cor3_max)
  logical debug
  integer ( kind = 4 ) edge(4,edge_max)
  integer ( kind = 4 ) face(order_max,face_max)
  real ( kind = 4 ) face_area(face_max)
  integer ( kind = 4 ) face_index
  integer ( kind = 4 ) face_material(face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_object(face_max)
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) face_rank(face_max)
  real ( kind = 4 ) face_tex_uv(2,face_max)
  integer ( kind = 4 ) face_tier(face_max)
  character ( len = * ) filein_name
  character ( len = 10 ) filein_type
  character ( len = * ) fileout_name
  character ( len = 10 ) fileout_type
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_prune
  integer ( kind = 4 ) list(cor3_max)
  character ( len = * ) material_name(material_max)
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) group_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) material_num
  integer ( kind = 4 ) object_num
  integer ( kind = 4 ) point(point_max)
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) texture_num
  character ( len = * ) object_name
  logical s_eqi
  character ( len = * ) texture_name(texture_max)
  real ( kind = 4 ) texture_temp(2,order_max*face_max)
  real ( kind = 4 ) transform_matrix(4,4)
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  real ( kind = 4 ) vertex_tex_uv(2,order_max,face_max)
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  real ( kind = 4 ) z
!
!  Print an introductory message.
!
  call hello ( cor3_max, face_max, line_max, material_max, order_max, &
    point_max, texture_max )
!
!  Get the next user command.
!
  do
 
    ierror = 0

    write ( *, '(a)' ) ' '
    write ( *, '(a)' )'Enter command ( or "Help" )'

    read ( *, '(a)', iostat = ios ) command

    if ( ios /= 0 ) then
      ierror = ios
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INTERACT - Fatal error!'
      write ( *, '(a)' ) '  Unexpected end of input.'
      return
    end if
!
!  << means read a new file, and add it to the current information.
!
    if ( command(1:2) == '<<' ) then
 
      filein_name = command(3:)
      call s_blank_delete ( filein_name )
 
      call infile ( filein_name, ierror, filein_type )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INTERACT - Warning!'
        write ( *, '(a)' ) '  The input file name was unacceptable!'
        cycle
      end if

      call data_read ( cor3, cor3_material, cor3_max, &
        cor3_normal, cor3_num, cor3_tex_uv, debug, face, face_area, &
        face_material, face_max, face_normal, face_num, face_order, &
        face_tex_uv, filein_name, filein_type, group_num, ierror, line_dex, &
        line_material, line_max, line_num, material_max, material_name, &
        material_num, material_rgba, object_num, order_max, &
        point, point_max, point_num, texture_max, texture_name, texture_num, &
        texture_temp, vertex_material, vertex_normal, vertex_tex_uv )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INTERACT - Warning!'
        write ( *, '(a)' ) '  The input data was not read properly.'
        ierror = 0
      end if
!
!  Get the next input file to be examined.
!
    else if ( command(1:1) == '<' ) then
 
      filein_name = command(2:)
      call s_blank_delete ( filein_name )
 
      call infile ( filein_name, ierror, filein_type )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INTERACT - Warning!'
        write ( *, '(a)' ) '  The input file name was unacceptable!'
        cycle
      end if

      call data_init ( cor3, cor3_material, cor3_max, cor3_normal, &
        cor3_num,  cor3_tex_uv, face, face_area, face_material, face_max, &
        face_normal, face_num, face_order, face_tex_uv, group_num, line_dex, &
        line_material, line_max, line_num, material_max, material_name, &
        material_num, material_rgba, object_name, object_num, order_max, &
        point, point_max, point_num, texture_max, texture_name, texture_num, &
        texture_temp, transform_matrix, vertex_material, vertex_normal, &
        vertex_tex_uv )

      call data_read ( cor3, cor3_material, cor3_max, &
        cor3_normal, cor3_num, cor3_tex_uv, debug, face, face_area, &
        face_material, face_max, face_normal, face_num, face_order, &
        face_tex_uv, filein_name, filein_type, group_num, ierror, line_dex, &
        line_material, line_max, line_num, material_max, material_name, &
        material_num, material_rgba, object_num, order_max, &
        point, point_max, point_num, texture_max, texture_name, texture_num, &
        texture_temp, vertex_material, vertex_normal, vertex_tex_uv )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' ' 
        write ( *, '(a)' ) 'INTERACT - Warning!'
        write ( *, '(a)' ) '  The input data was not read properly.'
        ierror = 0
      end if
!
!  The ">" command specifies output.
!
    else if ( command(1:1) == '>' ) then

      fileout_name = command(2:)

      call s_blank_delete ( fileout_name )
!
!  Check the output filename.
!
      call outfile ( filein_name, fileout_name, ierror, fileout_type )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INTERACT - Warning!'
        write ( *, '(a)' ) '  Improper output file name.'
        ierror = 0
        cycle
      end if
!
!  Write the output file.
!
      call data_write ( cor3, cor3_material, cor3_max, cor3_normal, cor3_num, &
        cor3_tex_uv, debug, face, face_material, face_max, face_normal, &
        face_num, face_order, face_tex_uv, filein_name, fileout_name, &
        fileout_type, ierror, line_dex, line_material, line_max, line_num, &
        line_prune, material_name, material_max, material_num, material_rgba, &
        object_name, order_max, point, point_max, point_num, texture_max, &
        texture_name, texture_num, vertex_material, vertex_normal, &
        vertex_tex_uv )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INTERACT - Warning!'
        write ( *, '(a)' ) '  An error occurred during output.'
        ierror = 0
        cycle
      end if
!
!  B: Switch byte swapping option:
!
    else if ( s_eqi ( command(1:1), 'B' ) ) then

      byte_swap = .not. byte_swap

      if ( byte_swap ) then
        write ( *, '(a)' ) 'Byte swapping set to TRUE.'
      else
        write ( *, '(a)' ) 'Byte swapping set to FALSE.'
      end if
!
!  D: Switch debug option:
!
    else if ( s_eqi ( command(1:1), 'D' ) ) then

      debug = .not. debug

      if ( debug ) then
        write ( *, '(a)' ) 'Debug option set to TRUE.'
      else
        write ( *, '(a)' ) 'Debug option set to FALSE.'
      end if
!
!  F: Check a face:
!
    else if ( s_eqi ( command(1:1), 'F' ) ) then

      if ( face_num <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Input a graphical object with faces first!'
        cycle
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'Enter a face between 1 and ', face_num
      read ( *, * ) face_index

      call face_print ( cor3, cor3_max, face, face_index, face_material, &
        face_max, face_normal, face_num, face_order, order_max, &
        vertex_material, vertex_normal )
!
!  HELP:
!
    else if ( s_eqi ( command(1:1), 'H' ) ) then
 
      call help
!
!  INFO:
!
    else if ( s_eqi ( command, 'INFO' ) ) then

      call news
!
!  INVERT:
!  Make an inverted copy of the object to give it thickness.
!
    else if ( s_eqi ( command(1:3), 'INV' ) ) then

      call object_invert ( cor3, cor3_material, cor3_max, cor3_normal, &
        cor3_num, face, face_material, face_max, face_normal, face_num, &
        face_order, material_max, material_name, material_num, material_rgba, &
        order_max, vertex_material, vertex_normal )
!
!  LINE_PRUNE:
!  Set line pruning option.
!
    else if ( s_eqi ( command, 'LINE_PRUNE' ) ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SET THE LINE PRUNING OPTION:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  0 means no pruning;'
      write ( *, '(a)' ) '  nonzero means only generate line (I,J)'
      write ( *, '(a)' ) '    if I < J.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The line pruning option is now ', line_prune
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Enter your line_pruning option:'

      read ( *, * ) line_prune
!
!  LINES:
!  Convert all faces to lines.
!
    else if ( s_eqi ( command, 'LINES' ) ) then
 
      if ( 0 < face_num ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INTERACT - Note:'
        write ( *, '(a)' ) '  Face information will be converted'
        write ( *, '(a)' ) '  to line information.'

        call face_to_line ( debug, face, face_max, face_num, face_order, &
          line_dex, line_material, line_max, line_num, line_prune, &
          order_max, vertex_material )

        if ( line_max < line_num ) then
 
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'INTERACT - Note:'
          write ( *, '(a)' ) '  Some face information was lost.'
          write ( *, '(a,i6)' ) '  The maximum number of lines is ', line_max
          write ( *, '(a,i6)' ) '  but we would need at least ', line_num

          line_num = line_max

        end if

        face_num = 0

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INTERACT - Note:'
        write ( *, '(a)' ) '  There were no faces to convert.'

      end if
!
!  NORMALS: recompute the normal vectors for vertices on faces,
!  and average these to get face normal vectors.
!
    else if ( s_eqi ( command(1:1), 'N' ) ) then

      do iface = 1, face_num
        do ivert = 1, face_order(iface)
          vertex_normal(1:3,ivert,iface) = 0.0E+00
        end do
      end do

      call vertex_normal_set ( cor3, cor3_max, face, face_max, &
        face_num, face_order, order_max, vertex_normal )

      cor3_normal(1:3,1:cor3_num) = 0.0E+00

      call cor3_normal_set ( cor3_max, cor3_normal, face, face_area, &
        face_max, face_num, face_order, order_max, vertex_normal )

      face_normal(1:3,1:face_num) = 0.0E+00

      call face_normal_ave ( face_max, face_normal, face_num, face_order, &
        order_max, vertex_normal )
!
!  OHELL: 
!    Use the node normals.
!    Set the vertex normals equal to the node normals.
!    Set the face normals by averaging vertex normals.
!
    else if ( s_eqi ( command(1:1), 'O' ) ) then
!
!  Recompute any zero vertex normals from vertex positions.
!
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Making sure all vertex normals are defined.'

      call vertex_normal_set ( cor3, cor3_max, face, face_max, &
        face_num, face_order, order_max, vertex_normal )
!
!  Compute node normals by averaging vertex normals.
!
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Averaging vertex normals at each node.'

      cor3_normal(1:3,1:cor3_num) = 0.0E+00

      call cor3_normal_set ( cor3_max, cor3_normal, face, face_area, &
        face_max, face_num, face_order, order_max, vertex_normal )
!
!  Copy node normals into vertex normals.
!
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Replacing vertex normals by average.'

      do iface = 1, face_num
        do ivert = 1, face_order(iface)
          icor3 = face(ivert,iface)
          vertex_normal(1:3,ivert,iface) = cor3_normal(1:3,icor3)
        end do
      end do
!
!  Recompute zero face normals by averaging vertex normals.
!
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Averaging vertex normals to get face normals.'

      call face_normal_ave ( face_max, face_normal, face_num, face_order, &
        order_max, vertex_normal )
!
!  QUIT:
!
    else if ( s_eqi ( command(1:1), 'Q' ) ) then
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INTERACT - End of interaction.'
      ierror = 0
      exit
!
!  REVERSE: Reverse normal vectors.
!
    else if ( s_eqi ( command(1:3), 'REV' ) ) then
 
      cor3_normal(1:3,1:cor3_num) = - cor3_normal(1:3,1:cor3_num)

      face_normal(1:3,1:face_num) = - face_normal(1:3,1:face_num)

      do iface = 1, face_num
        do ivert = 1, face_order(iface)
          vertex_normal(1:3,ivert,iface) = - vertex_normal(1:3,ivert,iface)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INTERACT - Note:'
      write ( *, '(a)' ) '  Reversed node, face and vertex normals.'
!
!  RELAX: Smooth the surface via relaxation.
!
    else if ( s_eqi ( command(1:3), 'REL' ) ) then

      call node_relax ( cor3, cor3_max, cor3_new, cor3_num, face, face_max, &
        face_num, face_order, order_max )

      cor3(1:3,1:cor3_num) = cor3_new(1:3,1:cor3_num)
!
!  S: Select a few faces, discard rest:
!
    else if ( s_eqi ( command(1:1), 'S' ) ) then
 
      call face_subset ( cor3, cor3_max, cor3_num, face, face_material, &
        face_max, face_normal, face_num, face_order, ierror, line_num, &
        list, order_max, vertex_material, vertex_normal )
!
!  T: Transform data.
!
    else if ( s_eqi ( command(1:1), 'T' ) ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'For now, we only offer point scaling.'
      write ( *, '(a)' ) 'Enter X, Y, Z scale factors:'
      read ( *, * ) x, y, z

      cor3(1,1:cor3_num) = cor3(1,1:cor3_num) * x
      cor3(2,1:cor3_num) = cor3(2,1:cor3_num) * y
      cor3(3,1:cor3_num) = cor3(3,1:cor3_num) * z

      call cor3_range ( cor3, cor3_max, cor3_num )

      do iface = 1, face_num
        do ivert = 1, face_order(iface)
          vertex_normal(1:3,ivert,iface) = 0.0E+00
        end do
      end do

      call vertex_normal_set ( cor3, cor3_max, face, face_max, &
        face_num, face_order, order_max, vertex_normal )

      cor3_normal(1:3,1:cor3_num) = 0.0E+00

      call cor3_normal_set ( cor3_max, cor3_normal, face, face_area, &
        face_max, face_num, face_order, order_max, vertex_normal )

      face_normal(1:3,1:face_num) = 0.0E+00

      call face_normal_ave ( face_max, face_normal, face_num, face_order, &
        order_max, vertex_normal )
!
!  U: Renumber faces, count objects:
!
    else if ( s_eqi ( command(1:1), 'U' ) ) then
 
      call face_check ( edge, edge_max, edge_num, face, face_material, &
        face_max, face_normal, face_num, face_object, face_order, face_rank, &
        face_tier, object_num, order_max, vertex_material, vertex_normal )
!
!  V: Convert polygons to triangles.
!
    else if ( s_eqi ( command(1:1), 'V' ) ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Convert polygonal faces to triangles.'
 
      call poly_2_tri ( face, face_material, face_max, face_num, face_order, &
        ierror, order_max, vertex_material )

      if ( ierror /= 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Conversion attempt abandoned.'

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) 'Number of faces is now ', face_num

        do iface = 1, face_num
          do ivert = 1, face_order(iface)
            vertex_normal(1:3,ivert,iface) = 0.0E+00
          end do
        end do

        call vertex_normal_set ( cor3, cor3_max, face, face_max, &
          face_num, face_order, order_max, vertex_normal )

        cor3_normal(1:3,1:cor3_num) = 0.0E+00

        call cor3_normal_set ( cor3_max, cor3_normal, face, face_area, &
          face_max, face_num, face_order, order_max, vertex_normal )

        face_normal(1:3,1:face_num) = 0.0E+00

        call face_normal_ave ( face_max, face_normal, face_num, face_order, &
          order_max, vertex_normal )

      end if
!
!  W: Reverse the order of the nodes that define each face.
!
    else if ( s_eqi ( command(1:1), 'W' ) ) then
 
      call face_reverse_order ( cor3_max, cor3_normal, cor3_num, face, &
        face_max, face_normal, face_num, face_order, order_max, &
        vertex_material, vertex_normal, vertex_tex_uv )
!
!  Unintelligible!
!
    else
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INTERACT - Warning!'
      write ( *, '(a)' ) '  Your command was not recognized:'
      write ( *, '(a)' ) trim ( command )
 
    end if
 
    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INTERACT - Warning!'
      write ( *, '(a)' ) '  An error occurred during this action.'
      ierror = 0
    end if

  end do

  return
end
subroutine intnex ( line, ival, done )

!*****************************************************************************80
!
!! INTNEX "reads" integers from a string, one at a time.
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
!    Input, character ( len = * ) LINE, a string, presumably containing
!    integer ( kind = 4 )s.  These may be separated by spaces or commas.
!
!    Output, integer ( kind = 4 ) IVAL.  If DONE is FALSE, then IVAL contains the
!    "next" integer read from LINE.  If DONE is TRUE, then
!    IVAL is zero.
!
!    Input/output, logical DONE.
!    On input with a fresh value of LINE, the user should set
!    DONE to TRUE.
!    On output, the routine sets DONE to FALSE if another integer
!    was read, or TRUE if no more integers could be read.
!
  implicit none

  logical done
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) lchar
  character ( len = * ) line
  integer ( kind = 4 ), save :: next = 1

  ival = 0

  if ( done ) then
    next = 1
    done = .false.
  end if

  if ( len(line) < next ) then
    done = .true.
    return
  end if

  call s_to_i4 ( line(next:), ival, ierror, lchar )

  if ( ierror /= 0 .or. lchar == 0 ) then
    done = .true.
    next = 1
  else
    done = .false.
    next = next + lchar
  end if

  return
end
subroutine iv_point_write ( cor3, cor3_max, cor3_num, filein_name, &
  fileout_name, line_dex, line_max, line_num, iunit )

!*****************************************************************************80
!
!! IV_POINT_WRITE writes point and line data to an Inventor file.
!
!  Discussion:
!
!    This routine is only called when there are no faces.  In that
!    case, extra effort is made to display the points and lines which
!    ordinarily are not displayed when the faces carry the visual information.
!
!  Example:
!
!     #Inventor V2.0 ascii
!
!     Separator {
!       Info {
!         string "Inventor file generated by IVREAD.
!         Original data in file cube.iv."
!       }
!       Separator {
!         LightModel {
!           model PHONG
!         }
!         Coordinate3 {
!           point [
!                8.59816       5.55317      -3.05561,
!                8.59816       2.49756      0.000000E+00,
!                ...etc...
!                2.48695       2.49756      -3.05561,
!           ]
!         }
!         Material { diffuseColor 0.0 1.0 0.0 }
!         IndexedLineSet {
!           coordIndex [
!              0,    1,    2,   -1,    3,    4,    5,   -1,    7,    8,    9,
!            ...etc...
!            191,   -1,
!           ]
!         }
!         Material { diffuseColor 1.0 0.0 0.0 }
!         Transform { translation 1.0 1.0 1.0 }
!         Sphere { radius 0.05 }
!         Transform { translation 0.1 0.0 0.0 }
!         Sphere { radius 0.05 }
!         Transform { translation 0.0 0.1 0.0 }
!         Sphere { radius 0.05 }
!         Transform { translation 0.0 0.0 0.1 }
!         Sphere { radius 0.05 }
!       }
!     }
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) line_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  character ( len = * ) filein_name
  character ( len = * ) fileout_name
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) length
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_num
  character ( len = 100 ) text
  integer ( kind = 4 ) text_num
  character ( len = 20 ) word

  text_num = 0

  write ( iunit, '(a)' ) '#Inventor V2.0 ascii'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) 'Separator {'
  write ( iunit, '(a)' ) '  Info {'
  write ( iunit, '(a)' ) '    string "' // trim ( fileout_name ) &
    // ' generated by IVREAD."'
  write ( iunit, '(a)' ) '    string "Original data in file ' // &
    trim ( filein_name ) // '."'
  write ( iunit, '(a)' ) '  }'
  write ( iunit, '(a)' ) '  Separator {'
  text_num = text_num + 8

  write ( iunit, '(a)' ) '    LightModel {'
  write ( iunit, '(a)' ) '      model PHONG'
  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 3
!
!  Point coordinates.
!
  write ( iunit, '(a)' ) '    Coordinate3 {'
  write ( iunit, '(a)' ) '      point ['
  text_num = text_num + 2
 
  do icor3 = 1, cor3_num
    write ( text, '(3f12.4,'','')' ) cor3(1:3,icor3)
    call s_blanks_delete ( text )
    write ( iunit, '(8x,a)' ) trim ( text )
    text_num = text_num + 1
  end do

  write ( iunit, '(a)' ) '      ]'
  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 2
!
!  IndexedLineSet
!
  if ( 0 < line_num ) then

    write ( iunit, '(a)' ) '    Material { diffuseColor 0.0 1.0 0.0}'
    write ( iunit, '(a)' ) '    IndexedLineSet {'
!
!  IndexedLineSet coordIndex
!
    write ( iunit, '(a)' ) '      coordIndex ['
    text_num = text_num + 2

    text = ' '
    length = 0

    do j = 1, line_num
 
      write ( word, '(i8,'','')' ) line_dex(j) - 1
      call s_cat ( text, word, text )
      length = length + 1

      if ( line_dex(j)-1 == -1 .or. length >= 10 .or. j >= line_num ) then

        call s_blanks_delete ( text )
        write ( iunit, '(8x,a)' ) trim ( text )
        text_num = text_num + 1
        text = ' '
        length = 0

      end if

    end do
 
    write ( iunit, '(a)' ) '      ]'
    write ( iunit, '(a)' ) '    }'
    text_num = text_num + 2

  end if
!
!  Draw Spheres.
!
  write ( iunit, '(a)' ) '    Material { diffuseColor 1.0 0.0 0.0}'
  text_num = text_num + 1

  do icor3 = 1, cor3_num

    if ( icor3 == 1 ) then
      write ( text, '(3f12.4)' ) cor3(1:3,icor3)
    else
      write ( text, '(3f12.4)' ) cor3(1:3,icor3) - cor3(1:3,icor3-1)
    end if

    call s_blanks_delete ( text )

    write ( iunit, '(a)' ) '    Transform { translation ' // &
      trim ( text ) // '}'
    write ( iunit, '(a)' ) '    Sphere { radius 0.03 }'
    text_num = text_num + 2
  end do
!
!  Close up the Separator node.
!
  write ( iunit, '(a)' ) '  }'
!
!  Close up the Separator node.
!
  write ( iunit, '(a)' ) '}'

  text_num = text_num + 2
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'IV_POINT_WRITE - Wrote ', text_num, &
    ' text lines to ' // trim ( fileout_name )

  return
end
subroutine iv_read ( bad_num, cor3, cor3_max, cor3_num, &
  debug, face, face_max, face_num, face_order, filein_name, ierror, &
  iunit, line_dex, line_material, line_max, line_num, material_max, &
  material_num, material_rgba, order_max, text_num, &
  texture_max, texture_name, texture_num, texture_temp, &
  vertex_material, vertex_normal, vertex_tex_uv )

!*****************************************************************************80
!
!! IV_READ reads graphics information from an Inventor file.
!
!  Diagnostics:
!
!    For now, we are going to apply the following kludge, which is an
!    improvement over the previous situation.  The transform matrix
!    will be initialized to the identity; every time a new transform
!    matrix is specified, it will OVERWRITE the old one.  Every point
!    and vector that is read in will be multiplied by the current
!    transform matrix.  That's it for now.  We need to start using
!    the transform matrix, and eventually, we need to start using
!    it more accurately than this.
!
!  Discussion:
!
!    It is intended that the information read from the file can
!    either start a whole new graphics object, or simply be added
!    to a current graphics object via the '<<' command.  
!
!    This is controlled by whether the input values have been zeroed 
!    out or not.  This routine simply tacks on the information it 
!    finds to the current graphics object.
!
!  Example:
!
!     #Inventor V2.0 ascii
!
!     Separator {
!       Info {
!         string "Inventor file generated by IVREAD.
!         Original data in file cube.iv."
!       }
!       Separator {
!         LightModel {
!           model PHONG
!         }
!         MatrixTransform { matrix
!           0.9  0.0  0.0  0.0E+00
!           0.0 -0.9  0.0  0.0E+00
!           0.0  0.0 -1.5  0.0E+00
!           0.0  0.0  0.0  1.0E+00
!         }
!         Material {
!           ambientColor  0.2 0.2 0.2
!           diffuseColor  [
!             0.8 0.8 0.8,
!             0.7 0.1 0.1,
!             0.1 0.8 0.2,
!           ]
!           emissiveColor 0.0 0.0 0.0E+00
!           specularColor 0.0 0.0 0.0E+00
!           shininess     0.2
!           transparency  [
!             0.0, 0.5, 1.0,
!           ]
!         }
!         Texture2 {
!           filename      "fred.rgb"
!           wrapS         REPEAT
!           wrapT         REPEAT
!           model         MODULATE
!           blendColor    0.0 0.0 0.0E+00
!         }
!         TextureCoordinateBinding {
!           value PER_VERTEX_INDEXED
!         }
!         MaterialBinding {
!           value PER_VERTEX_INDEXED
!         }
!         NormalBinding {
!           value PER_VERTEX_INDEXED
!         }
!         ShapeHints {
!           vertexOrdering COUNTERCLOCKWISE
!           shapeType UNKNOWN_SHAPE_TYPE
!           faceType CONVEX
!           creaseAngle 6.28319
!         }
!         Coordinate3 {
!           point [
!                8.59816       5.55317      -3.05561,
!                8.59816       2.49756      0.000000E+00,
!                ...etc...
!                2.48695       2.49756      -3.05561,
!           ]
!         }
!         TextureCoordinate2 {
!           point [
!                0.0  1.0,
!                0.1, 0.8,
!                ...etc...
!                0.4  0.7,
!           ]
!         }
!         Normal {
!           vector [
!             0.71 0.71 0.0,
!             ...etc...
!             0.32 0.32 0.41,
!           ]
!         }
!
!         IndexedLineSet {
!           coordIndex [
!              0,    1,    2,   -1,    3,    4,    5,   -1,    7,    8,    9,
!            ...etc...
!            191,   -1,
!           ]
!           materialIndex [
!              0,    0,    0,   -1,    1,    1,    1,   -1,    2,    2,    2,
!            ...etc...
!             64,   -1,
!           ]
!         }
!
!         IndexedFaceSet {
!           coordIndex [
!              0,    1,    2,   -1,    3,    4,    5,   -1,    7,    8,    9,
!            ...etc...
!            191,   -1,
!           ]
!           materialIndex [
!              0,    0,    0,   -1,    1,    1,    1,   -1,    2,    2,    2,
!            ...etc...
!             64,   -1,
!           ]
!           normalIndex [
!              0,    0,    0,   -1,    1,    1,    1,   -1,    2,    2,    2,
!            ...etc...
!             64,   -1,
!           ]
!           textureCoordIndex [
!              0,    0,    0,   -1,    1,    1,    1,   -1,    2,    2,    2,
!            ...etc...
!             64,   -1,
!           ]
!         }
!
!         IndexedTriangleStripSet {
!           vertexProperty VertexProperty {
!             vertex [ x y z,
!                      ...
!                      x y z ]
!             normal [ x y z,
!                      ...
!                      x y z ]
!             materialBinding OVERALL
!             normalBinding PER_VERTEX_INDEXED
!           }
!           coordIndex [ i, j, k, l, m, -1, n, o, p, q, r, s, t, u, -1,
!             v, w, x, -1 ..., -1 ]
!           normalIndex -1
!         }
!
!       }
!     }
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BAD_NUM, number of bad lines of text read.
!
!    Input/output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input/output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, logical DEBUG, debugging switch.
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input/output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices 
!    per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which data is read.
!
!    Input/output, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, 
!    terminated by -1.
!
!    Input/output, integer ( kind = 4 ) LINE_MATERIAL(LINE_MAX), material index for line.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer ( kind = 4 ) TEXTURE_MAX, the maximum number of textures.
!
!    Input/output, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Input/output, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Input/output, integer ( kind = 4 ) TEXTURE_NUM, the number of textures.
!
!    Output, integer ( kind = 4 ) TEXT_NUM, number of lines in the file.
!
!    Input/output, character ( len = * ) TEXTURE_NAME(TEXTURE_MAX), texture 
!    names.
!
!    Workspace, real TEXTURE_TEMP(2,ORDER_MAX*FACE_MAX), texture coordinates.
!
!    Output, real ( kind = 4 ) TRANSFORM_MATRIX(4,4), the transformation matrix.
!
!    Input/output, integer ( kind = 4 ) VERTEX_MATERIAL(ORDER_MAX,FACE_MAX), vertex
!    materials.
!
!    Output, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
!    Input/output, real ( kind = 4 ) VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture 
!    coordinates.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ), parameter :: level_max = 10
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) texture_max

  real ( kind = 4 ) b
  integer ( kind = 4 ) bad_num
  character ( len = 4 ) char4
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) cor3_num_old
  logical debug
  logical done
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_num2
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  real ( kind = 4 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) icolor
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) iface_num
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iprint
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) iuv
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) iword
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) ixyz
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) iz
  integer ( kind = 4 ) jval
  integer ( kind = 4 ) k
  integer ( kind = 4 ) level
  character ( len = 256 ) level_name(0:level_max)
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) line_num2
  logical lval
! character ( len = 30 ) material_binding
  integer ( kind = 4 ) material_default
  integer ( kind = 4 ) material_num
  integer ( kind = 4 ) material_num_old
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ) nlbrack
  integer ( kind = 4 ) normal_bad_num
! character ( len = 30 ) normal_binding
  real ( kind = 4 ) normal_temp(3,order_max*face_max)
  integer ( kind = 4 ) normal_temp_num
  integer ( kind = 4 ) nrbrack
  integer ( kind = 4 ) nu
  integer ( kind = 4 ) nv
  integer ( kind = 4 ), parameter :: OFFSET = 1
  real ( kind = 4 ) r
  real ( kind = 4 ) r3vec(3)
  real ( kind = 4 ) rgb(3)
  real ( kind = 4 ) rval
  logical s_eqi
  logical s_is_i4
  character ( len = 256 ) text
  integer ( kind = 4 ) text_num
! character ( len = 30 ) texture_coordinate_binding
  character ( len = * ) texture_name(texture_max)
  integer ( kind = 4 ) texture_num
  real ( kind = 4 ) texture_temp(2,order_max*face_max)
  integer ( kind = 4 ) texture_temp_num
  real ( kind = 4 ) transform_matrix(4,4)
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  real ( kind = 4 ) vertex_tex_uv(2,order_max,face_max)
  character ( len = 256 ) word
  character ( len = 256 ) wordm1

  cor3_num_old = cor3_num
  face_num2 = face_num
  icol = 0
  ierror = 0
  iface_num = face_num
  iprint = 1
  irow = 0
  iuv = 0
  ix = 0
  ixyz = 0
  iy = 0
  iz = 0
  jval = 0
  level = 0
  level_name(0) = 'Top'
  line_num2 = line_num
  material_num_old = material_num
! material_binding = 'PER_VERTEX_INDEXED'
  material_default = 1
  nlbrack = 0
  normal_bad_num = 0
! normal_binding = 'PER_VERTEX_INDEXED'
  normal_temp_num = 0
  nrbrack = 0
  nu = 0
  nv = 0
  rval = 0.0E+00
  texture_temp_num = 0
  call tmat_init ( transform_matrix )
  word = ' '
  wordm1 = ' '
!
!  Read a line of text from the file.
!
  do
 
    read ( iunit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      exit
    end if

    text_num = text_num + 1

    if ( debug ) then

      if ( text_num == iprint ) then
        write ( *, '(a,i6)' ) 'Line ', iprint
        iprint = 2 * iprint
      end if

    end if
 
    if ( text == ' ' ) then
      cycle
    end if
!
!  The first line of the file must be the header.
!
    if ( text_num == 1 ) then

      if ( .not. s_eqi ( text(1:9), '#Inventor' ) ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IV_READ - Fatal error!'
        write ( *, '(a)' ) '  The input file has a bad header.'
        write ( *, '(a)' ) trim ( text )
        return
      else
        cycle
      end if

    end if
 
    done = .true.
    iword = 0
!
!  Save the previous word read.  It helps when a word depends
!  on its context.
!
20    continue

    if ( word /= ' ' .and. word /= ',' ) then
      wordm1 = word
    end if
!
!  Read a word from the line.
!
    call word_next_read ( text, word, done )
!
!  If no more words in this line, read in a whole new line.
!
    if ( done ) then
      cycle
    end if
!
!  Skip over comments.
!
    if ( word(1:1) == '#' ) then
      cycle
    end if
!
!  Ignore blanks and commas.
!
    if ( word == ' ' .or. word == ',' ) then
      go to 20
    end if
!
!  Count the words in the current line, and the total.
!
    iword = iword + 1
!
!  If the word is a curly or square bracket, count it.
!
    if ( word == '{' .or. word == '[' ) then

      nlbrack = nlbrack + 1

      if ( debug ) then
        write ( *, '(a)' ) word(1:1)
      end if

    else if ( word .eq. '}' .or. word == ']' ) then

      if ( debug ) then
        write ( *, '(a)' ) word(1:1)
      end if

      nrbrack = nrbrack + 1

      if ( nlbrack < nrbrack ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IV_READ - Fatal error!'
        write ( *, '(a,i6)' ) '  Extraneous right bracket, line ', text_num
        write ( *, '(a)' ) trim ( text )
        write ( *, '(a)' ) 'Currently processing field:'
        write ( *, '(a)' ) trim ( level_name(level) )
        ierror = 1
        return
      end if

    end if
!
!  If the word is a left bracket, then the previous word
!  is the name of a node.
!
    if ( word == '{' .or. word == '[' ) then
 
      level = nlbrack - nrbrack
      if ( level < 0 ) then
        write ( *, '(a)' ) 'Too many right brackets!'
        level = 0
      else if ( level_max < level ) then
        write ( *, '(a)' ) 'Too many left brackets!'
        level = level_max
      end if

      level_name(level) = wordm1
 
      if ( debug ) then
        write ( *, '(a)' ) ' '
        do i = 0, level

          write ( *, '(i3,2x,a)' ) i, trim ( level_name(i) )
        end do
      end if
 
    end if
!
!  BASECOLOR
!
    if ( s_eqi ( level_name(level), 'BASECOLOR' ) ) then
 
      if ( word == '{' ) then

      else if ( word == '}' ) then

        level = nlbrack - nrbrack

      else if ( s_eqi ( word, 'RGB' ) ) then

      else

        go to 99

      end if
!
!  COORDINATE3
!
    else if ( s_eqi ( level_name(level), 'COORDINATE3' ) ) then
 
      if ( word == '{' ) then

      else if ( word == '}' ) then

        level = nlbrack - nrbrack

      else if ( s_eqi ( word, 'POINT' ) ) then

      else 
        go to 99
      end if
!
!  COORDINATE4
!
    else if ( s_eqi ( level_name(level), 'COORDINATE4' ) ) then
 
      if ( word == '{' ) then
 
      else if ( word == '}' ) then

        level = nlbrack - nrbrack

      else if ( s_eqi ( word, 'POINT' ) ) then

      else
        go to 99
      end if
!
!  COORDINDEX
!
    else if ( s_eqi ( level_name(level), 'COORDINDEX' ) ) then 
 
      if ( word == '[' ) then

        ivert = 0

      else if ( word == ']' ) then

        level = nlbrack - nrbrack
!
!  (indexedlineset) COORDINDEX
!
      else if ( s_eqi ( level_name(level-1), 'INDEXEDLINESET' ) ) then

        lval = s_is_i4 ( word, jval )
 
        if ( lval ) then

          if ( jval < -1 ) then

            bad_num = bad_num + 1

          else

            line_num = line_num + 1

            if ( line_num <= line_max ) then

              if ( jval == -1 ) then
                line_material(line_num) = -1 + OFFSET
              else
                line_material(line_num) = material_default
                jval = jval + cor3_num_old
              end if

              line_dex(line_num) = jval + OFFSET

            end if

          end if

        else
          go to 99
        end if
!
!  (indexedfaceset) COORDINDEX
!
      else if ( s_eqi ( level_name(level-1), 'INDEXEDFACESET' ) ) then
 
        lval = s_is_i4 ( word, jval )
 
        if ( lval ) then

          if ( jval == -1 ) then

            ivert = 0

          else

            ivert = ivert + 1

            if ( ivert == 1 ) then
              face_num = face_num + 1
              if ( face_num <= face_max ) then
                face_order(face_num) = 0
              end if
            end if

            if ( face_num <= face_max ) then
              face_order(face_num) = face_order(face_num) + 1
              if ( ivert <= order_max ) then
                face(ivert,face_num) = jval + cor3_num_old + OFFSET
              end if
            end if

            vertex_material(ivert,face_num) = material_default

          end if

        end if
!
!  (indexednurbssurface) COORDINDEX
!
      else if ( s_eqi ( level_name(level-1), 'INDEXEDNURBSSURFACE' ) ) then

        lval = s_is_i4 ( word, jval )

        if ( lval ) then

        else if ( word == '[' ) then

        else if ( word == ']' ) then
 
          do i = 1, nu-1
            line_num = line_num + 1
            if ( line_num <= line_max ) then
              line_dex(line_num) = i - 1 + cor3_num_old + OFFSET
            end if
          end do
 
          do i = nu, nu*(nv-1), nu
            line_num = line_num + 1
            if ( line_num <= line_max ) then
              line_dex(line_num) = i - 1 + cor3_num_old + OFFSET
            end if
          end do
 
          do i = nu*nv, nu*nv-nu+2, -1
            line_num = line_num + 1
            if ( line_num <= line_max ) then
              line_dex(line_num) = i - 1 + cor3_num_old + OFFSET
            end if
          end do
 
          do i = nu*(nv-1)+1, 1, -nu
            line_num = line_num + 1
            if ( line_num <= line_max ) then
              line_dex(line_num) = i - 1 + cor3_num_old + OFFSET
            end if
          end do
 
          line_num = line_num + 1
          if ( line_num <= line_max ) then
            line_dex(line_num) = -1 + OFFSET
          end if

        end if
!
!  (indexedtrianglestripset) COORDINDEX
!
!  First three coordinate indices I1, I2, I3 define a triangle.
!  Next triangle is defined by I2, I3, I4 (actually, I4, I3, I2
!  to stay with same counterclockwise sense).
!  Next triangle is defined by I3, I4, I5 ( don't need to reverse
!  odd numbered triangles) and so on.
!  List is terminated with -1.
!
      else if ( s_eqi ( level_name(level-1), 'INDEXEDTRIANGLESTRIPSET' ) ) then

        lval = s_is_i4 ( word, jval )
 
        if ( lval ) then

          if ( jval == -1 ) then

            ivert = 0

          else

            ivert = ivert + 1

            ix = iy
            iy = iz 
            iz = jval + cor3_num_old

            if ( ivert == 1 ) then

              face_num = face_num + 1
              if ( face_num <= face_max ) then
                face(ivert,face_num) = jval + cor3_num_old + OFFSET
                face_order(face_num) = 3
              end if

            else if ( ivert == 2 .or. ivert == 3 ) then

              if ( face_num <= face_max ) then
                face(ivert,face_num) = jval + cor3_num_old + OFFSET
              end if

            else

              face_num = face_num + 1
              if ( face_num <= face_max ) then
                face_order(face_num) = 3
                if ( mod ( ivert, 2 ) == 1 ) then
                  face(1,face_num) = ix + OFFSET
                  face(2,face_num) = iy + OFFSET
                  face(3,face_num) = iz + OFFSET
                else
                  face(1,face_num) = iz + OFFSET
                  face(2,face_num) = iy + OFFSET
                  face(3,face_num) = ix + OFFSET
                end if
              end if

            end if
!
!  ??? This can't be right.???
!
!  Very very tentative guess as to how indices into the normal
!  vector array are set up...
!
!            if ( face_num <= face_max .and. ivert >= 3 ) then
!               do j = 1, order_max
!                 vertex_normal(1:3,j,face_num) = normal(1:3,ix + OFFSET)
!               end do
!             end if

          end if

        end if

      end if
!
!  INDEXEDFACESET
!
    else if ( s_eqi ( level_name(level), 'INDEXEDFACESET' ) ) then
 
      if ( word == '{' ) then
 
      else if ( word == '}' ) then

        level = nlbrack - nrbrack

      else if ( s_eqi ( word, 'COORDINDEX' ) ) then

        ivert = 0
!
!  19 October 2001
!
!  If material binding is per node indexed...
! 
      else if ( s_eqi ( word, 'MATERIALINDEX' ) ) then

      else if ( s_eqi ( word, 'NORMALINDEX' ) ) then

      else if ( s_eqi ( word, 'TEXTURECOORDINDEX' ) ) then

        if ( texture_num <= 0 ) then
          texture_num = 1
          call i4_to_s_zero ( texture_num, char4 )
          texture_name(texture_num) = 'Texture_' // char4
        end if

      else
        go to 99
      end if
!
!  INDEXEDLINESET
!
    else if ( s_eqi ( level_name(level), 'INDEXEDLINESET' ) ) then 
 
      if ( word == '{' ) then

      else if ( word == '}' ) then

        level = nlbrack - nrbrack
 
      else if ( s_eqi ( word, 'COORDINDEX' ) ) then
!
!  19 October 2001:
!  If I start collecting this data, then I have to allow lines to
!  have materials by line, by segment, or by vertex...
!
      else if ( s_eqi ( word, 'MATERIALINDEX' ) ) then

      else
        go to 99
      end if
!
!  INDEXEDNURBSSURFACE
!
    else if ( s_eqi ( level_name(level), 'INDEXEDNURBSSURFACE' ) ) then
 
      if ( word == '{' ) then

      else if ( word == '}' ) then

        level = nlbrack - nrbrack
 
      else if ( s_eqi ( word, 'NUMUCONTROLPOINTS') ) then

        call word_next_read ( text, word, done )
        lval = s_is_i4 ( word, jval )

        if ( lval ) then
          nu = jval
        else
          nu = 0
          go to 99
        end if

      else if ( s_eqi ( word, 'NUMVCONTROLPOINTS' ) ) then

        call word_next_read ( text, word, done)
        lval = s_is_i4 ( word, jval ) 

        if ( lval ) then
          nv = jval
        else
          nv = 0
          go to 99
        end if

      else if ( s_eqi ( word, 'COORDINDEX' ) ) then

      else if ( s_eqi ( word, 'UKNOTVECTOR' ) ) then

      else if ( s_eqi ( word, 'VKNOTVECTOR' ) ) then

      else
        go to 99
      end if
!
!  INDEXEDTRIANGLESTRIPSET
!
    else if ( s_eqi ( level_name(level), 'INDEXEDTRIANGLESTRIPSET' ) ) then 
 
      if ( word == '{' ) then

      else if ( word == '}' ) then

        level = nlbrack - nrbrack
 
      else if ( s_eqi ( word, 'VERTEXPROPERTY' ) ) then

        call word_next_read ( text, word, done )

      else if ( s_eqi ( word, 'COORDINDEX' ) ) then

        ivert = 0

      else if ( s_eqi ( word, 'NORMALINDEX' ) ) then

        call word_next_read ( text, word, done )

        if ( word == '[' ) then

          nlbrack = nlbrack + 1
          level = level + 1
          level_name(level) = 'NORMALINDEX'

        else if ( word == '-1' ) then

          do iface = 1, face_num
            do ivert = 1, face_order(iface)
              k = face(ivert,iface)
              vertex_normal(1:3,ivert,iface) = normal_temp(1:3,k)
            end do
          end do

        end if

      else
        go to 99
      end if
!
!  INFO
!
    else if ( s_eqi ( level_name(level), 'INFO' ) ) then

      if ( word == '{' ) then
 
      else if ( word == '}' ) then

        level = nlbrack - nrbrack

      else if ( s_eqi ( word, 'STRING' ) ) then

      else if ( word == '"' ) then

      else

      end if
!
!  LIGHTMODEL
!
!  Read, but ignore.
!
    else if ( s_eqi ( level_name(level), 'LIGHTMODEL' ) ) then

      if ( word == '{' ) then

      else if ( word == '}' ) then

        level = nlbrack - nrbrack

      else if ( s_eqi ( word, 'MODEL' ) ) then

      else

      end if
!
!  MATERIAL
!
!  Read, but ignore for now, except that the ambient color information
!  is used to set (default) material #1.
!
    else if ( s_eqi ( level_name(level), 'MATERIAL' ) ) then

      if ( word == '{' ) then

      else if ( word == '}' ) then

        level = nlbrack - nrbrack

      else if ( s_eqi ( word, 'AMBIENTCOLOR' ) ) then

        call word_next_read ( text, word, done )
        call s_is_r4 ( word, r, lval )
        call word_next_read ( text, word, done )
        call s_is_r4 ( word, g, lval )
        call word_next_read ( text, word, done )
        call s_is_r4 ( word, b, lval )

        material_rgba(1:3,1) = (/ r, g, b /)

      else if ( s_eqi ( word, 'EMISSIVECOLOR' ) ) then

        call word_next_read ( text, word, done )
        call s_is_r4 ( word, r, lval )
        call word_next_read ( text, word, done )
        call s_is_r4 ( word, g, lval )
        call word_next_read ( text, word, done )
        call s_is_r4 ( word, b, lval )

      else if ( s_eqi ( word, 'DIFFUSECOLOR' ) ) then

        call word_next_read ( text, word, done )
        call s_is_r4 ( word, r, lval )
        call word_next_read ( text, word, done )
        call s_is_r4 ( word, g, lval )
        call word_next_read ( text, word, done )
        call s_is_r4 ( word, b, lval )

      else if ( s_eqi ( word, 'SHININESS' ) ) then

        call word_next_read ( text, word, done )
        call s_is_r4 ( word, r, lval )

      else if ( s_eqi ( word, 'SPECULARCOLOR' ) ) then

        call word_next_read ( text, word, done )
        call s_is_r4 ( word, r, lval )
        call word_next_read ( text, word, done )
        call s_is_r4 ( word, g, lval )
        call word_next_read ( text, word, done )
        call s_is_r4 ( word, b, lval )

      else if ( s_eqi ( word, 'TRANSPARENCY' ) ) then

        call word_next_read ( text, word, done )
        call s_is_r4 ( word, r, lval )

      else

      end if
!
!  MATERIALBINDING
!
!    OVERALL: Whole object has same material;
!    PER_FACE: One material for each face;
!    PER_FACE_INDEXED: One material for each face, indexed;
!    PER_PART: One material for each part;
!    PER_PART_INDEXED: One material for each part, indexed;
!    PER_VERTEX: One material for each vertex;
!    PER_VERTEX_INDEXED: One material for each vertex, indexed.
!
    else if ( s_eqi ( level_name(level), 'MATERIALBINDING' ) ) then
    
      if ( word == '{' ) then

      else if ( word == '}' ) then
        level = nlbrack - nrbrack
      else if ( s_eqi ( word, 'VALUE' ) ) then
        call word_next_read ( text, word, done )
!         material_binding = word
      else
        bad_num = bad_num + 1
        write ( *, '(a)' ) '  Bad data on level ' // trim ( level_name(level) )
        write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        write ( *, '(a,i6)' ) '  Line number: ', text_num
        write ( *, '(a)' ) trim ( text )
      end if
!
!  MATERIALINDEX
!
    else if ( s_eqi ( level_name(level), 'MATERIALINDEX' ) ) then
    
      if ( word == '[' ) then

        ivert = 0

      else if ( word == ']' ) then

        level = nlbrack - nrbrack
!
!  (indexedfaceset) MATERIALINDEX
!
      else if ( s_eqi ( level_name(level-1), 'INDEXEDFACESET' ) ) then

        lval = s_is_i4 ( word, jval )
 
        if ( lval ) then

          if ( jval == -1 ) then

            ivert = 0

          else

            ivert = ivert + 1

            if ( ivert == 1 ) then
              face_num2 = face_num2 + 1
            end if

            if ( face_num2 <= face_max ) then
              if ( jval /= -1 ) then
                jval = jval + cor3_num_old
              end if
              vertex_material(ivert,face_num2) = jval + OFFSET
            end if

          end if

        else
          go to 99
        end if
!
!  (indexedlineset) MATERIALINDEX
!
      else if ( s_eqi ( level_name(level-1), 'INDEXEDLINESET' ) ) then

        lval = s_is_i4 ( word, jval )
 
        if ( lval ) then

          line_num2 = line_num2 + 1
          if ( line_num2 <= line_max ) then
            if ( jval /= -1 ) then
!             jval = jval + cor3_num_old
              jval = jval + material_num_old
            end if
            line_material(line_num2) = jval + OFFSET
          end if

        else
          go to 99
        end if

      else
 
        lval = s_is_i4 ( word, jval )
 
        if ( lval ) then

        else
          go to 99
        end if

      end if
!
!  MATRIXTRANSFORM
!
    else if ( s_eqi ( level_name(level), 'MATRIXTRANSFORM' ) ) then

        if ( word == '{' ) then

        else if ( word == '}' ) then

          if ( irow /= 4 .or. icol .ne. 4 ) then
            write ( *, '(a)' ) 'Incomplete MatrixTransform!'
          end if

          level = nlbrack - nrbrack

        else if ( s_eqi ( word, 'MATRIX' ) ) then

          irow = 1
          icol = 0

        else

          call s_is_r4 ( word, rval, lval )

          if ( lval ) then

            icol = icol + 1
            if ( 4 < icol ) then
              icol = 1
              irow = irow + 1
              if ( 4 < irow ) then
                go to 99
              end if
            end if

            transform_matrix(irow,icol) = rval

          else
            go to 99
          end if

        end if
!
!  NORMAL
!
!  The field "VECTOR" may be followed by three numbers,
!  (handled here),  or by a square bracket, and sets of three numbers.
!
    else if ( s_eqi ( level_name(level), 'NORMAL' ) ) then
!
!  (vertexproperty) NORMAL
!  Just stick the normal vectors into NORMAL_TEMP for now,
!  retrieve them at COORDINDEX time.
!
      if ( s_eqi ( level_name(level-1), 'VERTEXPROPERTY' ) ) then

        if ( word == '[' ) then

          ixyz = 0

        else if ( word == ']' ) then

          level = nlbrack - nrbrack

        else

          call s_is_r4 ( word, rval, lval )
 
          if ( lval ) then
 
            ixyz = ixyz + 1
 
            if ( 3 < ixyz ) then
              ixyz = 1
            end if

            if ( ixyz == 1 ) then
              normal_temp_num = normal_temp_num + 1
            end if
 
            if ( normal_temp_num <= order_max*face_max ) then

              r3vec(ixyz) = rval

              if ( ixyz == 3 ) then
                call tmat_mxv ( transform_matrix, r3vec, r3vec )
                normal_temp(1:3,normal_temp_num) = r3vec(1:3)
              end if

            end if
 
          else
            go to 99
          end if

        end if
!
!  (anythingelse) NORMAL.
!
      else

        if ( word == '{' ) then

          ixyz = 0
          normal_temp_num = 0

        else if ( word == '}' ) then

          level = nlbrack - nrbrack

        else if ( s_eqi ( word, 'VECTOR' ) ) then

        else

          go to 99

        end if

      end if
!
!  NORMALBINDING 
!    OVERALL: Whole object has same normal.
!    PER_FACE: One normal per face;
!    PER_FACE_INDEXED: one normal per face, indexed;
!    PER_FACE: one normal per part;
!    PER_FACE_INDEXED: one normal per part, indexed.
!    PER_VERTEX: one normal per vertex;
!    PER_VERTEX_INDEXED: one normal per vertex, indexed.
!
    else if ( s_eqi ( level_name(level), 'NORMALBINDING' ) ) then
    
      if ( word == '{' ) then

      else if ( word == '}' ) then
        level = nlbrack - nrbrack
      else if ( s_eqi ( word, 'VALUE' ) ) then
        call word_next_read ( text, word, done )
!         normal_binding = word
      else
        bad_num = bad_num + 1
        write ( *, '(a)' ) '  Bad data on level ' // trim ( level_name(level) )
        write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        write ( *, '(a,i6)' ) '  Line number: ', text_num
        write ( *, '(a)' ) trim ( text )
      end if
!
!  NORMALINDEX
!
    else if ( s_eqi ( level_name(level), 'NORMALINDEX' ) ) then
!
!  (indexedtrianglestripset) NORMALINDEX
!
      if ( s_eqi ( level_name(level-1), 'INDEXEDTRIANGLESTRIPSET' ) ) then

        lval = s_is_i4 ( word, jval )

        if ( lval ) then

        else if ( word == '[' ) then

        else if ( word == ']' ) then
          level = nlbrack - nrbrack
        end if
!
!  (indexedFaceSet) NORMALINDEX
!
      else 

        if ( word == '[' ) then
          ivert = 0
        else if ( word == ']' ) then
          level = nlbrack - nrbrack
        else
 
          lval = s_is_i4 ( word, jval )
 
          if ( lval ) then

            if ( jval == -1 ) then

              ivert = 0

            else

              ivert = ivert + 1

              if ( ivert == 1 ) then
                iface_num = iface_num + 1
              end if

              if ( iface_num <= face_max .and. jval <= normal_temp_num ) then
                vertex_normal(1:3,ivert,iface_num) = normal_temp(1:3,jval + OFFSET)
              end if

            end if

          else
            go to 99
          end if

        end if

      end if
!
!  (coordinate3) POINT
!
    else if ( s_eqi ( level_name(level), 'POINT' ) ) then
 
      if ( s_eqi ( level_name(level-1), 'COORDINATE3' ) ) then

        if ( word == '[' ) then

          ixyz = 0
          cor3_num_old = cor3_num
 
        else if ( word == ']' ) then

          level = nlbrack - nrbrack

        else
 
          call s_is_r4 ( word, rval, lval )

          if ( lval ) then
 
            ixyz = ixyz + 1
 
            if ( 3 < ixyz ) then
              ixyz = 1
            end if
 
            r3vec(ixyz) = rval

            if ( ixyz == 3 ) then

              cor3_num = cor3_num + 1
              call tmat_mxp ( transform_matrix, r3vec, r3vec )

              if ( cor3_num <= cor3_max ) then
                cor3(1:3,cor3_num) = r3vec(1:3)
              end if

            end if
 
          else
            go to 99
          end if
 
        end if
!
!  (texturecoordinate2) POINT
!
      else if ( s_eqi ( level_name(level-1), 'TEXTURECOORDINATE2' ) ) then

        if ( word == '[' ) then

          iuv = 0
          texture_temp_num = 0
 
        else if ( word == ']' ) then

          level = nlbrack - nrbrack

        else
 
          call s_is_r4 ( word, rval, lval )

          if ( lval ) then
 
            iuv = iuv + 1
 
            if ( 2 < iuv ) then
              iuv = 1
            end if
 
            texture_temp(iuv,texture_temp_num+1) = rval

            if ( iuv == 2 ) then
              texture_temp_num = texture_temp_num + 1
            end if
 
          else
            go to 99
          end if
 
        end if

      end if
!
!  RGB
!
    else if ( s_eqi ( level_name(level), 'RGB' ) ) then

      if ( s_eqi ( level_name(level-1), 'BASECOLOR' ) ) then

        if ( word == '[' ) then
 
          icolor = 0

        else if ( word == ']' ) then

          level = nlbrack - nrbrack

        else 

          call s_is_r4 ( word, rval, lval )
 
          if ( lval ) then
 
            icolor = icolor + 1
 
            rgb(icolor) = rval

            if ( icolor == 3 ) then
              material_num = material_num + 1
              material_rgba(1:3,material_num) = rgb(1:3)
              icolor = 0
            end if

          else
            go to 99
          end if

        end if

      end if
!
!  SEPARATOR
!
    else if ( s_eqi ( level_name(level), 'SEPARATOR' ) ) then
 
      if ( word == '{' ) then
 
      else if ( word == '}' ) then

        level = nlbrack - nrbrack

      else
 
      end if
!
!  SHAPEHINTS
!
!  Read, but ignore.
!
    else if ( s_eqi ( level_name(level), 'SHAPEHINTS' ) ) then
    
      if ( word == '{' ) then

      else if ( word == '}' ) then

        level = nlbrack - nrbrack

      else if ( s_eqi ( word, 'CREASEANGLE' ) ) then

        call word_next_read ( text, word, done )
        call s_is_r4 ( word, rval, lval )
 
      else if ( s_eqi ( word, 'FACETYPE' ) ) then

        call word_next_read ( text, word, done )

      else if ( s_eqi ( word, 'SHAPETYPE' ) ) then

        call word_next_read ( text, word, done )

      else if ( s_eqi ( word, 'VERTEXORDERING' ) ) then

        call word_next_read ( text, word, done )

      else
        go to 99
      end if
!
!  TEXTURE2
!
    else if ( s_eqi ( level_name(level), 'TEXTURE2' ) ) then

      if ( word == '{' ) then

        texture_num = texture_num + 1

      else if ( word == '}' ) then

        level = nlbrack - nrbrack

      else if ( s_eqi ( word, 'BLENDCOLOR' ) ) then

        call word_next_read ( text, word, done )
        call s_is_r4 ( word, r, lval )
        call word_next_read ( text, word, done )
        call s_is_r4 ( word, g, lval )
        call word_next_read ( text, word, done )
        call s_is_r4 ( word, b, lval )

      else if ( s_eqi ( word, 'FILENAME' ) ) then

        call word_next_read ( text, word, done )
        call word_next_read ( text, word, done )
        texture_name(texture_num) = word
        call word_next_read ( text, word, done )

      else if ( s_eqi ( word, 'IMAGE' ) ) then

        go to 99

      else if ( s_eqi ( word, 'MODEL' ) ) then

        call word_next_read ( text, word, done )

      else if ( s_eqi ( word, 'WRAPS' ) ) then

        call word_next_read ( text, word, done )

      else if ( s_eqi ( word, 'WRAPT' ) ) then

        call word_next_read ( text, word, done )

      else

      end if
!
!  TEXTURECOORDINATE2
!
    else if ( s_eqi ( level_name(level), 'TEXTURECOORDINATE2' ) ) then
 
      if ( word == '{' ) then

      else if ( word == '}' ) then

        level = nlbrack - nrbrack

      else if ( s_eqi ( word, 'POINT' ) ) then

      else 
        go to 99
      end if
!
!  TEXTURECOORDINATEBINDING
!
!    PER_VERTEX: One texture coordinate for each vertex;
!    PER_VERTEX_INDEXED: One texture coordinate for each vertex, indexed.
!
    else if ( s_eqi ( level_name(level), 'TEXTURECOORDINATEBINDING' ) ) then
    
      if ( word == '{' ) then

      else if ( word == '}' ) then
        level = nlbrack - nrbrack
      else if ( s_eqi ( word, 'VALUE' ) ) then
        call word_next_read ( text, word, done )
!     texture_coordinate_binding = word
      else
        bad_num = bad_num + 1
        write ( *, '(a)' ) '  Bad data on level ' // trim ( level_name(level) )
        write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        write ( *, '(a,i6)' ) '  Line number: ', text_num
        write ( *, '(a)' ) trim ( text )
      end if
!
!  TEXTURECOORDINDEX
!
    else if ( s_eqi ( level_name(level), 'TEXTURECOORDINDEX' ) ) then

      if ( word == '[' ) then
        ivert = 0
        iface = 0
      else if ( word == ']' ) then
        level = nlbrack - nrbrack
      else
 
        lval = s_is_i4 ( word, jval )
 
        if ( lval ) then

          if ( jval == -1 ) then

            ivert = 0

          else

            ivert = ivert + 1

            if ( ivert == 1 ) then
              iface = iface + 1
            end if

            if ( iface <= face_max ) then
              vertex_tex_uv(1:2,ivert,iface) = texture_temp(1:2,jval + OFFSET)
            end if

          end if

        else
          go to 99
        end if

      end if
!
!  UKNOTVECTOR
!
    else if ( s_eqi ( level_name(level), 'UKNOTVECTOR' ) ) then

      if ( word == '[' ) then
        go to 20
      else if ( word == ']' ) then
        level = nlbrack - nrbrack
        go to 20
      else
        lval = s_is_i4 ( word, jval )
      end if
!
!  VECTOR
!
    else if ( s_eqi ( level_name(level), 'VECTOR' ) ) then

      if ( word == '[' ) then

      else if ( word == ']' ) then

        level = nlbrack - nrbrack
!
!  (normal) VECTOR.
!  For some reason, Joel's program is spewing out occasional
!  NAN normal vectors.  This should not halt the program.
!
      else if ( s_eqi ( level_name(level-1), 'NORMAL' ) ) then

        if ( word(1:13) == 'nan0x7ffffe00' ) then

          lval = .true.
          rval = 1.0E+00 / sqrt ( 3.0E+00 )
          normal_bad_num = normal_bad_num + 1

        else

          call s_is_r4 ( word, rval, lval )
 
        end if

        if ( lval ) then
 
          ixyz = ixyz + 1
 
          if ( 3 < ixyz ) then
            ixyz = 1
          end if

          if ( ixyz == 1 ) then
            normal_temp_num = normal_temp_num + 1
          end if
 
          r3vec(ixyz) = rval

          if ( ixyz == 3 ) then

            if ( normal_temp_num <= order_max * face_max ) then
              normal_temp(1:3,normal_temp_num) = r3vec(1:3)
            end if

          end if
 
        else
          go to 99
        end if

      end if
!
!  (vertexproperty) VERTEX
!
    else if ( s_eqi ( level_name(level), 'VERTEX' ) ) then

      if ( s_eqi ( level_name(level-1), 'VERTEXPROPERTY' ) ) then

        if ( word == '[' ) then

          ixyz = 0
          cor3_num_old = cor3_num

        else if ( word == ']' ) then

          level = nlbrack - nrbrack

        else

          call s_is_r4 ( word, rval, lval )
 
          if ( lval ) then
 
            ixyz = ixyz + 1
 
            if ( 3 < ixyz ) then
              ixyz = 1
            end if
 
            if ( cor3_num+1 <= cor3_max ) then
              cor3(ixyz,cor3_num+1) = rval
            end if
 
            if ( ixyz == 3 ) then
              cor3_num = cor3_num + 1
            end if
 
          else
            go to 99
          end if

        end if

      end if
!
!  (indexedtrianglestripset) VERTEXPROPERTY
!
    else if ( s_eqi ( level_name(level), 'VERTEXPROPERTY' ) ) then

      if ( word == '{' ) then

      else if ( word == '}' ) then

        level = nlbrack - nrbrack

      else if ( s_eqi ( word, 'VERTEX' ) ) then

      else if ( s_eqi ( word, 'MATERIALBINDING' ) ) then

        call word_next_read ( text, word, done )

!     material_binding = word

      else if ( s_eqi ( word, 'NORMAL' ) ) then

        ixyz = 0

      else if ( s_eqi ( word, 'NORMALBINDING' ) ) then

        call word_next_read ( text, word, done )

!     normal_binding = word

      else
        go to 99
      end if
!
!  VKNOTVECTOR
!
    else if ( s_eqi ( level_name(level), 'VKNOTVECTOR' ) ) then
      
      if ( word == '[' ) then
        go to 20
      else if ( word == ']' ) then
        level = nlbrack - nrbrack
        go to 20
      else
        lval = s_is_i4 ( word, jval )
      end if
!
!  Any other word:
!
    else
 
    end if
 
    go to 20
!
!  Bad data
!
99    continue

    bad_num = bad_num + 1

    if ( bad_num <= 10 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'IV_READ - Warning!'
      write ( *, '(a)' ) '  Bad data on level ' // trim ( level_name(level) )
      write ( *, '(a)' ) '  Bad word: ' // trim ( word )
      write ( *, '(a,i6)' ) '  Line number: ', text_num
      write ( *, '(a)' ) trim ( text )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'IV_READ - Fatal error!'
      write ( *, '(a)' ) '  Too many warnings!'
      return
    end if

  end do
!
!  End of information in file.
!
!
!  Reset the transformation matrix to the identity, cause we
!  went ahead and applied it.
!
  call tmat_init ( transform_matrix )
!
!  Check the "materials" defining a line.  
!
!  If COORDINDEX is -1, so should be the MATERIALINDEX.  
!  If COORDINDEX is not -1, then the MATERIALINDEX shouldn't be either.
!
  do i = 1, line_num

    if ( line_dex(i) == -1 + OFFSET ) then
      line_material(i) = -1 + OFFSET
    else if ( line_material(i) == -1 + OFFSET ) then
      line_material(i) = material_num
    end if

  end do
!
!  Complain once about bad entries in normal vectors.
!
  if ( 0 < normal_bad_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IV_READ - Warning!'
    write ( *, '(i6,a)' ) normal_bad_num, ' bad normal vector entries.'
  end if
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'IV_READ - Read ', text_num, ' text lines from ' &
    // trim ( filein_name )

  return
end
subroutine iv_write ( cor3, cor3_max, cor3_normal, cor3_num, face, &
  face_max, face_num, face_order, filein_name, fileout_name, &
  iunit, line_dex, line_material, line_max, line_num, material_max, &
  material_num, material_rgba, order_max, texture_max, &
  texture_num, texture_name, vertex_material, vertex_tex_uv )

!*****************************************************************************80
!
!! IV_WRITE writes graphics data to an Inventor file.
!
!  Example:
!
!     #Inventor V2.0 ascii
!
!     Separator {
!       Info {
!         string "Inventor file generated by IVREAD.
!         Original data in file cube.iv."
!       }
!       Separator {
!         LightModel {
!           model PHONG
!         }
!         MatrixTransform { matrix
!           0.9  0.0  0.0  0.0
!           0.0 -0.9  0.0  0.0
!           0.0  0.0 -1.5  0.0
!           0.0  0.0  0.0  1.0
!         }
!         Material {
!           ambientColor  0.2 0.2 0.2
!           diffuseColor  [
!             0.8 0.8 0.8,
!             0.7 0.1 0.1,
!             0.1 0.8 0.2,
!           ]
!           emissiveColor 0.0 0.0 0.0
!           specularColor 0.0 0.0 0.0
!           shininess     0.2
!           transparency  [
!             0.0, 0.5, 1.0,
!           ]
!         }
!         Texture2 {
!           filename      "fred.rgb"
!           wrapS         REPEAT
!           wrapT         REPEAT
!           model         MODULATE
!           blendColor    0.0 0.0 0.0
!         }
!         TextureCoordinateBinding {
!           value PER_VERTEX_INDEXED
!         }
!         MaterialBinding {
!           value PER_VERTEX_INDEXED
!         }
!         NormalBinding {
!           value PER_VERTEX_INDEXED
!         }
!         ShapeHints {
!           vertexOrdering COUNTERCLOCKWISE
!           shapeType UNKNOWN_SHAPE_TYPE
!           faceType CONVEX
!           creaseAngle 6.28319
!         }
!         Coordinate3 {
!           point [
!                8.59816       5.55317      -3.05561,
!                8.59816       2.49756      0.000000,
!                ...etc...
!                2.48695       2.49756      -3.05561,
!           ]
!         }
!         TextureCoordinate2 {
!           point [
!                0.0  1.0,
!                0.1, 0.8,
!                ...etc...
!                0.4  0.7,
!           ]
!         }
!         Normal {
!           vector [
!             0.71 0.71 0.0,
!             ...etc...
!             0.32 0.32 0.41,
!           ]
!         }
!
!         IndexedLineSet {
!           coordIndex [
!              0,    1,    2,   -1,    3,    4,    5,   -1,    7,    8,    9,
!            ...etc...
!            191,   -1,
!           ]
!           materialIndex [
!              0,    0,    0,   -1,    1,    1,    1,   -1,    2,    2,    2,
!            ...etc...
!             64,   -1,
!           ]
!         }
!
!         IndexedFaceSet {
!           coordIndex [
!              0,    1,    2,   -1,    3,    4,    5,   -1,    7,    8,    9,
!            ...etc...
!            191,   -1,
!           ]
!           materialIndex [
!              0,    0,    0,   -1,    1,    1,    1,   -1,    2,    2,    2,
!            ...etc...
!             64,   -1,
!           ]
!           normalIndex [
!              0,    0,    0,   -1,    1,    1,    1,   -1,    2,    2,    2,
!            ...etc...
!             64,   -1,
!           ]
!           textureCoordIndex [
!              0,    0,    0,   -1,    1,    1,    1,   -1,    2,    2,    2,
!            ...etc...
!             64,   -1,
!           ]
!         }
!
!         IndexedTriangleStripSet {
!           vertexProperty VertexProperty {
!             vertex [ x y z,
!                      ...
!                      x y z ]
!             normal [ x y z,
!                      ...
!                      x y z ]
!             materialBinding OVERALL
!             normalBinding PER_VERTEX_INDEXED
!           }
!           coordIndex [ i, j, k, l, m, -1, n, o, p, q, r, s, t, u, -1,
!             v, w, x, -1 ..., -1 ]
!           normalIndex -1
!         }
!
!       }
!     }
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, real ( kind = 4 ) COR3_NORMAL(3,COR3_MAX), normals at nodes.  
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer ( kind = 4 ) LINE_MATERIAL(LINE_MAX), the material of each line.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Input, real ( kind = 4 ) MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and A values.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer ( kind = 4 ) TEXTURE_NUM, the number of textures.
!
!    Input, real ( kind = 4 ) TRANSFORM_MATRIX(4,4), the transformation matrix.
!
!    Input, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input, real ( kind = 4 ) VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture
!    coordinates.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) texture_max

  real ( kind = 4 ) cor3(3,cor3_max)
  real ( kind = 4 ) cor3_normal(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  character ( len = * ) fileout_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jtemp
  integer ( kind = 4 ) length
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ) ndx
  integer ( kind = 4 ), parameter :: OFFSET = 1
  character ( len = 200 ) text
  integer ( kind = 4 ) text_num
  character ( len = * ) texture_name(texture_max)
  integer ( kind = 4 ) texture_num
  real ( kind = 4 ) transform_matrix(4,4)
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_tex_uv(2,order_max,face_max)
  character ( len = 20 ) word

  text_num = 0

  write ( iunit, '(a)' ) '#Inventor V2.0 ascii'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) 'Separator {'
  write ( iunit, '(a)' ) '  Info {'
  write ( iunit, '(a)' ) '    string "' // trim ( fileout_name ) &
    // ' generated by IVREAD."'
  write ( iunit, '(a)' ) '    string "Original data in file ' // &
    trim ( filein_name ) // '."'
  write ( iunit, '(a)' ) '  }'
  write ( iunit, '(a)' ) '  Separator {'
  text_num = text_num + 8
!
!  LightModel:
!
!    BASE_COLOR ignores light sources, and uses only diffuse color
!      and transparency.  Even without normal vector information,
!      the object will show up.  However, you won't get shadow
!      and lighting effects.
!
!    PHONG uses the Phong lighting model, accounting for light sources
!      and surface orientation.  This is the default.  I believe
!      you need accurate normal vector information in order for this
!      option to produce nice pictures.
!
!    DEPTH ignores light sources, and calculates lighting based on
!      the location of the object within the near and far planes
!      of the current camera's view volume.
!
  write ( iunit, '(a)' ) '    LightModel {'
  write ( iunit, '(a)' ) '      model PHONG'
  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 3
!
!  Transformation matrix.
!
  call tmat_init ( transform_matrix )

  write ( iunit, '(a)' ) '    MatrixTransform { matrix'
  do i = 1, 4
    write ( iunit, '(8x,4g14.6)' ) transform_matrix(i,1:4)
  end do
  write ( iunit, '(a)' ) '    }'
!
!  Material
!    I think you're only allowed one ambient color.  We take that from
!    material 1.
!
  write ( iunit, '(a)' ) '    Material {'
  write ( iunit, '(a,3f7.3)' ) '      ambientColor ', material_rgba(1:3,1)
  write ( iunit, '(a)' ) '      diffuseColor  ['
  text_num = text_num + 3

  do i = 1, material_num
    write ( iunit, '(8x,3f8.4,'','')' ) material_rgba(1:3,i)
    text_num = text_num + 1
  end do

  write ( iunit, '(a)' ) '      ]'
  write ( iunit, '(a)' ) '      emissiveColor 0.0 0.0 0.0'
  write ( iunit, '(a)' ) '      specularColor 0.0 0.0 0.0'
  write ( iunit, '(a)' ) '      shininess     0.2'
  write ( iunit, '(a)' ) '      transparency  ['
  text_num = text_num + 4

  do ilo = 1, material_num, 10
    ihi = min ( ilo + 9, material_num )
    write ( iunit, '(8x,10(f7.3,'',''))' ) 1.0E+00 - material_rgba(4,ilo:ihi)
    text_num = text_num + 1
  end do

  write ( iunit, '(a)' ) '      ]'
  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 2
!
!  TEXTURE2
!
!  FLAW: We can only handle one texture right now.
!
  if ( 0 < texture_num ) then
    write ( iunit, '(a)' ) '    Texture2 {'
    write ( iunit, '(a)' ) '      filename    "' // trim ( texture_name(1) ) &
      // '"'
    write ( iunit, '(a)' ) '      wrapS       REPEAT'
    write ( iunit, '(a)' ) '      wrapT       REPEAT'
    write ( iunit, '(a)' ) '      model       MODULATE'
    write ( iunit, '(a)' ) '      blendColor  0.0 0.0 0.0'
    write ( iunit, '(a)' ) '    }'
    text_num = text_num + 7
  end if
!
!  TextureCoordinateBinding
!
  write ( iunit, '(a)' ) '    TextureCoordinateBinding {'
  write ( iunit, '(a)' ) '      value PER_VERTEX_INDEXED'
  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 3
!
!  MaterialBinding
!
  write ( iunit, '(a)' ) '    MaterialBinding {'
  write ( iunit, '(a)' ) '      value PER_VERTEX_INDEXED'
  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 3
!
!  NormalBinding
!
!    PER_VERTEX promises that we will write a list of normal vectors
!    in a particular order, namely, the normal vectors for the vertices
!    of the first face, then the second face, and so on.
!
!    PER_VERTEX_INDEXED promises that we will write a list of normal vectors,
!    and then, as part of the IndexedFaceSet, we will give a list of
!    indices referencing this normal vector list.
!
  write ( iunit, '(a)' ) '    NormalBinding {'
  write ( iunit, '(a)' ) '      value PER_VERTEX_INDEXED'
  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 3
!
!  ShapeHints
!
  write ( iunit, '(a)' ) '    ShapeHints {'
  write ( iunit, '(a)' ) '      vertexOrdering COUNTERCLOCKWISE'
  write ( iunit, '(a)' ) '      shapeType UNKNOWN_SHAPE_TYPE'
  write ( iunit, '(a)' ) '      faceType CONVEX'
  write ( iunit, '(a)' ) '      creaseAngle 6.28319'
  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 6
!
!  Point coordinates.
!
  write ( iunit, '(a)' ) '    Coordinate3 {'
  write ( iunit, '(a)' ) '      point ['
  text_num = text_num + 2
 
  do icor3 = 1, cor3_num
    write ( text, '(3f12.4,'','')' ) cor3(1:3,icor3)
    call s_blanks_delete ( text )
    write ( iunit, '(8x,a)' ) trim ( text )
    text_num = text_num + 1
  end do

  write ( iunit, '(a)' ) '      ]'
  write ( iunit, '(a)' ) '    }'
  text_num = text_num + 2
!
!  Texture coordinates.
!
  if ( 0 < texture_num ) then

    write ( iunit, '(a)' ) '    TextureCoordinate2 {'
    write ( iunit, '(a)' ) '      point ['
    text_num = text_num + 2
 
    do iface = 1, face_num
      do ivert = 1, face_order(iface)
        write ( text, '(2f12.4,'','')' ) vertex_tex_uv(1:2,ivert,iface)
        call s_blanks_delete ( text )
        write ( iunit, '(8x,a)' ) trim ( text )
        text_num = text_num + 1
      end do
    end do

    write ( iunit, '(a)' ) '      ]'
    write ( iunit, '(a)' ) '    }'
    text_num = text_num + 2

  end if
!
!  Normal vectors.
!    Use the normal vectors associated with nodes.
!
  if ( 0 < face_num ) then

    write ( iunit, '(a)' ) '    Normal { '
    write ( iunit, '(a)' ) '      vector ['
    text_num = text_num + 2

    do icor3 = 1, cor3_num
      write ( text, '(3f12.4,'','')' ) cor3_normal(1:3,icor3)
      call s_blanks_delete ( text )
      write ( iunit, '(8x,a)' ) trim ( text )
      text_num = text_num + 1
    end do

    write ( iunit, '(a)' ) '      ]'
    write ( iunit, '(a)' ) '    }'
    text_num = text_num + 2

  end if
!
!  IndexedLineSet
!
  if ( 0 < line_num ) then

    write ( iunit, '(a)' ) '    IndexedLineSet {'
!
!  IndexedLineSet coordIndex
!
    write ( iunit, '(a)' ) '      coordIndex ['
    text_num = text_num + 2

    text = ' '
    length = 0

    do j = 1, line_num
 
      write ( word, '(i8,'','')' ) line_dex(j) - OFFSET
      call s_cat ( text, word, text )
      length = length + 1

      if ( line_dex(j)-OFFSET == -1 .or. length >= 10 .or. j >= line_num ) then

        call s_blanks_delete ( text )
        write ( iunit, '(8x,a)' ) trim ( text )
        text_num = text_num + 1
        text = ' '
        length = 0

      end if

    end do
 
    write ( iunit, '(a)' ) '      ]'
    text_num = text_num + 1
!
!  IndexedLineSet materialIndex.
!
    write ( iunit, '(a)' ) '      materialIndex ['
    text_num = text_num + 1

    text = ' '
    length = 0

    do j = 1, line_num
 
      write ( word, '(i8,'','')' ) line_material(j) - OFFSET

      call s_cat ( text, word, text )
      length = length + 1

      if ( line_material(j)-OFFSET == -1 .or. &
           length >= 10 .or. j >= line_num ) then

        call s_blanks_delete ( text )
        write ( iunit, '(8x,a)' ) trim ( text )
        text_num = text_num + 1
        text = ' '
        length = 0

      end if
 
    end do
 
    write ( iunit, '(a)' ) '      ]'
    write ( iunit, '(a)' ) '    }'
    text_num = text_num + 2

  end if
!
!  IndexedFaceSet.
!
  if ( 0 < face_num ) then

    write ( iunit, '(a)' ) '    IndexedFaceSet {'
!
!  IndexedFaceSet coordIndex
!
    write ( iunit, '(a)' ) '      coordIndex ['
    text_num = text_num + 2

    text = ' '
    length = 0
 
    do iface = 1, face_num

      do ivert = 1, face_order(iface) + 1

        if ( ivert <= face_order(iface) ) then
          itemp = face(ivert,iface) - OFFSET
        else
          itemp = 0 - OFFSET
        end if

        write ( word, '(i8,'','')' ) itemp
  
        call s_cat ( text, word, text )
        length = length + 1

        if ( itemp == -1 .or. length >= 10 .or. &
           ( iface == face_num .and. ivert == face_order(iface) + 1 )  ) then

          call s_blanks_delete ( text )
          write ( iunit, '(8x,a)' ) trim ( text )
          text_num = text_num + 1
          text = ' '
          length = 0

        end if

      end do

    end do

    write ( iunit, '(a)' ) '      ]'
    text_num = text_num + 1
!
!  IndexedFaceSet normalIndex
!
    if ( 0 < texture_num ) then

      write ( iunit, '(a)' ) '      normalIndex ['
      text_num = text_num + 1

      text = ' '
      length = 0

      do iface = 1, face_num
  
        do ivert = 1, face_order(iface) + 1

          if ( ivert <= face_order(iface) ) then
            itemp = face(ivert,iface) - OFFSET
          else
            itemp = 0 - OFFSET
          end if
  
          write ( word, '(i8,'','')' ) itemp

          call s_cat ( text, word, text )
          length = length + 1

          if ( itemp == -1 .or. length >= 10 .or. &
             ( iface == face_num .and. ivert == face_order(iface) + 1 )  ) then

            call s_blanks_delete ( text )
            write ( iunit, '(8x,a)' ) trim ( text )
            text_num = text_num + 1
            text = ' '
            length = 0

          end if

        end do

      end do

      write ( iunit, '(a)' ) '      ]'
      text_num = text_num + 1

    end if
!
!  IndexedFaceSet textureCoordIndex
!
    write ( iunit, '(a)' ) '      textureCoordIndex ['
    text_num = text_num + 1

    text = ' '
    length = 0
    itemp = 0

    do iface = 1, face_num

      do ivert = 1, face_order(iface) + 1

        if ( ivert <= face_order(iface) ) then
          jtemp = itemp
          itemp = itemp + 1
        else
          jtemp = - 1
        end if
  
        write ( word, '(i8,'','')' ) jtemp

        call s_cat ( text, word, text )
        length = length + 1

        if ( jtemp == -1 .or. length >= 10 .or. &
           ( iface == face_num .and. ivert == face_order(iface) + 1 )  ) then

          call s_blanks_delete ( text )
          write ( iunit, '(8x,a)' ) trim ( text )
          text_num = text_num + 1
          text = ' '
          length = 0

        end if

      end do

    end do

    write ( iunit, '(a)' ) '      ]'
    text_num = text_num + 1
!
!  IndexedFaceSet materialIndex
!
    write ( iunit, '(a)' ) '      materialIndex ['
    text_num = text_num + 1

    text = ' '
    length = 0
    ndx = 0
 
    do iface = 1, face_num
 
      do ivert = 1, face_order(iface) + 1

        if ( ivert <= face_order(iface) ) then
          itemp = vertex_material(ivert,iface) - OFFSET
          ndx = ndx + 1
        else
          itemp = 0 - OFFSET
        end if

        write ( word, '(i8,'','')' ) itemp

        call s_cat ( text, word, text )
        length = length + 1

        if ( itemp == -1 .or. length >= 10 .or. &
           ( iface == face_num .and. ivert == face_order(iface) + 1 )  ) then

          call s_blanks_delete ( text )
          write ( iunit, '(8x,a)' ) trim ( text )
          text_num = text_num + 1
          text = ' '
          length = 0
  
        end if
 
      end do

    end do

    write ( iunit, '(a)' ) '      ]'
    text_num = text_num + 1

    write ( iunit, '(a)' ) '    }'
    text_num = text_num + 1

  end if
!
!  Close up the Separator node.
!
  write ( iunit, '(a)' ) '  }'
!
!  Close up the Separator node.
!
  write ( iunit, '(a)' ) '}'
 
  text_num = text_num + 2
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'IV_WRITE - Wrote ', text_num, &
    ' text lines to ' // trim ( fileout_name )


  return
end
subroutine i4vec_max ( nval, iarray, imax )

!*****************************************************************************80
!
!! I4VEC_MAX computes the maximum element of an integer array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVAL, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) IARRAY(NVAL), the array.
!
!    Output, integer ( kind = 4 ) IMAX, the value of the largest entry in the array.
!
  implicit none

  integer ( kind = 4 ) nval

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarray(nval)
  integer ( kind = 4 ) imax

  if ( nval <= 0 ) then

    imax = 0

  else

    imax = iarray(1)

    do i = 2, nval
      imax = max ( imax, iarray(i) )
    end do

  end if
 
  return
end
subroutine i4vec_reverse ( n, x )

!*****************************************************************************80
!
!! I4VEC_REVERSE reverses the elements of an I4VEC.
!
!  Example:
!
!    Input:
!
!      N = 5, X = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      X = ( 15, 14, 13, 12, 11 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) X(N), the array to be reversed.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) temp
  integer ( kind = 4 ) x(n)

  do i = 1, n/2
    temp = x(i)
    x(i) = x(n+1-i)
    x(n+1-i) = temp
  end do

  return
end
subroutine i4vec_rotate ( n, m, a )

!*****************************************************************************80
!
!! I4VEC_ROTATE rotates an object in place.
!
!  Example:
!
!    Input:
!
!      N = 5, M = 2
!      A = ( 1, 2, 3, 4, 5 )
!
!    Output:
!
!      A = ( 4, 5, 1, 2, 3 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) M, the number of positions to the right that
!    each element should be moved.  Elements that shift pass position
!    N "wrap around" to the beginning of the array.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be rotated.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mcopy
  integer ( kind = 4 ) nset
  integer ( kind = 4 ) temp
!
!  Force M to be positive, between 0 and N-1.
!
  mcopy = i4_modp ( m, n )

  if ( mcopy == 0 ) then
    return
  end if

  istart = 0
  nset = 0

  do

    istart = istart + 1

    if ( n < istart ) then
      exit
    end if

    temp = a(istart)
    iget = istart
!
!  Copy the new value into the vacated entry.
!
    do

      iput = iget

      iget = iget - mcopy

      if ( iget < 1 ) then
        iget = iget + n
      end if

      if ( iget == istart ) then
        exit
      end if

      a(iput) = a(iget)
      nset = nset + 1

    end do

    a(iput) = temp
    nset = nset + 1

    if ( nset >= n ) then
      exit
    end if

  end do

  return
end
function lcon ( chr )

!*****************************************************************************80
!
!! LCON reports whether a character is a control character or not.
!
!  Discussion:
!
!    A "control character" has ASCII code <= 31 or ASCII code => 127.
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
!    Input, character CHR, is the character to be tested.
!
!    Output, logical LCON, TRUE if CHR is a control character, and
!    FALSE otherwise.
!
  implicit none

  character chr
  logical lcon

  if ( ichar ( chr ) <= 31 .or. ichar ( chr ) >= 127 ) then
    lcon = .true.
  else
    lcon = .false.
  end if

  return
end
subroutine mesh_t3 ( face, face_max, face_num, face_order, m, n, order_max )

!*****************************************************************************80
!
!! MESH_T3 produces a grid of pairs of 3 node triangles.
!
!  Example:
!
!    Input:
!
!      M = 4, N = 3
!
!    Output:
!
!      FACE =
!         1,  4,  2;
!         5,  2,  4;
!         2,  5,  3;
!         6,  3,  5;
!         4,  7,  5;
!         8,  5,  7
!         5,  8,  6;
!         9,  6,  8;
!         7, 10,  8;
!        11,  8, 10;
!         8, 11,  9;
!        12,  9, 11.
!
!  Diagram:
!
!    3----6----9---12
!    |\ 4 |\ 8 |\12 |
!    | \  | \  | \  |
!    |  \ |  \ |  \ |
!    |  3\|  7\| 11\|
!    2----5----8---11
!    |\ 2 |\ 6 |\10 |
!    | \  | \  | \  |
!    |  \ |  \ |  \ |
!    |  1\|  5\|  9\|
!    1----4----7---10
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX); FACE(I,J) contains the
!    index of the I-th vertex of the J-th face.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input/output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the order of each face.
!
!    Input, integer ( kind = 4 ) M, N, the number of nodes that form the mesh.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum order for any face.
!
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  if ( face_max < face_num + 2 * ( m - 1 ) * ( n - 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH_T3 - Fatal error!'
    write ( *, '(a)' ) '  Not enough storage for new faces.'
    write ( *, '(a,i6)' ) '  Current number of faces is ', face_num
    write ( *, '(a,i6)' ) '  Maximum number of faces is ', face_max
    write ( *, '(a,i6)' ) '  New faces needed: ', 2 * ( m - 1 ) * ( n - 1 )
    stop
  end if

  do j = 1, n - 1
    do i = 1, m - 1

      face_num = face_num + 1
      face_order(face_num) = 3

      face(1,face_num) = ( i - 1 ) * n + j
      face(2,face_num) =   i       * n + j
      face(3,face_num) = ( i - 1 ) * n + j + 1

      face_num = face_num + 1
      face_order(face_num) = 3

      face(1,face_num) =   i       * n + j + 1
      face(2,face_num) = ( i - 1 ) * n + j + 1
      face(3,face_num) =   i       * n + j

    end do
  end do

  return
end
subroutine news

!*****************************************************************************80
!
!! NEWS prints out news (old and new) about the program.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 June 2002
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NEWS:'
  write ( *, '(a)' ) '  This is a list of changes to the program:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  27 August 2003:'
  write ( *, '(a)' ) '    Added OFF_READ and OFF_WRITE.'
  write ( *, '(a)' ) '  13 June 2002:'
  write ( *, '(a)' ) '    Modifying IV_POINT_WRITE to add line data.'
  write ( *, '(a)' ) '    My goal is to display Voronoi diagrams on a sphere.'
  write ( *, '(a)' ) '  04 June 2002:'
  write ( *, '(a)' ) '    Added XYZ_READ.'
  write ( *, '(a)' ) '    Added IV_POINT_WRITE.'
  write ( *, '(a)' ) '  23 February 2002:'
  write ( *, '(a)' ) '    Modified the PostScript header.'
  write ( *, '(a)' ) '  13 October 2001:'
  write ( *, '(a)' ) '    Added TS_WRITE.'
  write ( *, '(a)' ) '  05 June 2001:'
  write ( *, '(a)' ) '    Restored XYZ_WRITE.'
  write ( *, '(a)' ) '  14 June 2000:'
  write ( *, '(a)' ) '    FORTRAN90 conversion.'
  write ( *, '(a)' ) '  06 October 1999:'
  write ( *, '(a)' ) '    Trying to retrieve PS_WRITE from old code.'
  write ( *, '(a)' ) '  03 August 1999:'
  write ( *, '(a)' ) '    Corrected TMAT_ROT_VECTOR.'
  write ( *, '(a)' ) '  29 June 1999:'
  write ( *, '(a)' ) '    IV_READ and IV_WRITE support UV textures.'
  write ( *, '(a)' ) '  24 June 1999:'
  write ( *, '(a)' ) '    HRC_WRITE and TXT_WRITE write texture info.'
  write ( *, '(a)' ) '    BYU_WRITE writes Movie.BYU format.'
  write ( *, '(a)' ) '  23 June 1999:'
  write ( *, '(a)' ) '    HRC_READ and HRC_WRITE use vertex textures.'
  write ( *, '(a)' ) '  08 June 1999:'
  write ( *, '(a)' ) '    Added simple TECPLOT output.'
  write ( *, '(a)' ) '    Added Greg Hood "TRI" triangle output.'
  write ( *, '(a)' ) '  02 June 1999:'
  write ( *, '(a)' ) '    Adding material names.'
  write ( *, '(a)' ) '  26 May 1999:'
  write ( *, '(a)' ) '    Internal LINE_PRUNE option added.'
  write ( *, '(a)' ) '  24 May 1999:'
  write ( *, '(a)' ) '    STLA_WRITE automatically decomposes any'
  write ( *, '(a)' ) '    non-triangular faces before writing them.'
  write ( *, '(a)' ) '  23 May 1999:'
  write ( *, '(a)' ) '    DXF_WRITE can output faces.'
  write ( *, '(a)' ) '  22 May 1999:'
  write ( *, '(a)' ) '    VLA output includes line versions of faces.'
  write ( *, '(a)' ) '  04 May 1999:'
  write ( *, '(a)' ) '    SMF_READ/WRITE support face/node material,'
  write ( *, '(a)' ) '    normal, texture coordinates.'
  write ( *, '(a)' ) '  01 May 1999:'
  write ( *, '(a)' ) '    Relaxation smoothing option added.'
  write ( *, '(a)' ) '  27 April 1999:'
  write ( *, '(a)' ) '    SMF_READ and SMF_WRITE handle SMF2.0 colors.'
  write ( *, '(a)' ) '  23 April 1999:'
  write ( *, '(a)' ) '    Better POV_WRITE material treatment.'
  write ( *, '(a)' ) '    FACE_MATERIAL vs VERTEX_MAT fixup.'
  write ( *, '(a)' ) '  21 April 1999:'
  write ( *, '(a)' ) '    Trying to get IV_READ to merge two files.'
  write ( *, '(a)' ) '  20 April 1999:'
  write ( *, '(a)' ) '    Added << option, trying to get OBJ_READ'
  write ( *, '(a)' ) '      to merge two files.'
  write ( *, '(a)' ) '  26 March 1999'
  write ( *, '(a)' ) '    Added RGB_TO_HUE routine;'
  write ( *, '(a)' ) '    Adding UCD_WRITE.'
  write ( *, '(a)' ) '  23 February 1999'
  write ( *, '(a)' ) '    In HRC_WRITE, made specular and ambient'
  write ( *, '(a)' ) '    material colors the same as diffuse.'
  write ( *, '(a)' ) '  15 February 1999'
  write ( *, '(a)' ) '    Trying to get grid lines in OOGL.'
  write ( *, '(a)' ) '  13 February 1999'
  write ( *, '(a)' ) '    Added face area calculation.'
  write ( *, '(a)' ) '  12 February 1999'
  write ( *, '(a)' ) '    Added new color scheme to IV_WRITE.'
  write ( *, '(a)' ) '  10 February 1999'
  write ( *, '(a)' ) '    HRC_READ should now be able to read '
  write ( *, '(a)' ) '    material data.'
  write ( *, '(a)' ) '  09 February 1999'
  write ( *, '(a)' ) '    HRC_WRITE now writes out material data.'
  write ( *, '(a)' ) '    Moving all RGB material information'
  write ( *, '(a)' ) '      into MATERIAL_RGBA, with other items'
  write ( *, '(a)' ) '      using pointers to index it.'
  write ( *, '(a)' ) '  08 February 1999'
  write ( *, '(a)' ) '    Adding "OOGL" file format for Greg Foss.'
  write ( *, '(a)' ) '  02 December 1998'
  write ( *, '(a)' ) '    Restored VRML write.'
  write ( *, '(a)' ) '    Set up simple hooks for texture map names.'
  write ( *, '(a)' ) '  19 November 1998'
  write ( *, '(a)' ) '    IV_WRITE uses PER_VERTEX normal binding.'
  write ( *, '(a)' ) '  18 November 1998'
  write ( *, '(a)' ) '    Added node-based normal vectors.'
  write ( *, '(a)' ) '  17 November 1998'
  write ( *, '(a)' ) '    Added face node ordering reversal option.'
  write ( *, '(a)' ) '  23 October 1998'
  write ( *, '(a)' ) '    Added polygon to triangle option.'
  write ( *, '(a)' ) '  20 October 1998'
  write ( *, '(a)' ) '    Added interactive scaling patch.'
  write ( *, '(a)' ) '    Inserted TMAT routines.'
  write ( *, '(a)' ) '  19 October 1998'
  write ( *, '(a)' ) '    SMF_READ and SMF_WRITE added.'
  write ( *, '(a)' ) '  12 October 1998'
  write ( *, '(a)' ) '    Added FACE_CHECK code.'
  write ( *, '(a)' ) '  08 October 1998'
  write ( *, '(a)' ) '    Added POV_WRITE;'
  write ( *, '(a)' ) '    Added SET_VERTEX_NORMAL;'
  write ( *, '(a)' ) '    Modified normal vector computation.'
  write ( *, '(a)' ) '  30 August 1998'
  write ( *, '(a)' ) '    Still trying to fix up normals, because'
  write ( *, '(a)' ) '    of OBJ_READ and OBJ_WRITE complications.'
  write ( *, '(a)' ) '  29 August 1998'
  write ( *, '(a)' ) '    OBJ_READ and OBJ_WRITE now handle normals,'
  write ( *, '(a)' ) '    and read and write normals to file.'
  write ( *, '(a)' ) '    OBJ_READ can handle // face format.'
  write ( *, '(a)' ) '  28 August 1998'
  write ( *, '(a)' ) '    STLA_READ and STLA_WRITE seem OK after'
  write ( *, '(a)' ) '    the normal changes.'
  write ( *, '(a)' ) '  27 August 1998'
  write ( *, '(a)' ) '    Trying better NORMAL storage approach.'
  write ( *, '(a)' ) '  21 August 1998'
  write ( *, '(a)' ) '    Trying to add HRC_READ.'
  write ( *, '(a)' ) '    TXT_WRITE improved.'
  write ( *, '(a)' ) '  20 August 1998'
  write ( *, '(a)' ) '    Adding linear splines to HRC_WRITE.'
  write ( *, '(a)' ) '  19 August 1998'
  write ( *, '(a)' ) '    Automatic normal computation for OBJ files.'
  write ( *, '(a)' ) '    SoftImage HRC output added.'
  write ( *, '(a)' ) '    Normal vector computation improved.'
  write ( *, '(a)' ) '  18 August 1998'
  write ( *, '(a)' ) '    Improved treatment of face materials and normals.'
  write ( *, '(a)' ) '  17 August 1998'
  write ( *, '(a)' ) '    The maximum number of vertices per face'
  write ( *, '(a)' ) '    was increased to 35.'
  write ( *, '(a)' ) '    The maximum input line length was increased'
  write ( *, '(a)' ) '    to 256 characters.'
  write ( *, '(a)' ) '  10 August 1998:'
  write ( *, '(a)' ) '    Output DXF files have a comment now.'
  write ( *, '(a)' ) '    OBJ_READ corrected line indexing problem.'
  write ( *, '(a)' ) '  24 July 1998:'
  write ( *, '(a)' ) '    INCHECK checks the input data.'
  write ( *, '(a)' ) '    DXF_READ suppresses duplicate points.'
  write ( *, '(a)' ) '    Removed grid routines.'
  write ( *, '(a)' ) '    LINES(2,*) -> LINE_DEX(), LINE_MATERIAL().'
  write ( *, '(a)' ) '    PS and VRML output dropped.'
  write ( *, '(a)' ) '    OBJ_WRITE line output tightened up.'
  write ( *, '(a)' ) '  22 July 1998:'
  write ( *, '(a)' ) '    STLA_READ suppresses duplicate points.'
  write ( *, '(a)' ) '  21 July 1998:'
  write ( *, '(a)' ) '    VLA_READ suppresses duplicate points.'
  write ( *, '(a)' ) '    OBJ_WRITE outputs line data now.'
  write ( *, '(a)' ) '  15 July 1998:'
  write ( *, '(a)' ) '    Added STLA_READ and STLA_WRITE.'
  write ( *, '(a)' ) '  11 July 1998:'
  write ( *, '(a)' ) '    DXF_READ and DXF_WRITE use IV data.'
  write ( *, '(a)' ) '    Dropped XYZ data structures.'
  write ( *, '(a)' ) '  10 July 1998:'
  write ( *, '(a)' ) '    Dropped XYZ input/output option.'
  write ( *, '(a)' ) '    VLA_READ and VLA_WRITE use IV data.'
  write ( *, '(a)' ) '  08 July 1998:'
  write ( *, '(a)' ) '    Added OBJ_READ and OBJ_WRITE.'
  write ( *, '(a)' ) '    Set ORDER_MAX=4, to allow for quad faces.'
  write ( *, '(a)' ) '  05 July 1998:'
  write ( *, '(a)' ) '    Added RF command to reverse faces.'
  write ( *, '(a)' ) '    Fixed 0/1 index based problem for FACE.'
  write ( *, '(a)' ) '  04 July 1998:'
  write ( *, '(a)' ) '    Added CHECK command to examine a face.'
  write ( *, '(a)' ) '  03 July 1998:'
  write ( *, '(a)' ) '    Only converting data when necessary.'
  write ( *, '(a)' ) '    Alternate IV triangles have opposite sense.'
  write ( *, '(a)' ) '    NORMALS command will recompute normals.'
  write ( *, '(a)' ) '  02 July 1998:'
  write ( *, '(a)' ) '    Trying to write simple ASE files.'
  write ( *, '(a)' ) '  01 July 1998:'
  write ( *, '(a)' ) '    Tentative attempts to read new IV data,'
  write ( *, '(a)' ) '    INDEXEDTRIANGLESTRIP.'
  write ( *, '(a)' ) '  03 June 1998:'
  write ( *, '(a)' ) '    MATERIALINDEX works for IV faces AND lines.'
  write ( *, '(a)' ) '  02 June 1998:'
  write ( *, '(a)' ) '    Trying to sort out -1/0 LINES confusion.'
  write ( *, '(a)' ) '    VRML_WRITE uses Inventor data.'
  write ( *, '(a)' ) '    I got VRML_WRITE to do color lines.'
  write ( *, '(a)' ) '    I need to reconcile RGBCOLOR/FACE/NODE.'
  write ( *, '(a)' ) '  15 May 1998:'
  write ( *, '(a)' ) '    Set up RGBCOLOR for new color handling.'
  write ( *, '(a)' ) '  08 May 1998:'
  write ( *, '(a)' ) '    Preparing for IV PER_VERTEX_INDEXED.'
  write ( *, '(a)' ) '  06 May 1998:'
  write ( *, '(a)' ) '    Added "reverse normal" option.'
  write ( *, '(a)' ) '    Added 3D projected plane grid.'
  write ( *, '(a)' ) '  05 May 1998:'
  write ( *, '(a)' ) '    Added projection into 3D plane.'
  write ( *, '(a)' ) '    Added spherical grid lines.'
  write ( *, '(a)' ) '  04 May 1998:'
  write ( *, '(a)' ) '    Sphere projection set up.'
  write ( *, '(a)' ) '    The PostScript output seems OK.'
  write ( *, '(a)' ) '  30 April 1998:'
  write ( *, '(a)' ) '    Adding 2D PostScript output option.'
  write ( *, '(a)' ) '  23 April 1998:'
  write ( *, '(a)' ) '    ASE->IV surface information sorta works.'
  write ( *, '(a)' ) '  22 April 1998:'
  write ( *, '(a)' ) '    IVREAD accepts command line arguments.'
  write ( *, '(a)' ) '  21 April 1998:'
  write ( *, '(a)' ) '    IV_WRITE now writes a default material.'
  write ( *, '(a)' ) '  20 April 1998:'
  write ( *, '(a)' ) '    ASE_READ tries to read vertex color.'
  write ( *, '(a)' ) '  17 April 1998:'
  write ( *, '(a)' ) '    ASE_READ reads transform matrix.'
  write ( *, '(a)' ) '    Overhauled ReadIV routine.'
  write ( *, '(a)' ) '  16 April 1998:'
  write ( *, '(a)' ) '    Increased big array limits.'
  write ( *, '(a)' ) '    Adding ASE_READ.'
  write ( *, '(a)' ) '  15 April 1998:'
  write ( *, '(a)' ) '    VLA intensities parameterized in INTENSE.'
  write ( *, '(a)' ) '    Got VRML option to work.'
  write ( *, '(a)' ) '  14 April 1998:'
  write ( *, '(a)' ) '    Added VRML_WRITE.'
  write ( *, '(a)' ) '  13 April 1998:'
  write ( *, '(a)' ) '    Made program command driven.'
  write ( *, '(a)' ) '    Started projection option.'
  write ( *, '(a)' ) '  10 April 1998:'
  write ( *, '(a)' ) '    Added Min/Max coordinate print.'
  write ( *, '(a)' ) '    Compressed IV output.'

  return
end
subroutine node_relax ( cor3, cor3_max, cor3_new, cor3_num, face, face_max, &
  face_num, face_order, order_max )

!*****************************************************************************80
!
!! NODE_RELAX smooths a shape by an averaging operation on the node positions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Output, real ( kind = 4 ) COR3_NEW(3,COR3_MAX), the new, averaged coordinates of 
!    the nodes.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces. 
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), is the number of nodes
!    making up each face.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, is the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!  Local variables:
!
!    Integer COR3_NUMBER(COR3_MAX), the number of node neighbors of node I.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3(3,cor3_max)
  real ( kind = 4 ) cor3_new(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) cor3_number(cor3_max)
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) inode
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) jnode
!
!  COR3_NEW will contain the new averaged coordinates.
!
  cor3_number(1:cor3_num) = 0
  cor3_new(1:3,1:cor3_num) = 0.0E+00
!
!  Consider each edge.  Essentially, the edge (I,J) is a signal to
!  add the old coordinates of I to the new J coordinates, and vice versa.
!
!  Because we are using a face representation, many, perhaps all the
!  edges, will show up repeatedly, probably twice.  To keep the algorithm
!  simple, for now we will simply use an edge every time it shows up
!  in a face, which means that edges that occur in multiple faces
!  will be weighted more.
!
  do iface = 1, face_num

    inode = face(face_order(iface),iface)

    do ivert = 1, face_order(iface)

      jnode = inode
      inode = face(ivert,iface)

      cor3_number(inode) = cor3_number(inode) + 1
      cor3_number(jnode) = cor3_number(jnode) + 1

      cor3_new(1:3,jnode) = cor3_new(1:3,jnode) + cor3(1:3,inode)
      cor3_new(1:3,inode) = cor3_new(1:3,inode) + cor3(1:3,jnode)

    end do

  end do
!
!  Copy the new into the old.
!
  do icor3 = 1, cor3_num

    if ( cor3_number(icor3) /= 0 ) then
      cor3_new(1:3,icor3) = cor3_new(1:3,icor3) / real ( cor3_number(icor3) )
    end if

  end do

  return
end
subroutine node_to_vertex_material ( cor3_material, cor3_max, face, &
  face_max, face_num, face_order, order_max, vertex_material )

!*****************************************************************************80
!
!! NODE_TO_VERTEX_MATERIAL extends node material definitions to vertices.
!
!  Discussion:
!
!    A NODE is a point in space.
!    A VERTEX is a node as used in a particular face.
!    One node may be used as a vertex in several faces, or none.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) COR3_MATERIAL(COR3_MAX), the material index of each node.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input/output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Output, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  integer ( kind = 4 ) cor3_material(cor3_max)
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) node
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) vertex_material(order_max,face_max)

  do iface = 1, face_num
    do ivert = 1, face_order(iface)
      node = face(ivert,iface)
      vertex_material(ivert,iface) = cor3_material(node)
    end do
  end do

  return
end
subroutine obj_read ( bad_num, cor3, cor3_max, cor3_num, face, &
  face_material, face_max, face_num, face_order, filein_name, &
  group_num, ierror, iunit, line_dex, line_material, line_max, &
  line_num, material_max, material_name, material_num, material_rgba, &
  object_num, order_max, text_num, vertex_material, vertex_normal )

!*****************************************************************************80
!
!! OBJ_READ reads graphics information from a Wavefront OBJ file.
!
!  Discussion:
!
!    It is intended that the information read from the file can
!    either start a whole new graphics object, or simply be added
!    to a current graphics object via the '<<' command.  
!
!    This is controlled by whether the input values have been zeroed 
!    out or not.  This routine simply tacks on the information it 
!    finds to the current graphics object.
!
!  Example:
!
!    #  magnolia.obj
!
!    mtllib ./vp.mtl
!
!    g
!    v -3.269770 -39.572201 0.876128
!    v -3.263720 -39.507999 2.160890
!    ...
!    v 0.000000 -9.988540 0.000000
!    vn 1.0 0.0 0.0E+00
!    ...
!    vn 0.0 1.0 0.0E+00
!    g stem
!    s 1
!    usemtl brownskn
!    f 8 9 11 10
!    f 12 13 15 14
!    ...
!    f 788 806 774
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BAD_NUM, the number of bad text lines encountered.
!
!    Input/output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input/output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input/output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices 
!    per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Input/output, integer ( kind = 4 ) GROUP_NUM, the number of groups.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which data is read.
!
!    Input/output, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated 
!    by -1.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input/output, integer ( kind = 4 ) LINE_MATERIAL(LINE_MAX), material index for 
!    each line.
!
!    Input/output, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Input/output, character ( len = * ) MATERIAL_NAME(MATERIAL_MAX), 
!    material names.
!
!    Input/output, integer ( kind = 4 ) OBJECT_NUM, the number of objects.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Output, integer ( kind = 4 ) TEXT_NUM, the number of lines of text read from
!    the file.
!
!    Input/output, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Output, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max

  integer ( kind = 4 ) bad_num
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) cor3_num_base
  logical done
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  integer ( kind = 4 ) group_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) iword
  integer ( kind = 4 ) lchar
  character ( len = 256 ) line
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_num
  character ( len = * ) material_name(material_max)
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  real ( kind = 4 ) normal_temp(3,order_max*face_max)
  integer ( kind = 4 ) object_num
  integer ( kind = 4 ), parameter :: OFFSET = 1
  real ( kind = 4 ) r
  real ( kind = 4 ) s
  logical s_eqi
  integer ( kind = 4 ) text_num
  real ( kind = 4 ) temp
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  integer ( kind = 4 ) vertex_normal_num
  character ( len = 256 ) word
  character ( len = 256 ) word1

  ierror = 0
!
!  Save a copy of the input value of COR3_NUM to use as a base.
!
  cor3_num_base = cor3_num

  bad_num = 0
  text_num = 0

  vertex_normal_num = 0
  word = ' '
!
!  Read a line of text from the file.
!
  do
 
    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    text_num = text_num + 1
!
!  Replace any control characters (in particular, TAB's) by blanks.
!
    call s_control_blank ( line )
 
    done = .true.
    iword = 0
!
!  Read a word from the line.
!
    call word_next_read ( line, word, done )
!
!  If no more words in this line, read a new line.
!
    if ( done ) then
      cycle
    end if
!
!  If this word begins with '#' or '$', then it's a comment.  Read a new line.
!
    if ( word(1:1) == '#' .or. word(1:1) == '$' ) then
      cycle
    end if

    iword = iword + 1
    if ( iword == 1 ) then
      word1 = word
    end if
!
!  BEVEL
!  Bevel interpolation.
!
    if ( s_eqi ( word1, 'BEVEL' ) ) then
!  
!  BMAT
!  Basis matrix.
!
    else if ( s_eqi ( word1, 'BMAT' ) ) then
!
!  C_INTERP
!  Color interpolation.
!
    else if ( s_eqi ( word1, 'C_INTERP' ) ) then
!
!  CON
!  Connectivity between free form surfaces.
!
    else if ( s_eqi ( word1, 'CON' ) ) then
!
!  CSTYPE
!  Curve or surface type.
!
    else if ( s_eqi ( word1, 'CSTYPE' ) ) then
!
!  CTECH
!  Curve approximation technique.
!
    else if ( s_eqi ( word1, 'CTECH' ) ) then
!
!  CURV
!  Curve.
!
    else if ( s_eqi ( word1, 'CURV' ) ) then
!
!  CURV2
!  2D curve.
!
    else if ( s_eqi ( word1, 'CURV2' ) ) then
!
!  D_INTERP
!  Dissolve interpolation.
!
    else if ( s_eqi ( word1, 'D_INTERP' ) ) then
!
!  DEG
!  Degree.
!
    else if ( s_eqi ( word1, 'DEG' ) ) then
!
!  END
!  End statement.
!
    else if ( s_eqi ( word1, 'END' ) ) then
!
!  F V1 V2 V3 ...
!    or
!  F V1/VT1/VN1 V2/VT2/VN2 ...
!    or
!  F V1//VN1 V2//VN2 ...
!
!  Face.
!  A face is defined by the vertices.
!  Optionally, slashes may be used to include the texture vertex
!  and vertex normal indices.
!
    else if ( s_eqi ( word1, 'F' ) ) then

      face_num = face_num + 1

      ivert = 0

      do

        ivert = ivert + 1
 
        call word_next_read ( line, word, done )

        if ( done ) then
          exit
        end if
!
!  Locate the slash characters in the word, if any.
!
        i1 = index ( word, '/' )
        if ( 0 < i1 ) then
          i2 = index ( word(i1+1:), '/' ) + i1
        else
          i2 = 0
        end if
!
!  Read the vertex index.
!
        call s_to_i4 ( word, itemp, ierror, lchar )

        if ( ierror /= 0 ) then
          itemp = -1
          ierror = 0
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'OBJ_READ - Error!'
          write ( *, '(a)' ) '  Bad FACE field.'
          write ( *, '(a)' ) trim ( word )
        end if

        if ( ivert <= order_max .and. face_num <= face_max ) then
          face(ivert,face_num) = itemp + cor3_num_base
          vertex_material(ivert,face_num) = material_num
        end if

        if ( face_num <= face_max ) then

          face_material(face_num) = material_num
          face_order(face_num) = ivert

        end if
!
!  If there are two slashes, then read the data following the second one.
!
        if ( 0 < i2 ) then

          call s_to_i4 ( word(i2+1:), itemp, ierror, lchar )

          if ( 1 <= itemp .and. itemp <= vertex_normal_num ) then
            vertex_normal(1:3,ivert,face_num) = normal_temp(1:3,itemp)
          end if

        end if

      end do
!
!  G
!  Group name.
!
    else if ( s_eqi ( word1, 'G' ) ) then

      group_num = group_num + 1
!
!  HOLE
!  Inner trimming loop.
!
    else if ( s_eqi ( word1, 'HOLE' ) ) then
!
!  L
!  A line, described by a sequence of vertex indices.
!  Are the vertex indices 0 based or 1 based?
!
    else if ( s_eqi ( word1, 'L' ) ) then

      do
 
      call word_next_read ( line, word, done )
!
!  If no more indices, tack a "-1" on the end.
!
      if ( done ) then

        line_num = line_num + 1

        if ( line_num <= line_max ) then
          line_dex(line_num) = -1 + OFFSET
          line_material(line_num) = -1 + OFFSET
        end if

        exit

      end if
!
!  Otherwise, extract the node index and add it to the line list.
!
      call s_to_i4 ( word, itemp, ierror, lchar )

      line_num = line_num + 1

      if ( line_num <= line_max ) then
        line_dex(line_num) = itemp
        line_material(line_num) = material_num
      end if

    end do
!
!  LOD
!  Level of detail.
!
    else if ( s_eqi ( word1, 'LOD' ) ) then
!
!  MG
!  Merging group.
!
    else if ( s_eqi ( word1, 'MG' ) ) then
!
!  MTLLIB
!  Material library.
!
    else if ( s_eqi ( word1, 'MTLLIB' ) ) then
!
!  O
!  Object name.
!
    else if ( s_eqi ( word1, 'O' ) ) then

      object_num = object_num + 1
!
!  P
!  Point.
!
    else if ( s_eqi ( word1, 'P' ) ) then
!
!  PARM
!  Parameter values.
!
    else if ( s_eqi ( word1, 'PARM' ) ) then
!
!  S
!  Smoothing group.
!
    else if ( s_eqi ( word1, 'S' ) ) then
!
!  SCRV
!  Special curve.
!
    else if ( s_eqi ( word1, 'SCRV' ) ) then
!
!  SHADOW_OBJ
!  Shadow casting.
!
    else if ( s_eqi ( word1, 'SHADOW_OBJ' ) ) then
!
!  SP
!  Special point.
!
    else if ( s_eqi ( word1, 'SP' ) ) then
!
!  STECH
!  Surface approximation technique.
!
    else if ( s_eqi ( word1, 'STECH' ) ) then
!
!  STEP
!  Stepsize.
!
    else if ( s_eqi ( word1, 'STEP' ) ) then
!
!  SURF
!  Surface.
!
    else if ( s_eqi ( word1, 'SURF' ) ) then
!
!  TRACE_OBJ
!  Ray tracing.
!
    else if ( s_eqi ( word1, 'TRACE_OBJ' ) ) then
!
!  TRIM
!  Outer trimming loop.
!
    else if ( s_eqi ( word1, 'TRIM' ) ) then
!
!  USEMTL
!  Material name.
!
    else if ( s_eqi ( word1, 'USEMTL' ) ) then

      call word_next_read ( line, word, done )

      material_num = material_num + 1

      if ( material_num <= material_max ) then
        material_name(material_num) = word
        call r4_random ( 0.0E+00, 1.0E+00, r )
        material_rgba(1,material_num) = r
        call r4_random ( 0.0E+00, 1.0E+00, r )
        material_rgba(2,material_num) = r
        call r4_random ( 0.0E+00, 1.0E+00, r )
        material_rgba(3,material_num) = r
        material_rgba(4,material_num) = 1.0E+00
      end if
!
!  V X Y Z W
!  Geometric vertex.
!
!  (X,Y,Z) is the coordinate of the vertex.
!  W is optional, a weight used for rational curves and surfaces.
!  The default for W is 1.
!
    else if ( s_eqi ( word1, 'V' ) ) then

      cor3_num = cor3_num + 1

      do i = 1, 3
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, temp, ierror, lchar )
        if ( cor3_num <= cor3_max ) then
          cor3(i,cor3_num) = temp
        end if
      end do
!
!  VN
!  Vertex normals.
!
    else if ( s_eqi ( word1, 'VN' ) ) then

      vertex_normal_num = vertex_normal_num + 1

      if ( vertex_normal_num <= order_max*face_max ) then

        do i = 1, 3
          call word_next_read ( line, word, done )
          call s_to_r4 ( word, temp, ierror, lchar )
          normal_temp(i,vertex_normal_num) = temp
        end do

      end if
!
!  VT
!  Vertex texture.
!
    else if ( s_eqi ( word1, 'VT' ) ) then
!
!  VP
!  Parameter space vertices.
!
    else if ( s_eqi ( word1, 'VP' ) ) then
!
!  Unrecognized keyword.
!
    else

      bad_num = bad_num + 1

      if ( bad_num <= 10 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) 'OBJ_READ: Bad data on line ', text_num
        write ( *, '(a)' ) '  Bad word: ' // trim ( word )
      end if

    end if

  end do
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'OBJ_READ - Read ', text_num, ' text lines from ' &
    // trim ( filein_name )

  return
end
subroutine obj_write ( cor3, cor3_max, cor3_num, face, face_max, face_num, &
  face_order, filein_name, fileout_name, iunit, line_dex, line_max, &
  line_num, order_max, vertex_normal )

!*****************************************************************************80
!
!! OBJ_WRITE writes graphics information to a WaveFront OBJ file.
!
!  Example:
!
!    #  magnolia.obj
!
!    mtllib ./vp.mtl
!
!    g Group001
!    v -3.269770 -39.572201 0.876128
!    v -3.263720 -39.507999 2.160890
!    ...
!    v 0.000000 -9.988540 0.000000
!    vn 0.0 1.0 0.0E+00
!    vn 1.0 0.0 0.0E+00
!    ...
!    vn 0.0 0.0 1.0E+00
!    s 1
!    usemtl brownskn
!    f 8//1 9//2 11//3 10//4
!    f 12//5 13//6 15//7 14//8
!    ...
!    f 788//800 806//803 774//807
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  character ( len = * ) fileout_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) indexvn
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) nl
  integer ( kind = 4 ), parameter :: OFFSET = 1
  character ( len = 256 ) text
  integer ( kind = 4 ) text_num
  character ( len = 256 ) text2
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  real ( kind = 4 ) w

  text_num = 0
  write ( iunit, '(a)' ) '# ' // trim ( fileout_name ) // ' created by IVREAD.'
  write ( iunit, '(a)' ) '# Original data in ' // trim ( filein_name ) // '.'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) 'g Group001'
  write ( iunit, '(a)' ) ' '

  text_num = text_num + 5
!
!  V: vertex coordinates.
!
  w = 1.0E+00
  do j = 1, cor3_num
    write ( text, '(a1,2x,4g14.6)' ) 'v', cor3(1:3,j), w
    call s_blanks_delete ( text )
    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1
  end do
!
!  VN: vertex face normal vectors.
!
  if ( 0 < face_num ) then
    write ( iunit, '(a)' ) ' '
    text_num = text_num + 1
  end if

  do iface = 1, face_num

    do ivert = 1, face_order(iface)

      write ( text, '(a2,2x,3f7.3)' ) 'vn', vertex_normal(1:3,ivert,iface)
      call s_blanks_delete ( text )
      write ( iunit, '(a)' ) trim ( text )
      text_num = text_num + 1

    end do

  end do
!
!  F: Faces, specified as 
!    vertex index/texture vertex index/normal index
!
  if ( 0 < face_num ) then
    write ( iunit, '(a)' ) ' '
    text_num = text_num + 1
  end if

  indexvn = 0

  do iface = 1, face_num

    text = 'f'

    do ivert = 1, face_order(iface)
      indexvn = indexvn + 1
      text2 = ' '     
      write ( text2(2:), '(i8, ''//'', i8 )' ) face(ivert,iface), indexvn
      call s_blank_delete ( text2(2:) )
      call s_cat ( text, text2, text )
    end do

    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1

  end do
!
!  L: lines, specified as a sequence of vertex indices.
!
  nl = 0

  if ( 0 < line_num ) then

    write ( iunit, '(a)' ) ' '
    text_num = text_num + 1
  
    ihi = 0

    ilo = ihi

    do

      do

        ihi = ihi + 1

        if ( line_num <= ihi ) then
          exit
        end if

        if ( line_dex(ihi) == -1 + OFFSET ) then
          exit
        end if

      end do

      write ( iunit, '(a,20i8)' ) 'l', line_dex(ilo+1:ihi-1)

      text_num = text_num + 1
      nl = nl + 1

      if ( line_num <= ihi ) then
        exit
      end if

      ilo = ihi

    end do

  end if
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'OBJ_WRITE - Wrote ', text_num, &
    ' text lines to ' // trim ( fileout_name )

  return
end
subroutine object_build ( face, face_num, face_object, face_order, &
  face_rank, face_tier, object_num, order_max )

!*****************************************************************************80
!
!! OBJECT_BUILD builds edge-connected "objects" out of polygonal faces.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_NUM), the nodes making faces.
!
!    Output, integer ( kind = 4 ) FACE_OBJECT(FACE_NUM), describes the objects.
!    FACE_OBJECT(I) is the index of the edge-connected "object" to 
!    which face I belongs.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_NUM), the number of nodes per face.
!
!    Output, integer ( kind = 4 ) FACE_RANK(FACE_NUM), is an ordered list of faces.
!    FACE_RANK(1) is the index of the face in the first tier of the 
!    first object, followed by second tier faces, and so on until
!    object one is complete.  Object two follows, and so on.
!
!    Output, integer ( kind = 4 ) FACE_TIER(FACE_NUM).  FACE_TIER(I) is the "tier"
!    of face I in its object.  The seed of the object has tier 1,
!    the neighbors of the seed have tier 2, and so on.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, is the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Output, integer ( kind = 4 ) OBJECT_NUM, the number of objects.
!
  implicit none

  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) face_num

  logical, parameter :: DEBUG = .false.
  integer ( kind = 4 ) face(order_max,face_num)
  integer ( kind = 4 ) face_object(face_num)
  integer ( kind = 4 ) face_order(face_num)
  integer ( kind = 4 ) face_rank(face_num)
  integer ( kind = 4 ) face_tier(face_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ihi_next
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ilo_next
  integer ( kind = 4 ) iprint
  integer ( kind = 4 ) irank
  integer ( kind = 4 ) jface
  integer ( kind = 4 ) jrank
  integer ( kind = 4 ) object_num
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) tier
  integer ( kind = 4 ) touch
!
!  Initialization.
!
  iprint = 1
  object_num = 0

  if ( face_num <= 0 ) then
    return
  end if

  face_object(1:face_num) = 0
  face_rank(1:face_num) = 0
  face_tier(1:face_num) = 0

  irank = 0

  seed = 1
!
!  Begin the next object, seeded with face SEED.
!
10    continue

  tier = 1

  object_num = object_num + 1
  irank = irank + 1
  jrank = irank
  if ( irank == iprint .or. irank == face_num ) then
    write ( *, '(2i6)' ) irank, seed
    iprint = 2 * iprint
  end if

  face_rank(irank) = seed
  face_tier(seed) = tier
  face_object(seed) = object_num

  ilo = irank
  ihi = irank
!
!  Begin the next tier of faces, which are neighbors of faces we
!  found in the previous tier.
!
20    continue

  tier = tier + 1

  ilo_next = ihi + 1
  ihi_next = ihi

  do jface = 1, face_num

    if ( face_tier(jface) == 0 ) then

      do i = ilo, ihi

        iface = face_rank(i)

        call face_touch ( face, face_order, order_max, face_num, iface, &
          jface, touch )

        if ( touch /= 0 ) then
          if ( DEBUG ) then
            write ( *, '(a,2i6)' ) 'Touching faces: ', iface, jface
          end if
          ihi_next = ihi_next + 1
          irank = irank + 1
          if ( irank == iprint .or. irank == face_num ) then
            write ( *, '(2i6)' ) irank, jface
            iprint = 2 * iprint
          end if
          face_rank(irank) = jface
          face_tier(jface) = tier
          face_object(jface) = object_num
          go to 30
        end if

      end do

    end if

30      continue

  end do

  if ( ilo_next <= ihi_next ) then
    ilo = ilo_next
    ihi = ihi_next
    go to 20
  end if

  write ( *, '(a,i6,a,i6,a)' ) 'Object ', object_num, ' uses ', &
    irank + 1 - jrank, ' faces.'
  jrank = irank
!
!  No neighbors were found, so this object is complete.  
!  Search for an unused face, which will be the seed of the next object.
!
  do iface = 1, face_num

    if ( face_tier(iface) == 0 ) then
      seed = iface
      go to 10
    end if

  end do

  return
end
subroutine object_invert ( cor3, cor3_material, cor3_max, cor3_normal, &
  cor3_num, face, face_material, face_max, face_normal, face_num, &
  face_order, material_max, material_name, material_num, material_rgba, &
  order_max, vertex_material, vertex_normal )

!*****************************************************************************80
!
!! OBJECT_INVERT makes an inverted duplicate of the object.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices 
!    per face.
!
!    Input/output, character ( len = * ) MATERIAL_NAME(MATERIAL_MAX), 
!    material names.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Output, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Input/output, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_material(cor3_max)
  real ( kind = 4 ) cor3_normal(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  real ( kind = 4 ), parameter :: EPS = 0.01E+00
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) icor32
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) iface2
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) ivert2
  integer ( kind = 4 ) j
  character ( len = * ) material_name(material_max)
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
!
!  Check.
!
  if ( face_max < 2 * face_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'OBJECT_INVERT - Fatal error!'
    write ( *, '(a)' ) '  2 * FACE_NUM exceeds FACE_MAX.'
    return
  end if

  if ( cor3_max < 2 * cor3_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'OBJECT_INVERT - Fatal error!'
    write ( *, '(a)' ) '  2 * COR3_NUM exceeds COR3_MAX.'
    return
  end if
!
!  If there aren't 5 materials, add some.
!
  if ( material_num < 5 ) then

    material_num = material_num + 1
    if ( material_num <= material_max ) then
      material_name(material_num) = 'Red_Material'
      material_rgba(1,material_num) = 1.0E+00
      material_rgba(2,material_num) = 0.0E+00
      material_rgba(3,material_num) = 0.0E+00
      material_rgba(4,material_num) = 1.0E+00
      write ( *, '(a)' ) 'OBJECT_INVERT - Adding dummy red material'
    end if
  end if

  if ( material_num < 5 ) then

    material_num = material_num + 1
    if ( material_num <= material_max ) then
      material_name(material_num) = 'Green_Material'
      material_rgba(1,material_num) = 0.0E+00
      material_rgba(2,material_num) = 1.0E+00
      material_rgba(3,material_num) = 0.0E+00
      material_rgba(4,material_num) = 1.0E+00
      write ( *, '(a)' ) 'OBJECT_INVERT - Adding dummy green material'
    end if
  end if
!
!  Generate new points, displaced by EPS in the negative direction.
!
  do icor3 = 1, cor3_num

    icor32 = icor3 + cor3_num

    if ( cor3_material(icor3) == 1 ) then
      cor3_material(icor32) = 2
    else if ( cor3_material(icor3) == 3 ) then
      cor3_material(icor32) = 4
    else if ( cor3_material(icor3) == 4 ) then
      cor3_material(icor3) = 3
      cor3_material(icor32) = 4
    else if ( cor3_material(icor3) == 5 ) then
      cor3_material(icor3) = 3
      cor3_material(icor32) = 5
    end if

    cor3(1:3,icor32) = cor3(1:3,icor3) - EPS * cor3_normal(1:3,icor3)
    cor3_normal(1:3,icor32) = - cor3_normal(1:3,icor3)

  end do
!
!  Generate new faces.
!
  do iface = 1, face_num

    iface2 = face_num + iface
    face_order(iface2) = face_order(iface)

    if ( face_material(iface) == 1 ) then
      face_material(iface2) = 2
    else if ( face_material(iface) == 3 ) then
      face_material(iface2) = 4
    else if ( face_material(iface) == 4 ) then
      face_material(iface) = 3
      face_material(iface2) = 4
    else if ( face_material(iface) == 5 ) then
      face_material(iface) = 3
      face_material(iface2) = 5
    end if

    do ivert = 1, face_order(iface)

      ivert2 = face_order(iface) + 1 - ivert
      face(ivert2,iface2) = face(ivert,iface) + cor3_num

      if ( vertex_material(ivert,iface) == 1 ) then
        vertex_material(ivert2,iface2) = 2
      else if ( vertex_material(ivert,iface) == 3 ) then
        vertex_material(ivert2,iface2) = 4
      else if ( vertex_material(ivert,iface) == 4 ) then
        vertex_material(ivert,iface) = 3
        vertex_material(ivert2,iface2) = 4
      else if ( vertex_material(ivert,iface) == 5 ) then
        vertex_material(ivert,iface) = 3
        vertex_material(ivert2,iface2) = 5
      end if

      vertex_normal(1:3,ivert2,iface2) = - vertex_normal(1:3,ivert,iface)

    end do

    face_normal(1:3,iface2) = - face_normal(1:3,iface)

  end do

  cor3_num = 2 * cor3_num
  face_num = 2 * face_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'OBJECT_INVERT: Information:'
  write ( *, '(a,i6)' ) '  Number of points = ', cor3_num
  write ( *, '(a,i6)' ) '  Number of faces =  ', face_num

  return
end
subroutine off_read ( cor3, cor3_max, cor3_num, face, face_max, &
  face_num, face_order, filein_name, ierror, iunit, order_max, text_num )

!*****************************************************************************80
!
!! OFF_READ reads graphics information from a GEOMVIEW OFF file.
!
!  Example:
!
!    OFF
!    8  6  12
!    0.0 0.0 0.0
!    0.0 0.0 1.0
!    0.0 1.0 0.0
!    0.0 1.0 1.0
!    1.0 0.0 0.0
!    1.0 0.0 1.0
!    1.0 1.0 0.0
!    1.0 1.0 1.0
!    4  0 2 3 1
!    4  4 5 7 6
!    4  0 1 5 4
!    4  2 6 7 3
!    4  0 4 6 2
!    4  1 3 7 5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input/output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input/output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices
!    per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which data is read.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Output, integer ( kind = 4 ) TEXT_NUM, the number of lines of text read from
!    the file.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) cor3_num_base
  logical done
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) lchar
  character ( len = 256 ) line
  integer ( kind = 4 ), parameter :: OFFSET = 1
  logical s_eqi
  integer ( kind = 4 ) text_num
  real ( kind = 4 ) temp
  character ( len = 256 ) word

  ierror = 0
!
!  Save a copy of the input value of COR3_NUM to use as a base.
!
  cor3_num_base = cor3_num

  text_num = 0

  word = ' '
!
!  Read a line of text from the file that is not "OFF" or blank or a comment.
!
  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = 1
      return
    end if

    text_num = text_num + 1

    if ( line(1:1) == '#' ) then
      cycle
    end if

    if ( s_eqi ( line(1:3), 'OFF' ) ) then
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    exit

  end do
!
!  Presumably, LINE now contains the text defining the problem size.
!
  done = .true.

  call word_next_read ( line, word, done )
  call s_to_i4 ( word, cor3_num, ierror, lchar )

  call word_next_read ( line, word, done )
  call s_to_i4 ( word, face_num, ierror, lchar )

  call word_next_read ( line, word, done )
  call s_to_i4 ( word, edge_num, ierror, lchar )
!
!  Now we expect the coordinates of each vertex.
!
  do j = 1, cor3_num

    done = .true.
    read ( iunit, '(a)', iostat = ios ) line
    if ( ios /= 0 ) then
      ierror = 1
      return
    end if
    text_num = text_num + 1

    do i = 1, 3
      call word_next_read ( line, word, done )
      call s_to_r4 ( word, cor3(i,j), ierror, lchar )
    end do

  end do
!
!  Now we expect the order of each face, and the vertex list.
!
  do j = 1, face_num

    done = .true.
    read ( iunit, '(a)', iostat = ios ) line
    if ( ios /= 0 ) then
      ierror = 1
      return
    end if
    text_num = text_num + 1

    call word_next_read ( line, word, done )
    call s_to_i4 ( word, face_order(j), ierror, lchar )

    do i = 1, face_order(j)
      call word_next_read ( line, word, done )
      call s_to_i4 ( word, face(i,j), ierror, lchar )
      face(i,j) = face(i,j) + OFFSET
    end do

  end do
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'OFF_READ - Read ', text_num, ' text lines from ' &
    // trim ( filein_name )

  return
end
subroutine off_write ( cor3, cor3_max, cor3_num, face, face_max, face_num, &
  face_order, fileout_name, iunit, order_max )

!*****************************************************************************80
!
!! OFF_WRITE writes graphics information to a GEOMVIEW OFF file.
!
!  Example:
!
!    OFF
!    8  6  12
!    0.0 0.0 0.0
!    0.0 0.0 1.0
!    0.0 1.0 0.0
!    0.0 1.0 1.0
!    1.0 0.0 0.0
!    1.0 0.0 1.0
!    1.0 1.0 0.0
!    1.0 1.0 1.0
!    4  0 2 3 1
!    4  4 5 7 6
!    4  0 1 5 4
!    4  2 6 7 3
!    4  0 4 6 2
!    4  1 3 7 5
!    
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) fileout_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: OFFSET = 1
  character ( len = 256 ) text
  integer ( kind = 4 ) text_num
  character ( len = 256 ) text2
!
!  "Magic Number"
!
  text_num = 0
  write ( iunit, '(a)' ) 'OFF'
  text_num = text_num + 1
!
!  Compute EDGE_NUM.
!
  call edge_count ( face_max, order_max, face_num, face, face_order, edge_num )
!
!  Counts.
!
  write ( text, '(3i6)' ) cor3_num, face_num, edge_num
  call s_blanks_delete ( text )
  write ( iunit, '(a)' ) trim ( text )
  text_num = text_num + 1
!
!  Vertex coordinates.
!
  do j = 1, cor3_num
    write ( text, '(2x,3g14.6)' ) cor3(1:3,j)
    call s_blanks_delete ( text )
    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1
  end do
!
!  Face count, and vertex indices.
!
  do iface = 1, face_num

    write ( text, '(2x,i4)' ) face_order(iface)

    do ivert = 1, face_order(iface)
      write ( text2, '(2x,i8)' ) face(ivert,iface) - OFFSET
      call s_blank_delete ( text2(3:) )
      call s_cat ( text, text2, text )
    end do

    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1

  end do
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'OFF_WRITE - Wrote ', text_num, &
    ' text lines to ' // trim ( fileout_name )

  return
end
subroutine oogl_read ( cor3, cor3_material, cor3_max, cor3_normal, cor3_num, &
  face, face_area, face_material, face_max, face_normal, face_num, &
  face_order, filein_name, ierror, iunit, material_max, material_name, &
  material_num, material_rgba, ncol_oogl, nrow_oogl, order_max, text_num, &
  vertex_material, vertex_normal )

!*****************************************************************************80
!
!! OOGL_READ reads graphics information from a OOGL file.
!
!  Diagnostics:
!
!    Note that raw READ statements are used.  As written, the
!    code can't handle a blank line, or deal with a case where
!    information runs over to a new line.
!
!  Discussion:
!
!    It is intended that the information read from the file can
!    either start a whole new graphics object, or simply be added
!    to a current graphics object via the '<<' command.  
!
!    This is controlled by whether the input values have been zeroed 
!    out or not.  This routine simply tacks on the information it 
!    finds to the current graphics object.
!
!  Example:
!
!    {CMESH
!      5 3
!      0.0  0.0  0.0  0.94  0.70  0.15  1.00
!      1.0  0.0  0.0  0.94  0.70  0.15  1.00
!      2.0  0.0  0.0  0.94  0.70  0.15  1.00
!      3.0  0.0  0.0  0.94  0.70  0.15  1.00
!      4.0  0.0  0.0  0.94  0.70  0.15  1.00
!      0.0  1.0  0.0  0.94  0.70  0.15  1.00
!      1.0  1.0  0.0  0.94  0.70  0.15  1.00
!      2.0  1.0  0.0  0.94  0.70  0.15  1.00
!      3.0  1.0  0.0  0.94  0.70  0.15  1.00
!      4.0  1.0  0.0  0.94  0.70  0.15  1.00
!      0.0  2.0  0.0  0.94  0.70  0.15  1.00
!      1.0  2.0  0.0  0.94  0.70  0.15  1.00
!      2.0  2.0  0.0  0.94  0.70  0.15  1.00
!      3.0  2.0  0.0  0.94  0.70  0.15  1.00
!      4.0  2.0  0.0  0.94  0.70  0.15  1.00
!    }
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input/output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Workspace, real FACE_AREA(FACE_MAX), the area of each face.
! 
!    Input/output, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input/output, real ( kind = 4 ) FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices 
!    per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which data is read.
!
!    Input/output, character ( len = * ) MATERIAL_NAME(MATERIAL_MAX), 
!    material names.
!
!    Input/output, real ( kind = 4 ) MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B 
!    and A values.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input/output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input/output, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Output, integer ( kind = 4 ) TEXT_NUM, the number of lines of text read from
!    the file.
!
!    ?, integer NCOL_OOGL, ?
!
!    ?, integer NROW_OOGL, ?
!
!    Input/output, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    ?, real VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), ?
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  integer ( kind = 4 ) b2l2
  integer ( kind = 4 ) b2r2
  integer ( kind = 4 ) black
  character ( len = 4 ) char4
  logical clipping
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_material(cor3_max)
  real ( kind = 4 ) cor3_normal(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  real ( kind = 4 ) dist
  integer ( kind = 4 ) face(order_max,face_max)
  real ( kind = 4 ) face_area(face_max)
  integer ( kind = 4 ) face_material(face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  real ( kind = 4 ) g
  logical grid
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icor3
  logical identify
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) imat
  logical invert
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jvert
  integer ( kind = 4 ) k
  character ( len = * ) material_name(material_max)
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ) nclip
  integer ( kind = 4 ) ncol_oogl
  integer ( kind = 4 ) ngold
  integer ( kind = 4 ) nrow_oogl
  real ( kind = 4 ) r
  real ( kind = 4 ) rgba(4)
  integer ( kind = 4 ) t2l2
  integer ( kind = 4 ) t2r2
  integer ( kind = 4 ) text_num
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  integer ( kind = 4 ) white
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  real ( kind = 4 ) z

  ierror = 0
  text_num = 0
!
!  30 March 1999
!  KLUDGE WARNING:
!  For our current purposes, we need to have the following materials defined.
!
!  Set up material #1 = GOLD.
!
  material_num = material_num + 1
  if ( material_num <= material_max ) then
    material_name(material_num) = 'Gold_Material'
    material_rgba(1,material_num) = 0.94E+00
    material_rgba(2,material_num) = 0.70E+00
    material_rgba(3,material_num) = 0.15E+00
    material_rgba(4,material_num) = 1.0E+00
  end if
!
!  Set up material #2 = DARK BLUE.
!
  material_num = material_num + 1
  if ( material_num <= material_max ) then
    material_name(material_num) = 'Dark_Blue_Material'
    material_rgba(1,material_num) = 0.24E+00
    material_rgba(2,material_num) = 0.00E+00
    material_rgba(3,material_num) = 0.85E+00
    material_rgba(4,material_num) = 1.0E+00
  end if
!
!  Set up material #3 = LIGHT BLUE.
!
  material_num = material_num + 1
  if ( material_num <= material_max ) then
    material_name(material_num) = 'Light_Blue_Material'
    material_rgba(1,material_num) = 0.24E+00
    material_rgba(2,material_num) = 0.70E+00
    material_rgba(3,material_num) = 0.85E+00
    material_rgba(4,material_num) = 1.0E+00
  end if
!
!  A node with grid coordinates (I,J) will be mapped to a
!  node number ( I - 1 ) * NCOL + J
!
  read ( iunit, *, iostat = ios )

  if ( ios /= 0 ) then
    go to 50
  end if

  text_num = text_num + 1

  read ( iunit, *, iostat = ios ) ncol_oogl, nrow_oogl

  if ( ios /= 0 ) then
    go to 50
  end if

  text_num = text_num + 1

  do i = 1, nrow_oogl
    do j = 1, ncol_oogl

      read ( iunit, *, iostat = ios ) x, y, z, r, g, b, a

      if ( ios /= 0 ) then
        go to 50
      end if

      text_num = text_num + 1

      cor3_num = cor3_num + 1
!
!  25 February 1999:
!  In the data we've seen, the first and last columns have almost
!  identical X,Y,Z coordinates and it might help if they were identical.
!
      if ( cor3_num <= cor3_max ) then

        cor3(1,cor3_num) = x
        cor3(2,cor3_num) = y
        cor3(3,cor3_num) = z

        if ( j == ncol_oogl ) then

          itemp = cor3_num + 1 - ncol_oogl

          dist = sqrt ( ( cor3(1,cor3_num) - cor3(1,itemp) )**2 &
                      + ( cor3(2,cor3_num) - cor3(2,itemp) )**2 &
                      + ( cor3(3,cor3_num) - cor3(3,itemp) )**2 )

          if ( 0.0E+00 < dist .and. dist < 0.000001E+00 ) then
            cor3(1,cor3_num) = cor3(1,itemp)
            cor3(2,cor3_num) = cor3(2,itemp)
            cor3(3,cor3_num) = cor3(3,itemp)
          end if

        end if

      end if
!
!  Some of the input data has had RGBA values that are negative.
!  Do not allow this.
!
      rgba(1) = min ( max ( 0.0E+00, r ), 1.0E+00 )
      rgba(2) = min ( max ( 0.0E+00, g ), 1.0E+00 )
      rgba(3) = min ( max ( 0.0E+00, b ), 1.0E+00 )
      rgba(4) = min ( max ( 0.0E+00, a ), 1.0E+00 )
!
!  See if the RGBA values of this material match those of a material
!  that has already been defined.
!
      if ( material_num <= 1000 ) then
        call r4col_find ( 4, 4, material_num, material_rgba, rgba, imat )
      else
        imat = 0
      end if

      if ( imat == 0 ) then

        material_num = material_num + 1

        if ( material_num <= material_max ) then

          call i4_to_s_zero ( material_num, char4 )

          material_name(material_num) = 'Material_' // char4              
          material_rgba(1:4,material_num) = rgba(1:4)
          imat = material_num

        else

          imat = 0

        end if

      end if

      if ( cor3_num <= cor3_max ) then
        cor3_material(cor3_num) = imat
      end if
  
    end do
  end do

  read ( iunit, *, iostat = ios )

  if ( ios /= 0 ) then
    go to 50
  end if

  text_num = text_num + 1
!
!  Now set up the faces from the grid of points that were defined.
!
  clipping = .false.
  nclip = 0

  do i = 2, nrow_oogl
    do j = 2, ncol_oogl

      b2l2 = ( i - 2 ) * ncol_oogl + j - 1
      b2r2 = ( i - 2 ) * ncol_oogl + j
      t2l2 = ( i - 1 ) * ncol_oogl + j - 1
      t2r2 = ( i - 1 ) * ncol_oogl + j

      if ( .not. clipping  ) then

        face_num = face_num + 1
        face_material(face_num) = cor3_material(t2r2)
        face_order(face_num) = 3

        k = t2r2
        face(1,face_num) = k
        vertex_material(1,face_num) = cor3_material(t2r2)

        k = b2r2
        face(2,face_num) = k
        vertex_material(2,face_num) = cor3_material(b2r2)

        k = b2l2
        face(3,face_num) = k
        vertex_material(3,face_num) = cor3_material(b2l2)

        face_num = face_num + 1
        face_material(face_num) = cor3_material(t2l2)
        face_order(face_num) = 3

        k = t2l2
        face(1,face_num) = k
        vertex_material(1,face_num) = cor3_material(t2l2)

        k = t2r2
        face(2,face_num) = k
        vertex_material(2,face_num) = cor3_material(t2r2)

        k = b2l2
        face(3,face_num) = k
        vertex_material(3,face_num) = cor3_material(b2l2)

      else 

        ngold = 0
        if ( cor3_material(t2r2) == 3 ) then
          ngold = ngold + 1
        end if
        if ( cor3_material(b2r2) == 3 ) then
          ngold = ngold + 1
        end if
        if ( cor3_material(b2l2) == 3 ) then
          ngold = ngold + 1
        end if
        if ( cor3_material(t2l2) == 3 ) then
          ngold = ngold + 1
        end if

        if ( 3 <= ngold ) then

          face_num = face_num + 1
          face_material(face_num) = 3
          face_order(face_num) = 3

          k = t2r2
          face(1,face_num) = k
          vertex_material(1,face_num) = 3

          k = b2r2
          face(2,face_num) = k
          vertex_material(2,face_num) = 3

          k = b2l2
          face(3,face_num) = k
          vertex_material(3,face_num) = 3

          face_num = face_num + 1
          face_material(face_num) = 3
          face_order(face_num) = 3

          k = t2l2
          face(1,face_num) = k
          vertex_material(1,face_num) = 3

          k = t2r2
          face(2,face_num) = k
          vertex_material(2,face_num) = 3

          k = b2l2
          face(3,face_num) = k
          vertex_material(3,face_num) = 3

        else

          nclip = nclip + 2

        end if

      end if

    end do
  end do

  if ( 0 < nclip ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'OOGL_READ:'
    write ( *, '(a,i6,a)' ) '  "Clipped" ', nclip, ' faces.'
  end if
!
!  Set "identify" TRUE to identify the repeated points at the first
!  and last, and along the sides.
!
  identify = .true.

  if ( identify ) then

    do iface = 1, face_num
      do jvert = 1, face_order(iface)
 
        i = face(jvert,iface)
        if ( i <= ncol_oogl ) then
          face(jvert,iface) = 1
        else if ( i <= nrow_oogl * ncol_oogl .and. &
          ( nrow_oogl - 1 ) * ncol_oogl < i ) then
          face(jvert,iface) = ( nrow_oogl - 1 ) * ncol_oogl + 1
        end if

      end do
    end do

    do iface = 1, face_num
      do jvert = 1, face_order(iface)
        i = face(jvert,iface)
        if ( i <= nrow_oogl * ncol_oogl .and. mod ( i, ncol_oogl ) == 0 ) then
          face(jvert,iface) = face(jvert,iface) - ncol_oogl + 1
        end if

      end do
    end do

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'OOGL_READ - Information:'
  write ( *, '(a,i6)' ) '  Initial number of points = ', cor3_num
  write ( *, '(a,i6)' ) '  Initial number of faces =  ', face_num

  invert = .false.
  grid = .false.

  if ( invert .or. grid) then
!
!  Set up the normal vector information.
!
    call vertex_normal_set ( cor3, cor3_max, face, face_max, &
      face_num, face_order, order_max, vertex_normal )

    call face_area_set ( cor3, cor3_max, face, face_area, face_max, &
      face_num, face_order, order_max ) 

    cor3_normal(1:3,1:cor3_num) = 0.0E+00

    call cor3_normal_set ( cor3_max, cor3_normal, face, face_area, &
      face_max, face_num, face_order, order_max, vertex_normal )

  end if
!
!  Make the other side of the surface.
!
  if ( invert ) then
    call object_invert ( cor3, cor3_material, cor3_max, cor3_normal, &
      cor3_num, face, face_material, face_max, face_normal, face_num, &
      face_order, material_max, material_name, material_num, material_rgba, &
      order_max, vertex_material, vertex_normal )
  end if
!
!  Make the grid.
!
  if ( grid ) then

    black = 1
    white = 2

    call oogl_grid ( cor3, cor3_material, cor3_normal, face, face_material, &
      face_order, cor3_max, face_max, order_max, cor3_num, face_num, &
      ncol_oogl, nrow_oogl, black, white, invert, vertex_material )

  end if
!
!  06 April 1999
!  KLUDGE #2 WARNING:
!  For our current purposes, we need to guarantee that there 
!  is at least one node of each material.
!
  cor3_num = cor3_num + 1
  cor3(1,cor3_num) = cor3(1,1)
  cor3(2,cor3_num) = cor3(1,1)
  cor3(3,cor3_num) = cor3(1,1)
  cor3_material(cor3_num) = 1
  cor3_normal(1,cor3_num) = 1.0E+00
  cor3_normal(2,cor3_num) = 0.0E+00
  cor3_normal(3,cor3_num) = 0.0E+00

  cor3_num = cor3_num + 1
  cor3(1,cor3_num) = cor3(1,1)
  cor3(2,cor3_num) = cor3(1,1)
  cor3(3,cor3_num) = cor3(1,1)
  cor3_material(cor3_num) = 2
  cor3_normal(1,cor3_num) = 1.0E+00
  cor3_normal(2,cor3_num) = 0.0E+00
  cor3_normal(3,cor3_num) = 0.0E+00

  cor3_num = cor3_num + 1
  cor3(1,cor3_num) = cor3(1,1)
  cor3(2,cor3_num) = cor3(1,1)
  cor3(3,cor3_num) = cor3(1,1)
  cor3_material(cor3_num) = 3
  cor3_normal(1,cor3_num) = 1.0E+00
  cor3_normal(2,cor3_num) = 0.0E+00
  cor3_normal(3,cor3_num) = 0.0E+00
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'OOGL_READ - Read ', text_num, ' text lines from ' &
    // trim ( filein_name )

  return
!
!  Unexpected end of information.
!
50    continue

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'OOGL_READ - Fatal error!'
  write ( *, '(a)' ) '  Unexpected end of information!'

  return
end
subroutine oogl_grid ( cor3, cor3_material, cor3_normal, face, face_material, &
  face_order, cor3_max, face_max, order_max, cor3_num, face_num, ncol_oogl, &
  nrow_oogl, black, white, invert, vertex_material )

!*****************************************************************************80
!
!! OOGL_GRID adds a grid to an OOGL data file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Output, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which data is read.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Output, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Output, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  integer ( kind = 4 ) base
  integer ( kind = 4 ) b2l2
  integer ( kind = 4 ) b2r1u
  integer ( kind = 4 ) b2r2
  integer ( kind = 4 ) b2r2u
  integer ( kind = 4 ) black
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_material(cor3_max)
  real ( kind = 4 ) cor3_normal(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  real ( kind = 4 ), parameter :: EPS = 0.0025E+00
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) grid_nx
  integer ( kind = 4 ), parameter :: GRID_NX_NUM = 20
  integer ( kind = 4 ) grid_ny
  integer ( kind = 4 ), parameter :: GRID_NY_NUM = 20
  real ( kind = 4 ), parameter :: GRID_WIDTH = 0.008E+00
  integer ( kind = 4 ) i
  logical invert
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ncol_oogl
  real ( kind = 4 ) norm
  integer ( kind = 4 ) nrow_oogl
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) t1l2u
  integer ( kind = 4 ) t1r2u
  integer ( kind = 4 ) t2l2
  integer ( kind = 4 ) t2l2u
  integer ( kind = 4 ) t2r1u
  integer ( kind = 4 ) t2r2
  integer ( kind = 4 ) t2r2u
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  integer ( kind = 4 ) white
!
!  Determine the number of grid lines.
!
  grid_nx = ncol_oogl / GRID_NX_NUM
  grid_nx = min ( grid_nx, ncol_oogl - 1 )
  grid_nx = max ( grid_nx, 1 )

  grid_ny = nrow_oogl / GRID_NY_NUM
  grid_ny = min ( grid_ny, nrow_oogl - 1 )
  grid_ny = max ( grid_ny, 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'OOGL_GRID:'
  write ( *, '(a,i6)' ) '  Grid spacing is ', grid_nx
  write ( *, '(a,i6)' ) '  by ', grid_ny
!
!  Do the grid lines that I think of as running along the "top" of
!  the affected faces.
!
!    T2L2U---- ---- T2R2U
!    T1L2U---- ---- T1R2U
!    |    .... ....    |
!    +--- ---- ---- ---+
!
!    T2L2 ---- ---- T2R2
!    |    .... ....    |
!    |    .... ....    |
!    B2L2 ---- ---- B2R2
!
  do i = grid_ny+1, nrow_oogl, grid_ny
    do j = 2, ncol_oogl

      b2l2 = ( i - 2 ) * ncol_oogl + j - 1
      b2r2 = ( i - 2 ) * ncol_oogl + j
      t2l2 = ( i - 1 ) * ncol_oogl + j - 1
      t2r2 = ( i - 1 ) * ncol_oogl + j

      cor3_num = cor3_num + 1
      t2l2u = cor3_num
      cor3(1:3,t2l2u) = cor3(1:3,t2l2) + EPS * cor3_normal(1:3,t2l2)
      cor3_material(cor3_num) = black

      cor3_num = cor3_num + 1
      t2r2u = cor3_num
      do k = 1, 3
        cor3(k,t2r2u) = cor3(k,t2r2) + EPS * cor3_normal(k,t2r2)
      end do
      cor3_material(cor3_num) = black

      cor3_num = cor3_num + 1
      t1l2u = cor3_num

      norm = 0.0E+00
      do k = 1, 3
        norm = norm + ( cor3(k,b2l2) + EPS * cor3_normal(k,b2l2) &
          - cor3(k,t2l2) )**2
      end do
      norm = sqrt ( norm )
      if ( norm /= 0.0E+00 ) then
        norm = 1.0E+00 / norm
      end if

      do k = 1, 3
        cor3(k,t1l2u) = cor3(k,t2l2u) + GRID_WIDTH * ( cor3(k,b2l2) &
          + EPS * cor3_normal(k,b2l2) - cor3(k,t2l2) ) * norm
      end do

      cor3_material(cor3_num) = black


      cor3_num = cor3_num + 1
      t1r2u = cor3_num

      norm = 0.0E+00
      do k = 1, 3
        norm = norm + ( cor3(k,b2r2) + EPS * cor3_normal(k,b2r2) &
          - cor3(k,t2r2) )**2
      end do
      norm = sqrt ( norm )
      if ( norm /= 0.0E+00 ) then
        norm = 1.0E+00 / norm
      end if

      do k = 1, 3
        cor3(k,t1r2u) = cor3(k,t2r2u) &
          + GRID_WIDTH * ( cor3(k,b2r2) + EPS * cor3_normal(k,b2r2) &
          - cor3(k,t2r2) ) * norm
      end do

      cor3_material(cor3_num) = black

      face_num = face_num + 1
      face_material(face_num) = black
      face_order(face_num) = 3
      nvert = 0

      nvert = nvert + 1
      k = t2l2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      nvert = nvert + 1
      k = t2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      nvert = nvert + 1
      k = t1l2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      face_num = face_num + 1
      face_material(face_num) = black
      face_order(face_num) = 3
      nvert = 0

      nvert = nvert + 1
      k = t1l2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      nvert = nvert + 1
      k = t2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      nvert = nvert + 1
      k = t1r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

    end do
  end do

  write ( *, '(a)' ) 'OOGL_GRID2: Done "top" grid.'
!
!  Do the grid lines that I think of as running along the "right" of
!  the affected faces.
!
!    +--- ---- T2R1U-T2R2U
!    |    ....    |    |
!    |    ....    |    |
!    +--- ---- B2R1U-B2R2U
!
!    T2L2 ---- ---- T2R2
!    |    .... ....    |
!    |    .... ....    |
!    B2L2 ---- ---- B2R2
!
  do j = grid_nx+1, ncol_oogl, grid_nx
    do i = 2, nrow_oogl

      b2l2 = ( i - 2 ) * ncol_oogl + j - 1
      b2r2 = ( i - 2 ) * ncol_oogl + j
      t2l2 = ( i - 1 ) * ncol_oogl + j - 1
      t2r2 = ( i - 1 ) * ncol_oogl + j

      cor3_num = cor3_num + 1
      b2r2u = cor3_num
      do k = 1, 3
        cor3(k,b2r2u) = cor3(k,b2r2) + EPS * cor3_normal(k,b2r2)
      end do
      cor3_material(cor3_num) = black

      cor3_num = cor3_num + 1
      t2r2u = cor3_num
      do k = 1, 3
        cor3(k,t2r2u) = cor3(k,t2r2) + EPS * cor3_normal(k,t2r2)
      end do
      cor3_material(cor3_num) = black


      cor3_num = cor3_num + 1
      t2r1u = cor3_num

      norm = 0.0E+00
      do k = 1, 3
        norm = norm + ( cor3(k,t2l2) + EPS * cor3_normal(k,t2l2) &
          - cor3(k,t2r2) )**2
      end do
      norm = sqrt ( norm )
      if ( norm /= 0.0E+00 ) then
        norm = 1.0E+00 / norm
      end if

      do k = 1, 3
        cor3(k,t2r1u) = cor3(k,t2r2u) + GRID_WIDTH * ( cor3(k,t2l2) + &
          EPS * cor3_normal(k,t2l2) - cor3(k,t2r2) ) * norm
      end do

      cor3_material(cor3_num) = black

      cor3_num = cor3_num + 1
      b2r1u = cor3_num

      norm = 0.0E+00
      do k = 1, 3
        norm = norm + ( cor3(k,b2l2) + EPS * cor3_normal(k,b2l2) &
          - cor3(k,b2r2) )**2
      end do
      norm = sqrt ( norm )
      if ( norm /= 0.0E+00 ) then
        norm = 1.0E+00 / norm
      end if

      do k = 1, 3
        cor3(k,b2r1u) = cor3(k,b2r2u) + GRID_WIDTH * ( cor3(k,b2l2) &
          + EPS * cor3_normal(k,b2l2) - cor3(k,b2r2) ) * norm
      end do

      cor3_material(cor3_num) = black

      face_num = face_num + 1
      face_material(face_num) = black
      face_order(face_num) = 3
      nvert = 0

      nvert = nvert + 1
      k = b2r1u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      nvert = nvert + 1
      k = t2r1u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      nvert = nvert + 1
      k = t2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      face_num = face_num + 1
      face_material(face_num) = black
      face_order(face_num) = 3
      nvert = 0

      nvert = nvert + 1
      k = b2r1u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      nvert = nvert + 1
      k = t2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

      nvert = nvert + 1
      k = b2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = black

    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'OOGL_GRID - Added black "grid lines".'
  write ( *, '(a,i6)' ) '  Number of points = ', cor3_num
  write ( *, '(a,i6)' ) '  Number of faces =  ', face_num
!
!  Do the grid lines that I think of as running along the "top" of
!  the affected faces.
!
!    T2L2U---- ---- T2R2U
!    T1L2U---- ---- T1R2U
!    |    .... ....    |
!    +--- ---- ---- ---+
!
!    T2L2 ---- ---- T2R2
!    |    .... ....    |
!    |    .... ....    |
!    B2L2 ---- ---- B2R2
!
  if ( .not. invert ) then
    return
  end if

  base = ncol_oogl * nrow_oogl
  write ( *, '(a,i6)' ) 'BASE = ', base

  do i = grid_ny+1, nrow_oogl, grid_ny
    do j = 2, ncol_oogl

      b2l2 = base + ( i - 2 ) * ncol_oogl + j - 1
      b2r2 = base + ( i - 2 ) * ncol_oogl + j
      t2l2 = base + ( i - 1 ) * ncol_oogl + j - 1
      t2r2 = base + ( i - 1 ) * ncol_oogl + j

      cor3_num = cor3_num + 1
      t2l2u = cor3_num
      do k = 1, 3
        cor3(k,t2l2u) = cor3(k,t2l2) + EPS * cor3_normal(k,t2l2)
      end do
      cor3_material(cor3_num) = white

      cor3_num = cor3_num + 1
      t2r2u = cor3_num
      do k = 1, 3
        cor3(k,t2r2u) = cor3(k,t2r2) + EPS * cor3_normal(k,t2r2)
      end do
      cor3_material(cor3_num) = white

      cor3_num = cor3_num + 1
      t1l2u = cor3_num

      norm = 0.0E+00
      do k = 1, 3
        norm = norm + ( cor3(k,b2l2) + EPS * cor3_normal(k,b2l2) &
          - cor3(k,t2l2) )**2
      end do
      norm = sqrt ( norm )
      if ( norm /= 0.0E+00 ) then
        norm = 1.0E+00 / norm
      end if

      do k = 1, 3
        cor3(k,t1l2u) = cor3(k,t2l2u) + GRID_WIDTH * ( cor3(k,b2l2) &
          + EPS * cor3_normal(k,b2l2) - cor3(k,t2l2) ) * norm
      end do

      cor3_material(cor3_num) = white

      cor3_num = cor3_num + 1
      t1r2u = cor3_num

      norm = 0.0E+00
      do k = 1, 3
        norm = norm + ( cor3(k,b2r2) + EPS * cor3_normal(k,b2r2) &
          - cor3(k,t2r2) )**2
      end do
      norm = sqrt ( norm )
      if ( norm /= 0.0E+00 ) then
        norm = 1.0E+00 / norm
      end if

      do k = 1, 3
        cor3(k,t1r2u) = cor3(k,t2r2u) + GRID_WIDTH * ( cor3(k,b2r2) &
          + EPS * cor3_normal(k,b2r2) - cor3(k,t2r2) ) * norm
      end do

      cor3_material(cor3_num) = white

      face_num = face_num + 1
      face_material(face_num) = white
      face_order(face_num) = 3
      nvert = 0

      nvert = nvert + 1
      k = t2l2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      nvert = nvert + 1
      k = t2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      nvert = nvert + 1
      k = t1l2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      face_num = face_num + 1
      face_material(face_num) = white
      face_order(face_num) = 3
      nvert = 0

      nvert = nvert + 1
      k = t1l2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      nvert = nvert + 1
      k = t2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      nvert = nvert + 1
      k = t1r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'OOGL_GRID - Added white top "grid lines".'
  write ( *, '(a,i6)' ) '  Number of points = ', cor3_num
  write ( *, '(a,i6)' ) '  Number of faces =  ', face_num
!
!  Do the grid lines that I think of as running along the "right" of
!  the affected faces.
!
!    +--- ---- T2R1U-T2R2U
!    |    ....    |    |
!    |    ....    |    |
!    +--- ---- B2R1U-B2R2U
!
!    T2L2 ---- ---- T2R2
!    |    .... ....    |
!    |    .... ....    |
!    B2L2 ---- ---- B2R2
!
  do j = grid_nx+1, ncol_oogl, grid_nx
    do i = 2, nrow_oogl

      b2l2 = base + ( i - 2 ) * ncol_oogl + j - 1
      b2r2 = base + ( i - 2 ) * ncol_oogl + j
      t2l2 = base + ( i - 1 ) * ncol_oogl + j - 1
      t2r2 = base + ( i - 1 ) * ncol_oogl + j

      cor3_num = cor3_num + 1
      b2r2u = cor3_num

      do k = 1, 3
        cor3(k,b2r2u) = cor3(k,b2r2) + EPS * cor3_normal(k,b2r2)
      end do
      cor3_material(cor3_num) = white

      cor3_num = cor3_num + 1
      t2r2u = cor3_num
      do k = 1, 3
        cor3(k,t2r2u) = cor3(k,t2r2) + EPS * cor3_normal(k,t2r2)
      end do
      cor3_material(cor3_num) = white

      cor3_num = cor3_num + 1
      t2r1u = cor3_num

      norm = 0.0E+00
      do k = 1, 3
        norm = norm + ( cor3(k,t2l2) + EPS * cor3_normal(k,t2l2) &
          - cor3(k,t2r2) )**2
      end do
      norm = sqrt ( norm )
      if ( norm /= 0.0E+00 ) then
        norm = 1.0E+00 / norm
      end if

      do k = 1, 3
        cor3(k,t2r1u) = cor3(k,t2r2u) + GRID_WIDTH * ( cor3(k,t2l2) &
          + EPS * cor3_normal(k,t2l2) - cor3(k,t2r2) ) * norm
      end do

      cor3_material(cor3_num) = white


      cor3_num = cor3_num + 1
      b2r1u = cor3_num

      norm = 0.0E+00
      do k = 1, 3
        norm = norm + ( cor3(k,b2l2) + EPS * cor3_normal(k,b2l2) &
          - cor3(k,b2r2) )**2
      end do
      norm = sqrt ( norm )
      if ( norm /= 0.0E+00 ) then
        norm = 1.0E+00 / norm
      end if

      do k = 1, 3
        cor3(k,b2r1u) = cor3(k,b2r2u) + GRID_WIDTH * ( cor3(k,b2l2) &
          + EPS * cor3_normal(k,b2l2) - cor3(k,b2r2) ) * norm
      end do

      cor3_material(cor3_num) = white



      face_num = face_num + 1
      face_material(face_num) = white
      face_order(face_num) = 3
      nvert = 0

      nvert = nvert + 1
      k = b2r1u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      nvert = nvert + 1
      k = t2r1u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      nvert = nvert + 1
      k = t2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      face_num = face_num + 1
      face_material(face_num) = white
      face_order(face_num) = 3
      nvert = 0

      nvert = nvert + 1
      k = b2r1u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      nvert = nvert + 1
      k = t2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

      nvert = nvert + 1
      k = b2r2u
      face(nvert,face_num) = k
      vertex_material(nvert,face_num) = white

    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'OOGL_GRID - Added white "grid lines".'
  write ( *, '(a,i6)' ) '  Number of points = ', cor3_num
  write ( *, '(a,i6)' ) '  Number of faces =  ', face_num

  return
end
subroutine outfile ( filein_name, fileout_name, ierror, fileout_type )

!*****************************************************************************80
!
!! OUTFILE determines the output filename and type.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Output, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Output, character ( len = 10 ) FILEOUT_TYPE, the type of the file, which is
!    set to the filename extension.  Typical values include
!    'ase', 'dxf', 'iv', 'obj', 'pov', 'ps', 'smf', 'stl', 'stla', 
!    'tec', 'tri', 'ts', 'txt', 'vla', 'wrl', 'xgl', or 'xyz'.
!
  implicit none

  character ( len = * ) filein_name
  character ( len = * ) fileout_name
  character ( len = 10 ) fileout_type
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  logical s_eqi

  ierror = 0

  if ( filein_name == ' ' ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'OUTFILE - Error!'
    write ( *, '(a)' ) '  You must read a file IN before you can'
    write ( *, '(a)' ) '  write a file OUT.'
    return
  end if
 
  if ( fileout_name == ' ' ) then
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'OUTFILE:'
    write ( *, '(a)' ) '  Enter the output filename to be created,'
    write ( *, '(a)' ) '  or hit return if done.'
 
    read ( *, '(a)', iostat = ios ) fileout_name
 
    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'OUTFILE - Error!'
      write ( *, '(a)' ) '  The output file was not specified correctly.'
      ierror = ios
      fileout_name = ' '
      return
    end if
 
  end if
!
!  Determine the output file type.
!
  call file_name_ext_get ( fileout_name, i1, i2 )

  if ( i1 /= 0 ) then
    fileout_type = fileout_name(i1:i2)
  else
    fileout_type = ' '
  end if
 
  if ( .not. ( &
    s_eqi ( fileout_type, 'ASE'  ) .or. &
    s_eqi ( fileout_type, 'BYU'  ) .or. &
    s_eqi ( fileout_type, 'DXF'  ) .or. &
    s_eqi ( fileout_type, 'HRC'  ) .or. &
    s_eqi ( fileout_type, 'IV'   ) .or. &
    s_eqi ( fileout_type, 'OBJ'  ) .or. &
    s_eqi ( fileout_type, 'OFF'  ) .or. &
    s_eqi ( fileout_type, 'POV'  ) .or. s_eqi ( fileout_type, 'PS'   ) .or. &
    s_eqi ( fileout_type, 'SMF'  ) .or. s_eqi ( fileout_type, 'STL'  ) .or. &
    s_eqi ( fileout_type, 'STLA' ) .or. s_eqi ( fileout_type, 'TEC'  ) .or. &
    s_eqi ( fileout_type, 'TRI'  ) .or. s_eqi ( fileout_type, 'TRIA' ) .or. &
    s_eqi ( fileout_type, 'TS'   ) .or. &
    s_eqi ( fileout_type, 'TXT'  ) .or. s_eqi ( fileout_type, 'UCD'  ) .or. &
    s_eqi ( fileout_type, 'VLA'  ) .or. s_eqi ( fileout_type, 'WRL'  ) .or. &
    s_eqi ( fileout_type, 'XGL'  ) .or. s_eqi ( fileout_type, 'XYZ' ) ) ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'OutFile could not determine the file type.'
    write ( *, '(a)' ) '  The output file name is:'
    write ( *, '(a)' ) trim ( fileout_name )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The file type should occur after the period.'
    write ( *, '(a)' ) '  Please specify the file type you are using:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  "ase", "byu", "dxf", "hrc", "iv",'
    write ( *, '(a)' ) '  "obj", "off", "pov", "ps",  "smf",'
    write ( *, '(a)' ) '  "stl", "tec", "tri", "ts",  "txt",'
    write ( *, '(a)' ) '  "ucd", "vla", "wrl", "xgl", "xyz":'

    read ( *, '(a)' ) fileout_type
    call s_cap ( fileout_type )

    if ( .not. ( &
      s_eqi ( fileout_type, 'ASE'  ) .or. s_eqi ( fileout_type, 'BYU'  ) .or. &
      s_eqi ( fileout_type, 'DXF'  ) .or. s_eqi ( fileout_type, 'HRC'  ) .or. &
      s_eqi ( fileout_type, 'IV'   ) .or. s_eqi ( fileout_type, 'OBJ'  ) .or. &
      s_eqi ( fileout_type, 'OFF'  ) .or. &
      s_eqi ( fileout_type, 'POV'  ) .or. s_eqi ( fileout_type, 'PS'   ) .or. &
      s_eqi ( fileout_type, 'SMF'  ) .or. s_eqi ( fileout_type, 'STL'  ) .or. &
      s_eqi ( fileout_type, 'STLA' ) .or. s_eqi ( fileout_type, 'TEC'  ) .or. &
      s_eqi ( fileout_type, 'TRI'  ) .or. s_eqi ( fileout_type, 'TRIA' ) .or. &
      s_eqi ( fileout_type, 'TS'   ) .or. &
      s_eqi ( fileout_type, 'TXT'  ) .or. s_eqi ( fileout_type, 'UCD'  ) .or. &
      s_eqi ( fileout_type, 'VLA'  ) .or. s_eqi ( fileout_type, 'WRL'  ) .or. &
      s_eqi ( fileout_type, 'XGL'  ) .or. s_eqi ( fileout_type, 'XYZ' ) ) ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'OUTFILE - Error!'
      write ( *, '(a)' ) '  The file type was not acceptable!'
      return
    end if

  end if

  return
end
function pi ( )

!*****************************************************************************80
!
!! PI returns the value of pi.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 4 ) PI, the value of pi.
!
  implicit none

  real ( kind = 4 ) pi

  pi = 3.141592653589793E+00

  return
end
subroutine plane_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d )

!*****************************************************************************80
!
!! PLANE_EXP2IMP_3D converts an explicit plane to implicit form in 3D.
!
!  Discussion:
!
!    The explicit form of a plane in 3D is
!
!      (X1,Y1,Z1), (X2,Y2,Z2), (X3,Y3,Z3).
!
!    The implicit form of a plane in 3D is
!
!      A * X + B * Y + C * Z + D = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X1, Y1, Z1, X2, Y2, X2, X3, Y3, Z3, are three points 
!    on the plane, which must be distinct, and not collinear.
!
!    Output, real ( kind = 4 ) A, B, C, D, coefficients which describe the plane.
!
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  real ( kind = 4 ) c
  real ( kind = 4 ) d
  real ( kind = 4 ) x1
  real ( kind = 4 ) x2
  real ( kind = 4 ) x3
  real ( kind = 4 ) y1
  real ( kind = 4 ) y2
  real ( kind = 4 ) y3
  real ( kind = 4 ) z1
  real ( kind = 4 ) z2
  real ( kind = 4 ) z3

  a = ( y2 - y1 ) * ( z3 - z1 ) - ( z2 - z1 ) * ( y3 - y1 )
  b = ( z2 - z1 ) * ( x3 - x1 ) - ( x2 - x1 ) * ( z3 - z1 )
  c = ( x2 - x1 ) * ( y3 - y1 ) - ( y2 - y1 ) * ( x3 - x1 ) 
  d = - x2 * a - y2 * b - z2 * c 

  return
end
subroutine plane_imp_point_nearest_3d ( a, b, c, d, x, y, z, xn, yn, zn )

!*****************************************************************************80
!
!! PLANE_IMP_POINT_NEAREST_3D: nearest point on a implicit plane to a point in 3D.
!
!  Discussion:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A, B, C, D, coefficients that define the plane as
!    the set of points for which A*X+B*Y+C*Z+D = 0.
!
!    Input, real ( kind = 4 ) X, Y, Z, the coordinates of the point.
!
!    Output, real ( kind = 4 ) XN, YN, ZN, the coordinates of the nearest point on
!    the plane.
! 
  implicit none

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  real ( kind = 4 ) c
  real ( kind = 4 ) d
  real ( kind = 4 ) t
  real ( kind = 4 ) x
  real ( kind = 4 ) xn
  real ( kind = 4 ) y
  real ( kind = 4 ) yn
  real ( kind = 4 ) z
  real ( kind = 4 ) zn

  if ( a == 0.0E+00 .and. b == 0.0E+00 .and. c == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_POINT_NEAREST_3D - Fatal error!'
    write ( *, '(a)' ) '  A = B = C = 0.'
    stop
  end if
!
!  The normal N to the plane is (A,B,C).
!
!  The line defined by (XN-X)/A = (YN-Y)/B = (ZN-Z)/C = T
!  goes through (X,Y,Z) and is parallel to N.
!
!  Solving for the point (XN,YN,ZN) we get
!
!    XN = A*T+X
!    YN = B*T+Y
!    ZN = C*T+Z
!
!  Now place these values in the equation for the plane:
!
!    A*(A*T+X) + B*(B*T+Y) + C*(C*T+Z) + D = 0
!
!  and solve for T:
!
!    T = (-A*X-B*Y-C*Z-D) / (A * A + B * B + C * C )
!
  t = - ( a * x + b * y + c * z + d ) / ( a * a + b * b + c * c )
 
  xn = x + a * t
  yn = y + b * t
  zn = z + c * t
 
  return
end
subroutine points_distance_3d ( dis, x1, y1, z1, x2, y2, z2 )

!*****************************************************************************80
!
!! POINTS_DISTANCE_3D finds the distance between two points in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 4 ) DIS, the distance between the points.
!
!    Input, real ( kind = 4 ) X1, Y1, Z1, X2, Y2, Z2, determines the pair of points
!    (X1,Y1,Z1) and (X2,Y2,Z2) whose distance apart is be determined.
!
  implicit none

  real ( kind = 4 ) dis
  real ( kind = 4 ) x1
  real ( kind = 4 ) x2
  real ( kind = 4 ) y1
  real ( kind = 4 ) y2
  real ( kind = 4 ) z1
  real ( kind = 4 ) z2

  dis = sqrt ( ( x1 - x2 )**2 + ( y1 - y2 )**2 + ( z1 - z2 )**2 )
 
  return
end
subroutine poly_2_tri ( face, face_material, face_max, face_num, face_order, &
  ierror, order_max, vertex_material )

!*****************************************************************************80
!
!! POLY_2_TRI converts a collection of polygons into a collection of triangles.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    1, the algorithm failed because FACE_MAX was too small.
!    2, the algorithm failed because there were faces of order < 3.
!    3, the algorithm failed because there were faces of order > ORDER_MAX.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces allowed.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices allowed per face.
!
!    Input/output, integer ( kind = 4 ) FACE_NUM, the number of faces.  This value is updated
!    on return.
!
!    Input/output, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  implicit none

  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_num2
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) iface_old
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) k
  integer ( kind = 4 ) vertex_material(order_max,face_max)

  ierror = 0
  face_num2 = 0

  do iface = 1, face_num

    if ( face_order(iface) < 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLY_2_TRI - Fatal error!'
      write ( *, '(a,i6,a)' ) '  Face ', iface, ' is illegal.'
      write ( *, '(a,i6)' ) '  Number of vertices is ', face_order(iface)
      ierror = 2
      return
    else if ( order_max < face_order(iface) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLY_2_TRI - Fatal error!'
      write ( *, '(a,i6,a)' ) '  Face ', iface, ' is illegal.'
      write ( *, '(a,i6)' ) '  Number of vertices is ', face_order(iface)
      write ( *, '(a,i6)' ) '  ORDER_MAX is ', order_max
      ierror = 3
      return
    end if

    do ivert = 3, face_order(iface)
      face_num2 = face_num2 + 1
    end do

  end do

  if ( face_max < face_num2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLY_2_TRI - Fatal error!'
    write ( *, '(a)' ) '  FACE_MAX is too small to replace all faces'
    write ( *, '(a)' ) '  by triangles.'
    write ( *, '(a,i6)' ) '  FACE_MAX = ', face_max
    write ( *, '(a,i6)' ) '  FACE_NUM2 = ', face_num2
    ierror = 1
    return
  end if

  iface_old = face_num
  k = face_order(iface_old)

  do iface = face_num2, 1, -1

    if ( k < 3 ) then
      iface_old = iface_old - 1
      k = face_order(iface_old)
    end if

    face_material(iface) = face_material(iface_old)
    face_order(iface) = 3
    face(1,iface) = face(1,iface_old)
    vertex_material(1,iface) = vertex_material(1,iface_old)
    do ivert = 2, 3
      face(ivert,iface) = face(k+ivert-3,iface_old)
      vertex_material(ivert,iface) = vertex_material(k+ivert-3,iface_old)
    end do

    k = k - 1

  end do

  face_num = face_num2

  return
end
subroutine pov_write ( cor3, cor3_max, face, face_material, face_max, &
  face_num, face_order, filein_name, fileout_name, iunit, &
  material_max, material_num, material_rgba, order_max, vertex_normal )

!*****************************************************************************80
!
!! POV_WRITE writes graphics information to a POV file.
!
!  Example:
!
!    // cone.pov created by IVREAD.
!    // Original data in cone.iv.
!
!    #version 3.0E+00
!    #include "colors.inc"
!    #include "shapes.inc"
!    global_settings { assumed_gamma 2.2 }
!
!    camera { 
!     right < 4/3, 0, 0>
!     up < 0, 1, 0 >
!     sky < 0, 1, 0 >
!     angle 20
!     location < 0, 0, -300 >
!     look_at < 0, 0, 0>
!    }
!
!    light_source { < 20, 50, -100 > color White }
!
!    background { color SkyBlue }
!
!    #declare Material001 = texture { 
!      pigment { color rgb < 0.8, 0.2, 0.2> } 
!      finish { ambient 0.2 diffuse 0.5 }
!    }
! 
!    #declare Material002 = texture { 
!      pigment { color rgb < 0.2, 0.2, 0.8> } 
!      finish { ambient 0.2 diffuse 0.5 }
!    }
!
!    mesh {
!      smooth_triangle { 
!        < 0.29, -0.29, 0.0>, < 0.0, 0.0, -1.0 >,
!        < 38.85, 10.03, 0.0>, < 0.0, 0.0, -1.0 >,
!        < 40.21, -0.29, 0.0>, <  0.0, 0.0, -1.0 > 
!        texture { Material002 } }
!        ...
!      smooth_triangle { 
!        <  0.29, -0.29, 70.4142 >, < 0.0,  0.0, 1.0 >,
!        <  8.56,  -2.51, 70.4142 >, < 0.0,  0.0, 1.0 >,
!        <  8.85, -0.29, 70.4142 >, < 0.0,  0.0, 1.0 > 
!        texture { Material001 } }
!    }
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), face materials.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, real ( kind = 4 ) MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and A values.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) b
  character ( len = 4 ) char4
  character comma
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  character ( len = * ) fileout_name
  real ( kind = 4 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_mat
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) r
  character ( len = 100 ) text
  integer ( kind = 4 ) text_num
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)

  text_num = 0

  write ( iunit, '(a)' ) '// ' // trim ( fileout_name ) // ' created by IVREAD.'
  write ( iunit, '(a)' ) '// Original data in ' // trim ( filein_name ) // '.'

  text_num = text_num + 2
!
!  Initial declarations.
!
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) '#version 3.0'
  write ( iunit, '(a)' ) '#include "colors.inc"'
  write ( iunit, '(a)' ) '#include "shapes.inc"'
  write ( iunit, '(a)' ) 'global_settings { assumed_gamma 2.2 }'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) 'camera { '
  write ( iunit, '(a)' ) ' right < 4/3, 0, 0>'
  write ( iunit, '(a)' ) ' up < 0, 1, 0 >'
  write ( iunit, '(a)' ) ' sky < 0, 1, 0 >'
  write ( iunit, '(a)' ) ' angle 20'
  write ( iunit, '(a)' ) ' location < 0, 0, -300 >'
  write ( iunit, '(a)' ) ' look_at < 0, 0, 0>'
  write ( iunit, '(a)' ) '}'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) 'light_source { < 20, 50, -100 > color White }'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) 'background { color SkyBlue }'

  text_num = text_num + 15
!
!  Declare RGB textures.
!
  do i = 1, material_num

    write ( iunit, '(a)' ) ' '

    call i4_to_s_zero ( i, char4 )
    write ( iunit, '(a)' ) '#declare Material_' // char4 // ' = texture { '

    r = material_rgba(1,i)
    g = material_rgba(2,i)
    b = material_rgba(3,i)
    write ( iunit, '(a,f4.2,a,f4.2,a,f4.2,a)' ) '  pigment { color rgb < ', &
      r, ',', g, ',', b ,' > } '

    write ( iunit, '(a)' ) '  finish { ambient 0.2 diffuse 0.5 }'

    write ( iunit, '(a)' ) '}'

  end do
!
!  Write one big object.
!
  write ( iunit, '(a)' ) 'mesh {'
  text_num = text_num + 1
!
!  Do the next face.
!
  do i = 1, face_num
!
!  Break the face up into triangles, anchored at node 1.
!
    do jlo = 1, face_order(i) - 2

      write ( iunit, '(a)' )  '  smooth_triangle { '
      text_num = text_num + 1

      do j = jlo, jlo + 2

        if ( j == jlo ) then
          jj = 1
        else
          jj = j
        end if

        k = face(jj,i)

        if ( j < jlo + 2 ) then
          comma = ','
        else
          comma = ' '
        end if
 
        write ( text, '(a,3(f10.3,a),3(f6.2,a),a )' ) &
          '<', cor3(1,k), ',', cor3(2,k), ',', cor3(3,k), '>, <', &
          vertex_normal(1,jj,i), ',', vertex_normal(2,jj,i), ',', &
          vertex_normal(3,jj,i), '>', comma

        call s_blanks_delete ( text )
        write ( iunit, '(a)' ) trim ( text )
        text_num = text_num + 1

      end do

      i_mat = face_material(i)
      call i4_to_s_zero ( i_mat, char4 )
      write ( iunit, '(a)' ) 'texture { Material_' // char4 // ' } }'
      text_num = text_num + 1

    end do

  end do

  write ( iunit, '(a)' ) '}'
  text_num = text_num + 1
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'POV_WRITE - Wrote ', text_num, ' text lines to ' &
    // trim ( fileout_name )

  return
end
subroutine project_2d ( cor2, cor3, ierror, cor2_max, cor3_max, cor2_num, &
  cor3_num )

!*****************************************************************************80
!
!! PROJECT_2D projects 3D data to 2D based on user choices.
!
!  Discussion:
!
!    Projections include:
!
!      drop X coordinate, display YZ;
!      drop Y coordinate, display XZ;
!      drop Z coordinate, display XY;
!      orthographic projection into a plane;
!      perspective projection into a plane through a focus point;
!      project X into YZ using an angle THETA;
!      project Y into XZ using an angle THETA;
!      project Z into XY using an angle THETA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 4 ) COR2(2,COR2_MAX), the projected 2D data.
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the data to project.
!
!    Input, integer ( kind = 4 ) COR2_MAX, COR3_MAX, the maximum number of 2D and
!    3D points that can be handled.
!
!    Output, integer ( kind = 4 ) COR2_NUM, the number of 2D points.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of 3D points.
!
  implicit none

  integer ( kind = 4 ) cor2_max
  integer ( kind = 4 ) cor3_max

  real ( kind = 4 ) cor2(2,cor2_max)
  integer ( kind = 4 ) cor2_num
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 20 ) isay
  logical s_eqi
  real ( kind = 4 ) theta
  real ( kind = 4 ) x1
  real ( kind = 4 ) x2
  real ( kind = 4 ) x3
  real ( kind = 4 ) xf
  real ( kind = 4 ) y1
  real ( kind = 4 ) y2
  real ( kind = 4 ) y3
  real ( kind = 4 ) yf
  real ( kind = 4 ) z1
  real ( kind = 4 ) z2
  real ( kind = 4 ) z3
  real ( kind = 4 ) zf

  ierror = 0

  if ( cor3_num <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PROJECT_2D - Fatal error!'
    write ( *, '(a)' ) '  Input COR3_NUM <= 0.'
    ierror = 1
    return
  end if

  cor2_num = cor3_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Choose a projection from 3D -> 2D:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '-X     drop X coordinate, display YZ;'
  write ( *, '(a)' ) '-Y     drop Y coordinate, display XZ;'
  write ( *, '(a)' ) '-Z     drop Z coordinate, display XY;'
  write ( *, '(a)' ) 'OPLANE orthographic projection into plane.'
  write ( *, '(a)' ) 'PPLANE perspective projection into plane.'
  write ( *, '(a)' ) 'PX     project X into YZ using THETA;'
  write ( *, '(a)' ) 'PY     project Y into XZ using THETA;'
  write ( *, '(a)' ) 'PZ     project Z into XY using THETA;'
 
  read ( *, '(a)', iostat = ios ) isay

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PROJECT_2D - Unexpected end of input.'
    ierror = ios
    return
  end if
 
  if ( s_eqi ( isay, '-X' ) ) then
 
    do i = 1, cor3_num
      cor2(1,i) = cor3(2,i)
      cor2(2,i) = cor3(3,i)
    end do
 
  else if ( s_eqi ( isay, '-Y' ) ) then
 
    do i = 1, cor3_num
      cor2(1,i) = cor3(1,i)
      cor2(2,i) = cor3(3,i)
    end do
 
  else if ( s_eqi ( isay, '-Z' ) ) then
 
    do i = 1, cor3_num
      cor2(1,i) = cor3(1,i)
      cor2(2,i) = cor3(2,i)
    end do
 
  else if ( s_eqi ( isay, 'OPLANE' ) ) then
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter 3 (X,Y,Z) points on the plane:'
    read ( *, *, iostat = ios ) x1, y1, z1, x2, y2, z2, x3, y3, z3

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PROJECT_2D - Unexpected end of input.'
      ierror = ios
      return
    end if
 
    call project_oplane ( x1, y1, z1, x2, y2, z2, x3, y3, z3, &
      cor2, cor3, cor2_max, cor3_max, cor3_num )

  else if ( s_eqi ( isay, 'PPLANE' ) ) then
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter 3 (X,Y,Z) points on the plane:'
    read ( *, *, iostat = ios ) x1, y1, z1, x2, y2, z2, x3, y3, z3

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PROJECT_2D - Unexpected end of input.'
      ierror = ios
      return
    end if

    write ( *, '(a)' ) 'Enter focus point (X,Y,Z):'
    read ( *, *, iostat = ios ) xf, yf, zf

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PROJECT_2D - Unexpected end of input.'
      ierror = ios
      return
    end if
 
    call project_pplane ( x1, y1, z1, x2, y2, z2, x3, y3, z3, xf, yf, zf, &
      cor2, cor3, cor2_max, cor3_max, cor2_num, cor3_num )
 
  else if ( s_eqi ( isay, 'PX' ) ) then
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter projection angle THETA in degrees:'
    read ( *, *, iostat = ios ) theta

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PROJECT_2D - Unexpected end of input.'
      ierror = ios
      return
    end if

    call project_angle ( 'X', cor2, cor3, cor2_max, cor3_max, cor2_num, theta )

  else if ( s_eqi ( isay, 'PY' ) ) then
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter projection angle THETA in degrees:'
    read ( *, *, iostat = ios ) theta

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PROJECT_2D - Unexpected end of input.'
      ierror = ios
      return
    end if

    call project_angle ( 'Y', cor2, cor3, cor2_max, cor3_max, cor2_num, theta )

  else if ( s_eqi ( isay, 'PZ' ) ) then
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter projection angle THETA in degrees:'
    read ( *, *, iostat = ios ) theta

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PROJECT_2D - Unexpected end of input.'
      ierror = ios
      return
    end if

    call project_angle ( 'Z', cor2, cor3, cor2_max, cor3_max, cor2_num, theta )

  else
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PROJECT_2D - Error!'
    write ( *, '(a)' ) '  Unrecognized projection option!'
    cor2_num = 0
    ierror = 1
    
  end if
 
  return
end
subroutine project_angle ( cor, cor2, cor3, cor2_max, cor3_max, cor2_num, &
  theta )

!*****************************************************************************80
!
!! PROJECT_ANGLE converts 3D data to 2D using a presentation angle.
!
!  Discussion:
!
!    A "presentation angle" THETA is used to project the 3D point
!    (X3D, Y3D, Z3D) to the 2D projection (XVAL,YVAL).
!
!    The formula used if COR = 'X' is
!
!      X2D = Y3D - sin ( THETA ) * X3D
!      Y2D = Z3D - sin ( THETA ) * X3D
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character COR, the coordinate to be projected.
!    COR should be 'X', 'Y', or 'Z'.
!
!    Output, real ( kind = 4 ) COR2(2,COR2_MAX), the 2D projections.
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the 3D points to be projected.
!
!    Input, integer ( kind = 4 ) COR2_MAX, COR3_MAX, the maximum number of 2D and
!    3D points allowed.
!
!    Input, integer ( kind = 4 ) COR2_NUM, the number of 2D values to be computed.
!
!    Input, real ( kind = 4 ) THETA, the presentation angle in degrees.
!
  implicit none

  integer ( kind = 4 ) cor2_max
  integer ( kind = 4 ) cor3_max

  character cor
  real ( kind = 4 ) cor2(2,cor2_max)
  integer ( kind = 4 ) cor2_num
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) i
  real ( kind = 4 ) pi
  logical s_eqi
  real ( kind = 4 ) stheta
  real ( kind = 4 ) theta

  stheta = sin ( pi() * theta / 180.0E+00 )

  if ( cor2_max < cor2_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PROJECT_ANGLE - Fatal error!'
    write ( *, '(a)' ) '  COR2_NUM is greater than COR2_MAX.'
    stop
  end if

  if (cor3_max < cor2_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PROJECT_ANGLE - Fatal error!'
    write ( *, '(a)' ) '  COR2_NUM is greater than COR3_MAX.'
    stop
  end if

  if ( s_eqi ( cor, 'X' ) ) then

    do i = 1, cor2_num
      cor2(1,i) = cor3(2,i) - stheta * cor3(1,i)
      cor2(2,i) = cor3(3,i) - stheta * cor3(1,i)
    end do

  else if ( s_eqi ( cor, 'Y' ) ) then

    do i = 1, cor2_num
      cor2(1,i) = cor3(1,i) - stheta * cor3(2,i)
      cor2(2,i) = cor3(3,i) - stheta * cor3(2,i)
    end do

  else if ( s_eqi ( cor, 'Z' ) ) then

    cor2(1,1:cor2_num) = cor3(1,1:cor2_num) - stheta * cor3(3,1:cor2_num)
    cor2(2,1:cor2_num) = cor3(2,1:cor2_num) - stheta * cor3(3,1:cor2_num)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PROJECT_ANGLE - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized axis.'
    stop

  end if

  return
end
subroutine project_oplane ( x1, y1, z1, x2, y2, z2, x3, y3, z3, cor2, cor3, &
  cor2_max, cor3_max, cor3_num )

!*****************************************************************************80
!
!! PROJECT_OPLANE projects 3D points onto an orthographic plane.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are the 
!    coordinates of three points on the plane.
!
!    Output, real ( kind = 4 ) COR2(2,COR2_MAX), the "local" in-plane coordinates
!    of the projections of the object points.
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the (X,Y,Z) coordinates of the object 
!     points.
!
!    Input, integer ( kind = 4 ) COR2_MAX, COR3_MAX, the maximum number of 2D and
!    3D points that can be handled.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points to project.
!
  implicit none

  integer ( kind = 4 ) cor2_max
  integer ( kind = 4 ) cor3_max

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  real ( kind = 4 ) c
  real ( kind = 4 ) cor2(2,cor2_max)
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  real ( kind = 4 ) d
  real ( kind = 4 ) dot
  integer ( kind = 4 ) i
  real ( kind = 4 ) v1(3)
  real ( kind = 4 ) v2(3)
  real ( kind = 4 ) x1
  real ( kind = 4 ) x2
  real ( kind = 4 ) x3
  real ( kind = 4 ) xn
  real ( kind = 4 ) y1
  real ( kind = 4 ) y2
  real ( kind = 4 ) y3
  real ( kind = 4 ) yn
  real ( kind = 4 ) z1
  real ( kind = 4 ) z2
  real ( kind = 4 ) z3
  real ( kind = 4 ) zn
!
!  Put the plane into ABCD form.
!
  call plane_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d )
!
!  Compute the two basis vectors for the affine plane.
!
  v1(1) = x2 - x1
  v1(2) = y2 - y1
  v1(3) = z2 - z1

  call vector_unit_nd ( 3, v1 )

  v2(1) = x3 - x1
  v2(2) = y3 - y1
  v2(3) = z3 - z1

  dot = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)

  v2(1:3) = v2(1:3) - dot * v1(1:3)

  call vector_unit_nd ( 3, v2 )
!
!  For each point, its image in the plane is the nearest point 
!  in the plane.
!
  do i = 1, min ( cor3_num, cor3_max )
 
    call plane_imp_point_nearest_3d ( a, b, c, d, cor3(1,i), cor3(2,i), &
      cor3(3,i), xn, yn, zn )

    cor2(1,i) = ( xn - x1 ) * v1(1) + ( yn - y1 ) * v1(2) + ( zn - z1 ) * v1(3)
    cor2(2,i) = ( xn - x1 ) * v2(1) + ( yn - y1 ) * v2(2) + ( zn - z1 ) * v2(3)
 
  end do

  return
end
subroutine project_pplane ( x1, y1, z1, x2, y2, z2, x3, y3, z3, xf, yf, zf, &
  cor2, cor3, cor2_max, cor3_max, cor2_num, cor3_num )

!*****************************************************************************80
!
!! PROJECT_PPLANE projects a point through a focus point onto a perspective plane.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, are the 
!    coordinates of three points on the plane.
!
!    Input, real ( kind = 4 ) XF, YF, ZF, are the coordinates of the focus point.
!
!    Output, real ( kind = 4 ) COR2(2,COR2_MAX), the "local" in-plane coordinates
!    of the projections of the object points.
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the (X,Y,Z) coordinates of points
!    to be projected.
!
!    Input, integer ( kind = 4 ) COR2_MAX, COR3_MAX, the maximum number of 2D and
!    3D points that can be handled.
!
!    Output, integer ( kind = 4 ) COR2_NUM, the number of projected points.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points to project.
!
  implicit none

  integer ( kind = 4 ) cor2_max
  integer ( kind = 4 ) cor3_max

  real ( kind = 4 ) a
  real ( kind = 4 ) alpha
  real ( kind = 4 ) angle_rad_3d
  real ( kind = 4 ) b
  real ( kind = 4 ) beta
  real ( kind = 4 ) c
  real ( kind = 4 ) cor2(2,cor2_max)
  integer ( kind = 4 ) cor2_num
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  real ( kind = 4 ) d
  real ( kind = 4 ) disfo
  real ( kind = 4 ) disfn
  real ( kind = 4 ) dot
  integer ( kind = 4 ) i
  real ( kind = 4 ) v1(3)
  real ( kind = 4 ) v2(3)
  real ( kind = 4 ) x1
  real ( kind = 4 ) x2
  real ( kind = 4 ) x3
  real ( kind = 4 ) xf
  real ( kind = 4 ) xn
  real ( kind = 4 ) xp
  real ( kind = 4 ) y1
  real ( kind = 4 ) y2
  real ( kind = 4 ) y3
  real ( kind = 4 ) yf
  real ( kind = 4 ) yn
  real ( kind = 4 ) yp
  real ( kind = 4 ) z1
  real ( kind = 4 ) z2
  real ( kind = 4 ) z3
  real ( kind = 4 ) zf
  real ( kind = 4 ) zn
  real ( kind = 4 ) zp
!
!  Put the plane into ABCD form.
!
  call plane_exp2imp_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, a, b, c, d )
!
!  Get the nearest point on the plane to the focus.
!
  call plane_imp_point_nearest_3d ( a, b, c, d, xf, yf, zf, xn, yn, zn )
!
!  Get the distance from the focus to the plane.
!
  call points_distance_3d ( disfn, xf, yf, zf, xn, yn, zn )
!
!  If the focus lies in the plane, this is bad.  We could still
!  project points that actually lie in the plane, but we'll
!  just bail out.
!
  if ( disfn == 0.0E+00 ) then

    cor2_num = 0
    return

  end if
!
!  Compute the two basis vectors for the affine plane.
!
  v1(1) = x2 - x1
  v1(2) = y2 - y1
  v1(3) = z2 - z1

  call vector_unit_nd ( 3, v1 )

  v2(1) = x3 - x1
  v2(2) = y3 - y1
  v2(3) = z3 - z1

  dot = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)

  v2(1:3) = v2(1:3) - dot * v1(1:3)

  call vector_unit_nd ( 3, v2 )
!
!  Process the points.
!
  do i = 1, min ( cor3_num, cor3_max )
!
!  Get the distance from the focus to the object.
!
    call points_distance_3d ( disfo, xf, yf, zf, cor3(1,i), cor3(2,i), &
      cor3(3,i) )
 
    if ( disfo == 0.0E+00 ) then
 
      xp = xn
      yp = yn
      zp = zn 

    else
!
!  Compute ALPHA, the angle between (OBJECT-FOCUS) and (NEAREST-FOCUS).
!
      alpha = angle_rad_3d ( cor3(1,i), cor3(2,i), cor3(3,i), &
        xf, yf, zf, xn, yn, zn )
 
      if ( cos ( alpha ) == 0.0E+00 ) then
 
        xp = xn
        yp = yn
        zp = zn
 
      else
!
!  Multiplier BETA is Dist(NEAREST-FOCUS) / ( Cos(ALPHA)*Dist(OBJECT-FOCUS) )
!
        beta = disfn / ( cos ( alpha ) * disfo )
!
!  Set the projected point.
!
        xp = xf + beta * ( cor3(1,i) - xf )
        yp = yf + beta * ( cor3(2,i) - yf )
        zp = zf + beta * ( cor3(3,i) - zf )
 
      end if
 
    end if
 
    cor2(1,i) = ( xp - x1 ) * v1(1) + ( yp - y1 ) * v1(2) + ( zp - z1 ) * v1(3)
    cor2(2,i) = ( xp - x1 ) * v2(1) + ( yp - y1 ) * v2(2) + ( zp - z1 ) * v2(3)

  end do
 
  return
end
subroutine ps_write ( cor2, cor2_max, cor2_num, face, face_material, &
  face_max, face_num, face_order, fileout_name, iunit, line_dex, &
  line_material, line_max, line_num, material_max, &
  material_rgba, order_max )

!*****************************************************************************80
!
!! PS_WRITE writes 2D face and line information to a PostScript file.
!
!  Discussion:
!
!    The intent is that a 3D model will be projected in some way
!    to a 2D model that can be printed out as a standard PostScript object.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR2(2,COR2_MAX), the X and Y components of 2D points.
!
!    Input, integer ( kind = 4 ) COR2_MAX, the maximum number of 2D points.
!
!    Input, integer ( kind = 4 ) COR2_NUM, the number of 2D points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces defined.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which data is written.
!
!    Input, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer ( kind = 4 ) LINE_MATERIAL(LINE_MAX), material index for line.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, real ( kind = 4 ) MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and A values.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
  implicit none

  integer ( kind = 4 ) cor2_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) alpha
  real ( kind = 4 ) blue
  real ( kind = 4 ) cor2(2,cor2_max)
  integer ( kind = 4 ) cor2_num
  character ( len = 8 ) date
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) fileout_name
  real ( kind = 4 ) green
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) imat
  integer ( kind = 4 ) imat_old
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) k
  logical lineopen
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) margin
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ) page_num
  integer ( kind = 4 ) pagexmax
  integer ( kind = 4 ) pagexmin
  integer ( kind = 4 ) pageymax
  integer ( kind = 4 ) pageymin
  integer ( kind = 4 ) plotxmax
  integer ( kind = 4 ) plotxmin
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymax
  integer ( kind = 4 ) plotymin
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px
  integer ( kind = 4 ) py
  real ( kind = 4 ) red
  integer ( kind = 4 ) text_num
  character ( len = 10 ) time
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin
  real ( kind = 4 ) ymax
  real ( kind = 4 ) ymin
!
!  Initialization
!
  imat_old = -999
  page_num = 1
  text_num = 0

  call date_and_time ( date, time )
!
!  Compute the bounding box.
!
  xmin = minval ( cor2(1,1:cor2_num) )
  xmax = maxval ( cor2(1,1:cor2_num) )
  ymin = minval ( cor2(2,1:cor2_num) )
  ymax = maxval ( cor2(2,1:cor2_num) )

  if ( xmin == xmax ) then
    xmin = cor2(1,1) - 0.5E+00
    xmax = cor2(1,1) + 0.5E+00
  end if

  if ( ymin == ymax ) then
    ymin = cor2(2,1) - 0.5E+00
    ymax = cor2(2,1) + 0.5E+00
  end if
!
!  Compute the scale factor.
!
  pagexmax = 612
  pagexmin = 0
  pageymax = 792
  pageymin = 0

  margin = 36

  plotxmax = pagexmax - margin
  plotxmin = pagexmin + margin
  plotymax = pageymax - margin
  plotymin = pageymin + margin

  alpha = min ( ( plotxmax - plotxmin ) / ( xmax - xmin ), &
                ( plotymax - plotymin ) / ( ymax - ymin ) )
!
!  Adjust PLOTXMIN and PLOTYMIN to center the image.
!
  plotxmin2 = 0.5E+00 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) )
  plotymin2 = 0.5E+00 * ( plotymin + plotymax - alpha * ( ymax - ymin ) )
!
!  Prolog
!
  write ( iunit, '(a)' ) '%!PS-Adobe-3.0'
  write ( iunit, '(a)' ) '%%Creator: ivread.f90'
  write ( iunit, '(a)' ) '%%Title: ' // trim ( fileout_name )
  write ( iunit, '(a,4i5)' ) '%%BoundingBox', plotxmin, plotymin, plotxmax, &
    plotymax
  write ( iunit, '(a)' ) '%%CreationDate: ' // date // '  ' // time
  write ( iunit, '(a)' ) '%%Pages: 1'
  write ( iunit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( iunit, '(a)' ) '%%LanguageLevel: 2'
  write ( iunit, '(a)' ) '%%EndComments'
  write ( iunit, '(a)' ) '%%BeginProlog'
  write ( iunit, '(a)' ) '%%EndProlog'
  write ( iunit, '(a)' ) '/inch {72 mul} def'
  write ( iunit, '(a)' ) '%%Page: 1 1'
  write ( iunit, '(a)' ) 'save'

  text_num = text_num + 14
!
!  Fill the faces.
!
  red = 0.7E+00
  green = 0.7E+00
  blue = 0.0E+00
  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'

  text_num = text_num + 1

  do iface = 1, face_num

    imat = face_material(iface)

    if ( imat /= imat_old ) then
      imat_old = imat
      red = material_rgba(1,imat)
      green = material_rgba(2,imat)
      blue = material_rgba(3,imat)
      write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
    end if

    jhi = face_order(iface)

    do j = 1, jhi + 1

      if ( j <= face_order(iface) ) then
        k = face(j,iface)
      else
        k = face(1,iface)
      end if

      px = plotxmin2 + nint ( alpha * ( cor2(1,k) - xmin ) )
      py = plotymin2 + nint ( alpha * ( cor2(2,k) - ymin ) )

      if ( j == 1 ) then
        write ( iunit, '(a)' ) ' newpath'
        write ( iunit, '(2i4,a)' ) px, py, ' moveto'
        text_num = text_num + 2
      else
        write ( iunit, '(2i4,a)' ) px, py, ' lineto'
        text_num = text_num + 1
      end if

    end do

    write ( iunit, '(a)' ) ' fill'
    text_num = text_num + 1

  end do
!
!  Draw the boundaries of the faces as black lines.
!
  red = 0.0E+00
  green = 0.0E+00
  blue = 0.0E+00
  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
  text_num = text_num + 1

  do iface = 1, face_num

    jhi = face_order(iface)

    do j = 1, jhi + 1

      if ( j <= face_order(iface) ) then
        k = face(j,iface)
      else
        k = face(1,iface)
      end if

      px = plotxmin2 + nint ( alpha * ( cor2(1,k) - xmin ) )
      py = plotymin2 + nint ( alpha * ( cor2(2,k) - ymin ) )

      if ( j == 1 ) then
        write ( iunit, '(a)' ) ' newpath'
        write ( iunit, '(2i4,a)' ) px, py, ' moveto'
        text_num = text_num + 2
      else
        write ( iunit, '(2i4,a)' ) px, py, ' lineto'
        text_num = text_num + 1
      end if

    end do

    write ( iunit, '(a)' ) ' stroke'
    text_num = text_num + 1

  end do
!
!  Draw other lines.
!
!  We need to set the color of the lines as specified by the user 
!  using LINE_MAT.
!
  lineopen = .false.

  do i = 1, line_num

    j = line_dex(i)

    if ( j <= 0 ) then

      write ( iunit, '(a)' ) ' stroke'
      text_num = text_num + 1
      lineopen = .false.

    else

      imat = line_material(i)

      if ( imat /= imat_old ) then
        imat_old = imat
        red = material_rgba(1,imat)
        green = material_rgba(2,imat)
        blue = material_rgba(3,imat)
        write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
      end if

      px = plotxmin2 + nint ( alpha * ( cor2(1,j) - xmin ) )
      py = plotymin2 + nint ( alpha * ( cor2(2,j) - ymin ) )

      if ( lineopen ) then
        write ( iunit, '(2i4,a)' ) px, py, ' lineto'
        text_num = text_num + 1
      else
        write ( iunit, '(a)' ) ' newpath'
        write ( iunit, '(2i4,a)' ) px, py, ' moveto'
        text_num = text_num + 2
        lineopen = .true.
      end if

    end if

  end do
!
!  Write the epilog.
!
  write ( iunit, '(a)' ) 'showpage'
  write ( iunit, '(a)' ) 'grestore'
  write ( iunit, '(a)' ) '%%Trailer'
  write ( iunit, '(a)' ) 'end'
  write ( iunit, '(a)' ) '%%EOF'
  text_num = text_num + 5
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'PS_WRITE - Wrote ', text_num, ' text lines to ' &
    // trim ( fileout_name )

  return
end
subroutine r4_random ( rlo, rhi, r )

!*****************************************************************************80
!
!! R4_RANDOM returns a random real in a given range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) RLO, RHI, the minimum and maximum values.
!
!    Output, real ( kind = 4 ) R, the randomly chosen value.
!
  implicit none

  logical, save :: seed = .false.
  real ( kind = 4 ) r
  real ( kind = 4 ) rhi
  real ( kind = 4 ) rlo
  real ( kind = 4 ) t

  if ( .not. seed ) then
    call random_seed
    seed = .true.
  end if
!
!  Pick a random number in (0,1).
!
  call random_number ( harvest = t )
!
!  Set R.
!
  r = ( 1.0E+00 - t ) * rlo + t * rhi

  return
end
subroutine r4_swap ( x, y )

!*****************************************************************************80
!
!! R4_SWAP switches two R4's.
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
!    Input/output, real ( kind = 4 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 4 ) x
  real ( kind = 4 ) y
  real ( kind = 4 ) z

  z = x
  x = y
  y = z

  return
end
subroutine r4col_find ( lda, m, n, a, x, icol )

!*****************************************************************************80
!
!! R4COL_FIND seeks a table column equal to a real vector.
!
!  Example:
!
!    Input:
!
!      M = 3,
!      N = 4,
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!      x = ( 3.,
!            7.,
!           11. )
!
!    Output:
!
!      ICOL = 3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 4 ) A(LDA,N), a table of numbers, regarded as
!    N columns of vectors of length M.
!
!    Input, real ( kind = 4 ) X(M), a vector to be matched with a column of A.
!
!    Output, integer ( kind = 4 ) ICOL, the index of the first column of A
!    which exactly matches every entry of X, or 0 if no match
!    could be found.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) j
  real ( kind = 4 ) x(m)

  icol = 0

  do j = 1, n

    icol = j

    do i = 1, m
      if ( x(i) /= a(i,j) ) then
        icol = 0
        exit
      end if
    end do

    if ( icol /= 0 ) then
      return
    end if

  end do

  return
end
subroutine rgb_to_hue ( r, g, b, h )

!*****************************************************************************80
!
!! RGB_TO_HUE converts (R,G,B) colors to a hue value between 0 and 1.
!
!  Discussion:
!
!    The hue computed here should be the same as the H value computed
!    for HLS and HSV, except that it ranges from 0 to 1 instead of
!    0 to 360.
!
!    A monochromatic color ( white, black, or a shade of gray) does not
!    have a hue.  This routine will return a special value of H = -1
!    for such cases.
!
!  Example:
!
!    Color    R    G    B     H
!
!    red      1.0  0.0  0.0   0.00
!    yellow   1.0  1.0  0.0   0.16
!    green    0.0  1.0  0.0   0.33
!    cyan     0.0  1.0  1.0   0.50
!    blue     0.0  0.0  1.0   0.67
!    magenta  1.0  0.0  1.0   0.83
!
!    black    0.0  0.0  0.0  -1.00
!    gray     0.5  0.5  0.5  -1.00
!    white    1.0  1.0  1.0  -1.00
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) R, G, B, the red, green and blue values of the color.
!    These values should be between 0 and 1.
!
!    Output, real ( kind = 4 ) H, the corresponding hue of the color, or -1.0 if
!    the color is monochromatic.
!
  implicit none

  real ( kind = 4 ) b
  real ( kind = 4 ) b2
  real ( kind = 4 ) g
  real ( kind = 4 ) g2
  real ( kind = 4 ) h
  real ( kind = 4 ) r
  real ( kind = 4 ) r2
  real ( kind = 4 ) rgbmax
  real ( kind = 4 ) rgbmin
!
!  Make sure the colors are between 0 and 1.
!
  r2 = min ( max ( r, 0.0E+00 ), 1.0E+00 )
  g2 = min ( max ( g, 0.0E+00 ), 1.0E+00 )
  b2 = min ( max ( b, 0.0E+00 ), 1.0E+00 )
!
!  Compute the minimum and maximum of R, G and B.
!
  rgbmax = r2
  rgbmax = max ( rgbmax, g2 )
  rgbmax = max ( rgbmax, b2 )

  rgbmin = r2
  rgbmin = min ( rgbmin, g2 )
  rgbmin = min ( rgbmin, b2 )
!
!  If RGBMAX = RGBMIN, then the color has no hue.
!
  if ( rgbmax == rgbmin ) then

    h = - 1.0E+00
!
!  Otherwise, we need to determine the dominant color.
!
  else

    if ( r2 == rgbmax ) then
      h = ( g2 - b2 ) / ( rgbmax - rgbmin )
    else if ( g2 == rgbmax ) then
      h = 2.0E+00 + ( b2 - r2 ) / ( rgbmax - rgbmin )
    else if ( b2 == rgbmax ) then
      h = 4.0E+00 + ( r2 - g2 ) / ( rgbmax - rgbmin )
    end if

    h = h / 6.0E+00
!
!  Make sure H lies between 0 and 1.0.
!
    if ( h < 0.0E+00 ) then
      h = h + 1.0E+00
    else if ( 1.0E+00 < h ) then
      h = h - 1.0E+00
    end if

  end if

  return
end
subroutine relnex ( line, rval, done )

!*****************************************************************************80
!
!! RELNEX "reads" real numbers from a string, one at a time.
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
!    Input, character ( len = * ) LINE, a string, presumably containing real
!    numbers.  These may be separated by spaces or commas.
!
!    Output, real ( kind = 4 ) RVAL.  If DONE is FALSE, then RVAL contains the
!    "next" real value read from LINE.  If DONE is TRUE, then
!    RVAL is zero.
!
!    Input/output, logical DONE.
!    On input with a fresh value of LINE, the user should set
!    DONE to TRUE.
!    On output, the routine sets DONE to FALSE if another real
!    value was read, or TRUE if no more reals could be read.
!
  implicit none

  logical done
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) lchar
  character ( len = * ) line
  integer ( kind = 4 ), save :: next = 1
  real ( kind = 4 ) rval

  rval = 0.0E+00

  if ( done ) then
    next = 1
    done = .false.
  end if

  if ( len_trim ( line ) < next ) then
    done = .true.
    return
  end if

  call s_to_r4 ( line(next:), rval, ierror, lchar )

  if ( ierror /= 0 .or. lchar == 0 ) then
    done = .true.
    next = 1
  else
    done = .false.
    next = next + lchar
  end if

  return
end
subroutine r4vec_to_s ( n, x, s )

!*****************************************************************************80
!
!! R4VEC_TO_S "writes" an R4VEC into a string.
!
!  Discussion:
!
!    The values will be separated by commas and a single space.
!    If the string is too short, then data will be lost.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of X.
!
!    Input, real ( kind = 4 ) X(N), a vector to be written to a string.
!
!    Output, character ( len = * ) S, a string to which the real vector
!    has been written.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  character ( len = * ) s
  character ( len = 14 ) s2
  real ( kind = 4 ) x(n)

  do i = 1, n

    if ( x(i) == 0.0E+00 ) then
      s2 = '0'
    else if ( 1.0E+10 <= abs ( x(i) ) ) then
      write ( s2, '(g14.6)' ) x(i)
      call s_trim_zeros ( s2 )
    else if ( real ( int ( x(i) ) ) == x(i) ) then
      write ( s2, '(i12)' ) int ( x(i) )
    else
      write ( s2, '(g14.6)' ) x(i)
      call s_trim_zeros ( s2 )
    end if

    if ( i == 1 ) then
      s = adjustl ( s2 )
    else
      s = trim ( s ) // ', ' // adjustl ( s2 )
    end if

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
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) nchar
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  iput = 0
  nchar = len_trim ( s )

  do iget = 1, nchar

    c = s(iget:iget)

    if ( c /= ' ' .and. c /= TAB ) then
      iput = iput + 1
      s(iput:iput) = c
    end if

  end do

  s(iput+1:) = ' '

  return
end
subroutine s_blanks_delete ( s )

!*****************************************************************************80
!
!! S_BLANKS_DELETE replaces consecutive blanks by one blank.
!
!  Discussion:
!
!    The remaining characters are left justified and right padded with blanks.
!    TAB characters are converted to spaces.
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

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nchar
  character newchr
  character oldchr
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  j = 0
  newchr = ' '
  nchar = len_trim ( s )

  do i = 1, nchar

    oldchr = newchr
    newchr = s(i:i)

    if ( newchr == TAB ) then
      newchr = ' '
    end if

    s(i:i) = ' '

    if ( oldchr /= ' ' .or. newchr /= ' ' ) then
      j = j + 1
      s(j:j) = newchr
    end if

  end do

  return
end
subroutine s_cap ( string )

!*****************************************************************************80
!
!! S_CAP replaces any lowercase letters by uppercase ones in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) STRING, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nchar
  character ( len = * ) string

  nchar = len_trim ( string )

  do i = 1, nchar

    c = string(i:i)
    call ch_cap ( c )
    string(i:i) = c

  end do

  return
end
subroutine s_cat ( s1, s2, s3 )

!*****************************************************************************80
!
!! S_CAT concatenates two strings to make a third string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the "prefix" string.
!
!    Input, character ( len = * ) S2, the "postfix" string.
!
!    Output, character ( len = * ) S3, the string made by
!    concatenating S1 and S2, ignoring any trailing blanks.
!
  implicit none

  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3

  s3 = trim ( s1 ) // trim ( s2 )

  return
end
subroutine s_control_blank ( s )

!*****************************************************************************80
!
!! S_CONTROL_BLANK replaces control characters with blanks.
!
!  Discussion:
!
!    A "control character" has ASCII code <= 31 or 127 <= ASCII code.
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
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  logical ch_is_control
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nchar
  character ( len = * ) s

  nchar = len_trim ( s )

  do i = 1, nchar
    if ( ch_is_control ( s(i:i) ) ) then
      s(i:i) = ' '
    end if
  end do

  return
end
function s_eqi ( strng1, strng2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
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
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( s1 /= s2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

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
function s_is_i4 ( string, ival )

!*****************************************************************************80
!
!! S_IS_I4 returns .TRUE. if STRING represents an integer.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the string to be checked.
!
!    Output, integer ( kind = 4 ) IVAL.  If S_IS_INT is TRUE, then IVAL is the
!    integer ( kind = 4 ) represented.  Otherwise IVAL is 0.
!
!    Output, logical S_IS_I4, .TRUE. if STRING represents an integer.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) lenc
  logical s_is_i4
  character ( len = * ) string

  lenc = len_trim ( string )

  call s_to_i4 ( string, ival, ierror, lchar )

  if ( ierror == 0 .and. lenc <= lchar ) then
    s_is_i4 = .true.
  else
    s_is_i4 = .false.
    ival = 0
  end if

  return
end
subroutine s_is_r4 ( string, rval, lval )

!*****************************************************************************80
!
!! S_IS_R4 returns .TRUE. if STRING represents a real number.
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
!    Input, character ( len = * ) STRING, the string to be checked.
!
!    Output, real ( kind = 4 ) RVAL.  If ISREAL is TRUE, then RVAL is the real
!    number represented.  Otherwise RVAL is 0.
!
!    Output, logical LVAL, .TRUE. if STRING represents a real number.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) lchar
  logical lval
  real ( kind = 4 ) rval
  character ( len = * ) string

  call s_to_r4 ( string, rval, ierror, lchar )

  if ( ierror == 0 .and. len_trim ( string ) <= lchar ) then
    lval = .true.
  else
    lval = .false.
    rval = 0.0E+00
  end if

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If STRING is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character of S used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        last = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    last = len_trim ( s )
  else
    ierror = 1
    last = 0
  end if

  return
end
subroutine s_to_r4 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R4 reads an R4 from a string.
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
!    '1'               1.0E+00
!    '     1   '       1.0E+00
!    '1A'              1.0E+00
!    '12,34,56'        12.0E+00
!    '  34 7'          34.0E+00
!    '-1E2ABCD'        -100.0E+00
!    '-1X2ABCD'        -1.0E+00
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0E+00
!    '17d2'            1700.0E+00
!    '-14e-2'         -0.14
!    'e2'              100.0E+00
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
!    Output, real ( kind = 4 ) R, the real value that was read from the string.
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
  real ( kind = 4 ) r
  real ( kind = 4 ) rbot
  real ( kind = 4 ) rexp
  real ( kind = 4 ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0E+00
  lchar = - 1
  isgn = 1
  rtop = 0.0E+00
  rbot = 1.0E+00
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
        rtop = 10.0E+00 * rtop + real ( ndig )
      else if ( ihave == 5 ) then
        rtop = 10.0E+00 * rtop + real ( ndig )
        rbot = 10.0E+00 * rbot
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
    rexp = 1.0E+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0E+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0E+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

  return
end
subroutine s_trim_zeros ( s )

!*****************************************************************************80
!
!! S_TRIM_ZEROS removes trailing zeros from a string.
!
!  Example:
!
!    Input:
!
!      S = '1401.072500'
!
!    Output:
!
!      S = '1401.0725'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be operated on.
!
  implicit none

  integer ( kind = 4 ) i
  character ( len = * ) s

  i = len_trim ( s )

  do while ( 0 < i .and. s(i:i) == '0' )
    s(i:i) = ' '
    i = i - 1
  end do

  return
end
subroutine smf_read ( bad_num, cor3, cor3_material, cor3_max, cor3_normal, &
  cor3_num, cor3_tex_uv, debug, face, face_material, face_max, &
  face_normal, face_num, face_order, face_tex_uv, filein_name, &
  group_num, ierror, iunit, material_max, material_name, material_num, &
  material_rgba, order_max, text_num, texture_max, texture_num, &
  texture_name, vertex_material )

!*****************************************************************************80
!
!! SMF_READ reads graphics information from an SMF file.
!
!  Discussion:
!
!    It is intended that the information read from the file can
!    either start a whole new graphics object, or simply be added
!    to a current graphics object via the '<<' command.  
!
!    This is controlled by whether the input values have been zeroed 
!    out or not.  This routine simply tacks on the information it 
!    finds to the current graphics object.
!
!  Example:
!
!    #SMF2.0
!    #  cube_face.smf
!    #  This example demonstrates how an RGB color can be assigned to
!    #  each face of an object.
!    #    
!    # First, define the geometry of the cube.
!    #
!    v 0.0  0.0  0.0
!    v 1.0  0.0  0.0
!    v 0.0  1.0  0.0
!    v 1.0  1.0  0.0
!    v 0.0  0.0  1.0
!    v 1.0  0.0  1.0
!    v 0.0  1.0  1.0
!    v 1.0  1.0  1.0
!    f 1 4 2
!    f 1 3 4
!    f 5 6 8
!    f 5 8 7
!    f 1 2 6
!    f 1 6 5
!    f 2 4 8
!    f 2 8 6
!    f 4 3 7
!    f 4 7 8
!    f 3 1 5
!    f 3 5 7
!    #
!    #  Colors will be bound 1 per face.
!    #
!    bind c face
!    c 1.0  0.0  0.0
!    c 1.0  0.0  0.0
!    c 0.0  1.0  0.0
!    c 0.0  1.0  0.0
!    c 0.0  0.0  1.0
!    c 0.0  0.0  1.0
!    c 1.0  1.0  0.0
!    c 1.0  1.0  0.0
!    c 0.0  1.0  1.0
!    c 0.0  1.0  1.0
!    c 1.0  0.0  1.0
!    c 1.0  0.0  1.0
!    #
!    #  Normal vectors will be bound 1 per face.
!    #
!    bind n face
!    n  0.0   0.0  -1.0
!    n  0.0   0.0  -1.0
!    n  0.0   0.0   1.0
!    n  0.0   0.0   1.0
!    n  0.0  -1.0   0.0
!    n  0.0  -1.0   0.0
!    n  1.0   0.0   0.0
!    n  1.0   0.0   0.0
!    n  0.0   1.0   0.0
!    n  0.0   1.0   0.0
!    n -1.0   0.0   0.0
!    n -1.0   0.0   0.0
!    #
!    #  Texture coordinate pairs will be bound 1 per face.
!    #
!    bind r face
!    r  0.0   0.0
!    r  0.0   0.1
!    r  0.0   0.2
!    r  0.0   0.3
!    r  0.1   0.0
!    r  0.1   0.1
!    r  0.1   0.2
!    r  0.1   0.3
!    r  0.2   0.0
!    r  0.2   0.1
!    r  0.2   0.2
!    r  0.2   0.3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input/output, integer ( kind = 4 ) COR3_MATERIAL(COR3_MAX), the material index of 
!    each node.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input/output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input/output, real ( kind = 4 ) COR3_TEX_UV(2,COR3_MAX), UV texture coordinates 
!    for nodes.
!
!    Input, logical DEBUG, debugging switch.
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input/output, real ( kind = 4 ) FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of 
!    vertices per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which data is read.
!
!    Input/output, character ( len = * ) MATERIAL_NAME(MATERIAL_MAX), 
!    material names.
!
!    Input/output, real ( kind = 4 ) MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and 
!    A values.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Output, integer ( kind = 4 ) BAD_NUM, the number of bad lines of text in the file.
!
!    Input/output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input/output, integer ( kind = 4 ) GROUP_NUM, the number of groups.
!
!    Input/output, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Output, integer ( kind = 4 ) TEXT_NUM, the number of lines of text in the file.
!
!    Input/output, character ( len = * ) TEXTURE_NAME(TEXTURE_MAX), 
!    texture names.
!
!    Output, real ( kind = 4 ) TRANSFORM_MATRIX(4,4), the transformation matrix.
!
!    Input/output, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) texture_max

  real ( kind = 4 ) angle
  character axis
  real ( kind = 4 ) b
  integer ( kind = 4 ) bad_num
  character ( len = 4 ) char4
  character cnr
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_material(cor3_max)
  real ( kind = 4 ) cor3_normal(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  real ( kind = 4 ) cor3_tex_uv(2,cor3_max)
  logical debug
  logical done
  real ( kind = 4 ) dx
  real ( kind = 4 ) dy
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_count
  integer ( kind = 4 ) face_material(face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  real ( kind = 4 ) face_tex_uv(2,face_max)
  character ( len = * ) filein_name
  real ( kind = 4 ) g
  integer ( kind = 4 ) group_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icor3_normal
  integer ( kind = 4 ) icor3_tex_uv
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface_normal
  integer ( kind = 4 ) iface_tex_uv
  integer ( kind = 4 ) imat
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) iword
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) level
  character ( len = 256 ) line
  character ( len = 30 ) material_binding
  character ( len = * ) material_name(material_max)
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ) node_count
  character ( len = 30 ) normal_binding
  real ( kind = 4 ) r
  real ( kind = 4 ) rgba(4)
  logical s_eqi
  real ( kind = 4 ) sx
  real ( kind = 4 ) sy
  real ( kind = 4 ) sz
  real ( kind = 4 ) temp
  integer ( kind = 4 ) text_num
  character ( len = 30 ) texture_binding
  character ( len = * ) texture_name(texture_max)
  integer ( kind = 4 ) texture_num
  real ( kind = 4 ) transform_matrix(4,4)
  character ( len = 30 ) type
  real ( kind = 4 ) u
  real ( kind = 4 ) v
  integer ( kind = 4 ) vertex_base
  integer ( kind = 4 ) vertex_correction
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  character ( len = 256 ) word
  character ( len = 256 ) word1
  real ( kind = 4 ) x
  real ( kind = 4 ) xvec(3)
  real ( kind = 4 ) y
  real ( kind = 4 ) z

  face_count = 0
  ierror = 0
  icor3_normal = 0
  icor3_tex_uv = 0
  iface_normal = 0
  iface_tex_uv = 0
  level = 0
  material_binding = 'UNDEFINED'
  normal_binding = 'UNDEFINED'
  node_count = 0
  texture_binding = 'UNDEFINED'
  vertex_base = cor3_num
  vertex_correction = 0
  word = ' '

  call tmat_init ( transform_matrix )
!
!  Read a line of text from the file.
!
10    continue
 
  read ( iunit, '(a)', iostat = ios ) line

  if ( ios /= 0 ) then
    go to 30
  end if

  if ( debug ) then
    write ( *, '(a)' ) trim ( line )
  end if

  text_num = text_num + 1
!
!  Replace any control characters (in particular, TAB's) by blanks.
!
  call s_control_blank ( line )
 
  done = .true.
  iword = 0
!
!  Read a word from the line.
!
  call word_next_read ( line, word, done )

  if ( debug ) then
    write ( *, '(a)' ) trim ( word )
  end if
!
!  If no more words in this line, read a new line.
!
  if ( done ) then
    go to 10
  end if
!
!  If this word begins with '#' or '$', then it's a comment.  Read a new line.
!
  if ( word(1:1) == '#' .or. word(1:1) == '$' ) then
    go to 10
  end if

  iword = iword + 1
  if ( iword == 1 ) then
    word1 = word
  end if
!
!  BEGIN
!  Reset the transformation matrix to identity.
!  Node numbering starts at zero again.  (Really, this is level based)
!  (Really should define a new transformation matrix, and concatenate.)
!  (Also, might need to keep track of level.)
!
  if ( s_eqi ( word1, 'BEGIN' ) ) then

    level = level + 1

    vertex_base = cor3_num
    group_num = group_num + 1
    call tmat_init ( transform_matrix )
!
!  BIND [c|n|r] [vertex|face]
!  Specify the binding for RGB color, Normal, or Texture.
!  Options are "vertex" or "face"
!
  else if ( s_eqi ( word1, 'BIND' ) ) then

    call word_next_read ( line, cnr, done )

    call word_next_read ( line, type, done )

    if ( s_eqi ( cnr, 'C' ) ) then

      if ( s_eqi ( type, 'VERTEX' ) ) then
        material_binding = 'PER_VERTEX'
      else if ( s_eqi ( type, 'FACE' ) ) then
        material_binding = 'PER_FACE'
      end if

    else if ( s_eqi ( cnr, 'N' ) ) then

      if ( s_eqi ( type, 'VERTEX' ) ) then
        normal_binding = 'PER_VERTEX'
      else if ( s_eqi ( type, 'FACE' ) ) then
        normal_binding = 'PER_FACE'
      end if

    else if ( s_eqi ( cnr, 'R' ) ) then

      if ( s_eqi ( type, 'VERTEX' ) ) then
        texture_binding = 'PER_VERTEX'
      else if ( s_eqi ( type, 'FACE' ) ) then
        texture_binding = 'PER_FACE'
      end if

    end if
!
!  C <r> <g> <b>
!  Specify an RGB color, with R, G, B between 0.0 and 1.0.
!
  else if ( s_eqi ( word1, 'C' ) ) then

    call word_next_read ( line, word, done )
    call s_to_r4 ( word, r, ierror, lchar )
    call word_next_read ( line, word, done )
    call s_to_r4 ( word, g, ierror, lchar )
    call word_next_read ( line, word, done )
    call s_to_r4 ( word, b, ierror, lchar )
!
!    Set up a temporary material (R,G,B,1.0).
!    Add the material to the material database, or find the index of
!      a matching material already in.
!    Assign the material of the node or face to this index.
!
    rgba(1) = r
    rgba(2) = g
    rgba(3) = b
    rgba(4) = 1.0E+00

    if ( material_num <= 1000 ) then
      call r4col_find ( 4, 4, material_num, material_rgba, rgba, imat )
    else
      imat = 0
    end if

    if ( imat == 0 ) then

      material_num = material_num + 1

      if ( material_num <= material_max ) then

        call i4_to_s_zero ( material_num, char4 )

        material_name(material_num) = 'Material_' // char4
        material_rgba(1:4,material_num) = rgba(1:4)
        imat = material_num

      else

        imat = 0

      end if

    end if

    if ( material_binding == 'PER_FACE' ) then

      face_count = face_count + 1
      face_material(face_count) = imat

    else if ( material_binding == 'PER_VERTEX' ) then

      node_count = node_count + 1
      cor3_material(node_count) = imat

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SMF_READ - Fatal error!'
      write ( *, '(a)' ) '  Material binding undefined!'
      stop

    end if
!
!  END
!  Drop down a level.
!
  else if ( s_eqi ( word1, 'END' ) ) then

    level = level - 1

    if ( level < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SMF_READ - Fatal error!'
      write ( *, '(a)' ) '  More END statements than BEGINs.'
      write ( *, '(a,i6)' ) '  Stopping on line ', text_num
      stop
    end if
!
!  F V1 V2 V3 ...
!  Face.
!  A face is defined by the vertices.
!
  else if ( s_eqi ( word1, 'F' ) ) then

    face_num = face_num + 1
    face_material(face_num) = material_num

    ivert = 0

    do

      ivert = ivert + 1
 
      call word_next_read ( line, word, done )

      if ( done ) then
        exit
      end if
!
!  Read the vertex index.
!  Note that vertex indices start back at 0 each time a BEGIN is entered.
!  The strategy here won't handle nested BEGIN's, just one at a time.
!
      call s_to_i4 ( word, itemp, ierror, lchar )

      if ( ierror /= 0 ) then
        itemp = -1
        ierror = 0
        write ( *, '(a)' ) 'SMF_READ - Error!'
        write ( *, '(a)' ) '  Bad FACE field.'
        write ( *, '(a)' ) trim ( word )
      end if

      if ( ivert <= order_max .and. face_num <= face_max ) then
        face(ivert,face_num) = itemp + vertex_base
        vertex_material(ivert,face_num) = material_num
      end if

      if ( face_num <= face_max ) then
        face_order(face_num) = ivert
      end if

    end do

    go to 10
!
!  N <x> <y> <z>
!  Specify a normal vector.
!
  else if ( s_eqi ( word1, 'N' ) ) then

    call word_next_read ( line, word, done )
    call s_to_r4 ( word, x, ierror, lchar )
    call word_next_read ( line, word, done )
    call s_to_r4 ( word, y, ierror, lchar )
    call word_next_read ( line, word, done )
    call s_to_r4 ( word, z, ierror, lchar )

    if ( normal_binding == 'PER_FACE' ) then

      iface_normal = iface_normal + 1

      face_normal(1,iface_normal) = x
      face_normal(2,iface_normal) = y
      face_normal(3,iface_normal) = z

    else if ( normal_binding == 'PER_VERTEX' ) then

      icor3_normal = icor3_normal + 1

      cor3_normal(1,icor3_normal) = x
      cor3_normal(2,icor3_normal) = y
      cor3_normal(3,icor3_normal) = z

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SMF_READ - Fatal error!'
      write ( *, '(a)' ) '  Normal binding undefined!'
      stop
                
    end if
!
!  R <u> <v>
!  Specify a texture coordinate.
!  (Read, but ignore for now.)
!
  else if ( s_eqi ( word1, 'R' ) ) then

    call word_next_read ( line, word, done )
    call s_to_r4 ( word, u, ierror, lchar )
    call word_next_read ( line, word, done )
    call s_to_r4 ( word, v, ierror, lchar )

    if ( texture_binding == 'PER_FACE' ) then

      iface_tex_uv = iface_tex_uv + 1

      face_tex_uv(1,iface_tex_uv) = u
      face_tex_uv(2,iface_tex_uv) = v

    else if ( texture_binding == 'PER_VERTEX' ) then

      icor3_tex_uv = icor3_tex_uv + 1

      cor3_tex_uv(1,icor3_tex_uv) = u
      cor3_tex_uv(2,icor3_tex_uv) = v

    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SMF_READ - Fatal error!'
      write ( *, '(a)' ) '  Texture binding undefined!'
      stop
    end if
!
!  ROT [x|y|z] <theta>
!  (Read, but ignore for now.)
!
  else if ( s_eqi ( word1, 'ROT' ) ) then

    call word_next_read ( line, axis, done )

    call word_next_read ( line, word, done )
    call s_to_r4 ( word, angle, ierror, lchar )

    call tmat_rot_axis ( transform_matrix, transform_matrix, angle, axis )
!
!  SCALE <sx> <sy> <sz>
!  (Read, but ignore for now.)
!
  else if ( s_eqi ( word1, 'SCALE' ) ) then

    call word_next_read ( line, word, done )
    call s_to_r4 ( word, sx, ierror, lchar )
    call word_next_read ( line, word, done )
    call s_to_r4 ( word, sy, ierror, lchar )
    call word_next_read ( line, word, done )
    call s_to_r4 ( word, sz, ierror, lchar )

    call tmat_scale ( transform_matrix, transform_matrix, sx, sy, sz )
!
!  SET VERTEX_CORRECTION <i>
!  Specify increment to add to vertex indices in file.
!  (Read, but ignore for now.)
!
  else if ( s_eqi ( word1, 'SET' ) ) then

    call word_next_read ( line, word, done )

    call word_next_read ( line, word, done )
    call s_to_i4 ( word, vertex_correction, ierror, lchar )
!
!  T V1 V2 V3
!  Triangle
!  (Added 30 October 2001, JVB)
!
  else if ( s_eqi ( word1, 'T' ) ) then

    face_num = face_num + 1
    face_material(face_num) = material_num

    do ivert = 1, 3
 
      call word_next_read ( line, word, done )

      if ( done ) then
        exit
      end if
!
!  Read the vertex index.
!  Note that vertex indices start back at 0 each time a BEGIN is entered.
!  The strategy here won't handle nested BEGIN's, just one at a time.
!
      call s_to_i4 ( word, itemp, ierror, lchar )

      if ( ierror /= 0 ) then
        itemp = -1
        ierror = 0
        write ( *, '(a)' ) 'SMF_READ - Error!'
        write ( *, '(a)' ) '  Bad TRIANGLE field.'
        write ( *, '(a)' ) trim ( word )
      end if

      if ( ivert <= order_max .and. face_num <= face_max ) then
        face(ivert,face_num) = itemp + vertex_base
        vertex_material(ivert,face_num) = material_num
      end if

      if ( face_num <= face_max ) then
        face_order(face_num) = ivert
      end if

    end do

    go to 10
!
!  T_SCALE <dx> <dy>
!  Specify a translation to texture coordinates.
!  (Read, but ignore for now.)
!
  else if ( s_eqi ( word1, 'T_SCALE' ) ) then

    call word_next_read ( line, word, done )
    call s_to_r4 ( word, dx, ierror, lchar )
    call word_next_read ( line, word, done )
    call s_to_r4 ( word, dy, ierror, lchar )
!
!  T_TRANS <dx> <dy>
!  Specify a translation to texture coordinates.
!  (Read, but ignore for now.)
!
  else if ( s_eqi ( word1, 'T_TRANS' ) ) then

    call word_next_read ( line, word, done )
    call s_to_r4 ( word, dx, ierror, lchar )
    call word_next_read ( line, word, done )
    call s_to_r4 ( word, dy, ierror, lchar )
!
!  TEX <filename>
!  Specify a filename containing the texture.
!
  else if ( s_eqi ( word1, 'TEX' ) ) then

    call word_next_read ( line, word, done )

    texture_num = texture_num + 1
    texture_name(texture_num) = word
!
!  TRANS <dx> <dy> <dz>
!
  else if ( s_eqi ( word1, 'TRANS' ) ) then

    call word_next_read ( line, word, done )
    call s_to_r4 ( word, x, ierror, lchar )
    call word_next_read ( line, word, done )
    call s_to_r4 ( word, y, ierror, lchar )
    call word_next_read ( line, word, done )
    call s_to_r4 ( word, z, ierror, lchar )

    call tmat_trans ( transform_matrix, transform_matrix, x, y, z )
!
!  V X Y Z
!  Geometric vertex.
!
  else if ( s_eqi ( word1, 'V' ) ) then

    cor3_num = cor3_num + 1

    do i = 1, 3
      call word_next_read ( line, word, done )
      call s_to_r4 ( word, temp, ierror, lchar )
      xvec(i) = temp
    end do
!
!  Apply current transformation matrix.
!  Right now, we can only handle one matrix, not a stack of
!  matrices representing nested BEGIN/END's.
!
    call tmat_mxp ( transform_matrix, xvec, xvec )

    if ( cor3_num <= cor3_max ) then
      cor3(1:3,cor3_num) = xvec(1:3)
    end if
!
!  Unrecognized keyword.
!
  else

    bad_num = bad_num + 1

    if ( bad_num <= 10 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'SMF_READ: Bad data on line ', text_num
      write ( *, '(a)' ) '  Bad word: ' // trim ( word )
    end if

  end if

  go to 10
!
!  End of information in file.
!
30    continue
!
!  Extend the material definition 
!  * from the face to the vertices and nodes, or
!  * from the vertices to the faces and nodes.
!
  if ( material_binding == 'PER_FACE' ) then

    do ivert = 1, order_max
      vertex_material(ivert,1:face_num) = face_material(1:face_num)
    end do

    call vertex_to_node_material ( cor3_material, cor3_max, face, &
      face_order, face_max, order_max, face_num, vertex_material )

  else if ( material_binding == 'PER_VERTEX' ) then

    call node_to_vertex_material ( cor3_material, cor3_max, face, &
      face_max, face_num, face_order, order_max, vertex_material )

    face_material(1:face_num) = vertex_material(1,1:face_num)

  end if
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'SMF_READ - Read ', text_num, ' text lines from ' &
    // trim ( filein_name )

  return
end
subroutine smf_write ( cor3, cor3_material, cor3_max, cor3_normal, cor3_num, &
  cor3_tex_uv, face, face_max, face_num, face_order, filein_name, &
  fileout_name, iunit, material_max, material_rgba, order_max, &
  texture_max, texture_name, texture_num )

!*****************************************************************************80
!
!! SMF_WRITE writes graphics information to an SMF file.
!
!  Example:
!
!    #SMF2.0E+00
!    #  cube_face.smf
!    #  This example demonstrates how an RGB color can be assigned to
!    #  each face of an object.
!    #    
!    # First, define the geometry of the cube.
!    #
!    v 0.0  0.0  0.0E+00
!    v 1.0  0.0  0.0E+00
!    v 0.0  1.0  0.0E+00
!    v 1.0  1.0  0.0E+00
!    v 0.0  0.0  1.0E+00
!    v 1.0  0.0  1.0E+00
!    v 0.0  1.0  1.0E+00
!    v 1.0  1.0  1.0E+00
!    f 1 4 2
!    f 1 3 4
!    f 5 6 8
!    f 5 8 7
!    f 1 2 6
!    f 1 6 5
!    f 2 4 8
!    f 2 8 6
!    f 4 3 7
!    f 4 7 8
!    f 3 1 5
!    f 3 5 7
!    #
!    #  Colors will be bound 1 per face.
!    #
!    bind c face
!    c 1.0  0.0  0.0
!    c 1.0  0.0  0.0
!    c 0.0  1.0  0.0
!    c 0.0  1.0  0.0
!    c 0.0  0.0  1.0
!    c 0.0  0.0  1.0
!    c 1.0  1.0  0.0
!    c 1.0  1.0  0.0
!    c 0.0  1.0  1.0
!    c 0.0  1.0  1.0
!    c 1.0  0.0  1.0
!    c 1.0  0.0  1.0
!    #
!    #  Normal vectors will be bound 1 per face.
!    #
!    bind n face
!    n  0.0   0.0  -1.0
!    n  0.0   0.0  -1.0
!    n  0.0   0.0   1.0
!    n  0.0   0.0   1.0
!    n  0.0  -1.0   0.0
!    n  0.0  -1.0   0.0
!    n  1.0   0.0   0.0
!    n  1.0   0.0   0.0
!    n  0.0   1.0   0.0
!    n  0.0   1.0   0.0
!    n -1.0   0.0   0.0
!    n -1.0   0.0   0.0
!    #
!    #  Texture coordinate pairs will be bound 1 per face.
!    #
!    bind r face
!    r  0.0   0.0
!    r  0.0   0.1
!    r  0.0   0.2
!    r  0.0   0.3
!    r  0.1   0.0
!    r  0.1   0.1
!    r  0.1   0.2
!    r  0.1   0.3
!    r  0.2   0.0
!    r  0.2   0.1
!    r  0.2   0.2
!    r  0.2   0.3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MATERIAL(COR3_MAX), the material index of each node.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, real ( kind = 4 ) COR3_NORMAL(3,COR3_MAX), normals at nodes.  
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, real ( kind = 4 ) COR3_TEX_UV(2,COR3_MAX), UV texture coordinates for nodes.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, real ( kind = 4 ) MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and A values.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer ( kind = 4 ) TEXTURE_MAX, the maximum number of textures.
!
!    Input, character ( len = * ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
!    Input, integer ( kind = 4 ) TEXTURE_NUM, the number of textures.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) texture_max

  real ( kind = 4 ) b
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_material(cor3_max)
  real ( kind = 4 ) cor3_normal(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  real ( kind = 4 ) cor3_tex_uv(2,cor3_max)
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  character ( len = * ) fileout_name
  real ( kind = 4 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) imat
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ivert
  real ( kind = 4 ) material_rgba(4,material_max)
  real ( kind = 4 ) r
  character ( len = 256 ) text
  integer ( kind = 4 ) text_num
  character ( len = * ) texture_name(texture_max)
  integer ( kind = 4 ) texture_num

  text_num = 0

  write ( iunit, '(a)' ) '#$SMF 1.0'
  write ( iunit, '(a, i8)' ) '#$vertices ', cor3_num
  write ( iunit, '(a, i8)' ) '#$faces ', face_num
  write ( iunit, '(a)' ) '#'
  write ( iunit, '(a)' ) '# ' // trim ( fileout_name ) // ' created by IVREAD.'
  write ( iunit, '(a)' ) '# Original data in ' // trim ( filein_name ) // '.'
  write ( iunit, '(a)' ) '#'
  text_num = text_num + 7
!
!  V: vertex coordinates.
!
  do icor3 = 1, cor3_num
    write ( text, '(a1,2x,3g14.6)' ) 'v', cor3(1:3,icor3)
    call s_blanks_delete ( text )
    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1
  end do
!
!  F: Faces.
!
  if ( 0 < face_num ) then
    write ( iunit, '(a)' ) ' '
    text_num = text_num + 1
  end if

  do iface = 1, face_num
    write ( text, '(a1,2x,10i8)' ) 'f', &
      ( face(ivert,iface), ivert = 1, face_order(iface) )
    call s_blanks_delete ( text )
    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1
  end do
!
!  Material binding.
!
  write ( iunit, '(a)' ) 'bind c vertex'
  text_num = text_num + 1
!
!  Material RGB values at each node.
!
  do icor3 = 1, cor3_num
    imat = cor3_material(icor3)
    r = material_rgba(1,imat)
    g = material_rgba(2,imat)
    b = material_rgba(3,imat)
    write ( iunit, '(a,1x,3f6.2)' ) 'c', r, g, b
    text_num = text_num + 1
  end do
!
!  Normal binding.
!
  write ( iunit, '(a)' ) 'bind n vertex'
  text_num = text_num + 1
!
!  Normal vector at each node.
!
  do icor3 = 1, cor3_num
    write ( iunit, '(a,1x,3f6.2)' ) 'n', cor3_normal(1:3,icor3)
    text_num = text_num + 1
  end do

  if ( 0 < texture_num ) then
!
!  Texture filename
!
    write ( iunit, '(a)' ) 'tex ' // trim ( texture_name(1) )
    text_num = text_num + 1
!
!  Texture binding.
!
    write ( iunit, '(a)' ) 'bind r vertex'
    text_num = text_num + 1
!
!  Texture coordinates at each node.
!
    do icor3 = 1, cor3_num
      write ( iunit, '(a,1x,3f6.2)' ) 'r', cor3_tex_uv(1:2,icor3)
      text_num = text_num + 1
    end do

  end if
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'SMF_WRITE - Wrote ', text_num, &
    ' text lines to ' // trim ( fileout_name )

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into linear order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 November 2000
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert WIlf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I precedes J, ISGN = +1 if J precedes I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I precedes J;
!    0 <= ISGN means J precedes I.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    n1 = n
    k = n / 2
    k1 = k
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i = i + 1
      end if

      j = k1
      k1 = i
      indx = - 1
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        indx = 0
      else
        i = n1
        n1 = n1 - 1
        j = 1
        indx = 1
      end if

      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i = 2 * k1

    if ( i == n1 ) then
      j = k1
      k1 = i
      indx = - 1
      return
    else if ( i <= n1 ) then
      j = i + 1
      indx = - 2
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    indx = 0
  else
    i = n1
    n1 = n1 - 1
    j = 1
    indx = 1
  end if

  return
end
subroutine stla_read ( bad_num, cor3, cor3_max, cor3_num, dup_num, face, &
  face_material, face_max, face_normal, face_num, face_order, &
  filein_name, ierror, iunit, material_num, object_num, order_max, &
  text_num, vertex_material, vertex_normal )

!*****************************************************************************80
!
!! STLA_READ reads graphics information from an ASCII StereoLithography file.
!
!  Discussion:
!
!    It is intended that the information read from the file can
!    either start a whole new graphics object, or simply be added
!    to a current graphics object via the '<<' command.  
!
!    This is controlled by whether the input values have been zeroed 
!    out or not.  This routine simply tacks on the information it 
!    finds to the current graphics object.
!
!  Example:
!
!    solid MYSOLID
!      facet normal 0.4 0.4 0.2
!        outerloop
!          vertex  1.0 2.1 3.2
!          vertex  2.1 3.7 4.5
!          vertex  3.1 4.5 6.7
!        endloop
!      endfacet
!      ...
!      facet normal 0.2 0.2 0.4
!        outerloop
!          vertex  2.0 2.3 3.4
!          vertex  3.1 3.2 6.5
!          vertex  4.1 5.5 9.0
!        endloop
!      endfacet
!    endsolid MYSOLID
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    3D Systems, Inc,
!    Stereolithography Interface Specification,
!    October 1989.
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input/output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input/output, real ( kind = 4 ) FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input/output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices 
!    per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which data is read.
!
!    Input/output, integer ( kind = 4 ) LINE_NUM, the number of line items.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input/output, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input/output, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals 
!    at vertices.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  integer ( kind = 4 ) bad_num
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  logical done
  integer ( kind = 4 ) dup_num
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) material_num
  integer ( kind = 4 ) object_num
  real ( kind = 4 ) rval
  logical s_eqi
  real ( kind = 4 ) temp(3)
  character ( len = 256 ) text
  integer ( kind = 4 ) text_num
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  character ( len = 256 ) word1
  character ( len = 256 ) word2

  ierror = 0
  istate = 0
!
!  Read the next line of text.
!
  do
 
    read ( iunit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      if ( istate /= 0 .and. istate /= 1 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning.'
        write ( *, '(a)' ) '  End-of-file, but model not finished.'
      end if
      exit
    end if

    text_num = text_num + 1
    done = .true.
!
!  Read the first word in the line.
!
    call word_next_read ( text, word1, done )
!
!  "Doctor" the text, changing a beginning occurrence of:
!
!      END FACET to ENDFACET
!      END LOOP to ENDLOOP
!      END SOLID to ENDSOLID
!      FACET NORMAL to FACETNORMAL
!      OUTER LOOP to OUTERLOOP
!
    if ( s_eqi ( word1, 'END' ) .or. &
         s_eqi ( word1, 'FACET' ) .or. &
         s_eqi ( word1, 'OUTER' ) ) then

      call word_next_read ( text, word2, done )
      call s_cat ( word1, word2, word1 )

    end if
!
!  This first word tells us what to do.
!
!  SOLID - begin a new solid.
!    Valid in state 0, moves to state 1.
!  ENDSOLID - end current solid.
!    Valid in state 1, moves to state 0.
!
!  FACETNORMAL - begin a new facet.
!    Valid in state 0 or 1, moves to state 2.
!  ENDFACET - end current facet.
!    Valid in state 2, moves to state 1.
!
!  OUTERLOOP - begin a list of vertices.
!    Valid in state 2, moves to state 3.
!  ENDLOOP - end vertex list.
!    Valid in state 3, moves to state 2.
!
!  VERTEX - give coordinates of next vertex.
!    Valid in state 3.
!
!  End of file - 
!    Valid in state 0 or 1.
!
    if ( s_eqi ( word1, 'SOLID' ) ) then

      if ( istate == 0 ) then
        istate = 1
        object_num = object_num + 1
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning!'
        write ( *, '(a)' ) '  Model not in right state for SOLID.'
        bad_num = bad_num + 1
        ierror = 1
        exit
      end if

    else if ( s_eqi ( word1, 'ENDSOLID' ) ) then

      if ( istate == 1 ) then
        istate = 0
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning!'
        write ( *, '(a)' ) '  Model not in right state for ENDSOLID.'
        bad_num = bad_num + 1
        ierror = 1
        exit
      end if

    else if ( s_eqi ( word1, 'FACETNORMAL' ) ) then

      if ( istate == 0 .or. istate == 1 ) then

        istate = 2
        face_num = face_num + 1

        if ( face_num <= face_max ) then

          face_material(face_num) = material_num
          face_order(face_num) = 0

          do i = 1, 3
            face_normal(i,face_num) = 0.0E+00
            call word_next_read ( text, word2, done )
            if ( .not. done ) then
              call s_to_r4 ( word2, rval, ierror, lchar )
              if ( ierror == 0 ) then
                face_normal(i,face_num) = rval
              end if
            end if
          end do

        end if

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning!'
        write ( *, '(a)' ) '  Model not in right state for FACET.'
        bad_num = bad_num + 1
        ierror = 1
        exit

      end if

    else if ( s_eqi ( word1, 'ENDFACET' ) ) then

      if ( istate == 2 ) then
        istate = 1
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning!'
        write ( *, '(a)' ) '  Model not in right state for ENDFACET.'
        bad_num = bad_num + 1
        ierror = 1
        exit
      end if

    else if ( s_eqi ( word1, 'OUTERLOOP' ) ) then

      if ( istate == 2 ) then
        istate = 3
        ivert = 0
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning!'
        write ( *, '(a)' ) '  Model not in right state for OUTERLOOP.'
        bad_num = bad_num + 1
        ierror = 1
        exit
      end if
 
    else if ( s_eqi ( word1, 'ENDLOOP' ) ) then

      if ( istate == 3 ) then
        istate = 2
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning!'
        write ( *, '(a)' ) '  Model not in right state for ENDLOOP.'
        bad_num = bad_num + 1
        ierror = 1
        exit
      end if
 
    else if ( s_eqi ( word1, 'VERTEX' ) ) then

      if ( istate == 3 ) then

        do i = 1, 3
          call word_next_read ( text, word2, done )
          call s_to_r4 ( word2, rval, ierror, lchar )
          temp(i) = rval
        end do
!
!  If the coordinate values already exist in COR3, then  
!  save space by using the index of a previous copy.
!
        if ( cor3_num <= 1000 ) then
          call r4col_find ( 3, 3, cor3_num, cor3, temp, icor3 )
        else
          icor3 = 0
        end if

        if ( icor3 == 0 ) then
          cor3_num = cor3_num + 1
          icor3 = cor3_num
          if ( cor3_num <= cor3_max ) then
            cor3(1:3,cor3_num) = temp(1:3)
          end if
        else
          dup_num = dup_num + 1
        end if

        ivert = ivert + 1
 
        if ( ivert <= order_max .and. face_num <= face_max ) then
          face(ivert,face_num) = icor3
          vertex_material(ivert,face_num) = material_num
          vertex_normal(1:3,ivert,face_num) = face_normal(1:3,face_num)
        end if

        if ( face_num <= face_max .and. face_order(face_num) < order_max ) then

          face_order(face_num) = face_order(face_num) + 1

        end if

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'STLA_READ - Warning!'
        write ( *, '(a)' ) '  Model not in right state for VERTEX.'
        bad_num = bad_num + 1
        ierror = 1
        exit

      end if

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STLA_READ - Warning!'
      write ( *, '(a)' ) '  Unrecognized line in file.'
      bad_num = bad_num + 1
      ierror = 1
      exit

    end if
 
  end do
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'STLA_READ - Read ', text_num, &
    ' text lines from ' // trim ( filein_name )

  return
end
subroutine stla_write ( cor3, cor3_max, face, face_max, face_normal, &
  face_num, face_order, filein_name, fileout_name, iunit, order_max )

!*****************************************************************************80
!
!! STLA_WRITE writes graphics information to an ASCII StereoLithography file.
!
!  Example:
!
!    solid MYSOLID
!      facet normal 0.4 0.4 0.2
!        outerloop
!          vertex  1.0 2.1 3.2
!          vertex  2.1 3.7 4.5
!          vertex  3.1 4.5 6.7
!        endloop
!      endfacet
!      ...
!      facet normal 0.2 0.2 0.4
!        outerloop
!          vertex  2.0 2.3 3.4
!          vertex  3.1 3.2 6.5
!          vertex  4.1 5.5 9.0
!        endloop
!      endfacet
!    endsolid MYSOLID
!
!  Discussion:
!
!    The polygons in an STL file should only be triangular.  This routine 
!    will try to automatically decompose higher-order polygonal faces into 
!    suitable triangles, without actually modifying the internal graphics 
!    data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    3D Systems, Inc,
!    Stereolithography Interface Specification,
!    October 1989.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, real ( kind = 4 ) FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) face(order_max,face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_num2
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  character ( len = * ) fileout_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jvert
  integer ( kind = 4 ) node
  integer ( kind = 4 ) text_num

  text_num = 0

  write ( iunit, '(a,a)' ) 'solid MYSOLID created by IVREAD, ' // &
    'original data in ', trim ( filein_name )

  text_num = text_num + 1
  face_num2 = 0

  do iface = 1, face_num

    do jvert = 3, face_order(iface)

      face_num2 = face_num2 + 1

      write ( iunit, '(''  facet normal '', 3g14.6)' ) face_normal(1:3,iface)
      text_num = text_num + 1

      write ( iunit, '(a)' ) '    outer loop'
      text_num = text_num + 1

      node = face(1,iface)
      write ( iunit, '(''      vertex '', 3g14.6)' ) cor3(1:3,node)
      text_num = text_num + 1

      node = face(jvert-1,iface)
      write ( iunit, '(''      vertex '', 3g14.6)' ) cor3(1:3,node)
      text_num = text_num + 1

      node = face(jvert,iface)
      write ( iunit, '(''      vertex '', 3g14.6)' ) cor3(1:3,node)
      text_num = text_num + 1

      write ( iunit, '(a)' ) '  endloop'
      write ( iunit, '(a)' ) 'endfacet'
      text_num = text_num + 2

    end do

  end do

  write ( iunit, '(a)' ) 'endsolid MYSOLID'
  text_num = text_num + 1
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'STLA_WRITE - Wrote ', text_num, &
    ' text lines to ' // trim ( fileout_name )

  if ( face_num2 /= face_num ) then
    write ( *, '(a,i6)' ) '  Number of faces in original data was ', face_num
    write ( *, '(a,i6)' ) &
      '  Number of triangular faces in decomposed data was ', face_num2
  end if

  return
end
subroutine tec_write ( cor3, cor3_material, cor3_max, cor3_num, face, &
  face_max, face_num, face_order, fileout_name, iunit, material_max, &
  material_rgba, order_max )

!*****************************************************************************80
!
!! TEC_WRITE writes graphics information to a TECPLOT file.
!
!  Discussion:
!
!    The file format used is appropriate for 3D finite element
!    surface zone data.  Polygons are decomposed into triangles where
!    necessary.
!
!  Example:
!
!    TITLE = "cube.tec created by IVREAD."
!    VARIABLES = "X", "Y", "Z", "R", "G", "B"
!    ZONE T="TRIANGLES", N=8, E=12, F=FEPOINT, ET=TRIANGLE
!    0.0 0.0 0.0 0.0 0.0 0.0
!    1.0 0.0 0.0 1.0 0.0 0.0
!    1.0 1.0 0.0 1.0 1.0 0.0
!    0.0 1.0 0.0 0.0 1.0 0.0
!    0.0 0.0 1.0 0.0 0.0 1.0
!    1.0 0.0 1.0 1.0 0.0 1.0
!    1.0 1.0 1.0 1.0 1.0 1.0
!    0.0 1.0 1.0 0.0 1.0 1.0
!    1 4 2
!    2 4 3
!    1 5 8
!    1 2 5
!    2 6 5
!    2 3 6
!    3 7 6
!    3 4 7
!    4 8 7
!    4 1 8
!    5 6 8
!    6 7 8
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MATERIAL(COR3_MAX), the material of each node.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, real ( kind = 4 ) MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and A values.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) b
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_material(cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face2(3)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_num2
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) fileout_name
  real ( kind = 4 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) imat
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jlo
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ) text_num
  real ( kind = 4 ) r
!
!  Determine the number of triangular faces.
!
  face_num2 = 0
  do i = 1, face_num
    do jlo = 1, face_order(i) - 2
      face_num2 = face_num2 + 1
    end do
  end do

  text_num = 0

  write ( iunit, '(a)' ) 'TITLE = "' // trim ( fileout_name ) &
    // ' created by IVREAD."'
  write ( iunit, '(a)' ) 'VARIABLES = "X", "Y", "Z", "R", "G", "B"'
  write ( iunit, '(a,i6,a,i6,a)' ) 'ZONE T="TRIANGLES", N=', cor3_num, &
    ', E=', face_num2, ', F=FEPOINT, ET=TRIANGLE'

  text_num = text_num + 3
!
!  Write out X, Y, Z, R, G, B per node.
!
  do icor3 = 1, cor3_num
    imat = cor3_material(icor3)
    r = material_rgba(1,imat)
    g = material_rgba(2,imat)
    b = material_rgba(3,imat)
    write ( iunit, '(6g11.3)' ) cor3(1,icor3), cor3(2,icor3), cor3(3,icor3), &
      r, g, b
    text_num = text_num + 1
  end do
!
!  Do the next face.
!
  do iface = 1, face_num
!
!  Break the face up into triangles, anchored at node 1.
!
    do jlo = 1, face_order(iface) - 2

      face2(1) = face(    1,iface)
      face2(2) = face(jlo+1,iface)
      face2(3) = face(jlo+2,iface)

      write ( iunit, '(3i6)' ) face2(1), face2(2), face2(3)
      text_num = text_num + 1

    end do

  end do
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'TEC_WRITE - Wrote ', text_num, ' text lines to ' &
    // trim ( fileout_name )

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2001
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
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine tmat_init ( a )

!*****************************************************************************80
!
!! TMAT_INIT initializes the geometric transformation matrix.
!
!  Discussion:
!
!    The geometric transformation matrix can be thought of as a 4 by 4
!    matrix "A" having components:
!
!      r11 r12 r13 t1
!      r21 r22 r23 t2
!      r31 r32 r33 t3
!        0   0   0  1
!
!    This matrix encodes the rotations, scalings and translations that
!    are applied to graphical objects.
!
!    A point P = (x,y,z) is rewritten in "homogeneous coordinates" as 
!    PH = (x,y,z,1).  Then to apply the transformations encoded in A to 
!    the point P, we simply compute A * PH.
!
!    Individual transformations, such as a scaling, can be represented
!    by simple versions of the transformation matrix.  If the matrix
!    A represents the current set of transformations, and we wish to 
!    apply a new transformation B, then the original points are
!    transformed twice:  B * ( A * PH ).  The new transformation B can
!    be combined with the original one A, to give a single matrix C that
!    encodes both transformations: C = B * A.
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
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(4,4), the geometric transformation matrix.
!
  implicit none

  real ( kind = 4 ) a(4,4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, 4
    do j = 1, 4
      if ( i == j ) then
        a(i,j) = 1.0E+00
      else
        a(i,j) = 0.0E+00
      end if
    end do
  end do

  return
end
subroutine tmat_mxm ( a, b, c )

!*****************************************************************************80
!
!! TMAT_MXM multiplies two geometric transformation matrices.
!
!  Discussion:
!
!    The product is accumulated in a temporary array, and then assigned
!    to the result.  Therefore, it is legal for any two, or all three,
!    of the arguments to share memory.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(4,4), the first geometric transformation matrix.
!
!    Input, real ( kind = 4 ) B(4,4), the second geometric transformation matrix.
!
!    Output, real ( kind = 4 ) C(4,4), the product A * B.
!
  implicit none

  real ( kind = 4 ) a(4,4)
  real ( kind = 4 ) b(4,4)
  real ( kind = 4 ) c(4,4)

  c(1:4,1:4) = matmul ( a(1:4,1:4), b(1:4,1:4) )

  return
end
subroutine tmat_mxp ( a, x, y )

!*****************************************************************************80
!
!! TMAT_MXP multiplies a geometric transformation matrix times a point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(4,4), the geometric transformation matrix.
!
!    Input, real ( kind = 4 ) X(3), the point to be multiplied.  The fourth component
!    of X is implicitly assigned the value of 1.
!
!    Output, real ( kind = 4 ) Y(3), the result of A*X.  The product is accumulated in 
!    a temporary vector, and then assigned to the result.  Therefore, it 
!    is legal for X and Y to share memory.
!
  implicit none

  real ( kind = 4 ) a(4,4)
  real ( kind = 4 ) x(3)
  real ( kind = 4 ) y(3)

  y(1:3) = a(1:3,4) + matmul ( a(1:3,1:3), x(1:3) )

  return
end
subroutine tmat_mxp2 ( a, x, y, n )

!*****************************************************************************80
!
!! TMAT_MXP2 multiplies a geometric transformation matrix times N points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(4,4), the geometric transformation matrix.
!
!    Input, real ( kind = 4 ) X(3,N), the points to be multiplied.  
!
!    Output, real ( kind = 4 ) Y(3,N), the transformed points.  Each product is 
!    accumulated in a temporary vector, and then assigned to the
!    result.  Therefore, it is legal for X and Y to share memory.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 4 ) a(4,4)
  integer ( kind = 4 ) k
  real ( kind = 4 ) x(3,n)
  real ( kind = 4 ) y(3,n)

  do k = 1, n

    y(1:3,k) = a(1:3,4) + matmul ( a(1:3,1:3), x(1:3,k) )

  end do

  return
end
subroutine tmat_mxv ( a, x, y )

!*****************************************************************************80
!
!! TMAT_MXV multiplies a geometric transformation matrix times a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(4,4), the geometric transformation matrix.
!
!    Input, real ( kind = 4 ) X(3), the vector to be multiplied.  The fourth component
!    of X is implicitly assigned the value of 1.
!
!    Output, real ( kind = 4 ) Y(3), the result of A*X.  The product is accumulated in 
!    a temporary vector, and then assigned to the result.  Therefore, it 
!    is legal for X and Y to share memory.
!
  implicit none

  real ( kind = 4 ) a(4,4)
  real ( kind = 4 ) x(3)
  real ( kind = 4 ) y(3)

  y(1:3) = a(1:3,4) + matmul ( a(1:3,1:3), x(1:3) )

  return
end
subroutine tmat_rot_axis ( a, b, angle, axis )

!*****************************************************************************80
!
!! TMAT_ROT_AXIS applies a coordinate axis rotation to the geometric transformation matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(4,4), the current geometric transformation matrix.
!
!    Output, real ( kind = 4 ) B(4,4), the modified geometric transformation matrix.
!    A and B may share the same memory.
!
!    Input, real ( kind = 4 ) ANGLE, the angle, in degrees, of the rotation.
!
!    Input, character AXIS, is 'X', 'Y' or 'Z', specifying the coordinate
!    axis about which the rotation occurs.
!
  implicit none

  real ( kind = 4 ) a(4,4)
  real ( kind = 4 ) angle
  real ( kind = 4 ) angle_rad
  character axis
  real ( kind = 4 ) b(4,4)
  real ( kind = 4 ) c(4,4)
  real ( kind = 4 ) d(4,4)
  real ( kind = 4 ) degrees_to_radians
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  angle_rad = degrees_to_radians ( angle )

  call tmat_init ( c )

  if ( axis == 'X' .or. axis == 'x' ) then
    c(2,2) =   cos ( angle_rad )
    c(2,3) = - sin ( angle_rad )
    c(3,2) =   sin ( angle_rad )
    c(3,3) =   cos ( angle_rad )
  else if ( axis == 'Y' .or. axis == 'y' ) then
    c(1,1) =   cos ( angle_rad )
    c(1,3) =   sin ( angle_rad )
    c(3,1) = - sin ( angle_rad )
    c(3,3) =   cos ( angle_rad )
  else if ( axis == 'Z' .or. axis == 'z' ) then
    c(1,1) =   cos ( angle_rad )
    c(1,2) = - sin ( angle_rad )
    c(2,1) =   sin ( angle_rad )
    c(2,2) =   cos ( angle_rad )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TMAT_ROT_AXIS - Fatal error!'
    write ( *, '(a)' ) '  Illegal rotation axis: ' // axis
    write ( *, '(a)' ) '  Legal choices are ''X'', ''Y'', or ''Z''.'
    return
  end if

  call tmat_mxm ( c, a, d )

  b(1:4,1:4) = d(1:4,1:4)

  return
end
subroutine tmat_rot_vector ( a, b, angle, axis )

!*****************************************************************************80
!
!! TMAT_ROT_VECTOR: arbitrary axis rotation to geometric transformation matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(4,4), the current geometric transformation matrix.
!
!    Output, real ( kind = 4 ) B(4,4), the modified geometric transformation matrix.
!    A and B may share the same memory.
!
!    Input, real ( kind = 4 ) ANGLE, the angle, in degrees, of the rotation.
!
!    Input, real ( kind = 4 ) AXIS(3), the axis vector about which rotation occurs.
!    AXIS may not be the zero vector.
!
  implicit none

  real ( kind = 4 ) a(4,4)
  real ( kind = 4 ) angle
  real ( kind = 4 ) angle_rad
  real ( kind = 4 ) axis(3)
  real ( kind = 4 ) b(4,4)
  real ( kind = 4 ) c(4,4)
  real ( kind = 4 ) ca
  real ( kind = 4 ) d(4,4)
  real ( kind = 4 ) degrees_to_radians
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 4 ) norm
  real ( kind = 4 ) sa
  real ( kind = 4 ) v1
  real ( kind = 4 ) v2
  real ( kind = 4 ) v3

  v1 = axis(1)
  v2 = axis(2)
  v3 = axis(3)

  norm = sqrt ( v1 * v1 + v2 * v2 + v3 * v3 )

  if ( norm == 0.0E+00 ) then
    return
  end if

  v1 = v1 / norm
  v2 = v2 / norm
  v3 = v3 / norm

  angle_rad = degrees_to_radians ( angle )
  ca = cos ( angle_rad )
  sa = sin ( angle_rad )

  call tmat_init ( c )

  c(1,1) =                    v1 * v1 + ca * ( 1.0E+00 - v1 * v1 )
  c(1,2) = ( 1.0E+00 - ca ) * v1 * v2 - sa * v3
  c(1,3) = ( 1.0E+00 - ca ) * v1 * v3 + sa * v2

  c(2,1) = ( 1.0E+00 - ca ) * v2 * v1 + sa * v3
  c(2,2) =                    v2 * v2 + ca * ( 1.0E+00 - v2 * v2 )
  c(2,3) = ( 1.0E+00 - ca ) * v2 * v3 - sa * v1

  c(3,1) = ( 1.0E+00 - ca ) * v3 * v1 - sa * v2
  c(3,2) = ( 1.0E+00 - ca ) * v3 * v2 + sa * v1
  c(3,3) =                    v3 * v3 + ca * ( 1.0E+00 - v3 * v3 )

  call tmat_mxm ( c, a, d )

  b(1:4,1:4) = d(1:4,1:4)

  return
end
subroutine tmat_scale ( a, b, sx, sy, sz )

!*****************************************************************************80
!
!! TMAT_SCALE applies a scaling to the geometric transformation matrix.
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
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(4,4), the current geometric transformation matrix.
!
!    Output, real ( kind = 4 ) B(4,4), the modified geometric transformation matrix.
!    A and B may share the same memory.
!
!    Input, real ( kind = 4 ) SX, SY, SZ, the scalings to be applied to the X, Y and
!    Z coordinates.
!
  implicit none

  real ( kind = 4 ) a(4,4)
  real ( kind = 4 ) b(4,4)
  real ( kind = 4 ) c(4,4)
  real ( kind = 4 ) d(4,4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 4 ) sx
  real ( kind = 4 ) sy
  real ( kind = 4 ) sz

  call tmat_init ( c )

  c(1,1) = sx
  c(2,2) = sy
  c(3,3) = sz

  call tmat_mxm ( c, a, d )

  b(1:4,1:4) = d(1:4,1:4)

  return
end
subroutine tmat_shear ( a, b, axis, s )

!*****************************************************************************80
!
!! TMAT_SHEAR applies a shear to the geometric transformation matrix.
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
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(4,4), the current geometric transformation matrix.
!
!    Output, real ( kind = 4 ) B(4,4), the modified geometric transformation matrix.
!    A and B may share the same memory.
!
!    Input, character ( len = 2 ) AXIS, is 'XY', 'XZ', 'YX', 'YZ', 'ZX' or 'ZY',
!    specifying the shear equation:
!
!      XY:  x' = x + s * y;
!      XZ:  x' = x + s * z;
!      YX:  y' = y + s * x;
!      YZ:  y' = y + s * z;
!      ZX:  z' = z + s * x;
!      ZY:  z' = z + s * y.
!
!    Input, real ( kind = 4 ) S, the shear coefficient.
!
  implicit none

  real ( kind = 4 ) a(4,4)
  character ( len = 2 ) axis
  real ( kind = 4 ) b(4,4)
  real ( kind = 4 ) c(4,4)
  real ( kind = 4 ) d(4,4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 4 ) s

  call tmat_init ( c )

  if ( axis == 'XY' .or. axis == 'xy' ) then
    c(1,2) = s
  else if ( axis == 'XZ' .or. axis == 'xz' ) then
    c(1,3) = s
  else if ( axis == 'YX' .or. axis == 'yx' ) then
    c(2,1) = s
  else if ( axis == 'YZ' .or. axis == 'yz' ) then
    c(2,3) = s
  else if ( axis == 'ZX' .or. axis == 'zx' ) then
    c(3,1) = s
  else if ( axis == 'ZY' .or. axis == 'zy' ) then
    c(3,2) = s
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TMAT_SHEAR - Fatal error!'
    write ( *, '(a)' ) '  Illegal shear axis: ' // axis
    write ( *, '(a)' ) '  Legal choices are XY, XZ, YX, YZ, ZX, or ZY.'
    return
  end if

  call tmat_mxm ( c, a, d )

  b(1:4,1:4) = d(1:4,1:4)

  return
end
subroutine tmat_trans ( a, b, x, y, z )

!*****************************************************************************80
!
!! TMAT_TRANS applies a translation to the geometric transformation matrix.
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
!  Reference:
!
!    Foley, van Dam, Feiner, Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) A(4,4), the current geometric transformation matrix.
!
!    Output, real ( kind = 4 ) B(4,4), the modified transformation matrix.
!    A and B may share the same memory.
!
!    Input, real ( kind = 4 ) X, Y, Z, the translation.  This may be thought of as the
!    point that the origin moves to under the translation.
!
  implicit none

  real ( kind = 4 ) a(4,4)
  real ( kind = 4 ) b(4,4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  real ( kind = 4 ) z

  do i = 1, 4
    do j = 1, 4

      if ( i == 1 .and. j == 4 ) then
        b(1,4) = a(1,4) + x
      else if ( i == 2 .and. j == 4 ) then
        b(2,4) = a(2,4) + y
      else if ( i == 3 .and. j == 4 ) then
        b(3,4) = a(3,4) + z
      else
        b(i,j) = a(i,j)
      end if

    end do
  end do

  return
end
subroutine tria_read ( cor3, cor3_max, cor3_num, dup_num, face, &
  face_material, face_max, face_num, face_order, filein_name, ierror, &
  iunit, order_max, text_num, vertex_material, vertex_normal )

!*****************************************************************************80
!
!! TRIA_READ reads graphics information from an ASCII triangle file.
!
!  Discussion:
!
!    It is intended that the information read from the file can
!    either start a whole new graphics object, or simply be added
!    to a current graphics object via the '<<' command.
!
!    This is controlled by whether the input values have been zeroed
!    out or not.  This routine simply tacks on the information it
!    finds to the current graphics object.
!
!  Example:
!
!    12                    <-- Number of triangles
!
!                          (x,y,z) and (nx,ny,nz) of normal vector at:
!
!    0.0 0.0 0.0 0.3 0.3 0.3   node 1 of triangle 1.
!    1.0 0.0 0.0 0.3 0.1 0.3   node 2 of triangle 1,
!    0.0 1.0 0.0 0.3 0.1 0.3   node 3 of triangle 1,
!    1.0 0.5 0.0 0.3 0.1 0.3   node 1 of triangle 2,
!    ...
!    0.0 0.5 0.5 0.3 0.1 0.3   node 3 of triangle 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input/output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices 
!    per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which data is read.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Output, integer ( kind = 4 ) DUP_NUM, the number of duplicate nodes discovered.
!
!    Input/output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Output, integer ( kind = 4 ) TEXT_NUM, the number of lines of input text.
!
!    Input/output, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input/output, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals 
!    at vertices.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  real ( kind = 4 ) cvec(3)
  integer ( kind = 4 ) dup_num
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_num2
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) iface_hi
  integer ( kind = 4 ) iface_lo
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ivert
  real ( kind = 4 ) r1
  real ( kind = 4 ) r2
  real ( kind = 4 ) r3
  real ( kind = 4 ) r4
  real ( kind = 4 ) r5
  real ( kind = 4 ) r6
  integer ( kind = 4 ) text_num
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
!
  ierror = 0
!
!  Read the number of (triangular) faces.
!  (This is added on to the current number, if any).
!
  read ( iunit, *, iostat = ios ) face_num2

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIA_READ - Warning.'
    write ( *, '(a)' ) '  End-of-file, but the model was not finished.'
    return
  end if

  text_num = text_num + 1
!
!  For each triangle.
!
  iface_lo = face_num + 1
  iface_hi = face_num + face_num2

  do iface = iface_lo, iface_hi

    if ( iface <= FACE_MAX ) then
      face_order(iface) = 3
      face_material(iface) = 1
    end if
!
!  For each face of a triangle:
!
    do ivert = 1, face_order(iface)

      read ( iunit, *, iostat = ios ) r1, r2, r3, r4, r5, r6

      if ( ios /= 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIA_READ - Warning.'
        write ( *, '(a)' ) '  End-of-file, but the model was not finished.'
        return
      end if

      text_num = text_num + 1

      cvec(1) = r1
      cvec(2) = r2
      cvec(3) = r3

      if ( cor3_num <= 1000 ) then
        call r4col_find ( 3, 3, cor3_num, cor3, cvec, icor3 )
      else
        icor3 = 0
      end if
 
      if ( icor3 == 0 ) then

        cor3_num = cor3_num + 1
        icor3 = cor3_num

        if ( cor3_num <= cor3_max ) then
          cor3(1:3,cor3_num) = cvec(1:3)
        end if

      else

        dup_num = dup_num + 1

      end if

      if ( iface <= FACE_MAX ) then

        face(ivert,iface) = icor3

        vertex_material(ivert,iface) = 1
        vertex_normal(1,ivert,iface) = r4
        vertex_normal(2,ivert,iface) = r5
        vertex_normal(3,ivert,iface) = r6

      end if

    end do

    face_num = iface

  end do
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'TRIA_READ - Read ', text_num, &
    ' text lines from ' // trim ( filein_name )

  return
end
subroutine tria_write ( cor3, cor3_max, cor3_normal, face, face_max, face_num, &
  face_order, fileout_name, iunit, order_max )

!*****************************************************************************80
!
!! TRIA_WRITE writes the graphics data to an ASCII "triangle" file.
!
!  Discussion:
!
!    This is just a private format that Greg Hood requested from me.
!
!  Example:
!
!    12                        <-- Number of triangles
!
!                              (x,y,z) and (nx,ny,nz) of normal vector at:
!
!    0.0 0.0 0.0 0.3 0.3 0.3   node 1 of triangle 1.
!    1.0 0.0 0.0 0.3 0.1 0.3   node 2 of triangle 1,
!    0.0 1.0 0.0 0.3 0.1 0.3   node 3 of triangle 1,
!    1.0 0.5 0.0 0.3 0.1 0.3   node 1 of triangle 2,
!    ...
!    0.0 0.5 0.5 0.3 0.1 0.3   node 3 of triangle 12.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, real ( kind = 4 ) COR3_NORMAL(3,COR3_MAX), the normal vector at each node.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3(3,cor3_max)
  real ( kind = 4 ) cor3_normal(3,cor3_max)
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_num2
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) face2(3)
  character ( len = * ) fileout_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  real ( kind = 4 ) nx
  real ( kind = 4 ) ny
  real ( kind = 4 ) nz
  integer ( kind = 4 ) text_num
  real ( kind = 4 ) x
  real ( kind = 4 ) y
  real ( kind = 4 ) z
!
  text_num = 0
!
!  Determine the number of triangular faces.
!
  face_num2 = 0
  do i = 1, face_num
    do jlo = 1, face_order(i) - 2
      face_num2 = face_num2 + 1
    end do
  end do

  write ( iunit, '(i6)' ) face_num2
  text_num = text_num + 1
!
!  Do the next face.
!
  do iface = 1, face_num
!
!  Break the face up into triangles, anchored at node 1.
!
    do jlo = 1, face_order(iface) - 2

      face2(1) = face(    1,iface)
      face2(2) = face(jlo+1,iface)
      face2(3) = face(jlo+2,iface)

      do k = 1, 3

        icor3 = face2(k)

        x = cor3(1,icor3)
        y = cor3(2,icor3)
        z = cor3(3,icor3)

        nx = cor3_normal(1,icor3)
        ny = cor3_normal(2,icor3)
        nz = cor3_normal(3,icor3)

        write ( iunit, '(6f10.4)' ) x, y, z, nx, ny, nz

        text_num = text_num + 1

      end do

    end do

  end do
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'TRIA_WRITE - Wrote ', text_num, &
    ' text lines to ' // trim ( fileout_name )

  return
end
subroutine ts_read ( bad_num, cor3, cor3_max, cor3_num, face, face_material, &
  face_max, face_num, face_order, filein_name, ierror, iunit, line_dex, &
  line_material, line_max, line_num, material_max, material_num, &
  material_rgba, order_max, point, point_max, point_num )

!*****************************************************************************80
!
!! TS_READ reads graphics information from a Mathematica TS file.
!
!  Example:
!
!    % Graphics3D objects
!    boundingbox
!    0.0 0.0 0.0
!    0.0 1.0 1.0
!    viewpoint
!    1.3 -2.4 2.0
!    ambientlight
!    0.0 0.0 0.0
!    lightsources <-- pairs of positions and RGB colors
!    1.0 0.0 1.0
!    1.0 0.0 0.0
!    1.0 1.0 1.0
!    0.0 1.0 0.0
!    0.0 1.0 1.0
!    0.0 0.0 0.1
!    polygon
!    0.0 0.0 0.0
!    0.0 1.0 0.0
!    0.0 1.0 1.0
!    point
!    0.1 0.4 0.3
!    line
!    0.1 0.1 0.1
!    0.1 0.9 0.1
!    0.9 0.9 0.1
!    0.9 0.1 0.1
!    0.1 0.1 0.1
!    line
!    0.1 0.1 0.9
!    0.1 0.9 0.9
!    0.9 0.9 0.9
!    0.9 0.1 0.9
!    0.1 0.1 0.9
!    mesh 3 2  <-- 3 by 2 mesh
!    0.0       <-- Z(1,1)
!    1.0       <-- Z(1,2)
!    0.5       <-- Z(2,1)
!    1.0       <-- Z(2,2)
!    1.0       <-- Z(3,1)
!    2.0       <-- Z(3,2)
!    colormesh 3 2 <-- 3 by 2 mesh, and 2 by 1 central colors
!    0.0
!    1.0
!    0.5
!    1.0
!    1.0
!    2.0
!    1.0 0.0 0.0 <-- RGB for Z(1,1),Z(1,2),Z(2,1),Z(2,2)
!    0.0 1.0 0.0 <-- RGB for Z(2,1),Z(2,2),Z(3,1),Z(3,2)
!    color
!    0.0 0.0 1.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) BAD_NUM, the number of bad text lines encountered.
!
!    Input/output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input/output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input/output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices 
!    per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which data is read.
!
!    Input/output, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated
!    by -1.
!
!    Input/output, integer ( kind = 4 ) LINE_MATERIAL(LINE_MAX), material index for 
!    each line.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input/output, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Input/output, character ( len = * ) MATERIAL_NAME(MATERIAL_MAX), 
!    material names.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Output, integer ( kind = 4 ) TEXT_NUM, the number of lines of text read from
!    the file.
!
!    Input/output, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) point_max

  real ( kind = 4 ) b
  integer ( kind = 4 ) bad_num
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  logical done
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  real ( kind = 4 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstep
  integer ( kind = 4 ) lchar
  character ( len = 256 ) line
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) m
  integer ( kind = 4 ) material_default
! character ( len = * ) material_name(material_max)
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  character ( len = 100 ) mode
  integer ( kind = 4 ) n
! integer object_num
  integer ( kind = 4 ), parameter :: offset = 1
  integer ( kind = 4 ) point(point_max)
  integer ( kind = 4 ) point_num
  real ( kind = 4 ) r
  logical s_eqi
  integer ( kind = 4 ) step
  integer ( kind = 4 ) text_num
  character ( len = 256 ) word
  real ( kind = 4 ) x
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin
  real ( kind = 4 ) y
  real ( kind = 4 ) ymax
  real ( kind = 4 ) ymin
  real ( kind = 4 ) z
  real ( kind = 4 ) zmax
  real ( kind = 4 ) zmin
!
  ierror = 0
!
!  Save a copy of the input value of COR3_NUM to use as a base.
!
  bad_num = 0
  material_default = 1
  text_num = 0
  word = ' '
  mode = 'UNKNOWN'
!
!  Read a line of text from the file.
!
  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TSREAD - Error!'
      write ( *, '(a)' ) '  IERROR nonzero!'
      write ( *, '(a)' ) trim ( line )
      write ( *, '(a,i6)' ) 'text_num = ', text_num
    end if

    text_num = text_num + 1
!
!  Replace any control characters (in particular, TAB's) by blanks.
!
    call s_control_blank ( line )

    done = .true.
!
!  Try to read the first word from the line.
!
    call word_next_read ( line, word, done )
!
!  If no words in this line, cycle.
!
    if ( done ) then
      cycle
    end if
!
!  If this word begins with '%', then it's a comment.  Cycle.
!
    if ( word(1:1) == '%' ) then

      mode = 'COMMENT'

    else if ( s_eqi ( word, 'AMBIENTLIGHT' ) ) then

      mode = word
      step = 0

    else if ( s_eqi ( word, 'BOUNDINGBOX' ) ) then

      mode = word
      step = 0

    else if ( s_eqi ( word, 'COLOR' ) ) then

      mode = word
      step = 0

    else if ( s_eqi ( word, 'COLORMESH' ) ) then

      mode = word
      call word_next_read ( line, word, done )
      call s_to_i4 ( word, m, ierror, lchar )
      call word_next_read ( line, word, done )
      call s_to_i4 ( word, n, ierror, lchar )
      istep = 1
      jstep = 1

    else if ( s_eqi ( word, 'LIGHTSOURCES' ) ) then

      mode = word
      step = 0

    else if ( s_eqi ( word, 'LINE' ) ) then

      mode = word
      step = 0

    else if ( s_eqi ( word, 'MESH' ) ) then

      mode = word
      call word_next_read ( line, word, done )
      call s_to_i4 ( word, m, ierror, lchar )
      call word_next_read ( line, word, done )
      call s_to_i4 ( word, n, ierror, lchar )

      step = 0

    else if ( s_eqi ( word, 'POINT' ) ) then

      mode = word
      step = 0

    else if ( s_eqi ( word, 'POLYGON' ) ) then

      mode = word
      step = 0

    else if ( s_eqi ( word, 'VIEWPOINT' ) ) then

      mode = word
      step = 0
!
!  Otherwise, we expect this to be numeric data.
!
    else if ( s_eqi ( mode, 'AMBIENTLIGHT' ) ) then

      if ( step == 0 ) then

        call s_to_r4 ( word, r, ierror, lchar )
        if ( ierror /= 0 ) then
          bad_num = bad_num + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) 'TS_READ: Bad AMBIENTLIGHT data on line ', &
            text_num
          write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        end if
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, g, ierror, lchar )
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, b, ierror, lchar )

        material_rgba(1:3,1) = (/ r, g, b /)

      else
        bad_num = bad_num + 1
        if ( bad_num <= 10 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) 'TS_READ: Bad AMBIENTLIGHT data on line ', &
            text_num
          write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        end if
      end if

      step = step + 1

    else if ( s_eqi ( mode, 'BOUNDINGBOX' ) ) then

      if ( step == 0 ) then

        call s_to_r4 ( word, xmin, ierror, lchar )
        if ( ierror /= 0 ) then
          bad_num = bad_num + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) 'TS_READ: Bad BOUNDINGBOX XMIN data on line ', &
            text_num
          write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        else
          write ( *, '(a,g14.6)' ) 'XMIN = ', xmin
        end if

        call word_next_read ( line, word, done )
        call s_to_r4 ( word, ymin, ierror, lchar )
        if ( ierror /= 0 ) then
          bad_num = bad_num + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) 'TS_READ: Bad BOUNDINGBOX YMIN data on line ', &
            text_num
          write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        else
          write ( *, '(a,g14.6)' ) 'YMIN = ', ymin
        end if

        call word_next_read ( line, word, done )
        call s_to_r4 ( word, zmin, ierror, lchar )
        if ( ierror /= 0 ) then
          bad_num = bad_num + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) 'TS_READ: Bad BOUNDINGBOX ZMIN data on line ', &
            text_num
          write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        else
          write ( *, '(a,g14.6)' ) 'ZMIN = ', zmin
        end if

      else if ( step == 1 ) then

        call s_to_r4 ( word, xmax, ierror, lchar )
        if ( ierror /= 0 ) then
          bad_num = bad_num + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) 'TS_READ: Bad BOUNDINGBOX XMAX data on line ', &
            text_num
          write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        end if
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, ymax, ierror, lchar )
        if ( ierror /= 0 ) then
          bad_num = bad_num + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) 'TS_READ: Bad BOUNDINGBOX YMAX data on line ', &
            text_num
          write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        end if
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, zmax, ierror, lchar )
        if ( ierror /= 0 ) then
          bad_num = bad_num + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) 'TS_READ: Bad BOUNDINGBOX ZMAX data on line ', &
            text_num
          write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        end if
      else
        bad_num = bad_num + 1
        if ( bad_num <= 10 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) 'TS_READ: Bad BOUNDINGBOX data on line ', &
            text_num
          write ( *, '(a,i6)' ) '  Step = ', step
          write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        end if
      end if

      step = step + 1
!
!  This command should generate a new material,
!  and it resets the material default.
!
    else if ( s_eqi ( mode, 'COLOR' ) ) then

      if ( step == 0 ) then

        call s_to_r4 ( word, r, ierror, lchar )
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, g, ierror, lchar )
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, b, ierror, lchar )

        material_num = material_num + 1
        material_rgba(1:3,material_num) = (/ r, g, b /)
        material_default = material_num

      else

        bad_num = bad_num + 1
        if ( bad_num <= 10 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) 'TS_READ: Bad COLOR data on line ', text_num
          write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        end if

      end if

      step = step + 1
!
!  Hold off on this until we've got MESH working.
!
!  The colors defined here should not affect the default material.
!
    else if ( s_eqi ( mode, 'COLORMESH' ) ) then

      if ( step == 0 ) then
        call mesh_t3 ( face, face_max, face_num, face_order, m, n, order_max )
      end if

      if ( step < m * n ) then

        j = j + 1

        if ( n < j ) then
          j = 1
          i = i + 1
        end if

        x = real ( i ) * xmax / real ( m )
        y = real ( j ) * ymax / real ( n )
        call s_to_r4 ( word, z, ierror, lchar )

        cor3_num = cor3_num + 1
        cor3(1:3,cor3_num) = (/ x, y, z /)

      else

        if ( step == m * n ) then
          j = 0
          i = 1
!
!         Assign faces 1,2 material_num + 1
!                      3,4 material_num + 2
!                      2*(m-1)*(n-1)-1,2*(m-1)*(n-1) material_num+(m-1)*(n-1)
!
        end if

        j = j + 1

        if ( n - 1 < j ) then
          j = 1
          i = i + 1
        end if

        call s_to_r4 ( word, r, ierror, lchar )
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, g, ierror, lchar )
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, b, ierror, lchar )

        material_num = material_num + 1
        material_rgba(1:3,material_num) = (/ r, g, b /)
        material_default = material_num

      end if

      step = step + 1

    else if ( s_eqi ( mode, 'LIGHTSOURCES' ) ) then

      if ( mod(step,2) == 0 ) then
        call s_to_r4 ( word, x, ierror, lchar )
        if ( ierror /= 0 ) then
          bad_num = bad_num + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) 'TS_READ: Bad LIGHTSOURCES data on line ', &
            text_num
          write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        end if
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, y, ierror, lchar )
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, z, ierror, lchar )
      else
        call s_to_r4 ( word, r, ierror, lchar )
        if ( ierror /= 0 ) then
          bad_num = bad_num + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) 'TS_READ: Bad LIGHTSOURCES RGB data on line ', &
            text_num
          write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        end if
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, g, ierror, lchar )
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, b, ierror, lchar )
      end if

      step = step + 1
!
!  Use the default material, whatever that is, for the line.
!
    else if ( s_eqi ( mode, 'LINE' ) ) then

      if ( step == 0 ) then
        line_num = line_num + 1
        line_dex(line_num) = -1 + offset
      end if

      call s_to_r4 ( word, x, ierror, lchar )
      call word_next_read ( line, word, done )
      call s_to_r4 ( word, y, ierror, lchar )
      call word_next_read ( line, word, done )
      call s_to_r4 ( word, z, ierror, lchar )

      cor3_num = cor3_num + 1
      cor3(1:3,cor3_num) = (/ x, y, z /)

      line_dex(line_num) = cor3_num
      line_material(line_num) = material_default
      line_num = line_num + 1
      line_dex(line_num) = -1
!
!  This is trick because we have to run through the set of indices.
!
    else if ( s_eqi ( mode, 'MESH' ) ) then

      if ( step == 0 ) then
        call mesh_t3 ( face, face_max, face_num, face_order, m, n, order_max )
      end if

      j = j + 1

      if ( n < j ) then
        j = 1
        i = i + 1
      end if

      x = real ( i ) * xmax / real ( m )
      y = real ( j ) * ymax / real ( n )
      call s_to_r4 ( word, z, ierror, lchar )

      cor3_num = cor3_num + 1
      cor3(1:3,cor3_num) = (/ x, y, z /)

      step = step + 1

    else if ( s_eqi ( mode, 'POINT' ) ) then

      if ( step == 0 ) then
        call s_to_r4 ( word, x, ierror, lchar )
        if ( ierror /= 0 ) then
          bad_num = bad_num + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) 'TS_READ: Bad POINT data on line ', text_num
          write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        end if
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, y, ierror, lchar )
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, z, ierror, lchar )

        cor3_num = cor3_num + 1
        cor3(1:3,cor3_num) = (/ x, y, z /)

        point_num = point_num + 1
        point(point_num) = cor3_num

      else
        bad_num = bad_num + 1
        if ( bad_num <= 10 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) 'TS_READ: Bad POINT data on line ', text_num
          write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        end if
      end if

      step = step + 1

    else if ( s_eqi ( mode, 'POLYGON' ) ) then

      if ( step == 0 ) then
        face_num = face_num + 1
        face_order(face_num) = 0
      end if

      call s_to_r4 ( word, x, ierror, lchar )
        if ( ierror /= 0 ) then
          bad_num = bad_num + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) 'TS_READ: Bad POLYGON data on line ', text_num
          write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        end if
      call word_next_read ( line, word, done )
      call s_to_r4 ( word, y, ierror, lchar )
      call word_next_read ( line, word, done )
      call s_to_r4 ( word, z, ierror, lchar )

      cor3_num = cor3_num + 1
      cor3(1:3,cor3_num) = (/ x, y, z /)

      face_order(face_num) = step + 1
      face(step+1,face_num) = cor3_num

      step = step + 1

    else if ( s_eqi ( mode, 'VIEWPOINT' ) ) then

      if ( step == 0 ) then
        call s_to_r4 ( word, x, ierror, lchar )
        if ( ierror /= 0 ) then
          bad_num = bad_num + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) 'TS_READ: Bad VIEWPOINT data on line ', text_num
          write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        end if
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, y, ierror, lchar )
        call word_next_read ( line, word, done )
        call s_to_r4 ( word, z, ierror, lchar )
      else
        bad_num = bad_num + 1
        if ( bad_num <= 10 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) 'TS_READ: Bad VIEWPOINT data on line ', text_num
          write ( *, '(a)' ) '  Bad word: ' // trim ( word )
        end if
      end if

      step = step + 1
!
!  Unrecognized keyword.
!
    else

      bad_num = bad_num + 1

      if ( bad_num <= 10 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) 'TS_READ: Bad data on line ', text_num
        write ( *, '(a)' ) '  Bad word: ' // trim ( word )
      end if

    end if

  end do
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'TS_READ - Read ', text_num, &
    ' text lines from ' // trim ( filein_name )

  return
end
subroutine ts_write ( cor3, cor3_max, cor3_num, face, face_material, &
  face_max, face_num, face_order, filein_name, fileout_name, iunit, &
  line_dex, line_material, line_max, line_num, material_max, material_num, &
  material_rgba, order_max, point, point_max, point_num )

!*****************************************************************************80
!
!! TS_WRITE writes graphics information to a Mathematica TS file.
!
!  Example:
!
!    % Graphics3D objects
!    boundingbox
!    0 0 0
!    0 1 1
!    viewpoint
!    1.3 -2.4 2.
!    ambientlight
!    0 0 0
!    lightsources
!    1. 0. 1.
!    1 0 0
!    1. 1. 1.
!    0 1 0
!    0. 1. 1.
!    0 0 1
!    polygon
!    0 0 0
!    0 1 0
!    0 1 1
!    point
!    0.1 0.4 0.2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, real ( kind = 4 ) POINT(POINT_MAX), the indices of points to display.
!
!    Input, integer ( kind = 4 ) POINT_MAX, the maximum number of points to display.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points to display.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) point_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_i
  integer ( kind = 4 ) face_material(face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  character ( len = * ) fileout_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_i
  integer ( kind = 4 ) line_line_i_old
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) material_default
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ), parameter :: OFFSET = 1
  integer ( kind = 4 ) point(point_max)
  integer ( kind = 4 ) point_num
  character ( len = 256 ) text
  integer ( kind = 4 ) text_num
  character ( len = 256 ) text2
  integer ( kind = 4 ) v
  integer ( kind = 4 ) vert_i
  real ( kind = 4 ) xmax
  real ( kind = 4 ) xmin
  real ( kind = 4 ) xview
  real ( kind = 4 ) ymax
  real ( kind = 4 ) ymin
  real ( kind = 4 ) yview
  real ( kind = 4 ) zmax
  real ( kind = 4 ) zmin
  real ( kind = 4 ) zview
!
  text_num = 0
  material_default = 1
!
!  Header.
!
  write ( iunit, '(a)' ) '% Graphics3D objects'
  text_num = text_num + 1
!
!  BOUNDINGBOX
!
  xmax = maxval ( cor3(1,1:cor3_num) )
  xmin = minval ( cor3(1,1:cor3_num) )
  ymax = maxval ( cor3(2,1:cor3_num) )
  ymin = minval ( cor3(2,1:cor3_num) )
  zmax = maxval ( cor3(3,1:cor3_num) )
  zmin = minval ( cor3(3,1:cor3_num) )

  write ( iunit, '(a)' ) 'boundingbox'
  write ( iunit, '(3g14.6)' ) xmin, ymin, zmin
  write ( iunit, '(3g14.6)' ) xmax, ymax, zmax
  text_num = text_num + 3
!
!  VIEWPOINT
!    x, y, z
!
  xview = 5.0 * xmax
  yview = 3.5 * ymax
  zview = 2.0 * zmax
  write ( iunit, '(a)' ) 'viewpoint'
  write ( iunit, '(3g14.6)' ) xview, yview, zview
  text_num = text_num + 2
!
!  AMBIENTLIGHT
!    r, g, b
!
  write ( iunit, '(a)' ) 'ambientlight'
  write ( iunit, '(3g14.6)' ) material_rgba(1:3,1)
  text_num = text_num + 2
!
!  LIGHTSOURCES
!    x, y, z
!    r, g, b
!
  write ( iunit, '(a)' ) 'lightsources'
  write ( iunit, '(3g14.6)' ) 2.0 * xmax, 0.0,        2.0 * zmax
  write ( iunit, '(3g14.6)' ) 1.0,        0.0,        0.0
  write ( iunit, '(3g14.6)' ) 2.0 * xmax, 2.0 * ymax, 2.0 * zmax
  write ( iunit, '(3g14.6)' ) 0.0,        1.0,        0.0
  write ( iunit, '(3g14.6)' ) 0.0,        2.0 * ymax, 2.0 * zmax
  write ( iunit, '(3g14.6)' ) 1.0,        0.0,        1.0
  text_num = text_num + 7
!
!  POLYGON
!  For each face, write the coordinates of the vertices.
!
  do face_i = 1, face_num

    if ( face_material(face_i) /= material_default ) then
      material_default = face_material(face_i)
      write ( iunit, '(a)' ) 'color'
      write ( iunit, '(3g14.6)' ) material_rgba(1:3,material_default)
      text_num = text_num + 1
    end if

    write ( iunit, '(a)' ) 'polygon'
    text_num = text_num + 1

    do vert_i = 1, face_order(face_i)
      v = face(vert_i,face_i)
      write ( iunit, '(3g14.6)' ) cor3(1:3,v)
      text_num = text_num + 1
    end do

  end do
!
!  LINE
!
  line_i = 0
  line_line_i_old = -1 + OFFSET

  do while ( line_i < line_num )

    line_i = line_i + 1

    if ( line_dex(line_i) /= -1 + OFFSET .and. &
      line_line_i_old == -1 + OFFSET ) then

      if ( line_material(line_i) /= material_default ) then
        if ( 0 < line_material(line_i) .and. &
          line_material(line_i) <= material_num ) then
          material_default = line_material(line_i)
          write ( iunit, '(a)' ) 'color'
          write ( iunit, '(3g14.6)') material_rgba(1:3,material_default)
          text_num = text_num + 2
        end if
      end if

      write ( iunit, '(a)' ) 'line'
      text_num = text_num + 1

    end if

    if ( line_dex(line_i) /= -1 + OFFSET ) then
      v = line_dex(line_i)
      write ( iunit, '(3g14.6)' ) cor3(1:3,v)
      text_num = text_num + 1
    end if

    line_line_i_old = line_dex(line_i)

  end do
!
!  POINT
!
  do i = 1, point_num
    write ( iunit, '(a)' ) 'point'
    v = point(i)
    write ( iunit, '(3g14.6)' ) cor3(1:3,v)
    text_num = text_num + 2
  end do
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'TS_WRITE - Wrote ', text_num, &
    ' text lines to ' // trim ( fileout_name )

  return
end
subroutine txt_write ( cor3, cor3_material, cor3_max, cor3_normal, cor3_num, &
  cor3_tex_uv, face, face_material, face_max, face_normal, face_num, &
  face_order, face_tex_uv, filein_name, fileout_name, iunit, line_dex, &
  line_material, line_max, line_num, material_max, material_name, &
  material_num, material_rgba, object_name, order_max, point, point_max, &
  point_num, texture_max, texture_name, texture_num, vertex_material, &
  vertex_normal, vertex_tex_uv )

!*****************************************************************************80
!
!! TXT_WRITE writes the graphics data to a text file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MATERIAL(COR3_MAX), the material of each node.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, real ( kind = 4 ) COR3_NORMAL(3,COR3_MAX), normals at nodes.  
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, real ( kind = 4 ) COR3_TEX_UV(2,COR3_MAX), UV texture coordinates for nodes.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, real ( kind = 4 ) FACE_NORMAL(3,FACE_MAX), the normal vector at each face.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer ( kind = 4 ) LINE_MATERIAL(LINE_MAX), material index for each line.
!
!    Input, character ( len = * ) MATERIAL_NAME(MATERIAL_MAX), material names.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer ( kind = 4 ) TEXTURE_MAX, the maximum number of textures.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Input, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Input, integer ( kind = 4 ) TEXTURE_NUM, the number of textures.
!
!    Input, character ( len = * ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
!    Input, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Input, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.  
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) point_max
  integer ( kind = 4 ) texture_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_material(cor3_max)
  real ( kind = 4 ) cor3_normal(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  real ( kind = 4 ) cor3_tex_uv(2,cor3_max)
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  real ( kind = 4 ) face_normal(3,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  real ( kind = 4 ) face_tex_uv(2,face_max)
  character ( len = * ) filein_name
  character ( len = * ) fileout_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) imat
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) length
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_num
  character ( len = * ) material_name(material_max)
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  character ( len = * ) object_name
  integer ( kind = 4 ), parameter :: OFFSET = 1
  integer ( kind = 4 ) point(point_max)
  integer ( kind = 4 ) point_num
  character ( len = 100 ) text
  integer ( kind = 4 ) text_num
  character ( len = * ) texture_name(texture_max)
  integer ( kind = 4 ) texture_num
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  real ( kind = 4 ) vertex_tex_uv(2,order_max,face_max)
  character ( len = 10 ) word
!
  text_num = 0

  write ( iunit, '(a)' ) trim ( fileout_name ) // ' created by IVREAD.'
  write ( iunit, '(a)' ) 'Original data in ' // trim ( filein_name ) // '.'
  write ( iunit, '(a)' ) 'Object name is ' // trim ( object_name ) // '.'
  text_num = text_num + 3
!
!  NODES.
!
  write ( iunit, * ) ' '
  write ( iunit, * ) cor3_num, ' nodes:'
  text_num = text_num + 2

  if ( 0 < cor3_num ) then

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Node  Coordinates:'
    write ( iunit, '(a)' ) '========  ============'
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do j = 1, cor3_num
      write ( iunit, '(i8,2x,3g12.4)' ) j-OFFSET, cor3(1:3,j)
      text_num = text_num + 1
    end do

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Node  Normal vector'
    write ( iunit, '(a)' ) '========  ============'
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do j = 1, cor3_num
      write ( iunit, '(i8,2x,3g12.4)' ) j-OFFSET, cor3_normal(1:3,j)
      text_num = text_num + 1
    end do

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Node  Material'
    write ( iunit, '(a)' ) '========  ============'
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do j = 1, cor3_num
      write ( iunit, '(i8,2x,i8)' ) j-OFFSET, cor3_material(j)-OFFSET
      text_num = text_num + 1
    end do

    if ( 0 < texture_num ) then

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Node  Texture coordinates:'
    write ( iunit, '(a)' ) '========  ===================='
    write ( iunit, * ) ' '
    text_num = text_num + 4

      do j = 1, cor3_num
        write ( iunit, '(i8,2x,3g12.4)' ) j-OFFSET, cor3_tex_uv(1:2,j)
        text_num = text_num + 1
      end do

    end if

  end if
!
!  POINTS
!
  write ( iunit, * ) ' '
  write ( iunit, * ) point_num, ' point data items.'
  text_num = text_num + 2
  do i = 1, point_num
    write ( iunit, * ) i, point(i)
  end do
!
!  LINES
!
  write ( iunit, * ) ' '
  write ( iunit, * ) line_num, ' line data items.'
  text_num = text_num + 2

  if ( 0 < line_num ) then

    write ( iunit, * ) ' '
    write ( iunit, * ) '  Line index data:'
    write ( iunit, * ) ' '
    text_num = text_num + 3

    text = ' '
    length = 0
      
    do j = 1, line_num
 
      write ( word, '(i8,'','')' ) line_dex(j) - OFFSET
      call s_cat ( text, word, text )
      length = length + 1

      if ( line_dex(j)-OFFSET == -1 .or. 10 <= length .or. line_num <= j ) then

        call s_blanks_delete ( text )
        write ( iunit, '(8x,a)' ) trim ( text )
        text_num = text_num + 1
        text = ' '
        length = 0

      end if

    end do

    write ( iunit, * ) ' '
    write ( iunit, * ) 'Line materials:'
    write ( iunit, * ) ' '
    text_num = text_num + 3

    text = ' '
    length = 0
     
    do j = 1, line_num
 
      write ( word, '(i8,'','')' ) line_material(j) - OFFSET
      call s_cat ( text, word, text )
      length = length + 1

      if ( line_material(j)-OFFSET == -1 .or. &
           10 <= length .or. line_num <= j ) then

        call s_blanks_delete ( text )
        write ( iunit, '(8x,a)' ) trim ( text )
        text_num = text_num + 1
        text = ' '
        length = 0

      end if

    end do

  end if
!
!  FACES
!
  write ( iunit, * ) ' '
  write ( iunit, * ) face_num, ' faces:'
  text_num = text_num + 2

  if ( 0 < face_num ) then

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Face  Material     Order'
    write ( iunit, '(a)' ) '========  ========  ========'
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do iface = 1, face_num
      write ( iunit, '(i8,2x,i8,2x,i8)' ) iface-OFFSET, &
        face_material(iface)-OFFSET, face_order(iface)
      text_num = text_num + 1
    end do

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Face  Vertices:'
    write ( iunit, '(a)' ) '========  ========================'
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do iface = 1, face_num

      do jlo = 1, face_order(iface), 10
        jhi = min ( jlo + 9, face_order(iface) )
        if ( jlo == 1 ) then
          write ( iunit, '(i8,2x,10i8)' ) iface-OFFSET, &
                ( face(ivert,iface)-OFFSET, ivert = jlo, jhi )
        else
          write ( iunit, '(10x,10i8)' ) &
                ( face(ivert,iface)-OFFSET, ivert = jlo, jhi )
        end if
        text_num = text_num + 1
      end do

    end do

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Face    Vertex Material Indices'
    write ( iunit, '(a)' ) '========  ========================='
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do iface = 1, face_num

      do jlo = 1, face_order(iface), 10
        jhi = min ( jlo + 9, face_order(iface) )
        if ( jlo == 1 ) then
          write ( iunit, '(i8,2x,10i8)' ) iface-OFFSET, &
                ( vertex_material(ivert,iface)-OFFSET, ivert = jlo, jhi )
        else
          write ( iunit, '(10x,10i8)' ) &
                ( vertex_material(ivert,iface)-OFFSET, ivert = jlo, jhi )
        end if
        text_num = text_num + 1
      end do

    end do

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Face  Normal vector'
    write ( iunit, '(a)' ) '========  ========================'
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do iface = 1, face_num
      write ( iunit, '(i8,3g12.4)' ) iface-OFFSET, face_normal(1:3,iface)
      text_num = text_num + 1
    end do

    if ( 0 < texture_num ) then

      write ( iunit, * ) ' '
      write ( iunit, '(a)' ) '    Face  Texture coordinates:'
      write ( iunit, '(a)' ) '========  ========================'
      write ( iunit, * ) ' '
      text_num = text_num + 4

      do iface = 1, face_num
        write ( iunit, '(i8,2x,3g12.4)' ) iface-OFFSET, face_tex_uv(1:2,iface)
        text_num = text_num + 1
      end do

    end if

  end if
!
!  VERTICES.
!
  if ( 0 < face_num ) then

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Face    Vertex  Normal vector:'
    write ( iunit, '(a)' ) '========  ========  =============='
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do iface = 1, face_num
      write ( iunit, '(a)' ) ' '
      text_num = text_num + 1
      do ivert = 1, face_order(iface)
        write ( iunit, '(i8,2x,i8,2x,3f12.4)' ) iface-OFFSET, ivert-OFFSET, &
          vertex_normal(1:3,ivert,iface)
        text_num = text_num + 1
      end do
    end do

    if ( 0 < texture_num ) then

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) '    Face    Vertex  Texture coordinates:'
    write ( iunit, '(a)' ) '========  ========  =================='
    write ( iunit, * ) ' '
    text_num = text_num + 4

      do iface = 1, face_num
        write ( iunit, '(a)' ) ' '
        text_num = text_num + 1
        do ivert = 1, face_order(iface)
          write ( iunit, '(i8,2x,i8,2x,3f12.4)' ) iface-OFFSET, ivert-OFFSET, &
            vertex_tex_uv(1:2,ivert,iface)
          text_num = text_num + 1
        end do
      end do

    end if

  end if
!
!  MATERIALS
!
  write ( iunit, * ) ' '
  write ( iunit, * ) material_num, ' materials:'
  write ( iunit, * ) ' '
  text_num = text_num + 3

  if ( 0 < material_num ) then
    write ( iunit, '(a)' ) 'Material        Name            ' // &
      '    R         G         B         A'
    write ( iunit, '(a)' ) '========  ====================  ' // &
      '========================================'
    write ( iunit, * ) ' '
    do imat = 1, material_num
      write ( iunit, '(i8,2x,a20,2x,4f10.6)' ) imat-OFFSET, &
        material_name(imat)(1:20), material_rgba(1:4,imat)
    end do
  end if
!
!  TEXTURES:
!
  write ( iunit, * ) ' '
  write ( iunit, * ) texture_num, ' textures:'
  text_num = text_num + 2

  if ( 0 < texture_num ) then

    write ( iunit, * ) ' '
    write ( iunit, '(a)' ) ' Texture    Name'
    write ( iunit, '(a)' ) '========  ========'
    write ( iunit, * ) ' '
    text_num = text_num + 4

    do i = 1, texture_num
      write ( iunit, '(i6,2x,a)' ) i-OFFSET, trim ( texture_name(i) )
      text_num = text_num + 1
    end do
  end if
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'TXT_WRITE - Wrote ', text_num, &
    ' text lines to ' // trim ( fileout_name )
 
  return
end
subroutine ucd_write ( cor3, cor3_material, cor3_max, cor3_num, face, &
  face_material, face_max, face_num, face_order, fileout_name, iunit, &
  material_max, material_num, material_rgba, order_max )

!*****************************************************************************80
!
!! UCD_WRITE writes graphics data to an AVS UCD file.
!
!  Example:
!
!    #  cube.ucd created by IVREAD.
!    #
!    #  Material RGB to hue map:
!    #
!    #  material    R    G      B   Alpha     Hue
!    #
!    #    0       0.94  0.70  0.15  1.000   0.116
!    #    1       0.24  0.70  0.85  1.000   0.541
!    #    2       0.24  0.00  0.85  1.000   0.666
!    #
!    #  The node data is
!    #    node # / material # / RGBA / Hue
!    #
!    8  6  6  0  0
!    0  0.0  0.0  0.0E+00
!    1  1.0  0.0  0.0E+00
!    2  1.0  1.0  0.0E+00
!    3  0.0  1.0  0.0E+00
!    4  0.0  0.0  1.0E+00
!    5  1.0  0.0  1.0E+00
!    6  1.0  1.0  1.0E+00
!    7  0.0  1.0  1.0E+00
!    0  0  quad  0  1  2  3
!    1  0  quad  0  4  5  1 
!    2  0  quad  1  5  6  2
!    3  0  quad  2  6  7  3
!    4  0  quad  3  7  4  0
!    5  0  quad  4  7  6  5
!    3  1 4 1
!    material, 0...2
!    RGBA, 0-1/0-1/0-1/0-1
!    Hue, 0-1
!    0  0  0.94  0.70  0.15  1.0  0.116
!    1  0  0.94  0.70  0.15  1.0  0.116
!    2  0  0.94  0.70  0.15  1.0  0.116
!    3  0  0.94  0.70  0.15  1.0  0.116
!    4  1  0.24  0.70  0.85  1.0  0.541
!    5  1  0.24  0.70  0.85  1.0  0.541
!    6  2  0.24  0.24  0.85  0.0  0.666
!    7  2  0.24  0.24  0.85  0.0  0.666
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Input, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Input, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) a
  real ( kind = 4 ) b
  character ( len = 4 ) cell_type
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_material(cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) fileout_name
  real ( kind = 4 ) g
  real ( kind = 4 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imat
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ), parameter :: OFFSET = 1
  real ( kind = 4 ) r
  character ( len = 100 ) text
  integer ( kind = 4 ) text_num
!
  text_num = 0

  write ( iunit, '(a)' ) '#  ' // trim ( fileout_name ) // ' created by IVREAD.'
  write ( iunit, '(a)' ) '#'
  write ( iunit, '(a)' ) '#  Material RGB to Hue map:'
  write ( iunit, '(a)' ) '#'
  write ( iunit, '(a)' ) '#  material    R    G      B     Alpha  Hue'
  write ( iunit, '(a)' ) '#'
  text_num = text_num + 6

  do j = 1, material_num
    r = material_rgba(1,j)
    g = material_rgba(2,j)
    b = material_rgba(3,j)
    a = material_rgba(4,j)
    call rgb_to_hue ( r, g, b, h )
    write ( text, '(a,2x,i4,5f7.3)' ) '#  ', j-OFFSET, r, g, b, a, h
    call s_blanks_delete ( text )
    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1
  end do

  write ( iunit, '(a)' ) '#'
  write ( iunit, '(a)' ) '#  The node data is'
  write ( iunit, '(a)' ) '#    node # / material # / RGBA / Hue'
  write ( iunit, '(a)' ) '#'
  text_num = text_num + 4

  write ( text, '(i6,2x,i6,2x,a)' ) cor3_num, face_num, '6  0  0'
  call s_blanks_delete ( text )
  write ( iunit, '(a)' ) trim ( text )
  text_num = text_num + 1

  do j = 1, cor3_num
    write ( text, '(i6,2x,3f7.3)' ) j - OFFSET, cor3(1:3,j)
    call s_blanks_delete ( text )
    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1
  end do
!
!  NOTE:
!  UCD only accepts triangles and quadrilaterals, not higher order
!  polygons.  We would need to break polygons up to proceed.
!
!  Also, we just use the material of vertex 1 of the face as the
!  material of the face.
!
  do j = 1, face_num

    if ( face_order(j) == 3 ) then
      cell_type = 'tri'
    else if ( face_order(j) .eq. 4 ) then
      cell_type = 'quad'
    else
      cell_type = '???'
    end if

    imat = face_material(j) - OFFSET
    write ( text, '(i6,2x,i6,2x,a4,2x,10i6)' ) j-OFFSET, imat, &
      cell_type, ( face(i,j) - OFFSET, i = 1, face_order(j) )
    call s_blanks_delete ( text )
    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1

  end do

  write ( iunit, '(a)' )  '3  1  4  1'
  write ( iunit, '(a,i6)' )  'material, 0...', material_num - OFFSET
  write ( iunit, '(a)' )  'RGBA, 0-1/0-1/0-1/0-1'
  write ( iunit, '(a)' )  'Hue, 0-1'
  text_num = text_num + 4

  do j = 1, cor3_num
    imat = cor3_material(j)
    r = material_rgba(1,imat)
    g = material_rgba(2,imat)
    b = material_rgba(3,imat)
    a = material_rgba(4,imat)
    call rgb_to_hue ( r, g, b, h )
    write ( text, '(i6,2x,i6,2x,5f7.3)' ) j-OFFSET, imat-OFFSET, r, g, b, a, h
    call s_blanks_delete ( text )
    write ( iunit, '(a)' ) trim ( text )
    text_num = text_num + 1
  end do
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'UCD_WRITE - Wrote ', text_num, &
    ' text lines to ' // trim ( fileout_name )


  return
end
subroutine vector_unit_nd ( n, v )

!*****************************************************************************80
!
!! VECTOR_UNIT_ND normalizes a vector in ND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, real ( kind = 4 ) V(N), the vector to be normalized.  On output,
!    V should have unit Euclidean norm.  However, if the input vector
!    has zero Euclidean norm, it is not altered.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 4 ) enorm_nd
  integer ( kind = 4 ) i
  real ( kind = 4 ) temp
  real ( kind = 4 ) v(n)

  temp = enorm_nd ( n, v )
 
  if ( temp /= 0.0E+00 ) then
    v(1:n) = v(1:n) / temp
  end if
 
  return
end
subroutine vertex_normal_set ( cor3, cor3_max, face, face_max, &
  face_num, face_order, order_max, vertex_normal )

!*****************************************************************************80
!
!! VERTEX_NORMAL_SET recomputes the face vertex normal vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input/output, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normal vectors
!    at vertices.  
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3(3,cor3_max)
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) jp2
  real ( kind = 4 ) norm
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) fix_num
  integer ( kind = 4 ) vec_num
  integer ( kind = 4 ) zero_num
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  real ( kind = 4 ) x0
  real ( kind = 4 ) x1
  real ( kind = 4 ) x2
  real ( kind = 4 ) xc
  real ( kind = 4 ) y0
  real ( kind = 4 ) y1
  real ( kind = 4 ) y2
  real ( kind = 4 ) yc
  real ( kind = 4 ) z0
  real ( kind = 4 ) z1
  real ( kind = 4 ) z2
  real ( kind = 4 ) zc

  if ( face_num <= 0 ) then
    return
  end if

  fix_num = 0
  vec_num = 0
  zero_num = 0
!
!  Consider each face.
!
  do iface = 1, face_num
!
!  Consider each vertex.
!
    do ivert = 1, face_order(iface)

      vec_num = vec_num + 1

      norm = sqrt ( sum ( vertex_normal(1:3,ivert,iface)**2 ) )

      if ( norm == 0.0E+00 ) then

        fix_num = fix_num + 1

        i0 = face(ivert,iface)
        x0 = cor3(1,i0)
        y0 = cor3(2,i0)
        z0 = cor3(3,i0)

        jp1 = ivert + 1
        if ( face_order(iface) < jp1 ) then
          jp1 = jp1 - face_order(iface)
        end if
        i1 = face(jp1,iface)
        x1 = cor3(1,i1)
        y1 = cor3(2,i1)
        z1 = cor3(3,i1)

        jp2 = ivert + 2
        if ( face_order(iface) < jp2 ) then
          jp2 = jp2 - face_order(iface)
        end if
        i2 = face(jp2,iface)
        x2 = cor3(1,i2)
        y2 = cor3(2,i2)
        z2 = cor3(3,i2)

        call cross0_3d ( x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc )

        norm = sqrt ( xc * xc + yc * yc + zc * zc )

        if ( norm == 0.0E+00 ) then
          zero_num = zero_num + 1
          xc = 1.0E+00 / sqrt ( 3.0E+00 )
          yc = 1.0E+00 / sqrt ( 3.0E+00 )
          zc = 1.0E+00 / sqrt ( 3.0E+00 )
        else
          xc = xc / norm
          yc = yc / norm
          zc = zc / norm
        end if

        vertex_normal(1:3,ivert,iface) = (/ xc, yc, zc /)

      end if

    end do

  end do

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VERTEX_NORMAL_SET - Note:'
    write ( *, '(a,i6)' ) '  Number of vertex normal vectors is ', vec_num
    write ( *, '(a,i6)' ) '  Zero vectors were ', fix_num
    write ( *, '(a,i6)' ) '  Zero areas were ', zero_num
  end if

  return
end
subroutine vertex_to_node_material ( cor3_material, cor3_max, face, &
  face_order, face_max, order_max, face_num, vertex_material )

!*****************************************************************************80
!
!! VERTEX_TO_NODE_MATERIAL extends vertex material definitions to nodes.
!
!  Discussion:
!
!    A NODE is a point in space.
!    A VERTEX is a node as used in a particular face.
!    One node may be used as a vertex in several faces, or none.
!    This routine simply runs through all the vertices, and assigns
!    the material of the vertex to the corresponding node.  If a
!    node appears as a vertex several times, then the node will 
!    end up having the material of the vertex that occurs "last".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) COR3_MATERIAL(COR3_MAX), the material index of each node.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) VERTEX_MAT(ORDER_MAX,FACE_MAX), vertex materials.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) order_max

  integer ( kind = 4 ) cor3_material(cor3_max)
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_order(face_max)
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) node
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) vertex_material(order_max,face_max)

  do iface = 1, face_num
    do ivert = 1, face_order(iface)
      node = face(ivert,iface)
      cor3_material(node) = vertex_material(ivert,iface)
    end do
  end do

  return
end
subroutine vla_read ( bad_num, cor3, cor3_max, cor3_num, dup_num, &
  filein_name, ierror, iunit, line_dex, line_material, line_max, &
  line_num, material_num, text_num )

!*****************************************************************************80
!
!! VLA_READ reads graphics information from a VLA file.
!
!  Comments:
!
!    Internal comments begin with a semicolon in column 1.
!
!    The X, Y, Z coordinates of points begin with a "P" to
!    denote the beginning of a line, and "L" to denote the
!    continuation of a line.  The fourth entry is intensity, which 
!    should be between 0.0 and 1.0.
!
!  Discussion:
!
!    It is intended that the information read from the file can
!    either start a whole new graphics object, or simply be added
!    to a current graphics object via the '<<' command.  
!
!    This is controlled by whether the input values have been zeroed 
!    out or not.  This routine simply tacks on the information it 
!    finds to the current graphics object.
!
!  Example:
!
!     set comment cube.vla created by IVREAD
!     set comment from data in file cube.iv
!     set comment
!     set intensity EXPLICIT
!     set parametric NON_PARAMETRIC
!     set filecontent LINES
!     set filetype NEW
!     set depthcue 0
!     set defaultdraw stellar
!     set coordsys RIGHT
!     set author IVREAD
!     set site Buhl Planetarium
!     set library_id UNKNOWN
!     ; DXF LINE entity
!     P   8.59816       5.55317      -3.05561       1.00000
!     L   8.59816       2.49756      0.000000E+00   1.00000
!     L   8.59816       2.49756      -3.05561       1.00000
!     L   8.59816       5.55317      -3.05561       1.00000
!     ; DXF LINE entity
!     P   8.59816       5.55317      0.000000E+00   1.00000
!     ...etc...
!     L   2.48695       2.49756      -3.05561       1.00000
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which data is read.
!
!    Output, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Output, integer ( kind = 4 ) LINE_MATERIAL(LINE_MAX), material index for each line.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line items.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) line_max

  integer ( kind = 4 ) bad_num
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  real ( kind = 4 ) cvec(3)
  logical done
  integer ( kind = 4 ) dup_num
  character ( len = * ) filein_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) iword
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) material_num
  integer ( kind = 4 ), parameter :: OFFSET = 1
  character prevcode
  real ( kind = 4 ) rval
  logical s_eqi
  character ( len = 256 ) text
  integer ( kind = 4 ) text_num
  character ( len = 256 ) word
  character ( len = 256 ) word1

  ierror = 0
  prevcode = 'P'

10    continue
 
  read ( iunit, '(a)', iostat = ios ) text

  if ( ios /= 0 ) then
    go to 30
  end if

  text_num = text_num + 1
 
  done = .true.
  iword = 0
 
20    continue
 
  call word_next_read ( text, word, done )
!
!  If no more words in this line, read a new line.
!
  if ( done ) then
    go to 10
  end if
 
  iword = iword + 1
!
!  The first word in the line tells us what's happening.
!
  if ( iword == 1 ) then
    word1 = word
  end if
!
!  If WORD1 is "SET", then we regard this line as comments.
!
  if ( s_eqi ( word1, 'set' ) ) then
!
!  If WORD1 is ";", then we regard this line as comments.
!
  else if ( word1 == ';' ) then
!
!  If WORD1 is "P", then this is the initial point on a line.
!  If WORD1 is "L", then this is a followup point on a line.
!
  else if ( s_eqi ( word1, 'P' ) .or. s_eqi ( word1, 'L' ) ) then
!
!  Terminate the current line if the new code is 'P'.
!
    if ( s_eqi ( prevcode, 'L' ) .and. s_eqi ( word1, 'P' ) ) then

      line_num = line_num + 1

      if ( line_num <= line_max ) then
        line_dex(line_num) = -1 + OFFSET
        line_material(line_num) = -1 + OFFSET
      end if

    end if

    prevcode = word1
!
!  Read in the coordinates of the point, stored temporarily in CVEC.
!
    do i = 1, 3

      call word_next_read ( text, word, done )

      if ( done ) then
        bad_num = bad_num + 1
        go to 10
      end if

      call s_to_r4 ( word, rval, ierror, lchar )

      if ( ierror /= 0 ) then
        bad_num = bad_num + 1
        go to 10
      end if

      if ( cor3_num <= cor3_max ) then
        cvec(i) = rval
      end if

    end do
!
!  If the values in CVEC already exist in COR3, then  
!  save space by using the index of a previous copy.
!
!  Otherwise, add CVEC to COR3, and increment COR3_NUM.
!
    if ( cor3_num <= 1000 ) then
      call r4col_find ( 3, 3, cor3_num, cor3, cvec, icor3 )
    else
      icor3 = 0
    end if

    if ( icor3 == 0 ) then
      cor3_num = cor3_num + 1
      icor3 = cor3_num
      if ( cor3_num <= cor3_max ) then
        cor3(1:3,cor3_num) = cvec(1:3)
      end if
    else
      dup_num = dup_num + 1
    end if
!
!  Now define the line.
!
    line_num = line_num + 1

    if ( line_num <= line_max ) then
      line_dex(line_num) = icor3 - 1 + OFFSET
      line_material(line_num) = material_num
    end if

    go to 10
!
!  If the first word is unrecognized, then skip the whole line.
!
  else

    bad_num = bad_num + 1
    go to 10

  end if
 
  go to 20

30    continue
!
!  Terminate the very last line.
!
  if ( 0 < line_num ) then

    line_num = line_num + 1

    if ( line_num <= line_max ) then
      line_dex(line_num) = -1 + OFFSET
      line_material(line_num) = -1 + OFFSET
    end if

  end if
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'VLA_READ - Read ', text_num, &
    ' text lines from ' // trim ( filein_name )

  return
end
subroutine vla_write ( cor3, cor3_max, filein_name, fileout_name, iunit, &
  line_dex, line_max, line_num )

!*****************************************************************************80
!
!! VLA_WRITE writes graphics data to a VLA file.
!
!  Discussion:
!
!    Comments begin with a semicolon in column 1.
!    The X, Y, Z coordinates of points begin with a "P" to
!    denote the beginning of a line, and "L" to denote the
!    continuation of a line.  The fourth entry is intensity, which 
!    should be between 0.0 and 1.0.
!
!  Example:
!
!     set comment cube.vla created by IVREAD
!     set comment Original data in cube.iv.
!     set comment
!     set intensity EXPLICIT
!     set parametric NON_PARAMETRIC
!     set filecontent LINES
!     set filetype NEW
!     set depthcue 0
!     set defaultdraw stellar
!     set coordsys RIGHT
!     set author IVREAD
!     set site Buhl Planetarium
!     set library_id UNKNOWN
!     P   8.59816       5.55317      -3.05561       1.00000
!     L   8.59816       2.49756      0.000000E+00   1.00000
!     L   8.59816       2.49756      -3.05561       1.00000
!     L   8.59816       5.55317      -3.05561       1.00000
!     P   8.59816       5.55317      0.000000E+00   1.00000
!     ...etc...
!     L   2.48695       2.49756      -3.05561       1.00000
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) line_max

  character code
  real ( kind = 4 ) cor3(3,cor3_max)
  character ( len = * ) filein_name
  character ( len = * ) fileout_name
  integer ( kind = 4 ) i
  real ( kind = 4 ) intense
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_num
  logical newline
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) np
  integer ( kind = 4 ), parameter :: OFFSET = 1
  integer ( kind = 4 ) text_num

  intense = 1.0E+00
  np = 0
  nl = 0
 
  write ( iunit, '(a)' ) 'set comment ' // trim ( fileout_name ) // &
    ' created by IVREAD'

  write ( iunit, '(a)' ) 'set comment Original data in ' &
    // trim ( filein_name ) // '.'

  write ( iunit, '(a)' ) 'set comment'
  write ( iunit, '(a)' ) 'set intensity EXPLICIT'
  write ( iunit, '(a)' ) 'set parametric NON_PARAMETRIC'
  write ( iunit, '(a)' ) 'set filecontent LINES'
  write ( iunit, '(a)' ) 'set filetype NEW'
  write ( iunit, '(a)' ) 'set depthcue 0'
  write ( iunit, '(a)' ) 'set defaultdraw stellar'
  write ( iunit, '(a)' ) 'set coordsys RIGHT'
  write ( iunit, '(a)' ) 'set author IVREAD'
  write ( iunit, '(a)' ) 'set site Buhl Planetarium'
  write ( iunit, '(a)' ) 'set library_id UNKNOWN'

  text_num = 13
!
!  Print the line data.
!
  newline = .TRUE.

  do j = 1, line_num
 
    k = line_dex(j) - OFFSET

    if ( k == - 1 ) then

      newline = .TRUE.

    else 

      if ( newline ) then
        code = 'P'
        np = np + 1
      else
        code = 'L'
        nl = nl + 1
      end if

      write ( iunit, '(a,4g12.4)' ) code, cor3(1:3,k+OFFSET), intense

      text_num = text_num + 1
      newline = .FALSE.

    end if

  end do
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'VLA_WRITE - Wrote ', text_num, &
    ' text lines to ' // trim ( fileout_name )
 
  return
end
subroutine vrml_read ( bad_num, cor3, cor3_material, cor3_max, cor3_num, &
  face, face_material, face_max, face_num, face_order, filein_name, &
  ierror, iunit, line_dex, line_material, line_max, line_num, &
  material_max, material_name, material_num, material_rgba, &
  order_max, text_num, texture_max, texture_name, texture_num, &
  vertex_material, vertex_normal, vertex_tex_uv )

!*****************************************************************************80
!
!! VRML_READ reads graphics information from a VRML file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input/output, integer ( kind = 4 ) COR3_MATERIAL(COR3_MAX), the material index 
!    of each node.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input/output, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input/output, integer ( kind = 4 ) FACE_MATERIAL(FACE_MAX), the material of each face.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices 
!    per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which data is read.
!
!    Output, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Output, integer ( kind = 4 ) LINE_MATERIAL(LINE_MAX), material index for line.
!
!    Output, character ( len = * ) MATERIAL_NAME(MATERIAL_MAX), 
!    material names.
!
!    Output, real ( kind = 4 ) MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and A values.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
!    Input, integer ( kind = 4 ) TEXTURE_MAX, the maximum number of textures.
!
!    Output, integer ( kind = 4 ) BAD_NUM, the number of "bad" lines of input text.
!
!    Output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Output, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Output, integer ( kind = 4 ) TEXTURE_NUM, the number of textures.
!
!    Output, integer ( kind = 4 ) TEXT_NUM, the number of lines of input text.
!
!    Output, character ( len = * ) TEXTURE_NAME(TEXTURE_MAX), texture names.
!
!    Output, integer ( kind = 4 ) VERTEX_MATERIAL(ORDER_MAX,FACE_MAX), vertex materials.
!
!    Output, real ( kind = 4 ) VERTEX_NORMAL(3,ORDER_MAX,FACE_MAX), normals at vertices.
!
!    Output, real ( kind = 4 ) VERTEX_TEX_UV(2,ORDER_MAX,FACE_MAX), vertex texture
!    coordinates.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ), parameter :: ivec_max = 20
  integer ( kind = 4 ), parameter :: level_max = 10
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ), parameter :: r4vec_max = 20
  integer ( kind = 4 ) texture_max

  real ( kind = 4 ) angle
  real ( kind = 4 ) axis(3)
  integer ( kind = 4 ) bad_num
  character ( len = 4 ) char4
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_material(cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) cor3_num_old
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_num_old
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) imat
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivec(ivec_max)
  integer ( kind = 4 ) ivec_num
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) level
  character ( len = 256 ) level_name(0:level_max)
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) line_num_old
  character ( len = 100 ) material_binding
  character ( len = * ) material_name(material_max)
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ) new_color
  integer ( kind = 4 ) new_color_index
  integer ( kind = 4 ) cor3_new
  integer ( kind = 4 ) new_face
  integer ( kind = 4 ) node
  integer ( kind = 4 ) nlbrack
  integer ( kind = 4 ) nrbrack
  integer ( kind = 4 ), parameter :: offset = 1
  integer ( kind = 4 ) overall_mat
  real ( kind = 4 ) r01
  real ( kind = 4 ) r02
  real ( kind = 4 ) r03
  real ( kind = 4 ) r04
  real ( kind = 4 ) r05
  real ( kind = 4 ) r06
  real ( kind = 4 ) r07
  real ( kind = 4 ) r08
  real ( kind = 4 ) r09
  real ( kind = 4 ) r10
  real ( kind = 4 ) r11
  real ( kind = 4 ) r12
  real ( kind = 4 ) rval
  real ( kind = 4 ) rvec(r4vec_max)
  integer ( kind = 4 ) r4vec_num
  real ( kind = 4 ) rx
  real ( kind = 4 ) ry
  real ( kind = 4 ) rz
  logical s_eqi
  real ( kind = 4 ) sx
  real ( kind = 4 ) sy
  real ( kind = 4 ) sz
  character ( len = 256 ) text
  integer ( kind = 4 ) text_num
  character ( len = * ) texture_name(texture_max)
  integer ( kind = 4 ) texture_num
  real ( kind = 4 ) transform_matrix(4,4)
  real ( kind = 4 ) tx
  real ( kind = 4 ) ty
  real ( kind = 4 ) tz
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  real ( kind = 4 ) vertex_normal(3,order_max,face_max)
  real ( kind = 4 ) vertex_tex_uv(2,order_max,face_max)
  character ( len = 256 ) word
  character ( len = 256 ) wordm1

  ierror = 0
  level = 0
  level_name(0) = 'Top'
  nlbrack = 0
  nrbrack = 0
  ivec_num = 0
  r4vec_num = 0
!
!  Save old counts in order to be able to add a new object to a set
!  of existing ones.
!
  cor3_num_old = cor3_num
  face_num_old = face_num
  line_num_old = line_num
!
!  Initialize the transformation matrix.
!
  call tmat_init ( transform_matrix )

  word = ' '
  text = ' '
!
!  Read the next word and its trailing context, from the input file.
!
10    continue

  wordm1 = word

11    continue

  call file_get_next_word ( iunit, word, text, text_num, ierror )

  if ( ierror /= 0 ) then
    go to 50
  end if

  if ( debug ) then
    write ( *, '(a)' ) trim ( word )
  end if
!
!  The first line of the file must be the header.
!
  if ( text_num == 1 ) then

    if ( s_eqi ( word, '#VRML' ) ) then
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      if ( s_eqi ( word, 'V2.0' ) ) then
        call file_get_next_word ( iunit, word, text, text_num, ierror )
        if ( s_eqi ( word, 'utf8' ) ) then
          go to 10
        end if
      end if
    end if

    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VRML_READ - Fatal error!'
    write ( *, '(a)' ) '  The input file has a bad header.'
    write ( *, '(a)' ) trim ( text )
    return

  end if
!
!  Skip a comment and all text following it on a line.
!
  if ( word(1:1) == '#' ) then
    text = ' '
    go to 10
  end if
!
!  Ignore a blank line.
!
  if ( word == ' ' ) then
    go to 10
  end if
!
!  Ignore an isolated comma.
!
  if ( word == ',' ) then
    go to 10
  end if
!
!  If the word is a curly or square bracket, count it.
!
  if ( word == '{' .or. word == '[' ) then

    nlbrack = nlbrack + 1

  else if ( word .eq. '}' .or. word == ']' ) then

    nrbrack = nrbrack + 1

    if ( nlbrack < nrbrack ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'VRML_READ - Fatal error!'
      write ( *, '(a,i6)' ) '  Extraneous right bracket, line ', text_num
      write ( *, '(a)' ) trim ( text )
      write ( *, '(a)' ) '  Currently processing field:'
      write ( *, '(a)' ) trim ( level_name(level) )
      ierror = 1
      return
    end if

  end if
!
!  If the word is DEF, then read the next word right now,
!  and DON'T copy WORD into WORDM1.  
!
  if ( s_eqi ( word, 'DEF' ) ) then
    call file_get_next_word ( iunit, word, text, text_num, ierror )
    write ( *, '(a)' ) 'Skipping DEF ' // trim ( word )
    go to 11
  end if
!
!  If the word is a left bracket, then the previous word
!  is the name of a node.
!
  if ( word == '{' .or. word == '[' ) then

    level = nlbrack - nrbrack
    if ( level < 0 ) then
      write ( *, '(a)' ) 'Too many right brackets!'
      level = 0
    else if ( level_max < level ) then
      write ( *, '(a)' ) 'Too many left brackets!'
      level = level_max
    end if

    level_name(level) = wordm1

    if ( debug ) then
      write ( *, '(a)' ) ' '
      do i = 0, level
        write ( *, '(i3,2x,a)' ) i, trim ( level_name(i) )
      end do
    end if

  end if
!
!  ANCHOR
!
  if ( s_eqi ( level_name(level), 'Anchor' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  APPEARANCE
!
  else if ( s_eqi ( level_name(level), 'Appearance' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'material' ) ) then
    else if ( s_eqi ( word, 'texture' ) ) then
    else if ( s_eqi ( word, 'textureTransform' ) ) then

    end if
!
!  AUDIOCLIP
!
  else if ( s_eqi ( level_name(level), 'AudioClip' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  BACKGROUND
!
  else if ( s_eqi ( level_name(level), 'Background' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'skyColor' ) ) then
    else if ( s_eqi ( word, 'skyAngle' ) ) then
    else if ( s_eqi ( word, 'groundColor' ) ) then
    else if ( s_eqi ( word, 'groundAngle' ) ) then

    end if
!
!  BILLBOARD
!
  else if ( s_eqi ( level_name(level), 'Billboard' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  BOX
!
  else if ( s_eqi ( level_name(level), 'Box' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  COLLISION
!
  else if ( s_eqi ( level_name(level), 'Collision' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  COLOR
!
  else if ( level_name(level) == 'Color' ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( word == 'color' ) then

      write ( *, '(a)' ) 'DEBUG: COLOR saw color.'

    end if
!
!  COLOR { COLOR [] }
!
  else if ( level_name(level) == 'color' ) then

    if ( level_name(level-1) == 'Color' ) then

      if ( word == '[' ) then

        r4vec_num = 0
        new_color = 0

      else if ( word == ']' ) then

        if ( r4vec_num == 3 ) then

          if ( material_num <= material_max ) then
            new_color = new_color + 1
            material_rgba(1,material_num) = rvec(1)
            material_rgba(2,material_num) = rvec(2)
            material_rgba(3,material_num) = rvec(3)
            material_rgba(4,material_num) = 1.0E+00
            call i4_to_s_zero ( material_num, char4 )
            material_name(material_num) = 'Material_' // char4
          end if

        end if

        r4vec_num = 0
        level = nlbrack - nrbrack

      else

        call s_to_r4 ( word, rval, ierror, lchar )
        r4vec_num = r4vec_num + 1
        rvec(r4vec_num) = rval

        if ( r4vec_num == 3 ) then

          material_num = material_num + 1

          if ( material_num <= material_max ) then
            new_color = new_color + 1
            material_rgba(1,material_num) = rvec(1)
            material_rgba(2,material_num) = rvec(2)
            material_rgba(3,material_num) = rvec(3)
            material_rgba(4,material_num) = 1.0E+00
            call i4_to_s_zero ( material_num, char4 )
            material_name(material_num) = 'Material_' // char4
          end if

          r4vec_num = 0

        end if

      end if

    end if
!
!  COLORINDEX
!
  else if ( s_eqi ( level_name(level), 'colorIndex' ) ) then

    if ( word == '[' ) then

      ivec_num = 0

    else if ( word == ']' ) then

      write ( *, '(a,i6)' ) 'Hey, IVEC_NUM is ', ivec_num
      new_color_index = ivec_num

!         if ( ivec_num /= 0 ) then

!           face_num = face_num + 1

!           if ( face_num <= face_max ) then

!             if ( order_max < ivec_num ) then
!               ivec_num = order_max
!             end if

!             do i = 1, ivec_num
!               face(i,face_num) = ivec(i)
!             end do
!             face_order(face_num) = ivec_num
!             ivec_num = 0

!           end if

!         end if

      level = nlbrack - nrbrack

    else

      call s_to_i4 ( word, ival, ierror, lchar )

      if ( ivec_num < ivec_max ) then
        ivec_num = ivec_num + 1
        ivec(ivec_num) = ival + OFFSET
      end if

    end if
!
!  COLORINTERPOLATOR
!
  else if ( s_eqi ( level_name(level), 'ColorInterpolator' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  CONE
!
  else if ( s_eqi ( level_name(level), 'CONE' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  COORDINATE
!
  else if ( s_eqi ( level_name(level), 'Coordinate' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'point' ) ) then

    end if
!
!  COORDINATE POINT
!
  else if ( s_eqi ( level_name(level), 'point' ) ) then

    if ( s_eqi ( level_name(level-1), 'Coordinate' ) ) then

      if ( word == '[' ) then

        r4vec_num = 0

      else if ( word == ']' ) then

        if ( r4vec_num == 3 ) then

          cor3_num = cor3_num + 1

          if ( cor3_num <= cor3_max ) then

            call tmat_mxv ( transform_matrix, rvec, rvec )

            cor3(1:3,cor3_num) = rvec(1:3)

          end if

        end if

        r4vec_num = 0
        level = nlbrack - nrbrack

      else

        call s_to_r4 ( word, rval, ierror, lchar )
        r4vec_num = r4vec_num + 1
        rvec(r4vec_num) = rval

        if ( r4vec_num == 3 ) then

          cor3_num = cor3_num + 1

          if ( cor3_num <= cor3_max ) then
            call tmat_mxv ( transform_matrix, rvec, rvec )
            cor3(1:3,cor3_num) = rvec(1:3)
          end if

          r4vec_num = 0

        end if

      end if

    end if
!
!  COORDINATEINTERPOLATOR
!
  else if ( s_eqi ( level_name(level), 'CoordinateInterpolator' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  COORDINDEX
!
  else if ( s_eqi ( level_name(level), 'coordIndex' ) ) then

    if ( word == '[' ) then

      ivec_num = 0

    else if ( word == ']' ) then

      if ( ivec_num /= 0 ) then

        face_num = face_num + 1

        if ( face_num <= face_max ) then

          if ( order_max < ivec_num ) then
            ivec_num = order_max
          end if

          do i = 1, ivec_num
            face(i,face_num) = ivec(i)
          end do
          face_order(face_num) = ivec_num
          ivec_num = 0

        end if

      end if

      level = nlbrack - nrbrack

    else

      call s_to_i4 ( word, ival, ierror, lchar )

      if ( ival /= -1 ) then

        if ( ivec_num < ivec_max ) then
          ivec_num = ivec_num + 1
          ivec(ivec_num) = ival + cor3_num_old + OFFSET
        end if

      else

        face_num = face_num + 1

        if ( face_num <= face_max ) then

          if ( order_max < ivec_num ) then
            ivec_num = order_max
          end if

          do i = 1, ivec_num
            face(i,face_num) = ivec(i)
          end do
          face_order(face_num) = ivec_num

          ivec_num = 0

        end if
      end if

    end if
!
!  CYLINDER
!
  else if ( s_eqi ( level_name(level), 'Cylinder' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  CYLINDERSENSOR
!
  else if ( s_eqi ( level_name(level), 'CylinderSensor' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  DIRECTIONALLIGHT
!
  else if ( s_eqi ( level_name(level), 'DirectionalLight' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  ELEVATIONGRID
!
  else if ( s_eqi ( level_name(level), 'ElevationGrid' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  EXTRUSION
!
  else if ( s_eqi ( level_name(level), 'Extrusion' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  FOG
!
  else if ( s_eqi ( level_name(level), 'Fog' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  FONTSTYLE
!
  else if ( s_eqi ( level_name(level), 'FontStyle' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  GROUP
!
  else if ( s_eqi ( level_name(level), 'Group' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  IMAGETEXTURE
!
  else if ( s_eqi ( level_name(level), 'ImageTexture' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  INDEXEDFACESET
!
  else if ( s_eqi ( level_name(level), 'IndexedFaceSet' ) ) then

    if ( word == '{' ) then

      material_binding = 'PerVertex'
      cor3_num_old = cor3_num
      face_num_old = face_num
      new_color = 0
      new_color_index = 0

      if ( material_num == 0 ) then

        material_num = material_num + 1

        material_rgba(1,material_num) = 1.0E+00
        material_rgba(2,material_num) = 0.0E+00
        material_rgba(3,material_num) = 0.0E+00
        material_rgba(4,material_num) = 1.0E+00
        call i4_to_s_zero ( material_num, char4 )
        material_name(material_num) = 'Material_' // char4

      end if

      overall_mat = material_num

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

      imat = overall_mat

      write ( *, '(a,i6)' ) 'New_Color = ', new_color
      write ( *, '(a,i6)' ) 'New_Color_Index = ', new_color_index
      write ( *, '(a)' ) 'Material binding is ' // material_binding(1:9)

      if ( material_binding == 'PerVertex' ) then

        cor3_new = min ( cor3_num, cor3_max ) - min ( cor3_num_old, cor3_max )

        do i = 1, cor3_new
          if ( new_color /= 0 ) then
            if ( new_color_index == 0 ) then
              j = mod ( i-1, new_color )
            else
              k = mod ( i-1, new_color_index ) + 1
              j = ivec(k)
            end if
            imat = overall_mat + j + 1
          end if
          icor3 = cor3_num_old + i
          cor3_material(icor3) = imat
        end do

        new_face = min ( face_num, face_max ) - min ( face_num_old, face_max )

        do i = 1, new_face
          iface = face_num_old + i
          do ivert = 1, face_order(iface)
            node = face(ivert,iface)
            vertex_material(ivert,iface) = cor3_material(node)
          end do
        end do

        ivert = 1
        do i = 1, new_face
          iface = face_num_old + i
          face_material(iface) = vertex_material(ivert,iface)
        end do

      else if ( material_binding == 'PerFace' ) then

        new_face = min ( face_num, face_max ) - min ( face_num_old, face_max )

        do i = 1, new_face

          if ( new_color /= 0 ) then
            if ( new_color_index == 0 ) then
              j = mod ( i-1, new_color )
            else
              k = mod ( i-1, new_color_index ) + 1
              j = ivec(k)
            end if
            imat = overall_mat + j + 1
          end if

          iface = face_num_old + i
          face_material(iface) = imat
        end do 

        do i = 1, new_face
          iface = face_num_old + i
          do ivert = 1, face_order(iface)
            vertex_material(ivert,iface) = face_material(iface)
          end do
        end do

        do i = 1, new_face
          iface = face_num_old + i
          do ivert = 1, face_order(iface)
            node = face(ivert,iface)
            cor3_material(node) = vertex_material(ivert,iface)
          end do
        end do

      else

        write ( *, '(a)' ) 'Cannot decide what material binding is...'

      end if

      cor3_num_old = cor3_num
      face_num_old = face_num

    else if ( s_eqi ( word, 'ccw' ) ) then
    else if ( s_eqi ( word, 'color' ) ) then
    else if ( s_eqi ( word, 'colorIndex' ) ) then
    else if ( s_eqi ( word, 'colorPerVertex' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )

      if ( s_eqi ( word, 'TRUE' ) ) then
        material_binding = 'PerVertex'
      else if ( s_eqi ( word, 'FALSE' ) ) then
        material_binding = 'PerFace'
      end if

    else if ( s_eqi ( word, 'convex' ) ) then
    else if ( s_eqi ( word, 'coord' ) ) then
    else if ( s_eqi ( word, 'coordIndex' ) ) then
    else if ( s_eqi ( word, 'creaseAngle' ) ) then
    else if ( s_eqi ( word, 'normal' ) ) then
    else if ( s_eqi ( word, 'normalIndex' ) ) then
    else if ( s_eqi ( word, 'normalPerVertex' ) ) then
    else if ( s_eqi ( word, 'solid' ) ) then
    else if ( s_eqi ( word, 'texCoord' ) ) then
    else if ( s_eqi ( word, 'texCoordIndex' ) ) then

    end if
!
!  INDEXEDLINESET
!
  else if ( s_eqi ( level_name(level), 'IndexedLineSet' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  INLINE
!
  else if ( s_eqi ( level_name(level), 'Inline' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  LOD
!
  else if ( s_eqi ( level_name(level), 'LOD' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  MATERIAL
!
  else if ( s_eqi ( level_name(level), 'Material' ) ) then

    if ( word == '{' ) then

      r01 = 0.2

      r02 = 0.8
      r03 = 0.8
      r04 = 0.8

      r05 = 0.0E+00
      r06 = 0.0E+00
      r07 = 0.0E+00

      r08 = 0.2

      r09 = 0.0E+00
      r10 = 0.0E+00
      r11 = 0.0E+00

      r12 = 0.0E+00

    else if ( word == '}' ) then

      material_num = material_num + 1

      if ( material_num <= material_max ) then
        material_rgba(1,material_num) = r02
        material_rgba(2,material_num) = r03
        material_rgba(3,material_num) = r04
        material_rgba(4,material_num) = 1.0E+00 - r12
        call i4_to_s_zero ( material_num, char4 )
        material_name(material_num) = 'Material_' // char4
      end if

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'ambientIntensity' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, r01, ierror, lchar )

    else if ( s_eqi ( word, 'diffuseColor' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, r02, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, r03, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, r04, ierror, lchar )

    else if ( s_eqi ( word, 'emissiveColor' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, r05, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, r06, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, r07, ierror, lchar )

    else if ( s_eqi ( word, 'shininess' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, r08, ierror, lchar )

    else if ( s_eqi ( word, 'specularColor' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, r09, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, r10, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, r11, ierror, lchar )

    else if ( s_eqi ( word, 'transparency' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, r12, ierror, lchar )

    end if
!
!  MOVIETEXTURE
!
  else if ( s_eqi ( level_name(level), 'MovieTexture' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  NAVIGATIONINFO
!
  else if ( s_eqi ( level_name(level), 'NavigationInfo' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else

    end if
!
!  NORMAL
!
  else if ( s_eqi ( level_name(level), 'Normal' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  NORMALINTERPOLATOR
!
  else if ( s_eqi ( level_name(level), 'NormalInterpolator' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  ORIENTATIONINTERPOLATOR
!
  else if ( s_eqi ( level_name(level), 'OrientationInterpolator' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  PIXELTEXTURE
!
  else if ( s_eqi ( level_name(level), 'PixelTexture' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  PLANESENSOR
!
  else if ( s_eqi ( level_name(level), 'PlaneSensor' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  POINTLIGHT
!
  else if ( s_eqi ( level_name(level), 'PointLight' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  POINTSET
!
  else if ( s_eqi ( level_name(level), 'PointSet' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  POSITIONINTERPOLATOR
!
  else if ( s_eqi ( level_name(level), 'PositionInterpolator' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  PROXIMITYSENSOR
!
  else if ( s_eqi ( level_name(level), 'ProximitySensor' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  SCALARINTERPOLATOR
!
  else if ( s_eqi ( level_name(level), 'ScalarInterpolator' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  SCRIPT
!
  else if ( s_eqi ( level_name(level), 'Script' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  SHAPE
!
  else if ( s_eqi ( level_name(level), 'Shape' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'appearance' ) ) then
    else if ( s_eqi ( word, 'geometry' ) ) then

    end if
!
!  SOUND
!
  else if ( s_eqi ( level_name(level), 'Sound' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  SPHERE
!
  else if ( s_eqi ( level_name(level), 'Sphere' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  SPHERESENSOR
!
  else if ( s_eqi ( level_name(level), 'SphereSensor' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  SPOTLIGHT
!
  else if ( s_eqi ( level_name(level), 'SpotLight' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  SWITCH
!
  else if ( s_eqi ( level_name(level), 'Switch' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  TEXT
!
  else if ( s_eqi ( level_name(level), 'Text' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  TEXTURECOORDINATE
!
  else if ( s_eqi ( level_name(level), 'TextureCoordinate' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  TEXTURETRANSFORM
!
  else if ( s_eqi ( level_name(level), 'TextureTransform' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  TIMESENSOR
!
  else if ( s_eqi ( level_name(level), 'TimeSensor' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  TOP
!
  else if ( s_eqi ( level_name(level), 'Top' ) ) then

    if ( word == 'Anchor' ) then
    else if ( word == 'Appearance' ) then
    else if ( word == 'AudioClip' ) then
    else if ( word == 'Background' ) then
    else if ( word == 'Billboard' ) then
    else if ( word == 'Box' ) then
    else if ( word == 'Collision' ) then
    else if ( word == 'Color' ) then
    else if ( word == 'ColorInterpolator' ) then
    else if ( word == 'Cone' ) then
    else if ( word == 'CoordinateInterpolator' ) then
    else if ( word == 'Cylinder' ) then
    else if ( word == 'CylinderSensor' ) then

    else if ( word == 'DirectionalLight' ) then
    else if ( word == 'ElevationGrid' ) then
    else if ( word == 'Extrusion' ) then
    else if ( word == 'Fog' ) then
    else if ( word == 'Fontstyle' ) then
    else if ( word == 'Group' ) then
    else if ( word == 'ImageTexture' ) then
    else if ( word == 'IndexedFaceSet' ) then
    else if ( word == 'IndexedLineSet' ) then
    else if ( word == 'Inline' ) then
    else if ( word == 'LOD' ) then
    else if ( word == 'Material' ) then
    else if ( word == 'MovieTexture' ) then
    else if ( word == 'NavigationInfo' ) then
    else if ( word == 'Normal' ) then
    else if ( word == 'NormalInterpolator' ) then
    else if ( word == 'OrientationInterpolator' ) then
    else if ( word == 'PixelTexture' ) then
    else if ( word == 'PlaneSensor' ) then
    else if ( word == 'PointLight' ) then
    else if ( word == 'PointSet' ) then
    else if ( word == 'PositionInterpolator' ) then
    else if ( word == 'ProximitySensor' ) then
    else if ( word == 'ScalarInterpolator' ) then
    else if ( word == 'Script' ) then
    else if ( word == 'Sound' ) then
    else if ( word == 'Sphere' ) then
    else if ( word == 'SphereSensor' ) then
    else if ( word == 'SpotLight' ) then
    else if ( word == 'Switch' ) then
    else if ( word == 'Text' ) then
    else if ( word == 'TextureCoordinate' ) then
    else if ( word == 'TextureTransform' ) then
    else if ( word == 'TimeSensor' ) then
    else if ( word == 'TouchSensor' ) then
    else if ( word == 'Transform' ) then
    else if ( word == 'Viewpoint' ) then
    else if ( word == 'VisibilitySensor' ) then
    else if ( word == 'WorldInfo' ) then
    end if
!
!  TOUCHSENSOR
!
  else if ( s_eqi ( level_name(level), 'TouchSensor' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  TRANSFORM
!
  else if ( s_eqi ( level_name(level), 'Transform' ) ) then

    if ( word == '{' ) then

      angle = 0.0E+00

      rx = 1.0E+00
      ry = 1.0E+00
      rz = 1.0E+00

      sx = 1.0E+00
      sy = 1.0E+00
      sz = 1.0E+00

      tx = 1.0E+00
      ty = 1.0E+00
      tz = 1.0E+00

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    else if ( s_eqi ( word, 'rotation' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, rx, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, ry, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, rz, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, angle, ierror, lchar )

    else if ( s_eqi ( word, 'scale' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, sx, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, sy, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, sz, ierror, lchar )

    else if ( s_eqi ( word, 'translation' ) ) then

      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, tx, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )
      call s_to_r4 ( word, ty, ierror, lchar )
      call file_get_next_word ( iunit, word, text, text_num, ierror )

      call s_to_r4 ( word, tz, ierror, lchar )

    else if ( s_eqi ( word, 'children' ) ) then

      call tmat_init ( transform_matrix )

      call tmat_scale ( transform_matrix, transform_matrix, sx, sy, sz )

      axis(1) = rx
      axis(2) = ry
      axis(3) = rz

      call tmat_rot_vector ( transform_matrix, transform_matrix, angle, axis )

      call tmat_trans ( transform_matrix, transform_matrix, tx, ty, tz )

    end if
!
!  VIEWPOINT
!
  else if ( s_eqi ( level_name(level), 'Viewpoint' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  VISIBILITYSENSOR
!
  else if ( s_eqi ( level_name(level), 'VisibilitySensor' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  WORLDINFO
!
  else if ( s_eqi ( level_name(level), 'WorldInfo' ) ) then

    if ( word == '{' ) then

    else if ( word == '}' ) then

      level = nlbrack - nrbrack

    end if
!
!  Any other word:
!
  else

  end if

  go to 10
!
!  Bad data
!
99    continue

  bad_num = bad_num + 1

  if ( bad_num <= 10 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VRML_READ - Warning!'
    write ( *, '(a)' ) '  Bad data on level ' // trim ( level_name(level) )
    write ( *, '(a)' ) '  Bad word: ' // trim ( word )
    write ( *, '(a,i6)' ) '  Line number: ', text_num
    write ( *, '(a)' ) trim ( text )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VRML_READ - Fatal error!'
    write ( *, '(a)' ) '  Too many warnings!'
    return
  end if

  go to 10
!
!  End of information in file.
!
50    continue

  ierror = 0
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'VRML_READ - Read ', text_num, &
    ' text lines from ' // trim ( filein_name )

  return
end
subroutine vrml_write ( cor3, cor3_max, cor3_num, face, face_max, face_num, &
  face_order, filein_name, fileout_name, iunit, line_dex, line_material, &
  line_max, line_num, material_max, material_num, material_rgba, order_max, &
  vertex_material )

!*****************************************************************************80
!
!! VRML_WRITE writes graphics data to a VRML file.
!
!  The VRML files written by this routine have the form:
!
!
!     #VRML V2.0 utf8
!
!       WorldInfo {
!         title "cube.iv."
!         string "VRML file generated by IVREAD.
!       }
!
!       Group {
!         children [
!
!           Shape {
!
!             appearance Appearance {
!               material Material {
!                 diffuseColor   0.0 0.0 0.0E+00
!                 emissiveColor  0.0 0.0 0.0E+00
!                 shininess      1.0E+00
!               }
!             } #end of appearance
!
!             geometry IndexedLineSet {
!
!               coord Coordinate {
!                 point [
!                   8.59816       5.55317      -3.05561
!                   8.59816       2.49756      0.000000E+00
!                   ...etc...
!                   2.48695       2.49756      -3.05561
!                 ]
!               }
!
!               coordIndex [
!                   0    1     2    -1     3     4     5     6     7     8    -1
!                   9   10    -1    11    12    -1    13    14    15    -1    16
!                 ...etc...
!                 191    -1
!               ]
!
!               colorPerVertex TRUE
!
!               colorIndex [
!                   0    0     0    -1     2     3     1     1     4     7    -1
!                  10    9    -1     7     7    -1     3     2     2    -1    12
!                 ...etc...
!                 180    -1
!               ]
!
!             }  #end of geometry
!
!           }  #end of Shape
!
!         ]  #end of children
!       }  #end of Group
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer ( kind = 4 ) LINE_MATERIAL(LINE_MAX), material index for line.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) line_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) filein_name
  character ( len = * ) fileout_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ) length
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_material(line_max)
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ) ndx
  integer ( kind = 4 ), parameter :: OFFSET = 1
  character ( len = 100 ) text
  integer ( kind = 4 ) text_num
  integer ( kind = 4 ) vertex_material(order_max,face_max)
  character ( len = 10 ) word
!
  write ( iunit, '(a)' ) '#VRML V2.0 utf8'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) '  WorldInfo {'
  write ( iunit, '(a)' ) '    title "' // trim ( fileout_name ) //'"'
  write ( iunit, '(a)' ) '    info "VRML file generated by IVREAD."'
  write ( iunit, '(a)' ) '    info "Original data in file ' // &
    trim ( filein_name ) // '."'
  write ( iunit, '(a)' ) '  }'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) '  Group {'
  write ( iunit, '(a)' ) '    children ['
  write ( iunit, '(a)' ) '      Shape {'
  write ( iunit, '(a)' ) '        appearance Appearance {'
  write ( iunit, '(a)' ) '          material Material {'
  write ( iunit, '(a)' ) '            diffuseColor   0.0 0.0 0.0'
  write ( iunit, '(a)' ) '            emissiveColor  0.0 0.0 0.0'
  write ( iunit, '(a)' ) '            shininess      1.0'
  write ( iunit, '(a)' ) '          }'
  write ( iunit, '(a)' ) '        } '

  text_num = 18
!
!  IndexedLineSet
!
  if ( 0 < line_num ) then

    write ( iunit, '(a)' ) '        geometry IndexedLineSet {'
!
!  IndexedLineSet coord
!
    write ( iunit, '(a)' ) '          coord Coordinate {'
    write ( iunit, '(a)' ) '            point ['
 
    text_num = text_num + 3
 
    do icor3 = 1, cor3_num
      write ( text, '(3f12.4,'','')' ) cor3(1:3,icor3)
      call s_blanks_delete ( text )
      write ( iunit, '(8x,a)' ) trim ( text )
      text_num = text_num + 1
    end do

    write ( iunit, '(a)' ) '            ]'
    write ( iunit, '(a)' ) '          }' 
    text_num = text_num + 2
!
!  IndexedLineSet coordIndex.
!
    write ( iunit, '(a)' ) '          coordIndex ['

    text_num = text_num + 1

    text = ' '
    length = 0
    
    do j = 1, line_num
 
      write ( word, '(i8,'','')' ) line_dex(j) - OFFSET
      call s_cat ( text, word, text )
      length = length + 1

      if ( line_dex(j)-OFFSET == -1 .or. &
           10 <= length .or. &
           line_num <= j ) then

        call s_blanks_delete ( text )
        write ( iunit, '(8x,a)' ) trim ( text )
        text_num = text_num + 1
        text = ' '
        length = 0

      end if

    end do
 
    write ( iunit, '(a)' ) '          ]'
    text_num = text_num + 1
!
!  Colors. (materials)
!
    write ( iunit, '(a)' ) '          color Color {'
    write ( iunit, '(a)' ) '            color ['
    text_num = text_num + 2

    do j = 1, material_num
      write ( iunit, '(3f12.4,'','')' ) material_rgba(1:3,j)
      text_num = text_num + 1
    end do

    write ( iunit, '(a)' ) '            ]'
    write ( iunit, '(a)' ) '          }'
    write ( iunit, '(a)' ) '          colorPerVertex TRUE'
!
!  IndexedLineset colorIndex
!
    write ( iunit, '(a)' ) '          colorIndex ['

    text_num = text_num + 4

    text = ' '
    length = 0
    
    do j = 1, line_num
 
      write ( word, '(i8,'','')' ) line_material(j) - OFFSET
      call s_cat ( text, word, text )
      length = length + 1

      if ( line_dex(j)-OFFSET == -1 .or. length >= 10 .or. j >= line_num ) then

        call s_blanks_delete ( text )
        write ( iunit, '(8x,a)' ) trim ( text )
        text_num = text_num + 1
        text = ' '
        length = 0

      end if

    end do

    if ( text /= ' ' ) then
      write ( iunit, '(a)' ) trim ( text )
      text_num = text_num + 1
      text = ' '
    end if

    write ( iunit, '(a)' ) '          ]'
    write ( iunit, '(a)' ) '        }'
    text_num = text_num + 2

  end if
!
!  End of IndexedLineSet
!
!
!  IndexedFaceSet
!
  if ( 0 < face_num ) then

    write ( iunit, '(a)' ) '        geometry IndexedFaceSet {'
!
!  IndexedFaceSet coord
!
    write ( iunit, '(a)' ) '          coord Coordinate {'
    write ( iunit, '(a)' ) '            point ['
 
    text_num = text_num + 3
 
    do icor3 = 1, cor3_num
      write ( text, '(3f12.4,'','')' ) cor3(1:3,icor3)
      call s_blanks_delete ( text )
      write ( iunit, '(8x,a)' ) trim ( text )
      text_num = text_num + 1
    end do

    write ( iunit, '(a)' ) '            ]'
    write ( iunit, '(a)' ) '          }' 
!
!  IndexedFaceSet coordIndex.
!
    write ( iunit, '(a)' ) '          coordIndex ['

    text_num = text_num + 3

    text = ' '
    length = 0
 
    do iface = 1, face_num

      do ivert = 1, face_order(iface) + 1

        if ( ivert <= face_order(iface) ) then
          itemp = face(ivert,iface) - OFFSET
        else
          itemp = 0 - OFFSET
        end if

        write ( word, '(i8,'','')' ) itemp

        call s_cat ( text, word, text )
        length = length + 1

        if ( itemp == -1 .or. length >= 10 .or. &
           ( iface == face_num .and. ivert == face_order(iface) + 1 ) ) then

          call s_blanks_delete ( text )
          write ( iunit, '(8x,a)' ) trim ( text )
          text_num = text_num + 1
          text = ' '
          length = 0

        end if
 
      end do

    end do
 
    write ( iunit, '(a)' ) '          ]'
    text_num = text_num + 1
!
!  IndexedFaceset colorIndex
!
    write ( iunit, '(a)' ) '          colorIndex ['

    text_num = text_num + 4

    text = ' '
    length = 0
    ndx = 0
 
    do iface = 1, face_num
 
      do ivert = 1, face_order(iface) + 1

        if ( ivert <= face_order(iface) ) then
          itemp = vertex_material(ivert,iface) - OFFSET
          ndx = ndx + 1
        else
          itemp = 0 - OFFSET
        end if

        write ( word, '(i8,'','')' ) itemp

        call s_cat ( text, word, text )
        length = length + 1

        if ( itemp == -1 .or.  length >= 10 .or. &
           ( iface == face_num .and.  ivert == face_order(iface) + 1 )  ) then

          call s_blanks_delete ( text )
          write ( iunit, '(8x,a)' ) trim ( text )
          text_num = text_num + 1
          text = ' '
          length = 0
  
        end if

      end do
 
    end do

    if ( text /= ' ' ) then
      write ( iunit, '(a)' ) trim ( text )
      text_num = text_num + 1
      text = ' '
    end if

    write ( iunit, '(a)' ) '          ]'
    write ( iunit, '(a)' ) '        }'
    text_num = text_num + 2

  end if
!
!  End of IndexedFaceSet
!
!  End of:
!        Shape
!      children
!    Group
!
  write ( iunit, '(a)' ) '      }'
  write ( iunit, '(a)' ) '    ]'
  write ( iunit, '(a)' ) '  }'
 
  text_num = text_num + 3
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'VRML_WRITE - Wrote ', text_num, &
    ' text lines to ' // trim ( fileout_name )

  return
end
subroutine word_next_read ( s, word, done )

!*****************************************************************************80
!
!! WORD_NEXT_READ "reads" words from a string, one at a time.
!
!  Discussion:
!
!    The following characters are considered to be a single word,
!    whether surrounded by spaces or not:
!
!      " ( ) { } [ ]
!
!    Also, if there is a trailing comma on the word, it is stripped off.
!    This is to facilitate the reading of lists.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string, presumably containing words
!    separated by spaces.
!
!    Output, character ( len = * ) WORD.
!    If DONE is FALSE, then WORD contains the "next" word read.
!    If DONE is TRUE, then WORD is blank, because there was no more to read.
!
!    Input/output, logical DONE.
!    On input with a fresh string, set DONE to TRUE.
!    On output, the routine sets DONE:
!      FALSE if another word was read,
!      TRUE if no more words could be read.
!
  implicit none

  logical done
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ), save :: lenc = 0
  integer ( kind = 4 ), save :: next = 1
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
  character ( len = * ) word
!
!  We "remember" LENC and NEXT from the previous call.
!
!  An input value of DONE = TRUE signals a new line of text to examine.
!
  if ( done ) then

    next = 1
    done = .false.
    lenc = len_trim ( s )

    if ( lenc <= 0 ) then
      done = .true.
      word = ' '
      return
    end if

  end if
!
!  Beginning at index NEXT, search the string for the next nonblank,
!  which signals the beginning of a word.
!
  ilo = next
!
!  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
!
  do

    if ( lenc < ilo ) then
      word = ' '
      done = .true.
      next = lenc + 1
      return
    end if
!
!  If the current character is blank, skip to the next one.
!
    if ( s(ilo:ilo) /= ' ' .and. s(ilo:ilo) /= TAB ) then
      exit
    end if

    ilo = ilo + 1

  end do
!
!  ILO is the index of the next nonblank character in the string.
!
!  If this initial nonblank is a special character,
!  then that's the whole word as far as we're concerned,
!  so return immediately.
!
  if ( s(ilo:ilo) == '"' .or. &
       s(ilo:ilo) == '(' .or. &
       s(ilo:ilo) == ')' .or. &
       s(ilo:ilo) == '{' .or. &
       s(ilo:ilo) == '}' .or. &
       s(ilo:ilo) == '[' .or. &
       s(ilo:ilo) == ']' ) then

    word = s(ilo:ilo)
    next = ilo + 1
    return

  end if
!
!  Now search for the last contiguous character that is not a
!  blank, TAB, or special character.
!
  next = ilo + 1

  do while ( next <= lenc )

    if ( s(next:next) == ' ' ) then
      exit
    else if ( s(next:next) == TAB ) then
      exit
    else if ( s(next:next) == '"' ) then
      exit
    else if ( s(next:next) == '(' ) then
      exit
    else if ( s(next:next) == ')' ) then
      exit
    else if ( s(next:next) == '{' ) then
      exit
    else if ( s(next:next) == '}' ) then
      exit
    else if ( s(next:next) == '[' ) then
      exit
    else if ( s(next:next) == ']' ) then
      exit
    end if

    next = next + 1

  end do

  if ( s(next-1:next-1) == ',' ) then
    word = s(ilo:next-2)
  else
    word = s(ilo:next-1)
  end if

  return
end
subroutine xgl_write ( cor3, cor3_max, cor3_normal, cor3_num, face, &
  face_material, face_max, face_num, face_order, fileout_name, iunit, &
  material_max, material_num, material_rgba, order_max )

!*****************************************************************************80
!
!! XGL_WRITE writes graphics data to an XGL file.
!
!  Discussion:
!
!    Only triangular faces are allowed.
!
!  Example:
!
!    <WORLD>
!
!      <BACKGROUND>
!        <BACKCOLOR>r,g,b</BACKCOLOR>
!      </BACKGROUND>
! 
!      <LIGHTING>
!        <AMBIENT>r,g,b</AMBIENT>
!        <DIRECTIONALLIGHT>
!          <DIFFUSE>r,g,b</DIFFUSE>
!          <DIRECTION>x,y,z</DIRECTION>
!          <SPECULAR>r,g,b</SPECULAR>
!        </DIRECTIONALLIGHT>
!      </LIGHTING>
!
!      <MESH ID="0">
!
!        <P ID = "0"> x, y, z </P>
!        ...
!        <P ID = "99"> x, y, z </P>
!
!        <N ID = "0"> nx, ny, nz </N>
!        ...
!        <N ID = "55"> nx, ny, nz </N>
!
!        <MAT ID = "0">
!          <ALPHA></ALPHA>
!          <AMB>r,g,b</AMB>
!          <DIFF>r,g,b</DIFF>
!          <EMISS>r,g,b</EMISS>
!          <SHINE></SHINE>
!          <SPEC>r,g,b</SPEC>
!        </MAT>
!
!        <F>
!          <MATREF>0</MATREF>
!          <FV1><PREF>0</PREF><NREF>0</NREF></FV1>
!          <FV2><PREF>1</PREF><NREF>1</NREF></FV2>
!          <FV3><PREF>2</PREF><NREF>2</NREF></FV3>
!        </F>
!        ...
!        <F>
!          <MATREF>0</MATREF>
!          <FV1><PREF>12</PREF><NREF>12</NREF></FV1>
!          <FV2><PREF>13</PREF><NREF>13</NREF></FV2>
!          <FV3><PREF>14</PREF><NREF>14</NREF></FV3>
!        </F>
!
!       </MESH>
!
!       <OBJECT>
!         <TRANSFORM>
!           <FORWARD>x,y,z</FORWARD>
!           <POSITION>x,y,z</POSITION>
!           <SCALE>x,y,z</SCALE>
!           <UP>x,y,z</UP>
!         </TRANSFORM>
!         <MESHREF>0</MESHREF>
!       </OBJECT>
!
!       <OBJECT>
!         <TRANSFORM>
!           <FORWARD>x,y,z</FORWARD>
!           <POSITION>x,y,z</POSITION>
!           <SCALE>x,y,z</SCALE>
!           <UP>x,y,z</UP>
!         </TRANSFORM>
!         <MESHREF>0</MESHREF>
!       </OBJECT>
!
!     </WORLD>
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, real ( kind = 4 ) COR3_NORMAL(3,COR3_MAX), normals at nodes.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE(ORDER_MAX,FACE_MAX), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_MAX, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(FACE_MAX), the number of vertices per face.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer ( kind = 4 ) MATERIAL_MAX, the maximum number of materials.
!
!    Input, integer ( kind = 4 ) MATERIAL_NUM, the number of materials.
!
!    Input, real ( kind = 4 ) MATERIAL_RGBA(4,MATERIAL_MAX), material R, G, B and A values.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices per face.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) face_max
  integer ( kind = 4 ) material_max
  integer ( kind = 4 ) order_max

  real ( kind = 4 ) background_rgb(3)
  real ( kind = 4 ) cor3(3,cor3_max)
  real ( kind = 4 ) cor3_normal(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  integer ( kind = 4 ) face(order_max,face_max)
  integer ( kind = 4 ) face_material(face_max)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) face_order(face_max)
  character ( len = * ) fileout_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  real ( kind = 4 ) light_ambient_rgb(3)
  real ( kind = 4 ) light_diffuse_rgb(3)
  real ( kind = 4 ) light_direction(3)
  real ( kind = 4 ) light_specular_rgb(3)
  real ( kind = 4 ) material_alpha(1)
  real ( kind = 4 ) material_amb_rgb(3)
  real ( kind = 4 ) material_diff_rgb(3)
  real ( kind = 4 ) material_emiss_rgb(3)
  real ( kind = 4 ) material_shine(1)
  real ( kind = 4 ) material_spec_rgb(3)
  integer ( kind = 4 ) material
  integer ( kind = 4 ) material_num
  real ( kind = 4 ) material_rgba(4,material_max)
  integer ( kind = 4 ) mesh
  integer ( kind = 4 ), parameter :: mesh_num = 1
  integer ( kind = 4 ), parameter :: object_num = 1
  integer ( kind = 4 ), parameter :: OFFSET = 1
  character ( len = 100 ) s
  character s1
  character ( len = 6 ) s2
  integer ( kind = 4 ) text_num
  real ( kind = 4 ) transform_forward(3)
  real ( kind = 4 ) transform_position(3)
  real ( kind = 4 ) transform_scale(3)
  real ( kind = 4 ) transform_up(3)
!
  text_num = 0
!
!  Define some dummy values.
!
  background_rgb(1) = 0.1E+00
  background_rgb(2) = 0.1E+00
  background_rgb(3) = 0.1E+00

  light_ambient_rgb(1) = 0.2E+00
  light_ambient_rgb(2) = 0.1E+00
  light_ambient_rgb(3) = 0.1E+00

  light_diffuse_rgb(1) = 0.1E+00
  light_diffuse_rgb(2) = 0.2E+00
  light_diffuse_rgb(3) = 0.1E+00

  light_direction(1) =   0.0E+00
  light_direction(2) =   0.0E+00
  light_direction(3) = 100.0E+00

  light_specular_rgb(1) = 0.1E+00
  light_specular_rgb(2) = 0.1E+00
  light_specular_rgb(3) = 0.2E+00

  material_alpha = 0.9E+00

  material_amb_rgb(1) = 0.1E+00
  material_amb_rgb(2) = 0.1E+00
  material_amb_rgb(3) = 0.1E+00

  material_diff_rgb(1) = 0.2E+00
  material_diff_rgb(2) = 0.1E+00
  material_diff_rgb(3) = 0.1E+00

  material_emiss_rgb(1) = 0.1E+00
  material_emiss_rgb(2) = 0.2E+00
  material_emiss_rgb(3) = 0.1E+00

  material_shine = 0.8E+00

  material_spec_rgb(1) = 0.1E+00
  material_spec_rgb(2) = 0.1E+00
  material_spec_rgb(3) = 0.2E+00

  transform_forward(1) = 0.0E+00
  transform_forward(2) = 0.0E+00
  transform_forward(3) = 0.0E+00

  transform_position(1) = 0.0E+00
  transform_position(2) = 0.0E+00
  transform_position(3) = 0.0E+00

  transform_scale(1) = 1.0E+00
  transform_scale(2) = 1.0E+00
  transform_scale(3) = 1.0E+00

  transform_up(1) = 1.0E+00
  transform_up(2) = 1.0E+00
  transform_up(3) = 1.0E+00

  write ( iunit, '(a)' ) '<WORLD>'
  write ( iunit, '(a)' ) ' '

  write ( iunit, '(a)' ) '  <BACKGROUND>'
  call r4vec_to_s ( 3, background_rgb, s )
  write ( iunit, '(a)' ) '    <BACKCOLOR> ' // trim ( s ) // ' </BACKCOLOR>'
  write ( iunit, '(a)' ) '  </BACKGROUND>'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) '  <LIGHTING>'
  call r4vec_to_s ( 3, light_ambient_rgb, s )
  write ( iunit, '(a)' ) '    <AMBIENT> ' // trim ( s ) // ' </AMBIENT>'
  write ( iunit, '(a)' ) '    <DIRECTIONALLIGHT>'
  call r4vec_to_s ( 3, light_diffuse_rgb, s )
  write ( iunit, '(a)' ) '      <DIFFUSE> ' // trim ( s ) // ' </DIFFUSE>'
  call r4vec_to_s ( 3, light_direction, s )
  write ( iunit, '(a)' ) '      <DIRECTION> ' // trim ( s ) // ' </DIRECTION>'
  call r4vec_to_s ( 3, light_specular_rgb, s )
  write ( iunit, '(a)' ) '      <SPECULAR> ' // trim ( s ) // ' </SPECULAR>'
  write ( iunit, '(a)' ) '    </DIRECTIONALLIGHT>'
  write ( iunit, '(a)' ) '  </LIGHTING>'

  text_num = text_num + 14

  do mesh = 1, mesh_num

    write ( iunit, '(a)' ) ' '
    write ( s2, '(i6)' ) mesh - OFFSET
    s2 = adjustl ( s2 )
    write ( iunit, '(a)' ) '  <MESH ID = "' // trim ( s2 ) // '">'

    write ( iunit, '(a)' ) ' '
    text_num = text_num + 3

    do j = 1, cor3_num
      write ( s2, '(i6)' ) j - OFFSET
      s2 = adjustl ( s2 )
      call r4vec_to_s ( 3, cor3(1,j), s )
      write ( iunit, '(a)' ) '    <P ID="' // trim ( s2 ) // '"> ' // &
        trim ( s ) // ' </P>'
      text_num = text_num + 1
    end do

    write ( iunit, '(a)' ) ' '
    text_num = text_num + 1

    do j = 1, cor3_num
      write ( s2, '(i6)' ) j - OFFSET
      s2 = adjustl ( s2 )
      call r4vec_to_s ( 3, cor3_normal(1,j), s )
      write ( iunit, '(a)' ) '    <N ID="' // trim ( s2 ) // '"> ' // &
        trim ( s ) // ' </N>'
      text_num = text_num + 1
    end do

    do material = 1, material_num
      write ( iunit, '(a)' ) ' '
      write ( s2, '(i6)' ) material - OFFSET
      s2 = adjustl ( s2 )
      write ( iunit, '(a)' ) '    <MAT ID="' // trim ( s2 ) // '">'
      call r4vec_to_s ( 1, material_alpha, s )
      write ( iunit, '(a)' ) '      <ALPHA> ' // trim ( s ) // ' </ALPHA>'
      call r4vec_to_s ( 3, material_amb_rgb, s )
      write ( iunit, '(a)' ) '      <AMB> ' // trim ( s ) // ' </AMB>'
      call r4vec_to_s ( 3, material_diff_rgb, s )
      write ( iunit, '(a)' ) '      <DIFF> ' // trim ( s ) // ' </DIFF>'
      call r4vec_to_s ( 3, material_emiss_rgb, s )
      write ( iunit, '(a)' ) '      <EMISS> ' // trim ( s ) // ' </EMISS>'
      call r4vec_to_s ( 1, material_shine, s )
      write ( iunit, '(a)' ) '      <SHINE> ' // trim ( s ) // ' </SHINE>'
      call r4vec_to_s ( 3, material_spec_rgb, s )
      write ( iunit, '(a)' ) '      <SPEC> ' // trim ( s ) // ' </SPEC>'
      write ( iunit, '(a)' ) '    </MAT>'
      text_num = text_num + 9
    end do

    write ( iunit, '(a)' ) ' '
    text_num = text_num + 1

    do iface = 1, face_num

      write ( iunit, '(a)' ) '    <F>'
      write ( s2, '(i6)' ) face_material(iface) - OFFSET
      s2 = adjustl ( s2 )
      write ( iunit, '(a)' ) '      <MATREF> ' // trim ( s2 ) // ' </MATREF>'
      text_num = text_num + 2

      do ivert = 1, face_order(iface)
        write ( s1, '(i1)' ) ivert
        write ( s2, '(i6)' ) face(ivert,iface) - OFFSET
        s2 = adjustl ( s2 )
        write ( iunit, '(a)' ) '      <FV' // trim ( s1 ) // '><PREF> ' // &
          trim ( s2 ) // ' </PREF><NREF> ' // &
          trim ( s2 ) // ' </NREF></FV' // trim ( s1 ) // '>'
        text_num = text_num + 1
      end do
      write ( iunit, '(a)' ) '    </F>'
      text_num = text_num + 1
    end do

    write ( iunit, '(a)' ) '  </MESH>'
    text_num = text_num + 1

  end do

  write ( iunit, '(a)' ) ' '
  text_num = text_num + 1

  do i = 1, object_num

    write ( iunit, '(a)' ) '  <OBJECT>'
    write ( iunit, '(a)' ) '    <TRANSFORM>'
    call r4vec_to_s ( 3, transform_forward, s )
    write ( iunit, '(a)' ) '      <FORWARD> ' // trim ( s ) // ' </FORWARD>'
    call r4vec_to_s ( 3, transform_position, s )
    write ( iunit, '(a)' ) '      <POSITION> ' // trim ( s ) // ' </POSITION>'
    call r4vec_to_s ( 3, transform_scale, s )
    write ( iunit, '(a)' ) '      <SCALE> ' // trim ( s ) // ' </SCALE>'
    call r4vec_to_s ( 3, transform_up, s )
    write ( iunit, '(a)' ) '      <UP> ' // trim ( s ) // ' </UP>'
    write ( iunit, '(a)' ) '    </TRANSFORM>'
    mesh = 1
    write ( s2, '(i6)' ) mesh - OFFSET
    s2 = adjustl ( s2 )
    write ( iunit, '(a)' ) '    <MESHREF> ' // trim ( s2 ) // ' </MESHREF>'
    write ( iunit, '(a)' ) '  </OBJECT>'
    text_num = text_num + 9

  end do

  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) '</WORLD>'
  text_num = text_num + 2
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'XGL_WRITE - Wrote ', text_num, &
    ' text lines to ' // trim ( fileout_name )

  return
end
subroutine xyz_read ( cor3, cor3_max, cor3_num, filein_name, ierror, iunit, &
  line_dex, line_max, line_num, text_num )

!*****************************************************************************80
!
!! XYZ_READ reads graphics information from an XYZ file.
!
!  Discussion:
!
!    Comment lines begin with '#";
!
!    The XYZ coordinates of a point are written on a single line;
!
!    If the next line is also a point, then the two points are to be
!    joined by a line segment.
!
!  Example:
!
!     # cube.xyz
!     #
!     #  First the points:
!
!     0 0 0
!
!     0 0 1
!
!     0 1 0
!
!     0 1 1
!
!     1 0 0
!
!     1 0 1
!
!     1 1 0
!
!     1 1 1
!
!     #
!     #  Now the lines
!     #
!     0 0 0
!     0 0 1
!     0 1 1
!     0 1 0
!     0 0 0
!
!     1 0 0
!     1 0 1
!     1 1 0
!     1 1 1
!     1 0 0
!
!     0 0 0
!     1 0 0
!
!     0 0 1
!     1 0 1
!
!     0 1 0
!     1 1 0
!
!     0 1 1
!     1 1 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input/output, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit from which data is read.
!
!    Input/output, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated
!    by -1.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Input/output, integer ( kind = 4 ) LINE_MATERIAL(LINE_MAX), material index for
!    each line.
!
!    Input/output, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
!    Output, integer ( kind = 4 ) TEXT_NUM, the number of lines read from the file.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) line_max

  logical add_on
  integer ( kind = 4 ), parameter :: BLANK = 0
  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  logical done
  character ( len = * ) filein_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lchar
  character ( len = 256 ) line
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ), parameter :: NONBLANK = 1
  integer ( kind = 4 ), parameter :: OFFSET = 1
  real ( kind = 4 ) temp
  integer ( kind = 4 ) text_num
  character ( len = 100 ) word

  ierror = 0
  text_num = 0
  word = ' '
  last = BLANK
!
!  Read a line of text from the file.
!
  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    text_num = text_num + 1
!
!  If this line begins with '#' , then it's a comment.  Read a new line.
!
    if ( line(1:1) == '#' ) then
      cycle
    end if
!
!  If this line is blank, then record that information.
!
    if ( len_trim ( line ) == 0 ) then
      last = BLANK
      cycle
    end if
!
!  This line records a node's coordinates.
!
    cor3_num = cor3_num + 1

    if ( cor3_num <= cor3_max ) then

      done = .true.

      do i = 1, 3

        call word_next_read ( line, word, done )

        call s_to_r4 ( word, temp, ierror, lchar )

        if ( ierror /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'XYZ_READ - Fatal error!'
          write ( *, '(a,i6)' ) '  S_TO_R4 returned IERROR = ', ierror
          write ( *, '(a,i6)' ) '  Reading (X,Y,Z) component ', i
          write ( *, '(a,i6)' ) '  Line number ', text_num
          exit
        end if

        cor3(i,cor3_num) = temp

      end do

      if ( ierror /= 0 ) then
        exit
      end if

    end if
!
!  If the previous input line was not blank, then this point was connected
!  to the previous one.
!
    if ( last == NONBLANK ) then

      add_on = .false.

      if ( line_num >= 2 ) then
        if ( line_dex(line_num-1) == cor3_num - 1 .and. &
             line_dex(line_num) == 0 ) then
          add_on = .true.
        end if
      end if

      if ( add_on ) then
        line_dex(line_num) = cor3_num
        line_num = line_num + 1
        line_dex(line_num) = 0
      else
        line_num = line_num + 1
        line_dex(line_num) = cor3_num - 1
        line_num = line_num + 1
        line_dex(line_num) = cor3_num
        line_num = line_num + 1
        line_dex(line_num) = 0
      end if

    end if

    last = NONBLANK

  end do
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) 'XYZ_READ - Read ', text_num, ' text lines from ' &
    // trim ( filein_name )

  return
end
subroutine xyz_write ( cor3, cor3_max, cor3_num, filein_name, fileout_name, &
  iunit, line_dex, line_max, line_num )

!*****************************************************************************80
!
!! XYZ_WRITE writes graphics data to an XYZ file.
!
!  Discussion:
!
!    Comments begin with a "#" in column 1.
!    Points listed consecutively are to be joined by lines.
!    Isolated points are separated by blanks.
!
!  Example:
!
!     # cube.xyz
!     #
!     #  First the points:
!
!     0 0 0
!
!     0 0 1
!
!     0 1 0
!
!     0 1 1
!
!     1 0 0
!
!     1 0 1
!
!     1 1 0
!
!     1 1 1
!
!     #
!     #  Now the lines
!     #
!     0 0 0
!     0 0 1
!     0 1 1
!     0 1 0
!     0 0 0
!
!     1 0 0
!     1 0 1
!     1 1 0
!     1 1 1
!     1 0 0
!
!     0 0 0
!     1 0 0
!
!     0 0 1
!     1 0 1
!
!     0 1 0
!     1 1 0
!
!     0 1 1
!     1 1 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) COR3(3,COR3_MAX), the coordinates of points.
!
!    Input, integer ( kind = 4 ) COR3_MAX, the maximum number of points.
!
!    Input, integer ( kind = 4 ) COR3_NUM, the number of points.
!
!    Input, character ( len = * ) FILEIN_NAME, the name of the input file.
!
!    Input, character ( len = * ) FILEOUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit to which output is written.
!
!    Input, integer ( kind = 4 ) LINE_DEX(LINE_MAX), nodes forming a line, terminated by -1.
!
!    Input, integer ( kind = 4 ) LINE_MAX, the maximum number of line definition items.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of line definition items.
!
  implicit none

  integer ( kind = 4 ) cor3_max
  integer ( kind = 4 ) line_max

  real ( kind = 4 ) cor3(3,cor3_max)
  integer ( kind = 4 ) cor3_num
  character ( len = * ) filein_name
  character ( len = * ) fileout_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) line_dex(line_max)
  integer ( kind = 4 ) line_num
  logical newline
  integer ( kind = 4 ), parameter :: OFFSET = 1
  integer ( kind = 4 ) text_num

  write ( iunit, '(a)' ) '# ' // trim ( fileout_name ) // ' created by IVREAD.'
  write ( iunit, '(a)' ) '#'
  text_num = 2
!
!  Print the points.
!
  do j = 1, cor3_num
    write ( iunit, '(a)' ) ' '
    write ( iunit, '(3g12.4)' ) cor3(1:3,j)
    text_num = text_num + 2
  end do
!
!  Print the line data.
!
  newline = .TRUE.

  do j = 1, line_num
 
    k = line_dex(j) - OFFSET

    if ( k == -1 ) then

      newline = .TRUE.

    else 

      if ( newline ) then
        write ( iunit, '(a)' ) ' '
        text_num = text_num + 1
        newline = .FALSE.
      end if

      write ( iunit, '(3g12.4)' ) cor3(1:3,k+OFFSET)

      text_num = text_num + 1

    end if

  end do
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XYZ_WRITE - Wrote ', text_num, ' text lines to ' &
    // trim ( fileout_name )
 
  return
end
