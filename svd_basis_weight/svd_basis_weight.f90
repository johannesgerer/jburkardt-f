program main

!*****************************************************************************80
!
!! MAIN is the main program for SVD_BASIS_WRITE.
!
!  Discussion:
!
!    SVD_BASIS_WEIGHT forms a basis from the SVD of weighted data vectors.
!
!    This program uses the singular value decomposition (SVD) to analyze
!    a set of weighted data, and extract a number of dominant modes.
!
!    This program is intended as an intermediate application, in
!    the following situation:
!
!    A) a "high fidelity" or "high resolution" PDE solver is used
!       to determine many (say N = 500) solutions of a discretized
!       PDE at various times, or parameter values.  Each solution
!       may be regarded as an M vector.  Typically, each solution
!       involves an M by M linear system, greatly reduced in
!       complexity because of bandedness or sparsity.
!
!    B) The user determines a weight vector W, with one value assigned
!       to each solution or vector.  Depending on the problem, the
!       weights might be chosen beforehand, or computed automatically
!       by some natural system, perhaps related to a varying time step
!       size, or other reasons.
!
!    C) This program is applied to extract L dominant modes from
!       the N weighted solutions.  This is done using the singular value
!       decomposition of the M by N matrix, each of whose columns
!       is one of the original solution vectors scaled by the weight.
!
!    D) a "reduced order model" program may then attempt to solve
!       a discretized version of the PDE, using the L dominant
!       modes as basis vectors.  Typically, this means that a dense
!       L by L linear system will be involved.
!
!    Thus, the program might read in 500 solution files, and a
!    weight file, and write out 5 or 10 files of the corresponding size
!    and "shape", representing the dominant solution modes.
!
!    An option has been added to compute the average of the vectors
!    and subtract it before SVD processing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gal Berkooz, Philip Holmes, John Lumley,
!    The proper orthogonal decomposition in the analysis of turbulent flows,
!    Annual Review of Fluid Mechanics,
!    Volume 25, 1993, pages 539-575.
!
!    John Burkardt, Max Gunzburger, Hyung-Chun Lee,
!    Centroidal Voronoi Tessellation-Based Reduced-Order
!    Modelling of Complex Systems,
!    SIAM Journal on Scientific Computing,
!    Volume 28, Number 2, 2006, pages 459-484.
!
!    Lawrence Sirovich,
!    Turbulence and the dynamics of coherent structures, Parts I-III,
!    Quarterly of Applied Mathematics,
!    Volume 45, Number 3, 1987, pages 561-590.
!
  implicit none

  integer ( kind = 4 ), parameter :: data_file_base_max = 20

  character average_char
  logical :: average_normalization
  real ( kind = 8 ) average_value
  character ( len = 255 ) basis_file
  integer ( kind = 4 ) basis_num
  logical, parameter :: clean = .true.
  logical comment
  character comment_char
  integer ( kind = 4 ) comp_num
  character ( len = 255 ) data_file
  integer ( kind = 4 ) data_file_num
  character ( len = 255 ) data_file_base(data_file_base_max)
  integer ( kind = 4 ) data_file_base_num
  integer ( kind = 4 ) dim_num
  logical file_exist
  character ( len = 255 ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) node_num
  real ( kind = 8 ), allocatable, dimension ( :, :) :: point
  real ( kind = 8 ), allocatable, dimension ( : ) :: point_average
  integer ( kind = 4 ) point_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: sval
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: table
  real ( kind = 8 ) tol
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SVD_BASIS_WEIGHT'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Given a PDE for which:'
  write ( *, '(a)' ) '    C is the number of components of the solution '
  write ( *, '(a)' ) '      at any single point,'
  write ( *, '(a)' ) '    P is the number of points where a solution is given,'
  write ( *, '(a)' ) '    N is the number of solution vectors,'
  write ( *, '(a)' ) '    L is the number of modes to be extracted.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Then we let M = C*P be the abstract spatial dimension.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  If requested, we compute the average solution,'
  write ( *, '(a)' ) '  subtract it from each solution, and save that'
  write ( *, '(a)' ) '  as mode #0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Set up A, the M by N matrix of solution vectors,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Get A = U * S * V'', the singular value decomposition.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The first L columns of U are our dominant modes.'
  write ( *, '(a)' ) ' '
!
!  What is the basis size?
!
  call i4_input ( '  How many basis vectors (L) are to be extracted?',&
    basis_num, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVD_BASIS_WEIGHT - Fatal error!'
    write ( *, '(a)' ) '  Input error reading the basis size.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVD_BASIS_WEIGHT:'
    write ( *, '(a)' ) '  Abnormal end of execution.'
    write ( *, '(a)' ) ' '
    call timestamp ( )
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  L = ', basis_num
!
!  Gather one or more "base" file names.
!
  data_file_base_num = 0

  do

    if ( data_file_base_max <= data_file_base_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  No more base file names can be entered.'
      exit
    end if
!
!  Get the next base file name.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  You specify a consecutive sequence of file names'
    write ( *, '(a)' ) '  by giving the first "base" file name.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  If there are no more sequences to enter,'
    write ( *, '(a)' ) '  just hit RETURN.'

    file_name = ' '

    call s_input ( '  Enter a new base file name, or RETURN.', &
      file_name, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SVD_BASIS_WEIGHT - Fatal error!'
      write ( *, '(a)' ) '  Input error reading the base file name.'
      stop
    end if

    if ( len_trim ( file_name ) <= 0 .or. file_name == ' ' ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  RETURN was entered.'
      write ( *, '(a)' ) '  Presumably, there are no more file sequences.'
      exit
    end if

    data_file_base_num = data_file_base_num + 1
    data_file_base(data_file_base_num) = file_name

    write ( *, '(a)' ) ' '
    write ( *, '(i8,a)' ) data_file_base_num, ':  "' &
      // trim ( file_name ) // '".'
!
!  For the very first base file, get the data sizes.
!
    if ( data_file_base_num == 1 ) then

      call r8table_header_read ( file_name, comp_num, node_num )

      dim_num = comp_num * node_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  According to the first base file,'
      write ( *, '(a,i8)' ) &
        '  The number of solution components C =   ', comp_num
      write ( *, '(a,i8)' ) &
        '  The number of solution points P =       ', node_num
      write ( *, '(a,i8)' ) &
        '  The "size" of each solution M = (C*P) = ', dim_num
!
!  Idiocy check.  L must be less than or equal to M.
!
      if ( dim_num < basis_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SVD_BASIS_WEIGHT - Fatal error!'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  M < L.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) &
          '  That is, the number of modes requested (L) is greater'
        write ( *, '(a)' ) '  than the spatial dimension (M).'
        write ( *, '(a)' ) '  Technically, the program could pad out the answer'
        write ( *, '(a)' ) '  with L-M zero vectors, but instead, we will stop'
        write ( *, '(a)' ) '  assuming you made an error, or a misapprehension.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SVD_BASIS_WEIGHT:'
        write ( *, '(a)' ) '  Abnormal end of execution.'
        write ( *, '(a)' ) ' '
        call timestamp ( )
        stop
      end if

    end if

  end do
!
!  Count all the data files.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Counting the data files for each set.'
  write ( *, '(a)' ) ' '

  data_file_num = 0

  do i = 1, data_file_base_num

    data_file = data_file_base(i)

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Data set # ', i
    write ( *, '(a)' ) '  begins with file "' // trim ( data_file ) // '".'

    do

      if ( .not. file_exist ( data_file ) ) then
        exit
      end if

      data_file_num = data_file_num + 1

      call file_name_inc ( data_file )

    end do

    write ( *, '(a)' ) '  and terminates because there is no file'
    write ( *, '(a)' ) '  "' // trim ( data_file ) // '".'

  end do

  if ( data_file_num == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVD_BASIS_WEIGHT - Fatal error!'
    write ( *, '(a)' ) '  There do not seem to be any solution files;'
    write ( *, '(a)' ) '  that is, files whose names are "incremented"'
    write ( *, '(a)' ) '  versions of the first file name.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The first file we looked for was "' // &
      trim ( data_file ) // '".'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVD_BASIS_WEIGHT:'
    write ( *, '(a)' ) '  Abnormal end of execution.'
    write ( *, '(a)' ) ' '
    call timestamp ( )
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of data files N = ', data_file_num
!
!  Set up an array to hold all the data.
!
  point_num = data_file_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data is stored in an M by N matrix A.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The "spatial" dimension M is   ', dim_num
  write ( *, '(a,i8)' ) '  The number of data points N is ', point_num
!
!  Allocate space for the POINT array.
!
  allocate ( point(1:dim_num,1:point_num) )
  allocate ( table(1:comp_num,1:node_num) )
!
!  Read the data.
!
  l = 0

  do ii = 1, data_file_base_num

    data_file = data_file_base(ii)

    do

      if ( .not. file_exist ( data_file ) ) then
        exit
      end if

      l = l + 1

      call r8table_data_read ( data_file, comp_num, node_num, table )

      k = 0
      do j = 1, node_num
        do i = 1, comp_num
          k = k + 1
          point(k,l) = table(i,j)
        end do
      end do

      call file_name_inc ( data_file )

    end do

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data has been read into the matrix A.'
!
!----------------------------------------------------------------------------
!
!  Get the weights and scale the data.
!
!----------------------------------------------------------------------------
!
  call s_input ( '  Enter the weight file name:', file_name, ierror )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Weight file name is "' // trim ( file_name ) // '".'

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVD_BASIS_WEIGHT - Fatal error!'
    write ( *, '(a)' ) '  Input error reading the weight file name.'
    stop
  end if

  allocate ( weight(1:point_num) )

  call r8table_data_read ( file_name, 1, point_num, weight )

  do j = 1, point_num
    point(1:dim_num,j) = point(1:dim_num,j) * weight(j)
  end do
!
!----------------------------------------------------------------------------
!  Optionally, average the data, subtract the average from each entry,
!  and later save the average as vector #0.
!----------------------------------------------------------------------------
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SVD_BASIS_WEIGHT:'
  write ( *, '(a)' ) '  Averaging of data is optional.'
  write ( *, '(a)' ) '  The program can average the data vectors,'
  write ( *, '(a)' ) '  subtract it from each data vector,'
  write ( *, '(a)' ) '  and write out the data average vector as an'
  write ( *, '(a)' ) '  extra "mode 0" vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Do you want to compute and use the average? (Y/N)'

  call s_input ( '  Enter Y or N:', average_char, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVD_BASIS_WEIGHT - Warning!'
    write ( *, '(a)' ) '  Input error reading the average option.'
    write ( *, '(a)' ) '  We will assume averaging is NOT used.'
    average_char = 'N'
  end if

  if ( average_char == 'Y' .or. average_char == 'y' ) then
    average_normalization = .true.
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The user has requested the average vector.'
  else
    average_normalization = .false.
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The user does not want the average vector.'
  end if

  if ( average_normalization ) then

    allocate ( point_average(1:dim_num) )

    do i = 1, dim_num
      point_average(i) = sum ( point(i,1:point_num) )
    end do

    point_average(1:dim_num) = point_average(1:dim_num) &
      / real ( point_num, kind = 8 )

    do i = 1, dim_num
      point(i,1:point_num) = point(i,1:point_num) - point_average(i)
    end do

  end if
!
!----------------------------------------------------------------------------
!
!  Compute the SVD of A.
!
!----------------------------------------------------------------------------
!
  allocate ( sval(1:basis_num) )

  call singular_vectors ( dim_num, point_num, basis_num, point, sval )
!
!----------------------------------------------------------------------------
!
!  "Clean" the output data.  We are having problems with some vectors
!  containing a few very tiny (and meaningless) values.
!
!----------------------------------------------------------------------------
!
  if ( clean ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Because the CLEAN option is on,'
    write ( *, '(a)' ) '  we will set very tiny vector entries to 0.'

    tol = epsilon ( tol )

    do j = 1, basis_num
      do i = 1, dim_num
        if ( abs ( point(i,j) ) < tol ) then
          point(i,j) = 0.0D+00;
        end if
      end do
    end do
  end if
!
!----------------------------------------------------------------------------
!
!  Write the first L left singular vectors (columns of U) to files.
!
!----------------------------------------------------------------------------
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SVD_BASIS_WEIGHT:'
  write ( *, '(a)' ) '  Ready to write the left singular vectors to files.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Do you want comments in the header of the file?'
  write ( *, '(a)' ) '  (These begin with the "#" character.) (Y/N)'

  call s_input ( '  Enter Y or N:', comment_char, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SVD_BASIS_WEIGHT - Warning!'
    write ( *, '(a)' ) '  Input error reading the comment option.'
    write ( *, '(a)' ) '  We will assume comments are acceptable.'
    comment_char = 'Y'
  end if

  if ( comment_char == 'Y' .or. comment_char == 'y' ) then
    comment = .true.
  else
    comment = .false.
  end if

  basis_file = 'svd_000.txt'

  if ( average_normalization ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Writing average file "' // trim ( basis_file ) // '".'

    average_value = 0.0D+00

    call basis_write ( basis_file, comp_num, node_num, average_value, &
      point_average(1:dim_num), comment )

  end if

  do j = 1, basis_num

    call file_name_inc ( basis_file )

    if ( j == 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Writing first file "' // trim ( basis_file ) // '".'
    end if

    if ( j == basis_num ) then
      write ( *, '(a)' ) '  Writing last file  "' // trim ( basis_file ) // '".'
    end if

    call basis_write ( basis_file, comp_num, node_num, sval(j), &
      point(1:dim_num,j), comment )

  end do

  deallocate ( point )
  if ( allocated ( point_average ) ) then
    deallocate ( point_average )
  end if
  deallocate ( sval )
  deallocate ( table )
  deallocate ( weight )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SVD_BASIS_WEIGHT'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine basis_write ( file_out_name, m, n, s, u, comment )

!*****************************************************************************80
!
!! BASIS_WRITE writes a basis vector to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the file to write.
!
!    Input, integer ( kind = 4 ) M, the number of data components.
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, real ( kind = 8 ) S, the associated singular value.
!
!    Input, real ( kind = 8 ) U(M,N), the data values.
!
!    Input, logical COMMENT, is TRUE if comments are to be included.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  logical comment
  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  character ( len = 10 ) form
  integer ( kind = 4 ) j
  real ( kind = 8 ) s
  character ( len = 40 ) string
  real ( kind = 8 ) u(m,n)

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace' )

  if ( comment ) then

    call timestring ( string )

    write ( file_out_unit, '(a)'       ) '#  ' // trim ( file_out_name )
    write ( file_out_unit, '(a)'       ) '#  created by BASIS_WRITE.F90,'
    write ( file_out_unit, '(a)'       ) '#  part of SVD_BASIS.F90,'
    write ( file_out_unit, '(a)'       ) '#'
    write ( file_out_unit, '(a)'       ) '#  Created on ' // trim ( string )
    write ( file_out_unit, '(a)'       ) '#'
    write ( file_out_unit, '(a,i8)'    ) '#  Number of components M = ', m
    write ( file_out_unit, '(a,i8)'    ) '#  Number of items N      = ', n
    write ( file_out_unit, '(a,g15.6)' ) '#  Singular value S = ', s
    write ( file_out_unit, '(a)'       ) '#'

  end if

  write ( form, '( ''('',i2,''g15.6)'' )' ) m

  do j = 1, n

    write ( file_out_unit, form ) u(1:m,j)

  end do

  close ( unit = file_out_unit )

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

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical ch_eqi

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
function ch_is_digit ( c )

!*****************************************************************************80
!
!! CH_IS_DIGIT is TRUE if a character is a decimal digit.
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
!    Input, character C, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, .TRUE. if C is a digit, .FALSE. otherwise.
!
  implicit none

  character c
  logical ch_is_digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
    ch_is_digit = .true.
  else
    ch_is_digit = .false.
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
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.
!    If C was 'illegal', then DIGIT is -1.
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
subroutine digit_inc ( c )

!*****************************************************************************80
!
!! DIGIT_INC increments a decimal digit.
!
!  Example:
!
!    Input  Output
!    -----  ------
!    '0'    '1'
!    '1'    '2'
!    ...
!    '8'    '9'
!    '9'    '0'
!    'A'    'A'
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
!    Input/output, character C, a digit to be incremented.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  call ch_to_digit ( c, digit )

  if ( digit == -1 ) then
    return
  end if

  digit = digit + 1

  if ( digit == 10 ) then
    digit = 0
  end if

  call digit_to_ch ( digit, c )

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
subroutine file_column_count ( input_filename, column_num )

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
!    Input, character ( len = * ) INPUT_FILENAME, the name of the file.
!
!    Output, integer ( kind = 4 ) COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer ( kind = 4 ) column_num
  logical got_one
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 256 ) line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    ' // trim ( input_filename )
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
function file_exist ( file_name )

!*****************************************************************************80
!
!! FILE_EXIST reports whether a file exists.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, logical FILE_EXIST, is TRUE if the file exists.
!
  implicit none

  character ( len = * ) file_name
  logical file_exist

  inquire ( file = file_name, exist = file_exist )

  return
end
subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC generates the next filename in a series.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected, and
!    if the name contains no digits, then nothing is done.
!
!  Example:
!
!      Input          Output
!      -----          ------
!      a7to11.txt     a7to12.txt
!      a7to99.txt     a8to00.txt
!      a9to99.txt     a0to00.txt
!      cat.txt        cat.txt
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
!    Input/output, character ( len = * ) FILE_NAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  logical ch_is_digit
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( file_name )

  do i = lens, 1, -1

    c = file_name(i:i)

    if ( ch_is_digit ( c ) ) then

      call digit_inc ( c )

      file_name(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do

  return
end
subroutine file_row_count ( input_filename, row_num )

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
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
!
  implicit none

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) comment_num
  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 100 ) line
  integer ( kind = 4 ) record_num
  integer ( kind = 4 ) row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( input_filename )
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
subroutine i4_input ( string, value, ierror )

!*****************************************************************************80
!
!! I4_INPUT prints a prompt string and reads an I4 from the user.
!
!  Discussion:
!
!    If the input line starts with a comment character ('#') or is
!    blank, the routine ignores that line, and tries to read the next one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the prompt string.
!
!    Output, integer ( kind = 4 ) VALUE, the value input by the user.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero
!    if no error occurred.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) last
  character ( len = 80 ) line
  character ( len = * ) string
  integer ( kind = 4 ) value

  ierror = 0
  value = huge ( value )
!
!  Write the prompt.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  do

    read ( *, '(a)', iostat = ierror ) line

    if ( ierror /= 0 ) then
      return
    end if
!
!  If the line begins with a comment character, go back and read the next line.
!
    if ( line(1:1) == '#' ) then
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if
!
!  Extract integer information from the string.
!
    call s_to_i4 ( line, value, ierror, last )

    if ( ierror /= 0 ) then
      value = huge ( value )
      return
    end if

    exit

  end do

  return
end
subroutine r8table_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! R8TABLE_HEADER_READ reads the header from a double precision table file.
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

  character ( len = * ) input_filename
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  call file_column_count ( input_filename, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8TABLE_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8TABLE_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  return
end
subroutine r8table_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! R8TABLE_DATA_READ reads data from a real table file.
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
!    08 October 2003
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

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  character ( len = 255 ) line
  real ( kind = 8 ) table(m,n)
  real ( kind = 8 ) x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8TABLE_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( input_filename )
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = 1
      exit
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
subroutine s_input ( string, value, ierror )

!*****************************************************************************80
!
!! S_INPUT prints a prompt string and reads a string from the user.
!
!  Discussion:
!
!    If the input line starts with a comment character ('#'),
!    the routine ignores that line, and tries to read the next one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the prompt string.
!
!    Output, character ( len = * ) VALUE, the value input by the user.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero
!    if no error occurred.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = 80 ) line
  character ( len = * ) string
  character ( len = * ) value

  ierror = 0
  value = ' '
!
!  Write the prompt.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  do

    line = ' '

    read ( *, '(a)', iostat = ierror ) line

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'S_INPUT: Fatal error!'
      write ( *, '(a)' ) '  Input error!'
      stop
    end if
!
!  If the line begins with a comment character, go back and read the next line.
!
    if ( line(1:1) == '#' ) then
      cycle
    end if

    if ( len_trim ( line ) <= 0 .or. line == ' ' ) then
      value = ' '
      exit
    end if

    value = line

    exit

  end do

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an integer value from a string.
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
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character of S used
!    to make IVAL.
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
subroutine s_to_r8 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 value from a string.
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
!    ' 2D-1'           0.2
!    '23.45'           23.45
!    '-4.2D+2'         -420.0
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
!    characters can be read to form a legal number.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) R, the value that was read from the string.
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

  character c
  logical ch_eqi
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
    if ( iterm == 1 .or. nchar <= lchar + 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar + 1 == nchar ) then
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
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) lchar
  real ( kind = 8 ) rvec(n)
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1
  ihi = len ( s )

  if ( ihi < 1 ) then
    ierror = 1
    return
  end if

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:ihi), rvec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = 2
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
subroutine singular_vectors ( m, n, basis_num, a, sval )

!*****************************************************************************80
!
!! SINGULAR_VECTORS computes the desired singular values.
!
!  Discussion:
!
!    The LAPACK SVD routine DGESVD is used to compute the singular
!    value decomposition:
!
!      A = U * S * V'
!
!    The specification of LWORK was recently (18 July 2007) corrected, from
!
!      3 * min ( m, n ) + max ( max ( m, n ), 2 * min ( m, n ) )
!
!    to
!
!      3 * min ( m, n ) + max ( max ( m, n ), 5 * min ( m, n ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford,
!    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum,
!    Sven Hammarling, Alan McKenney, Danny Sorensen,
!    LAPACK User's Guide,
!    Third Edition,
!    SIAM, 1999,
!    ISBN: 0898714478,
!    LC: QA76.73.F25L36
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input, integer ( kind = 4 ) BASIS_NUM, the number of basis vectors
!    to be extracted.
!
!    Input/output, real ( kind = 8 ) A(M,N); on input, the matrix whose
!    singular values are to be computed.  On output, A(M,1:BASIS_NUM)
!    contains the first BASIS_NUM left singular vectors.
!
!    Output, real ( kind = 8 ) SVAL(BASIS_NUM), the first BASIS_NUM
!    singular values.
!
  implicit none

  integer ( kind = 4 ) basis_num
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldvt
  character jobu
  character jobvt
  integer ( kind = 4 ) lwork
  real ( kind = 8 ) s(min(m,n))
  real ( kind = 8 ) sval(basis_num)
  real ( kind = 8 ) u(1,1)
  real ( kind = 8 ) vt(1,1)
  real ( kind = 8 ) work(3*min(m,n)+max(max(m,n),5*min(m,n)))

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SINGULAR_VECTORS'
  write ( *, '(a)' ) '  For an MxN matrix A in general storage,'
  write ( *, '(a)' ) '  The LAPACK routine DGESVD computes the '
  write ( *, '(a)' ) '  singular value decomposition:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    A = U * S * V'''
  write ( *, '(a)' ) ' '
!
!  Compute the eigenvalues and eigenvectors.
!
!  JOBU = 'O' means that the first min ( M, N ) columns of U
!  are to be computed and stored in the first min ( M, N ) columns of A.
!
  jobu = 'O'
!
!  JOBVT = 'N' means that V is not to be computed.
!
  jobvt = 'N'
  lda = m
  ldu = m
  ldvt = n
  lwork = 3 * min ( m, n ) + max ( max ( m, n ), 5 * min ( m, n ) )

  call dgesvd ( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, &
    lwork, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SINGULAR_VECTORS - Warning:'
    write ( *, '(a,i8)' ) '  DGESVD returned nonzero INFO = ', info
    return
  end if
!
!  Copy out the first BASIS_NUM singular values.
!
  sval(1:basis_num) = s(1:basis_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The leading singular values:'
  write ( *, '(a)' ) ' '

  do i = 1, basis_num
    write ( *, '(2x,i4,2x,g16.8)' ) i, sval(i)
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
subroutine timestring ( string )

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = 'May 31 2001   9:45:54.872 AM'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  character ( len = * ) string
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

  write ( string, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
