program main

!*****************************************************************************80
!
!! MAIN is the main program for the meshless basis function routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 2006
!
!  Author:
!
!    Lili Ju
!
  implicit none

  real      ( kind = 8 ), allocatable :: a(:)
  real      ( kind = 8 ), dimension ( 2 ) :: alpha = (/ 0.0D+00, 1.0D+00 /)
  real      ( kind = 8 ), allocatable :: b(:)
  integer   ( kind = 4 ), allocatable :: base(:)
  real      ( kind = 8 ), allocatable :: basis_center(:,:)
  character ( len = 80 ) :: basis_file = 'basis.dat'
  integer   ( kind = 4 ) :: basis_m = 0
  real      ( kind = 8 ), allocatable :: basis_radius(:)
  real      ( kind = 8 ) :: basis_radius_factor = 1.0D+00
  real      ( kind = 8 ), dimension ( 2 ) :: beta = (/ 0.0D+00, 1.0D+00 /)
  integer   ( kind = 4 ) :: center_choice = 0
  character ( len = 80 ) command
  real      ( kind = 8 ), allocatable :: cvt_center(:,:)
  character ( len = 80 ) :: cvt_file = 'cvt.dat'
  integer   ( kind = 4 ) cvt_m
  real      ( kind = 8 ), allocatable :: cvt_radius(:)
  real      ( kind = 8 ) :: cvt_radius_factor = 1.0D+00
  integer   ( kind = 4 ) :: density_function = 0
  character ( len = 80 ) :: halton_file = 'halton.dat'
  real      ( kind = 8 ), allocatable :: halton_points(:,:)
  integer   ( kind = 4 ) i_temp
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ios
  integer   ( kind = 4 ) last
  integer   ( kind = 4 ) :: maxit = 0
  integer   ( kind = 4 ) :: myrank = 0
  integer   ( kind = 4 ) :: n = 0
  integer   ( kind = 4 ) ncolumn
  integer   ( kind = 4 ) :: m = 2
  integer   ( kind = 4 ) :: ns = 0
  integer   ( kind = 4 ) :: radius_choice = 0
  logical                s_eqi
  integer   ( kind = 4 ) :: skip = 0
  character ( len = 80 ) :: uniform_file = 'uniform.dat'
  real      ( kind = 8 ), allocatable :: uniform_points(:,:)
  character ( len = 80 ) what

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MESHLESS'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Meshless basis function generation.'
!
!  Initialize the random number generator.
!
  call set_random_seed ( myrank )
!
!  Read a command, do a command.
!
  do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter command (or HELP for a list)'

    do

      read ( *, '(a)', iostat = ios ) command
      write ( *, '(a)' ) '  Command: " ' // trim ( command ) // '".'

      if ( command(1:1) /= '#' ) then
        exit
      end if

    end do

    call s_blank_delete ( command )

    if ( ios /= 0 ) then
      exit
    end if

    if ( s_eqi ( command, 'BASIS_MAKE' ) ) then

      if ( m <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BASIS_MAKE - Error!'
        write ( *, '(a)' ) '  The spatial dimension must be positive!'
        write ( *, '(a)' ) '  Enter the spatial dimension M with the command'
        write ( *, '(a)' ) '    "M = ???"'
        cycle
      end if

      if ( n <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BASIS_MAKE - Error!'
        write ( *, '(a)' ) '  The number of points N must be positive!'
        write ( *, '(a)' ) '  Enter the number of points with the command'
        write ( *, '(a)' ) '    "N = ???"'
        cycle
      end if

      basis_m = 7 * int ( sqrt ( real ( n, kind = 8 ) ) )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'Mysterious parameter BASIS_M = ', basis_m

      do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Initial center point determination:'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  1: use a uniform distribution;'
        write ( *, '(a)' ) '  2: use a Halton distribution.'

        read ( *, * ) center_choice

        if ( center_choice == 1 ) then

          call uniform_make ( m, n, a, b, density_function, &
            basis_center )
          exit

        else if ( center_choice == 2 ) then

          call halton_generate ( m, n, skip, base, basis_center )
          exit
        end if

      end do

      do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Initial radius determination:'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  1: radius method 1;'
        write ( *, '(a)' ) '  2: radius method 2.'

        read ( *, * ) radius_choice

        if ( radius_choice == 1 ) then

          if ( m <= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'RADIUS_MAKE_1 - Error!'
            write ( *, '(a)' ) '  The spatial dimension must be positive!'
            write ( *, '(a)' ) &
              '  Enter the spatial dimension M with the command'
            write ( *, '(a)' ) '    "M = ???"'
            exit
          end if

          call radius_make_1 ( m, a, b, n, basis_m, density_function, &
            basis_center, basis_radius )

          exit

        else if ( radius_choice == 2 ) then

          if ( m <= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'RADIUS_MAKE_2 - Error!'
            write ( *, '(a)' ) '  The spatial dimension must be positive!'
            write ( *, '(a)' ) &
              '  Enter the spatial dimension M with the command'
            write ( *, '(a)' ) '    "M = ???"'
            exit
          end if

          call radius_make_2 ( m, a, b, n, basis_m, basis_center, &
            basis_radius )

          exit
        end if

      end do

    else if ( s_eqi ( command, 'BASIS_OVERLAP' ) ) then

      call basis_overlap ( m, n, basis_center, basis_radius )

    else if ( s_eqi ( command, 'BASIS_PLOT' ) ) then

      call basis_plot ( m, n, basis_center, basis_radius )

    else if ( s_eqi ( command, 'BASIS_READ' ) ) then

      call file_column_count ( basis_file, ncolumn )

      m = ncolumn - 1

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'Setting M to ', m

      call file_line_count ( basis_file, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'Setting number of center points to ', n

      call basis_read ( m, n, basis_center, basis_radius, basis_file )

      call basis_box ( m, n, basis_center, basis_radius, a, b )

    else if ( s_eqi ( command, 'BASIS_SCALE' ) ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Enter a scale factor to apply to the radius'
      write ( *, '(a)' ) 'of the current basis function set:'

      read ( *, *, iostat = ios ) basis_radius_factor
      if ( ios /= 0 ) then
        exit
      end if

      basis_radius(1:n) = basis_radius_factor * basis_radius(1:n)

    else if ( s_eqi ( command, 'BASIS_WRITE' ) ) then

      if ( n <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Use the command radius_make first!'
        cycle
      end if

      call basis_write ( m, n, basis_center, basis_radius, basis_file )

    else if ( s_eqi ( command, 'CENTER_PLOT' ) ) then

      call center_plot ( m, n, basis_center, a, b )

    else if ( s_eqi ( command, 'CVT_MAKE' ) ) then

      if ( m <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CVT_MAKE - Error!'
        write ( *, '(a)' ) '  The spatial dimension must be positive!'
        write ( *, '(a)' ) '  Enter the spatial dimension M with the command'
        write ( *, '(a)' ) '    "M = ???"'
        cycle
      end if

      if ( n <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CVT_MAKE - Error!'
        write ( *, '(a)' ) '  The number of points N must be positive!'
        write ( *, '(a)' ) '  Enter the number of points with the command'
        write ( *, '(a)' ) '    "N = ???"'
        cycle
      end if

      cvt_m = 7 * int ( sqrt ( real ( n, kind = 8 ) ) )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'Mysterious parameter CVT_M = ', cvt_m

      do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Choose a basis center initialization:'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  1: use a uniform distribution;'
        write ( *, '(a)' ) '  2: use a Halton distribution.'

        read ( *, * ) center_choice

        if ( center_choice == 1 ) then
          call uniform_make ( m, n, a, b, density_function, cvt_center )
          exit
        else if ( center_choice == 2 ) then

          call halton_generate ( m, n, skip, base, cvt_center )
          exit
        end if

      end do

      do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Initial radius determination:'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  1: radius method 1;'
        write ( *, '(a)' ) '  2: radius method 2.'

        read ( *, * ) radius_choice

        if ( radius_choice == 1 ) then

          if ( m <= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'RADIUS_MAKE_2 - Error!'
            write ( *, '(a)' ) '  The spatial dimension must be positive!'
            write ( *, '(a)' ) &
              '  Enter the spatial dimension M with the command'
            write ( *, '(a)' ) '    "M = ???"'
            exit
          end if

          call radius_make_1 ( m, a, b, n, cvt_m, density_function, &
            cvt_center, cvt_radius )
          exit

        else if ( radius_choice == 2 ) then

          if ( m <= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'RADIUS_MAKE_2 - Error!'
            write ( *, '(a)' ) '  The spatial dimension must be positive!'
            write ( *, '(a)' ) &
              '  Enter the spatial dimension M with the command'
            write ( *, '(a)' ) '    "M = ???"'
            exit
          end if

          call radius_make_2 ( m, a, b, n, cvt_m, cvt_center, &
            cvt_radius )
          exit
        end if

      end do

      call cvt_make ( alpha, beta, m, a, b, maxit, ns, density_function, &
        n, cvt_center )

    else if ( s_eqi ( command, 'CVT_OVERLAP' ) ) then

      call basis_overlap ( m, n, cvt_center, cvt_radius )

    else if ( s_eqi ( command, 'CVT_PLOT' ) ) then

      call basis_plot ( m, n, cvt_center, cvt_radius )

    else if ( s_eqi ( command, 'CVT_READ' ) ) then

      call cvt_read ( m, n, cvt_center, cvt_radius, cvt_file )

    else if ( s_eqi ( command, 'CVT_SCALE' ) ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Enter a scale factor to apply to the radius'
      write ( *, '(a)' ) 'of the current CVT function set:'

      read ( *, *, iostat = ios ) cvt_radius_factor
      if ( ios /= 0 ) then
        exit
      end if

      cvt_radius(1:n) = cvt_radius_factor * cvt_radius(1:n)

    else if ( s_eqi ( command, 'CVT_WRITE' ) ) then

      if ( n <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Use the command cvt_make first!'
        cycle
      end if

      call cvt_write ( m, n, cvt_center, cvt_radius, cvt_file )

    else if ( s_eqi ( command(1:17), 'DENSITY_FUNCTION=' ) ) then

      call s_to_i4 ( command(18:), i_temp, ierror, last )

      if ( ierror == 0 ) then

        density_function = i_temp

        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  DENSITY_FUNCTION has been set to ', &
          density_function

      end if

    else if ( s_eqi ( command, 'HALTON_MAKE' ) ) then

      if ( m <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HALTON_MAKE - Error!'
        write ( *, '(a)' ) '  The spatial dimension must be positive!'
        write ( *, '(a)' ) '  Enter the spatial dimension M with the command'
        write ( *, '(a)' ) '    "M = ???"'
        cycle
      end if

      if ( n <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HALTON_MAKE - Error!'
        write ( *, '(a)' ) '  The number of points N must be positive!'
        write ( *, '(a)' ) '  Enter the number of points with the command'
        write ( *, '(a)' ) '    "N = ???"'
        cycle
      end if

      call halton_generate ( m, n, skip, base, halton_points )

    else if ( s_eqi ( command, 'HALTON_READ' ) ) then

      call halton_read ( m, n, halton_points, halton_file )

    else if ( s_eqi ( command, 'HALTON_WRITE' ) ) then

      if ( n <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MESHLESS - Error!'
        write ( *, '(a)' ) '  No Halton basis functions have been defined.'
        write ( *, '(a)' ) '  Use the command HALTON_MAKE or HALTON_READ first!'
        cycle
      end if

      call halton_write ( m, n, skip, base, halton_points, halton_file )

    else if ( s_eqi ( command, 'HELP' ) ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Commands:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DENSITY_FUNCTION = ... specify density.'
      write ( *, '(a)' ) '  MAXIT = ...   specify number of iterations.'
      write ( *, '(a)' ) '  N =     ...   specify the number of basis points.'
      write ( *, '(a)' ) '  M =     ...   specify the spatial dimension.'
      write ( *, '(a)' ) '  NS =    ...   specify the sampling point density'
      write ( *, '(a)' ) '                for estimating areas.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  BASIS_MAKE    create a set of basis'
      write ( *, '(a)' ) '                functions.'
      write ( *, '(a)' ) '  BASIS_OVERLAP count the basis function overlaps.'
      write ( *, '(a)' ) '  BASIS_PLOT    plot a set of basis functions.'
      write ( *, '(a)' ) '  BASIS_READ    read a set of basis'
      write ( *, '(a)' ) '                functions from a file.'
      write ( *, '(a)' ) '  BASIS_SCALE   make the basis functions "wider"'
      write ( *, '(a)' ) '  BASIS_WRITE   write the current basis functions '
      write ( *, '(a)' ) '                to a file.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  CENTER_PLOT'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  CVT_MAKE      create a set of CVT basis'
      write ( *, '(a)' ) '                functions.'
      write ( *, '(a)' ) '  CVT_OVERLAP   count the basis function overlaps.'
      write ( *, '(a)' ) '  CVT_PLOT      plot a set of CVT basis functions.'
      write ( *, '(a)' ) '  CVT_READ      read a set of CVT basis'
      write ( *, '(a)' ) '                functions from a file.'
      write ( *, '(a)' ) '  CVT_SCALE     make CVT functions "wider")'
      write ( *, '(a)' ) '  CVT_WRITE     write the current CVT basis'
      write ( *, '(a)' ) '                functions to a file.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  HALTON_MAKE   create a set of Halton basis'
      write ( *, '(a)' ) '                functions.'
      write ( *, '(a)' ) '                Specify M and N first!'
      write ( *, '(a)' ) '  HALTON_READ   read a set of Halton basis'
      write ( *, '(a)' ) '                functions from a file.'
      write ( *, '(a)' ) '  HALTON_WRITE  write the current Halton basis '
      write ( *, '(a)' ) '                functions to a file.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  UNIFORM_MAKE  create a set of uniform basis'
      write ( *, '(a)' ) '                functions.'
      write ( *, '(a)' ) '  UNIFORM_READ  read a set of uniform basis'
      write ( *, '(a)' ) '                functions from a file.'
      write ( *, '(a)' ) '  UNIFORM_WRITE write the current uniform basis'
      write ( *, '(a)' ) '                functions to a file.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Help          print the list of commands.'
      write ( *, '(a)' ) '  Print         print the value of some quantity.'
      write ( *, '(a)' ) '  Quit          terminate the program.'

    else if ( s_eqi ( command(1:2), 'M=' ) ) then

      call s_to_i4 ( command(3:), i_temp, ierror, last )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MESHLESS - Error!'
        write ( *, '(a)' ) '  Could not set M!'
        cycle
      end if

      m = i_temp
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  M has been set to ', m

      if ( allocated ( a ) ) then
        deallocate ( a ) 
      end if
      allocate ( a(m) )
      a(1:m) = -1.0D+00

      if ( allocated ( b ) ) then
        deallocate ( b ) 
      end if
      allocate ( b(m) )
      b(1:m) = +1.0D+00

      if ( allocated ( base ) ) then
        deallocate ( base )
      end if
      allocate ( base(m) )

      if ( allocated ( basis_center ) ) then
        deallocate ( basis_center )
      end if
      if ( 0 < n ) then
        allocate ( basis_center(m,n) )
      end if

      if ( allocated ( cvt_center ) ) then
        deallocate ( cvt_center )
      end if
      if ( 0 < n ) then
        allocate ( cvt_center(m,n) )
      end if

      if ( allocated ( halton_points ) ) then
        deallocate ( halton_points )
      end if
      if ( 0 < n ) then
        allocate ( halton_points(m,n) )
      end if

      if ( allocated ( uniform_points ) ) then
        deallocate ( uniform_points )
      end if
      if ( 0 < n ) then
        allocate ( uniform_points(m,n) )
      end if

    else if ( s_eqi ( command(1:6), 'MAXIT=' ) ) then

      call s_to_i4 ( command(7:), i_temp, ierror, last )

      if ( ierror == 0 ) then

        maxit = i_temp
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  MAXIT has been set to ', maxit

      end if

    else if ( s_eqi ( command(1:2), 'N=' ) ) then

      call s_to_i4 ( command(3:), i_temp, ierror, last )

      if ( ierror == 0 ) then

        n = i_temp
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  N has been set to ', n

      end if

      if ( allocated ( basis_center ) ) then
        deallocate ( basis_center )
      end if
      if ( 0 < m ) then
        allocate ( basis_center(m,n) )
      end if

      if ( allocated ( basis_radius ) ) then
        deallocate ( basis_radius )
      end if
      allocate ( basis_radius(n) )

      if ( allocated ( cvt_center ) ) then
        deallocate ( cvt_center )
      end if
      if ( 0 < m ) then
        allocate ( cvt_center(m,n) )
      end if

      if ( allocated ( cvt_radius ) ) then
        deallocate ( cvt_radius )
      end if
      allocate ( cvt_radius(n) )

      if ( allocated ( halton_points ) ) then
        deallocate ( halton_points )
      end if
      if ( 0 < m ) then
        allocate ( halton_points(m,n) )
      end if

      if ( allocated ( uniform_points ) ) then
        deallocate ( uniform_points )
      end if
      if ( 0 < m ) then
        allocate ( uniform_points(m,n) )
      end if

    else if ( s_eqi ( command(1:3), 'NS=' ) ) then

      call s_to_i4 ( command(4:), i_temp, ierror, last )

      if ( ierror == 0 ) then

        ns = i_temp
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  NS has been set to ', ns

      end if

    else if ( s_eqi ( command(1:5), 'PRINT' ) ) then

      what = adjustl ( command(6:) )
      if ( len_trim ( what ) == 0 ) then
        what = 'ALL'
      end if

      if ( s_eqi ( what, 'N' ) .or. s_eqi ( what, 'ALL' ) ) then
        write ( *, '(a,i6)' ) '  The number of points, N = ', n
      end if

      if ( s_eqi ( what, 'M' ) .or. s_eqi ( what, 'ALL' ) ) then
        write ( *, '(a,i6)' ) '  The spatial dimension, M = ', m
      end if

      if ( s_eqi ( what, 'MAXIT' ) .or. s_eqi ( what, 'ALL' ) ) then
        write ( *, '(a,i6)' ) &
          '  The maximum number of iterations, MAXIT = ', maxit
      end if

      if ( s_eqi ( what, 'NS' ) .or. s_eqi ( what, 'ALL' ) ) then
        write ( *, '(a,i6)' ) '  I don''t know what this is, NS = ', ns
      end if

      if ( s_eqi ( what, 'DENSITY_FUNCTION' ) .or. s_eqi ( what, 'ALL' ) ) then
        write ( *, '(a,i6)' ) '  The density function, DENSITY_FUNCTION = ', &
          density_function
      end if

    else if ( s_eqi ( command, 'QUIT' ) ) then

      exit

    else if ( s_eqi ( command, 'UNIFORM_MAKE' ) ) then

      if ( m <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'UNIFORM_MAKE - Error!'
        write ( *, '(a)' ) '  The spatial dimension must be positive!'
        write ( *, '(a)' ) '  Enter the spatial dimension M with the command'
        write ( *, '(a)' ) '    "M = ???"'
        cycle
      end if

      if ( n <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'UNIFORM_MAKE - Error!'
        write ( *, '(a)' ) '  The number of points N must be positive!'
        write ( *, '(a)' ) '  Enter the number of points with the command'
        write ( *, '(a)' ) '    "N = ???"'
        cycle
      end if

      call uniform_make ( m, n, a, b, density_function, &
        uniform_points )

    else if ( s_eqi ( command, 'UNIFORM_READ' ) ) then

      call uniform_read ( m, n, uniform_points, uniform_file )

    else if ( s_eqi ( command, 'UNIFORM_WRITE' ) ) then

      if ( n <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Use the command UNIFORM_MAKE first!'
        cycle
      end if

      call uniform_write ( m, n, uniform_points, uniform_file )

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MESHLESS - Warning'
      write ( *, '(a)' ) '  Your command was not recognized.'

    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MESHLESS'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine basis_box ( m, n, center, radius, a, b )

!*****************************************************************************80
!
!! BASIS_BOX finds a box that contains all the basis functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) CENTER(M,N), the basis function centers.
!
!    Input, real ( kind = 8 ) RADIUS(N), the basis function radii.
!
!    Output, real ( kind = 8 ) A(M), B(M), the minimum and maximum
!    values in each dimension.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ), dimension ( m, n ) :: center
  integer ( kind = 4 ) i
  real ( kind = 8 ) radius(n)

  do i = 1, m
    a(i) = minval ( center(i,1:n) - radius(1:n) )
    b(i) = maxval ( center(i,1:n) + radius(1:n) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_BOX'
  write ( *, '(a)' ) '  Compute the bounding box.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Box limits for center points:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Dimension  Lower  Upper'
  write ( *, '(a)' ) ' '

  do i = 1, m
    write ( *, '(i4,2g14.6)' ) &
      i, minval ( center(i,1:n) ), maxval ( center(i,1:n) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Box limits for radius:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Lower  Upper'
  write ( *, '(a)' ) ' '
  write ( *, '(2g14.6)' ) minval ( radius(1:n) ), maxval ( radius(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Box limits for basis support:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Dimension  Lower  Upper'
  write ( *, '(a)' ) ' '

  do i = 1, m
    write ( *, '(i4,2g14.6)' ) i, a(i), b(i)
  end do

  return
end
subroutine basis_overlap ( m, n, center, radius )

!*****************************************************************************80
!
!! BASIS_OVERLAP counts the number of basis functions with overlapping support.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 January 2001
!
!  Author:
!
!    Lili Ju
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of basis functions.
!
!    Input, real ( kind = 8 ) CENTER(M,N), the basis function centers.
!
!    Input, real ( kind = 8 ) RADIUS(N), the basis function radii.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) m

  real ( kind = 8 ), dimension ( m, n ) :: center
  real ( kind = 8 ) dist
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mel
  integer ( kind = 4 ) nel
  real ( kind = 8 ), dimension ( n ) :: radius

  if ( n <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASIS_OVERLAP - Warning!'
    write ( *, '(a)' ) '  N <= 1, so no overlap possible.'
    return
  end if

  mel = ( n * ( n - 1 ) ) / 2
  nel = 0

  do i = 1, n
    do j = i+1, n

      dist = sqrt ( sum ( ( center(1:m,i) - center(1:m,j) )**2 ) )

      if ( dist <= ( radius(i) + radius(j) ) ) then
        nel = nel + 1
      end if

    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_OVERLAP:'
  write ( *, '(a)' ) '  Count the pairs of basis functions '
  write ( *, '(a)' ) '  with overlapping support.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum possible count = ', mel
  write ( *, '(a,i6)' ) '  Actual count =           ', nel
  write ( *, '(a,f10.4)' ) &
    '  Percentage =             ', &
    real ( 100 * nel, kind = 8 ) / real ( mel, kind = 8 )

  return
end
subroutine basis_plot ( m, n, center, radius )

!*****************************************************************************80
!
!! BASIS_PLOT plots the basis functions in the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of basis functions.
!
!    Input, real ( kind = 8 ) CENTER(M,N), the basis function centers.
!
!    Input, real ( kind = 8 ) RADIUS(N), the basis function radii.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) m

  real ( kind = 8 ) blue
  real ( kind = 8 ) center(m,n)
  character ( len = 80 ) file_name
  logical, parameter :: filled = .false.
  real ( kind = 8 ) green
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) r
  real ( kind = 8 ) radius(n)
  real ( kind = 8 ) red
  real ( kind = 8 ) rxmax
  real ( kind = 8 ) rxmin
  real ( kind = 8 ) rymax
  real ( kind = 8 ) rymin
  real ( kind = 8 ) x
  integer ( kind = 4 ), parameter :: x_ps_max = 576
  integer ( kind = 4 ), parameter :: x_ps_min = 36
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  integer ( kind = 4 ), parameter :: y_ps_max = 756
  integer ( kind = 4 ), parameter :: y_ps_min = 36
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  do i = 1, n
    if ( i == 1 ) then
      rxmin = center(1,i) - radius(i)
      rxmax = center(1,i) + radius(i)
      rymin = center(2,i) - radius(i)
      rymax = center(2,i) + radius(i)
    else
      rxmin = min ( rxmin, center(1,i) - radius(i) )
      rxmax = max ( rxmax, center(1,i) + radius(i) )
      rymin = min ( rymin, center(2,i) - radius(i) )
      rymax = max ( rymax, center(2,i) + radius(i) )
    end if
  end do

  write ( *, * ) 'RXMIN, RXMAX = ', rxmin, rxmax
  write ( *, * ) 'RYMIN, RYMAX = ', rymin, rymax

  call get_unit ( iunit )

  file_name = 'basis.eps'

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASIS_PLOT'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    return
  end if

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
    y_ps_max )

  xmin = rxmin - 0.1D+00 * ( rxmax - rxmin )
  xmax = rxmax + 0.1D+00 * ( rxmax - rxmin )
  ymin = rymin - 0.1D+00 * ( rymax - rymin )
  ymax = rymax + 0.1D+00 * ( rymax - rymin )

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Draw the outline of the region.
!
  red = 0.0D+00
  green = 0.0D+00
  blue = 1.0D+00

  call ps_color_line_set ( red, green, blue )

  call ps_moveto ( rxmin, rymin )
  call ps_lineto ( rxmax, rymin )
  call ps_lineto ( rxmax, rymax )
  call ps_lineto ( rxmin, rymax )
  call ps_lineto ( rxmin, rymin )
!
!  Draw a grid in the region.
!
  nx = 11
  ny = 11

  red = 0.5D+00
  green = 0.5D+00
  blue = 0.5D+00

  call ps_grid_cartesian ( rxmin, rxmax, nx, rymin, rymax, ny )
!
!  Draw the basis function support disks.
!
  red = 0.2D+00
  green = 0.2D+00
  blue = 0.2D+00

  call ps_color_line_set ( red, green, blue )

  red = 0.9D+00
  green = 0.6D+00
  blue = 0.6D+00

  call ps_color_fill_set ( red, green, blue )

  do i = 1, n

    x = center(1,i)
    y = center(2,i)
    r = radius(i)

    call ps_mark_disk ( x, y )

    if ( filled ) then
      call ps_circle_fill ( x, y, r )
    end if

    call ps_circle ( x, y, r )

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_PLOT'
  write ( *, '(a)' ) &
    '  Created a basis function plot file ' // trim ( file_name )

  return
end
subroutine basis_read ( m, n, center, radius, input_file )

!*****************************************************************************80
!
!! BASIS_READ reads basis data from a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Output, integer ( kind = 4 ) N, the number of points.
!
!    Output, real ( kind = 8 ) CENTER(M,N), the basis function centers.
!
!    Output, real ( kind = 8 ) RADIUS(N), the basis function radii.
!
!    Input, character ( len = * ) INPUT_FILE, the name of the file.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ), dimension ( m ) :: c_temp
  real ( kind = 8 ), dimension ( m, * ) :: center
  character ( len = * ) input_file
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) n
  real ( kind = 8 ) r_temp
  real ( kind = 8 ), dimension ( * ) :: radius

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_READ:'
  write ( *, '(a)' ) '  Reading basis data file: ' // trim ( input_file )

  open ( unit = 12, file = input_file, form = 'formatted', status = 'old' )

  n = 0

  do

    read ( 12, *, iostat = ios ) c_temp(1:m), r_temp

    if ( ios /= 0 ) then
      exit
    end if

    n = n + 1
    center(1:m,n) = c_temp(1:m)
    radius(n) = r_temp

  end do

  close ( unit = 12 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_READ:'
  write ( *, '(a,i6)' ) '  Number of data points read was ', n

  return
end
subroutine basis_write ( m, n, center, radius, output_file )

!*****************************************************************************80
!
!! BASIS_WRITE writes the basis data to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) CENTER(M,N), the basis function centers.
!
!    Input, real ( kind = 8 ) RADIUS(N), the basis function radii.
!
!    Input, character ( len = * ) OUTPUT_FILE, the name of the output file.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  character ( len = * ) output_file
  real ( kind = 8 ), dimension ( m, n ) :: center
  real ( kind = 8 ), dimension ( n ) :: radius

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_WRITE:'
  write ( *, '(a)' ) '  Write basis data file: ' // trim ( output_file )
  write ( *, '(a,i6)' ) '  Number of data points is ', n

  open ( unit = 12, file = output_file, form = 'formatted', status = 'replace' )

  do i = 1, n
    write ( 12, * ) center(1:m,i), radius(i)
  end do

  close ( unit = 12 )

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
subroutine center_plot ( m, n, center, a, b )

!*****************************************************************************80
!
!! CENTER_PLOT plots the center in the region.
!
!  Discussion:
!
!    My idea of plotting the basis functions is OK if there are just a
!    few, but quickly becomes useless.  Here, at least, I can plot just
!    the center points, and get an idea of the distribution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of basis functions.
!
!    Input, real ( kind = 8 ) CENTER(M,N), the basis function centers.
!
!    Input, real ( kind = 8 ) A(M), B(M), the minimum and maximum
!    values in each dimension.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) blue
  real ( kind = 8 ) center(m,n)
  character ( len = 80 ) file_name
  logical, parameter :: filled = .false.
  real ( kind = 8 ) green
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) marker_size
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) red
  real ( kind = 8 ) rxmax
  real ( kind = 8 ) rxmin
  real ( kind = 8 ) rymax
  real ( kind = 8 ) rymin
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  NOT CORRECT FOR 2 < M.
!
  rxmin = a(1)
  rxmax = b(1)
  rymin = a(2)
  rymax = b(2)

  iunit = 1
  file_name = 'center_plot.ps'

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CENTER_PLOT'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    return
  end if

  xmin = ( rxmin - 0.1D+00 * ( rxmax - rxmin ) )
  xmax = ( rxmax + 0.1D+00 * ( rxmax - rxmin ) )
  ymin = ( rymin - 0.1D+00 * ( rymax - rymin ) )
  ymax = ( rymax + 0.1D+00 * ( rymax - rymin ) )

  call ps_file_head ( file_name )

  call ps_page_head ( xmax, xmin, ymax, ymin )
!
!  Draw the outline of the region.
!
  red = 0.0D+00
  green = 0.0D+00
  blue = 1.0D+00

  call ps_color_line_set ( red, green, blue )

  call ps_moveto ( rxmin, rymin )
  call ps_lineto ( rxmax, rymin )
  call ps_lineto ( rxmax, rymax )
  call ps_lineto ( rxmin, rymax )
  call ps_lineto ( rxmin, rymin )
!
!  Draw a grid in the region.
!
  red = 0.2D+00
  green = 0.3D+00
  blue = 0.4D+00

  call ps_color_line_set ( red, green, blue )

  nx = 11
  ny = 11

  call ps_grid_cartesian ( rxmin, rxmax, nx, rymin, rymax, ny )

  marker_size = 3
  call ps_marker_size ( marker_size )
!
!  Draw the center points.
!
  x(1:n) = center(1,1:n)
  y(1:n) = center(2,1:n)

  call ps_mark_circles ( n, x, y )

  call ps_page_tail

  call ps_file_tail

  call ps_file_close ( iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CENTER_PLOT'
  write ( *, '(a)' ) &
    '  Created a basis function plot file ' // trim ( file_name )

  return
end
subroutine cvt_make ( alpha, beta, m, a, b, maxit, ns, density_function, &
  n, center )

!*****************************************************************************80
!
!! CVT_MAKE computes centroidal Voronoi tessellation generators.
!
!  Discussion:
!
!    The routine is given an initial set of points that define a 
!    Voronoi tessellation.  It computes new points that define a Voronoi
!    tessellation, and which have the property that they are the centroids
!    of their Voronoi regions.
!
!    No stopping criterion is used, for now, except to iterate to MAXIT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2001
!
!  Author:
!
!    Lili Ju
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA(2), BETA(2), coefficients used in the
!    generation program.  ALPHA(1) and ALPHA(2) must lie between 0 and 1
!    and sum to 1.  The same is true for BETA.  Common values are
!    ALPHA(1) = BETA(1) = 0,
!    ALPHA(2) = BETA(2) = 1.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, real ( kind = 8 ) A(M), B(M), the coordinates of the
!    two extreme corners of the box that defines the region.
!
!    Input, integer ( kind = 4 ) MAXIT, the maximum number of iterations.
!
!    Input, integer ( kind = 4 ) NS, the average number of sampling points 
!    per generator, on each step of the iteration.
!
!    Input, integer ( kind = 4 ) DENSITY_FUNCTION, specifies the density function.
!    1: d(x) = 1.0;
!    2: d(x) = exp ( - 4.0 * ( sum(x(1:n)**2) ) )
!    3: d(x) = exp ( - 3.0 * ( 1.0 - sum(x(1:n)**2) ) )
!
!    Input, integer ( kind = 4 ) N, the number of generators.
!
!    Input/output, real ( kind = 8 ) CENTER(M,N).
!    On input, initial points that generate the Voronoi regions.
!    On output, points that generate the Voronoi regions, which are
!    also the centroids of those regions.
!
!  Local variables:
!
!    Integer UPDATES(N), counts the number of times a generator 
!    has been updated.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(1:m)
  real ( kind = 8 ) alpha(2)
  real ( kind = 8 ) b(1:m)
  real ( kind = 8 ) beta(2)
  real ( kind = 8 ) center(m,n)
  real ( kind = 8 ) center2(m,n)
  integer ( kind = 4 ) count(n)
  integer ( kind = 4 ) density_function
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) iloop
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) maxit
  integer ( kind = 4 ) ns
  real ( kind = 8 ) p1(m)
  real ( kind = 8 ) p2(m)
  real ( kind = 8 ) s
  integer ( kind = 4 ) updates(n)
  real ( kind = 8 ) y(m)

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CVT_MAKE - Error!'
    write ( *, '(a)' ) '  The spatial dimension must be positive!'
    write ( *, '(a)' ) '  Enter a spatial dimension with the "M = " command.'
    return
  end if
!
!  UPDATES starts at 1.
!
  updates(1:n) = 1
!
!  Iterate MAXIT times.
!
  do iloop = 1, maxit

    center2(1:m,1:n) = 0.0D+00
    count(1:n) = 0

    do j = 1, n * ns
!
!  Generate a random sampling point Y.
!
      call random_generator ( m, a, b, density_function, y )
!
!  Find the nearest generator CENTER.  Its index is IC.
!
      call find_closest ( m, y, n, center, ic, s )
!
!  Add Y to the averaging data for CENTER(*,IC).
!
      center2(1:m,ic) = center2(1:m,ic) + y(1:m)
      count(ic) = count(ic) + 1

    end do
!
!  For each cell J, average the estimated centroid CENTER with the 
!  averaged hit points CENTER2.  (There are COUNT(J) of these new points).
!
    do j = 1, n

      if ( count(j) /= 0.0 ) then
!
!  This loop and the next could be combined, eliminating the
!  temporary variables P1 and P2.
!
        do k = 1, m
          
          p1(k) = ( alpha(1) * updates(j) + beta(1) ) * center(k,j)
          p2(k) = ( alpha(2) * updates(j) + beta(2) ) * center2(k,j) &
            / dble ( count(j) )
        end do

        do k = 1, m
          center(k,j) = ( p1(k) + p2(k) ) / dble ( updates(j) + 1 )
        end do

        updates(j) = updates(j) + 1

      end if

    end do

  end do

  return
end
subroutine cvt_read ( m, n, center, radius, input_file )

!*****************************************************************************80
!
!! CVT_READ reads the CVT data from a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, real ( kind = 8 ) CENTER(M,N), the centroidal points.
!
!    Input, character ( len = * ) INPUT_FILE, the name of the file.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  character ( len = * ) input_file
  real ( kind = 8 ), dimension ( m, n ) :: center
  real ( kind = 8 ), dimension ( n ) :: radius

  write ( *, '(a)' ) ' '
  write ( *, '(a)')  'CVT_READ:'
  write ( *, '(a)' ) '  Reading Centroidal Voronoi data file: ' &
    // trim ( input_file )
  write ( *, '(a,i6)' ) '  Number of data points is ', n

  open ( unit = 12, file = input_file, form = 'formatted', status = 'old' )

  do i = 1, n
    read ( 12, * ) center(1:m,i), radius(i)
  end do

  close ( unit = 12 )

  return
end
subroutine cvt_write ( m, n, center, radius, output_file )

!*****************************************************************************80
!
!! CVT_WRITE writes the CVT data to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) CENTER(M,N), the centroidal points.
!
!    Input, character ( len = * ) OUTPUT_FILE, the name of the output file.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  character ( len = * ) output_file
  real ( kind = 8 ), dimension ( m, n ) :: center
  real ( kind = 8 ), dimension ( n ) :: radius

  write ( *, '(a)' ) ' '
  write ( *, '(a)' )  'CVT_WRITE:'
  write ( *, '(a)' ) '  Write Centroidal Voronoi data file: ' &
    // trim ( output_file )
  write ( *, '(a,i6)' ) '  Number of data points is ', n

  open ( unit = 12, file = output_file, form = 'formatted', status = 'replace' )

  do i = 1, n
    write ( 12, * ) center(1:m,i), radius(i)
  end do

  close ( unit = 12 )

  return
end
function density ( m, x, density_function )

!*****************************************************************************80
!
!! DENSITY evaluates the density function.
!
!  Discussion:
!
!    The density function is used to control the generation of random points
!    that sample the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 January 2001
!
!  Author:
!
!    Lili Ju
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, real ( kind = 8 ) X(M), the location of the point.
!
!    Input, integer ( kind = 4 ) DENSITY_FUNCTION, chooses the density function.
!    1: d(x) = 1.0;
!    2: d(x) = exp ( - 4.0 * ( sum(x(1:n)**2) ) )
!    3: d(x) = exp ( - 3.0 * ( 1.0 - sum(x(1:n)**2) ) )
!
!    Output, real ( kind = 8 ) DENSITY, the value of the density 
!    function, which should be between 0 and 1 in the region.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) density
  integer ( kind = 4 ) density_function
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(m)
!
!  Here is an odd feature that sets density to zero if the point
!  is more than 1 away from the origin.  This is not what I documented,
!  and should probably be controllable by the user!
!
  if ( 1.0D+00 < sum ( x(1:m)**2 ) ) then
    density = 0.0D+00
    return
  end if

  if ( density_function == 1 ) then
    density = 1.0D+00
  else if ( density_function == 2 ) then
    density = exp ( - 4.0D+00 * sum ( x(1:m)**2 ) )
  else if ( density_function == 3 ) then
    density = exp ( - 3.0D+00 * ( 1.0D+00 - sum ( x(1:m)**2 ) ) )
  else
    density = 1.0D+00
  end if

  return
end
subroutine eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
  y_ps_max )

!*****************************************************************************80
!
!! EPS_FILE_HEAD writes header information to an encapsulated PostScript file.
!
!  Discussion:
!
!    The file should contain the description of only one page, but this
!    is not currently checked.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) X_PS_MIN, Y_PS_MIN, X_PS_MAX, Y_PS_MAX, the minimum 
!    and maximum X and Y values of the data, in PostScript units.  Any data
!    that lies outside this range will not show up properly.  A reasonable
!    set of values might be 0, 0, 612, 792, or, for a half inch margin,
!    36, 36, 576, 756.
!
  implicit none

  character ( len = 8 ) date
  character ( len = * ) file_name
  real ( kind = 8 ) line_blue
  real ( kind = 8 ) line_green
  real ( kind = 8 ) line_red
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  integer ( kind = 4 ) x_ps_max
  integer ( kind = 4 ) x_ps_min
  integer ( kind = 4 ) y_ps_max
  integer ( kind = 4 ) y_ps_min
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPS_FILE_HEAD - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 1 is required.'
    return
  end if
!
!  Initialization
!
  call ps_default ( )
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call date_and_time ( date )
!
!  Write the prolog.
!
  write ( unit, '(a)' )     '%!PS-Adobe-3.0 EPSF-3.0'
  write ( unit, '(a)' )     '%%Creator: ps_write.f90'
  write ( unit, '(a)' )     '%%Title: ' // trim ( file_name )
  write ( unit, '(a)' )     '%%CreationDate: '// trim ( date )
  write ( unit, '(a)' )     '%%Pages: 1'
  write ( unit, '(a,4i6)' ) '%%BoundingBox:', &
    x_ps_min, y_ps_min, x_ps_max, y_ps_max
  write ( unit, '(a)' )     '%%Document-Fonts: Times-Roman'
  write ( unit, '(a)' )     '%%LanguageLevel: 1'
  write ( unit, '(a)' )     '%%EndComments'
  write ( unit, '(a)' )     '%%BeginProlog'
  write ( unit, '(a)' )     '/inch {72 mul} def'
  write ( unit, '(a)' )     '%%EndProlog'
!
!  Set the font.
!
  write ( unit, '(a)' ) '/Times-Roman findfont'
  write ( unit, '(a)' ) '1.00 inch scalefont'
  write ( unit, '(a)' ) 'setfont'
!
!  Set the line color.
!
  line_red = 0.0D+00
  line_green = 0.0D+00
  line_blue = 0.0D+00

  call ps_color_line ( 'SET', line_red, line_green, line_blue )
!
!  Reset the state.
!
  state = 2

  call ps_setting_int ( 'SET', 'STATE', state )

  return
end
subroutine eps_file_tail ( )

!*****************************************************************************80
!
!! EPS_FILE_TAIL writes trailer information to an encapsulated PostScript file.
!
!  Discussion:
!
!    Looks like that penultimate 'end' line is not wanted, so I commented
!    it out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    None
!
  implicit none

  integer ( kind = 4 ) num_pages
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state == 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPS_FILE_TAIL - Warning!'
    write ( *, '(a)' ) '  A page was open.  It is being forced closed.'
    state = 2
    call ps_setting_int ( 'SET', 'STATE', state )
  end if

  if ( state /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPS_FILE_TAIL - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 2 is required.'
    return
  end if
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )
!
!  Retrieve the number of pages.
!
  call ps_setting_int ( 'GET', 'NUM_PAGES', num_pages )

  if ( 1 < num_pages ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPS_FILE_TAIL - Warning!'
    write ( *, '(a)' ) '  An encapsulated PostScript file describes ONE page.'
    write ( *, '(a,i9,a)' ) '  This file describes ', num_pages, ' pages.'
    write ( *, '(a)' ) '  It is not a legal EPS file.'
  end if
!
!  Write the epilog.
!
  write ( unit, '(a)' ) '%%Trailer'
! write ( unit, '(a)' ) 'end'
  write ( unit, '(a)' ) '%%EOF'
!
!  Zero out the number of pages.
!
  num_pages = 0

  call ps_setting_int ( 'SET', 'NUM_PAGES', num_pages )
!
!  Reset the state.
!
  state = 4

  call ps_setting_int ( 'SET', 'STATE', state )

  return
end
subroutine file_column_count ( file_name, ncolumn )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Each line of the file is presumed to consist of NCOLUMN words, separated
!    by spaces.
!
!    The routine reads the first line and counts the number of words.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) NCOLUMN, the number of columns assumed to be in the file.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = 256 ) line
  integer ( kind = 4 ) ncolumn
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    ncolumn = - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_LINE_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    ' // trim ( file_name )
    return
  end if
!
!  Read one line.
!
  read ( iunit, '(a)', iostat = ios ) line

  if ( ios /= 0 ) then
    ncolumn = -2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_LINE_COUNT - Fatal error!'
    write ( *, '(a)' ) '  The initial line of the file could not be read.'
    write ( *, '(a)' ) '    ' // trim ( file_name )
    return
  end if

  close ( unit = iunit )

  call s_word_count ( line, ncolumn )

  return
end
subroutine file_line_count ( file_name, nline )

!*****************************************************************************80
!
!! FILE_LINE_COUNT counts the number of lines in a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) NLINE, the number of lines found in the file.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) nline

  nline = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    nline = - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_LINE_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    ' // trim ( file_name )
    return
  end if
!
!  Count the lines.
!
  do

    read ( iunit, '(a)', iostat = ios )

    if ( ios /= 0 ) then
      exit
    end if

    nline = nline + 1

  end do

  close ( unit = iunit )

  return
end
subroutine find_closest ( m, y, n, center, ic, s )

!*****************************************************************************80
!
!! FIND_CLOSEST finds the center point CENTER(1:M,IC) closest to Y(1:M).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2001
!
!  Author:
!
!    Lili Ju
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, real ( kind = 8 ) Y(M), the point to be checked.
!
!    Input, integer ( kind = 4 ) N, the number of center points.
!
!    Input, real ( kind = 8 ) CENTER(M,N), the center points.
!
!    Output, integer ( kind = 4 ) IC, the index of the nearest center point.
!
!    Output, real ( kind = 8 ) S, the distance.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) m

  real ( kind = 8 ) center(m,n)
  real ( kind = 8 ) dist_sq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  real ( kind = 8 ) s
  real ( kind = 8 ) y(m)

  s = huge ( s )

  do i = 1, n

    dist_sq = sum ( ( center(1:m,i) - y(1:m) )**2 )

    if ( dist_sq < s ) then
      s = dist_sq
      ic = i
    end if

  end do

  s = sqrt ( s )

  return
end
subroutine find_re ( m, pt, r, n, center, np, dist, ord, ierr )

!*****************************************************************************80
!
!! FIND_RE finds the number of center points near a given point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2001
!
!  Author:
!
!    Lili Ju
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, real ( kind = 8 ) PT(M), the base point.
!
!    Input, real ( kind = 8 ) R, the maximum distance from PT 
!    to be considered.
!
!    Input, integer ( kind = 4 ) N, the number of center points.
!
!    Input, real ( kind = 8 ) CENTER(M,N), the center points.
!
!    Output, integer ( kind = 4 ) NP, the number of neighboring center points found.
!
!    Output, real ( kind = 8 ) DIST(1:NP), lists the Euclidean distance 
!    between PT and each of the neighboring center points found.
!
!    Output, integer ( kind = 4 ) ORD(1:NP), lists the indices in CENTER of the points found.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    -1, if no points at all were found.
!    0, if at least one point was found.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) m

  real ( kind = 8 ) center(m,n)
  real ( kind = 8 ) dist(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) np
  integer ( kind = 4 ) ord(n)
  real ( kind = 8 ) pt(m)
  real ( kind = 8 ) r
!
!  Search for a center point CENTER(I) which is no more than R away in X or Y,
!  and which is not identical to PT.  Count the number of such points found,
!  their distance, and indexes.
!
  np = 0

  do i = 1, n
!
!  We do not allow the case where a point in CENTER is exactly equal to PT.
!
    if ( all ( pt(1:m) == center(1:m,i) ) ) then
      cycle
    end if
!
!  A center point is a suitable neighbor if its L-Infinity distance
!  from PT is no more than R.  But for some reason, we then record
!  its L2 distance...Is this an inconsistency?
!
    if ( all ( abs ( pt(1:m) - center(1:m,i) ) <= r ) ) then
      np = np + 1
      dist(np) = sqrt ( sum ( ( pt(1:m) - center(1:m,i) )**2 ) )
      ord(np) = i
    end if

  end do

  if ( np == 0 ) then
    ierr = -1
  else
    ierr = 0
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
!    Output, integer ( kind = 4 ) IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 ) then

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
subroutine halton_generate ( m, n, skip, base, r )

!*****************************************************************************80
!
!! HALTON_GENERATE generates a Halton dataset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points to generate.
!
!    Input, integer ( kind = 4 ) SKIP, the number of initial points to skip.
!
!    Output, integer ( kind = 4 ) BASE(M), the bases used for the sequence.
!
!    Output, real ( kind = 8 ) R(M,N), the points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ), dimension ( m ) :: base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) prime
  real ( kind = 8 ), dimension ( m, n ) :: r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) skip

  do i = 1, m
    base(i) = prime(i)
  end do

  do j = 1, n
    seed = skip + j - 1
    call i4_to_halton ( seed, base, m, r(1:m,j) )
  end do

  return
end
subroutine halton_read ( m, n, points, input_file )

!*****************************************************************************80
!
!! HALTON_READ reads Halton data from a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, real ( kind = 8 ) POINTS(M,N), the Halton points.
!
!    Input, character ( len = * ) INPUT_FILE, the name of the file.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  character ( len = * ) input_file
  real ( kind = 8 ), dimension ( m, n ) :: points

  write ( *, '(a)' ) ' '
  write ( *, '(a)' )  'HALTON_READ:'
  write ( *, '(a)' ) '  Reading Halton data file: ' // trim ( input_file )
  write ( *, '(a,i6)' ) '  Number of data points is ', n

  open ( unit = 12, file = input_file, form = 'formatted', status = 'old' )

  do i = 1, n
    read ( 12, * ) points(1:m,i)
  end do

  close ( unit = 12 )

  return
end
subroutine halton_write ( m, n, skip, base, r, file_out_name )

!*****************************************************************************80
!
!! HALTON_WRITE writes a Halton dataset to a file.
!
!  Discussion:
!
!    The initial lines of the file are comments, which begin with a
!    '#' character.
!
!    Thereafter, each line of the file contains the M-dimensional
!    components of the SKIP+I-1 entry of the Halton sequence.
!
!    For the Halton sequence, the value of SKIP is the same
!    as the value of SEED used to generate the first point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of (successive) points.
!
!    Input, integer ( kind = 4 ) SKIP, the number of skipped points.
!
!    Input, integer ( kind = 4 ) BASE(M), the bases used for the sequence.
!
!    Input, real ( kind = 8 ) R(M,N), the points.
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of
!    the output file.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) base(m)
  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mhi
  integer ( kind = 4 ) mlo
  real ( kind = 8 ) r(m,n)
  integer ( kind = 4 ) skip
  character ( len = 40 ) string

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HALTON_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if

  call timestring ( string )

  write ( file_out_unit, '(a)'      ) '#  ' // trim ( file_out_name )
  write ( file_out_unit, '(a)'      ) '#  created by MESHLESS.F90'
  write ( file_out_unit, '(a)'      ) '#'
  write ( file_out_unit, '(a)'      ) '#  File generated on ' // trim ( string )
  write ( file_out_unit, '(a)'      ) '#'
  write ( file_out_unit, '(a,i6)'   ) '#  Spatial dimension M = ', m
  write ( file_out_unit, '(a,i6)'   ) '#  Number of points N = ', n
  do mlo = 1, m, 10
    mhi = min ( mlo + 9, m )
    if ( mlo == 1 ) then
      write ( file_out_unit, '(a,20i6)' ) '#  Bases: ', base(mlo:mhi)
    else
      write ( file_out_unit, '(a,20i6)' ) '#         ', base(mlo:mhi)
    end if
  end do
  write ( file_out_unit, '(a,i6)'   ) '#  Initial values skipped = ', skip
  write ( file_out_unit, '(a)'      ) '#'

  write ( string, '(a,i3,a)' ) '(', m, 'f10.6)'

  do j = 1, n
    write ( file_out_unit, string ) r(1:m,j)
  end do

  close ( unit = file_out_unit )

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two integer values.
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
subroutine i4_to_halton ( seed, base, ndim, r )

!*****************************************************************************80
!
!! I4_TO_HALTON computes an element of a Halton sequence.
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
!  Reference:
!
!    John Halton,
!    On the efficiency of certain quasi-random sequences of points
!    in evaluating multi-dimensional integrals,
!    Numerische Mathematik,
!    Volume 2, pages 84-90, 1960.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED, the index of the desired element.
!    Only the absolute value of SEED is considered.  SEED = 0 is allowed,
!    and returns R = 0.
!
!    Input, integer ( kind = 4 ) BASE(NDIM), the Halton bases, which should be
!    distinct prime numbers.  This routine only checks that each base
!    is greater than 1.
!
!    Input, integer ( kind = 4 ) NDIM, the dimension of the sequence.
!
!    Output, real ( kind = 8 ) R(NDIM), the SEED-th element of the
!    Halton sequence for the given bases.
!
  implicit none

  integer ( kind = 4 ) ndim

  integer ( kind = 4 ) base(ndim)
  real ( kind = 8 ) base_inv(ndim)
  integer ( kind = 4 ) digit(ndim)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r(ndim)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed2(ndim)

  seed2(1:ndim) = abs ( seed )

  r(1:ndim) = 0.0D+00

  if ( any ( base(1:ndim) <= 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_HALTON - Fatal error!'
    write ( *, '(a)' ) '  An input base BASE is <= 1!'
    do i = 1, ndim
      write ( *, '(i6,i6)' ) i, base(i)
    end do
    stop
  end if

  base_inv(1:ndim) = 1.0D+00 / real ( base(1:ndim), kind = 8 )

  do while ( any ( seed2(1:ndim) /= 0 ) )
    digit(1:ndim) = mod ( seed2(1:ndim), base(1:ndim) )
    r(1:ndim) = r(1:ndim) + real ( digit(1:ndim), kind = 8 ) * base_inv(1:ndim)
    base_inv(1:ndim) = base_inv(1:ndim) / real ( base(1:ndim), kind = 8 )
    seed2(1:ndim) = seed2(1:ndim) / base(1:ndim)
  end do

  return
end
function point_inside_box_2d ( x1, y1, x2, y2, x, y )

!*****************************************************************************80
!
!! POINT_INSIDE_BOX_2D determines if a point is inside a box in 2D.
!
!  Discussion:
!
!    A "box" is defined by its "left down" corner and its
!    "right up" corner, and all the points between.  It is
!    assumed that the sides of the box align with coordinate directions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, the two corners of the box.
!
!    Input, real ( kind = 8 ) X, Y, the point to be checked.
!
!    Output, logical POINT_INSIDE_BOX_2D, is .TRUE. if (X,Y) is inside the
!    box, or on its boundary, and .FALSE. otherwise.
!
  implicit none

  logical point_inside_box_2d
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  if ( x1 <= x .and. x <= x2 .and. &
       y1 <= y .and. y <= y2 ) then
    point_inside_box_2d = .true.
  else
    point_inside_box_2d = .false.
  end if

  return
end
function prime ( n )

!*****************************************************************************80
!
!! PRIME returns any of the first PRIME_MAX prime numbers.
!
!  Discussion:
!
!    PRIME_MAX is 1500, and the largest prime stored is 12553.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz and Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964, pages 870-873.
!
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, pages 95-98.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the index of the desired prime number.
!    N = -1 returns PRIME_MAX, the index of the largest prime available.
!    N = 0 is legal, returning PRIME = 1.
!    It should generally be true that 0 <= N <= PRIME_MAX.
!
!    Output, integer ( kind = 4 ) PRIME, the N-th prime.  If N is out of range, PRIME
!    is returned as 0.
!
  implicit none

  integer ( kind = 4 ), parameter :: prime_max = 1500

  integer ( kind = 4 ), save :: icall = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save, dimension ( prime_max ) :: npvec
  integer ( kind = 4 ) prime

  if ( icall == 0 ) then

    icall = 1

    npvec(1:100) = (/ &
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29, &
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71, &
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113, &
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173, &
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229, &
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281, &
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349, &
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409, &
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463, &
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541 /)

    npvec(101:200) = (/ &
      547,  557,  563,  569,  571,  577,  587,  593,  599,  601, &
      607,  613,  617,  619,  631,  641,  643,  647,  653,  659, &
      661,  673,  677,  683,  691,  701,  709,  719,  727,  733, &
      739,  743,  751,  757,  761,  769,  773,  787,  797,  809, &
      811,  821,  823,  827,  829,  839,  853,  857,  859,  863, &
      877,  881,  883,  887,  907,  911,  919,  929,  937,  941, &
      947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013, &
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, &
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, &
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 /)

    npvec(201:300) = (/ &
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, &
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, &
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, &
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, &
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, &
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, &
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, &
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, &
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987 /)

    npvec(301:400) = (/ &
     1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, &
     2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, &
     2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, &
     2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, &
     2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, &
     2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, &
     2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, &
     2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, &
     2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, &
     2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741 /)

    npvec(401:500) = (/ &
     2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, &
     2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, &
     2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, &
     3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, &
     3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, &
     3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, &
     3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, &
     3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, &
     3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, &
     3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571 /)

    npvec(501:600) = (/ &
     3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, &
     3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, &
     3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, &
     3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, &
     3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, &
     4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, &
     4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, &
     4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, &
     4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, &
     4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409 /)

    npvec(601:700) = (/ &
     4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, &
     4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, &
     4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, &
     4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, &
     4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, &
     4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, &
     4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, &
     5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, &
     5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, &
     5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279 /)

    npvec(701:800) = (/ &
     5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, &
     5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, &
     5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, &
     5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, &
     5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, &
     5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, &
     5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, &
     5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, &
     5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, &
     6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133 /)

    npvec(801:900) = (/ &
     6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, &
     6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, &
     6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, &
     6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, &
     6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, &
     6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, &
     6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, &
     6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, &
     6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, &
     6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997 /)

    npvec(901:1000) = (/ &
     7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, &
     7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, &
     7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, &
     7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, &
     7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, &
     7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, &
     7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, &
     7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, &
     7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, &
     7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /)

    npvec(1001:1100) = (/ &
     7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, &
     8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, &
     8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, &
     8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, &
     8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, &
     8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, &
     8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, &
     8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, &
     8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, &
     8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831 /)

    npvec(1101:1200) = (/ &
     8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, &
     8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, &
     9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, &
     9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, &
     9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, &
     9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, &
     9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, &
     9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, &
     9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, &
     9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733 /)

    npvec(1201:1300) = (/ &
     9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, &
     9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, &
     9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007, &
    10009,10037,10039,10061,10067,10069,10079,10091,10093,10099, &
    10103,10111,10133,10139,10141,10151,10159,10163,10169,10177, &
    10181,10193,10211,10223,10243,10247,10253,10259,10267,10271, &
    10273,10289,10301,10303,10313,10321,10331,10333,10337,10343, &
    10357,10369,10391,10399,10427,10429,10433,10453,10457,10459, &
    10463,10477,10487,10499,10501,10513,10529,10531,10559,10567, &
    10589,10597,10601,10607,10613,10627,10631,10639,10651,10657 /)

    npvec(1301:1400) = (/ &
    10663,10667,10687,10691,10709,10711,10723,10729,10733,10739, &
    10753,10771,10781,10789,10799,10831,10837,10847,10853,10859, &
    10861,10867,10883,10889,10891,10903,10909,19037,10939,10949, &
    10957,10973,10979,10987,10993,11003,11027,11047,11057,11059, &
    11069,11071,11083,11087,11093,11113,11117,11119,11131,11149, &
    11159,11161,11171,11173,11177,11197,11213,11239,11243,11251, &
    11257,11261,11273,11279,11287,11299,11311,11317,11321,11329, &
    11351,11353,11369,11383,11393,11399,11411,11423,11437,11443, &
    11447,11467,11471,11483,11489,11491,11497,11503,11519,11527, &
    11549,11551,11579,11587,11593,11597,11617,11621,11633,11657 /)

    npvec(1401:1500) = (/ &
    11677,11681,11689,11699,11701,11717,11719,11731,11743,11777, &
    11779,11783,11789,11801,11807,11813,11821,11827,11831,11833, &
    11839,11863,11867,11887,11897,11903,11909,11923,11927,11933, &
    11939,11941,11953,11959,11969,11971,11981,11987,12007,12011, &
    12037,12041,12043,12049,12071,12073,12097,12101,12107,12109, &
    12113,12119,12143,12149,12157,12161,12163,12197,12203,12211, &
    12227,12239,12241,12251,12253,12263,12269,12277,12281,12289, &
    12301,12323,12329,12343,12347,12373,12377,12379,12391,12401, &
    12409,12413,12421,12433,12437,12451,12457,12473,12479,12487, &
    12491,12497,12503,12511,12517,12527,12539,12541,12547,12553 /)

  end if

  if ( n == -1 ) then
    prime = prime_max
  else if ( n == 0 ) then
    prime = 1
  else if ( n <= prime_max ) then
    prime = npvec(n)
  else
    prime = 0
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRIME - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal prime index N = ', n
    write ( *, '(a,i6)' ) '  N must be between 0 and PRIME_MAX =', prime_max
    stop
  end if

  return
end
subroutine ps_circle ( x0, y0, r )

!*****************************************************************************80
!
!! PS_CIRCLE draws a circle.
!
!  Discussion:
!
!    As a side effect, the current point is set to the center of the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center
!    of the circle.
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ), parameter :: angle_max = 360
  integer ( kind = 4 ), parameter :: angle_min = 0
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) pr
  integer ( kind = 4 ) pxcen
  integer ( kind = 4 ) pycen
  real ( kind = 8 ) r
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  real ( kind = 8 ) x0
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y0
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_CIRCLE - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )
  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'YMIN', ymin )

  pxcen = plotxmin2 + nint ( alpha * ( x0 - xmin ) )
  pycen = plotymin2 + nint ( alpha * ( y0 - ymin ) )
  pr = nint ( alpha * r )

  write ( unit, '(a)' ) 'newpath'
  write ( unit, '(5i6,a)' ) pxcen, pycen, pr, angle_min, angle_max, ' arc'
!
!  Draw the circle.
!
  write ( unit, '(a)' ) 'closepath stroke'

  call ps_setting_real ( 'SET', 'XCUR', x0 )
  call ps_setting_real ( 'SET', 'YCUR', y0 )

  return
end
subroutine ps_color_fill_set ( r, g, b )

!*****************************************************************************80
!
!! PS_COLOR_FILL_SET sets the fill color.
!
!  Discussion:
!
!    By calling this routine, you guarantee that a check will be made
!    of the current fill color.  If the current and new fill colors are
!    the same, then we skip the extraneous action of setting the color.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB values for the new fill color.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b_old
  real ( kind = 8 ) g
  real ( kind = 8 ) g_old
  real ( kind = 8 ) r
  real ( kind = 8 ) r_old
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Check the state.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_COLOR_FILL_SET - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  A PostScript state of 1 or more is required.'
    return
  end if
!
!  Get the current colors.
!
  call ps_setting_real ( 'GET', 'FILL_RED', r_old )
  call ps_setting_real ( 'GET', 'FILL_GREEN', g_old )
  call ps_setting_real ( 'GET', 'FILL_BLUE', b_old )
!
!  If any color has changed, we need to reset them.
!
  if ( r_old /= r .or. g_old /= g .or. b_old /= b ) then

    call ps_setting_int ( 'GET', 'UNIT', unit )

    call ps_comment ( 'Set RGB line color.' )

    write ( unit, '(3f7.4,a)' ) r, g, b, ' setrgbcolor'

    call ps_setting_real ( 'SET', 'FILL_RED', r )
    call ps_setting_real ( 'SET', 'FILL_GREEN', g )
    call ps_setting_real ( 'SET', 'FILL_BLUE', b )

  end if

  return
end
subroutine ps_color_line ( action, r, g, b )

!*****************************************************************************80
!
!! PS_COLOR_LINE handles the line color.
!
!  Discussion:
!
!    By calling this routine, you can temporarily set the line color,
!    draw some lines, and then restore it to whatever it was.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the desired action.
!    'SET', set the line color to RGB.
!    'GET', set RGB to the current line color.
!    'PUSH', push a value onto the RGB stack.
!    'POP', pop the RGB stack.
!
!    Input, real ( kind = 8 ) R, G, B, the RGB values for the new line color.
!
  implicit none

  integer ( kind = 4 ), parameter :: nstack = 10

  character ( len = * ) action
  real ( kind = 8 ) b
  real ( kind = 8 ) b_old
  real ( kind = 8 ) b_stack(nstack)
  real ( kind = 8 ) g
  real ( kind = 8 ) g_old
  real ( kind = 8 ) g_stack(nstack)
  integer ( kind = 4 ), save :: istack = 0
  real ( kind = 8 ) r
  real ( kind = 8 ) r_old
  real ( kind = 8 ) r_stack(nstack)
  logical s_eqi

  if ( s_eqi ( action, 'SET' ) ) then

    call ps_color_line_set ( r, g, b )

  else if ( s_eqi ( action, 'GET' ) ) then

    call ps_setting_real ( 'GET', 'LINE_RED', r )
    call ps_setting_real ( 'GET', 'LINE_GREEN', g )
    call ps_setting_real ( 'GET', 'LINE_BLUE', b )

  else if ( s_eqi ( action, 'POP' ) ) then

    if ( 0 < istack ) then
      r = r_stack(istack)
      g = g_stack(istack)
      b = b_stack(istack)
      istack = istack - 1
    end if

    call ps_color_line_set ( r, g, b )

  else if ( s_eqi ( action, 'PUSH' ) ) then

    call ps_setting_real ( 'GET', 'LINE_RED', r_old )
    call ps_setting_real ( 'GET', 'LINE_GREEN', g_old )
    call ps_setting_real ( 'GET', 'LINE_BLUE', b_old )

    if ( istack <= nstack ) then
      istack = istack + 1
      r_stack(istack) = r_old
      g_stack(istack) = g_old
      b_stack(istack) = b_old
    end if

    call ps_color_line_set ( r, g, b )

  end if

  return
end
subroutine ps_color_line_set ( r, g, b )

!*****************************************************************************80
!
!! PS_COLOR_LINE_SET sets the line color.
!
!  Discussion:
!
!    By calling this routine, you guarantee that a check will be made
!    of the current line color.  If the current and new line colors are
!    the same, then we skip the extraneous action of setting the color.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB values for the new line color.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b_old
  real ( kind = 8 ) g
  real ( kind = 8 ) g_old
  real ( kind = 8 ) r
  real ( kind = 8 ) r_old
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Check the state.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_COLOR_LINE_SET - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  A PostScript state of at least 1 is required.'
    return
  end if
!
!  Get the current colors.
!
  call ps_setting_real ( 'GET', 'LINE_RED', r_old )
  call ps_setting_real ( 'GET', 'LINE_GREEN', g_old )
  call ps_setting_real ( 'GET', 'LINE_BLUE', b_old )
!
!  If any color has changed, we need to reset them.
!
  if ( r_old /= r .or. g_old /= g .or. b_old /= b ) then

    call ps_setting_int ( 'GET', 'UNIT', unit )

    call ps_comment ( 'Set RGB line color.' )

    write ( unit, '(3f7.4,a)' ) r, g, b, ' setrgbcolor'

    call ps_setting_real ( 'SET', 'LINE_RED', r )
    call ps_setting_real ( 'SET', 'LINE_GREEN', g )
    call ps_setting_real ( 'SET', 'LINE_BLUE', b )

  end if

  return
end
subroutine ps_comment ( string )

!*****************************************************************************80
!
!! PS_COMMENT inserts a comment into the PostScript file.
!
!  Discussion:
!
!    A comment begins with a percent sign in column 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the comment.
!
  implicit none

  character ( len = * ) string
  integer ( kind = 4 ) unit
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )
!
!  Write the comment.
!
  if ( len_trim ( string ) == 0 ) then
    write ( unit, '(a)' ) '%'
  else
    write ( unit, '(a)' ) '%'
    write ( unit, '(a2,a)' ) '% ', trim ( string )
    write ( unit, '(a)' ) '%'
  end if

  return
end
subroutine ps_default ( )

!*****************************************************************************80
!
!! PS_DEFAULT sets the internal settings to their default values
!
!  Discussion:
!
!    Certain variables are not reset, including the number of pages,
!    the unit number, the internal state, and variables relating to
!    the size and shape of the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    None
!
  implicit none

  real ( kind = 8 ) fill_blue
  real ( kind = 8 ) fill_green
  real ( kind = 8 ) fill_red
  real ( kind = 8 ) font_size
  real ( kind = 8 ) line_blue
  real ( kind = 8 ) line_green
  real ( kind = 8 ) line_red
  integer ( kind = 4 ) line_width
  integer ( kind = 4 ) marker_size

  line_width = 1
  marker_size = 5

  call ps_setting_int ( 'SET', 'LINE_WIDTH', line_width )
  call ps_setting_int ( 'SET', 'MARKER_SIZE', marker_size )

  fill_blue = 0.7D+00
  fill_green = 0.7D+00
  fill_red = 0.7D+00
  font_size = 0.1D+00
  line_blue = 0.0D+00
  line_green = 0.0D+00
  line_red = 0.0D+00

  call ps_setting_real ( 'SET', 'FILL_BLUE', fill_blue )
  call ps_setting_real ( 'SET', 'FILL_GREEN', fill_green )
  call ps_setting_real ( 'SET', 'FILL_RED', fill_red )
  call ps_setting_real ( 'SET', 'FONT_SIZE', font_size )
  call ps_setting_real ( 'SET', 'LINE_BLUE', line_blue )
  call ps_setting_real ( 'SET', 'LINE_GREEN', line_green )
  call ps_setting_real ( 'SET', 'LINE_RED', line_red )

  return
end
subroutine ps_file_close ( unit )

!*****************************************************************************80
!
!! PS_FILE_CLOSE closes a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) UNIT, the FORTRAN unit to which output was written.
!
  implicit none

  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state < 1 .or. 4 < state ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_CLOSE - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 1, 2, 3 or 4 is required.'
    return
  end if

  close ( unit = unit )

  state = 0
  call ps_setting_int ( 'SET', 'STATE', state )

  unit = 0
  call ps_setting_int ( 'SET', 'UNIT', unit )

  return
end
subroutine ps_file_head ( file_name )

!*****************************************************************************80
!
!! PS_FILE_HEAD writes header information to a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
  implicit none

  character ( len = 8 ) date
  character ( len = * ) file_name
  real ( kind = 8 ) line_blue
  real ( kind = 8 ) line_green
  real ( kind = 8 ) line_red
  integer ( kind = 4 ) margin
  integer ( kind = 4 ) pagexmax
  integer ( kind = 4 ) pagexmin
  integer ( kind = 4 ) pageymax
  integer ( kind = 4 ) pageymin
  integer ( kind = 4 ) plotxmax
  integer ( kind = 4 ) plotxmin
  integer ( kind = 4 ) plotymax
  integer ( kind = 4 ) plotymin
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_HEAD - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 1 is required.'
    return
  end if
!
!  Initialization
!
  call ps_default ( )
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
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call date_and_time ( date )
!
!  Write the prolog.
!
  write ( unit, '(a)' )     '%!PS-Adobe-1.0'
  write ( unit, '(a)' )     '%%Creator: ps_write.f90'
  write ( unit, '(a)' )     '%%Title: ' // trim ( file_name )
  write ( unit, '(a)' )     '%%CreationDate: ' // trim ( date )
  write ( unit, '(a)' )     '%%Pages: (atend)'
  write ( unit, '(a,4i6)' ) '%%BoundingBox:', plotxmin, plotymin, plotxmax, &
    plotymax
  write ( unit, '(a)' )     '%%Document-Fonts: Times-Roman'
  write ( unit, '(a)' )     '%%LanguageLevel: 1'
  write ( unit, '(a)' )     '%%EndComments'
  write ( unit, '(a)' )     '%%BeginProlog'
  write ( unit, '(a)' )     '/inch {72 mul} def'
  write ( unit, '(a)' )     '%%EndProlog'
!
!  Set the font.
!
  call ps_comment ( 'Set the font:' )

  write ( unit, '(a)' ) '/Times-Roman findfont'
  write ( unit, '(a)' ) '1.00 inch scalefont'
  write ( unit, '(a)' ) 'setfont'
!
!  Set the line color.
!
  line_red = 0.0D+00
  line_green = 0.0D+00
  line_blue = 0.0D+00

  call ps_color_line ( 'SET', line_red, line_green, line_blue )
!
!  Reset the state.
!
  state = 2

  call ps_setting_int ( 'SET', 'STATE', state )

  return
end
subroutine ps_file_open ( file_name, unit, ierror )

!*****************************************************************************80
!
!! PS_FILE_OPEN opens a new version of a PostScript file with a given name.
!
!  Discussion:
!
!    If a file of the given name already exists, it is deleted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) UNIT, the FORTRAN unit to which output should
!    be written.
!
!    Input, character ( len = 80 ) FILE_NAME, the name of the output file.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    nonzero, the file could not be created.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_OPEN - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 0 is required.'
    write ( *, '(a)' ) '  Call PS_FILE_CLOSE first!'
    return
  end if

  ierror = 0
!
!  Now create a new empty file of the given name.
!
  open ( unit = unit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    ierror = ios
    return
  end if

  state = 1
  call ps_setting_int ( 'SET', 'STATE', state )

  call ps_setting_int ( 'SET', 'UNIT', unit )

  return
end
subroutine ps_file_tail ( )

!*****************************************************************************80
!
!! PS_FILE_TAIL writes trailer information to a PostScript file.
!
!  Discussion:
!
!    Looks like that penultimate 'end' line is not wanted, so
!    I commented it out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    None
!
  implicit none

  integer ( kind = 4 ) num_pages
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state == 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_TAIL - Warning!'
    write ( *, '(a)' ) '  A page was open.  It is being forced closed.'
    state = 2
    call ps_setting_int ( 'SET', 'STATE', state )
  end if

  if ( state /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_TAIL - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 2 is required.'
    return
  end if
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )
!
!  Retrieve the number of pages.
!
  call ps_setting_int ( 'GET', 'NUM_PAGES', num_pages )
!
!  Write the epilog.
!
  write ( unit, '(a)' ) '%%Trailer'
  write ( unit, '(a,i6)' ) '%%Pages: ', num_pages
! write ( unit, '(a)' ) 'end'
  write ( unit, '(a)' ) '%%EOF'
!
!  Zero out the number of pages.
!
  num_pages = 0

  call ps_setting_int ( 'SET', 'NUM_PAGES', num_pages )
!
!  Reset the state.
!
  state = 4

  call ps_setting_int ( 'SET', 'STATE', state )

  return
end
subroutine ps_grid_cartesian ( xmin, xmax, nx, ymin, ymax, ny )

!*****************************************************************************80
!
!! PS_GRID_CARTESIAN draws a cartesian grid.
!
!  Discussion:
!
!    The current point is not modified by this call.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XMIN, XMAX, the minimum and maximum values
!    at which X grid lines should be drawn.
!
!    Input, integer ( kind = 4 ) NX, the number of X grid lines.
!    If NX is not positive, no X grid lines are drawn.
!    If NX is 1, a single grid line is drawn midway.
!
!    Input, real ( kind = 8 ) YMIN, YMAX, the minimum and maximum values 
!    at which Y grid lines should be drawn.
!
!    Input, integer ( kind = 4 ) NY, the number of Y grid lines.
!    If NY is not positive, no Y grid lines are drawn.
!    If NY is 1, a single grid line is drawn midway.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px
  integer ( kind = 4 ) py
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xmin2
  real ( kind = 8 ) y
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) ymin2
!
!  At least one of NX and NY must be positive.
!
  if ( nx < 1 .and. ny < 1 ) then
    return
  end if
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_GRID_CARTESIAN - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )
!
!  Get settings.
!
  alpha = 0.0D+00
  xmin2 = 0.0D+00
  ymin2 = 0.0D+00

  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )

  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin2 )
  call ps_setting_real ( 'GET', 'YMIN', ymin2 )
!
!  Draw the vertical (X) grid lines.
!
  do i = 1, nx

    if ( 1 < nx ) then

      x = ( real ( nx - i,     kind = 8 ) * xmin   &
          + real (      i - 1, kind = 8 ) * xmax ) &
          / real ( nx     - 1, kind = 8 )

    else if ( nx == 1 ) then

      x = 0.5D+00 * ( xmin + xmax )

    end if

    px = plotxmin2 + nint ( alpha * ( x - xmin2 ) )

    write ( unit, '(a)' ) 'newpath'

    py = plotymin2 + nint ( alpha * ( ymin - ymin2 ) )
    write ( unit, '(2i6,a)' ) px, py, ' moveto'

    py = plotymin2 + nint ( alpha * ( ymax - ymin2 ) )
    write ( unit, '(2i6,a)' ) px, py, ' lineto'

    write ( unit, '(a)' ) 'stroke'

  end do
!
!  Draw the horizontal (Y) grid lines.
!
  do i = 1, ny

    if ( 1 < ny ) then

      y = ( real ( ny - i,     kind = 8 ) * ymin   &
          + real (      i - 1, kind = 8 ) * ymax ) &
          / real ( ny     - 1, kind = 8 )

    else if ( ny == 1 ) then

      y = 0.5D+00 * ( ymin + ymax )

    end if

    py = plotymin2 + nint ( alpha * ( y - ymin2 ) )

    write ( unit, '(a)' ) 'newpath'

    px = plotxmin2 + nint ( alpha * ( xmin - xmin2 ) )
    write ( unit, '(2i6,a)' ) px, py, ' moveto'

    px = plotxmin2 + nint ( alpha * ( xmax - xmin2 ) )
    write ( unit, '(2i6,a)' ) px, py, ' lineto'

    write ( unit, '(a)' ) 'stroke'

  end do

  return
end
subroutine ps_line_closed ( npoint, x, y )

!*****************************************************************************80
!
!! PS_LINE_CLOSED adds the graph of a closed line to a PostScript file.
!
!  Discussion:
!
!    A "closed" line is one in which the last point is connected back
!    to the first one.
!
!    The current point is set to the first (and logically last) point
!    in the list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPOINT, the number of points in the line.
!
!    Input, real ( kind = 8 ) X(NPOINT), Y(NPOINT), the X and Y components
!    of the points.
!
  implicit none

  integer ( kind = 4 ) npoint

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px
  integer ( kind = 4 ) py
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  real ( kind = 8 ) x(npoint)
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(npoint)
  real ( kind = 8 ) ymin
!
!  Refuse to handle fewer than 2 points.
!
  if ( npoint < 2 ) then
    return
  end if
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_LINE_CLOSED - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )

  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )

  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
!
!  Draw lines.
!
  write ( unit, '(a)' ) 'newpath'

  px = plotxmin2 + nint ( alpha * ( x(1) - xmin ) )
  py = plotymin2 + nint ( alpha * ( y(1) - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

  do i = 2, npoint
    px = plotxmin2 + nint ( alpha * ( x(i) - xmin ) )
    py = plotymin2 + nint ( alpha * ( y(i) - ymin ) )
    write ( unit, '(2i6,a)' ) px, py, ' lineto'
  end do
!
!  Add the final extra segment to the initial point.
!
  px = plotxmin2 + nint ( alpha * ( x(1) - xmin ) )
  py = plotymin2 + nint ( alpha * ( y(1) - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' lineto'
!
!  Draw the line.
!
  write ( unit, '(a)' ) 'stroke'

  call ps_setting_real ( 'SET', 'XCUR', x(1) )
  call ps_setting_real ( 'SET', 'YCUR', y(1) )

  return
end
subroutine ps_lineto ( x, y )

!*****************************************************************************80
!
!! PS_LINETO draws a line from the current point to the given point.
!
!  Discussion:
!
!    The current point is updated to the given point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the X and Y components of the new point.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px
  integer ( kind = 4 ) py
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  real ( kind = 8 ) x
  real ( kind = 8 ) xcur
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) ycur
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_LINE - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )

  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )

  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XCUR', xcur )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'YCUR', ycur )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
!
!  Draw the line.
!
  write ( unit, '(a)' ) 'newpath'

  px = plotxmin2 + nint ( alpha * ( xcur - xmin ) )
  py = plotymin2 + nint ( alpha * ( ycur - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

  px = plotxmin2 + nint ( alpha * ( x - xmin ) )
  py = plotymin2 + nint ( alpha * ( y - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' lineto'
!
!  Draw the line.
!
  write ( unit, '(a)' ) 'stroke'

  call ps_setting_real ( 'SET', 'XCUR', x )
  call ps_setting_real ( 'SET', 'YCUR', y )

  return
end
subroutine ps_mark_circles ( n, x, y )

!*****************************************************************************80
!
!! PS_MARK_CIRCLES marks points with a small open circle.
!
!  Discussion:
!
!    The current point is set to the center of the last circle.
!
!    The circles are drawn with the current RGB line colors.
!
!    The circles are drawn the current marker size.
!
!    Points outside the region are not marked.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the points to mark.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) marker_size
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  logical point_inside_box_2d
  integer ( kind = 4 ) pxcen
  integer ( kind = 4 ) pycen
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_MARK_CIRCLE - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'MARKER_SIZE', marker_size )
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )
  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'XMAX', xmax )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
  call ps_setting_real ( 'GET', 'YMAX', ymax )

  write ( unit, '(a)' ) 'newpath'

  do i = 1, n

    if ( .not. point_inside_box_2d ( xmin, ymin, xmax, ymax, x(i), y(i) ) ) then
      cycle
    end if

    pxcen = plotxmin2 + nint ( alpha * ( x(i) - xmin ) )
    pycen = plotymin2 + nint ( alpha * ( y(i) - ymin ) )
    write ( unit, '(3i6,a)' ) pxcen, pycen, marker_size, &
      ' 0 360 arc closepath stroke'

  end do

  call ps_setting_real ( 'SET', 'XCUR', x )
  call ps_setting_real ( 'SET', 'YCUR', y )

  return
end
subroutine ps_mark_disk ( x, y )

!*****************************************************************************80
!
!! PS_MARK_DISK marks a point with a small filled disk.
!
!  Discussion:
!
!    The current point is set to the center of the disk.
!
!    The circle is drawn with the current RGB fill colors.
!
!    The circle is drawn the current marker size.
!
!    Points outside the region are not marked.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of the point to mark.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) marker_size
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  logical point_inside_box_2d
  integer ( kind = 4 ) pxcen
  integer ( kind = 4 ) pycen
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_MARK_DISK - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'MARKER_SIZE', marker_size )
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )
  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'XMAX', xmax )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
  call ps_setting_real ( 'GET', 'YMAX', ymax )
!
!  If the point is outside the plot box, don't draw it.
!
  if ( .not. point_inside_box_2d ( xmin, ymin, xmax, ymax, x, y ) ) then
    return
  end if

  write ( unit, '(a)' ) 'newpath'

  pxcen = plotxmin2 + nint ( alpha * ( x - xmin ) )
  pycen = plotymin2 + nint ( alpha * ( y - ymin ) )

  write ( unit, '(3i6,a)' ) pxcen, pycen, marker_size, &
    ' 0 360 arc closepath fill'

  call ps_setting_real ( 'SET', 'XCUR', x )
  call ps_setting_real ( 'SET', 'YCUR', y )

  return
end
subroutine ps_marker_size ( marker_size )

!*****************************************************************************80
!
!! PS_MARKER_SIZE sets the marker size.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MARKER_SIZE, the marker size.
!    0 is invisible, 1 is a single point.
!    A typical value is 3, 5 or 8.
!
  implicit none

  integer ( kind = 4 ) marker_size
  integer ( kind = 4 ) state
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 2 .and. state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_MARKER_SIZE - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 2 or 3 is required.'
    return
  end if

  call ps_setting_int ( 'SET', 'MARKER_SIZE', marker_size )

  return
end
subroutine ps_moveto ( x, y )

!*****************************************************************************80
!
!! PS_MOVETO "moves to" a new point, which becomes the current point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the X and Y components of the current point.
!
  implicit none
!
  real ( kind = 8 ) alpha
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px
  integer ( kind = 4 ) py
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  real ( kind = 8 ) x
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_MOVETO - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )

  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )

  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
!
!  Move to the new point.
!
  px = plotxmin2 + nint ( alpha * ( x - xmin ) )
  py = plotymin2 + nint ( alpha * ( y - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

  call ps_setting_real ( 'SET', 'XCUR', x )
  call ps_setting_real ( 'SET', 'YCUR', y )

  return
end
subroutine ps_page_head ( xmin, ymin, xmax, ymax )

!*****************************************************************************80
!
!! PS_PAGE_HEAD writes header information on a new page.
!
!  Discussion:
!
!    I think an earlier version of this code, which wrote
!    "%% Page:" rather than "%%Page:" may have caused problems
!    for some interpreters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XMIN, YMIN, XMAX, YMAX, the minimum and maximum X
!    and Y values of the data to be drawn on this page.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) num_pages
  integer ( kind = 4 ) state
  real ( kind = 8 ) line_blue
  real ( kind = 8 ) line_green
  real ( kind = 8 ) line_red
  integer ( kind = 4 ) margin
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
  integer ( kind = 4 ) unit
  real ( kind = 8 ) xcur
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmax2
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xmin2
  real ( kind = 8 ) xvec(4)
  real ( kind = 8 ) ycur
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymax2
  real ( kind = 8 ) ymin
  real ( kind = 8 ) ymin2
  real ( kind = 8 ) yvec(4)
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state == 3 ) then
    state = 2
    call ps_setting_int ( 'SET', 'STATE', state )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_PAGE_HEAD - Warning!'
    write ( *, '(a)' ) '  The current open page is forced closed.'
  end if

  if ( state /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_PAGE_HEAD - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 2 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'NUM_PAGES', num_pages )

  num_pages = num_pages + 1

  call ps_setting_int ( 'SET', 'NUM_PAGES', num_pages )

  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(a,i6,i6)' ) '%%Page: ', num_pages, num_pages
  write ( unit, '(a)' ) 'save'
!
!  Reset the state.
!
  state = 3

  call ps_setting_int ( 'SET', 'STATE', state )
!
!  Determine and store parameters.
!
  if ( xmax == xmin ) then
    xmax2 = xmax + 1.0D+00
    xmin2 = xmax - 1.0D+00
  else
    xmax2 = xmax
    xmin2 = xmin
  end if

  if ( ymax == ymin ) then
    ymax2 = ymax + 1.0D+00
    ymin2 = ymax - 1.0D+00
  else
    ymax2 = ymax
    ymin2 = ymin
  end if
!
!  Set the value of "current point".
!
  xcur = xmin
  ycur = ymin
!
!  Set the conversion factors.
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

  alpha = min ( real ( plotxmax - plotxmin, kind = 8 ) / ( xmax2 - xmin2 ), &
                real ( plotymax - plotymin, kind = 8 ) / ( ymax2 - ymin2 ) )
!
!  Adjust PLOTXMIN and PLOTYMIN to center the image.
!
  plotxmin2 = nint ( 0.5D+00 * &
    ( real ( plotxmin + plotxmax, kind = 8 ) - alpha * ( xmax2 - xmin2 ) ) )

  plotymin2 = nint ( 0.5D+00 * &
    ( real ( plotymin + plotymax, kind = 8 ) - alpha * ( ymax2 - ymin2 ) ) )
!
!  Store data.
!
  call ps_setting_int ( 'SET', 'PXMIN', plotxmin2 )
  call ps_setting_int ( 'SET', 'PYMIN', plotymin2 )

  call ps_setting_real ( 'SET', 'ALPHA', alpha )
  call ps_setting_real ( 'SET', 'XCUR', xcur )
  call ps_setting_real ( 'SET', 'XMIN', xmin )
  call ps_setting_real ( 'SET', 'XMAX', xmax )
  call ps_setting_real ( 'SET', 'YCUR', ycur )
  call ps_setting_real ( 'SET', 'YMIN', ymin )
  call ps_setting_real ( 'SET', 'YMAX', ymax )
!
!  Draw a gray border around the page.
!
  line_red = 0.9D+00
  line_green = 0.9D+00
  line_blue = 0.9D+00

  call ps_color_line ( 'PUSH', line_red, line_green, line_blue )

  call ps_comment ( 'Draw a gray border around the page.' )

  xvec(1:4) = (/ xmin, xmax, xmax, xmin /)
  yvec(1:4) = (/ ymin, ymin, ymax, ymax /)

  call ps_line_closed ( 4, xvec, yvec )

  call ps_color_line ( 'POP', line_red, line_green, line_blue )

  return
end
subroutine ps_page_tail ( )

!*****************************************************************************80
!
!! PS_PAGE_TAIL writes tail information at the end of a page.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    None
!
  implicit none

  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_PAGE_TAIL - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(a)' ) 'restore showpage'

  call ps_comment ( 'End of page' )
!
!  Reset the state.
!
  state = 2

  call ps_setting_int ( 'SET', 'STATE', state )

  return
end
subroutine ps_setting_int ( action, variable, value )

!*****************************************************************************80
!
!! PS_SETTING_INT sets, gets, or prints integer internal PS_WRITE parameters.
!
!  Discussion:
!
!    Normally, the user does not call this routine.  It is a utility
!    used by the package.
!
!    I'd like a more sophisticated pop and push.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the desired action:
!    'GET' to get the current value of VARIABLE, or
!    'POP' to return the current value and set a new value;
!    'SET' to set a new value of VARIABLE, or
!    'PUSH' to return the current value and set a new value;
!    'PRINT' to print the current value of VARIABLE.
!
!    Input, character ( len = * ) VARIABLE, the variable to get or set:
!    'LINE_WIDTH', the line width.
!      0 is the very thinnest line possible,
!      1 is more usual, 2 is thicker, and so on.
!    'MARKER_SIZE', the size of marker circles and disks, in PostScript points;
!    'NUM_PAGES', the number of pages begun or completed;
!    'PXMIN', the location of the left hand margin of the region
!       in PostScript points;
!    'PYMIN', the location of the lower margin of the region
!       in PostScript points;
!    'STATE', the current internal state,
!      0, file not open,
!      1, file open, no header written, no page open,
!      2, file open, header written, no page open,
!      3, file open, header written, page open.
!      4, file open, header written, trailer written.
!    'UNIT', the FORTRAN output unit associated with the PostScript file.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    If ACTION = 'GET', then VALUE is an output quantity, and is the
!    current internal value of the variable.
!
!    If ACTION = 'SET', then VALUE is an input quantity, and the
!    current internal value of the variable is set to this value.
!
!    If ACTION = 'PRINT', then VALUE is ignored.
!
  implicit none

  character ( len = * ) action
  integer ( kind = 4 ), save :: line_width = 1
  integer ( kind = 4 ), save :: marker_size = 0
  integer ( kind = 4 ), save :: num_pages = 0
  integer ( kind = 4 ), save :: pxmin = 0
  integer ( kind = 4 ), save :: pymin = 0
  integer ( kind = 4 ), save :: state = 0
  integer ( kind = 4 ), save :: unit = 0
  integer ( kind = 4 ) value
  character ( len = * ) variable

  if ( variable == 'LINE_WIDTH' ) then

    if ( action == 'GET' ) then
      value = line_width
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Line width, LINE_WIDTH = ', line_width
    else if ( action == 'SET' ) then
      line_width = value
    else if ( action == 'POP' ) then
      call i4_swap ( line_width, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( line_width, value )
    end if

  else if ( variable == 'MARKER_SIZE' ) then

    if ( action == 'GET' ) then
      value = marker_size
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Marker size, MARKER_SIZE = ', marker_size
    else if ( action == 'SET' ) then
      marker_size = value
    else if ( action == 'POP' ) then
      call i4_swap ( marker_size, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( marker_size, value )
    end if

  else if ( variable == 'NUM_PAGES' ) then

    if ( action == 'GET' ) then
      value = num_pages
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Number of pages, NUM_PAGES = ', num_pages
    else if ( action == 'SET' ) then
      num_pages = value
    end if

  else if ( variable == 'PXMIN' ) then

    if ( action == 'GET' ) then
      value = pxmin
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'PostScript minimum X point, PXMIN = ', pxmin
    else if ( action == 'SET' ) then
      pxmin = value
    else if ( action == 'POP' ) then
      call i4_swap ( pxmin, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( pxmin, value )
    end if

  else if ( variable == 'PYMIN' ) then

    if ( action == 'GET' ) then
      value = pymin
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'PostScript minimum Y point, PYMIN = ', pymin
    else if ( action == 'SET' ) then
      pymin = value
    else if ( action == 'POP' ) then
      call i4_swap ( pymin, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( pymin, value )
    end if

  else if ( variable == 'STATE' ) then

    if ( action == 'GET' ) then
      value = state
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Current internal state, STATE = ', state
    else if ( action == 'SET' ) then
      state = value
    else if ( action == 'POP' ) then
      call i4_swap ( state, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( state, value )
    end if

  else if ( variable == 'UNIT' ) then

    if ( action == 'GET' ) then
      value = unit
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Current FORTRAN unit, UNIT = ', unit
    else if ( action == 'SET' ) then
      unit = value
    else if ( action == 'POP' ) then
      call i4_swap ( unit, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( unit, value )
    end if

  end if

  return
end
subroutine ps_setting_real ( action, variable, value )

!*****************************************************************************80
!
!! PS_SETTING_REAL sets, gets, or prints real internal PS_WRITE parameters.
!
!  Discussion:
!
!    I'd like a more sophisticated pop and push.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is either:
!    'GET' to get the current value, or
!    'POP' to return the current value and set a new one;
!    'PRINT' to print the current value, or
!    'SET' to set the current value or
!    'PUSH' to set a new value and return the current one.
!
!    Input, character ( len = * ) VARIABLE, the variable to get or set:
!    'ALPHA', the scale factor from XY user space to PostScript points;
!    'FILL_BLUE', the intensity of the blue fill color, between 0.0 and 1.0.
!    'FILL_GREEN', the intensity of the green fill color, between 0.0 and 1.0.
!    'FILL_RED', the intensity of the red fill color, between 0.0 and 1.0.
!    'FONT_SIZE', the font size, in inches.
!    'LINE_BLUE', the blue component of the line color, between 0.0 and 1.0.
!    'LINE_GREEN', the green component of the line color, between 0.0 and 1.0.
!    'LINE_RED', the red component of the line color, between 0.0 and 1.0.
!    'XCUR', the current X location.
!    'XMAX', maximum X value of the data.
!    'XMIN', minimum X value of the data.
!    'YCUR', the current Y location.
!    'YMAX', maximum Y value of the data.
!    'YMIN', minimum Y value of the data.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    If ACTION = 'GET', then VALUE is an output quantity, and is the
!    current internal value of the variable.
!
!    If ACTION = 'SET', then VALUE is an input quantity, and the
!    current internal value of the variable is set to this value.
!
!    If ACTION = 'PRINT', then VALUE is ignored.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: alpha = 0.0D+00
  real ( kind = 8 ), save :: fill_blue = 0.7D+00
  real ( kind = 8 ), save :: fill_green = 0.7D+00
  real ( kind = 8 ), save :: fill_red = 0.7D+00
  real ( kind = 8 ), save :: font_size = 0.1D+00
  real ( kind = 8 ), save :: line_blue = 0.0D+00
  real ( kind = 8 ), save :: line_green = 0.0D+00
  real ( kind = 8 ), save :: line_red = 0.0D+00
  real ( kind = 8 ) value
  character ( len = * ) variable
  real ( kind = 8 ), save :: xcur = 0.0D+00
  real ( kind = 8 ), save :: xmax = 1.0D+00
  real ( kind = 8 ), save :: xmin = 0.0D+00
  real ( kind = 8 ), save :: ycur = 0.0D+00
  real ( kind = 8 ), save :: ymax = 0.0D+00
  real ( kind = 8 ), save :: ymin = 0.0D+00

  if ( variable == 'ALPHA' ) then

    if ( action == 'GET' ) then
      value = alpha
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Scale factor from user to PS, ALPHA = ', alpha
    else if ( action == 'SET' ) then
      alpha = value
    else if ( action == 'POP' ) then
      call r8_swap ( alpha, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( alpha, value )
    end if

  else if ( variable == 'FILL_BLUE' ) then

    if ( action == 'GET' ) then
      value = fill_blue
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Blue fill RGB value, FILL_BLUE = ', fill_blue
    else if ( action == 'SET' ) then
      fill_blue = value
    else if ( action == 'POP' ) then
      call r8_swap ( fill_blue, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( fill_blue, value )
    end if

  else if ( variable == 'FILL_GREEN' ) then

    if ( action == 'GET' ) then
      value = fill_green
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Green fill RGB value, FILL_GREEN = ', fill_green
    else if ( action == 'SET' ) then
      fill_green = value
    else if ( action == 'POP' ) then
      call r8_swap ( fill_green, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( fill_green, value )
    end if

  else if ( variable == 'FILL_RED' ) then

    if ( action == 'GET' ) then
      value = fill_red
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'RED fill RGB value, FILL_RED = ', fill_red
    else if ( action == 'SET' ) then
      fill_red = value
    else if ( action == 'POP' ) then
      call r8_swap ( fill_red, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( fill_red, value )
    end if

  else if ( variable == 'FONT_SIZE' ) then

    if ( action == 'GET' ) then
      value = font_size
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Font size, FONT_SIZE = ', font_size
    else if ( action == 'SET' ) then
      font_size = value
    else if ( action == 'POP' ) then
      call r8_swap ( font_size, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( font_size, value )
    end if

  else if ( variable == 'LINE_BLUE' ) then

    if ( action == 'GET' ) then
      value = line_blue
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Blue line RGB value, LINE_BLUE = ', line_blue
    else if ( action == 'SET' ) then
      line_blue = value
    else if ( action == 'POP' ) then
      call r8_swap ( line_blue, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( line_blue, value )
    end if

  else if ( variable == 'LINE_GREEN' ) then

    if ( action == 'GET' ) then
      value = line_green
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Green line RGB value, LINE_GREEN = ', line_green
    else if ( action == 'SET' ) then
      line_green = value
    else if ( action == 'POP' ) then
      call r8_swap ( line_green, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( line_green, value )
    end if

  else if ( variable == 'LINE_RED' ) then

    if ( action == 'GET' ) then
      value = line_red
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Red line RGB value, LINE_RED = ', line_red
    else if ( action == 'SET' ) then
      line_red = value
    else if ( action == 'POP' ) then
      call r8_swap ( line_red, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( line_red, value )
    end if

  else if ( variable == 'XCUR' ) then

    if ( action == 'GET' ) then
      value = xcur
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Current X location, XCUR = ', xcur
    else if ( action == 'SET' ) then
      xcur = value
    else if ( action == 'POP' ) then
      call r8_swap ( xcur, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( xcur, value )
    end if

  else if ( variable == 'XMAX' ) then

    if ( action == 'GET' ) then
      value = xmax
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Maximum X value, XMAX = ', xmax
    else if ( action == 'SET' ) then
      xmax = value
    else if ( action == 'POP' ) then
      call r8_swap ( xmax, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( xmax, value )
    end if

  else if ( variable == 'XMIN' ) then

    if ( action == 'GET' ) then
      value = xmin
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Minimum X value, XMIN = ', xmin
    else if ( action == 'SET' ) then
      xmin = value
    else if ( action == 'POP' ) then
      call r8_swap ( xmin, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( xmin, value )
    end if

  else if ( variable == 'YCUR' ) then

    if ( action == 'GET' ) then
      value = ycur
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Current Y location, YCUR = ', ycur
    else if ( action == 'SET' ) then
      ycur = value
    else if ( action == 'POP' ) then
      call r8_swap ( ycur, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( ycur, value )
    end if

  else if ( variable == 'YMAX' ) then

    if ( action == 'GET' ) then
      value = ymax
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Maximum Y value, YMAX = ', ymax
    else if ( action == 'SET' ) then
      ymax = value
    else if ( action == 'POP' ) then
      call r8_swap ( ymax, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( ymax, value )
    end if

  else if ( variable == 'YMIN' ) then

    if ( action == 'GET' ) then
      value = ymin
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Minimum Y value, YMIN = ', ymin
    else if ( action == 'SET' ) then
      ymin = value
    else if ( action == 'POP' ) then
      call r8_swap ( ymin, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( ymin, value )
    end if

  end if

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two real values.
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
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
subroutine radius_make_1 ( m, a, b, n, basis_m, density_function, &
  center, radius )

!*****************************************************************************80
!
!! RADIUS_MAKE_1 uses algorithm 1 to determine a suitable radius for each center.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2001
!
!  Author:
!
!    Lili Ju
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, real ( kind = 8 ) A(M), B(M), the coordinates of the
!    two extreme corners of the box that defines the region.
!
!    Input, integer ( kind = 4 ) N, the number of Voronoi regions.
!
!    Input, integer ( kind = 4 ) BASIS_M, ?
!
!    Input, integer ( kind = 4 ) DENSITY_FUNCTION, specifies the density function.
!    1: d(x) = 1.0;
!    2: d(x) = exp ( - 4.0 * ( sum(x(1:n)**2) ) )
!    3: d(x) = exp ( - 3.0 * ( 1.0 - sum(x(1:n)**2) ) )
!
!    Input, real ( kind = 8 ) CENTER(M,N), the centers of Voronoi regions.
!
!    Output, real ( kind = 8 ) RADIUS(N), the radii.
!
  implicit none

  integer ( kind = 4 ) basis_m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) asspts(m,basis_m*basis_m)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) cc
  real ( kind = 8 ) center(m,n)
  real ( kind = 8 ) dd
  integer ( kind = 4 ) density_function
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) mm
  real ( kind = 8 ) pt(m)
  real ( kind = 8 ) radius(n)
  real ( kind = 8 ) s
!
!  If 2 < M, do we want MM = BASIS_M**M?
!
  mm = basis_m * basis_m

  call random_number ( asspts(1:m,1:mm/2) )
!
!  This CAN'T be right for 2 < M:
!  And how are we being assured that ASSPTS(2,I) is between A(2) and B(2)???
!
  do i = 1, mm/2
    asspts(1,i) = a(1) + ( b(1) - a(1) ) * asspts(1,i)
    cc = - sqrt ( 1.0D+00 - asspts(1,i)**2 )
    dd = - cc
    asspts(2,i) = cc + ( dd - cc ) * asspts(2,i)
  end do

  do i = mm/2+1, mm
    call random_generator ( m, a, b, density_function, asspts(1:m,i) )
  end do

  radius(1:n) = 0.0

  do i = 1, mm
    pt(1:m) = asspts(1:m,i)
    call find_closest ( m, pt, n, center, ic, s )

    if ( radius(ic) < s ) then
      radius(ic) = s
    end if

  end do

  return
end
subroutine radius_make_2 ( m, a, b, n, basis_m, center, radius )

!*****************************************************************************80
!
!! RADIUS_MAKE_2 uses algorithm 2 to determine a suitable radius for each center.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2001
!
!  Author:
!
!    Lili Ju
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, real ( kind = 8 ) A(M), B(M), the coordinates of the
!    two extreme corners of the box that defines the region.
!
!    Input, integer ( kind = 4 ) N, the number of Voronoi regions.
!
!    Input, integer ( kind = 4 ) BASIS_M, ?
!
!    Input, real ( kind = 8 ) CENTER(M,N), the centers of Voronoi regions.
!
!    Output, real ( kind = 8 ) RADIUS(N), the radii.
!
  implicit none

  integer ( kind = 4 ) basis_m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) center(m,n)
  real ( kind = 8 ) dd
  real ( kind = 8 ) dist(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) np
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) mmn
  integer ( kind = 4 ) ord(n)
  real ( kind = 8 ) pps(m,n+basis_m*basis_m)
  real ( kind = 8 ) pt(m)
  real ( kind = 8 ) r
  real ( kind = 8 ) radius(n)

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RADIUS_MAKE_2 - Error!'
    write ( *, '(a)' ) '  The spatial dimension must be positive!'
    write ( *, '(a)' ) '  Enter a spatial dimension with the "M = " command.'
    return
  end if

  if ( basis_m <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RADIUS_MAKE_2 - Fatal error!'
    write ( *, '(a)' ) '  Input parameter BASIS_M <= 1.'
    write ( *, '(a,i6)' ) '  BASIS_M = ', basis_m
    stop
  end if

  mm = basis_m * basis_m
  mmn = mm + n

  pps(1:m,1:n) = center(1:m,1:n)
!
!  Select BASIS_M**2 evenly spaced points in [A(1:M),B(1:M)]
!
!  Need to modify this for 2 < M.
!  This is hard-wired for 2D right now!
!
  do i = 1, basis_m
    do j = 1, basis_m
      pps(1,n+(i-1)*basis_m+j) = &
        a(1) + dble (i-1) * ( b(1) - a(1) ) / dble ( basis_m - 1 )
      pps(2,n+(i-1)*basis_m+j) = &
        a(2) + dble (j-1) * ( b(2) - a(2) ) / dble ( basis_m - 1 )
    end do
  end do
!
!  This is hard-wired for 2D right now!
!
  r = 0.2 * sqrt ( ( b(1) - a(1) ) * ( b(2) - a(2) ) / dble ( n ) )

  radius(1:n) = 0.0

  do i = 1, mm + n

    pt(1:m) = pps(1:m,i)
    ierr = -1
    dd = 0.0D+00
    kk = 0

    do while ( ierr /= 0 )

      call find_re ( m, pt, r, n, center, np, dist, ord, ierr )

      if ( ierr /= 0 ) then
        r = 1.5 * r
      else
        dd = 2.0
        do j = 1, np
          if ( dist(j) < dd ) then
            dd = dist(j)
            kk = ord(j)
          end if
        end do
      end if

    end do

    if ( kk /= 0 ) then
      radius(kk) = max ( radius(kk), dd )
    end if

  end do

  return
end
subroutine random_generator ( m, a, b, density_function, z )

!*****************************************************************************80
!
!! RANDOM_GENERATOR returns a random point Z(1:M) in [A(1:M),B(1:M)].
!
!  Discussion:
!
!    The region is an M dimensional box.
!
!    A density function DENSITY(X(1:M) controls whether the points are 
!    uniformly distributed in the region, or are clustered or scattered.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2001
!
!  Author:
!
!    Lili Ju
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper
!    ranges for the variable.
!
!    Input, integer ( kind = 4 ) DENSITY_FUNCTION, specifies the density function.
!    1: density(x) = 1.0
!    2: density(x) = exp ( - 4.0 * ( sum(x(1:n)**2) ) )
!    3: density(x) = exp ( - 3.0 * ( 1.0 - sum(x(1:n)**2) ) )
!
!    Output, real ( kind = 8 ) Z(M), the random point.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) density
  integer ( kind = 4 ) density_function
  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) z(m)

  do
!
!  Generate a point at random.
!
    do i = 1, m
      call random_number ( r )
      x(i) = ( ( 1.0D+00 - r ) * a(i) + r * b(i) )
    end do
!
!  The density function determines if we should accept or reject the point.
!
    call random_number ( r )

    if ( r < density ( m, x, density_function ) ) then
      z(1:m) = x(1:m)
      exit
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
function s_eqi ( s1, s2 )

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
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character ( len = * ) s1
  character ( len = * ) s2

  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

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
!    If the string is blank, then IVAL will be returned 0.
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
subroutine set_random_seed ( myrank )

!*****************************************************************************80
!
!! SET_RANDOM_SEED initializes the FORTRAN90 random number generator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2001
!
!  Author:
!
!    Lili Ju
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MYRANK, an identifier for each processor.  (This
!    is only used in multi-processor applications).
!
  implicit none

  character ( len = 10 ) big_ben(3)
  integer ( kind = 4 ) date_time(8)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) myrank
  integer ( kind = 4 ), allocatable :: seed(:)
!
!  Initialize the random seed routine.
!
  call random_seed ( )
!
!  Request the size of a typical seed.
!  (It's probably just 1!)
!
  call random_seed ( size = k )
!
!  Set up space for a seed.
!
  allocate ( seed(k) )
!
!  Get the date and time.
!
  call date_and_time ( big_ben(1), big_ben(2), big_ben(3), date_time )
!
!  Make up a "random" value based on date and time information.
!
  do i = 1, k
    seed(i) = date_time(8-mod(i-1,8)) + i * ( myrank + 1 ) * 100
  end do
!
!  Send this random value back to the RANDOM_SEED routine, to be
!  used as the seed of the random number generator.
!
  call random_seed ( put = seed(1:k) )

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
subroutine uniform_make ( m, n, a, b, density_function, points )

!*****************************************************************************80
!
!! UNIFORM_MAKE sets up the uniform random data points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 February 2004
!
!  Author:
!
!    Lili Ju
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of data points to generate.
!
!    Input, real ( kind = 8 ) A(M), B(M), the coordinates of the
!    two extreme corners of the box that defines the region.
!
!    Input, integer ( kind = 4 ) DENSITY_FUNCTION, specifies the density function.
!    1: d(x) = 1.0
!    2: d(x) = exp ( - 4.0 * ( sum(x(1:n)**2) ) )
!    3: d(x) = exp ( - 3.0 * ( 1.0 - sum(x(1:n)**2) ) )
!
!    Output, real ( kind = 8 ) POINTS(M,N), the data points.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) m

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) density_function
  integer ( kind = 4 ) i
  real ( kind = 8 ) points(m,n)

  do i = 1, n
    call random_generator ( m, a, b, density_function, points(1:m,i) )
  end do

  write ( *, '(a,i6)' ) '  Number of Uniform random data points created was ', n

  return
end
subroutine uniform_read ( m, n, points, input_file )

!*****************************************************************************80
!
!! UNIFORM_READ reads uniform random data from a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, real ( kind = 8 ) POINTS(M,N), the uniform random points.
!
!    Input, character ( len = * ) INPUT_FILE, the name of the file.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  character ( len = * ) input_file
  real ( kind = 8 ), dimension ( m, n ) :: points

  write ( *, '(a)' ) ' '
  write ( *, '(a)' )  'UNIFORM_READ:'
  write ( *, '(a)' ) '  Reading Uniform Random data file: ' &
    // trim ( input_file )
  write ( *, '(a,i6)' ) '  Number of data points is ', n

  open ( unit = 12, file = input_file, form = 'formatted', status = 'old' )

  do i = 1, n
    read ( 12, * ) points(1:m,i)
  end do

  close ( unit = 12 )

  return
end
subroutine uniform_write ( m, n, points, output_file )

!*****************************************************************************80
!
!! UNIFORM_WRITE writes the uniform random data to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) POINTS(M,N), the uniform random points.
!
!    Input, character ( len = * ) OUTPUT_FILE, the name of the output file.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  character ( len = * ) output_file
  real ( kind = 8 ), dimension ( m, n ) :: points

  write ( *, '(a)' ) ' '
  write ( *, '(a)' )  'UNIFORM_WRITE:'
  write ( *, '(a)' ) '  Write Uniform Random data file: ' &
    // trim ( output_file )
  write ( *, '(a,i6)' ) '  Number of data points is ', n

  open ( unit = 12, file = output_file, form = 'formatted', status = 'replace' )

  do i = 1, n
    write ( 12, * ) points(1:m,i)
  end do

  close ( unit = 12 )

  return
end
