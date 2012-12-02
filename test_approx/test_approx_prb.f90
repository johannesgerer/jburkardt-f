program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_APPROX_PRB.
!
!  Discussion:
!
!    TEST_APPROX_PRB calls the TEST_APPROX tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp (  )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_APPROX_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_APPROX library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )

  call test10 ( )
  call test11 ( )
  call test12 ( )
  call test13 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_APPROX_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 shows how P00_TITLE can be called.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Demonstrate some of the bookkeeping routines.'
  write ( *, '(a)' ) '  P00_PROB_NUM returns the number of problems.'
  write ( *, '(a)' ) '  P00_TITLE returns the problem title.'
  write ( *, '(a)' ) '  P00_LIMIT returns the problem limits.'

  call p00_prob_num ( prob_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of problems = ', prob_num
  write ( *, '(a)' ) ' '

  do prob = 1, prob_num

    call p00_title ( prob, title )
    write ( *, '(2x,i2,2x,a)' ) prob, '"' // trim ( title ) // '".'

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 shows how P00_STORY can be called.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  P00_STORY prints the problem "story".'

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', prob

    call p00_story ( prob )

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 uses polynomial interpolation on data vector problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_tab = 12

  real ( kind = 8 ) diftab(max_tab)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  character mark
  integer ( kind = 4 ) ntab
  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  character ( len = 80 ) title
  real ( kind = 8 ) x
  real ( kind = 8 ) xdata(max_tab)
  real ( kind = 8 ) yapprox
  real ( kind = 8 ) ydata(max_tab)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Polynomial interpolation to a vector of data.'

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num

    call p00_title ( prob, title )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', prob
    write ( *, '(2x,a)' ) trim ( title )

    call p00_data_num ( prob, data_num )

    write ( *, '(2x,a,i8)' ) '  DATA_NUM = ', data_num

    if ( max_tab < data_num ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Skipped problem ', prob
      write ( *, '(a)' ) '  Too big.'

    else

      call p00_dat ( prob, data_num, xdata, ydata )

      ntab = data_num

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Interpolating polynomial order = ', ntab
      write ( *, '(a)' ) ' '
!
!  Construct the interpolating polynomial via finite differences.
!
      call data_to_dif ( ntab, xdata, ydata, diftab )
!
!  Print out the approximation, including midpoints of the intervals.
!
      do i = 1, ntab

        if ( i < ntab ) then
          jhi = 2
        else
          jhi = 1
        end if

        do j = 1, jhi

          if ( i < ntab ) then
            x = ( real ( jhi - j + 1, kind = 8 ) * xdata(i)     &
                + real (       j - 1, kind = 8 ) * xdata(i+1) ) &
                / real ( jhi,         kind = 8 )
          else
            x = xdata(ntab)
          end if

          if ( j == 1 ) then
            mark = '*'
          else
            mark = ' '
          end if

          call dif_val ( ntab, xdata, diftab, x, yapprox )

          write ( *, '(2x,a,2g14.6)' ) mark, x, yapprox

        end do

      end do

    end if

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 uses linear spline interpolation on all problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  character mark
  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  character ( len = 80 ) title
  real ( kind = 8 ), allocatable, dimension ( : ) :: xdata
  real ( kind = 8 ) xval
  real ( kind = 8 ), allocatable, dimension ( : ) :: ydata
  real ( kind = 8 ) ypval
  real ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Linear spline interpolation.'

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num

    call p00_title ( prob, title )

    call p00_data_num ( prob, data_num )
    allocate ( xdata(1:data_num) )
    allocate ( ydata(1:data_num) )

    call p00_dat ( prob, data_num, xdata, ydata )

    a = xdata(1)
    b = xdata(data_num)

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', prob
    write ( *, '(2x,a)' ) trim ( title )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '       X          Y          Y'' '
    write ( *, '(a)' ) ' '
!
!  Evaluate the interpolation function.
!
    imax = 2 * data_num - 1

    do i = 1, imax

      xval = ( real ( imax - i,     kind = 8 ) * a   &
             + real (        i - 1, kind = 8 ) * b ) &
             / real ( imax     - 1, kind = 8 )

      call spline_linear_val ( data_num, xdata, ydata, xval, yval, ypval )

      if ( mod ( i, 2 ) == 1 ) then
        mark = '*'
      else
        mark = ' '
      end if

      write ( *, '(2x,a,3g14.6)' ) mark, xval, yval, ypval

    end do

    deallocate ( xdata )
    deallocate ( ydata )

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 uses Overhauser spline interpolation on all problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: num_dim = 1

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jmax
  character mark
  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  character ( len = 80 ) title
  real ( kind = 8 ), allocatable, dimension ( : ) :: xdata
  real ( kind = 8 ) xval
  real ( kind = 8 ), allocatable, dimension ( : ) :: ydata
  real ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Overhauser spline interpolation.'

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num

    call p00_title ( prob, title )

    call p00_data_num ( prob, data_num )

    allocate ( xdata(1:data_num) )
    allocate ( ydata(1:data_num) )

    call p00_dat ( prob, data_num, xdata, ydata )

    a = xdata(1)
    b = xdata(data_num)

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', prob
    write ( *, '(2x,a)' ) trim ( title )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  X   Y'
    write ( *, '(a)' ) ' '
!
!  Evaluate the interpolation function.
!
    do i = 1, data_num - 1

      jmax = 3

      if ( i == data_num - 1 ) then
        jhi = jmax
      else
        jhi = jmax - 1
      end if

      do j = 1, jhi

        xval = ( real ( jmax - j,     kind = 8 ) * xdata(i)     &
               + real (        j - 1, kind = 8 ) * xdata(i+1) ) &
               / real ( jmax     - 1, kind = 8 )

        call spline_overhauser_val ( num_dim, data_num, xdata, ydata, xval, &
          yval )

        if ( j == 1 .or. j == 3 ) then
          mark = '*'
        else
          mark = ' '
        end if

        write ( *, '(2x,a,2g14.6)' ) mark, xval, yval

      end do

    end do

    deallocate ( xdata )
    deallocate ( ydata )

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 uses cubic spline interpolation on all problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jmax
  character mark
  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  character ( len = 80 ) title
  real ( kind = 8 ), allocatable, dimension ( : ) :: xdata
  real ( kind = 8 ) xval
  real ( kind = 8 ) ybcbeg
  real ( kind = 8 ) ybcend
  real ( kind = 8 ), allocatable, dimension ( : ) :: ydata
  real ( kind = 8 ), allocatable, dimension ( : ) :: ypp
  real ( kind = 8 ) yppval
  real ( kind = 8 ) ypval
  real ( kind = 8 ) yval

  ibcbeg = 0
  ibcend = 0
  ybcbeg = 0.0D+00
  ybcend = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Cubic spline interpolation.'

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num

    call p00_title ( prob, title )

    call p00_data_num ( prob, data_num )
    allocate ( xdata(1:data_num) )
    allocate ( ydata(1:data_num) )
    allocate ( ypp(1:data_num) )

    call p00_dat ( prob, data_num, xdata, ydata )

    a = xdata(1)
    b = xdata(data_num)
!
!  Set up the interpolation function.
!
    call spline_cubic_set ( data_num, xdata, ydata, ibcbeg, ybcbeg, &
      ibcend, ybcend, ypp )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', prob
    write ( *, '(2x,a)' ) trim ( title )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    X   Y'
    write ( *, '(a)' ) ' '
!
!  Evaluate the interpolation function.
!
    do i = 1, data_num - 1

      jmax = 3

      if ( i == data_num - 1 ) then
        jhi = jmax
      else
        jhi = jmax - 1
      end if

      do j = 1, jhi

        xval = ( real ( jmax - j,     kind = 8 ) * xdata(i)     &
               + real (        j - 1, kind = 8 ) * xdata(i+1) ) &
               / real ( jmax     - 1, kind = 8 )

        call spline_cubic_val ( data_num, xdata, ydata, ypp, xval, yval, ypval, &
          yppval )

        if ( j == 1 .or. j == 3 ) then
          mark = '*'
        else
          mark = ' '
        end if

        write ( *, '(2x,a,2g14.6)' ) mark, xval, yval

      end do

    end do

    deallocate ( xdata )
    deallocate ( ydata )
    deallocate ( ypp )

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 plots an Overhauser spline interpolant for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) :: approx_filename = 'test07_approx.txt'
  character ( len = 80 ) :: data_filename = 'test07_data.txt'
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) :: jmax = 7
  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) nplot
  integer ( kind = 4 ), parameter :: num_dim = 1
  integer ( kind = 4 ) plot
  integer ( kind = 4 ) prob
  real ( kind = 8 ), allocatable, dimension ( : ) :: xdata
  real ( kind = 8 ), allocatable, dimension ( : ) :: xplot
  real ( kind = 8 ) xval
  real ( kind = 8 ), allocatable, dimension ( : ) :: ydata
  real ( kind = 8 ), allocatable, dimension ( : ) :: yplot
  real ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  Plot an Overhauser spline interpolant for problem 7.'
!
!  Get the problem data.
!
  prob = 7

  call p00_data_num ( prob, data_num )

  allocate ( xdata(1:data_num) )
  allocate ( ydata(1:data_num) )

  call p00_dat ( prob, data_num, xdata, ydata )

  call r8vec2_write ( data_filename, data_num, xdata, ydata )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data values stored in "' &
    // trim ( data_filename ) // '".'
!
!  Evaluate the approximating function.
!
  nplot = ( jmax - 1 ) * ( data_num - 1 ) + 1

  allocate ( xplot(1:nplot) )
  allocate ( yplot(1:nplot) )

  plot = 0

  do i = 1, data_num - 1

    if ( i == data_num - 1 ) then
      jhi = jmax
    else
      jhi = jmax - 1
    end if

    do j = 1, jhi

      xval = ( real ( jmax - j,     kind = 8 ) * xdata(i)     &
             + real (        j - 1, kind = 8 ) * xdata(i+1) ) &
             / real ( jmax     - 1, kind = 8 )

      call spline_overhauser_val ( num_dim, data_num, xdata, ydata, xval, yval )

      plot = plot + 1
      xplot(plot) = xval
      yplot(plot) = yval

    end do

  end do

  call r8vec2_write ( approx_filename, nplot, xplot, yplot )

  write ( *, '(a)' ) '  Approximant values stored in "' &
    // trim ( approx_filename ) // '".'

  deallocate ( xdata )
  deallocate ( xplot )
  deallocate ( ydata )
  deallocate ( yplot )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 plots a cubic spline interpolant for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) :: approx_filename = 'test08_approx.txt'
  character ( len = 80 ) :: data_filename = 'test08_data.txt'
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ), parameter :: jmax = 7
  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) nplot
  integer ( kind = 4 ) plot
  integer ( kind = 4 ) prob
  real ( kind = 8 ), allocatable, dimension ( : ) :: xdata
  real ( kind = 8 ), allocatable, dimension ( : ) :: xplot
  real ( kind = 8 ) xval
  real ( kind = 8 ) ybcbeg
  real ( kind = 8 ) ybcend
  real ( kind = 8 ), allocatable, dimension ( : ) :: ydata
  real ( kind = 8 ), allocatable, dimension ( : ) :: yplot
  real ( kind = 8 ), allocatable, dimension ( : ) :: ypp
  real ( kind = 8 ) yppval
  real ( kind = 8 ) ypval
  real ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  Plot a cubic spline interpolant for problem 7.'

  prob = 7
!
!  Get the data.
!
  call p00_data_num ( prob, data_num )

  allocate ( xdata(1:data_num) )
  allocate ( ydata(1:data_num) )
  allocate ( ypp(1:data_num) )

  call p00_dat ( prob, data_num, xdata, ydata )

  call r8vec2_write ( data_filename, data_num, xdata, ydata )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data values stored in "' &
    // trim ( data_filename ) // '".'
!
!  Set up the interpolation function.
!
  ibcbeg = 0
  ibcend = 0
  ybcbeg = 0.0D+00
  ybcend = 0.0D+00

  call spline_cubic_set ( data_num, xdata, ydata, ibcbeg, ybcbeg, ibcend, &
    ybcend, ypp )
!
!  Evaluate the interpolation function.
!
  plot = 0
  nplot = ( jmax - 1 ) * ( data_num - 1 ) + 1

  allocate ( xplot(1:nplot) )
  allocate ( yplot(1:nplot) )

  do i = 1, data_num - 1

    if ( i == data_num - 1 ) then
      jhi = jmax
    else
      jhi = jmax - 1
    end if

    do j = 1, jhi

      xval = ( real ( jmax - j,     kind = 8 ) * xdata(i)     &
             + real (        j - 1, kind = 8 ) * xdata(i+1) ) &
             / real ( jmax     - 1, kind = 8 )

      call spline_cubic_val ( data_num, xdata, ydata, ypp, xval, yval, ypval, &
        yppval )

      plot = plot + 1
      xplot(plot) = xval
      yplot(plot) = yval

    end do

  end do

  call r8vec2_write ( approx_filename, nplot, xplot, yplot )

  write ( *, '(a)' ) '  Approximant values stored in "' &
    // trim ( approx_filename ) // '".'

  deallocate ( xdata )
  deallocate ( xplot )
  deallocate ( ydata )
  deallocate ( yplot )
  deallocate ( ypp )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 uses B spline approximation on all problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jmax
  character mark
  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) prob
  integer ( kind = 4 ) prob_num
  character ( len = 80 ) title
  real ( kind = 8 ), allocatable, dimension ( : ) :: xdata
  real ( kind = 8 ) xval
  real ( kind = 8 ), allocatable, dimension ( : ) :: ydata
  real ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  B spline approximation.'

  call p00_prob_num ( prob_num )

  do prob = 1, prob_num

    call p00_title ( prob, title )

    call p00_data_num ( prob, data_num )
    allocate ( xdata(1:data_num) )
    allocate ( ydata(1:data_num) )

    call p00_dat ( prob, data_num, xdata, ydata )

    a = xdata(1)
    b = xdata(data_num)

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem ', prob
    write ( *, '(2x,a)' ) trim ( title )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '       X        Y'
    write ( *, '(a)' ) ' '
!
!  Evaluate the interpolation function.
!
    do i = 1, data_num - 1

      jmax = 3

      if ( i == data_num - 1 ) then
        jhi = jmax
      else
        jhi = jmax - 1
      end if

      do j = 1, jhi

        xval = ( real ( jmax - j,     kind = 8 ) * xdata(i)     &
               + real (        j - 1, kind = 8 ) * xdata(i+1) ) &
               / real ( jmax     - 1, kind = 8 )

        call spline_b_val ( data_num, xdata, ydata, xval, yval )

        if ( j == 1 .or. j == 3 ) then
          mark = '*'
        else
          mark = ' '
        end if

        write ( *, '(2x,a,2g14.6)' ) mark, xval, yval

      end do

    end do

    deallocate ( xdata )
    deallocate ( ydata )

  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 plots a B spline approximant for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) :: approx_filename = 'test10_approx.txt'
  character ( len = 80 ) :: data_filename = 'test10_data.txt'
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ), parameter :: jmax = 7
  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) nplot
  integer ( kind = 4 ) plot
  integer ( kind = 4 ) prob
  character ( len = 80 ) title
  real ( kind = 8 ), allocatable, dimension ( : ) :: xdata
  real ( kind = 8 ), allocatable, dimension ( : ) :: xplot
  real ( kind = 8 ) xval
  real ( kind = 8 ), allocatable, dimension ( : ) :: ydata
  real ( kind = 8 ), allocatable, dimension ( : ) :: yplot
  real ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  Plot a B spline approximant for problem 7'

  prob = 7

  call p00_title ( prob, title )
!
!  Get the data.
!
  call p00_data_num ( prob, data_num )

  allocate ( xdata(1:data_num) )
  allocate ( ydata(1:data_num) )

  call p00_dat ( prob, data_num, xdata, ydata )

  call r8vec2_write ( data_filename, data_num, xdata, ydata )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data values stored in "' &
    // trim ( data_filename ) // '".'
!
!  Evaluate the approximation function.
!
  plot = 0
  nplot = ( jmax - 1 ) * ( data_num - 1 ) + 1

  allocate ( xplot(1:nplot) )
  allocate ( yplot(1:nplot) )

  do i = 1, data_num - 1

    if ( i == data_num - 1 ) then
      jhi = jmax
    else
      jhi = jmax - 1
    end if

    do j = 1, jhi

      xval = ( real ( jmax - j,     kind = 8 ) * xdata(i)     &
             + real (        j - 1, kind = 8 ) * xdata(i+1) ) &
             / real ( jmax     - 1, kind = 8 )

      call spline_b_val ( data_num, xdata, ydata, xval, yval )

      plot = plot + 1
      xplot(plot) = xval
      yplot(plot) = yval

    end do

  end do

  call r8vec2_write ( approx_filename, nplot, xplot, yplot )

  write ( *, '(a)' ) '  Approximant values stored in "' &
    // trim ( approx_filename ) // '".'

  deallocate ( xdata )
  deallocate ( xplot )
  deallocate ( ydata )
  deallocate ( yplot )

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 plots a beta spline approximant for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) :: approx_filename = 'test11_approx.txt'
  real ( kind = 8 ) beta1
  real ( kind = 8 ) beta2
  character ( len = 80 ) :: data_filename = 'test11_data.txt'
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ), parameter :: jmax = 7
  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) nplot
  integer ( kind = 4 ) plot
  integer ( kind = 4 ) prob
  character ( len = 80 ) title
  real ( kind = 8 ), allocatable, dimension ( : ) :: xdata
  real ( kind = 8 ), allocatable, dimension ( : ) :: xplot
  real ( kind = 8 ) xval
  real ( kind = 8 ), allocatable, dimension ( : ) :: ydata
  real ( kind = 8 ), allocatable, dimension ( : ) :: yplot
  real ( kind = 8 ) yval

  beta1 = 100.0D+00
  beta2 = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  Plot a beta spline approximant for problem 7'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
  write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2

  prob = 7

  call p00_title ( prob, title )
!
!  Get the data.
!
  call p00_data_num ( prob, data_num )

  allocate ( xdata(1:data_num) )
  allocate ( ydata(1:data_num) )

  call p00_dat ( prob, data_num, xdata, ydata )

  call r8vec2_write ( data_filename, data_num, xdata, ydata )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data values stored in "' &
    // trim ( data_filename ) // '".'
!
!  Evaluate the interpolation function.
!
  plot = 0
  nplot = ( jmax - 1 ) * ( data_num - 1 ) + 1

  allocate ( xplot(1:nplot) )
  allocate ( yplot(1:nplot) )

  do i = 1, data_num - 1

    if ( i == data_num - 1 ) then
      jhi = jmax
    else
      jhi = jmax - 1
    end if

    do j = 1, jhi

      xval = ( real ( jmax - j,     kind = 8 ) * xdata(i)     &
             + real (        j - 1, kind = 8 ) * xdata(i+1) ) &
             / real ( jmax     - 1, kind = 8 )

      call spline_beta_val ( beta1, beta2, data_num, xdata, ydata, xval, yval )

      plot = plot + 1
      xplot(plot) = xval
      yplot(plot) = yval

    end do

  end do

  call r8vec2_write ( approx_filename, nplot, xplot, yplot )

  write ( *, '(a)' ) '  Approximant values stored in "' &
    // trim ( approx_filename ) // '".'

  deallocate ( xdata )
  deallocate ( xplot )
  deallocate ( ydata )
  deallocate ( yplot )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 plots a Bernstein spline approximant for problem 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  character ( len = 80 ) :: approx_filename = 'test12_approx.txt'
  real ( kind = 8 ) b
  character ( len = 80 ) :: data_filename = 'test12_data.txt'
  integer ( kind = 4 ) i
  integer ( kind = 4 ) data_num
  integer ( kind = 4 ), parameter :: nplot = 101
  integer ( kind = 4 ) plot
  integer ( kind = 4 ) prob
  real ( kind = 8 ), allocatable, dimension ( : ) :: xdata
  real ( kind = 8 ), allocatable, dimension ( : ) :: xplot
  real ( kind = 8 ) xval
  real ( kind = 8 ), allocatable, dimension ( : ) :: ydata
  real ( kind = 8 ), allocatable, dimension ( : ) :: yplot
  real ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  Plot a Bernstein approximant for problem 5.'
  write ( *, '(a)' ) '  Note that the Bernstein approximant requires equally'
  write ( *, '(a)' ) '  spaced data!'

  prob = 5
!
!  Get the data.
!
  call p00_data_num ( prob, data_num )

  allocate ( xdata(1:data_num) )
  allocate ( ydata(1:data_num) )

  call p00_dat ( prob, data_num, xdata, ydata )

  call r8vec2_write ( data_filename, data_num, xdata, ydata )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data values stored in "' &
    // trim ( data_filename ) // '".'
!
!  Evaluate the approximant function.
!
  allocate ( xplot(1:nplot) )
  allocate ( yplot(1:nplot) )

  a = xdata(1)
  b = xdata(data_num)

  do plot = 1, nplot

    xval = ( real ( nplot - plot,     kind = 8 ) * a     &
           + real (         plot - 1, kind = 8 ) * b ) &
           / real ( nplot        - 1, kind = 8 )

    call bpab_approx ( data_num - 1, a, b, ydata, xval, yval )

    xplot(plot) = xval
    yplot(plot) = yval

  end do

  call r8vec2_write ( approx_filename, nplot, xplot, yplot )

  write ( *, '(a)' ) '  Approximant values stored in "' &
    // trim ( approx_filename ) // '".'

  deallocate ( xdata )
  deallocate ( xplot )
  deallocate ( ydata )
  deallocate ( yplot )

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 plots a cubic spline interpolant for problem 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nplot = 101

  character ( len = 80 ) :: approx_filename = 'test13_approx.txt'
  character ( len = 80 ) :: data_filename = 'test13_data.txt'
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  integer ( kind = 4 ) j
  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) prob
  real ( kind = 8 ), allocatable, dimension ( : ) :: xdata
  real ( kind = 8 ) xplot(nplot)
  real ( kind = 8 ) xval
  real ( kind = 8 ) ybcbeg
  real ( kind = 8 ) ybcend
  real ( kind = 8 ), allocatable, dimension ( : ) :: ydata
  real ( kind = 8 ) yplot(nplot)
  real ( kind = 8 ), allocatable, dimension ( : ) :: ypp
  real ( kind = 8 ) yppval
  real ( kind = 8 ) ypval
  real ( kind = 8 ) yval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  Plot a cubic spline interpolant for problem 5'

  prob = 5

  call p00_data_num ( prob, data_num )

  allocate ( xdata(1:data_num) )
  allocate ( ydata(1:data_num) )
  allocate ( ypp(1:data_num) )

  call p00_dat ( prob, data_num, xdata, ydata )

  call r8vec2_write ( data_filename, data_num, xdata, ydata )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data values stored in "' &
    // trim ( data_filename ) // '".'
!
!  Set up the interpolation function.
!
  ibcbeg = 0
  ibcend = 0
  ybcbeg = 0.0D+00
  ybcend = 0.0D+00

  call spline_cubic_set ( data_num, xdata, ydata, ibcbeg, ybcbeg, ibcend, &
    ybcend, ypp )
!
!  Evaluate the interpolation function.
!
  do j = 1, nplot

    xval = ( real ( nplot - j,     kind = 8 ) * xdata(1)    &
           + real (         j - 1, kind = 8 ) * xdata(data_num) ) &
           / real ( nplot     - 1, kind = 8 )

    call spline_cubic_val ( data_num, xdata, ydata, ypp, xval, yval, ypval, &
      yppval )

    xplot(j) = xval
    yplot(j) = yval

  end do

  call r8vec2_write ( approx_filename, nplot, xplot, yplot )

  write ( *, '(a)' ) '  Interpolant values stored in "' &
    // trim ( approx_filename ) // '".'

  deallocate ( xdata )
  deallocate ( ydata )
  deallocate ( ypp )

  return
end
