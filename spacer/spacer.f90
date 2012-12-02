program main

!*****************************************************************************80
!
!! MAIN is the main program for SPACER.
!
!  Discussion:
!
!    SPACER reads a euclidean distance matrix file and outputs a
!    principal coordinates analysis.
!
!    SEQUENCE_MAX is the maximum number of sequences (or other objects)
!    that can be ordinated.  To work with a greater value, change the
!    definition of SEQUENCE_MAX in EVERY subroutine that mentions it!
!
!  Modified:
!
!    12 February 2003
!
!  Author:
!
!    Des Higgins,
!    EMBL Data Library,
!    April 1992
!
!  Reference:
!
!    John Gower,
!    Some distance properties of latent root and vector methods
!    used in multivariate analysis,
!    Biometrika,
!    Volume 53, 1966, pages 325-328.
!
!    Des Higgins,
!    Sequence ordinations: a multivariate analysis approach
!    to analysing large sequence data sets,
!    CABIOS,
!    Volume 8, 1992, pages 15-22.
!
  implicit none

  integer, parameter :: sequence_max = 650

  real dist(sequence_max,sequence_max)
  integer i
  character ( len = 100 ) input_name
  integer input_unit
  integer ios
  integer j
  character ( len = 120 ) line
  integer eig_number
  character ( len = 100 ) output_name
  integer output_unit
  integer sequence_number
  character symbol(sequence_max)
  character ( len = 94 ) :: symbol_internal = &
    'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz12345678' // &
    '90!@#$%^&*()`-=[] ; ,./~_+{}|:"<>?'
  character ( len = 100 ) symbol_name
  integer symbol_number
  integer symbol_unit

  symbol_internal(88:88) = char ( 92 )
  symbol_internal(90:90) = char ( 39 )

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPACER'
  write ( *, '(a)' ) '  Principal coordinate analysis program.'
!
!  Get the distance matrix data.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter the distance matrix file name:'
  read ( *, '(a)' ) input_name

  input_unit = 66

  open ( unit = input_unit, file = input_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPACER - Fatal error!'
    write ( *, '(a)' ) '  Could not open the distance matrix file.'
    stop
  end if

  call dist_mat_read ( input_unit, sequence_max, sequence_number, dist )

  close ( unit = input_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of objects  = ', sequence_number
!
!  Get the symbols.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  You can use a file with one symbol per seq. '
  write ( *, '(a)' ) '  Ignore this option by pressing <RETURN>'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter the symbol input filename:'
  read ( *, '(a)' ) symbol_name

  symbol_unit = 69
  symbol_number = 0

  do i = 1, sequence_max
    symbol(i) = '*'
  end do

  if ( symbol_name /= ' ' ) then

    open ( unit = symbol_unit, file = symbol_name, iostat = ios, &
      status = 'old' )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPACER - Fatal error!'
      write ( *, '(a)' ) '  Could not open the symbol file.'
      stop
    end if

    do

      read ( symbol_unit, '(a)', iostat = ios ) line

      if ( ios /= 0 ) then
        exit
      end if

      do i = 1, len ( line )

        if ( line(i:i) /= ' ' ) then

          symbol_number = symbol_number + 1

          if ( symbol_number <= sequence_max ) then
            symbol(symbol_number) = line(i:i)
          end if

        end if

      end do


    end do

    if ( sequence_max < symbol_number ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPACER - Warning!'
      write ( *, '(a)' ) '  Number of symbols in file was ', symbol_number
      write ( *, '(a)' ) '  Number of symbols stored was  ', sequence_max
      symbol_number = sequence_max
    end if

    close ( unit = symbol_unit )

  else

    symbol_number = 94
    do i = 1, symbol_number
      symbol(i) = symbol_internal(i:i)
    end do

  end if
!
!  The number of eigenvectors to be calculated is NEIG.
!
  eig_number = min ( 10, sequence_number-1 )
!
!  Set the euclidean distance.
!
  dist(1:sequence_number,1:sequence_number) = &
    - 0.5E+00 * dist(1:sequence_number,1:sequence_number)**2
!
!  Get the eigenvectors and print the results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter the output filename:'
  read ( *, '(a)' ) output_name

  open ( unit = output_unit, file = output_name, iostat = ios, &
    status = 'replace' )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPACER - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if

  call pcoord ( output_unit, sequence_max, sequence_number, dist, symbol, &
    eig_number )

  close ( unit = output_unit )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPACER:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine dist_mat_read ( input_unit, sequence_max, sequence_number, dist )

!*****************************************************************************80
!
!! DIST_MAT_READ reads the distance matrix from a file.
!
!  Modified:
!
!    12 February 2003
!
!  Author:
!
!    Des Higgins,
!    EMBL Data Library,
!    April 1992
!
!  Parameters:
!
!    Input, integer INPUT_UNIT, the FORTRAN unit from which the data
!    should be read.
!
!    Input, integer SEQUENCE_MAX, the maximum size of the distance matrix.
!
!    Output, integer SEQUENCE_NUMBER, the number of sequences.
!
!    Output, real DIST(SEQUENCE_MAX,SEQUENCE_MAX), the distance matrix.
!
  implicit none

  integer sequence_max

  real dist(sequence_max,sequence_max)
  integer i
  integer input_unit
  integer j
  integer sequence_number

  read ( input_unit, * ) sequence_number
  read ( input_unit, * ) ( ( dist(i,j), j = 1, i ), i = 1, sequence_number )
!
!  Copy the symmetric part.
!
  do i = 1, sequence_number
    dist(1:i-1,i) = dist(i,1:i-1)
  end do

  return
end
subroutine eigen ( n_max, n, a, d )

!*****************************************************************************80
!
!! EIGEN computes the eigenvalues and eigenvectors of a symmetric matrix.
!
!  Discussion:
!
!    EIGEN proceeds through two steps:
!    1) reduce the matrix to tridiagonal form.
!    2) compute the eigenvalues and eigenvectors.
!
!  Modified:
!
!    12 February 2003
!
!  Author:
!
!    Des Higgins,
!    EMBL Data Library,
!    April 1992
!
!  Parameters:
!
!    Input, integer N_MAX, the leading dimension of the matrix.
!    N_MAX must be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real A(N_MAX,N), on input, the N by N symmetric matrix.
!    On output, the N by N matrix of eigenvectors.  The I-th eigenvector
!    is stored in column I of A.
!
!    Output, real D(N), the N eigenvalues of the matrix.
!
  implicit none

  integer n
  integer n_max

  real a(n_max,n)
  real d(n)
  real e(n)

  call tred2 ( n_max, n, a, d, e )

  call tqli ( n_max, n, d, e, a )

  return
end
subroutine eigensort ( n_max, n, evectors, evalues )

!*****************************************************************************80
!
!! EIGENSORT sorts the eigenvalues and adjusts the eigenvectors.
!
!  Discussion:
!
!    Uses a QUICKSORT adapted from Kernighan and Plauger.
!
!  Modified:
!
!    12 February 2003
!
!  Author:
!
!    Des Higgins,
!    EMBL Data Library,
!    April 1992
!
!  Reference:
!
!    Kernighan, Plauger,
!    Software Tools,
!    1986, page 115.
!
!  Parameters:
!
!    Input, integer N_MAX, the leading dimension of the matrix.
!    N_MAX must be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real EVECTORS(N_MAX,N), the N by N eigenvector matrix.
!
!    Input/output, real EVALUES(N), the eigenvalues.
!
  implicit none

  integer n
  integer n_max

  real evalues(n)
  real evectors(n_max,n)
  integer i
  integer j
  integer lst(50)
  integer m
  integer p
  real pivlin
  integer ust(50)

  lst(1) = 1
  ust(1) = n
  p = 1

  do

    if ( p <= 0 ) then
      exit
    end if

    if ( lst(p) >= ust(p) ) then

      p = p - 1

    else

      i = lst(p) - 1
      j = ust(p)
      pivlin = evalues(j)

      do

        if ( j <= i ) then
          exit
        end if

        do

          i = i + 1

          if ( evalues(i) >= pivlin ) then
            exit
          end if

        end do

        do

          j = j - 1

          if ( evalues(j) <= pivlin ) then
            exit
          end if

          if ( j <= i ) then
            exit
          end if

        end do

        if ( i < j ) then

          call r4_swap ( evalues(i), evalues(j) )

          do m = 1, n
            call r4_swap ( evectors(m,i), evectors(m,j) )
          end do

        end if

      end do

      j = ust(p)

      call r4_swap ( evalues(i), evalues(j) )

      do m = 1, n
        call r4_swap ( evectors(m,i), evectors(m,j) )
      end do

      if ( i-lst(p) < ust(p)-i ) then
        lst(p+1) = lst(p)
        ust(p+1) = i - 1
        lst(p) = i + 1
      else
        lst(p+1) = i + 1
        ust(p+1) = ust(p)
        ust(p) = i - 1
      end if

      p = p + 1

    end if

  end do

  return

end
subroutine pcoord ( output_unit, n_max, n, a, symbol, eig_number )

!*****************************************************************************80
!
!! PCOORD performs a principal coordinate analysis of distance matrix data.
!
!  Modified:
!
!    12 February 2003
!
!  Author:
!
!    Des Higgins,
!    EMBL Data Library,
!    April 1992
!
!  Parameters:
!
!    Input, integer OUTPUT_UNIT, the FORTRAN unit number to which
!    the plot information should be written.
!
!    Input, integer N_MAX, the leading dimension of the matrix.
!    N_MAX must be at least N.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real A(N_MAX,N), on input, the N by N symmetric matrix.
!    On output, the N by N matrix of eigenvectors.  The I-th eigenvector
!    is stored in column I of A.
!
!    Input, character SYMBOL(N_MAX), identification symbols to be used in
!    plotting.
!
!    Input, integer EIG_NUMBER, the number of eigenvectors to be analyzed.
!
  implicit none

  integer, parameter :: eig_max = 12

  integer n_max

  real a(n_max,n_max)
  real b(n_max,eig_max)
  real e(n_max)
  integer i
  integer j
  integer jhi
  integer jlo
  integer n
  integer eig_number
  integer output_unit
  real percy
  real r(n_max)
  real sei
  real sum2
  real sumall
  real sumneg
  real sumpos
  character symbol(n_max)
  real xaxis(n_max)
  real yaxis(n_max)

  do i = 1, n
    e(i) = sum ( a(i,1:n) )
  end do

  sum2 = sum ( e(1:n) )
!
!  The sum is the mean for all values.
!
!  Vector E contains the column averages.
!
  e(1:n) = e(1:n) / real ( n )

  sum2 = sum2 / real ( n**2 )

  do i = 1, n
    do j = 1, n
      a(i,j) = a(i,j) - e(i) - e(j) + sum2
    end do
  end do

  write ( output_unit, '(a)' ) ' '
  write ( output_unit, '(a)' ) 'SPACER'
  write ( output_unit, '(a)' ) '  Principal coordinates analysis'
  write ( output_unit, '(a)' ) ' '
  write ( output_unit, '(a)' ) '  The method was first described in:'
  write ( output_unit, '(a)' ) ' '
  write ( output_unit, '(a)' ) '    Gower, J. C. (1966)'
  write ( output_unit, '(a)' ) '    Some distance properties of latent'
  write ( output_unit, '(a)' ) '    root and vector methods used in '
  write ( output_unit, '(a)' ) '    multivariate analysis.'
  write ( output_unit, '(a)' ) '    Biometrika, volume 53, pages 325-328.'
  write ( output_unit, '(a)' ) ' '
  write ( output_unit, '(a)' ) '  This software was described in:'
  write ( output_unit, '(a)' ) ' '
  write ( output_unit, '(a)' ) '    Higgins, D. G. (1992)'
  write ( output_unit, '(a)' ) '    Sequence ordinations: a multivariate'
  write ( output_unit, '(a)' ) '    analysis approach to analysing large'
  write ( output_unit, '(a)' ) '    sequence data sets.'
  write ( output_unit, '(a)' ) '    CABIOS, volume 8, pages 15-22.'
  write ( output_unit, '(a)' ) ' '
  write ( output_unit, '(a)' ) ' '
  write ( output_unit, '(a)' ) '  DIAGONAL ELEMENTS OF TRANSFORMED MATRIX'
  write ( output_unit, '(a)' ) ' '
  write ( output_unit, '(a)' ) '    Squared distance of each point from centroid:'
  write ( output_unit, '(a)' ) ' '

  write ( output_unit, '(1x,5f14.6)' ) ( a(i,i), i = 1, n )

  sum2 = 0.0
  do i = 1, n
    sum2 = sum2 + a(i,i)
  end do

  write ( output_unit, '(a)') ' '
  write ( output_unit, '(a,f14.6)' ) '  TRACE (total variation ) = ', sum2
!
!  Extract the eigenvalues (E) and eigenvectors (A after the call) of A
!
  call eigen ( n_max, n, a, e )
!
!  Sort the eigenvalues (E) in increasing order,
!  and the columns of A accordingly.
!
  call eigensort ( n_max, n, a, e )

  write ( output_unit, '(a)' ) ' '
  write ( output_unit, '(a)' ) '  EIGENVALUES OF TRANSFORMED DISTANCE MATRIX'
  write ( output_unit, '(a)' ) ' '
  write ( output_unit, '(1x,5f14.6)' ) ( e(j), j = n, 1, - 1 )

  sumall = 0.0E+00
  sumpos = 0.0E+00
  sumneg = 0.0E+00

  do i = 1, n

    sumall = sumall + e(i)

    if ( e(i) > 0.0E+00 ) then
      sumpos = sumpos + e(i)
    else if ( e(i) < 0.0E+00 ) then
      sumneg = sumneg + e(i)
    end if

  end do

  write ( output_unit, '(a)' ) ' '
  write ( output_unit, '(a)' ) '  Eigenvalue sums:'
  write ( output_unit, '(a)' ) ' '
  write ( output_unit, '(a,g14.6)' ) '    Positive = ', sumpos
  write ( output_unit, '(a,g14.6)' ) '    Negative = ', sumneg
  write ( output_unit, '(a,g14.6)' ) '    All      = ', sumall

  do i = n, n-eig_number+1, -1
    if ( e(i) < 0.0E+00 ) then
      eig_number = n - i
      exit
    else if ( e(i) == 0.0E+00 ) then
      r(i) = 0.0E+00
    else
      r(i) = sqrt ( e(i) )
    end if
  end do
!
!  Compute B, the scaled eigenvectors.
!
  do j = 1, eig_number
    do i = 1, n
      b(i,j) = a(i,n-j+1) * r(n-j+1)
    end do
  end do

  write ( output_unit, '(a)' ) ' '
  write ( output_unit, '(a)' ) '  Specimen coordinates'
  write ( output_unit, '(a)' ) ' '

  do i = 1, n
    jlo = 1
    jhi = min ( eig_number, 6 )
    write ( output_unit, '(i5,2x,6f12.5)' ) i, ( b(i,j), j = jlo, jhi )
    jlo = jhi + 1
    jhi = eig_number
    write ( output_unit, '(7x,6f12.5)' ) ( b(i,j), j = jlo, jhi )
  end do

  if ( eig_number >= 2 ) then
    write ( output_unit, '(a)' ) ' '
    write ( output_unit, '(a)' ) '  Plot of first (horizontal) and second ' // &
      '(vertical) axes'
    write ( output_unit, '(a)' ) ' '
    do i = 1, n
      xaxis(i) = b(i,1)
      yaxis(i) = b(i,2)
    end do
    call plotxy ( output_unit, n, xaxis, yaxis, symbol )
  end if

  if ( eig_number >= 3 ) then
    write ( output_unit, '(a)' ) ' '
    write ( output_unit, '(a)' ) '  Plot of first (horizontal) and third ' // &
      '(vertical) axes'
    write ( output_unit, '(a)' ) ' '
    do i = 1, n
      yaxis(i) = b(i,3)
    end do
    call plotxy ( output_unit, n, xaxis, yaxis, symbol )
  end if

  write ( output_unit, '(a)' ) ' '
  write ( output_unit, '(a)' ) '         AXIS     PERCENTAGE VARIATION' // &
    '      CUMULATIVE'
  write ( output_unit, '(a)' ) ' '

  sei = 0.0E+00
  do i = n, n-eig_number+1, -1
    percy = e(i) * 100.0E+00 / sumpos
    sei = sei + percy
    write ( output_unit, '(a,i4,7x,f10.2,15x,f10.2)' ) '       ', n-i+1, &
      percy, sei
  end do

  return
end
subroutine plotxy ( output_unit, n, x, y, symbol )

!*****************************************************************************80
!
!! PLOTXY creates a scatterplot of pairs of X, Y data.
!
!  Modified:
!
!    12 February 2003
!
!  Author:
!
!    Des Higgins,
!    EMBL Data Library,
!    April 1992
!
!  Parameters:
!
!    Input, integer OUTPUT_UNIT, the FORTRAN unit number to which
!    the plot information should be written.
!
!    Input, integer N, the number of data items.
!
!    Input, real X(N), Y(N), the points to be plotted.
!
  implicit none

  integer, parameter :: line_len = 79
  integer, parameter :: page_len = 50

  integer n

  character graph(0:line_len,0:page_len)
  integer i
  integer j
  integer output_unit
  character symbol(n)
  real x(n)
  integer xc
  real xmax
  real xmin
  real xrange
  real xscale
  integer xzero
  real y(n)
  integer yc
  real ymax
  real ymin
  real yrange
  real yscale
  integer yzero

  xmax = maxval ( x(1:n) )
  xmin = minval ( x(1:n) )

  xrange = xmax - xmin
  xscale = real ( line_len - 1 ) / xrange

  if ( xmin <= 0.0E+00 ) then
    xzero = int ( - xmin * xscale )
  else
    xzero = 0
  end if

  ymax = maxval ( y(1:n) )
  ymin = minval ( y(1:n) )

  yrange = ymax - ymin
  yscale = real ( page_len - 1 ) / yrange

  if ( ymin <= 0.0E+00 ) then
    yzero = int ( - ymin * yscale )
  else
    yzero = 0
  end if

  do i = 0, line_len
    do j = 0, page_len
      graph(i,j) = ' '
    end do
  end do

  do i = 0, line_len
    graph(i,0) = '-'
    graph(i,yzero) = '-'
    graph(i,page_len) = '-'
  end do

  do i = 0, page_len
    graph(0,i) = '|'
    graph(xzero,i) = '|'
    graph(line_len,i) = '|'
  end do
!
!  For each object I, try to place the symbol ALPH(I) at the location
!  appropriate for (X(I), Y(I)).
!
  do i = 1, n

    xc = 1 + int ( ( x(i) - xmin ) * xscale )
    yc = 1 + int ( ( y(i) - ymin ) * yscale )

    if ( graph(xc,yc) == ' ' .or. graph(xc,yc) == '-' .or. &
         graph(xc,yc) == '|' .or. graph(xc,yc) == symbol(i) ) then
      graph(xc,yc) = symbol(i)
    else
      graph(xc,yc) = '*'
    end if

  end do

  write ( output_unit, '(a)' ) ' '
  do i = page_len, 0, -1
    write ( output_unit, '(200a1)' ) ( graph(j,i), j = 0, line_len )
  end do

  return
end
subroutine r4_swap ( x, y )

!*****************************************************************************80
!
!! R4_SWAP switches two R4's.
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
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real x
  real y
  real z

  z = x
  x = y
  y = z

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
  character ( len = 10 )  time
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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine tqli ( n_max, n, d, e, z )

!*****************************************************************************80
!
!! TQLI finds eigenvalues and eigenvectors of a real symmetric tridiagonal matrix.
!
!  Parameters:
!
!    Input, integer N_MAX, the maximum dimension of the Z array.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real D(N), E(N).  On input, the diagonal and off-diagonal
!    elements of the matrix.  On output, D contains the eigenvalues.
!
!    Output, real Z(N_MAX,N), the eigenvector matrix.
!
  implicit none

  integer n
  integer n_max

  real b
  real c
  real d(n)
  real dd
  real e(n)
  real f
  real g
  integer i
  integer im
  integer iter
  integer k
  integer l
  integer m
  real p
  real r
  real s
  real z(n_max,n)

  do i = 2, n
    e(i-1) = e(i)
  end do

  e(n) = 0.0

  do l = 1, n

    iter = 0

10  continue

    do im = l, n-1
      dd = abs ( d(im) ) + abs ( d(im+1) )
      if ( abs ( e(im) ) + dd == dd ) then
        m = im
        go to 20
      end if
    end do

    m = n

20  continue

    if ( m /= l ) then

      iter = iter + 1

      if ( iter >= 30 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TQLI - Fatal error!'
        write ( *, '(a)' ) '  No convergence.'
        write ( *, '(a,i6)' ) '  Number of iterations = ', iter
        stop
      end if

      g = ( d(l+1) - d(l) ) / ( 2.0 * e(l) )
      r = sqrt ( ( g * g ) + 1.0 )
      g = d(m) - d(l) + e(l) / ( g + sign ( r, g ) )
      c = 1.0
      s = 1.0
      p = 0.0

      do i = m-1, l, -1

        f = s * e(i)
        b = c * e(i)

        if ( abs ( f ) >= abs ( g ) ) then
          c = g / f
          r = sqrt ( ( c * c ) + 1.0 )
          e(i+1) = f * r
          s = 1.0 / r
          c = c * s
        else
          s = f / g
          r = sqrt ( ( s * s ) + 1.0 )
          e(i+1) = g * r
          c = 1.0 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b

        do k = 1, n
          f = z(k,i+1)
          z(k,i+1) = s * z(k,i) + c * f
          z(k,i) = c * z(k,i) - s * f
        end do

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0

    end if

    if ( m /= l ) then
      go to 10
    end if

  end do

  return
end
subroutine tred2 ( n_max, n, a, d, e )

!*****************************************************************************80
!
!! TRED2 tridiagonalizes a real symmetric matrix.
!
!  Parameters:
!
!    Input, integer N_MAX, the maximum dimension of the array A.
!
!    Input, integer N, the order of the array A.
!
!    Input, real A(N_MAX,N), the N by N array.
!
!    Output, real D(N), E(N), the diagonal and offdiagonal entries of
!    a symmetric matrix that is similar to A.
!
  implicit none

  integer n
  integer n_max

  real a(n_max,n)
  real d(n)
  real e(n)
  real f
  real g
  real h
  real hh
  integer i
  integer j
  integer k
  integer l
  real scale

  do i = n, 2, -1

    l = i - 1
    scale = 0.0
    h     = 0.0

    if ( l > 1 ) then

      do k = 1, l
        scale = scale + abs ( a(i,k) )
      end do

      if ( scale == 0.0 ) then

        e(i) = a(i,l)

      else

        do k = 1, l
          a(i,k) = a(i,k) / scale
          h = h + a(i,k) * a(i,k)
        end do

        f = a(i,l)

        if ( f > 0.0 ) then
          g = - sqrt ( h )
        else
          g = sqrt ( h )
        end if

        e(i) = scale * g
        h = h - f * g
        a(i,l) = f - g

        f = 0.0

        do j = 1, l

          a(j,i) = a(i,j) / h

          g = 0.0
          do k = 1, j
            g = g + a(j,k) * a(i,k)
          end do

          do k = j+1, l
            g = g + a(k,j) * a(i,k)
          end do

          e(j) = g / h
          f = f + a(i,j) * e(j)

        end do

        hh = f / ( h + h )

        do j = 1, l

          f = a(i,j)
          g = e(j) - hh * f
          e(j) = g
          do k = 1, j
            a(j,k) = a(j,k) - ( f * e(k) + g * a(i,k) )
          end do

        end do

      end if

    else

      e(i) = a(i,l)

    end if

    d(i) = h

  end do

  d(1) = 0.0
  e(1) = 0.0

  do i = 1, n

    l = i - 1

    if ( d(i) > 0.0 ) then

      do j = 1, l

        g = 0.0
        do k = 1, l
          g = g + a(i,k) * a(k,j)
        end do

        do k = 1, l
          a(k,j) = a(k,j) - g * a(k,i)
        end do

      end do

    end if

    d(i) = a(i,i)

    a(i,i) = 1.0
    do j = 1, l
      a(i,j) = 0.0
      a(j,i) = 0.0
    end do

  end do

  return
end
