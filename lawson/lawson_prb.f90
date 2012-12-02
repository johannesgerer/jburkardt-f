program main

!*****************************************************************************80
!
!! MAIN is the main program for LAWSON_PRB.
!
!  Discussion:
!
!    LAWSON_PRB calls tests for the LAWSON package.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAWSON_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Tests the LAWSON library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAWSON_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests HFT, HS1 and COV.
!
!  Discussion:
!
!    HFT and HS1 solve least squares problems.
!    COV computes the associated covariance matrix.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
  implicit none

  integer, parameter :: mda = 8

  real ( kind = 8 ) a(mda,mda)
  real ( kind = 8 ) anoise
  real ( kind = 8 ) b(mda)
  real ( kind = 8 ) dummy
  real ( kind = 8 ) gen
  real ( kind = 8 ) h(mda)
  integer i
  integer j
  integer jm1
  integer k
  integer l
  integer m
  integer mmn
  integer mn1
  integer mn2
  integer n
  integer noise
  integer np1
  integer npj
  real ( kind = 8 ) sm
  real ( kind = 8 ) srsmsq

  300 format (/'  Estimated parameters,  x = a**(+)*b,', &
      ' computed by ''hft,hs1'''// &
      (9x,i6,e16.8,i6,e16.8,i6,e16.8,i6,e16.8,i6,e16.8))
  305 format (/'  Residual length  = ',e12.4)
  310 format (/'  Covariance matrix (unscaled) of', &
      ' estimated parameters'/9x,'computed by ''cov''.'/1x)
  320 format (1x,2i3,e16.8,2i3,e16.8,2i3,e16.8,2i3,e16.8,2i3,e16.8)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  HFT factors a least squares problem;'
  write ( *, '(a)' ) '  HS1 solves a factored least squares problem;'
  write ( *, '(a)' ) '  COV computes the associated covariance matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  No checking will be made for rank deficiency'
  write ( *, '(a)' ) '  in this test.  Such checks are made in later tests.'

  do noise = 1, 2
!
!  Initialize the data generation function.
!
    anoise = -1.0D+00
    dummy = gen ( anoise )

    if ( noise == 1 ) then
      anoise = 0.0D+00
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  No noise used in matrix generation.'
    else
      anoise = 0.0001D+00
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Matrix generation noise level = ', anoise
    end if

    do mn1 = 1, 6, 5

      mn2 = mn1+2

      do m = mn1, mn2

        do n = mn1, m

          np1 = n+1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '   M   N'
          write ( *, '(2i4)' ) m, n

          do i = 1, m
            do j = 1, n
              a(i,j) = gen(anoise)
            end do
          end do

          do i = 1,m
            b(i) = gen(anoise)
          end do

          if ( m < n ) go to 180
!
!  Apply algorithm HFT to factor the matrix.
!
          do j = 1, n
            call h12 ( 1, j, j+1, m, a(1,j), 1, h(j), a(1,min(j+1,n)), 1, &
              mda, n-j )
          end do
!
!  Apply algorithm hs1
!  Apply the transformations  q(n)...q(1) = q to b
!  replacing the previous contents of the array, b .
!
              do j = 1, n
                call h12 ( 2, j, j+1, m, a(1,j), 1, h(j), b, 1, 1, 1 )
              end do
!
!  Solve the triangular system for the solution X.
!  Store X in the array B.
!
              do k = 1, n

                 i = np1-k

                 sm = dot_product ( a(i,i+1:n), b(i+1:n) )

                 if ( a(i,i) == 0.0D+00 ) then
                   write ( *, '(a)' ) ' '
                   write ( *, '(a)' ) '  Terminating this case.'
                   write ( *, '(a)' ) '  A divisor is exactly zero.'
                   go to 180
                 end if

                 b(i) = ( b(i) - sm ) / a(i,i)

              end do
!
!  Compute the norm of the residual.
!
              srsmsq = 0.0
              do j = 1, m-n
                srsmsq = srsmsq+b(n+j)**2
              end do
              srsmsq = sqrt ( srsmsq )
!
!  Begin algorithm COV.
!  Compute unscaled covariance matrix   ((a**t)*a)**(-1)
!
              do j = 1, n
                a(j,j) = 1.0D+00 / a(j,j)
              end do

              do i = 1, n-1
                do j = i+1,n
                  sm = 0.0D+00
                  do l = i,j-1
                    sm = sm+a(i,l)*dble(a(l,j))
                  end do
                  a(i,j) = -sm*a(j,j)
                end do
              end do
!
!  The upper triangle of A has been inverted upon itself.
!
              do i = 1, n
                do j = i, n
                  sm = 0.0D+00
                  do l = j, n
                    sm = sm + a(i,l) * dble(a(j,l))
                  end do
                  a(i,j) = sm
                end do
              end do
!
!  the upper triangular part of the
!  symmetric matrix (a**t*a)**(-1) has
!  replaced the upper triangular part of
!  the a array.
!
              write (*,300) (i,b(i),i = 1,n)
              write (*,305) srsmsq
              write (*,310)
              do i = 1, n
                write (*,320) (i,j,a(i,j),j = i,n)
              end do

180           continue

        end do

      end do

    end do

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 demonstrates algorithms HFTI and COV.
!
!  Discussion:
!
!    Algorithm HFTI solves least squares problems.
!    Algorithm COV computes the associated unscaled covariance matrix.
!
!  Modified:
!
!    20 April 2004
!
!  Author:
!
!    Charles Lawson, Richard Hanson
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
  implicit none

  integer, parameter :: mda = 8
  integer, parameter :: mdb = 8
  integer, parameter :: nb = 1

  real ( kind = 8 ) a(mda,mda)
  real ( kind = 8 ) anoise
  real ( kind = 8 ) anorm
  real ( kind = 8 ) b(mdb,nb)
  real ( kind = 8 ) dummy
  real ( kind = 8 ) g(mda)
  real ( kind = 8 ) gen
  real ( kind = 8 ) h(mda)
  integer i
  integer ii
  integer ip(mda)
  integer j
  integer jm1
  integer k
  integer km1
  integer kp1
  integer krank
  integer l
  integer m
  integer mn1
  integer n
  integer noise
  real ( kind = 8 ) sm
  real ( kind = 8 ) srsmsq(nb)
  real ( kind = 8 ) tau
  real ( kind = 8 ) temp

  220 format (9x,2i3,e16.8,2i3,e16.8,2i3,e16.8,2i3,e16.8,2i3,e16.8)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Demonstrate the algorithms HFTI and COV.'

  do noise = 1, 2
!
!  Initialize the data generation function
!
    anoise = -1.0D+00
    dummy = gen ( anoise )

    anorm =  500.0D+00

    if ( noise == 1 ) then
      anoise = 0.0D+00
      tau = 0.5D+00
    else if ( noise == 2 ) then
      anoise = 0.0001D+00
      tau = anorm * anoise * 10.0D+00
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Use a relative noise level of        ', anoise
    write ( *, '(a,g14.6)' ) '  The matrix norm is approximately     ', anorm
    write ( *, '(a,g14.6)' ) '  The absolute pseudorank tolerance is ', tau

    do mn1 = 1, 6, 5

      do m = mn1, mn1 + 2

        do n = mn1, mn1 + 2

          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) '  M = ', m
          write ( *, '(a,i6)' ) '  N = ', n

          do i = 1, m
            do j = 1, n
              a(i,j) = gen ( anoise )
            end do
          end do

          do i = 1, m
            b(i,1) = gen ( anoise )
          end do

          call hfti ( a, mda, m, n, b, mdb, nb, tau, krank, srsmsq, h, g, ip )

          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' ) '  Pseudorank = ', krank
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Estimated parameters X = A**(+)*B from HFTI:'
          write ( *, '(a)' ) ' '
          write ( *, '(5g14.6)' ) b(1:n,1)
          write ( *, '(a)' ) ' '
          write ( *, '(a,g14.6)' ) '  Residual norm = ', srsmsq

          if ( n <= krank ) then
!
!  Algorithm COV begins here
!
          do j = 1, n
            a(j,j) = 1.0D+00 / a(j,j)
          end do

          do i = 1, n - 1
            do j = i + 1, n
              sm =  0.0D+00
              do l = i, j-1
                sm = sm + a(i,l) * a(l,j)
              end do
              a(i,j) = - sm * a(j,j)
            end do
          end do
!
!  The upper triangle of A has been inverted upon itself.
!
          do i = 1,n

            do j = i,n
              sm = 0.0D+00
              do l = j,n
                sm = sm + a(i,l) * dble ( a(j,l) )
              end do
              a(i,j) = sm
              end do
            end do

            do ii = 2,n

              i = n+1-ii

              if ( ip(i) /= i ) then

                k = ip(i)

                temp   = a(i,i)
                a(i,i) = a(k,k)
                a(k,k) = temp

                do l = 2, i
                  temp     = a(l-1,i)
                  a(l-1,i) = a(l-1,k)
                  a(l-1,k) = temp
                end do

                do l = i+1, k-1
                  temp   = a(i,l)
                  a(i,l) = a(l,k)
                  a(l,k) = temp
                end do

                do l = k+1, n
                  temp   = a(i,l)
                  a(i,l) = a(k,l)
                  a(k,l) = temp
                end do

              end if

            end do
!
!  The covariance has been computed and repermuted.
!  The upper triangular part of the symmetric matrix (a**t*a)**(-1) has
!  replaced the upper triangular part of the A array.
!
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Unscaled covariance matrix of estimated '
          write ( *, '(a)' ) '  parameters computed by COV'
          write ( *, '(a)' ) ' '
          do i = 1,n
            write (*,220) (i,j,a(i,j),j = i,n)
          end do

          end if

        end do
      end do
    end do
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 demonstrates the use of SVDRS,
!
!  Discussion:
!
!    SVDRS computes the singular value decomposition of a matrix A and
!    solves the least squares problem A*X = B.
!
!
!     the SVD  a =  u*(s,0)**t*v**t is computed so that..
!     (1)  u**t*b replaces the m+1 st. column of  b.
!
!     (2)  u**t   replaces the m by m identity in
!     the first m columns of B.
!
!     (3) v replaces the first n rows and columns of a.
!
!     (4) the diagonal entries of the S matrix
!     replace the first n entries of the array S.
!
!     The array S must be dimensioned at least 3*N.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
  implicit none

  integer, parameter :: mda = 8
  integer, parameter :: mdb = 8
  integer, parameter :: mx = 8

  real ( kind = 8 ) a(mda,mx)
  real ( kind = 8 ) aa(mda,mx)
  real ( kind = 8 ) anoise
  real ( kind = 8 ) b(mdb,mdb+1)
  real ( kind = 8 ) dn
  real ( kind = 8 ) dummy
  real ( kind = 8 ) floatn
  real ( kind = 8 ) gen
  integer i
  integer j
  integer kp1
  integer krank
  integer l
  integer m
  integer minmn
  integer mn1
  integer mn2
  integer n
  integer noise
  real ( kind = 8 ) rho
  real ( kind = 8 ) s(mx)
  real ( kind = 8 ) sm
  real ( kind = 8 ) srsmsq
  real ( kind = 8 ) t
  real ( kind = 8 ) tau
  real ( kind = 8 ) work(mx,2)
  real ( kind = 8 ) x(mx)

  220 format (/'  frobenius norm(a-u*(s,0)**t*v**t)'/ &
      '  ---------------------------------  = ',g12.4/ &
     '    sqrt(n) * spectral norm of a')
  230 format (2x,i5,e16.8,i5,e16.8,i5,e16.8)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Demonstrate the use of SVDRS.'

  do noise = 1, 2

    anoise = -1.0D+00
    dummy = gen ( anoise )

     if ( noise == 1 ) then
       anoise = 0.0D+00
       rho = 0.001D+00
     else
       anoise = 0.0001D+00
       rho = 10.0D+00 * anoise
     end if

     write ( *, '(a)' ) ' '
     write ( *, * ) '  The relative noise level in the data will be ', anoise
     write ( *, * ) '  The relative pseudorank tolerance will be ', rho

     do mn1 = 1, 6, 5
        mn2 = mn1+2
        do m = mn1, mn1 + 2
           do n = mn1, mn1 + 2

              write ( *, '(a)' ) ' '
              write ( *, '(a,i6)' ) '  M = ', m
              write ( *, '(a,i6)' ) '  N = ', n

              b(1:m,1:m) = 0.0D+00
              do i = 1,m
                b(i,i) = 1.0D+00
              end do

              do i = 1, m
                do j = 1,n
                  a(i,j) = gen(anoise)
                end do
              end do

              aa(1:m,1:n) = a(1:m,1:n)

              do i = 1, m
                b(i,m+1) = gen(anoise)
              end do
!
!  Compute the singular value decomposition.
!
              call svdrs ( a, mda, m, n, b, mdb, m+1, s )

              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) '  The singular value vector S:'
              write ( *, '(a)' ) ' '
              write ( *, '(4x,5f14.6)' ) s(1:n)
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) '  Transformed right side, U''*B'
              write ( *, '(a)' ) ' '
              write (*,230) (i,b(i,m+1),i = 1,m)
!
!  Test for disparity of ratio of singular values.
!
              krank = n
              tau = rho * s(1)
              do i = 1, n
                if ( s(i) <= tau ) then
                  krank = i - 1
                  exit
                end if
              end do

              write ( *, '(a)' ) ' '
              write ( *, * ) '  Absolute pseudorank tolerance TAU = ', tau
              write ( *, * ) '  Pseudorank = ', krank
!
!  Compute solution vector assuming pseudorank is KRANK.
!
              do i = 1, krank
                b(i,m+1) = b(i,m+1) / s(i)
              end do

              do i = 1, n
                x(i) = dot_product ( a(i,1:krank), b(1:krank,m+1) )
              end do
!
!  Compute predicted norm of residual vector.
!
              srsmsq = 0.0D+00
              if ( krank /= m ) then
                 do i = krank+1, m
                   srsmsq = srsmsq + b(i,m+1)**2
                 end do
                 srsmsq = sqrt ( srsmsq )
              end if

              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) '  Estimated X = A**(+) * B, computed by SVDRS:'
              write ( *, '(a)' ) ' '
              write ( *, '(5g14.6)' ) x(1:n)
              write ( *, '(a)' ) ' '
              write ( *, * ) '  Residual norm = ', srsmsq
!
!  Compute the frobenius norm of  a' - v * (s,0) * u'.
!
!  Compute v*s first.
!
              minmn = min(m,n)
              do j = 1,minmn
                a(1:n,j) = a(1:n,j) * s(j)
              end do

              dn = 0.0D+00
              do j = 1,m
                 do i = 1,n
                    sm = 0.0D+00
                    do l = 1,minmn
                      sm = sm+a(i,l)*b(l,j)
                    end do
!
!  Computed difference of (i,j) th entry of  a' - v * (s,0) * u'.
!
                    t = aa(j,i)-sm
                    dn = dn + t * t
                  end do
              end do

!             write ( *, * ) '  S(1) = ', s(1)

              dn = sqrt ( dn ) / ( sqrt ( dble ( n ) ) * s(1) )

              write (*,220) dn
        end do
      end do
    end do
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 demonstrates SVA, which carries out singular value analysis.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
  implicit none

  integer, parameter :: mda = 15
  integer, parameter :: mx = 5

  real ( kind = 8 ) a(mda,mx)
  real ( kind = 8 ) b(mda)
  real ( kind = 8 ) d(mx)
  character ( len = 80 ) :: file_name = 'lawson_prb_input.txt'
  integer i
  integer j
  integer, save, dimension ( 4 ) :: kpvec = (/ 1, 111111, -1, 76 /)
  integer m
  integer n
  character ( len = 8 ), dimension ( mx ) :: names = (/ &
    '  fire', ' water', ' earth', '   air', 'cosmos' /)
  real ( kind = 8 ) sing(mx)
  real ( kind = 8 ) work(2*mx)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Demonstrate singular value analysis.'
  write ( *, '(a)' ) '  Read data from file ' // trim ( file_name )
  write ( *, '(a)' ) '  Then call SVA.'

  m = mda
  n = mx

  open ( unit = 10, file = file_name, status= 'old')

  read (10,'(6f12.0)') ((a(i,j),j = 1,n),b(i),i=1,m)

  close ( unit = 10 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Listing of input matrix, a, and vector, b, follows..'
  write ( *, '(/(1x,f11.7,4f12.7,f13.4))') ((a(i,j),j = 1,n),b(i),i=1,m)
  write ( *, '(a)' ) ' '

  call sva ( a, mda, m, n, mda, b, sing, kpvec, names, 1, d, work )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests BNDACC and BNDSOL.
!
!  Discussion:
!
!    BNDACC and BNDSOL solve sequentially the banded least squares problem
!    that arises in spline curve fitting.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
  implicit none

  integer, parameter :: mdg = 12
  integer, parameter :: mxy = 12
  integer, parameter :: nband = 4

  real ( kind = 8 ) b(10)
  real ( kind = 8 ) c(mxy)
  real ( kind = 8 ) cov(mxy)
  real ( kind = 8 ) g(mdg,nband+1)
  real ( kind = 8 ) h
  integer i
  integer ic
  integer ig
  integer ip
  integer ir
  integer j
  integer jt
  integer l
  integer m
  integer mt
  integer nbp
  integer nc
  real ( kind = 8 ) p1
  real ( kind = 8 ) p2
  real ( kind = 8 ) q(4)
  real ( kind = 8 ) r
  real ( kind = 8 ) rdummy
  real ( kind = 8 ) rnorm
  real ( kind = 8 ) sigfac
  real ( kind = 8 ) sigsq
  real ( kind = 8 ) u
  real ( kind = 8 ) x(mxy)
  real ( kind = 8 ), dimension ( mxy ) :: y = (/ &
    2.2D+00, 4.0D+00, 5.0D+00, 4.6D+00, 2.8D+00, 2.7D+00, &
    3.8D+00, 5.1D+00, 6.1D+00, 6.3D+00, 5.0D+00, 2.0D+00 /)
  real ( kind = 8 ) yfit

  160 format (/'   i',8x,'x',10x,'y',6x,'yfit',4x,'r = y-yfit/1x')
  170 format (1x,i3,4x,f6.0,4x,f6.2,4x,f6.2,4x,f8.4)

  190 format (3(2x,2i4,g15.7))
  200 format (/' covariance matrix of the spline coefficients.')

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  BNDACC accumulates a banded matrix.'
  write ( *, '(a)' ) '  BNDSOL solves an associated banded least'
  write ( *, '(a)' ) '    squares problem.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Execute a sequence of cubic spline fits'
  write ( *, '(a)' ) '  to a discrete set of data.'

  m = mdg
!
!  Set the data points.
!
  do i = 1, m
    x(i) = dble ( 2 * i )
  end do
!
!  Loop thru cases using increasing numbers of breakpoints.
!
  do nbp = 5, 10

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The number of breakpoints is ', nbp

    nc = nbp + 2
!
!  Set the breakpoints.
!
    do i = 1, nbp
      b(i) = ( dble ( nbp - i ) * x(1) + dble ( i - 1 ) * x(m) ) &
        / dble ( nbp - 1 )
    end do

    h = ( x(m) - x(1) ) / dble ( nbp - 1 )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Breakpoints:'
    write ( *, '(a)' ) ' '
    write ( *, '(5g14.6)' ) b(1:nbp)
!
!  Initialize IR and IP before first call to BNDACC.
!
    ir = 1
    ip = 1
    i = 1
    jt = 1
!
!  Set row for I-th data point.
!
    do

      mt = 0

      do while ( x(i) <= b(jt+1) )

        u = ( x(i) - b(jt) ) / h
        ig = ir + mt

        g(ig,1) = p1 ( 1.0D+00 - u )
        g(ig,2) = p2 ( 1.0D+00 - u )
        g(ig,3) = p2 ( u )
        g(ig,4) = p1 ( u )
        g(ig,5) = y(i)

        mt = mt + 1

        if ( i == m ) then
          exit
        end if

        i = i + 1

      end do
!
!  Send the block of data to the band matrix accumulator.
!
      call bndacc ( g, mdg, nband, ip, ir, mt, jt )

      if ( i == m ) then
        exit
      end if

      jt = jt + 1

    end do
!
!  Compute the solution C.
!
    call bndsol ( 1, g, mdg, nband, ip, ir, c, nc, rnorm )
!
!  Print the solution.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C:'
    write ( *, '(a)' ) ' '
    write ( *, '(5(2x,g14.6))' ) c(1:nc)
!
!  Print the norm of the residual.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  RNORM = ', rnorm
!
!  Compute and print x,y,yfit,r = y-yfit
!
    write (*,160)
    jt = 1
    do i = 1, m

      do while ( x(i) > b(jt+1) )
        jt = jt + 1
      end do

      u = ( x(i) - b(jt) ) / h
      q(1) = p1 ( 1.0D+00 - u )
      q(2) = p2 ( 1.0D+00 - u )
      q(3) = p2(u)
      q(4) = p1(u)
      yfit = 0.0D+00

      do l = 1, 4
        ic = jt - 1 + l
        yfit = yfit + c(ic) * q(l)
      end do

      r = y(i) - yfit
      write (*,170) i,x(i),y(i),yfit,r

    end do
!
!  Compute the residual vector norm.
!
    if ( m > nc ) then
      sigsq = rnorm**2 / ( m - nc )
      sigfac = sqrt ( sigsq )
      write ( *, '(a,g14.6)' ) '  SIGFAC = ', sigfac
      write (*,200)
!
!  Compute and print columns of covariance.
!
      do j = 1, nc

        cov(1:nc) = 0.0D+00
        cov(j) = 1.0D+00

        call bndsol ( 2, g, mdg, nband, ip, ir, cov, nc, rdummy )

        call bndsol ( 3, g, mdg, nband, ip, ir, cov, nc, rdummy )
!
!  Compute the J-th column of the covariance matrix.
!
        cov(1:nc) = cov(1:nc) * sigsq

        write (*,190) (l,j,cov(l),l = 1,nc)

      end do

    end if

  end do

  return
end
function p1 ( t )

!*****************************************************************************80
!
!! P1 ???
!
!  Modified:
!
!    22 October 2008
!
  implicit none

  real ( kind = 8 ) p1
  real ( kind = 8 ) t

  p1 = 0.25D+00 * t**3

  return
end
function p2 ( t )

!*****************************************************************************80
!
!! P2 ???
!
!  Modified:
!
!    22 October 2008
!
  implicit none

  real ( kind = 8 ) p2
  real ( kind = 8 ) t

  p2 = -(1.0D+00 - t)**2 * ( 1.0D+00 + t ) * 0.75D+00 + 1.0D+00

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 demonstrate the use of LDP.
!
!  Discussion:
!
!    LDP is used for least distance programming.
!
!    This test example solves the constrained line data fitting problem
!    of chapter 23.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
  implicit none

  integer, parameter :: mde = 4
  integer, parameter :: me = 4
  integer, parameter :: mdgh = 3
  integer, parameter :: mg = 3
  integer, parameter :: mx = 2
  integer, parameter :: n = 2

  real ( kind = 8 ) e(mde,n)
  real ( kind = 8 ) f(me,1)
  real ( kind = 8 ) g(mdgh,mx)
  real ( kind = 8 ) g2(mdgh,mx)
  real ( kind = 8 ), dimension ( mg ) :: h = (/ 0.0D+00, 0.0D+00, -1.0D+00 /)
  real ( kind = 8 ) h2(mg)
  integer i
  integer indx(mg)
  integer j
  integer l
  integer mode
  integer np1
  real ( kind = 8 ) res
  real ( kind = 8 ) s(mx)
  real ( kind = 8 ) sm
  real ( kind = 8 ), dimension ( me ) :: t = (/ &
    0.25D+00, 0.5D+00, 0.5D+00, 0.8D+00 /)
  real ( kind = 8 ), dimension ( me ) :: w = (/ &
    0.50D+00, 0.6D+00, 0.7D+00, 1.2D+00 /)
  real ( kind = 8 ) wldp((n+1)*(mg+2)+2*mg)
  real ( kind = 8 ) work(2*mx)
  real ( kind = 8 ) x(mx)
  real ( kind = 8 ) z(mx)
  real ( kind = 8 ) znorm

  200 format (' mode (from ldp) = ',i3,',  znorm = ',f10.5)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  LDP carries out least distance programming.'
!
!  Define the least squares and constraint matrices.
!
  e(1:me,1) = t(1:me)
  e(1:me,2) = 1.0D+00
  f(1:me,1) = w(1:me)

  g(1,1) = 1.0D+00
  g(1,2) = 0.0D+00
  g(2,1) = 0.0D+00
  g(2,2) = 1.0D+00
  g(3,1) = -1.0D+00
  g(3,2) = -1.0D+00
!
!  Compute the singular value decomposition of the matrix E.
!
  call svdrs ( e, mde, me, n, f, me, 1, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  V:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,2f10.5)' ) (e(i,1:n),i=1,n)
  write ( *, '(a)' ) ' '
  write ( *, '(4x,2f10.5)' ) f(1:me,1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Singular value vector S:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,2f10.5)' ) s(1:n)
!
!  Generally, rank determination and Levenberg-Marquardt
!  stabilization could be inserted here.
!
!  Define the constraint matrix for the Z coordinate system.
!
  do i = 1, mg
    do j = 1, n
      sm = 0.0D+00
      do l = 1, n
        sm = sm + g(i,l) * e(l,j)
      end do
      g2(i,j) = sm / s(j)
    end do
  end do
!
!  Define the constraint right hand side for the Z coordinate system.
!
  do i = 1, mg
    h2(i) = h(i) - dot_product ( g2(i,1:n), f(1:n,1) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  G tilde = '
  write ( *, '(2f10.5)') (g2(i,1:n),i=1,mg)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  H tilde = '
  write ( *, '(3f10.5)' ) h2(1:mg)
!
!  Solve the constrained problem in Z coordinates.
!
  call ldp ( g2, mdgh, mg, n, h2, z, znorm, wldp, indx, mode )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  LDP MODE = ', mode
  write ( *, '(a,g14.6)' ) '  ZNORM = ', znorm

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Z = '
  write ( *, '(a)' ) ' '
  write ( *, '(2g14.6)' ) z(1:mx)
!
!  Transform back from Z coordinates to X coordinates.
!
  z(1:n) = ( z(1:n) + f(1:n,1) ) / s(1:n)

  x(1:n) = matmul ( e(1:n,1:n), z(1:n) )

  res = sqrt ( znorm**2 + sum ( f(n+1:me,1)**2 ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The coefficients of the fitted line'
  write ( *, '(a)' ) '    F(T) = X(1) * T + X(2)'
  write ( *, '(a)' ) '  are:'
  write ( *, '(5g14.6)' ) x(1:n)
!
!  Compute the residuals.
!
  f(1:me,1) = w(1:me) - x(1) * t(1:me) - x(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The residual vector:'
  write ( *, '(a)' ) ' '
  do i = 1, me
    write ( *, '(g14.6)' ) f(i,1)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The residual norm: '
  write ( *, '(a)' ) ' '
  write ( *, '(g14.6)' ) res

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 demonstrates the use of QRBD.
!
!  Modified:
!
!    22 October 2008
!
!  Author:
!
!    Charles Lawson, Richard Hanson
!
!  Reference:
!
!    Charles Lawson, Richard Hanson,
!    Solving Least Squares Problems,
!    SIAM, 1995,
!    ISBN: 0898713560,
!    LC: QA275.L38.
!
  implicit none

  integer, parameter :: n = 3

  real ( kind = 8 ) bd(n,n)
  real ( kind = 8 ) e(n)
  integer i
  integer ipass
  real ( kind = 8 ) q(n)
  real ( kind = 8 ) s(n,n)
  real ( kind = 8 ) u(n,n)
  real ( kind = 8 ) v(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  QRBD computes the singular values S of a bidiagonal'
  write ( *, '(a)' ) '    matrix BD, and can also compute the decomposition'
  write ( *, '(a)' ) '    factors U and V, so that'
  write ( *, '(a)' ) '      S = U * BD * V.'

  call random_number ( harvest = q(1:n) )
  call random_number ( harvest = e(2:n) )

  bd(1:n,1:n) = 0.0D+00
  do i = 1, n
    bd(i,i) = q(i)
  end do
  do i = 2, n
    bd(i-1,i) = e(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The bidiagonal matrix BD:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(5f12.4)' ) bd(i,1:n)
  end do

  u(1:n,1:n) = 0.0D+00
  v(1:n,1:n) = 0.0D+00
  do i = 1, n
    u(i,i) = 1.0D+00
    v(i,i) = 1.0D+00
  end do

  call qrbd ( ipass, q, e, n, v, n, n, u, n, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Error flag IPASS = ', ipass
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The singular values of BD:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, q(i)
  end do

  s(1:n,1:n) = 0.0D+00
  do i = 1, n
    s(i,i) = q(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The factor U:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(5f12.4)' ) u(i,1:n)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The factor V:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(5f12.4)' ) v(i,1:n)
  end do

  bd(1:n,1:n) = matmul ( matmul ( transpose ( u(1:n,1:n) ), s(1:n,1:n) ), &
    transpose ( v(1:n,1:n) ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The product U'' * S * V'' = BD:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(5f12.4)' ) bd(i,1:n)
  end do

  return
end
