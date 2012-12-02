!!
!                                                                         !
!        N  A  S     P A R A L L E L     B E N C H M A R K S  3.3         !
!                                                                         !
!                      S E R I A L     V E R S I O N                      !
!                                                                         !
!                                   C G                                   !
!                                                                         !
!!
!                                                                         !
!    This benchmark is a serial version of the NPB CG code.               !
!    Refer to NAS Technical Reports 95-020 for details.                   !
!                                                                         !
!    Permission to use, copy, distribute and modify this software         !
!    for any purpose with or without fee is hereby granted.  We           !
!    request, however, that all derived work reference the NAS            !
!    Parallel Benchmarks 3.3. This software is provided "as is"           !
!    without express or implied warranty.                                 !
!                                                                         !
!    Information on NPB 3.3, including the technical report, the          !
!    original specifications, source code, results and information        !
!    on how to submit new results, is available at:                       !
!                                                                         !
!           http://www.nas.nasa.gov/Software/NPB/                         !
!                                                                         !
!    Send comments or suggestions to  npb@nas.nasa.gov                    !
!                                                                         !
!          NAS Parallel Benchmarks Group                                  !
!          NASA Ames Research Center                                      !
!          Mail Stop: T27A-1                                              !
!          Moffett Field, CA   94035-1000                                 !
!                                                                         !
!          E-mail:  npb@nas.nasa.gov                                      !
!          Fax:     (650) 604-3957                                        !
!                                                                         !
!!


!
!      NPB CG serial version
!
! Authors: M. Yarrow
!          C. Kuszmaul
!
!
!
!
      program cg

!*****************************************************************************80


      implicit none

      include 'globals.h'


      common / main_int_mem /  colidx,     rowstr, &
                               iv,         arow,     acol
      integer                  colidx(nz), rowstr(na+1), &
                               iv(na),  arow(na), acol(naz)


      common / main_flt_mem /  aelt,     a, &
                               x, &
                               z, &
                               p, &
                               q, &
                               r
      double precision         aelt(naz), a(nz), &
                               x(na+2), &
                               z(na+2), &
                               p(na+2), &
                               q(na+2), &
                               r(na+2)



      integer            i, j, k, it

      double precision   zeta, randlc
      external           randlc
      double precision   rnorm
      double precision   norm_temp1,norm_temp2

      double precision   t, mflops, tmax
      character          class
      logical            verified
      double precision   zeta_verify_value, epsilon, err

      integer   fstatus
      character t_names(t_last)*8

      do i = 1, T_last
         call timer_clear( i )
      end do

      open(unit=2, file='timer.flag', status='old', iostat=fstatus)
      if (fstatus .eq. 0) then
         timeron = .true.
         t_names(t_init) = 'init'
         t_names(t_bench) = 'benchmk'
         t_names(t_conj_grad) = 'conjgd'
         close(2)
      else
         timeron = .false.
      endif

      call timer_start( T_init )

      firstrow = 1
      lastrow  = na
      firstcol = 1
      lastcol  = na


      if( na .eq. 1400 .and. &
          nonzer .eq. 7 .and. &
          niter .eq. 15 .and. &
          shift .eq. 10.d0 ) then
         class = 'S'
         zeta_verify_value = 8.5971775078648d0
      else if( na .eq. 7000 .and. &
               nonzer .eq. 8 .and. &
               niter .eq. 15 .and. &
               shift .eq. 12.d0 ) then
         class = 'W'
         zeta_verify_value = 10.362595087124d0
      else if( na .eq. 14000 .and. &
               nonzer .eq. 11 .and. &
               niter .eq. 15 .and. &
               shift .eq. 20.d0 ) then
         class = 'A'
         zeta_verify_value = 17.130235054029d0
      else if( na .eq. 75000 .and. &
               nonzer .eq. 13 .and. &
               niter .eq. 75 .and. &
               shift .eq. 60.d0 ) then
         class = 'B'
         zeta_verify_value = 22.712745482631d0
      else if( na .eq. 150000 .and. &
               nonzer .eq. 15 .and. &
               niter .eq. 75 .and. &
               shift .eq. 110.d0 ) then
         class = 'C'
         zeta_verify_value = 28.973605592845d0
      else if( na .eq. 1500000 .and. &
               nonzer .eq. 21 .and. &
               niter .eq. 100 .and. &
               shift .eq. 500.d0 ) then
         class = 'D'
         zeta_verify_value = 52.514532105794d0
      else if( na .eq. 9000000 .and. &
               nonzer .eq. 26 .and. &
               niter .eq. 100 .and. &
               shift .eq. 1.5d3 ) then
         class = 'E'
         zeta_verify_value = 77.522164599383d0
      else
         class = 'U'
      endif

      write( *,1000 )
      write( *,1001 ) na
      write( *,1002 ) niter
      write( *,* )
 1000 format(//,' NAS Parallel Benchmarks (NPB3.3-SER)', &
                ' - CG Benchmark', /)
 1001 format(' Size: ', i11 )
 1002 format(' Iterations: ', i5 )



      naa = na
      nzz = nz


!
!  Inialize random number generator
!
      tran    = 314159265.0D0
      amult   = 1220703125.0D0
      zeta    = randlc( tran, amult )

!
!
!
      call makea(naa, nzz, a, colidx, rowstr, &
                 firstrow, lastrow, firstcol, lastcol, &
                 arow, acol, aelt, iv)



!
!  Note: as a result of the above call to makea:
!        values of j used in indexing rowstr go from 1 --> lastrow-firstrow+1
!        values of colidx which are col indexes go from firstcol --> lastcol
!        So:
!        Shift the col index vals from actual (firstcol --> lastcol )
!        to local, i.e., (1 --> lastcol-firstcol+1)
!
      do j=1,lastrow-firstrow+1
         do k=rowstr(j),rowstr(j+1)-1
            colidx(k) = colidx(k) - firstcol + 1
         end do
      end do

!
!  set starting vector to (1, 1, .... 1)
!
      do i = 1, na+1
         x(i) = 1.0D0
      end do
      do j=1, lastcol-firstcol+1
         q(j) = 0.0d0
         z(j) = 0.0d0
         r(j) = 0.0d0
         p(j) = 0.0d0
      end do

      zeta  = 0.0d0

!
!
!  Do one iteration untimed to init all code and data page tables
!                    (then reinit, start timing, to niter its)
!
      do it = 1, 1

!
!  The call to the conjugate gradient routine:
!
         call conj_grad ( colidx, &
                          rowstr, &
                          x, &
                          z, &
                          a, &
                          p, &
                          q, &
                          r, &
                          rnorm )

!
!  zeta = shift + 1/(x.z)
!  So, first: (x.z)
!  Also, find norm of z
!  So, first: (z.z)
!
         norm_temp1 = 0.0d0
         norm_temp2 = 0.0d0
         do j=1, lastcol-firstcol+1
            norm_temp1 = norm_temp1 + x(j)*z(j)
            norm_temp2 = norm_temp2 + z(j)*z(j)
         end do

         norm_temp2 = 1.0d0 / sqrt( norm_temp2 )


!
!  Normalize z to obtain x
!
         do j=1, lastcol-firstcol+1
            x(j) = norm_temp2*z(j)
         end do


      end do                              ! end of do one iteration untim


!
!  set starting vector to (1, 1, .... 1)
!
!
!
!
      do i = 1, na+1
         x(i) = 1.0D0
      end do

      zeta  = 0.0d0

      call timer_stop( T_init )

      write (*, 2000) timer_read(T_init)
 2000 format(' Initialization time = ',f15.3,' seconds')

      call timer_start( T_bench )

!
!
!  Main Iteration for inverse power method
!
!
      do it = 1, niter

!
!  The call to the conjugate gradient routine:
!
         if ( timeron ) call timer_start( T_conj_grad )
         call conj_grad ( colidx, &
                          rowstr, &
                          x, &
                          z, &
                          a, &
                          p, &
                          q, &
                          r, &
                          rnorm )
         if ( timeron ) call timer_stop( T_conj_grad )


!
!  zeta = shift + 1/(x.z)
!  So, first: (x.z)
!  Also, find norm of z
!  So, first: (z.z)
!
         norm_temp1 = 0.0d0
         norm_temp2 = 0.0d0
         do j=1, lastcol-firstcol+1
            norm_temp1 = norm_temp1 + x(j)*z(j)
            norm_temp2 = norm_temp2 + z(j)*z(j)
         end do


         norm_temp2 = 1.0d0 / sqrt( norm_temp2 )


         zeta = shift + 1.0d0 / norm_temp1
         if( it .eq. 1 ) write( *,9000 )
         write( *,9001 ) it, rnorm, zeta

 9000    format( /,'   iteration           ||r||                 zeta' )
 9001    format( 4x, i5, 7x, e20.14, f20.13 )

!
!  Normalize z to obtain x
!
         do j=1, lastcol-firstcol+1
            x(j) = norm_temp2*z(j)
         end do


      end do                              ! end of main iter inv pow meth

      call timer_stop( T_bench )

!
!  End of timed section
!

      t = timer_read( T_bench )


      write(*,100)
 100  format(' Benchmark completed ')

      epsilon = 1.d-10
      if (class .ne. 'U') then

!         err = abs( zeta - zeta_verify_value)
         err = abs( zeta - zeta_verify_value )/zeta_verify_value
         if( err .le. epsilon ) then
            verified = .TRUE.
            write(*, 200)
            write(*, 201) zeta
            write(*, 202) err
 200        format(' VERIFICATION SUCCESSFUL ')
 201        format(' Zeta is    ', E20.13)
 202        format(' Error is   ', E20.13)
         else
            verified = .FALSE.
            write(*, 300)
            write(*, 301) zeta
            write(*, 302) zeta_verify_value
 300        format(' VERIFICATION FAILED')
 301        format(' Zeta                ', E20.13)
 302        format(' The correct zeta is ', E20.13)
         endif
      else
         verified = .FALSE.
         write (*, 400)
         write (*, 401)
         write (*, 201) zeta
 400     format(' Problem size unknown')
 401     format(' NO VERIFICATION PERFORMED')
      endif


      if( t .ne. 0. ) then
         mflops = float( 2*niter*na ) &
                     * ( 3.+float( nonzer*(nonzer+1) ) &
                       + 25.*(5.+float( nonzer*(nonzer+1) )) &
                       + 3. ) / t / 1000000.0
      else
         mflops = 0.0
      endif


         call print_results('CG', class, na, 0, 0, &
                            niter, t, &
                            mflops, '          floating point', &
                            verified, npbversion, compiletime, &
                            cs1, cs2, cs3, cs4, cs5, cs6, cs7)



 600  format( i4, 2e19.12)


!
!      More timers
!
      if (.not.timeron) goto 999

      tmax = timer_read(T_bench)
      if (tmax .eq. 0.0) tmax = 1.0

      write(*,800)
 800  format('  SECTION   Time (secs)')
      do i=1, t_last
         t = timer_read(i)
         if (i.eq.t_init) then
            write(*,810) t_names(i), t
         else
            write(*,810) t_names(i), t, t*100./tmax
            if (i.eq.t_conj_grad) then
               t = tmax - t
               write(*,820) 'rest', t, t*100./tmax
            endif
         endif
 810     format(2x,a8,':',f9.3:'  (',f6.2,'%)')
 820     format('    --> total ',a8,':',f9.3,'  (',f6.2,'%)')
      end do

 999  continue

      stop
      end
      subroutine conj_grad ( colidx, &
                             rowstr, &
                             x, &
                             z, &
                             a, &
                             p, &
                             q, &
                             r, &
                             rnorm )

!*****************************************************************************80
!
!  Floating point arrays here are named as in NPB1 spec discussion of
!  CG algorithm
!

      implicit none


      include 'globals.h'


      double precision   x(*), &
                         z(*), &
                         a(nzz)
      integer            colidx(nzz), rowstr(naa+1)

      double precision   p(*), &
                         q(*), &
                         r(*)


      integer   j, k
      integer   cgit, cgitmax

      double precision   d, sum, rho, rho0, alpha, beta, rnorm

      data      cgitmax / 25 /


      rho = 0.0d0

!
!  Initialize the CG algorithm:
!
      do j=1,naa+1
         q(j) = 0.0d0
         z(j) = 0.0d0
         r(j) = x(j)
         p(j) = r(j)
      end do


!
!  rho = r.r
!  Now, obtain the norm of r: First, sum squares of r elements locally...
!
      do j=1, lastcol-firstcol+1
         rho = rho + r(j)*r(j)
      end do

!
!  The conj grad iteration loop
!
      do cgit = 1, cgitmax

!
!  q = A.p
!  The partition submatrix-vector multiply: use workspace w
!
!
!  NOTE: this version of the multiply is actually (slightly: maybe %5)
!        faster on the sp2 on 16 nodes than is the unrolled-by-2 version
!        below.   On the Cray t3d, the reverse is true, i.e., the
!        unrolled-by-two version is some 10% faster.
!        The unrolled-by-8 version below is significantly faster
!        on the Cray t3d - overall speed of code is 1.5 times faster.
!
         do j=1,lastrow-firstrow+1
            sum = 0.d0
            do k=rowstr(j),rowstr(j+1)-1
               sum = sum + a(k)*p(colidx(k))
            end do
            q(j) = sum
         end do

!          do j=1,lastrow-firstrow+1
!             i = rowstr(j)
!             iresidue = mod( rowstr(j+1)-i, 2 )
!             sum1 = 0.d0
!             sum2 = 0.d0
!             if( iresidue .eq. 1 )
!      &          sum1 = sum1 + a(i)*p(colidx(i))
!             do k=i+iresidue, rowstr(j+1)-2, 2
!                sum1 = sum1 + a(k)  *p(colidx(k))
!                sum2 = sum2 + a(k+1)*p(colidx(k+1))
!             end do
!             q(j) = sum1 + sum2
!          end do

!          do j=1,lastrow-firstrow+1
!             i = rowstr(j)
!             iresidue = mod( rowstr(j+1)-i, 8 )
!             sum = 0.d0
!             do k=i,i+iresidue-1
!                sum = sum +  a(k)*p(colidx(k))
!             end do
!             do k=i+iresidue, rowstr(j+1)-8, 8
!                sum = sum + a(k  )*p(colidx(k  ))
!      &                   + a(k+1)*p(colidx(k+1))
!      &                   + a(k+2)*p(colidx(k+2))
!      &                   + a(k+3)*p(colidx(k+3))
!      &                   + a(k+4)*p(colidx(k+4))
!      &                   + a(k+5)*p(colidx(k+5))
!      &                   + a(k+6)*p(colidx(k+6))
!      &                   + a(k+7)*p(colidx(k+7))
!             end do
!             q(j) = sum
!          end do



!
!  Obtain p.q
!
         d = 0.0d0
         do j=1, lastcol-firstcol+1
            d = d + p(j)*q(j)
         end do


!
!  Obtain alpha = rho / (p.q)
!
         alpha = rho / d

!
!  Save a temporary of rho
!
         rho0 = rho

!
!  Obtain z = z + alpha*p
!  and    r = r - alpha*q
!
         rho = 0.0d0
         do j=1, lastcol-firstcol+1
            z(j) = z(j) + alpha*p(j)
            r(j) = r(j) - alpha*q(j)
         end do

!
!  rho = r.r
!  Now, obtain the norm of r: First, sum squares of r elements locally...
!
         do j=1, lastcol-firstcol+1
            rho = rho + r(j)*r(j)
         end do

!
!  Obtain beta:
!
         beta = rho / rho0

!
!  p = r + beta*p
!
         do j=1, lastcol-firstcol+1
            p(j) = r(j) + beta*p(j)
         end do


      end do                             ! end of do cgit=1,cgitmax


!
!  Compute residual norm explicitly:  ||r|| = ||x - A.z||
!  First, form A.z
!  The partition submatrix-vector multiply
!
      sum = 0.0d0
      do j=1,lastrow-firstrow+1
         d = 0.d0
         do k=rowstr(j),rowstr(j+1)-1
            d = d + a(k)*z(colidx(k))
         end do
         r(j) = d
      end do
!
!  At this point, r contains A.z
!
      do j=1, lastcol-firstcol+1
         d   = x(j) - r(j)
         sum = sum + d*d
      end do

      rnorm = sqrt( sum )

      return
      end
      subroutine makea( n, nz, a, colidx, rowstr, &
                        firstrow, lastrow, firstcol, lastcol, &
                        arow, acol, aelt, iv )

!*****************************************************************************80

      include             'npbparams.h'

      integer             n, nz
      integer             firstrow, lastrow, firstcol, lastcol
      integer             colidx(nz), rowstr(n+1)
      integer             iv(n), arow(n), acol(nonzer+1,n)
      double precision    aelt(nonzer+1,n)
      double precision    a(nz)
!
!       generate the test problem for benchmark 6
!       makea generates a sparse matrix with a
!       prescribed sparsity distribution
!
!       parameter    type        usage
!
!       input
!
!       n            i           number of cols/rows of matrix
!       nz           i           nonzeros as declared array size
!       rcond        r*8         condition number
!       shift        r*8         main diagonal shift
!
!       output
!
!       a            r*8         array for nonzeros
!       colidx       i           col indices
!       rowstr       i           row pointers
!
!       workspace
!
!       iv, arow, acol i
!       aelt           r*8
!

      integer          i, iouter, ivelt, nzv, nn1
      integer          ivc(nonzer+1)
      double precision vc(nonzer+1)

!
!      nonzer is approximately  (int(sqrt(nnza /n)));
!

      external          sparse, sprnvc, vecset

!
!    nn1 is the smallest power of two not less than n
!

      nn1 = 1
 50   continue
        nn1 = 2 * nn1
        if (nn1 .lt. n) goto 50

!
!  Generate nonzero positions and save for the use in sparse.
!

      do iouter = 1, n
         nzv = nonzer
         call sprnvc( n, nzv, nn1, vc, ivc )
         call vecset( n, vc, ivc, nzv, iouter, .5D0 )
         arow(iouter) = nzv
         do ivelt = 1, nzv
            acol(ivelt, iouter) = ivc(ivelt)
            aelt(ivelt, iouter) = vc(ivelt)
         end do
      end do

!
!       ... make the sparse matrix from list of elements with duplicates
!           (iv is used as  workspace)
!
      call sparse( a, colidx, rowstr, n, nz, nonzer, arow, acol, &
                   aelt, firstrow, lastrow, &
                   iv, rcond, shift )
      return

      end
      subroutine sparse( a, colidx, rowstr, n, nz, nonzer, arow, acol, &
                         aelt, firstrow, lastrow, &
                         nzloc, rcond, shift )

!*****************************************************************************80

      implicit           logical (a-z)
      integer            colidx(*), rowstr(*)
      integer            firstrow, lastrow
      integer            n, nz, nonzer, arow(*), acol(nonzer+1,*)
      double precision   a(*), aelt(nonzer+1,*), rcond, shift

!
!       rows range from firstrow to lastrow
!       the rowstr pointers are defined for nrows = lastrow-firstrow+1 values
!
      integer            nzloc(n), nrows

!
!       generate a sparse matrix from a list of
!       [col, row, element] tri
!
      integer            i, j, j1, j2, nza, k, kk, nzrow, jcol
      double precision   xi, size, scale, ratio, va

!
!    how many rows of result
!
      nrows = lastrow - firstrow + 1

!
!     ...count the number of triples in each row
!
      do j = 1, nrows+1
         rowstr(j) = 0
      end do

      do i = 1, n
         do nza = 1, arow(i)
            j = acol(nza, i) + 1
            rowstr(j) = rowstr(j) + arow(i)
         end do
      end do

      rowstr(1) = 1
      do j = 2, nrows+1
         rowstr(j) = rowstr(j) + rowstr(j-1)
      end do
      nza = rowstr(nrows+1) - 1

!
!     ... rowstr(j) now is the location of the first nonzero
!           of row j of a
!

      if (nza .gt. nz) then
         write(*,*) 'Space for matrix elements exceeded in sparse'
         write(*,*) 'nza, nzmax = ',nza, nz
         stop
      endif


!
!     ... preload data pages
!
      do j = 1, nrows
         do k = rowstr(j), rowstr(j+1)-1
             a(k) = 0.d0
             colidx(k) = 0
         end do
         nzloc(j) = 0
      end do

!
!     ... generate actual values by summing duplicates
!

      size = 1.0D0
      ratio = rcond ** (1.0D0 / dfloat(n))

      do i = 1, n
         do nza = 1, arow(i)
            j = acol(nza, i)

            scale = size * aelt(nza, i)
            do nzrow = 1, arow(i)
               jcol = acol(nzrow, i)
               va = aelt(nzrow, i) * scale

!
!       ... add the identity * rcond to the generated matrix to bound
!           the smallest eigenvalue from below by rcond
!
               if (jcol .eq. j .and. j .eq. i) then
                  va = va + rcond - shift
               endif

               do k = rowstr(j), rowstr(j+1)-1
                  if (colidx(k) .gt. jcol) then
!
!       ... insert colidx here orderly
!
                     do kk = rowstr(j+1)-2, k, -1
                        if (colidx(kk) .gt. 0) then
                           a(kk+1)  = a(kk)
                           colidx(kk+1) = colidx(kk)
                        endif
                     end do
                     colidx(k) = jcol
                     a(k)  = 0.d0
                     goto 40
                  else if (colidx(k) .eq. 0) then
                     colidx(k) = jcol
                     goto 40
                  else if (colidx(k) .eq. jcol) then
!
!       ... mark the duplicated entry
!
                     nzloc(j) = nzloc(j) + 1
                     goto 40
                  endif
               end do
               print *,'internal error in sparse: i=',i
               stop
   40          continue
               a(k) = a(k) + va
            end do
   60       continue
         end do
         size = size * ratio
      end do


!
!       ... remove empty entries and generate final results
!
      do j = 2, nrows
         nzloc(j) = nzloc(j) + nzloc(j-1)
      end do

      do j = 1, nrows
         if (j .gt. 1) then
            j1 = rowstr(j) - nzloc(j-1)
         else
            j1 = 1
         endif
         j2 = rowstr(j+1) - nzloc(j) - 1
         nza = rowstr(j)
         do k = j1, j2
            a(k) = a(nza)
            colidx(k) = colidx(nza)
            nza = nza + 1
         end do
      end do
      do j = 2, nrows+1
         rowstr(j) = rowstr(j) - nzloc(j-1)
      end do
      nza = rowstr(nrows+1) - 1


!C       write (*, 11000) nza
      return
11000   format ( //,'final nonzero count in sparse ', &
                  /,'number of nonzeros       = ', i16 )
      end
      subroutine sprnvc( n, nz, nn1, v, iv )

!*****************************************************************************80

      implicit           none
      double precision   v(*)
      integer            n, nz, nn1, iv(*)
      common /urando/    amult, tran
      double precision   amult, tran


!
!       generate a sparse n-vector (v, iv)
!       having nzv nonzeros
!
!       mark(i) is set to 1 if position i is nonzero.
!       mark is all zero on entry and is reset to all zero before exit
!       this corrects a performance bug found by John G. Lewis, caused by
!       reinitialization of mark on every one of the n calls to sprnvc
!

        integer            nzv, ii, i, icnvrt

        external           randlc, icnvrt
        double precision   randlc, vecelt, vecloc


        nzv = 0

100     continue
        if (nzv .ge. nz) goto 110

         vecelt = randlc( tran, amult )

!
!   generate an integer between 1 and n in a portable manner
!
         vecloc = randlc(tran, amult)
         i = icnvrt(vecloc, nn1) + 1
         if (i .gt. n) goto 100

!
!  was this integer generated already?
!
         do ii = 1, nzv
            if (iv(ii) .eq. i) goto 100
         end do
         nzv = nzv + 1
         v(nzv) = vecelt
         iv(nzv) = i
         goto 100
110     continue

      return
      end
      function icnvrt(x, ipwr2)

!*****************************************************************************80

      implicit           logical (a-z)
      double precision   x
      integer            ipwr2, icnvrt

!
!    scale a double precision number x in (0,1) by a power of 2 and chop it
!
      icnvrt = int(ipwr2 * x)

      return
      end
      subroutine vecset(n, v, iv, nzv, i, val)

!*****************************************************************************80

      implicit           logical (a-z)
      integer            n, iv(*), nzv, i, k
      double precision   v(*), val

!
!       set ith element of sparse vector (v, iv) with
!       nzv nonzeros to val
!

      logical set

      set = .false.
      do k = 1, nzv
         if (iv(k) .eq. i) then
            v(k) = val
            set  = .true.
         endif
      end do
      if (.not. set) then
         nzv     = nzv + 1
         v(nzv)  = val
         iv(nzv) = i
      endif
      return
      end
      subroutine print_results(name, class, n1, n2, n3, niter, &
                     t, mops, optype, verified, npbversion, &
                     compiletime, cs1, cs2, cs3, cs4, cs5, cs6, cs7)

!*****************************************************************************80
!
!! PRINT_RESULTS prints the results.
!
      implicit none

      character name*(*)
      character class*1
      integer   n1, n2, n3, niter, j
      double precision t, mops
      character optype*24, size*15
      logical   verified
      character*(*) npbversion, compiletime, &
                    cs1, cs2, cs3, cs4, cs5, cs6, cs7

         write (*, 2) name
 2       format(//, ' ', A, ' Benchmark Completed.')

         write (*, 3) Class
 3       format(' Class           = ', 12x, a12)
!
!   If this is not a grid-based problem (EP, FT, CG), then
!   we only print n1, which contains some measure of the
!   problem size. In that case, n2 and n3 are both zero.
!   Otherwise, we print the grid size n1xn2xn3
!
         if ((n2 .eq. 0) .and. (n3 .eq. 0)) then
            if (name(1:2) .eq. 'EP') then
               write(size, '(f15.0)' ) 2.d0**n1
               j = 15
               if (size(j:j) .eq. '.') then
                  size(j:j) = ' '
                  j = j - 1
               endif
               write (*,42) size(1:j)
 42            format(' Size            = ',9x, a15)
            else
               write (*,44) n1
 44            format(' Size            = ',12x, i12)
            endif
         else
            write (*, 4) n1,n2,n3
 4          format(' Size            =  ',9x, i4,'x',i4,'x',i4)
         endif

         write (*, 5) niter
 5       format(' Iterations      = ', 12x, i12)

         write (*, 6) t
 6       format(' Time in seconds = ',12x, f12.2)

         write (*,9) mops
 9       format(' Mop/s total     = ',12x, f12.2)

         write(*, 11) optype
 11      format(' Operation type  = ', a24)

         if (verified) then
            write(*,12) '  SUCCESSFUL'
         else
            write(*,12) 'UNSUCCESSFUL'
         endif
 12      format(' Verification    = ', 12x, a)

         write(*,13) npbversion
 13      format(' Version         = ', 12x, a12)

         write(*,14) compiletime
 14      format(' Compile date    = ', 12x, a12)


         write (*,121) cs1
 121     format(/, ' Compile options:', /, &
                '    F77          = ', A)

         write (*,122) cs2
 122     format('    FLINK        = ', A)

         write (*,123) cs3
 123     format('    F_LIB        = ', A)

         write (*,124) cs4
 124     format('    F_INC        = ', A)

         write (*,125) cs5
 125     format('    FFLAGS       = ', A)

         write (*,126) cs6
 126     format('    FLINKFLAGS   = ', A)

         write(*, 127) cs7
 127     format('    RAND         = ', A)

         write (*,130)
 130     format(//' Please send all errors/feedbacks to:'// &
                  ' NPB Development Team'/ &
                  ' npb@nas.nasa.gov'//)

         return
         end
      double precision function randlc(x, a)

!*****************************************************************************80
!
!! RANDLC returns a uniform pseudorandom double precision number.
!
!  Discussion:
!
!    The number returned is in the range (0, 1).
!
!    The algorithm uses the linear congruential generator
!
!   x_{k+1} = a x_k  (mod 2^46)
!
!   where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
!   before repeating.  The argument A is the same as 'a' in the above formula,
!   and X is the same as x_0.  A and X must be odd double precision integers
!   in the range (1, 2^46).  The returned value RANDLC is normalized to be
!   between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
!   the new seed x_1, so that subsequent calls to RANDLC using the same
!   arguments will generate a continuous sequence.
!
      implicit none

      double precision x, a
      integer*8 i246m1, Lx, La
      double precision d2m46

      parameter(d2m46=0.5d0**46)

      save i246m1
      data i246m1/X'00003FFFFFFFFFFF'/

      Lx = X
      La = A

      Lx   = iand(Lx*La,i246m1)
      randlc = d2m46*dble(Lx)
      x    = dble(Lx)
      return
      end
      subroutine timer_clear ( n )

!*****************************************************************************80
!
!! TIMER_CLEAR clears the timer.
!
      implicit none

      integer n

      double precision start(64), elapsed(64)
      common /tt/ start, elapsed

      elapsed(n) = 0.0
      return
      end
      subroutine timer_start ( n )

!*****************************************************************************80
!
!! TIMER_START starts the timer.
!
      implicit none

      external         elapsed_time
      double precision elapsed_time
      integer n
      double precision start(64), elapsed(64)
      common /tt/ start, elapsed

      start(n) = elapsed_time()

      return
      end
      subroutine timer_stop ( n )

!*****************************************************************************80
!
!! TIMER_STOP stops the timer.
!
      implicit none

      external         elapsed_time
      double precision elapsed_time
      integer n
      double precision start(64), elapsed(64)
      common /tt/ start, elapsed
      double precision t, now
      now = elapsed_time()
      t = now - start(n)
      elapsed(n) = elapsed(n) + t

      return
      end
      double precision function timer_read ( n )

!*****************************************************************************80
!
!! TIMER_READ reads the timer.
!
      implicit none

      integer n
      double precision start(64), elapsed(64)
      common /tt/ start, elapsed

      timer_read = elapsed(n)
      return
      end
      double precision function elapsed_time ( )

!*****************************************************************************80
!
!! ELAPSED_TIME measures wall clock time.
!
      implicit none

      double precision t
!
!  This function must measure wall clock time, not CPU time.
!  Since there is no portable timer in Fortran (77)
!  we call a routine compiled in C (though the C source may have
!  to be tweaked).
!
      call wtime(t)
!
!  The following is not ok for "official" results because it reports
!  CPU time not wall clock time. It may be useful for developing/testing
!  on timeshared Crays, though.
!     call second(t)
!
      elapsed_time = t

      return
      end


