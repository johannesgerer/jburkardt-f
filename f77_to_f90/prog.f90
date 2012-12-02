!  sblas2.f  02 June 2000                                               
!                                                                       
      subroutine smxpy ( n1, y, n2, lda, x, a ) 
!                                                                       
!***********************************************************************
!                                                                       
!! SMXPY computes Y = Y + A * X, where A is a matrix.                   
!                                                                       
!                                                                       
!  Parameters:                                                          
!                                                                       
!    Input, integer N1, the number of entries in Y.                     
!                                                                       
!    Input/output, real Y(N1), the vector to which A*X is to be added.  
!                                                                       
!    Input, integer N2, the number of entries in X.                     
!                                                                       
!    Input, integer LDA, the leading dimension of the array A.          
!    LDA must be at least N1.                                           
!                                                                       
!    Input, real X(N2), the vector to be multiplied by A.               
!                                                                       
!    Input, real A(LDA,N2), the N1 by N2 matrix which is to multiply X. 
!                                                                       
      integer lda 
      integer n2 
!                                                                       
      real a(lda,n2) 
      integer i 
      integer j 
      integer n1 
      real x(n2) 
      real y(n1) 
!                                                                       
      do j = 1, n2 
        do i = 1, n1 
          y(i) = y(i) + a(i,j) * x(j) 
        end do 
      end do 
                                                                        
      return 
      END                                           
      subroutine sxmpy ( n1, ldy, y, n2, ldx, x, lda, a ) 
!                                                                       
!***********************************************************************
!                                                                       
!! SXMPY computes y(1,j)=y(1,j) + sum(i= 1 to n2) x(1,i) * a(i,j)       
!                                                                       
!                                                                       
!  Parameters:                                                          
!                                                                       
!    Input, integer N1, the number of entries in the row vector Y.      
!    This is usually the number of columns in the two dimensional array 
!                                                                       
!    Input, integer LDY, the leading dimension of the array Y.          
!                                                                       
!    Input/output, real Y(LDY,N1), an array containing the row          
!    vector Y as its first row.  On output, X * A has been added to the 
!    first row of y.                                                    
!                                                                       
!    Input, integer N2, the number of entries in the row vector X,      
!    This is usually the number of columns in the two dimensional array 
!                                                                       
!    Input, integer LDX, the leading dimension of the array X.          
!                                                                       
!    Input, real X(LDX,N2), an array containing the row vector X as its 
!    first row.                                                         
!                                                                       
!    Input, integer LDA, the leading dimension of the array A.          
!                                                                       
!    Input, real A(LDA,N1), an array which is to be multiplied by X.    
!                                                                       
      integer lda 
      integer ldx 
      integer ldy 
      integer n1 
      integer n2 
!                                                                       
      real a(lda,n1) 
      integer i 
      integer j 
      real x(ldx,n2) 
      real y(ldy,n1) 
!                                                                       
      do j = 1, n2 
        do i = 1, n1 
          y(1,i) = y(1,i) + x(1,j) * a(j,i) 
        end do 
      end do 
                                                                        
      return 
      END                                           
      subroutine sgbmv ( trans, m, n, kl, ku, alpha, a, lda, x, incx,   &
     &                   beta, y, incy )                                
!                                                                       
!***********************************************************************
!                                                                       
!! SGBMV computes Y := ALPHA * A * X + BETA * Y.                        
!                                                                       
!                                                                       
!  Discussion:                                                          
!                                                                       
!    SGBMV can also compute Y := ALPHA * A' * X + BETA * Y.             
!                                                                       
!    ALPHA and BETA are scalars, X and Y are vectors and A is an        
!    M by N band matrix, with KL sub-diagonals and KU super-diagonals.  
!                                                                       
!  Parameters:                                                          
!                                                                       
!  trans - character.                                                   
!           on entry, trans specifies the operation to be performed as  
!           follows:                                                    
!                                                                       
!              trans = 'N' or 'N'   y := alpha * a * x + beta*y.        
!                                                                       
!              trans = 'T' or 'T'   y := alpha * a' * x + beta*y.       
!                                                                       
!              trans = 'C' or 'C'   y := alpha * a' * x + beta*y.       
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  m - integer.                                                         
!           on entry, m specifies the number of rows of the matrix a.   
!           m must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  n - integer.                                                         
!           on entry, n specifies the number of columns of the matrix a.
!           n must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  kl - integer.                                                        
!           on entry, kl specifies the number of sub-diagonals of the   
!           matrix a. kl must satisfy  0  <= kl.                        
!           unchanged on exit.                                          
!                                                                       
!  ku - integer.                                                        
!           on entry, ku specifies the number of super-diagonals of the 
!           matrix a. ku must satisfy  0  <= ku.                        
!           unchanged on exit.                                          
!                                                                       
!  alpha - real            .                                            
!           on entry, alpha specifies the scalar alpha.                 
!           unchanged on exit.                                          
!                                                                       
!  a - real             array of dimension ( lda, n ).                  
!           before entry, the leading ( kl + ku + 1 ) by n part of the  
!           array a must contain the matrix of coefficients, supplied   
!           column by column, with the leading diagonal of the matrix in
!           row ( ku + 1 ) of the array, the first super-diagonal       
!           starting at position 2 in row ku, the first sub-diagonal    
!           starting at position 1 in row ( ku + 2 ), and so on.        
!           elements in the array a that do not correspond to elements  
!           in the band matrix (such as the top left ku by ku triangle) 
!           are not referenced.                                         
!           the following program segment will transfer a band matrix   
!           from conventional full matrix storage to band storage:      
!                                                                       
!                 do j = 1, n                                           
!                    k = ku + 1 - j                                     
!                    do i = max ( 1, j - ku ), min ( m, j + kl )        
!                       a( k + i,j) = matrix(i,j)                       
!                    end do                                             
!                 end do                                                
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  lda - integer.                                                       
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least          
!           ( kl + ku + 1 ).                                            
!           unchanged on exit.                                          
!                                                                       
!  x - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incx ) ) when trans = 'N' or 'N'     
!           and at least                                                
!           ( 1 + ( m - 1 ) * abs( incx ) ) otherwise.                  
!           before entry, the incremented array x must contain the      
!           vector x.                                                   
!           unchanged on exit.                                          
!                                                                       
!  incx - integer.                                                      
!           on entry, incx specifies the increment for the elements of  
!           x. incx must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!  beta - real            .                                             
!           on entry, beta specifies the scalar beta. when beta is      
!           supplied as zero then y need not be set on input.           
!           unchanged on exit.                                          
!                                                                       
!  y - real             array of dimension at least                     
!           ( 1 + ( m - 1 ) * abs( incy ) ) when trans = 'N' or 'N'     
!           and at least                                                
!           ( 1 + ( n - 1 ) * abs( incy ) ) otherwise.                  
!           before entry, the incremented array y must contain the      
!           vector y. on exit, y is overwritten by the updated vector y.
!                                                                       
!  incy - integer.                                                      
!           on entry, incy specifies the increment for the elements of  
!           y. incy must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!                                                                       
!  level 2 blas routine.                                                
!                                                                       
! -- written on 22-october-1986.                                        
!     jack dongarra, argonne national lab.                              
!     jeremy du croz, nag central office.                               
!     sven hammarling, nag central office.                              
!     richard hanson, sandia national labs.                             
!                                                                       
      integer lda 
!                                                                       
      real alpha 
      real beta 
      integer incx 
      integer incy 
      integer kl 
      integer ku 
      integer m 
      integer n 
      character trans 
      real               a( lda, * ), x( * ), y( * ) 
      real               temp 
      integer            i, info, ix, iy, j, jx, jy, k, kup1, kx, ky 
      integer lenx, leny 
                                                                        
      logical            lsame 
      external           lsame 
                                                                        
      external           xerbla 
!                                                                       
!  test the input parameters.                                           
!                                                                       
      info = 0 
                                                                        
      if     ( .not.lsame ( trans, 'N' ).and.                           &
     &         .not.lsame ( trans, 'T' ).and.                           &
     &         .not.lsame ( trans, 'C' )     ) then                     
         info = 1 
      else if ( m < 0 ) then 
         info = 2 
      else if ( n < 0 ) then 
         info = 3 
      else if ( kl < 0 ) then 
         info = 4 
      else if ( ku < 0 ) then 
         info = 5 
      else if ( lda < ( kl + ku + 1 ) ) then 
         info = 8 
      else if ( incx == 0 ) then 
         info = 10 
      else if ( incy == 0 ) then 
         info = 13 
      end if 
                                                                        
      if ( info /= 0 ) then 
         call xerbla ( 'sgbmv ', info ) 
         return 
      end if 
!                                                                       
!  Quick return if possible.                                            
!                                                                       
      if ( m == 0 ) then 
        return 
      else if ( n == 0 ) then 
        return 
      else if ( alpha == 0.0 .and. beta == 1.0 ) then 
        return 
      end if 
!                                                                       
!  set  lenx  and  leny, the lengths of the vectors x and y, and set    
!  up the start points in  x  and  y.                                   
!                                                                       
      if ( lsame ( trans, 'N' ) ) then 
         lenx = n 
         leny = m 
      else 
         lenx = m 
         leny = n 
      end if 
                                                                        
      if ( incx > 0 ) then 
         kx = 1 
      else 
         kx = 1 - ( lenx - 1 ) * incx 
      end if 
                                                                        
      if ( incy > 0 ) then 
         ky = 1 
      else 
         ky = 1 - ( leny - 1 ) * incy 
      end if 
!                                                                       
!  Start the operations. in this version the elements of A are          
!  accessed sequentially with one pass through the band part of A.      
!                                                                       
!  First form  y := beta*y.                                             
!                                                                       
      if ( beta /= 1.0 ) then 
         if ( incy == 1 ) then 
            if ( beta == 0.0 ) then 
               do i = 1, leny 
                  y(i) = 0.0 
               end do 
            else 
               do i = 1, leny 
                  y(i) = beta * y(i) 
               end do 
            end if 
         else 
            iy = ky 
            if ( beta == 0.0 ) then 
               do i = 1, leny 
                  y(iy) = 0.0 
                  iy = iy + incy 
               end do 
            else 
               do i = 1, leny 
                  y(iy) = beta * y(iy) 
                  iy = iy + incy 
               end do 
            end if 
         end if 
      end if 
                                                                        
      if ( alpha == 0.0 ) then 
        return 
      end if 
                                                                        
      kup1 = ku + 1 
      if ( lsame ( trans, 'N' ) ) then 
!                                                                       
!  Form  y := alpha * a * x + y.                                        
!                                                                       
         jx = kx 
         if ( incy == 1 ) then 
            do j = 1, n 
               if ( x(jx) /= 0.0 ) then 
                  temp = alpha * x(jx) 
                  k = kup1 - j 
                  do i = max ( 1, j - ku ), min ( m, j + kl ) 
                     y(i) = y(i) + temp * a( k + i,j) 
                  end do 
               end if 
               jx = jx + incx 
            end do 
         else 
            do j = 1, n 
               if ( x(jx) /= 0.0 ) then 
                  temp = alpha * x(jx) 
                  iy = ky 
                  k = kup1 - j 
                  do i = max ( 1, j - ku ), min ( m, j + kl ) 
                     y(iy) = y(iy) + temp * a( k + i,j) 
                     iy = iy + incy 
                  end do 
               end if 
               jx = jx + incx 
               if ( j>ku ) then 
                  ky = ky + incy 
               end if 
            end do 
         end if 
      else 
!                                                                       
!  Form  y := alpha * a' * x + y.                                       
!                                                                       
         jy = ky 
         if ( incx == 1 ) then 
            do j = 1, n 
               temp = 0.0 
               k = kup1 - j 
               do i = max ( 1, j - ku ), min ( m, j + kl ) 
                  temp = temp + a( k + i,j) * x(i) 
               end do 
               y(jy) = y(jy) + alpha * temp 
               jy = jy + incy 
            end do 
         else 
            do j = 1, n 
               temp = 0.0 
               ix = kx 
               k = kup1 - j 
               do i = max ( 1, j - ku ), min ( m, j + kl ) 
                  temp = temp + a( k + i,j) * x(ix) 
                  ix = ix + incx 
               end do 
               y(jy) = y(jy) + alpha * temp 
               jy = jy + incy 
               if ( j>ku ) then 
                 kx = kx + incx 
               end if 
            end do 
         end if 
      end if 
                                                                        
      return 
      END                                           
      subroutine sgemv ( trans, m, n, alpha, a, lda, x, incx,           &
     &                   beta, y, incy )                                
!                                                                       
!***********************************************************************
!                                                                       
!! SGEMV performs one of the matrix-vector operations                   
!                                                                       
!     y := alpha * a * x + beta*y,   or   y := alpha * a'*x + beta*y,   
!                                                                       
!  where alpha and beta are scalars, x and y are vectors and a is an    
!  m by n matrix.                                                       
!                                                                       
!  Parameters:                                                          
!                                                                       
!  trans - character.                                                   
!           on entry, trans specifies the operation to be performed as  
!           follows:                                                    
!                                                                       
!              trans = 'N' or 'N'   y := alpha * a * x + beta*y.        
!                                                                       
!              trans = 'T' or 'T'   y := alpha * a' * x + beta*y.       
!                                                                       
!              trans = 'C' or 'C'   y := alpha * a' * x + beta*y.       
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  m - integer.                                                         
!           on entry, m specifies the number of rows of the matrix a.   
!           m must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  n - integer.                                                         
!           on entry, n specifies the number of columns of the matrix a.
!           n must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  alpha - real            .                                            
!           on entry, alpha specifies the scalar alpha.                 
!           unchanged on exit.                                          
!                                                                       
!  a - real             array of dimension ( lda, n ).                  
!           before entry, the leading m by n part of the array a must   
!           contain the matrix of coefficients.                         
!           unchanged on exit.                                          
!                                                                       
!  lda - integer.                                                       
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least          
!           max ( 1, m ).                                               
!           unchanged on exit.                                          
!                                                                       
!  x - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incx ) ) when trans = 'N' or 'N'     
!           and at least                                                
!           ( 1 + ( m - 1 ) * abs( incx ) ) otherwise.                  
!           before entry, the incremented array x must contain the      
!           vector x.                                                   
!           unchanged on exit.                                          
!                                                                       
!  incx - integer.                                                      
!           on entry, incx specifies the increment for the elements of  
!           x. incx must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!  beta - real            .                                             
!           on entry, beta specifies the scalar beta. when beta is      
!           supplied as zero then y need not be set on input.           
!           unchanged on exit.                                          
!                                                                       
!  y - real             array of dimension at least                     
!           ( 1 + ( m - 1 ) * abs( incy ) ) when trans = 'N' or 'N'     
!           and at least                                                
!           ( 1 + ( n - 1 ) * abs( incy ) ) otherwise.                  
!           before entry with beta non-zero, the incremented array y    
!           must contain the vector y. on exit, y is overwritten by the 
!           updated vector y.                                           
!                                                                       
!  incy - integer.                                                      
!           on entry, incy specifies the increment for the elements of  
!           y. incy must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!                                                                       
!  level 2 blas routine.                                                
!                                                                       
! -- written on 22-october-1986.                                        
!     jack dongarra, argonne national lab.                              
!     jeremy du croz, nag central office.                               
!     sven hammarling, nag central office.                              
!     richard hanson, sandia national labs.                             
!                                                                       
!                                                                       
      real               alpha, beta 
      integer            incx, incy, lda, m, n 
      character        trans 
!     .. array arguments ..                                             
      real               a( lda, * ), x( * ), y( * ) 
      real               temp 
      integer            i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny 
!     .. external functions ..                                          
      logical            lsame 
      external           lsame 
!     .. external subroutines ..                                        
      external           xerbla 
!     .. intrinsic functions ..                                         
      intrinsic          max 
!     ..                                                                
!     .. executable statements ..                                       
!                                                                       
!     test the input parameters.                                        
!                                                                       
      info = 0 
      if     ( .not.lsame ( trans, 'N' ).and.                           &
     &         .not.lsame ( trans, 'T' ).and.                           &
     &         .not.lsame ( trans, 'C' )     ) then                     
         info = 1 
      else if ( m < 0 ) then 
         info = 2 
      else if ( n < 0 ) then 
         info = 3 
      else if ( lda < max ( 1, m ) ) then 
         info = 6 
      else if ( incx == 0 ) then 
         info = 8 
      else if ( incy == 0 ) then 
         info = 11 
      end if 
                                                                        
      if ( info /= 0 ) then 
         call xerbla ( 'sgemv ', info ) 
         return 
      end if 
!                                                                       
!  Quick return if possible.                                            
!                                                                       
      if ( ( m == 0 ).or.( n== 0 ).or.                                  &
     &    ( ( alpha == 0.0 ).and.( beta== 1.0 ) ) )                     &
     &   return                                                         
!                                                                       
!  set  lenx  and  leny, the lengths of the vectors x and y, and set    
!  up the start points in  x  and  y.                                   
!                                                                       
      if ( lsame ( trans, 'N' ) ) then 
         lenx = n 
         leny = m 
      else 
         lenx = m 
         leny = n 
      end if 
                                                                        
      if ( incx > 0 ) then 
         kx = 1 
      else 
         kx = 1 - ( lenx - 1 ) * incx 
      end if 
                                                                        
      if ( incy > 0 ) then 
         ky = 1 
      else 
         ky = 1 - ( leny - 1 ) * incy 
      end if 
!                                                                       
!  Start the operations. in this version the elements of a are          
!  accessed sequentially with one pass through a.                       
!                                                                       
!  First form  y := beta*y.                                             
!                                                                       
      if ( beta /= 1.0 ) then 
         if ( incy == 1 ) then 
            if ( beta == 0.0 ) then 
               do i = 1, leny 
                  y(i) = 0.0 
               end do 
            else 
               do i = 1, leny 
                  y(i) = beta * y(i) 
               end do 
            end if 
         else 
            iy = ky 
            if ( beta == 0.0 ) then 
               do i = 1, leny 
                  y(iy) = 0.0 
                  iy = iy + incy 
               end do 
            else 
               do i = 1, leny 
                  y(iy) = beta * y(iy) 
                  iy = iy + incy 
               end do 
            end if 
         end if 
      end if 
                                                                        
      if ( alpha == 0.0 ) then 
        return 
      end if 
                                                                        
      if ( lsame ( trans, 'N' ) ) then 
!                                                                       
!  Form  y := alpha * a * x + y.                                        
!                                                                       
         jx = kx 
         if ( incy == 1 ) then 
            do j = 1, n 
               if ( x(jx) /= 0.0 ) then 
                  temp = alpha * x(jx) 
                  do i = 1, m 
                     y(i) = y(i) + temp * a(i,j) 
                  end do 
               end if 
               jx = jx + incx 
            end do 
         else 
            do j = 1, n 
               if ( x(jx) /= 0.0 ) then 
                  temp = alpha * x(jx) 
                  iy = ky 
                  do i = 1, m 
                     y(iy) = y(iy) + temp * a(i,j) 
                     iy = iy + incy 
                  end do 
               end if 
               jx = jx + incx 
            end do 
         end if 
      else 
!                                                                       
!  Form  y := alpha * a' * x + y.                                       
!                                                                       
         jy = ky 
         if ( incx == 1 ) then 
            do j = 1, n 
               temp = 0.0 
               do i = 1, m 
                  temp = temp + a(i,j) * x(i) 
               end do 
               y(jy) = y(jy) + alpha * temp 
               jy = jy + incy 
            end do 
         else 
            do j = 1, n 
               temp = 0.0 
               ix = kx 
               do i = 1, m 
                  temp = temp + a(i,j) * x(ix) 
                  ix = ix + incx 
               end do 
               y(jy) = y(jy) + alpha * temp 
               jy = jy + incy 
            end do 
         end if 
      end if 
                                                                        
      return 
      END                                           
      subroutine sger  ( m, n, alpha, x, incx, y, incy, a, lda ) 
!                                                                       
!***********************************************************************
!                                                                       
!! SGER performs the rank 1 operation                                   
!                                                                       
!     a := alpha * x*y' + a,                                            
!                                                                       
!  where alpha is a scalar, x is an m element vector, y is an n element 
!  vector and a is an m by n matrix.                                    
!                                                                       
!  Parameters:                                                          
!                                                                       
!  m - integer.                                                         
!           on entry, m specifies the number of rows of the matrix a.   
!           m must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  n - integer.                                                         
!           on entry, n specifies the number of columns of the matrix a.
!           n must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  alpha - real            .                                            
!           on entry, alpha specifies the scalar alpha.                 
!           unchanged on exit.                                          
!                                                                       
!  x - real             array of dimension at least                     
!           ( 1 + ( m - 1 ) * abs( incx ) ).                            
!           before entry, the incremented array x must contain the m    
!           element vector x.                                           
!           unchanged on exit.                                          
!                                                                       
!  incx - integer.                                                      
!           on entry, incx specifies the increment for the elements of  
!           x. incx must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!  y - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incy ) ).                            
!           before entry, the incremented array y must contain the n    
!           element vector y.                                           
!           unchanged on exit.                                          
!                                                                       
!  incy - integer.                                                      
!           on entry, incy specifies the increment for the elements of  
!           y. incy must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!  a - real             array of dimension ( lda, n ).                  
!           before entry, the leading m by n part of the array a must   
!           contain the matrix of coefficients. on exit, a is           
!           overwritten by the updated matrix.                          
!                                                                       
!  lda - integer.                                                       
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least          
!           max ( 1, m ).                                               
!           unchanged on exit.                                          
!                                                                       
!                                                                       
!  level 2 blas routine.                                                
!                                                                       
! -- written on 22-october-1986.                                        
!     jack dongarra, argonne national lab.                              
!     jeremy du croz, nag central office.                               
!     sven hammarling, nag central office.                              
!     richard hanson, sandia national labs.                             
!                                                                       
      real               alpha 
      integer            incx, incy, lda, m, n 
!     .. array arguments ..                                             
      real               a( lda, * ), x( * ), y( * ) 
!     ..                                                                
!                                                                       
!     .. local scalars ..                                               
      real               temp 
      integer            i, info, ix, j, jy, kx 
!     .. external subroutines ..                                        
      external           xerbla 
!     .. intrinsic functions ..                                         
      intrinsic          max 
!     ..                                                                
!     .. executable statements ..                                       
!                                                                       
!     test the input parameters.                                        
!                                                                       
      info = 0 
      if     ( m < 0 ) then 
         info = 1 
      else if ( n < 0 ) then 
         info = 2 
      else if ( incx == 0 ) then 
         info = 5 
      else if ( incy == 0 ) then 
         info = 7 
      else if ( lda < max ( 1, m ) ) then 
         info = 9 
      end if 
      if ( info /= 0 ) then 
         call xerbla ( 'sger  ', info ) 
         return 
      end if 
!                                                                       
!  Quick return if possible.                                            
!                                                                       
      if ( ( m == 0 ).or.( n== 0 ).or.( alpha==0.0 ) ) then 
        return 
      end if 
!                                                                       
!  Start the operations. in this version the elements of a are          
!  accessed sequentially with one pass through a.                       
!                                                                       
      if ( incy > 0 ) then 
         jy = 1 
      else 
         jy = 1 - ( n - 1 ) * incy 
      end if 
                                                                        
      if ( incx == 1 ) then 
                                                                        
         do j = 1, n 
            if ( y(jy) /= 0.0 ) then 
               temp = alpha * y(jy) 
               do i = 1, m 
                  a(i,j) = a(i,j) + x(i) * temp 
              end do 
            end if 
            jy = jy + incy 
         end do 
                                                                        
      else 
                                                                        
         if ( incx > 0 ) then 
            kx = 1 
         else 
            kx = 1 - ( m - 1 ) * incx 
         end if 
                                                                        
         do j = 1, n 
            if ( y(jy) /= 0.0 ) then 
               temp = alpha * y(jy) 
               ix = kx 
               do i = 1, m 
                  a(i,j) = a(i,j) + x(ix) * temp 
                  ix = ix + incx 
               end do 
            end if 
            jy = jy + incy 
         end do 
                                                                        
      end if 
                                                                        
      return 
      END                                           
      subroutine ssbmv ( uplo, n, k, alpha, a, lda, x, incx,            &
     &                   beta, y, incy )                                
!                                                                       
!***********************************************************************
!                                                                       
!! SSBMV performs the matrix-vector  operation                          
!                                                                       
!     y := alpha * a * x + beta*y,                                      
!                                                                       
!  where alpha and beta are scalars, x and y are n element vectors and  
!  a is an n by n symmetric band matrix, with k super-diagonals.        
!                                                                       
!  Parameters:                                                          
!                                                                       
!  uplo - character.                                                    
!           on entry, uplo specifies whether the upper or lower         
!           triangular part of the band matrix a is being supplied as   
!           follows:                                                    
!                                                                       
!              uplo = 'U' or 'U'   the upper triangular part of a is    
!                                  being supplied.                      
!                                                                       
!              uplo = 'L' or 'L'   the lower triangular part of a is    
!                                  being supplied.                      
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  n - integer.                                                         
!           on entry, n specifies the order of the matrix a.            
!           n must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  k - integer.                                                         
!           on entry, k specifies the number of super-diagonals of the  
!           matrix a. k must satisfy  0  <= k.                          
!           unchanged on exit.                                          
!                                                                       
!  alpha - real            .                                            
!           on entry, alpha specifies the scalar alpha.                 
!           unchanged on exit.                                          
!                                                                       
!  a - real             array of dimension ( lda, n ).                  
!           before entry with uplo = 'U' or 'U', the leading ( k + 1 )  
!           by n part of the array a must contain the upper triangular  
!           band part of the symmetric matrix, supplied column by       
!           column, with the leading diagonal of the matrix in row      
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. the top left k by k triangle
!           of the array a is not referenced.                           
!           the following program segment will transfer the upper       
!           triangular part of a symmetric band matrix from conventional
!           full matrix storage to band storage:                        
!                                                                       
!                 do j = 1, n                                           
!                    m = k + 1 - j                                      
!                    do i = max ( 1, j - k ), j                         
!                       a( m + i,j) = matrix(i,j)                       
!                    end do                                             
!                 end do                                                
!                                                                       
!           before entry with uplo = 'L' or 'L', the leading ( k + 1 )  
!           by n part of the array a must contain the lower triangular  
!           band part of the symmetric matrix, supplied column by       
!           column, with the leading diagonal of the matrix in row 1 of 
!           the array, the first sub-diagonal starting at position 1 in 
!           row 2, and so on. the bottom right k by k triangle of the   
!           array a is not referenced.                                  
!           the following program segment will transfer the lower       
!           triangular part of a symmetric band matrix from conventional
!           full matrix storage to band storage:                        
!                                                                       
!                 do j = 1, n                                           
!                    m = 1 - j                                          
!                    do i = j, min ( n, j + k )                         
!                       a( m + i,j) = matrix(i,j)                       
!                    end do                                             
!                 end do                                                
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  lda - integer.                                                       
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least          
!           ( k + 1 ).                                                  
!           unchanged on exit.                                          
!                                                                       
!  x - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incx ) ).                            
!           before entry, the incremented array x must contain the      
!           vector x.                                                   
!           unchanged on exit.                                          
!                                                                       
!  incx - integer.                                                      
!           on entry, incx specifies the increment for the elements of  
!           x. incx must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!  beta - real            .                                             
!           on entry, beta specifies the scalar beta.                   
!           unchanged on exit.                                          
!                                                                       
!  y - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incy ) ).                            
!           before entry, the incremented array y must contain the      
!           vector y. on exit, y is overwritten by the updated vector y.
!                                                                       
!  incy - integer.                                                      
!           on entry, incy specifies the increment for the elements of  
!           y. incy must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!                                                                       
!  level 2 blas routine.                                                
!                                                                       
! -- written on 22-october-1986.                                        
!     jack dongarra, argonne national lab.                              
!     jeremy du croz, nag central office.                               
!     sven hammarling, nag central office.                              
!     richard hanson, sandia national labs.                             
!                                                                       
      real               alpha, beta 
      integer            incx, incy, k, lda, n 
      character        uplo 
!     .. array arguments ..                                             
      real               a( lda, * ), x( * ), y( * ) 
!     ..                                                                
!                                                                       
!     .. local scalars ..                                               
      real               temp1, temp2 
      integer            i, info, ix, iy, j, jx, jy, kplus1, kx, ky, l 
!     .. external functions ..                                          
      logical            lsame 
      external           lsame 
!     .. external subroutines ..                                        
      external           xerbla 
!     .. intrinsic functions ..                                         
      intrinsic          max, min 
!     ..                                                                
!     .. executable statements ..                                       
!                                                                       
!  test the input parameters.                                           
!                                                                       
      info = 0 
      if     ( .not.lsame ( uplo, 'U' ).and.                            &
     &         .not.lsame ( uplo, 'L' )     ) then                      
         info = 1 
      else if ( n < 0 ) then 
         info = 2 
      else if ( k < 0 ) then 
         info = 3 
      else if ( lda < ( k + 1 ) ) then 
         info = 6 
      else if ( incx == 0 ) then 
         info = 8 
      else if ( incy == 0 ) then 
         info = 11 
      end if 
      if ( info /= 0 ) then 
         call xerbla ( 'ssbmv ', info ) 
         return 
      end if 
!                                                                       
!  Quick return if possible.                                            
!                                                                       
      if ( ( n == 0 ).or.( ( alpha== 0.0 ).and.( beta== 1.0 ) ) ) then 
        return 
      end if 
!                                                                       
!  set up the start points in  x  and  y.                               
!                                                                       
      if ( incx>0 ) then 
         kx = 1 
      else 
         kx = 1 - ( n - 1 ) * incx 
      end if 
                                                                        
      if ( incy>0 ) then 
         ky = 1 
      else 
         ky = 1 - ( n - 1 ) * incy 
      end if 
!                                                                       
!  start the operations. in this version the elements of the array a    
!  are accessed sequentially with one pass through a.                   
!                                                                       
!     first form  y := beta*y.                                          
!                                                                       
      if ( beta /= 1.0 ) then 
         if ( incy == 1 ) then 
            if ( beta == 0.0 ) then 
               do i = 1, n 
                  y(i) = 0.0 
               end do 
            else 
               do i = 1, n 
                  y(i) = beta*y(i) 
               end do 
            end if 
         else 
            iy = ky 
            if ( beta == 0.0 ) then 
               do i = 1, n 
                  y(iy) = 0.0 
                  iy = iy + incy 
               end do 
            else 
               do i = 1, n 
                  y(iy) = beta*y(iy) 
                  iy = iy + incy 
               end do 
            end if 
         end if 
      end if 
                                                                        
      if ( alpha == 0.0 ) then 
        return 
      end if 
                                                                        
      if ( lsame ( uplo, 'U' ) ) then 
!                                                                       
!  form  y  when upper triangle of a is stored.                         
!                                                                       
         kplus1 = k + 1 
         if ( ( incx == 1 ).and.( incy== 1 ) ) then 
            do j = 1, n 
               temp1 = alpha * x(j) 
               temp2 = 0.0 
               l = kplus1 - j 
               do i = max ( 1, j - k ), j - 1 
                  y(i) = y(i) + temp1 * a( l + i,j) 
                  temp2 = temp2 + a( l + i,j) * x(i) 
               end do 
               y(j) = y(j) + temp1 * a( kplus1,j) + alpha * temp2 
            end do 
         else 
            jx = kx 
            jy = ky 
            do j = 1, n 
               temp1 = alpha * x(jx) 
               temp2 = 0.0 
               ix = kx 
               iy = ky 
               l = kplus1 - j 
               do i = max ( 1, j - k ), j - 1 
                  y(iy) = y(iy) + temp1 * a( l + i,j) 
                  temp2 = temp2 + a( l + i,j) * x(ix) 
                  ix = ix + incx 
                  iy = iy + incy 
               end do 
               y(jy) = y(jy) + temp1 * a( kplus1,j) + alpha * temp2 
               jx = jx + incx 
               jy = jy + incy 
               if ( j>k ) then 
                  kx = kx + incx 
                  ky = ky + incy 
               end if 
            end do 
         end if 
      else 
!                                                                       
!  form  y  when lower triangle of a is stored.                         
!                                                                       
         if ( ( incx == 1 ).and.( incy== 1 ) ) then 
            do j = 1, n 
               temp1 = alpha * x(j) 
               temp2 = 0.0 
               y(j) = y(j) + temp1 * a( 1,j) 
               l = 1 - j 
               do i = j + 1, min ( n, j + k ) 
                  y(i) = y(i) + temp1 * a( l + i,j) 
                  temp2 = temp2 + a( l + i,j) * x(i) 
               end do 
               y(j) = y(j) + alpha * temp2 
            end do 
         else 
            jx = kx 
            jy = ky 
            do j = 1, n 
               temp1 = alpha * x(jx) 
               temp2 = 0.0 
               y(jy) = y(jy) + temp1 * a( 1,j) 
               l = 1 - j 
               ix = jx 
               iy = jy 
               do i = j + 1, min ( n, j + k ) 
                  ix = ix + incx 
                  iy = iy + incy 
                  y(iy) = y(iy) + temp1 * a( l + i,j) 
                  temp2 = temp2 + a( l + i,j) * x(ix) 
               end do 
               y(jy) = y(jy) + alpha * temp2 
               jx = jx + incx 
               jy = jy + incy 
            end do 
         end if 
      end if 
                                                                        
      return 
      END                                           
      subroutine sspmv ( uplo, n, alpha, ap, x, incx, beta, y, incy ) 
!                                                                       
!***********************************************************************
!                                                                       
!! SSPMV performs the matrix-vector operation                           
!                                                                       
!     y := alpha * a * x + beta*y,                                      
!                                                                       
!  where alpha and beta are scalars, x and y are n element vectors and  
!  a is an n by n symmetric matrix, supplied in packed form.            
!                                                                       
!  Parameters:                                                          
!                                                                       
!  uplo - character.                                                    
!           on entry, uplo specifies whether the upper or lower         
!           triangular part of the matrix a is supplied in the packed   
!           array ap as follows:                                        
!                                                                       
!              uplo = 'U' or 'U'   the upper triangular part of a is    
!                                  supplied in ap.                      
!                                                                       
!              uplo = 'L' or 'L'   the lower triangular part of a is    
!                                  supplied in ap.                      
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  n - integer.                                                         
!           on entry, n specifies the order of the matrix a.            
!           n must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  alpha - real            .                                            
!           on entry, alpha specifies the scalar alpha.                 
!           unchanged on exit.                                          
!                                                                       
!  ap - real             array of dimension at least                    
!           ( ( n*( n + 1 ) )/2 ).                                      
!           before entry with uplo = 'U' or 'U', the array ap must      
!           contain the upper triangular part of the symmetric matrix   
!           packed sequentially, column by column, so that ap( 1 )      
!           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )   
!           and a( 2, 2 ) respectively, and so on.                      
!           before entry with uplo = 'L' or 'L', the array ap must      
!           contain the lower triangular part of the symmetric matrix   
!           packed sequentially, column by column, so that ap( 1 )      
!           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )   
!           and a( 3, 1 ) respectively, and so on.                      
!           unchanged on exit.                                          
!                                                                       
!  x - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incx ) ).                            
!           before entry, the incremented array x must contain the n    
!           element vector x.                                           
!           unchanged on exit.                                          
!                                                                       
!  incx - integer.                                                      
!           on entry, incx specifies the increment for the elements of  
!           x. incx must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!  beta - real            .                                             
!           on entry, beta specifies the scalar beta. when beta is      
!           supplied as zero then y need not be set on input.           
!           unchanged on exit.                                          
!                                                                       
!  y - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incy ) ).                            
!           before entry, the incremented array y must contain the n    
!           element vector y. on exit, y is overwritten by the updated  
!           vector y.                                                   
!                                                                       
!  incy - integer.                                                      
!           on entry, incy specifies the increment for the elements of  
!           y. incy must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!                                                                       
!  level 2 blas routine.                                                
!                                                                       
! -- written on 22-october-1986.                                        
!     jack dongarra, argonne national lab.                              
!     jeremy du croz, nag central office.                               
!     sven hammarling, nag central office.                              
!     richard hanson, sandia national labs.                             
!                                                                       
      real               alpha, beta 
      integer            incx, incy, n 
      character        uplo 
!     .. array arguments ..                                             
      real               ap( * ), x( * ), y( * ) 
!     ..                                                                
!                                                                       
!     .. local scalars ..                                               
      real               temp1, temp2 
      integer            i, info, ix, iy, j, jx, jy, k, kk, kx, ky 
!     .. external functions ..                                          
      logical            lsame 
      external           lsame 
!     .. external subroutines ..                                        
      external           xerbla 
!     ..                                                                
!     .. executable statements ..                                       
!                                                                       
!     test the input parameters.                                        
!                                                                       
      info = 0 
      if ( .not.lsame ( uplo, 'U' ).and..not.lsame( uplo, 'L' )) then 
         info = 1 
      else if ( n < 0 ) then 
         info = 2 
      else if ( incx == 0 ) then 
         info = 6 
      else if ( incy == 0 ) then 
         info = 9 
      end if 
      if ( info /= 0 ) then 
         call xerbla ( 'sspmv ', info ) 
         return 
      end if 
!                                                                       
!  Quick return if possible.                                            
!                                                                       
      if ( ( n == 0 ).or.( ( alpha== 0.0 ).and.( beta== 1.0 ) ) ) then 
        return 
      end if 
!                                                                       
!  Set up the start points in  x  and  y.                               
!                                                                       
      if ( incx > 0 ) then 
         kx = 1 
      else 
         kx = 1 - ( n - 1 ) * incx 
      end if 
                                                                        
      if ( incy > 0 ) then 
         ky = 1 
      else 
         ky = 1 - ( n - 1 ) * incy 
      end if 
!                                                                       
!  Start the operations. in this version the elements of the array ap   
!  are accessed sequentially with one pass through ap.                  
!                                                                       
!  First form  y := beta*y.                                             
!                                                                       
      if ( beta /= 1.0 ) then 
         if ( incy == 1 ) then 
            if ( beta == 0.0 ) then 
               do i = 1, n 
                  y(i) = 0.0 
               end do 
            else 
               do i = 1, n 
                  y(i) = beta * y(i) 
               end do 
            end if 
         else 
            iy = ky 
            if ( beta == 0.0 ) then 
               do i = 1, n 
                  y(iy) = 0.0 
                  iy = iy + incy 
               end do 
            else 
               do i = 1, n 
                  y(iy) = beta * y(iy) 
                  iy = iy + incy 
               end do 
            end if 
         end if 
      end if 
                                                                        
      if ( alpha == 0.0 ) then 
        return 
      end if 
                                                                        
      kk = 1 
      if ( lsame ( uplo, 'U' ) ) then 
!                                                                       
!  form  y  when ap contains the upper triangle.                        
!                                                                       
         if ( ( incx == 1 ).and.( incy== 1 ) ) then 
            do j = 1, n 
               temp1 = alpha * x(j) 
               temp2 = 0.0 
               k = kk 
               do i = 1, j - 1 
                  y(i) = y(i) + temp1 * ap(k) 
                  temp2 = temp2 + ap(k) * x(i) 
                  k = k + 1 
               end do 
               y(j) = y(j) + temp1 * ap( kk + j - 1 ) + alpha * temp2 
               kk = kk + j 
            end do 
         else 
            jx = kx 
            jy = ky 
            do j = 1, n 
               temp1 = alpha * x(jx) 
               temp2 = 0.0 
               ix = kx 
               iy = ky 
               do k = kk, kk + j - 2 
                  y(iy) = y(iy) + temp1 * ap(k) 
                  temp2 = temp2 + ap(k) * x(ix) 
                  ix = ix + incx 
                  iy = iy + incy 
               end do 
               y(jy) = y(jy) + temp1 * ap( kk + j - 1 ) + alpha * temp2 
               jx = jx + incx 
               jy = jy + incy 
               kk = kk + j 
            end do 
         end if 
      else 
!                                                                       
!  form  y  when ap contains the lower triangle.                        
!                                                                       
         if ( ( incx == 1 ).and.( incy== 1 ) ) then 
            do j = 1, n 
               temp1 = alpha * x(j) 
               temp2 = 0.0 
               y(j) = y(j) + temp1 * ap(kk) 
               k = kk + 1 
               do i = j + 1, n 
                  y(i) = y(i) + temp1 * ap(k) 
                  temp2 = temp2 + ap(k) * x(i) 
                  k = k + 1 
               end do 
               y(j) = y(j) + alpha * temp2 
               kk = kk + ( n - j + 1 ) 
            end do 
         else 
            jx = kx 
            jy = ky 
            do j = 1, n 
               temp1 = alpha * x(jx) 
               temp2 = 0.0 
               y(jy) = y(jy) + temp1 * ap(kk) 
               ix = jx 
               iy = jy 
               do k = kk + 1, kk + n - j 
                  ix = ix + incx 
                  iy = iy + incy 
                  y(iy) = y(iy) + temp1 * ap(k) 
                  temp2 = temp2 + ap(k) * x(ix) 
               end do 
               y(jy) = y(jy) + alpha * temp2 
               jx = jx + incx 
               jy = jy + incy 
               kk = kk + ( n - j + 1 ) 
            end do 
         end if 
      end if 
                                                                        
      return 
      END                                           
      subroutine sspr  ( uplo, n, alpha, x, incx, ap ) 
!                                                                       
!***********************************************************************
!                                                                       
!! SSPR performs the symmetric rank 1 operation                         
!                                                                       
!     a := alpha * x*x' + a,                                            
!                                                                       
!  where alpha is a real scalar, x is an n element vector and a is an   
!  n by n symmetric matrix, supplied in packed form.                    
!                                                                       
!  Parameters:                                                          
!                                                                       
!  uplo - character.                                                    
!           on entry, uplo specifies whether the upper or lower         
!           triangular part of the matrix a is supplied in the packed   
!           array ap as follows:                                        
!                                                                       
!              uplo = 'U' or 'U'   the upper triangular part of a is    
!                                  supplied in ap.                      
!                                                                       
!              uplo = 'L' or 'L'   the lower triangular part of a is    
!                                  supplied in ap.                      
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  n - integer.                                                         
!           on entry, n specifies the order of the matrix a.            
!           n must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  alpha - real            .                                            
!           on entry, alpha specifies the scalar alpha.                 
!           unchanged on exit.                                          
!                                                                       
!  x - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incx ) ).                            
!           before entry, the incremented array x must contain the n    
!           element vector x.                                           
!           unchanged on exit.                                          
!                                                                       
!  incx - integer.                                                      
!           on entry, incx specifies the increment for the elements of  
!           x. incx must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!  ap - real             array of dimension at least                    
!           ( ( n*( n + 1 ) )/2 ).                                      
!           before entry with  uplo = 'U' or 'U', the array ap must     
!           contain the upper triangular part of the symmetric matrix   
!           packed sequentially, column by column, so that ap( 1 )      
!           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )   
!           and a( 2, 2 ) respectively, and so on. on exit, the array   
!           ap is overwritten by the upper triangular part of the       
!           updated matrix.                                             
!           before entry with uplo = 'L' or 'L', the array ap must      
!           contain the lower triangular part of the symmetric matrix   
!           packed sequentially, column by column, so that ap( 1 )      
!           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )   
!           and a( 3, 1 ) respectively, and so on. on exit, the array   
!           ap is overwritten by the lower triangular part of the       
!           updated matrix.                                             
!                                                                       
!                                                                       
!  level 2 blas routine.                                                
!                                                                       
! -- written on 22-october-1986.                                        
!     jack dongarra, argonne national lab.                              
!     jeremy du croz, nag central office.                               
!     sven hammarling, nag central office.                              
!     richard hanson, sandia national labs.                             
!                                                                       
      real               alpha 
      integer            incx, n 
      character        uplo 
!     .. array arguments ..                                             
      real               ap( * ), x( * ) 
!     ..                                                                
!                                                                       
!     .. local scalars ..                                               
      real               temp 
      integer            i, info, ix, j, jx, k, kk, kx 
!     .. external functions ..                                          
      logical            lsame 
      external           lsame 
!     .. external subroutines ..                                        
      external           xerbla 
!     ..                                                                
!     .. executable statements ..                                       
!                                                                       
!     test the input parameters.                                        
!                                                                       
      info = 0 
      if ( .not.lsame ( uplo, 'U' ).and..not.lsame( uplo, 'L' ) ) then 
         info = 1 
      else if ( n < 0 ) then 
         info = 2 
      else if ( incx == 0 ) then 
         info = 5 
      end if 
      if ( info /= 0 ) then 
         call xerbla ( 'sspr  ', info ) 
         return 
      end if 
!                                                                       
!  Quick return if possible.                                            
!                                                                       
      if ( ( n == 0 ).or.( alpha== 0.0 ) ) then 
        return 
      end if 
!                                                                       
!  set the start point in x if the increment is not unity.              
!                                                                       
      if ( incx <= 0 ) then 
         kx = 1 - ( n - 1 ) * incx 
      else if ( incx /= 1 ) then 
         kx = 1 
      end if 
!                                                                       
!  start the operations. in this version the elements of the array ap   
!  are accessed sequentially with one pass through ap.                  
!                                                                       
      kk = 1 
      if ( lsame ( uplo, 'U' ) ) then 
!                                                                       
!  form  a  when upper triangle is stored in ap.                        
!                                                                       
         if ( incx == 1 ) then 
            do j = 1, n 
               if ( x(j) /= 0.0 ) then 
                  temp = alpha * x(j) 
                  k = kk 
                  do i = 1, j 
                     ap(k) = ap( k ) + x(i) * temp 
                     k = k + 1 
                  end do 
               end if 
               kk = kk + j 
            end do 
         else 
            jx = kx 
            do j = 1, n 
               if ( x(jx) /= 0.0 ) then 
                  temp = alpha * x(jx) 
                  ix = kx 
                  do k = kk, kk + j - 1 
                     ap(k) = ap( k ) + x(ix) * temp 
                     ix = ix + incx 
                  end do 
               end if 
               jx = jx + incx 
               kk = kk + j 
            end do 
         end if 
      else 
!                                                                       
!  Form  a  when lower triangle is stored in ap.                        
!                                                                       
         if ( incx == 1 ) then 
            do j = 1, n 
               if ( x(j) /= 0.0 ) then 
                  temp = alpha * x(j) 
                  k = kk 
                  do i = j, n 
                     ap(k) = ap( k ) + x(i) * temp 
                     k = k + 1 
                  end do 
               end if 
               kk = kk + n - j + 1 
            end do 
         else 
            jx = kx 
            do j = 1, n 
               if ( x(jx) /= 0.0 ) then 
                  temp = alpha * x(jx) 
                  ix = jx 
                  do k = kk, kk + n - j 
                     ap(k) = ap( k ) + x(ix) * temp 
                     ix = ix + incx 
                  end do 
               end if 
               jx = jx + incx 
               kk = kk + n - j + 1 
            end do 
         end if 
      end if 
                                                                        
      return 
      END                                           
      subroutine sspr2 ( uplo, n, alpha, x, incx, y, incy, ap ) 
!                                                                       
!***********************************************************************
!                                                                       
!! SSPR2 performs the symmetric rank 2 operation                        
!                                                                       
!     a := alpha * x*y' + alpha * y*x' + a,                             
!                                                                       
!  where alpha is a scalar, x and y are n element vectors and a is an   
!  n by n symmetric matrix, supplied in packed form.                    
!                                                                       
!  Parameters:                                                          
!                                                                       
!  uplo - character.                                                    
!           on entry, uplo specifies whether the upper or lower         
!           triangular part of the matrix a is supplied in the packed   
!           array ap as follows:                                        
!                                                                       
!              uplo = 'U' or 'U'   the upper triangular part of a is    
!                                  supplied in ap.                      
!                                                                       
!              uplo = 'L' or 'L'   the lower triangular part of a is    
!                                  supplied in ap.                      
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  n - integer.                                                         
!           on entry, n specifies the order of the matrix a.            
!           n must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  alpha - real            .                                            
!           on entry, alpha specifies the scalar alpha.                 
!           unchanged on exit.                                          
!                                                                       
!  x - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incx ) ).                            
!           before entry, the incremented array x must contain the n    
!           element vector x.                                           
!           unchanged on exit.                                          
!                                                                       
!  incx - integer.                                                      
!           on entry, incx specifies the increment for the elements of  
!           x. incx must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!  y - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incy ) ).                            
!           before entry, the incremented array y must contain the n    
!           element vector y.                                           
!           unchanged on exit.                                          
!                                                                       
!  incy - integer.                                                      
!           on entry, incy specifies the increment for the elements of  
!           y. incy must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!  ap - real             array of dimension at least                    
!           ( ( n*( n + 1 ) )/2 ).                                      
!           before entry with  uplo = 'U' or 'U', the array ap must     
!           contain the upper triangular part of the symmetric matrix   
!           packed sequentially, column by column, so that ap( 1 )      
!           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )   
!           and a( 2, 2 ) respectively, and so on. on exit, the array   
!           ap is overwritten by the upper triangular part of the       
!           updated matrix.                                             
!           before entry with uplo = 'L' or 'L', the array ap must      
!           contain the lower triangular part of the symmetric matrix   
!           packed sequentially, column by column, so that ap( 1 )      
!           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )   
!           and a( 3, 1 ) respectively, and so on. on exit, the array   
!           ap is overwritten by the lower triangular part of the       
!           updated matrix.                                             
!                                                                       
!                                                                       
!  level 2 blas routine.                                                
!                                                                       
! -- written on 22-october-1986.                                        
!     jack dongarra, argonne national lab.                              
!     jeremy du croz, nag central office.                               
!     sven hammarling, nag central office.                              
!     richard hanson, sandia national labs.                             
!                                                                       
      real               alpha 
      integer            incx, incy, n 
      character        uplo 
!     .. array arguments ..                                             
      real               ap( * ), x( * ), y( * ) 
!     ..                                                                
!                                                                       
!                                                                       
!     .. local scalars ..                                               
      real               temp1, temp2 
      integer            i, info, ix, iy, j, jx, jy, k, kk, kx, ky 
!     .. external functions ..                                          
      logical            lsame 
      external           lsame 
!     .. external subroutines ..                                        
      external           xerbla 
!     ..                                                                
!     .. executable statements ..                                       
!                                                                       
!     test the input parameters.                                        
!                                                                       
      info = 0 
      if     ( .not.lsame ( uplo, 'U' ).and.                            &
     &         .not.lsame ( uplo, 'L' )     ) then                      
         info = 1 
      else if ( n < 0 ) then 
         info = 2 
      else if ( incx == 0 ) then 
         info = 5 
      else if ( incy == 0 ) then 
         info = 7 
      end if 
      if ( info /= 0 ) then 
         call xerbla ( 'sspr2 ', info ) 
         return 
      end if 
!                                                                       
!  Quick return if possible.                                            
!                                                                       
      if ( ( n == 0 ).or.( alpha== 0.0 ) ) then 
        return 
      end if 
!                                                                       
!  Set up the start points in x and y if the increments are not both    
!     unity.                                                            
!                                                                       
      if ( ( incx /= 1 ).or.( incy/=1 ) ) then 
         if ( incx>0 ) then 
            kx = 1 
         else 
            kx = 1 - ( n - 1 ) * incx 
         end if 
         if ( incy>0 ) then 
            ky = 1 
         else 
            ky = 1 - ( n - 1 ) * incy 
         end if 
         jx = kx 
         jy = ky 
      end if 
!                                                                       
!  Start the operations. in this version the elements of the array ap   
!     are accessed sequentially with one pass through ap.               
!                                                                       
      kk = 1 
      if ( lsame ( uplo, 'U' ) ) then 
!                                                                       
!  Form  a  when upper triangle is stored in ap.                        
!                                                                       
         if ( ( incx == 1 ).and.( incy== 1 ) ) then 
            do j = 1, n 
               if ( ( x(j) /= 0.0 ).or.( y(j)/=0.0 ) ) then 
                  temp1 = alpha * y(j) 
                  temp2 = alpha * x(j) 
                  k = kk 
                  do i = 1, j 
                     ap(k) = ap( k ) + x(i) * temp1 + y(i) * temp2 
                     k = k + 1 
                  end do 
               end if 
               kk = kk + j 
            end do 
         else 
            do j = 1, n 
               if ( ( x(jx) /= 0.0 ).or.( y(jy)/=0.0 ) ) then 
                  temp1 = alpha * y(jy) 
                  temp2 = alpha * x(jx) 
                  ix = kx 
                  iy = ky 
                  do k = kk, kk + j - 1 
                     ap(k) = ap( k ) + x(ix) * temp1 + y(iy) * temp2 
                     ix = ix + incx 
                     iy = iy + incy 
                  end do 
               end if 
               jx = jx + incx 
               jy = jy + incy 
               kk = kk + j 
            end do 
         end if 
      else 
!                                                                       
!  Form  a  when lower triangle is stored in ap.                        
!                                                                       
         if ( ( incx == 1 ).and.( incy== 1 ) ) then 
            do j = 1, n 
               if ( ( x(j) /= 0.0 ).or.( y(j)/= 0.0 ) ) then 
                  temp1 = alpha * y(j) 
                  temp2 = alpha * x(j) 
                  k = kk 
                  do i = j, n 
                     ap(k) = ap( k ) + x(i) * temp1 + y(i) * temp2 
                     k = k + 1 
                  end do 
               end if 
               kk = kk + n - j + 1 
            end do 
         else 
            do j = 1, n 
               if ( ( x(jx) /= 0.0 ).or.( y(jy)/= 0.0  ) ) then 
                  temp1 = alpha * y(jy) 
                  temp2 = alpha * x(jx) 
                  ix = jx 
                  iy = jy 
                  do k = kk, kk + n - j 
                     ap(k) = ap( k ) + x(ix) * temp1 + y(iy) * temp2 
                     ix = ix + incx 
                     iy = iy + incy 
                  end do 
               end if 
               jx = jx + incx 
               jy = jy + incy 
               kk = kk + n - j + 1 
            end do 
         end if 
      end if 
                                                                        
      return 
      END                                           
      subroutine ssymv ( uplo, n, alpha, a, lda, x, incx,               &
     &                   beta, y, incy )                                
!                                                                       
!***********************************************************************
!                                                                       
!! SSYMV performs the matrix-vector  operation                          
!                                                                       
!     y := alpha * a * x + beta*y,                                      
!                                                                       
!  where alpha and beta are scalars, x and y are n element vectors and  
!  a is an n by n symmetric matrix.                                     
!                                                                       
!  Parameters:                                                          
!                                                                       
!  uplo - character.                                                    
!           on entry, uplo specifies whether the upper or lower         
!           triangular part of the array a is to be referenced as       
!           follows:                                                    
!                                                                       
!              uplo = 'U' or 'U'   only the upper triangular part of a  
!                                  is to be referenced.                 
!                                                                       
!              uplo = 'L' or 'L'   only the lower triangular part of a  
!                                  is to be referenced.                 
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  n - integer.                                                         
!           on entry, n specifies the order of the matrix a.            
!           n must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  alpha - real            .                                            
!           on entry, alpha specifies the scalar alpha.                 
!           unchanged on exit.                                          
!                                                                       
!  a - real             array of dimension ( lda, n ).                  
!           before entry with  uplo = 'U' or 'U', the leading n by n    
!           upper triangular part of the array a must contain the upper 
!           triangular part of the symmetric matrix and the strictly    
!           lower triangular part of a is not referenced.               
!           before entry with uplo = 'L' or 'L', the leading n by n     
!           lower triangular part of the array a must contain the lower 
!           triangular part of the symmetric matrix and the strictly    
!           upper triangular part of a is not referenced.               
!           unchanged on exit.                                          
!                                                                       
!  lda - integer.                                                       
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least          
!           max ( 1, n ).                                               
!           unchanged on exit.                                          
!                                                                       
!  x - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incx ) ).                            
!           before entry, the incremented array x must contain the n    
!           element vector x.                                           
!           unchanged on exit.                                          
!                                                                       
!  incx - integer.                                                      
!           on entry, incx specifies the increment for the elements of  
!           x. incx must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!  beta - real            .                                             
!           on entry, beta specifies the scalar beta. when beta is      
!           supplied as zero then y need not be set on input.           
!           unchanged on exit.                                          
!                                                                       
!  y - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incy ) ).                            
!           before entry, the incremented array y must contain the n    
!           element vector y. on exit, y is overwritten by the updated  
!           vector y.                                                   
!                                                                       
!  incy - integer.                                                      
!           on entry, incy specifies the increment for the elements of  
!           y. incy must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!                                                                       
!  level 2 blas routine.                                                
!                                                                       
! -- written on 22-october-1986.                                        
!     jack dongarra, argonne national lab.                              
!     jeremy du croz, nag central office.                               
!     sven hammarling, nag central office.                              
!     richard hanson, sandia national labs.                             
!                                                                       
      real               alpha, beta 
      integer            incx, incy, lda, n 
      character        uplo 
!     .. array arguments ..                                             
      real               a( lda, * ), x( * ), y( * ) 
!     ..                                                                
!                                                                       
!     .. parameters ..                                                  
!     .. local scalars ..                                               
      real               temp1, temp2 
      integer            i, info, ix, iy, j, jx, jy, kx, ky 
!     .. external functions ..                                          
      logical            lsame 
      external           lsame 
!     .. external subroutines ..                                        
      external           xerbla 
!     .. intrinsic functions ..                                         
      intrinsic          max 
!     ..                                                                
!     .. executable statements ..                                       
!                                                                       
!     test the input parameters.                                        
!                                                                       
      info = 0 
      if     ( .not.lsame ( uplo, 'U' ).and.                            &
     &         .not.lsame ( uplo, 'L' )     ) then                      
         info = 1 
      else if ( n < 0 ) then 
         info = 2 
      else if ( lda < max ( 1, n ) ) then 
         info = 5 
      else if ( incx == 0 ) then 
         info = 7 
      else if ( incy == 0 ) then 
         info = 10 
      end if 
      if ( info /= 0 ) then 
         call xerbla ( 'ssymv ', info ) 
         return 
      end if 
!                                                                       
!  Quick return if possible.                                            
!                                                                       
      if ( n == 0 ) then 
        return 
      else if ( alpha == 0.0 .and. beta== 1.0 ) then 
        return 
      end if 
!                                                                       
!  Set up the start points in  x  and  y.                               
!                                                                       
      if ( incx>0 ) then 
         kx = 1 
      else 
         kx = 1 - ( n - 1 ) * incx 
      end if 
      if ( incy>0 ) then 
         ky = 1 
      else 
         ky = 1 - ( n - 1 ) * incy 
      end if 
!                                                                       
!  Start the operations. in this version the elements of a are          
!     accessed sequentially with one pass through the triangular part   
!     of a.                                                             
!                                                                       
!     first form  y := beta*y.                                          
!                                                                       
      if ( beta /= 1.0 ) then 
         if ( incy == 1 ) then 
            if ( beta == 0.0 ) then 
               do i = 1, n 
                  y(i) = 0.0 
               end do 
            else 
               do i = 1, n 
                  y(i) = beta*y(i) 
               end do 
            end if 
         else 
            iy = ky 
            if ( beta == 0.0 ) then 
               do i = 1, n 
                  y(iy) = 0.0 
                  iy = iy + incy 
               end do 
            else 
               do i = 1, n 
                  y(iy) = beta*y(iy) 
                  iy = iy + incy 
               end do 
            end if 
         end if 
      end if 
                                                                        
      if ( alpha == 0.0 ) then 
        return 
      end if 
                                                                        
      if ( lsame ( uplo, 'U' ) ) then 
!                                                                       
!  Form  y  when a is stored in upper triangle.                         
!                                                                       
         if ( ( incx == 1 ).and.( incy== 1 ) ) then 
            do j = 1, n 
               temp1 = alpha * x(j) 
               temp2 = 0.0 
               do i = 1, j - 1 
                  y(i) = y(i) + temp1 * a(i,j) 
                  temp2 = temp2 + a(i,j) * x(i) 
               end do 
               y(j) = y(j) + temp1 * a(j,j) + alpha * temp2 
            end do 
         else 
            jx = kx 
            jy = ky 
            do j = 1, n 
               temp1 = alpha * x(jx) 
               temp2 = 0.0 
               ix = kx 
               iy = ky 
               do i = 1, j - 1 
                  y(iy) = y(iy) + temp1 * a(i,j) 
                  temp2 = temp2 + a(i,j) * x(ix) 
                  ix = ix + incx 
                  iy = iy + incy 
               end do 
               y(jy) = y(jy) + temp1 * a(j,j) + alpha * temp2 
               jx = jx + incx 
               jy = jy + incy 
            end do 
         end if 
      else 
!                                                                       
!  Form  y  when a is stored in lower triangle.                         
!                                                                       
         if ( ( incx == 1 ).and.( incy== 1 ) ) then 
            do j = 1, n 
               temp1 = alpha * x(j) 
               temp2 = 0.0 
               y(j) = y(j) + temp1 * a(j,j) 
               do i = j + 1, n 
                  y(i) = y(i) + temp1 * a(i,j) 
                  temp2 = temp2 + a(i,j) * x(i) 
               end do 
               y(j) = y(j) + alpha * temp2 
            end do 
         else 
            jx = kx 
            jy = ky 
            do j = 1, n 
               temp1 = alpha * x(jx) 
               temp2 = 0.0 
               y(jy) = y(jy) + temp1 * a(j,j) 
               ix = jx 
               iy = jy 
               do i = j + 1, n 
                  ix = ix + incx 
                  iy = iy + incy 
                  y(iy) = y(iy) + temp1 * a(i,j) 
                  temp2 = temp2 + a(i,j) * x(ix) 
               end do 
               y(jy) = y(jy) + alpha * temp2 
               jx = jx + incx 
               jy = jy + incy 
            end do 
         end if 
      end if 
                                                                        
      return 
      END                                           
      subroutine ssyr ( uplo, n, alpha, x, incx, a, lda ) 
!                                                                       
!***********************************************************************
!                                                                       
!! SSYR performs the symmetric rank 1 operation                         
!                                                                       
!     a := alpha * x*x' + a,                                            
!                                                                       
!  where alpha is a real scalar, x is an n element vector and a is an   
!  n by n symmetric matrix.                                             
!                                                                       
!  Parameters:                                                          
!                                                                       
!  uplo - character.                                                    
!           on entry, uplo specifies whether the upper or lower         
!           triangular part of the array a is to be referenced as       
!           follows:                                                    
!                                                                       
!              uplo = 'U' or 'U'   only the upper triangular part of a  
!                                  is to be referenced.                 
!                                                                       
!              uplo = 'L' or 'L'   only the lower triangular part of a  
!                                  is to be referenced.                 
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  n - integer.                                                         
!           on entry, n specifies the order of the matrix a.            
!           n must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  alpha - real            .                                            
!           on entry, alpha specifies the scalar alpha.                 
!           unchanged on exit.                                          
!                                                                       
!  x - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incx ) ).                            
!           before entry, the incremented array x must contain the n    
!           element vector x.                                           
!           unchanged on exit.                                          
!                                                                       
!  incx - integer.                                                      
!           on entry, incx specifies the increment for the elements of  
!           x. incx must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!  a - real             array of dimension ( lda, n ).                  
!           before entry with  uplo = 'U' or 'U', the leading n by n    
!           upper triangular part of the array a must contain the upper 
!           triangular part of the symmetric matrix and the strictly    
!           lower triangular part of a is not referenced. on exit, the  
!           upper triangular part of the array a is overwritten by the  
!           upper triangular part of the updated matrix.                
!           before entry with uplo = 'L' or 'L', the leading n by n     
!           lower triangular part of the array a must contain the lower 
!           triangular part of the symmetric matrix and the strictly    
!           upper triangular part of a is not referenced. on exit, the  
!           lower triangular part of the array a is overwritten by the  
!           lower triangular part of the updated matrix.                
!                                                                       
!  lda - integer.                                                       
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least          
!           max ( 1, n ).                                               
!           unchanged on exit.                                          
!                                                                       
!                                                                       
!  level 2 blas routine.                                                
!                                                                       
! -- written on 22-october-1986.                                        
!     jack dongarra, argonne national lab.                              
!     jeremy du croz, nag central office.                               
!     sven hammarling, nag central office.                              
!     richard hanson, sandia national labs.                             
!                                                                       
      real               alpha 
      integer            incx, lda, n 
      character        uplo 
!     .. array arguments ..                                             
      real               a( lda, * ), x( * ) 
!     ..                                                                
!                                                                       
!     .. parameters ..                                                  
!     .. local scalars ..                                               
      real               temp 
      integer            i, info, ix, j, jx, kx 
!     .. external functions ..                                          
      logical            lsame 
      external           lsame 
!     .. external subroutines ..                                        
      external           xerbla 
!     .. intrinsic functions ..                                         
      intrinsic          max 
!     ..                                                                
!     .. executable statements ..                                       
!                                                                       
!     test the input parameters.                                        
!                                                                       
      info = 0 
      if     ( .not.lsame ( uplo, 'U' ).and.                            &
     &         .not.lsame ( uplo, 'L' )     ) then                      
         info = 1 
      else if ( n < 0 ) then 
         info = 2 
      else if ( incx == 0 ) then 
         info = 5 
      else if ( lda < max ( 1, n ) ) then 
         info = 7 
      end if 
      if ( info /= 0 ) then 
         call xerbla ( 'ssyr  ', info ) 
         return 
      end if 
!                                                                       
!  Quick return if possible.                                            
!                                                                       
      if ( ( n == 0 ).or.( alpha== 0.0 ) ) then 
        return 
      end if 
!                                                                       
!     set the start point in x if the increment is not unity.           
!                                                                       
      if ( incx <= 0 ) then 
         kx = 1 - ( n - 1 ) * incx 
      else if ( incx /= 1 ) then 
         kx = 1 
      end if 
!                                                                       
!  Start the operations. in this version the elements of a are          
!     accessed sequentially with one pass through the triangular part   
!     of a.                                                             
!                                                                       
      if ( lsame ( uplo, 'U' ) ) then 
!                                                                       
!  Form  a  when a is stored in upper triangle.                         
!                                                                       
         if ( incx == 1 ) then 
            do j = 1, n 
               if ( x(j) /= 0.0 ) then 
                  temp = alpha * x(j) 
                  do i = 1, j 
                     a(i,j) = a(i,j) + x(i) * temp 
                  end do 
               end if 
            end do 
         else 
            jx = kx 
            do j = 1, n 
               if ( x(jx) /= 0.0 ) then 
                  temp = alpha * x(jx) 
                  ix = kx 
                  do i = 1, j 
                     a(i,j) = a(i,j) + x(ix) * temp 
                     ix = ix + incx 
                  end do 
               end if 
               jx = jx + incx 
            end do 
         end if 
      else 
!                                                                       
!  Form  a  when a is stored in lower triangle.                         
!                                                                       
         if ( incx == 1 ) then 
            do j = 1, n 
               if ( x(j) /= 0.0 ) then 
                  temp = alpha * x(j) 
                  do i = j, n 
                     a(i,j) = a(i,j) + x(i) * temp 
                  end do 
               end if 
            end do 
         else 
            jx = kx 
            do j = 1, n 
               if ( x(jx) /= 0.0 ) then 
                  temp = alpha * x(jx) 
                  ix = jx 
                  do i = j, n 
                     a(i,j) = a(i,j) + x(ix) * temp 
                     ix = ix + incx 
                  end do 
               end if 
               jx = jx + incx 
            end do 
         end if 
      end if 
                                                                        
      return 
      END                                           
      subroutine ssyr2 ( uplo, n, alpha, x, incx, y, incy, a, lda ) 
!                                                                       
!***********************************************************************
!                                                                       
!! SSYR2 performs the symmetric rank 2 operation                        
!                                                                       
!     a := alpha * x*y' + alpha * y*x' + a,                             
!                                                                       
!  where alpha is a scalar, x and y are n element vectors and a is an n 
!  by n symmetric matrix.                                               
!                                                                       
!  Parameters:                                                          
!                                                                       
!  uplo - character.                                                    
!           on entry, uplo specifies whether the upper or lower         
!           triangular part of the array a is to be referenced as       
!           follows:                                                    
!                                                                       
!              uplo = 'U' or 'U'   only the upper triangular part of a  
!                                  is to be referenced.                 
!                                                                       
!              uplo = 'L' or 'L'   only the lower triangular part of a  
!                                  is to be referenced.                 
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  n - integer.                                                         
!           on entry, n specifies the order of the matrix a.            
!           n must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  alpha - real            .                                            
!           on entry, alpha specifies the scalar alpha.                 
!           unchanged on exit.                                          
!                                                                       
!  x - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incx ) ).                            
!           before entry, the incremented array x must contain the n    
!           element vector x.                                           
!           unchanged on exit.                                          
!                                                                       
!  incx - integer.                                                      
!           on entry, incx specifies the increment for the elements of  
!           x. incx must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!  y - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incy ) ).                            
!           before entry, the incremented array y must contain the n    
!           element vector y.                                           
!           unchanged on exit.                                          
!                                                                       
!  incy - integer.                                                      
!           on entry, incy specifies the increment for the elements of  
!           y. incy must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!  a - real             array of dimension ( lda, n ).                  
!           before entry with  uplo = 'U' or 'U', the leading n by n    
!           upper triangular part of the array a must contain the upper 
!           triangular part of the symmetric matrix and the strictly    
!           lower triangular part of a is not referenced. on exit, the  
!           upper triangular part of the array a is overwritten by the  
!           upper triangular part of the updated matrix.                
!           before entry with uplo = 'L' or 'L', the leading n by n     
!           lower triangular part of the array a must contain the lower 
!           triangular part of the symmetric matrix and the strictly    
!           upper triangular part of a is not referenced. on exit, the  
!           lower triangular part of the array a is overwritten by the  
!           lower triangular part of the updated matrix.                
!                                                                       
!  lda - integer.                                                       
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least          
!           max ( 1, n ).                                               
!           unchanged on exit.                                          
!                                                                       
!                                                                       
!  level 2 blas routine.                                                
!                                                                       
! -- written on 22-october-1986.                                        
!     jack dongarra, argonne national lab.                              
!     jeremy du croz, nag central office.                               
!     sven hammarling, nag central office.                              
!     richard hanson, sandia national labs.                             
!                                                                       
      real               alpha 
      integer            incx, incy, lda, n 
      character        uplo 
!     .. array arguments ..                                             
      real               a( lda, * ), x( * ), y( * ) 
!     .. parameters ..                                                  
!     .. local scalars ..                                               
      real               temp1, temp2 
      integer            i, info, ix, iy, j, jx, jy, kx, ky 
!     .. external functions ..                                          
      logical            lsame 
      external           lsame 
!     .. external subroutines ..                                        
      external           xerbla 
!     .. intrinsic functions ..                                         
      intrinsic          max 
!     ..                                                                
!     .. executable statements ..                                       
!                                                                       
!     test the input parameters.                                        
!                                                                       
      info = 0 
      if     ( .not.lsame ( uplo, 'U' ).and.                            &
     &         .not.lsame ( uplo, 'L' )     ) then                      
         info = 1 
      else if ( n < 0 ) then 
         info = 2 
      else if ( incx == 0 ) then 
         info = 5 
      else if ( incy == 0 ) then 
         info = 7 
      else if ( lda < max ( 1, n ) ) then 
         info = 9 
      end if 
      if ( info /= 0 ) then 
         call xerbla ( 'ssyr2 ', info ) 
         return 
      end if 
!                                                                       
!  Quick return if possible.                                            
!                                                                       
      if ( ( n == 0 ).or.( alpha== 0.0 ) ) then 
        return 
      end if 
!                                                                       
!  Set up the start points in x and y if the increments are not both    
!     unity.                                                            
!                                                                       
      if ( ( incx /= 1 ).or.( incy/=1 ) ) then 
         if ( incx>0 ) then 
            kx = 1 
         else 
            kx = 1 - ( n - 1 ) * incx 
         end if 
         if ( incy>0 ) then 
            ky = 1 
         else 
            ky = 1 - ( n - 1 ) * incy 
         end if 
         jx = kx 
         jy = ky 
      end if 
!                                                                       
!  Start the operations. in this version the elements of a are          
!     accessed sequentially with one pass through the triangular part   
!     of a.                                                             
!                                                                       
      if ( lsame ( uplo, 'U' ) ) then 
!                                                                       
!  Form  a  when a is stored in the upper triangle.                     
!                                                                       
         if ( ( incx == 1 ).and.( incy== 1 ) ) then 
            do j = 1, n 
               if ( ( x(j) /= 0.0 ).or.( y(j)/= 0.0  ) ) then 
                  temp1 = alpha * y(j) 
                  temp2 = alpha * x(j) 
                  do i = 1, j 
                     a(i,j) = a(i,j) + x(i) * temp1 + y(i) * temp2 
                  end do 
               end if 
            end do 
         else 
            do j = 1, n 
               if ( ( x(jx) /= 0.0 ).or.( y(jy)/= 0.0  ) ) then 
                  temp1 = alpha * y(jy) 
                  temp2 = alpha * x(jx) 
                  ix = kx 
                  iy = ky 
                  do i = 1, j 
                     a(i,j) = a(i,j) + x(ix) * temp1 + y(iy) * temp2 
                     ix = ix + incx 
                     iy = iy + incy 
                  end do 
               end if 
               jx = jx + incx 
               jy = jy + incy 
            end do 
         end if 
      else 
!                                                                       
!  Form  a  when a is stored in the lower triangle.                     
!                                                                       
         if ( ( incx == 1 ).and.( incy== 1 ) ) then 
            do j = 1, n 
               if ( ( x(j) /= 0.0 ).or.( y(j)/= 0.0  ) ) then 
                  temp1 = alpha * y(j) 
                  temp2 = alpha * x(j) 
                  do i = j, n 
                     a(i,j) = a(i,j) + x(i) * temp1 + y(i) * temp2 
                  end do 
               end if 
            end do 
         else 
            do j = 1, n 
               if ( ( x(jx) /= 0.0 ).or.( y(jy)/= 0.0  ) ) then 
                  temp1 = alpha * y(jy) 
                  temp2 = alpha * x(jx) 
                  ix = jx 
                  iy = jy 
                  do i = j, n 
                     a(i,j) = a(i,j) + x(ix) * temp1 + y(iy) * temp2 
                     ix = ix + incx 
                     iy = iy + incy 
                  end do 
               end if 
               jx = jx + incx 
               jy = jy + incy 
            end do 
         end if 
      end if 
                                                                        
      return 
      END                                           
      subroutine stbmv ( uplo, trans, diag, n, k, a, lda, x, incx ) 
!                                                                       
!***********************************************************************
!                                                                       
!! STBMV performs one of the matrix-vector operations                   
!                                                                       
!     x := a * x,   or   x := a'*x,                                     
!                                                                       
!  where x is an n element vector and  a is an n by n unit, or non-unit,
!  upper or lower triangular band matrix, with ( k + 1 ) diagonals.     
!                                                                       
!  Parameters:                                                          
!                                                                       
!  uplo - character.                                                    
!           on entry, uplo specifies whether the matrix is an upper or  
!           lower triangular matrix as follows:                         
!                                                                       
!              uplo = 'U' or 'U'   a is an upper triangular matrix.     
!                                                                       
!              uplo = 'L' or 'L'   a is a lower triangular matrix.      
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  trans - character.                                                   
!           on entry, trans specifies the operation to be performed as  
!           follows:                                                    
!                                                                       
!              trans = 'N' or 'N'   x := a * x.                         
!                                                                       
!              trans = 'T' or 'T'   x := a' * x.                        
!                                                                       
!              trans = 'C' or 'C'   x := a' * x.                        
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  diag - character.                                                    
!           on entry, diag specifies whether or not a is unit           
!           triangular as follows:                                      
!                                                                       
!              diag = 'U' or 'U'   a is assumed to be unit triangular.  
!                                                                       
!              diag = 'N' or 'N'   a is not assumed to be unit          
!                                  triangular.                          
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  n - integer.                                                         
!           on entry, n specifies the order of the matrix a.            
!           n must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  k - integer.                                                         
!           on entry with uplo = 'U' or 'U', k specifies the number of  
!           super-diagonals of the matrix a.                            
!           on entry with uplo = 'L' or 'L', k specifies the number of  
!           sub-diagonals of the matrix a.                              
!           k must satisfy  0  <= k.                                    
!           unchanged on exit.                                          
!                                                                       
!  a - real             array of dimension ( lda, n ).                  
!           before entry with uplo = 'U' or 'U', the leading ( k + 1 )  
!           by n part of the array a must contain the upper triangular  
!           band part of the matrix of coefficients, supplied column by 
!           column, with the leading diagonal of the matrix in row      
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. the top left k by k triangle
!           of the array a is not referenced.                           
!           the following program segment will transfer an upper        
!           triangular band matrix from conventional full matrix storage
!           to band storage:                                            
!                                                                       
!                 do j = 1, n                                           
!                    m = k + 1 - j                                      
!                    do i = max ( 1, j - k ), j                         
!                       a( m + i,j) = matrix(i,j)                       
!                    end do                                             
!                 end do                                                
!                                                                       
!           before entry with uplo = 'L' or 'L', the leading ( k + 1 )  
!           by n part of the array a must contain the lower triangular  
!           band part of the matrix of coefficients, supplied column by 
!           column, with the leading diagonal of the matrix in row 1 of 
!           the array, the first sub-diagonal starting at position 1 in 
!           row 2, and so on. the bottom right k by k triangle of the   
!           array a is not referenced.                                  
!           the following program segment will transfer a lower         
!           triangular band matrix from conventional full matrix storage
!           to band storage:                                            
!                                                                       
!                 do j = 1, n                                           
!                    m = 1 - j                                          
!                    do i = j, min ( n, j + k )                         
!                       a( m + i,j) = matrix(i,j)                       
!                    end do                                             
!                 end do                                                
!                                                                       
!           note that when diag = 'U' or 'U' the elements of the array a
!           corresponding to the diagonal elements of the matrix are not
!           referenced, but are assumed to be unity.                    
!           unchanged on exit.                                          
!                                                                       
!  lda - integer.                                                       
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least          
!           ( k + 1 ).                                                  
!           unchanged on exit.                                          
!                                                                       
!  x - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incx ) ).                            
!           before entry, the incremented array x must contain the n    
!           element vector x. on exit, x is overwritten with the        
!           tranformed vector x.                                        
!                                                                       
!  incx - integer.                                                      
!           on entry, incx specifies the increment for the elements of  
!           x. incx must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!                                                                       
!  level 2 blas routine.                                                
!                                                                       
! -- written on 22-october-1986.                                        
!     jack dongarra, argonne national lab.                              
!     jeremy du croz, nag central office.                               
!     sven hammarling, nag central office.                              
!     richard hanson, sandia national labs.                             
!                                                                       
      integer            incx, k, lda, n 
      character        diag, trans, uplo 
!     .. array arguments ..                                             
      real               a( lda, * ), x( * ) 
!     ..                                                                
!     .. parameters ..                                                  
!     .. local scalars ..                                               
      real               temp 
      integer            i, info, ix, j, jx, kplus1, kx, l 
      logical            nounit 
!     .. external functions ..                                          
      logical            lsame 
      external           lsame 
!     .. external subroutines ..                                        
      external           xerbla 
!     .. intrinsic functions ..                                         
      intrinsic          max, min 
!     ..                                                                
!     .. executable statements ..                                       
!                                                                       
!     test the input parameters.                                        
!                                                                       
      info = 0 
      if     ( .not.lsame ( uplo , 'U' ).and.                           &
     &         .not.lsame ( uplo , 'L' )     ) then                     
         info = 1 
      else if ( .not.lsame ( trans, 'N' ).and.                          &
     &         .not.lsame ( trans, 'T' ).and.                           &
     &         .not.lsame ( trans, 'C' )     ) then                     
         info = 2 
      else if ( .not.lsame ( diag , 'U' ).and.                          &
     &         .not.lsame ( diag , 'N' )     ) then                     
         info = 3 
      else if ( n < 0 ) then 
         info = 4 
      else if ( k < 0 ) then 
         info = 5 
      else if ( lda < ( k + 1 ) ) then 
         info = 7 
      else if ( incx == 0 ) then 
         info = 9 
      end if 
      if ( info /= 0 ) then 
         call xerbla ( 'stbmv ', info ) 
         return 
      end if 
!                                                                       
!  Quick return if possible.                                            
!                                                                       
      if ( n == 0 ) then 
        return 
      end if 
!                                                                       
      nounit = lsame ( diag, 'N' ) 
!                                                                       
!  Set up the start point in x if the increment is not unity. this      
!     will be  ( n - 1 ) * incx   too small for descending loops.       
!                                                                       
      if ( incx <= 0 ) then 
         kx = 1 - ( n - 1 ) * incx 
      else if ( incx /= 1 ) then 
         kx = 1 
      end if 
!                                                                       
!  Start the operations. in this version the elements of a are          
!     accessed sequentially with one pass through a.                    
!                                                                       
      if ( lsame ( trans, 'N' ) ) then 
!                                                                       
!         form  x := a * x.                                             
!                                                                       
         if ( lsame ( uplo, 'U' ) ) then 
            kplus1 = k + 1 
            if ( incx == 1 ) then 
               do j = 1, n 
                  if ( x(j) /= 0.0 ) then 
                     temp = x(j) 
                     l = kplus1 - j 
                     do i = max ( 1, j - k ), j - 1 
                        x(i) = x(i) + temp * a( l + i,j) 
                     end do 
                                                                        
                     if ( nounit ) then 
                       x(j) = x(j) * a( kplus1,j) 
                     end if 
                                                                        
                  end if 
               end do 
            else 
               jx = kx 
               do j = 1, n 
                  if ( x(jx) /= 0.0 ) then 
                     temp = x(jx) 
                     ix = kx 
                     l = kplus1 - j 
                     do i = max ( 1, j - k ), j - 1 
                        x(ix) = x(ix) + temp * a( l + i,j) 
                        ix = ix + incx 
                     end do 
                                                                        
                     if ( nounit ) then 
                       x(jx) = x(jx) * a( kplus1,j) 
                     end if 
                                                                        
                  end if 
                  jx = jx + incx 
                                                                        
                  if ( j>k ) then 
                    kx = kx + incx 
                  end if 
                                                                        
               end do 
            end if 
         else 
            if ( incx == 1 ) then 
               do j = n, 1, -1 
                  if ( x(j) /= 0.0 ) then 
                     temp = x(j) 
                     l = 1 - j 
                     do i = min ( n, j + k ), j + 1, -1 
                        x(i) = x(i) + temp * a( l + i,j) 
                     end do 
                                                                        
                     if ( nounit ) then 
                       x(j) = x(j) * a( 1,j) 
                     end if 
                                                                        
                  end if 
               end do 
            else 
               kx = kx + ( n - 1 ) * incx 
               jx = kx 
               do j = n, 1, -1 
                  if ( x(jx) /= 0.0 ) then 
                     temp = x(jx) 
                     ix = kx 
                     l = 1 - j 
                     do i = min ( n, j + k ), j + 1, -1 
                        x(ix) = x(ix) + temp * a( l + i,j) 
                        ix = ix - incx 
                     end do 
                                                                        
                     if ( nounit ) then 
                       x(jx) = x(jx) * a(1,j) 
                     end if 
                                                                        
                  end if 
                  jx = jx - incx 
                                                                        
                  if ( ( n - j ) >= k ) then 
                    kx = kx - incx 
                  end if 
                                                                        
               end do 
            end if 
         end if 
      else 
!                                                                       
!  Form  x := a' * x.                                                   
!                                                                       
         if ( lsame ( uplo, 'U' ) ) then 
            kplus1 = k + 1 
            if ( incx == 1 ) then 
               do j = n, 1, -1 
                  temp = x(j) 
                  l = kplus1 - j 
                                                                        
                  if ( nounit ) then 
                    temp = temp * a( kplus1,j) 
                  end if 
                                                                        
                  do i = j - 1, max ( 1, j - k ), -1 
                     temp = temp + a( l + i,j) * x(i) 
                  end do 
                  x(j) = temp 
               end do 
            else 
               kx = kx + ( n - 1 ) * incx 
               jx = kx 
               do j = n, 1, -1 
                  temp = x(jx) 
                  kx = kx - incx 
                  ix = kx 
                  l = kplus1 - j 
                                                                        
                  if ( nounit ) then 
                    temp = temp * a( kplus1,j) 
                  end if 
                                                                        
                  do i = j - 1, max ( 1, j - k ), -1 
                     temp = temp + a( l + i,j) * x(ix) 
                     ix = ix - incx 
                  end do 
                  x(jx) = temp 
                  jx = jx - incx 
               end do 
            end if 
         else 
            if ( incx == 1 ) then 
               do j = 1, n 
                  temp = x(j) 
                  l = 1 - j 
                                                                        
                  if ( nounit ) then 
                    temp = temp * a( 1,j) 
                  end if 
                                                                        
                  do i = j + 1, min ( n, j + k ) 
                     temp = temp + a( l + i,j) * x(i) 
                  end do 
                  x(j) = temp 
               end do 
            else 
               jx = kx 
               do j = 1, n 
                  temp = x(jx) 
                  kx = kx + incx 
                  ix = kx 
                  l = 1 - j 
                                                                        
                  if ( nounit ) then 
                    temp = temp * a( 1,j) 
                  end if 
                                                                        
                  do i = j + 1, min ( n, j + k ) 
                     temp = temp + a( l + i,j) * x(ix) 
                     ix = ix + incx 
                  end do 
                  x(jx) = temp 
                  jx = jx + incx 
               end do 
            end if 
         end if 
      end if 
                                                                        
      return 
      END                                           
      subroutine stbsv ( uplo, trans, diag, n, k, a, lda, x, incx ) 
!                                                                       
!***********************************************************************
!                                                                       
!! STBSV solves one of the systems of equations                         
!                                                                       
!     a * x = b,   or   a'*x = b,                                       
!                                                                       
!  where b and x are n element vectors and a is an n by n unit, or      
!  non-unit, upper or lower triangular band matrix, with ( k + 1 )      
!  diagonals.                                                           
!                                                                       
!  no test for singularity or near-singularity is included in this      
!  routine. such tests must be performed before calling this routine.   
!                                                                       
!  Parameters:                                                          
!                                                                       
!  uplo - character.                                                    
!           on entry, uplo specifies whether the matrix is an upper or  
!           lower triangular matrix as follows:                         
!                                                                       
!              uplo = 'U' or 'U'   a is an upper triangular matrix.     
!                                                                       
!              uplo = 'L' or 'L'   a is a lower triangular matrix.      
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  trans - character.                                                   
!           on entry, trans specifies the equations to be solved as     
!           follows:                                                    
!                                                                       
!              trans = 'N' or 'N'   a * x = b.                          
!                                                                       
!              trans = 'T' or 'T'   a' * x = b.                         
!                                                                       
!              trans = 'C' or 'C'   a' * x = b.                         
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  diag - character.                                                    
!           on entry, diag specifies whether or not a is unit           
!           triangular as follows:                                      
!                                                                       
!              diag = 'U' or 'U'   a is assumed to be unit triangular.  
!                                                                       
!              diag = 'N' or 'N'   a is not assumed to be unit          
!                                  triangular.                          
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  n - integer.                                                         
!           on entry, n specifies the order of the matrix a.            
!           n must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  k - integer.                                                         
!           on entry with uplo = 'U' or 'U', k specifies the number of  
!           super-diagonals of the matrix a.                            
!           on entry with uplo = 'L' or 'L', k specifies the number of  
!           sub-diagonals of the matrix a.                              
!           k must satisfy  0  <= k.                                    
!           unchanged on exit.                                          
!                                                                       
!  a - real             array of dimension ( lda, n ).                  
!           before entry with uplo = 'U' or 'U', the leading ( k + 1 )  
!           by n part of the array a must contain the upper triangular  
!           band part of the matrix of coefficients, supplied column by 
!           column, with the leading diagonal of the matrix in row      
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. the top left k by k triangle
!           of the array a is not referenced.                           
!           the following program segment will transfer an upper        
!           triangular band matrix from conventional full matrix storage
!           to band storage:                                            
!                                                                       
!                 do j = 1, n                                           
!                    m = k + 1 - j                                      
!                    do i = max ( 1, j - k ), j                         
!                       a( m + i,j) = matrix(i,j)                       
!                    end do                                             
!                 end do                                                
!                                                                       
!           before entry with uplo = 'L' or 'L', the leading ( k + 1 )  
!           by n part of the array a must contain the lower triangular  
!           band part of the matrix of coefficients, supplied column by 
!           column, with the leading diagonal of the matrix in row 1 of 
!           the array, the first sub-diagonal starting at position 1 in 
!           row 2, and so on. the bottom right k by k triangle of the   
!           array a is not referenced.                                  
!           the following program segment will transfer a lower         
!           triangular band matrix from conventional full matrix storage
!           to band storage:                                            
!                                                                       
!                 do j = 1, n                                           
!                    m = 1 - j                                          
!                    do i = j, min ( n, j + k )                         
!                       a( m + i,j) = matrix(i,j)                       
!                    end do                                             
!                 end do                                                
!                                                                       
!           note that when diag = 'U' or 'U' the elements of the array a
!           corresponding to the diagonal elements of the matrix are not
!           referenced, but are assumed to be unity.                    
!           unchanged on exit.                                          
!                                                                       
!  lda - integer.                                                       
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least          
!           ( k + 1 ).                                                  
!           unchanged on exit.                                          
!                                                                       
!  x - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incx ) ).                            
!           before entry, the incremented array x must contain the n    
!           element right-hand side vector b. on exit, x is overwritten 
!           with the solution vector x.                                 
!                                                                       
!  incx - integer.                                                      
!           on entry, incx specifies the increment for the elements of  
!           x. incx must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!                                                                       
!  level 2 blas routine.                                                
!                                                                       
! -- written on 22-october-1986.                                        
!     jack dongarra, argonne national lab.                              
!     jeremy du croz, nag central office.                               
!     sven hammarling, nag central office.                              
!     richard hanson, sandia national labs.                             
!                                                                       
      integer            incx, k, lda, n 
      character        diag, trans, uplo 
!     .. array arguments ..                                             
      real               a( lda, * ), x( * ) 
!     ..                                                                
!                                                                       
!     .. parameters ..                                                  
!     .. local scalars ..                                               
      real               temp 
      integer            i, info, ix, j, jx, kplus1, kx, l 
      logical            nounit 
!     .. external functions ..                                          
      logical            lsame 
      external           lsame 
!     .. external subroutines ..                                        
      external           xerbla 
!     .. intrinsic functions ..                                         
      intrinsic          max, min 
!     ..                                                                
!     .. executable statements ..                                       
!                                                                       
!     test the input parameters.                                        
!                                                                       
      info = 0 
      if     ( .not.lsame ( uplo , 'U' ).and.                           &
     &         .not.lsame ( uplo , 'L' )     ) then                     
         info = 1 
      else if ( .not.lsame ( trans, 'N' ).and.                          &
     &         .not.lsame ( trans, 'T' ).and.                           &
     &         .not.lsame ( trans, 'C' )     ) then                     
         info = 2 
      else if ( .not.lsame ( diag , 'U' ).and.                          &
     &         .not.lsame ( diag , 'N' )     ) then                     
         info = 3 
      else if ( n < 0 ) then 
         info = 4 
      else if ( k < 0 ) then 
         info = 5 
      else if ( lda < ( k + 1 ) ) then 
         info = 7 
      else if ( incx == 0 ) then 
         info = 9 
      end if 
      if ( info /= 0 ) then 
         call xerbla ( 'stbsv ', info ) 
         return 
      end if 
!                                                                       
!  Quick return if possible.                                            
!                                                                       
      if ( n == 0 ) then 
        return 
      end if 
!                                                                       
      nounit = lsame ( diag, 'N' ) 
!                                                                       
!  Set up the start point in x if the increment is not unity. this      
!     will be  ( n - 1 ) * incx  too small for descending loops.        
!                                                                       
      if ( incx <= 0 ) then 
         kx = 1 - ( n - 1 ) * incx 
      else if ( incx /= 1 ) then 
         kx = 1 
      end if 
!                                                                       
!  Start the operations. in this version the elements of a are          
!     accessed by sequentially with one pass through a.                 
!                                                                       
      if ( lsame ( trans, 'N' ) ) then 
!                                                                       
!  Form  x := inv( a ) * x.                                             
!                                                                       
         if ( lsame ( uplo, 'U' ) ) then 
            kplus1 = k + 1 
            if ( incx == 1 ) then 
               do j = n, 1, -1 
                  if ( x(j) /= 0.0 ) then 
                     l = kplus1 - j 
                                                                        
                     if ( nounit ) then 
                       x(j) = x(j) / a(kplus1,j) 
                     end if 
                                                                        
                     temp = x(j) 
                     do i = j - 1, max ( 1, j - k ), -1 
                        x(i) = x(i) - temp * a(l+i,j) 
                     end do 
                  end if 
               end do 
            else 
               kx = kx + ( n - 1 ) * incx 
               jx = kx 
               do j = n, 1, -1 
                  kx = kx - incx 
                  if ( x(jx) /= 0.0 ) then 
                     ix = kx 
                     l = kplus1 - j 
                     if ( nounit ) then 
                       x(jx) = x(jx) / a(kplus1,j) 
                     end if 
                     temp = x(jx) 
                     do i = j - 1, max ( 1, j - k ), -1 
                        x(ix) = x(ix) - temp * a(l+i,j) 
                        ix = ix - incx 
                     end do 
                  end if 
                  jx = jx - incx 
               end do 
            end if 
         else 
            if ( incx == 1 ) then 
               do j = 1, n 
                  if ( x(j) /= 0.0 ) then 
                     l = 1 - j 
                     if ( nounit ) then 
                       x(j) = x(j) / a(1,j) 
                     end if 
                     temp = x(j) 
                     do i = j + 1, min ( n, j + k ) 
                        x(i) = x(i) - temp * a(l+i,j) 
                     end do 
                  end if 
               end do 
            else 
               jx = kx 
               do j = 1, n 
                  kx = kx + incx 
                  if ( x(jx) /= 0.0 ) then 
                     ix = kx 
                     l = 1 - j 
                     if ( nounit ) then 
                       x(jx) = x(jx) / a( 1,j) 
                     end if 
                     temp = x(jx) 
                     do i = j + 1, min ( n, j + k ) 
                        x(ix) = x(ix) - temp * a( l + i,j) 
                        ix = ix + incx 
                     end do 
                  end if 
                  jx = jx + incx 
               end do 
            end if 
         end if 
      else 
!                                                                       
!  Form  x := inv( a') * x.                                             
!                                                                       
         if ( lsame ( uplo, 'U' ) ) then 
            kplus1 = k + 1 
            if ( incx == 1 ) then 
               do j = 1, n 
                  temp = x(j) 
                  l = kplus1 - j 
                  do i = max ( 1, j - k ), j - 1 
                     temp = temp - a( l + i,j) * x(i) 
                  end do 
                  if ( nounit ) then 
                    temp = temp/a( kplus1,j) 
                  end if 
                  x(j) = temp 
               end do 
            else 
               jx = kx 
               do j = 1, n 
                  temp = x(jx) 
                  ix = kx 
                  l = kplus1 - j 
                  do i = max ( 1, j - k ), j - 1 
                     temp = temp - a( l + i,j) * x(ix) 
                     ix = ix + incx 
                  end do 
                  if ( nounit ) then 
                    temp = temp/a( kplus1,j) 
                  end if 
                  x(jx) = temp 
                  jx = jx + incx 
                  if ( j>k ) then 
                    kx = kx + incx 
                  end if 
               end do 
            end if 
         else 
            if ( incx == 1 ) then 
               do j = n, 1, -1 
                  temp = x(j) 
                  l = 1 - j 
                  do i = min ( n, j + k ), j + 1, -1 
                     temp = temp - a( l + i,j) * x(i) 
                  end do 
                  if ( nounit ) then 
                    temp = temp/a( 1,j) 
                  end if 
                  x(j) = temp 
               end do 
            else 
               kx = kx + ( n - 1 ) * incx 
               jx = kx 
               do j = n, 1, -1 
                  temp = x(jx) 
                  ix = kx 
                  l = 1 - j 
                  do i = min ( n, j + k ), j + 1, -1 
                     temp = temp - a( l + i,j) * x(ix) 
                     ix = ix - incx 
                  end do 
                  if ( nounit ) then 
                    temp = temp/a( 1,j) 
                  end if 
                  x(jx) = temp 
                  jx = jx - incx 
                  if ( ( n - j ) >= k ) then 
                    kx = kx - incx 
                  end if 
               end do 
            end if 
         end if 
      end if 
                                                                        
      return 
      END                                           
      subroutine stpmv ( uplo, trans, diag, n, ap, x, incx ) 
!                                                                       
!***********************************************************************
!                                                                       
!! STPMV performs one of the matrix-vector operations                   
!                                                                       
!     x := a * x,   or   x := a'*x,                                     
!                                                                       
!  where x is an n element vector and  a is an n by n unit, or non-unit,
!  upper or lower triangular matrix, supplied in packed form.           
!                                                                       
!  Parameters:                                                          
!                                                                       
!  uplo - character.                                                    
!           on entry, uplo specifies whether the matrix is an upper or  
!           lower triangular matrix as follows:                         
!                                                                       
!              uplo = 'U' or 'U'   a is an upper triangular matrix.     
!                                                                       
!              uplo = 'L' or 'L'   a is a lower triangular matrix.      
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  trans - character.                                                   
!           on entry, trans specifies the operation to be performed as  
!           follows:                                                    
!                                                                       
!              trans = 'N' or 'N'   x := a * x.                         
!                                                                       
!              trans = 'T' or 'T'   x := a' * x.                        
!                                                                       
!              trans = 'C' or 'C'   x := a' * x.                        
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  diag - character.                                                    
!           on entry, diag specifies whether or not a is unit           
!           triangular as follows:                                      
!                                                                       
!              diag = 'U' or 'U'   a is assumed to be unit triangular.  
!                                                                       
!              diag = 'N' or 'N'   a is not assumed to be unit          
!                                  triangular.                          
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  n - integer.                                                         
!           on entry, n specifies the order of the matrix a.            
!           n must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  ap - real             array of dimension at least                    
!           ( ( n*( n + 1 ) )/2 ).                                      
!           before entry with  uplo = 'U' or 'U', the array ap must     
!           contain the upper triangular matrix packed sequentially,    
!           column by column, so that ap( 1 ) contains a( 1, 1 ),       
!           ap( 2 ) and ap( 3 ) contain a( 1, 2 ) and a( 2, 2 )         
!           respectively, and so on.                                    
!           before entry with uplo = 'L' or 'L', the array ap must      
!           contain the lower triangular matrix packed sequentially,    
!           column by column, so that ap( 1 ) contains a( 1, 1 ),       
!           ap( 2 ) and ap( 3 ) contain a( 2, 1 ) and a( 3, 1 )         
!           respectively, and so on.                                    
!           note that when  diag = 'U' or 'U', the diagonal elements of 
!           a are not referenced, but are assumed to be unity.          
!           unchanged on exit.                                          
!                                                                       
!  x - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incx ) ).                            
!           before entry, the incremented array x must contain the n    
!           element vector x. on exit, x is overwritten with the        
!           tranformed vector x.                                        
!                                                                       
!  incx - integer.                                                      
!           on entry, incx specifies the increment for the elements of  
!           x. incx must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!                                                                       
!  level 2 blas routine.                                                
!                                                                       
! -- written on 22-october-1986.                                        
!     jack dongarra, argonne national lab.                              
!     jeremy du croz, nag central office.                               
!     sven hammarling, nag central office.                              
!     richard hanson, sandia national labs.                             
!                                                                       
      integer            incx, n 
      character        diag, trans, uplo 
!     .. array arguments ..                                             
      real               ap( * ), x( * ) 
!     ..                                                                
!                                                                       
!     .. parameters ..                                                  
!     .. local scalars ..                                               
      real               temp 
      integer            i, info, ix, j, jx, k, kk, kx 
      logical            nounit 
!     .. external functions ..                                          
      logical            lsame 
      external           lsame 
!     .. external subroutines ..                                        
      external           xerbla 
!     ..                                                                
!     .. executable statements ..                                       
!                                                                       
!     test the input parameters.                                        
!                                                                       
      info = 0 
      if     ( .not.lsame ( uplo , 'U' ).and.                           &
     &         .not.lsame ( uplo , 'L' )     ) then                     
         info = 1 
      else if ( .not.lsame ( trans, 'N' ).and.                          &
     &         .not.lsame ( trans, 'T' ).and.                           &
     &         .not.lsame ( trans, 'C' )     ) then                     
         info = 2 
      else if ( .not.lsame ( diag , 'U' ).and.                          &
     &         .not.lsame ( diag , 'N' )     ) then                     
         info = 3 
      else if ( n < 0 ) then 
         info = 4 
      else if ( incx == 0 ) then 
         info = 7 
      end if 
      if ( info /= 0 ) then 
         call xerbla ( 'stpmv ', info ) 
         return 
      end if 
!                                                                       
!  Quick return if possible.                                            
!                                                                       
      if ( n == 0 ) then 
        return 
      end if 
!                                                                       
      nounit = lsame ( diag, 'N' ) 
!                                                                       
!  Set up the start point in x if the increment is not unity. this      
!     will be  ( n - 1 ) * incx  too small for descending loops.        
!                                                                       
      if ( incx <= 0 ) then 
         kx = 1 - ( n - 1 ) * incx 
      else if ( incx /= 1 ) then 
         kx = 1 
      end if 
!                                                                       
!  Start the operations. in this version the elements of ap are         
!     accessed sequentially with one pass through ap.                   
!                                                                       
      if ( lsame ( trans, 'N' ) ) then 
!                                                                       
!  Form  x:= a * x.                                                     
!                                                                       
         if ( lsame ( uplo, 'U' ) ) then 
            kk = 1 
            if ( incx == 1 ) then 
               do j = 1, n 
                  if ( x(j) /= 0.0 ) then 
                     temp = x(j) 
                     k = kk 
                     do i = 1, j - 1 
                        x(i) = x(i) + temp * ap(k) 
                        k = k + 1 
                     end do 
                     if ( nounit ) then 
                       x(j) = x(j) * ap( kk + j - 1 ) 
                     end if 
                  end if 
                  kk = kk + j 
               end do 
            else 
               jx = kx 
               do j = 1, n 
                  if ( x(jx) /= 0.0 ) then 
                     temp = x(jx) 
                     ix = kx 
                     do k = kk, kk + j - 2 
                        x(ix) = x(ix) + temp * ap(k) 
                        ix = ix + incx 
                     end do 
                     if ( nounit ) then 
                       x(jx) = x(jx) * ap( kk + j - 1 ) 
                     end if 
                  end if 
                  jx = jx + incx 
                  kk = kk + j 
               end do 
            end if 
         else 
            kk = ( n*( n + 1 ) )/2 
            if ( incx == 1 ) then 
               do j = n, 1, -1 
                  if ( x(j) /= 0.0 ) then 
                     temp = x(j) 
                     k = kk 
                     do i = n, j + 1, -1 
                        x(i) = x(i) + temp * ap(k) 
                        k = k - 1 
                     end do 
                     if ( nounit ) then 
                       x(j) = x(j) * ap( kk - n + j ) 
                     end if 
                  end if 
                  kk = kk - ( n - j + 1 ) 
               end do 
            else 
               kx = kx + ( n - 1 ) * incx 
               jx = kx 
               do j = n, 1, -1 
                  if ( x(jx) /= 0.0 ) then 
                     temp = x(jx) 
                     ix = kx 
                     do k = kk, kk - ( n - ( j + 1 ) ), -1 
                        x(ix) = x(ix) + temp * ap(k) 
                        ix = ix - incx 
                     end do 
                     if ( nounit ) then 
                       x(jx) = x(jx) * ap( kk - n + j ) 
                     end if 
                  end if 
                  jx = jx - incx 
                  kk = kk - ( n - j + 1 ) 
               end do 
            end if 
         end if 
      else 
!                                                                       
!  Form  x := a' * x.                                                   
!                                                                       
         if ( lsame ( uplo, 'U' ) ) then 
            kk = ( n*( n + 1 ) )/2 
            if ( incx == 1 ) then 
               do j = n, 1, -1 
                  temp = x(j) 
                  if ( nounit ) then 
                    temp = temp * ap(kk) 
                  end if 
                  k = kk - 1 
                  do i = j - 1, 1, -1 
                     temp = temp + ap(k) * x(i) 
                     k = k - 1 
                  end do 
                  x(j) = temp 
                  kk = kk - j 
               end do 
            else 
               jx = kx + ( n - 1 ) * incx 
               do j = n, 1, -1 
                  temp = x(jx) 
                  ix = jx 
                  if ( nounit ) then 
                    temp = temp * ap(kk) 
                  end if 
                  do k = kk - 1, kk - j + 1, -1 
                     ix = ix - incx 
                     temp = temp + ap(k) * x(ix) 
                  end do 
                  x(jx) = temp 
                  jx = jx - incx 
                  kk = kk - j 
               end do 
            end if 
         else 
            kk = 1 
            if ( incx == 1 ) then 
               do j = 1, n 
                  temp = x(j) 
                  if ( nounit ) then 
                    temp = temp * ap(kk) 
                  end if 
                  k = kk + 1 
                  do i = j + 1, n 
                     temp = temp + ap(k) * x(i) 
                     k = k + 1 
                  end do 
                  x(j) = temp 
                  kk = kk + ( n - j + 1 ) 
               end do 
            else 
               jx = kx 
               do j = 1, n 
                  temp = x(jx) 
                  ix = jx 
                  if ( nounit ) then 
                    temp = temp * ap(kk) 
                  end if 
                  do k = kk + 1, kk + n - j 
                     ix = ix + incx 
                     temp = temp + ap(k) * x(ix) 
                  end do 
                  x(jx) = temp 
                  jx = jx + incx 
                  kk = kk + ( n - j + 1 ) 
               end do 
            end if 
         end if 
      end if 
                                                                        
      return 
      END                                           
      subroutine stpsv ( uplo, trans, diag, n, ap, x, incx ) 
!                                                                       
!***********************************************************************
!                                                                       
!! STPSV solves one of the systems of equations                         
!                                                                       
!     a * x = b,   or   a'*x = b,                                       
!                                                                       
!  where b and x are n element vectors and a is an n by n unit, or      
!  non-unit, upper or lower triangular matrix, supplied in packed form. 
!                                                                       
!  no test for singularity or near-singularity is included in this      
!  routine. such tests must be performed before calling this routine.   
!                                                                       
!  Parameters:                                                          
!                                                                       
!  uplo - character.                                                    
!           on entry, uplo specifies whether the matrix is an upper or  
!           lower triangular matrix as follows:                         
!                                                                       
!              uplo = 'U' or 'U'   a is an upper triangular matrix.     
!                                                                       
!              uplo = 'L' or 'L'   a is a lower triangular matrix.      
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  trans - character.                                                   
!           on entry, trans specifies the equations to be solved as     
!           follows:                                                    
!                                                                       
!              trans = 'N' or 'N'   a * x = b.                          
!                                                                       
!              trans = 'T' or 'T'   a' * x = b.                         
!                                                                       
!              trans = 'C' or 'C'   a' * x = b.                         
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  diag - character.                                                    
!           on entry, diag specifies whether or not a is unit           
!           triangular as follows:                                      
!                                                                       
!              diag = 'U' or 'U'   a is assumed to be unit triangular.  
!                                                                       
!              diag = 'N' or 'N'   a is not assumed to be unit          
!                                  triangular.                          
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  n - integer.                                                         
!           on entry, n specifies the order of the matrix a.            
!           n must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  ap - real             array of dimension at least                    
!           ( ( n*( n + 1 ) )/2 ).                                      
!           before entry with  uplo = 'U' or 'U', the array ap must     
!           contain the upper triangular matrix packed sequentially,    
!           column by column, so that ap( 1 ) contains a( 1, 1 ),       
!           ap( 2 ) and ap( 3 ) contain a( 1, 2 ) and a( 2, 2 )         
!           respectively, and so on.                                    
!           before entry with uplo = 'L' or 'L', the array ap must      
!           contain the lower triangular matrix packed sequentially,    
!           column by column, so that ap( 1 ) contains a( 1, 1 ),       
!           ap( 2 ) and ap( 3 ) contain a( 2, 1 ) and a( 3, 1 )         
!           respectively, and so on.                                    
!           note that when  diag = 'U' or 'U', the diagonal elements of 
!           a are not referenced, but are assumed to be unity.          
!           unchanged on exit.                                          
!                                                                       
!  x - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incx ) ).                            
!           before entry, the incremented array x must contain the n    
!           element right-hand side vector b. on exit, x is overwritten 
!           with the solution vector x.                                 
!                                                                       
!  incx - integer.                                                      
!           on entry, incx specifies the increment for the elements of  
!           x. incx must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!                                                                       
!  level 2 blas routine.                                                
!                                                                       
! -- written on 22-october-1986.                                        
!     jack dongarra, argonne national lab.                              
!     jeremy du croz, nag central office.                               
!     sven hammarling, nag central office.                              
!     richard hanson, sandia national labs.                             
!                                                                       
      integer            incx, n 
      character        diag, trans, uplo 
!     .. array arguments ..                                             
      real               ap( * ), x( * ) 
!     ..                                                                
!                                                                       
!     .. parameters ..                                                  
!     .. local scalars ..                                               
      real               temp 
      integer            i, info, ix, j, jx, k, kk, kx 
      logical            nounit 
!     .. external functions ..                                          
      logical            lsame 
      external           lsame 
!     .. external subroutines ..                                        
      external           xerbla 
!     ..                                                                
!     .. executable statements ..                                       
!                                                                       
!     test the input parameters.                                        
!                                                                       
      info = 0 
      if     ( .not.lsame ( uplo , 'U' ).and.                           &
     &         .not.lsame ( uplo , 'L' )     ) then                     
         info = 1 
      else if ( .not.lsame ( trans, 'N' ).and.                          &
     &         .not.lsame ( trans, 'T' ).and.                           &
     &         .not.lsame ( trans, 'C' )     ) then                     
         info = 2 
      else if ( .not.lsame ( diag , 'U' ).and.                          &
     &         .not.lsame ( diag , 'N' )     ) then                     
         info = 3 
      else if ( n < 0 ) then 
         info = 4 
      else if ( incx == 0 ) then 
         info = 7 
      end if 
      if ( info /= 0 ) then 
         call xerbla ( 'stpsv ', info ) 
         return 
      end if 
!                                                                       
!  Quick return if possible.                                            
!                                                                       
      if ( n == 0 ) then 
        return 
      end if 
                                                                        
      nounit = lsame ( diag, 'N' ) 
!                                                                       
!  Set up the start point in x if the increment is not unity. this      
!     will be  ( n - 1 ) * incx  too small for descending loops.        
!                                                                       
      if ( incx <= 0 ) then 
         kx = 1 - ( n - 1 ) * incx 
      else if ( incx /= 1 ) then 
         kx = 1 
      end if 
!                                                                       
!  Start the operations. in this version the elements of ap are         
!     accessed sequentially with one pass through ap.                   
!                                                                       
      if ( lsame ( trans, 'N' ) ) then 
!                                                                       
!  Form  x := inv( a ) * x.                                             
!                                                                       
         if ( lsame ( uplo, 'U' ) ) then 
            kk = ( n*( n + 1 ) )/2 
            if ( incx == 1 ) then 
               do j = n, 1, -1 
                  if ( x(j) /= 0.0 ) then 
                     if ( nounit ) then 
                       x(j) = x(j)/ap(kk) 
                     end if 
                     temp = x(j) 
                     k = kk - 1 
                     do i = j - 1, 1, -1 
                        x(i) = x(i) - temp * ap(k) 
                        k = k - 1 
                     end do 
                  end if 
                  kk = kk - j 
               end do 
            else 
               jx = kx + ( n - 1 ) * incx 
               do j = n, 1, -1 
                  if ( x(jx) /= 0.0 ) then 
                     if ( nounit ) then 
                       x(jx) = x(jx) / ap(kk) 
                     end if 
                     temp = x(jx) 
                     ix = jx 
                     do k = kk - 1, kk - j + 1, -1 
                        ix = ix - incx 
                        x(ix) = x(ix) - temp * ap(k) 
                     end do 
                  end if 
                  jx = jx - incx 
                  kk = kk - j 
               end do 
            end if 
         else 
            kk = 1 
            if ( incx == 1 ) then 
               do j = 1, n 
                  if ( x(j) /= 0.0 ) then 
                     if ( nounit ) then 
                       x(j) = x(j) / ap(kk) 
                     end if 
                     temp = x(j) 
                     k = kk + 1 
                     do i = j + 1, n 
                        x(i) = x(i) - temp * ap(k) 
                        k = k + 1 
                     end do 
                  end if 
                  kk = kk + ( n - j + 1 ) 
               end do 
            else 
               jx = kx 
               do j = 1, n 
                  if ( x(jx) /= 0.0 ) then 
                     if ( nounit ) then 
                       x(jx) = x(jx) / ap(kk) 
                     end if 
                     temp = x(jx) 
                     ix = jx 
                     do k = kk + 1, kk + n - j 
                        ix = ix + incx 
                        x(ix) = x(ix) - temp * ap(k) 
                     end do 
                  end if 
                  jx = jx + incx 
                  kk = kk + ( n - j + 1 ) 
               end do 
            end if 
         end if 
      else 
!                                                                       
!  Form  x := inv( a' ) * x.                                            
!                                                                       
         if ( lsame ( uplo, 'U' ) ) then 
            kk = 1 
            if ( incx == 1 ) then 
               do j = 1, n 
                  temp = x(j) 
                  k = kk 
                  do i = 1, j - 1 
                     temp = temp - ap(k) * x(i) 
                     k = k + 1 
                  end do 
                  if ( nounit ) then 
                    temp = temp / ap( kk + j - 1 ) 
                  end if 
                  x(j) = temp 
                  kk = kk + j 
               end do 
            else 
               jx = kx 
               do j = 1, n 
                  temp = x(jx) 
                  ix = kx 
                  do k = kk, kk + j - 2 
                     temp = temp - ap(k) * x(ix) 
                     ix = ix + incx 
                  end do 
                  if ( nounit ) then 
                    temp = temp/ap( kk + j - 1 ) 
                  end if 
                  x(jx) = temp 
                  jx = jx + incx 
                  kk = kk + j 
               end do 
            end if 
         else 
            kk = ( n*( n + 1 ) )/2 
            if ( incx == 1 ) then 
               do j = n, 1, -1 
                  temp = x(j) 
                  k = kk 
                  do i = n, j + 1, -1 
                     temp = temp - ap(k) * x(i) 
                     k = k - 1 
                  end do 
                  if ( nounit ) then 
                    temp = temp / ap( kk - n + j ) 
                  end if 
                  x(j) = temp 
                  kk = kk - ( n - j + 1 ) 
               end do 
            else 
               kx = kx + ( n - 1 ) * incx 
               jx = kx 
               do j = n, 1, -1 
                  temp = x(jx) 
                  ix = kx 
                  do k = kk, kk - ( n - ( j + 1 ) ), -1 
                     temp = temp - ap(k) * x(ix) 
                     ix = ix - incx 
                  end do 
                  if ( nounit ) then 
                    temp = temp / ap( kk - n + j ) 
                  end if 
                  x(jx) = temp 
                  jx = jx - incx 
                  kk = kk - ( n - j + 1 ) 
               end do 
            end if 
         end if 
      end if 
                                                                        
      return 
      END                                           
      subroutine strmv ( uplo, trans, diag, n, a, lda, x, incx ) 
!                                                                       
!***********************************************************************
!                                                                       
!! STRMV performs one of the matrix-vector operations                   
!                                                                       
!     x := a * x,   or   x := a'*x,                                     
!                                                                       
!  where x is an n element vector and  a is an n by n unit, or non-unit,
!  upper or lower triangular matrix.                                    
!                                                                       
!  Parameters:                                                          
!                                                                       
!  uplo - character.                                                    
!           on entry, uplo specifies whether the matrix is an upper or  
!           lower triangular matrix as follows:                         
!                                                                       
!              uplo = 'U' or 'U'   a is an upper triangular matrix.     
!                                                                       
!              uplo = 'L' or 'L'   a is a lower triangular matrix.      
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  trans - character.                                                   
!           on entry, trans specifies the operation to be performed as  
!           follows:                                                    
!                                                                       
!              trans = 'N' or 'N'   x := a * x.                         
!                                                                       
!              trans = 'T' or 'T'   x := a' * x.                        
!                                                                       
!              trans = 'C' or 'C'   x := a' * x.                        
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  diag - character.                                                    
!           on entry, diag specifies whether or not a is unit           
!           triangular as follows:                                      
!                                                                       
!              diag = 'U' or 'U'   a is assumed to be unit triangular.  
!                                                                       
!              diag = 'N' or 'N'   a is not assumed to be unit          
!                                  triangular.                          
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  n - integer.                                                         
!           on entry, n specifies the order of the matrix a.            
!           n must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  a - real             array of dimension ( lda, n ).                  
!           before entry with  uplo = 'U' or 'U', the leading n by n    
!           upper triangular part of the array a must contain the upper 
!           triangular matrix and the strictly lower triangular part of 
!           a is not referenced.                                        
!           before entry with uplo = 'L' or 'L', the leading n by n     
!           lower triangular part of the array a must contain the lower 
!           triangular matrix and the strictly upper triangular part of 
!           a is not referenced.                                        
!           note that when  diag = 'U' or 'U', the diagonal elements of 
!           a are not referenced either, but are assumed to be unity.   
!           unchanged on exit.                                          
!                                                                       
!  lda - integer.                                                       
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least          
!           max ( 1, n ).                                               
!           unchanged on exit.                                          
!                                                                       
!  x - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incx ) ).                            
!           before entry, the incremented array x must contain the n    
!           element vector x. on exit, x is overwritten with the        
!           tranformed vector x.                                        
!                                                                       
!  incx - integer.                                                      
!           on entry, incx specifies the increment for the elements of  
!           x. incx must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!                                                                       
!  level 2 blas routine.                                                
!                                                                       
! -- written on 22-october-1986.                                        
!     jack dongarra, argonne national lab.                              
!     jeremy du croz, nag central office.                               
!     sven hammarling, nag central office.                              
!     richard hanson, sandia national labs.                             
!                                                                       
      integer            incx, lda, n 
      character        diag, trans, uplo 
!     .. array arguments ..                                             
      real               a( lda, * ), x( * ) 
!     ..                                                                
!                                                                       
!     .. local scalars ..                                               
      real               temp 
      integer            i, info, ix, j, jx, kx 
      logical            nounit 
!     .. external functions ..                                          
      logical            lsame 
      external           lsame 
!     .. external subroutines ..                                        
      external           xerbla 
!     .. intrinsic functions ..                                         
      intrinsic          max 
!     ..                                                                
!     .. executable statements ..                                       
!                                                                       
!     test the input parameters.                                        
!                                                                       
      info = 0 
      if     ( .not.lsame ( uplo , 'U' ).and.                           &
     &         .not.lsame ( uplo , 'L' )     ) then                     
         info = 1 
      else if ( .not.lsame ( trans, 'N' ).and.                          &
     &         .not.lsame ( trans, 'T' ).and.                           &
     &         .not.lsame ( trans, 'C' )     ) then                     
         info = 2 
      else if ( .not.lsame ( diag , 'U' ).and.                          &
     &         .not.lsame ( diag , 'N' )     ) then                     
         info = 3 
      else if ( n < 0 ) then 
         info = 4 
      else if ( lda < max ( 1, n ) ) then 
         info = 6 
      else if ( incx == 0 ) then 
         info = 8 
      end if 
      if ( info /= 0 ) then 
         call xerbla ( 'strmv ', info ) 
         return 
      end if 
!                                                                       
!  Quick return if possible.                                            
!                                                                       
      if ( n == 0 ) then 
        return 
      end if 
!                                                                       
      nounit = lsame ( diag, 'N' ) 
!                                                                       
!  Set up the start point in x if the increment is not unity. this      
!     will be  ( n - 1 ) * incx  too small for descending loops.        
!                                                                       
      if ( incx <= 0 ) then 
         kx = 1 - ( n - 1 ) * incx 
      else if ( incx /= 1 ) then 
         kx = 1 
      end if 
!                                                                       
!  Start the operations. in this version the elements of a are          
!     accessed sequentially with one pass through a.                    
!                                                                       
      if ( lsame ( trans, 'N' ) ) then 
!                                                                       
!  Form  x := a * x.                                                    
!                                                                       
         if ( lsame ( uplo, 'U' ) ) then 
            if ( incx == 1 ) then 
               do j = 1, n 
                  if ( x(j) /= 0.0 ) then 
                     temp = x(j) 
                     do i = 1, j - 1 
                        x(i) = x(i) + temp * a(i,j) 
                     end do 
                     if ( nounit ) then 
                       x(j) = x(j) * a(j,j) 
                     end if 
                  end if 
               end do 
            else 
               jx = kx 
               do j = 1, n 
                  if ( x(jx) /= 0.0 ) then 
                     temp = x(jx) 
                     ix = kx 
                     do i = 1, j - 1 
                        x(ix) = x(ix) + temp * a(i,j) 
                        ix = ix + incx 
                     end do 
                     if ( nounit ) then 
                       x(jx) = x(jx) * a(j,j) 
                     end if 
                  end if 
                  jx = jx + incx 
               end do 
            end if 
         else 
            if ( incx == 1 ) then 
               do j = n, 1, -1 
                  if ( x(j) /= 0.0 ) then 
                     temp = x(j) 
                     do i = n, j + 1, -1 
                        x(i) = x(i) + temp * a(i,j) 
                     end do 
                     if ( nounit ) then 
                       x(j) = x(j) * a(j,j) 
                     end if 
                  end if 
               end do 
            else 
               kx = kx + ( n - 1 ) * incx 
               jx = kx 
               do j = n, 1, -1 
                  if ( x(jx) /= 0.0 ) then 
                     temp = x(jx) 
                     ix = kx 
                     do i = n, j + 1, -1 
                        x(ix) = x(ix) + temp * a(i,j) 
                        ix = ix - incx 
                     end do 
                     if ( nounit ) then 
                       x(jx) = x(jx) * a(j,j) 
                     end if 
                  end if 
                  jx = jx - incx 
               end do 
            end if 
         end if 
      else 
!                                                                       
!  Form  x := a' * x.                                                   
!                                                                       
         if ( lsame ( uplo, 'U' ) ) then 
            if ( incx == 1 ) then 
               do j = n, 1, -1 
                  temp = x(j) 
                  if ( nounit ) then 
                    temp = temp * a(j,j) 
                  end if 
                  do i = j - 1, 1, -1 
                     temp = temp + a(i,j) * x(i) 
                  end do 
                  x(j) = temp 
               end do 
            else 
               jx = kx + ( n - 1 ) * incx 
               do j = n, 1, -1 
                  temp = x(jx) 
                  ix = jx 
                  if ( nounit ) then 
                    temp = temp * a(j,j) 
                  end if 
                  do i = j - 1, 1, -1 
                     ix = ix - incx 
                     temp = temp + a(i,j) * x(ix) 
                  end do 
                  x(jx) = temp 
                  jx = jx - incx 
               end do 
            end if 
         else 
            if ( incx == 1 ) then 
               do j = 1, n 
                  temp = x(j) 
                  if ( nounit ) then 
                    temp = temp * a(j,j) 
                  end if 
                  do i = j + 1, n 
                     temp = temp + a(i,j) * x(i) 
                  end do 
                  x(j) = temp 
               end do 
            else 
               jx = kx 
               do j = 1, n 
                  temp = x(jx) 
                  ix = jx 
                  if ( nounit ) then 
                    temp = temp * a(j,j) 
                  end if 
                  do i = j + 1, n 
                     ix = ix + incx 
                     temp = temp + a(i,j) * x(ix) 
                  end do 
                  x(jx) = temp 
                  jx = jx + incx 
               end do 
            end if 
         end if 
      end if 
                                                                        
      return 
      END                                           
      subroutine strsv ( uplo, trans, diag, n, a, lda, x, incx ) 
!                                                                       
!***********************************************************************
!                                                                       
!! STRSV solves one of the systems of equations                         
!                                                                       
!     a * x = b,   or   a'*x = b,                                       
!                                                                       
!  where b and x are n element vectors and a is an n by n unit, or      
!  non-unit, upper or lower triangular matrix.                          
!                                                                       
!  no test for singularity or near-singularity is included in this      
!  routine. such tests must be performed before calling this routine.   
!                                                                       
!  Parameters:                                                          
!                                                                       
!  uplo - character.                                                    
!           on entry, uplo specifies whether the matrix is an upper or  
!           lower triangular matrix as follows:                         
!                                                                       
!              uplo = 'U' or 'U'   a is an upper triangular matrix.     
!                                                                       
!              uplo = 'L' or 'L'   a is a lower triangular matrix.      
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  trans - character.                                                   
!           on entry, trans specifies the equations to be solved as     
!           follows:                                                    
!                                                                       
!              trans = 'N' or 'N'   a * x = b.                          
!                                                                       
!              trans = 'T' or 'T'   a' * x = b.                         
!                                                                       
!              trans = 'C' or 'C'   a' * x = b.                         
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  diag - character.                                                    
!           on entry, diag specifies whether or not a is unit           
!           triangular as follows:                                      
!                                                                       
!              diag = 'U' or 'U'   a is assumed to be unit triangular.  
!                                                                       
!              diag = 'N' or 'N'   a is not assumed to be unit          
!                                  triangular.                          
!                                                                       
!           unchanged on exit.                                          
!                                                                       
!  n - integer.                                                         
!           on entry, n specifies the order of the matrix a.            
!           n must be at least 0.                                       
!           unchanged on exit.                                          
!                                                                       
!  a - real             array of dimension ( lda, n ).                  
!           before entry with  uplo = 'U' or 'U', the leading n by n    
!           upper triangular part of the array a must contain the upper 
!           triangular matrix and the strictly lower triangular part of 
!           a is not referenced.                                        
!           before entry with uplo = 'L' or 'L', the leading n by n     
!           lower triangular part of the array a must contain the lower 
!           triangular matrix and the strictly upper triangular part of 
!           a is not referenced.                                        
!           note that when  diag = 'U' or 'U', the diagonal elements of 
!           a are not referenced either, but are assumed to be unity.   
!           unchanged on exit.                                          
!                                                                       
!  lda - integer.                                                       
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least          
!           max ( 1, n ).                                               
!           unchanged on exit.                                          
!                                                                       
!  x - real             array of dimension at least                     
!           ( 1 + ( n - 1 ) * abs( incx ) ).                            
!           before entry, the incremented array x must contain the n    
!           element right-hand side vector b. on exit, x is overwritten 
!           with the solution vector x.                                 
!                                                                       
!  incx - integer.                                                      
!           on entry, incx specifies the increment for the elements of  
!           x. incx must not be 0.                                      
!           unchanged on exit.                                          
!                                                                       
!                                                                       
!  level 2 blas routine.                                                
!                                                                       
! -- written on 22-october-1986.                                        
!     jack dongarra, argonne national lab.                              
!     jeremy du croz, nag central office.                               
!     sven hammarling, nag central office.                              
!     richard hanson, sandia national labs.                             
!                                                                       
      integer            incx, lda, n 
      character        diag, trans, uplo 
!     .. array arguments ..                                             
      real               a( lda, * ), x( * ) 
!     ..                                                                
!     .. local scalars ..                                               
      real               temp 
      integer            i, info, ix, j, jx, kx 
      logical            nounit 
!     .. external functions ..                                          
      logical            lsame 
      external           lsame 
!     .. external subroutines ..                                        
      external           xerbla 
!     .. intrinsic functions ..                                         
      intrinsic          max 
!     ..                                                                
!     .. executable statements ..                                       
!                                                                       
!     test the input parameters.                                        
!                                                                       
      info = 0 
      if     ( .not.lsame ( uplo , 'U' ).and.                           &
     &         .not.lsame ( uplo , 'L' )     ) then                     
         info = 1 
      else if ( .not.lsame ( trans, 'N' ).and.                          &
     &         .not.lsame ( trans, 'T' ).and.                           &
     &         .not.lsame ( trans, 'C' )     ) then                     
         info = 2 
      else if ( .not.lsame ( diag , 'U' ).and.                          &
     &         .not.lsame ( diag , 'N' )     ) then                     
         info = 3 
      else if ( n < 0 ) then 
         info = 4 
      else if ( lda < max ( 1, n ) ) then 
         info = 6 
      else if ( incx == 0 ) then 
         info = 8 
      end if 
      if ( info /= 0 ) then 
         call xerbla ( 'strsv ', info ) 
         return 
      end if 
!                                                                       
!  Quick return if possible.                                            
!                                                                       
      if ( n == 0 ) then 
        return 
      end if 
!                                                                       
      nounit = lsame ( diag, 'N' ) 
!                                                                       
!  Set up the start point in x if the increment is not unity. this      
!     will be  ( n - 1 ) * incx  too small for descending loops.        
!                                                                       
      if ( incx <= 0 ) then 
         kx = 1 - ( n - 1 ) * incx 
      else if ( incx /= 1 ) then 
         kx = 1 
      end if 
!                                                                       
!  Start the operations. in this version the elements of a are          
!     accessed sequentially with one pass through a.                    
!                                                                       
      if ( lsame ( trans, 'N' ) ) then 
!                                                                       
!  Form  x := inv( a ) * x.                                             
!                                                                       
         if ( lsame ( uplo, 'U' ) ) then 
            if ( incx == 1 ) then 
               do j = n, 1, -1 
                  if ( x(j) /= 0.0 ) then 
                     if ( nounit ) then 
                       x(j) = x(j)/a(j,j) 
                     end if 
                     temp = x(j) 
                     do i = j - 1, 1, -1 
                        x(i) = x(i) - temp * a(i,j) 
                     end do 
                  end if 
               end do 
            else 
               jx = kx + ( n - 1 ) * incx 
               do j = n, 1, -1 
                  if ( x(jx) /= 0.0 ) then 
                     if ( nounit ) then 
                       x(jx) = x(jx)/a(j,j) 
                     end if 
                     temp = x(jx) 
                     ix = jx 
                     do i = j - 1, 1, -1 
                        ix = ix - incx 
                        x(ix) = x(ix) - temp * a(i,j) 
                     end do 
                  end if 
                  jx = jx - incx 
               end do 
            end if 
         else 
            if ( incx == 1 ) then 
               do j = 1, n 
                  if ( x(j) /= 0.0 ) then 
                     if ( nounit ) then 
                       x(j) = x(j)/a(j,j) 
                     end if 
                     temp = x(j) 
                     do i = j + 1, n 
                        x(i) = x(i) - temp * a(i,j) 
                     end do 
                  end if 
               end do 
            else 
               jx = kx 
               do j = 1, n 
                  if ( x(jx) /= 0.0 ) then 
                     if ( nounit ) then 
                       x(jx) = x(jx)/a(j,j) 
                     end if 
                     temp = x(jx) 
                     ix = jx 
                     do i = j + 1, n 
                        ix = ix + incx 
                        x(ix) = x(ix) - temp * a(i,j) 
                     end do 
                  end if 
                  jx = jx + incx 
               end do 
            end if 
         end if 
      else 
!                                                                       
!  Form  x := inv( a' ) * x.                                            
!                                                                       
         if ( lsame ( uplo, 'U' ) ) then 
            if ( incx == 1 ) then 
               do j = 1, n 
                  temp = x(j) 
                  do i = 1, j - 1 
                     temp = temp - a(i,j) * x(i) 
                  end do 
                  if ( nounit ) then 
                    temp = temp/a(j,j) 
                  end if 
                  x(j) = temp 
               end do 
            else 
               jx = kx 
               do j = 1, n 
                  temp = x(jx) 
                  ix = kx 
                  do i = 1, j - 1 
                     temp = temp - a(i,j) * x(ix) 
                     ix = ix + incx 
                  end do 
                  if ( nounit ) then 
                    temp = temp/a(j,j) 
                  end if 
                  x(jx) = temp 
                  jx = jx + incx 
               end do 
            end if 
         else 
            if ( incx == 1 ) then 
               do j = n, 1, -1 
                  temp = x(j) 
                  do i = n, j + 1, -1 
                     temp = temp - a(i,j) * x(i) 
                  end do 
                  if ( nounit ) then 
                    temp = temp/a(j,j) 
                  end if 
                  x(j) = temp 
               end do 
            else 
               kx = kx + ( n - 1 ) * incx 
               jx = kx 
               do j = n, 1, -1 
                  temp = x(jx) 
                  ix = kx 
                  do i = n, j + 1, -1 
                     temp = temp - a(i,j) * x(ix) 
                     ix = ix - incx 
                  end do 
                  if ( nounit ) then 
                    temp = temp/a(j,j) 
                  end if 
                  x(jx) = temp 
                  jx = jx - incx 
               end do 
            end if 
         end if 
      end if 
                                                                        
      return 
      END                                           
      subroutine xerbla ( srname, info ) 
!                                                                       
!***********************************************************************
!                                                                       
!! XERBLA is an error handler for the level 2 blas routines.            
!                                                                       
!  it is called by the level 2 blas routines if an input parameter is   
!  invalid.                                                             
!                                                                       
!  installers should consider modifying the stop statement in order to  
!  call system-specific exception-handling facilities.                  
!                                                                       
!  Parameters:                                                          
!                                                                       
!  srname - character*6.                                                
!           on entry, srname specifies the name of the routine which    
!           called xerbla.                                              
!                                                                       
!  info - integer.                                                      
!           on entry, info specifies the position of the invalid        
!           parameter in the parameter-list of the calling routine.     
!                                                                       
!                                                                       
!  auxiliary routine for level 2 blas.                                  
!                                                                       
!  written on 20-july-1986.                                             
!                                                                       
      integer            info 
      character*6        srname 
!                                                                       
      write (*,99999) srname, info 
                                                                        
      stop 
!                                                                       
99999 format ( ' ** on entry to ', a6, ' parameter number ', i2,        &
     &         ' had an illegal value' )                                
                                                                        
      END                                           
      function lsame ( ca, cb ) 
!                                                                       
!***********************************************************************
!                                                                       
!! LSAME tests if ca is the same letter as cb regardless of case.       
!  cb is assumed to be an upper case letter. lsame returns .true. if    
!  ca is either the same as cb or the equivalent lower case letter.     
!                                                                       
!  n.b. this version of the routine is only correct for ascii code.     
!       installers must modify the routine for other character-codes.   
!                                                                       
!       for ebcdic systems the constant ioff must be changed to -64.    
!       for cdc systems using 6-12 bit representations, the system-     
!       specific code in comments must be activated.                    
!                                                                       
!  Parameters:                                                          
!                                                                       
!  ca - character                                                       
!  cb - character                                                       
!           on entry, ca and cb specify characters to be compared.      
!           unchanged on exit.                                          
!                                                                       
!                                                                       
!  auxiliary routine for level 2 blas.                                  
!                                                                       
! -- written on 20-july-1986                                            
!     richard hanson, sandia national labs.                             
!     jeremy du croz, nag central office.                               
!                                                                       
      character            ca, cb 
      logical lsame 
      integer                ioff 
      parameter            ( ioff=32 ) 
!                                                                       
!  Test for equality.                                                   
!                                                                       
      lsame = ca == cb 
!                                                                       
!  Now test for equivalence                                             
!                                                                       
      if ( .not.lsame ) then 
         lsame = ichar(ca) - ioff == ichar(cb) 
      end if 
                                                                        
      return 
      END                                           
