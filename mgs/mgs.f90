subroutine mgs ( m, n, c, r, q )

  integer m
  integer n

  real c(m,n)
  real q(m,n)
  real r(n,n)
  real x(m)
  real xn

  do k = 1, n
  do j = 1, m
  x(j) = c(j,k)
  end do
  xn = 0
  do j = 1, m
  xn = xn + x(j)*x(j)
  end do
  r(k,k) = sqrt(xn)
  if ( 0.0 < r(k,k) ) then
  do j = 1, m
  q(j,k) = c(j,k)/r(k,k)
  end do
  else
  do j = 1, m
  q(j,k) = 0.0
  end do
  end if
  do j = k+1, n
  r(k,j) = 0
  do i=1, m
  r(k,j) = r(k,j) + q(i,k)*c(i,j)
  end do
  do i=1,m
  c(i,j) = c(i,j) - q(i,k)*r(k,j)
  end do
  end do
  end do

  return
end
