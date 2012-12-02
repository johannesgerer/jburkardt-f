!  Here are some floating comments.
!
program alpha

!*****************************************************************************80
!
!! ALPHA is a main program.
!
  call beta ( y )
  x = gamma ( y )

  stop
end
subroutine beta ( y, delta, gamma )

!*****************************************************************************80
!
!! BETA is a subroutine
!
!  real function zeta ( x ) is commented out.
!
  real delta
  real gamma
  character ( len = 50 ) s
  real y

  gamma = 7.0
  y = 3.14159265
!
!  This line might confuse the extractor!
!
  s = 'function gamma ( 17 )'

  return
end


function gamma ( alpha, beta, delta )

!*****************************************************************************80
!
!! GAMMA is a function with no type statement.
!
!  Here's another commented out line for suckers.
!
!     integer function epsilon ( x )
!
  real alpha
  real beta
  real delta
  real gamma
  real :: x = 17.0

  gamma = sqrt ( x )

  return
end
complex function delta ( x )

!*****************************************************************************80
!
!! DELTA is a complex function with type statement.
!
  complex x

  delta = x * x

  return
end

  write ( *, * ) ' '
  write ( *, * ) 'This is a main program with no PROGRAM statement.'
  stop
end
integer function epsilon ( x )

!*****************************************************************************80
!
!! EPSILON is an integer function with type statement.
!
  integer x

  epsilon = x + 1

  return
end
real function zeta ( x )

!*****************************************************************************80
!
!! ZETA is a real function with type statement.
!
  real x

  zeta = 1.0 / x

  return
end
module eta

!*****************************************************************************80
!
!! ETA is a module.
!
  real :: x = 17.0

end
function beta ( y, delta, gamma )

!*****************************************************************************80
!
!! BETA is a function, but this is the second module with this name!
!
  beta = 17.0

  return
end
block data

!*****************************************************************************80
!
!! This is a blank block data routine.
!
  real :: x = 1.0
  common / morley / x

end
block data theta

!*****************************************************************************80
!
!! THETA is a block data routine.
!
  real x
  common / samantha / x
  save / samantha /
  data x / 1.0 /

end
recursive function iota ( x ) result ( value )

!*****************************************************************************80
!
!! IOTA is a recursive function.
!
  if ( x <= 0.0 ) then
    value = 0.0
  else
    value = x + iota ( x - 1 )
  end if

  return
end
block data beta

!*****************************************************************************80
!
!! BETA is a blockdata routine, but this is the third module with this name.
!
  real x
  common / smith / x
  save / smith /
  data x / 1.0 /

end
blockdata enid

!*****************************************************************************80
!
!! ENID is a blockdata routine.
!
  real x
  common / smith / x
  save / smith /
  data x / 1.0 /

end
