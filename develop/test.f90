program test_mltbc_weights
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none

  integer, parameter :: nstep = 8
  real(real64)       :: weights(nstep), cum(nstep)
  integer            :: i
  character(len=32)  :: method

  interface
    subroutine mltbc_compute_weights(nstep, method, weights)
      use, intrinsic :: iso_fortran_env, only: real64
      integer, intent(in)                :: nstep
      character(len=*), intent(in)       :: method
      real(real64), intent(inout)        :: weights(nstep)
    end subroutine mltbc_compute_weights
  end interface

  method = 'DolphChebyshev'

  ! Call the subroutine
  call mltbc_compute_weights(nstep, method, weights)

  ! Print weights
  print *, '=== Weights (method: ', trim(method), ') ==='
  do i = 1, nstep
     print '(i3,2x,f12.6)', i, weights(i)
  end do

  ! Check sum
  print *, 'Sum of weights = ', sum(weights)

  ! Check symmetry
  print *, 'Symmetry check (weight(i) vs weight(n-i+1)):'
  do i = 1, nstep / 2
     print '(i3,2x,f8.5,2x,f8.5)', i, weights(i), weights(nstep - i + 1)
  end do

  ! Cumulative sum
  cum(1) = weights(1)
  do i = 2, nstep
     cum(i) = cum(i - 1) + weights(i)
  end do
  print *, 'Cumulative application:'
  do i = 1, nstep
     print '(i3,2x,f8.5)', i, cum(i)
  end do

end program test_mltbc_weights


