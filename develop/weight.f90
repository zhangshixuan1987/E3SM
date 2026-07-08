subroutine mltbc_compute_weights(nstep, method, weights)
  !
  ! mltbc_compute_weights:
  !    Construct weighting function for ML-predicted tendency application.
  !
  ! Input:
  !   - nstep : number of steps to apply nudging
  !   - method: type of weighting method
  !
  ! Output:
  !   - weights(nstep): computed weight array (normalized)
  !
  !===============================================================
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none

  integer, intent(in)                :: nstep
  character(len=*), intent(in)       :: method
  real(real64), intent(inout)            :: weights(nstep)

  ! Local variables
  integer :: m, i, n, center, start
  real(real64) :: pi, norm
  real(real64), allocatable :: xwgt(:)

  pi = 4.0_real64 * atan(1.0_real64)
  weights(:) = 0.0_real64

  select case (trim(method))

    case ('STEP')
       weights(:) = 1.0_real64

    case ('IMT')
       weights(1) = 1.0_real64

    case ('Linear')
       weights(:) = 1.0_real64 / real(nstep, real64)

    case ('TopHat')
       m = int(nstep / 2)
       norm = 0.0_real64
       do i = 1, nstep
          if (i <= m) then
             weights(i) = real(i, real64)
          else
             weights(i) = real(nstep - i + 1, real64)
          end if
          norm = norm + weights(i)
       end do
       if (norm > 0.0_real64) weights(:) = weights(:) / norm

    case ('DolphChebyshev', 'DolphChebyshev1')
       m = (nstep + 1) / 2
       allocate(xwgt(2 * m + 1))
       do i = 1, 2 * m + 1
          n = i - m - 1
          if (n == 0) then
             xwgt(i) = 1.0_real64 / real(m, real64)
          else if (abs(n) > m) then
             xwgt(i) = 0.0_real64
          else
             xwgt(i) = sin(n * pi / real(m + 1, real64)) * (m + 1) / (n * pi) * &
                       sin(n * pi / real(m, real64)) / (n * pi)
          end if
       end do
       norm = sum(xwgt)
       ! Extract centered nstep weights
       center = (2 * m + 1 + 1) / 2          ! central index
       start = center - (nstep - 1) / 2      ! ensure symmetric center
       weights(:) = xwgt(start : start + nstep - 1) / norm
       deallocate(xwgt)
  end select

  ! Optional rescaling for special method variants
  if (trim(method) == 'IMT1' .or. trim(method) == 'DolphChebyshev1') then
     weights(:) = weights(:) * nstep
  end if

end subroutine mltbc_compute_weights
