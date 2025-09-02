program fftw_verification_test
  !
  ! A simple program to verify the behavior of FFTW3 forward and backward transforms.
  ! We define a plane wave in real space and check if its Fourier transform
  ! results in a single peak in reciprocal space, as predicted by theory.
  !
  use, intrinsic :: iso_fortran_env
  implicit none

  include 'fftw3.f'

  ! Define double precision
  integer, parameter :: dp = REAL64

  ! Grid dimensions
  integer, parameter :: nx = 8
  integer, parameter :: ny = 8
  integer, parameter :: nz = 8
  integer, parameter :: N_total = nx * ny * nz

  ! Test plane wave frequency indices (qx, qy, qz)
  ! Should be within [0, nx-1], [0, ny-1], [0, nz-1]
  integer, parameter :: qx = 2
  integer, parameter :: qy = 3
  integer, parameter :: qz = 4

  ! Arrays for the function and its transform
  complex(dp), allocatable :: test_function(:,:,:)
  complex(dp), allocatable :: ft_result(:,:,:)
  complex(dp), allocatable :: inverse_ft_result(:,:,:)

  ! FFTW specific variables
  complex(dp), allocatable :: fftw_in(:,:,:), fftw_out(:,:,:)
  integer(8) :: plan_f, plan_b

  ! Loop counters and constants
  integer :: i, j, k
  real(dp), parameter :: pi = 3.141592653589793238_dp
  complex(dp), parameter :: iimag = cmplx(0.0_dp, 1.0_dp, kind=dp)

  ! Allocate arrays
  allocate(test_function(nx, ny, nz))
  allocate(ft_result(nx, ny, nz))
  allocate(inverse_ft_result(nx, ny, nz))
  allocate(fftw_in(nx, ny, nz), fftw_out(nx, ny, nz))

  ! 1. === DEFINE THE TEST FUNCTION (a plane wave) ===
  ! Fortran arrays are 1-based, so we use (i-1), (j-1), (k-1) to match
  ! the theoretical formula which is 0-based.
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        test_function(i,j,k) = exp( 2.0_dp * pi * iimag * &
             ( real(qx, dp) * real(i-1, dp) / real(nx, dp) + &
               real(qy, dp) * real(j-1, dp) / real(ny, dp) + &
               real(qz, dp) * real(k-1, dp) / real(nz, dp) ) )
      end do
    end do
  end do

  ! 2. === SETUP FFTW PLANS ===
  ! Create plans for forward and backward 3D DFTs.
  call dfftw_plan_dft_3d(plan_f, nx, ny, nz, fftw_in, fftw_out, FFTW_FORWARD, FFTW_MEASURE)
  call dfftw_plan_dft_3d(plan_b, nx, ny, nz, fftw_in, fftw_out, FFTW_BACKWARD, FFTW_MEASURE)

  ! 3. === EXECUTE FORWARD TRANSFORM ===
  fftw_in = test_function
  call dfftw_execute(plan_f)
  ft_result = fftw_out

  ! 4. === VERIFY THE FORWARD TRANSFORM RESULT ===
  write(*,*) '--- Verification of Forward Fourier Transform ---'
  write(*,*) 'Grid size: ', nx, 'x', ny, 'x', nz, ' = ', N_total, ' points.'
  write(*,*) 'Input plane wave frequency indices (qx,qy,qz): ', qx, qy, qz
  write(*,*) 'Theoretical peak value should be close to: ', real(N_total, dp)
  write(*,*)

  ! The peak should be at (qx+1, qy+1, qz+1) due to 1-based indexing
  write(*,'(A,3I4,A,2F12.6)') 'Peak found at index (', qx+1, qy+1, qz+1, &
       ') with value: ', ft_result(qx+1, qy+1, qz+1)

  ! Check a few other points, which should be close to zero.
  write(*,'(A,3I4,A,2F12.6)') 'Value at index      (', 1, 1, 1, &
       ') is:           ', ft_result(1, 1, 1)
  write(*,'(A,3I4,A,2F12.6)') 'Value at index      (', 2, 2, 2, &
       ') is:           ', ft_result(2, 2, 2)
  write(*,*)

  ! 5. === EXECUTE BACKWARD TRANSFORM to recover original function ===
  fftw_in = ft_result
  call dfftw_execute(plan_b)
  ! IMPORTANT: Normalize the result of the backward transform
  inverse_ft_result = fftw_out / real(N_total, dp)

  ! 6. === VERIFY THE BACKWARD TRANSFORM RESULT ===
  write(*,*) '--- Verification of Backward Fourier Transform ---'
  write(*,*) 'Comparing original function with the back-transformed one.'
  write(*,*) 'The difference should be close to zero.'

  ! Check the value at a sample point (e.g., i=5, j=5, k=5)
  i = 4; j = 6; k = 5;
  write(*,'(A,2F12.6)') 'Original value at (5,5,5):      ', test_function(i,j,k)
  write(*,'(A,2F12.6)') 'Back-transformed value at (5,5,5):', inverse_ft_result(i,j,k)
  write(*,'(A,E12.4)')  'Absolute difference:              ', abs(test_function(i,j,k) - inverse_ft_result(i,j,k))
  write(*,*)

  ! 7. === CLEANUP ===
  call dfftw_destroy_plan(plan_f)
  call dfftw_destroy_plan(plan_b)

  deallocate(test_function, ft_result, inverse_ft_result, fftw_in, fftw_out)

end program fftw_verification_test
