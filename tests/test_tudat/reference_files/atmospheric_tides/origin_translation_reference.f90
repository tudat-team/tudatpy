! Independent reference calculation of published IERS2010 Eq20/Table7.6.
! Newly written optional test generator; no upstream software is copied.
! Complex Fourier synthesis and Cartesian cross-product basis, independent of the production C++ trigonometric path.
program atmospheric_origin_reference
  use iso_fortran_env, only: real64
  implicit none
  real(real64) :: c(3,4), epochs(7), pi, theta, radial(3), east(3), north(3), cart(3), local(3)
  complex(real64) :: wave1, wave2
  integer :: i,j
  c(1,:) = [2.1188e-4_real64,-7.6861e-4_real64,1.4472e-4_real64,-1.7844e-4_real64]
  c(2,:) = [-7.2766e-4_real64,-2.3582e-4_real64,-3.2691e-4_real64,-1.5878e-4_real64]
  c(3,:) = [-1.2176e-5_real64,3.2243e-5_real64,-9.6271e-5_real64,1.6976e-5_real64]
  epochs = [-43200._real64,-21600._real64,0._real64,21600._real64,12345.6789_real64,-78583.25_real64,203726.875_real64]
  pi=acos(-1._real64)
  ! Explicit Cartesian test position; it establishes no provider grid latitude convention.
  radial = [4500000._real64,-3000000._real64,2800000._real64]
  radial = radial/sqrt(sum(radial**2))
  east = [-radial(2),radial(1),0._real64]
  east = east/sqrt(sum(east**2))
  north(1)=radial(2)*east(3)-radial(3)*east(2)
  north(2)=radial(3)*east(1)-radial(1)*east(3)
  north(3)=radial(1)*east(2)-radial(2)*east(1)
  do i=1,size(epochs)
    theta=2*pi*(epochs(i)/86400._real64+0.5_real64)
    wave1=exp(cmplx(0._real64,theta,kind=real64))
    wave2=wave1**2
    do j=1,3
      cart(j)=real(cmplx(c(j,1),-c(j,2),kind=real64)*wave1 + &
                   cmplx(c(j,3),-c(j,4),kind=real64)*wave2,kind=real64)
    enddo
    local = [dot_product(radial,cart),dot_product(north,cart),dot_product(east,cart)]
    write(*,'(7ES25.16)') epochs(i),cart,local
  enddo
end program atmospheric_origin_reference
