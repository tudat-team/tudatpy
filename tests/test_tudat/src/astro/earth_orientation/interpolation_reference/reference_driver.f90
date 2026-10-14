! Original Tudat driver; external reference routines remain unmodified.
program lagrange_reference_driver
  implicit none
  integer, parameter :: n=8
  integer :: i, j, count, year, month, day, status
  real(8) :: mjd(n), atomic_offset(n), fraction, input(n,6), x(n), y(n), q, qdat, elapsed, result(6)
  do i=1,n
     read(*,*) year, month, day, fraction, mjd(i), input(i,:)
     call IAU_DAT(year,month,day,fraction,atomic_offset(i),status)
     if(status.lt.0) stop 1
  enddo
  x=(mjd-mjd(1))*86400d0+atomic_offset-atomic_offset(1)
  input(:,3)=input(:,3)-atomic_offset
  write(*,'(ES26.17E3)') atomic_offset(1)
  read(*,*) count
  do i=1,count
     read(*,*) year,month,day,fraction,q
     call IAU_DAT(year,month,day,fraction,qdat,status)
     if(status.lt.0) stop 2
     elapsed=(q-mjd(1))*86400d0+qdat-atomic_offset(1)
     do j=1,6
        y=input(:,j)
        call LAGINT(x,y,n,elapsed,result(j))
     enddo
     write(*,'(8ES26.17E3)') q,elapsed,result
  enddo
end program lagrange_reference_driver
