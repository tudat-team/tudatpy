! Original Tudat test driver, neither distributed nor endorsed by IERS.
! Reference routines are compiled unmodified in an external cache.
program reference_driver
  implicit none
  integer :: n,i,year,month,day,status
  double precision :: mjd,t,pm(2),dut1,dlod,orthogonal(3),ocean(4),args(6),gmst,delta_at,tt_seconds
  double precision, parameter :: pi=3.141592653589793238462643d0
  read(*,*) n
  do i=1,n
    read(*,*) mjd,year,month,day
    call iau_dat(year,month,day,mjd-floor(mjd),delta_at,status)
    if (status < 0) stop 1
    tt_seconds=(mjd-51544.5d0)*86400d0+delta_at+32.184d0
    call pmsdnut2(mjd,pm)
    call utlibr(mjd,dut1,dlod)
    call ortho_eop(mjd,orthogonal)
    call pmut1_oceans(mjd,ocean(1),ocean(2),ocean(3),ocean(4))
    t=(mjd-51544.5d0)/36525d0
    call fundarg(t,args(2),args(3),args(4),args(5),args(6))
    gmst=mod(67310.54841d0+t*((8640184.812866d0+3155760000d0)+t*(0.093104d0+t*(-0.0000062d0))),86400d0)
    args(1)=mod(gmst/(86400d0/(2*pi))+pi,2*pi)
    write(*,'(18(es26.17e3,1x))') tt_seconds,pm,dut1,dlod,orthogonal,ocean,args
  end do
end program
