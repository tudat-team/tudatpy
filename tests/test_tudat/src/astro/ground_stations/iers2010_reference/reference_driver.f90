! Original Tudat fixture driver; not distributed or endorsed by the IERS Center.
! Calls the unmodified IERS reference routines; no Tudat implementation is linked.
program reference_driver
  implicit none
  integer :: year, month, day, status, count, i
  double precision :: hour, station(3), sun(3), moon(3), displacement(3)
  double precision :: mjd0, mjd, delta_at, utc_seconds, tt_seconds, step2diurnal(3), step2long(3), centuries
  read (*,*) count
  do i = 1, count
    read (*,*) year, month, day, hour, station, sun, moon
    call dehanttideinel(station, year, month, day, hour, sun, moon, displacement)
    call iau_cal2jd(year, month, day, mjd0, mjd, status)
    if (status /= 0) stop 1
    call iau_dat(year, month, day, hour/24d0, delta_at, status)
    if (status /= 0) stop 2
    utc_seconds = ((mjd0-2451545d0)+mjd)*86400d0+hour*3600d0
    tt_seconds = utc_seconds+delta_at+32.184d0
    centuries = tt_seconds/(86400d0*36525d0)
    call step2diu(station, hour, centuries, step2diurnal)
    call step2lon(station, centuries, step2long)
    write (*,'(11(es26.17e3,1x))') utc_seconds, tt_seconds, displacement, step2diurnal, step2long
  end do
end program

! Link adapters for the IERS main routine's DAT/CAL2JD names and the
! unmodified SOFA-supplied iau_DAT/iau_CAL2JD definitions. No algorithm changes.
subroutine dat(year, month, day, fraction, delta_at, status)
  implicit none
  integer :: year, month, day, status
  double precision :: fraction, delta_at
  call iau_dat(year, month, day, fraction, delta_at, status)
end subroutine
subroutine cal2jd(year, month, day, mjd0, mjd, status)
  implicit none
  integer :: year, month, day, status
  double precision :: mjd0, mjd
  call iau_cal2jd(year, month, day, mjd0, mjd, status)
end subroutine
