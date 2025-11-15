      subroutine call_mcd4idl(zkey,xz,xlon,xlat,hireskey,datekey,xdate,localtime,dset,scena, &
            perturkey,seedin,gwlength,extvarkeys, pres,dens,temp,zonwind,merwind,            &
            meanvar,extvar,seedout,ier)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use MCD, only: call_mcd
      implicit none
! -------------------------------------------
      integer,parameter :: nextvar=101  ! size of extra variable array extvar()
      integer,parameter :: nmeanvar=6   ! size of meanvar() array
!
!     inputs
!     ******      
      integer,intent(in) :: zkey      ! flag for vertical coordinate type
      real   ,intent(in) :: xz        ! vertical coordinate (m or Pa, depending on zkey)
      real   ,intent(in) :: xlon      ! east longitude (degrees) 
      real   ,intent(in) :: xlat      ! latitude (degrees)
      integer,intent(in) :: hireskey  ! flag: high resolution if =1
      integer,intent(in) :: datekey   ! flag: 0=earth date 1= Mars date
      real*8 ,intent(in) :: xdate     ! earth julian date (or Ls)
      real   ,intent(in) :: localtime ! true solar time at longitude lon
      integer,intent(in) :: perturkey ! flag: perturbation type
      real   ,intent(in) :: seedin    ! seed for random number generation
      real   ,intent(in) :: gwlength  ! gravity wave perturbation (vertical wavelength of)
      integer,intent(in) :: scena     ! dust scenario
      integer,intent(in) :: extvarkeys(nextvar) ! flag: compute extra variable i
                                      !       if i==1
      character*(*),intent(in) :: dset      ! path to datafiles
!
!     outputs
!     *******      
      real   ,intent(out) :: meanvar(nmeanvar) ! array for mean values
      real   ,intent(out) :: extvar(nextvar)   ! array for extra outputs
      real   ,intent(out) :: seedout ! current value of 'seed' (used by ran1)
      real   ,intent(out) :: pres    ! atmospheric pressure
      real   ,intent(out) :: dens    ! atmospheric density
      real   ,intent(out) :: temp    ! atmospheric temperature
      real   ,intent(out) :: zonwind ! zonal (eastward) wind
      real   ,intent(out) :: merwind ! meridional (northward) wind
      integer,intent(out) :: ier     ! status (error) code
!
!     local variables
!     *************** 
!
      integer,parameter :: nextvar1=100 ! size of extra variable array extvar()
      integer,parameter :: nmeanvar1=5  ! size of meanvar() array
      integer extvarkeys1(nextvar1)     ! flag: compute extra variable i
      real    meanvar1(nmeanvar1)       ! array for mean values
      real    extvar1(nextvar1)         ! array for extra outputs
      integer,parameter :: str_max=50   ! max size of dset
      character(len=str_max) :: dset1
      integer i,n

      do i=1,nextvar1
       extvarkeys1(i) = extvarkeys(i+1)
      enddo

      n  = ichar(dset(1:1))

      do i=1,str_max
       if(i.le.n) then
        dset1(i:i)=dset(i+1:i+1)
       else
        dset1(i:i)= ' '
       endif 
      enddo

      call call_mcd(zkey,xz,xlon,xlat,hireskey, datekey,xdate,localtime,dset1,scena, &
            perturkey,seedin,gwlength,extvarkeys1,pres,dens,temp,zonwind,merwind,    &
            meanvar1,extvar1,seedout,ier)

      do i=1,nextvar1
       extvar(i+1)     = extvar1(i)
      enddo
      do i=1,nmeanvar1
       meanvar(i+1) = meanvar1(i)
      enddo

      end ! end of call_mcd

