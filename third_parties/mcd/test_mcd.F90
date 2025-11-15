      program test_mcd

!     ---------------------------------------------------------
!     This is a simple program used to test and give an example
!     of the call to the main mcd subroutine "call_mcd".
!     ---------------------------------------------------------

      use MCD
      

      implicit none



!ccccccccccccc CALL_MCD arguments ccccccccccccccccccccccccccccccccccccccccccccc
     
!     Inputs:

      integer          zkey            ! flag to choose the type of z coordinates
      real             z               ! value of the z coordinate (m or Pa)
      real             lon             ! east longitude (degrees)
      real             lat             ! north latitude (degrees)
      integer          hireskey        ! high resolution flag (0: off, 1: on) 
      integer          datekey         ! date flag (0: Earth date 1: Mars date)
      double precision date            ! Julian date (if datekey=0) or
                                       ! solar longitude Ls (if datekey=1) (degrees)
      real             localtime       ! local time at longitude lon (only if datekey=1) (Martian hours)
      character(len=100) :: dset=" "   ! path to MCD datasets; unset here (ie: defaults to "MCD_DATA/")
      integer          scena           ! dust and solar EUV scenario
      integer          perturkey       ! perturbation type
      real             seedin          ! random generator seed and flag (if perturkey=1,2,3 or 4)
                                       ! coefficient to multiply std. dev. by (if perturkey=5)
      real             gwlength        ! Gravity wave wavelength (needed if perturkey=3 or 4)
      integer          extvarkeys(100) ! extra output variables (1: yes, 0: no)

!     Outputs:

      real pres                        ! atmospheric pressure (Pa)
      real dens                        ! atmospheric density (kg.m-3)
      real temp                        ! atmospheric temperature (K)
      real zonwind                     ! zonal wind (m/s)
      real merwind                     ! meridional wind (m/s)
      real meanvar(5)                  ! unperturbed values of main meteorological variables
      real extvar(100)                 ! extra output variables
      real seedout                     ! current value of random generator seed index
      integer ier                      ! call_mcd status (=0 if all went well)
      
!     Local variables
      
      integer month,day,year,hour,minute,second   ! for Earth date input
      real ls                                     ! for user input of Ls (if using "martian time")
      character choice_date*1
      integer i 

      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! TIME    
      write(*,*) 'Do you want to use Earth date (e) or Mars date (m) ?'
 888  continue
      read(*,fmt='(a)') choice_date
      write(*,*) choice_date

      if (choice_date.eq.'e') then       ! Earth date
         datekey=0
         localtime=0. !compulsary with earth date
         write(*,*)'enter date : day month year hour minute second'
         read (*,*) day,month,year,hour,minute,second
         write(*,'(t2,i2,t4,a1,t5,i2,t7,a1,t8,i4,t14,i2,t16,a1,t17,i2,t19,a1,t20,i2)')  &
               day,'/',month,'/',year,hour,':',minute,':',second
         call julian(month,day,year,hour,minute,second,ier,date)
         write(*,'(''date, in Julian days : '',f16.8)') date
      else if (choice_date.eq.'m')  then     ! Mars date
         datekey=1
         write(*,*)'enter solar longitude Ls (deg)'
         read(*,*) ls
         write(*,*)'mars local time (0 < local time < 24) ?'
         read(*,*) localtime
         date=ls
      else
         goto 888
      end if

! LOCATION
      write(*,*)'choose your type for vertical coordinates:'
      write(*,*)'1 radius from center of planet (in meters )'
      write(*,*)'2 height above areoid (in meters )'
      write(*,*)'3 height above the surface (in meters )'
      write(*,*)'4 pressure level (in Pa )'
      read(*,*)  zkey
      write(*,'(i1)') zkey
      if (zkey.eq.1) then
       write(*,*)  'please enter the radius from center of planet (m)' 
      elseif (zkey.eq.2) then
       write(*,*)  'please enter height above areoid (m)' 
      elseif (zkey.eq.3) then
       write(*,*)  'please enter height above surface (m)' 
      elseif (zkey.eq.4) then
       write(*,*)  'please enter the pressure level in Pa' 
      endif
      
      read(*,*) z
      write(*,'(f10.3)') z

!      set hires flag
      write(*,*) ' high resolution? (1: yes, 0: no)'
      read(*,*) hireskey
      write(*,'(i1)') hireskey

      write(*,*)'latitude in deg ?'
      read(*,*)  lat
      write(*,'(f10.3)') lat
      write(*,*)'EAST longitude in deg ?'
      read(*,*)lon
      write(*,'(f10.3)') lon

! DUST/EUV scenario
      write(*,*)'dust/EUV scenarios ?'
      write(*,*)'1= Climatology     typical Mars year dust scenario'
      write(*,*)'                   average solar EUV conditions'  
      write(*,*)'2= Climatology     typical Mars year dust scenario'
      write(*,*)'                   minimum solar EUV conditions'  
      write(*,*)'3= Climatology     typical Mars year dust scenario '
      write(*,*)'                   maximum solar EUV conditions'
      write(*,*)'4= dust storm      constant dust opacity = 5'
      write(*,*)'                   min solar EUV conditions'
      write(*,*)'5= dust storm      constant dust opacity = 5' 
      write(*,*)'                   ave solar EUV conditions'  
      write(*,*)'6= dust storm      constant dust opacity = 5' 
      write(*,*)'                   max solar EUV conditions'  
      write(*,*)'7= warm scenario   warm scenario: dustier conditions'
      write(*,*)'                   max solar EUV conditions'  
      write(*,*)'8= cold scenario   cold scenario: clearer conditions'
      write(*,*)'                   min solar EUV conditions'  
      read(*,*)scena
      write(*,'(i1)') scena

      write(*,*)' perturbation : none = 1 ;  large scale = 2 ;',' small scale= 3 ; small+large = 4 ; n sig =5'
      read(*,*)perturkey
      write(*,'(i1)') perturkey
      write(*,*)' seedin and gwlength:'
      read(*,*)seedin,gwlength
      write(*,'(f10.0,f10.2)') seedin,gwlength
      write(*,*)'extra variables : yes = 1, no = 0'
      read(*,*)extvarkeys(1)
      write(*,'(i1)')extvarkeys(1)
      write(*,'(a1)') ' '
      do i=2,100
        ! propagate extvarkeys(1) to extvarkeys(:)
        extvarkeys(i)=extvarkeys(1)
      enddo

!     You can call call_mcd in a loop for a vertical profile, 
!     by just changing the value of xz at each iteration before
!     the call call_mcd() below.

      call call_mcd(zkey,z,lon,lat,hireskey,datekey,date,localtime,dset,scena,   &  
                    perturkey,seedin,gwlength,extvarkeys,                        &
                    pres,dens,temp,zonwind,merwind,meanvar,extvar,seedout,ier)

      if (ier.eq.0) then
         write(*,'(''p   =             '',1pe12.2,'' Pa'')') pres
         write(*,'(''rho =             '',1pe12.2,'' kg/m**3'')') dens
         write(*,'(''T   =             '',1pe12.2,'' K'')') temp
         write(*,'(''Zonal wind      = '',1pe12.2,'' m/s'')') zonwind
         write(*,'(''Meridional wind = '',1pe12.2,'' m/s'')') merwind
         write(*,'(a1)') ' '
         do i=1,5
            write(*,'(''meanvar('',i2,'') = '',1pe12.2)') i,meanvar(i)
         end do
         write(*,'(a1)') ' '
         if(extvarkeys(1).ne.0) then 
          do i = 1,85 ! write all 85 extra variables
           write(*,'(''extvar('',i2,'')  = '',1pe12.2)') i,extvar(i)
          enddo
         else ! write the first 13 extvar()
          do i=1,13
           write(*,'(''extvar('',i2,'')  = '',1pe12.2)') i,extvar(i)
          enddo
         end if
      else
         write(*,*)'CALL_MCD ERROR !!'
         write(*,*)'         returned error code: ', ier
      end if

      end

