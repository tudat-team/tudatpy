#include "fintrf.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Gateway routine for access to call_mcd in MCD.F90
!
      subroutine mexFunction(nlhs,plhs,nrhs,prhs)
      use MCD
      IMPLICIT NONE
      mwpointer prhs(*),plhs(*),x_pr,y_pr
      mwpointer mxGetPr
      double precision mxGetScalar
      integer   nlhs,nrhs
      
!!!!!!!!!!!!!!! CALL_MCD arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! inputs:
      integer zkey            
      real    xz              
      real    xlon            
      real    xlat            
      integer hireskey        
      integer datekey         
      double  precision xdate 
      real    localtime       
      character*50 dset       
      integer dust            
      integer perturkey       
      real    seedin          
      real    gwlength        
      integer extvarkeys(100) 
! outputs:
      real    pres            
      real    ro              
      real    temp            
      real    u               
      real    v               
      real    meanvar(5)      
      real    extvar(100)     
      real    seedout         
      integer ier             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !integer*4        i4buffer(100)
      double precision dbuffer(100)
      double precision dbuffer5(5)
      mwpointer        dset_pr
      mwpointer        extvarkeys_pr
      real*8           r8buffer(100)
!
      zkey            = IDINT(mxGetScalar(prhs( 1)))
      xz              = SNGL( mxGetScalar(prhs( 2)))
      xlon            = SNGL( mxGetScalar(prhs( 3)))
      xlat            = SNGL( mxGetScalar(prhs( 4)))
      hireskey        = IDINT(mxGetScalar(prhs( 5)))
      datekey         = IDINT(mxGetScalar(prhs( 6)))
      xdate           =       mxGetScalar(prhs( 7))
      localtime       = SNGL( mxGetScalar(prhs( 8)))
      dset_pr         =                   prhs( 9)
      dust            = IDINT(mxGetScalar(prhs(10)))
      perturkey       = IDINT(mxGetScalar(prhs(11)))
      seedin          = SNGL( mxGetScalar(prhs(12)))
      gwlength        = SNGL( mxGetScalar(prhs(13)))
      extvarkeys_pr   =           mxGetPr(prhs(14))

      call mxArrayToCharacter(dset_pr,dset)
      !call mxCopyPtrToInteger4(extvarkeys_pr,i4buffer,100)
      call mxCopyPtrToReal8(extvarkeys_pr,r8buffer,100)
      !extvarkeys = INT(i4buffer)
      extvarkeys = INT(r8buffer)
      !
      ! call the calculation routine
      !
      call call_mcd(zkey,xz,xlon,xlat,hireskey,datekey,xdate,localtime,dset,dust,   &
                    perturkey,seedin,gwlength,extvarkeys,pres,ro,temp,u,v,          &
                    meanvar,extvar,seedout,ier)

      call mxCopyReal8ToPtr(DBLE(pres),mxGetPr(prhs(15)),1)
      call mxCopyReal8ToPtr(DBLE(ro  ),mxGetPr(prhs(16)),1)
      call mxCopyReal8ToPtr(DBLE(temp),mxGetPr(prhs(17)),1)
      call mxCopyReal8ToPtr(DBLE(u   ),mxGetPr(prhs(18)),1)
      call mxCopyReal8ToPtr(DBLE(v   ),mxGetPr(prhs(19)),1)

      dbuffer5 = DBLE(meanvar)
      call mxCopyReal8ToPtr(dbuffer5,mxGetPr(prhs(20)),5)
      dbuffer =  DBLE(extvar)
      call mxCopyReal8ToPtr(dbuffer,mxGetPr(prhs(21)),100)

      call mxCopyReal8ToPtr(DBLE(seedout),mxGetPr(prhs(22)),1)
      call mxCopyReal8ToPtr(DBLE(ier    ),mxGetPr(prhs(23)),1)

      return 
      end
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copies MATLAB character string to Fortran character string
! Input:   mx  -->  Address of mxArray containing char string
! Outputs: s   -->  Fortran character string (blank padded at end)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine mxArrayToCharacter( mx, s )
      implicit none
      mwPointer mx
      character*(*) s
      mwPointer, external :: mxGetData
      integer*4, external :: mxIsChar
      mwPointer, external :: mxGetNumberOfElements
      mwPointer pr
      mwPointer n
!
! Check for NULL input
!
      if(mx==0) then
       s = ' '
       return
      endif
!
! Check for char string input
!
      if( mxIsChar(mx) == 0 ) then
       s = ' '
       return
      endif
!
! Get char data pointer and number of characters
!
      pr = mxGetData(mx)
      n  = mxGetNumberOfElements(mx)
!
! Call utility routine to treat MATLAB char data as integer*2
!
      call mxCopyI2toCharacter( %VAL(pr), n, s )
      return
      end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copy 1 character at a time, using the 2-byte MATLAB char value as
! an integer value containing the ASCII character code. If there are
! not enough characters to fill s, blank pad the end.
!
      subroutine mxCopyI2toCharacter( I2, n, s )
      implicit none
      mwPointer n
      integer*2 I2(n)
      character*(*) s
      mwPointer i, m
!
      m = len(s)
      do i=1,min(m,n)
       s(i:i) = achar(I2(i))
      enddo
      if( n < m )then
       s(n+1:) = ' '
      end if
      return
      end subroutine
