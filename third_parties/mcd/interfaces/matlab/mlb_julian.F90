#include "fintrf.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Gateway routine for access to julian.F
!
      subroutine mexFunction(nlhs,plhs,nrhs,prhs)
      use MCD
      IMPLICIT NONE
      mwpointer prhs(*),plhs(*),x_pr,y_pr
      mwpointer mxGetPr
      double precision mxGetScalar
      integer   nlhs,nrhs

!!!!!!!!!!!!!!! JULIAN arguments   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! inputs:
      integer month   
      integer day
      integer year
      integer hour
      integer minute
      integer second
! outputs:
      integer ierr
      real*8  date
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      month           = IDINT(mxGetScalar(prhs( 1)))
      day             = IDINT(mxGetScalar(prhs( 2)))
      year            = IDINT(mxGetScalar(prhs( 3)))
      hour            = IDINT(mxGetScalar(prhs( 4)))
      minute          = IDINT(mxGetScalar(prhs( 5)))
      second          = IDINT(mxGetScalar(prhs( 6)))
      !
      ! call the calculation routine
      !
      call julian(month,day,year,hour,minute,second,ierr,date)

      call mxCopyReal8ToPtr(DBLE(ierr   ),mxGetPr(prhs(7)),1)
      call mxCopyReal8ToPtr(DBLE(date   ),mxGetPr(prhs(8)),1)

      return 
      end
