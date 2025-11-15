      
      subroutine julian(month,day,year,hour,minute,second,ierr,date)


      implicit none
!
!     Given Earth date and time compute and local time on Mars
!     Updated version by B. Dolla and F. Forget, 2005
!     Inputs
!
      integer,intent(in)  :: month   
      integer,intent(in)  :: day
      integer,intent(in)  :: year
      integer,intent(in)  :: hour       !All Earth GMT values
      integer,intent(in)  :: minute
      integer,intent(in)  :: second
!
!     Output
!
      integer,intent(out) :: ierr       !0 if ok >0 if there is a problem
      real*8, intent(out) :: date       !Julian date
!
!     Local
!
      integer nday
      integer daynumber(12)   !days for months of the year
      data    daynumber/0,31,59,90,120,151,181,212,243,273,304,334/
      integer jul
      integer j
!
!     Check ranges
!
      ierr=0
      if ((month.lt.1).or.(month.gt.12)) then
        ierr=1
        return
      endif
      if ((day.lt.1).or.(day.gt.31)) then
        ierr=2
        return
      endif
      if (year.lt.1) then
        ierr=3
        return
      endif
      if ((hour.lt.0).or.(hour.gt.23)) then
        ierr=4
        return
      endif
      if ((minute.lt.0).or.(minute.gt.59)) then
        ierr=5
        return
      endif
      if ((second.lt.0).or.(second.gt.60)) then
        ierr=6
        return
      endif
!
!     Calculate Julian date
!
      nday=daynumber(month)+day-1
!
!     Correct for leap year
!     We use the followings conventions
!     GREGORIAN CALENDAR: a year is bissextil if it is a multiple
!     of 4 but not of 100 or if it is a multiple of 400.
!     JULIAN CALENDAR: a year is bissextil if it is a multiple of 4
!
!     The JULIAN calendar ends on the 4th october 1582
!     The GREGORIAN calendar begins on the 15th october 1582
!     Hence there are 10 days missing... e.g. the 10th october 1582 does not exist!!
!

      jul=0
      IF (year.LT.1582) jul=1
      IF ((year.EQ.1582).AND.(month.LT.10)) jul=1
      IF ((year.EQ.1582).AND.(month.EQ.10).AND.(day.LT.15)) jul=1
      
      IF (jul.EQ.0) THEN
         IF ((mod(year,4).EQ.0).AND.(mod(year,100).NE.0).AND.(month.GT.2)) nday=nday+1
         IF ((mod(year,400).EQ.0).AND.(month.GT.2)) nday=nday+1
      ENDIF
      IF (jul.EQ.1) THEN
         IF ((mod(year,4).EQ.0).AND.(month.GT.2)) nday=nday+1
         nday=nday+10
      ENDIF
!
!     We use 1968 as a year of reference, the julian date being 2.4398565d6
!
      IF (year.GT.1968) THEN
         DO j=1968,year-1,1
            nday=nday+365
            IF ((mod(j,4).EQ.0).AND.(mod(j,100).NE.0)) nday=nday+1
            IF (mod(j,400).EQ.0) nday=nday+1
         ENDDO
      ENDIF

      IF (year.LT.1968) THEN
         jul=1
         DO j=year,1967,1
            IF (j.GT.1581) jul=0
            IF (jul.EQ.0) THEN
               nday=nday-365
               IF ((mod(j,4).EQ.0).AND.(mod(j,100).NE.0)) nday=nday-1
               IF (mod(j,400).EQ.0) nday=nday-1
            ENDIF
            IF (jul.EQ.1) THEN
               nday=nday-365
               IF (mod(j,4).EQ.0) nday=nday-1
            ENDIF
         ENDDO
      ENDIF
!
!     Compute Julian date
!
      date=2.4398565d6+nday
      date=date+hour/24.0d0+minute/1.440d3+second/8.6400d4
!
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
