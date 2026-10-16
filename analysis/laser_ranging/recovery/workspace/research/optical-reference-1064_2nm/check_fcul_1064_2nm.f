C     Actual Grasse vacuum wavelength; original McDonald reference weather.
      PROGRAM CHECKFCULIR
      IMPLICIT NONE
      DOUBLE PRECISION TOTAL, HYDRO, WET
      CALL FCULZD_HPA(30.67166667D0, 2010.344D0,
     &     798.4188D0, 14.322D0, 1.0642D0, TOTAL, HYDRO, WET)
      WRITE(*,'(A,ES26.17)') 'FCUL_ZTD = ', TOTAL
      WRITE(*,'(A,ES26.17)') 'FCUL_ZHD = ', HYDRO
      WRITE(*,'(A,ES26.17)') 'FCUL_ZWD = ', WET
      END
