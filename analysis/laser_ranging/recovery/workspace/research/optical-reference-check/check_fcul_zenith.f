C     Evaluate the unmodified IERS routine with its printed test inputs.
C     This driver is original Tudat reference-generation code.
      PROGRAM CHECKFCUL
      IMPLICIT NONE
      DOUBLE PRECISION TOTAL, HYDRO, WET

      CALL FCULZD_HPA(30.67166667D0, 2010.344D0,
     &     798.4188D0, 14.322D0, 0.532D0, TOTAL, HYDRO, WET)

      WRITE(*,'(A,ES26.17)') 'FCUL_ZTD = ', TOTAL
      WRITE(*,'(A,ES26.17)') 'FCUL_ZHD = ', HYDRO
      WRITE(*,'(A,ES26.17)') 'FCUL_ZWD = ', WET
      END
