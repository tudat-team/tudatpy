; julian.pro
;
; This function performs checking of variables before a call
; is made to an external routine via CALL_EXTERNAL.
;
; The data types need to be the same as those expected by the
; external routine.
;
FUNCTION JULIAN, mt, dy, yr, hr, mn, sc, ie, dt

month = long(mt)
day   = long(dy)
year  = long(yr)
hour  = long(hr)
minute= long(mn)
second= long(sc)
ierr  = long(ie)
date  = double(dt)
;
; Call the external routine:
result = CALL_EXTERNAL('call_mcd.so', 'julian_wrapf_', $
         month, day, year, hour, minute, second, ierr, date)
IF (result NE 1) THEN MESSAGE, 'Error calling external routine julian'

dt = date
; Return result:
RETURN, ierr
END
