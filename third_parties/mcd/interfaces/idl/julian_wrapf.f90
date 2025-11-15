! julian_wrapf.f
!
! This function is the wrapper routine called by IDL
!
      INTEGER*4 FUNCTION julian_wrapf(argc, argv)
      IMPLICIT NONE
!
! The declaration below is integer*8 for machines with
! 64-bit address spaces and integer*4 for 32-bit
!
      INTEGER*8 argc, argv(*)
      
!
! Call The Fortran Routine:
! julian(month,day,year,hour,minute,second,ierr,date)
!
      CALL julian(%val(argv(1)), %val(argv(2)), %val(argv(3)),    &
       %val(argv(4)), %val(argv(5)), %val(argv(6)), %val(argv(7)),&
       %val(argv(8)))
!
! Give the function vecadd a value for return to IDL, this
! facilitates the checking of whether the routine executes.
!
      julian_wrapf=1
!
      END
