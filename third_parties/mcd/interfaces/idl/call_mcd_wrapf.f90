! call_mcd_wrapf.f
!
! This function is the wrapper routine called by IDL
!
      INTEGER*4 FUNCTION call_mcd_wrapf(argc, argv)
      IMPLICIT NONE
!
! The declaration below is integer*8 for machines with
! 64-bit address spaces and integer*4 for 32-bit
!
      INTEGER*8, target :: argc, argv(*)
!
! Call The Fortran Routine:
!
      CALL call_mcd4idl(%val(argv( 1)), %val(argv( 2)), %val(argv( 3)),  &
       %val(argv( 4)), %val(argv( 5)), %val(argv( 6)), %val(argv( 7)),   &
       %val(argv( 8)), %val(argv( 9)), %val(argv(10)), %val(argv(11)),   &
       %val(argv(12)), %val(argv(13)), %val(argv(14)), %val(argv(15)),   &
       %val(argv(16)), %val(argv(17)), %val(argv(18)), %val(argv(19)),   &
       %val(argv(20)), %val(argv(21)), %val(argv(22)), %val(argv(23)))
!
! Give the function a value for return to IDL, this
! facilitates the checking of whether the routine executes.
!
      call_mcd_wrapf=1
!
      END
