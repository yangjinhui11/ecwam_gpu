interface
SUBROUTINE READWIND (CDTWIR, FILNM, LLNOTOPENED, IREAD, &
 & NXS, NXE, NYS, NYE, FIELDG)
 use parkind_wave, only:&
 & jwim
 USE YOWDRVTYPE , ONLY : FORCING_FIELDS
INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      CHARACTER(LEN=14), INTENT(INOUT) :: CDTWIR
      CHARACTER(LEN=24), INTENT(INOUT) :: FILNM
      LOGICAL, INTENT(INOUT) :: LLNOTOPENED
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FIELDG
END SUBROUTINE READWIND
end interface