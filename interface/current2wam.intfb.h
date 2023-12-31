interface
SUBROUTINE CURRENT2WAM (FILNM, IREAD, CDATEIN,        &
     &                        BLK2LOC,                      &
     &                        NXS, NXE, NYS, NYE, FIELDG,   &
     &                        WVENVI)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 USE YOWDRVTYPE , ONLY : FORCING_FIELDS, WVGRIDLOC, ENVIRONMENT
 USE YOWGRID , ONLY : NPROMA_WAM, NCHNK
    INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      CHARACTER(LEN=24), INTENT(IN) :: FILNM
      CHARACTER(LEN=14), INTENT(INOUT) :: CDATEIN
      TYPE(WVGRIDLOC), INTENT(IN) :: BLK2LOC
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(IN) :: FIELDG
      TYPE(ENVIRONMENT), INTENT(OUT) :: WVENVI
 END SUBROUTINE CURRENT2WAM
end interface
