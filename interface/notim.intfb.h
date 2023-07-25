interface
SUBROUTINE NOTIM (CDTWIS, CDTWIE,              &
 &                NXS, NXE, NYS, NYE, FIELDG,  &
 &                BLK2LOC, WVENVI, FF_NEXT,    &
 &                IREAD, LWCUR, NEMO2WAM)
 use parkind_wave, only:&
 & jwim
 USE YOWDRVTYPE , ONLY : WVGRIDLOC, ENVIRONMENT, FORCING_FIELDS, OCEAN2WAVE
 USE YOWGRID , ONLY : NPROMA_WAM, NCHNK
CHARACTER(LEN=14), INTENT(IN) :: CDTWIS, CDTWIE
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FIELDG
      TYPE(WVGRIDLOC), INTENT(IN) :: BLK2LOC
      TYPE(ENVIRONMENT), INTENT(INOUT) :: WVENVI
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NEXT
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      LOGICAL, INTENT(IN) :: LWCUR
      TYPE(OCEAN2WAVE), INTENT(INOUT) :: NEMO2WAM
END SUBROUTINE NOTIM
end interface
