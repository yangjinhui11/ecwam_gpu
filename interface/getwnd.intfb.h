interface
SUBROUTINE GETWND (BLK2LOC,                               &
 &                 NXS, NXE, NYS, NYE, FIELDG,            &
 &                 WVENVI,                                &
 &                 FF_NOW,                                &
 &                 CDTWIS, LWNDFILE, LCLOSEWND, IREAD,    &
 &                 LWCUR, NEMO2WAM,                       &
 &                 ICODE_WND)
 use parkind_wave, only:&
 & jwim,&
 & jwrb,&
 & jwro
 USE YOWDRVTYPE , ONLY :WVGRIDLOC, FORCING_FIELDS,ENVIRONMENT,OCEAN2WAVE
 USE YOWGRID , ONLY : NPROMA_WAM, NCHNK
 
      TYPE(WVGRIDLOC), INTENT(IN) :: BLK2LOC
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FIELDG
      TYPE(ENVIRONMENT), INTENT(IN) :: WVENVI

      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NOW
      CHARACTER(LEN=14), INTENT(IN) :: CDTWIS
      LOGICAL, INTENT(IN) :: LWNDFILE, LCLOSEWND
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      LOGICAL, INTENT(IN) :: LWCUR
      TYPE(OCEAN2WAVE), INTENT(IN) :: NEMO2WAM
      INTEGER(KIND=JWIM), INTENT(OUT) :: ICODE_WND
END SUBROUTINE GETWND
end interface
