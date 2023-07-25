interface
SUBROUTINE PREWIND (BLK2LOC, WVENVI, FF_NOW, FF_NEXT,       &
 &                  NXS, NXE, NYS, NYE, LLINIT_FIELDG,      &
 &                  LLINIT, IREAD,                          &
 &                  NFIELDS, NGPTOTG, NC, NR,               &
 &                  FIELDS, LWCUR, MASK_IN,                 &
 &                  NEMO2WAM)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
USE YOWDRVTYPE  , ONLY : WVGRIDLOC, ENVIRONMENT, FORCING_FIELDS, OCEAN2WAVE
 USE YOWGRID , ONLY : NPROMA_WAM, NCHNK
TYPE(WVGRIDLOC), INTENT(IN)         :: BLK2LOC
      TYPE(ENVIRONMENT), INTENT(INOUT)    :: WVENVI
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NOW
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NEXT
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      LOGICAL, INTENT(IN) :: LLINIT_FIELDG
      LOGICAL, INTENT(IN) :: LLINIT
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      INTEGER(KIND=JWIM), INTENT(IN) :: NFIELDS
      INTEGER(KIND=JWIM), INTENT(IN) :: NGPTOTG
      INTEGER(KIND=JWIM), INTENT(IN) :: NC
      INTEGER(KIND=JWIM), INTENT(IN) :: NR
      REAL(KIND=JWRB),DIMENSION(NGPTOTG, NFIELDS), INTENT(IN) :: FIELDS
      LOGICAL, INTENT(IN) :: LWCUR
      INTEGER(KIND=JWIM),DIMENSION(NGPTOTG), INTENT(INOUT)  :: MASK_IN
      TYPE(OCEAN2WAVE), INTENT(INOUT) :: NEMO2WAM
END SUBROUTINE PREWIND
end interface
