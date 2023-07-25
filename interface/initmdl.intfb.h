interface
SUBROUTINE INITMDL (NADV,                                 &
 &                  IREAD,                                &
 &                  BLK2GLO,  BLK2LOC,                    &
 &                  WVENVI, WVPRPT, FF_NOW,               &
 &                  FL1,                                  &
 &                  NFIELDS, NGPTOTG, NC, NR,             &
 &                  FIELDS, LWCUR, MASK_IN, PRPLRADI,     &
 &                  NEMO2WAM)
 use parkind_wave, only:&
 & jwim,&
 & jwrb,JWRU
 use yowdrvtype , only:&
 & wvgridglo,&
 & wvgridloc,&
 & environment,&
 & frequency,&
 & forcing_fields,&
 & ocean2wave
 use yowgrid , only:&
 & nproma_wam,&
 & nchnk
 use yowparam , only:&
 & nang,&
 & nfre,&
 & niblo
    INTEGER(KIND=JWIM), INTENT(OUT) :: NADV
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      TYPE(WVGRIDLOC), INTENT(IN) :: BLK2LOC
      TYPE(ENVIRONMENT), INTENT(INOUT) :: WVENVI
      TYPE(FREQUENCY), INTENT(INOUT) :: WVPRPT
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NOW
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1
      INTEGER(KIND=JWIM), INTENT(IN) :: NFIELDS
      INTEGER(KIND=JWIM), INTENT(IN) :: NGPTOTG
      INTEGER(KIND=JWIM), INTENT(IN) :: NC
      INTEGER(KIND=JWIM), INTENT(IN) :: NR
      REAL(KIND=JWRB),DIMENSION(NGPTOTG,NFIELDS), INTENT(IN) :: FIELDS
      LOGICAL, INTENT(IN) :: LWCUR
      INTEGER(KIND=JWIM),DIMENSION(NGPTOTG), INTENT(INOUT) :: MASK_IN
      REAL(KIND=JWRB), INTENT(IN) :: PRPLRADI
      TYPE(OCEAN2WAVE), INTENT(INOUT) :: NEMO2WAM
END SUBROUTINE INITMDL
end interface
