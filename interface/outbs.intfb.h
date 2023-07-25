interface
SUBROUTINE OUTBS (MIJ, FL1, XLLWS, &
 & WVPRPT, WVENVI, FF_NOW, INTFLDS, NEMO2WAM, &
 & BOUT)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY, FORCING_FIELDS,  &
     &                         INTGT_PARAM_FIELDS, OCEAN2WAVE
 use yowcout , only:&
 & niprmout
 USE YOWGRID , ONLY : NPROMA_WAM, NCHNK
 USE YOWPARAM , ONLY : NANG ,NFRE
INTEGER(KIND=JWIM), DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN)          :: MIJ
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(IN) :: XLLWS
      TYPE(FREQUENCY), INTENT(IN)                                           :: WVPRPT
      TYPE(ENVIRONMENT), INTENT(IN)                                         :: WVENVI
      TYPE(FORCING_FIELDS), INTENT(INOUT)                                   :: FF_NOW
      TYPE(INTGT_PARAM_FIELDS), INTENT(IN)                                  :: INTFLDS
      TYPE(OCEAN2WAVE), INTENT(IN)                                          :: NEMO2WAM
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NIPRMOUT, NCHNK), INTENT(OUT)  :: BOUT
END SUBROUTINE OUTBS
end interface
