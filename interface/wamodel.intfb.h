interface
SUBROUTINE WAMODEL (NADV, LDSTOP, LDWRRE, BLK2GLO,             &
 &                  WVENVI, WVPRPT, FF_NOW, FF_NEXT, INTFLDS,  &
 &                  WAM2NEMO, NEMO2WAM, FL1, TIME1)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 USE YOWDRVTYPE , ONLY : WVGRIDGLO, ENVIRONMENT, FREQUENCY, FORCING_FIELDS, &
 & INTGT_PARAM_FIELDS, WAVE2OCEAN, OCEAN2WAVE
 USE YOWGRID , ONLY : NPROMA_WAM, NCHNK
 USE YOWPARAM , ONLY : NIBLO ,NANG ,NFRE
INTEGER(KIND=JWIM), INTENT(IN)                                           :: NADV
      LOGICAL, INTENT(INOUT)                                                   :: LDSTOP, LDWRRE
      TYPE(WVGRIDGLO), INTENT(IN)                                              :: BLK2GLO
      TYPE(ENVIRONMENT), INTENT(INOUT)                                         :: WVENVI
      TYPE(FREQUENCY), INTENT(INOUT)                                           :: WVPRPT
      TYPE(FORCING_FIELDS), INTENT(INOUT)                                      :: FF_NOW
      TYPE(FORCING_FIELDS), INTENT(IN)                                         :: FF_NEXT
      TYPE(INTGT_PARAM_FIELDS), INTENT(INOUT)                                  :: INTFLDS
      TYPE(WAVE2OCEAN), INTENT(INOUT)                                          :: WAM2NEMO
      TYPE(OCEAN2WAVE), INTENT(IN)                                             :: NEMO2WAM
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1

      REAL(KIND=JWRB), INTENT(INOUT) :: TIME1(3)
END SUBROUTINE WAMODEL
end interface
