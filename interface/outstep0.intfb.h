interface
SUBROUTINE OUTSTEP0 (WVENVI, WVPRPT, FF_NOW, INTFLDS,  &
 &                   WAM2NEMO, NEMO2WAM, FL1)
 use parkind_wave, only:&
 & jwrb
 USE YOWDRVTYPE , ONLY : ENVIRONMENT, FREQUENCY, FORCING_FIELDS, &
 & INTGT_PARAM_FIELDS, WAVE2OCEAN, OCEAN2WAVE
 USE YOWGRID , ONLY : NPROMA_WAM, NCHNK
 use yowparam , only:&
 & nang,&
 & nfre
TYPE(ENVIRONMENT), INTENT(INOUT)                                     :: WVENVI
      TYPE(FREQUENCY), INTENT(INOUT)                                 :: WVPRPT
      TYPE(FORCING_FIELDS), INTENT(INOUT)                            :: FF_NOW
      TYPE(INTGT_PARAM_FIELDS), INTENT(INOUT)                        :: INTFLDS
      TYPE(WAVE2OCEAN), INTENT(INOUT)                                :: WAM2NEMO
      TYPE(OCEAN2WAVE), INTENT(IN)                                   :: NEMO2WAM
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1
END SUBROUTINE OUTSTEP0
end interface
