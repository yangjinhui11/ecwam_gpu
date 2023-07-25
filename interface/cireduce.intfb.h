interface
SUBROUTINE CIREDUCE (WVPRPT, FF_NOW)
 use parkind_wave, only:&
 & jwrb
 USE YOWGRID , ONLY : NPROMA_WAM, NCHNK
 USE YOWPARAM , ONLY : NFRE
 USE YOWDRVTYPE ,ONLY: FREQUENCY,FORCING_FIELDS
    TYPE(FREQUENCY), INTENT(INOUT)            :: WVPRPT
      TYPE(FORCING_FIELDS), INTENT(IN)          :: FF_NOW
END SUBROUTINE CIREDUCE
end interface
