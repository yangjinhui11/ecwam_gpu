interface
SUBROUTINE INITDPTHFLDS(WVENVI, WVPRPT, WVPRPT_LAND)
 USE YOWDRVTYPE , ONLY : ENVIRONMENT, FREQUENCY,FREQUENCY_LAND
 USE YOWGRID , ONLY : NPROMA_WAM, NCHNK
 USE YOWPARAM , ONLY : NFRE
    TYPE(ENVIRONMENT), INTENT(INOUT) :: WVENVI
      TYPE(FREQUENCY), INTENT(INOUT) :: WVPRPT
      TYPE(FREQUENCY_LAND), INTENT(INOUT) :: WVPRPT_LAND
END SUBROUTINE INITDPTHFLDS
end interface
