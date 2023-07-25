interface
SUBROUTINE PROPAG_WAM (BLK2GLO, WVENVI, WVPRPT, FL1)
 use parkind_wave, only:&
 & jwrb
 USE YOWDRVTYPE , ONLY : WVGRIDGLO, ENVIRONMENT, FREQUENCY
 use yowgrid , only:&
 & nproma_wam,&
 & nchnk
 use yowparam , only:&
 & nang,&
 & nfre,&
 & niblo
TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      TYPE(ENVIRONMENT), INTENT(IN) :: WVENVI
      TYPE(FREQUENCY), INTENT(IN) :: WVPRPT
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1
END SUBROUTINE PROPAG_WAM
end interface
