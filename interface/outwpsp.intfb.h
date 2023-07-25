interface
SUBROUTINE OUTWPSP (FL1, FF_NOW)
 use parkind_wave, only:&
 & jwrb
 USE YOWDRVTYPE , ONLY : FORCING_FIELDS
 use yowgrid , only:&
 & nproma_wam,&
 & nchnk
 use yowparam , only:&
 & nang,&
 & nfre
REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(IN) :: FL1
      TYPE(FORCING_FIELDS), INTENT(IN) :: FF_NOW
END SUBROUTINE OUTWPSP
end interface
