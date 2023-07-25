interface
SUBROUTINE OUTBC (FL1, BLK2GLO, IU19)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 USE YOWDRVTYPE , ONLY : WVGRIDGLO
 use yowparam , only:&
 & niblo,&
 & nang,&
 & nfre
 use yowcpbo , only:&
 & gbounc
 USE YOWGRID , ONLY : NPROMA_WAM, NCHNK
REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(IN) :: FL1
      TYPE(WVGRIDGLO), INTENT(IN)                         :: BLK2GLO
      INTEGER(KIND=JWIM),DIMENSION(GBOUNC), INTENT(IN) :: IU19
END SUBROUTINE OUTBC
end interface
