interface
SUBROUTINE BOUINPT (IU02, FL1, NSTART, NEND)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 use yowgrid , only:&
 & nproma_wam,&
 & nchnk
 use yowmpp , only:&
 & nproc
 USE YOWPARAM , ONLY : NANG ,NFRE
 INTEGER(KIND=JWIM), INTENT(INOUT) :: IU02
 REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1
 INTEGER(KIND=JWIM), DIMENSION(NPROC), INTENT(IN) :: NSTART, NEND
END SUBROUTINE BOUINPT
end interface
