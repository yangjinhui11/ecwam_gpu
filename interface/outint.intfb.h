interface
 SUBROUTINE OUTINT(CDATE, CDATED, IFCST, BOUT)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 use yowcout , only:&
 & niprmout
 use yowgrid , only:&
 & nproma_wam,&
 & nchnk
 CHARACTER(LEN=14), INTENT(IN) :: CDATE, CDATED
 INTEGER(KIND=JWIM), INTENT(IN) :: IFCST
 REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NIPRMOUT, NCHNK), INTENT(IN) :: BOUT
 END SUBROUTINE OUTINT
end interface
