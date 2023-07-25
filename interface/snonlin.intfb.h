interface
 SUBROUTINE SNONLIN (KIJS, KIJL, FL1, FLD, SL, WAVNUM, DEPTH, AKMEAN)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 USE YOWPARAM , ONLY : NANG ,NFRE
 INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT):: FLD, SL
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: DEPTH, AKMEAN
 END SUBROUTINE SNONLIN
end interface
