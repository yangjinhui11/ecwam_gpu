interface
 SUBROUTINE SDISSIP (KIJS, KIJL, FL1, FLD, SL, &
 & INDEP, WAVNUM, XK2CG, &
 & EMEAN, F1MEAN, XKMEAN, &
 & UFRIC, COSWDIF, RAORW)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 USE YOWPARAM , ONLY : NANG ,NFRE
 INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
 INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: INDEP
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM, XK2CG
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: EMEAN, F1MEAN, XKMEAN
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: UFRIC, RAORW
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG), INTENT(IN) :: COSWDIF
 END SUBROUTINE SDISSIP
end interface
