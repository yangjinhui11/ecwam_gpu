interface
 SUBROUTINE MEANSQS_LF(NFRE_EFF, KIJS, KIJL, F, WAVNUM, XMSS)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 USE YOWPARAM , ONLY : NANG ,NFRE
 INTEGER(KIND=JWIM), INTENT(IN) :: NFRE_EFF
 INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: F
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: XMSS
 END SUBROUTINE MEANSQS_LF
end interface
