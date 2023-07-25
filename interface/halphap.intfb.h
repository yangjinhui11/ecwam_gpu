interface
SUBROUTINE HALPHAP(KIJS, KIJL, WAVNUM, COSWDIF, FL1, HALP)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 USE YOWPARAM , ONLY : NANG , NFRE
 INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG), INTENT(IN) :: COSWDIF
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: HALP
END SUBROUTINE HALPHAP
end interface
