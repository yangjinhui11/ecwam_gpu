interface
 SUBROUTINE FKMEAN (KIJS, KIJL, FL1, WAVNUM, &
 & EM, FM1, F1, AK, XK)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 USE YOWPARAM , ONLY : NANG ,NFRE
 INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: WAVNUM
 REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(OUT) :: EM, FM1, F1, AK, XK
 END SUBROUTINE FKMEAN
end interface
