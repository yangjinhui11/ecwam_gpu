interface
 SUBROUTINE SEBTMEAN (KIJS, KIJL, FL1, TB, TT, EBT)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 USE YOWPARAM , ONLY : NANG, NFRE
 INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG, NFRE), INTENT(IN) :: FL1
 REAL(KIND=JWRB), INTENT(IN) :: TB, TT
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: EBT
 END SUBROUTINE SEBTMEAN
end interface
