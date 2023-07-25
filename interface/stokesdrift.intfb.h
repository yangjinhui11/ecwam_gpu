interface
 SUBROUTINE STOKESDRIFT(KIJS, KIJL, FL1, STOKFAC, WSWAVE, WDWAVE, CICOVER, USTOKES, VSTOKES)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 use yowparam , only:&
 & nang,&
 & nfre
 INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: STOKFAC
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: WSWAVE, WDWAVE, CICOVER
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: USTOKES, VSTOKES
 END SUBROUTINE STOKESDRIFT
end interface