interface
 SUBROUTINE CIWAF (KIJS, KIJL, CGROUP, CICOVER, CITHICK, CIWA)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 USE YOWPARAM , ONLY : NFRE
 INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
 REAL(KIND=JWRB),DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: CGROUP
 REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: CICOVER
 REAL(KIND=JWRB),DIMENSION(KIJS:KIJL), INTENT(IN) :: CITHICK
 REAL(KIND=JWRB),DIMENSION(KIJS:KIJL,NFRE), INTENT(OUT) :: CIWA
 END SUBROUTINE CIWAF
end interface
