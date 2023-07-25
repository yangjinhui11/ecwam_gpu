interface
 SUBROUTINE WEFLUX (KIJS, KIJL, FL1, CGROUP, &
 & NFRE, NANG, DFIM, DELTH, &
 & COSTH, SINTH, &
 & WEFMAG, WEFDIR)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: CGROUP
 INTEGER(KIND=JWIM), INTENT(IN) :: NFRE, NANG
 REAL(KIND=JWRB), INTENT(IN) :: DELTH
 REAL(KIND=JWRB), DIMENSION(NFRE), INTENT(IN) :: DFIM
 REAL(KIND=JWRB), DIMENSION(NANG), INTENT(IN) :: COSTH, SINTH
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: WEFMAG, WEFDIR
 END SUBROUTINE WEFLUX
end interface
