interface
 SUBROUTINE STRESSO (KIJS, KIJL, MIJ, RHOWGDFTH, &
 & FL1, SL, SPOS, &
 & CINV, &
 & WDWAVE, UFRIC, Z0M, AIRD, RNFAC, &
 & COSWDIF, SINWDIF2, &
 & TAUW, TAUWDIR, PHIWA, LLPHIWA)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 USE YOWPARAM , ONLY : NANG ,NFRE
 INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
 INTEGER(KIND=JWIM), INTENT(IN) :: MIJ(KIJS:KIJL)
 REAL(KIND=JWRB),DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: RHOWGDFTH
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1, SL, SPOS
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: CINV
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: WDWAVE, UFRIC, Z0M, AIRD, RNFAC
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG), INTENT(IN) :: COSWDIF, SINWDIF2
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: TAUW, TAUWDIR, PHIWA
 LOGICAL, INTENT(IN) :: LLPHIWA
 END SUBROUTINE STRESSO
end interface
