interface
SUBROUTINE TAU_PHI_HF(KIJS, KIJL, MIJ, LTAUWSHELTER, UFRIC, Z0M, &
 & FL1, AIRD, RNFAC, &
 & COSWDIF, SINWDIF2, &
 & UST, TAUHF, PHIHF, LLPHIHF)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 USE YOWPARAM , ONLY : NANG ,NFRE
 INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
 INTEGER(KIND=JWIM), INTENT(IN) :: MIJ(KIJS:KIJL)
 LOGICAL, INTENT(IN) :: LTAUWSHELTER
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: UFRIC, Z0M
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: FL1
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: AIRD, RNFAC
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL, NANG), INTENT(IN) :: COSWDIF, SINWDIF2
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: UST
 REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(OUT) :: TAUHF, PHIHF
 LOGICAL, INTENT(IN) :: LLPHIHF
END SUBROUTINE TAU_PHI_HF
end interface
