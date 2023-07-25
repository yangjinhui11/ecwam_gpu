interface
SUBROUTINE WNFLUXES (KIJS, KIJL,                       &
 &                   MIJ, RHOWGDFTH,                   &
 &                   CINV,                             &
 &                   SSURF, CICOVER,                   &
 &                   PHIWA,                            &
 &                   EM, F1, WSWAVE, WDWAVE,           &
 &                   UFRIC, AIRD,   &
 &                   NPHIEPS, NTAUOC, NSWH, NMWP, NEMOTAUX, &
 &                   NEMOTAUY, NEMOWSWAVE, NEMOPHIF, &
 &                   TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUOC, &
 &                   PHIOCD, PHIEPS, PHIAW, &
 &                   LNUPD)
 use parkind_wave, only:&
 & jwim,&
 & jwrb,JWRO
 use yowdrvtype , only:&
 & intgt_param_fields,&
 & wave2ocean
 USE YOWPARAM , ONLY : NANG ,NFRE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: MIJ
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: RHOWGDFTH
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NFRE), INTENT(IN) :: CINV
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL,NANG,NFRE), INTENT(IN) :: SSURF
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: CICOVER
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: PHIWA
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: EM, F1, WSWAVE, WDWAVE
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(IN) :: UFRIC, AIRD
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUOC
      REAL(KIND=JWRB), DIMENSION(KIJS:KIJL), INTENT(INOUT) :: PHIOCD, PHIEPS, PHIAW
      REAL(KIND=JWRO), DIMENSION(KIJS:KIJL), INTENT(INOUT) ::  NPHIEPS, NTAUOC, NSWH, NMWP, NEMOTAUX
      REAL(KIND=JWRO), DIMENSION(KIJS:KIJL), INTENT(INOUT) ::  NEMOTAUY, NEMOWSWAVE, NEMOPHIF
      LOGICAL, INTENT(IN) :: LNUPD
END SUBROUTINE WNFLUXES
end interface
