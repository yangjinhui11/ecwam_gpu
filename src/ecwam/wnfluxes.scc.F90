! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE WNFLUXES_SCC_MOD
  CONTAINS
  SUBROUTINE WNFLUXES_SCC (KIJS, KIJL, MIJ, RHOWGDFTH, CINV, SSURF, CICOVER, PHIWA, EM, F1, WSWAVE, WDWAVE, UFRIC, AIRD,  &
  & NPHIEPS, NTAUOC, NSWH, NMWP, NEMOTAUX, NEMOTAUY, NEMOWSWAVE, NEMOPHIF, TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUOC, PHIOCD,  &
  & PHIEPS, PHIAW, LNUPD)
    
    ! ----------------------------------------------------------------------
    
    !**** *WNFLUXES* - WAVE FLUXES CALCULATION
    
    !*    PURPOSE.
    !     --------
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *WNFLUXES* (KIJS, KIJL,
    !    &                     MIJ, RHOWGDFTH,
    !    &                     CINV,
    !    &                     SSURF, CICOVER,
    !    &                     PHIWA,
    !    &                     EM, F1, WSWAVE, WDWAVE,
    !    &                     UFRIC, AIRD, INTFLDS, WAM2NEMO,
    !    &                     LNUPD)
    !          *KIJS*    - INDEX OF FIRST GRIDPOINT.
    !          *KIJL*    - INDEX OF LAST GRIDPOINT.
    !          *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE
    !         *RHOWGDFTH    - WATER DENSITY * G * DF * DTHETA
    !                         FOR TRAPEZOIDAL INTEGRATION BETWEEN FR(1) and FR(MIJ)
    !                         !!!!!!!!  RHOWGDFTH=0 FOR FR > FR(MIJ)
    !          *CINV*   - INVERSE PHASE SPEED.
    !          *SSURF*  - CONTRIBUTION OF ALL SOURCE TERMS ACTING ON
    !                     THE SURFACE MOMENTUM AND ENERGY FLUXES.
    !          *CICOVER*- SEA ICE COVER.
    !          *PHIWA*  - ENERGY FLUX FROM WIND INTO WAVES INTEGRATED
    !                     OVER THE FULL FREQUENCY RANGE.
    !          *EM*     - MEAN WAVE VARIANCE.
    !          *F1*     - MEAN WAVE FREQUENCY BASED ON f*F INTEGRATION.
    !          *WSWAVE* - WIND SPEED IN M/S.
    !          *WDWAVE* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC CONVENTION
    !          *UFRIC*  - FRICTION VELOCITY IN M/S.
    !          *AIRD*   - AIR DENSITY IN KG/M3.
    !          *INTFLDS*-  INTEGRATED/DERIVED PARAMETERS
    !          WAM2NEMO*- WAVE FIELDS PASSED TO NEMO
    !          *LNUPD*  - UPDATE NEMO FIELDS.
    
    ! ----------------------------------------------------------------------
    
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU, JWRO
    USE YOWDRVTYPE, ONLY: FORCING_FIELDS, INTGT_PARAM_FIELDS, WAVE2OCEAN
    
    USE YOWALTAS, ONLY: EGRCRV, AFCRV, BFCRV
    USE YOWCOUP, ONLY: LWNEMOCOU, LWNEMOTAUOC
    USE YOWFRED, ONLY: FR, COSTH, SINTH
    USE YOWICE, ONLY: LICERUN, LWAMRSETCI, CITHRSH, CIBLOCK
    
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: TAUOCMIN, TAUOCMAX, PHIEPSMIN, PHIEPSMAX, EPSUS, EPSU10, G, ZPI
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    INTEGER(KIND=JWIM), INTENT(IN), DIMENSION(KIJL) :: MIJ
    
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: RHOWGDFTH
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: CINV
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: SSURF
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: CICOVER
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: PHIWA
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: EM, F1, WSWAVE, WDWAVE
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: UFRIC, AIRD
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUOC
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: PHIOCD, PHIEPS, PHIAW
    REAL(KIND=JWRO), INTENT(INOUT), DIMENSION(KIJL) :: NPHIEPS, NTAUOC, NSWH, NMWP, NEMOTAUX
    REAL(KIND=JWRO), INTENT(INOUT), DIMENSION(KIJL) :: NEMOTAUY, NEMOWSWAVE, NEMOPHIF
    LOGICAL, INTENT(IN) :: LNUPD
    
    
    INTEGER(KIND=JWIM) :: IJ, K, M, MAXIJ
    
    !     FICTITIOUS VALUE OF THE NORMALISED WAVE ENERGY FLUX UNDER THE SEA ICE
    !     (negative because it is defined as leaving the waves)
    REAL(KIND=JWRB), PARAMETER :: PHIOC_ICE = -3.75_JWRB
    REAL(KIND=JWRB), PARAMETER :: PHIAW_ICE = 3.75_JWRB
    
    !     USE HERSBACH 2011 FOR CD(U10) (SEE ALSO EDSON et al. 2013)
    REAL(KIND=JWRB), PARAMETER :: C1 = 1.03E-3_JWRB
    REAL(KIND=JWRB), PARAMETER :: C2 = 0.04E-3_JWRB
    REAL(KIND=JWRB), PARAMETER :: P1 = 1.48_JWRB
    REAL(KIND=JWRB), PARAMETER :: P2 = -0.21_JWRB
    REAL(KIND=JWRB), PARAMETER :: CDMAX = 0.003_JWRB
    
    REAL(KIND=JWRB), PARAMETER :: EFD_MIN = 0.0625_JWRB    ! corresponds to min Hs=1m under sea ice
    REAL(KIND=JWRB), PARAMETER :: EFD_MAX = 6.25_JWRB    ! corresponds to max Hs=10m under sea ice
    
    REAL(KIND=JWRB) :: TAU, XN, TAUO
    REAL(KIND=JWRB) :: U10P, CD_BULK, CD_WAVE, CD_ICE
    REAL(KIND=JWRB) :: CNST
    REAL(KIND=JWRB) :: EPSUS3
    REAL(KIND=JWRB) :: CITHRSH_INV
    REAL(KIND=JWRB) :: EFD, FFD, EFD_FAC, FFD_FAC
    
    REAL(KIND=JWRB) :: XSTRESS, YSTRESS
    REAL(KIND=JWRB) :: USTAR
    REAL(KIND=JWRB) :: PHILF
    REAL(KIND=JWRB) :: OOVAL
    REAL(KIND=JWRB) :: EM_OC, F1_OC
    REAL(KIND=JWRB) :: CMRHOWGDFTH
    REAL(KIND=JWRB) :: SUMT, SUMX, SUMY
!$acc routine vector
!$acc data present( MIJ, RHOWGDFTH, CINV, SSURF, CICOVER, PHIWA, EM, F1, WSWAVE, WDWAVE, UFRIC, AIRD, NPHIEPS, NTAUOC, NSWH,  &
!$acc & NMWP, NEMOTAUX, NEMOTAUY, NEMOWSWAVE, NEMOPHIF, TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUOC, PHIOCD, PHIEPS, PHIAW )
    
    ! ----------------------------------------------------------------------
    
    
    EPSUS3 = EPSUS*SQRT(EPSUS)
    
    CITHRSH_INV = 1._JWRB / MAX(CITHRSH, 0.01_JWRB)
    
    EFD_FAC = 4.0_JWRB*EGRCRV / G**2
    FFD_FAC = (EGRCRV / AFCRV)**(1.0_JWRB / BFCRV)*G
    
!$loki separator
    MAXIJ = -1
    
!$acc loop vector reduction( max: maxij )
    DO IJ=KIJS,KIJL
      MAXIJ = MAX(MAXIJ, MIJ(IJ))
    END DO
!$acc loop vector
    DO IJ=KIJS,KIJL
      
      !*    DETERMINE NORMALIZED FLUXES FROM AIR TO WAVE AND FROM WAVE TO OCEAN.
      !     -------------------------------------------------------------------
      
      !     ENERGY FLUX from SSURF
      !     MOMENTUM FLUX FROM SSURF
      PHILF = 0.0_JWRB
      XSTRESS = 0.0_JWRB
      YSTRESS = 0.0_JWRB
      
      !     THE INTEGRATION ONLY UP TO FR=MIJ
!$acc loop seq
      DO M=1,MAXIJ
        K = 1
        SUMT = SSURF(IJ, K, M)
        SUMX = SINTH(K)*SSURF(IJ, K, M)
        SUMY = COSTH(K)*SSURF(IJ, K, M)
!$acc loop seq
        DO K=2,NANG
          SUMT = SUMT + SSURF(IJ, K, M)
          SUMX = SUMX + SINTH(K)*SSURF(IJ, K, M)
          SUMY = SUMY + COSTH(K)*SSURF(IJ, K, M)
        END DO
        PHILF = PHILF + SUMT*RHOWGDFTH(IJ, M)
        CMRHOWGDFTH = CINV(IJ, M)*RHOWGDFTH(IJ, M)
        XSTRESS = XSTRESS + SUMX*CMRHOWGDFTH
        YSTRESS = YSTRESS + SUMY*CMRHOWGDFTH
      END DO
      
      IF (LICERUN .and. LWAMRSETCI) THEN
        IF (CICOVER(IJ) > CIBLOCK) THEN
          OOVAL = EXP(-MIN((CICOVER(IJ)*CITHRSH_INV)**4, 10._JWRB))
          !           ADJUST USTAR FOR THE PRESENCE OF SEA ICE
          U10P = MAX(WSWAVE(IJ), EPSU10)
          CD_BULK = MIN((C1 + C2*U10P**P1)*U10P**P2, CDMAX)
          CD_WAVE = (UFRIC(IJ) / U10P)**2
          CD_ICE = OOVAL*CD_WAVE + (1.0_JWRB - OOVAL)*CD_BULK
          USTAR = MAX(SQRT(CD_ICE)*U10P, EPSUS)
          
          ! EM_OC and F1_OC with fully developed model ENERGY
          ! The significant wave height derived from EM_OC will be used
          ! by NEMO as a scaling factor as if it was open ocean
          EFD = MIN(EFD_FAC*USTAR**4, EFD_MAX)
          EM_OC = MAX(OOVAL*EM(IJ) + (1.0_JWRB - OOVAL)*EFD, EFD_MIN)
          FFD = FFD_FAC / USTAR
          F1_OC = OOVAL*F1(IJ) + (1.0_JWRB - OOVAL)*FFD
          F1_OC = MIN(MAX(F1_OC, FR(2)), FR(NFRE))
        ELSE
          OOVAL = 1.0_JWRB
          USTAR = UFRIC(IJ)
          EM_OC = EM(IJ)
          F1_OC = F1(IJ)
        END IF
      ELSE
        OOVAL = 1.0_JWRB
        USTAR = UFRIC(IJ)
        EM_OC = EM(IJ)
        F1_OC = F1(IJ)
      END IF
      
      
      TAU = AIRD(IJ)*MAX(USTAR**2, EPSUS)
      TAUXD(IJ) = TAU*SIN(WDWAVE(IJ))
      TAUYD(IJ) = TAU*COS(WDWAVE(IJ))
      
      TAUOCXD(IJ) = TAUXD(IJ) - OOVAL*XSTRESS
      TAUOCYD(IJ) = TAUYD(IJ) - OOVAL*YSTRESS
      TAUO = SQRT(TAUOCXD(IJ)**2 + TAUOCYD(IJ)**2)
      TAUOC(IJ) = MIN(MAX(TAUO / TAU, TAUOCMIN), TAUOCMAX)
      
      XN = AIRD(IJ)*MAX(USTAR**3, EPSUS3)
      PHIOCD(IJ) = OOVAL*(PHILF - PHIWA(IJ)) + (1.0_JWRB - OOVAL)*PHIOC_ICE*XN
      
      PHIEPS(IJ) = PHIOCD(IJ) / XN
      PHIEPS(IJ) = MIN(MAX(PHIEPS(IJ), PHIEPSMIN), PHIEPSMAX)
      
      PHIOCD(IJ) = PHIEPS(IJ)*XN
      
      PHIAW(IJ) = PHIWA(IJ) / XN
      PHIAW(IJ) = OOVAL*PHIWA(IJ) / XN + (1.0_JWRB - OOVAL)*PHIAW_ICE
      
      IF (LWNEMOCOU .and. LNUPD) THEN
        NPHIEPS(IJ) = PHIEPS(IJ)
        NTAUOC(IJ) = TAUOC(IJ)
        IF (EM_OC /= 0.0_JWRB) THEN
          NSWH(IJ) = 4.0_JWRO*SQRT(EM_OC)
        ELSE
          NSWH(IJ) = 0.0_JWRO
        END IF
        IF (F1_OC /= 0.0_JWRB) THEN
          NMWP(IJ) = 1.0_JWRO / F1_OC
        ELSE
          NMWP(IJ) = 0.0_JWRO
        END IF
        
        IF (LWNEMOTAUOC) THEN
          NEMOTAUX(IJ) = NEMOTAUX(IJ) + TAUOCXD(IJ)
          NEMOTAUY(IJ) = NEMOTAUY(IJ) + TAUOCYD(IJ)
        ELSE
          NEMOTAUX(IJ) = NEMOTAUX(IJ) + TAUXD(IJ)
          NEMOTAUY(IJ) = NEMOTAUY(IJ) + TAUYD(IJ)
        END IF
        NEMOWSWAVE(IJ) = NEMOWSWAVE(IJ) + WSWAVE(IJ)
        NEMOPHIF(IJ) = NEMOPHIF(IJ) + PHIOCD(IJ)
      END IF
      
      ! ----------------------------------------------------------------------
      
    END DO
!$acc end data
  END SUBROUTINE WNFLUXES_SCC
END MODULE WNFLUXES_SCC_MOD
