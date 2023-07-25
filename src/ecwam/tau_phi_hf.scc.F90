! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE TAU_PHI_HF_SCC_MOD
  !CONTAINED SUBROUTINES:
  ! - OMEGAGC
  ! - TAU_PHI_HF
  ! - MEANSQS_GC
  CONTAINS
  SUBROUTINE OMEGAGC_SCC (UST, NS, XKS, OMS)
    
    !***  DETERMINE THE CUT-OFF ANGULAR FREQUENCY FOR THE GRAV-CAPILLARY WAVES
    !     !!!! rounded to the closest index of XK_GC  !!!!!
    
    !     AUTHOR: PETER JANSSEN
    !     ------
    
    !     REFERENCES:
    !     ----------
    
    !     VIERS PAPER EQ.(29)
    
    !----------------------------------------------------------------------
    
    USE NS_GC_SCC_MOD, ONLY: NS_GC_SCC
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: OMEGA_GC, XK_GC
    
    
    !----------------------------------------------------------------------
    
    IMPLICIT NONE
    REAL(KIND=JWRB), INTENT(IN) :: UST
    INTEGER(KIND=JWIM), INTENT(OUT) :: NS    ! index in array XK_GC corresponding to XKS and OMS
    REAL(KIND=JWRB), INTENT(OUT) :: XKS    ! cut-off wave number
    REAL(KIND=JWRB), INTENT(OUT) :: OMS    ! cut-off angular frequency
    
    
    
    
    ! ----------------------------------------------------------------------
!$acc routine seq
    
    
    NS = NS_GC_SCC(UST)
    XKS = XK_GC(NS)
    OMS = OMEGA_GC(NS)
    
    
  END SUBROUTINE OMEGAGC_SCC
  
  SUBROUTINE TAU_PHI_HF_SCC (KIJS, KIJL, MIJ, LTAUWSHELTER, UFRIC, Z0M, FL1, AIRD, RNFAC, COSWDIF, SINWDIF2, UST, TAUHF, PHIHF,  &
  & LLPHIHF)
    
    ! ----------------------------------------------------------------------
    
    !**** *TAU_PHI_HF* - COMPUTATION OF HIGH-FREQUENCY STRESS.
    !                                   HIGH-FREQUENCY ENERGY FLUX.
    
    !     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90
    !     JEAN BIDLOT  ECMWF  JANUARY 2017
    
    !*    PURPOSE.
    !     ---------
    
    !       COMPUTE HIGH-FREQUENCY WAVE STRESS AND ENERGY FLUX
    
    !**   INTERFACE.
    !     ---------
    
    !       *CALL* *TAU_PHI_HF(KIJS, KIJL, MIJ, LTAUWSHELTER, UFRIC, UST, Z0M,
    !                          FL1, AIRD, RNFAC,
    !                          COSWDIF, SINWDIF2,
    !                          UST, TAUHF, PHIHF, LLPHIHF)
    !          *KIJS*         - INDEX OF FIRST GRIDPOINT
    !          *KIJL*         - INDEX OF LAST GRIDPOINT
    !          *MIJ*          - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
    !          *LTAUWSHELTER* - if true then TAUWSHELTER
    !          *FL1*          - WAVE SPECTRUM.
    !          *AIRD*         - AIR DENSITY IN KG/M**3.
    !          *RNFAC*        - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
    !          *UFRIC*        - FRICTION VELOCITY
    !          *COSWDIF*      - COS(TH(K)-WDWAVE(IJ))
    !          *SINWDIF2*     - SIN(TH(K)-WDWAVE(IJ))**2
    !          *UST*          - REDUCED FRICTION VELOCITY DUE TO SHELTERING
    !          *Z0M*          - ROUGHNESS LENGTH
    !          *TAUHF*        - HIGH-FREQUENCY STRESS
    !          *PHIHF*        - HIGH-FREQUENCY ENERGY FLUX INTO OCEAN
    !          *LLPHIHF*      - TRUE IF PHIHF NEEDS TO COMPUTED
    
    
    !     METHOD.
    !     -------
    
    !       IT NEEDS A CALL TO INIT_X0TAUHF TO INITIALISE
    !       SEE REFERENCE FOR WAVE STRESS CALCULATION.
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    !       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
    
    ! ----------------------------------------------------------------------
    
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWCOUP, ONLY: X0TAUHF, JTOT_TAUHF, WTAUHF, LLGCBZ0, LLNORMAGAM
    USE YOWFRED, ONLY: ZPIFR, FR5, TH, DELTH
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: G, GM1, ZPI, ZPI4GM1, ZPI4GM2
    USE YOWPHYS, ONLY: ZALP, XKAPPA, TAUWSHELTER, GAMNCONST
    USE YOWTEST, ONLY: IU06
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    INTEGER(KIND=JWIM), INTENT(IN) :: MIJ(KIJL)
    LOGICAL, INTENT(IN) :: LTAUWSHELTER
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: UFRIC, Z0M
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: AIRD, RNFAC
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: COSWDIF, SINWDIF2
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: UST
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: TAUHF, PHIHF
    LOGICAL, INTENT(IN) :: LLPHIHF
    
    
    INTEGER(KIND=JWIM) :: J, IJ, K
    INTEGER(KIND=JWIM) :: NS
    
    REAL(KIND=JWRB), PARAMETER :: ZSUPMAX = 0.0_JWRB    !  LOG(1.)
    REAL(KIND=JWRB) :: OMEGA, OMEGACC
    REAL(KIND=JWRB) :: X0G
    REAL(KIND=JWRB) :: YC, Y, CM1, ZX, ZARG, ZLOG, ZBETA
    REAL(KIND=JWRB) :: FNC, FNC2
    REAL(KIND=JWRB) :: GAMNORMA    ! RENORMALISATION FACTOR OF THE GROWTH RATE
    REAL(KIND=JWRB) :: ZNZ, CONFG
    REAL(KIND=JWRB) :: COSW, FCOSW2
    
    REAL(KIND=JWRB) :: XKS, OMS
    REAL(KIND=JWRB) :: SQRTZ0OG, ZSUP, ZINF, DELZ
    REAL(KIND=JWRB) :: TAUL, XLOGGZ0, SQRTGZ0
    REAL(KIND=JWRB) :: USTPH
    REAL(KIND=JWRB) :: CONST1, CONST2, CONSTTAU, CONSTPHI
    REAL(KIND=JWRB) :: F1DCOS2, F1DCOS3
    REAL(KIND=JWRB) :: F1D, F1DSIN2
!$acc routine vector
!$acc data present( MIJ, UFRIC, Z0M, FL1, AIRD, RNFAC, COSWDIF, SINWDIF2, UST, TAUHF, PHIHF )
    
!$acc loop vector
    DO IJ=KIJS,KIJL
      
      ! ----------------------------------------------------------------------
      
      
      IF (LLGCBZ0) THEN
!$loki inline
        CALL OMEGAGC_SCC(UFRIC(IJ), NS, XKS, OMS)
      END IF
      
      !     See INIT_X0TAUHF
      X0G = X0TAUHF*G
      
      IF (LLPHIHF)       USTPH = UST(IJ)
      
      !*    COMPUTE THE INTEGRALS
      !     ---------------------
      
      XLOGGZ0 = LOG(G*Z0M(IJ))
      OMEGACC = MAX(ZPIFR(MIJ(IJ)), X0G / UST(IJ))
      SQRTZ0OG = SQRT(Z0M(IJ)*GM1)
      SQRTGZ0 = 1.0_JWRB / SQRTZ0OG
      YC = OMEGACC*SQRTZ0OG
      ZINF = LOG(YC)
      
      CONSTTAU = ZPI4GM2*FR5(MIJ(IJ))
      
      K = 1
      COSW = MAX(COSWDIF(IJ, K), 0.0_JWRB)
      FCOSW2 = FL1(IJ, K, MIJ(IJ))*COSW**2
      F1DCOS3 = FCOSW2*COSW
      F1DCOS2 = FCOSW2
      F1DSIN2 = FL1(IJ, K, MIJ(IJ))*SINWDIF2(IJ, K)
      F1D = FL1(IJ, K, MIJ(IJ))
!$acc loop seq
      DO K=2,NANG
        COSW = MAX(COSWDIF(IJ, K), 0.0_JWRB)
        FCOSW2 = FL1(IJ, K, MIJ(IJ))*COSW**2
        F1DCOS3 = F1DCOS3 + FCOSW2*COSW
        F1DCOS2 = F1DCOS2 + FCOSW2
        F1DSIN2 = F1DSIN2 + FL1(IJ, K, MIJ(IJ))*SINWDIF2(IJ, K)
        F1D = F1D + FL1(IJ, K, MIJ(IJ))
      END DO
      F1DCOS3 = DELTH*F1DCOS3
      F1DCOS2 = DELTH*F1DCOS2
      F1DSIN2 = DELTH*F1DSIN2
      F1D = DELTH*F1D
      
      IF (LLNORMAGAM) THEN
        CONFG = GAMNCONST*FR5(MIJ(IJ))*RNFAC(IJ)*SQRTGZ0
        CONST1 = CONFG*F1DSIN2
        CONST2 = CONFG*F1D
      ELSE
        CONST1 = 0.0_JWRB
        CONST2 = 0.0_JWRB
      END IF
      
      
      !     TAUHF :
      IF (LLGCBZ0) THEN
        ZSUP = MIN(LOG(OMS*SQRTZ0OG), ZSUPMAX)
      ELSE
        ZSUP = ZSUPMAX
      END IF
      
      TAUL = UST(IJ)**2
      DELZ = MAX((ZSUP - ZINF) / REAL(JTOT_TAUHF - 1, kind=JWRB), 0.0_JWRB)
      TAUHF(IJ) = 0.0_JWRB
      
      ! Intergrals are integrated following a change of variable : Z=LOG(Y)
      IF (LTAUWSHELTER) THEN
!$acc loop seq
        DO J=1,JTOT_TAUHF
          Y = EXP(ZINF + REAL(J - 1, kind=JWRB)*DELZ)
          OMEGA = Y*SQRTGZ0
          CM1 = OMEGA*GM1
          ZX = UST(IJ)*CM1 + ZALP
          ZARG = XKAPPA / ZX
          ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
          ZLOG = MIN(ZLOG, 0.0_JWRB)
          ZBETA = ZLOG**4*EXP(ZLOG)
          ZNZ = ZBETA*UST(IJ)*Y
          GAMNORMA = (1.0_JWRB + CONST1*ZNZ) / (1.0_JWRB + CONST2*ZNZ)
          FNC2 = F1DCOS3*CONSTTAU*ZBETA*TAUL*WTAUHF(J)*DELZ*GAMNORMA
          TAUL = MAX(TAUL - TAUWSHELTER*FNC2, 0.0_JWRB)
          
          UST(IJ) = SQRT(TAUL)
          TAUHF(IJ) = TAUHF(IJ) + FNC2
        END DO
      ELSE
!$acc loop seq
        DO J=1,JTOT_TAUHF
          Y = EXP(ZINF + REAL(J - 1, kind=JWRB)*DELZ)
          OMEGA = Y*SQRTGZ0
          CM1 = OMEGA*GM1
          ZX = UST(IJ)*CM1 + ZALP
          ZARG = XKAPPA / ZX
          ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
          ZLOG = MIN(ZLOG, 0.0_JWRB)
          ZBETA = ZLOG**4*EXP(ZLOG)
          FNC2 = ZBETA*WTAUHF(J)
          ZNZ = ZBETA*UST(IJ)*Y
          GAMNORMA = (1.0_JWRB + CONST1*ZNZ) / (1.0_JWRB + CONST2*ZNZ)
          TAUHF(IJ) = TAUHF(IJ) + FNC2*GAMNORMA
        END DO
        TAUHF(IJ) = F1DCOS3*CONSTTAU*TAUL*TAUHF(IJ)*DELZ
      END IF
      
      
      PHIHF(IJ) = 0.0_JWRB
      IF (LLPHIHF) THEN
        !       PHIHF:
        !       We are neglecting the gravity-capillary contribution
        !       Recompute DELZ over the full interval
        TAUL = USTPH**2
        ZSUP = ZSUPMAX
        DELZ = MAX((ZSUP - ZINF) / REAL(JTOT_TAUHF - 1, kind=JWRB), 0.0_JWRB)
        
        CONSTPHI = AIRD(IJ)*ZPI4GM1*FR5(MIJ(IJ))
        
        ! Intergrals are integrated following a change of variable : Z=LOG(Y)
        IF (LTAUWSHELTER) THEN
!$acc loop seq
          DO J=1,JTOT_TAUHF
            Y = EXP(ZINF + REAL(J - 1, kind=JWRB)*DELZ)
            OMEGA = Y*SQRTGZ0
            CM1 = OMEGA*GM1
            ZX = USTPH*CM1 + ZALP
            ZARG = XKAPPA / ZX
            ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
            ZLOG = MIN(ZLOG, 0.0_JWRB)
            ZBETA = ZLOG**4*EXP(ZLOG)
            ZNZ = ZBETA*UST(IJ)*Y
            GAMNORMA = (1.0_JWRB + CONST1*ZNZ) / (1.0_JWRB + CONST2*ZNZ)
            FNC2 = ZBETA*TAUL*WTAUHF(J)*DELZ*GAMNORMA
            TAUL = MAX(TAUL - TAUWSHELTER*F1DCOS3*CONSTTAU*FNC2, 0.0_JWRB)
            USTPH = SQRT(TAUL)
            PHIHF(IJ) = PHIHF(IJ) + FNC2 / Y
          END DO
          PHIHF(IJ) = F1DCOS2*CONSTPHI*SQRTZ0OG*PHIHF(IJ)
        ELSE
!$acc loop seq
          DO J=1,JTOT_TAUHF
            Y = EXP(ZINF + REAL(J - 1, kind=JWRB)*DELZ)
            OMEGA = Y*SQRTGZ0
            CM1 = OMEGA*GM1
            ZX = USTPH*CM1 + ZALP
            ZARG = XKAPPA / ZX
            ZLOG = XLOGGZ0 + 2.0_JWRB*LOG(CM1) + ZARG
            ZLOG = MIN(ZLOG, 0.0_JWRB)
            ZBETA = ZLOG**4*EXP(ZLOG)
            ZNZ = ZBETA*UST(IJ)*Y
            GAMNORMA = (1.0_JWRB + CONST1*ZNZ) / (1.0_JWRB + CONST2*ZNZ)
            FNC2 = ZBETA*WTAUHF(J)*GAMNORMA
            PHIHF(IJ) = PHIHF(IJ) + FNC2 / Y
          END DO
          PHIHF(IJ) = F1DCOS2*CONSTPHI*SQRTZ0OG*TAUL*PHIHF(IJ)*DELZ
        END IF
      END IF
      
      
    END DO
!$acc end data
  END SUBROUTINE TAU_PHI_HF_SCC
  
  SUBROUTINE MEANSQS_GC (XKMSS, KIJS, KIJL, HALP, USTAR, XMSSCG, FRGC)
    
    !***  DETERMINE MSS FOR GRAV-CAP WAVES UP TO WAVE NUMBER XKMSS
    
    !     AUTHOR: PETER JANSSEN
    !     ------
    
    !     REFERENCES:
    !     ----------
    
    !     VIERS PAPER EQ.(29)
    
    !----------------------------------------------------------------------
    
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: NWAV_GC, XLOGKRATIOM1_GC, XKM_GC, VG_GC, C2OSQRTVG_GC, DELKCC_GC, DELKCC_GC_NS
    USE YOWPCONS, ONLY: G, ZPI, SURFT
    
    
    !----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    
    REAL(KIND=JWRB), INTENT(IN) :: XKMSS    ! WAVE NUMBER CUT-OFF
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: HALP    ! 1/2 Phillips parameter
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: USTAR    ! friction velocity
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: XMSSCG    ! mean square slope for gravity-capillary waves
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: FRGC    ! Frequency from which the gravity-capillary spectrum is approximated
    
    
    INTEGER(KIND=JWIM) :: IJ, I, NE
    INTEGER(KIND=JWIM), DIMENSION(KIJL) :: NS
    REAL(KIND=JWRB), DIMENSION(KIJL) :: XKS, OMS, COEF
    
    !     INCLUDE FUNCTIONS FROM GRAVITY-CAPILLARY DISPERSION REALTIONS
#include "gc_dispersion.h"
    
    ! ----------------------------------------------------------------------
    
    
    NE = MIN(MAX(NINT(LOG(XKMSS*XKM_GC(1))*XLOGKRATIOM1_GC), 1), NWAV_GC)
    
    DO IJ=KIJS,KIJL
      CALL OMEGAGC(USTAR(IJ), NS(IJ), XKS(IJ), OMS(IJ))
      FRGC(IJ) = OMS(IJ) / ZPI
      IF (XKS(IJ) > XKMSS) THEN
        NS(IJ) = NE
        XMSSCG(IJ) = 0.0_JWRB
      ELSE
        XMSSCG(IJ) = DELKCC_GC_NS(NS(IJ))*XKM_GC(NS(IJ))
      END IF
    END DO
    
    DO IJ=KIJS,KIJL
      DO I=NS(IJ) + 1,NE
        !         ANALYTICAL FORM INERTIAL SUB RANGE F(k) = k**(-4)*BB
        !         BB = COEF(IJ)*SQRT(VG_GC(I))/C_GC(I)**2
        !         mss :  integral of k**2 F(k)  k dk
        XMSSCG(IJ) = XMSSCG(IJ) + DELKCC_GC(I)*XKM_GC(I)
      END DO
      COEF(IJ) = C2OSQRTVG_GC(NS(IJ))*HALP(IJ)
      XMSSCG(IJ) = XMSSCG(IJ)*COEF(IJ)
    END DO
    
    
  END SUBROUTINE MEANSQS_GC
END MODULE TAU_PHI_HF_SCC_MOD
