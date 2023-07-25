! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE STRESSO_SCC_MOD
  CONTAINS
  SUBROUTINE STRESSO_SCC (KIJS, KIJL, MIJ, RHOWGDFTH, FL1, SL, SPOS, CINV, WDWAVE, UFRIC, Z0M, AIRD, RNFAC, COSWDIF, SINWDIF2,  &
  & TAUW, TAUWDIR, PHIWA, LLPHIWA)
    
    ! ----------------------------------------------------------------------
    
    !**** *STRESSO* - COMPUTATION OF WAVE STRESS.
    
    !     H. GUNTHER      GKSS/ECMWF NOVEMBER   1989 CODE MOVED FROM SINPUT.
    !     P.A.E.M. JANSSEN     KNMI  AUGUST     1990
    !     J. BIDLOT            ECMWF FEBRUARY   1996-97
    !     S. ABDALLA           ECMWF OCTOBER    2001 INTRODUCTION OF VARIABLE
    !                                                AIR DENSITY
    !     P.A.E.M. JANSSEN     ECMWF            2011  ADD FLUX CALULATIONS
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTE NORMALIZED WAVE STRESS FROM INPUT SOURCE FUNCTION
    
    !**   INTERFACE.
    !     ----------
    
    !        *CALL* *STRESSO (KIJS, KIJL, MIJ, RHOWGDFTH,
    !                         FL1, SL, SPOS,
    !    &                    CINV,
    !    &                    WDWAVE, UFRIC, Z0M, AIRD, RNFAC,
    !    &                    COSWDIF, SINWDIF2,
    !    &                    TAUW, TAUWDIR, PHIWA)*
    !         *KIJS*        - INDEX OF FIRST GRIDPOINT.
    !         *KIJL*        - INDEX OF LAST GRIDPOINT.
    !         *MIJ*         - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
    !         *RHOWGDFTH    - WATER DENSITY * G * DF * DTHETA
    !         *FL1*         - WAVE SPECTRUM.
    !         *SL*          - WIND INPUT SOURCE FUNCTION ARRAY (positive and negative contributions).
    !         *SPOS*        - POSITIVE WIND INPUT SOURCE FUNCTION ARRAY.
    !         *CINV*        - INVERSE PHASE VELOCITY.
    !         *WDWAVE*      - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
    !                         NOTATION (POINTING ANGLE OF WIND VECTOR,
    !                         CLOCKWISE FROM NORTH).
    !         *UFRIC*       - FRICTION VELOCITY IN M/S.
    !         *Z0M*         - ROUGHNESS LENGTH IN M.
    !         *AIRD*        - AIR DENSITY IN KG/M**3.
    !         *RNFAC*       - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
    !         *COSWDIF*     - COS(TH(K)-WDWAVE(IJ))
    !         *SINWDIF2*    - SIN(TH(K)-WDWAVE(IJ))**2
    !         *TAUW*        - KINEMATIC WAVE STRESS IN (M/S)**2
    !         *TAUWDIR*     - KINEMATIC WAVE STRESS DIRECTION
    !         *PHIWA*       - ENERGY FLUX FROM WIND INTO WAVES INTEGRATED
    !                         OVER THE FULL FREQUENCY RANGE.
    !         *LLPHIWA*     - TRUE IF PHIWA NEEDS TO BE COMPUTED
    
    !     METHOD.
    !     -------
    
    !       THE INPUT SOURCE FUNCTION IS INTEGRATED OVER FREQUENCY
    !       AND DIRECTIONS.
    !       BECAUSE ARRAY *SPOS* IS USED, ONLY THE INPUT SOURCE
    !       HAS TO BE STORED IN *SPOS* (CALL FIRST SINPUT, THEN
    !       STRESSO, AND THEN THE REST OF THE SOURCE FUNCTIONS)
    
    !     REFERENCE.
    !     ----------
    !       P. JANSSEN,
    
    ! ----------------------------------------------------------------------
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWCOUP, ONLY: LLGCBZ0
    USE YOWFRED, ONLY: FR, RHOWG_DFIM, DELTH, TH, COSTH, SINTH, FR5
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPHYS, ONLY: TAUWSHELTER
    USE YOWTABL, ONLY: EPS1
    USE YOWSTAT, ONLY: IPHYS
    
    USE TAU_PHI_HF_SCC_MOD, ONLY: TAU_PHI_HF_SCC
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    INTEGER(KIND=JWIM), INTENT(IN) :: MIJ(KIJL)
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: RHOWGDFTH
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1, SL, SPOS
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: CINV
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: WDWAVE, UFRIC, Z0M, AIRD, RNFAC
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: COSWDIF, SINWDIF2
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: TAUW, TAUWDIR, PHIWA
    LOGICAL, INTENT(IN) :: LLPHIWA
    
    
    INTEGER(KIND=JWIM) :: IJ, M, K, I, J, II, MAXIJ
    
    REAL(KIND=JWRB) :: TAUTOUS2
    REAL(KIND=JWRB) :: COSW, FCOSW2
    REAL(KIND=JWRB), DIMENSION(KIJL) :: XSTRESS, YSTRESS
    REAL(KIND=JWRB), DIMENSION(KIJL) :: TAUHF, PHIHF
    REAL(KIND=JWRB), DIMENSION(KIJL) :: USDIRP, UST
    
    REAL(KIND=JWRB) :: CMRHOWGDFTH
    REAL(KIND=JWRB) :: TAUX, TAUY, TAUPX, TAUPY
    REAL(KIND=JWRB) :: SUMT, SUMX, SUMY
    
    LOGICAL :: LTAUWSHELTER
!$acc routine vector
!$acc data present( MIJ, RHOWGDFTH, FL1, SL, SPOS, CINV, WDWAVE, UFRIC, Z0M, AIRD, RNFAC, COSWDIF, SINWDIF2, TAUW, TAUWDIR,  &
!$acc & PHIWA )
    
    ! ----------------------------------------------------------------------
    
    
    MAXIJ = -1
    
!$acc loop vector reduction( max: MAXIJ )
    DO IJ=KIJS,KIJL
      MAXIJ = MAX(MAXIJ, MIJ(IJ))
    END DO
!$acc loop vector
    DO IJ=KIJS,KIJL
      
      PHIWA(IJ) = 0.0_JWRB
      XSTRESS(IJ) = 0.0_JWRB
      YSTRESS(IJ) = 0.0_JWRB
      
      !*    CONTRIBUTION TO THE WAVE STRESS FROM THE NEGATIVE PART OF THE WIND INPUT
      !     ------------------------------------------------------------------------
      
      IF (LLPHIWA) THEN
        !     full energy flux due to negative Sinput (SL-SPOS)
        !     we assume that above NFRE, the contibutions can be neglected
!$acc loop seq
        DO M=1,NFRE
!$acc loop seq
          DO K=1,NANG
            PHIWA(IJ) = PHIWA(IJ) + (SL(IJ, K, M) - SPOS(IJ, K, M))*RHOWG_DFIM(M)
          END DO
        END DO
      END IF
      
      !*    CALCULATE LOW-FREQUENCY CONTRIBUTION TO STRESS AND ENERGY FLUX (positive sinput).
      !     ---------------------------------------------------------------------------------
!$acc loop seq
      DO M=1,MAXIJ
        !     THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
        K = 1
        SUMX = SPOS(IJ, K, M)*SINTH(K)
        SUMY = SPOS(IJ, K, M)*COSTH(K)
!$acc loop seq
        DO K=2,NANG
          SUMX = SUMX + SPOS(IJ, K, M)*SINTH(K)
          SUMY = SUMY + SPOS(IJ, K, M)*COSTH(K)
        END DO
        CMRHOWGDFTH = RHOWGDFTH(IJ, M)*CINV(IJ, M)
        XSTRESS(IJ) = XSTRESS(IJ) + CMRHOWGDFTH*SUMX
        YSTRESS(IJ) = YSTRESS(IJ) + CMRHOWGDFTH*SUMY
      END DO
      
      !     TAUW is the kinematic wave stress !
      XSTRESS(IJ) = XSTRESS(IJ) / MAX(AIRD(IJ), 1.0_JWRB)
      YSTRESS(IJ) = YSTRESS(IJ) / MAX(AIRD(IJ), 1.0_JWRB)
      
      IF (LLPHIWA) THEN
!$acc loop seq
        DO M=1,MAXIJ
          !       THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
          K = 1
          SUMT = SPOS(IJ, K, M)
!$acc loop seq
          DO K=2,NANG
            SUMT = SUMT + SPOS(IJ, K, M)
          END DO
          PHIWA(IJ) = PHIWA(IJ) + RHOWGDFTH(IJ, M)*SUMT
        END DO
      END IF
      
      !*    CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS and energy flux (positive sinput).
      !     ----------------------------------------------------------------------------------
      
      IF (IPHYS == 0 .or. TAUWSHELTER == 0.0_JWRB) THEN
        LTAUWSHELTER = .false.
        USDIRP(IJ) = WDWAVE(IJ)
        UST(IJ) = UFRIC(IJ)
      ELSE
        LTAUWSHELTER = .true.
        TAUX = UFRIC(IJ)**2*SIN(WDWAVE(IJ))
        TAUY = UFRIC(IJ)**2*COS(WDWAVE(IJ))
        TAUPX = TAUX - TAUWSHELTER*XSTRESS(IJ)
        TAUPY = TAUY - TAUWSHELTER*YSTRESS(IJ)
        USDIRP(IJ) = ATAN2(TAUPX, TAUPY)
        UST(IJ) = (TAUPX**2 + TAUPY**2)**0.25_JWRB
      END IF
      
    END DO
    CALL TAU_PHI_HF_SCC(KIJS, KIJL, MIJ, LTAUWSHELTER, UFRIC, Z0M, FL1, AIRD, RNFAC, COSWDIF, SINWDIF2, UST, TAUHF, PHIHF,  &
    & LLPHIWA)
!$acc loop vector
    DO IJ=KIJS,KIJL
      
      XSTRESS(IJ) = XSTRESS(IJ) + TAUHF(IJ)*SIN(USDIRP(IJ))
      YSTRESS(IJ) = YSTRESS(IJ) + TAUHF(IJ)*COS(USDIRP(IJ))
      TAUW(IJ) = SQRT(XSTRESS(IJ)**2 + YSTRESS(IJ)**2)
      TAUW(IJ) = MAX(TAUW(IJ), 0.0_JWRB)
      TAUWDIR(IJ) = ATAN2(XSTRESS(IJ), YSTRESS(IJ))
      
      IF (.not.LLGCBZ0) THEN
        TAUTOUS2 = 1.0_JWRB / (1.0_JWRB + EPS1)
        TAUW(IJ) = MIN(TAUW(IJ), UFRIC(IJ)**2*TAUTOUS2)
      END IF
      
      IF (LLPHIWA) THEN
        PHIWA(IJ) = PHIWA(IJ) + PHIHF(IJ)
      END IF
      
      
    END DO
!$acc end data
  END SUBROUTINE STRESSO_SCC
END MODULE STRESSO_SCC_MOD
