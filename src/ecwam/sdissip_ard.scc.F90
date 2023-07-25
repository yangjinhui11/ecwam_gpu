! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SDISSIP_ARD_SCC_MOD
  CONTAINS
  SUBROUTINE SDISSIP_ARD_SCC (KIJS, KIJL, FL1, FLD, SL, INDEP, WAVNUM, XK2CG, UFRIC, COSWDIF, RAORW)
    ! ----------------------------------------------------------------------
    
    !**** *SDISSIP_ARD* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.
    
    !     LOTFI AOUF       METEO FRANCE 2013
    !     FABRICE ARDHUIN  IFREMER  2013
    
    
    !*    PURPOSE.
    !     --------
    !       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
    !       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
    !       OF DISSIPATION SOURCE FUNCTION.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *SDISSIP_ARD (KIJS, KIJL, FL1, FLD,SL,*
    !                            INDEP, WAVNUM, XK2CG,
    !                            UFRIC, COSWDIF, RAORW)*
    !          *KIJS*   - INDEX OF FIRST GRIDPOINT
    !          *KIJL*   - INDEX OF LAST GRIDPOINT
    !          *FL1*    - SPECTRUM.
    !          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
    !          *SL*     - TOTAL SOURCE FUNCTION ARRAY
    !          *INDEP*  - DEPTH INDEX
    !          *WAVNUM* - WAVE NUMBER
    !          *XK2CG*  - (WAVE NUMBER)**2 * GROUP SPEED
    !          *UFRIC*  - FRICTION VELOCITY IN M/S.
    !          *RAORW*  - RATIO AIR DENSITY TO WATER DENSITY
    !          *COSWDIF*-  COS(TH(K)-WDWAVE(IJ))
    
    
    !     METHOD.
    !     -------
    
    !       SEE REFERENCES.
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    !       ARDHUIN et AL. JPO DOI:10.1175/20110JPO4324.1
    
    
    ! ----------------------------------------------------------------------
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWFRED, ONLY: FR, TH, ZPIFR
    USE YOWPCONS, ONLY: G, ZPI
    USE YOWPARAM, ONLY: NANG, NFRE, NANG_PARAM
    USE YOWPHYS, ONLY: SDSBR, ISDSDTH, ISB, IPSAT, SSDSC2, SSDSC4, SSDSC6, MICHE, SSDSC3, SSDSBRF1, BRKPBCOEF, SSDSC5, NSDSNTH,  &
    & NDIKCUMUL, INDICESSAT, SATWEIGHTS, CUMULW
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL, NANG, NFRE) :: FLD, SL
    INTEGER(KIND=JWIM), INTENT(IN), DIMENSION(KIJL) :: INDEP
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM, XK2CG
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: UFRIC, RAORW
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: COSWDIF
    
    
    INTEGER(KIND=JWIM) :: IJ, K, M, I, J, M2, K2, KK
    
    REAL(KIND=JWRB) :: TPIINV, TPIINVH, TMP01, TMP03
    REAL(KIND=JWRB) :: EPSR, SSDSC6M1, ZCOEF, ZCOEFM1
    
    
    REAL(KIND=JWRB) :: SSDSC2_SIG
    REAL(KIND=JWRB) :: FACTURB, BTH, BTH0
    REAL(KIND=JWRB), DIMENSION(NANG_PARAM) :: SCUMUL, D
    
    REAL(KIND=JWRB) :: RENEWALFREQ
!$acc routine vector
!$acc data present( FL1, FLD, SL, INDEP, WAVNUM, XK2CG, UFRIC, COSWDIF, RAORW )
    
!$acc loop vector private( SCUMUL, D )
    DO IJ=KIJS,KIJL
      
      ! ----------------------------------------------------------------------
      
      
      ! INITIALISATION
      
      EPSR = SQRT(SDSBR)
      
      TPIINV = 1.0_JWRB / ZPI
      TPIINVH = 0.5_JWRB*TPIINV
      TMP03 = 1.0_JWRB / (SDSBR*MICHE)
      SSDSC6M1 = 1._JWRB - SSDSC6
      
!$acc loop seq
      DO M=1,NFRE
        
        ! SATURATION TERM
        SSDSC2_SIG = SSDSC2*ZPIFR(M)
        ZCOEF = SSDSC2_SIG*SSDSC6
        ZCOEFM1 = SSDSC2_SIG*SSDSC6M1
        
        ! COMPUTE SATURATION SPECTRUM
        BTH0 = 0.0_JWRB
        
!$acc loop seq
        DO K=1,NANG
          BTH = 0.0_JWRB
          ! integrates in directional sector
!$acc loop seq
          DO K2=1,NSDSNTH*2 + 1
            KK = INDICESSAT(K, K2)
            BTH = BTH + SATWEIGHTS(K, K2)*FL1(IJ, KK, M)
          END DO
          BTH = BTH*WAVNUM(IJ, M)*TPIINV*XK2CG(IJ, M)
          BTH0 = MAX(BTH0, BTH)
          
          D(K) = ZCOEFM1*MAX(0._JWRB, BTH*TMP03 - SSDSC4)**IPSAT
          
          SCUMUL(K) = MAX(SQRT(ABS(BTH)) - EPSR, 0._JWRB)**2
        END DO
        
!$acc loop seq
        DO K=1,NANG
          ! cumulative term
          D(K) = D(K) + ZCOEF*MAX(0._JWRB, BTH0*TMP03 - SSDSC4)**IPSAT
          IF (BTH0 <= SDSBR) THEN
            SCUMUL(K) = 0._JWRB
          END IF
          
        END DO
        
        IF (M > NDIKCUMUL) THEN
          ! CUMULATIVE TERM
          IF (SSDSC3 /= 0.0_JWRB) THEN
            
!$acc loop seq
            DO K=1,NANG
              ! Correction of saturation level for shallow-water kinematics
              ! Cumulative effect based on lambda   (breaking probability is
              ! the expected rate of sweeping by larger breaking waves)
              
              RENEWALFREQ = 0.0_JWRB
              
!$acc loop seq
              DO M2=1,M - NDIKCUMUL
!$acc loop seq
                DO K2=1,NANG
                  KK = ABS(K2 - K)
                  IF (KK > NANG / 2)                   KK = KK - NANG / 2
                  ! Integrates over frequencies M2 and directions K2 to
                  ! Integration is performed from M2=1 to a frequency lower than M: IK-NDIKCUMUL
                  RENEWALFREQ = RENEWALFREQ + CUMULW(INDEP(IJ), KK, M2, M)*SCUMUL(K2)
                END DO
              END DO
              
              D(K) = D(K) + RENEWALFREQ
            END DO
          END IF
        END IF
        
        !       WAVE-TURBULENCE INTERACTION TERM
        IF (SSDSC5 /= 0.0_JWRB) THEN
          TMP01 = 2._JWRB*SSDSC5 / G
          FACTURB = TMP01*RAORW(IJ)*UFRIC(IJ)*UFRIC(IJ)
!$acc loop seq
          DO K=1,NANG
            D(K) = D(K) - ZPIFR(M)*WAVNUM(IJ, M)*FACTURB*COSWDIF(IJ, K)
          END DO
        END IF
        
        
        ! ADD ALL CONTRIBUTIONS TO SOURCE TERM
!$acc loop seq
        DO K=1,NANG
          SL(IJ, K, M) = SL(IJ, K, M) + D(K)*FL1(IJ, K, M)
          FLD(IJ, K, M) = FLD(IJ, K, M) + D(K)
        END DO
      END DO
      
      
    END DO
!$acc end data
  END SUBROUTINE SDISSIP_ARD_SCC
END MODULE SDISSIP_ARD_SCC_MOD
