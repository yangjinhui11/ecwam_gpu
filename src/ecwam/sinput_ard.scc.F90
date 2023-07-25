! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE SINPUT_ARD_SCC_MOD
  !CONTAINED SUBROUTINES:
  ! - WSIGSTAR
  ! - SINPUT_ARD
  ! - SINPUT_JAN
  CONTAINS
  SUBROUTINE WSIGSTAR_SCC (WSWAVE, UFRIC, Z0M, WSTAR, SIG_N)
    ! ----------------------------------------------------------------------
    
    !**** *WSIGSTAR* - COMPUTATION OF THE RELATIVE STANDARD DEVIATION OF USTAR.
    
    !*    PURPOSE.
    !     ---------
    
    !     COMPUTES THE STANDARD DEVIATION OF USTAR DUE TO SMALL SCALE GUSTINESS
    !     RELATIVE TO USTAR
    
    !**   INTERFACE.
    !     ----------
    
    !     *CALL* *WSIGSTAR (KIJS, KIJL, WSWAVE, UFRIC, Z0M, WSTAR, SIG_N)
    !             *KIJS*   - INDEX OF FIRST GRIDPOINT.
    !             *KIJL*   - INDEX OF LAST GRIDPOINT.
    !             *WSWAVE* - 10M WIND SPEED (m/s).
    !             *UFRIC*  - NEW FRICTION VELOCITY IN M/S.
    !             *Z0M*    - ROUGHNESS LENGTH IN M.
    !             *WSTAR*  - FREE CONVECTION VELOCITY SCALE (M/S).
    !             *SIG_N*  - ESTINATED RELATIVE STANDARD DEVIATION OF USTAR.
    
    !     METHOD.
    !     -------
    
    !     USE PANOFSKY (1991) TO EXPRESS THE STANDARD DEVIATION OF U10 IN TERMS
    !     USTAR AND  w* THE CONVECTIVE VELOCITY SCALE.
    !     (but with the background gustiness set to 0.)
    !     and USTAR=SQRT(Cd)*U10 to DERIVE THE STANDARD DEVIATION OF USTAR.
    !     WITH CD=A+B*U10 (see below).
    
    !     REFERENCE.
    !     ----------
    
    !     SEE SECTION 3.2.1 OF THE WAM DOCUMENTATION.
    
    ! ----------------------------------------------------------------------
    
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWCOUP, ONLY: LLGCBZ0
    USE YOWPCONS, ONLY: G, EPSUS, ACDLIN, BCDLIN
    USE YOWPHYS, ONLY: XKAPPA, RNUM, ALPHAMIN, ALPHAMAX
    USE YOWWIND, ONLY: WSPMIN
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    REAL(KIND=JWRB), INTENT(IN) :: WSWAVE, UFRIC, Z0M, WSTAR
    REAL(KIND=JWRB), INTENT(OUT) :: SIG_N
    
    REAL(KIND=JWRB), PARAMETER :: BG_GUST = 0.0_JWRB    ! NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)
    REAL(KIND=JWRB), PARAMETER :: ONETHIRD = 1.0_JWRB / 3.0_JWRB
    REAL(KIND=JWRB), PARAMETER :: SIG_NMAX = 0.9_JWRB    ! MAX OF RELATIVE STANDARD DEVIATION OF USTAR
    
    REAL(KIND=JWRB), PARAMETER :: LOG10 = LOG(10.0_JWRB)
    REAL(KIND=JWRB), PARAMETER :: C1 = 1.03E-3_JWRB
    REAL(KIND=JWRB), PARAMETER :: C2 = 0.04E-3_JWRB
    REAL(KIND=JWRB), PARAMETER :: P1 = 1.48_JWRB
    REAL(KIND=JWRB), PARAMETER :: P2 = -0.21_JWRB
    
    REAL(KIND=JWRB) :: ZCHAR, C_D, DC_DDU, SIG_CONV
    REAL(KIND=JWRB) :: XKAPPAD, U10, C2U10P1, U10P2
    REAL(KIND=JWRB) :: BCD, U10M1, ZN, Z0VIS
    
    
    ! ----------------------------------------------------------------------
!$acc routine seq
    
    
    IF (LLGCBZ0) THEN
      ZN = RNUM
      
      U10M1 = 1.0_JWRB / MAX(WSWAVE, WSPMIN)
      ! CHARNOCK:
      Z0VIS = ZN / MAX(UFRIC, EPSUS)
      ZCHAR = G*(Z0M - Z0VIS) / MAX(UFRIC**2, EPSUS)
      ZCHAR = MAX(MIN(ZCHAR, ALPHAMAX), ALPHAMIN)
      
      BCD = BCDLIN*SQRT(ZCHAR)
      C_D = ACDLIN + BCD*WSWAVE
      DC_DDU = BCD
      SIG_CONV = 1.0_JWRB + 0.5_JWRB*WSWAVE / C_D*DC_DDU
      SIG_N = MIN(SIG_NMAX, SIG_CONV*U10M1*(BG_GUST*UFRIC**3 + 0.5_JWRB*XKAPPA*WSTAR**3)**ONETHIRD)
    ELSE
      ZN = 0.0_JWRB
      
      !!! for consistency I have kept the old method, even though the new method above could be used,
      !!! but until LLGCBZ0 is the default, keep the old scheme whe it is not...
      !
      !       IN THE FOLLOWING U10 IS ESTIMATED ASSUMING EVERYTHING IS
      !       BASED ON U*
      !
      XKAPPAD = 1.0_JWRB / XKAPPA
      U10 = UFRIC*XKAPPAD*(LOG10 - LOG(Z0M))
      U10 = MAX(U10, WSPMIN)
      U10M1 = 1.0_JWRB / U10
      C2U10P1 = C2*U10**P1
      U10P2 = U10**P2
      C_D = (C1 + C2U10P1)*U10P2
      DC_DDU = (P2*C1 + (P1 + P2)*C2U10P1)*U10P2*U10M1
      SIG_CONV = 1.0_JWRB + 0.5_JWRB*U10 / C_D*DC_DDU
      SIG_N = MIN(SIG_NMAX, SIG_CONV*U10M1*(BG_GUST*UFRIC**3 + 0.5_JWRB*XKAPPA*WSTAR**3)**ONETHIRD)
    END IF
    
    
  END SUBROUTINE WSIGSTAR_SCC
  SUBROUTINE SINPUT_ARD_SCC (NGST, LLSNEG, KIJS, KIJL, FL1, WAVNUM, CINV, XK2CG, WDWAVE, WSWAVE, UFRIC, Z0M, COSWDIF, SINWDIF2,  &
  & RAORW, WSTAR, RNFAC, FLD, SL, SPOS, XLLWS)
    ! ----------------------------------------------------------------------
    
    !**** *SINPUT_ARD* - COMPUTATION OF INPUT SOURCE FUNCTION.
    
    
    !*    PURPOSE.
    !     ---------
    
    !       COMPUTE THE WIND INPUT SOURCE TRERM BASED ON ARDHUIN ET AL. 2010.
    
    !       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
    !       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
    !       INPUT SOURCE FUNCTION.
    !
    !       GUSTINESS IS INTRODUCED FOLL0WING THE APPROACH OF JANSSEN(1986),
    !       USING A GAUSS-HERMITE APPROXIMATION SUGGESTED BY MILES(1997).
    !       IN THE PRESENT VERSION ONLY TWO HERMITE POLYNOMIALS ARE UTILISED
    !       IN THE EVALUATION OF THE PROBABILITY INTEGRAL. EXPLICITELY ONE THEN
    !       FINDS:
    !
    !             <GAMMA(X)> = 0.5*( GAMMA(X(1+SIG)) + GAMMA(X(1-SIG)) )
    !
    !       WHERE X IS THE FRICTION VELOCITY AND SIG IS THE RELATIVE GUSTINESS
    !       LEVEL.
    
    !**   INTERFACE.
    !     ----------
    
    !     *CALL* *SINPUT_ARD (NGST, LLSNEG, KIJS, KIJL, FL1,
    !    &                    WAVNUM, CINV, XK2CG,
    !    &                    WSWAVE, WDWAVE, UFRIC, Z0M,
    !    &                    COSWDIF, SINWDIF2,
    !    &                    RAORW, WSTAR, RNFAC,
    !    &                    FLD, SL, SPOS, XLLWS)
    !         *NGST* - IF = 1 THEN NO GUSTINESS PARAMETERISATION
    !                - IF = 2 THEN GUSTINESS PARAMETERISATION
    !         *LLSNEG- IF TRUE THEN THE NEGATIVE SINPUT (SWELL DAMPING) WILL BE COMPUTED
    !         *KIJS* - INDEX OF FIRST GRIDPOINT.
    !         *KIJL* - INDEX OF LAST GRIDPOINT.
    !          *FL1* - SPECTRUM.
    !       *WAVNUM* - WAVE NUMBER.
    !         *CINV* - INVERSE PHASE VELOCITY.
    !       *XK2CG*  - (WAVNUM)**2 * GROUP SPPED.
    !       *WDWAVE* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
    !                  NOTATION (POINTING ANGLE OF WIND VECTOR,
    !                  CLOCKWISE FROM NORTH).
    !        *UFRIC* - NEW FRICTION VELOCITY IN M/S.
    !        *Z0M* - ROUGHNESS LENGTH IN M.
    !      *COSWDIF* - COS(TH(K)-WDWAVE(IJ))
    !     *SINWDIF2* - SIN(TH(K)-WDWAVE(IJ))**2
    !        *RAORW* - RATIO AIR DENSITY TO WATER DENSITY.
    !        *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
    !        *RNFAC* - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
    !          *FLD* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
    !           *SL* - TOTAL SOURCE FUNCTION ARRAY.
    !         *SPOS* - POSITIVE SOURCE FUNCTION ARRAY.
    !       *XLLWS*  - = 1 WHERE SINPUT IS POSITIVE
    
    !     METHOD.
    !     -------
    
    !       SEE REFERENCE.
    
    
    ! ----------------------------------------------------------------------
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWCOUP, ONLY: LLCAPCHNK, LLNORMAGAM
    USE YOWFRED, ONLY: FR, TH, DFIM, COSTH, SINTH, ZPIFR, DELTH
    USE YOWPARAM, ONLY: NANG, NFRE, NANG_PARAM
    USE YOWPCONS, ONLY: G, GM1, EPSMIN, EPSUS, ZPI
    USE YOWPHYS, ONLY: ZALP, TAUWSHELTER, XKAPPA, BETAMAXOXKAPPA2, RN1_RN, RNU, RNUM, SWELLF, SWELLF2, SWELLF3, SWELLF4,  &
    & SWELLF5, SWELLF6, SWELLF7, SWELLF7M1, Z0RAT, Z0TUBMAX, ABMIN, ABMAX
    USE YOWTEST, ONLY: IU06
    USE YOWTABL, ONLY: IAB, SWELLFT
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: NGST
    LOGICAL, INTENT(IN) :: LLSNEG
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM, CINV, XK2CG
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: WDWAVE, WSWAVE, UFRIC, Z0M
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: RAORW, WSTAR, RNFAC
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: COSWDIF, SINWDIF2
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: FLD, SL, SPOS
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: XLLWS
    
    
    INTEGER(KIND=JWIM) :: IJ, K, M, IND, IGST
    
    REAL(KIND=JWRB) :: CONSTN
    REAL(KIND=JWRB) :: AVG_GST, ABS_TAUWSHELTER
    REAL(KIND=JWRB) :: CONST1
    REAL(KIND=JWRB) :: ZNZ
    REAL(KIND=JWRB) :: X1, X2, ZLOG, ZLOG1, ZLOG2, ZLOG2X, XV1, XV2, ZBETA1, ZBETA2
    REAL(KIND=JWRB) :: XI, X, DELI1, DELI2
    REAL(KIND=JWRB) :: FU, FUD, NU_AIR, SMOOTH, HFTSWELLF6, Z0TUB
    REAL(KIND=JWRB) :: FAC_NU_AIR, FACM1_NU_AIR
    REAL(KIND=JWRB) :: ARG, DELABM1
    REAL(KIND=JWRB) :: TAUPX, TAUPY
    REAL(KIND=JWRB) :: DSTAB2
    
    REAL(KIND=JWRB) :: SIG2, COEF, COEF5, DFIM_SIG2, COSLP
    
    REAL(KIND=JWRB) :: XNGAMCONST
    REAL(KIND=JWRB) :: CONSTF, CONST11, CONST22
    REAL(KIND=JWRB) :: Z0VIS, Z0NOZ, FWW
    REAL(KIND=JWRB) :: PVISC, PTURB
    REAL(KIND=JWRB) :: ZCN
    REAL(KIND=JWRB) :: SIG_N, UORBT, AORB, TEMP, RE, RE_C, ZORB
    REAL(KIND=JWRB) :: CNSN, SUMF, SUMFSIN2
    REAL(KIND=JWRB) :: CSTRNFAC
    REAL(KIND=JWRB) :: FLP_AVG, SLP_AVG
    REAL(KIND=JWRB) :: ROGOROAIR, AIRD_PVISC
    REAL(KIND=JWRB) :: DSTAB1, TEMP1, TEMP2
    
    REAL(KIND=JWRB), DIMENSION(2) :: XSTRESS, YSTRESS, FLP, SLP
    REAL(KIND=JWRB), DIMENSION(2) :: USG2, TAUX, TAUY, USTP, USTPM1, USDIRP, UCN
    REAL(KIND=JWRB), DIMENSION(2) :: UCNZALPD
    REAL(KIND=JWRB), DIMENSION(2) :: GAMNORMA    ! ! RENORMALISATION FACTOR OF THE GROWTH RATE
    REAL(KIND=JWRB), DIMENSION(2, NANG_PARAM) :: GAM0, DSTAB
    
    LOGICAL :: LTAUWSHELTER
!$acc routine vector
!$acc data present( FL1, WAVNUM, CINV, XK2CG, WDWAVE, WSWAVE, UFRIC, Z0M, COSWDIF, SINWDIF2, RAORW, WSTAR, RNFAC, FLD, SL, SPOS,  &
!$acc & XLLWS )
    ! ----------------------------------------------------------------------
    
    
    AVG_GST = 1.0_JWRB / NGST
    CONST1 = BETAMAXOXKAPPA2
    CONSTN = DELTH / (XKAPPA*ZPI)
    
    ABS_TAUWSHELTER = ABS(TAUWSHELTER)
    IF (ABS_TAUWSHELTER == 0.0_JWRB) THEN
      LTAUWSHELTER = .false.
    ELSE
      LTAUWSHELTER = .true.
    END IF
!$loki separator
    
!$acc loop vector private( XSTRESS, YSTRESS, FLP, SLP, USG2, TAUX, TAUY, USTP, USTPM1, USDIRP, UCN, UCNZALPD, GAMNORMA, GAM0,  &
!$acc & DSTAB )
    DO IJ=KIJS,KIJL
      
      IF (NGST > 1) THEN
!$loki inline
        CALL WSIGSTAR_SCC(WSWAVE(IJ), UFRIC(IJ), Z0M(IJ), WSTAR(IJ), SIG_N)
      END IF
      
      
      IF (LLNORMAGAM) THEN
        CSTRNFAC = CONSTN*RNFAC(IJ) / RAORW(IJ)
      END IF
      
      
      !     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
      
      ! ----------------------------------------------------------------------
      IF (LLSNEG) THEN
        !!!!  only for the negative sinput
        NU_AIR = RNU
        FACM1_NU_AIR = 4.0_JWRB / NU_AIR
        
        FAC_NU_AIR = RNUM
        
        FU = ABS(SWELLF3)
        FUD = SWELLF2
        DELABM1 = REAL(IAB) / (ABMAX - ABMIN)
        
        
        !       computation of Uorb and Aorb
        UORBT = EPSMIN
        AORB = EPSMIN
        
!$acc loop seq
        DO M=1,NFRE
          SIG2 = ZPIFR(M)**2
          DFIM_SIG2 = DFIM(M)*SIG2
          
          K = 1
          TEMP = FL1(IJ, K, M)
!$acc loop seq
          DO K=2,NANG
            TEMP = TEMP + FL1(IJ, K, M)
          END DO
          
          UORBT = UORBT + DFIM_SIG2*TEMP
          AORB = AORB + DFIM(M)*TEMP
        END DO
        
        UORBT = 2.0_JWRB*SQRT(UORBT)          ! this is the significant orbital amplitude
        AORB = 2.0_JWRB*SQRT(AORB)          ! this 1/2 Hs
        RE = FACM1_NU_AIR*UORBT*AORB          ! this is the Reynolds number
        Z0VIS = FAC_NU_AIR / MAX(UFRIC(IJ), 0.0001_JWRB)
        Z0TUB = Z0RAT*MIN(Z0TUBMAX, Z0M(IJ))
        Z0NOZ = MAX(Z0VIS, Z0TUB)
        ZORB = AORB / Z0NOZ
        
        !         compute fww
        XI = (LOG10(MAX(ZORB, 3.0_JWRB)) - ABMIN)*DELABM1
        IND = MIN(IAB - 1, INT(XI))
        DELI1 = MIN(1.0_JWRB, XI - REAL(IND, kind=JWRB))
        DELI2 = 1.0_JWRB - DELI1
        FWW = SWELLFT(IND)*DELI2 + SWELLFT(IND + 1)*DELI1
        TEMP2 = FWW*UORBT
        
        !       Define the critical Reynolds number
        IF (SWELLF6 == 1.0_JWRB) THEN
          RE_C = SWELLF4
        ELSE
          HFTSWELLF6 = 1.0_JWRB - SWELLF6
          RE_C = SWELLF4*(2.0_JWRB / AORB)**HFTSWELLF6
        END IF
        
        !       Swell damping weight between viscous and turbulent boundary layer
        IF (SWELLF7 > 0.0_JWRB) THEN
          SMOOTH = 0.5_JWRB*TANH((RE - RE_C)*SWELLF7M1)
          PTURB = 0.5_JWRB + SMOOTH
          PVISC = 0.5_JWRB - SMOOTH
        ELSE
          IF (RE <= RE_C) THEN
            PTURB = 0.0_JWRB
            PVISC = 0.5_JWRB
          ELSE
            PTURB = 0.5_JWRB
            PVISC = 0.0_JWRB
          END IF
        END IF
        
        AIRD_PVISC = PVISC*RAORW(IJ)
        
      END IF
      
      
      
      ! Initialisation
      
      IF (NGST == 1) THEN
        USTP(1) = UFRIC(IJ)
      ELSE
        USTP(1) = UFRIC(IJ)*(1.0_JWRB + SIG_N)
        USTP(2) = UFRIC(IJ)*(1.0_JWRB - SIG_N)
      END IF
      
!$acc loop seq
      DO IGST=1,NGST
        USTPM1(IGST) = 1.0_JWRB / MAX(USTP(IGST), EPSUS)
      END DO
      
      IF (LTAUWSHELTER) THEN
!$acc loop seq
        DO IGST=1,NGST
          XSTRESS(IGST) = 0.0_JWRB
          YSTRESS(IGST) = 0.0_JWRB
          USG2(IGST) = USTP(IGST)**2
          TAUX(IGST) = USG2(IGST)*SIN(WDWAVE(IJ))
          TAUY(IGST) = USG2(IGST)*COS(WDWAVE(IJ))
        END DO
        
        ROGOROAIR = G / RAORW(IJ)
      END IF
      
      
      !*    2. MAIN LOOP OVER FREQUENCIES.
      !        ---------------------------
      
      IF (.not.LLNORMAGAM) THEN
!$acc loop seq
        DO IGST=1,NGST
          GAMNORMA(IGST) = 1.0_JWRB
        END DO
      END IF
      
      IF (.not.LLSNEG) THEN
!$acc loop seq
        DO K=1,NANG
!$acc loop seq
          DO IGST=1,NGST
            DSTAB(IGST, K) = 0.0_JWRB
          END DO
        END DO
      END IF
      
!$acc loop seq
      DO M=1,NFRE
        
        IF (LTAUWSHELTER) THEN
!$acc loop seq
          DO IGST=1,NGST
            TAUPX = TAUX(IGST) - ABS_TAUWSHELTER*XSTRESS(IGST)
            TAUPY = TAUY(IGST) - ABS_TAUWSHELTER*YSTRESS(IGST)
            USDIRP(IGST) = ATAN2(TAUPX, TAUPY)
            USTP(IGST) = (TAUPX**2 + TAUPY**2)**0.25_JWRB
            USTPM1(IGST) = 1.0_JWRB / MAX(USTP(IGST), EPSUS)
          END DO
          
          CONSTF = ROGOROAIR*CINV(IJ, M)*DFIM(M)
        END IF
        
        
        !*      PRECALCULATE FREQUENCY DEPENDENCE.
        !       ----------------------------------
        
!$acc loop seq
        DO IGST=1,NGST
          UCN(IGST) = USTP(IGST)*CINV(IJ, M)
          UCNZALPD(IGST) = XKAPPA / (UCN(IGST) + ZALP)
        END DO
        ZCN = LOG(WAVNUM(IJ, M)*Z0M(IJ))
        CNSN = ZPIFR(M)*CONST1*RAORW(IJ)
        
        !*    2.1 LOOP OVER DIRECTIONS.
        !         ---------------------
        
!$acc loop seq
        DO K=1,NANG
          XLLWS(IJ, K, M) = 0.0_JWRB
        END DO
        
        IF (LLSNEG) THEN
          !       SWELL DAMPING:
          
          SIG2 = ZPIFR(M)**2
          DFIM_SIG2 = DFIM(M)*SIG2
          
          COEF = -SWELLF*16._JWRB*SIG2 / G
          COEF5 = -SWELLF5*2._JWRB*SQRT(2._JWRB*NU_AIR*ZPIFR(M))
          
          DSTAB1 = COEF5*AIRD_PVISC*WAVNUM(IJ, M)
          TEMP1 = COEF*RAORW(IJ)
        END IF
        
!$acc loop seq
        DO K=1,NANG
!$acc loop seq
          DO IGST=1,NGST
            
            SUMF = 0.0_JWRB
            SUMFSIN2 = 0.0_JWRB
            
            IF (LTAUWSHELTER) THEN
              COSLP = COS(TH(K) - USDIRP(IGST))
            ELSE
              COSLP = COSWDIF(IJ, K)
            END IF
            
            GAM0(IGST, K) = 0._JWRB
            IF (COSLP > 0.01_JWRB) THEN
              X = COSLP*UCN(IGST)
              ZLOG = ZCN + UCNZALPD(IGST) / COSLP
              IF (ZLOG < 0.0_JWRB) THEN
                ZLOG2X = ZLOG*ZLOG*X
                GAM0(IGST, K) = EXP(ZLOG)*ZLOG2X*ZLOG2X*CNSN
                XLLWS(IJ, K, M) = 1.0_JWRB
              END IF
            END IF
            
            IF (LLSNEG) THEN
              DSTAB2 = TEMP1*(TEMP2 + (FU + FUD*COSLP)*USTP(IGST))
              DSTAB(IGST, K) = DSTAB1 + PTURB*DSTAB2
            END IF
            
            SUMF = SUMF + GAM0(IGST, K)*FL1(IJ, K, M)
            SUMFSIN2 = SUMFSIN2 + GAM0(IGST, K)*FL1(IJ, K, M)*SINWDIF2(IJ, K)
          END DO
        END DO
        
        IF (LLNORMAGAM) THEN
          
          XNGAMCONST = CSTRNFAC*XK2CG(IJ, M)
!$acc loop seq
          DO IGST=1,NGST
            ZNZ = XNGAMCONST*USTPM1(IGST)
            GAMNORMA(IGST) = (1.0_JWRB + ZNZ*SUMFSIN2) / (1.0_JWRB + ZNZ*SUMF)
          END DO
          
        END IF
        
        
        
        !*    2.2 UPDATE THE SHELTERING STRESS (in any),
        !         AND THEN ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
        !         ---------------------------------------------------------
        
!$acc loop seq
        DO K=1,NANG
          
!$acc loop seq
          DO IGST=1,NGST
            ! SLP: only the positive contributions
            SLP(IGST) = GAM0(IGST, K)*GAMNORMA(IGST)
            FLP(IGST) = SLP(IGST) + DSTAB(IGST, K)
          END DO
          
!$acc loop seq
          DO IGST=1,NGST
            SLP(IGST) = SLP(IGST)*FL1(IJ, K, M)
          END DO
          
          IF (LTAUWSHELTER) THEN
            CONST11 = CONSTF*SINTH(K)
            CONST22 = CONSTF*COSTH(K)
!$acc loop seq
            DO IGST=1,NGST
              XSTRESS(IGST) = XSTRESS(IGST) + SLP(IGST)*CONST11
              YSTRESS(IGST) = YSTRESS(IGST) + SLP(IGST)*CONST22
            END DO
          END IF
          
          IGST = 1
          SLP_AVG = SLP(IGST)
          FLP_AVG = FLP(IGST)
!$acc loop seq
          DO IGST=2,NGST
            SLP_AVG = SLP_AVG + SLP(IGST)
            FLP_AVG = FLP_AVG + FLP(IGST)
          END DO
          
          SPOS(IJ, K, M) = AVG_GST*SLP_AVG
          FLD(IJ, K, M) = AVG_GST*FLP_AVG
          SL(IJ, K, M) = FLD(IJ, K, M)*FL1(IJ, K, M)
          
        END DO
        
      END DO
      ! END LOOP OVER FREQUENCIES
      
      
    END DO
!$acc end data
  END SUBROUTINE SINPUT_ARD_SCC
  SUBROUTINE SINPUT_JAN_SCC (NGST, LLSNEG, KIJS, KIJL, FL1, WAVNUM, CINV, XK2CG, WSWAVE, UFRIC, Z0M, COSWDIF, SINWDIF2, RAORW,  &
  & WSTAR, RNFAC, FLD, SL, SPOS, XLLWS)
    ! ----------------------------------------------------------------------
    
    !**** *SINPUT_JAN* - COMPUTATION OF INPUT SOURCE FUNCTION.
    
    !     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
    
    !     OPTIMIZED BY : H. GUENTHER
    
    !     MODIFIED BY :
    !       J-R BIDLOT NOVEMBER 1995
    !       J-R BIDLOT FEBRUARY 1996-97
    !       J-R BIDLOT FEBRUARY 1999 : INTRODUCE ICALL AND NCALL
    !       P.A.E.M. JANSSEN MAY 2000 : INTRODUCE GUSTINESS
    !       J-R BIDLOT FEBRUARY 2001 : MAKE IT FULLY IMPLICIT BY ONLY
    !                                  USING NEW STRESS AND ROUGHNESS.
    !       S. ABDALLA OCTOBER 2001:  INTRODUCTION OF VARIABLE AIR
    !                                 DENSITY AND STABILITY-DEPENDENT
    !                                 WIND GUSTINESS
    !       P.A.E.M. JANSSEN OCTOBER 2008: INTRODUCE DAMPING WHEN WAVES ARE
    !                                      RUNNING FASTER THAN THE WIND.
    !       J-R BIDLOT JANUARY 2013: SHALLOW WATER FORMULATION.
    
    !*    PURPOSE.
    !     ---------
    
    !       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
    !       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
    !       INPUT SOURCE FUNCTION.
    !
    !       GUSTINESS IS INTRODUCED FOLL0WING THE APPROACH OF JANSSEN(1986),
    !       USING A GAUSS-HERMITE APPROXIMATION SUGGESTED BY MILES(1997).
    !       IN THE PRESENT VERSION ONLY TWO HERMITE POLYNOMIALS ARE UTILISED
    !       IN THE EVALUATION OF THE PROBABILITY INTEGRAL. EXPLICITELY ONE THEN
    !       FINDS:
    !
    !             <GAMMA(X)> = 0.5*( GAMMA(X(1+SIG)) + GAMMA(X(1-SIG)) )
    !
    !       WHERE X IS THE FRICTION VELOCITY AND SIG IS THE RELATIVE GUSTINESS
    !       LEVEL.
    
    !**   INTERFACE.
    !     ----------
    
    !     *CALL* *SINPUT_JAN (NGST, LLSNEG, KIJS, KIJL, FL1,
    !    &                    WAVNUM, CINV, XK2CG,
    !    &                    WDWAVE, WSWAVE, UFRIC, Z0M,
    !    &                    COSWDIF, SINWDIF2,
    !    &                    RAORW, WSTAR, RNFAC,
    !    &                    FLD, SL, SPOS, XLLWS)
    !         *NGST* - IF = 1 THEN NO GUSTINESS PARAMETERISATION
    !                - IF = 2 THEN GUSTINESS PARAMETERISATION
    !         *LLSNEG- IF TRUE THEN THE NEGATIVE SINPUT (SWELL DAMPING) WILL BE COMPUTED
    !         *KIJS* - INDEX OF FIRST GRIDPOINT.
    !         *KIJL* - INDEX OF LAST GRIDPOINT.
    !          *FL1* - SPECTRUM.
    !       *WAVNUM* - WAVE NUMBER.
    !         *CINV* - INVERSE PHASE VELOCITY.
    !       *XK2CG*  - (WAVNUM)**2 * GROUP SPPED.
    !       *WDWAVE* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
    !                  NOTATION (POINTING ANGLE OF WIND VECTOR,
    !                  CLOCKWISE FROM NORTH).
    !        *UFRIC* - FRICTION VELOCITY IN M/S.
    !        *Z0M*   - ROUGHNESS LENGTH IN M.
    !      *COSWDIF* - COS(TH(K)-WDWAVE(IJ))
    !     *SINWDIF2* - SIN(TH(K)-WDWAVE(IJ))**2
    !        *RAORW* - RATIO AIR DENSITY TO WATER DENSITY
    !        *RNFAC* - WIND DEPENDENT FACTOR USED IN THE GROWTH RENORMALISATION.
    !        *WSTAR* - FREE CONVECTION VELOCITY SCALE (M/S).
    !          *FLD* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
    !           *SL* - TOTAL SOURCE FUNCTION ARRAY.
    !         *SPOS* - ONLY POSITIVE PART OF INPUT SOURCE FUNCTION ARRAY.
    !        *XLLWS* - 1 WHERE SINPUT IS POSITIVE.
    
    
    !     METHOD.
    !     -------
    
    !       SEE REFERENCE.
    
    !     EXTERNALS.
    !     ----------
    
    !       WSIGSTAR.
    
    !     MODIFICATIONS
    !     -------------
    
    !     - REMOVAL OF CALL TO CRAY SPECIFIC FUNCTIONS EXPHF AND ALOGHF
    !       BY THEIR STANDARD FORTRAN EQUIVALENT EXP and ALOGHF
    !     - MODIFIED TO MAKE INTEGRATION SCHEME FULLY IMPLICIT
    !     - INTRODUCTION OF VARIABLE AIR DENSITY
    !     - INTRODUCTION OF WIND GUSTINESS
    
    !     REFERENCE.
    !     ----------
    
    !       P. JANSSEN, J.P.O., 1989.
    !       P. JANSSEN, J.P.O., 1991
    
    ! ----------------------------------------------------------------------
    
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWCOUP, ONLY: LLNORMAGAM
    USE YOWFRED, ONLY: ZPIFR, DELTH, TH
    USE YOWFRED, ONLY: FR, TH, ZPIFR
    USE YOWPARAM, ONLY: NANG, NFRE, NANG_PARAM
    USE YOWPCONS, ONLY: G, GM1, ZPI, EPSUS
    USE YOWPHYS, ONLY: ZALP, XKAPPA, BETAMAXOXKAPPA2
    USE YOWSTAT, ONLY: IDAMPING
    USE YOWTEST, ONLY: IU06
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: NGST
    LOGICAL, INTENT(IN) :: LLSNEG
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM, CINV, XK2CG
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: WSWAVE, UFRIC, Z0M
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG) :: COSWDIF, SINWDIF2
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: RAORW, WSTAR, RNFAC
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: FLD, SL, SPOS
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL, NANG, NFRE) :: XLLWS
    
    
    INTEGER(KIND=JWIM) :: IJ, IG, K, M
    INTEGER(KIND=JWIM) :: IGST
    
    REAL(KIND=JWRB) :: CONST1, CONST3, XKAPPAD
    REAL(KIND=JWRB) :: CONSTN
    REAL(KIND=JWRB) :: ZNZ
    REAL(KIND=JWRB) :: X, ZLOG, ZLOG2X, ZBETA
    REAL(KIND=JWRB) :: TEMPD
    
    REAL(KIND=JWRB), DIMENSION(2) :: WSIN
    REAL(KIND=JWRB) :: ZTANHKD
    REAL(KIND=JWRB) :: SIG_N
    REAL(KIND=JWRB) :: CNSN
    REAL(KIND=JWRB) :: SUMF, SUMFSIN2
    REAL(KIND=JWRB) :: CSTRNFAC
    REAL(KIND=JWRB) :: UFAC1, UFAC2
    REAL(KIND=JWRB), DIMENSION(2) :: GAMNORMA    ! ! RENORMALISATION FACTOR OF THE GROWTH RATE
    REAL(KIND=JWRB), DIMENSION(2) :: SIGDEV, US, Z0, UCN, ZCN
    REAL(KIND=JWRB), DIMENSION(2) :: USTPM1
    REAL(KIND=JWRB), DIMENSION(2) :: XVD, UCND, CONST3_UCN2
    REAL(KIND=JWRB), DIMENSION(2, NANG_PARAM) :: GAM0
    
    LOGICAL :: LZ
!$acc routine vector
!$acc data present( FL1, WAVNUM, CINV, XK2CG, WSWAVE, UFRIC, Z0M, COSWDIF, SINWDIF2, RAORW, WSTAR, RNFAC, FLD, SL, SPOS, XLLWS )
    
    ! ----------------------------------------------------------------------
    
    
    CONST1 = BETAMAXOXKAPPA2
    CONST3 = 2.0_JWRB*XKAPPA / CONST1      ! SEE IDAMPING
    XKAPPAD = 1.E0_JWRB / XKAPPA
    
    CONST3 = IDAMPING*CONST3
    
    CONSTN = DELTH / (XKAPPA*ZPI)
!$loki separator
    
!$acc loop vector private( WSIN, GAMNORMA, SIGDEV, US, Z0, UCN, ZCN, USTPM1, XVD, UCND, CONST3_UCN2, GAM0 )
    DO IJ=KIJS,KIJL
      
      !     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
      IF (NGST > 1) THEN
!$loki inline
        CALL WSIGSTAR_SCC(WSWAVE(IJ), UFRIC(IJ), Z0M(IJ), WSTAR(IJ), SIG_N)
      END IF
      
      !     DEFINE WHERE SINPUT WILL BE EVALUATED IN RELATIVE TERM WRT USTAR
      !     DEFINE ALSO THE RELATIVE WEIGHT OF EACH.
      
      IF (NGST == 1) THEN
        WSIN(1) = 1.0_JWRB
        SIGDEV(1) = 1.0_JWRB
      ELSE
        WSIN(1) = 0.5_JWRB
        WSIN(2) = 0.5_JWRB
        SIGDEV(1) = 1.0_JWRB - SIG_N
        SIGDEV(2) = 1.0_JWRB + SIG_N
      END IF
      
      
      IF (NGST == 1) THEN
        US(1) = UFRIC(IJ)
        Z0(1) = Z0M(IJ)
      ELSE
!$acc loop seq
        DO IGST=1,NGST
          US(IGST) = UFRIC(IJ)*SIGDEV(IGST)
          Z0(IGST) = Z0M(IJ)
        END DO
      END IF
      
!$acc loop seq
      DO IGST=1,NGST
        USTPM1(IGST) = 1.0_JWRB / MAX(US(IGST), EPSUS)
      END DO
      
      ! ----------------------------------------------------------------------
      
      !*    2. LOOP OVER FREQUENCIES.
      !        ----------------------
      
!$acc loop seq
      DO M=1,NFRE
        
        !*      PRECALCULATE FREQUENCY DEPENDENCE.
        !       ----------------------------------
        
        ZTANHKD = ZPIFR(M)**2 / (G*WAVNUM(IJ, M))
        CNSN = CONST1*ZPIFR(M)*ZTANHKD*RAORW(IJ)
        
!$acc loop seq
        DO IGST=1,NGST
          UCN(IGST) = US(IGST)*CINV(IJ, M) + ZALP
          CONST3_UCN2(IGST) = CONST3*UCN(IGST)**2
          UCND(IGST) = 1.0_JWRB / UCN(IGST)
          ZCN(IGST) = LOG(WAVNUM(IJ, M)*Z0(IGST))
          XVD(IGST) = 1.0_JWRB / (-US(IGST)*XKAPPAD*ZCN(IGST)*CINV(IJ, M))
        END DO
        
        !*    2.1 LOOP OVER DIRECTIONS.
        !         ---------------------
        
        !       WIND INPUT:
!$acc loop seq
        DO K=1,NANG
          XLLWS(IJ, K, M) = 0.0_JWRB
          
!$acc loop seq
          DO IGST=1,NGST
            
            IF (COSWDIF(IJ, K) > 0.01_JWRB) THEN
              LZ = .true.
              TEMPD = XKAPPA / COSWDIF(IJ, K)
            ELSE
              LZ = .false.
              TEMPD = XKAPPA
            END IF
            
            GAM0(IGST, K) = 0.0_JWRB
            IF (LZ) THEN
              ZLOG = ZCN(IGST) + TEMPD*UCND(IGST)
              IF (ZLOG < 0.0_JWRB) THEN
                X = COSWDIF(IJ, K)*UCN(IGST)
                ZLOG2X = ZLOG*ZLOG*X
                GAM0(IGST, K) = ZLOG2X*ZLOG2X*EXP(ZLOG)*CNSN
                XLLWS(IJ, K, M) = 1.0_JWRB
              END IF
            END IF
          END DO
          
        END DO
        
        
        IF (LLNORMAGAM) THEN
          
          SUMF = 0.0_JWRB
          SUMFSIN2 = 0.0_JWRB
!$acc loop seq
          DO K=1,NANG
!$acc loop seq
            DO IGST=1,NGST
              SUMF = SUMF + GAM0(IGST, K)*FL1(IJ, K, M)
              SUMFSIN2 = SUMFSIN2 + GAM0(IGST, K)*FL1(IJ, K, M)*SINWDIF2(IJ, K)
            END DO
            
            CSTRNFAC = CONSTN*RNFAC(IJ) / RAORW(IJ)
            ZNZ = CSTRNFAC*XK2CG(IJ, M)*USTPM1(IGST)
            GAMNORMA(IGST) = (1.0_JWRB + ZNZ*SUMFSIN2) / (1.0_JWRB + ZNZ*SUMF)
            
          END DO
        ELSE
!$acc loop seq
          DO IGST=1,NGST
            GAMNORMA(IGST) = 1.0_JWRB
          END DO
        END IF
        
!$acc loop seq
        DO K=1,NANG
          UFAC1 = WSIN(1)*GAM0(1, K)*GAMNORMA(1)
!$acc loop seq
          DO IGST=2,NGST
            UFAC1 = UFAC1 + WSIN(IGST)*GAM0(IGST, K)*GAMNORMA(IGST)
          END DO
          
          UFAC2 = 0.0_JWRB
          IF (LLSNEG) THEN
            !         SWELL DAMPING:
            ZBETA = CONST3_UCN2(1)*(COSWDIF(IJ, K) - XVD(1))
            UFAC2 = WSIN(1)*ZBETA
!$acc loop seq
            DO IGST=2,NGST
              ZBETA = CONST3_UCN2(IGST)*(COSWDIF(IJ, K) - XVD(IGST))
              UFAC2 = UFAC2 + WSIN(IGST)*ZBETA
            END DO
          END IF
          
          FLD(IJ, K, M) = UFAC1 + UFAC2*CNSN
          SPOS(IJ, K, M) = UFAC1*FL1(IJ, K, M)
          SL(IJ, K, M) = FLD(IJ, K, M)*FL1(IJ, K, M)
        END DO
        
        !*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
        !         ------------------------------------------------
        
      END DO
      
      
    END DO
!$acc end data
  END SUBROUTINE SINPUT_JAN_SCC
END MODULE SINPUT_ARD_SCC_MOD
