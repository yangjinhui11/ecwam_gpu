! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE STOKESDRIFT_SCC_MOD
  CONTAINS
  SUBROUTINE STOKESDRIFT_SCC (KIJS, KIJL, FL1, STOKFAC, WSWAVE, WDWAVE, CICOVER, USTOKES, VSTOKES)
    
    !
    !***  *STOKESDRIFT*   DETERMINES THE STOKES DRIFT
    !
    !     PETER JANSSEN MARCH 2009
    !
    !     PURPOSE.
    !     --------
    !
    !              DETERMINATION OF STOKES DRIFT VECTOR
    !
    !     INTERFACE.
    !     ----------
    !              *CALL*  *STOKESDRIFT(KIJS, KIJL, FL1, STOKFAC, WSWAVE,WDWAVE,CICOVER,USTOKES,VSTOKES)*
    !
    !                       INPUT:
    !                            *KIJS*   - FIRST GRIDPOINT
    !                            *KIJL*   - LAST GRIDPOINT
    !                            *FL1*    - 2-D SPECTRUM
    !                            *STOKFAC*- FACTOR TO COMPUTE THE STOKES DRIFT
    !                            Auxilliary fields to specify Stokes when model sea ice cover the blocking threshold
    !                            as 0.016*WSWAVE, aligned in the wind direction
    !                            *WSWAVE* - WIND SPEED IN M/S.
    !                            *WDWAVE* - WIND DIRECTION IN RADIANS.
    !                            *CICOVER*- SEA ICE COVER.
    !
    !                       OUTPUT:
    !                            *USTOKES*   - U-COMPONENT STOKES DRIFT
    !                            *VSTOKES*   - V-COMPONENT STOKES DRIFT
    !
    !     METHOD.
    !     -------
    !              DETERMINE U- AND V-COMPONENT OF STOKES DRIFT FOLLOWING
    !              K.E. KENYON, J.G.R., 74, 6991-6994
    !
    !     EXTERNALS.
    !     ----------
    !              NONE
    !
    !
    !-----------------------------------------------------------------------
    
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWPCONS, ONLY: G, ZPI
    USE YOWFRED, ONLY: FR, DFIM, DELTH, TH, DFIM_SIM, FRATIO, COSTH, SINTH
    USE YOWICE, ONLY: LICERUN, LWAMRSETCI, CITHRSH
    USE YOWPARAM, ONLY: NANG, NFRE, NFRE_ODD
    
    
    ! ----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: STOKFAC
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: WSWAVE, WDWAVE, CICOVER
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: USTOKES, VSTOKES
    
    
    INTEGER(KIND=JWIM) :: IJ, M, K
    
    REAL(KIND=JWRB), PARAMETER :: STMAX = 1.5_JWRB    ! maximum magnitude (this is for safety when coupled)
    REAL(KIND=JWRB) :: CONST, FAC, FAC1, FAC2, FAC3
    REAL(KIND=JWRB) :: STFAC
!$acc routine vector
!$acc data present( FL1, STOKFAC, WSWAVE, WDWAVE, CICOVER, USTOKES, VSTOKES )
    
!$acc loop vector
    DO IJ=KIJS,KIJL
      
      ! ----------------------------------------------------------------------
      
      
      
      !***  1. DETERMINE STOKE DRIFT VECTOR.
      !     --------------------------------
      
      CONST = 2.0_JWRB*DELTH*ZPI**3 / G*FR(NFRE_ODD)**4
      
      !***  1.1 PERFORM INTEGRATION.
      !     ------------------------
      
      USTOKES(IJ) = 0.0_JWRB
      VSTOKES(IJ) = 0.0_JWRB
      
!$acc loop seq
      DO M=1,NFRE_ODD
        STFAC = STOKFAC(IJ, M)*DFIM_SIM(M)
!$acc loop seq
        DO K=1,NANG
          FAC3 = STFAC*FL1(IJ, K, M)
          USTOKES(IJ) = USTOKES(IJ) + FAC3*SINTH(K)
          VSTOKES(IJ) = VSTOKES(IJ) + FAC3*COSTH(K)
        END DO
      END DO
      
      !***  1.2 ADD CONTRIBUTION OF UNRESOLVED WAVES.
      !     -----------------------------------------
      
!$acc loop seq
      DO K=1,NANG
        FAC1 = CONST*SINTH(K)
        FAC2 = CONST*COSTH(K)
        USTOKES(IJ) = USTOKES(IJ) + FAC1*FL1(IJ, K, NFRE_ODD)
        VSTOKES(IJ) = VSTOKES(IJ) + FAC2*FL1(IJ, K, NFRE_ODD)
      END DO
      
      
      !***  1.3 Sea Ice exception
      !     ---------------------
      IF (LICERUN .and. LWAMRSETCI) THEN
        IF (CICOVER(IJ) > CITHRSH) THEN
          USTOKES(IJ) = 0.016_JWRB*WSWAVE(IJ)*SIN(WDWAVE(IJ))*(1.0_JWRB - CICOVER(IJ))
          VSTOKES(IJ) = 0.016_JWRB*WSWAVE(IJ)*COS(WDWAVE(IJ))*(1.0_JWRB - CICOVER(IJ))
        END IF
      END IF
      
      !***  1.4 Protection
      !     --------------
      
      USTOKES(IJ) = MIN(MAX(USTOKES(IJ), -STMAX), STMAX)
      VSTOKES(IJ) = MIN(MAX(VSTOKES(IJ), -STMAX), STMAX)
      
      
    END DO
!$acc end data
  END SUBROUTINE STOKESDRIFT_SCC
END MODULE STOKESDRIFT_SCC_MOD
