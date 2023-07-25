! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE CIMSSTRN_SCC_MOD
  CONTAINS
  SUBROUTINE CIMSSTRN_SCC (KIJS, KIJL, FL1, WAVNUM, DEPTH, CITHICK, STRN)
    
    ! ----------------------------------------------------------------------
    
    !**** *CIMSSTRN* - COMPUTATION OF THE MEAN SQUARE WAVE STRAIN IN SEA ICE.
    
    !     J. BIDLOT  ECMWF  JANUARY 2013.
    
    !*    PURPOSE.
    !     --------
    
    !       COMPUTES MEAN SQUARE WAVE STRAIN AT EACH GRID POINT.
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *CIMSSTRN (KIJS, KIJL, FL1, WAVNUM, DEPTH, CITHICK, STRN)*
    !              *KIJS*    - INDEX OF FIRST GRIDPOINT
    !              *KIJL*    - INDEX OF LAST GRIDPOINT
    !              *FL1*     - SPECTRUM.
    !              *WAVNUM*  - OPEN WATER WAVE NUMBER
    !              *DEPTH*   - WATER DEPTH
    !              *CITHICK* - SEA ICE THICKNESS
    !              *STRN*    - MEAN SQUARE WAVE STRAIN IN ICE (OUTPUT).
    
    !     METHOD.
    !     -------
    
    !      !!! IT ASSUMES SO DEFAULT SETTING FOR THE MECHANICAL PROPERTIES OF
    !          THE SEA ICE (SEE AKI_ICE) !!!!!!!
    
    !       NONE.
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    !       NONE.
    
    ! ----------------------------------------------------------------------
    
    USE AKI_ICE_SCC_MOD, ONLY: AKI_ICE_SCC
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWICE, ONLY: FLMIN
    USE YOWFRED, ONLY: FR, DFIM, DELTH
    USE YOWPARAM, ONLY: NANG, NFRE
    USE YOWPCONS, ONLY: G, ZPI, ROWATER
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: DEPTH
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: CITHICK
    REAL(KIND=JWRB), INTENT(OUT), DIMENSION(KIJL) :: STRN
    
    
    INTEGER(KIND=JWIM) :: IJ, M, K
    REAL(KIND=JWRB) :: F1LIM
    REAL(KIND=JWRB) :: XKI, E, SUME
!$acc routine vector
!$acc data present( FL1, WAVNUM, DEPTH, CITHICK, STRN )
    
!$acc loop vector
    DO IJ=KIJS,KIJL
      
      ! ----------------------------------------------------------------------
      
      
      !*    1. INITIALISE
      !        ----------
      
      F1LIM = FLMIN / DELTH
      
      STRN(IJ) = 0.0_JWRB
      
      ! ----------------------------------------------------------------------
      
      !*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
      !        ------------------------------------------
      
!$acc loop seq
      DO M=1,NFRE
        XKI = AKI_ICE_SCC(G, WAVNUM(IJ, M), DEPTH(IJ), ROWATER, CITHICK(IJ))
        E = 0.5_JWRB*CITHICK(IJ)*XKI**3 / WAVNUM(IJ, M)
        
        SUME = 0.0_JWRB
!$acc loop seq
        DO K=1,NANG
          SUME = SUME + FL1(IJ, K, M)
        END DO
        
        IF (SUME > F1LIM) THEN
          STRN(IJ) = STRN(IJ) + E**2*SUME*DFIM(M)
        END IF
        
      END DO
      
      
    END DO
!$acc end data
  END SUBROUTINE CIMSSTRN_SCC
END MODULE CIMSSTRN_SCC_MOD
