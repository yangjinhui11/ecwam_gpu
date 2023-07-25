! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE STOKESTRN_SCC_MOD
  CONTAINS
  SUBROUTINE STOKESTRN_SCC (KIJS, KIJL, FL1, WAVNUM, STOKFAC, DEPTH, WSWAVE, WDWAVE, CICOVER, CITHICK, USTOKES, VSTOKES, STRNMS,  &
  & NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN)
    
    ! ----------------------------------------------------------------------
    
    !**** *STOKESTRN* - WRAPPER TO CALL STOKESDRIFT and CIMSSTRN
    
    !*    PURPOSE.
    !     --------
    
    !**   INTERFACE.
    !     ----------
    
    !       *CALL* *STOKESTRN (KIJS, KIJL, FL1, WAVNUM, STOKFAC, DEPTH, FF_NOW, INTFLDS, WAM2NEMO)
    
    !          *KIJS*    - INDEX OF FIRST GRIDPOINT.
    !          *KIJL*    - INDEX OF LAST GRIDPOINT.
    !          *FL1*     - SPECTRUM(INPUT).
    !          *WAVNUM*  - WAVE NUMBER.
    !          *STOKFAC* - STOKES DRIFT FACTOR.
    !          *DEPTH*   - WATER DEPTH.
    !          *FF_NOW*  - FORCING FIELDS AT CURRENT TIME.
    !          *INTFLDS* - INTEGRATED/DERIVED PARAMETERS
    !          *WAM2NEMO*- WAVE FIELDS PASSED TO NEMO
    
    ! ----------------------------------------------------------------------
    
    USE STOKESDRIFT_SCC_MOD, ONLY: STOKESDRIFT_SCC
    USE CIMSSTRN_SCC_MOD, ONLY: CIMSSTRN_SCC
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU, JWRO
    USE YOWDRVTYPE, ONLY: FORCING_FIELDS, INTGT_PARAM_FIELDS, WAVE2OCEAN
    
    USE YOWCOUP, ONLY: LWCOU, LWNEMOCOU, LWNEMOCOUSEND, LWNEMOCOUSTK, LWNEMOCOUSTRN
    USE YOWPARAM, ONLY: NANG, NFRE
    
    
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
    
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NANG, NFRE) :: FL1
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL, NFRE) :: WAVNUM, STOKFAC
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: DEPTH
    REAL(KIND=JWRB), INTENT(IN), DIMENSION(KIJL) :: WSWAVE, WDWAVE, CICOVER, CITHICK
    REAL(KIND=JWRB), INTENT(INOUT), DIMENSION(KIJL) :: USTOKES, VSTOKES, STRNMS
    REAL(KIND=JWRO), INTENT(INOUT), DIMENSION(KIJL) :: NEMOUSTOKES, NEMOVSTOKES, NEMOSTRN
    
    
    INTEGER :: IJ
!$acc routine vector
!$acc data present( FL1, WAVNUM, STOKFAC, DEPTH, WSWAVE, WDWAVE, CICOVER, CITHICK, USTOKES, VSTOKES, STRNMS, NEMOUSTOKES,  &
!$acc & NEMOVSTOKES, NEMOSTRN )
    
    ! ----------------------------------------------------------------------
    
    
    CALL STOKESDRIFT_SCC(KIJS, KIJL, FL1, STOKFAC, WSWAVE, WDWAVE, CICOVER, USTOKES, VSTOKES)
    
    IF (LWNEMOCOUSTRN)     CALL CIMSSTRN_SCC(KIJS, KIJL, FL1, WAVNUM, DEPTH, CITHICK, STRNMS)
    
!$acc loop vector
    DO IJ=KIJS,KIJL
      
      
      IF (LWNEMOCOU .and. (LWNEMOCOUSEND .and. LWCOU .or. .not.LWCOU)) THEN
        IF (LWNEMOCOUSTK) THEN
          NEMOUSTOKES(IJ) = USTOKES(IJ)
          NEMOVSTOKES(IJ) = VSTOKES(IJ)
        ELSE
          NEMOUSTOKES(IJ) = 0.0_JWRO
          NEMOVSTOKES(IJ) = 0.0_JWRO
        END IF
        
        IF (LWNEMOCOUSTRN)         NEMOSTRN(IJ) = STRNMS(IJ)
      END IF
      
      
      ! ----------------------------------------------------------------------
      
    END DO
!$acc end data
  END SUBROUTINE STOKESTRN_SCC
END MODULE STOKESTRN_SCC_MOD
