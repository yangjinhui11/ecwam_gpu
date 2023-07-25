! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE CHNKMIN_SCC_MOD
  CONTAINS
  FUNCTION CHNKMIN_SCC (U10)
    
    ! ----------------------------------------------------------------------
    
    !**** *CHNKMIN* - FUNCTION TO COMPUTE THE MINMUM CHARNOCK
    
    !*    PURPOSE.
    !     -------
    
    
    !**   INTERFACE.
    !     ----------
    
    !       *FUNCTION* *CHNKMIN (U10)*
    
    !     METHOD.
    !     -------
    
    !     CHNKMIN = ALPHAMIN + (ALPHA-ALPHAMIN)*0.5_JWRB*(1.0_JWRB-TANH(U10-A))
    
    !     EXTERNALS.
    !     ----------
    
    !       NONE.
    
    !     REFERENCE.
    !     ----------
    
    !       NONE.
    
    ! ----------------------------------------------------------------------
    
    USE PARKIND_WAVE, ONLY: JWIM, JWRB, JWRU
    
    USE YOWPHYS, ONLY: ALPHA, ALPHAMIN, CHNKMIN_U
    ! ----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    REAL(KIND=JWRB) :: CHNKMIN_SCC
    REAL(KIND=JWRB), INTENT(IN) :: U10
    
    ! ----------------------------------------------------------------------
    
!$acc routine seq
    
    CHNKMIN_SCC = ALPHAMIN + (ALPHA - ALPHAMIN)*0.5_JWRB*(1.0_JWRB - TANH(U10 - CHNKMIN_U))
    
    
  END FUNCTION CHNKMIN_SCC
END MODULE CHNKMIN_SCC_MOD
