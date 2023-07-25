! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE CIREDUCE_MOD

CONTAINS

SUBROUTINE CIREDUCE_SCC (CIWA_DPTR,CGROUP_DPTR, CICOVER_DPTR,CITHICK_DPTR)

! ----------------------------------------------------------------------

!**** *CIREDUCE* - COMPUTE SEA ICE REDUCTION FACTOR FOR SOURCE TERMS 
!                  AND THE SEA ICE WAVE ATTENUATION FACTORS

!           IF THERE IS NO SEA ICE INFORMATION OR
!           ALL SEA ICE COVER POINTS WILL BE MASKED
!           THEN CIWA WILL BE SET ON THE FIRST CALL. NOTHING WILL BE DONE
!           IN ALL FOLLOWING CALLS

!!!! currently also setting parametric sea ice thickness !!!!

!*    PURPOSE.
!     --------

!       CIREDUCE COMPUTES SEA ICE SOURCE TERM REDUCTION FACTOR.

!**   INTERFACE.
!     ----------

!       *CALL* *CIREDUCE (CGROUP, CICOVER, CITHICK, CIWA)

!          *CGROUP*  - GROUP SPEED.
!          *CICOVER* - SEA ICE COVER.
!          *CITHICK* - SEA ICE THICKNESS. 
!          *CIWA*-     SEA ICE WAVE ATTENUATION FACTOR. 

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE 
      USE YOWPARAM , ONLY : NFRE

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOWDRVTYPE ,ONLY: FREQUENCY, FORCING_FIELDS
      
      USE CIWAF_MOD  ,ONLY: CIWAF_SCC
! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "ciwaf.intfb.h"

      REAL(KIND=JWRB),INTENT(INOUT):: CIWA_DPTR(:,:,:)
      REAL(KIND=JWRB),INTENT(IN):: CGROUP_DPTR(:,:,:)
      REAL(KIND=JWRB),INTENT(IN):: CICOVER_DPTR(:,:)
      REAL(KIND=JWRB),INTENT(IN):: CITHICK_DPTR(:,:)

      INTEGER(KIND=JWIM) :: IJ, M 
      INTEGER(KIND=JWIM) :: ICHNK

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      LOGICAL, SAVE :: LLFRST

      DATA LLFRST / .TRUE. /

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CIREDUCE_SCC',0,ZHOOK_HANDLE)


        IF( .NOT. LICERUN .OR. LMASKICE ) THEN

          IF (LLFRST) THEN
            LLFRST=.FALSE.
!           NO REDUCTION, EITHER THERE IS NO SEA ICE INFORMATION OR
!           ALL SEA ICE COVER POINTS WILL BE MASKED
            CALL GSTATS(1493,0)
!$acc kernels present(CIWA_DPTR) 
            DO ICHNK = 1, NCHNK
               CIWA_DPTR(:,:,ICHNK) = 1.0_JWRB
            ENDDO
!$acc end kernels
            CALL GSTATS(1493,1)
          ENDIF

        ELSE

          CALL GSTATS(1493,0)
!         DETERMINE THE WAVE ATTENUATION FACTOR
!$acc parallel loop gang present(CGROUP_DPTR,CICOVER_DPTR,CITHICK_DPTR,CIWA_DPTR)
          DO ICHNK = 1, NCHNK
            CALL CIWAF_SCC(1, NPROMA_WAM, CGROUP_DPTR(:,:,ICHNK), CICOVER_DPTR(:,ICHNK), &
&                      CITHICK_DPTR(:,ICHNK), CIWA_DPTR(:,:,ICHNK))
          ENDDO
          CALL GSTATS(1493,1)
        ENDIF

IF (LHOOK) CALL DR_HOOK('CIREDUCE_SCC',1,ZHOOK_HANDLE)

END SUBROUTINE CIREDUCE_SCC

END MODULE
