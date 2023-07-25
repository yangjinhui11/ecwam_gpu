! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE NEWWIND_MOD

CONTAINS
SUBROUTINE NEWWIND_SCC (CDATE, CDATEWH, LLNEWFILE,           &
 &                  CIWA_DPTR,CGROUP_DPTR, FF_NOW_FIELD, FF_NEXT_FIELD)
! ----------------------------------------------------------------------

!**** *NEWWIND* - HANDLING OF WINDS OUTSIDE MULTITASKED AREA    

!     P.A.E.M. JANSSEN  KNMI/ECMWF  SEPTEMBER 1994
!     J. BIDLOT         ECMWF       FEBRUARY 1996  MESSAGE PASSING 
!     S. ABDALLA        ECMWF       OCTOBER 2001   INCLUSION OF AIR
!                                                  DENSITY AND Zi/L
!     J. BIDLOT         ECMWF       AUGUST 2008: CLEAN-UP 

!*    PURPOSE.
!     --------

!       GETS FORCING FIELDS WHEN NEEDED.                                   

!**   INTERFACE.
!     ----------

!       *CALL* *NEWWIND (CDATE, CDATEWH, LLNEWFILE,
!                        WVPRPT, FF_NOW, FF_NEXT,
!      *CDATE*    - START DATE OF SOURCE FUNCTION INTEGRATION
!      *CDATEWH*  - DATE OF THE NEXT FORCING FIELDS
!      *LLNEWFILE*- TRUE IF NEW WIND FILE HAS BEEN OPENED
!      *WVPRPT*   - WAVE PROPERTIES FIELDS
!      *FF_NOW*   - DATA STRUCTURE WITH THE CURRENT FORCING FIELDS
!      *FF_NEXT*  - DATA STRUCTURE WITH THE NEXT FORCING FIELDS

!     METHOD.
!     -------

!       READ NEW WINDS AND ASSOCIATED FIELDS WHEN NEEDED.

!     EXTERNALS.
!     ---------

!       *INCDATE*   - UPDATE DATE TIME GROUP.

!     REFERENCE.
!     ----------

!       NONE

! ----------------------------------------------------------------------

!*    *PARAMETER*  FOR ARRAY DIMENSIONS.

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FREQUENCY, FORCING_FIELDS
      USE YOWFIELD_MOD, ONLY : FORCING_FIELDS_FIELD
      USE YOWPCONS , ONLY : ACD      ,BCD      ,EPSMIN
      USE YOWCOUP  , ONLY : LWCOU
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWPARAM , ONLY : NFRE
      USE YOWPHYS  , ONLY : ALPHA
      USE YOWSTAT  , ONLY : IDELWO
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : CDATEWL  ,CDAWIFL  ,CDATEFL  ,CDTNEXT  ,         &
     &                      NSTORE   ,WSPMIN_RESET_TAUW  ,USTMIN_RESET_TAUW
      USE YOWWNDG  , ONLY : ICODE    ,ICODE_CPL

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      
      USE CIREDUCE_MOD, ONLY: CIREDUCE_SCC
! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "cireduce.intfb.h"
#include "incdate.intfb.h"

      CHARACTER(LEN=14), INTENT(IN)    :: CDATE
      CHARACTER(LEN=14), INTENT(INOUT) :: CDATEWH
      LOGICAL, INTENT(INOUT)           :: LLNEWFILE
      REAL(KIND=JWRB), INTENT(INOUT) :: CIWA_DPTR(:,:,:)  ! CIWA on device Yangjinhui
      REAL(KIND=JWRB), INTENT(INOUT) :: CGROUP_DPTR(:,:,:)  ! point on device Yangjinhui
      TYPE(FORCING_FIELDS_FIELD), INTENT(INOUT)  :: FF_NOW_FIELD
      TYPE(FORCING_FIELDS_FIELD), INTENT(INOUT)     :: FF_NEXT_FIELD


      INTEGER(KIND=JWIM) :: ICODE_WND
      INTEGER(KIND=JWIM) :: ICHNK, KIJS, KIJL, IJ

      REAL(KIND=JWRB) :: WGHT, TLWMAX
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: AIRD_NOW(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WDWAVE_NOW(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CICOVER_NOW(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSWAVE_NOW(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSTAR_NOW(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: UFRIC_NOW(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUW_NOW(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUWDIR_NOW(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: Z0M_NOW(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: Z0B_NOW(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CHRNCK_NOW(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CITHICK_NOW(:,:) => NULL()
  
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: AIRD_NEXT(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WDWAVE_NEXT(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CICOVER_NEXT(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSWAVE_NEXT(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSTAR_NEXT(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: UFRIC_NEXT(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUW_NEXT(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUWDIR_NEXT(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: Z0M_NEXT(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: Z0B_NEXT(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CHRNCK_NEXT(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CITHICK_NEXT(:,:) => NULL()
! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('NEWWIND_SCC',0,ZHOOK_HANDLE)

      IF (LWCOU) THEN
        ICODE_WND = ICODE_CPL
      ELSE
        ICODE_WND = ICODE
      ENDIF

      WGHT = 1.0_JWRB/MAX(WSPMIN_RESET_TAUW,EPSMIN)

      IF (CDATE >= CDATEWH) THEN

!*    2. NEW WIND INPUT.
!        ---------------
        IF (CDATE >= CDATEFL) THEN
            LLNEWFILE = .TRUE.
        ENDIF

!*    2.2 NEW WINDS ARE READ IN.
!         ----------------------

        CDATEWL = CDTNEXT

        CALL GSTATS(1492,0)
    CALL FF_NEXT_FIELD%DEVICE_POINTER(AIRD=AIRD_NEXT,WDWAVE=WDWAVE_NEXT,CICOVER=CICOVER_NEXT, WSWAVE=WSWAVE_NEXT, &
&                          WSTAR=WSTAR_NEXT,UFRIC=UFRIC_NEXT,TAUW=TAUW_NEXT, TAUWDIR=TAUWDIR_NEXT, &
&                          Z0M=Z0M_NEXT, Z0B=Z0B_NEXT,CHRNCK=CHRNCK_NEXT,CITHICK=CITHICK_NEXT)
    CALL FF_NOW_FIELD%DEVICE_POINTER(AIRD=AIRD_NOW,WDWAVE=WDWAVE_NOW,CICOVER=CICOVER_NOW, WSWAVE=WSWAVE_NOW, &
&                          WSTAR=WSTAR_NOW,UFRIC=UFRIC_NOW,TAUW=TAUW_NOW, TAUWDIR=TAUWDIR_NOW, &
&                          Z0M=Z0M_NOW, Z0B=Z0B_NOW,CHRNCK=CHRNCK_NOW,CITHICK=CITHICK_NOW)
!if(associated(wswave_now,target=wswave_next))then
!    print*,"newwind_scc now and next point to the same target"
!endif

!$acc parallel present(AIRD_NEXT,WDWAVE_NEXT,CICOVER_NEXT,WSWAVE_NEXT,WSTAR_NEXT,UFRIC_NEXT,TAUW_NEXT,TAUWDIR_NEXT,Z0M_NEXT,Z0B_NEXT,CHRNCK_NEXT,CITHICK_NEXT,&
!$acc &AIRD_NOW,WDWAVE_NOW,CICOVER_NOW,WSWAVE_NOW,WSTAR_NOW,UFRIC_NOW,TAUW_NOW,TAUWDIR_NOW,Z0M_NOW,Z0B_NOW,CHRNCK_NOW,CITHICK_NOW)
!$acc loop gang
        DO ICHNK = 1, NCHNK
          KIJS = 1
          KIJL = NPROMA_WAM
          IF (ICODE_WND == 3 ) THEN
!$acc loop vector independent
            DO IJ = KIJS, KIJL
              WSWAVE_NOW(IJ,ICHNK) = WSWAVE_NEXT(IJ,ICHNK)
! adapt first estimate of wave induced stress for low winds
! to a fraction of the simple relation u*^2 = Cd(U10) * U10^2
! where this fraction varies from 0 for U10=0 to 1 for U10=WSPMIN_RESET_TAUW
              IF (WSWAVE_NOW(IJ,ICHNK) < WSPMIN_RESET_TAUW) THEN
                TLWMAX = WGHT * (ACD+BCD*WSWAVE_NOW(IJ,ICHNK)) * WSWAVE_NOW(IJ,ICHNK)**3
                TAUW_NOW(IJ,ICHNK) = MIN(TAUW_NOW(IJ,ICHNK), TLWMAX)
              ENDIF
            ENDDO
          ELSE
!$acc loop vector
            DO IJ = KIJS, KIJL
              UFRIC_NOW(IJ,ICHNK) = UFRIC_NEXT(IJ,ICHNK)
! update the estimate of TAUW
              TAUW_NOW(IJ,ICHNK) = UFRIC_NOW(IJ,ICHNK)**2 * &
&(1.0_JWRB-(ALPHA/CHRNCK_NOW(IJ,ICHNK))**2)
! adapt first estimate of wave induced stress for low winds
              IF (UFRIC_NOW(IJ,ICHNK) < USTMIN_RESET_TAUW) TAUW_NOW(IJ,ICHNK) = 0.0_JWRB
            ENDDO
          ENDIF

          WDWAVE_NOW(KIJS:KIJL,ICHNK)  = WDWAVE_NEXT(KIJS:KIJL,ICHNK)
          AIRD_NOW(KIJS:KIJL,ICHNK)    = AIRD_NEXT(KIJS:KIJL,ICHNK)
          WSTAR_NOW(KIJS:KIJL,ICHNK)   = WSTAR_NEXT(KIJS:KIJL,ICHNK)
          CICOVER_NOW(KIJS:KIJL,ICHNK) = CICOVER_NEXT(KIJS:KIJL,ICHNK)
          CITHICK_NOW(KIJS:KIJL,ICHNK) = CITHICK_NEXT(KIJS:KIJL,ICHNK)

        ENDDO
!$acc end parallel 
        CALL GSTATS(1492,1)


        CALL INCDATE(CDATEWH, IDELWO)

!       UPDATE THE SEA ICE REDUCTION FACTOR
        CALL CIREDUCE_SCC (CIWA_DPTR,CGROUP_DPTR,CICOVER_NOW,CITHICK_NOW)

      ENDIF

IF (LHOOK) CALL DR_HOOK('NEWWIND_SCC',1,ZHOOK_HANDLE)

END SUBROUTINE NEWWIND_SCC
END MODULE
