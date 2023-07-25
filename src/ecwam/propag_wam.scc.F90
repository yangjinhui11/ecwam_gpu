! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
MODULE PROPAG_WAM_MOD

#include "wam_user_clock.intfb.h"
CONTAINS
SUBROUTINE PROPAG_WAM_SCC (BLK2GLO, WVENVI, WVPRPT, FL1)

! ----------------------------------------------------------------------

!**** *PROPAG_WAM* - WAVE PROPGATION

!*    PURPOSE.
!     --------

!     PROPAGATION

!**   INTERFACE.
!     ----------

!     *CALL* *PROPAG_WAM (BLK2GLO, WVENVI, WVPRPT, FL1)
!          *BLK2GLO*   - BLOCK TO GRID TRANSFORMATION
!          *WVENVI*    - WAVE ENVIRONMENT FIELDS
!          *WVPRPT*    - WAVE PROPERTIES FIELDS
!          *FL1*       - SPECTRUM


!     METHOD.
!     -------

!     EXTERNALS.
!     ----------


! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO!, ENVIRONMENT, FREQUENCY
      USE YOWFIELD_MOD, ONLY : FREQUENCY_FIELD, ENVIRONMENT_FIELD
      USE YOWCURR  , ONLY : LLCHKCFL ,LLCHKCFLA
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, KIJL4CHNK, IJFROMCHNK
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED ,NIBLO , LLUNSTR
      USE YOWREFD  , ONLY : LLUPDTTD ,THDD     ,THDC     ,SDOT
      USE YOWSTAT  , ONLY : IPROPAGS ,IFRELFMAX, DELPRO_LF, IDELPRO
      USE YOWUBUF  , ONLY : LUPDTWGHT
#ifdef WAM_HAVE_UNWAM
      USE UNWAM    , ONLY : PROPAG_UNWAM
#endif
      USE YOWTEST, ONLY: IU06
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
    
! Yangjinhui 
      USE MPEXCHNG_MOD, ONLY: MPEXCHNG_SCC
      USE PROENVHALO_MOD, ONLY: PROENVHALO_SCC
      USE PROPDOT_MOD   , ONLY: PROPDOT_SCC
      USE CTUWUPDT_MOD  , ONLY:CTUWUPDT_SCC
      USE PROPAGS2_MOD  , ONLY:PROPAGS2_SCC
      USE PROPAGS1_MOD  , ONLY:PROPAGS1_SCC
      USE PROPAGS_MOD   , ONLY:PROPAGS_SCC
      USE GPUOBJ_MOD   ,  ONLY:get_mem_info
! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"

      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      TYPE(ENVIRONMENT_FIELD), INTENT(INOUT) :: WVENVI
      TYPE(FREQUENCY_FIELD), INTENT(INOUT) :: WVPRPT
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1


      INTEGER(KIND=JWIM) :: IJ, K, M, J,TIME0
      INTEGER(KIND=JWIM) :: JKGLO, NPROMA, MTHREADS
      INTEGER(KIND=JWIM) :: NSTEP_LF, ISUBST
!$    INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS
      INTEGER(KIND=JWIM) :: IJSG, IJLG, ICHNK, KIJS, KIJL, IJSB, IJLB

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     Spectra extended with the halo exchange for the propagation
!     But limited to NFRE_RED frequencies
!!! the advection schemes are still written in block structure
      REAL(KIND=JWRB),ALLOCATABLE :: FL1_EXT(:, :, :), FL3_EXT(:,:,:)
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED) :: WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1) :: DELLAM1_EXT, COSPHM1_EXT 
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1) :: DEPTH_EXT, UCUR_EXT, VCUR_EXT

      LOGICAL :: L1STCALL

      DATA L1STCALL / .TRUE. /

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PROPAG_WAM',0,ZHOOK_HANDLE)

ALLOCATE(FL1_EXT(NINF:NSUP+1, NANG, NFRE_RED))
ALLOCATE(FL3_EXT(NINF:NSUP+1, NANG, NFRE_RED))
      IF (NIBLO > 1) THEN

        IJSG = IJFROMCHNK(1,1)
        IJLG = IJSG + SUM(KIJL4CHNK) - 1

        MTHREADS=1
!$      MTHREADS=OMP_GET_MAX_THREADS()
        NPROMA=(IJLG-IJSG+1)/MTHREADS + 1
        NPROMA=NPROMA_WAM ! Yangjinhui

!!! the advection schemes are still written in block structure
!!! mapping chuncks to block ONLY for actual grid points !!!!
!$acc enter data create(FL1_EXT,FL3_EXT)
!$acc parallel loop gang present(IJFROMCHNK,KIJL4CHNK,FL1_EXT,FL1)
        DO ICHNK = 1, NCHNK
          KIJS = 1
          IJSB = IJFROMCHNK(KIJS, ICHNK)
          KIJL = KIJL4CHNK(ICHNK)
          IJLB = IJFROMCHNK(KIJL, ICHNK)
!$acc loop vector collapse(2)
          DO M = 1, NFRE_RED
            DO K = 1, NANG
              FL1_EXT(IJSB:IJLB, K, M) = FL1(KIJS:KIJL, K, M, ICHNK)
            ENDDO
          ENDDO
        ENDDO
!       SET THE DUMMY LAND POINT TO 0.
!$acc kernels present(FL1_EXT)
        FL1_EXT(NSUP+1,:,:) = 0.0_JWRB 
!$acc end kernels

        IF (LLUNSTR) THEN 

        WRITE(0,*) '!!! ********************************* !!'
        WRITE(0,*) '!!! in PROPAG_WAM Not yet ready !!!' 
        WRITE(0,*) '!!! PROPAG_UNWAM will need to be adapted to new FL1 structure !!!' 
        WRITE(0,*) '!!! ********************************* !!'
        CALL ABORT1

!!!!!           CALL PROPAG_UNWAM(FL1_EXT, FL1)

        ELSE


!          OBTAIN INFORMATION AT NEIGHBORING GRID POINTS (HALO)
!          ----------------------------------------------------
           CALL MPEXCHNG_SCC(FL1_EXT, NANG, 1, NFRE_RED)


           CALL GSTATS(1430,0)

!          UPDATE DOT THETA TERM
!          ---------------------
!$acc  enter data create(WAVNUM_EXT,CGROUP_EXT,OMOSNH2KD_EXT,DELLAM1_EXT,COSPHM1_EXT,DEPTH_EXT,UCUR_EXT,VCUR_EXT) 
           IF (LLUPDTTD) THEN
             IF (.NOT.ALLOCATED(THDC)) THEN
                    ALLOCATE(THDC(IJSG:IJLG, NANG))
                    ALLOCATE(THDD(IJSG:IJLG, NANG))
                    ALLOCATE(SDOT(IJSG:IJLG, NANG, NFRE_RED))
             ENDIF
!            NEED HALO VALUES
             CALL  PROENVHALO_SCC (NINF, NSUP,                            &
&                              WVPRPT,                                &
&                              WVENVI,                                &
&                              WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
&                              DELLAM1_EXT, COSPHM1_EXT,              & 
&                              DEPTH_EXT, UCUR_EXT, VCUR_EXT )


!            DOT THETA TERM:

!$acc parallel loop gang &
!$acc &present(BLK2GLO,WAVNUM_EXT,CGROUP_EXT,OMOSNH2KD_EXT,COSPHM1_EXT,DEPTH_EXT,&
!$acc &UCUR_EXT,VCUR_EXT,THDC,THDD,SDOT)
             DO JKGLO = IJSG, IJLG, NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1, IJLG)
               CALL PROPDOT_SCC(KIJS, KIJL, NINF, NSUP,                                     &
     &                      BLK2GLO,                                                    &
     &                      WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT,                      &
     &                      COSPHM1_EXT, DEPTH_EXT, UCUR_EXT, VCUR_EXT,                 &
     &                      THDC(KIJS:KIJL,:), THDD(KIJS:KIJL,:), SDOT(KIJS:KIJL,:,:))
             ENDDO

             LLUPDTTD = .FALSE.
           ENDIF

!          IPROPAGS = 2 is the default option (the other options are kept but usually not used)
!          -----------------------------------------------------------------------------------
           SELECT CASE (IPROPAGS)

           CASE(2)

             IF (LUPDTWGHT) THEN
!              NEED HALO VALUES
               CALL  PROENVHALO_SCC (NINF, NSUP,                            &
&                                WVPRPT,                                &
&                                WVENVI,                                &
&                                WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
&                                DELLAM1_EXT, COSPHM1_EXT,              & 
&                                DEPTH_EXT, UCUR_EXT, VCUR_EXT )

!              COMPUTES ADVECTION WEIGTHS AND CHECK CFL CRITERIA
               CALL CTUWUPDT_SCC(IJSG, IJLG, NINF, NSUP,                      &
&                            BLK2GLO,                                     &
&                            CGROUP_EXT, OMOSNH2KD_EXT,                   &
&                            COSPHM1_EXT, DEPTH_EXT, UCUR_EXT, VCUR_EXT )

               LUPDTWGHT=.FALSE.
             ENDIF
!$acc parallel loop gang present(FL1_EXT,FL3_EXT)
             DO JKGLO = IJSG, IJLG, NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1, IJLG)
               CALL PROPAGS2_SCC(FL1_EXT, FL3_EXT, NINF, NSUP, KIJS, KIJL, NANG, 1, NFRE_RED)
             ENDDO

!            SUB TIME STEPPING FOR FAST WAVES (only if IFRELFMAX > 0)
             IF (IFRELFMAX > 0 ) THEN
               NSTEP_LF = NINT(REAL(IDELPRO, JWRB)/DELPRO_LF)
               ISUBST = 2  ! The first step was done as part of the previous call to PROPAGS2

               DO WHILE (ISUBST <= NSTEP_LF)
!$acc parallel loop gang present(FL1_EXT,FL3_EXT)
                 DO JKGLO = IJSG, IJLG, NPROMA
                   KIJS=JKGLO
                   KIJL=MIN(KIJS+NPROMA-1, IJLG)
!$acc loop vector collapse(3)
                   DO M = 1, IFRELFMAX 
                     DO K = 1, NANG
                       DO IJ = KIJS, KIJL
                         FL1_EXT(IJ, K, M) = FL3_EXT(IJ, K, M)
                       ENDDO
                     ENDDO
                   ENDDO
                 ENDDO
                 CALL MPEXCHNG_SCC(FL1_EXT(:,:,1:IFRELFMAX), NANG, 1, IFRELFMAX)
!$acc parallel loop gang present(FL1_EXT,FL3_EXT)
                 DO JKGLO = IJSG, IJLG, NPROMA
                   KIJS=JKGLO
                   KIJL=MIN(KIJS+NPROMA-1, IJLG)
                   CALL PROPAGS2_SCC(FL1_EXT(:,:,1:IFRELFMAX), FL3_EXT(:,:,1:IFRELFMAX), NINF, NSUP, KIJS, KIJL, NANG, 1, IFRELFMAX)
                 ENDDO

                 ISUBST = ISUBST + 1

               ENDDO
             ENDIF  ! end sub time steps (if needed)

           CASE(1)
             IF (L1STCALL .OR. LLCHKCFLA) LLCHKCFL=.TRUE.

!            NEED HALO VALUES
             CALL  PROENVHALO_SCC (NINF, NSUP,                            &
&                              WVPRPT,                                &
&                              WVENVI,                                &
&                              WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
&                              DELLAM1_EXT, COSPHM1_EXT,              & 
&                              DEPTH_EXT, UCUR_EXT, VCUR_EXT )

!$acc parallel loop gang &
!$acc &present(FL1_EXT,FL3_EXT,BLK2GLO,DEPTH_EXT,CGROUP_EXT,OMOSNH2KD_EXT,&
!$acc &DELLAM1_EXT,COSPHM1_EXT,UCUR_EXT,VCUR_EXT)
             DO JKGLO = IJSG, IJLG, NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1, IJLG)
               CALL PROPAGS1_SCC(FL1_EXT, FL3_EXT, NINF, NSUP, KIJS, KIJL,       &
&                            BLK2GLO,                                        &
&                            DEPTH_EXT,                                      &
&                            CGROUP_EXT, OMOSNH2KD_EXT,                      &
&                            DELLAM1_EXT, COSPHM1_EXT,                       &
&                            UCUR_EXT, VCUR_EXT,                             &
&                            L1STCALL)
             ENDDO

           CASE(0)
             IF (L1STCALL .OR. LLCHKCFLA) LLCHKCFL=.TRUE.

!            NEED HALO VALUES
             CALL  PROENVHALO_SCC (NINF, NSUP,                            &
&                              WVPRPT,                                &
&                              WVENVI,                                &
&                              WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
&                              DELLAM1_EXT, COSPHM1_EXT,              & 
&                              DEPTH_EXT, UCUR_EXT, VCUR_EXT )

!$acc parallel loop gang present(&
!$acc& FL1_EXT,FL3_EXT,BLK2GLO,DEPTH_EXT,CGROUP_EXT,OMOSNH2KD_EXT,&
!$acc& DELLAM1_EXT,COSPHM1_EXT,UCUR_EXT,VCUR_EXT)
             DO JKGLO = IJSG, IJLG, NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1, IJLG)
               CALL PROPAGS_SCC(FL1_EXT, FL3_EXT, NINF, NSUP, KIJS, KIJL,       &
&                           BLK2GLO,                                        &
&                           DEPTH_EXT,                                      &
&                           CGROUP_EXT, OMOSNH2KD_EXT,                      &
&                           DELLAM1_EXT, COSPHM1_EXT,                       &
&                           UCUR_EXT, VCUR_EXT,                             &
&                           L1STCALL)
             ENDDO
           END SELECT 

!!! the advection schemes are still written in block structure
!!!  So need to convert back to the nproma_wam chuncks
!$acc parallel present(FL1,FL3_EXT,IJFROMCHNK,KIJL4CHNK)
!$acc loop gang
          DO ICHNK = 1, NCHNK
            KIJS = 1
            IJSB = IJFROMCHNK(KIJS, ICHNK)
            KIJL = KIJL4CHNK(ICHNK)
            IJLB = IJFROMCHNK(KIJL, ICHNK)
!$acc loop vector collapse(2)
            DO M = 1, NFRE_RED
              DO K = 1, NANG
                FL1(KIJS:KIJL, K, M, ICHNK) = FL3_EXT(IJSB:IJLB, K, M)
              ENDDO
            ENDDO
          ENDDO
! Yangjinhui

!$acc loop gang
          DO ICHNK = 1, NCHNK
            KIJS = 1
            IJSB = IJFROMCHNK(KIJS, ICHNK)
            KIJL = KIJL4CHNK(ICHNK)
            IJLB = IJFROMCHNK(KIJL, ICHNK)
            
            IF (KIJL < NPROMA_WAM) THEN
              !!! make sure fictious points keep values of the first point in the chunk
!$acc loop vector collapse(2)
              DO M = 1, NFRE_RED
                DO K = 1, NANG
                  FL1(KIJL+1:NPROMA_WAM, K, M, ICHNK) = FL1(1, K, M, ICHNK)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
!$acc end parallel


           CALL GSTATS(1430,1)

        ENDIF  ! end propagation

      ENDIF ! more than one grid point

      L1STCALL=.FALSE.
      LLCHKCFL=.FALSE.
! Yangjinhui
!$acc exit data delete(FL1_EXT,FL3_EXT)
!$acc exit data delete(WAVNUM_EXT,CGROUP_EXT,OMOSNH2KD_EXT,DELLAM1_EXT,COSPHM1_EXT,DEPTH_EXT,UCUR_EXT,VCUR_EXT)
    !IF (ALLOCATED(THDD)) THEN
    !   DEALLOCATE(THDD)
    !    DEALLOCATE(THDC)
    !    DEALLOCATE(SDOT)
    !ENDIF
IF (LHOOK) CALL DR_HOOK('PROPAG_WAM',1,ZHOOK_HANDLE)

END SUBROUTINE PROPAG_WAM_SCC
END MODULE
