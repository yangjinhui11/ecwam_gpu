! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE WAMINTGR (CDTPRA, CDATE, CDATEWH, CDTIMP, CDTIMPNEXT,  &
 &                   BLK2GLO,                                     &
 &                   WVENVI, WVPRPT, FF_NOW, FF_NEXT, INTFLDS,    &
 &                   WAM2NEMO, MIJ, FL1, XLLWS, TIME1)


! ----------------------------------------------------------------------

!**** *WAMINTGR* - 3-G WAM MODEL - TIME INTEGRATION OF WAVE FIELDS.

!*    PURPOSE.
!     --------

!       COMPUTATION OF THE 2-D FREQUENCY-DIRECTION WAVE SPECTRUM AT ALL
!       GRID POINTS FOR A GIVEN INITIAL SPECTRUM AND FORCING SURFACE
!       STRESS FIELD.

!     REFERENCE.
!     ----------

!         IFS DOCUMENTATION, part VII 

! -------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
USE YOWDRVTYPE  , ONLY : WVGRIDGLO, ENVIRONMENT, FREQUENCY, FORCING_FIELDS,  &
 &                       INTGT_PARAM_FIELDS, WAVE2OCEAN

USE YOWCOUP  , ONLY : LWNEMOCOU, NEMONTAU
USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
USE YOWPARAM , ONLY : NIBLO, NANG, NFRE
USE YOWPCONS , ONLY : EPSMIN
USE YOWSTAT  , ONLY : CDTPRO, IDELPRO, IDELT, IDELWI, LLSOURCE
USE YOWWIND  , ONLY : CDAWIFL, CDATEWO, CDATEFL
USE YOWFIELD_MOD, ONLY : FREQUENCY_FIELD, ENVIRONMENT_FIELD, FORCING_FIELDS_FIELD, &
 &  WAVE2OCEAN_FIELD, INTGT_PARAM_FIELDS_FIELD, SOURCE_CONTRIBS_FIELD
USE YOWMPP   , ONLY : IRANK
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
#ifdef WAM_PHYS_GPU
USE IMPLSCH_SCC_MOD, ONLY: implsch_SCC,INIT_IMPLSCH
#endif

!.....Adding global variables to driver symbol table for offload instructions
USE YOWINDN, ONLY: dal2, k11w, inlcoef, k21w, fklam, fklam1, mfrstlw, ikp1, kfrh, af11, ikm1, dal1, mlsthg, k1w, ikm, fklap,  &
& fklap1, k2w, ikp, rnlcoef
USE YOWCOUP, ONLY: lwnemotauoc, lwnemocoustrn, lwnemocoustk, lwcou, wtauhf, x0tauhf, llgcbz0, lwnemocousend, lwvflx_snl,  &
& lwflux, llnormagam, llcapchnk
USE YOWPHYS, ONLY: satweights, cdisvis, tauwshelter, rnum, chnkmin_u, z0rat, dthrn_u, nsdsnth, z0tubmax, alpha, alphapmax,  &
& alphamin, bmaxokap, gamnconst, indicessat, delta_sdis, swellf5, ndikcumul, cdis, ang_gc_c, cumulw, ang_gc_a, tailfactor,  &
& ang_gc_b, dthrn_a, rnu, betamaxoxkappa2, zalp, tailfactor_pm, rn1_rn
USE YOWPCONS, ONLY: zpi4gm2, zpi4gm1, gm1, sqrtgosurft, zpi, g
USE YOWFRED, ONLY: omxkm3_gc, zpifr, dfim_sim, sinth, delkcc_gc_ns, dfimfr2, flogsprdm1, costh, xkm_gc, c2osqrtvg_gc, xk_gc,  &
& xkmsqrtvgoc2_gc, delkcc_omxkm3_gc, fr5, fr, om3gmkm_gc, cofrm4, omega_gc, rhowg_dfim, cm_gc, flmax, th, nwav_gc, dfim,  &
& delth, dfimofr, dfimfr
USE YOWPARAM, ONLY: nfre_odd, llunstr, nfre_red
USE YOWICE, ONLY: lmaskice, lwamrsetci, ciblock, cdicwa, cithrsh, licerun, cithrsh_tail, lciwabr
USE YOWWNDG, ONLY: icode, icode_cpl
USE YOWSTAT, ONLY: lbiwbk, idamping, isnonlin, iphys
USE YOWALTAS, ONLY: bfcrv, egrcrv, afcrv
USE YOWCOUT, ONLY: lwfluxout
USE YOWWIND, ONLY: wspmin
USE YOWTABL, ONLY: swellft
USE YOWTEST, ONLY: IU06      
USE GPUOBJ_MOD, ONLY: WVPRPT_FIELD,WVENVI_FIELD,FF_NOW_FIELD,WAM2NEMO_FIELD,INTFLDS_FIELD,SRC_CONTRIBS ,&
&get_device_info, FF_NEXT_FIELD,get_mem_info

USE PROPAG_WAM_MOD, ONLY: PROPAG_WAM_SCC
USE NEWWIND_MOD , ONLY: NEWWIND_SCC
! ----------------------------------------------------------------------

IMPLICIT NONE

#ifndef WAM_PHYS_GPU
#include "implsch.intfb.h"
#endif
#include "incdate.intfb.h"
#include "newwind.intfb.h"
#include "propag_wam.intfb.h"
#include "wam_user_clock.intfb.h"

CHARACTER(LEN=14), INTENT(IN)                                            :: CDTPRA  ! DATE FOR CALL PROPAGATION
CHARACTER(LEN=14), INTENT(INOUT)                                         :: CDATE  ! CURRENT DATE
CHARACTER(LEN=14), INTENT(INOUT)                                         :: CDATEWH ! DATE OF THE NEXT FORCING FIELDS
CHARACTER(LEN=14), INTENT(INOUT)                                         :: CDTIMP  ! START DATE OF SOURCE FUNCTION INTEGRATION
CHARACTER(LEN=14), INTENT(INOUT)                                         :: CDTIMPNEXT  ! NEXT START DATE OF SOURCE FUNCTION INTEGRATION
TYPE(WVGRIDGLO), INTENT(IN)                                              :: BLK2GLO  ! BLOCK TO GRID TRANSFORMATION
TYPE(ENVIRONMENT), INTENT(INOUT)                                         :: WVENVI !  WAVE ENVIRONMENT FIELDS
TYPE(FREQUENCY), INTENT(INOUT)                                           :: WVPRPT  ! WAVE PROPERTIES FIELDS
TYPE(FORCING_FIELDS), INTENT(INOUT)                                      :: FF_NOW  ! FORCING FIELDS AT CURRENT TIME
TYPE(FORCING_FIELDS), INTENT(IN)                                         :: FF_NEXT  !  DATA STRUCTURE WITH THE NEXT FORCING FIELDS
TYPE(INTGT_PARAM_FIELDS), INTENT(INOUT)                                  :: INTFLDS  ! INTEGRATED/DERIVED PARAMETERS
TYPE(WAVE2OCEAN), INTENT(INOUT)                                          :: WAM2NEMO  ! WAVE FIELDS PASSED TO NEMO
INTEGER(KIND=JWIM), DIMENSION(NPROMA_WAM, NCHNK), INTENT(INOUT)          :: MIJ  ! LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE
REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1
REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: XLLWS ! TOTAL WINDSEA MASK FROM INPUT SOURCE TERM

REAL(KIND=JWRB), INTENT(INOUT) :: TIME1(3)
REAL(KIND=JWRB) :: TIME0, TIME2, T1


INTEGER(KIND=JWIM) :: IJ, K, M
INTEGER(KIND=JWIM) :: ICHNK
INTEGER(KIND=JWIM) :: IDELWH

! Objects to store fields
!TYPE(FREQUENCY_FIELD) :: WVPRPT_FIELD
!TYPE(ENVIRONMENT_FIELD) :: WVENVI_FIELD
!TYPE(FORCING_FIELDS_FIELD) :: FF_NOW_FIELD
!TYPE(WAVE2OCEAN_FIELD) :: WAM2NEMO_FIELD
!TYPE(INTGT_PARAM_FIELDS_FIELD) :: INTFLDS_FIELD
!TYPE(SOURCE_CONTRIBS_FIELD) :: SRC_CONTRIBS

#ifdef WAM_PHYS_GPU
  ! DEVICE POINTERS
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: FL1_DPTR(:,:,:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: XLLWS_DPTR(:,:,:,:) => NULL()
  INTEGER(KIND=JWIM), POINTER, CONTIGUOUS :: MIJ_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WAVNUM_DPTR(:,:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CGROUP_DPTR(:,:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CIWA_DPTR(:,:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CINV_DPTR(:,:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: XK2CG_DPTR(:,:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: STOKFAC_DPTR(:,:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: EMAXDPT_DPTR(:,:) => NULL()
  INTEGER(KIND=JWIM), POINTER, CONTIGUOUS :: INDEP_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: DEPTH_DPTR(:,:) => NULL()
  INTEGER(KIND=JWIM), POINTER, CONTIGUOUS :: IOBND_DPTR(:,:) => NULL()
  INTEGER(KIND=JWIM), POINTER, CONTIGUOUS :: IODP_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CICOVER_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSWAVE_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WDWAVE_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: AIRD_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSTAR_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: UFRIC_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUW_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUWDIR_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: Z0M_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: Z0B_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CHRNCK_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CITHICK_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NEMOUSTOKES_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NEMOVSTOKES_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NEMOSTRN_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NPHIEPS_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NTAUOC_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NSWH_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NMWP_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NEMOTAUX_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NEMOTAUY_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NEMOWSWAVE_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: NEMOPHIF_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSEMEAN_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSFMEAN_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: USTOKES_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: VSTOKES_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: STRNMS_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUXD_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUYD_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUOCXD_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUOCYD_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUOC_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: PHIOCD_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: PHIEPS_DPTR(:,:) => NULL()
  REAL(KIND=JWRB), POINTER, CONTIGUOUS :: PHIAW_DPTR(:,:) => NULL()
#endif

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


LOGICAL, SAVE :: LLNEWFILE,FIRST_TIME

DATA LLNEWFILE / .FALSE. /
DATA FIRST_TIME / .TRUE. /

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WAMINTGR',0,ZHOOK_HANDLE)

!*     PROPAGATION TIME
!      ----------------

TIME0=-WAM_USER_CLOCK()
CALL WVPRPT_FIELD%DEVICE_POINTER(WAVNUM=WAVNUM_DPTR, CGROUP=CGROUP_DPTR, CIWA=CIWA_DPTR, CINV=CINV_DPTR, XK2CG=XK2CG_DPTR,  &
      & STOKFAC=STOKFAC_DPTR)
CALL WVENVI_FIELD%DEVICE_POINTER(EMAXDPT=EMAXDPT_DPTR, INDEP=INDEP_DPTR, DEPTH=DEPTH_DPTR, IOBND=IOBND_DPTR, IODP=IODP_DPTR)
CALL FF_NOW_FIELD%DEVICE_POINTER(AIRD=AIRD_DPTR, WDWAVE=WDWAVE_DPTR, CICOVER=CICOVER_DPTR, WSWAVE=WSWAVE_DPTR,  &
      & WSTAR=WSTAR_DPTR, UFRIC=UFRIC_DPTR, TAUW=TAUW_DPTR, TAUWDIR=TAUWDIR_DPTR, Z0M=Z0M_DPTR, Z0B=Z0B_DPTR,  &
      & CHRNCK=CHRNCK_DPTR, CITHICK=CITHICK_DPTR)
CALL SRC_CONTRIBS%DEVICE_POINTER(FL1=FL1_DPTR, XLLWS=XLLWS_DPTR, MIJ=MIJ_DPTR)
      
CALL WAM2NEMO_FIELD%DEVICE_POINTER(NEMOUSTOKES=NEMOUSTOKES_DPTR, NEMOVSTOKES=NEMOVSTOKES_DPTR, NEMOSTRN=NEMOSTRN_DPTR,  &
      & NPHIEPS=NPHIEPS_DPTR, NTAUOC=NTAUOC_DPTR, NSWH=NSWH_DPTR, NMWP=NMWP_DPTR, NEMOTAUX=NEMOTAUX_DPTR,  &
      & NEMOTAUY=NEMOTAUY_DPTR, NEMOWSWAVE=NEMOWSWAVE_DPTR, NEMOPHIF=NEMOPHIF_DPTR)
CALL INTFLDS_FIELD%DEVICE_POINTER(WSEMEAN=WSEMEAN_DPTR, WSFMEAN=WSFMEAN_DPTR, USTOKES=USTOKES_DPTR,  &
      & VSTOKES=VSTOKES_DPTR, STRNMS=STRNMS_DPTR, TAUXD=TAUXD_DPTR, TAUYD=TAUYD_DPTR, TAUOCXD=TAUOCXD_DPTR,  &
      & TAUOCYD=TAUOCYD_DPTR, TAUOC=TAUOC_DPTR, PHIOCD=PHIOCD_DPTR, PHIEPS=PHIEPS_DPTR, PHIAW=PHIAW_DPTR)
!$acc update device(  &
!$acc & dal2,k11w,lwnemotauoc,satweights,zpi4gm2,omxkm3_gc,rnum,zpifr,dfim_sim,chnkmin_u,z0rat,dthrn_u,ciblock,fklam1,ikp1,sinth,zpi4gm1,x0tauhf,af11,nfre,llgcbz0,bmaxokap,llunstr,flogsprdm1,gm1,lbiwbk,lwfluxout,xkm_gc,idamping,ikm1,cdicwa,dal1,indicessat,ndikcumul,fr,lwnemocousend,egrcrv,lwvflx_snl,sqrtgosurft,cithrsh,cdis,mlsthg,om3gmkm_gc,cofrm4,k1w,isnonlin,afcrv,licerun,cm_gc,fklap,zpi,llnormagam,lwnemocou,rnu,llcapchnk,cithrsh_tail,g,icode_cpl,nwav_gc,dfim,delth,k2w,dfimfr,dfimofr,swellft,ikp,rnlcoef,cdisvis,inlcoef,lwnemocoustrn,tauwshelter,k21w,nfre_odd,fklam,lwcou,lmaskice,lwamrsetci,nsdsnth,mfrstlw,z0tubmax,icode,wtauhf,kfrh,alpha,alphapmax,delkcc_gc_ns,alphamin,gamnconst,dfimfr2,nang,costh,bfcrv,delta_sdis,c2osqrtvg_gc,swellf5,xk_gc,xkmsqrtvgoc2_gc,delkcc_omxkm3_gc,nfre_red,lwflux,ang_gc_c,ikm,omega_gc,rhowg_dfim,cumulw,ang_gc_a,flmax,wspmin,th,tailfactor,idelt,ang_gc_b,dthrn_a,fklap1,lciwabr,iphys,betamaxoxkappa2,zalp,tailfactor_pm,lwnemocoustk,rn1_rn,fr5 &
!$acc &  )
IF (CDATE == CDTPRA) THEN
#ifdef WAM_PHYS_GPU
  CALL PROPAG_WAM_SCC(BLK2GLO, WVENVI_FIELD, WVPRPT_FIELD, FL1_DPTR)
!CALL SRC_CONTRIBS%ENSURE_HOST(LFL1=.TRUE.)
WRITE(IU06,*)"PROPAG_WAM_SCC TIME:",(TIME0+WAM_USER_CLOCK())*1.E-06
T1=-WAM_USER_CLOCK()
#else
  CALL PROPAG_WAM(BLK2GLO, WVENVI, WVPRPT, FL1)
#endif
  TIME1(1) = TIME1(1) + (TIME0+WAM_USER_CLOCK())*1.E-06
  CDATE = CDTPRO
  
ENDIF

!* RETRIEVING NEW FORCING FIELDS IF NEEDED.
!  ----------------------------------------
#ifdef WAM_PHYS_GPU
CALL NEWWIND_SCC(CDTIMP, CDATEWH, LLNEWFILE,      &
 &           CIWA_DPTR,CGROUP_DPTR, FF_NOW_FIELD, FF_NEXT_FIELD)
#else
CALL NEWWIND(CDTIMP, CDATEWH, LLNEWFILE,      &
 &           WVPRPT, FF_NOW, FF_NEXT)
#endif
! IT IS TIME TO INTEGRATE THE SOURCE TERMS
! ----------------------------------------
IF (CDATE >= CDTIMPNEXT) THEN
! COMPUTE UPDATE DUE TO SOURCE TERMS
  CALL GSTATS(1431,0)
  IF (LLSOURCE) THEN
    
    IF(FIRST_TIME)THEN
        CALL GPU_OFFLOAD_SIZE()
        FIRST_TIME=.FALSE.
    ENDIF
    TIME2=-WAM_USER_CLOCK()
#ifdef WAM_PHYS_GPU
CALL INIT_IMPLSCH(1,NPROMA_WAM,NCHNK,NANG, NFRE) ! only allocate local array first time

!$acc data present(FL1_DPTR,XLLWS_DPTR,MIJ_DPTR,WAVNUM_DPTR,CGROUP_DPTR,CIWA_DPTR,CINV_DPTR,XK2CG_DPTR,STOKFAC_DPTR,&
!$acc &            EMAXDPT_DPTR,INDEP_DPTR,DEPTH_DPTR,IOBND_DPTR,IODP_DPTR,CICOVER_DPTR,WSWAVE_DPTR,WDWAVE_DPTR,AIRD_DPTR,&
!$acc &            WSTAR_DPTR,UFRIC_DPTR,TAUW_DPTR,TAUWDIR_DPTR,Z0M_DPTR,Z0B_DPTR,CHRNCK_DPTR,CITHICK_DPTR,NEMOUSTOKES_DPTR,&
!$acc &            NEMOVSTOKES_DPTR,NEMOSTRN_DPTR,NPHIEPS_DPTR,NTAUOC_DPTR,NSWH_DPTR,NMWP_DPTR,NEMOTAUX_DPTR,NEMOTAUY_DPTR,&
!$acc &            NEMOWSWAVE_DPTR,NEMOPHIF_DPTR,WSEMEAN_DPTR,WSFMEAN_DPTR,USTOKES_DPTR,VSTOKES_DPTR,STRNMS_DPTR,TAUXD_DPTR,&
!$acc &            TAUYD_DPTR,TAUOCXD_DPTR,TAUOCYD_DPTR,TAUOC_DPTR,PHIOCD_DPTR,PHIEPS_DPTR,PHIAW_DPTR)
!$acc parallel loop gang vector_length(NPROMA_WAM)
      DO ICHNK=1,NCHNK
        CALL IMPLSCH_SCC(1, NPROMA_WAM, FL1_DPTR(:,:,:,ICHNK), WAVNUM_DPTR(:,:,ICHNK), &
        & CGROUP_DPTR(:,:,ICHNK), CIWA_DPTR(:,:,ICHNK), CINV_DPTR(:,:,ICHNK), &
        & XK2CG_DPTR(:,:,ICHNK), STOKFAC_DPTR(:,:,ICHNK), EMAXDPT_DPTR(:,ICHNK), &
        & INDEP_DPTR(:,ICHNK), DEPTH_DPTR(:,ICHNK), IOBND_DPTR(:,ICHNK), &
        & IODP_DPTR(:,ICHNK), AIRD_DPTR(:,ICHNK), WDWAVE_DPTR(:,ICHNK),  &
        & CICOVER_DPTR(:,ICHNK), WSWAVE_DPTR(:,ICHNK), WSTAR_DPTR(:,ICHNK), &
        & UFRIC_DPTR(:,ICHNK), TAUW_DPTR(:,ICHNK), TAUWDIR_DPTR(:,ICHNK), &
        & Z0M_DPTR(:,ICHNK), Z0B_DPTR(:,ICHNK), CHRNCK_DPTR(:,ICHNK), &
        & CITHICK_DPTR(:,ICHNK), NEMOUSTOKES_DPTR(:,ICHNK), NEMOVSTOKES_DPTR(:,ICHNK), &
        & NEMOSTRN_DPTR(:,ICHNK), NPHIEPS_DPTR(:,ICHNK), NTAUOC_DPTR(:,ICHNK), &
        & NSWH_DPTR(:,ICHNK), NMWP_DPTR(:,ICHNK), NEMOTAUX_DPTR(:,ICHNK), &
        & NEMOTAUY_DPTR(:,ICHNK), NEMOWSWAVE_DPTR(:,ICHNK), NEMOPHIF_DPTR(:,ICHNK), &
        & WSEMEAN_DPTR(:,ICHNK), WSFMEAN_DPTR(:,ICHNK), USTOKES_DPTR(:,ICHNK), &
        & VSTOKES_DPTR(:,ICHNK), STRNMS_DPTR(:,ICHNK), TAUXD_DPTR(:,ICHNK), &
        & TAUYD_DPTR(:,ICHNK), TAUOCXD_DPTR(:,ICHNK), TAUOCYD_DPTR(:,ICHNK), &
        & TAUOC_DPTR(:,ICHNK), PHIOCD_DPTR(:,ICHNK), PHIEPS_DPTR(:,ICHNK),  &
        & PHIAW_DPTR(:,ICHNK), MIJ_DPTR(:,ICHNK), XLLWS_DPTR(:,:,:,ICHNK),ICHNK)
      END DO
!$acc end parallel loop
!$acc end data

!CALL SRC_CONTRIBS%ENSURE_HOST(LFL1=.TRUE.)
!WRITE(IU06,*)"wamintgr output of IMPLSCH_SCC FL1(1:10) ",FL1(1:10, 1, 1,1),"max",maxval(FL1),minval(FL1)
#else
!$OMP  PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK) FIRSTPRIVATE(WVPRPT_FIELD, WVENVI_FIELD, FF_NOW_FIELD, WAM2NEMO_FIELD, &
!$OMP& INTFLDS_FIELD, SRC_CONTRIBS)
      DO ICHNK=1,NCHNK
        CALL WVPRPT_FIELD%UPDATE_VIEW(ICHNK)
        CALL WVENVI_FIELD%UPDATE_VIEW(ICHNK)
        CALL FF_NOW_FIELD%UPDATE_VIEW(ICHNK)
        CALL WAM2NEMO_FIELD%UPDATE_VIEW(ICHNK)
        CALL INTFLDS_FIELD%UPDATE_VIEW(ICHNK)
        CALL SRC_CONTRIBS%UPDATE_VIEW(ICHNK)
        
        CALL IMPLSCH(1, NPROMA_WAM, SRC_CONTRIBS%FL1, WVPRPT_FIELD%WAVNUM, WVPRPT_FIELD%CGROUP, WVPRPT_FIELD%CIWA,  &
        & WVPRPT_FIELD%CINV, WVPRPT_FIELD%XK2CG, WVPRPT_FIELD%STOKFAC, WVENVI_FIELD%EMAXDPT, WVENVI_FIELD%INDEP,  &
        & WVENVI_FIELD%DEPTH, WVENVI_FIELD%IOBND, WVENVI_FIELD%IODP, FF_NOW_FIELD%AIRD, FF_NOW_FIELD%WDWAVE,  &
        & FF_NOW_FIELD%CICOVER, FF_NOW_FIELD%WSWAVE, FF_NOW_FIELD%WSTAR, FF_NOW_FIELD%UFRIC, FF_NOW_FIELD%TAUW,  &
        & FF_NOW_FIELD%TAUWDIR, FF_NOW_FIELD%Z0M, FF_NOW_FIELD%Z0B, FF_NOW_FIELD%CHRNCK, FF_NOW_FIELD%CITHICK,  &
        & WAM2NEMO_FIELD%NEMOUSTOKES, WAM2NEMO_FIELD%NEMOVSTOKES, WAM2NEMO_FIELD%NEMOSTRN, WAM2NEMO_FIELD%NPHIEPS,  &
        & WAM2NEMO_FIELD%NTAUOC, WAM2NEMO_FIELD%NSWH, WAM2NEMO_FIELD%NMWP, WAM2NEMO_FIELD%NEMOTAUX, WAM2NEMO_FIELD%NEMOTAUY,  &
        & WAM2NEMO_FIELD%NEMOWSWAVE, WAM2NEMO_FIELD%NEMOPHIF, INTFLDS_FIELD%WSEMEAN, INTFLDS_FIELD%WSFMEAN,  &
        & INTFLDS_FIELD%USTOKES, INTFLDS_FIELD%VSTOKES, INTFLDS_FIELD%STRNMS, INTFLDS_FIELD%TAUXD, INTFLDS_FIELD%TAUYD,  &
        & INTFLDS_FIELD%TAUOCXD, INTFLDS_FIELD%TAUOCYD, INTFLDS_FIELD%TAUOC, INTFLDS_FIELD%PHIOCD, INTFLDS_FIELD%PHIEPS,  &
        & INTFLDS_FIELD%PHIAW, SRC_CONTRIBS%MIJ, SRC_CONTRIBS%XLLWS)
      END DO
!$OMP END PARALLEL DO
#endif
    TIME1(3) = TIME1(3) + (TIME2+WAM_USER_CLOCK())*1.E-06

    IF (LWNEMOCOU) NEMONTAU = NEMONTAU + 1

  ELSE
!   NO SOURCE TERM CONTRIBUTION
!$OMP      PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK)  
    DO ICHNK = 1, NCHNK
      MIJ(:,ICHNK) = NFRE
      FL1(:,:,:,ICHNK) = MAX(FL1(:,:,:,ICHNK), EPSMIN)
      XLLWS(:,:,:,ICHNK) = 0.0_JWRB
    ENDDO
!$OMP      END PARALLEL DO
  ENDIF
  CALL GSTATS(1431,1)


!*       UPDATE FORCING FIELDS TIME COUNTER 
!        ----------------------------------
  IF (LLNEWFILE) THEN
    LLNEWFILE = .FALSE.
    IDELWH = MAX(IDELWI, IDELPRO)
    CALL INCDATE(CDAWIFL, IDELWH)
    CALL INCDATE(CDATEFL, IDELWH)
  ENDIF

  CDATEWO=CDATEWH
  CDTIMP=CDTIMPNEXT
  CALL INCDATE(CDTIMPNEXT, IDELT)   

ENDIF

    WRITE(IU06,*)"WAMINTGR implsch time",(T1+WAM_USER_CLOCK())*1.E-06,"s"
    WRITE(IU06,*)"WAMINTGR time ",(TIME0+WAM_USER_CLOCK())*1.E-06,"s"
IF (LHOOK) CALL DR_HOOK('WAMINTGR',1,ZHOOK_HANDLE)

  CONTAINS
    SUBROUTINE GPU_OFFLOAD_SIZE()
      INTEGER, PARAMETER :: LargeInt = SELECTED_INT_KIND(18)
      INTEGER(KIND=LargeInt) :: size0, size1
  
      size0 = sizeof(WVPRPT%WAVNUM)
      size0 = size0 + sizeof(WVPRPT%CGROUP)
      size0 = size0 + sizeof(WVPRPT%CIWA)
      size0 = size0 + sizeof(WVPRPT%CINV)
      size0 = size0 + sizeof(WVPRPT%XK2CG)
      size0 = size0 + sizeof(WVPRPT%STOKFAC)
 
      size0 = size0 + sizeof(WVENVI%EMAXDPT)
      size0 = size0 + sizeof(WVENVI%INDEP)
      size0 = size0 + sizeof(WVENVI%DEPTH)
      size0 = size0 + sizeof(WVENVI%IOBND)
      size0 = size0 + sizeof(WVENVI%IODP)
 
      size0 = size0 + sizeof(FF_NOW%AIRD)
      size0 = size0 + sizeof(FF_NOW%WDWAVE)
      size0 = size0 + sizeof(FF_NOW%CICOVER)
      size0 = size0 + sizeof(FF_NOW%WSWAVE)
      size0 = size0 + sizeof(FF_NOW%WSTAR)
      size0 = size0 + sizeof(FF_NOW%UFRIC)
      size0 = size0 + sizeof(FF_NOW%TAUW)
      size0 = size0 + sizeof(FF_NOW%TAUWDIR)
      size0 = size0 + sizeof(FF_NOW%Z0M)
      size0 = size0 + sizeof(FF_NOW%Z0B)
      size0 = size0 + sizeof(FF_NOW%CHRNCK)
      size0 = size0 + sizeof(FF_NOW%CITHICK)
 
      size0 = size0 + sizeof(WAM2NEMO%NEMOUSTOKES)
      size0 = size0 + sizeof(WAM2NEMO%NEMOVSTOKES)
      size0 = size0 + sizeof(WAM2NEMO%NEMOSTRN)
      size0 = size0 + sizeof(WAM2NEMO%NPHIEPS)
      size0 = size0 + sizeof(WAM2NEMO%NTAUOC)
      size0 = size0 + sizeof(WAM2NEMO%NSWH)
      size0 = size0 + sizeof(WAM2NEMO%NMWP)
      size0 = size0 + sizeof(WAM2NEMO%NEMOTAUX)
      size0 = size0 + sizeof(WAM2NEMO%NEMOTAUY)
      size0 = size0 + sizeof(WAM2NEMO%NEMOWSWAVE)
      size0 = size0 + sizeof(WAM2NEMO%NEMOPHIF)
 
      size0 = size0 + sizeof(INTFLDS%WSEMEAN)
      size0 = size0 + sizeof(INTFLDS%WSFMEAN)
      size0 = size0 + sizeof(INTFLDS%USTOKES)
      size0 = size0 + sizeof(INTFLDS%VSTOKES)
      size0 = size0 + sizeof(INTFLDS%STRNMS)
      size0 = size0 + sizeof(INTFLDS%TAUXD)
      size0 = size0 + sizeof(INTFLDS%TAUYD)
      size0 = size0 + sizeof(INTFLDS%TAUOCXD)
      size0 = size0 + sizeof(INTFLDS%TAUOCYD)
      size0 = size0 + sizeof(INTFLDS%TAUOC)
      size0 = size0 + sizeof(INTFLDS%PHIOCD)
      size0 = size0 + sizeof(INTFLDS%PHIEPS)
      size0 = size0 + sizeof(INTFLDS%PHIAW)
 
      size0 = size0 + sizeof(FL1)
      size0 = size0 + sizeof(XLLWS)
      size0 = size0 + sizeof(MIJ)
 
      size1 = sizeof(dal2) +sizeof(k11w) +sizeof(lwnemotauoc) +sizeof(satweights) +sizeof(zpi4gm2) +&
 &            sizeof(omxkm3_gc) +sizeof(rnum) +sizeof(zpifr) +sizeof(dfim_sim)
      size1 = size1 + sizeof(chnkmin_u) +sizeof(z0rat) +sizeof(dthrn_u) +sizeof(ciblock) +sizeof(fklam1) +&
 &            sizeof(ikp1) +sizeof(sinth) +sizeof(zpi4gm1) +sizeof(x0tauhf)
      size1 = size1 + sizeof(af11) +sizeof(nfre) +sizeof(llgcbz0) +sizeof(bmaxokap) +sizeof(llunstr) +&
 &            sizeof(flogsprdm1) +sizeof(gm1) +sizeof(lbiwbk) +sizeof(lwfluxout)
      size1 = size1 + sizeof(xkm_gc) +sizeof(idamping) +sizeof(ikm1) +sizeof(cdicwa) +sizeof(dal1) +&
 &            sizeof(indicessat) +sizeof(ndikcumul) +sizeof(fr) +sizeof(lwnemocousend)
      size1 = size1 + sizeof(egrcrv) +sizeof(lwvflx_snl) +sizeof(sqrtgosurft) +sizeof(cithrsh) +sizeof(cdis) +&
 &            sizeof(mlsthg) +sizeof(om3gmkm_gc) +sizeof(cofrm4) +sizeof(k1w)
      size1 = size1 + sizeof(isnonlin) +sizeof(afcrv) +sizeof(licerun) +sizeof(cm_gc) +sizeof(fklap) +sizeof(zpi) +&
 &            sizeof(llnormagam) +sizeof(lwnemocou) +sizeof(rnu) +sizeof(llcapchnk)
      size1 = size1 + sizeof(cithrsh_tail) +sizeof(g) +sizeof(icode_cpl) +sizeof(nwav_gc) +sizeof(dfim) +sizeof(delth) +&
 &            sizeof(k2w) +sizeof(dfimfr) +sizeof(dfimofr) +sizeof(swellft)
      size1 = size1 + sizeof(ikp) +sizeof(rnlcoef) +sizeof(cdisvis) +sizeof(inlcoef) +sizeof(lwnemocoustrn) +&
 &            sizeof(tauwshelter) +sizeof(k21w) +sizeof(nfre_odd) +sizeof(fklam)
      size1 = size1 + sizeof(lwcou) +sizeof(lmaskice) +sizeof(lwamrsetci) +sizeof(nsdsnth) +sizeof(mfrstlw) +&
 &            sizeof(z0tubmax) +sizeof(icode) +sizeof(wtauhf) +sizeof(kfrh) +sizeof(alpha)
      size1 = size1 + sizeof(alphapmax) +sizeof(delkcc_gc_ns) +sizeof(alphamin) +sizeof(gamnconst) +sizeof(dfimfr2) +&
 &            sizeof(nang) +sizeof(costh) +sizeof(bfcrv) +sizeof(delta_sdis)
      size1 = size1 + sizeof(c2osqrtvg_gc) +sizeof(swellf5) +sizeof(xk_gc) +sizeof(xkmsqrtvgoc2_gc) +&
 &            sizeof(delkcc_omxkm3_gc) +sizeof(nfre_red) +sizeof(lwflux) +sizeof(ang_gc_c)
      size1 = size1 + sizeof(ikm) +sizeof(omega_gc) +sizeof(rhowg_dfim) +sizeof(cumulw) +sizeof(ang_gc_a) +&
 &            sizeof(flmax) +sizeof(wspmin) +sizeof(th) +sizeof(tailfactor) +sizeof(idelt)
      size1 = size1 + sizeof(ang_gc_b) +sizeof(dthrn_a) +sizeof(fklap1) +sizeof(lciwabr) +sizeof(iphys) +&
 &            sizeof(betamaxoxkappa2) +sizeof(zalp) +sizeof(tailfactor_pm)
      size1 = size1 + sizeof(lwnemocoustk) +sizeof(rn1_rn) +sizeof(fr5)

      write(IU06,*)"GPU offload size arguments:", dble(size0/1024./1024./1024.)," GB"
      write(IU06,*) "GPU offload size global vars:",dble(size1/1024./1024./1024.)," GB"
      write(IU06,*) "GPU offload size total:",dble((size0+size1)/1024./1024./1024.)," GB"
    END SUBROUTINE GPU_OFFLOAD_SIZE
END SUBROUTINE WAMINTGR
