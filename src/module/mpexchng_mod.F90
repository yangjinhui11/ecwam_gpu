module MPEXCHNG_MOD

contains

subroutine MPEXCHNG_SCC(FLD,NDIM2,ND3S,ND3E)

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,    &
     &          KTAG
      USE YOWSPEC, ONLY   : NTOPE    ,NTOPEMAX ,IJTOPE   ,NGBTOPE  ,    &
     &          NTOPELST   ,NFROMPE  ,NFROMPEMAX,NIJSTART,NGBFROMPE,    &
     &          NFROMPELST

      USE MPL_MODULE, ONLY : MPL_RECV, MPL_SEND, MPL_WAIT=>MPL_WAITS, &
                           & JP_NON_BLOCKING_STANDARD
      USE MPL_MODULE, ONLY : MPL_COMM
      USE YOWTEST, ONLY: IU06
      USE MPI
!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: NDIM2, ND3S, ND3E
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NDIM2, ND3S:ND3E), INTENT(INOUT):: FLD

      INTEGER(KIND=JWIM) :: INGB, M, K, IH, IJ, KCOUNT, IPROC, IR
      INTEGER(KIND=JWIM) :: NBUFMAX
      INTEGER(KIND=JWIM) :: NDIM3,IERROR
      INTEGER(KIND=JWIM), DIMENSION(NGBTOPE+NGBFROMPE) :: IREQ

      REAL(KIND=JWRB) :: ZDUM(2)
      REAL(KIND=JWRB), ALLOCATABLE :: ZCOMBUFS(:,:)
      REAL(KIND=JWRB), ALLOCATABLE :: ZCOMBUFR(:,:)

      LOGICAL :: LLOK
      real::start,finish
!----------------------------------------------------------------------


      IF (NPROC <= 1) THEN
        RETURN
      ENDIF

      CALL GSTATS_BARRIER(736)

      NDIM3 = ND3E-ND3S+1

      NBUFMAX=MAX(NTOPEMAX,NFROMPEMAX)*NDIM2*NDIM3
      ALLOCATE(ZCOMBUFS(NBUFMAX,NGBTOPE))
      ALLOCATE(ZCOMBUFR(NBUFMAX,NGBFROMPE))
!$acc enter data create(ZCOMBUFS,ZCOMBUFR)
!     PACK SEND BUFFERS FOR NGBTOPE NEIGHBOURING PE's
!     -------------------------------------------------
!$acc parallel present(FLD,NTOPELST,IJTOPE,NTOPE,ZCOMBUFS)
!$acc loop gang 
      DO INGB=1,NGBTOPE
        IPROC=NTOPELST(INGB)
!$acc loop vector collapse(3) private(kcount) independent
        DO M = ND3S, ND3E
          DO K = 1, NDIM2
            DO IH = 1, NTOPE(IPROC)
              IJ=IJTOPE(IH,IPROC)
              KCOUNT= (M-ND3S)*NDIM2*NTOPE(IPROC)+(K-1)*NTOPE(IPROC)+IH
              ZCOMBUFS(KCOUNT,INGB)=FLD(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$acc end parallel

      IR=0
      CALL GSTATS(676,0)

      DO INGB=1,NGBFROMPE
        IR=IR+1
        IPROC=NFROMPELST(INGB)
        KCOUNT=NDIM3*NDIM2*NFROMPE(IPROC)
        CALL MPL_RECV(ZCOMBUFR(1:KCOUNT,INGB),KSOURCE=IPROC,KTAG=KTAG,  &
     &     KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ(IR),         &
     &     CDSTRING='MPEXCHNG:')
!!$acc host_data use_device(ZCOMBUFR(1:KCOUNT,INGB))
!    call MPI_IRECV(ZCOMBUFR(1:KCOUNT,INGB),KCOUNT,MPI_REAL8,IPROC-1,KTAG,MPL_COMM,IREQ(IR),IERROR)
!!$acc end host_data
      ENDDO

      DO INGB=1,NGBTOPE
        IR=IR+1
        IPROC=NTOPELST(INGB)
        KCOUNT=NDIM3*NDIM2*NTOPE(IPROC)
!$acc update host(ZCOMBUFS(1:KCOUNT,INGB))
        CALL MPL_SEND(ZCOMBUFS(1:KCOUNT,INGB),KDEST=IPROC,KTAG=KTAG,    &
     &     KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ(IR),         &
     &     CDSTRING='MPEXCHNG:')
!!$acc host_data use_device(ZCOMBUFS(1:KCOUNT,INGB))
!    CALL MPI_ISEND(ZCOMBUFS(1:KCOUNT,INGB),KCOUNT,MPI_REAL8,IPROC-1,KTAG,MPL_COMM,IREQ(IR),IERROR)
!!$acc end host_data
      ENDDO

!     NOW WAIT FOR ALL TO COMPLETE

      CALL MPL_WAIT(KREQUEST=IREQ(1:IR),CDSTRING='MPEXCHNG:')
      
      DO INGB=1,NGBFROMPE
        IPROC=NFROMPELST(INGB)
        KCOUNT=NDIM3*NDIM2*NFROMPE(IPROC)
        !$acc update device(ZCOMBUFR(1:KCOUNT,INGB))
      ENDDO


      CALL GSTATS(676,1)

!     DECODE THE RECEIVED BUFFERS

      CALL GSTATS(1893,0)
!$acc kernels present(FLD,NFROMPELST,NFROMPE,NIJSTART,ZCOMBUFR)
!$acc loop gang
      DO INGB=1,NGBFROMPE
        IPROC=NFROMPELST(INGB)
!$acc loop vector collapse(3) private(KCOUNT) independent
        DO M = ND3S, ND3E
          DO K = 1, NDIM2
            DO IH = 1, NFROMPE(IPROC)
              IJ=NIJSTART(IPROC)+IH-1
              KCOUNT= (M-ND3S)*NDIM2*NFROMPE(IPROC)+(K-1)*NFROMPE(IPROC)+IH
              FLD(IJ,K,M)=ZCOMBUFR(KCOUNT,INGB)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$acc end kernels
      CALL GSTATS(1893,1)

      KTAG=KTAG+1

!$acc exit data delete(ZCOMBUFS,ZCOMBUFR)
      DEALLOCATE(ZCOMBUFS)
      DEALLOCATE(ZCOMBUFR)


      END SUBROUTINE MPEXCHNG_SCC


endmodule
