SUBROUTINE GSTATS_PRINT(KULOUT,PAVEAVE,KLEN)

!**** *GSTATS_PRINT* - print timing statistics

!     PURPOSE.
!     --------
!       To print out timings gathered by GSTATS


!**   INTERFACE.
!     ----------
!       *CALL* *GSTATS_PRINT*

!        EXPLICIT ARGUMENTS     None
!        --------------------


!        IMPLICIT ARGUMENTS
!        --------------------
!        Module YOMSTATS

!     METHOD.
!     -------


!     EXTERNALS.   
!     ----------   

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!        Mats Hamrud ECMWF

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 98-11-15
!        D.Salmond: 99-09-21 : Timer for SLCOMM2
!        G.Mozdzynski 05-09-25 : fix master ncalls overwrite for nproc>1
!        C.Larsson    8-May-2006 : Added xml file output
!        G.Mozdzynski 16-Oct-2007 : xml file output under switch LXML_STATS
!        P.Towers     11-May-2011 : mpl comms statistics output
!     ------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE YOMGSTATS
USE YOMMPI   , ONLY : MREALT
USE MPL_MODULE
USE MPL_STATS_MOD

IMPLICIT NONE

INTEGER(KIND=JPIM) :: KULOUT,KLEN
REAL(KIND=JPRB) :: PAVEAVE(0:KLEN)
CHARACTER*7  CLACTION(0:3)
CHARACTER(LEN=JPMAXDELAYS*10) CLTEMP
INTEGER(KIND=JPIM),PARAMETER :: JPARRAYS=8
REAL(KIND=JPRB) :: ZREABUF(JPARRAYS*(JPMAXSTAT+1))
REAL(KIND=JPRB) :: ZAVEAVE(0:JPMAXSTAT),ZAVEMAX(0:JPMAXSTAT),ZTIMELCALL(0:JPMAXSTAT),&
         &ZTHISTIME(0:JPMAXSTAT),ZFRACMAX(0:JPMAXSTAT),&
         &ZSUMMAX(0:JPMAXSTAT),ZSUMTOT(0:JPMAXSTAT)
REAL(KIND=JPRB) :: ZT_SUM,ZT_SUM2,ZT_SUM3,ZT_SUMIO,ZT_SUM4,ZT_SUM5,ZT_SUMB
REAL(KIND=JPRB) :: ZDELAY,ZDELAY_MAX

REAL(KIND=JPRB) :: ZMPL(NPROC_STATS)
REAL(KIND=JPRB) :: ZBAR(NPROC_STATS)
REAL(KIND=JPRB) :: ZGBR(NPROC_STATS)
REAL(KIND=JPRB) :: ZGB2(NPROC_STATS)
REAL(KIND=JPRB) :: ZOMP(NPROC_STATS)
REAL(KIND=JPRB) :: ZIO (NPROC_STATS)
REAL(KIND=JPRB) :: ZSER(NPROC_STATS)
REAL(KIND=JPRB) :: ZMXD(NPROC_STATS)

INTEGER(KIND=JPIM) :: ICALLSX(0:JPMAXSTAT)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: ICALLS, ILBUF, ILSEND, ILRECV, &
             &ISEND, ITAG, JJ, JNUM, JROC, JCALL, ICALLER,IACTION
INTEGER(KIND=JPIM) :: IMEM, INUM, JMEM
INTEGER(KIND=JPIM) :: JDELAY, IDELAY
INTEGER(KIND=JPIM) :: NSEND,NRECV

!     LOCAL REAL SCALARS
REAL(KIND=JPRB) :: ZAVE, ZAVETCPU, ZAVEVCPU, ZCOMTIM, ZDETAIL,&
          &ZFRAC, ZMAX, ZMEAN, ZSTDDEV, ZSUM, ZSUMB, &
          &ZTOTAL, ZTOTCPU, ZTOTUNBAL, ZTOTVCPU, &
          &ZUNBAL, ZMEANT, ZMAXT

REAL(KIND=JPRB)    :: SBYTES,RBYTES,SENDRATE,RECVRATE
REAL(KIND=JPRB)    :: AVGSENDLEN,AVGRECVLEN
REAL(KIND=JPRB)    :: MAXCOMMTIME(501:1000)
REAL(KIND=JPRB)    :: TOTSENDBYTES(501:1000)
REAL(KIND=JPRB)    :: TOTRECVBYTES(501:1000)

INTEGER(KIND=JPIM) :: IXMLLUN  

!     ------------------------------------------------------------------

! write(0,*) "JPMAXSTAT,NPRNT_STATS =",JPMAXSTAT,NPRNT_STATS

ILBUF = JPARRAYS*(JPMAXSTAT+1)
ILRECV = 500*4
ZAVEAVE(:) = 0.0_JPRB
ZAVEMAX(:) = 0.0_JPRB
ZFRACMAX(:)= 0.0_JPRB
ZSUMMAX(:)= 0.0_JPRB
ZSUMTOT(:)= 0.0_JPRB

! OPEN GSTATS.XML for xml statistics
IF(LXML_STATS)THEN
  IXMLLUN=40
  OPEN (UNIT=IXMLLUN, FILE='gstats.xml',ACTION='write')
  WRITE(IXMLLUN,'(A)')'<?xml version="1.0" encoding="UTF-8"?>'
  WRITE(IXMLLUN,'(A)')'<gstats '
  WRITE(IXMLLUN,'(A)')' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"&
  & xsi:schemaLocation="http://www.ecmwf.int/services/prepifs/gstats   ./gstats.xsd">'
ENDIF

WRITE(KULOUT,'(A)')'===-=== START OF TIMING STATISTICS ===-==='
IF(LSYNCSTATS.AND.NPROC_STATS>1) THEN
  WRITE(KULOUT,'(A)')'START OF TIMINGS SYNCRONIZED'
ENDIF
IF(LDETAILED_STATS) THEN
  LSTATS_COMMS=.TRUE.
  LSTATS_OMP=.TRUE.
ENDIF
IF(LBARRIER_STATS2) THEN
  DO JNUM=1,JBMAXBASE
    IF(NBAR_PTR(JNUM) > 0) THEN
      INUM=NBAR_PTR(JNUM)
      CCDESC(INUM)=CCDESC(JNUM)
      CCTYPE(INUM)='GB2'
    ENDIF
  ENDDO
ENDIF
IF (LSTATS .AND. MYPROC_STATS /= 1) THEN
  JJ = 1
  DO JNUM=0,JPMAXSTAT
    ZREABUF(JJ  ) = TIMESUM(JNUM)
    ZREABUF(JJ+1) = TIMESQSUM(JNUM)
    ZREABUF(JJ+2) = TIMEMAX(JNUM)
    ZREABUF(JJ+3) = TIMESUMB(JNUM)
    ZREABUF(JJ+4) = TIMELCALL(JNUM)
    ZREABUF(JJ+5) = NCALLS(JNUM)
    ZREABUF(JJ+6) = TTCPUSUM(JNUM)
    ZREABUF(JJ+7) = TVCPUSUM(JNUM)
    JJ = JJ+JPARRAYS
  ENDDO
  ILSEND = ILBUF
  ISEND =1
  ITAG = JPTAGSTAT
  CALL MPL_SEND(ZREABUF(1:ILSEND),KDEST=NPRCIDS_STATS(ISEND), &
       & KTAG=ITAG,CDSTRING='GSTATS_PRINT:')

  IF(LSTATS_MPL) THEN
    JJ=1
    DO JNUM=501,1000
      ZREABUF(JJ  ) = NUMSEND(JNUM)
      ZREABUF(JJ+1) = SENDBYTES(JNUM)
      ZREABUF(JJ+2) = NUMRECV(JNUM)
      ZREABUF(JJ+3) = RECVBYTES(JNUM)
      JJ=JJ+4
    ENDDO

    ILSEND = JJ-1
    ITAG = JPTAGSTAT + 1

    CALL MPL_SEND(ZREABUF(1:ILSEND),KDEST=NPRCIDS_STATS(ISEND), &
         & KTAG=ITAG,CDSTRING='GSTATS_PRINT:')

  ENDIF
ELSEIF(LSTATS) THEN
  IF(LSTATS_MPL) THEN
    DO JNUM=501,1000
      MAXCOMMTIME(JNUM)=0.0_JPRB
      TOTSENDBYTES(JNUM)=0.0_JPRB
      TOTRECVBYTES(JNUM)=0.0_JPRB
    ENDDO
  ENDIF
  DO JROC=1,NPROC_STATS
    IF (JROC /= 1) THEN
      ITAG = JPTAGSTAT
      CALL MPL_RECV(ZREABUF(1:ILBUF),KSOURCE=NPRCIDS_STATS(JROC), &
       & KTAG=ITAG,CDSTRING='GSTATS_PRINT:')
      JJ = 1
      DO JNUM=0,JPMAXSTAT
        TIMESUM(JNUM)   = ZREABUF(JJ  )
        TIMESQSUM(JNUM) = ZREABUF(JJ+1)
        TIMEMAX(JNUM)   = ZREABUF(JJ+2)
        TIMESUMB(JNUM)  = ZREABUF(JJ+3)
        TIMELCALL(JNUM) = ZREABUF(JJ+4)
        ICALLSX(JNUM)    = NINT(ZREABUF(JJ+5))
        TTCPUSUM(JNUM)  = ZREABUF(JJ+6)
        TVCPUSUM(JNUM)  = ZREABUF(JJ+7)
        JJ = JJ+JPARRAYS
      ENDDO

      IF(LSTATS_MPL) THEN
        ITAG = JPTAGSTAT+1
        CALL MPL_RECV(ZREABUF(1:ILRECV),KSOURCE=NPRCIDS_STATS(JROC), &
         & KTAG=ITAG,CDSTRING='GSTATS_PRINT:')
        JJ = 1
        DO JNUM=501,1000
          NUMSEND(JNUM)   = NINT(ZREABUF(JJ))
          SENDBYTES(JNUM) = ZREABUF(JJ+1)
          NUMRECV(JNUM)   = NINT(ZREABUF(JJ+2))
          RECVBYTES(JNUM) = ZREABUF(JJ+3)
          JJ=JJ+4
        ENDDO
      ENDIF

    ELSE
      ICALLSX(:)=NCALLS(:)

    ENDIF
    IF (JROC == 1) THEN
      ZTOTAL=TIMESUM(0)
      ZTOTCPU = TTCPUSUM(0)
      ZTOTVCPU = TVCPUSUM(0)
    ENDIF
    IF(.NOT. LSTATSCPU) THEN
      TTCPUSUM(1:JPMAXSTAT) = -0.0_JPRB
      TVCPUSUM(1:JPMAXSTAT) = -0.0_JPRB
    ENDIF
    ZT_SUM=0.0_JPRB
    ZT_SUM2=0.0_JPRB
    ZT_SUM3=0.0_JPRB
    ZT_SUM4=0.0_JPRB
    ZT_SUM5=0.0_JPRB
    ZT_SUMIO=0.0_JPRB
    ZT_SUMB=0.0_JPRB
    IF( LDETAILED_STATS .AND. JROC <= NPRNT_STATS ) THEN
      WRITE(KULOUT,'(A,I4)') 'TIMING STATISTICS:PROCESSOR=',JROC
      IF(LXML_STATS)THEN
        WRITE(IXMLLUN,'(A,I4,A)')'<timing processor="',JROC,'">'
      ENDIF
      IF(NPROC_STATS > 1) THEN
        WRITE(KULOUT,'(A,F6.1,A)')'STARTUP COST ',&
         &TIME_START(JROC),' SECONDS'
      ENDIF
      WRITE(KULOUT,'(A)')&
       &' NUM     ROUTINE                                     '//&
       &'CALLS  SUM(s)   AVE(ms)   CPUAVE(ms) VAVE(ms) '//&
       &'STDDEV(ms)  MAX(ms) '//&
       &'SUMB(s) FRAC(%)'
    ENDIF
    DO JNUM=0,JPMAXSTAT
      IF(ICALLSX(JNUM) > 1) THEN
        ICALLS = ICALLSX(JNUM)/2
        ZSUM = TIMESUM(JNUM)
        ZAVE = TIMESUM(JNUM)/ICALLS*1000._JPRB
        ZMAX = TIMEMAX(JNUM)*1000._JPRB
        ZSUMB = TIMESUMB(JNUM)
        ZFRAC = TIMESUM(JNUM)/ZTOTAL*100.0_JPRB
        ZFRACMAX(JNUM)=MAX(ZFRACMAX(JNUM),ZFRAC)
        ZSUMMAX(JNUM)=MAX(ZSUMMAX(JNUM),TIMESUM(JNUM))
        ZSUMTOT(JNUM)=ZSUMTOT(JNUM)+ZSUM
        ZAVEAVE(JNUM)=ZAVEAVE(JNUM)+ZAVE
        ZAVEMAX(JNUM)=MAX(ZAVEMAX(JNUM),ZAVE)
        ZAVETCPU = TTCPUSUM(JNUM)/ICALLS*1000._JPRB
        ZAVEVCPU = TVCPUSUM(JNUM)/ICALLS*1000._JPRB
        IF(ICALLS > 1 ) THEN
          ZSTDDEV = 1000._JPRB*&
           &SQRT(MAX((TIMESQSUM(JNUM)-TIMESUM(JNUM)**2/ICALLS)&
           &/(ICALLS-1),0.0_JPRB))
        ELSE
          ZSTDDEV = 0.0_JPRB
        ENDIF
        IF(CCTYPE(JNUM).EQ.'MPL') THEN
          ZT_SUM=ZT_SUM+ZSUM
          ZT_SUMB=ZT_SUMB+ZSUMB
        ENDIF
        IF(CCTYPE(JNUM).EQ.'BAR' .OR. CCTYPE(JNUM).EQ.'GBR' .OR. CCTYPE(JNUM).EQ.'GB2' ) THEN
          ZT_SUM4=ZT_SUM4+ZSUM
          ZT_SUMB=ZT_SUMB+ZSUMB
        ENDIF
        IF(CCTYPE(JNUM).EQ.'OMP') THEN
          ZT_SUM2=ZT_SUM2+ZSUM
          ZT_SUMB=ZT_SUMB+ZSUMB
        ENDIF
        IF(CCTYPE(JNUM).EQ.'IO-') THEN
          ZT_SUMIO=ZT_SUMIO+ZSUM
          ZT_SUMB=ZT_SUMB+ZSUMB
        ENDIF
        IF(CCTYPE(JNUM).EQ.'SER') THEN
          ZT_SUM3=ZT_SUM3+ZSUM
          ZT_SUMB=ZT_SUMB+ZSUMB
        ENDIF
        IF(CCTYPE(JNUM).EQ.'MXD') THEN
          ZT_SUM5=ZT_SUM5+ZSUM
          ZT_SUMB=ZT_SUMB+ZSUMB
        ENDIF
        IF( LDETAILED_STATS .AND. JROC <= NPRNT_STATS ) THEN
!         IF(JNUM < 501 .OR. LSTATS_COMMS .OR. LSTATS_OMP) THEN 
            WRITE(KULOUT,'(I4,1X,A3,1X,A40,1X,I5,6(1X,F9.1),1X,F5.1,1X,F8.2)')&
             &JNUM,CCTYPE(JNUM),CCDESC(JNUM),ICALLS,ZSUM,ZAVE,ZAVETCPU,ZAVEVCPU,&
             &ZSTDDEV,ZMAX,ZSUMB,ZFRAC
            IF(LXML_STATS)THEN
              WRITE(IXMLLUN,&
               & '(A,I4,A,//,A,A40,A,//,A,I5,A,//,6(A,F9.1,A,//),A,F5.1,A,//,A,F8.2,A,//,A)')&
               & '<num id="',JNUM,'">',&
               & '<description>',CCDESC(JNUM),'</description>',&
               & '<sum unit="seconds">',ICALLS,'</sum>',&
               & '<sumcall unit="ms">',ZSUM,'</sumcall>',&
               & '<average unit="ms">',ZAVE,'</average>',&
               & '<cpuaverage unit="ms">',ZAVETCPU,'</cpuaverage>',&
               & '<vave unit="ms">',ZAVEVCPU,'</vave>',&
               & '<stddev unit="ms">',ZSTDDEV,'</stddev>',&
               & '<max unit="ms">',ZMAX,'</max>',&
               & '<sumb unit="seconds">',ZSUMB,'</sumb>',&
               & '<frac unit="percent">',ZFRAC,'</frac>',&
               & '</num>'
             ENDIF
!        ENDIF
        ENDIF
      ENDIF
    ENDDO
!    ZCOMTIM = SUM(TIMESUM(51:99))+SUM(TIMESUM(151:199))
!    ZDETAIL = SUM(TIMESUM(:))-TIMESUM(0)-TIMESUM(1)-SUM(TIMESUM(20:23))
    IF( LDETAILED_STATS .AND. JROC <= NPRNT_STATS ) THEN
      WRITE(KULOUT,*) ''
      WRITE(KULOUT,'((A,2F8.1))')&
       &'CPU-TIME AND VECTOR CPU-TIME AS PERCENT OF TOTAL ',&
       &TTCPUSUM(0)/TIMESUM(0)*100.0_JPRB,TVCPUSUM(0)/TIMESUM(0)*100.0_JPRB
      IF(LXML_STATS)THEN
        WRITE(IXMLLUN,'((A,F8.1,A,//,A,F8.1,A))')&
         &'<cpufraction>',&
         &TTCPUSUM(0)/TIMESUM(0)*100.0_JPRB,&
         &'</cpufraction>',&
         &'<cpuvectorfraction>',&
         &TVCPUSUM(0)/TIMESUM(0)*100.0_JPRB,&
         &'</cpuvectorfraction>'
      ENDIF


      IF(ZT_SUM > 0.0_JPRB) THEN
        WRITE(KULOUT,'(A,F10.1,A,F6.2,A)')'SUMMED TIME IN COMMUNICATIONS   = '&
         & ,ZT_SUM, ' SECONDS ',ZT_SUM/TIMESUM(0)*100.0_JPRB,&
         &' PERCENT OF TOTAL'
        IF(LXML_STATS)THEN
          WRITE(IXMLLUN,'(A,F10.1,A,//,A,F6.2,A)')'<zcom unit="seconds">',&
           &ZT_SUM,'</zcom>',&
           &'<fraczcom unit="percent">',ZT_SUM/TIMESUM(0)*100.0_JPRB,&
           &'</fraczcom>'
        ENDIF
      ENDIF
      IF(ZT_SUM2 > 0.0_JPRB) THEN
        WRITE(KULOUT,'(A,F10.1,A,F6.2,A)')'SUMMED TIME IN PARALLEL REGIONS = '&
         & ,ZT_SUM2, ' SECONDS ',ZT_SUM2/TIMESUM(0)*100.0_JPRB,&
         &' PERCENT OF TOTAL'
        IF(LXML_STATS)THEN
          WRITE(IXMLLUN,'(A,F10.1,A,//,A,F6.2,A)') &
           &'<parallelztime unit="seconds">',&
           &ZT_SUM2, '</parallelztime>',&
           &'<fracparallelztime unit="percent">',ZT_SUM2/TIMESUM(0)*100.0_JPRB,&
           &'</fracparallelztime>'
        ENDIF
      ENDIF
      IF(ZT_SUMIO > 0.0_JPRB) THEN
        WRITE(KULOUT,'(A,F10.1,A,F6.2,A)')'SUMMED TIME IN I/O SECTIONS     = '&
         & ,ZT_SUMIO, ' SECONDS ',ZT_SUMIO/TIMESUM(0)*100.0_JPRB,&
         &' PERCENT OF TOTAL'
        IF(LXML_STATS)THEN
          WRITE(IXMLLUN,'(A,F10.1,A,//,A,F6.2,A)')'<ioztime unit="seconds">',&
           &ZT_SUMIO, '</ioztime>',&
           &'<fracioztime unit="percent">',ZT_SUMIO/TIMESUM(0)*100.0_JPRB,&
           &'</fracioztime>'
        ENDIF
      ENDIF
      IF(ZT_SUM3 > 0.0_JPRB) THEN
        WRITE(KULOUT,'(A,F10.1,A,F6.2,A)')'SUMMED TIME IN SERIAL SECTIONS  = '&
        & ,ZT_SUM3, ' SECONDS ',ZT_SUM3/TIMESUM(0)*100.0_JPRB,&
         &' PERCENT OF TOTAL'
        IF(LXML_STATS)THEN
          WRITE(IXMLLUN,'(A,F10.1,A,//,A,F6.2,A)')'<serialztime unit="seconds">',&
           & ZT_SUM3,'</serialztime>',&
           &'<fracserialztime unit="percent">',&
           &ZT_SUM3/TIMESUM(0)*100.0_JPRB,&
           &'</fracserialztime>'
        ENDIF
      ENDIF
      IF(ZT_SUM4 > 0.0_JPRB) THEN
        WRITE(KULOUT,'(A,F10.1,A,F6.2,A)')'SUMMED TIME IN BARRIERS         = '&
         & ,ZT_SUM4, ' SECONDS ',ZT_SUM4/TIMESUM(0)*100.0_JPRB,&
         &' PERCENT OF TOTAL'
        IF(LXML_STATS)THEN
          WRITE(IXMLLUN,'(A,F10.1,A,//,A,F6.2,A)')&
           &'<barrierztime unit="seconds">',&
           &ZT_SUM4,'</barrierztime>',&
           & '<fracbarrierztime unit="percent">',&
           & ZT_SUM4/TIMESUM(0)*100.0_JPRB,'</fracbarrierztime>'
        ENDIF
      ENDIF
      IF(ZT_SUM5 > 0.0_JPRB) THEN
        WRITE(KULOUT,'(A,F10.1,A,F6.2,A)')'SUMMED TIME IN MIXED SECTIONS   = '&
         & ,ZT_SUM5, ' SECONDS ',ZT_SUM5/TIMESUM(0)*100.0_JPRB,&
         &' PERCENT OF TOTAL'
        WRITE(IXMLLUN,'(A,F10.1,A,//,A,F6.2,A)')&
         &'<mixedztime unit="seconds">',&
         &ZT_SUM5,'</mixedztime>',&
         & '<fracmixedztime unit="percent">',&
         & ZT_SUM5/TIMESUM(0)*100.0_JPRB,'</fracmixedztime>'
      ENDIF
      IF(LSTATS_COMMS.AND.LSTATS_OMP)THEN
        WRITE(KULOUT,'(A,F8.2)')'FRACTION OF TOTAL TIME ACCOUNTED FOR ',&
         & (ZT_SUM+ZT_SUM2+ZT_SUMIO+ZT_SUM3+ZT_SUM4+ZT_SUM5)/TIMESUM(0)*100.0_JPRB
        WRITE(KULOUT,'(A,F8.2)')'FRACTION OF TOTAL TIME ACCOUNTED FOR INCLUDING SUMB ',&
         & (ZT_SUM+ZT_SUM2+ZT_SUMIO+ZT_SUM3+ZT_SUM4+ZT_SUM5+ZT_SUMB)/TIMESUM(0)*100.0_JPRB
        WRITE(KULOUT,'(" ")')
        IF(LXML_STATS)THEN
          WRITE(IXMLLUN,'(A,F8.2,A)')'<fractotal unit="percent">',&
           &(ZT_SUM+ZT_SUM2+ZT_SUMIO+ZT_SUM3+ZT_SUM4+ZT_SUM5)/TIMESUM(0)*100.0_JPRB,&
           &'</fractotal>'
        ENDIF
      ENDIF
    ENDIF
    IF( LDETAILED_STATS .AND. JROC < 3 ) THEN
      IF(LXML_STATS)THEN
        WRITE(IXMLLUN,'(A)')'</timing>'
      ENDIF
    ENDIF
    IF( LSTATS_MPL .AND. JROC <= NPRNT_STATS ) THEN
      WRITE(KULOUT,'(/,A,I4,/)') 'COMMUNICATIONS STATISTICS:PROCESSOR=',JROC
      WRITE(KULOUT,'(A)') &
       &' NUM     ROUTINE                              '//&
       &'  SUM(s)   SENDS  AVG(kb) TOTAL(MB) MB/s '//&
       &'  RECVS  AVG(kb) TOTAL(MB) MB/s '

      do jnum=501,1000
        if(numsend(jnum) /= 0 .or. numrecv(jnum) /= 0 ) then
           sendrate=sendbytes(jnum)*1.e-6/timesum(jnum)
           recvrate=recvbytes(jnum)*1.e-6/timesum(jnum)
           if(numsend(jnum) /= 0) then
             avgsendlen=sendbytes(jnum)*1.e-3/numsend(jnum)
           else
             avgsendlen=0.0_JPRB
           endif
           if(numrecv(jnum) /= 0) then
             avgrecvlen=recvbytes(jnum)*1.e-3/numrecv(jnum)
           else
             avgsendlen=0.0_JPRB
           endif
           if(numrecv(jnum) /= 0) then
             avgrecvlen=recvbytes(jnum)*1.e-3/numrecv(jnum)
           else
             avgrecvlen=0.0_JPRB
           endif
           write(kulout,'(I6,1X,A40,f6.1,2(I8,3F8.1))') &
            &  jnum,ccdesc(jnum),timesum(jnum),&
            &  numsend(jnum),avgsendlen,sendbytes(jnum)*1.e-6, sendrate, &
            &  numrecv(jnum),avgrecvlen,recvbytes(jnum)*1.e-6, recvrate
        endif
      enddo
      WRITE(KULOUT,'(/,A,I4,/)') 'UNKNOWN COMMUNICATIONS STATISTICS:PROCESSOR=', JROC
      WRITE(KULOUT,'(A)') &
       &' NUM     BEFORE ROUTINE                        '//&
       &'    SENDS TOTAL(MB) '//&
       &'RECVS TOTAL(MB)  '
      do jnum=501,1000
        if(unknown_numsend(jnum) /= 0 .or. unknown_numrecv(jnum) /= 0 ) then
           write(kulout,'(I6,1X,A40,2(I8,F8.1))') &
            &  jnum,ccdesc(jnum),&
            &  unknown_numsend(jnum),unknown_sendbytes(jnum)*1.e-6, &
            &  unknown_numrecv(jnum),unknown_recvbytes(jnum)*1.e-6
        endif
      enddo
      write(kulout,'(7x,"TOTAL",35x,2(I8,F8.1),//)') &
       & sum(unknown_numsend(:)),sum(unknown_sendbytes(:))*1.e-6 , &
       & sum(unknown_numrecv(:)),sum(unknown_recvbytes(:))*1.e-6

     ENDIF

     if(lstats_mpl) then
      do jnum=501,1000
         totsendbytes(jnum) = totsendbytes(jnum) + sendbytes(jnum)
         totrecvbytes(jnum) = totrecvbytes(jnum) + recvbytes(jnum)
         if(sendbytes(jnum).gt.0.0_JPRB.or. &
         &  recvbytes(jnum).gt.0.0_JPRB) then
           maxcommtime(jnum)  = MAX(maxcommtime(jnum),timesum(jnum))
         endif
      enddo
     endif

  ENDDO
  IF(LXML_STATS)THEN
    WRITE(IXMLLUN,'(A)')'<timing_all_processors>'
  ENDIF
  WRITE(KULOUT,*) ''
  WRITE(KULOUT,'(A)') 'STATS FOR ALL PROCESSORS'
  WRITE(KULOUT,'(A)') &
  &' NUM ROUTINE                                     CALLS  MEAN(ms)   MAX(ms)   FRAC(%)  UNBAL(%)'
  ZTOTUNBAL = 0.0_JPRB
  DO JNUM=0,500
    IF(NCALLS(JNUM) > 1) THEN
      ICALLS = NCALLS(JNUM)/2
      ZMEAN = ZAVEAVE(JNUM)/NPROC_STATS
      ZMAX  = ZAVEMAX(JNUM)
      ZMEANT = ZSUMTOT(JNUM)/NPROC_STATS
      ZMAXT  = ZSUMMAX(JNUM)
      IF(ZMEANT .NE. 0.0)THEN
        ZUNBAL= (ZMAXT-ZMEANT)/ZTOTAL*100._JPRB
      ELSE
        ZUNBAL=0.0
      ENDIF
      ZFRAC=ZFRACMAX(JNUM)
      WRITE(KULOUT,'(I4,1X,A40,1X,I8,2(1X,F9.1),2(1X,F9.2))')&
       &JNUM,CCDESC(JNUM),ICALLS,ZMEAN,ZMAX,ZFRAC,ZUNBAL

      IF(LXML_STATS)THEN
        WRITE(IXMLLUN,'(A,I4,A,//,A,A40,A,//,A,I8,A,2(A,F9.1,A,//),2(A,F9.2,A,//),A)')&
         &'<item id="',JNUM,'">',&
         &'<description>',CCDESC(JNUM),'</description>',&
         &'<calls>',ICALLS,'</calls>',&
         &'<mean unit="ms">',ZMEAN,'</mean>',&
         &'<max unit="ms">',ZMAX,'</max>',&
         &'<fraction unit="percent">',ZFRAC,'</fraction>',&
         &'<unbalanced unit="percent">',ZUNBAL,'</unbalanced>','</item>'
      ENDIF
    ENDIF
  ENDDO
  IF(LXML_STATS)THEN
    WRITE(IXMLLUN,'(A)')'</timing_all_processors>'
  ENDIF

IF(LSTATS_COMMS)THEN
  IF(LXML_STATS)THEN
    WRITE(IXMLLUN,'(A)')'<timing_communications>'
  ENDIF
  WRITE(KULOUT,*) ''
  WRITE(KULOUT,'(A)') 'STATS FOR COMMUNICATIONS'
  WRITE(KULOUT,'(A)')  &
  &' NUM ROUTINE                                     CALLS    MEAN(ms)  MAX(ms)   FRAC(%)  UNBAL(%)'
  ZT_SUM=0.0
  DO JNUM=500,JPMAXSTAT
    IF((CCTYPE(JNUM).eq."MPL".OR.CCTYPE(JNUM).eq."BAR".OR.CCTYPE(JNUM).eq."GBR".OR.CCTYPE(JNUM).eq."GB2") &
     & .AND.NCALLS(JNUM) > 1) THEN
      ICALLS = NCALLS(JNUM)/2
      ZMEAN = ZAVEAVE(JNUM)/NPROC_STATS
      ZMAX  = ZAVEMAX(JNUM)
      ZMEANT = ZSUMTOT(JNUM)/NPROC_STATS
      ZMAXT  = ZSUMMAX(JNUM)
      IF(ZMEANT .NE. 0.0)THEN
        ZUNBAL= (ZMAXT-ZMEANT)/ZTOTAL*100._JPRB
      ELSE
        ZUNBAL=0.0
      ENDIF
      ZFRAC=ZFRACMAX(JNUM)
      ZTOTUNBAL = ZTOTUNBAL+(ZMAXT-ZMEANT)
      WRITE(KULOUT,'(I4,1X,A40,1X,I8,2(1X,F9.1),2(1X,F9.2))')&
       &JNUM,CCDESC(JNUM),ICALLS,ZMEAN,ZMAX,ZFRAC,ZUNBAL
      IF(LXML_STATS)THEN
        WRITE(IXMLLUN,'(A,I4,A,//,A,A40,A,//,A,I8,A,//,2(A,F9.1,A,//),2(A,F9.2,A,//),A)')&
         &'<comitem id="',JNUM,'">',&
         &'<description>',CCDESC(JNUM),'</description>',&
         &'<calls>',ICALLS,'</calls>',&
         &'<mean unit="ms">',ZMEAN,'</mean>',&
         &'<max unit="ms">',ZMAX,'</max>',&
         &'<fraction unit="percent">',ZFRAC,'</fraction>',&
         &'<unbalanced unit="percent">',ZUNBAL,'</unbalanced>','</comitem>'
      ENDIF

      ZT_SUM=ZT_SUM+ZMEANT
    ENDIF
  ENDDO

  WRITE(KULOUT,*) ''
  WRITE(KULOUT,'(A,F10.1,A)')'SUMMED TIME IN COMMUNICATIONS   = ',ZT_SUM, ' SECONDS '
  if(lstats_mpl) then
    WRITE(KULOUT,'(/,A,/)') 'TOTAL COMMUNICATIONS VOLUMES AND BANDWIDTH'
    WRITE(KULOUT,'(A)') &
     &'   NUM   ROUTINE                              '//&
     &'  SUM(s)  SEND(GB)  RECV(GB)    GB/s'
    do jnum=501,1000
       if(totsendbytes(jnum).gt.0.0_JPRB.or.totrecvbytes(jnum).gt.0.0_JPRB) then
          write(kulout,'(I6,1X,A40,f6.1,2F10.1,F8.1)') &
           &  jnum,ccdesc(jnum),maxcommtime(jnum),&
           &  totsendbytes(jnum)*1.e-9, &
           &  totrecvbytes(jnum)*1.e-9 , &
           &  (totsendbytes(jnum)*1.e-9)/maxcommtime(jnum)
       endif
    enddo
    write(kulout,'(/,A,42x,f6.1,2F10.1,F8.1)') &
     &   'TOTAL', &
     &   SUM(maxcommtime) , &
     &   SUM(totsendbytes)*1.e-9, &
     &   SUM(totrecvbytes)*1.e-9, &
     &   (SUM(totsendbytes)*1.e-9)/SUM(maxcommtime)
  endif

  IF(LXML_STATS)THEN
    WRITE(IXMLLUN,'(A,F10.1,A)')'<zcom unit="seconds">',ZT_SUM, '</zcom>'
  ENDIF
  WRITE(KULOUT,*) ''
  IF(LXML_STATS)THEN
    WRITE(IXMLLUN,'(A)')'</timing_communications>'
  ENDIF

ENDIF
IF(LSTATS_OMP)THEN
  IF(LXML_STATS)THEN
    WRITE(IXMLLUN,'(A)')'<timing_parallel>'
  ENDIF
  WRITE(KULOUT,*) ''
  WRITE(KULOUT,'(A)') 'STATS FOR PARALLEL REGIONS'
  WRITE(KULOUT,'(A)')  &
  &' NUM ROUTINE                                     CALLS    MEAN(ms)  MAX(ms)   FRAC(%)  UNBAL(%)'
  ZT_SUM=0.0
  DO JNUM=500,JPMAXSTAT
    IF(CCTYPE(JNUM).eq."OMP".AND.NCALLS(JNUM) > 1) THEN
      ICALLS = NCALLS(JNUM)/2
      ZMEAN = ZAVEAVE(JNUM)/NPROC_STATS
      ZMAX  = ZAVEMAX(JNUM)
      ZMEANT = ZSUMTOT(JNUM)/NPROC_STATS
      ZMAXT  = ZSUMMAX(JNUM)
      IF(ZMEANT .NE. 0.0)THEN
        ZUNBAL= (ZMAXT-ZMEANT)/ZTOTAL*100._JPRB
      ELSE
        ZUNBAL=0.0
      ENDIF
      ZFRAC=ZFRACMAX(JNUM)
      ZTOTUNBAL = ZTOTUNBAL+(ZMAXT-ZMEANT)
      WRITE(KULOUT,'(I4,1X,A40,1X,I8,2(1X,F9.1),2(1X,F9.2))')&
       &JNUM,CCDESC(JNUM),ICALLS,ZMEAN,ZMAX,ZFRAC,ZUNBAL
      IF(LXML_STATS)THEN
        WRITE(IXMLLUN,'(A,I4,A,//,A,A40,A,//,A,I8,A,//,2(A,F9.1,A,//),2(A,F9.2,A,//),A)')&
         &'<parallelitem id="',JNUM,'">',&
         &'<description>',CCDESC(JNUM),'</description>',&
         &'<calls>',ICALLS,'</calls>',&
         &'<mean unit="seconds">',ZMEAN,'</mean>',&
         &'<max unit="seconds">',ZMAX,'</max>',&
         &'<fraction unit="percent">',ZFRAC,'</fraction>',&
         &'<unbalanced unit="percent">',ZUNBAL,'</unbalanced>',&
         &'</parallelitem>'
      ENDIF

      ZT_SUM=ZT_SUM+ZMEANT
    ENDIF
  ENDDO

  WRITE(KULOUT,*) ''
  WRITE(KULOUT,'(A,F10.1,A)')'SUMMED TIME IN PARALLEL REGIONS = ',ZT_SUM, ' SECONDS '
  WRITE(KULOUT,*) ''

  IF(LXML_STATS)THEN
    WRITE(IXMLLUN,'(A,F10.1,A)')'<zpar unit="seconds">',ZT_SUM, '</zpar>'
  ENDIF

  WRITE(KULOUT,*) ''
  WRITE(KULOUT,'(A)') 'STATS FOR I/O REGIONS'
  WRITE(KULOUT,'(A)')  &
  &' NUM ROUTINE                                     CALLS    MEAN(ms)  MAX(ms)   FRAC(%)  UNBAL(%)'
  ZT_SUM=0.0
  DO JNUM=500,JPMAXSTAT
    IF(CCTYPE(JNUM).eq."IO-".AND.NCALLS(JNUM) > 1) THEN
      ICALLS = NCALLS(JNUM)/2
      ZMEAN = ZAVEAVE(JNUM)/NPROC_STATS
      ZMAX  = ZAVEMAX(JNUM)
      ZMEANT = ZSUMTOT(JNUM)/NPROC_STATS
      ZMAXT  = ZSUMMAX(JNUM)
      IF(ZMEANT .NE. 0.0)THEN
        ZUNBAL= (ZMAXT-ZMEANT)/ZTOTAL*100._JPRB
      ELSE
        ZUNBAL=0.0
      ENDIF
      ZFRAC=ZFRACMAX(JNUM)
      ZTOTUNBAL = ZTOTUNBAL+(ZMAXT-ZMEANT)
      WRITE(KULOUT,'(I4,1X,A40,1X,I8,2(1X,F9.1),2(1X,F9.2))')&
       &JNUM,CCDESC(JNUM),ICALLS,ZMEAN,ZMAX,ZFRAC,ZUNBAL
      IF(LXML_STATS)THEN
        WRITE(IXMLLUN,'(A,I4,A,//,A,A40,A,//,A,I8,A,//,2(A,F9.1,A,//),2(A,F9.2,A,//),A)')&
         &'<para_io_item id="',JNUM,'">',&
         &'<description>',CCDESC(JNUM),'</description>',&
         &'<calls>',ICALLS,'</calls>',&
         &'<mean unit="seconds">',ZMEAN,'</mean>','<max unit="seconds">',&
         & ZMAX,'</max>',&
         &'<fraction unit="percent">',ZFRAC,'</fraction>',&
         &'<unbalanced unit="percent">',ZUNBAL,'</unbalanced>',&
         &'</para_io_item>'
      ENDIF

      ZT_SUM=ZT_SUM+ZMEANT
    ENDIF
  ENDDO

  WRITE(KULOUT,*) ''
  WRITE(KULOUT,'(A,F10.1,A)')'SUMMED TIME IN I/O REGIONS = ',&
       &ZT_SUM, ' SECONDS '
  WRITE(KULOUT,*) ''

  IF(LXML_STATS)THEN
    WRITE(IXMLLUN,'(A,F10.1,A)')'<zio unit="seconds">',ZT_SUM,'</zio>'
  ENDIF

  WRITE(KULOUT,*) ''
  WRITE(KULOUT,'(A)') 'STATS FOR SERIAL(no OMP) REGIONS'
  WRITE(KULOUT,'(A)')  &
  &' NUM ROUTINE                                     CALLS    MEAN(ms)  MAX(ms)   FRAC(%)  UNBAL(%)'
  ZT_SUM=0.0
  DO JNUM=500,JPMAXSTAT
    IF(CCTYPE(JNUM).eq."SER".AND.NCALLS(JNUM) > 1) THEN
      ICALLS = NCALLS(JNUM)/2
      ZMEAN = ZAVEAVE(JNUM)/NPROC_STATS
      ZMAX  = ZAVEMAX(JNUM)
      ZMEANT = ZSUMTOT(JNUM)/NPROC_STATS
      ZMAXT  = ZSUMMAX(JNUM)
      IF(ZMEANT .NE. 0.0)THEN
        ZUNBAL= (ZMAXT-ZMEANT)/ZTOTAL*100._JPRB
      ELSE
        ZUNBAL=0.0
      ENDIF
      ZFRAC=ZFRACMAX(JNUM)
      ZTOTUNBAL = ZTOTUNBAL+(ZMAXT-ZMEANT)
      WRITE(KULOUT,'(I4,1X,A40,1X,I8,2(1X,F9.1),2(1X,F9.2))')&
       &JNUM,CCDESC(JNUM),ICALLS,ZMEAN,ZMAX,ZFRAC,ZUNBAL
      IF(LXML_STATS)THEN
        WRITE(IXMLLUN,'(A,I4,A,A,A40,A,A,I8,A,2(A,F9.1,A),2(A,F9.2,A,//),A)')&
         &'<serialitem id="',JNUM,'">',&
         &'<description>',CCDESC(JNUM),'</description>',&
         &'<calls>',ICALLS,'</calls>',&
         &'<mean unit="ms">',ZMEAN,'</mean>','<max unit="ms">',ZMAX,'</max>',&
         &'<fraction unit="percent">',ZFRAC,'</fraction>',&
         &'<unbalanced unit="percent">',ZUNBAL,'</unbalanced>','</serialitem>'
      ENDIF

      ZT_SUM=ZT_SUM+ZMEANT
    ENDIF
  ENDDO

  WRITE(KULOUT,*) ''
  WRITE(KULOUT,'(A,F10.1,A)')'SUMMED TIME IN SERIAL REGIONS = ',ZT_SUM, ' SECONDS '
  WRITE(KULOUT,*) ''

  IF(LXML_STATS)THEN
    WRITE(IXMLLUN,'(A,F10.1,A)')'<zserial unit="seconds">',&
     &ZT_SUM, '</zserial>'
  ENDIF

  WRITE(KULOUT,*) ''
  WRITE(KULOUT,'(A)') 'STATS FOR MIXED SECTIONS'
  WRITE(KULOUT,'(A)')  &
  &' NUM ROUTINE                                     CALLS    MEAN(ms)  MAX(ms)   FRAC(%)  UNBAL(%)'
  ZT_SUM=0.0
  DO JNUM=500,JPMAXSTAT
    IF(CCTYPE(JNUM).eq."MXD".AND.NCALLS(JNUM) > 1) THEN
      ICALLS = NCALLS(JNUM)/2
      ZMEAN = ZAVEAVE(JNUM)/NPROC_STATS
      ZMAX  = ZAVEMAX(JNUM)
      ZMEANT = ZSUMTOT(JNUM)/NPROC_STATS
      ZMAXT  = ZSUMMAX(JNUM)
      IF(ZMEANT .NE. 0.0)THEN
        ZUNBAL= (ZMAXT-ZMEANT)/ZTOTAL*100._JPRB
      ELSE
        ZUNBAL=0.0
      ENDIF
      ZFRAC=ZFRACMAX(JNUM)
      ZTOTUNBAL = ZTOTUNBAL+(ZMAXT-ZMEANT)
      WRITE(KULOUT,'(I4,1X,A40,1X,I8,2(1X,F9.1),2(1X,F9.2))')&
       &JNUM,CCDESC(JNUM),ICALLS,ZMEAN,ZMAX,ZFRAC,ZUNBAL
      IF(LXML_STATS)THEN
        WRITE(IXMLLUN,'(A,I4,A,A,A40,A,A,I8,A,2(A,F9.1,A),2(A,F9.2,A,//),A)')&
         &'<mixeditem id="',JNUM,'">',&
         &'<description>',CCDESC(JNUM),'</description>',&
         &'<calls>',ICALLS,'</calls>',&
         &'<mean unit="ms">',ZMEAN,'</mean>','<max unit="ms">',ZMAX,'</max>',&
         &'<fraction unit="percent">',ZFRAC,'</fraction>',&
         &'<unbalanced unit="percent">',ZUNBAL,'</unbalanced>','</mixeditem>'
      ENDIF

      ZT_SUM=ZT_SUM+ZMEANT
    ENDIF
  ENDDO

  WRITE(KULOUT,*) ''
  WRITE(KULOUT,'(A,F10.1,A)')'SUMMED TIME IN MIXED SECTIONS = ',ZT_SUM, ' SECONDS '
  WRITE(KULOUT,*) ''

  IF(LXML_STATS)THEN
    WRITE(IXMLLUN,'(A,F10.1,A)')'<zmixed unit="seconds">',&
     &ZT_SUM, '</zmixed>'
  ENDIF

  IF(LXML_STATS)THEN
    WRITE(IXMLLUN,'(A)')'</timing_parallel>'
  ENDIF

ENDIF

  WRITE(KULOUT,'(A,F10.1,A,F4.1,A)')&
   &'TOTAL MEASURED IMBALANCE =',ZTOTUNBAL,&
   &' SECONDS, ',ZTOTUNBAL/ZTOTAL*100._JPRB,' PERCENT'
ELSE
  ZTOTAL=TIMESUM(0)
  ZTOTCPU = TTCPUSUM(0)
  ZTOTVCPU = TVCPUSUM(0)
ENDIF

IF ( MYPROC_STATS == 1) THEN
  WRITE(KULOUT,'(3(A,F10.1)/)')'TOTAL WALLCLOCK TIME ',ZTOTAL,&
   &' CPU TIME',ZTOTCPU,' VECTOR TIME ',ZTOTVCPU
  IF(LXML_STATS)THEN
    WRITE(IXMLLUN,'(3(A,F10.1,A,//)/)')'<totalwallclocktime>',ZTOTAL,&
     &'</totalwallclocktime>',&
     &'<cputime>',ZTOTCPU,'</cputime>',&
     & '<vectortime>',ZTOTVCPU,'</vectortime>'
  ENDIF
ENDIF
IF( LDETAILED_STATS )THEN
  ITAG = JPTAGSTAT
  ZDELAY_MAX=0.0_JPRB
  ZMPL(:)=0.0_JPRB
  ZBAR(:)=0.0_JPRB
  ZGBR(:)=0.0_JPRB
  ZGB2(:)=0.0_JPRB
  ZOMP(:)=0.0_JPRB
  ZIO (:)=0.0_JPRB
  ZSER(:)=0.0_JPRB
  ZMXD(:)=0.0_JPRB
  DO JROC=1,NPROC_STATS
    IF( JROC > 1 )THEN
      IF( MYPROC_STATS == JROC )THEN
        ISEND=1
        CALL MPL_SEND(NDELAY_INDEX,KDEST=NPRCIDS_STATS(ISEND), &
         & KTAG=ITAG,CDSTRING='GSTATS_PRINT:')
        IF( NDELAY_INDEX > 0 )THEN
          CALL MPL_SEND(NDELAY_COUNTER(1:NDELAY_INDEX),KDEST=NPRCIDS_STATS(ISEND), &
           & KTAG=ITAG+1,CDSTRING='GSTATS_PRINT:')
          CALL MPL_SEND(TDELAY_VALUE(1:NDELAY_INDEX),KDEST=NPRCIDS_STATS(ISEND), &
           & KTAG=ITAG+2,CDSTRING='GSTATS_PRINT:')
          DO JDELAY=1,NDELAY_INDEX
            CLTEMP((JDELAY-1)*10+1:JDELAY*10)=CDELAY_TIME(JDELAY)
          ENDDO
          CALL MPL_SEND(CLTEMP(1:NDELAY_INDEX*10),KDEST=NPRCIDS_STATS(ISEND), &
           & KTAG=ITAG+3,CDSTRING='GSTATS_PRINT:')
        ENDIF
      ENDIF
      IF( MYPROC_STATS == 1 )THEN
        CALL MPL_RECV(NDELAY_INDEX,KSOURCE=NPRCIDS_STATS(JROC), &
         & KTAG=ITAG,CDSTRING='GSTATS_PRINT:')
        IF( NDELAY_INDEX > 0 )THEN
          CALL MPL_RECV(NDELAY_COUNTER(1:NDELAY_INDEX),KSOURCE=NPRCIDS_STATS(JROC), &
           & KTAG=ITAG+1,CDSTRING='GSTATS_PRINT:')
          CALL MPL_RECV(TDELAY_VALUE(1:NDELAY_INDEX),KSOURCE=NPRCIDS_STATS(JROC), &
           & KTAG=ITAG+2,CDSTRING='GSTATS_PRINT:')
          CALL MPL_RECV(CLTEMP(1:NDELAY_INDEX*10),KSOURCE=NPRCIDS_STATS(JROC), &
            & KTAG=ITAG+3,CDSTRING='GSTATS_PRINT:')
          DO JDELAY=1,NDELAY_INDEX
            CDELAY_TIME(JDELAY)=CLTEMP((JDELAY-1)*10+1:JDELAY*10)
          ENDDO
        ENDIF
      ENDIF
    ENDIF
    IF( MYPROC_STATS == 1 .AND. NDELAY_INDEX > 0 )THEN
      WRITE(KULOUT,'("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")')
      WRITE(KULOUT,'("PROC=",I6," NUMBER OF UNEXPECTED DELAYS=",I4)') JROC,NDELAY_INDEX
      IF( NDELAY_INDEX == JPMAXDELAYS )THEN
        WRITE(KULOUT,'(" NOTE THAT THE MAXIMUM NUMBER OF DELAYS HAS BEEN REACHED =",I6)')JPMAXDELAYS
      ENDIF
      ZDELAY=SUM(TDELAY_VALUE(1:NDELAY_INDEX))
      WRITE(KULOUT,'("TOTAL UNEXPECTED DELAY TIME (SECS) =",F9.1)') ZDELAY
      ZDELAY_MAX=MAX(ZDELAY_MAX,ZDELAY)
      DO JDELAY=1,NDELAY_INDEX
        WRITE(KULOUT,'(A,":",A,":",A,1X,I4,1X,A3,1X,A40,1X,F6.1)')&
         &CDELAY_TIME(JDELAY)(1:2),CDELAY_TIME(JDELAY)(3:4),CDELAY_TIME(JDELAY)(5:6),&
         &NDELAY_COUNTER(JDELAY),CCTYPE(NDELAY_COUNTER(JDELAY)),&
         &CCDESC(NDELAY_COUNTER(JDELAY)),TDELAY_VALUE(JDELAY)
         IF( CCTYPE(NDELAY_COUNTER(JDELAY)) .EQ. 'MPL' ) ZMPL(JROC)=ZMPL(JROC)+TDELAY_VALUE(JDELAY)
         IF( CCTYPE(NDELAY_COUNTER(JDELAY)) .EQ. 'BAR' ) ZBAR(JROC)=ZBAR(JROC)+TDELAY_VALUE(JDELAY)
         IF( CCTYPE(NDELAY_COUNTER(JDELAY)) .EQ. 'GBR' ) ZGBR(JROC)=ZGBR(JROC)+TDELAY_VALUE(JDELAY)
         IF( CCTYPE(NDELAY_COUNTER(JDELAY)) .EQ. 'GB2' ) ZGB2(JROC)=ZGB2(JROC)+TDELAY_VALUE(JDELAY)
         IF( CCTYPE(NDELAY_COUNTER(JDELAY)) .EQ. 'OMP' ) ZOMP(JROC)=ZOMP(JROC)+TDELAY_VALUE(JDELAY)
         IF( CCTYPE(NDELAY_COUNTER(JDELAY)) .EQ. 'IO-' ) ZIO (JROC)=ZIO (JROC)+TDELAY_VALUE(JDELAY)
         IF( CCTYPE(NDELAY_COUNTER(JDELAY)) .EQ. 'SER' ) ZSER(JROC)=ZSER(JROC)+TDELAY_VALUE(JDELAY)
         IF( CCTYPE(NDELAY_COUNTER(JDELAY)) .EQ. 'MXD' ) ZMXD(JROC)=ZMXD(JROC)+TDELAY_VALUE(JDELAY)
      ENDDO
      WRITE(KULOUT,'(" ")')
      WRITE(KULOUT,'("PROC=",I6," UNEXPECTED DELAYS SORTED BY COUNTER")') JROC
      DO JNUM=500,JPMAXSTAT
        IDELAY=0
        ZDELAY=0.0_JPRB
        DO JDELAY=1,NDELAY_INDEX
          IF( NDELAY_COUNTER(JDELAY) == JNUM )THEN
            IDELAY=IDELAY+1
            ZDELAY=ZDELAY+TDELAY_VALUE(JDELAY)
          ENDIF
        ENDDO
        IF( IDELAY /= 0 )THEN
          WRITE(KULOUT,'(I4,1X,A3,1X,A40,1X,I4,3X,F6.1)')&
           &JNUM,CCTYPE(JNUM),CCDESC(JNUM),IDELAY,ZDELAY
        ENDIF
      ENDDO
      WRITE(KULOUT,'(" ")')
      WRITE(KULOUT,'(" ")')
    ENDIF
    CALL MPL_BARRIER(CDSTRING='GSTATS_PRINT')
  ENDDO
  IF( MYPROC_STATS == 1 )THEN
    WRITE(KULOUT,'("MAXIMUM TOTAL UNEXPECTED DELAY TIME (SECS) =",F9.1)') ZDELAY_MAX
    WRITE(KULOUT,'(" ")')
    WRITE(KULOUT,'(" ")')
    WRITE(KULOUT,'("  PROC   ","     MPL   ","     BAR   ","     GBR   ","     GB2   ","     OMP   ",&
     &"     IO-   ","     SER   ","     MXD   ")')
    DO JROC=1,NPROC_STATS
      WRITE(KULOUT,'(I6,8(2X,F9.1))') JROC,ZMPL(JROC),ZBAR(JROC),ZGBR(JROC),ZGB2(JROC),&
      &ZOMP(JROC),ZIO (JROC),ZSER(JROC),ZMXD(JROC)
    ENDDO
    WRITE(KULOUT,'(" ")')
    WRITE(KULOUT,'(" ")')
  ENDIF
ENDIF

!   Trace stats

IF (LTRACE_STATS) THEN
  WRITE(KULOUT,'(A)') '=== TRACE OF CALLS TO GSTATS'
  IF (NCALLS_TOTAL > NTRACE_STATS) THEN
    WRITE(KULOUT,'(A,2I10)') ' ONLY PART OF TRACE STORED AS BUFFER TO SMALL ',&
     & NCALLS_TOTAL,NTRACE_STATS
  ENDIF
  WRITE(KULOUT,'(A)') '==='
  CLACTION(0)='ON'
  CLACTION(1)='OFF'
  CLACTION(2)='SUSPEND'
  CLACTION(3)='RESUME'
  DO JCALL=1,MIN(NCALLS_TOTAL,NTRACE_STATS)
    ICALLER = MOD(NCALL_TRACE(JCALL),(JPMAXSTAT+1))
    IACTION   = NCALL_TRACE(JCALL)/(JPMAXSTAT+1)
    IF (IACTION == 0) THEN
      ZTIMELCALL(ICALLER) = TIME_TRACE(JCALL)
      ZTHISTIME(ICALLER) = 0.0_JPRB
    ELSEIF (IACTION == 2) THEN
      ZTHISTIME(ICALLER) = TIME_TRACE(JCALL)-ZTIMELCALL(ICALLER)
    ELSEIF (IACTION == 3) THEN
      ZTIMELCALL(ICALLER) = TIME_TRACE(JCALL)
    ENDIF
    IF (IACTION == 1) THEN
      WRITE(KULOUT,'(1X,F10.3,1X,A,1X,A,1X,F10.3)') &
       &TIME_TRACE(JCALL),CCDESC(ICALLER),CLACTION(IACTION),&
       &ZTHISTIME(ICALLER)+(TIME_TRACE(JCALL)-ZTIMELCALL(ICALLER))
    ELSE
      WRITE(KULOUT,'(1X,F10.3,1X,A,1X,A)') TIME_TRACE(JCALL),CCDESC(ICALLER),&
       & CLACTION(IACTION)
    ENDIF
  ENDDO
ENDIF

IF(LSTATS .AND. MYPROC_STATS == 1) THEN
  PAVEAVE(0:KLEN) = ZAVEAVE(0:KLEN)
ELSE
  PAVEAVE(0:KLEN) = 0.0_JPRB
ENDIF

WRITE(KULOUT,'(/A)')'===-=== END   OF TIMING STATISTICS ===-==='


IF(LSTATS_MEM)THEN
  IF(LXML_STATS)THEN
    WRITE(IXMLLUN,'(A)')'<memory>'
  ENDIF
  WRITE(KULOUT,*) ''
  WRITE(KULOUT,*) 'STATS FOR MEMORY'
  WRITE(KULOUT,*)  &
 &' NUM ROUTINE               CALLS  CALL    MAXINCR   TOTINCR   MININCR'
  WRITE(KULOUT,*)  &
 &'                                   NO      (KB)      (KB)      (KB)'
  DO JNUM=0,JPMAXSTAT
    IF(NCALLS(JNUM) > 1) THEN
      ICALLS = NCALLS(JNUM)/2
      IMEM=NTMEM(JNUM,1)
      INUM=NTMEM(JNUM,3)/2
      JMEM=NTMEM(JNUM,4)
      WRITE(KULOUT,'(I4,1X,A20,1X,I8,1X,I6,3(1X,I9))')&
       &JNUM,CCDESC(JNUM),ICALLS,INUM,IMEM,JMEM,NTMEM(JNUM,5)

      IF(LXML_STATS)THEN
        WRITE(IXMLLUN,'(A,I4,A,//,A,A20,A,//,A,I8,A,//,A,I6,A,//,3(A,I9,A,//))')&
         &'<memitem id="',JNUM,'"/>',&
         &'<description>',CCDESC(JNUM),'</description>',&
         &'<calls>',ICALLS,'</calls>',&
         &'<callnum>',INUM,'</callnum>','<maxincr unit="kb">',IMEM,'</maxincr>',&
         &'<totincr unit="kb">',JMEM,'</totincr>',&
         &'<minincr unit="kb">',NTMEM(JNUM,5),'</minincr>'
      ENDIF
    ENDIF
  ENDDO

  WRITE(KULOUT,*) ''
  WRITE(KULOUT,'(/A)')'===-=== END   OF MEMORY STATISTICS ===-==='
  WRITE(KULOUT,*) ''
  IF(LXML_STATS)THEN
    WRITE(IXMLLUN,'(A)')'</memory>'
  ENDIF
ENDIF
IF(LXML_STATS)THEN
  WRITE(IXMLLUN,'(A)')'</gstats>'
  CLOSE(IXMLLUN)
ENDIF

RETURN
END SUBROUTINE GSTATS_PRINT
