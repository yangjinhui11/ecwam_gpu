interface
 SUBROUTINE MPDISTRIBFL(ISEND, ITAG, NBLKS, NBLKE, KINF, KSUP, &
 & MINF, MSUP, FL)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 use yowmpp , only:&
 & nproc
 use yowparam , only:&
 & niblo
 INTEGER(KIND=JWIM), INTENT(IN) :: ISEND, ITAG, KINF, KSUP, MINF, MSUP
 INTEGER(KIND=JWIM), DIMENSION(NPROC), INTENT(IN) :: NBLKS, NBLKE
 REAL(KIND=JWRB), DIMENSION(NIBLO,KINF:KSUP,MINF:MSUP), INTENT(INOUT) :: FL
 END SUBROUTINE MPDISTRIBFL
end interface
