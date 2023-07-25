interface
 SUBROUTINE SAVSPEC(FL1, NBLKS, NBLKE, CDTPRO, CDATEF, CDATER)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 use yowgrid , only:&
 & nproma_wam,&
 & nchnk
 use yowmpp , only:&
 & nproc
 use yowparam , only:&
 & nang,&
 & nfre
 REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(IN) :: FL1
 INTEGER(KIND=JWIM), DIMENSION(NPROC), INTENT(IN) :: NBLKS, NBLKE
 CHARACTER(LEN=14), INTENT(IN) :: CDTPRO, CDATEF, CDATER
 END SUBROUTINE SAVSPEC
end interface
