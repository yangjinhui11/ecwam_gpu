interface
 SUBROUTINE MPDISTRIBSCFLD(ISEND, ITAG, NBLKS, NBLKE, FIELD)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 use yowmpp , only:&
 & nproc
 use yowparam , only:&
 & niblo
 INTEGER(KIND=JWIM), INTENT(IN) :: ISEND, ITAG
 INTEGER(KIND=JWIM), DIMENSION(NPROC), INTENT(IN) :: NBLKS, NBLKE
 REAL(KIND=JWRB), DIMENSION(NIBLO), INTENT(INOUT) :: FIELD
 END SUBROUTINE MPDISTRIBSCFLD
end interface
