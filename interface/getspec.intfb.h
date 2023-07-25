interface
SUBROUTINE GETSPEC(FL1, BLK2GLO, BLK2LOC, WVENVI, NBLKS, NBLKE, IREAD)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 use yowdrvtype , only:&
 & wvgridglo,&
 & wvgridloc,&
 & environment
 use yowgrid , only:&
 & nproma_wam,&
 & nchnk
 use yowmpp , only:&
 & nproc
 use yowparam , only:&
 & nang,&
 & nfre,&
 & niblo
    REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(OUT) :: FL1
      TYPE(WVGRIDGLO), INTENT(IN)                                            :: BLK2GLO
      TYPE(WVGRIDLOC), INTENT(IN)                                            :: BLK2LOC
      TYPE(ENVIRONMENT), INTENT(IN)                                          :: WVENVI
      INTEGER(KIND=JWIM), DIMENSION(NPROC), INTENT(IN)                       :: NBLKS, NBLKE
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
END SUBROUTINE GETSPEC
end interface
