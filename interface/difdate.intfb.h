interface
 SUBROUTINE DIFDATE (CDATE1, CDATE2, KSHIFT)
 use parkind_wave, only:&
 & jwim
 INTEGER(KIND=JWIM) :: M, IL, ISI, IDUM, IRET, ISHIFT, KSHIFT
 CHARACTER(LEN=*) :: CDATE1, CDATE2
 END SUBROUTINE DIFDATE
end interface
