interface
 SUBROUTINE GSFILE (IU06, IUNIT, ICH, CDATE, CDATEF, FILEID, OPT)
 use parkind_wave, only:&
 & jwim
 INTEGER(KIND=JWIM), INTENT(IN) :: IU06, ICH
 INTEGER(KIND=JWIM), INTENT(INOUT) :: IUNIT
 CHARACTER(LEN=1), INTENT(IN) :: OPT
 CHARACTER(LEN=3), INTENT(IN) :: FILEID
 CHARACTER(LEN=14), INTENT(IN) :: CDATE, CDATEF
 END SUBROUTINE GSFILE
end interface