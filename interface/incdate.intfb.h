interface
 SUBROUTINE INCDATE (CDATE, ISHIFT)
 use parkind_wave, only:&
 & jwim
 INTEGER(KIND=JWIM), INTENT(IN) :: ISHIFT
 CHARACTER(LEN=*), INTENT(INOUT) :: CDATE
 END SUBROUTINE INCDATE
end interface
