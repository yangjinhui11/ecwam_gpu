interface
 SUBROUTINE WRITEFL(FL, IJINF, IJSUP, KINF, KSUP, MINF, MSUP, &
 & FILENAME, IUNIT, LOUNIT, LRSTPARAL)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 INTEGER(KIND=JWIM), INTENT(IN) :: IJINF, IJSUP, KINF, KSUP, MINF, MSUP
 INTEGER(KIND=JWIM), INTENT(INOUT) :: IUNIT
 REAL(KIND=JWRB), DIMENSION(IJINF:IJSUP,KINF:KSUP,MINF:MSUP), INTENT(INOUT) :: FL
 CHARACTER(LEN=296), INTENT(IN) :: FILENAME
 LOGICAL, INTENT(IN) :: LOUNIT, LRSTPARAL
 END SUBROUTINE WRITEFL
end interface
