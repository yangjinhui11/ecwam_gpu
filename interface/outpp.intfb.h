interface
 SUBROUTINE OUTPP (CDATE, IUOUT, ID1, ID2, NGX, NGY, TITL, CONST, &
 & ARRAY, AMOWEP, AMOSOP, AMOEAP, AMONOP)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 INTEGER(KIND=JWIM), INTENT(IN) :: IUOUT, ID1, ID2, NGX, NGY
 REAL(KIND=JWRB), INTENT(IN) :: CONST
 REAL(KIND=JWRB), INTENT(IN) :: AMOWEP, AMOSOP, AMOEAP, AMONOP
 REAL(KIND=JWRB), DIMENSION(ID1,ID2), INTENT(IN) :: ARRAY
 CHARACTER(LEN=14), INTENT(IN) :: CDATE
 CHARACTER(LEN=100), INTENT(IN) :: TITL
 END SUBROUTINE OUTPP
end interface
