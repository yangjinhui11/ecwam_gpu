interface
SUBROUTINE WGRIBENCODE_MODEL (IU06, ITEST, I1, I2, FIELD, &
 & ITABLE, IPARAM, KLEV, IK, IM, &
 & CDATE, IFCST, MARSTYPE, &
 & IGRIB_HANDLE)
 use parkind_wave, only:&
 & jwim,&
 & jwrb
 INTEGER(KIND=JWIM), INTENT(IN) :: IU06, ITEST, I1, I2
 INTEGER(KIND=JWIM), INTENT(IN) :: ITABLE, IPARAM, KLEV, IK, IM, IFCST
 INTEGER(KIND=JWIM), INTENT(OUT) :: IGRIB_HANDLE
 CHARACTER(LEN=2), INTENT(IN) :: MARSTYPE
 CHARACTER(LEN=14), INTENT(IN) :: CDATE
 REAL(KIND=JWRB), INTENT(INOUT) :: FIELD(I1,I2)
END SUBROUTINE WGRIBENCODE_MODEL
end interface
