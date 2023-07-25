interface
SUBROUTINE MICEP (IPARAM, KIJS, KIJL, IFROMIJ, JFROMIJ,    &
 &                NXS, NXE, NYS, NYE, FIELDG,              &
 &                CICVR, CITH, NEMOCICOVER, NEMOCITHICK)
 use parkind_wave, only:&
 & jwim,&
 & jwrb,&
 & jwro
 USE YOWDRVTYPE , ONLY : FORCING_FIELDS
INTEGER(KIND=JWIM), INTENT(IN) :: IPARAM, KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: IFROMIJ  ,JFROMIJ
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(IN) :: FIELDG
      REAL(KIND=JWRB), DIMENSION (KIJS:KIJL), INTENT(INOUT) :: CICVR, CITH
      REAL(KIND=JWRO), DIMENSION (KIJS:KIJL), INTENT(IN) :: NEMOCICOVER, NEMOCITHICK
END SUBROUTINE MICEP
end interface
